#' @title Cross-sample segmentation
#'
#' @description SCOPE offers a cross-sample Poisson likelihood-based recursive segmentation,
#'  enabling shared breakpoints across cells from the same genetic background.
#'
#' @param Y raw read depth matrix after quality control procedure
#' @param Yhat normalized read depth matrix
#' @param sampname vector of sample names
#' @param ref GRanges object after quality control procedure
#' @param lmax maximum CNV length in number of bins returned
#' @param mode format of returned copy numbers. Only integer mode is supported for scDNA-seq data.
#' @param segment.CODEX2 logical, whether to perform individual segmentation. Default is \code{FALSE}.
#'
#' @return A list with components
#'   \item{poolcall}{Cross-sample CNV callings indicating shared breakpoints}
#'   \item{finalcall}{Final cross-sample segmented callset of CNVs with genotyping results}
#'   \item{finalcall_CODEX2}{Final individual segmented callset of CNVs with genotyping results}
#'   \item{image.orig}{A matrix giving logarithm of normalized z-scores}
#'   \item{image.seg}{A matrix of logarithm of estimated copy number over 2}
#'   \item{iCN}{A matrix of inferred integer copy number profiles}
#'
#' @author Rujin Wang \email{rujin@email.unc.edu}
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges RangesList Views countOverlaps
#' @export
segmentCBScs = function(Y, Yhat, sampname, ref, lmax, mode = "integer", segment.CODEX2 = FALSE) {
  lmax = min(length(ref), lmax)

  poolcall = NULL
  lmax = lmax - 1
  message("Cross-sample segmenting for ", ncol(Y), " samples.")
  Z0 = NULL

  for(sampno in 1:ncol(Y)) {
    message("Calculating scan statistic for sample ", sampno, ": ", sampname[sampno], ".")
    y = Y[, sampno]
    ysum = sum(y)
    yhat = Yhat[, sampno]
    yhatsum = sum(yhat)
    num = length(y)
    y = c(y, rep(0, lmax))
    yhat = c(yhat, rep(0, lmax))
    i = rep(1:num, rep(lmax, num))
    j = rep(1:lmax, num) + i
    yact = rep(0, length(i))
    lambda = rep(0, length(i))
    for (k in 1:num) {
      yact[(lmax * k - (lmax - 1)):(lmax * k)] <- cumsum(y[k:(k + lmax)])[-1]
      lambda[(lmax * k - (lmax - 1)):(lmax * k)] <- cumsum(yhat[k:(k + lmax)])[-1]
    }
    i <- i[j <= num]
    yact <- yact[j <= num]
    lambda <- lambda[j <= num]
    j <- j[j <= num]

    yact[lambda < 20] <- 20
    lambda[lambda < 20] <- 20
    if (mode == "integer") {
      chat <- round(2 * (yact/lambda))
      chat.res <- round(2 * ((ysum - yact)/(yhatsum - lambda)))
      chat.res[which(is.na(chat.res))] <- 2 # when lmax = the length of entire chr: yhatsum - lambda = 0
    }
    lratio <- (1 - chat/2) * lambda + log((chat + 1e-04)/2.0001) * yact
    lratio <- lratio + (1 - chat.res/2) * (yhatsum - lambda) + log((chat.res + 1e-04)/2.0001) * (ysum - yact)
    chat[chat > 14] <- 14
    Z0 = cbind(Z0, lratio)
  }
  Z = rowSums(Z0)

  if (sum(Z > 0) > 0) {
    if (sum(Z > 0) >= 2) {
      finalmat <- (cbind(i, j, Z))[Z > 0, ]
      finalmat <- finalmat[order(-finalmat[, 3]), ]
      s <- 1
      while (s <= (nrow(finalmat))) {
        rowstart <- finalmat[s, 1]
        rowend <- finalmat[s, 2]
        rowsel <- (finalmat[, 1] <= rowend & finalmat[, 2] >= rowstart)
        rowsel[s] <- FALSE
        finalmat <- finalmat[!rowsel, ]
        if (is.vector(finalmat)) {
          finalmat <- t(as.matrix(finalmat))
        }
        s <- s + 1
      }
    }
    if (sum(Z > 0) == 1) {
      finalmat <- (cbind(i, j, Z))[Z > 0, ]
      finalmat <- t(as.matrix(finalmat))
    }
    finalmat <- round(finalmat, digits = 3)
    loglikeij <- rep(NA, nrow(finalmat))
    mBIC <- rep(NA, nrow(finalmat))

    kappa1 = 3/2
    kappa2 = 2.27
    N = ncol(Y)
    T = num

    for (s in 1:nrow(finalmat)){
      tau <- sort(unique(c(as.vector(finalmat[1:s, 1:2]), 1, num)))
      m <- length(tau) - 2
      J = matrix(data = NA, nrow = m+2, ncol = N)

      # deltas
      Y0 = Y
      Yhat0 = Yhat
      Y0[Y0 <= 20] = 20
      Yhat0[Yhat0 <= 20] = 20

      muhat = matrix(data = NA, nrow = m+1, ncol = N)
      rhat = matrix(data = NA, nrow = m+1, ncol = N)
      muhat[1,] = round(apply(Y0[1:tau[2], ,drop = F], 2, sum) / tau[2])
      rhat[1,] = round((apply(Y0[1:tau[2], ,drop = F], 2, sum)/apply(Yhat0[1:tau[2], ,drop = F], 2, sum)))
      for(r in 1:m){
        muhat[r+1,] = round(apply(Yhat0[(tau[r+1]+1):tau[r+2], ,drop = F], 2, sum) / (tau[r+2] - tau[r+1] + 1))
        rhat[r+1,] = round((apply(Y0[(tau[r+1]+1):tau[r+2], ,drop = F], 2, sum)/apply(Yhat0[(tau[r+1]+1):tau[r+2], ,drop = F], 2, sum)))
      }
      carriershat = rhat[2:(m+1),] - rhat[1:m,]
      carriershat[carriershat!=0] = 1
      J[2:(m+1),] = carriershat
      deltahat = muhat[2:(m+1),] - muhat[1:m,]
      deltahatJ = deltahat * J[2:(m+1),]
      if(is.null(dim(deltahatJ))){
        deltahatJ = matrix(data=deltahat,nrow=1,ncol=length(deltahat))
      }
      deltahatJ.sq.sum = apply(deltahatJ^2,1,sum)
      deltahatJ.sq.sum = max(deltahatJ.sq.sum,1)    # lower bound at 1, so that log(...) would at least be non-negative.

      # M and pi
      M = sum(carriershat)
      pihat = M/(N*m)

      temp = matrix(c(tau[1:(length(tau)-1)], c(tau[2:(length(tau)-1)] - 1, tau[length(tau)])), ncol = 2, byrow = FALSE)
      L = matrix(data = 0, ncol = N, nrow = (m+1))
      for (r in 1:nrow(L)) {
        yact.temp = apply(Y[temp[r,1]:temp[r,2], ,drop = F], 2, sum)
        lambda.temp = apply(Yhat[temp[r,1]:temp[r,2], ,drop = F], 2, sum)
        yact.temp[lambda.temp < 20] <- 20
        lambda.temp[lambda.temp < 20] <- 20
        L[r,] = (1-round(2*yact.temp/lambda.temp)/2)*lambda.temp + log((round(2 * (yact.temp/lambda.temp)) + 1e-04)/2.0001) * yact.temp
      }
      loglikeij[s] = sum(L)

      term1 = loglikeij[s]
      if(M == 0){
        term2 = 0
      } else{
        term2 = - M/2 * log(2 * loglikeij[s] / M)
      }
      term3 = - log(choose(T, m))
      # term3 = - m * log(choose(T, m))
      term4 = - M/2
      term5 = - sum(log(deltahatJ.sq.sum))
      term6 = - m * (kappa1 - kappa2)
      if(pihat == 0 || pihat ==1){
        term7 = 0
      } else {
        term7 = (M * log(pihat) + (N*m - M) * log(1-pihat))
      }

      mbic = term1 + term2 + term3 + term4 + term5 + term6 + term7
      mBIC[s] <- mbic
    }
    mBIC <- round(mBIC, digits = 3)

    if (mBIC[1] > 0) {
      # finalmat <- cbind(rep(sampname_qc[sampno], nrow(finalmat)), rep(chr, nrow(finalmat)), finalmat)
      finalmat <- (cbind(finalmat, mBIC)[1:which.max(mBIC), ,drop = F])
      finalmat.temp = finalmat[order(finalmat[,1]),]
      if(is.null(dim(finalmat.temp))){
        finalmat.temp = matrix(finalmat[order(finalmat[,1]),], nrow = 1, byrow = TRUE)
      }
      if(nrow(finalmat.temp) != 1){
        for (k in 1:(nrow(finalmat.temp)-1)) {
          if(finalmat.temp[k,2] + 1 != finalmat.temp[k+1, 1]){
            finalmat.temp = rbind(finalmat.temp, c(finalmat.temp[k,2] + 1, finalmat.temp[k+1,1] - 1, NA, NA))
          }
        }
        finalmat.temp = finalmat.temp[order(finalmat.temp[,1]),]
      }
      if(finalmat.temp[1,1] != 1){
        finalmat.temp = rbind(c(1, finalmat.temp[1,1] - 1, NA, NA), finalmat.temp)
      }
      if(finalmat.temp[nrow(finalmat.temp),2] != length(ref)){
        finalmat.temp = rbind(finalmat.temp, c(finalmat.temp[nrow(finalmat.temp),2] + 1, length(ref), NA, NA))
      }
      poolcall <- finalmat.temp[,c(1,2)]
      # tau <- sort(unique(c(as.vector(finalmat[1:which.max(mBIC), 1:2]), 1, num)))
      # if(tau[1] == 1 & tau[2] == 2){
      #   poolcall = matrix(c(1, tau[3:(length(tau)-1)], 2, c(tau[4:(length(tau)-1)] - 1, tau[length(tau)])), ncol = 2, byrow = FALSE)
      # } else{
      #   poolcall = matrix(c(tau[1:(length(tau)-1)], c(tau[2:(length(tau)-1)] - 1, tau[length(tau)])), ncol = 2, byrow = FALSE)
      # }
      colnames(poolcall) = c("st", "end")
    }
  }

  if (is.vector(poolcall)) {
    poolcall <- t(as.matrix(poolcall))
  }

  finalcall = NULL
  image.orig=log(pmax(0.001,Y)/Yhat)
  image.orig[image.orig<=-2]=-2
  image.orig[image.orig>2]=2

  image.seg = matrix(data=NA, nrow = nrow(image.orig), ncol = ncol(image.orig))
  iCN = matrix(data=NA, nrow = nrow(image.orig), ncol = ncol(image.orig))
  colnames(image.seg) = colnames(Y)
  colnames(iCN) = colnames(Y)

  for (i in 1:nrow(poolcall)) {
    st_bin = poolcall[i,"st"]
    ed_bin = poolcall[i,"end"]
    st_bp <- start(ref)[as.numeric(poolcall[i,"st"])]
    ed_bp <- end(ref)[as.numeric(poolcall[i,"end"])]

    yact.ind = colSums(Y[st_bin:ed_bin, ,drop = F])
    lambda.ind = colSums(Yhat[st_bin:ed_bin, ,drop = F])
    if (mode == "integer") {
      chat.ind <- round(2 * (yact.ind/lambda.ind))
    } else if (mode == "fraction") {
      chat.ind <- 2 * (yact/lambda)
    }
    chat.ind[chat.ind > 14] = 14

    image.seg[poolcall[i,1]:poolcall[i,2], ] = matrix(nrow = (poolcall[i,2] - poolcall[i,1] + 1),
                                                      ncol = ncol(image.seg), data = log(chat.ind / 2), byrow = TRUE)
    iCN[poolcall[i,1]:poolcall[i,2], ] = matrix(nrow = (poolcall[i,2] - poolcall[i,1] + 1),
                                                ncol = ncol(iCN), data = chat.ind, byrow = TRUE)

    temp = cbind(colnames(Y), rep(st_bin, ncol(Y)), rep(ed_bin, ncol(Y)), rep(st_bp, ncol(Y)),
                 rep(ed_bp, ncol(Y)), chat.ind)
    temp = temp[chat.ind != 2, ]
    finalcall = rbind(finalcall, temp)
    rownames(finalcall) = NULL
    colnames(finalcall) = c("sample_name", "st_bin", "ed_bin", "st_bp", "ed_bp", "cnv_no")
  }

  finalcall = as.data.frame(finalcall)

  finalcall$st_bin = as.numeric(paste(finalcall$st_bin))
  finalcall$ed_bin = as.numeric(paste(finalcall$ed_bin))
  finalcall$st_bp = as.numeric(paste(finalcall$st_bp))
  finalcall$ed_bp = as.numeric(paste(finalcall$ed_bp))
  finalcall$cnv_no = as.numeric(paste(finalcall$cnv_no))


  finalcall_CODEX2 = NULL
  if(segment.CODEX2){
    # CODEX2 segmenting by individuals
    for(sampno in 1:ncol(Y)) {
      message("Segmenting sample ", sampno, ": ", sampname[sampno], ".")
      # y = Y[, sampno]
      # yhat = Yhat[, sampno]
      # num = length(y)
      # y = c(y, rep(0, lmax))
      # yhat = c(yhat, rep(0, lmax))
      # i = rep(1:num, rep(lmax, num))
      # j = rep(1:lmax, num) + i
      # yact = rep(0, length(i))
      # lambda = rep(0, length(i))
      # for (k in 1:num) {
      # yact[(lmax * k - (lmax - 1)):(lmax * k)] <- cumsum(y[k:(k + lmax)])[-1]
      # lambda[(lmax * k - (lmax - 1)):(lmax * k)] <- cumsum(yhat[k:(k + lmax)])[-1]
      # }
      # i <- i[j <= num]
      # yact <- yact[j <= num]
      # lambda <- lambda[j <= num]
      # j <- j[j <= num]
      # yact[lambda < 20] <- 20
      # lambda[lambda < 20] <- 20
      # if (mode == "integer") {
      # chat <- round(2 * (yact/lambda))
      # } else if (mode == "fraction") {
      # chat <- 2 * (yact/lambda)
      # }
      # lratio <- (1 - chat/2) * lambda + log((chat + 1e-04)/2.0001) * yact
      # chat[chat > 14] <- 14

      y = Y[, sampno]
      ysum = sum(y)
      yhat = Yhat[, sampno]
      yhatsum = sum(yhat)
      num = length(y)
      y = c(y, rep(0, lmax))
      yhat = c(yhat, rep(0, lmax))
      i = rep(1:num, rep(lmax, num))
      j = rep(1:lmax, num) + i
      yact = rep(0, length(i))
      lambda = rep(0, length(i))
      for (k in 1:num) {
        yact[(lmax * k - (lmax - 1)):(lmax * k)] <- cumsum(y[k:(k + lmax)])[-1]
        lambda[(lmax * k - (lmax - 1)):(lmax * k)] <- cumsum(yhat[k:(k + lmax)])[-1]
      }
      i <- i[j <= num]
      yact <- yact[j <= num]
      lambda <- lambda[j <= num]
      j <- j[j <= num]

      yact[lambda < 20] <- 20
      lambda[lambda < 20] <- 20
      if (mode == "integer") {
        chat <- round(2 * (yact/lambda))
        chat.res <- round(2 * ((ysum - yact)/(yhatsum - lambda)))
        chat.res[which(is.na(chat.res))] <- 2 # when lmax = the length of entire chr: yhatsum - lambda = 0
      }
      lratio <- (1 - chat/2) * lambda + log((chat + 1e-04)/2.0001) * yact
      lratio <- lratio + (1 - chat.res/2) * (yhatsum - lambda) + log((chat.res + 1e-04)/2.0001) * (ysum - yact)
      chat[chat > 14] <- 14

      if (sum(lratio > 0) > 0) {
        if (sum(lratio > 0) >= 2) {
          finalmat <- (cbind(i, j, yact, lambda, chat, lratio))[lratio > 0, ]
          finalmat <- finalmat[order(-finalmat[, 6]), ]
          s <- 1
          while (s <= (nrow(finalmat))) {
            rowstart <- finalmat[s, 1]
            rowend <- finalmat[s, 2]
            rowsel <- (finalmat[, 1] <= rowend & finalmat[, 2] >= rowstart)
            rowsel[s] <- FALSE
            finalmat <- finalmat[!rowsel, ]
            if (is.vector(finalmat)) {
              finalmat <- t(as.matrix(finalmat))
            }
            s <- s + 1
          }
        }
        if (sum(lratio > 0) == 1) {
          finalmat <- (cbind(i, j, yact, lambda, chat, lratio))[lratio > 0, ]
          finalmat <- t(as.matrix(finalmat))
        }
        finalmat <- round(finalmat, digits = 3)
        loglikeij <- cumsum(finalmat[, 6])
        mBIC <- rep(NA, length(loglikeij))
        for (s in 1:nrow(finalmat)) {
          tau <- sort(unique(c(as.vector(finalmat[1:s, 1:2]), 1, num)))
          P <- length(tau) - 2
          mbic <- loglikeij[s]
          mbic <- mbic - 0.5 * sum(log(tau[2:length(tau)] - tau[1:(length(tau) - 1)]))
          mbic <- mbic + (0.5 - P) * log(num)
          mBIC[s] <- mbic
        }
        mBIC <- round(mBIC, digits = 3)
        if (mBIC[1] > 0) {
          finalmat <- cbind(rep(sampname[sampno], nrow(finalmat)), finalmat)
          finalmat <- (cbind(finalmat, mBIC)[1:which.max(mBIC), ])
          finalcall_CODEX2 <- rbind(finalcall_CODEX2, finalmat)
        }
      }
    }

    if (is.vector(finalcall_CODEX2)) {
      finalcall_CODEX2 <- t(as.matrix(finalcall_CODEX2))
    }

    colnames(finalcall_CODEX2) <- c("sample_name", "st_bin", "ed_bin", "raw_cov",
                                    "norm_cov", "copy_no", "lratio", "mBIC")
    rownames(finalcall_CODEX2) <- rep("", nrow(finalcall_CODEX2))
    lratio = as.numeric(finalcall_CODEX2[, "lratio"])
    rownames(finalcall_CODEX2) = paste("cnv", 1:nrow(finalcall_CODEX2), sep = "")
    finalcall_CODEX2 = as.data.frame(finalcall_CODEX2)
    finalcall_CODEX2$copy_no = as.numeric(paste(finalcall_CODEX2$copy_no))
    finalcall_CODEX2$raw_cov = as.numeric(paste(finalcall_CODEX2$raw_cov))
    finalcall_CODEX2$norm_cov = as.numeric(paste(finalcall_CODEX2$norm_cov))
    finalcall_CODEX2$st_bin = as.numeric(paste(finalcall_CODEX2$st_bin))
    finalcall_CODEX2$ed_bin = as.numeric(paste(finalcall_CODEX2$ed_bin))
    finalcall_CODEX2$lratio = as.numeric(paste(finalcall_CODEX2$lratio))
    finalcall_CODEX2$mBIC = as.numeric(paste(finalcall_CODEX2$mBIC))
  }

  return(list(poolcall = poolcall, finalcall = finalcall, finalcall_CODEX2 = finalcall_CODEX2, image.orig = image.orig, image.seg = image.seg, iCN = iCN))
}
