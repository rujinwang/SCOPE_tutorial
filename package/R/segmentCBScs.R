#' @title Cross-sample segmentation
#'
#' @description SCOPE offers a cross-sample Poisson likelihood-based recursive segmentation,
#'  enabling shared breakpoints across cells from the same genetic background.
#'
#' @param Y raw read depth matrix after quality control procedure
#' @param Yhat normalized read depth matrix
#' @param sampname vector of sample names
#' @param ref GRanges object after quality control procedure
#' @param mode format of returned copy numbers. Only integer mode is supported for scDNA-seq data.
#' @param segment.idv logical, whether to perform individual segmentation. Default is \code{FALSE}.
#' @param max.ns a number specifying how many rounds of nested structure searching would be performed. Defalut is \code{0}.
#'
#' @return A list with components
#'   \item{poolcall}{Cross-sample CNV callings indicating shared breakpoints}
#'   \item{finalcall}{Final cross-sample segmented callset of CNVs with genotyping results}
#'   \item{finalcall.idv}{Final individual segmented callset of CNVs with genotyping results}
#'   \item{image.orig}{A matrix giving logarithm of normalized z-scores}
#'   \item{image.seg}{A matrix of logarithm of estimated copy number over 2}
#'   \item{iCN}{A matrix of inferred integer copy number profiles}
#'   \item{iCN.idv}{A matrix of inferred integer copy number profiles by individual segmentation}
#'
#' @author Rujin Wang \email{rujin@email.unc.edu}
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges RangesList Views countOverlaps
#' @export
segmentCBScs = function(Y, Yhat, sampname, ref, mode = "integer", segment.idv = FALSE, max.ns = 0) {
  create_chptsmat = function(mat1, st_end){
    st = st_end[1]
    end = st_end[2]
    mat1 = mat1[order(mat1[,1]), , drop = F]
    if(mat1[1,1] != st){
      newchptsmat = matrix(data = c(st, mat1[1,1]-1, mat1[1,]), ncol = 2, byrow = TRUE)
    } else{
      newchptsmat = t(as.matrix(mat1[1,]))
    }
    if(nrow(mat1) > 1){
      for (r in 2:nrow(mat1)) {
        if(mat1[r,1] != mat1[r-1,2] + 1){
          newchptsmat = rbind(newchptsmat, matrix(data = c(mat1[r-1,2] + 1, mat1[r,1] - 1, mat1[r,]), ncol = 2, byrow = TRUE))
        }else{
          newchptsmat = rbind(newchptsmat, matrix(data = mat1[r,], ncol = 2, byrow = TRUE))
        }
      }
    }
    if(mat1[nrow(mat1),2] != end){
      newchptsmat = rbind(newchptsmat, matrix(data = c(mat1[nrow(mat1),2]+1, end), nrow = 1, ncol = 2))
    }
    return(newchptsmat)
  }

  compute_cs_lratio = function(Y, Yhat, sampname, this.chpts, msgprint = TRUE){
    Z0 = NULL
    if(this.chpts[2] - this.chpts[1] < 2){
      return(list(i = NA, j = NA, Z = NA, finalmat = t(as.matrix(this.chpts)), is.nested = 0))
    } else{
      lmax = nrow(Y) - 1
      for(sampno in 1:ncol(Y)) {
        if(msgprint){
          message("Calculating scan statistic for sample ", sampno, ": ", sampname[sampno], ".")
        }

        y = Y[, sampno]
        yhat = Yhat[, sampno]

        if(any(yhat < 10)){
          if(any(y[yhat < 10] < 10)){
            y[yhat < 10] <- 10
          }
          yhat[yhat < 10] <- 10
        }

        ysum = sum(y)
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

        if (mode == "integer") {
          chat <- round(2 * (yact/lambda))
          idx.noLeft = which(i == this.chpts[1])
          idx.noRight = which(j == this.chpts[2])
          chat.L = rep(NA, length(chat))
          chat.R = rep(NA, length(chat))
          yact.L = rep(NA, length(chat))
          yact.R = rep(NA, length(chat))
          lambda.L = rep(NA, length(chat))
          lambda.R = rep(NA, length(chat))
          for (x in 1:length(chat)) {
            if(x %in% idx.noLeft){
              yact.L[x] = 0
              lambda.L[x] = 0
            }else{
              yact.L[x] = sum(y[1:(i[x]-1)])
              lambda.L[x] = sum(yhat[1:(i[x]-1)])
            }

            if(x %in% idx.noRight){
              yact.R[x] = 0
              lambda.R[x] = 0
            }else{
              yact.R[x] = sum(y[(j[x]+1):length(y)])
              lambda.R[x] = sum(yhat[(j[x]+1):length(y)])
            }
          }
          chat.L <- round(2 * (yact.L/lambda.L))
          chat.R <- round(2 * (yact.R/lambda.R))
          chat.L[which(is.na(chat.L))] <- 2 # when lmax = the length of entire chr: yhatsum - lambda = 0
          chat.R[which(is.na(chat.R))] <- 2 # when lmax = the length of entire chr: yhatsum - lambda = 0
        }

        lratio.C = (1 - chat/2) * lambda + log((chat + 1e-04)/2.0001) * yact
        lratio.L = (1 - chat.L/2) * lambda.L + log((chat.L + 1e-04)/2.0001) * yact.L
        lratio.R = (1 - chat.R/2) * lambda.R + log((chat.R + 1e-04)/2.0001) * yact.R
        lratio = lratio.C + lratio.L + lratio.R
        chat[chat > 14] <- 14

        Z0 = cbind(Z0, lratio)
      }
      Z = rowSums(Z0)

      winlag = this.chpts[1] - 1
      i = i + winlag
      j = j + winlag

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
        finalmat <- finalmat[finalmat[,2] - finalmat[,1] > 1, ,drop = F]
        finalmat <- finalmat[finalmat[,3] > 10, , drop = F]
        finalmat <- round(finalmat, digits = 3)

        if(nrow(finalmat) == 0){
          is.nested = 0
        } else if(nrow(finalmat) == 1 & finalmat[1,1] == this.chpts[1] & finalmat[1,2] == this.chpts[2]){
          is.nested = 0
        } else if(nrow(finalmat) > 2 & length(unique(finalmat[,3]))==1){
          is.nested = 0
        } else if(max(finalmat[,2] - finalmat[,1]) <= 5){
          is.nested = 0
        } else if(nrow(finalmat) == 1 & (this.chpts[2] - this.chpts[1]) - (finalmat[1,2] - finalmat[1,1]) == 1){
          is.nested = 0
        } else{
          is.nested = 1
        }
      } else{
        finalmat = NULL
        is.nested = 0
      }

      # Avoid losing the boundaries
      if(!is.null(finalmat) && nrow(finalmat)>0){
        min.st = min(finalmat[,1:2])
        if(this.chpts[1] < min.st){
          idx = which(i == this.chpts[1] & j == (min.st - 1))
          idx2 = which(finalmat[,1] == min.st)
          if(length(idx)!=0){
            Z.temp = Z[idx]
            if(Z.temp > 0){
              finalmat[idx2,1] = this.chpts[1]
            } else{
              finalmat = rbind(c(this.chpts[1], min.st-1, 0), finalmat)
            }
          } else{
            finalmat[idx2,1] = this.chpts[1]
          }
        }
        max.ed = max(finalmat[,1:2])
        if(this.chpts[2] > max.ed){
          idx = which(i ==  (max.ed + 1) & j == this.chpts[2])
          idx2 = which(finalmat[,2] == max.ed)
          if(length(idx)!=0){
            Z.temp = Z[idx]
            if(Z.temp > 0){
              finalmat[idx2,2] = this.chpts[2]
            } else{
              finalmat = rbind(finalmat, c(max.ed+1, this.chpts[2], 0))
            }
          } else{
            finalmat[idx2,2] = this.chpts[2]
          }
        }
      }

      # Avoid only one-point fragment: merge into segment with smaller Z/lratio
      if(!is.null(finalmat) && nrow(finalmat)>1){
        for (r in 1:nrow(finalmat)) {
          singlepoint = finalmat[r,1] - 1
          sg.idx = which(finalmat[,2] == singlepoint - 1)
          if(length(sg.idx)!=0 && r < sg.idx){
            finalmat[sg.idx,2] = singlepoint
          }
        }
        for (r in 1:nrow(finalmat)) {
          singlepoint = finalmat[r,2] + 1
          sg.idx = which(finalmat[,1] == singlepoint + 1)
          if(length(sg.idx)!=0 && r < sg.idx){
            finalmat[sg.idx,1] = singlepoint
          }
        }
      }

      if(!is.null(finalmat) && nrow(finalmat)>1){
        finalmat = create_chptsmat(finalmat[,1:2,drop = F], this.chpts)
      }
      return(list(i = i, j = j, Z = Z, finalmat = finalmat, is.nested = is.nested))
    }
  }

  search_cs_nested = function(Y, Yhat, sampname, chptsmat, ...){
    message("Performing cross-sample nested search... \n")
    nest.list = vector("list", nrow(chptsmat))
    for (r in 1:nrow(chptsmat)) {
      nest.list[[r]] = compute_cs_lratio(Y[chptsmat[r,1]:chptsmat[r,2], , drop = F], Yhat[chptsmat[r,1]:chptsmat[r,2], , drop = F], sampname, chptsmat[r,], msgprint = FALSE)
    }
    return(nest.list)
  }

  compute_lratio = function(this.y, this.yhat, this.chpts){
    if(this.chpts[2] - this.chpts[1] < 2){
      return(list(i = NA, j = NA, yact = NA, lambda = NA, chat = NA, lratio = NA, finalmat = t(as.matrix(this.chpts)), is.nested = 0))
    } else{
      if(any(this.yhat < 10)){
        if(any(this.y[this.yhat < 10] < 10)){
          this.y[this.yhat < 10] <- 10
        }
        this.yhat[this.yhat < 10] <- 10
      }
      lmax = length(this.y) - 1
      num = length(this.y)
      this.y = c(this.y, rep(0, lmax))
      this.yhat = c(this.yhat, rep(0, lmax))
      i = rep(1:num, rep(lmax, num))
      j = rep(1:lmax, num) + i
      this.yact = rep(0, length(i))
      this.lambda = rep(0, length(i))
      for (k in 1:num) {
        this.yact[(lmax * k - (lmax - 1)):(lmax * k)] <- cumsum(this.y[k:(k + lmax)])[-1]
        this.lambda[(lmax * k - (lmax - 1)):(lmax * k)] <- cumsum(this.yhat[k:(k + lmax)])[-1]
      }
      i <- i[j <= num]
      this.yact <- this.yact[j <= num]
      this.lambda <- this.lambda[j <= num]
      j <- j[j <= num]

      if (mode == "integer") {
        chat <- round(2 * (this.yact/this.lambda))
        idx.noLeft = which(i == this.chpts[1])
        idx.noRight = which(j == this.chpts[2])
        chat.L = rep(NA, length(chat))
        chat.R = rep(NA, length(chat))
        this.yact.L = rep(NA, length(chat))
        this.yact.R = rep(NA, length(chat))
        this.lambda.L = rep(NA, length(chat))
        this.lambda.R = rep(NA, length(chat))
        for (x in 1:length(chat)) {
          if(x %in% idx.noLeft){
            this.yact.L[x] = 0
            this.lambda.L[x] = 0
          }else{
            this.yact.L[x] = sum(this.y[1:(i[x]-1)])
            this.lambda.L[x] = sum(this.yhat[1:(i[x]-1)])
          }

          if(x %in% idx.noRight){
            this.yact.R[x] = 0
            this.lambda.R[x] = 0
          }else{
            this.yact.R[x] = sum(this.y[(j[x]+1):length(this.y)])
            this.lambda.R[x] = sum(this.yhat[(j[x]+1):length(this.y)])
          }

        }
        chat.L <- round(2 * (this.yact.L/this.lambda.L))
        chat.R <- round(2 * (this.yact.R/this.lambda.R))
        chat.L[which(is.na(chat.L))] <- 2 # when lmax = the length of entire chr: yhatsum - lambda = 0
        chat.R[which(is.na(chat.R))] <- 2 # when lmax = the length of entire chr: yhatsum - lambda = 0
      }
      lratio.C = (1 - chat/2) * this.lambda + log((chat + 1e-04)/2.0001) * this.yact
      lratio.L = (1 - chat.L/2) * this.lambda.L + log((chat.L + 1e-04)/2.0001) * this.yact.L
      lratio.R = (1 - chat.R/2) * this.lambda.R + log((chat.R + 1e-04)/2.0001) * this.yact.R
      lratio = lratio.C + lratio.L + lratio.R
      chat[chat > 14] <- 14

      winlag = this.chpts[1] - 1
      i = i + winlag
      j = j + winlag

      if (sum(lratio > 0) > 0) {
        if (sum(lratio > 0) >= 2) {
          finalmat <- (cbind(i, j, this.yact, this.lambda, chat, lratio, lratio.C, lratio.L, lratio.R))[lratio > 0, ]
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
          finalmat <- (cbind(i, j, this.yact, this.lambda, chat, lratio, lratio.C, lratio.L, lratio.R))[lratio > 0, ]
          finalmat <- t(as.matrix(finalmat))
        }
        finalmat <- finalmat[finalmat[,2] - finalmat[,1] > 1, ,drop = F]
        finalmat <- finalmat[finalmat[,6] > 10, , drop = F]
        finalmat <- round(finalmat, digits = 3)
        if(nrow(finalmat) == 0){
          is.nested = 0
        } else if(nrow(finalmat) == 1 & finalmat[1,1] == this.chpts[1] & finalmat[1,2] == this.chpts[2]){
          is.nested = 0
        } else if(nrow(finalmat) > 2 & length(unique(finalmat[,6]))==1){
          is.nested = 0
        } else if(nrow(finalmat) == 1 & (this.chpts[2] - this.chpts[1]) - (finalmat[1,2] - finalmat[1,1]) == 1){
          is.nested = 0
        } else{
          is.nested = 1
        }
      }else{
        finalmat = NULL
        is.nested = 0
      }

      # Avoid losing the boundaries
      if(!is.null(finalmat) && nrow(finalmat)>0){
        min.st = min(finalmat[,1:2])
        if(this.chpts[1] < min.st){
          idx = which(i == this.chpts[1] & j == (min.st - 1))
          idx2 = which(finalmat[,1] == min.st)
          if(length(idx)!=0){
            lratio.temp = lratio[idx]
            if(lratio.temp > 0){
              finalmat[idx2,1] = this.chpts[1]
            } else{
              finalmat = rbind(c(this.chpts[1], min.st-1, 0), finalmat)
            }
          } else{
            finalmat[idx2,1] = this.chpts[1]
          }
        }
        max.ed = max(finalmat[,1:2])
        if(this.chpts[2] > max.ed){
          idx = which(i ==  (max.ed + 1) & j == this.chpts[2])
          idx2 = which(finalmat[,2] == max.ed)
          if(length(idx)!=0){
            lratio.temp = lratio[idx]
            if(lratio.temp > 0){
              finalmat[idx2,2] = this.chpts[2]
            } else{
              finalmat = rbind(finalmat, c(max.ed+1, this.chpts[2], 0))
            }
          } else{
            finalmat[idx2,2] = this.chpts[2]
          }
        }
      }

      # Avoid only one-point fragment: merge into segment with smaller Z/lratio
      if(!is.null(finalmat) && nrow(finalmat)>1){
        for (r in 1:nrow(finalmat)) {
          singlepoint = finalmat[r,1] - 1
          sg.idx = which(finalmat[,2] == singlepoint - 1)
          if(length(sg.idx)!=0 && r < sg.idx){
            finalmat[sg.idx,2] = singlepoint
          }
        }
        for (r in 1:nrow(finalmat)) {
          singlepoint = finalmat[r,2] + 1
          sg.idx = which(finalmat[,1] == singlepoint + 1)
          if(length(sg.idx)!=0 && r < sg.idx){
            finalmat[sg.idx,1] = singlepoint
          }
        }
      }

      if(!is.null(finalmat) && nrow(finalmat)>1){
        finalmat = create_chptsmat(finalmat[,1:2,drop = F], this.chpts)
      }
      return(list(i = i, j = j, yact = this.yact, lambda = this.lambda, chat = chat, lratio = lratio, lratio.C = lratio.C, lratio.L = lratio.L, lratio.R = lratio.R, finalmat = finalmat, is.nested = is.nested))
    }
  }


  search_nested = function(this.y, this.yhat, chptsmat){
    message("Performing individual nested search... \n")
    nest.list = vector("list", nrow(chptsmat))
    for (r in 1:nrow(chptsmat)) {
      nest.list[[r]] = compute_lratio(this.y[chptsmat[r,1]:chptsmat[r,2]], this.yhat[chptsmat[r,1]:chptsmat[r,2]], chptsmat[r,])
    }
    return(nest.list)
  }

  poolcall = NULL
  message("Cross-sample segmenting for ", ncol(Y), " samples.")

  chpts0 = c(1, nrow(Y))
  cs.scan = compute_cs_lratio(Y, Yhat, sampname, chpts0)
  i = cs.scan$i
  j = cs.scan$j
  Z = cs.scan$Z
  init.cs.finalmat = cs.scan$finalmat

  if(!is.null(init.cs.finalmat) && nrow(init.cs.finalmat)>0){
    # Further cross-sample nested searching
    chpts = init.cs.finalmat[,1:2]
    if(is.null(dim(chpts))){
      chpts = t(as.matrix(chpts))
    }

    keep_going = 1
    # number of nested searching
    num_ns = 1
    while (!all(keep_going==0) & num_ns <= max.ns) {
      nested.output = search_cs_nested(Y, Yhat, sampname, chpts)
      keep_going = sapply(nested.output, function(z) z$is.nested)

      newchpts = NULL
      for (r in 1:length(keep_going)) {
        if(keep_going[r]){
          newchpts = rbind(newchpts, nested.output[[r]]$finalmat[,1:2])
        } else{
          newchpts = rbind(newchpts, chpts[r,])
        }
      }
      chpts = newchpts
      num_ns = num_ns + 1
    }

    # backward compute Z after nested searching
    temp = NULL
    for(r in 1:nrow(chpts)){
      idx = which(i == chpts[r,1] & j == chpts[r,2])
      if(length(idx)!=0){
        temp = c(temp, Z[idx])
      } else{
        temp = c(temp, rep(NA, 4))
      }
    }
    cs.finalmat = cbind(chpts, temp)
    cs.finalmat = cs.finalmat[!is.na(cs.finalmat[,3]), , drop = F]
    cs.finalmat = cs.finalmat[order(-cs.finalmat[,3]), , drop = F]
  } else{ # if init.cs.cs.finalmat = NULL, the whole chromosome is normal
    cs.finalmat = matrix(data = c(chpts0, 0), nrow = 1, ncol = 3, byrow = TRUE)
  }

  cs.loglikeij <- rep(NA, nrow(cs.finalmat))
  cs.mBIC <- rep(NA, nrow(cs.finalmat))

  kappa1 = 3/2
  kappa2 = 2.27
  N = ncol(Y)
  T = nrow(Y)

  for (s in 1:nrow(cs.finalmat)){
    tau <- sort(unique(c(as.vector(cs.finalmat[1:s, 1:2]), 1, nrow(Y))))
    m <- length(tau) - 2
    if(m > 0){
      J = matrix(data = NA, nrow = m+2, ncol = N)

      # deltas
      Y0 = Y
      Yhat0 = Yhat
      Y0[Y0 <= 20] = 20
      Yhat0[Yhat0 <= 20] = 20

      muhat = matrix(data = NA, nrow = m+1, ncol = N)
      rhat = matrix(data = NA, nrow = m+1, ncol = N)
      muhat[1,] = round(apply(Y0[1:tau[2], ,drop = F], 2, sum) / tau[2])
      rhat[1,] = round(2*(apply(Y0[1:tau[2], ,drop = F], 2, sum)/apply(Yhat0[1:tau[2], ,drop = F], 2, sum)))
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

      cs.temp = create_chptsmat(cs.finalmat[1:s,1:2,drop = F], chpts0)
      L = matrix(data = NA, ncol = N, nrow = nrow(cs.temp))
      for (r in 1:nrow(L)) {
        yact.temp = apply(Y[cs.temp[r,1]:cs.temp[r,2], ,drop = F], 2, sum)
        lambda.temp = apply(Yhat[cs.temp[r,1]:cs.temp[r,2], ,drop = F], 2, sum)
        yact.temp[lambda.temp < 20] <- 20
        lambda.temp[lambda.temp < 20] <- 20
        L[r,] = (1-round(2*yact.temp/lambda.temp)/2)*lambda.temp + log((round(2 * (yact.temp/lambda.temp)) + 1e-04)/2.0001) * yact.temp
      }
      cs.loglikeij[s] = sum(L)

      term1 = cs.loglikeij[s]
      if(M == 0){
        term2 = 0
      } else{
        term2 = - M/2 * log(2 * cs.loglikeij[s] / M)
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
      cs.mBIC[s] <- mbic
    } else{
      cs.mBIC[s] = 0
    }
  }
  cs.mBIC <- round(cs.mBIC, digits = 3)

  cs.finalmat <- (cbind(cs.finalmat, cs.mBIC)[1:which.max(cs.mBIC), ,drop = F])
  poolcall <- create_chptsmat(cs.finalmat[,1:2,drop = F], chpts0)
  colnames(poolcall) = c("st", "end")

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
    # st_bp <- start(ref)[as.numeric(poolcall[i,"st"])]
    # ed_bp <- end(ref)[as.numeric(poolcall[i,"end"])]

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

    # temp = cbind(colnames(Y), rep(st_bin, ncol(Y)), rep(ed_bin, ncol(Y)), rep(st_bp, ncol(Y)),
    #              rep(ed_bp, ncol(Y)), chat.ind)
    temp = cbind(colnames(Y), rep(st_bin, ncol(Y)), rep(ed_bin, ncol(Y)), chat.ind)
    temp = temp[chat.ind != 2, ]
    finalcall = rbind(finalcall, temp)
    rownames(finalcall) = NULL
    # colnames(finalcall) = c("sample_name", "st_bin", "ed_bin", "st_bp", "ed_bp", "cnv_no")
    colnames(finalcall) = c("sample_name", "st_bin", "ed_bin", "cnv_no")
  }

  finalcall = as.data.frame(finalcall)

  finalcall$st_bin = as.numeric(paste(finalcall$st_bin))
  finalcall$ed_bin = as.numeric(paste(finalcall$ed_bin))
  # finalcall$st_bp = as.numeric(paste(finalcall$st_bp))
  # finalcall$ed_bp = as.numeric(paste(finalcall$ed_bp))
  finalcall$cnv_no = as.numeric(paste(finalcall$cnv_no))


  finalcall.idv = NULL
  iCN.idv = matrix(data=NA, nrow = nrow(Y), ncol = ncol(Y))

  if(segment.idv){
    # CBS segmenting by individuals
    for(sampno in 1:ncol(Y)) {
      message("Segmenting sample ", sampno, ": ", sampname[sampno], ".")
      chpts0 = c(1, nrow(Y))
      idv.scan = compute_lratio(Y[, sampno], Yhat[, sampno], chpts0)
      i = idv.scan$i
      j = idv.scan$j
      yact = idv.scan$yact
      lambda = idv.scan$lambda
      chat = idv.scan$chat
      lratio = idv.scan$lratio
      init.finalmat = idv.scan$finalmat

      if(!is.null(init.finalmat) && nrow(init.finalmat)>0){
        # Further nested searching
        chpts = init.finalmat[,1:2]
        if(is.null(dim(chpts))){
          chpts = t(as.matrix(chpts))
        }

        keep_going = 1
        # number of nested searching
        num_ns = 1
        while (!all(keep_going==0) & num_ns <= 2) {
          nested.output = search_nested(Y[, sampno], Yhat[, sampno], chpts)
          keep_going = sapply(nested.output, function(z) z$is.nested)

          newchpts = NULL
          for (r in 1:length(keep_going)) {
            if(keep_going[r]){
              newchpts = rbind(newchpts, nested.output[[r]]$finalmat[,1:2])
            } else{
              newchpts = rbind(newchpts, chpts[r,])
            }
          }
          chpts = newchpts
          num_ns = num_ns + 1
        }

        # backward compute lratio after nested searching
        temp = NULL
        for(r in 1:nrow(chpts)){
          idx = which(i == chpts[r,1] & j == chpts[r,2])
          if(length(idx)!=0){
            temp = rbind(temp, c(yact[idx], lambda[idx], chat[idx], lratio[idx]))
          } else{
            temp = rbind(temp, rep(NA, 4))
          }
        }
        finalmat = cbind(chpts, temp)
        # finalmat = finalmat[!is.na(finalmat[,6]) & finalmat[,6] > 0,]
        finalmat = finalmat[!is.na(finalmat[,6]), , drop = F]
        finalmat = finalmat[order(-finalmat[,6]), , drop = F]
      } else{ # if init.finalmat = NULL, the whole chromosome is normal
        finalmat = matrix(data = c(chpts0, sum(Y[,sampno]), sum(Yhat[,sampno]), 2, 0), nrow = 1, ncol = 6, byrow = TRUE)
      }

      loglikeij <- rep(NA, nrow(finalmat))
      mBIC <- rep(NA, length(loglikeij))

      for (s in 1:nrow(finalmat)) {
        tau <- sort(unique(c(as.vector(finalmat[1:s, 1:2]), chpts0)))
        P <- length(tau) - 2

        temp = create_chptsmat(finalmat[1:s,1:2,drop = F], chpts0)
        L = rep(NA, nrow(temp))
        for (r in 1:length(L)) {
          yact.temp = sum(Y[temp[r,1]:temp[r,2], sampno])
          lambda.temp = sum(Yhat[temp[r,1]:temp[r,2], sampno])
          yact.temp[lambda.temp < 20] <- 20
          lambda.temp[lambda.temp < 20] <- 20
          L[r] = (1-round(2*yact.temp/lambda.temp)/2)*lambda.temp + log((round(2 * (yact.temp/lambda.temp)) + 1e-04)/2.0001) * yact.temp
        }
        loglikeij[s] = sum(L)

        mbic <- loglikeij[s]
        mbic <- mbic - 0.5 * sum(log(tau[2:length(tau)] - tau[1:(length(tau) - 1)]))
        mbic <- mbic + (0.5 - P) * log(nrow(Y))
        mBIC[s] <- mbic
      }
      mBIC <- round(mBIC, digits = 3)

      if (mBIC[1] > 0) {
        finalmat <- (cbind(finalmat, mBIC)[1:which.max(mBIC), , drop = F])
        idv.seg.out <- create_chptsmat(finalmat[,1:2,drop = F], chpts0)
        idv.seg.out = cbind(idv.seg.out, rep(NA, nrow(idv.seg.out)))
        # output copy_no for the entire segments
        for(r in 1:nrow(idv.seg.out)){
          idx = which(i == idv.seg.out[r,1] & j == idv.seg.out[r,2])
          if(length(idx)!=0){
            idv.seg.out[r,3] = chat[idx]
          } else{
            # this segment consists of only one point
            idx2 = idv.seg.out[r,1]
            y.idv = Y[idx2,sampno]
            yhat.idv = Yhat[idx2,sampno]
            if(yhat.idv < 20){
              yhat.idv = 20
            }
            idv.seg.out[r,3] = round(2 * (y.idv/yhat.idv))
          }
          iCN.idv[idv.seg.out[r,1]:idv.seg.out[r,2],sampno] = idv.seg.out[r,3]
        }
        finalmat <- cbind(rep(sampname[sampno], nrow(finalmat)), finalmat)
        finalcall.idv <- rbind(finalcall.idv, finalmat)
      }
    }

    if (is.vector(finalcall.idv)) {
      finalcall.idv <- t(as.matrix(finalcall.idv))
    }

    colnames(finalcall.idv) <- c("sample_name", "st_bin", "ed_bin", "raw_cov", "norm_cov", "copy_no", "lratio", "mBIC")
    rownames(finalcall.idv) <- rep("", nrow(finalcall.idv))
    lratio = as.numeric(finalcall.idv[, "lratio"])
    rownames(finalcall.idv) = paste("cnv", 1:nrow(finalcall.idv), sep = "")
    finalcall.idv = as.data.frame(finalcall.idv)
    finalcall.idv$copy_no = as.numeric(paste(finalcall.idv$copy_no))
    finalcall.idv$raw_cov = as.numeric(paste(finalcall.idv$raw_cov))
    finalcall.idv$norm_cov = as.numeric(paste(finalcall.idv$norm_cov))
    finalcall.idv$st_bin = as.numeric(paste(finalcall.idv$st_bin))
    finalcall.idv$ed_bin = as.numeric(paste(finalcall.idv$ed_bin))
    finalcall.idv$lratio = as.numeric(paste(finalcall.idv$lratio))
    finalcall.idv$mBIC = as.numeric(paste(finalcall.idv$mBIC))
  }

  return(list(poolcall = poolcall, finalcall = finalcall, finalcall.idv = finalcall.idv, image.orig = image.orig, image.seg = image.seg, iCN = iCN, iCN.idv = iCN.idv))
}
