#' @title Normalization of read depth with latent factors using Expectation-Maximization algorithm
#'  under the case-control setting
#'
#' @description Fit a Poisson generalized linear model to normalize the raw read depth data from
#'  single-cell DNA sequencing, with latent factors under the case-control setting. Model GC
#'  content bias using an expectation-maximization algorithm, which accounts for the different
#'  copy number states.
#'
#' @param Y_qc read depth matrix after quality control
#' @param gc_qc vector of GC content for each bin after quality control
#' @param K Number of latent Poisson factors
#' @param norm_index indices of normal/diploid cells
#' @param N library size factor, which is computed from the genome-wide read depth data
#' @param T a vector of integers indicating number of CNV groups. Use BIC to select optimal
#'  number of CNV groups. If \code{T = 1}, assume all reads are from normal regions so that
#'  EM algorithm is not implemented. Otherwise, we assume there is always a CNV group of heterozygous
#'  deletion and a group of null region. The rest groups are representative of different duplication
#'  states.
#' @param alpha0 initialized matrix of absolute copy number
#' @param beta0 a vector of initialized bin-specific biases returned from CODEX2 without latent factors
#' @param minCountQC the minimum read coverage required for normalization and EM fitting. Defalut is \code{20}
#'
#' @return A list with components
#'   \item{Yhat}{A list of normalized read depth matrix with EM}
#'   \item{alpha.hat}{A list of absolute copy number matrix}
#'   \item{fGC.hat}{A list of EM estimated GC content bias matrix}
#'   \item{beta.hat}{A list of EM estimated bin-specific bias vector}
#'   \item{g.hat}{A list of estimated Poisson latent factor}
#'   \item{h.hat}{A list of estimated Poisson latent factor}
#'   \item{AIC}{AIC for model selection}
#'   \item{BIC}{BIC for model selection}
#'   \item{RSS}{RSS for model selection}
#'   \item{K}{Number of latent Poisson factors}
#'
#' @author Rujin Wang \email{rujin@email.unc.edu}
#' @import stats
#' @export
normalize_scope = function (Y_qc, gc_qc, K, norm_index, N, T, alpha0, beta0, minCountQC = 20) {
  if (max(K) > length(norm_index))
    stop("Number of latent Poisson factors K cannot exceed the number of
         normal samples!")
  Ntotal <- N
  N <- round(N/median(N)*median(colSums(Y_qc)))
  Nmat <- matrix(nrow = nrow(Y_qc), ncol = ncol(Y_qc), data = N, byrow = TRUE)

  Yhat = vector('list', length(K))
  fGC.hat <- vector('list', length(K))
  alpha.hat <- vector('list', length(K))
  beta.hat <- vector('list', length(K))
  g.hat <- vector('list', length(K))
  h.hat <- vector('list', length(K))
  AIC <- rep(NA, length = length(K))
  BIC <- rep(NA, length = length(K))
  RSS <- rep(NA, length = length(K))

  # Initialization
  message("Initialization ...")
  gcfit.temp=Y_qc/Nmat/beta0
  offset=Nmat*matrix(nrow = nrow(Y_qc), ncol = ncol(Y_qc), data = beta0, byrow = FALSE)
  fhat.temp=getfGC(gcfit.temp = gcfit.temp, gctemp = gc_qc, Y = Y_qc, norm_index = norm_index, offset = offset, T = T, alpha = alpha0, minCountQC = minCountQC)
  fhat0=fhat.temp$fGC.hat
  alpha0=fhat.temp$alpha

  for (ki in 1:length(K)) {
    k = K[ki]
    message("Computing normalization with k = ", k, " latent factors ...",sep="")
    message("k = ", k)
    maxiter = 10
    maxhiter = 50
    BHTHRESH = 1e-04
    HHTHRESH = 1e-05
    iter = 1
    fhat = matrix(nrow = nrow(Y_qc), ncol = ncol(Y_qc), data = 0)
    betahat = beta0
    betahatmat = matrix(nrow = nrow(Y_qc), ncol = ncol(Y_qc),
                        data = betahat, byrow = FALSE)
    ghat = matrix(0, nrow = nrow(Y_qc), ncol = k)
    hhat = matrix(0, nrow = ncol(Y_qc), ncol = k)
    bhdiff = rep(Inf, maxiter)
    fhdiff = rep(Inf, maxiter)

    betahatlist = vector('list', maxiter)
    fhatlist = vector('list', maxiter)
    ghatlist = vector('list', maxiter)
    hhatlist = vector('list', maxiter)
    alphahatlist = vector('list', maxiter)

    while (iter <= maxiter) {
      if(iter==1){fhatnew=fhat0; alpha=alpha0}
      if(iter>1){
        gcfit.temp=Y_qc/Nmat/betahat/exp(ghat %*% t(hhat))
        offset=Nmat*matrix(nrow = nrow(Y_qc), ncol = ncol(Y_qc),
                           data = betahat, byrow = FALSE)*exp(ghat %*% t(hhat))
        fhat.temp=getfGC(gcfit.temp=gcfit.temp, gctemp=gc_qc, Y=Y_qc, norm_index=norm_index, offset=offset, T = T, alpha=alpha0, minCountQC = minCountQC)
        fhatnew=fhat.temp$fGC.hat
        alpha=fhat.temp$alpha
      }
      fhatnew[fhatnew<quantile(fhatnew, 0.005)]=quantile(fhatnew, 0.005)
      betahatnew = apply((Y_qc/(fhatnew * Nmat * exp(ghat %*% t(hhat))))[,norm_index], 1, median)
      betahatnew[betahatnew <= 0] = min(betahatnew[betahatnew > 0])
      bhdiff[iter] = sum((betahatnew - betahat)^2)/length(betahat)
      fhdiff[iter] = sum((fhatnew - fhat)^2)/length(fhat)
      if (fhdiff[iter] > min(fhdiff))
        break
      message("Iteration ", iter, "\t", "beta diff =",
              signif(bhdiff[iter], 3), "\t", "f(GC) diff =", signif(fhdiff[iter], 3))
      fhat = fhatnew
      betahat = betahatnew
      betahatmat = matrix(nrow = nrow(Y_qc), ncol = ncol(Y_qc),
                          data = betahat, byrow = FALSE)
      L = log(Nmat * fhat * betahatmat * alpha / 2)
      logmat = log(pmax(Y_qc, 1)) - L
      logmat = logmat - matrix(nrow = nrow(Y_qc), ncol = ncol(Y_qc), data = apply(logmat, 1, mean), byrow = FALSE)
      hhat = svd(logmat, nu = k, nv = k)$v
      hhatnew = hhat
      hiter = 1
      hhdiff = rep(Inf, maxhiter)
      while (hiter <= maxhiter) {
        for (s in 1:nrow(Y_qc)) {
          temp=try(glm(formula = Y_qc[s,norm_index] ~ hhat[norm_index,] -
                         1, offset = L[s,norm_index], family = poisson)$coefficients,silent=TRUE)
          if(is.character(temp)){
            temp=lm(log(pmax(Y_qc[s,norm_index],1)) ~ hhat[norm_index,] -
                      1, offset = log(L[s,norm_index]))$coefficients
          }
          ghat[s, ] = temp
        }
        # avoid overflow or underflow of the g latent factors
        ghat[is.na(ghat)]=0
        if(max(ghat) >= 30){
          ghat=apply(ghat,2,function(z){
            z[z>quantile(z,0.995)] = min(quantile(z,0.995),30)
            z})
        }
        if(min(ghat) <= -30){
          ghat=apply(ghat,2,function(z){
            z[z<quantile(z,0.005)] = max(quantile(z,0.005),-30)
            z})
        }
        for (t in 1:ncol(Y_qc)) {
          hhatnew[t, ] = glm(formula = Y_qc[, t] ~ ghat -
                               1, offset = L[, t], family = poisson)$coefficients
        }
        gh = ghat %*% t(hhatnew)
        gh <- scale(gh, center = TRUE, scale = FALSE)
        hhatnew = svd(gh, nu = k, nv = k)$v
        hhdiff[hiter] = sum((hhatnew - hhat)^2)/length(hhat)
        message("\t\t\t", "hhat diff =", signif(hhdiff[hiter], 3))
        hhat = hhatnew
        if (hhdiff[hiter] < HHTHRESH)
          break
        if (hiter > 10 & (rank(hhdiff))[hiter] <= 3)
          break
        hiter = hiter + 1
      }
      alphahatlist[[iter]] = alpha
      fhatlist[[iter]] = fhat
      betahatlist[[iter]] = betahat
      ghatlist[[iter]] = ghat
      hhatlist[[iter]] = hhat
      if (bhdiff[iter] < BHTHRESH)
        break
      if (iter > 5 & bhdiff[iter] > 1)
        break
      iter = iter + 1
    }
    optIter = which.min(fhdiff)
    message(paste("Stop at Iteration ", optIter, ".", sep = ""))
    alpha = alphahatlist[[optIter]]
    fhat = fhatlist[[optIter]]
    betahat = betahatlist[[optIter]]
    ghat = ghatlist[[optIter]]
    hhat = hhatlist[[optIter]]
    betahatmat = matrix(nrow = nrow(Y_qc), ncol = ncol(Y_qc),
                        data = betahat, byrow = FALSE)
    Yhat[[ki]] = pmax(round(fhat * Nmat * betahatmat * exp(ghat %*% t(hhat)),0),1)
    alpha.hat[[ki]] <- alpha
    fGC.hat[[ki]] <- signif(fhat, 3)
    beta.hat[[ki]] <- signif(betahat, 3)
    h.hat[[ki]] <- signif(hhat, 3)
    g.hat[[ki]] <- signif(ghat, 3)
    Yhat.temp=Yhat[[ki]]*alpha/2
    AIC[ki] = 2 * sum(Y_qc * log(pmax(Yhat.temp,1)) - Yhat.temp) -
      2 * (length(ghat) + length(hhat))
    BIC[ki] = 2 * sum(Y_qc * log(pmax(Yhat.temp,1)) - Yhat.temp) -
      (length(ghat) + length(hhat)) * log(length(Y_qc))
    RSS[ki] = sum((Y_qc - Yhat.temp)^2/length(Y_qc))
    message("AIC", k, " = ", round(AIC[ki], 3))
    message("BIC", k, " = ", round(BIC[ki], 3))
    message("RSS", k, " = ", round(RSS[ki], 3), "\n")
  }
  list(Yhat = Yhat, alpha.hat = alpha.hat, fGC.hat = fGC.hat, beta.hat = beta.hat,
       g.hat = g.hat, h.hat = h.hat, AIC = AIC, BIC = BIC, RSS = RSS,
       K = K)
}
