#' @title Normalization of read depth without latent factors under the case-control setting
#'
#' @description Assuming that all reads are from diploid regions, fit a Poisson generalized
#'  linear model to normalize the raw read depth data from
#'  single-cell DNA sequencing, without latent factors under the case-control setting.
#'
#' @param Y_qc read depth matrix after quality control
#' @param gc_qc vector of GC content for each bin after quality control
#' @param K default is \code{K = 1}
#' @param norm_index indices of normal/diploid cells
#' @param N library size factor, which is computed from the genome-wide read depth data
#'
#' @return A list with components
#'   \item{Yhat}{A list of normalized read depth matrix}
#'   \item{fGC.hat}{A list of estimated GC content bias matrix}
#'   \item{beta.hat}{A list of estimated bin-specific bias vector}
#'
#' @author Rujin Wang \email{rujin@email.unc.edu}
#' @import stats
#' @export
normalize_codex2_ns_noK = function (Y_qc, gc_qc, K = 1, norm_index, N) {
  Ntotal <- N
  N <- round(N/median(N)*median(colSums(Y_qc)))
  Nmat <- matrix(nrow = nrow(Y_qc), ncol = ncol(Y_qc), data = N, byrow = TRUE)
  Yhat = list(length = length(K))
  fGC.hat <- list(length = length(K))
  beta.hat <- list(length = length(K))
  for (ki in 1:length(K)) {
    k = K[ki]
    message("Computing normalization with no latent factors")
    message("k = ", k)
    maxiter = 10
    maxhiter = 50
    BHTHRESH = 1e-04
    HHTHRESH = 1e-05
    iter = 1
    fhat = matrix(nrow = nrow(Y_qc), ncol = ncol(Y_qc), data = 0)
    fhatnew = matrix(nrow = nrow(Y_qc), ncol = ncol(Y_qc))
    betahat = rep(1, nrow(Y_qc))
    betahatmat = matrix(nrow = nrow(Y_qc), ncol = ncol(Y_qc),
                        data = betahat, byrow = FALSE)
    bhdiff = rep(Inf, maxiter)
    fhdiff = rep(Inf, maxiter)
    betahatlist = list(length = maxiter)
    fhatlist = list(length = maxiter)
    while (iter <= maxiter) {
      gcfit = Y_qc/Nmat/betahatmat
      fhatnew <- apply(gcfit, 2, function(z) {
        spl <- smooth.spline(gc_qc, z)
        temp <- predict(spl, gc_qc)$y
        temp[temp <= 0] <- min(temp[temp > 0])
        temp
      })
      fhatnew[fhatnew<quantile(fhatnew, 0.005)]=quantile(fhatnew, 0.005)
      betahatnew = apply((Y_qc/(fhatnew * Nmat))[,norm_index], 1, median)
      betahatnew[betahatnew <= 0] = min(betahatnew[betahatnew > 0])
      bhdiff[iter] = sum((betahatnew - betahat)^2)/length(betahat)
      fhdiff[iter] = sum((fhatnew - fhat)^2)/length(fhat)
      if (fhdiff[iter] > min(fhdiff))
        break
      message("Iteration ", iter, "\t", "beta diff =", signif(bhdiff[iter],
                                                              3), "\t", "f(GC) diff =", signif(fhdiff[iter], 3))
      fhat = fhatnew
      betahat = betahatnew
      betahatmat = matrix(nrow = nrow(Y_qc), ncol = ncol(Y_qc),
                          data = betahat, byrow = FALSE)
      fhatlist[[iter]] = fhat
      betahatlist[[iter]] = betahat
      if (bhdiff[iter] < BHTHRESH)
        break
      if (iter > 5 & bhdiff[iter] > 1)
        break
      iter = iter + 1
    }
    optIter = which.min(fhdiff)
    message(paste("Stop at Iteration ", optIter, ".", sep = ""))
    fhat = fhatlist[[optIter]]
    betahat = betahatlist[[optIter]]
    betahatmat = matrix(nrow = nrow(Y_qc), ncol = ncol(Y_qc),
                        data = betahat, byrow = FALSE)
    Yhat[[ki]] = pmax(round(fhat * Nmat * betahatmat,0),1)
    fGC.hat[[ki]] <- signif(fhat, 3)
    beta.hat[[ki]] <- signif(betahat, 3)
  }
  list(Yhat = Yhat, fGC.hat = fGC.hat, beta.hat = beta.hat)
}
