#' @title Normalization of read depth without latent factors using Expectation-Maximization algorithm
#'  under the case-control setting
#'
#' @description Fit a Poisson generalized linear model to normalize the raw read depth data from
#'  single-cell DNA sequencing, without latent factors under the case-control setting. Model GC
#'  content bias using an expectation-maximization algorithm, which accounts for the different
#'  copy number states.
#'
#' @param Y_qc read depth matrix after quality control
#' @param gc_qc vector of GC content for each bin after quality control
#' @param K default is \code{K = 1}
#' @param norm_index indices of normal/diploid cells
#' @param N library size factor, which is computed from the genome-wide read depth data
#' @param numGroup a vector of integers indicating number of CNV groups. Use BIC to select optimal
#'  number of CNV groups. If \code{numGroup = 1}, assume all reads are from normal regions so that
#'  EM algorithm is not implemented. Otherwise, we assume there is always a CNV group of heterozygous
#'  deletion and a group of null region. The rest groups are representative of different duplication
#'  states.
#' @param qc.thres the lower bound of bin-specific GC content bias threshold
#' @param min.prop the minmimum of mixture proportion for candidate CNV groups,
#' which serves as a stopping metric for EM algorithm
#'
#' @return A list with components
#'   \item{Yhat}{A list of normalized read depth matrix with EM}
#'   \item{fGC.hat}{A list of EM estimated GC content bias matrix}
#'   \item{beta.hat}{A list of EM estimated bin-specific bias vector}
#'
#' @author Rujin Wang \email{rujin@email.unc.edu}
#' @import stats
#' @export
normalize_codex2_ns_noK_EM_random = function (Y_qc, gc_qc, K = 1, norm_index, N, numGroup, qc.thres = 5e-5, min.prop = 2e-3) {
  Ntotal <- N
  N <- round(N/median(N)*median(colSums(Y_qc)))
  Nmat <- matrix(nrow = nrow(Y_qc), ncol = ncol(Y_qc), data = N, byrow = TRUE)
  Yhat = list(length = length(K))
  fGC.hat <- list(length = length(K))
  beta.hat <- list(length = length(K))
  for (ki in 1:length(K)) {
    k = K[ki]
    message("Computing normalization with no latent factors ...")
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
      if(iter>2){
        for(j in 1:ncol(gcfit)){
          cat(j,'\t')
		  if (j %in% norm_index){
			fGC.obj = multi_run.Pois(gcfitj = gcfit[,j], gctemp=gc_qc, Yj=Y_qc[,j], Nj=N[j], betatemp=betahatmat[,j], numGroup=1, rerun=10, qc.thres = qc.thres, min.prop = min.prop)
		  } else{
			fGC.obj = multi_run.Pois(gcfitj = gcfit[,j], gctemp=gc_qc, Yj=Y_qc[,j], Nj=N[j], betatemp=betahatmat[,j], numGroup=numGroup, rerun=10, qc.thres = qc.thres, min.prop = min.prop)
		  }
		  temp=fGC.obj$fGCi.keep
          temp[temp <= 0] <- min(temp[temp > 0])
          fhatnew[,j]=temp
          cat('\n')
        }
      } else{
        fhatnew <- apply(gcfit, 2, function(gcfitj) {
          spl <- smooth.spline(gc_qc, gcfitj)
          temp <- predict(spl, gc_qc)$y
          temp[temp <= 0] <- min(temp[temp > 0])
          temp})
      }
      betahatnew = apply((Y_qc/(fhatnew * Nmat))[,norm_index], 1, median)
      betahatnew[betahatnew <= 0] = min(betahatnew[betahatnew > 0])
      bhdiff[iter] = sum((betahatnew - betahat)^2)/length(betahat)
      fhdiff[iter] = sum((fhatnew - fhat)^2)/length(fhat)
      message("Iteration ", iter, "\t", "beta diff =", signif(bhdiff[iter],
                                                              3), "\t", "f(GC) diff =", signif(fhdiff[iter], 3))
      fhat = fhatnew
      betahat = betahatnew
      betahatmat = matrix(nrow = nrow(Y_qc), ncol = ncol(Y_qc),
                          data = betahat, byrow = FALSE)
      fhatlist[[iter]] = fhat
      betahatlist[[iter]] = betahat
      if (bhdiff[iter] < BHTHRESH & iter >2)
        break
      if (iter > 5 & bhdiff[iter] > 1)
        break
      iter = iter + 1
    }
    optIter = which.min(fhdiff)
    message(paste("Stop at Iteration ", optIter, ".", sep = ""))
    if(iter > maxiter){
      iter = iter - 1
    }
    fhat = fhatlist[[iter]]
    betahat = betahatlist[[iter]]
    betahatmat = matrix(nrow = nrow(Y_qc), ncol = ncol(Y_qc),
                        data = betahat, byrow = FALSE)
    Yhat[[ki]] = pmax(round(fhat * Nmat * betahatmat,0),1)
    fGC.hat[[ki]] <- signif(fhat, 3)
    beta.hat[[ki]] <- signif(betahat, 3)
  }
  list(Yhat = Yhat, fGC.hat = fGC.hat, beta.hat = beta.hat)
}
