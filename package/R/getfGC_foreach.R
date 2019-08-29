if(getRversion() >= "2.15.1"){
  utils::globalVariables(c("Ti"))
}
#' @title Get GC content bias fitting using Expectation-Maximization algorithm across all cells in parallel
#'
#' @description Fit a Poisson generalized linear model to normalize the raw read depth data from
#'  single-cell DNA sequencing, under the case-control setting. SCOPE
#'  implements an EM algorithm with random initialization to unmask null regions. This is a
#'  bulit-in function within \code{\link{normalize_scope_foreach}}.
#'
#' @param gcfit.temp bin-specific GC content bias matrix across all single cells
#' @param gctemp a vector giving values of bin-specific GC content after \code{gcfit.tempj}
#'  outlier removal (for better fitting)
#' @param Y read depth matrix across all cells
#' @param norm_index indices of normal/diploid cells
#' @param offset offset matrix across all cells
#' @param T a vector of integers indicating number of CNV groups. Use BIC to select optimal
#'  number of CNV groups. If \code{T = 1}, assume all reads are from normal regions so that
#'  EM algorithm is not implemented. Otherwise, we assume there is always a CNV group of heterozygous
#'  deletion and a group of null region. The rest groups are representative of different duplication
#'  states.
#' @param alpha adjusted/absolute copy number matrix across all cells
#' @param minCountQC the minimum read coverage required for EM fitting. Defalut is \code{20}
#'
#' @return A list with components
#'   \item{fGC.hat}{EM estimated GC content bias matrix given optimal CNV group selection}
#'   \item{alpha}{Absolute copy number matrix adjusted by EM fitting}
#'
#' @author Rujin Wang \email{rujin@email.unc.edu}
#' @import stats foreach parallel doParallel
#' @export
getfGC_foreach = function(gcfit.temp, gctemp, Y, norm_index, offset, T, alpha, minCountQC){
	fGC.hat = matrix(ncol = ncol(Y), nrow = nrow(Y))
	for(j in 1:ncol(Y)){
		cat(j, '\t')
		if( j %in% norm_index){
			alpha[,j] = 2
			loe.fit = loess(gcfit.temp[,j]~gctemp)
			fGC.hat[,j] = loe.fit$fitted
		} else{
			# begin foreach
			nCores <- detectCores() - 1
			registerDoParallel(nCores)

			TList <- foreach(Ti = T, .export = c("fitGC", "Estep", "Mstep", "getfGCj_foreach")) %dopar% {
				getfGCj_foreach(gcfit.tempj = gcfit.temp[,j], gctemp = gctemp, Yj = Y[,j], offsetj = offset[,j], T = Ti, draw.plot = F, alphaj = alpha[,j], minCountQC = minCountQC)
			}

			# When you're done, clean up the cluster
			stopImplicitCluster()
			#end foreach

			fGCi.obj = lapply(TList, function(x) x[["fGCi.obj"]])
			Z.obj = lapply(TList, function(x) x[["Z.obj"]])
			vec_pi.obj = lapply(TList, function(x) x[["vec_pi.obj"]])
			loglik = sapply(TList, function(x) x[["loglik"]])
			BIC = sapply(TList, function(x) x[["BIC"]])

			fGCj = list(fGCi.obj = fGCi.obj, Z.obj = Z.obj, vec_pi.obj = vec_pi.obj, loglik = loglik, BIC = BIC)
			if(which.max(fGCj$BIC)==1){
				alpha[,j] = 2
			} else{
				alpha[,j] = apply(fGCj$Z.obj[[which.max(fGCj$BIC)]],1,which.max) # This is to update alpha
			}
			fGC.hat[,j] = fGCj$fGCi.obj[[which.max(fGCj$BIC)]]
		}
	}
	return(list(fGC.hat = fGC.hat, alpha = alpha))
}
