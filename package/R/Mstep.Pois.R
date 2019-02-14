#' @title Maximization step for SCOPE normalization
#'
#' @description Maximization step of GC content fitting for SCOPE normalization.
#'
#' @param Z matrix of optimized missing data from the E-step \code{\link{Estep.Pois}}
#' @param gcfitj vector of bin-specific GC content biases for each single cell
#' @param gctemp a vector giving values of GC content for each bin after quality control
#'
#' @return A list with components
#'   \item{vec_pi}{Vector of incident rate for CNV events that span bin \emph{i},
#'    which is to be fed into E-step \code{\link{Estep.Pois}}}
#'   \item{fGCi}{Estimated vector of GC content bias fitting}
#'
#' @author Rujin Wang \email{rujin@email.unc.edu}
#' @importFrom stats smooth.spline predict
#' @export
Mstep.Pois = function(Z, gcfitj, gctemp){
	vec_pi = apply(Z, 2, mean)
	gcfit.temp = gcfitj / (Z %*% as.matrix((1:ncol(Z))/2))
	spl <- smooth.spline(gctemp[gcfitj>0], gcfit.temp[gcfitj>0], spar = 0.9)
	fGCi = pmax(1e-20, predict(spl, gctemp)$y)
	return(list(vec_pi = vec_pi, fGCi = fGCi))
}
