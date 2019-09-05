#' @title Maximization step for SCOPE normalization
#'
#' @description Maximization step of GC content fitting for SCOPE normalization.
#'
#' @param Z matrix of optimized missing data from the E-step \code{\link{Estep}}
#' @param gcfit.tempj vector of bin-specific GC content biases for each single cell
#' @param gctemp a vector giving values of GC content for each bin after quality control
#'
#' @return A list with components
#'   \item{vec_pi}{Vector of incident rate for CNV events that span bin \emph{i},
#'    which is to be fed into E-step \code{\link{Estep}}}
#'   \item{fGCi}{Estimated vector of GC content bias fitting}
#'
#' @author Rujin Wang \email{rujin@email.unc.edu}
#' @export
Mstep = function(Z, gcfit.tempj, gctemp){
  gcfit.temp = gcfit.tempj/(Z%*%as.matrix((1:ncol(Z))/2))
  fGCi = fitGC(gctemp, gcfit.temp)
  vec_pi = colSums(Z)/nrow(Z)
  return(list(vec_pi = vec_pi, fGCi = fGCi))
}
