#' @title Fit a loess curve of f(GC)
#'
#' @description Fit a loess curve of f(GC)
#'
#' @param gctemp a vector giving values of GC content for each bin after quality control
#' @param gcfit.temp vector of bin-specific GC content biases for each single cell
#'
#' @return
#'   \item{fGCi}{Estimated vector of GC content bias fitting}
#'
#' @author Rujin Wang \email{rujin@email.unc.edu}
#' @importFrom stats loess
#' @export
fitGC = function(gctemp, gcfit.temp){
  loe.fit = loess(gcfit.temp~gctemp)
  fGCi = loe.fit$fitted
  temp = min(fGCi[fGCi>0])
  fGCi[fGCi <= 0] = temp
  return(fGCi)
}
