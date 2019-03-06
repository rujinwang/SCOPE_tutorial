#' @title Expectation step for SCOPE normalization
#'
#' @description Expectation step of GC content fitting for SCOPE normalization.
#'
#' @param fGCi estimated vector of GC content bias fitting from the M-step \code{\link{Mstep}}
#' @param vec_pi vector of incident rates for CNV events that span bin \emph{i}
#' from the M-step \code{\link{Mstep}}
#' @param Yj read depth vector for each single cell
#' @param offsetj vector of offset terms for each cell
#'
#' @return
#'   \item{Z}{Matrix of optimized missing data to be fed into M-step \code{\link{Mstep}}}
#'
#' @author Rujin Wang \email{rujin@email.unc.edu}
#' @importFrom stats dpois
#' @export
Estep = function(fGCi, vec_pi, Yj, offsetj){
  P = matrix(nrow = length(fGCi), ncol = length(vec_pi))
  vec_pi[vec_pi==0] = 1e-100
  lambda =  matrix(nrow = nrow(P), ncol = ncol(P), data = offsetj*fGCi) *
    matrix(nrow = nrow(P),ncol = ncol(P), data = (1:length(vec_pi))/2, byrow = TRUE)
  P = dpois(matrix(nrow = nrow(P), ncol = ncol(P), data = Yj), lambda = lambda, log = T) +
    matrix(nrow = nrow(P),ncol = ncol(P), data = log(vec_pi), byrow = TRUE)
  Z = matrix(nrow = length(fGCi), ncol = length(vec_pi))
  for(k in 1:length(vec_pi)){
    Z[,k] = 1/(1+apply(exp(apply(P, 2, function(x){x-P[,k]})[,-k, drop = F]), 1, sum))
  }
  return(Z)
}
