#' @title Expectation step for SCOPE normalization
#'
#' @description Expectation step of GC content fitting for SCOPE normalization.
#'
#' @param Yj read depth vector for each single cell
#' @param Nj a numerirc value of total number of reads per cell
#' @param betatemp bin-specific bias estimated using negative control samples
#' @param fGCi estimated vector of GC content bias fitting from the M-step \code{\link{Mstep.Pois}}
#' @param vec_pi vector of incident rates for CNV events that span bin \emph{i}
#' from the M-step \code{\link{Mstep.Pois}}
#' @param min.prop the minmimum of mixture proportion for candidate CNV groups,
#' which serves as a stopping metric for EM algorithm
#'
#' @return A list with components
#'   \item{Z}{Matrix of optimized missing data to be fed into M-step \code{\link{Mstep.Pois}}}
#'   \item{obs_LL}{Observed log-Likelihood}
#'   \item{keep_going}{If normal proportion or a mixture proportion reached the minimum value,
#'    \code{keep_going = FALSE}, and output \code{NULL}, aka this case won't be involved in optimal number of CNV
#'    group selection based upon BIC.}
#'
#' @author Rujin Wang \email{rujin@email.unc.edu}
#' @importFrom stats dpois
#' @export
Estep.Pois = function(Yj, Nj, betatemp, fGCi, vec_pi, min.prop){
	keep_going = TRUE
	if (any(vec_pi<=min.prop)){
		message("A mixture proportion has gone to zero!!")
		keep_going = FALSE
	}
	if (vec_pi[2]<5e-2){
		message("Normal proportion has gone to zero!!")
		keep_going = FALSE
	}
	Z = matrix(nrow = length(fGCi), ncol = length(vec_pi))
	for (k in 1:length(vec_pi)){
		Z[,k] = log(vec_pi[k]) + dpois(Yj, lambda = Nj * betatemp * fGCi * (k/2), log = T)
	}
	obs_LL = 0
	for (ii in 1:length(fGCi)){
		tmp_num = logsumexp(Z[ii,])
		obs_LL = obs_LL + tmp_num
		Z[ii,] = exp(Z[ii,] - tmp_num)
	}
	return(list(Z = Z, obs_LL = obs_LL, keep_going = keep_going))
}
