#' @title Logarithm of summation of exponentials
#'
#' @description Computes the logarithm of summation of numeric exponentials.
#'
#' @param xx a numeric value or vector
#'
#' @return A numeric value giving natural logarithm of summation of exponentials
#'
#' @author Rujin Wang \email{rujin@email.unc.edu}
#'
#' @export
logsumexp = function(xx){
	if(length(xx)==1){
		xx
	}else{
		max_val = max(xx)
		xx2 = xx - max_val
		log(sum(exp(xx2))) + max_val
	}
}
