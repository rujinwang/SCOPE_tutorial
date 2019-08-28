#' @title Ploidy pre-initialization
#'
#' @description Pre-estimate ploidies across all cells
#'
#' @param Y raw read depth matrix after quality control procedure
#' @param Yhat normalized read depth matrix
#' @param ref GRanges object after quality control procedure
#' @param maxPloidy maximum ploidy candidate. Defalut is \code{6}
#' @param minPloidy minimum ploidy candidate. Defalut is \code{1.5}
#' @param minBinWidth the minimum number of bins for a changed segment. Defalut is \code{5}
#'
#' @return
#'   \item{ploidy.SoS}{Vector of pre-estimated ploidies for each cell}
#'
#' @author Rujin Wang \email{rujin@email.unc.edu}
#' @importFrom DNAcopy CNA smooth.CNA segment
#' @importFrom GenomeInfoDb seqnames
#' @importFrom BSgenome start
#' @importFrom IRanges end
#' @export
PreEst_ploidy = function(Y, Yhat, ref, maxPloidy = 6, minPloidy = 1.5, minBinWidth = 5){
	ploidy.SoS=rep(NA, ncol(Y))

	breaks = matrix(0,nrow(Y),ncol(Y))
	RCNP = matrix(0,nrow(Y),ncol(Y))
	final = matrix(0,nrow(Y),ncol(Y))
	X = seq(minPloidy, maxPloidy, by=0.05)
	n_ploidy = length(X)
	SoS = matrix(0,n_ploidy,ncol(Y))

	normal = (Y+1) / (Yhat + 1)

	for(k in 1:ncol(Y)){
		if(k%%5 == 1){
			cat(k, "\t")
		}

		lr = log(normal[,k])
		loc = data.frame(seq = as.character(seqnames(ref)), start = start(ref), end = end(ref))
		CNA.object   = CNA(genomdat = lr, chrom = loc[,1], maploc = as.numeric(loc[,2]), data.type = 'logratio')
		CNA.smoothed = smooth.CNA(CNA.object)
		segs         = segment(CNA.smoothed, verbose=0, min.width=minBinWidth)
		frag         = segs$output[,2:3]
		len = dim(frag)[1]
		bps = array(0, len)
		for (j in 1:len){
			bps[j]=which((loc[,1]==frag[j,1]) & (as.numeric(loc[,2])==frag[j,2]))
		}
		bps = sort(bps)
		bps[(len=len+1)] = nrow(Y)
		breaks[bps,k] = 1
  	RCNP[,k][1:bps[2]] = median(normal[,k][1:bps[2]])
		for(i in 2:(len-1)){
			RCNP[,k][bps[i]:(bps[i+1]-1)] = median(normal[,k][bps[i]:(bps[i+1]-1)])
		}
		RCNP[,k] = RCNP[,k]/mean(RCNP[,k])

		SCNP = RCNP[,k] %o% X
		FSCP = round(SCNP)
		Diff2 = (SCNP - FSCP) ^ 2
		SoS[,k] = colSums(Diff2, na.rm = FALSE, dims = 1)
		ploidy.SoS[k] = X[which.min(SoS[,k])]
	}
	return(ploidy.SoS)
}
