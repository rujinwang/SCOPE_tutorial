#' @title Visualize EM fitting for each cell. BIC is used to choose optimal CNV group.
#'
#' @description Determine the number of CNV group via BIC for each. A pdf file containing EM fitting results and plots is generated.
#'
#' @param Y read depth matrix across all cells
#' @param gcfit.temp bin-specific GC content bias matrix across all single cells
#' @param gctemp a vector giving values of bin-specific GC content after \code{gcfit.tempj}
#'  outlier removal (for better fitting)
#' @param norm_index indices of normal/diploid cells
#' @param offset offset matrix across all cells
#' @param T a vector of integers indicating number of CNV groups. Use BIC to select optimal
#'  number of CNV groups. If \code{T = 1}, assume all reads are from normal regions so that
#'  EM algorithm is not implemented. Otherwise, we assume there is always a CNV group of heterozygous
#'  deletion and a group of null region. The rest groups are representative of different duplication
#'  states.
#' @param alpha adjusted/absolute copy number matrix across all cells
#' @param minCountQC the minimum read coverage required for EM fitting. Defalut is \code{20}
#' @param filename the name of output pdf file
#'
#' @return pdf file with EM fitting results and two plots: log likelihood, and BIC versus the number of CNV groups.
#'
#' @author Rujin Wang \email{rujin@email.unc.edu}
#' @import grDevices
#' @export
choiceofT = function (Y, gcfit.temp, gctemp, norm_index, offset, T, alpha, minCountQC = 20, filename){
	pdf(file = filename, width = 8, height = 10)
	for(j in 1:ncol(Y)){
		cat(j, '\t')
		if( j %in% norm_index){
			fGCj = getfGCj(gcfit.tempj = gcfit.temp[,j], gctemp = gctemp, Yj = Y[,j], offsetj = offset[,j], T = 1, draw.plot = TRUE, alphaj = alpha[,j], minCountQC)
		} else{
			fGCj = getfGCj(gcfit.tempj = gcfit.temp[,j], gctemp = gctemp, Yj = Y[,j], offsetj = offset[,j], T = T, draw.plot = TRUE, alphaj = alpha[,j], minCountQC)
		}
	}
	dev.off()
}
