#' @title Get bam file directories, sample names, and whole genomic bins
#'
#' @description Get bam file directories, sample names, and whole genomic bins from .bed file
#'
#' @param bamdir vector of the directory of a bam file. Should be in the same order as sample names in \code{sampname}.
#' @param bedFile path to the whole genome bed file
#' @param sampname vector of sample names. Should be in the same order as bam directories in \code{bamdir}.
#'
#' @return A list with components
#'   \item{bamdir}{A vector of bam directories}
#'   \item{sampname}{A vector of sample names}
#'   \item{ref}{A GRanges object specifying whole genomic bin positions}
#'
#' @author Rujin Wang \email{rujin@email.unc.edu}
#' @import utils
#' @import GenomeInfoDb
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges
#' @export
getbambed_scope = function(bamdir, bedFile, sampname){
  if (!file.exists(bedFile))
    stop("Please check the bed file directory provided. File could not be \nfound!")
  exomtarg <- read.table(bedFile, sep = "\t")
  ref <- GRanges(seqnames = exomtarg[,1], ranges = IRanges(start = exomtarg[,2], end = exomtarg[,3]))
  if(!any(grepl("chr", seqlevels(ref)))){
    seqlevels(ref) = paste(c(1:22, "X", "Y"), sep = "")
    ref <- sort(ref)
  } else{
    seqlevels(ref) = paste("chr", c(1:22, "X", "Y"), sep = "")
    ref <- sort(ref)
  }
  list(bamdir = bamdir, sampname = sampname, ref = ref)
}

