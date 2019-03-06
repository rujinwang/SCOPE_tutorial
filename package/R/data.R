#' @name mapp_hg19
#' @aliases mapp_hg19
#' @docType data
#' @title GRanges with mappability scores for hg19
#'
#' @description GRanges object specifying target positions with mappabilities across the
#'  whole genome.
#'
#' @return GRanges object with mappabilities for hg19
#'
#' @details GRanges of mappability track for 100-mers on the GRCh37/hg19 human reference
#'  genome from ENCODE.
#'
#' @usage data(mapp_hg19)
#'
#' @author Rujin Wang \email{rujin@email.unc.edu}
#'
#' @format A GRanges object with 21591667 ranges and 1 metadata column of mappability scores
#'
#' @keywords datasets
#'
#' @references \url{http://rohsdb.cmb.usc.edu/GBshape/cgi-bin/hgFileUi?db=hg19&g=wgEncodeMapability}
#'
NULL


#' @name mapp_hg38
#' @aliases mapp_hg38
#' @docType data
#' @title GRanges with mappability scores for hg38
#'
#' @description GRanges object specifying target positions with mappabilities across the
#'  whole genome.
#'
#' @return GRanges object with mappabilities for hg38
#'
#' @details Use liftOver utility to convert hg19 coordinates to hg38
#'
#' @usage data(mapp_hg38)
#'
#' @author Rujin Wang \email{rujin@email.unc.edu}
#'
#' @format A GRanges object with 21584930 ranges and 1 metadata column of mappability scores
#'
#' @keywords datasets
#'
#' @references \url{http://hgdownload.cse.ucsc.edu/goldenPath/hg19/liftOver/}
#'
NULL
