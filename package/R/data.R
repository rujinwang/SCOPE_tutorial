#' GRanges with mappability scores for hg19
#'
#' GRanges of mappability track for 100-mers on the GRCh37/hg19 human reference
#'  genome from ENCODE.
#'
#' @format A GRanges object with 21591667 ranges and 1 metadata column of mappability scores
#'
#' @references \url{http://rohsdb.cmb.usc.edu/GBshape/cgi-bin/hgFileUi?db=hg19&g=wgEncodeMapability}
"mapp_hg19"


#' GRanges with mappability scores for hg38
#'
#' Use liftOver utility to convert hg19 coordinates to hg38
#'
#' @format A GRanges object with 21584930 ranges and 1 metadata column of mappability scores
#'
#' @references \url{http://rohsdb.cmb.usc.edu/GBshape/cgi-bin/hgFileUi?db=hg19&g=wgEncodeMapability}
"mapp_hg38"
