% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/getmapp.R
\name{getmapp}
\alias{getmapp}
\title{Compute mappability}
\usage{
getmapp(ref, genome = NULL)
}
\arguments{
\item{ref}{GRanges object returned from \code{getbambed_scope}}

\item{genome}{by default, \code{genome = BSgenome.Hsapiens.UCSC.hg19}.
To calculate mappability for hg38, specify
\code{genome = BSgenome.Hsapiens.UCSC.hg38}}
}
\value{
\item{mapp}{Vector of mappability for each bin/target}
}
\description{
Compute mappability for each bin. Note that scDNA
sequencing is whole-genome amplification and the mappability
score is essential to determine variable binning method.
Mappability track for 100-mers on the GRCh37/hg19 human
reference genome from ENCODE is pre-saved. Compute the mean
of mappability scores that overlapped reads map to bins,
weighted by the width of mappability tracks on the genome
reference. Use liftOver utility to calculate mappability
for hg38, which is pre-saved as well.
}
\examples{
library(WGSmapp)
bedFile <- system.file('extdata',
                        'scWGA500kbsort.bed',
                        package = 'SCOPE')
bamfolder <- system.file('extdata', package = 'WGSmapp')
bamFile <- list.files(bamfolder, pattern = '*.dedup.bam$')
bamdir <- file.path(bamfolder, bamFile)
sampname_raw = sapply(strsplit(bamFile, '.', fixed = TRUE), '[', 1)
bambedObj <- getbambed_scope(bamdir = bamdir,
                            bedFile = bedFile,
                            sampname = sampname_raw)
bamdir <- bambedObj$bamdir
sampname_raw <- bambedObj$sampname
ref_raw <- bambedObj$ref

data('mapp_hg19')
mapp <- getmapp(ref_raw, genome = BSgenome.Hsapiens.UCSC.hg19)

}
\author{
Rujin Wang \email{rujin@email.unc.edu}
}
