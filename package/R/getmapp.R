if(getRversion() >= "2.15.1"){
  utils::globalVariables(c("BSgenome.Hsapiens.UCSC.hg19", "mapp_hg19", "mapp_hg38", "seqlevelsStyle<-"))
}
#' @title Compute mappability
#' @name getmapp
#'
#' @description Compute mappability for each bin. Note that scDNA sequencing is
#' whole-genome amplification and the mappability score is essential to
#' determine variable binning method. Mappability track for 100-mers on the
#' GRCh37/hg19 human reference genome from ENCODE is pre-saved. Compute the mean
#' of mappability scores that overlapped reads map to bins, weighted by the width
#' of mappability tracks on the genome reference. Use liftOver utility to calculate
#' mappability for hg38, which is pre-saved as well.
#'
#' @param ref GRanges object returned from \code{\link[CODEX2]{getbambed}}
#' @param genome by default, \code{genome = BSgenome.Hsapiens.UCSC.hg19}.
#' To calculate mappability for hg38, specify \code{genome = BSgenome.Hsapiens.UCSC.hg38}
#'
#' @return
#'   \item{mapp}{Vector of mappability for each bin/target}
#'
#' @author Rujin Wang \email{rujin@email.unc.edu}
#' @import BSgenome.Hsapiens.UCSC.hg19 utils
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges RangesList Views countOverlaps findOverlaps width
#' @importFrom GenomeInfoDb mapSeqlevels seqlevelsStyle seqnames
#' @export
getmapp = function (ref, genome = NULL) {
  if(is.null(genome)){genome = BSgenome.Hsapiens.UCSC.hg19}
  if(genome@provider_version == 'hg19'){mapp_gref = mapp_hg19}
  if(genome@provider_version == 'hg38'){mapp_gref = mapp_hg38}
  mapp <- rep(1, length(ref))
  seqlevelsStyle(ref)='UCSC'
  for(chr in unique(seqnames(ref))){
    message("Getting mappability for ", chr, sep = "")
    chr.index=which(as.matrix(seqnames(ref))==chr)
    ref.chr=GRanges(seqnames=chr, ranges=IRanges(start= start(ref)[chr.index] , end = end(ref)[chr.index]))
    mapp.chr=rep(1, length(ref.chr))
    overlap=as.matrix(findOverlaps(ref.chr, mapp_gref))
    for(i in unique(overlap[,1])){
      index.temp=overlap[which(overlap[,1]==i),2]
      mapp.chr[i]=sum((mapp_gref$score[index.temp])*(width(mapp_gref)[index.temp]))/
        sum(width(mapp_gref)[index.temp])
      #mapp.chr[i]=mean(mapp_gref$score[index.temp])
    }
  mapp[chr.index]=mapp.chr
  }
  mapp
}
