getcoverage.scDNA=function (bambedObj, mapqthres, mask.ref, seq) 
{
	ref <- bambedObj$ref
	bamdir <- bambedObj$bamdir
	sampname <- bambedObj$sampname
  
	Y=matrix(nrow=length(ref), ncol=length(sampname))
	rownames(Y) = paste(seqnames(ref), ":", start(ref), "-", end(ref), sep = "")
	colnames(Y) = sampname
	for (i in 1:length(sampname)) {
		bamurl <- bamdir[i]
		what <- c("rname","pos", "mapq", "qwidth")
		if(seq=='paired-end'){
			flag <- scanBamFlag(isPaired=TRUE, isDuplicate = FALSE, isUnmappedQuery = FALSE, 
								isNotPassingQualityControls = FALSE, isFirstMateRead = TRUE)
			param <- ScanBamParam(what = what, flag = flag)
			bam <- scanBam(bamurl, param = param)[[1]]
		} else if (seq=='single-end'){
			flag <- scanBamFlag(isPaired=FALSE, isDuplicate = FALSE, isUnmappedQuery = FALSE, 
								isNotPassingQualityControls = FALSE)
			param <- ScanBamParam(what = what, 
								  flag = flag)
			bam <- scanBam(bamurl, param = param)[[1]]
		}
		message("Getting coverage for sample ",i,': ', sampname[i], "...", sep = "")
		if(any(grepl("chr", bam$rname)==TRUE)){
			bam.ref=GRanges(seqnames=bam$rname, ranges=IRanges(start=bam[["pos"]], width=bam[["qwidth"]]))
		} else{
			bam.ref=GRanges(seqnames=paste0("chr", bam$rname), ranges=IRanges(start=bam[["pos"]], width=bam[["qwidth"]]))
		}
		# remove reads with low mapping quality
		bam.ref=bam.ref[bam$mapq>=mapqthres]
		# # only keep reads from chr1 to chr22 plus chrX chrY
		# bam.ref=bam.ref[!is.na(match(seqnames(bam.ref),paste('chr',c(1:22, 'X', 'Y'),sep=''))),]
		# remove reads mapped to masked regions
		bam.ref=bam.ref[countOverlaps(bam.ref, mask.ref)==0]
		Y[, i] <- countOverlaps(ref, bam.ref)
	}
	list(Y = Y)
}