getsampQC=function(bambedObj){
  ref <- bambedObj$ref
  bamdir <- bambedObj$bamdir
  sampname <- bambedObj$sampname
  QCmetric=matrix(ncol=8, nrow=length(sampname))
  rownames(QCmetric)=sampname
  colnames(QCmetric)=c('readlength','total','mapped','mapped_prop','non_dup','non_dup_prop','mapq20','mapq20_prop')
  for (i in 1:length(sampname)) {
    cat('Getting sample QC metric for sample',i,'\n')
    what <- c("rname","pos","strand","mapq","qwidth","flag")
    param <- ScanBamParam(what = what)
    aln <- scanBam(bamdir[i], param=param)
    aln <- aln[[1]]
    temp0=round(mean(aln$qwidth, na.rm=TRUE))
    temp1=length(aln$mapq) # total number of reads
    temp2=sum(!is.na(aln$rname)) # total mapped reads
    temp3=sum(!is.na(aln$rname) & aln$flag < 1024) # total mapped non-duplicate reads
    # https://broadinstitute.github.io/picard/explain-flags.html
    temp4=sum(aln$flag < 1024 & aln$mapq>=20,na.rm=T)
    QCmetric[i,]=c(temp0, temp1, temp2, round(temp2/temp1,3), 
                   temp3, round(temp3/temp1,3), 
                   temp4, round(temp4/temp1,3))
  }
  return(QCmetric)
}