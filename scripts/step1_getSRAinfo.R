rm(list = ls())

curr_dir = "/pine/scr/r/u/rujin"
out_dir = file.path(curr_dir,"Kim_Navin_et_al_Cell_2018")
setwd(out_dir)

# get sra session info using R
source('http://bioconductor.org/biocLite.R')
biocLite('SRAdb')
library(SRAdb)
srafile = getSRAdbFile()
library(DBI)
con = dbConnect(RSQLite::SQLite(), srafile)

sra.files=listSRAfile('SRP114962',con)
write.table(sra.files,file='breast_cancer.sra.files.txt',sep='\t',row.names=F,col.names=T,quote=F)

table.kim = read.table("SraRunTable.txt", header = TRUE, stringsAsFactors = FALSE, sep = '\t')
scDNA.SRA = table.kim[grep("CNV", table.kim$isolate),]

q("no")