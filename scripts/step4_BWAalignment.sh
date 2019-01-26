#!/bin/bash
pine=/pine/scr/r/u/rujin
kim=$pine/Kim_Navin_et_al_Cell_2018
sra_dir=$kim/sra
fastq_dir=$kim/fastq
align_dir=$kim/align

SRR=blah_blah

cd $fastq_dir

# 2) BWA
/proj/yuchaojlab/bin/bwa mem -M -t 16 \
	/proj/yuchaojlab/lib/ucsc.hg19.fasta `ls | grep "$SRR" | tr '\n' ' '` > $align_dir/"$SRR".sam

