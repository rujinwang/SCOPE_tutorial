#!/bin/bash
pine=/pine/scr/r/u/rujin
kim=$pine/Kim_Navin_et_al_Cell_2018
sra_dir=$kim/sra
fastq_dir=$kim/fastq
align_dir=$kim/align

SRR=blah_blah

cd $align_dir

java -Xmx40G -jar /proj/yuchaojlab/bin/picard.jar AddOrReplaceReadGroups \
	I="$SRR".sorted.bam O="$SRR".sorted.rg.bam RGID="$SRR" \
	RGLB=Chung_Et_Al RGPL=ILLUMINA RGPU=machine RGSM="$SRR"
	
/proj/yuchaojlab/bin/samtools index "$SRR".sorted.rg.bam
