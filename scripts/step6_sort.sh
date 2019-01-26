#!/bin/bash
pine=/pine/scr/r/u/rujin
kim=$pine/Kim_Navin_et_al_Cell_2018
sra_dir=$kim/sra
fastq_dir=$kim/fastq
align_dir=$kim/align

SRR=blah_blah

cd $align_dir

java -Xmx30G -jar /proj/yuchaojlab/bin/picard.jar SortSam \
	INPUT="$SRR".bam OUTPUT="$SRR".sorted.bam \
	SORT_ORDER=coordinate

