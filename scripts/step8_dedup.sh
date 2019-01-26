#!/bin/bash
pine=/pine/scr/r/u/rujin
kim=$pine/Kim_Navin_et_al_Cell_2018
sra_dir=$kim/sra
fastq_dir=$kim/fastq
align_dir=$kim/align

SRR=blah_blah

cd $align_dir

java -Xmx40G -jar /proj/yuchaojlab/bin/picard.jar MarkDuplicates \
	REMOVE_DUPLICATES=true \
	I="$SRR".sorted.rg.bam O="$SRR".sorted.rg.dedup.bam \
	METRICS_FILE="$SRR".sorted.rg.dedup.metrics.txt \
	PROGRAM_RECORD_ID= MarkDuplicates PROGRAM_GROUP_VERSION=null \
	PROGRAM_GROUP_NAME=MarkDuplicates
	
java -jar /proj/yuchaojlab/bin/picard.jar BuildBamIndex I="$SRR".sorted.rg.dedup.bam
