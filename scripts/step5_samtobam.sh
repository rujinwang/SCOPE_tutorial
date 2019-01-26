#!/bin/bash
pine=/pine/scr/r/u/rujin
kim=$pine/Kim_Navin_et_al_Cell_2018
sra_dir=$kim/sra
fastq_dir=$kim/fastq
align_dir=$kim/align

SRR=blah_blah

cd $align_dir
/proj/yuchaojlab/bin/samtools view -bS "$SRR".sam > "$SRR".bam
