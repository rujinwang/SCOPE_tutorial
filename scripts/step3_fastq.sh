#!/bin/bash
pine=/pine/scr/r/u/rujin
kim=$pine/Kim_Navin_et_al_Cell_2018
sra_dir=$kim/sra
fastq_dir=$kim/fastq
fastq_bin=/nas/longleaf/apps/sratoolkit/2.9.0/sra-tools/bin

sra_name=blah_blah

$fastq_bin/fastq-dump.2.9.0 -I -O $fastq_dir --split-files $sra_dir/$sra_name

