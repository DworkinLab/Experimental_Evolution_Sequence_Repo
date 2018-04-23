#!/bin/bash

#Variable for project:
project_dir=/home/paul/episodicData/novoalign

#Path to input directory
novo_rmd=${project_dir}/novo_rmd

#Path to output directory
novo_final=${project_dir}/novo_final

files=(${novo_rmd}/*_novo_rmd.bam)
for file in ${files[@]}
do
name=${file}
base=`basename ${name} _novo_rmd.bam`
samtools view -q 20 -F 0x0004 -b ${novo_rmd}/${base}_novo_rmd.bam > ${novo_final}/${base}_novo_final.bam
done
