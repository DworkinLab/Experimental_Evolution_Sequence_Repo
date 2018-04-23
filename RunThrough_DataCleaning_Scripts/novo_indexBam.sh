#! /bin/bash

#Variable for project:
project_dir=/home/paul/episodicData/novoalign

#Path to .bam files
novo_final=${project_dir}/novo_final

files=(${novo_final}/*_RG.bam)
for file in ${files[@]}
do
name=${file}
base=`basename ${name} _RG.bam`
samtools index ${novo_final}/${base}_RG.bam
done
