#! /bin/bash

#Variable for project:
project_dir=/home/paul/episodicData/novoalign

#Path to input directory
novo_dir=${project_dir}/novo_dir

#Path to output directory
novo_bam=${project_dir}/novo_bam

files=(${novo_dir}/*.sam)

for file in ${files[@]}
do
name=${file}
base=`basename ${name} .sam`
samtools view -b -S -q 20 ${novo_dir}/${base}.sam | samtools sort -o ${novo_bam}/${base}.bam
done
