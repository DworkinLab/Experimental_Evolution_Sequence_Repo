#! /bin/bash

# Variable for project:
project_dir=/home/paul/episodicData/novoalign

# Path to input directory
input=${project_dir}/novo_GATK

# Path to output novoalign pileup files
output=${project_dir}/novo_pileup

index_dir=/home/paul/episodicData/index_dir
ref_genome=${index_dir}/dmel-all-chromosome-r5.57_2.fasta

files=(${input}/*_merge_novo_final_realigned.bam)

for file in ${files[@]}

do

name=${file}

base=`basename ${name} _merge_novo_final_realigned.bam`

samtools mpileup -B -Q 0 -f ${ref_genome} ${input}/${base}_merge_novo_final_realigned.bam > ${output}/${base}.pileup

done
