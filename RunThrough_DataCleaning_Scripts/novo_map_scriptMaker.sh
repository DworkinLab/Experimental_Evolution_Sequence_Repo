#! /bin/bash

#Variable for project:
project_dir=/home/paul/episodicData/novoalign

#Create variable for reference genome
novo_index=${project_dir}/novo_index/dmel-all-chromosome-r5.57_2.nix

#Variable for path to Novoalign
novoalign=/usr/local/novoalign

#Path the trim outputs to be mapped
trim_dir=/home/paul/episodicData/trim_dir

#Path to output directory
novo_dir=${project_dir}/novo_dir

files=(${trim_dir}/*_R1_PE.fastq)

for file in ${files[@]}
do
name=${file}
base=`basename ${name} _R1_PE.fastq`
echo "${novoalign}/novoalign -d ${novo_index} -f ${trim_dir}/${base}_R1_PE.fastq ${trim_dir}/${base}_R2_PE.fastq -i 500,150 -o SAM > ${novo_dir}/${base}_novo.sam" > ./split_mappingScripts/${base}.sh

done
