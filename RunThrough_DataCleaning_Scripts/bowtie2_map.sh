#! /bin/bash

project_name=episodic_data_bowtie
project_dir1=/home/paul/episodicData
project_dir=/home/paul/episodicData/bowtie
index_dir=${project_dir1}/index_dir
ref_genome=${index_dir}/dmel-all-chromosome-r5.57.fasta.gz
ref_genome_base=${project_dir}/bowtie_indexes/dmel-all-chromosome-r5.57_2
trim_dir=${project_dir1}/trim_dir
bowtie2_dir=/usr/local/bowtie2/2.2.2
sam_dir=${project_dir}/sam_dir

files=(${trim_dir}/*_R1_PE.fastq.gz)
for file in ${files[@]}
do
name=${file}
base=`basename ${name} _R1_PE.fastq.gz`
${bowtie2_dir}/bowtie2 -x ${ref_genome_base} -1 ${trim_dir}/${base}_R1_PE.fastq.gz -2 ${trim_dir}/${base}_R2_PE.fastq.gz -S ${sam_dir}/${base}_bowtie_pe.sam
done
