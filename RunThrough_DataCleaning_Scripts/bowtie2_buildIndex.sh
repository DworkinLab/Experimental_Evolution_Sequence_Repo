#! /bin/bash

project_dir1=/home/paul/episodicData
project_dir=/home/paul/episodicData/bowtie

index_dir=${project_dir1}/index_dir
ref_genome=${index_dir}/dmel-all-chromosome-r5.57_2.fasta
trim_dir=${project_dir}/trim_dir
bowtie2_dir=/usr/local/bowtie2/2.2.2


${bowtie2_dir}/bowtie2-build ${ref_genome} ${project_dir}/bowtie_indexes/dmel-all-chromosome-r5.57_2
