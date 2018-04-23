#!/bin/bash

#Variable for project name (title of mpileup file)
project_name=novo_episodic

#Variable for project:
project_dir=/home/paul/episodicData/novoalign

#Path to input directory
novo_final=${project_dir}/novo_GATK

#Path to output directory
novo_mpileup=${project_dir}/novo_mpileup

#Variable for reference genome (non-indexed for novoAlign)
index_dir=/home/paul/episodicData/index_dir
ref_genome=${index_dir}/dmel-all-chromosome-r5.57_2.fasta

samtools mpileup -B -Q 0 -f ${ref_genome} ${novo_final}/*.bam > ${novo_mpileup}/${project_name}.mpileup
