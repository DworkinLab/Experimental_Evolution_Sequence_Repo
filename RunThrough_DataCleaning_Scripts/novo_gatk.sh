#!/bin/bash

#Variable for project name (file name)
project_name=novo_episodic

#Variable for project:
project_dir=/home/paul/episodicData/novoalign

#Path to input directory
novo_final=${project_dir}/novo_final

#Path to output directory
novo_GATK=${project_dir}/novo_GATK

#Variable for reference genome (non-zipped)
index_dir=/home/paul/episodicData/index_dir
ref_genome=${index_dir}/dmel-all-chromosome-r5.57_2.fasta

#Path to GATK
gatk=/usr/local/gatk/GenomeAnalysisTK.jar

files=(${novo_final}/*_RG.bam)
for file in ${files[@]}
do
name=${file}
base=`basename ${name} _RG.bam`

java -Xmx8g -jar ${gatk} -I ${novo_final}/${base}_RG.bam -R ${ref_genome} -T RealignerTargetCreator -o ${novo_GATK}/${base}.intervals

java -Xmx8g -jar ${gatk} -I ${novo_final}/${base}_RG.bam -R ${ref_genome} -T IndelRealigner -targetIntervals ${novo_GATK}/${base}.intervals -o ${novo_GATK}/${base}_realigned.bam

done
