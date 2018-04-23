#!/bin/bash

#Variable for project:
project_dir=/home/paul/episodicData/novoalign

#Path to input directory
novo_pic=${project_dir}/novo_pic

#Path to Picard
pic=/usr/local/picard-tools-1.131/picard.jar

#Path to output directory
novo_rmd=${project_dir}/novo_rmd

files=(${novo_pic}/*)
for file in ${files[@]}
do
name=${file}
base=`basename ${name} _novo_sort.bam`
java -Xmx2g -jar ${pic} MarkDuplicates I= ${novo_pic}/${base}_novo_sort.bam O= ${novo_rmd}/${base}_novo_rmd.bam M= ${novo_rmd}/dupstat.txt VALIDATION_STRINGENCY=SILENT REMOVE_DUPLICATES= true
done
