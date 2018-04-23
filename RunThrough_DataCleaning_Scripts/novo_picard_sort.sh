#!/bin/bash

#Variable for project:
project_dir=/home/paul/episodicData/novoalign

#Path to input directory
novo_merge=${project_dir}/novo_merge

#Path to Picard
pic=/usr/local/picard-tools-1.131/picard.jar

#Path to output directory
novo_pic=${project_dir}/novo_pic

#Path to tmp
novo_tmp=${project_dir}/novo_tmp

files=(${novo_merge}/*.bam)
for file in ${files[@]}
do
name=${file}

base=`basename ${name} .bam`
java -Xmx2g -Djava.io.tmpdir=${novo_tmp} -jar ${pic} SortSam \
I= ${novo_merge}/${base}.bam \
O= ${novo_pic}/${base}_novo_sort.bam \
VALIDATION_STRINGENCY=SILENT SO=coordinate TMP_DIR=${novo_tmp}

done
