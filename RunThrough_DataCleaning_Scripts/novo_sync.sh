#!/bin/bash

#Variable for project name (title of mpileup file)
project_name=novo_episodic

#Variable for project:
project_dir=/home/paul/episodicData/novoalign

#Path to input/output directory
novo_mpileup=${project_dir}/novo_mpileup

#Path and variable for script from PoPoolation to create .sync files
sync=/usr/local/popoolation/mpileup2sync.jar

java -ea -Xmx7g -jar ${sync} --input ${novo_mpileup}/${project_name}.mpileup --output ${novo_mpileup}/${project_name}.sync --fastq-type sanger --min-qual 20 --threads 2
