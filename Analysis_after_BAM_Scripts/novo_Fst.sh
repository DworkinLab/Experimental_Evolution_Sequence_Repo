#!/bin/bash

#Variable for project name (title of mpileup file)
project_name=novo_episodic

#Variable for project:
project_dir=/home/paul/episodicData/novoalign

#Path to .sync files
novo_mpileup=${project_dir}/novo_mpileup

#Path and variable for script from PoPoolation to create .sync files
fst=/usr/local/popoolation/fst-sliding.pl

mkdir ${novo_mpileup}/novo_fst
novo_fst=${novo_mpileup}/novo_fst

perl ${fst} --input ${novo_mpileup}/novo_episodic_main.sync --output ${novo_fst}/novo_episodic_main.fst --min-count 6 --min-coverage 10 --max-coverage 250 --min-covered-fraction 1 --window-size 500 --step-size 500 --pool-size 120
