#! /bin/bash

#Create variable for location of reference genome (fasta vs. fasta.gz?)
ref_genome=/home/paul/episodicData/index_dir/dmel-all-chromosome-r5.57_2.fasta

#Variable for project
project_dir=/home/paul/episodicData/novoalign

#Variable for novoalign
novoalign=/usr/local/novoalign

#Variable for output directory
novo_index=${project_dir}/novo_index

#Index the reference with novoindex

${novoalign}/novoindex ${novo_index}/dmel-all-chromosome-r5.57_2.nix  ${ref_genome}
