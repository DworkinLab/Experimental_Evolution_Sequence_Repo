#! /bin/bash

#Variable for project:
project_dir=/home/paul/episodicData/novoalign

# Variable for each script location
map_scripts=${project_dir}/novo_scripts/split_mappingScripts

#Variable for location of whole script (novo_scripts)
scripts=${project_dir}/novo_scripts


files=(${map_scripts}/*.sh)

for file in ${files[@]}
do
name=${file}
base=`basename ${name} .sh`
echo "${map_scripts}/${base}.sh &" >> ${scripts}/novo_parallel_map.sh

done
