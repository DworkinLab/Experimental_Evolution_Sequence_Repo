#! /bin/bash

#Variable for project name (title of mpileup file)
project_name=novo_episodic

#Variable for project:
project_dir=/home/paul/episodicData/novoalign

#Path to .sync files
SyncFiles=${project_dir}/novo_mpileup

#Output dir:
mkdir ${project_dir}/ChromoSubsets
subsets=${project_dir}/ChromoSubsets

# Need to copy three R scripts and add to a new directory (i.e. novo_Rscripts)
Rscripts=${project_dir}/novo_Rscripts


# The seperated .sync files
sync[0]=${SyncFiles}/novo_episodic_3R.sync
sync[1]=${SyncFiles}/novo_episodic_2R.sync
sync[2]=${SyncFiles}/novo_episodic_3L.sync
sync[3]=${SyncFiles}/novo_episodic_2L.sync
sync[4]=${SyncFiles}/novo_episodic_X.sync 

#Removing 4 b/c complete and short Run on its own (Does not need to be split)

for file in ${sync[@]}
do
name=${file}
base=`basename ${name} .sync`

mkdir ${subsets}/${base}_dir
basedir=${subsets}/${base}_dir

length=($(wc -l ${SyncFiles}/${base}.sync))
echo ${length}

#Split length into 11 segements (11th == length) (can extend this if to large)
cut=$((${length}/11))
cut_2=$((${cut}*2))
cut_3=$((${cut}*3))
cut_4=$((${cut}*4))
cut_5=$((${cut}*5))
cut_6=$((${cut}*6))
cut_7=$((${cut}*7))
cut_8=$((${cut}*8))
cut_9=$((${cut}*9))
cut_10=$((${cut}*10))

###

sed -n " 1, ${cut} p" ${SyncFiles}/${base}.sync > ${basedir}/${base}_1.sync
sed -n " $((${cut} + 1)), ${cut_2} p" ${SyncFiles}/${base}.sync > ${basedir}/${base}_2.sync
sed -n " $((${cut_2} + 1)), ${cut_3} p" ${SyncFiles}/${base}.sync > ${basedir}/${base}_3.sync
sed -n " $((${cut_3} + 1)), ${cut_4} p" ${SyncFiles}/${base}.sync > ${basedir}/${base}_4.sync
sed -n " $((${cut_4} + 1)), ${cut_5} p" ${SyncFiles}/${base}.sync > ${basedir}/${base}_5.sync
sed -n " $((${cut_5} + 1)), ${cut_6} p" ${SyncFiles}/${base}.sync > ${basedir}/${base}_6.sync
sed -n " $((${cut_6} + 1)), ${cut_7} p" ${SyncFiles}/${base}.sync > ${basedir}/${base}_7.sync
sed -n " $((${cut_7} + 1)), ${cut_8} p" ${SyncFiles}/${base}.sync > ${basedir}/${base}_8.sync
sed -n " $((${cut_8} + 1)), ${cut_9} p" ${SyncFiles}/${base}.sync > ${basedir}/${base}_9.sync
sed -n " $((${cut_9} + 1)), ${cut_10} p" ${SyncFiles}/${base}.sync > ${basedir}/${base}_10.sync
sed -n " $((${cut_10} + 1)), ${length} p" ${SyncFiles}/${base}.sync > ${basedir}/${base}_11.sync
done


mkdir ${subsets}/novo_episodic_4_dir
cp ${SyncFiles}/novo_episodic_4.sync ${subsets}/novo_episodic_4_dir

echo 'Done Splitting Files'

# Should now have 11 different .sync files to work with


###################

# Running Rscript below to go through each subset file created and convert to counts (see Sync_to_counts.R)

echo 'Running R script Sync To Counts'

Rscript ${Rscripts}/Sync_to_counts.R ${subsets}

echo ' Done R script Sync to Counts'

# can remove the sync files? (need a loop to enter each ${subsets} dir to remove *.sync

sync[0]=${SyncFiles}/novo_episodic_3R.sync
sync[1]=${SyncFiles}/novo_episodic_2R.sync
sync[2]=${SyncFiles}/novo_episodic_3L.sync
sync[3]=${SyncFiles}/novo_episodic_2L.sync
sync[4]=${SyncFiles}/novo_episodic_X.sync 
sync[5]=${SyncFiles}/novo_episodic_4.sync 

for file in ${sync[@]}
	do
	name=${file}
	base=`basename ${name} .sync`
	basedir=${subsets}/${base}_dir
	rm -f ${basedir}/*.sync
done

echo 'Removed Sync Files'

########################

# Run Model On each seperate file (creating a coeffs file)

echo 'Rscript Model'

#Rscript ${Rscripts}/Counts_to_model.R ${subsets}
# Left with the ${subsets} directory full of .csv files and .coeffs.csv files
# Note: this step takes the longest time, possibly search for method to parallelize this script (per chromosome, per file etc..)
# This R script is not the same as below (needs to set large overall directory), I don't recombend using this.
#

# Rscript ${Rscripts}/Counts_to_model.R ${subsets}/${basedir}

# Loop will not work because I need to call the "main" dir holding all the directories: Write seperate script for individual directory. (removes direcory option, args[1] == setwd() and ready to go.)

sync[0]=${SyncFiles}/novo_episodic_3R.sync
sync[1]=${SyncFiles}/novo_episodic_2R.sync
sync[2]=${SyncFiles}/novo_episodic_3L.sync
sync[3]=${SyncFiles}/novo_episodic_2L.sync
sync[4]=${SyncFiles}/novo_episodic_X.sync 
sync[5]=${SyncFiles}/novo_episodic_4.sync 

for file in ${sync[@]}
	do
	(name=${file}
	base=`basename ${name} .sync`
	basedir=${subsets}/${base}_dir
	Rscript ${Rscripts}/Counts_to_model_2.R ${basedir}) &	
done
wait

echo 'Done Model'
#########################

#Combine Coeffs file into one large chromosome file

#mkdir in main project location for combined coeffs 

# Takes into assumption the files are named "episodic_data_2L_11.sync.csv.coeffs.csv" to create a final .csv named "episodic_data_2L_chromo.csv" (removes all past the last _)

mkdir ${project_dir}/novo_coeffs
coeff_dir=${project_dir}/novo_coeffs

Rscript ${Rscripts}/Combine_chromo.R ${subsets} ${coeff_dir}

########################

#Remove all .csv intermediates from random files:

for file in ${sync[@]}
	do
	name=${file}
	base=`basename ${name} .sync`
	basedir=${subsets}/${base}_dir
	rm -f ${basedir}/*.csv
	rmdir ${base}_dir
done

rmdir ${subsets}

#left with the final output files (coeffs) for novoalign files

echo 'DONE'
