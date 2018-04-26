#! /bin/bash

### Set all variables (need to make an output directory):

#Variable for project name (title of mpileup file)
project_name=novo_episodic

#Variable for project:
project_dir=/home/paul/episodicData/novoalign

#Path to .sync files
SyncFiles=${project_dir}/novo_mpileup
	
mkdir ${SyncFiles}/splitsync_dir
splitSync=${SyncFiles}/splitsync_dir

#Output dir:
poolSeq=${project_dir}/novo_PoolSeq

# Need to copy three R scripts and add to a new directory (i.e. novo_Rscripts)
Rscripts=${project_dir}/novo_Rscripts
	
# The seperated .sync files
sync[0]=${SyncFiles}/novo_episodic_3R.sync
sync[1]=${SyncFiles}/novo_episodic_2R.sync
sync[2]=${SyncFiles}/novo_episodic_3L.sync
sync[3]=${SyncFiles}/novo_episodic_2L.sync
sync[4]=${SyncFiles}/novo_episodic_X.sync 
sync[5]=${SyncFiles}/novo_episodic_4.sync 

##-----------------------------------------------##

### Split into treatment vs. control

for file in ${sync[@]}
	do
	name=${file}
	base=`basename ${name} .sync`
	
	cat ${SyncFiles}/${base}.sync | awk '{print $1,$2,$3,$6,$7,$10, $11, $14, $15, $16, $16}' > ${splitSync}/${base}_Sel.sync
	
	cat ${SyncFiles}/${base}.sync | awk '{print $1,$2,$3,$4,$5,$8, $9, $12, $13, $16, $16}' > ${splitSync}/${base}_Con.sync

done


##------------------------------------------------##

### Split the sync files into many sized files (12):

files=(${splitSync}/*.sync)

for file in ${files[@]}
	do
	name=${file}
	base=`basename ${name} .sync`
	
	mkdir ${splitSync}/${base}_Split
	split_sync=${splitSync}/${base}_Split
	
	length=($(wc -l ${splitSync}/${base}.sync))
	#echo ${length}
		
	#Split length into 12 segements (12th == length) (can extend this if to large)
	cut=$((${length}/12))
	cut_2=$((${cut}*2))
	cut_3=$((${cut}*3))
	cut_4=$((${cut}*4))
	cut_5=$((${cut}*5))
	cut_6=$((${cut}*6))
	cut_7=$((${cut}*7))
	cut_8=$((${cut}*8))
	cut_9=$((${cut}*9))
	cut_10=$((${cut}*10))
	cut_11=$((${cut}*11))
	
	sed -n " 1, ${cut} p"  ${splitSync}/${base}.sync > ${split_sync}/${base}_1.sync

	sed -n " $((${cut} + 1)), ${cut_2} p"  ${splitSync}/${base}.sync >  ${split_sync}/${base}_2.sync

	sed -n " $((${cut_2} + 1)), ${cut_3} p"  ${splitSync}/${base}.sync > ${split_sync}/${base}_3.sync
	
	sed -n " $((${cut_3} + 1)), ${cut_4} p"  ${splitSync}/${base}.sync > ${split_sync}/${base}_4.sync

	sed -n " $((${cut_4} + 1)), ${cut_5} p"  ${splitSync}/${base}.sync > ${split_sync}/${base}_5.sync

	sed -n " $((${cut_5} + 1)), ${cut_6} p"  ${splitSync}/${base}.sync > ${split_sync}/${base}_6.sync

	sed -n " $((${cut_6} + 1)), ${cut_7} p"  ${splitSync}/${base}.sync > ${split_sync}/${base}_7.sync
	
	sed -n " $((${cut_7} + 1)), ${cut_8} p"  ${splitSync}/${base}.sync > ${split_sync}/${base}_8.sync
	
	sed -n " $((${cut_8} + 1)), ${cut_9} p"  ${splitSync}/${base}.sync > ${split_sync}/${base}_9.sync
	
	sed -n " $((${cut_9} + 1)), ${cut_10} p"  ${splitSync}/${base}.sync > ${split_sync}/${base}_10.sync
	
	sed -n " $((${cut_10} + 1)), ${cut_11} p"  ${splitSync}/${base}.sync > ${split_sync}/${base}_11.sync
	
	sed -n " $((${cut_11} + 1)), ${length} p"  ${splitSync}/${base}.sync > ${split_sync}/${base}_12.sync
	
	syncs=(${split_sync}/*.sync)
 	
	for file in ${syncs[@]}
	  	do
	  	(Chromo=$(cat ${file} | awk '{print $1; exit}')
	  	Rscript ${Rscripts}/poolSeq_selectionCoeff.R ${file} ${Chromo} ${split_sync}) &
	done 
	wait
	rm -f ${split_sync}/*.sync
done
wait

Rscript ${Rscripts}/combinePoolseqCSV.R ${splitSync}

##------------------------------------------------##
