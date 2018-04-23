# Mapping with Novoalign and the subsequent steps to clean data:
_________________________________________________________________________________________________________________________
_________________________________________________________________________________________________________________________

## Running Novoalign mapper

Starting with sequence files that have already been inspected with md5sum and fastqc and have been trimmed using trimmomatic (See other files for running trimmomatic)

[Novoalign tutorial](http://www.novocraft.com/documentation/novoalign-2/novoalign-ngs-quick-start-tutorial/basic-short-read-mapping/)

### Need to create directory for project and to house mapping outputs

```
# Project Directory
mkdir novoalign 
```

```
# mapping outputs
mkdir novo_dir
```

### Find path to call Novoalign variable

```
#Variable for novoalign
novoalign=/usr/local/novoalign
```


### Make a scripts directory to house scripts be ran 

```
mkdir novo_scripts
```

### Novoindex reference

Need index dir for novoindex files

```
mkdir novo_index
```

The reference genome needs to be indexed for novoalign mapping (with novoindex)

Script: novo_index.sh

```
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

```


### Unzip Files

Note: Compressed read files are not supported in unlicensed versions.

 - The unlicensed version of Novoalign (used here) does not support the zipped files, so need to unzip trimmomatic outputs

```
#From trim_dir

gunzip *.gz
```

### Novoalign 

Flags: 

- d -- Full pathname of indexed reference sequence from novoindex

- f -- Files containing the read sequences to be aligned  

- o -- Specifies output report format and options (SAM)  

- i ###,## -- Sets fragment orientation and approximate fragment length for proper pairs.
    ex. -i 250 50  Defaults to paired end Illumina or Mate Pair ABI with 250bp insert and 50bp standard deviation (possible check below)
     
     - 400, 100 found based on initial mapping with novoalign first run through
     
     - using 500, 150 (below)

### Checking insert Size

Using output from Picard Sort (if available) from Bowtie2 or BWA mem previously (before removing duplicates) use Picard CollectInsertSizeMetrics.jar to get summary statistics on file for insert size and other information

If other mappers not available, map using defaults or extreme insert / SD values (i.e. 0, 500) and then run this script on .bam files (needs to be sorted with Picard)

```
java -jar /usr/local/picard-tools-1.131/picard.jar CollectInsertSizeMetrics \
    I=/home/paul/episodicData/novoalign/novo_rmd/F115ConR1_TAGCTT_novo_merge_novo_rmd.bam \
    O=/home/paul/episodicData/novoalign/novo_rmd/insert_size_metrics.txt \
    H=/home/paul/episodicData/novoalign/novo_rmd/insert_size_histogram.pdf
```

Example output:

```
MEDIAN_INSERT_SIZE|MEDIAN_ABSOLUTE_DEVIATION|MIN_INSERT_SIZE|MAX_INSERT_SIZE|MEAN_INSERT_SIZE|STANDARD_DEVIATION|READ_PAIRS|PAIR_ORIENTATION|WIDTH_OF_10_PERCENT|WIDTH_OF_20_PERCENT|WIDTH_OF_30_PERCENT|WIDTH_OF_40_PERCENT|WIDTH_OF_50_PERCENT|WIDTH_OF_60_PERCENT|WIDTH_OF_70_PERCENT|WIDTH_OF_80_PERCENT|WIDTH_OF_90_PERCENT|WIDTH_OF_99_PERCENT|SAMPLE|LIBRARY|READ_GROUP

542|94|30|28389293|542.112611|156.897954|25505882|FR|35|69|107|145|189|241|305|391|533|893  
     
For generation 115 alone: mean = 542, SD = 156
```

### Running Novoalign (File By File)

This process will run one file at a time

The script: novo_map.sh
```
#! /bin/bash

#Variable for project:
project_dir=/home/paul/episodicData/novoalign

#Create variable for reference genome
novo_index=${project_dir}/novo_index/dmel-all-chromosome-r5.57_2.nix

#Variable for path to Novoalign
novoalign=/usr/local/novoalign

#Path the trim outputs to be mapped
trim_dir=/home/paul/episodicData/trim_dir

#Path to output directory for mapped files
novo_dir=${project_dir}/novo_dir

files=(${trim_dir}/*_R1_PE.fastq)

for file in ${files[@]}
do
name=${file}
base=`basename ${name} _R1_PE.fastq`

${novoalign}/novoalign -d ${novo_index} \
    -f ${trim_dir}/${base}_R1_PE.fastq ${trim_dir}/${base}_R2_PE.fastq \ 
    -i 500,150 -o SAM > ${novo_dir}/${base}_novo.sam

done
```

This takes a long time, as the unlicensed version can only uses 1 thread (100% computer)

*** From novoalign reference manual: -c 99 Sets the number of threads to be used. On licensed versions it defaults 
to the number of CPUs as reported by sysinfo(). On free version the option is disabled ***

### Running Novoalign (Running In Parallel)

A solution to run each file seperatly in a simple splitting method

*** alternative that may be an option: add the & after the code in the for Loop (before the done) and it should push that "for" to the background and run the next file in sequence *** 


__1) Make the script to make multiple scripts__

Make dir for all output scripts: 

```
mkdir split_mappingScripts
```

Script to create many scripts

- each different file normally looped through one by one above are put into a seperate script

Script: novo_map_scriptMaker.sh
```
#! /bin/bash

#Variable for project:
project_dir=/home/paul/episodicData/novoalign

#Create variable for reference genome
novo_index=${project_dir}/novo_index/dmel-all-chromosome-r5.57_2.nix

#Variable for path to Novoalign
novoalign=/usr/local/novoalign

#Path the trim outputs to be mapped
trim_dir=/home/paul/episodicData/trim_dir

#Path to output directory
novo_dir=${project_dir}/novo_dir

files=(${trim_dir}/*_R1_PE.fastq)

for file in ${files[@]}
do
name=${file}
base=`basename ${name} _R1_PE.fastq`
echo "${novoalign}/novoalign -d ${novo_index} -f ${trim_dir}/${base}_R1_PE.fastq ${trim_dir}/${base}_R2_PE.fastq -i 500,150 -o SAM > ${novo_dir}/${base}_novo.sam" > ./split_mappingScripts/${base}.sh

done
```

__2) Create script to call all and run in parallel (use "&" which puts job in background then multiple can run at a time)__

This creates a file that has all the scripts made in step 1) in a list with ''&'' at the end to run in parrallel

One method to run only half is make files/basename based on lane (i.e 001 or 002)

Script: novo_createParallel_scipt.sh
```
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
```

__3) Change permissions and run novo_parallel_map.sh__

If needed based on the computer space available, change input parametes to run subsets on different days (one option above)

Run on screen

Screen can be named with -S (ex. screen -S IDENTIFIERTITLE)

Can save all outputs of screen using script (ex. script LOGTITLE.log) and finish script with "exit"

```
novo_parallel_map.sh
```

This should map each file seperate in unison


### Save space again with trimmed files: Rezip

Rezip files in trim_dir (saves space)

From the trim_dir:
```
gzip *.fastq
```

### Left with mapped sequences with Novoalign: can continue with cleaning the sequence data to final .bam files
____________________________________________________________________________________________________________________________________
_________________________________________________________________________________________________________________________

## Cleaning the aligned Data

__Steps used for all mapping outputs: changing parameters for input/output directories__

### Change SAM files to BAM files: novo_samTobam.sh

-- Saves space (BAM files are binary compressed versions of SAM files)

-- need bam directory for .bam files (mkdir novo_bam)

Flags:

- b -- output is .bam

- S -- input is .sam

- q 20 -- quality mapping score of 20 (standard throughout all experiments)
     

Script: novo_samTobam.sh

```
#! /bin/bash

#Variable for project:
project_dir=/home/paul/episodicData/novoalign

#Path to input directory
novo_dir=${project_dir}/novo_dir

#Path to output directory
novo_bam=${project_dir}/novo_bam

files=(${novo_dir}/*.sam)

for file in ${files[@]}
do
name=${file}
base=`basename ${name} .sam`
samtools view -b -S -q 20 ${novo_dir}/${base}.sam | samtools sort -o ${novo_bam}/${base}.bam
done
```

### Merge 
```
mkdir novo_merge
```

no flags, just merging the two lanes of the illumina sequencing run

The script: novo_merge.sh
```    
#!/bin/bash

#Variable for project:
project_dir=/home/paul/episodicData/novoalign

#Path to input directory
novo_bam=${project_dir}/novo_bam

#Path to output directory
novo_merge=${project_dir}/novo_merge


files=(${novo_bam}/*_L001_novo.bam)
for file in ${files[@]}
do
name=${file}
base=`basename ${name} _L001_novo.bam`
samtools merge ${novo_merge}/${base}_novo_merge.bam ${novo_bam}/${base}_L001_novo.bam ${novo_bam}/${base}_L002_novo.bam
done
```

### Picard Sort 

Need to sord with Picard to mark and remove duplicates (needs to be sorted with Picard in order for downstream analysis)

Need a directory for outputs, and a temporary directory for space allocation
```
mkdir novo_pic

mkdir novo_tmp
```

Flags:

 - Xmx2g -- 2 Gb of memory allocated
 
 - Djava.io.tmpdir=${tmp} -- using my temporary directory due to errors in space allocation to avoid errors while running (not necessary but useful)

 - I -- input
 
 - O -- output

 - VALIDATION_STRINGENCY=SILENT -- stops Picard from reporting every issue that would ultimately be displayed

 - SO=coordinate -- sort order based on coordinate
  
 
Script: novo_picard_sort.sh 
```
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
```

### Remove Duplicates
```
mkdir novo_rmd
```
Flags: Similar to above

- Xmx2g -- ""
        
- MarkDuplicates -- ""
        
- I -- ""
        
- O -- ""
        
- M -- creates an output file of statistics of duplicates found

- VALIDATION_STRINGENCY=SILENT -- ""

- REMOVE_DUPLICATES= true -- get rid of any found duplicated regions


Script: novo_rmd.sh
```
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
```

### More QC and make final bam files

```
mkdir novo_final
```

Flags:

- q 20 -- ""

- F 0x0004 -- remove any unmapped reads (hexidecimal value for unmapped = 0x0004)

- b -- ""

Script: novo_final_qc.sh
```
#!/bin/bash

#Variable for project:
project_dir=/home/paul/episodicData/novoalign

#Path to input directory
novo_rmd=${project_dir}/novo_rmd

#Path to output directory
novo_final=${project_dir}/novo_final


files=(${novo_rmd}/*_novo_rmd.bam)
for file in ${files[@]}
do
name=${file}
base=`basename ${name} _novo_rmd.bam`
samtools view -q 20 -F 0x0004 -b ${novo_rmd}/${base}_novo_rmd.bam > ${novo_final}/${base}_novo_final.bam
done
```

### Merge the 2 ancestors

Need to merge the base generation additionaly (two sequence runs for ancestor need to merge: MGD2 and MGD)

Script: 
```
#!/bin/bash

#Variable for project:
project_dir=/home/paul/episodicData/novoalign

#Path to final directory
novo_final=${project_dir}/novo_final

samtools merge ${novo_final}/MGD3_SO_CAGATC_novo_merge_novo_final.bam \
${novo_final}/MGD2_SO_CAGATC_novo_merge_novo_final.bam \
${novo_final}/MGD_SO_CAGATC_novo_merge_novo_final.bam

mkdir ${novo_final}/Anc_unmerged
mv ${novo_final}/ MGD2_SO_CAGATC_novo_merge_novo_final.bam ${novo_final}/Anc_unmerged
mv ${novo_final}/ MGD_SO_CAGATC_novo_merge_novo_final.bam ${novo_final}/Anc_unmerged
```

### Indel Realigner with GATK


__1) Need an unzipped copy of the reference genome__

Made a second copy and unzipping the second copy
```
gunzip dmel-all-chromosome-r5.57_2.fasta.gz
```

__2) make a gatk directory__
```
mkdir novo_GATK
```

__3) Set up reference genome__

a) Create .dict file for reference: This creates a dictionary file for the ref genome with a header but no sam records (the header is only sequence records)

script: novo_dict_index.sh
```
#! /bin/bash

pic=/usr/local/picard-tools-1.131/picard.jar
index_dir=/home/paul/episodicData/index_dir
ref_genome=${index_dir}/dmel-all-chromosome-r5.57_2.fasta

java -jar ${pic} CreateSequenceDictionary R=${ref_genome} O=${index_dir}/dmel-all-chromosome-r5.57_2.dict
```

b) Create a .fai to reference genome: 
```
samtools faidx dmel-all-chromosome-r5.57_2.fasta
```

__4) Need Read Group headers__

For GATK to run, read group headers are needed to read for Indel Realignment, __BUT__ the details are not necessary and don't need to be accurate, just need to be there.

Flags:

- RGID --Read Group Identifier; for Illumina, are composed using the flowcell + lane name and number [using Lanes L001_L002 for now]

- RGLB -- DNA Preperation Library Identifier [library1 as place holder]

- RGPL -- platform/technology used to produce the read [Illumina]

- RGPU -- Platform Unit; details on the sequencing unit (i.e run barcode) [None, used for practice]

- RGSM -- Sample [Using the basename which is each unique sequence]

Script: novo_readgroups.sh
```
#! /bin/bash

#Variable for project:
project_dir=/home/paul/episodicData/novoalign

#Path to Picard
pic=/usr/local/picard-tools-1.131/picard.jar

#Path to .bam files
novo_final=${project_dir}/novo_final

files=(${novo_final}/*.bam)
for file in ${files[@]}
do
name=${file}
base=`basename ${name} .bam`

java -jar ${pic} AddOrReplaceReadGroups I=${novo_final}/${base}.bam O=${novo_final}/${base}_RG.bam RGID=L001_L002 RGLB=library1 RGPL=illumina RGPU=None RGSM=${base}

done
```

__5) Index the read group Bam files__

Need to have the .bam files indexed prior to the indel realignment

Script: novo_indexBam.sh
```
#! /bin/bash

#Variable for project:
project_dir=/home/paul/episodicData/novoalign

#Path to .bam files
novo_final=${project_dir}/novo_final

files=(${novo_final}/*_RG.bam)
for file in ${files[@]}
do
name=${file}
base=`basename ${name} _RG.bam`
samtools index ${novo_final}/${base}_RG.bam
done
```

__6) Run GATK indel realigner__

GATK indel realigner takes two steps, 1) target the indels to be raligned (.intervals file) and 2) realign the indels (.realigned files)

Script: novo_gatk.sh
```
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
```
