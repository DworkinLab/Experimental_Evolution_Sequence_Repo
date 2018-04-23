# Full run through of all sequence data - Trimming to Final .bam file format:
____________________________________________________________________________________________________

## Sections of file

**Sections 1:** Quality checks on data and trimming

**Section 2:** Mapping with BWA mem, Bowtie2 and Novoalign (3 mappers to ensure no false positive positions)

**Section 3:** Methods for cleaning data for one mapper (Novoalign) that mirrors methods followed for BWA mem and Bowtie2


____________________________________________________________________________________________________
____________________________________________________________________________________________________
____________________________________________________________________________________________________
____________________________________________________________________________________________________
____________________________________________________________________________________________________

## **Sections 1:** Quality checks on data and trimming

### md5sum:

md5sum all raw files: changes depending on the file name

```
md5sum - c md5.txt
```

### Fastqc; run as a quality control and view

```
fastqc -o ${output} ${input}/*.fastq.gz
```

### Index sequence for mapping:

Version 5.57 used as the reference sequence for *Drosophila melanogaster* obtained from [flybase.org](http://flybase.org/)
```
curl -O ftp://ftp.flybase.net/genomes/Drosophila_melanogaster/dmel_r5.57_FB2014_03/fasta/dmel-all-chromosome-r5.57.fasta.gz

bwa index dmel-all-chromosome-r5.57.fasta.gz
```
__________________________________________________________________________________

### Trimmomatic

Flags:

    -phred33 = may not need to be specified
    -trimlog = log of trim outputs
    -IlluminaClip = adapter removal (${adapter})
    -LEADING & TRAILING = 3; removal at start end end if below quality
    -MINLEN = minimum length of 36
    -MAXINFO = adaptive quality (balance b/w length and quality) = 0.5

Script:  
```
#! /bin/bash

#Variables for script:
project_name=episodic_data
project_dir=/home/paul/episodicData
raw_dir=${project_dir}/raw_dir

trimmomatic=/usr/local/trimmomatic
trim=${trimmomatic}/trimmomatic-0.33.jar

adapt_path=/usr/local/trimmomatic/adapters
adapter=${adapt_path}/TruSeq3-PE.fa:2:30:10

trim_dir=${project_dir}/trim_dir

#Loop to run on all raw (pairded) data files)
files=(${raw_dir}/*_R1_001.fastq.gz)
for file in ${files[@]} 
do
name=${file}
base=`basename ${name} _R1_001.fastq.gz`
java -jar ${trim} PE -phred33 -trimlog ${trim_dir}/trimlog.txt ${raw_dir}/${base}_R1_001.fastq.gz ${raw_dir}/${base}_R2_001.fastq.gz ${trim_dir}/${base}_R1_PE.fastq.gz ${trim_dir}/${base}_R1_SE.fastq.gz ${trim_dir}/${base}_R2_PE.fastq.gz ${trim_dir}/${base}_R2_SE.fastq.gz ILLUMINACLIP:${adapter} LEADING:3 TRAILING:3 MAXINFO:40:0.5 MINLEN:36

done
```

Have trimmed raw data files that can be read into different mappers.

____________________________________________________________________________________________________
____________________________________________________________________________________________________
____________________________________________________________________________________________________
____________________________________________________________________________________________________
____________________________________________________________________________________________________

## **Section 2:** Mapping with BWA mem, Bowtie2 and Novoalign (3 mappers to ensure no false positive positions)
_________________________________________________________________________________________________________________________

### BWA mapping

__1) Running bwa mem mapping__
Flags:

      -t 8 = number of processors
      
      -M = Mark shorter split hits as secondary (for Picard compatibility)

Script:  
```
#!/bin/bash

project_name=episodic_data
project_dir=/home/paul/episodicData
bwa_path=/usr/local/bwa/0.7.8
index_dir=${project_dir}/index_dir
ref_genome=${index_dir}/dmel-all-chromosome-r5.57.fasta.gz
trim_dir=${project_dir}/trim_dir
sam_dir=${project_dir}/sam_dir

cd ${bwa_path}
files=(${trim_dir}/*_R1_PE.fastq.gz)

for file in ${files[@]}
do
name=${file}
base=`basename ${name} _R1_PE.fastq.gz`
bwa mem -t 8 -M ${ref_genome} ${trim_dir}/${base}_R1_PE.fastq.gz ${trim_dir}/${base}_R2_PE.fastq.gz > ${sam_dir}/${base}_aligned_pe.SAM
done
```


### Bowtie2 mapping: Steps needed to index reference

__1) Build proper index:__

__1.1) Need unzipped index sequence__

```
gunzip /home/paul/episodicData/index_dir/dmel-all-chromosome-r5.57_2.fasta.gz
```

__1.2) build bowtie index__

```
#! /bin/bash

project_dir1=/home/paul/episodicData
project_dir=/home/paul/episodicData/bowtie

index_dir=${project_dir1}/index_dir
ref_genome=${index_dir}/dmel-all-chromosome-r5.57_2.fasta
trim_dir=${project_dir}/trim_dir
bowtie2_dir=/usr/local/bowtie2/2.2.2


${bowtie2_dir}/bowtie2-build ${ref_genome} ${project_dir}/bowtie_indexes/dmel-all-chromosome-r5.57_2
```
__2) Run Bowtie2 mapping__

Flags:
 - x = indexed reference
 
 - 1 = forward end
 
 - 2 = reverse end
 
 - S = same output
```
#! /bin/bash

project_name=episodic_data_bowtie
project_dir1=/home/paul/episodicData
project_dir=/home/paul/episodicData/bowtie
index_dir=${project_dir1}/index_dir
ref_genome=${index_dir}/dmel-all-chromosome-r5.57.fasta.gz
ref_genome_base=${project_dir}/bowtie_indexes/dmel-all-chromosome-r5.57_2
trim_dir=${project_dir1}/trim_dir
bowtie2_dir=/usr/local/bowtie2/2.2.2
sam_dir=${project_dir}/sam_dir

files=(${trim_dir}/*_R1_PE.fastq.gz)
for file in ${files[@]}
do
name=${file}
base=`basename ${name} _R1_PE.fastq.gz`
${bowtie2_dir}/bowtie2 -x ${ref_genome_base} -1 ${trim_dir}/${base}_R1_PE.fastq.gz -2 ${trim_dir}/${base}_R2_PE.fastq.gz -S ${sam_dir}/${base}_bowtie_pe.sam
done
```

_________________________________________________________________________________________________________________________
### Mapping with Novoalign and the subsequent steps to clean data: Many steps for set up

**Running Novoalign mapper**

[Novoalign tutorial](http://www.novocraft.com/documentation/novoalign-2/novoalign-ngs-quick-start-tutorial/basic-short-read-mapping/)


__1) Novoindex reference__
The reference genome needs to be indexed for novoalign mapping (with novoindex)

*Script: novo_index.sh*

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
__2) Unzip trimmed Files__

Note: Compressed read files are not supported in unlicensed versions.

 - The unlicensed version of Novoalign (used here) does not support the zipped files, so need to unzip trimmomatic outputs

```
gunzip *.gz
```


__3) Novoalign Flags:__

- d -- Full pathname of indexed reference sequence from novoindex

- f -- Files containing the read sequences to be aligned  

- o -- Specifies output report format and options (SAM)  

- i ###,## -- Sets fragment orientation and approximate fragment length for proper pairs.
    ex. -i 250 50  Defaults to paired end Illumina or Mate Pair ABI with 250bp insert and 50bp standard deviation (possible check below)
     
     - 400, 100 found based on initial mapping with novoalign first run through
     
     - using 500, 150 (below)

__4) Checking insert Size__

Using output from Picard Sort (if available) from Bowtie2 or BWA mem previously (before removing duplicates) use Picard CollectInsertSizeMetrics.jar to get summary statistics on file for insert size and other information

If other mappers not available, map using defaults or extreme insert / SD values (i.e. 0, 500) and then run this script on .bam files (needs to be sorted with Picard)

```
java -jar /usr/local/picard-tools-1.131/picard.jar CollectInsertSizeMetrics \
    I=/home/paul/episodicData/novoalign/novo_rmd/F115ConR1_TAGCTT_novo_merge_novo_rmd.bam \
    O=/home/paul/episodicData/novoalign/novo_rmd/insert_size_metrics.txt \
    H=/home/paul/episodicData/novoalign/novo_rmd/insert_size_histogram.pdf
```

__5) Running Novoalign (Running In Parallel)__
 

__5.1) Make the script to make multiple scripts__

Script to create many scripts (that run in parallel)

*Script: novo_map_scriptMaker.sh*
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

__5.2) Create script to call all and run in parallel (use "&" which puts job in background then multiple can run at a time)__

This creates a file that has all the scripts made in step 5.1) in a list with ''&'' at the end to run in parrallel

*Script: novo_createParallel_scipt.sh*
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

__5.3) Run novo_parallel_map.sh__

```
novo_parallel_map.sh
```

__6) Save space again with trimmed files: Rezip__

Rezip files in trim_dir (saves space)

```
gzip *.fastq
```

____________________________________________________________________________________________________
____________________________________________________________________________________________________
____________________________________________________________________________________________________
____________________________________________________________________________________________________
____________________________________________________________________________________________________

## **Section 3:** Methods for cleaning data for one mapper (Novoalign) that mirrors methods followed for BWA mem and Bowtie2

Have three sets of mapped data: the following methods (for Novoalign) is an example of next steps, but are done identically from here to final bams.

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

### Indel Realigner with GATK (6 steps)


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
