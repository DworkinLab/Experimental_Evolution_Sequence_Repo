# Full run through of all sequence data - Trimming to Final .bam file format and creating .sync file:
____________________________________________________________________________________________________

## Sections of file

**Sections 1:** Quality checks on data and trimming

**Section 2:** Mapping with BWA mem, Bowtie2 and Novoalign (3 mappers to ensure no false positive positions)

**Section 3:** Methods for cleaning data for one mapper (Novoalign) that mirrors methods followed for BWA mem and Bowtie2

**Section 4:** Create mpileup and .sync file for analysis

____________________________________________________________________________________________________
____________________________________________________________________________________________________
____________________________________________________________________________________________________
____________________________________________________________________________________________________
____________________________________________________________________________________________________

## **Sections 1:** Quality checks on data and trimming
____________________________________________________________________________________________________

### md5sum:

md5sum all raw files: changes depending on the file name

Check to ensure the files match (are correct) after moving around between people and computers

```
md5sum - c md5.txt
```
____________________________________________________________________________________________________

### Fastqc; run as a quality control

```
fastqc -o ${output} ${input}/*.fastq.gz
```
____________________________________________________________________________________________________

### Index sequence for mapping:

Version 5.57 used as the reference sequence for *Drosophila melanogaster* obtained from [flybase.org](http://flybase.org/)
```
curl -O ftp://ftp.flybase.net/genomes/Drosophila_melanogaster/dmel_r5.57_FB2014_03/fasta/dmel-all-chromosome-r5.57.fasta.gz

bwa index dmel-all-chromosome-r5.57.fasta.gz
```
__________________________________________________________________________________

### Trimmomatic

Trim raw reads for minimum length, quality, adapters, etc.

Flags:

    -phred33 = may not need to be specified
    -trimlog = log of trim outputs
    -IlluminaClip = adapter removal (${adapter})
    -LEADING & TRAILING = 3; removal at start end end if below quality
    -MINLEN = minimum length of 36
    -MAXINFO = adaptive quality (balance b/w length and quality) = 0.5

*Script: [trim.sh](https://github.com/PaulKnoops/Experimental_Evolution_Sequence_Repo/blob/master/RunThrough_DataCleaning_Scripts/trim.sh)*

Ex.
```
java -jar trimmomatic-0.33.jar PE -phred33 -trimlog ${trim_dir}/trimlog.txt ${raw_dir}/${base}_R1_001.fastq.gz ${raw_dir}/${base}_R2_001.fastq.gz ${trim_dir}/${base}_R1_PE.fastq.gz ${trim_dir}/${base}_R1_SE.fastq.gz ${trim_dir}/${base}_R2_PE.fastq.gz ${trim_dir}/${base}_R2_SE.fastq.gz ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 MAXINFO:40:0.5 MINLEN:36
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

*Script: [bwa_map.sh](https://github.com/PaulKnoops/Experimental_Evolution_Sequence_Repo/blob/master/RunThrough_DataCleaning_Scripts/bwa_map.sh)* 

Ex.
```
bwa mem -t 8 -M dmel-all-chromosome-r5.57.fasta.gz ${trim_dir}/${base}_R1_PE.fastq.gz ${trim_dir}/${base}_R2_PE.fastq.gz > ${ouput}/${base}_aligned_pe.SAM
```
____________________________________________________________________________________________________

### Bowtie2 mapping: 

__1) Build proper Bowtie2 index:__

__1.1) Need unzipped index sequence__

```
gunzip /home/paul/episodicData/index_dir/dmel-all-chromosome-r5.57_2.fasta.gz
```

__1.2) build bowtie index__

*Script: [bowtie2_buildIndex.sh](https://github.com/PaulKnoops/Experimental_Evolution_Sequence_Repo/blob/master/RunThrough_DataCleaning_Scripts/bowtie2_buildIndex.sh)*

Ex.
```
${bowtie2_dir}/bowtie2-build ${ref_genome_dir}/dmel-all-chromosome-r5.57_2.fasta ${ref_genome_dir}/bowtie_indexes/dmel-all-chromosome-r5.57_2
```
__2) Run Bowtie2 mapping__

Flags:

    - x = indexed reference
    - 1 = forward end
    - 2 = reverse end
    - S = sam output

*Script: [bowtie2_map.sh](https://github.com/PaulKnoops/Experimental_Evolution_Sequence_Repo/blob/master/RunThrough_DataCleaning_Scripts/bowtie2_map.sh)*

Ex.
```
${bowtie2_dir}/bowtie2 -x ${bowtie2-build_dir} -1 ${trim_dir}/${base}_R1_PE.fastq.gz -2 ${trim_dir}/${base}_R2_PE.fastq.gz -S ${output}/${base}_bowtie_pe.sam
```

_________________________________________________________________________________________________________________________

### Mapping with Novoalign and the subsequent steps to clean data: Many steps for set up

**Running Novoalign mapper example**: [Novoalign tutorial](http://www.novocraft.com/documentation/novoalign-2/novoalign-ngs-quick-start-tutorial/basic-short-read-mapping/)

__1) Novoindex reference (index reference genome)__

The reference genome needs to be indexed for novoalign mapping (with novoindex)

*Script: [novo_index.sh](https://github.com/PaulKnoops/Experimental_Evolution_Sequence_Repo/blob/master/RunThrough_DataCleaning_Scripts/novo_index.sh)*

Ex.
```
${novoalign_die}/novoindex ${novo_index}/dmel-all-chromosome-r5.57_2.nix  ${reference_genome}
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
    ex. -i 250 50  Defaults to paired end Illumina or Mate Pair ABI with 250bp insert and 50bp standard deviation (400, 100 found based on initial mapping with novoalign first run through)

__3.1) Checking insert Size__

Using output from Picard Sort (if available) from Bowtie2 or BWA mem previously (before removing duplicates) use Picard CollectInsertSizeMetrics.jar to get summary statistics on file for insert size and other information

If other mappers not available, map using defaults or extreme insert / SD values (i.e. 0, 500) and then run this script on .bam files (needs to be sorted with Picard)

```
java -jar /usr/local/picard-tools-1.131/picard.jar CollectInsertSizeMetrics \
    I=/home/paul/episodicData/novoalign/novo_rmd/F115ConR1_TAGCTT_novo_merge_novo_rmd.bam \
    O=/home/paul/episodicData/novoalign/novo_rmd/insert_size_metrics.txt \
    H=/home/paul/episodicData/novoalign/novo_rmd/insert_size_histogram.pdf
```

__3.2) Basis of novoalign script:__

Running novoalign for paired files run like below

Ex.
```
${novoalign_dir}/novoalign -d ${novo_index} \
    -f ${trim_dir}/${base}_R1_PE.fastq ${trim_dir}/${base}_R2_PE.fastq \ 
    -i 500,150 -o SAM > ${novo_dir}/${base}_novo.sam
```

But due to the high computational demand and limitations with free version of novoalign, running this in parallel more beneficial

__4) Running Novoalign (Running In Parallel)__
 
__4.1) Make the script to make multiple scripts__

Script to create many scripts (that run in parallel) (echo each line for each file)

*Script: [novo_map_scriptMaker.sh](https://github.com/PaulKnoops/Experimental_Evolution_Sequence_Repo/blob/master/RunThrough_DataCleaning_Scripts/novo_map_scriptMaker.sh)*

Ex.
```
echo "${novoalign_dir}/novoalign -d ${novo_index} -f ${trim_dir}/${base}_R1_PE.fastq ${trim_dir}/${base}_R2_PE.fastq -i 500,150 -o SAM > ${output}/${base}_novo.sam" > ./split_mappingScripts/${base}.sh
```

__4.2) Create script to call all and run in parallel (use "&" which puts job in background then multiple can run at a time)__

This creates a file that has all the scripts made in step 5.1) in a list with ''&'' at the end to run in parrallel

*Script: [novo_createParallel_scipt.sh](https://github.com/PaulKnoops/Experimental_Evolution_Sequence_Repo/blob/master/RunThrough_DataCleaning_Scripts/novo_createParallel_scipt.sh)*

Ex.
```
echo "${map_scripts}/${base}.sh &" >> ${scripts}/novo_parallel_map.sh
```

__4.3) Run novo_parallel_map.sh__

Running this line runs each sequence through novoalign in parallel

```
novo_parallel_map.sh
```

__5) Save space again with trimmed files: Rezip__

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
____________________________________________________________________________________________________

### Change SAM files to BAM files: novo_samTobam.sh

-- Saves space (BAM files are binary compressed versions of SAM files)

-- need bam directory for .bam files (mkdir novo_bam)

Flags:

    - b -- output is .bam
    - S -- input is .sam
    - q 20 -- quality mapping score of 20 (standard throughout all experiments)
     

*Script: [novo_samTobam.sh](https://github.com/PaulKnoops/Experimental_Evolution_Sequence_Repo/blob/master/RunThrough_DataCleaning_Scripts/novo_samTobam.sh)*

Ex.
```
samtools view -b -S -q 20 ${novo_dir}/${base}.sam | samtools sort -o ${novo_bam}/${base}.bam
```

____________________________________________________________________________________________________
### Merge 

No flags, just merging the two lanes of the illumina sequencing run into one file:

*The script: [novo_merge.sh](https://github.com/PaulKnoops/Experimental_Evolution_Sequence_Repo/blob/master/RunThrough_DataCleaning_Scripts/novo_merge.sh)*


Ex.
```    
samtools merge ${novo_merge}/${base}_novo_merge.bam ${novo_bam}/${base}_L001_novo.bam ${novo_bam}/${base}_L002_novo.bam
```
____________________________________________________________________________________________________

### Picard Sort 

Need to sord with Picard to mark and remove duplicates (needs to be sorted with Picard in order for downstream analysis)

Need a directory for outputs, and a temporary directory for space allocation (novo_tmp)

Flags:

    - Xmx2g -- 2 Gb of memory allocated
    - Djava.io.tmpdir=${tmp} -- using my temporary directory due to errors in space allocation to avoid errors while running (not necessary but useful)
    - I -- input
    - O -- output
    - VALIDATION_STRINGENCY=SILENT -- stops Picard from reporting every issue that would ultimately be displayed
    - SO=coordinate -- sort order based on coordinate
    - TMP_DIR -- temporary directory
  
 
*Script: [novo_picard_sort.sh](https://github.com/PaulKnoops/Experimental_Evolution_Sequence_Repo/blob/master/RunThrough_DataCleaning_Scripts/novo_picard_sort.sh)*

Ex.
```
java -Xmx2g -Djava.io.tmpdir=${novo_tmp} -jar ${pic} SortSam \
I= ${novo_merge}/${base}.bam \
O= ${novo_pic}/${base}_novo_sort.bam \
VALIDATION_STRINGENCY=SILENT SO=coordinate TMP_DIR=${novo_tmp}
```
____________________________________________________________________________________________________

### Remove Duplicates

Flags: Similar to above

    - Xmx2g -- ""
    - MarkDuplicates -- ""
    - I -- ""
    - O -- ""
    - M -- creates an output file of statistics of duplicates found
    - VALIDATION_STRINGENCY=SILENT -- ""
    - REMOVE_DUPLICATES= true -- get rid of any found duplicated regions


*Script: [novo_rmd.sh](https://github.com/PaulKnoops/Experimental_Evolution_Sequence_Repo/blob/master/RunThrough_DataCleaning_Scripts/novo_rmd.sh)*

Ex.
```
java -Xmx2g -jar ${picard_dir}/picard.jar MarkDuplicates I= ${novo_pic}/${base}_novo_sort.bam O= ${novo_rmd}/${base}_novo_rmd.bam M= ${novo_rmd}/dupstat.txt VALIDATION_STRINGENCY=SILENT REMOVE_DUPLICATES= true

```
____________________________________________________________________________________________________

### More QC and make final bam files

Flags:

    - q 20 -- ""
    - F 0x0004 -- remove any unmapped reads (hexidecimal value for unmapped = 0x0004)
    - b -- ""

*Script: [novo_final_qc.sh](https://github.com/PaulKnoops/Experimental_Evolution_Sequence_Repo/blob/master/RunThrough_DataCleaning_Scripts/novo_final_qc.sh)*

Ex.
```
samtools view -q 20 -F 0x0004 -b ${novo_rmd}/${base}_novo_rmd.bam > ${novo_final}/${base}_novo_final.bam
```
____________________________________________________________________________________________________

### Merge the 2 ancestors

Need to merge the base generation additionaly (two sequence runs for ancestor need to merge: MGD2 and MGD)

*Script: [novo_merge_anc.sh](https://github.com/PaulKnoops/Experimental_Evolution_Sequence_Repo/blob/master/RunThrough_DataCleaning_Scripts/novo_merge_anc.sh)*

____________________________________________________________________________________________________

### Indel Realigner with GATK (6 steps)

__1) Need an unzipped copy of the reference genome__

Made a second copy and unzipping the second copy (may already be present from Bowtie2 mapping)
```
gunzip dmel-all-chromosome-r5.57_2.fasta.gz
```

__2) Set up reference genome__

__2.1) Create .dict file for reference:__ 

This creates a dictionary file for the ref genome with a header but no sam records (the header is only sequence records)

*script: [novo_dict_index.sh](https://github.com/PaulKnoops/Experimental_Evolution_Sequence_Repo/blob/master/RunThrough_DataCleaning_Scripts/novo_dict_index.sh)*
```
java -jar ${picard_dir}/picard.jar CreateSequenceDictionary R=${reference_genome} O=${index_dir}/dmel-all-chromosome-r5.57_2.dict
```

__2.2) Create a .fai to reference genome: __
```
samtools faidx dmel-all-chromosome-r5.57_2.fasta
```

__3) Need Read Group headers__

For GATK to run, read group headers are needed to read for Indel Realignment, __BUT__ the details are not necessary and don't need to be accurate, just need to be there.

Flags:

    - RGID --Read Group Identifier; for Illumina, are composed using the flowcell + lane name and number [using Lanes L001_L002 for now]
    - RGLB -- DNA Preperation Library Identifier [library1 as place holder]
    - RGPL -- platform/technology used to produce the read [Illumina]
    - RGPU -- Platform Unit; details on the sequencing unit (i.e run barcode) [None, used for practice]
    - RGSM -- Sample [Using the basename which is each unique sequence]

*Script: [novo_readgroups.sh](https://github.com/PaulKnoops/Experimental_Evolution_Sequence_Repo/blob/master/RunThrough_DataCleaning_Scripts/novo_readgroups.sh)*
```
java -jar picard.jar AddOrReplaceReadGroups I=${novo_final}/${base}.bam O=${novo_final}/${base}_RG.bam RGID=L001_L002 RGLB=library1 RGPL=illumina RGPU=None RGSM=${base}
```

__4) Index the read group Bam files__

Need to have the .bam files indexed prior to the indel realignment

*Script: [novo_indexBam.sh](https://github.com/PaulKnoops/Experimental_Evolution_Sequence_Repo/blob/master/RunThrough_DataCleaning_Scripts/novo_indexBam.sh)*
```
samtools index ${novo_final}/${base}_RG.bam
```

__5) Run GATK indel realigner__

GATK indel realigner takes two steps, __1)__ target the indels to be raligned (.intervals file) and __2)__ realign the indels (.realigned files)

*Script: [novo_gatk.sh](https://github.com/PaulKnoops/Experimental_Evolution_Sequence_Repo/blob/master/RunThrough_DataCleaning_Scripts/novo_gatk.sh)*
```
java -Xmx8g -jar ${gatk_dir}/GenomeAnalysisTK.jar -I ${novo_final}/${base}_RG.bam -R ${ref_genome} -T RealignerTargetCreator -o ${novo_GATK}/${base}.intervals

java -Xmx8g -jar ${gatk_dir}/GenomeAnalysisTK.jar -I ${novo_final}/${base}_RG.bam -R ${ref_genome} -T IndelRealigner -targetIntervals ${novo_GATK}/${base}.intervals -o ${novo_GATK}/${base}_realigned.bam
```

The above steps for section 3 can be repeated for other mappers (with some changes to parameters) to have three sets of data files for 13 populations with final .bam files for further analysis

____________________________________________________________________________________________________
____________________________________________________________________________________________________
____________________________________________________________________________________________________
____________________________________________________________________________________________________
____________________________________________________________________________________________________


## **Section 4:** Create mpileup and .sync file for analysis

General here for novoalign: but similar for other mappers again

### Create mpileup

Flags:

    - B -- disable BAQ (base alignment quality) computation, helps to stop false SNPs passing through due to misalignment
    - Q -- minimum base quality (already filtered for 20, default is 13, just set to 0 and not worry about it)
    - f -- path to reference sequence
    

*Script: [novo_mpileup.sh](https://github.com/PaulKnoops/Experimental_Evolution_Sequence_Repo/blob/master/RunThrough_DataCleaning_Scripts/novo_mpileup.sh)*

Ex.
```
samtools mpileup -B -Q 0 -f ${reference_genome} ${novo_final}/*.bam > ${novo_mpileup}/${project_name}.mpileup
```

### Create .sync file 

Flags:

    - Xmx7g -- 7 Gb of memory allocated
    - input -- Input
    - output -- Output
    - fastq-type -- needed for base encoding
    - min-qual 20 -- already set to 20 before, but since custome script, use 20 as safer assumption
    - threads 2 -- Threads used

*Script: [novo_sync.sh](https://github.com/PaulKnoops/Experimental_Evolution_Sequence_Repo/blob/master/RunThrough_DataCleaning_Scripts/novo_sync.sh)*

Ex.
```
java -ea -Xmx7g -jar ${popoolation_dir}/mpileup2sync.jar --input ${novo_mpileup}/${project_name}.mpileup --output ${novo_mpileup}/${project_name}.sync --fastq-type sanger --min-qual 20 --threads 2
```
