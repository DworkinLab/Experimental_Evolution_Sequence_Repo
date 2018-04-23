
## Analysis of the variation

### Tajima's Pi

1) Create Directories

Will need two new directories, one for pileup files, and one for the output Pi files: Be sure they are in project dir

```
mkdir novo_pileup
mkdir novo_pi
```

2) Create pileup files of every .bam file

Each generated .bam file needs to be in mpileup format for the use in Popoolation Scripts:

Flags:

- B -- disable BAQ (base alignment quality) computation, helps to stop false SNPs passing through due to misalignment

- Q -- minimum base quality (already filtered for 20, default is 13, just set to 0 and not worry about it)

- f -- path to reference sequence
        
Script: novo_Pi_pileups.sh
```
#! /bin/bash

# Variable for project:
project_dir=/home/paul/episodicData/novoalign

# Path to input directory
input=${project_dir}/novo_GATK

# Path to output novoalign pileup files
output=${project_dir}/novo_pileup

index_dir=/home/paul/episodicData/index_dir
ref_genome=${index_dir}/dmel-all-chromosome-r5.57_2.fasta

files=(${input}/*_merge_novo_final_realigned.bam)

for file in ${files[@]}

do

name=${file}

base=`basename ${name} _merge_novo_final_realigned.bam`

samtools mpileup -B -Q 0 -f ${ref_genome} ${input}/${base}_merge_novo_final_realigned.bam > ${output}/${base}.pileup

done
```

3) Run script to calcualte Tajima's Pi

Using the Variance-sliding.pl script from Popoolation1

Flags:

- input -- input pileup file

- output -- output file with Tajima's Pi calculated

- measure [pi] -- Options include Tajima's Pi or Wattersons Theta or Tajima's D along chromosomes using a sliding window approach

- window-size [10000] -- size of the sliding window 

- step-size [10000] -- how far to move along with chromosome (if step size smaller, windows will overlap)
 
- min-count [2] -- minimum allele count 

- min-coverage [4] -- minimum coverage (not important if subsampling done..)

- max-coverage [400] --maximum coverage

- min-qual [20] -- minimum base quality (already filtered for 20 multiple times)

- pool-size [120] -- number of chromosomes (So double the number of individuals per pool)

- fastq-type [sanger] -- depending on the encoding of the fastq files

- min-covered-fraction [0.5] -- minimum percentage of sites having sufficient coverage in the given window -- 0.5 from example


Script: novo_tajima_pi.sh

```
#! /bin/bash

# Path to PoPoolation1 (Currently in Paul's Home directory)
popoolation=/home/paul/popoolation_1.2.2

# Variable for project:
project_dir=/home/paul/episodicData/novoalign

# Path to input directory
input=${project_dir}/novo_pileup

# Path to output Tajima Pi files
output=${project_dir}/novo_pi

files=(${input}/*.pileup)

for file in ${files[@]}

do

name=${file}

base=`basename ${name} .pileup`

perl ${popoolation}/Variance-sliding.pl \
	--input ${input}/${base}.pileup \
	--output ${output}/${base}.pi \
	--measure pi \
	--window-size 10000 \
	--step-size 10000 \
	--min-count 2 \
	--min-coverage 4 \
	--max-coverage 400 \
	--min-qual 20 \
	--pool-size 120 \
	--fastq-type sanger \
	--snp-output ${output}/${base}.snps \
	--min-covered-fraction 0.5
done
```
Need to check the files were created properly (check that there are no empty files).


4) Bring to local machine
```
scp paul@info.mcmaster.ca:/home/paul/episodicData/novoalign/novo_pi/*.pi /Users/paulknoops/Bioinformatics/R-projects_git/episodicSequenceData/R_scripts/Pi_Analysis_Novo
```

5) Function to create plots of tajima Pi data

R Function to run and create plots

See: R_scripts/Pi_Analysis_Novo/Pi_PlotFunction.R
```
Pi_PlotFunction <- function(x, y) {
  require(ggplot2)
  x2 <- gsub("\\_.*","",x)
  y2 <- y
  
  #bowtie.pi
  #novo.pi
  #Read in the data:
  Datt <- read.table(x)
  colnames(Datt) <- c('chr', 'window', 'windowCount', ' propInwindow', 'Pi')
  
  #Remove unnecessary regions: Not necessary based on later steps
  Datt$chr <- as.character(Datt$chr)
  Datt2 <- Datt

  #Remove "na" pi values
  Datt2 <- Datt2[-which(Datt2$Pi=="na"),]
  
  #Need the numbers for chromosomes for labelling and colours:
  DattX <- Datt2[which(Datt2$chr=="X"),]
  a <- dim(DattX)[1]
  DattX$number <- 1:a
  
  Datt2L <- Datt2[which(Datt2$chr=="2L"),]
  b <- dim(Datt2L)[1]
  Datt2L$number <- (a+1):(a+b)
  
  Datt2R <- Datt2[which(Datt2$chr=="2R"),]
  c <- dim(Datt2R)[1]
  Datt2R$number <- (a+b+1):(a+b+c)
  
  Datt3L <- Datt2[which(Datt2$chr=="3L"),]
  d <- dim(Datt3L)[1]
  Datt3L$number <- (a+b+c+1):(a+b+c+d)
  
  Datt3R <- Datt2[which(Datt2$chr=="3R"),]
  e <- dim(Datt3R)[1]
  Datt3R$number <- (a+b+c+d+1):(a+b+c+d+e)
  
  Datt4 <- Datt2[which(Datt2$chr=="4"),]
  f <- dim(Datt4)[1]
  Datt4$number <- (a+b+c+d+e+1):(a+b+c+d+e+f)
  
  #Full data frame of necessary chromosomes
  DattFull <- rbind(DattX, Datt2L, Datt2R, Datt3L, Datt3R, Datt4)
  
  #Pi as numeric
  DattFull$Pi=as.numeric(levels(DattFull$Pi))[DattFull$Pi]
  
  #Title:
  z2 <- paste(x2, y2, sep="_")
  
  # The plots: 
  Pi_plot <- ggplot(DattFull, aes(x = number, y= Pi, colour = chr)) 
  
  Pi_plot_2 <- Pi_plot + 
    geom_point(size=0.3, show.legend = F) +
    scale_y_continuous(limits=c(0, 0.02), breaks=seq(0, 0.02, 0.005)) + 
    xlab("") +
    scale_x_discrete(limits=c(1049, 3185, 5277, 7443, 9952, 11359), labels = c("X", "2L", '2R', '3L', '3R', "4")) +
    theme(text = element_text(size=20), 
          axis.text.x= element_text(size=15), axis.text.y= element_text(size=15)) +
    scale_colour_manual(values=c("#56B4E9", "#E69F00", 'grey30', 'grey46', 'wheat3', 'lemonchiffon4')) +
    ggtitle(z2)
  
  return(Pi_plot_2)
}
```

6) Run the function for each .pi file

```
Pi_PlotFunction('F115ConR1_TAGCTT_novo.pi', "Novoalign")
```














_________________________________________________________________________________________________________________________
_________________________________________________________________________________________________________________________
_________________________________________________________________________________________________________________________
_________________________________________________________________________________________________________________________
_________________________________________________________________________________________________________________________
_________________________________________________________________________________________________________________________
_________________________________________________________________________________________________________________________
_________________________________________________________________________________________________________________________
_________________________________________________________________________________________________________________________
_________________________________________________________________________________________________________________________
_________________________________________________________________________________________________________________________
_________________________________________________________________________________________________________________________
_________________________________________________________________________________________________________________________
_________________________________________________________________________________________________________________________
## Combining files together for analysis


### Create mpileup

```
mkdir novo_mpileup
```

Flags:

- B -- disable BAQ (base alignment quality) computation, helps to stop false SNPs passing through due to misalignment

- Q -- minimum base quality (already filtered for 20, default is 13, just set to 0 and not worry about it)

- f -- path to reference sequence
    

Script: novo_mpileup.sh
```
#!/bin/bash

#Variable for project name (title of mpileup file)
project_name=novo_episodic

#Variable for project:
project_dir=/home/paul/episodicData/novoalign

#Path to input directory
novo_final=${project_dir}/novo_GATK

#Path to output directory
novo_mpileup=${project_dir}/novo_mpileup

#Variable for reference genome (non-indexed for novoAlign)
index_dir=/home/paul/episodicData/index_dir
ref_genome=${index_dir}/dmel-all-chromosome-r5.57_2.fasta

samtools mpileup -B -Q 0 -f ${ref_genome} ${novo_final}/*.bam > ${novo_mpileup}/${project_name}.mpileup
```

### Create .sync file 

--use mpileup dir

Flags:

- Xmx7g -- ""

- input -- ""

- output -- ""

- fastq-type -- needed for base encoding

- min-qual 20 -- already set to 20 before, but since custome script, use 20 as safer assumption

- threads 2 -- ""

Script: 
```
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
```

### Split .sync file into chromosomes (easier to work with)
```
#!/bin/bash

#Variable for project name (title of mpileup file)
project_name=novo_episodic

#Variable for project:
project_dir=/home/paul/episodicData/novoalign

#Path to .sync files
novo_mpileup=${project_dir}/novo_mpileup

grep -v 'Het' ${novo_mpileup}/${project_name}.sync > ${novo_mpileup}/${project_name}_less_het.sync

wait

grep -v 'U' ${novo_mpileup}/${project_name}_less_het.sync > ${novo_mpileup}/${project_name}_removed_U_Het.sync

wait

grep -v 'dmel_mitochondrion_genome' ${novo_mpileup}/${project_name}_removed_U_Het.sync > ${novo_mpileup}/${project_name}_main.sync

wait

rm -f ${novo_mpileup}/${project_name}_less_het.sync

rm -f ${novo_mpileup}/${project_name}_removed_U_Het.sync

grep '3R' ${novo_mpileup}/${project_name}_main.sync > ${novo_mpileup}/${project_name}_3R.sync &
grep '2R' ${novo_mpileup}/${project_name}_main.sync > ${novo_mpileup}/${project_name}_2R.sync &
grep '3L' ${novo_mpileup}/${project_name}_main.sync > ${novo_mpileup}/${project_name}_3L.sync &
grep '2L' ${novo_mpileup}/${project_name}_main.sync > ${novo_mpileup}/${project_name}_2L.sync &
grep '^4' ${novo_mpileup}/${project_name}_main.sync > ${novo_mpileup}/${project_name}_4.sync &
grep 'X' ${novo_mpileup}/${project_name}_main.sync > ${novo_mpileup}/${project_name}_X.sync 
```

Sync File now set up to be easily worked with; with one version (being `_main.sync`) that has the 6 chromosomes of interest, and 6 files for these chromosomes seperate

# Easier layout for Fst, Poolseq, and Model scripts (& Pi) can be found [here](https://github.com/PaulKnoops/episodicSequenceData/blob/master/Analysis_after_sync_2018.md)

## Fst:

```
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

```


Split the file into each compasison: Run in R
```
# Script to read in one .fst file and split this into many .csv files based on comparisons
# .fst files generated from Popoolation2 fst-sliding.pl script
# Need to change the working directory, the input, and the number of comparisons present (i.e 6:ncol for .fst file)

### Packages Required (tidyverse, but more specifically tidyr)
  require(data.table)
  require(tidyverse)

### Set working directory to location of .fst file
  setwd("~/Bioinformatics/episodic_practice/FST")
  
### Read in the .fst file into R (requires data.table)
  Fst_novo <- fread('novo_episodic_main.fst')

### Make into long format
  XCC  <- gather(Fst_novo, Comparison, Fst_measure, 6:83, factor_key=TRUE)

### Remove intermediate:
  rm(Fst_novo)

### Seperate the Fst (ex. 1:2=na) into a comparison column and a fst column
  novo_Fst <- XCC %>%
    separate(Fst_measure, "=", into = c('Comp', 'Fst'))

### Remove intermediate:
  rm(XCC)

### Remove unnecessary column (column 6 has no value)
  novo_Fst <- novo_Fst[,c(1,2,3,4,5,7,8)]

### Rename columns:
  colnames(novo_Fst) <- c('chr', 'window', "num", 'frac', 'meanCov','Comp', 'Fst')

### Create list of all unique comparisons:
  X_compLIST <- unique(novo_Fst$Comp)

### For loop that will create a .csv file for every comparison:
  for (vxcx in X_compLIST){

    CXV_Comp <- novo_Fst[which(novo_Fst$Comp==vxcx),]
    
    write.csv(CXV_Comp, file=paste("Novo_fst_", vxcx, '.csv', sep = ""), row.names = FALSE)
    }

```
Which comparisons are needed? Find the order of files in the sync file (order in directory that .sync file created from)

i.e 1 == F115ConR1_TAGCTT_novo_merge_novo_final_realigned.bam, 13 == MGD3_SO_CAGATC_novo_merge_novo_final_realigned.bam
```
F115ConR1_TAGCTT_novo_merge_novo_final_realigned.bam
F115ConR2_GGCTAC_novo_merge_novo_final_realigned.bam
F115SelR1_GTTTCG_novo_merge_novo_final_realigned.bam
F115SelR2_GTGGCC_novo_merge_novo_final_realigned.bam
F38ConR1_ATCACG_novo_merge_novo_final_realigned.bam
F38ConR2_TTAGGC_novo_merge_novo_final_realigned.bam
F38SelR1_ACTTGA_novo_merge_novo_final_realigned.bam
F38SelR2_GATCAG_novo_merge_novo_final_realigned.bam
F77ConR1_ATGTCA_novo_merge_novo_final_realigned.bam
F77ConR2_ATTCCT_novo_merge_novo_final_realigned.bam
F77SelR1_TTAGGC_novo_merge_novo_final_realigned.bam
F77SelR2_GATCAG_novo_merge_novo_final_realigned.bam
MGD3_SO_CAGATC_novo_merge_novo_final_realigned.bam
```
## Running Custom R script: Logistic Regression

**LONG SCRIPT: (Change input and output, add Rscripts, etc.)**

Some assumptions are made in R scripts that may not work with this data: suggestion is to attempt a trail run with R scripts using 4th chrmosome (the smallest) to make sure all scripts work

The Basic model: glm of the major and minor allele count at each position for Treatment by Generation effects (ignoring cage)
```
      modlist_2[[i]] <- 
        glm(cbind(Major_count, Minor_count) ~ Treatment*Generation, 
            data = tmp2, family = "binomial")
```

The long script: basic outline == Breaks up sync files (per chromo) into small peices, converts to major/minor allele counts, runs model per position, and then brings back into one chromosome. 

Things need to be changed to run more efficiently

```
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
```

Now have 6 .csv files with output from model for each chromosome (check naming of 4th, might not include 4). 

### R Scripts for model: sync to counts, model itself etc.

**Script: Sync_to_counts.R**
To run on own:
```
Rscript Sync_to_counts.R '/home/paul/episodicData/novoalign/novo_mpileup'
```
```
# For loop sync to counts

## need next line to call arguments:

args <- commandArgs(trailingOnly = TRUE)


## Convert a .sync file into long format, filter somewhat, and have only position, treatment, Cage, Generation and Maj/Min counts

## Packages source code: only need these two for this script (need to be this order)

require('tidyr')
require('dplyr')

#1) Need to change details as needed above and below string of #####

#2) Needs a .sync file made by popoolation2

#3) Need to change most importantly for analysis the read in and read out names 

# Read in Data: Big Data Sets

#pwd a direcotry containing only the directories of interest (made with other sed -n script)


mydirs <- list.dirs(path = args[1], recursive = FALSE)

#includes that actual dir.. not with recursive = FALSE

for (dir in mydirs){

    setwd(dir)
  
  mysyncs <- list.files(pattern=".sync")
  
  for (sync in mysyncs){
  
      episodic_counts <- read.table(sync)
    
    #adjust colnames
    print("data read in")
    name.Columns <- c("Chromosome", "Position", "ref", "ConR1_115", "ConR2_115", "SelR2_115", "SelR1_115", "ConR1_38", "ConR2_38", "SelR1_38", "SelR2_38", "ConR1_77", "ConR2_77", "SelR1_77", "SelR2_77", "SelR1_0")
    colnames(episodic_counts) <- name.Columns
    
    #Add "replicates" of ancestor -- all are equal
    episodic_counts$SelR2_0 <- episodic_counts$SelR1_0
    episodic_counts$ConR1_0 <- episodic_counts$SelR1_0
    episodic_counts$ConR2_0 <- episodic_counts$SelR1_0
    
    
    #Need the ancestor to stay (after making long) to call major/minor alleles later
    episodic_counts$Ancestor <- episodic_counts$SelR1_0
    
    # Make long by bring populations down
    print("making long")
    long_episodic <- gather(episodic_counts, Population, Allele_Freq , ConR1_115:ConR2_0, factor_key=TRUE)
    
    rm(episodic_counts)
    print("removed counts")
    
    #Error???
    
    # All geneneric below for sync files (only real issue through file is population naming convention)
    ###################################################
    
    #Seperate the allele counts into independent columns for each base
    print("splitting allele freq")
    Episodic_split_2 <- long_episodic %>% 
      separate(Allele_Freq, c("A","T","C","G","N","del"), ":")
    
    rm(long_episodic)
    
    print("removed long")
    
    #Seperate the ancestor to base later things on
    Episodic_split_2 <- Episodic_split_2 %>% 
      separate(Ancestor, c("A_0","T_0","C_0","G_0","N_0","del_0"), ":")
    
    # as.numeric to multiple columns:
    cols.num <- c("A_0", "T_0", "C_0", "G_0", "N_0", "del_0", "A", "T", "C", "G", "N", "del")
    
    #Seems to take a long time for this step?
    Episodic_split_2[cols.num] <- sapply(Episodic_split_2[cols.num],as.numeric)
    
    #Get the sum of all the rows (all the different bases) for each population position:
    
    print("getting row sums")
    Episodic_split_2$sum <- (rowSums(Episodic_split_2[,11:16]))
    
    #Ancestor Major_Allele and minor allele:
    
    # Major allele of ancestor == the maximum positional count
    Episodic_split_2$anc_max <- apply(Episodic_split_2[,4:9], 1, max)
    # Minor is the ancestor second highest count
    Episodic_split_2$anc_min <- apply(Episodic_split_2[,4:9], 1, 
                                      function(x)max(x[x!=max(x)]))
    
    #Major / Minor Base name: match the number of anc_max with the column to call the correct base:
    
    Episodic_split_2 <- within(Episodic_split_2, {
      MajorAllele = ifelse(anc_max== Episodic_split_2[,4], "A", ifelse(anc_max== Episodic_split_2[,5],  "T", ifelse(anc_max== Episodic_split_2[,6],  "C",ifelse(anc_max== Episodic_split_2[,7],  "G", ifelse(anc_max== Episodic_split_2[,8],  "N", ifelse(anc_max== Episodic_split_2[,9],  "del", "N/A" ))))))})
    
    #Major Allele Count of evolved populations; match the Major allele with the count of certain columns for each population 
    
    Episodic_split_2 <- within(Episodic_split_2, {
      Maj_count = ifelse (MajorAllele == "A", Episodic_split_2[,11], ifelse (MajorAllele == "T", Episodic_split_2[,12], ifelse (MajorAllele == "C", Episodic_split_2[,13], ifelse (MajorAllele == "G", Episodic_split_2[,14], ifelse (MajorAllele == "N", Episodic_split_2[,15], ifelse (MajorAllele == "del", Episodic_split_2[,16], "N/A"))))))})
    
    
    # Same thing for minor allele: first ensure that if the sum of all counts == the Major coutn and the ancestor had no minor allele, their is no minor allele (N/A), then follow the same match of anc_min to a certain base
    
    Episodic_split_2 <- within(Episodic_split_2, {
      MinorAllele = ifelse(Maj_count==Episodic_split_2[,17] & anc_min==0, "N/A", ifelse(anc_min== Episodic_split_2[,4], "A", ifelse(anc_min== Episodic_split_2[,5],  "T", ifelse(anc_min== Episodic_split_2[,6],  "C",ifelse(anc_min== Episodic_split_2[,7],  "G", ifelse(anc_min== Episodic_split_2[,8],  "N", ifelse(anc_min== Episodic_split_2[,9],  "del", "Z") ))))))})
    
    
    #Minor Allele Count of the ancestreal minor allele count
    Episodic_split_2 <- within(Episodic_split_2, {
      Min_count = ifelse (MinorAllele == "A", Episodic_split_2[,11], ifelse (MinorAllele == "T", Episodic_split_2[,12], ifelse (MinorAllele == "C", Episodic_split_2[,13], ifelse (MinorAllele == "G", Episodic_split_2[,14], ifelse (MinorAllele == "N", Episodic_split_2[,15],ifelse (MinorAllele == "del", Episodic_split_2[,16],"N/A"))))))})
    
    print("called major and minor alleles and counts")
    # To determine the minor allele base if not specified by the ancestor (new allele brough up etc.)
    
    #max for the population (could be the minor allele)
    Episodic_split_2$maj_all <- apply(Episodic_split_2[,11:16], 1, max)
    
    #alt== second highest count for populations
    Episodic_split_2$alt_allele <- apply(Episodic_split_2[,11:16], 1, 
                                         function(x)max(x[x!=max(x)]))
    
    print("define unknown alleles")
    Episodic_split_2 <- within(Episodic_split_2, {
      Min_count_2 = ifelse (Maj_count == sum, 0, ifelse(Maj_count==maj_all, alt_allele, maj_all))})
    
    Episodic_split_2 <- within(Episodic_split_2, {
      MinorAllele_base = ifelse(Min_count_2==0, "N/A", ifelse(Min_count_2== Episodic_split_2[,11], "A", ifelse(Min_count_2== Episodic_split_2[,12],  "T", ifelse(Min_count_2== Episodic_split_2[,13],  "C",ifelse(Min_count_2== Episodic_split_2[,14],  "G", ifelse(Min_count_2== Episodic_split_2[,15],  "N", ifelse(Min_count_2== Episodic_split_2[,16],  "del", "Z") ))))))})
    
    # Remove unneeded columns (6,7,8,9,10,11,13,14,15)
    Episodic_split_2 <- subset(Episodic_split_2, select = -c(A_0,T_0,C_0,G_0,N_0,del_0,A,T,C,G,N,del,anc_max,anc_min, MinorAllele, Min_count, maj_all, alt_allele))
    
    print("removed unneeded columns")
    nam.col <- c("chr", "pos", "ref", "Population", "sum", "MajorAllele", "Major_count", "Minor_count", "MinorAllele")
    colnames(Episodic_split_2) <- nam.col
    
    
    #Remove unneccessary Columns (as needed)
    #Keep them all for now (except sum) as may be needed later
    #Episodic_split_2 <- subset( Episodic_split_2, select = -ref )
    #Episodic_split_2 <- subset( Episodic_split_2, select = -chr)
    #Episodic_split_2 <- subset( Episodic_split_2, select = -MajorAllele )
    #Episodic_split_2 <- subset( Episodic_split_2, select = -MinorAllele)
    Episodic_split_2<- subset( Episodic_split_2, select = -sum)
    
    
    ## Depends on the filter method:
    print("begin filtering")
    #Filter method: take the sum of each position, and must have at least 5 counts called (i.e over the 16 populations, the total of alleles called for the minor allele must be over 5)
    grp <- Episodic_split_2 %>%
      group_by(pos) %>%
      summarise(sum=sum(Minor_count))
    grp2 <- grp[which(grp$sum<=5),]
    Episodic_split_2 <- Episodic_split_2[!(Episodic_split_2$pos %in% grp2$pos),]
    
    
    #check that the number of obs for episodic_long2 == obs for those without 0's sum (*16 for number of "populations") (needed later as well == grp3)
    
    #grp3 <- grp[-which(grp$sum<=5),]
    rm(grp)
    rm(grp2)
    
    print("remove filter inermediates")
    #################################################
    #Should be all genetic above (from start specificed)
    
    ## Below depends on the population name layout etc. made above
    
    
    #Split Population into Treatment, Rep, and Generation - need to do twice, different seperators (change above??)
    
    print("seperate population to Treatment, Generation and Cage")
    episodic_long <- Episodic_split_2 %>%
      separate(Population, c("Treatment", "Generation"), "_")
    
    rm(Episodic_split_2)
    
    episodic_long <- episodic_long %>%
      separate(Treatment, c("Treatment", "Cage"), "R")
    
    cols.num <- c("Cage", "Generation", "Major_count", "Minor_count")
    episodic_long[cols.num] <- sapply(episodic_long[cols.num],as.numeric) 
    
    print("Have final episodic long; now write a csv")
  
    #will need to rename .csv files  
    write.csv(episodic_long, file=paste(sync, ".csv", sep=""))
  
    
    print("wrote csv and now done this .sync file")
  }
}

```


**Script: Counts_to_model_2.R**
	- Script for running each chromosome in unison
To run on own:
```
Rscript Counts_to_model.R 'DIRECTORY'
```
Script to run model for each chromosome seperatly
```
#Episodic data analysis: loop .csv files to run model:

args <- commandArgs(trailingOnly = TRUE)

  setwd(args[1])
  
  mycsvs <- list.files(pattern=".csv")
  
  for (file in mycsvs){
    
    episodic_long <- read.csv(file, h=T)
    
    #The Data: in long format, each position with Treatment, Cage and Generation, along with the Major and Mnor allele counts correponding to the ancestral major/minor allele
    
    #The full model:
    
    #Call each position
    
    position <- unique(episodic_long$pos)
    
    no.pos <- length(position)
    
    #Remove N/A's -- possibly not needed so hashed out.
    
    episodic_long <- na.omit(episodic_long)
    
    
    #Make list to store model
    
    modlist_2 <- as.list(1:no.pos)
    
    #Each model of a position, named for each mod will be position
    
    names(modlist_2) <- position
    
    #Empty Data Frame to store all coeffecients of model
    
    coeffs_df <- data.frame(NULL)
    
    #Run the  model for each position
    
    for(i in position){
      print(paste("Running entity:", i, "which is", which(position==i), "out of", no.pos, "file=", file))
      
      #Temporary data frame for only the one position
      
      tmp2 <- episodic_long[episodic_long$pos == i,]
      
      #The model: major vs. minor counts by Treatment, Generation and Treatment:Generation
      
      modlist_2[[i]] <- 
        glm(cbind(Major_count, Minor_count) ~ Treatment*Generation, 
            data = tmp2, family = "binomial")
      
      #Turn this model into a data frame of coefficients
      
      x <- as.data.frame(summary(modlist_2[[i]] )$coefficients)
      
      #Name the position of this model results with i
      
      x$position <- i
      x$chr <- tmp2$chr[1]
      #Add to data frame (total for the whole data set == coeffs_df + the newly made X)
      
      coeffs_df <- rbind(coeffs_df, x)
      
      #Remove i for safety and it starts over
      
      rm(i)
    }
    
    #Change column names to workable
    
    colnames(coeffs_df) <- c("Estimate", "Standard_error", "z-value", "p-value", "position", "chr")
    
    coeffs_df$Effects<-rownames(coeffs_df)
    
    coeffs_df$Effects_2 <- ifelse(grepl("TreatmentSel:Generation",coeffs_df$Effects),'T_Sel:Gen', ifelse(grepl("Intercept",coeffs_df$Effects),'Int', coeffs_df$Effects ))
    
    coeffs_df$Effects_2 <- ifelse(grepl("TreatmentSel",coeffs_df$Effects_2),'T_Sel', ifelse(grepl("Generation",coeffs_df$Effects_2),'Gen', coeffs_df$Effects_2))
    
    rownames(coeffs_df) <- c()
    
    #Make the p-values into -log10 p values
    coeffs_df$log_p <- -log10(coeffs_df$`p-value`)
    
    coeffs_df <- subset(coeffs_df, select = -c(Effects))
    
    coeffs_df$Effects <- ifelse(coeffs_df$Effects_2=='T_Sel', 'TreatmentSel', ifelse(coeffs_df$Effects_2=='Gen', 'Generation', ifelse(coeffs_df$Effects_2=='Int', 'Intercept', 'TreatmentSel:Generation')))
    
    coeffs_df <- subset(coeffs_df, select = -c(Effects_2))
    
    coeffs_df <- coeffs_df[-which(coeffs_df$log_p==0),]
    
    write.csv(coeffs_df, file=paste(file,".coeffs.csv", sep=""))
        rm(coeffs_df)
    rm(tmp2)
    rm(x)
    rm(modlist_2)
    rm(episodic_long)
    rm(no.pos)
    rm(position)
  }
```


**Script: Combine_chromo.R:**
To run on own:
```
Rscript Combine_chromo.R 'DIRECTORY' 'OutputDIRECTORY'
```
```

# Combine .coeffs.csv for two mappers
require(dplyr)

args <- commandArgs(trailingOnly = TRUE)

# Change to directory holding all directories:

mydirs <- list.dirs(path = args[1], recursive = FALSE)

for (dir in mydirs){

  setwd(dir)
  
  print("Read coeffs.csv files")

  mycsvs <- list.files(pattern=".coeffs.csv")

  Novoalign_Chromosome <- NULL
  
  for (file in mycsvs){
    print(file)
    coeffs1 <- read.csv(file, h=T)
    Novoalign_Chromosome  <- rbind(Novoalign_Chromosome , coeffs1)
    rm(coeffs1)
    J2 <- gsub("*_","", file)
    J3 <- gsub("\\..*","",J2)
}

x3 <- gsub("\\..*","",file)
J3 <- gsub('(.*)_\\w+', '\\1', x3)

X <- args[2]

write.csv(Novoalign_Chromosome , file=paste(X,"/",J3,"_chromo.csv", sep=""), row.names = FALSE)
rm(x3)
rm(J3)
rm(Novoalign_Chromosome)

}
```

## Next steps: combine chromosomes and effects of interest
__________________________________________________________

### Can combine with other mappers

Need a directory to hold chromosomal outputs:
```
mkdir /home/paul/coeffs_fullChromo
```


To run: have args[1] be the output directory created just above (which can hold script to avoid confusion with other mappers)
```
Rscript combine_threeMappers.R '/home/paul/coeffs_fullChromo'
```

```
#To run as Rscript and have output directory argument:
args <- commandArgs(trailingOnly = TRUE)

# combine the coeffs csv files for three mappers
require(dplyr)

# Will call these 6 Chromosoms (but this depends on the files all ending with [name2]_chromo.csv (ex. 2L_chromo.csv)
name2 <- c('2L', '2R', '3L', '3R', '4','X')

#Loop through all chromosomes: each read.csv has the directory holding all the chromosome subsets from the output of the model run just above. 

# Will combine each chromosome into one file with three mappers (defined in loops) and written to a new large .csv file

for (i in name2){

  print(i)
  
  file=paste(i,'_chromo.csv', sep = "")
  
  bwa_coeffs <- read.csv(paste('/home/paul/episodicData/R_dir/bwa_coeffs/episodic_data_', file, sep=""), header = TRUE)
  bwa_coeffs$mapper <- "bwa"
  
  bowtie_coeffs <- read.csv(paste('/home/paul/episodicData/bowtie/R_bowtie/bowtie_coeffs/episodic_data_bowtie_', file, sep=""), header = TRUE)
  bowtie_coeffs$mapper <- "bowtie"
  
  novo_coeffs <- read.csv(paste('/home/paul/episodicData/novoalign/novo_coeffs/novo_episodic_', file, sep=""), header = TRUE)
  novo_coeffs$mapper <- "novoalign"
  
  X <- rbind(bwa_coeffs, bowtie_coeffs, novo_coeffs)
  rm(coeffs_bowtie)
  rm(coeffs_bwa)
  rm(novo_coeffs)
  
  outputDir <- args[1]
 
  write.csv(X , file=paste(outputDir,"/",'Chromosome_', i, "_Full.csv", sep=""), row.names = FALSE)
}
```
### Choose an effect and only keep those with positions mapped on by all three mappers:

```
Rscript TxG_allThreeMappedPos.R '/home/paul/coeffs_fullChromo'
```
The script:
```
#To run as Rscript and have output directory argument:
args <- commandArgs(trailingOnly = TRUE)
#Can input manually:

setwd(args[1])
#setwd('/home/paul/coeffs_fullChromo')

#Can change the effect of interest by changing hashes at X_Effect (#EFFECT OF INTEREST#)

require(dplyr)
mycsvs <- list.files(pattern='_Full.csv')


for (file in mycsvs){
  #file = 'Chromosome_4_Full.csv'
  print(file)
  name <- gsub("^.*?_","",file)
  name2 <- gsub("_.*","",name)
  rm(name)
  X <- read.csv(file, h=T)

  
#EFFECT OF INTEREST#
#X_Effect <- X[which(X$Effects=="TreatmentSel"),]
X_Effect <- X[which(X$Effects=="TreatmentSel:Generation"),]
#X_Effect <- X[which(X$Effects=="Intercept"),]
#X_Effect <- X[which(X$Effects=="Generation"),]

title <- as.character(X_Effect$Effects[1])
title2 <- ifelse(title=='TreatmentSel', 'Treat', ifelse(title=='Generation', 'Gen', ifelse(title=='TreatmentSel:Generation', 'TxG', '
                                                                                           Int')))

rm(title)
rm(X)

# Only keep those with a position mapped by all three (create table of posiitons and all those >=3 kept)
tt <- table(X_Effect$position)
Effects_Final <- subset(X_Effect, position %in% names(tt[tt >= 3]))

rm(tt)
rm(X_Effect)

print('writing CSV')
write.csv(Effects_Final, file=paste0(name2, '_Chromosome_', title2, '.csv'), row.names = FALSE)


rm(Effects_Final)

print('Done and everything gone')
}
```
## Plotting and p.adjust:




_____________________________________________________________
_____________________________________________________________

# Running [poolSeq](https://github.com/ThomasTaus/poolSeq) R package:
**Note: running in 15 sections (more then below) took from March 20th -- April 17th!**

```
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
	  	Rscript ${Rscripts}/PoolSeq_SelCoeff.R ${file} ${Chromo} ${split_sync}) &
	done 
	wait
	rm -f ${split_sync}/*.sync
done
wait

Rscript ${Rscripts}/combinePoolseqCSV.R ${splitSync}

##------------------------------------------------##
```

## Pool Seq R Scripts: running poolseq, combining csv

**Script: PoolSeq_SelCoeff.R:**
To run on own:
```
Rscript PoolSeq_SelCoeff.R ".SYNCFILE" "CHROMOSOME_OF_SYNC" "INPUTDIR"
```
```
### Running the PoolSeq Package to run on .sync files for estimates of selection coefficeints per position
### Requires: R (>= 3.3.1), data.table (>= 1.9.4), foreach (>= 1.4.2), stringi (>= 0.4-1), matrixStats (>= 0.14.2)

  args <- commandArgs(trailingOnly = TRUE)
  
### Required Packages:

  #install.packages("/home/paul/poolSeq_0.3.2.tar.gz", repos=NULL, type="source")
  #install.packages("/home/paul/matrixStats_0.53.0.tar.gz", repos=NULL, type="source")

### Not available: so source seperate:
  #require(poolSeq)
  
### These are part of the dependencies for poolSeq
  
  require(methods)
  require(data.table)
  require(foreach)
  require(stringi)
  require(matrixStats)
  
### Source the scripts (Copied) for Pool-Seq (only one fails and is not needed)
  source('/home/paul/episodicData/novoalign/novo_Rscripts/Taus_Scripts/testTaus/loadaf.R')  
  #estne.R Fails.
  #source('/home/paul/episodicData/novoalign/novo_Rscripts/Taus_Scripts/testTaus/estne.R')
  source('/home/paul/episodicData/novoalign/novo_Rscripts/Taus_Scripts/testTaus/estsh.R')
  source('/home/paul/episodicData/novoalign/novo_Rscripts/Taus_Scripts/testTaus/idsel.R')
  source('/home/paul/episodicData/novoalign/novo_Rscripts/Taus_Scripts/testTaus/simaf.R')

### Possibly need custom function to read in manipulated .sync files:
	### Needed for manipulated .sync files (one basic change labeled at top of changed script:
  
  source("/home/paul/episodicData/novoalign/novo_Rscripts/Taus_ReadSync.R")

### Read in the data file for args[1]

  setwd(args[3])
  
  mySync <- read.sync_Personal(file=args[1], gen=c(115, 115, 38, 38, 77, 77, 0, 0), repl=c(1,2,1,2,1,2,1,2), polarization = "rising")

### Make data.frame of just alleles information to sort out relevent positions:

  ff <- as.data.frame(mySync@alleles)
  pst <- as.numeric(ff$pos)
  pst2 <- sort(pst)
  rm(pst)
  rm(ff)

### Create empty data frame to read into for estiamting S:
  
  DF <- data.frame(NULL)
  ccc <- c(0,38,77,115)

  for (i in pst2) {
  	Traj115 <- af.traj(mySync, args[2], repl=c(1,2), pos=i)
  	Bfsf <- estimateSH(Traj115, Ne=150, t=ccc, h=0.5, haploid = FALSE, simulate.p.value=TRUE)
  	Fd <- data.frame(Bfsf$s, Bfsf$p0, Bfsf$p.value)
 	 Fd$pos <- i
  	DF <- rbind(DF, Fd)
  	DF <- na.omit(DF)
  	#print(paste("Running entity:", i, "of", END))
  	rm(i)
	
	}
 
  x2 <- args[1]
  x3 <- gsub("\\..*","", x2)
  write.csv(DF, file=paste(args[3], "/", x3, ".csv", sep=""), row.names=FALSE)
  
  rm(DF)
  rm(mySync)
  rm(ccc)
  rm(pst2)
```


**Script: combinePoolseqCSV.R:**
To run on own:
```
Rscript combinePoolseqCSV.R "Dir_Holding_all_SPLIT_DIR's"
``` 
```
args <- commandArgs(trailingOnly = TRUE)
mydirs <- list.dirs(path=args[1], recursive = FALSE) 
XX <- args[1]
for (dir in mydirs){
	setwd(dir)
	J4 <- gsub('(.*)_\\w+', '\\1', dir)
	mycsvs <- list.files(pattern='.csv')
	X <- NULL
	for (file in mycsvs){
		X2 <- read.csv(file, h=T)
		X <- rbind(X, X2)
   	}
	write.csv(X, file=paste(J4,'.csv',sep=""), row.names = FALSE)
	rm(X)
	rm(J4)
}
```

 -- Note: may want to move combine into big (dir) loop and change script to not loop through all dir's but just the current one dir?
 
 ## Combine mappers and plot / position extract:




# Positions of interest: finding positions overlapping in model output and selection coefficients with high Con:Sel Fst values
