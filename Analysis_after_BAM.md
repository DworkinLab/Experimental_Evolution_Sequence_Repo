# The analysis of pooled sequence data after the final .bam files and .sync files have been created. 

- The analaysis here is for data mapped and finalized with three mappers: **bwa mem**, **bowtie2** and **novoalign** 

- See [associated script](https://github.com/PaulKnoops/Experimental_Evolution_Sequence_Repo/blob/master/RunThrough_DataCleaning.md) for steps up to final bams  
_______________________________________________________________________________________

## Outline of analysis:

### 1) Tajima's Pi of non-overlapping windows for each sequence
/
/
### 2) Run Fst on windows for each pairwise comparision of sequenced data and calculate average Fst across three mappers
/
/
### 3) per SNP logistic regression for each treatment by generation averaged for Novoalign, Bwa-mem and bowtie2
/
/
### 4) Average estimates of selection coefficient at each position for selection and control lineages for three mappers
/
/
### 5) Positions of interest for Fst, poolseq and model output
/
/
### 6) Trajectory of regions of interest based on model, Fst and selection coefficients
/
/

_______________________________________________________________________________________

## Notes / additional set up:

- For one mapper: Have a diretory with all the .final.bam files created and .mpileup /.sync files created using these .bam files

- Some analysis will be shown for one mapper but completed similarily for other mappers

- need to know the order of the .sync files: will be based on the order of the .bam files read into .sync 

- Some scripts require the .sync file split into chromosomes: script below shows method of splitting Novoalign sync

**Script:** [novo_split_sync_Chromosome.sh](https://github.com/PaulKnoops/Experimental_Evolution_Sequence_Repo/blob/master/Analysis_after_BAM_Scripts/novo_split_sync_Chromosome.sh)
_______________________________________________________________________________________

## 1) Tajima's Pi of non-overlapping windows for each sequence

### Create single pileup files for every .bam file

To run Pi function for popoolation1: each bam file has its own pileup format (created with mpileup)

Flags:

    - B -- disable BAQ (base alignment quality) computation, helps to stop false SNPs passing through due to misalignment
    - Q -- minimum base quality (already filtered for 20, default is 13, just set to 0 and not worry about it)
    - f -- path to reference sequence
    
*Script: [novo_PI_pileups.sh](https://github.com/PaulKnoops/Experimental_Evolution_Sequence_Repo/blob/master/Analysis_after_BAM_Scripts/novo_PI_pileups.sh)*

Ex.
```
samtools mpileup -B -Q 0 -f ${ref_genome} ${input}/${base}_merge_novo_final_realigned.bam > ${output}/${base}.pileup
```

### Run script to calcualte Tajima's Pi using the Variance-sliding.pl script from Popoolation1

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


*Script: [novo_tajima_pi.sh](https://github.com/PaulKnoops/Experimental_Evolution_Sequence_Repo/blob/master/Analysis_after_BAM_Scripts/novo_tajima_pi.sh)*

Ex. 
```
perl ${popoolation}/Variance-sliding.pl --input ${input}/${base}.pileup --output ${output}/${base}.pi --measure pi --window-size 10000 --step-size 10000 --min-count 2 --min-coverage 4 --max-coverage 400 --min-qual 20 --pool-size 120 --fastq-type sanger --snp-output ${output}/${base}.snps --min-covered-fraction 0.5
```
### Create plots of tajima Pi data

On local machine, this R function can run each .pi file to output a plot

**Script:** [Pi_PlotFunction.R](https://github.com/PaulKnoops/episodicSequenceData/blob/master/Analysis_after_sync_2018_scripts/Pi_PlotFunction.R)

### In R, run the function for each .pi file

ex. 
```
Pi_PlotFunction('FILE.pi', "Plot Title Details")
```

### Example Output: Ancestral Pi for Novoalign

![Ancestral Pi Plot for Novoalign](https://github.com/PaulKnoops/episodicSequenceData/blob/master/Analysis_after_sync_2018_plots/Ancestral_Pi.png)


**Questions**

Currently, Pi files for Novoalign and Bowtie2 generated.

- Should these be averaged (with bwa-mem that can be generated easily enough)

- is one mapper enough for this

- Do we want to show trends in nucleotide diversity (other Pi plots) or just focus on ancestoral Pi (overlay plots)

_______________________________________________________________________________________

## 2) Fst on windows of each pairwise comparision of sequences

### Running Fst

**Script:** [Novo_Fst.sh](https://github.com/PaulKnoops/episodicSequenceData/blob/master/Analysis_after_sync_2018_scripts/Novo_Fst.sh)

ex.
```
perl ${fst} --input ${novo_mpileup}/novo_episodic_main.sync --output ${novo_fst}/novo_episodic_main.fst --min-count 6 --min-coverage 10 --max-coverage 250 --min-covered-fraction 1 --window-size 500 --step-size 500 --pool-size 120
```
Flags:

- input -- input sync file
- output -- output file with Fst calculated 
- window-size [500] -- size of the window 
- step-size [500] -- distance to move along chromosome
- min-count [6] -- minimum allele count 
- min-coverage [10] -- minimum coverage
- max-coverage [250] --maximum coverage
- pool-size [120] -- double pooled size (diploid)
- min-covered-fraction [1] -- minimum percentage of sites having sufficient coverage in the given window


### In R, split the file into each compasison

**Script:** [novo_Fst_Split_Comparisons.R](https://github.com/PaulKnoops/episodicSequenceData/blob/master/Analysis_after_sync_2018_scripts/novo_Fst_Split_Comparisons.R)

R script that will split the .fst file into many .csv files with each comparison of interest (can choose the necessary ones from here)

### Combining three mappers:

**Script**[FST_combine3mappers.R](https://github.com/PaulKnoops/episodicSequenceData/blob/master/Analysis_after_sync_2018_scripts/FST_combine3mappers.R)

Will take the split comparisons, and combine those specified into one FST file with the average Fst between the three mappers

### Plotting Fst files for comparisons of interest

 - average Fst of three mappers **and** average between replicates
 
 - comparison betweeen control and selection lines
 
Generation 38:
[meanFst for F38](https://github.com/PaulKnoops/episodicSequenceData/blob/master/Analysis_after_sync_2018_plots/F38_meanFstPlot.png)

Generation 77: 
[meanFst for F77](https://github.com/PaulKnoops/episodicSequenceData/blob/master/Analysis_after_sync_2018_plots/F77_meanFstPlot.png)

Generation 115: 
![meanFst for F115](https://github.com/PaulKnoops/episodicSequenceData/blob/master/Analysis_after_sync_2018_plots/F115_meanFstPlot.png)




_______________________________________________________________________________________

## 3) per SNP logistic regression for each treatment by generation

**Long Script:*** [novo_regression_model_LONGSCRIPT.sh](https://github.com/PaulKnoops/episodicSequenceData/blob/master/Analysis_after_sync_2018_scripts/novo_regression_model_LONGSCRIPT.sh)

This script will break the chromosomal .sync files into smaller managable pieces and run through multiple R scripts while removing intermediates:

**R script to covert sync to Count data:** [Sync_to_counts.R](https://github.com/PaulKnoops/episodicSequenceData/blob/master/Analysis_after_sync_2018_scripts/Sync_to_counts.R)

Creates a file with the counts for the major and minor frequency (based on ancestor) that can run through the model

**R script for running the model for each position along the chromosome:** [Counts_to_model.R](https://github.com/PaulKnoops/episodicSequenceData/blob/master/Analysis_after_sync_2018_scripts/Counts_to_model.R)

In long script: this is set up to work in parallel, having each chromosome running at the same time (6 instances running over 11 sections)

NOTE: Script needs to be changed to run faster/ more efficiently (not done here b/c already completed)

Basic Model at each positon (tmp2):
```
modlist_2[[i]] <- 
        glm(cbind(Major_count, Minor_count) ~ Treatment*Generation, 
            data = tmp2, family = "binomial")
```

**R script to combine all the split chromosome pieces back into one chromosome:** [Combine_chromo.R](https://github.com/PaulKnoops/episodicSequenceData/blob/master/Analysis_after_sync_2018_scripts/Combine_chromo.R)

Recreates one chromosomal file

**R script to combine three mappers into one file** [model_combine3mappers.R](https://github.com/PaulKnoops/episodicSequenceData/blob/master/Analysis_after_sync_2018_scripts/model_combine3mappers.R)

Combines each of BWA-mem, Bowtie2 and Novoalign files into one file (keeping all information)

**R script to write files with coeffefficent of interest** [model_3mappersTxG.R](https://github.com/PaulKnoops/episodicSequenceData/blob/master/Analysis_after_sync_2018_scripts/model_3mappersTxG.R)

This script (choosing Treatment by Generation effect) keeps positions that are present in all three files (i.e position needs to be mapped three times)

**P.adjust?, Positions?,etc.**

**Plots**

Treatment x Generation -log10(meanP-value) for model output: NEED TO CHANGE FOR P.adjust!
![FullGenomeTxGPlot](https://github.com/PaulKnoops/episodicSequenceData/blob/master/Analysis_after_sync_2018_plots/CHROMO_meanP.png)

_______________________________________________________________________________________

## 4) estimates of selection coefficient at each position for selection and control lineages using [poolSeq](https://github.com/ThomasTaus/poolSeq) R package:

**Script:** [poolseq_SelectionCoefficientEstimate.sh](https://github.com/PaulKnoops/episodicSequenceData/blob/master/Analysis_after_sync_2018_scripts/poolseq_SelectionCoefficientEstimate.sh)

This script will break the file into two treatment .sync files, break apart these .sync files (smaller sized files), and run through a R script to run poolSeq Package. 

**Rscript: Running poolseq** [poolSeq_selectionCoeff.R](https://github.com/PaulKnoops/episodicSequenceData/blob/master/Analysis_after_sync_2018_scripts/poolSeq_selectionCoeff.R)

Note: to run, check poolseq is available, if not, source all PoolSeq scripts available from Taus git page. 

Also, to run with modified Sync files (Change is spacing), need a personal read.sync function: --> modifed script == [read.sync_personal_function.R](https://github.com/PaulKnoops/episodicSequenceData/blob/master/Analysis_after_sync_2018_scripts/read.sync_personal_function.R)

**Rscript: combining CSV files:** [combinePoolseqCSV.R](https://github.com/PaulKnoops/episodicSequenceData/blob/master/Analysis_after_sync_2018_scripts/combinePoolseqCSV.R)

Ex. For individual Chromosomes: [combine_poolseq_individual_Chromo.R](https://github.com/PaulKnoops/episodicSequenceData/blob/master/Analysis_after_sync_2018_scripts/combine_poolseq_individual_Chromo.R)

Notes: 

- need to edit to make more efficient to run (not taking 3 days per section)

- poolSeq may not work on certain versions of R; but can bring in Taus poolSeq scripts in manually and source (done in this script)

- may need to ensure updated packages and install if necessary (i.e matrixStats_0.53.0 installed for this reason)

- breaking the .sync file causes changes in structure, so a modified read_sync function is used (in R script; Taus_ReadSync.R))

_______________________________________________________________________________________

## 5) Comparison of Fst (of generation 115, control:Selection) and Model Outputs:

**First draft of plots** 

**Full Chromosome:**

 -- Numbers/ peaks do not line up perfectly (product of how they are plotted/ measured (windows vs. positions)
 
![FstV.model_full](https://github.com/PaulKnoops/episodicSequenceData/blob/master/Analysis_after_sync_2018_plots/log10p_fst_comboFUll.png)

 
_______________________________________________________________________________________

## 6) Trajectory of regions of interest based on model output

### Finding positions: 

**Fst window and ranges in data frame:**
```
### Adjust the FST output (FDR) and keep positons that have an Fst value after adjusting
 require(data.table)
 require(tidyverse)
### Read in the data (final generation) with comparison between Control and Selection lines:
 XC2 <- fread('combined_fst_1:3.csv')
 CX2 <- fread('combined_fst_2:4.csv')
### One data frame:
 ddat <- rbind(CX2, XC2)
### Mean Fst b/w two replicates:
 ddat2 <- ddat %>%
   group_by(chr, window) %>%
   summarise(meanFst = (mean(meanFst)))
### Fdr correction: Adjust the calculated Fst values with a false discovery rate:
 ddat2$adjustFst <- p.adjust(ddat2$meanFst, method = 'fdr')
### Remove positions with an Fst of 0
 ddat2 <-  ddat2[-which(ddat2$adjustFst==0),]
### The range of the 500 bp window (window == center)
 ddat2$minWindow <- ddat2$window -250
 ddat2$maxWindow <- ddat2$window +250

### need to sort positions based on chromosome (for trajectories):
 fst_X <- ddat2[which(ddat2$chr=='X'),]
 fst_2L <- ddat2[which(ddat2$chr=='2L'),]
 fst_2R <- ddat2[which(ddat2$chr=='2R'),]
 fst_3L <- ddat2[which(ddat2$chr=='3L'),]
 fst_3R <- ddat2[which(ddat2$chr=='3R'),]
 fst_4 <- ddat2[which(ddat2$chr=='4'),]

### create positional data frame (to be sources later:
 posfst_X <- as.data.frame(cbind(fst_X$window, fst_X$minWindow, fst_X$maxWindow))
 posfst_2L <- as.data.frame(cbind(fst_2L$window, fst_2L$minWindow, fst_2L$maxWindow))
 posfst_2R <- as.data.frame(cbind(fst_2R$window, fst_2R$minWindow, fst_2R$maxWindow))
 posfst_3L <- as.data.frame(cbind(fst_3L$window, fst_3L$minWindow, fst_3L$maxWindow))
 posfst_3R <- as.data.frame(cbind(fst_3R$window, fst_3R$minWindow, fst_3R$maxWindow))
 posfst_4 <- as.data.frame(cbind(fst_4$window, fst_4$minWindow, fst_4$maxWindow))

### Keep only position data frame for sourcing:
 rm(list=ls()[! ls() %in% c('posfst_X', "posfst_2L", 'posfst_2R', 'posfst_3R', 'posfst_3L','posfst_4')])
```
**Positions with significant change after FDR adjustments:**
```
### Finding Positions of interest from Model:

### Packages:
 require(dplyr)
 require(ggplot2)
 require(data.table)

### Read in each Chromosomal Data:  
 ddatX <- fread('../Data/X_Chromosome_TxG.csv', h=T)
 chr_X <- ddatX %>%
   group_by(position, chr, Effects) %>%
   summarise(mean_p = (mean(p.value)))
 rm(ddatX)
 a <- nrow(chr_X)
 chr_X$number <- 1:a

 ddat2L <- fread('../Data/2L_Chromosome_TxG.csv', h=T)
 chr_2L <- ddat2L %>%
   group_by(position, chr, Effects) %>%
   summarise(mean_p = (mean(p.value)))
 rm(ddat2L)
 b <- nrow(chr_2L)
 chr_2L$number <- (a+1):(a+b) 

 ddat2R <- fread('../Data/2R_Chromosome_TxG.csv', h=T)
 chr_2R <- ddat2R %>%
   group_by(position, chr, Effects) %>%
   summarise(mean_p = (mean(p.value)))
 rm(ddat2R)
 c <- nrow(chr_2R)
 chr_2R$number <- (a+b+1):(a+b+c)

 ddat3L <- fread('../Data/3L_Chromosome_TxG.csv', h=T)
 chr_3L <- ddat3L %>%
   group_by(position, chr, Effects) %>%
   summarise(mean_p = (mean(p.value)))
 rm(ddat3L)
 d <- nrow(chr_3L)
 chr_3L$number <- (a+b+c+1):(a+b+c+d)

 ddat3R <- fread('../Data/3R_Chromosome_TxG.csv', h=T)
 chr_3R <- ddat3R %>%
   group_by(position, chr, Effects) %>%
   summarise(mean_p = (mean(p.value)))
  rm(ddat3R)
 e <- nrow(chr_3R)
 chr_3R$number <- (a+b+c+d+1):(a+b+c+d+e)

 ddat4 <- fread('../Data/4_Chromosome_TxG.csv', h=T)
 chr_4 <- ddat4 %>%
   group_by(position, chr, Effects) %>%
   summarise(mean_p = (mean(p.value)))
 rm(ddat4)
 f <- nrow(chr_4)
 chr_4$number <- (a+b+c+d+e+1):(a+b+c+d+e+f)
 chr_4$chr <- as.character(chr_4$chr)

### One large data frame:
 CHROMOs <- rbind(chr_X, chr_2L, chr_2R, chr_3L, chr_3R, chr_4)

### Fdr correction: adjust p values with false discrovey rate correction:
 CHROMOs$adjustP <- p.adjust(CHROMOs$mean_p, method = 'fdr')
### Filter slightly to make it easier to wotk with:
 CHROMOs_3 <-  CHROMOs[which(CHROMOs$adjustP<0.9),]

# Split into chromosomes and keep significant postions:
 ddat2_X <- CHROMOs_3[which(CHROMOs_3$chr=='X' & CHROMOs_3$adjustP<0.05),]
 pos_X <- ddat2_X$position
 ddat2_2L <- CHROMOs_3[which(CHROMOs_3$chr=='2L' & CHROMOs_3$adjustP<0.05),]
 pos_2L <- ddat2_2L$position
 ddat2_2R <- CHROMOs_3[which(CHROMOs_3$chr=='2R' & CHROMOs_3$adjustP<0.05),]
 pos_2R <- ddat2_2R$position
 ddat2_3L <- CHROMOs_3[which(CHROMOs_3$chr=='3L' & CHROMOs_3$adjustP<0.05),]
 pos_3L <- ddat2_3L$position
 ddat2_3R <- CHROMOs_3[which(CHROMOs_3$chr=='3R' & CHROMOs_3$adjustP<0.05),]
 pos_3R <- ddat2_3R$position
 ddat2_4 <- CHROMOs_3[which(CHROMOs_3$chr=='4' & CHROMOs_3$adjustP<0.05),]
 pos_4 <- ddat2_4$position

### Remove all but positions:
 rm(list=ls()[! ls() %in% c('pos_X','pos_2L','pos_2R','pos_3L','pos_3R','pos_4')])
```

**Positions showing significant selection coefficents in selection but not controls:**
```
### Must be done for each chromosome seperatly: 
###-----    2R     -----###
 novo_sel_2R <- fread('novo_episodic_2R_Sel.csv')
 novo_con_2R <- fread('novo_episodic_2R_Con.csv')
 bwa_sel_2R <- fread('bwa_episodic_2R_Sel.csv')
 bwa_con_2R <- fread('bwa_episodic_2R_Con.csv')

 column.names <- c('selcoef', 'pval', 'pos', 'chr')
 colnames(novo_sel_2R) <- column.names
 colnames(novo_con_2R) <- column.names
 colnames(bwa_sel_2R) <- column.names
 colnames(bwa_con_2R) <- column.names

 novo_sel_2R <- na.omit(novo_sel_2R)
 novo_con_2R <- na.omit(novo_con_2R)
 bwa_sel_2R <- na.omit(bwa_sel_2R)
 bwa_con_2R <- na.omit(bwa_con_2R)

 novo_sel_2R$map <- 'Novo'
 novo_con_2R$map <- 'Novo'
 bwa_sel_2R$map <- 'bwa'
 bwa_con_2R$map <- 'bwa'

 novo_sel_2R$Treatment <- 'Sel'
 novo_con_2R$Treatment <- 'Con'
 bwa_sel_2R$Treatment <- 'Sel'
 bwa_con_2R$Treatment <- 'Con'

 Xcx_2R <- rbind(novo_con_2R, novo_sel_2R, bwa_con_2R, bwa_sel_2R)
###-----    2R     -----###

### Repeat this for all chromosomes (changing instances of 2R to chr) to create one large data frame:
 Xcx <- rbind(Xcx_2R, Xcx_2L)

### Make sure that all 2 (3) mappers have the position for chr/treatment:
 CXC <- Xcx %>%
   group_by(chr, Treatment, pos) %>%
   mutate(count = n())
 XCV <- CXC[which(CXC$count==2),]

### Get the mean selection coefficent and the least significant pvalue (max(pval)) between mappers:
 Zxc <- XCV %>%
   group_by(chr, Treatment, pos) %>%
   summarise(meanSelCoef = (mean(selcoef)),
             pval_max=max(pval))

### Adjust the p-value left with false discovery rate and keep significant postions:
 Zxc$adjustP <- p.adjust(Zxc$pval_max, method = 'fdr')
 Zxc$sig <- ifelse(Zxc$adjustP<0.05, "<0.05", ">0.05")
 Zxc_sig <- Zxc[which(Zxc$sig=='<0.05'),]

### Check if there are positons that are significant and present for both control and selection (and keep only positions that are only significant for selection:

 Zxc_count <- Zxc_sig %>%
   group_by(chr, pos) %>%
   mutate(count = n())
 Zxc_count2 <- Zxc_count[which(Zxc_count$count==1),]
 Zxc_count2 <- Zxc_count2[which(Zxc_count2$Treatment=='Sel'),]

### Positions by chromosome:

 poolseq_X <- Zxc_count2[which(Zxc_count2$chr=='X'),]
 pool_pos_X <- poolseq_X$pos
 poolseq_2L <- Zxc_count2[which(Zxc_count2$chr=='2L'),]
 pool_pos_2L <- poolseq_2L$pos
 poolseq_2R <- Zxc_count2[which(Zxc_count2$chr=='2R'),]
 pool_pos_2R <- poolseq_2R$pos
 poolseq_3L <- Zxc_count2[which(Zxc_count2$chr=='3L'),]
 pool_pos_3L <- poolseq_3L$pos
 poolseq_3R <- Zxc_count2[which(Zxc_count2$chr=='3R'),]
 pool_pos_3R <- poolseq_3R$pos
 poolseq_4 <- Zxc_count2[which(Zxc_count2$chr=='4'),]
 pool_pos_4 <- poolseq_4$pos
 
 ### keep only positions:
  rm(list=ls()[! ls() %in% c('pool_pos_2L', 'pool_pos_2R', 'pool_pos_3L', 'pool_pos_3R', 'pool_pos_4', 'pool_pos_X')])
```



Ex. with BWA -mem output sync files

**Rscript:** [extract_sig_Chromo_positions.R](https://github.com/PaulKnoops/episodicSequenceData/blob/master/Analysis_after_sync_2018_scripts/extract_sig_Chromo_positions.R)

running rscript:
```
Rscript extract_sig_Chromo_positions.R
```
