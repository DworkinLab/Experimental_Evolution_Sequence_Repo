# The analysis of pooled sequence data after the final .bam files and .sync files have been created. 

- The analaysis here is for data mapped and finalized with three mappers: **bwa mem**, **bowtie2** and **novoalign** 

- See [associated script](https://github.com/PaulKnoops/Experimental_Evolution_Sequence_Repo/blob/master/RunThrough_DataCleaning.md) for steps up to final bams  
_______________________________________________________________________________________

## Outline of analysis:

### 1) Tajima's Pi of non-overlapping windows for each sequence
/
/
### 2) Fst on windows for each pairwise comparision of sequenced data and calculate average Fst across three mappers
/
/
### 3) per SNP logistic regression for each treatment by generation averaged for Novoalign, Bwa-mem and bowtie2
/
/
### 4) Average estimates of selection coefficient at each position for selection and control lineages for two mappers
/
/
### 5) Positions of interest for Fst, poolseq and model output (overlap)
/
/
### 6) Running Gowinda for gene analysis from positions
/
/
### 7) Trajectory of regions of interest based on model, Fst and selection coefficients
/
/

_______________________________________________________________________________________

## Notes / additional set up:

- For one mapper: Have a diretory with all the .final.bam files created and .mpileup /.sync files created using these .bam files

- Some analysis will be shown for one mapper but completed similarily for other mappers

- need to know the order of the .sync files: will be based on the order of the .bam files read into .sync 

- Some scripts require the .sync file split into chromosomes: Ex. ***Script:** [novo_split_sync_Chromosome.sh](https://github.com/PaulKnoops/Experimental_Evolution_Sequence_Repo/blob/master/Analysis_after_BAM_Scripts/novo_split_sync_Chromosome.sh)*
_______________________________________________________________________________________

## 1) Tajima's Pi of non-overlapping windows for each sequence

### Create single pileup files for every .bam file

To run Pi function for popoolation1: each bam file has its own pileup format (created with mpileup)

Flags:

    - B -- disable BAQ (base alignment quality) computation, helps to stop false SNPs passing through due to misalignment
    - Q -- minimum base quality (already filtered for 20, default is 13, just set to 0 and not worry about it)
    - f -- path to reference sequence
    
***Script:** [novo_PI_pileups.sh](https://github.com/PaulKnoops/Experimental_Evolution_Sequence_Repo/blob/master/Analysis_after_BAM_Scripts/novo_PI_pileups.sh)*

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


***Script:** [novo_tajima_pi.sh](https://github.com/PaulKnoops/Experimental_Evolution_Sequence_Repo/blob/master/Analysis_after_BAM_Scripts/novo_tajima_pi.sh)*

Ex. 
```
perl ${popoolation}/Variance-sliding.pl --input ${input}/${base}.pileup --output ${output}/${base}.pi --measure pi --window-size 10000 --step-size 10000 --min-count 2 --min-coverage 4 --max-coverage 400 --min-qual 20 --pool-size 120 --fastq-type sanger --snp-output ${output}/${base}.snps --min-covered-fraction 0.5
```

Outputs of data were able to be loaded here:

***Novoalign Pi data:** [Pi_Novoalign](https://github.com/PaulKnoops/Experimental_Evolution_Sequence_Repo/tree/master/Data/Pi_Novoalign)*

***Bowtie Pi data:** [Pi_Bowtie](https://github.com/PaulKnoops/Experimental_Evolution_Sequence_Repo/tree/master/Data/Pi_Bowtie)*

### Create plots of tajima Pi data

This R function can run each .pi file to output a plot

***Script:** [Pi_plot_function.R](https://github.com/PaulKnoops/Experimental_Evolution_Sequence_Repo/blob/master/Analysis_after_BAM_Scripts/Pi_plot_function.R)*

This script can be updated and modified for different details on the plots

### In R, run the function for each .pi file

Ex. 
```
Pi_PlotFunction('FILE.pi', "Plot Title Details")
```
_______________________________________________________________________________________

## 2) Fst on windows of each pairwise comparision of sequences

### Running Fst

In 500 bp windows: calculates Fst values for each pairwise comparison between sequences (1-13) within the 500 bp window

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

**Script:** [novo_Fst.sh](https://github.com/PaulKnoops/Experimental_Evolution_Sequence_Repo/blob/master/Analysis_after_BAM_Scripts/novo_Fst.sh)

ex.
```
perl ${fst} --input ${novo_mpileup}/novo_episodic_main.sync --output ${novo_fst}/novo_episodic_main.fst --min-count 6 --min-coverage 10 --max-coverage 250 --min-covered-fraction 1 --window-size 500 --step-size 500 --pool-size 120
```

### In R, split the file into each compasison

**Script:** [novo_Fst_Split_Comparisons.R](https://github.com/PaulKnoops/Experimental_Evolution_Sequence_Repo/blob/master/Analysis_after_BAM_Scripts/novo_Fst_Split_Comparisons.R)

R script that will split the .fst file into many .csv files with each comparison (can choose the necessary ones from here)

### Combining three mappers output:

**Script:** [FST_combine3mappers.R](https://github.com/PaulKnoops/Experimental_Evolution_Sequence_Repo/blob/master/Analysis_after_BAM_Scripts/FST_combine3mappers.R)

Will take the split comparisons, and combine those specified into one FST file with the average Fst between the three mappers 

By changing the "patty" variable, the comparisons of interest can be combined and evaluated by taking the matching comparisons from three directories holding the different mappers output from previous script spliting FST output.

Some combined data files of interest can be found here: [Fst_combinedComparisons](https://github.com/PaulKnoops/Experimental_Evolution_Sequence_Repo/tree/master/Data/Fst_combinedComparisons)

### Plotting Fst files for Con:Sel comparison at three generations

 - average Fst of three mappers **and** average between replicates
 
 - comparison betweeen control and selection lines
 
 **Script:** [Fst_Plots.R](https://github.com/PaulKnoops/Experimental_Evolution_Sequence_Repo/blob/master/Analysis_after_BAM_Scripts/Fst_Plots.R)
 
Generation 38:
[meanFst for F38](https://github.com/PaulKnoops/episodicSequenceData/blob/master/Analysis_after_sync_2018_plots/F38_meanFstPlot.png)

Generation 77: 
[meanFst for F77](https://github.com/PaulKnoops/episodicSequenceData/blob/master/Analysis_after_sync_2018_plots/F77_meanFstPlot.png)

Generation 115: 
![meanFst for F115](https://github.com/PaulKnoops/episodicSequenceData/blob/master/Analysis_after_sync_2018_plots/F115_meanFstPlot.png)

_______________________________________________________________________________________

## 3) per SNP logistic regression for each treatment by generation

**Long Script:*** [novo_regression_model_LONGSCRIPT.sh](https://github.com/PaulKnoops/Experimental_Evolution_Sequence_Repo/blob/master/Analysis_after_BAM_Scripts/novo_regression_model_LONGSCRIPT.sh)

This script will break the chromosomal .sync files (i.e split per chromosome) into smaller managable pieces and run through multiple R scripts while removing intermediates:

The R script below are within the long script:

**R script to covert sync to Count data:** [Sync_to_counts.R](https://github.com/PaulKnoops/Experimental_Evolution_Sequence_Repo/blob/master/Analysis_after_BAM_Scripts/Sync_to_counts.R)

Creates a file with the counts for the major and minor frequency (based on ancestor) that can run through the model

**R script for running the model for each position along the chromosome:** [Counts_to_model.R](https://github.com/PaulKnoops/Experimental_Evolution_Sequence_Repo/blob/master/Analysis_after_BAM_Scripts/Counts_to_model.R)

In long script: this is set up to work in parallel, having each chromosome running at the same time (6 instances running over 11 sections)

NOTE: Script needs to be changed to run faster/ more efficiently (not done here b/c already completed)

Basic Model at each positon (tmp2):
```
modlist_2[[i]] <- 
        glm(cbind(Major_count, Minor_count) ~ Treatment*Generation, 
            data = tmp2, family = "binomial")
```


**After** long script complete:

**R script to combine all the split chromosome pieces back into one chromosome:** [Combine_model_Chromo.R](https://github.com/PaulKnoops/Experimental_Evolution_Sequence_Repo/blob/master/Analysis_after_BAM_Scripts/Combine_model_Chromo.R)


Recreates one chromosomal file

**R script to combine three mappers into one file** [model_combine3mappers.R](https://github.com/PaulKnoops/Experimental_Evolution_Sequence_Repo/blob/master/Analysis_after_BAM_Scripts/model_combine3mappers.R)

Combines each of BWA-mem, Bowtie2 and Novoalign files into one file (keeping all information)


**R script to write files with coeffefficent of interest** [model_3mappersTxG.R](https://github.com/PaulKnoops/Experimental_Evolution_Sequence_Repo/blob/master/Analysis_after_BAM_Scripts/model_3mappersTxG.R)

This script (choosing Treatment by Generation effect) keeps positions that are present in all three files (i.e position needs to be mapped three times)

**Rscript: Combine into one genome:** [combinemodelCHROMO.R](https://github.com/PaulKnoops/Experimental_Evolution_Sequence_Repo/blob/master/Analysis_after_BAM_Scripts/combinemodelCHROMO.R)

**Rscript: P.adjust:** [model_p.adjustFDR.R](https://github.com/PaulKnoops/Experimental_Evolution_Sequence_Repo/blob/master/Analysis_after_BAM_Scripts/model_p.adjustFDR.R)

Adjust the p-values found for multiple comparisons: adjusting with FDR 

The [output](https://github.com/PaulKnoops/Experimental_Evolution_Sequence_Repo/tree/master/Data/model_output_p.Adjust) results in a file containing all chromosomes with p.values (and log10(p)) for each significant position after adjustments. (Note: Bonferroni adjustments kept as well as smaller file to practice further analysis with).



**Plots**

Treatment x Generation -log10(meanP-value) for model output: NEED TO CHANGE FOR P.adjust!
![FullGenomeTxGPlot](https://github.com/PaulKnoops/episodicSequenceData/blob/master/Analysis_after_sync_2018_plots/CHROMO_meanP.png)

_______________________________________________________________________________________

## 4) estimates of selection coefficient at each position for selection and control lineages using [poolSeq](https://github.com/ThomasTaus/poolSeq) R package:

**Script:** [poolseq_SelectionCoefficientEstimate.sh](https://github.com/PaulKnoops/Experimental_Evolution_Sequence_Repo/blob/master/Analysis_after_BAM_Scripts/poolseq_SelectionCoefficientEstimate.sh)

This script will break the sync files into two treatment .sync files, break apart these .sync files (smaller sized files), and run through a R script to run poolSeq Package (poolSeq_selectonCoeff.R)

**Rscript: Running poolseq** [poolSeq_selectionCoeff.R](https://github.com/PaulKnoops/Experimental_Evolution_Sequence_Repo/blob/master/Analysis_after_BAM_Scripts/poolSeq_selectionCoeff.R)

Note: to run, check poolseq is available, if not, source all PoolSeq scripts available from Taus git page. 

Also, to run with modified Sync files which changes the spacing,  a personal read.sync function is needed: [read.sync_personal_function.R](https://github.com/PaulKnoops/Experimental_Evolution_Sequence_Repo/blob/master/Analysis_after_BAM_Scripts/read.sync_personal_function.R)

This function is taken from the poolseq scripts from [poolSeq](https://github.com/ThomasTaus/poolSeq) with slight modifications.

**Rscript: combining CSV files:** [combinePoolseqCSV.R](https://github.com/PaulKnoops/Experimental_Evolution_Sequence_Repo/blob/master/Analysis_after_BAM_Scripts/combinePoolseqCSV.R)

(I was impatiant and did this individually: ex. [combine_poolseq_individual_Chromo.R](https://github.com/PaulKnoops/episodicSequenceData/blob/master/Analysis_after_sync_2018_scripts/combine_poolseq_individual_Chromo.R)
 
**Rscript: To combine the mapppers:** [poolseq_combinemappers.R](https://github.com/PaulKnoops/Experimental_Evolution_Sequence_Repo/blob/master/Analysis_after_BAM_Scripts/poolseq_combinemappers.R)

This will combine mappers (only Novoalign and BWA mem) and keep the mean selection coefficient (mean diff <0.001) and the less significant p-value (i.e max P-value). (does write a second .csv used for positions in plotting)


Plot: control vs. selection lines selection coefficients. Red == control and other colours correspond to chromosomes of selection lines selection coefficients 

![poolseqoverlaycontrols](https://github.com/PaulKnoops/episodicSequenceData/blob/master/Analysis_after_sync_2018_plots/overlay_CON:SEL_poolseq.png)



Notes: 

- should be edited to make more efficient to run (not taking <3 days per section)

- poolSeq may not work on certain versions of R; but can bring in Taus poolSeq scripts in manually and source (done in this script)

- may need to ensure updated packages and install if necessary (i.e matrixStats_0.53.0 installed for this reason)

- breaking the .sync file causes changes in structure, so a modified read_sync function is used (in R script; Taus_ReadSync.R))

_______________________________________________________________________________________

## 5) Positions of interest for Fst, poolseq and model output (overlap)

Finding positions overlapping with the significant model output (after adjustments), the signifciant selection coefficients (only found in selection lines and not in controls) and that are found within a Fst window with a sufficiently high value.

**Selection Coeffients positions:** [positions_selCoef.R](https://github.com/PaulKnoops/Experimental_Evolution_Sequence_Repo/blob/master/Analysis_after_BAM_Scripts/positions_selCoef.R)
 
 **Fst Window Ranges:** [positions_FST.R](https://github.com/PaulKnoops/Experimental_Evolution_Sequence_Repo/blob/master/Analysis_after_BAM_Scripts/positions_FST.R)
 
 **Positions from model:** [positions_Model.R](https://github.com/PaulKnoops/Experimental_Evolution_Sequence_Repo/blob/master/Analysis_after_BAM_Scripts/positions_Model.R)
 
 The above three scripts are sourced in the position extract script
 
 **Extract positions:** [positions_Extract.R](https://github.com/PaulKnoops/Experimental_Evolution_Sequence_Repo/blob/master/Analysis_after_BAM_Scripts/positions_Extract.R)
 
 The positions that are found in both the model output and with a significant selection coefficient are first found, then checked if they are present within the 500 bp window from FST.
 
 Ends with a number of positons for each chromosome: have this for keeping the mean p-value from selcoef and model as well as less signifciant (max) p-value.
 

 ```
 Number of positions
  Chr - Max - Mean
   2L - 71 - 80
   2R - 51 - 52
   3L - 19 - 19
   3R - 23 - 24
   4  - 0  - 0
   X  - 73 - 86
 ```
_______________________________________________________________________________________

## 6) Running Gowinda for gene analysis from positions

### To Run Gowinda: [Gowinda source forge tutorial](https://sourceforge.net/p/gowinda/wiki/Tutorial/)

1) **GTF file**
Need a Gtf file: converted from a gff file from [FlyBase homepage](http://flybase.org/) *dmel-all-r5.57.gff.gz* (matching current index)

Converted gff to gtf with Gowinda script Gff2Gtf.py:

Ex.
```
python /home/paul/Gowinda/Gff2Gtf.py --input /home/paul/episodicData/index_dir/dmel-all-r5.57.gff > /home/paul/Gowinda/dmel-all-r5.57.gtf
```

2) **Gene Association File**

Gene associations with [FuncAssociate](http://llama.mshri.on.ca/funcassociate/download_go_associations)

3) **Candidate Positions** found [here](https://github.com/PaulKnoops/Experimental_Evolution_Sequence_Repo/blob/master/Analysis_after_BAM_Scripts/positions_Extract.R))

Need to write .txt like this:
```
write.table(FILE, file='candidatePos.txt', sep ="\t", col.names = F, row.names = F)
```
4) **Full genome and positions**

```
cat novo_episodic.sync  | awk '{print $1,$2}' > /home/paul/Gowinda/positions_1.txt
```
Read into R and re-write in tab deliminated format
```
Xc <- fread('positions_1.txt')
write.table(Xc, file='positions.txt', sep ="\t", col.names = F, row.names = F)
```

5) **Running Gowinda**

Note: already FDR adjust!

For details: See Gowinda source forge tutorial Example 1: Basic Example
```
java -Xmx32g -jar /home/paul/Gowinda/Gowinda-1.12.jar --snp-file /home/paul/episodicData/novoalign/novo_mpileup/novo_episodic.sync --candidate-snp-file /home/paul/Gowinda/candidatePos.txt --gene-set-file /home/paul/Gowinda/funcassociate_go_associations.txt --annotation-file /home/paul/Gowinda/dmel-all-r5.57.gtf --simulations 100000 --min-significance 1 --gene-definition gene --threads 8 --output-file results_gene_gene.txt --mode gene --min-genes 1
```

For details: See Gowinda source forge tutorial Example 3: high resolution GO term enrichment

```
java -Xmx32g -jar /home/paul/Gowinda/Gowinda-1.12.jar \
	--snp-file /home/paul/episodicData/novoalign/novo_mpileup/novo_episodic.sync \
	--candidate-snp-file /home/paul/Gowinda/candidatePos.txt \
	--gene-set-file /home/paul/Gowinda/funcassociate_go_associations.txt \
	--annotation-file /home/paul/Gowinda/dmel-all-r5.57.gtf \
	--simulations 100000 \
	--min-significance 1 \
	--gene-definition updownstream2000 \
	--threads 8 \
	--output-file /home/paul/Gowinda/results_snp_2000ud.txt \
	--mode snp \
	--min-genes 1
```

Found .... nothing?

_______________________________________________________________________________________
## 7) Trajectory of regions of interest based on model output

Take positions of each chromosome and take a .sync file to only keep the positions of interest. 

With BWA -mem output sync files

**Rscript:** [extract_sig_Chromo_positions.R](https://github.com/PaulKnoops/Experimental_Evolution_Sequence_Repo/blob/master/Analysis_after_BAM_Scripts/extract_sig_Chromo_positions.R)

Use these .sync files and can look at trajectories:

```
## Packages source code:
require('tidyverse')
require('tidyr')
require('dplyr')
require('data.table')

## Read in data: 
#episodic_data <- fread('bwa_2L_positions.sync', header = T)
#episodic_data <- fread('bwa_2R_positions.sync', header = T)
#episodic_data <- fread('bwa_3L_positions.sync', header = T)
#episodic_data <- fread('bwa_3R_positions.sync', header = T)
#episodic_data <- fread('bwa_4_positions.sync', header = T) #empty
episodic_data <- fread('bwa_X_positions.sync', header = T)
Chromosome <- 'X'

## Below is basically the same as script sync_to_counts.R
#episodic_data <- episodic_counts
episodic_data  <- episodic_data[,c(2:17)]
name.Columns <- c("Chromosome", "Position", "ref", "ConR1_115", "ConR2_115", "SelR1_115", "SelR2_115", "ConR1_38", "ConR2_38", "SelR1_38", "SelR2_38", "ConR1_77", "ConR2_77", "SelR1_77", "SelR2_77", "SelR1_0")
colnames(episodic_data) <- name.Columns

## Generic script for sync_to_counts basically:
#Add "replicates" of ancestor -- all are equal
episodic_data$SelR2_0 <- episodic_data$SelR1_0
episodic_data$ConR1_0 <- episodic_data$SelR1_0
episodic_data$ConR2_0 <- episodic_data$SelR1_0

#Need the ancestor to stay (after making long) to call major/minor alleles later
episodic_data$Ancestor <- episodic_data$SelR1_0

# Make long by bring populations down
long_episodic <- gather(episodic_data, Population, Allele_Freq , ConR1_115:ConR2_0, factor_key=TRUE)

#rm(episodic_data)

Episodic_split_2 <- long_episodic %>% 
  separate(Allele_Freq, c("A","T","C","G","N","del"), ":")

#rm(long_episodic)
#Seperate the ancestor to base later things on
Episodic_split_2 <- Episodic_split_2 %>% 
  separate(Ancestor, c("A_0","T_0","C_0","G_0","N_0","del_0"), ":")

# as.numeric to multiple columns:
cols.num <- c("A_0", "T_0", "C_0", "G_0", "N_0", "del_0", "A", "T", "C", "G", "N", "del")

#Seems to take a long time for this step?
Episodic_split_2[cols.num] <- sapply(Episodic_split_2[cols.num],as.numeric)

#Get the sum of all the rows (all the different bases) for each population position:

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

# To determine the minor allele base if not specified by the ancestor (new allele brough up etc.)

#max for the population (could be the minor allele)
Episodic_split_2$maj_all <- apply(Episodic_split_2[,11:16], 1, max)

#alt== second highest count for populations
Episodic_split_2$alt_allele <- apply(Episodic_split_2[,11:16], 1, 
                                     function(x)max(x[x!=max(x)]))

Episodic_split_2 <- within(Episodic_split_2, {
  Min_count_2 = ifelse (Maj_count == sum, 0, ifelse(Maj_count==maj_all, alt_allele, maj_all))})

Episodic_split_2 <- within(Episodic_split_2, {
  MinorAllele_base = ifelse(Min_count_2==0, "N/A", ifelse(Min_count_2== Episodic_split_2[,11], "A", ifelse(Min_count_2== Episodic_split_2[,12],  "T", ifelse(Min_count_2== Episodic_split_2[,13],  "C",ifelse(Min_count_2== Episodic_split_2[,14],  "G", ifelse(Min_count_2== Episodic_split_2[,15],  "N", ifelse(Min_count_2== Episodic_split_2[,16],  "del", "Z") ))))))})

#Episodic_split_2 <- subset(Episodic_split_2, select = -c(A_0,T_0,C_0,G_0,N_0,del_0,A,T,C,G,N,del,anc_max,anc_min, MinorAllele, Min_count, maj_all, alt_allele))
#keeping ancestor info
Episodic_split_2 <- subset(Episodic_split_2, select = -c(A_0,T_0,C_0,G_0,N_0,del_0,A,T,C,G,N,del, MinorAllele, Min_count, maj_all, alt_allele))
nam.col <- c("chr", "pos", "ref", "Population", "sum", 'AncMax', 'AncMin', "MajorAllele", "Major_count", "Minor_count", "MinorAllele")

#nam.col <- c("chr", "pos", "ref", "Population", "sum", "MajorAllele", "Major_count", "Minor_count", "MinorAllele")

colnames(Episodic_split_2) <- nam.col

#Split Population into Treatment, Rep, and Generation - need to do twice, different seperators (change above??)

episodic_long <- Episodic_split_2 %>%
  separate(Population, c("Treatment", "Generation"), "_")

#rm(Episodic_split_2)

cols.num <- c("Generation", "Major_count", "Minor_count", 'sum')
episodic_long[cols.num] <- sapply(episodic_long[cols.num],as.numeric) 

head(episodic_long)
episodic_long$majFreq <- as.numeric(episodic_long$Major_count/episodic_long$sum)
episodic_long$minFreq <- as.numeric(episodic_long$Minor_count/episodic_long$sum)
episodic_long$AncMinFreq <- as.numeric(episodic_long$AncMin/(episodic_long$AncMax + episodic_long$AncMin))

head(episodic_long)
mean(na.omit(episodic_long$majFreq))
episodic_group <- episodic_long %>%
  group_by(chr, Treatment, Generation) %>%
  summarise(mean_minFreq=mean(na.omit(minFreq)), 
            mean_majFreq=mean(na.omit(majFreq)))

#Plotsing trajectories:
###--###########
xc_mean <- mean(episodic_group$mean_minFreq)
ggplot(data = episodic_group, aes(x=Generation, y=mean_minFreq, color=Treatment)) +  geom_point(size=2, alpha = 0.95, position=position_dodge(width=10)) + scale_colour_manual(values=c("#56B4E9", "#E69F00", 'grey30', 'firebrick3')) #+ ylim(0,1) #+ geom_hline(yintercept = xc_mean)

#diff:
episodic_group_2 <- episodic_long %>%
  group_by(chr, Treatment, Generation) %>%
  summarise(mean_minFreq=mean(na.omit(minFreq)), 
            mean_majFreq=mean(na.omit(majFreq)), 
            ddif_min = mean(abs(AncMinFreq-minFreq)))

ggplot(data = episodic_group_2, aes(x=Generation, y=ddif_min, color=Treatment)) +  geom_point(size=2, alpha = 0.95, position=position_dodge(width=10)) + scale_colour_manual(values=c("#56B4E9", "#E69F00", 'grey30', 'firebrick3'))

head(episodic_group_2)
ggplot(data = episodic_group_2, aes(x=as.factor(Generation), y=ddif_min, color=Treatment, group=Treatment)) + geom_line() + geom_point(size=2, alpha = 0.95, position=position_dodge(width=0)) + scale_colour_manual(values=c("#56B4E9", "#E69F00", 'grey30', 'firebrick3')) + ggtitle(Chromosome) + xlab('Generation')

#plot positions (either random position (randPos) or signify a position (pos)
randPos <- sample(episodic_long$pos, 1)
episodic_pos <-episodic_long[ which(episodic_long$pos==randPos), ]
#pos <- 4961389
#episodic_pos <- episodic_long[which(episodic_long$pos==pos),]
ggplot(data = episodic_pos, aes(x=as.factor(Generation), y=minFreq, color=Treatment, group=Treatment, shape=MinorAllele)) + 
  scale_colour_manual(values=c("#56B4E9", "#E69F00", 'grey30', 'firebrick3')) +
  geom_line() + xlab("Generation") +
  geom_point(size=3, alpha=0.5, position = position_dodge(width=0.1))

```

