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

Source [Pi_plot_function.R](https://github.com/PaulKnoops/Experimental_Evolution_Sequence_Repo/blob/master/Analysis_after_BAM_Scripts/Pi_plot_function.R) and run files located in [data](https://github.com/PaulKnoops/Experimental_Evolution_Sequence_Repo/tree/master/Data) directory.

Ex. (currently have title details removed in function)
```
Pi_PlotFunction('FILE.pi', "Plot Title Details")
```

Example outputs:

**Ancestral Pi Plot with Bowtie2 mapping**
![BowtiePiPlot Ancestor](https://github.com/PaulKnoops/Experimental_Evolution_Sequence_Repo/blob/master/plots/PiPlot_Ancestor_Bowtie.png)

**Ancestral Pi Plot with Novoalign mapping**
![Novo_PiPlot_anc](https://github.com/PaulKnoops/Experimental_Evolution_Sequence_Repo/blob/master/plots/PiPlot_Ancestor_novoalign.png)
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

***Script:** [novo_Fst.sh](https://github.com/PaulKnoops/Experimental_Evolution_Sequence_Repo/blob/master/Analysis_after_BAM_Scripts/novo_Fst.sh)*

ex.
```
perl ${fst} --input ${novo_mpileup}/novo_episodic_main.sync --output ${novo_fst}/novo_episodic_main.fst --min-count 6 --min-coverage 10 --max-coverage 250 --min-covered-fraction 1 --window-size 500 --step-size 500 --pool-size 120
```

### In R, split the file into each compasison

***Script:** [novo_Fst_Split_Comparisons.R](https://github.com/PaulKnoops/Experimental_Evolution_Sequence_Repo/blob/master/Analysis_after_BAM_Scripts/novo_Fst_Split_Comparisons.R)*

R script that will split the .fst file into many .csv files with each comparison (can choose the necessary ones from here)

### Combining three mappers output:

***Script:** [FST_combine3mappers.R](https://github.com/PaulKnoops/Experimental_Evolution_Sequence_Repo/blob/master/Analysis_after_BAM_Scripts/FST_combine3mappers.R)*

Will take the split comparisons, and combine those specified into one FST file with the average Fst between the three mappers 

By changing the "patty" variable, the comparisons of interest can be combined and evaluated by taking the matching comparisons from three directories holding the different mappers output from previous script spliting FST output.

Some combined data files of interest can be found here: [Fst_combinedComparisons](https://github.com/PaulKnoops/Experimental_Evolution_Sequence_Repo/tree/master/Data/Fst_combinedComparisons)

### Plotting Fst files for Con:Sel comparison at three generations

 - average Fst of three mappers **and** average between replicates
 
 - comparison betweeen control and selection lines
 
***Script:** [Fst_Plots.R](https://github.com/PaulKnoops/Experimental_Evolution_Sequence_Repo/blob/master/Analysis_after_BAM_Scripts/Fst_Plots.R)*
 
 Can create plots by running script with data available from [Fst Data Repo](https://github.com/PaulKnoops/Experimental_Evolution_Sequence_Repo/tree/master/Data/Fst_combinedComparisons)
 
Generation 38:
![meanFst for F38](https://github.com/PaulKnoops/Experimental_Evolution_Sequence_Repo/blob/master/plots/meanFst_gen38_Sel:Con.png)

Generation 77: 
![meanFst for F77](https://github.com/PaulKnoops/Experimental_Evolution_Sequence_Repo/blob/master/plots/meanFst_gen77_Sel:Con.png)

Generation 115: 
![meanFst for F115](https://github.com/PaulKnoops/Experimental_Evolution_Sequence_Repo/blob/master/plots/meanFst_gen115_Sel:Con.png)

### Additional Plots:

If necessary: method to create sauron plots

Also includes details on quantiles used for positional cut offs later

**Sauron Plots:** [Fst_SauronPlots.R](https://github.com/PaulKnoops/Experimental_Evolution_Sequence_Repo/blob/master/Analysis_after_BAM_Scripts/Fst_SauronPlots.R)

_______________________________________________________________________________________

## 3) per SNP logistic regression for each treatment by generation

***Long Script:*** [novo_regression_model_LONGSCRIPT.sh](https://github.com/PaulKnoops/Experimental_Evolution_Sequence_Repo/blob/master/Analysis_after_BAM_Scripts/novo_regression_model_LONGSCRIPT.sh)*

This script will break the chromosomal .sync files (i.e split per chromosome) into smaller managable pieces and run through multiple R scripts while removing intermediates:

The R script below are within the long script:

***R script to covert sync to Count data:** [Sync_to_counts.R](https://github.com/PaulKnoops/Experimental_Evolution_Sequence_Repo/blob/master/Analysis_after_BAM_Scripts/Sync_to_counts.R)*

Creates a file with the counts for the major and minor frequency (based on ancestor) that can run through the model

***R script for running the model for each position along the chromosome:** [Counts_to_model.R](https://github.com/PaulKnoops/Experimental_Evolution_Sequence_Repo/blob/master/Analysis_after_BAM_Scripts/Counts_to_model.R)*

In long script: this is set up to work in parallel, having each chromosome running at the same time (6 instances running over 11 sections)

NOTE: Script needs to be changed to run faster/ more efficiently (not done here b/c already completed)

Basic Model at each positon (tmp2):
```
modlist_2[[i]] <- 
        glm(cbind(Major_count, Minor_count) ~ Treatment*Generation, 
            data = tmp2, family = "binomial")
```


**After** long script complete:

***R script to combine all the split chromosome pieces back into one chromosome:** [Combine_model_Chromo.R](https://github.com/PaulKnoops/Experimental_Evolution_Sequence_Repo/blob/master/Analysis_after_BAM_Scripts/Combine_model_Chromo.R)*


Recreates one chromosomal file

***R script to combine three mappers into one file** [model_combine3mappers.R](https://github.com/PaulKnoops/Experimental_Evolution_Sequence_Repo/blob/master/Analysis_after_BAM_Scripts/model_combine3mappers.R)*

Combines each of BWA-mem, Bowtie2 and Novoalign files into one file (keeping all information)


***R script to write files with coeffefficent of interest** [model_3mappersTxG.R](https://github.com/PaulKnoops/Experimental_Evolution_Sequence_Repo/blob/master/Analysis_after_BAM_Scripts/model_3mappersTxG.R)*

This script (choosing Treatment by Generation effect) keeps positions that are present in all three files (i.e position needs to be mapped three times)

***Rscript: Combine into one genome:** [combinemodelCHROMO.R](https://github.com/PaulKnoops/Experimental_Evolution_Sequence_Repo/blob/master/Analysis_after_BAM_Scripts/combinemodelCHROMO.R)*

The output results in a file containing all chromosomes together with meanPvalue and maxPvalue (less sig. pos), but is to large for Github (available elsewhere).

In the [output](https://github.com/PaulKnoops/Experimental_Evolution_Sequence_Repo/tree/master/Data/model_output_p.Adjust) directory is the chromosomes after FDR adjustments and a number file for plotting (and with Bonferroni for smaller easier to test out file)

***Rscript: P.adjust:** [model_p.adjustFDR.R](https://github.com/PaulKnoops/Experimental_Evolution_Sequence_Repo/blob/master/Analysis_after_BAM_Scripts/model_p.adjustFDR.R)*

Adjust the p-values found for multiple comparisons: adjusting with FDR (bonferroni and other p-value method available in script)


**Plots**

Treatment x Generation -log10(meanP-value) for model output with FDR Adjustments
![FullGenomeTxGPlot](https://github.com/PaulKnoops/Experimental_Evolution_Sequence_Repo/blob/master/plots/Model_FDR_p.adjustPlot.png)


Treatment x Generation -log10(meanP-value) for model output with Bonferroni Adjustments
![FullGenomeTxGPlot_Bonf](https://github.com/PaulKnoops/Experimental_Evolution_Sequence_Repo/blob/master/plots/Model_Bonf_p.adjust_Plot.png)

_______________________________________________________________________________________

## 4) estimates of selection coefficient at each position for selection and control lineages using [poolSeq](https://github.com/ThomasTaus/poolSeq) R package:

***Script:** [poolseq_SelectionCoefficientEstimate.sh](https://github.com/PaulKnoops/Experimental_Evolution_Sequence_Repo/blob/master/Analysis_after_BAM_Scripts/poolseq_SelectionCoefficientEstimate.sh)*

This script will break the sync files into two treatment .sync files, break apart these .sync files (smaller sized files), and run through a R script to run poolSeq Package (poolSeq_selectonCoeff.R)

***Rscript: Running poolseq** [poolSeq_selectionCoeff.R](https://github.com/PaulKnoops/Experimental_Evolution_Sequence_Repo/blob/master/Analysis_after_BAM_Scripts/poolSeq_selectionCoeff.R)*

Note: to run, check poolseq is available, if not, source all PoolSeq scripts available from Taus git page. 

Also, to run with modified Sync files which changes the spacing,  a personal read.sync function is needed: [read.sync_personal_function.R](https://github.com/PaulKnoops/Experimental_Evolution_Sequence_Repo/blob/master/Analysis_after_BAM_Scripts/read.sync_personal_function.R)

This function is taken from the poolseq scripts from [poolSeq](https://github.com/ThomasTaus/poolSeq) with slight modifications.

***Rscript: combining CSV files:** [combinePoolseqCSV.R](https://github.com/PaulKnoops/Experimental_Evolution_Sequence_Repo/blob/master/Analysis_after_BAM_Scripts/combinePoolseqCSV.R)*

(I was impatiant and did this individually: ex. [combine_poolseq_individual_Chromo.R](https://github.com/PaulKnoops/episodicSequenceData/blob/master/Analysis_after_sync_2018_scripts/combine_poolseq_individual_Chromo.R)
 
***Rscript: To combine the mapppers:** [poolseq_combinemappers.R](https://github.com/PaulKnoops/Experimental_Evolution_Sequence_Repo/blob/master/Analysis_after_BAM_Scripts/poolseq_combinemappers.R)*

This will combine mappers (only Novoalign and BWA mem) and keep the mean selection coefficient (mean difference between mappers < 0.001) and the less significant p-value (i.e max P-value). (does write a second .csv used for positions in plotting).


The output is to large for github, so the beginning of [poolseq_PlotSelCoef.R](https://github.com/PaulKnoops/Experimental_Evolution_Sequence_Repo/blob/master/Analysis_after_BAM_Scripts/poolseq_PlotSelCoef.R) script shows how this was broken down into files based on the p.adjust method (found in [Data](https://github.com/PaulKnoops/Experimental_Evolution_Sequence_Repo/tree/master/Data/poolseq_outputs) directory) as well as lengths for plotting.

To plot the data: [poolseq_PlotSelCoef.R](https://github.com/PaulKnoops/Experimental_Evolution_Sequence_Repo/blob/master/Analysis_after_BAM_Scripts/poolseq_PlotSelCoef.R)

Plot: control vs. selection lines selection coefficients. Red == control and other colours correspond to chromosomes of selection lines selection coefficients 

![poolseqoverlaycontrols](https://github.com/PaulKnoops/Experimental_Evolution_Sequence_Repo/blob/master/plots/poolseq_SelCoef_overlay.png)

Notes: 

- should be edited to make more efficient to run (not taking <3 days per section)

- poolSeq may not work on certain versions of R; but can bring in Taus poolSeq scripts in manually and source (done in this script)

- may need to ensure updated packages and install if necessary (i.e matrixStats_0.53.0 installed for this reason)

- breaking the .sync file causes changes in structure, so a modified read_sync function is used (in R script; Taus_ReadSync.R))

_______________________________________________________________________________________

## 5) Positions of interest for Fst, poolseq and model output (overlap)

Finding positions overlapping with the significant model output (after adjustments), the signifciant selection coefficients (only found in selection lines and not in controls) and that are found within a Fst window with a sufficiently high value.

***Selection Coeffients positions:** [positions_selCoef.R](https://github.com/PaulKnoops/Experimental_Evolution_Sequence_Repo/blob/master/Analysis_after_BAM_Scripts/positions_selCoef.R)*
 
 ***Fst Window Ranges:** [positions_FST.R](https://github.com/PaulKnoops/Experimental_Evolution_Sequence_Repo/blob/master/Analysis_after_BAM_Scripts/positions_FST.R)*
 
 ***Positions from model:** [positions_Model.R](https://github.com/PaulKnoops/Experimental_Evolution_Sequence_Repo/blob/master/Analysis_after_BAM_Scripts/positions_Model.R)*
 
 The above three scripts are sourced in the position extract script
 
 ***Extract positions:** [positions_Extract.R](https://github.com/PaulKnoops/Experimental_Evolution_Sequence_Repo/blob/master/Analysis_after_BAM_Scripts/positions_Extract.R)*
 
 The positions that are found in both the model output and with a significant selection coefficient are first found, then checked if they are present within the 500 bp window from FST.
 
 Ends with a number of positons for each chromosome and for a full chromosome candidate positons [found here](https://github.com/PaulKnoops/Experimental_Evolution_Sequence_Repo/tree/master/Data/Positions): have this for less signifciant (max) p-value (used for both selcoef and model) as well as for the mean p. The overlap with max_p and mean_p is shown below: went with the more stringent method (maxP) for candidate positions. Can combine the chromosomal meanP.csv files to create all candidates with mean P values.

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

3) ***Candidate Positions** found [here](https://github.com/PaulKnoops/Experimental_Evolution_Sequence_Repo/blob/master/Analysis_after_BAM_Scripts/positions_candidates_Full.R)*

Need to write .txt like this: tab deliminated text file
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

Note: already FDR adjusts within Gowinda!

For details: See [Gowinda source forge tutorial](https://sourceforge.net/p/gowinda/wiki/Tutorial/) Example 1: Basic Example
```
java -Xmx32g -jar /home/paul/Gowinda/Gowinda-1.12.jar --snp-file --snp-file /home/paul/Gowinda/positions.txt --candidate-snp-file /home/paul/Gowinda/candidatePos.txt --gene-set-file /home/paul/Gowinda/funcassociate_go_associations.txt --annotation-file /home/paul/Gowinda/dmel-all-r5.57.gtf --simulations 100000 --min-significance 1 --gene-definition gene --threads 8 --output-file results_gene_gene.txt --mode gene --min-genes 1
```

For details: See [Gowinda source forge tutorial](https://sourceforge.net/p/gowinda/wiki/Tutorial/) Example 3: high resolution GO term enrichment

```
java -Xmx32g -jar /home/paul/Gowinda/Gowinda-1.12.jar \
	--snp-file /home/paul/Gowinda/positions.txt \
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

Found .... nothing for both.... ??? Need to rerun possibly or try SNPeff

End of output
```
04/30/18 13:54:24: Finished - Candidate SNPs are a subset of all SNPs
04/30/18 13:54:24: Computing the number of genes and gene sets available for testing
04/30/18 13:56:25: Total genes in annotation: 12243; Annotated genes with SNP: 0 (genes without SNP are ignored).
04/30/18 13:56:25: Total gene sets: 9998; Gene sets with at least one gene having a SNP: 0  (gene sets without a single gene having a SNP are ignored)
04/30/18 13:56:25: Starting 100000 simulations for 261 candidate SNPs using 8 threads
04/30/18 13:56:25: Will randomly draw SNPs unless the number of random SNPs equals the number of candidate SNPs
04/30/18 13:56:25: This may take a while. Switch to the detailed log mode if you want to see the progress
04/30/18 13:56:28: Finished simulations
04/30/18 13:56:28: Estimating significance of the candidate SNPs
04/30/18 13:56:28: Starting FDR correction
04/30/18 13:56:28: Simulations detected SNPs in genes corresponding to 0 gene sets
04/30/18 13:56:28: Candidate SNPs show an association with 0 gene sets
04/30/18 13:56:45: Starting to write results to file: /home/paul/Gowinda/results_snp_2000ud.txt
04/30/18 13:56:45: Finished writing to file
04/30/18 13:56:45: FINISHED - Thank you for using Gowinda
```

_______________________________________________________________________________________
## 7) Trajectory of regions of interest based on model output

Take positions of [each chromosome](https://github.com/PaulKnoops/Experimental_Evolution_Sequence_Repo/tree/master/Data/Positions) and take a .sync file to only keep the positions of interest. 


***Rscript:** [extract_sig_Chromo_positions.R](https://github.com/PaulKnoops/Experimental_Evolution_Sequence_Repo/blob/master/Analysis_after_BAM_Scripts/extract_sig_Chromo_positions.R)*

Use these .sync files and can look at trajectories


***Rscript:** [trajectories_plots.R](https://github.com/PaulKnoops/Experimental_Evolution_Sequence_Repo/blob/master/Analysis_after_BAM_Scripts/trajectories_plots.R)*

Script will output different varients of trajectory plots of positions.

If the different mappers want to be combined, can be done with this ***Rscript:**[trajectoryPlots_combineMappers.R](https://github.com/PaulKnoops/Experimental_Evolution_Sequence_Repo/blob/master/Analysis_after_BAM_Scripts/trajectoryPlots_combineMappers.R)*


**Mean Difference (of three mappers minor frequencies) from the mean ancestral frequency averaged for all candidate positions** (plots made in [above script](https://github.com/PaulKnoops/Experimental_Evolution_Sequence_Repo/blob/master/Analysis_after_BAM_Scripts/trajectoryPlots_combineMappers.R).

![Difference_2](https://github.com/PaulKnoops/Experimental_Evolution_Sequence_Repo/blob/master/plots/Trajectory_Chromo2_diff.png)

![Difference_3](https://github.com/PaulKnoops/Experimental_Evolution_Sequence_Repo/blob/master/plots/Trajectories_Chromo3_Diff.png)

![Difference_X](https://github.com/PaulKnoops/Experimental_Evolution_Sequence_Repo/blob/master/plots/Trajectories_X_Diff.png)

**Positional Trajectories**

The minor frequency changes over generations of a random sampled position for each chromosome: can recreate at bottom of [script](https://github.com/PaulKnoops/Experimental_Evolution_Sequence_Repo/blob/master/Analysis_after_BAM_Scripts/trajectoryPlots_combineMappers.R)

![4positionson4chromosomes](https://github.com/PaulKnoops/Experimental_Evolution_Sequence_Repo/blob/master/plots/Trajectories_Rand_pos_2_3.png)

