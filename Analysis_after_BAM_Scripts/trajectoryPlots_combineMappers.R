## Combine three mappers: trajectories
require(data.table)
require(tidyverse)
# Based on starting from scripts:
source('multiplotFunction.R')
setwd("../Data/positions_syncfiles")
dir <- getwd()
patty <- c('2L_positions.sync','2R_positions.sync','3L_positions.sync','3R_positions.sync','X_positions.sync')

#Loop if all want to be done quickly (better to do individually (within loop))
for (patt in patty){
  
  # To run individual Chromosomes: unhash wanted chromosome and hash loop start and end (#for (patt in patty){) & (})
  # patt <- '2L_positions.sync'
  # patt <- '2R_positions.sync'
  #patt <- '3L_positions.sync'
  # patt <- '3R_positions.sync'
  # patt <- 'X_positions.sync'
  
  Xcs <- patt
  Xcs2 <- gsub("\\..*","", Xcs)
  Xcs3 <- gsub("_.*","", Xcs2)
  Chromosome <- Xcs3
  
  mypos <- list.files(path = dir, pattern=patt, full.names = F)
  
  DFS <- NULL
  #Loop to combine them:
  for (posy in mypos) {
  episodic_data <- fread(posy, header = T)
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
  
  #keeping ancestor info
  Episodic_split_2 <- subset(Episodic_split_2, select = -c(A_0,T_0,C_0,G_0,N_0,del_0,A,T,C,G,N,del, MinorAllele, Min_count, maj_all, alt_allele))
  nam.col <- c("chr", "pos", "ref", "Population", "sum", 'AncMax', 'AncMin', "MajorAllele", "Major_count", "Minor_count", "MinorAllele")
  
  colnames(Episodic_split_2) <- nam.col
  
  #Split Population into Treatment, Rep, and Generation - need to do twice, different seperators (change above??)
  
  episodic_long <- Episodic_split_2 %>%
    separate(Population, c("Treatment", "Generation"), "_")
  
  
  cols.num <- c("Generation", "Major_count", "Minor_count", 'sum')
  episodic_long[cols.num] <- sapply(episodic_long[cols.num],as.numeric) 
  
  episodic_long$majFreq <- as.numeric(episodic_long$Major_count/episodic_long$sum)
  episodic_long$minFreq <- as.numeric(episodic_long$Minor_count/episodic_long$sum)
  episodic_long$AncMinFreq <- as.numeric(episodic_long$AncMin/(episodic_long$AncMax + episodic_long$AncMin))

  DFS <- rbind(DFS, episodic_long)
  
  }

  DFS_2 <- DFS %>%
    group_by(chr,pos,Treatment, Generation) %>%
    summarise(mean_ancminFreq=mean(na.omit(AncMinFreq)), 
              mean_minFreq=mean(na.omit(minFreq)), 
              mean_majFreq=mean(na.omit(majFreq)),
              minallele_1=(MinorAllele[1]), 
              minallele_2=(MinorAllele[2]),
              minallele_3=(MinorAllele[3])
              )
  
  episodic_group_2 <- DFS_2 %>%
    group_by(chr, Treatment, Generation) %>%
    summarise(mean_minFreq=mean(na.omit(mean_minFreq)), 
              mean_majFreq=mean(na.omit(mean_majFreq)), 
              ddif_min = mean(abs(mean_ancminFreq-mean_minFreq)))
  
  assign(paste("Chromosome", Chromosome, 'meanDiff', sep = "_"), episodic_group_2)
  assign(paste("Chromosome", Chromosome, 'positions', sep = "_"), DFS_2)
  }

# List all Data frames created:
myDIFF <- ls(pattern = "*meanDiff")
mypos <- ls(pattern='*positions')

  # 2nd plot is better version of this one:
  #Ppplt <- ggplot(data = episodic_group_2, aes(x=Generation, y=ddif_min, color=Treatment)) +  geom_point(size=2, alpha = 0.95, position=position_dodge(width=10)) + scale_colour_manual(values=c("#56B4E9", "#E69F00", 'grey30', 'firebrick3')) + ggtitle(Chromosome)
  #print(Ppplt)


for (Ddif in myDIFF) {
  ddat <- eval(parse(text = Ddif))
  Chromo <- ddat$chr[1]
 Ppplt_2 <-  ggplot(data = ddat, aes(x=as.factor(Generation), y=ddif_min, color=Treatment, group=Treatment)) + geom_line(show.legend = T) + geom_point(size=2, alpha = 0.95, position=position_dodge(width=0), show.legend = F) + scale_colour_manual(values=c("#56B4E9", "#E69F00", 'grey30', 'firebrick3')) + ggtitle(Chromosome) + xlab('Generation') +
   ylab('Mean Diff from Gen0') + theme(legend.position="bottom") +
   ggtitle(Chromo)
 
  print(Ppplt_2)
  assign(paste("Plot", Chromo, 'meanDiff', sep = "_"), Ppplt_2)
}

multiplot(Plot_2L_meanDiff, Plot_2R_meanDiff, Plot_3L_meanDiff, Plot_3R_meanDiff, cols=2)  

multiplot(Plot_2L_meanDiff, Plot_2R_meanDiff, cols=2 )
multiplot(Plot_3L_meanDiff, Plot_3R_meanDiff, cols=2)
Plot_empty <- NULL
multiplot(Plot_X_meanDiff, Plot_empty,  cols=2)

# Do indivdually:
poosition <- Chromosome_2L_positions
#poosition <- Chromosome_2R_positions
#poosition <- Chromosome_3L_positions
#poosition <- Chromosome_3R_positions
#poosition <- Chromosome_X_positions

  randPos <- sample(poosition$pos, 1)
  episodic_pos <-poosition[ which(poosition$pos==randPos), ]
  #pos <- 4961389
  #episodic_pos <- poosition[which(poosition$pos==pos),]
  pplot_3 <- ggplot(data = episodic_pos, aes(x=as.factor(Generation), y=mean_minFreq, color=Treatment, group=Treatment)) + 
    scale_colour_manual(values=c("#56B4E9", "#E69F00", 'grey30', 'firebrick3')) +
    ggtitle(paste("Chromosome=", Chromosome, '  ', "Position=", randPos, sep = "")) +
    geom_line() + xlab("Generation") +
    geom_point(size=3, alpha=0.5, position = position_dodge(width=0.1))
  
  print(pplot_3)



