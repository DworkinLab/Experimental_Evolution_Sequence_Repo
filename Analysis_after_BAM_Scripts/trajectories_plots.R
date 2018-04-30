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

# Change based on the chromosome:
Chromosome <- 'X'

## Below is basically the same as script sync_to_counts.R b/w #####-----####
#####-------------------------------------------------------####
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

#####-------------------------------------------------------####

#Ploting trajectories:
#####-------------------------------------------------------####

# The mean minor frequency of chromosome (based on Ancestor) -- not super informative as freq can rise and fall.
xc_mean <- mean(episodic_group$mean_minFreq)
ggplot(data = episodic_group, aes(x=Generation, y=mean_minFreq, color=Treatment)) +  geom_point(size=2, alpha = 0.95, position=position_dodge(width=10)) + scale_colour_manual(values=c("#56B4E9", "#E69F00", 'grey30', 'firebrick3')) #+ ylim(0,1) #+ geom_hline(yintercept = xc_mean)

# Next two = absolute difference in frequency from ancestor.

#diff:
episodic_group_2 <- episodic_long %>%
  group_by(chr, Treatment, Generation) %>%
  summarise(mean_minFreq=mean(na.omit(minFreq)), 
            mean_majFreq=mean(na.omit(majFreq)), 
            ddif_min = mean(abs(AncMinFreq-minFreq)))

ggplot(data = episodic_group_2, aes(x=Generation, y=ddif_min, color=Treatment)) +  geom_point(size=2, alpha = 0.95, position=position_dodge(width=10)) + scale_colour_manual(values=c("#56B4E9", "#E69F00", 'grey30', 'firebrick3'))

head(episodic_group_2)
ggplot(data = episodic_group_2, aes(x=as.factor(Generation), y=ddif_min, color=Treatment, group=Treatment)) + geom_line() + geom_point(size=2, alpha = 0.95, position=position_dodge(width=0)) + scale_colour_manual(values=c("#56B4E9", "#E69F00", 'grey30', 'firebrick3')) + ggtitle(Chromosome) + xlab('Generation')



# Plot positions (either random position (randPos) or signify a position (pos)
randPos <- sample(episodic_long$pos, 1)
episodic_pos <-episodic_long[ which(episodic_long$pos==randPos), ]
#pos <- 4961389
#episodic_pos <- episodic_long[which(episodic_long$pos==pos),]
ggplot(data = episodic_pos, aes(x=as.factor(Generation), y=minFreq, color=Treatment, group=Treatment, shape=MinorAllele)) + 
  scale_colour_manual(values=c("#56B4E9", "#E69F00", 'grey30', 'firebrick3')) +
  geom_line() + xlab("Generation") +
  ggtitle(paste("Chromosome=", Chromosome, '  ', "Position=", randPos, sep = "")) +
  geom_point(size=3, alpha=0.5, position = position_dodge(width=0.1))
