
##-----------------------------------------------##
#Combines Chromosomes: while filtering mappers by keeping both mean pvalue and the less significant (max) p value
### Packages:

require(dplyr)
require(ggplot2)
require(data.table)

### Read in each Chromosomal Data:  
# Will read in chromosome
# based on position (x3), both the average p_value and the least significant value (max P) are calculated)

# wd = location of chromosome .csv files
setwd("/Users/paulknoops/Bioinformatics/Analysis_EpisodicSequenceData_2018/modelOutput_2018/Scripts")

ddatX <- fread('../Data/X_Chromosome_TxG.csv', h=T)
chr_X <- ddatX %>%
  group_by(position, chr, Effects) %>%
  summarise(mean_p = (mean(p.value)), max_p = (max(p.value)))
rm(ddatX)
a <- nrow(chr_X)
chr_X$number <- 1:a


ddat2L <- fread('../Data/2L_Chromosome_TxG.csv', h=T)
chr_2L <- ddat2L %>%
  group_by(position, chr, Effects) %>%
  summarise(mean_p = (mean(p.value)), max_p = (max(p.value)))
rm(ddat2L)
b <- nrow(chr_2L)
chr_2L$number <- (a+1):(a+b)


ddat2R <- fread('../Data/2R_Chromosome_TxG.csv', h=T)
chr_2R <- ddat2R %>%
  group_by(position, chr, Effects) %>%
  summarise(mean_p = (mean(p.value)), max_p = (max(p.value)))
rm(ddat2R)
c <- nrow(chr_2R)
chr_2R$number <- (a+b+1):(a+b+c)


ddat3L <- fread('../Data/3L_Chromosome_TxG.csv', h=T)
chr_3L <- ddat3L %>%
  group_by(position, chr, Effects) %>%
  summarise(mean_p = (mean(p.value)), max_p = (max(p.value)))
rm(ddat3L)
d <- nrow(chr_3L)
chr_3L$number <- (a+b+c+1):(a+b+c+d)

ddat3R <- fread('../Data/3R_Chromosome_TxG.csv', h=T)
chr_3R <- ddat3R %>%
  group_by(position, chr, Effects) %>%
  summarise(mean_p = (mean(p.value)), max_p = (max(p.value)))
rm(ddat3R)
e <- nrow(chr_3R)
chr_3R$number <- (a+b+c+d+1):(a+b+c+d+e)

ddat4 <- fread('../Data/4_Chromosome_TxG.csv', h=T)
chr_4 <- ddat4 %>%
  group_by(position, chr, Effects) %>%
  summarise(mean_p = (mean(p.value)), max_p = (max(p.value)))
rm(ddat4)
f <- nrow(chr_4)
chr_4$number <- (a+b+c+d+e+1):(a+b+c+d+e+f)
chr_4$chr <- as.character(chr_4$chr)

CHROMOs <- rbind(chr_X, chr_2L, chr_2R, chr_3L, chr_3R, chr_4)
CHROMOs <- rbind(chr_X, chr_2L, chr_2R, chr_3L, chr_3R, chr_4)
#rm(chr_X, chr_2L, chr_2R, chr_3L, chr_3R, chr_4)

#write.csv(CHROMOs, file = '../Data/CHROMO_model_Full.csv', row.names = F)

#Needed for plots:
CHROMO_NUMB <- c(a,b,c,d,e,f)
#write.csv(CHROMO_NUMB, file = '../Data/CHROMO_model_Full_LENGTH.csv', row.names = F)


