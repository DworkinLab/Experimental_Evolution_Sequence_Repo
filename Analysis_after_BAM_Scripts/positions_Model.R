### Finding Positions of interest from Model:
##-----------------------------------------------##
### Packages:
require(dplyr)
require(ggplot2)
require(data.table)

### Read in each Chromosomal Data: 
# Created from https://github.com/PaulKnoops/Experimental_Evolution_Sequence_Repo/blob/master/Analysis_after_BAM_Scripts/combinemodelCHROMO.R

CHROMOs <- fread('CHROMO_model_Full.csv')

###----------------------------------------##
#CHROMOs$adjustP <- p.adjust(CHROMOs$max_p, method = 'fdr')
#CHROMOs$adjustP <- p.adjust(CHROMOs$max_p, method = 'bonferroni')

CHROMOs$adjustP <- p.adjust(CHROMOs$mean_p, method = 'fdr')
#CHROMOs$adjustP <- p.adjust(CHROMOs$mean_p, method = 'bonferroni')


CHROMOs$logAdjustP <- -log10(CHROMOs$adjustP)
CHROMOs_4 <-  CHROMOs[-which(CHROMOs$adjustP>0.05),]
CHROMOs_3 <-  CHROMOs_4[-which(CHROMOs_4$logAdjustP==Inf),]

# split based on chromosome and significant positions:
# need to sort positions based on chromosome (for trajectories):
ddat2_X <- CHROMOs_3[which(CHROMOs_3$chr=='X' & CHROMOs_3$adjustP<0.05),]
ddat2_2L <- CHROMOs_3[which(CHROMOs_3$chr=='2L' & CHROMOs_3$adjustP<0.05),]
ddat2_2R <- CHROMOs_3[which(CHROMOs_3$chr=='2R' & CHROMOs_3$adjustP<0.05),]
ddat2_3L <- CHROMOs_3[which(CHROMOs_3$chr=='3L' & CHROMOs_3$adjustP<0.05),]
ddat2_3R <- CHROMOs_3[which(CHROMOs_3$chr=='3R' & CHROMOs_3$adjustP<0.05),]
ddat2_4 <- CHROMOs_3[which(CHROMOs_3$chr=='4' & CHROMOs_3$adjustP<0.05),]

# create and write a positional .csv file
pos_X <- ddat2_X$position
pos_2L <- ddat2_2L$position
pos_2R <- ddat2_2R$position
pos_3L <- ddat2_3L$position
pos_3R <- ddat2_3R$position
pos_4 <- ddat2_4$position


