# To be sourced for positions:

setwd('/Users/paulknoops/Bioinformatics/Analysis_EpisodicSequenceData_2018/poolSeq')
require(tidyverse)
require(data.table)
#Read in data:
## source(poolseq_CombineMappers.R)

#Created with https://github.com/PaulKnoops/Experimental_Evolution_Sequence_Repo/blob/master/Analysis_after_BAM_Scripts/poolseq_combinemappers.R

Zxc <- fread('SelCoef_Full.csv')

# adjust p values
#Zxc$adjustP <- p.adjust(Zxc$pval_max, method = 'fdr')
#Zxc$adjustP <- p.adjust(Zxc$pval_max, method = 'bonferroni')

Zxc$adjustP <- p.adjust(Zxc$pval_mean, method = 'fdr')
#Zxc$adjustP <- p.adjust(Zxc$pval_mean, method = 'bonferroni')

# Lable for significance 
Zxc$sig <- ifelse(Zxc$adjustP<0.05, "<0.05", ">0.05")
 
# Only significant positions:
Zxc_sig <- Zxc[which(Zxc$sig=='<0.05'),]
# 
Zxc_count <- Zxc_sig %>%
  group_by(chr, pos) %>%
  mutate(count = n())
Zxc_count2 <- Zxc_count[which(Zxc_count$count==1),]
Zxc_count2 <- Zxc_count2[which(Zxc_count2$Treatment=='Sel'),]

#positions: under selection (not removing any control seleted regions?

poolseq_X <- Zxc_count2[which(Zxc_count2$chr=='X' & Zxc_count2$Treatment=='Sel'),]
poolseq_2L <- Zxc_count2[which(Zxc_count2$chr=='2L' & Zxc_count2$Treatment=='Sel'),]
poolseq_2R <- Zxc_count2[which(Zxc_count2$chr=='2R' & Zxc_count2$Treatment=='Sel'),]
poolseq_3L <- Zxc_count2[which(Zxc_count2$chr=='3L' & Zxc_count2$Treatment=='Sel'),]
poolseq_3R <- Zxc_count2[which(Zxc_count2$chr=='3R' & Zxc_count2$Treatment=='Sel'),]
poolseq_4 <- Zxc_count2[which(Zxc_count2$chr=='4' & Zxc_count2$Treatment=='Sel'),]


# create and write a positional .csv file
pool_pos_X <- poolseq_X$pos
pool_pos_2L <- poolseq_2L$pos
pool_pos_2R <- poolseq_2R$pos
pool_pos_3L <- poolseq_3L$pos
pool_pos_3R <- poolseq_3R$pos
pool_pos_4 <- poolseq_4$pos
pool_pos_X <- poolseq_X$pos

