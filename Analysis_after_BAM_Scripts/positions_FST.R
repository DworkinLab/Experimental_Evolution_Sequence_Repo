# Fst: 
require(data.table)
require(tidyverse)

# Generation 115 comparisons: Sel v. Con (combined version (in DATA)
XC2 <- fread('combined_fst_1:3.csv')
CX2 <- fread('combined_fst_2:4.csv')

ddat <- rbind(CX2, XC2)
rm(CX2)
rm(XC2)
# Mead FST b/w replicates
ddat2 <- ddat %>%
  group_by(chr, window) %>%
  summarise(meanFst = (mean(meanFst)))

#hist(ddat2$meanFst, breaks = 50)
#mean(ddat2$meanFst)*2 == 0.25
# Using 0.2
#based on histogram and sauron interesct of "curve" and quantile (~75% among all positions below)
ddat2 <- ddat2[which(ddat2$meanFst>0.2),]

#Use positon list and IFELSE statement with thing between window max/min
ddat2$minWindow <- ddat2$window -250
ddat2$maxWindow <- ddat2$window +250

# need to sort positions based on chromosome (for trajectories):
fst_X <- ddat2[which(ddat2$chr=='X'),]
fst_2L <- ddat2[which(ddat2$chr=='2L'),]
fst_2R <- ddat2[which(ddat2$chr=='2R'),]
fst_3L <- ddat2[which(ddat2$chr=='3L'),]
fst_3R <- ddat2[which(ddat2$chr=='3R'),]
fst_4 <- ddat2[which(ddat2$chr=='4'),]

# create and write a positional .csv file
posfst_X <- as.data.frame(cbind(fst_X$window, fst_X$minWindow, fst_X$maxWindow))
posfst_2L <- as.data.frame(cbind(fst_2L$window, fst_2L$minWindow, fst_2L$maxWindow))
posfst_2R <- as.data.frame(cbind(fst_2R$window, fst_2R$minWindow, fst_2R$maxWindow))
posfst_3L <- as.data.frame(cbind(fst_3L$window, fst_3L$minWindow, fst_3L$maxWindow))
posfst_3R <- as.data.frame(cbind(fst_3R$window, fst_3R$minWindow, fst_3R$maxWindow))
posfst_4 <- as.data.frame(cbind(fst_4$window, fst_4$minWindow, fst_4$maxWindow))

