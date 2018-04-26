#Adjust P -> False Discovery Rate
#TxG_p.adjust:
#Read in data (Copied from other script)

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

CHROMOs <- rbind(chr_X, chr_2L, chr_2R, chr_3L, chr_3R, chr_4)
rm(chr_X, chr_2L, chr_2R, chr_3L, chr_3R, chr_4)

##-----------------------------------------------##
CHROMOs$adjustP <- p.adjust(CHROMOs$mean_p, method = 'fdr')
CHROMOs$logAdjustP <- -log10(CHROMOs$adjustP)
CHROMOs_4 <-  CHROMOs[-which(CHROMOs$adjustP>0.05),]
CHROMOs_4 <-  CHROMOs_4[-which(CHROMOs_4$logAdjustP==Inf),]

#write.csv(CHROMOs_4, file='CHROMO_FDR_Sig.csv', row.names = FALSE)


#Plot:
CHROMOs_5 <- fread("CHROMO_FDR_Sig.csv")

gdat <- ggplot(data = CHROMOs_5, aes(x=number, y=logAdjustP, colour = chr)) + 
geom_point(size = 0.5, show.legend = F) + 
theme(panel.background = element_blank()) + 
xlab("Chromosome") + 
scale_x_discrete(limits=c(a/2, a+(b/2), (a+b+ (c/2)), (a+b+c+(d/2)), (a+b+c+d+(e/2)), (a+b+c+d+e+(f/2))), 
labels = c("X", "2L", '2R', '3L', "3R", "4")) + 
scale_colour_manual(values=c("#56B4E9", "#E69F00", 'grey30', 'grey46', 'wheat3', 'lemonchiffon4')) + 
theme(text = element_text(size=20), axis.text.x= element_text(size=15), axis.text.y= element_text(size=15)) + 
ylab("-log10(p.adjust)")

#print(gdat)
