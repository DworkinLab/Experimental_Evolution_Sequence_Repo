setwd('/Users/paulknoops/Bioinformatics/Analysis_EpisodicSequenceData_2018/poolSeq')
require(tidyverse)
require(data.table)
#Read in data:
#How large data created:
## source(poolseq_CombineMappers.R)

# How subset is created (after adjustments)
# Zxc <- fread('SelCoef_Full.csv')

# adjust p values
# Zxc$adjustP <- p.adjust(Zxc$pval_max, method = 'fdr')
# Zxc$adjustP <- p.adjust(Zxc$pval_max, method = 'bonferroni')

# Lable for significance 
#Zxc$sig <- ifelse(Zxc$adjustP<0.05, "<0.05", ">0.05")
#ggplot(Zxc, aes(number, meanSelCoef, colour=sig)) + geom_point() 
# Only significant positions:
#Zxc_sig <- Zxc[which(Zxc$sig=='<0.05'),]

#Write the above and read in as smaller data set:
#write.csv(Zxc_sig, file='poolseq_FDR_Sigpos.csv', row.names = F)
#write.csv(Zxc_sig, file='poolseq_bonf_Sigpos.csv', row.names = F)

rm(list=ls())

## Read in made data: Three plots created (choose the nicest, 3rd one in the .md)
Zxc_sig <- fread('poolseq_FDR_Sigpos.csv')
#Zxc_sig <- fread('poolseq_bonf_Sigpos.csv')

Zxc_length <- fread('SelCoef_Full_LENGTHS.csv')

aa <- Zxc_length$x[1]
bb <- Zxc_length$x[2]
cc <- Zxc_length$x[3]
dd <- Zxc_length$x[4]
ee <- Zxc_length$x[5]
ff <- Zxc_length$x[6]

mmmax <- max(Zxc_sig$meanSelCoef)


pplot_4 <- ggplot(data = Zxc_sig, aes(x=number, y=meanSelCoef, colour=Treatment, shape=Treatment)) + 
  geom_point(size = 0.75, alpha=0.5) + 
  theme(panel.background = element_blank()) + 
  xlab("Chromosome") + 
  scale_x_discrete(limits=c(aa/2, aa+(bb/2), (aa+bb+(cc/2)), (aa+bb+cc+(dd/2)), (aa+bb+cc+dd+(ee/2)), (aa+bb+cc+dd+ee+(ff/2))), labels = c('X', "2L", "2R", '3L', '3R', '4')) + 
  scale_colour_manual(values=c("grey50", "firebrick3")) + 
  theme(text = element_text(size=20), 
        axis.text.x= element_text(size=15), 
        axis.text.y= element_text(size=15)) + 
  ylab("SelCoef") + geom_vline(xintercept = aa) + 
  geom_vline(xintercept = (aa+bb)) + 
  geom_vline(xintercept = (aa+bb+cc)) + 
  geom_vline(xintercept = (aa+bb+cc+dd)) + 
  geom_vline(xintercept = (aa+bb+cc+dd+ee)) + 
  ylim(0,mmmax)
print(pplot_4)

# 
Zxc_count <- Zxc_sig %>%
  group_by(chr, pos) %>%
  mutate(count = n())
Zxc_count2 <- Zxc_count[which(Zxc_count$count==1),]
Zxc_count2 <- Zxc_count2[which(Zxc_count2$Treatment=='Sel'),]

mmaaxx <- max(Zxc_count2$meanSelCoef)

pplot_5 <- ggplot(data = Zxc_count2, aes(x=number, y=meanSelCoef, colour=chr)) + 
  geom_point(size = 0.75, alpha=0.75, show.legend = F) + 
  theme(panel.background = element_blank()) + 
  xlab("Chromosome") + 
  scale_x_discrete(limits=c(aa/2, aa+(bb/2), (aa+bb+(cc/2)), (aa+bb+cc+(dd/2)), (aa+bb+cc+dd+(ee/2)), (aa+bb+cc+dd+ee+(ff/2))), labels = c('X', "2L", "2R", '3L', '3R', '4')) + 
  scale_colour_manual(values=c("#56B4E9", "#E69F00", 'grey30', 'grey46', 'wheat3', 'lemonchiffon4')) + 
  theme(text = element_text(size=20), 
        axis.text.x= element_text(size=15), 
        axis.text.y= element_text(size=15)) + 
  ylab("SelCoef") + 
  ylim(0,mmaaxx)
print(pplot_5)


Zxc_Sel <- Zxc_sig[which(Zxc_sig$Treatment=='Sel'),]
Zxc_Con <- Zxc_sig[which(Zxc_sig$Treatment=='Con'),]

pplot_6 <- ggplot(data = Zxc_Sel, aes(x=number, y=meanSelCoef, colour=chr)) + 
  geom_point(size = 0.66, alpha=1, show.legend = F) + 
  geom_point(data = Zxc_Con, size=0.33, alpha=0.5,  aes(x=number, y=meanSelCoef, colour=Treatment), show.legend = F) +
  theme(panel.background = element_blank()) + 
  xlab("Chromosome") + 
  scale_x_discrete(limits=c(aa/2, aa+(bb/2), (aa+bb+(cc/2)), (aa+bb+cc+(dd/2)), (aa+bb+cc+dd+(ee/2)), (aa+bb+cc+dd+ee+(ff/2))), labels = c('X', "2L", "2R", '3L', '3R', '4')) + 
  scale_colour_manual(values=c("#56B4E9", "#E69F00", 'grey30', 'grey46', 'wheat3','firebrick3', 'lemonchiffon4')) + 
  theme(text = element_text(size=20), 
        axis.text.x= element_text(size=15), 
        axis.text.y= element_text(size=15)) + 
  ylab("SelCoef") + 
  ylim(0,mmaaxx)
print(pplot_6)
