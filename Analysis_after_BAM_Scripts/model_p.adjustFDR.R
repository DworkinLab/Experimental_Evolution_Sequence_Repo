#Adjust P -> False Discovery Rate
#TxG_p.adjust:
#Read in data (Copied from other script)

### Packages:

require(dplyr)
require(ggplot2)
require(data.table)

#### Read in each Chromosomal Data: 
CHROMOs <- fread('../Data/CHROMO_model_Full.csv')
Zxc_length <- fread('../Data/CHROMO_model_Full_LENGTH.csv')

a <- Zxc_length$x[1]
b <- Zxc_length$x[2]
c <- Zxc_length$x[3]
d <- Zxc_length$x[4]
e <- Zxc_length$x[5]
f <- Zxc_length$x[6]

##-----------------------------------------------##
# 4 Adjustment methods: based on desired pvalue filtering.
#CHROMOs$adjustP <- p.adjust(CHROMOs$max_p, method = 'fdr')
#CHROMOs$adjustP <- p.adjust(CHROMOs$max_p, method = 'bonferroni')

CHROMOs$adjustP <- p.adjust(CHROMOs$mean_p, method = 'fdr')
#CHROMOs$adjustP <- p.adjust(CHROMOs$mean_p, method = 'bonferroni')

CHROMOs$logAdjustP <- -log10(CHROMOs$adjustP)
CHROMOs_4 <-  CHROMOs[-which(CHROMOs$adjustP>0.05),]
CHROMOs_4 <-  CHROMOs_4[-which(CHROMOs_4$logAdjustP==Inf),]

#Plot:
gdat <- ggplot(data = CHROMOs_4, aes(x=number, y=logAdjustP, colour = chr)) + 
geom_point(size = 0.5, show.legend = F) + 
theme(panel.background = element_blank()) + 
xlab("Chromosome") + 
scale_x_discrete(limits=c(a/2, a+(b/2), (a+b+ (c/2)), (a+b+c+(d/2)), (a+b+c+d+(e/2)), (a+b+c+d+e+(f/2))), 
labels = c("X", "2L", '2R', '3L', "3R", "4")) + 
scale_colour_manual(values=c("#56B4E9", "#E69F00", 'grey30', 'grey46', 'wheat3', 'lemonchiffon4')) + 
theme(text = element_text(size=20), axis.text.x= element_text(size=15), axis.text.y= element_text(size=15)) + 
ylab("-log10(p.adjust)")

#print(gdat)


CHROMOs_4  <- fread('../Data/Model_Chromo_FDRp.adjust.csv')
#CHROMOs_4  <- fread('../Data/Model_Chromo_Bonfp.adjust.csv')
Zxc_length <- fread('../Data/CHROMO_model_Full_LENGTH.csv')

a <- Zxc_length$x[1]
b <- Zxc_length$x[2]
c <- Zxc_length$x[3]
d <- Zxc_length$x[4]
e <- Zxc_length$x[5]
f <- Zxc_length$x[6]

gdat <- ggplot(data = CHROMOs_4, aes(x=number, y=logAdjustP, colour = chr)) + 
geom_point(size = 0.5, show.legend = F) + 
theme(panel.background = element_blank()) + 
xlab("Chromosome") + 
scale_x_discrete(limits=c(a/2, a+(b/2), (a+b+ (c/2)), (a+b+c+(d/2)), (a+b+c+d+(e/2)), (a+b+c+d+e+(f/2))), 
                 labels = c("X", "2L", '2R', '3L', "3R", "4")) + 
scale_colour_manual(values=c("#56B4E9", "#E69F00", 'grey30', 'grey46', 'wheat3', 'lemonchiffon4')) + 
theme(text = element_text(size=20), axis.text.x= element_text(size=15), 
      axis.text.y= element_text(size=15)) + 
ylab("-log10(p.adjust)")

print(gdat)

