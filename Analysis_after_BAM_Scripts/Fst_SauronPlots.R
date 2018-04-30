# Sauron Plots and quantiles:
require(tidyverse)

# Read in data:
# 115
# 1:2 = ConR1:ConR2, 1:3 = SelR1:ConR1, 3:4 = SelR1:SelR2, 2:4 = SelR2:ConR2
#CompR1 <- fread('../Data/Fst_combinedComparisons/combined_fst_1:3.csv')
#CompR2 <- fread('../Data/Fst_combinedComparisons/combined_fst_2:4.csv')

#Controls <- fread('../Data/Fst_combinedComparisons/combined_fst_1:2.csv')
#Selections <- fread('../Data/Fst_combinedComparisons/combined_fst_3:4.csv')



#38: 5:6 = ConR1:ConR2, 7:5 = SelR1:ConR1, 7:8 = SelR1:SelR2, 6:8 = SelR2:ConR2

CompR1 <- fread('../Data/Fst_combinedComparisons/combined_fst_5:7.csv')
CompR2 <- fread('../Data/Fst_combinedComparisons/combined_fst_6:8.csv')

Controls <- fread('../Data/Fst_combinedComparisons/combined_fst_5:6.csv')
Selections <- fread('../Data/Fst_combinedComparisons/combined_fst_7:8.csv')



datComp <- merge(CompR1, CompR2,by=c("window","chr"))
datComp$Thing <- "Comparison"
datNoncomp <- merge(Controls, Selections,by=c("window","chr"))
datNoncomp$Thing <- "WithinTreatment"
head(datComp)
head(datNoncomp)
#ggplot(datComp, aes(x=meanFst.x, y=meanFst.y)) + geom_point(size=0.5, alpha=0.5, colour='firebrick3')
#ggplot(datNoncomp, aes(x=meanFst.x, y=meanFst.y)) + geom_point(size=0.5, alpha=0.5, colour='grey30')

ppplt <- ggplot(datComp, aes(x=meanFst.x, y=meanFst.y)) +
  geom_point(size=0.5, alpha=0.5, colour='firebrick3') + 
  geom_point(data=datNoncomp, 
             aes(x=meanFst.x, y=meanFst.y), 
             size=0.5, alpha=0.5, 
             colour='grey30') + 
  ggtitle("Mean Fst Distribution") +
  xlab(expression(atop("ConR1:SelR1[Red]", 'ConR1:ConR2[Grey]'))) + 
  ylab(expression(atop("ConR2:SelR2[Red]", 'SelR1:SelR2[Grey]')))
print(ppplt)


#Can put as one plot if wanted
#source('multiplotFunction.R')
#ppl_115 <- ppplt
#ppl_38 <- ppplt
#multiplot(ppl_115, ppl_38, cols=1)

#Quantiles for interest sake:

with(datComp, quantile(meanFst.x, 
                       probs = c(0, 0.025, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.975, 0.99, 0.999)))

with(datComp, quantile(meanFst.y, 
                       probs = c(0, 0.025, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.975, 0.99, 0.999)))

with(datNoncomp, quantile(meanFst.x, 
                          probs = c(0, 0.025, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.975, 0.99, 0.999)))

with(datNoncomp, quantile(meanFst.y, 
                          probs = c(0, 0.025, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.975, 0.99, 0.999)))

