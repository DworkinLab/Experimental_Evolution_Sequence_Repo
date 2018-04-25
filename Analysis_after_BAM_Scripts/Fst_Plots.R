#Plot Treatment Vs. Control Lineages:

## Packages:
require(data.table)
require(tidyverse)

## Main directory holding fst directories for mappers:
setwd("../Data/Fst_combinedComparisons/")

# Control v. Treatment
#layout of list: two files that will be combined (averaged) and the title of the plot
patty <- c('_fst_1:3.csv', '_fst_2:4.csv', '115', '_fst_5:7.csv', '_fst_6:8.csv', '38', '_fst_9:11.csv', '_fst_10:12.csv', '77')
dir <- getwd()
bdbd <- 1

#assuming that patty list has joint replicates side by side (and an even number)
while (bdbd < length(patty)){
  ddat <- NULL
  patt <- patty[bdbd]
  mycomp <- list.files(path = dir, pattern=patt, full.names = T)
  XC2 <- fread(mycomp)
  ddat <- rbind(ddat, XC2)
  patt <- (patty[bdbd+1])
  mycomp <- list.files(path = dir, pattern=patt, full.names = T)
  XC2 <- fread(mycomp)
  ddat <- rbind(ddat, XC2)
  
  rm(XC2)
  tttle <- patty[bdbd+2]
  
  ddat2 <- ddat %>%
    group_by(chr, window) %>%
    summarise(meanFst = (mean(meanFst)))
  rm(ddat)
  
  ddat2$num <- 1:length(ddat2$meanFst)
 
  xcs <- mean(ddat2$meanFst)
  
  g <- nrow(ddat2[which(ddat2$chr=='2L'),])
  h <- nrow(ddat2[which(ddat2$chr=='2R'),])
  i <- nrow(ddat2[which(ddat2$chr=='3L'),])
  j <- nrow(ddat2[which(ddat2$chr=='3R'),])
  k <- nrow(ddat2[which(ddat2$chr=='4'),])
  l <- nrow(ddat2[which(ddat2$chr=='X'),])
  
  # To change the order for X to be first:
  ddat2$number <- c((l+1):(l+g), 
                    (l+g+1):(l+g+h), 
                    (l+g+h+1):(l+g+h+i),
                    (l+g+h+i+1):(l+g+h+i+j),
                    (l+g+h+i+j+1):(l+g+h+i+j+k), 
                    (1:l))
  
  ggxcv <-  ggplot(data = ddat2, aes(x=number, y=meanFst, color=chr))
  ggxcv2 <- ggxcv + 
    geom_point(size=0.5, show.legend = F, alpha = 0.75) + 
    theme(panel.background = element_blank()) +
    scale_y_continuous(limits=c(0, 1), breaks=seq(0, 1, 0.1)) +
    geom_hline(yintercept = xcs) + 
    xlab("Chromosome") +
    ggtitle(tttle) +
    scale_x_discrete(limits=c(l/2, l+(g/2), (l+g+(h/2)), (l+g+h+(i/2)), (l+g+h+i+(j/2)), (l+g+h+i+j+(k/2))), labels = c("X","2L", "2R", '3L', '3R', "4")) +
    scale_colour_manual(values=c("#56B4E9", "#E69F00", 'grey30', 'grey46', 'wheat3', 'lemonchiffon4')) +
    theme(text = element_text(size=20),
          axis.text.x= element_text(size=15), 
          axis.text.y= element_text(size=15))
  
  print(ggxcv2)

  bdbd <- bdbd+3
}
