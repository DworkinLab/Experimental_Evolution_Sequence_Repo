# Script to combine mappers into one average Fst value for each comparison of interest.
# Assumptions: 
    # one main directory (working dir) holds independent directories for each mapper (i.e novo_fst, bwa_fst, and bowtie_fst)
    # Comparison files split previously
    # have the known comprisons of interest (change in patty <- )
   
## Packages:
  require(data.table)
  require(tidyverse)

## Main directory holding fst directories for mappers:
  setwd("/Users/paulknoops/Bioinformatics/episodic_practice/FST")
  
## List directories in working directories:
  workdir <- getwd()
  mydirs <- list.dirs(path=workdir, recursive = FALSE)
  
## Make list of compasisons to combine (no use doing them all?)
  patty <- c('_fst_1:3.csv', '_fst_2:4.csv', '_fst_5:7.csv', '_fst_6:8.csv', '_fst_9:11.csv', '_fst_10:12.csv')

for (patt in patty){
  mtdf <- NULL
  for (dir in mydirs){
    mycomp <- list.files(path = dir, pattern=patt, full.names = T)
    for (file in mycomp){
      Xc <- fread(file)
    }
    Xc <- Xc[-which(Xc$Fst=='na'),]
    Xc$Fst <- as.numeric(Xc$Fst)
    mtdf <- rbind(mtdf, Xc)
  }

  comp_final <- mtdf %>%
    group_by(chr, window, Comp) %>%
    mutate(count = n())
  
  compcomb <- comp_final[which(comp_final$count==3),]
  
  comp_final_2 <- mtdf %>%
    group_by(chr, window, Comp) %>%
    summarise(meanFst = (mean(Fst)))
  
  write.csv(comp_final_2, file=paste("combined", patt, sep = ""), row.names = FALSE)
  }
