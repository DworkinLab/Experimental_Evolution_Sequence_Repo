# Script to read in one .fst file and split this into many .csv files based on comparisons
# .fst files generated from Popoolation2 fst-sliding.pl script
# Need to change the working directory, the input, and the number of comparisons present (i.e 6:ncol for .fst file)

### Packages Required (tidyverse, but more specifically tidyr)
  require(data.table)
  require(tidyverse)

### Set working directory to location of .fst file
  setwd("~/Bioinformatics/episodic_practice/FST")
  
### Read in the .fst file into R (requires data.table)
  Fst_novo <- fread('novo_episodic_main.fst')

### Make into long format
  ccol <- ncol(Fst_novo)
  XCC  <- gather(Fst_novo, Comparison, Fst_measure, 6:ccol, factor_key=TRUE)

### Remove intermediate:
  rm(Fst_novo)

### Seperate the Fst (ex. 1:2=na) into a comparison column and a fst column
  novo_Fst <- XCC %>%
    separate(Fst_measure, "=", into = c('Comp', 'Fst'))

### Remove intermediate:
  rm(XCC)

### Remove unnecessary column (column 6 has no value)
  novo_Fst <- novo_Fst[,c(1,2,3,4,5,7,8)]

### Rename columns:
  colnames(novo_Fst) <- c('chr', 'window', "num", 'frac', 'meanCov','Comp', 'Fst')

### Create list of all unique comparisons:
  X_compLIST <- unique(novo_Fst$Comp)

### For loop that will create a .csv file for every comparison:
  for (vxcx in X_compLIST){

    CXV_Comp <- novo_Fst[which(novo_Fst$Comp==vxcx),]
    
    write.csv(CXV_Comp, file=paste("Novo_fst_", vxcx, '.csv', sep = ""), row.names = FALSE)
    }
