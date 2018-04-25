#To run as Rscript and have output directory argument:
args <- commandArgs(trailingOnly = TRUE)
#Can input manually:

setwd(args[1])
#setwd('/home/paul/coeffs_fullChromo')

#Can change the effect of interest by changing hashes at X_Effect (#EFFECT OF INTEREST#)

require(dplyr)
mycsvs <- list.files(pattern='_Full.csv')


for (file in mycsvs){
  #file = 'Chromosome_4_Full.csv'
  print(file)
  name <- gsub("^.*?_","",file)
  name2 <- gsub("_.*","",name)
  rm(name)
  X <- read.csv(file, h=T)

  
#EFFECT OF INTEREST#
#X_Effect <- X[which(X$Effects=="TreatmentSel"),]
X_Effect <- X[which(X$Effects=="TreatmentSel:Generation"),]
#X_Effect <- X[which(X$Effects=="Intercept"),]
#X_Effect <- X[which(X$Effects=="Generation"),]

title <- as.character(X_Effect$Effects[1])
title2 <- ifelse(title=='TreatmentSel', 'Treat', ifelse(title=='Generation', 'Gen', ifelse(title=='TreatmentSel:Generation', 'TxG', '
                                                                                           Int')))

rm(title)
rm(X)

# Only keep those with a position mapped by all three (create table of posiitons and all those >=3 kept)
tt <- table(X_Effect$position)
Effects_Final <- subset(X_Effect, position %in% names(tt[tt >= 3]))

rm(tt)
rm(X_Effect)

print('writing CSV')
write.csv(Effects_Final, file=paste0(name2, '_Chromosome_', title2, '.csv'), row.names = FALSE)


rm(Effects_Final)

print('Done and everything gone')
}
