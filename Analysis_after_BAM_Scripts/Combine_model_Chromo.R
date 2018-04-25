# Combine .coeffs.csv for two mappers
require(dplyr)

args <- commandArgs(trailingOnly = TRUE)

# Change to directory holding all directories:

mydirs <- list.dirs(path = args[1], recursive = FALSE)

for (dir in mydirs){

  setwd(dir)
  
  print("Read coeffs.csv files")

  mycsvs <- list.files(pattern=".coeffs.csv")

  Novoalign_Chromosome <- NULL
  
  for (file in mycsvs){
    print(file)
    coeffs1 <- read.csv(file, h=T)
    Novoalign_Chromosome  <- rbind(Novoalign_Chromosome , coeffs1)
    rm(coeffs1)
    J2 <- gsub("*_","", file)
    J3 <- gsub("\\..*","",J2)
}

x3 <- gsub("\\..*","",file)
J3 <- gsub('(.*)_\\w+', '\\1', x3)

X <- args[2]

write.csv(Novoalign_Chromosome , file=paste(X,"/",J3,"_chromo.csv", sep=""), row.names = FALSE)
rm(x3)
rm(J3)
rm(Novoalign_Chromosome)

}
