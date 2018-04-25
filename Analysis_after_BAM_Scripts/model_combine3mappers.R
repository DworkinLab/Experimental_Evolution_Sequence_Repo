#To run as Rscript and have output directory argument:
# need to change the files directories:
args <- commandArgs(trailingOnly = TRUE)

# combine the coeffs csv files for three mappers
require(dplyr)

# Will call these 6 Chromosoms (but this depends on the files all ending with [name2]_chromo.csv (ex. 2L_chromo.csv)
name2 <- c('2L', '2R', '3L', '3R', '4','X')

#Loop through all chromosomes: each read.csv has the directory holding all the chromosome subsets from the output of the model run just above. 

# Will combine each chromosome into one file with three mappers (defined in loops) and written to a new large .csv file

for (i in name2){

  print(i)
  
  file=paste(i,'_chromo.csv', sep = "")
  
  bwa_coeffs <- read.csv(paste('/home/paul/episodicData/R_dir/bwa_coeffs/episodic_data_', file, sep=""), header = TRUE)
  bwa_coeffs$mapper <- "bwa"
  
  bowtie_coeffs <- read.csv(paste('/home/paul/episodicData/bowtie/R_bowtie/bowtie_coeffs/episodic_data_bowtie_', file, sep=""), header = TRUE)
  bowtie_coeffs$mapper <- "bowtie"
  
  novo_coeffs <- read.csv(paste('/home/paul/episodicData/novoalign/novo_coeffs/novo_episodic_', file, sep=""), header = TRUE)
  novo_coeffs$mapper <- "novoalign"
  
  X <- rbind(bwa_coeffs, bowtie_coeffs, novo_coeffs)
  rm(coeffs_bowtie)
  rm(coeffs_bwa)
  rm(novo_coeffs)
  
  outputDir <- args[1]
 
  write.csv(X , file=paste(outputDir,"/",'Chromosome_', i, "_Full.csv", sep=""), row.names = FALSE)
}
