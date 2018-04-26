################
### Running the PoolSeq Package to run on .sync files for estimates of selection coefficeints per position
### Requires: R (>= 3.3.1), data.table (>= 1.9.4), foreach (>= 1.4.2), stringi (>= 0.4-1), matrixStats (>= 0.14.2)

  args <- commandArgs(trailingOnly = TRUE)
  
### Required Packages:

  #install.packages("/home/paul/poolSeq_0.3.2.tar.gz", repos=NULL, type="source")
  #install.packages("/home/paul/matrixStats_0.53.0.tar.gz", repos=NULL, type="source")

### Not available: so source seperate:
  #require(poolSeq)
  
### These are part of the dependencies for poolSeq
  
  require(methods)
  require(data.table)
  require(foreach)
  require(stringi)
  require(matrixStats)
  
### Source the scripts (Copied) for Pool-Seq (only one fails and is not needed)
  source('/home/paul/episodicData/novoalign/novo_Rscripts/Taus_Scripts/testTaus/loadaf.R')  
  #estne.R Fails.
  #source('/home/paul/episodicData/novoalign/novo_Rscripts/Taus_Scripts/testTaus/estne.R')
  source('/home/paul/episodicData/novoalign/novo_Rscripts/Taus_Scripts/testTaus/estsh.R')
  source('/home/paul/episodicData/novoalign/novo_Rscripts/Taus_Scripts/testTaus/idsel.R')
  source('/home/paul/episodicData/novoalign/novo_Rscripts/Taus_Scripts/testTaus/simaf.R')

### Possibly need custom function to read in manipulated .sync files:
	### Needed for manipulated .sync files (one basic change labeled at top of changed script:
  
  source("/home/paul/episodicData/novoalign/novo_Rscripts/Taus_ReadSync.R")

### Read in the data file for args[1]

  setwd(args[3])
  
  mySync <- read.sync_Personal(file=args[1], gen=c(115, 115, 38, 38, 77, 77, 0, 0), repl=c(1,2,1,2,1,2,1,2), polarization = "rising")

### Make data.frame of just alleles information to sort out relevent positions:
# Turn alleles to data frame:
ff <- as.data.frame(mySync@alleles)

# Keep only positions:
pst <- as.numeric(ff$pos)
pst2 <- sort(pst)

# Generations:
ccc <- c(0,38,77,115)

rm(pst)
rm(ff)


### Create empty matrix to read into for estiamting S:
pbj <- matrix(NA,length(pst2), 3)

 for (i in 1:length(pst2)) {
    b_b <- pst2[i]
    TrajTEST <- af.traj(mySync, args[2], repl=c(1,2), pos=b_b)
    BfsfTEST <- estimateSH(TrajTEST, Ne=150, t=ccc, h=0.5, haploid = FALSE, simulate.p.value=TRUE)
    pbj[i,] <- c(BfsfTEST$s, BfsfTEST$p.value, b_b)
    rm(TrajTEST)
    rm(BfsfTEST)
    rm(b_b)
  }

x2 <- args[1]
x3 <- gsub("\\..*","", x2)
write.csv(pbj, file=paste(x3, ".csv", sep=""), row.names=FALSE)

rm(pbj)
rm(mySync)
rm(ccc)
rm(pst2)
