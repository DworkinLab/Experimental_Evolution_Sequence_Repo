# Modified function from Taus Pool Seq package 
# Changes the sep='' in the fread to modified .sync file which changed spacing

read.sync_Personal <- function(file, gen, repl, polarization=c("minor", "rising", "reference")) {
polarization <- match.arg(polarization)  
     
  cat("Reading sync file ...\n")
  # load sync-file
  syncDt <- fread(file, sep=" " , header=FALSE, stringsAsFactors=FALSE)
  setDT(syncDt)
  # if either 'gen' or 'repl' are not of propper length then stop
  if((ncol(syncDt)-3) %% length(gen) != 0 || (ncol(syncDt)-3) %% length(repl) != 0)
    stop("Either 'gen' (", length(gen), ") or 'repl' (", length(repl), ") is not a multiple of the number of populations (", ncol(syncDt)-3, ") specified in the sync-file.")
  gc()
  
  cat("Extracting biallelic counts ...\n")
  # extract numeric allele counts for A:T:C:G
  syncCnts <- lapply(syncDt[,-1:-3,with=FALSE], function(col) {
    matrix(as.numeric(stri_split_fixed(sub("(.*):[0-9]+:[0-9]+$", "\\1", col, perl=TRUE), pattern=":", simplify=TRUE)), ncol=4)
  } ) # <- bottleneck
  
  # get rid of data that is no longer needed
  snpCnt <- nrow(syncDt)
  chr <- as.character(syncDt$V1)
  pos <- syncDt$V2
  ref <- syncDt$V3
  popCnt <- ncol(syncDt)-3
  rm(syncDt)
  gc()
  
  # sum up counts across populations (time points and replicates) and add e (random uniform >=0 & <= 0.99) to make each count value unique
  sumCnts <- Reduce("+", syncCnts)
  sumCnts <- sumCnts + runif(nrow(sumCnts)*ncol(sumCnts), min=0, max=0.99)
  # deterime allele ranks for each position
  alleleRank <- rowRanks(sumCnts, ties.method="min")
  # extract 2 most common alleles (considering all populations)
  alleleCnts <- lapply(syncCnts, function(pop) {
    cbind(major=t(pop)[t(alleleRank) == 4], minor=t(pop)[t(alleleRank) == 3])
  } ) # <- bottleneck
  
  cat("Creating result object ...\n")
  # compute chromosome IDs (to later replace character vector by numeric one)
  chrNames <- unique(chr)
  chrID <- 1:length(chrNames)
  names(chrID) <- chrNames
  
  # extract major and minor allele for each position
  syncCntCol <- 1:4
  names(syncCntCol) <- c("A", "T", "C", "G")
  alleles <- data.table(chr=chr, pos=pos, ref=ref,
                        major=names(syncCntCol)[which(t(alleleRank) == 4) - 4*seq(0, nrow(alleleRank)-1)],
                        minor=names(syncCntCol)[which(t(alleleRank) == 3) - 4*seq(0, nrow(alleleRank)-1)],
                        rising=NA_character_)
  
  # if polarization is 'reference' -> check if ref-allele is either major or minor -> warning if not and polarization for minor instead
  if(polarization == "reference" && any(alleles$ref != alleles$minor & alleles$ref != alleles$major)) {
    warning("Cannot polarize for reference allele, because it is not among the two most common alleles for some SNPs. Changing polarization to 'minor'.")
    polarization <- "minor"
  }
  
  # combine generation and replicate info
  popInfo <- data.table(pop=1:popCnt, gen=gen, repl=repl)
  
  # add allele frequency and sequence coverage for each population
  for(r in unique(repl)) {
    for(i in seq(1:nrow(popInfo))[popInfo$repl == r]) {
      seqCov <- rowSums(alleleCnts[[i]])
      # compute allele frequencies according to 'polarization'
      if(polarization == "minor" || polarization == "rising") {
        alleles[,paste0("F", popInfo$gen[i], ".R", r, ".freq"):=alleleCnts[[i]][,"minor"]/seqCov]
      } else {
        alleles[,paste0("F", popInfo$gen[i], ".R", r, ".freq"):=ifelse(minor == ref, alleleCnts[[i]][,"minor"]/seqCov, alleleCnts[[i]][,"major"]/seqCov)]
      }
      alleles[,paste0("F", popInfo$gen[i], ".R", r, ".cov"):=seqCov]
    }
  }
  
  # if required then polarize allele counts for the rising allele
  if(polarization == "rising" && length(unique(popInfo$gen)) > 1) {
    ugens <- unique(popInfo$gen)
    minGen <- min(popInfo$gen)
    
    # if minGen allele frequency column is not available for all replicates then stop execution
    #if(sum(grepl(paste0("F", minGen, "\\.[R0-9]+\\.freq"), colnames(alleles))) != length(unique(repl)))
    #  stop("Cannot polarize for rising allele because not all replicates provide allele frequency estimates at generation F", minGen)
    
    # calculate mean allele frequency change per SNP and replicate
    meanAF <- foreach(r=unique(repl), .combine=cbind, .final=function(x) { if(is.matrix(x)) return(rowMeans(x)) else return(x) }) %do% {
      allCols <- grep(paste0("F[0-9]+\\.R", r, "\\.freq"), colnames(alleles), value=TRUE)
      allCols <- allCols[order(as.numeric(sub("F([0-9]+)\\..*", "\\1", allCols)))]
      rowMeans(alleles[,allCols[-1],with=FALSE]-alleles[[allCols[1]]])
    }
    
    # polarize allele frequencies
    needsPolarization <- meanAF < 0
    for(pop in grep("F[0-9]+\\.R[0-9]+\\.freq", colnames(alleles), value=TRUE)) {
      alleles[,eval(pop):=ifelse(needsPolarization, 1-alleles[[pop]], alleles[[pop]])]
    }
    
    # set column with rising allele
    alleles[,rising:=ifelse(needsPolarization, alleles$major, alleles$minor)]
  }
  
  # return sync-object for loaded sync-file
  return(new(Class="sync",
             gen=as.numeric(sub("F([0-9]+)\\.R[0-9]+.*", "\\1", colnames(alleles)[-1:-6])),
             repl=as.numeric(sub(".*\\.R([0-9]+)\\..*", "\\1", colnames(alleles)[-1:-6])),
             isAF=grepl(".*\\.freq$", colnames(alleles)[-1:-6]),
             polarization=polarization,
             alleles=alleles))
}
