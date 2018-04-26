
#Combine files from poolseq output!

args <- commandArgs(trailingOnly = TRUE)
mydirs <- list.dirs(path=args[1], recursive = FALSE) 
XX <- args[1]
for (dir in mydirs){
	setwd(dir)
	J4 <- gsub('(.*)_\\w+', '\\1', dir)
	mycsvs <- list.files(pattern='.csv')
	X <- NULL
	for (file in mycsvs){
		X2 <- read.csv(file, h=T)
		X <- rbind(X, X2)
   	}
	write.csv(X, file=paste(J4,'.csv',sep=""), row.names = FALSE)
	rm(X)
	rm(J4)
}
