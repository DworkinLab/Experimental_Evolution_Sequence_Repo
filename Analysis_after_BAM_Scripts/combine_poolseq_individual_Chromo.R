# set directory holding all .csv files to combine
setwd('/home/paul/episodicData/novoalign/novo_mpileup/splitsync_dir/novo_episodic_2L_Sel_Split')

#list Csvs
mycsvs <- list.files(pattern='.csv')
X <- NULL
for (file in mycsvs){
  X2 <- read.csv(file, h=T)
  X <- rbind(X, X2)
}
# change based on chromo!
X$chr <- '2L'

#write the CSV file !!! EASY PEASY
write.csv(X, file='/home/paul/episodicData/novoalign/novo_mpileup/novo_episodic_2L_Sel.csv', row.names = FALSE)
