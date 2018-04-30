# Filter order: does pval + ppool overlap: then is the position in FST window

# Only one option (but can change cut off point)
source('Fst_positons.R')

# Can change b/w FDR and Bonferroni in script, as well as p-value used (either max p-val or mean p-val)
source('model_positions.R')

# Can change b/w FDR and Bonferroni in script, as well as p-value used (either max p-val or mean p-val)
source('selCoef_positions.R')

#Remove all but positions
rm(list=ls()[! ls() %in% c('posfst_X', "posfst_2L", 'posfst_2R', 'posfst_3R', 'posfst_3L','posfst_4','pos_X','pos_2L','pos_2R','pos_3L','pos_3R','pos_4', 'pool_pos_2L', 'pool_pos_2R', 'pool_pos_3L','pool_pos_3R','pool_pos_4','pool_pos_X')])

CHR <- c('X', '2L', '2R', '4', '3L', '3R')
## 2R:
for (chr in CHR){
  posFst <- eval(parse(text = paste("posfst", chr, sep = '_')))
  pos_model <- eval(parse(text = paste("pos", chr, sep = '_')))
  pool_pos <- eval(parse(text = paste("pool_pos", chr, sep = '_')))
  Xcx <- as.data.frame(pos_model)
  Xcx$SIG <- ifelse(Xcx$pos_model  %in% pool_pos, 'SIG', "NOPE")
  Dfd <- Xcx[-which(Xcx$SIG=='NOPE'),]
  Xcs <- as.numeric(Dfd$pos_model)
  int3 <- NULL
  
  for (i in 1:length(posFst$V1)){
    tft <- as.list(posFst$V2[i]:posFst$V3[i])
    int <- intersect(tft, Xcs)
    int2 <- as.data.frame(int)
    int3 <- rbind(int3, int2)
  }
  xc <- unique(int3)
  assign(paste("int", chr, sep = "_"), xc)
  rm(int)
  rm(int2)
  rm(int3)
  rm(Xcs)
  rm(Dfd)
  rm(Xcx)
}

rm(list=ls()[! ls() %in% c('posfst_X', "posfst_2L", 'posfst_2R', 'posfst_3R', 'posfst_3L','posfst_4','pos_X','pos_2L','pos_2R','pos_3L','pos_3R','pos_4', 'pool_pos_2L', 'pool_pos_2R', 'pool_pos_3L','pool_pos_3R','pool_pos_4','pool_pos_X', 'int_2L','int_2R', 'int_3L','int_3R','int_4', 'int_X')]) 

#Depending on p-val method: can change name of output:
#write.csv(int_X, file='X_positions_maxP.csv', row.names = FALSE)
#write.csv(int_2L, file='2L_positions_maxP.csv', row.names = FALSE)
#write.csv(int_2R, file='2R_positions_maxP.csv', row.names = FALSE)
#write.csv(int_3L, file='3L_positions_maxP.csv', row.names = FALSE)
#write.csv(int_3R, file='3R_positions_maxP.csv', row.names = FALSE)
#write.csv(int_4, file='4_positions_maxP.csv', row.names = FALSE)
