#Combine chromosome into one data file of candidates:

Pos_X <- fread('X_positions_maxP.csv')
#Pos_X <- fread('X_positions_meanP.csv')
Pos_X$chr <- 'X'

Pos_2L <- fread('2L_positions_maxP.csv')
#Pos_2L <- fread('2L_positions_meanP.csv')
Pos_2L$chr <- '2L'

Pos_2R <- fread('2R_positions_maxP.csv')
#Pos_2R <- fread('2R_positions_meanP.csv')
Pos_2R$chr <- '2R'

Pos_3L <- fread('3L_positions_maxP.csv')
#Pos_3L <- fread('3L_positions_meanP.csv')
Pos_3L$chr <- '3L'

Pos_3R <- fread('3R_positions_maxP.csv')
#Pos_3R <- fread('3R_positions_meanP.csv')
Pos_3R$chr <- '3R'

Pos_4 <- fread('4_positions_maxP.csv')
#Pos_4 <- fread('4_positions_meanP.csv')
Pos_4$chr <- '4'

#  CHR-Max-Mean
#   2L-71-80
#   2R-51-52
#   3L-19-19
#   3R-23-24
#   4-0-0
#   X-73-86

Ppos <- rbind(Pos_2L, Pos_2R, Pos_3L, Pos_3R, Pos_4, Pos_X)
Ppos$pos <- Ppos$int
Ppos_2 <- subset( Ppos, select = -int )

write.csv(Ppos_2, file = "candidatePos.csv", row.names = F)

Xc <- fread('candidatePos.csv')
write.table(Xc, file='candidatePos.txt', sep ="\t", col.names = F, row.names = F)
cx <- fread('candidatePos.txt')
