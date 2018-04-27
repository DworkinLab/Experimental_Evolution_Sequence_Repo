# Combine mappers and keep mean selection coefficient and 2 pvalues: less significant p-values (max P) and mean p value

# Combine pool-seq outputs:
setwd('/Users/paulknoops/Bioinformatics/Analysis_EpisodicSequenceData_2018/poolSeq')
require(tidyverse)
require(data.table)
#Read in data:
###-----      X       -----###
column.names <- c('selcoef', 'pval', 'pos', 'chr')

# Selection
novo_sel_X <- fread('novo_episodic_X_Sel.csv')
bwa_sel_X <- fread('bwa_episodic_X_Sel.csv')
colnames(novo_sel_X) <- column.names
colnames(bwa_sel_X) <- column.names
novo_sel_X <- na.omit(novo_sel_X)
bwa_sel_X <- na.omit(bwa_sel_X)
novo_sel_X$Treatment <- 'Sel'
bwa_sel_X$Treatment <- 'Sel'
SelX <- rbind(novo_sel_X, bwa_sel_X)

Sel_X <- SelX %>%
  group_by(chr, Treatment, pos) %>%
  mutate(count = n())
Sel_X <- Sel_X[which(Sel_X$count==2),]

Sel_X <- Sel_X %>%
  group_by(chr, Treatment, pos) %>%
  summarise(meanSelCoef = (mean(selcoef)),
            pval_max=max(pval), pval_mean=mean(pval))

# Controls:
bwa_con_X <- fread('bwa_episodic_X_Con.csv')
novo_con_X <- fread('novo_episodic_X_Con.csv')
colnames(novo_con_X) <- column.names
colnames(bwa_con_X) <- column.names
bwa_con_X <- na.omit(bwa_con_X)
novo_con_X <- na.omit(novo_con_X)
novo_con_X$Treatment <- 'Con'
bwa_con_X$Treatment <- 'Con'
ConX <- rbind(novo_con_X, bwa_con_X)
Con_X <- ConX %>%
  group_by(chr, Treatment, pos) %>%
  mutate(count = n())
Con_X <- Con_X[which(Con_X$count==2),]

Con_X <- Con_X %>%
  group_by(chr, Treatment, pos) %>%
  summarise(meanSelCoef = (mean(selcoef)),
            pval_max=max(pval), pval_mean=mean(pval))


datX <- merge(Con_X, Sel_X, by = 'pos', all =T)

aa <- length(datX$pos)
datX$number <- 1:aa
#dat2L$number <- (1+aa):(aa+bb)
#dat2R$number <- (1+aa+bb):(aa+bb+cc)
#dat3L$number <- (1+aa+bb+cc):(aa+bb+cc+dd)
#dat3R$number <-  (aa+bb+cc+dd+1):(aa+bb+cc+dd+ee)
#dat4$number <- (aa+bb+cc+dd+ee+1):(aa+bb+cc+dd+ee+ff)

Con_X <- datX[c(1,2,3,4,5,6,12)]
Sel_X <- datX[c(1,7,8,9,10,11,12)]
column.names_2 <- c('pos','chr', 'Treatment', 'meanSelCoef', 'pval_max','pval_mean', 'number')
colnames(Sel_X) <- column.names_2
colnames(Con_X) <- column.names_2

Xcx_X <- rbind(Con_X, Sel_X)
rm(novo_con_X, novo_sel_X, bwa_con_X, bwa_sel_X)

#rm(list=ls()[! ls() %in% c('Xcx_2L', 'Xcx_2R', 'Xcx_3L','Xcx_3R','Xcx_4','Xcx_X', 'aa', 'bb', 'cc','dd','ee','ff')])


###-----      2L       -----###
column.names <- c('selcoef', 'pval', 'pos', 'chr')

# Selection
novo_sel_2L <- fread('novo_episodic_2L_Sel.csv')
bwa_sel_2L <- fread('bwa_episodic_2L_Sel.csv')
colnames(novo_sel_2L) <- column.names
colnames(bwa_sel_2L) <- column.names
novo_sel_2L <- na.omit(novo_sel_2L)
bwa_sel_2L <- na.omit(bwa_sel_2L)
novo_sel_2L$Treatment <- 'Sel'
bwa_sel_2L$Treatment <- 'Sel'
Sel2L <- rbind(novo_sel_2L, bwa_sel_2L)
Sel_2L <- Sel2L %>%
  group_by(chr, Treatment, pos) %>%
  mutate(count = n())
Sel_2L <- Sel_2L[which(Sel_2L$count==2),]
Sel_2L <- Sel_2L %>%
  group_by(chr, Treatment, pos) %>%
  summarise(meanSelCoef = (mean(selcoef)),
            pval_max=max(pval), pval_mean=mean(pval))

# Controls:
bwa_con_2L <- fread('bwa_episodic_2L_Con.csv')
novo_con_2L <- fread('novo_episodic_2L_Con.csv')
colnames(novo_con_2L) <- column.names
colnames(bwa_con_2L) <- column.names
bwa_con_2L <- na.omit(bwa_con_2L)
novo_con_2L <- na.omit(novo_con_2L)
novo_con_2L$Treatment <- 'Con'
bwa_con_2L$Treatment <- 'Con'
Con2L <- rbind(novo_con_2L, bwa_con_2L)
Con_2L <- Con2L %>%
  group_by(chr, Treatment, pos) %>%
  mutate(count = n())
Con_2L <- Con_2L[which(Con_2L$count==2),]
Con_2L <- Con_2L %>%
  group_by(chr, Treatment, pos) %>%
  summarise(meanSelCoef = (mean(selcoef)),
            pval_max=max(pval), pval_mean=mean(pval))

dat2L <- merge(Con_2L, Sel_2L, by = 'pos', all =T)

bb <- length(dat2L$pos)
#datX$number <- 1:aa
dat2L$number <- (1+aa):(aa+bb)
#dat2R$number <- (1+aa+bb):(aa+bb+cc)
#dat3L$number <- (1+aa+bb+cc):(aa+bb+cc+dd)
#dat3R$number <-  (aa+bb+cc+dd+1):(aa+bb+cc+dd+ee)
#dat4$number <- (aa+bb+cc+dd+ee+1):(aa+bb+cc+dd+ee+ff)
Con_2L <- dat2L[c(1,2,3,4,5,6,12)]
Sel_2L <- dat2L[c(1,7,8,9,10,11,12)]
column.names_2 <- c('pos','chr', 'Treatment', 'meanSelCoef', 'pval_max','pval_mean', 'number')
colnames(Sel_2L) <- column.names_2
colnames(Con_2L) <- column.names_2

Xcx_2L <- rbind(Con_2L, Sel_2L)
rm(novo_con_2L, novo_sel_2L, bwa_con_2L, bwa_sel_2L)

#rm(list=ls()[! ls() %in% c('Xcx_2L', 'Xcx_2R', 'Xcx_3L','Xcx_3R','Xcx_4','Xcx_X', 'aa', 'bb', 'cc','dd','ee','ff')])

###-----      2R       -----###
column.names <- c('selcoef', 'pval', 'pos', 'chr')

# Selection
novo_sel_2R <- fread('novo_episodic_2R_Sel.csv')
bwa_sel_2R <- fread('bwa_episodic_2R_Sel.csv')
colnames(novo_sel_2R) <- column.names
colnames(bwa_sel_2R) <- column.names
novo_sel_2R <- na.omit(novo_sel_2R)
bwa_sel_2R <- na.omit(bwa_sel_2R)
novo_sel_2R$Treatment <- 'Sel'
bwa_sel_2R$Treatment <- 'Sel'
Sel2R <- rbind(novo_sel_2R, bwa_sel_2R)
Sel_2R <- Sel2R %>%
  group_by(chr, Treatment, pos) %>%
  mutate(count = n())
Sel_2R <- Sel_2R[which(Sel_2R$count==2),]
Sel_2R <- Sel_2R %>%
  group_by(chr, Treatment, pos) %>%
  summarise(meanSelCoef = (mean(selcoef)),
            pval_max=max(pval), pval_mean=mean(pval))

# Controls:
bwa_con_2R <- fread('bwa_episodic_2R_Con.csv')
novo_con_2R <- fread('novo_episodic_2R_Con.csv')
colnames(novo_con_2R) <- column.names
colnames(bwa_con_2R) <- column.names
bwa_con_2R <- na.omit(bwa_con_2R)
novo_con_2R <- na.omit(novo_con_2R)
novo_con_2R$Treatment <- 'Con'
bwa_con_2R$Treatment <- 'Con'
Con2R <- rbind(novo_con_2R, bwa_con_2R)
Con_2R <- Con2R %>%
  group_by(chr, Treatment, pos) %>%
  mutate(count = n())
Con_2R <- Con_2R[which(Con_2R$count==2),]
Con_2R <- Con_2R %>%
  group_by(chr, Treatment, pos) %>%
  summarise(meanSelCoef = (mean(selcoef)),
            pval_max=max(pval), pval_mean=mean(pval))

dat2R <- merge(Con_2R, Sel_2R, by = 'pos', all =T)
cc <- length(dat2R$pos)
#datX$number <- 1:aa
#dat2L$number <- (1+aa):(aa+bb)
dat2R$number <- (1+aa+bb):(aa+bb+cc)
#dat3L$number <- (1+aa+bb+cc):(aa+bb+cc+dd)
#dat3R$number <-  (aa+bb+cc+dd+1):(aa+bb+cc+dd+ee)
#dat4$number <- (aa+bb+cc+dd+ee+1):(aa+bb+cc+dd+ee+ff)
Con_2R <- dat2R[c(1,2,3,4,5,6,12)]
Sel_2R <- dat2R[c(1,7,8,9,10,11,12)]
column.names_2 <- c('pos','chr', 'Treatment', 'meanSelCoef', 'pval_max','pval_mean', 'number')
colnames(Sel_2R) <- column.names_2
colnames(Con_2R) <- column.names_2

Xcx_2R <- rbind(Con_2R, Sel_2R)

#rm(list=ls()[! ls() %in% c('Xcx_2L', 'Xcx_2R', 'Xcx_3L','Xcx_3R','Xcx_4','Xcx_X', 'aa', 'bb', 'cc','dd','ee','ff')])


###-----      3L       -----###
column.names <- c('selcoef', 'pval', 'pos', 'chr')

# Selection
novo_sel_3L <- fread('novo_episodic_3L_Sel.csv')
bwa_sel_3L <- fread('bwa_episodic_3L_Sel.csv')
colnames(novo_sel_3L) <- column.names
colnames(bwa_sel_3L) <- column.names
novo_sel_3L <- na.omit(novo_sel_3L)
bwa_sel_3L <- na.omit(bwa_sel_3L)
novo_sel_3L$Treatment <- 'Sel'
bwa_sel_3L$Treatment <- 'Sel'
Sel3L <- rbind(novo_sel_3L, bwa_sel_3L)
Sel_3L <- Sel3L %>%
  group_by(chr, Treatment, pos) %>%
  mutate(count = n())
Sel_3L <- Sel_3L[which(Sel_3L$count==2),]
Sel_3L <- Sel_3L %>%
  group_by(chr, Treatment, pos) %>%
  summarise(meanSelCoef = (mean(selcoef)),
            pval_max=max(pval), pval_mean=mean(pval))

# Controls:
bwa_con_3L <- fread('bwa_episodic_3L_Con.csv')
novo_con_3L <- fread('novo_episodic_3L_Con.csv')
colnames(novo_con_3L) <- column.names
colnames(bwa_con_3L) <- column.names
bwa_con_3L <- na.omit(bwa_con_3L)
novo_con_3L <- na.omit(novo_con_3L)
novo_con_3L$Treatment <- 'Con'
bwa_con_3L$Treatment <- 'Con'
Con3L <- rbind(novo_con_3L, bwa_con_3L)
Con_3L <- Con3L %>%
  group_by(chr, Treatment, pos) %>%
  mutate(count = n())
Con_3L <- Con_3L[which(Con_3L$count==2),]
Con_3L <- Con_3L %>%
  group_by(chr, Treatment, pos) %>%
  summarise(meanSelCoef = (mean(selcoef)),
            pval_max=max(pval), pval_mean=mean(pval))

dat3L <- merge(Con_3L, Sel_3L, by = 'pos', all =T)
dd <- length(dat3L$pos)
#datX$number <- 1:aa
#dat2L$number <- (1+aa):(aa+bb)
#dat2R$number <- (1+aa+bb):(aa+bb+cc)
dat3L$number <- (1+aa+bb+cc):(aa+bb+cc+dd)
#dat3R$number <-  (aa+bb+cc+dd+1):(aa+bb+cc+dd+ee)
#dat4$number <- (aa+bb+cc+dd+ee+1):(aa+bb+cc+dd+ee+ff)
Con_3L <- dat3L[c(1,2,3,4,5,6,12)]
Sel_3L <- dat3L[c(1,7,8,9,10,11,12)]
column.names_2 <- c('pos','chr', 'Treatment', 'meanSelCoef', 'pval_max','pval_mean', 'number')
colnames(Sel_3L) <- column.names_2
colnames(Con_3L) <- column.names_2

Xcx_3L <- rbind(Con_3L, Sel_3L)
rm(novo_con_3L, novo_sel_3L, bwa_con_3L, bwa_sel_3L)

#rm(list=ls()[! ls() %in% c('Xcx_2L', 'Xcx_2R', 'Xcx_3L','Xcx_3R','Xcx_4','Xcx_X', 'aa', 'bb', 'cc','dd','ee','ff')])

###-----      3R       -----###
column.names <- c('selcoef', 'pval', 'pos', 'chr')

# Selection
novo_sel_3R <- fread('novo_episodic_3R_Sel.csv')
bwa_sel_3R <- fread('bwa_episodic_3R_Sel.csv')
colnames(novo_sel_3R) <- column.names
colnames(bwa_sel_3R) <- column.names
novo_sel_3R <- na.omit(novo_sel_3R)
bwa_sel_3R <- na.omit(bwa_sel_3R)
novo_sel_3R$Treatment <- 'Sel'
bwa_sel_3R$Treatment <- 'Sel'
Sel3R <- rbind(novo_sel_3R, bwa_sel_3R)
Sel_3R <- Sel3R %>%
  group_by(chr, Treatment, pos) %>%
  mutate(count = n())
Sel_3R <- Sel_3R[which(Sel_3R$count==2),]

Sel_3R <- Sel_3R %>%
  group_by(chr, Treatment, pos) %>%
  summarise(meanSelCoef = (mean(selcoef)),
            pval_max=max(pval), pval_mean=mean(pval))

# Controls:
bwa_con_3R <- fread('bwa_episodic_3R_Con.csv')
novo_con_3R <- fread('novo_episodic_3R_Con.csv')
colnames(novo_con_3R) <- column.names
colnames(bwa_con_3R) <- column.names
bwa_con_3R <- na.omit(bwa_con_3R)
novo_con_3R <- na.omit(novo_con_3R)
novo_con_3R$Treatment <- 'Con'
bwa_con_3R$Treatment <- 'Con'
Con3R <- rbind(novo_con_3R, bwa_con_3R)
Con_3R <- Con3R %>%
  group_by(chr, Treatment, pos) %>%
  mutate(count = n())
Con_3R <- Con_3R[which(Con_3R$count==2),]
Con_3R <- Con_3R %>%
  group_by(chr, Treatment, pos) %>%
  summarise(meanSelCoef = (mean(selcoef)),
            pval_max=max(pval), pval_mean=mean(pval))

dat3R <- merge(Con_3R, Sel_3R, by = 'pos', all =T)
ee <- length(dat3R$pos)
#datX$number <- 1:aa
#dat2L$number <- (1+aa):(aa+bb)
#dat2R$number <- (1+aa+bb):(aa+bb+cc)
#dat3L$number <- (1+aa+bb+cc):(aa+bb+cc+dd)
dat3R$number <-  (aa+bb+cc+dd+1):(aa+bb+cc+dd+ee)
#dat4$number <- (aa+bb+cc+dd+ee+1):(aa+bb+cc+dd+ee+ff)
Con_3R <- dat3R[c(1,2,3,4,5,6,12)]
Sel_3R <- dat3R[c(1,7,8,9,10,11,12)]
column.names_2 <- c('pos','chr', 'Treatment', 'meanSelCoef', 'pval_max','pval_mean', 'number')
colnames(Sel_3R) <- column.names_2
colnames(Con_3R) <- column.names_2

Xcx_3R <- rbind(Con_3R, Sel_3R)
rm(novo_con_3R, novo_sel_3R, bwa_con_3R, bwa_sel_3R)

#rm(list=ls()[! ls() %in% c('Xcx_2L', 'Xcx_2R', 'Xcx_3L','Xcx_3R','Xcx_4','Xcx_X', 'aa', 'bb', 'cc','dd','ee','ff')])

###-----      4       -----###
column.names <- c('selcoef', 'pval', 'pos', 'chr')

# Selection
novo_sel_4 <- fread('novo_episodic_4_Sel.csv')
bwa_sel_4 <- fread('bwa_episodic_4_Sel.csv')
colnames(novo_sel_4) <- column.names
colnames(bwa_sel_4) <- column.names
novo_sel_4 <- na.omit(novo_sel_4)
bwa_sel_4 <- na.omit(bwa_sel_4)
novo_sel_4$Treatment <- 'Sel'
bwa_sel_4$Treatment <- 'Sel'
Sel4 <- rbind(novo_sel_4, bwa_sel_4)

Sel_4 <- Sel4 %>%
  group_by(chr, Treatment, pos) %>%
  mutate(count = n())
Sel_4 <- Sel_4[which(Sel_4$count==2),]

Sel_4 <- Sel_4 %>%
  group_by(chr, Treatment, pos) %>%
  summarise(meanSelCoef = (mean(selcoef)),
            pval_max=max(pval), pval_mean=mean(pval))

# Controls:
bwa_con_4 <- fread('bwa_episodic_4_Con.csv')
novo_con_4 <- fread('novo_episodic_4_Con.csv')
colnames(novo_con_4) <- column.names
colnames(bwa_con_4) <- column.names
bwa_con_4 <- na.omit(bwa_con_4)
novo_con_4 <- na.omit(novo_con_4)
novo_con_4$Treatment <- 'Con'
bwa_con_4$Treatment <- 'Con'
Con4 <- rbind(novo_con_4, bwa_con_4)
Con_4 <- Con4 %>%
  group_by(chr, Treatment, pos) %>%
  mutate(count = n())
Con_4 <- Con_4[which(Con_4$count==2),]

Con_4 <- Con_4 %>%
  group_by(chr, Treatment, pos) %>%
  summarise(meanSelCoef = (mean(selcoef)),
            pval_max=max(pval), pval_mean=mean(pval))

dat4 <- merge(Con_4, Sel_4, by = 'pos', all =T)
ff <- length(dat4$pos)
#datX$number <- 1:aa
#dat2L$number <- (1+aa):(aa+bb)
#dat2R$number <- (1+aa+bb):(aa+bb+cc)
#dat3L$number <- (1+aa+bb+cc):(aa+bb+cc+dd)
#dat3R$number <-  (aa+bb+cc+dd+1):(aa+bb+cc+dd+ee)
dat4$number <- (aa+bb+cc+dd+ee+1):(aa+bb+cc+dd+ee+ff)
Con_4 <- dat4[c(1,2,3,4,5,6,12)]
Sel_4 <- dat4[c(1,7,8,9,10,11,12)]
column.names_2 <- c('pos','chr', 'Treatment', 'meanSelCoef', 'pval_max','pval_mean', 'number')
colnames(Sel_4) <- column.names_2
colnames(Con_4) <- column.names_2

Xcx_4 <- rbind(Con_4, Sel_4)
rm(novo_con_4, novo_sel_4, bwa_con_4, bwa_sel_4)


#rm(list=ls()[! ls() %in% c('Xcx_2L', 'Xcx_2R', 'Xcx_3L','Xcx_3R','Xcx_4','Xcx_X', 'aa', 'bb', 'cc','dd','ee','ff')])

###-----      Combine      -----###

Zxc <- rbind(Xcx_2L, Xcx_2R, Xcx_3L, Xcx_3R, Xcx_4, Xcx_X)
Zxc_lengths <- c(aa,bb,cc,dd,ee,ff)
#File with lengths for plots:
write.csv(Zxc_lengths, file="SelCoef_Full_LENGTHS.csv", row.names=FALSE)
#Full File:
write.csv(Zxc, file="SelCoef_Full.csv", row.names=FALSE)
