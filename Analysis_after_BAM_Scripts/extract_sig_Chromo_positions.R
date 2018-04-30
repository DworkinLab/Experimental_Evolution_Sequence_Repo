# Extract positions of interest:
## Assuming a list or data frame is available with all positions: 4 not done (no positions)

require('data.table')

name.Columns <- c("Chromosome", "Position", "ref", 
                  "ConR1_115", "ConR2_115", "SelR1_115", "SelR2_115", 
                  "ConR1_38", "ConR2_38", "SelR1_38", "SelR2_38", 
                  "ConR1_77", "ConR2_77", "SelR1_77", "SelR2_77", 
                  "SelR1_0")

### X Chromosome:
siglist <- fread('/home/paul/Positions/X_positions_maxP.csv', h=T)
siglist_2 <- as.list(siglist$int)
## Read in data:
episodic_data_bwa <- fread('/home/paul/episodicData/mpileup_dir/episodic_data_X.sync')
episodic_data_novo <- fread('/home/paul/episodicData/novoalign/novo_mpileup/novo_episodic_X.sync')
episodic_data_bowtie <- fread('/home/paul/episodicData/bowtie/mpileup_dir/bowtie2_episodic_X.sync')

colnames(episodic_data_bwa) <- name.Columns
colnames(episodic_data_novo) <- name.Columns
colnames(episodic_data_bowtie) <- name.Columns
## Keep rows with Position in list

episodic_counts_bwa <- episodic_data_bwa[episodic_data_bwa$Position %in% siglist_2 ,]
episodic_counts_novo <- episodic_data_novo[episodic_data_novo$Position %in% siglist_2 ,]
episodic_counts_bowtie <- episodic_data_bowtie[episodic_data_bowtie$Position %in% siglist_2 ,]

## write csv file with only positions
write.csv(episodic_counts_bwa, file="/home/paul/Positions/bwa_X_positions.sync")
write.csv(episodic_counts_novo, file="/home/paul/Positions/novo_X_positions.sync")
write.csv(episodic_counts_bowtie, file="/home/paul/Positions/bowtie_X_positions.sync")


### 2R Chromosome:
siglist <- fread('/home/paul/Positions/2R_positions_maxP.csv', h=T)
siglist_2 <- as.list(siglist$int)

episodic_data_bwa <- fread('/home/paul/episodicData/mpileup_dir/episodic_data_2R.sync')
episodic_data_novo <- fread('/home/paul/episodicData/novoalign/novo_mpileup/novo_episodic_2R.sync') # Does not exist
episodic_data_bowtie <- fread('/home/paul/episodicData/bowtie/mpileup_dir/bowtie2_episodic_2R.sync') # Does not exist

colnames(episodic_data_bwa) <- name.Columns
colnames(episodic_data_novo) <- name.Columns
colnames(episodic_data_bowtie) <- name.Columns
## Keep rows with Position in list

episodic_counts_bwa <- episodic_data_bwa[episodic_data_bwa$Position %in% siglist_2 ,]
episodic_counts_novo <- episodic_data_novo[episodic_data_novo$Position %in% siglist_2 ,]
episodic_counts_bowtie <- episodic_data_bowtie[episodic_data_bowtie$Position %in% siglist_2 ,]

## write csv file with only positions
write.csv(episodic_counts_bwa, file="/home/paul/Positions/bwa_2R_positions.sync")
write.csv(episodic_counts_novo, file="/home/paul/Positions/novo_2R_positions.sync")
write.csv(episodic_counts_bowtie, file="/home/paul/Positions/bowtie_2R_positions.sync")


### 2L Chromosome:
siglist <- fread('/home/paul/Positions/2L_positions_maxP.csv', h=T)
siglist_2 <- as.list(siglist$int)

episodic_data_bwa <- fread('/home/paul/episodicData/mpileup_dir/episodic_data_2L.sync')
episodic_data_novo <- fread('/home/paul/episodicData/novoalign/novo_mpileup/novo_episodic_2L.sync') # Does not exist
episodic_data_bowtie <- fread('/home/paul/episodicData/bowtie/mpileup_dir/bowtie2_episodic_2L.sync') # Does not exist

colnames(episodic_data_bwa) <- name.Columns
colnames(episodic_data_novo) <- name.Columns
colnames(episodic_data_bowtie) <- name.Columns

episodic_counts_bwa <- episodic_data_bwa[episodic_data_bwa$Position %in% siglist_2 ,]
episodic_counts_novo <- episodic_data_novo[episodic_data_novo$Position %in% siglist_2 ,]
episodic_counts_bowtie <- episodic_data_bowtie[episodic_data_bowtie$Position %in% siglist_2 ,]

write.csv(episodic_counts_bwa, file="/home/paul/Positions/bwa_2L_positions.sync")
write.csv(episodic_counts_novo, file="/home/paul/Positions/novo_2L_positions.sync")
write.csv(episodic_counts_bowtie, file="/home/paul/Positions/bowtie_2L_positions.sync")




### 3R Chromosome:
siglist <- fread('/home/paul/Positions/3R_positions_maxP.csv', h=T)
siglist_2 <- as.list(siglist$int)

episodic_data_bwa <- fread('/home/paul/episodicData/mpileup_dir/episodic_data_3R.sync')
episodic_data_novo <- fread('/home/paul/episodicData/novoalign/novo_mpileup/novo_episodic_3R.sync') # Does not exist
episodic_data_bowtie <- fread('/home/paul/episodicData/bowtie/mpileup_dir/bowtie2_episodic_3R.sync') # Does not exist

colnames(episodic_data_bwa) <- name.Columns
colnames(episodic_data_novo) <- name.Columns
colnames(episodic_data_bowtie) <- name.Columns

episodic_counts_bwa <- episodic_data_bwa[episodic_data_bwa$Position %in% siglist_2 ,]
episodic_counts_novo <- episodic_data_novo[episodic_data_novo$Position %in% siglist_2 ,]
episodic_counts_bowtie <- episodic_data_bowtie[episodic_data_bowtie$Position %in% siglist_2 ,]

write.csv(episodic_counts_bwa, file="/home/paul/Positions/bwa_3R_positions.sync")
write.csv(episodic_counts_novo, file="/home/paul/Positions/novo_3R_positions.sync")
write.csv(episodic_counts_bowtie, file="/home/paul/Positions/bowtie_3R_positions.sync")

### 3L Chromosome:
siglist <- fread('/home/paul/Positions/3L_positions_maxP.csv', h=T)
siglist_2 <- as.list(siglist$int)

episodic_data_bwa <- fread('/home/paul/episodicData/mpileup_dir/episodic_data_3L.sync')
episodic_data_novo <- fread('/home/paul/episodicData/novoalign/novo_mpileup/novo_episodic_3L.sync') # Does not exist
episodic_data_bowtie <- fread('/home/paul/episodicData/bowtie/mpileup_dir/bowtie2_episodic_3L.sync') # Does not exist

colnames(episodic_data_bwa) <- name.Columns
colnames(episodic_data_novo) <- name.Columns
colnames(episodic_data_bowtie) <- name.Columns

episodic_counts_bwa <- episodic_data_bwa[episodic_data_bwa$Position %in% siglist_2 ,]
episodic_counts_novo <- episodic_data_novo[episodic_data_novo$Position %in% siglist_2 ,]
episodic_counts_bowtie <- episodic_data_bowtie[episodic_data_bowtie$Position %in% siglist_2 ,]

write.csv(episodic_counts_bwa, file="/home/paul/Positions/bwa_3L_positions.sync")
write.csv(episodic_counts_novo, file="/home/paul/Positions/novo_3L_positions.sync")
write.csv(episodic_counts_bowtie, file="/home/paul/Positions/bowtie_3L_positions.sync")

### 4 Chromosome: Empty!!!
