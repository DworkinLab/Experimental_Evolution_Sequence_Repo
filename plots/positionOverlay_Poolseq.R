#Overlay Plots
require(data.table)
require(tidyverse)

xcs <- fread('../Data/Positions/candidatePos.csv')
#Positions

require(tidyverse)
require(data.table)

Zxc_sig <- fread('../Data/poolseq_outputs/poolseq_FDR_Sigpos.csv')
Zxc_length <- fread('../Data/poolseq_outputs/SelCoef_Full_LENGTHS.csv')
aa <- Zxc_length$x[1]
bb <- Zxc_length$x[2]
cc <- Zxc_length$x[3]
dd <- Zxc_length$x[4]
ee <- Zxc_length$x[5]
ff <- Zxc_length$x[6]

# 
Zxc_count <- Zxc_sig %>%
  group_by(chr, pos) %>%
  mutate(count = n())
Zxc_count2 <- Zxc_count[which(Zxc_count$count==1),]
Zxc_count2 <- Zxc_count2[which(Zxc_count2$Treatment=='Sel'),]

mmaaxx <- max(Zxc_count2$meanSelCoef)

Zxc_count2$sig <- ifelse(Zxc_count2$pos %in% xcs$pos & Zxc_count2$chr %in% xcs$chr, "Yes", 'No')


nnonSig_3R <- Zxc_count2[which(Zxc_count2$sig=="No"),]
ssig_3R <- Zxc_count2[which(Zxc_count2$sig=="Yes"),]

Chr_plot_sig <- ggplot(data = nnonSig_3R, aes(x=number, y=meanSelCoef, colour=chr)) + 
  geom_point(size = 0.8, alpha=0.5, show.legend = F) +
  scale_x_discrete(limits=c(aa/2, aa+(bb/2), (aa+bb+(cc/2)), (aa+bb+cc+(dd/2)), (aa+bb+cc+dd+(ee/2)), (aa+bb+cc+dd+ee+(ff/2))), labels = c('X', "2L", "2R", '3L', '3R', '4')) + 
  scale_colour_manual(values=c("#56B4E9", "#E69F00", 'grey30', 'grey46', 'wheat3', 'lemonchiffon4')) + 
  theme(text = element_text(size=20), 
        axis.text.x= element_text(size=15), 
        axis.text.y= element_text(size=15)) +
  ylab("SelCoef") +
  theme(panel.background = element_blank()) +
  geom_point(size=2, colour='black', alpha=1, show.legend = F, data=ssig_3R, aes(x=number, y=meanSelCoef)) +
  ylim(0,mmaaxx)

Chr_plot_sig
