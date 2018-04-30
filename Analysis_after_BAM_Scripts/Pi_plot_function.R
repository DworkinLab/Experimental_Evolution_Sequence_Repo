Pi_PlotFunction <- function(x, y) {
  require(ggplot2)
  x2 <- gsub("\\_.*","",x)
  y2 <- y
  
  #bowtie.pi
  #novo.pi
  #Read in the data:
  Datt <- read.table(x)
  colnames(Datt) <- c('chr', 'window', 'windowCount', ' propInwindow', 'Pi')
  
  #Remove unnecessary regions: Not necessary based on later steps
  Datt$chr <- as.character(Datt$chr)
  Datt2 <- Datt
  
  #Remove "na" pi values
  Datt2 <- Datt2[-which(Datt2$Pi=="na"),]
  
  #Need the numbers for chromosomes for labelling and colours:
  DattX <- Datt2[which(Datt2$chr=="X"),]
  a <- dim(DattX)[1]
  DattX$number <- 1:a
  
  Datt2L <- Datt2[which(Datt2$chr=="2L"),]
  b <- dim(Datt2L)[1]
  Datt2L$number <- (a+1):(a+b)
  
  Datt2R <- Datt2[which(Datt2$chr=="2R"),]
  c <- dim(Datt2R)[1]
  Datt2R$number <- (a+b+1):(a+b+c)
  
  Datt3L <- Datt2[which(Datt2$chr=="3L"),]
  d <- dim(Datt3L)[1]
  Datt3L$number <- (a+b+c+1):(a+b+c+d)
  
  Datt3R <- Datt2[which(Datt2$chr=="3R"),]
  e <- dim(Datt3R)[1]
  Datt3R$number <- (a+b+c+d+1):(a+b+c+d+e)
  
  Datt4 <- Datt2[which(Datt2$chr=="4"),]
  f <- dim(Datt4)[1]
  Datt4$number <- (a+b+c+d+e+1):(a+b+c+d+e+f)
  
  #Full data frame of necessary chromosomes
  DattFull <- rbind(DattX, Datt2L, Datt2R, Datt3L, Datt3R, Datt4)
  
  #Pi as numeric
  DattFull$Pi=as.numeric(levels(DattFull$Pi))[DattFull$Pi]
  
  #Title:
  z2 <- paste(x2, y2, sep="_")
  
  # The plots: 
  Pi_plot <- ggplot(DattFull, aes(x = number, y= Pi, colour = chr)) 
  
  Pi_plot_2 <- Pi_plot + 
    geom_point(size=0.3, show.legend = F) +
    scale_y_continuous(limits=c(0, 0.02), breaks=seq(0, 0.02, 0.005)) + 
    xlab("") +
    scale_x_discrete(limits=c(1049, 3185, 5277, 7443, 9952, 11359), labels = c("X", "2L", '2R', '3L', '3R', "4")) +
    theme(text = element_text(size=20), 
          axis.text.x= element_text(size=15), axis.text.y= element_text(size=15)) +
    scale_colour_manual(values=c("#56B4E9", "#E69F00", 'grey30', 'grey46', 'wheat3', 'lemonchiffon4')) #+ggtitle(z2)
  
  return(Pi_plot_2)
}
