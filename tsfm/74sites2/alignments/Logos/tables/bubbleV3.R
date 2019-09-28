# provide path to the folder that contains the tables
TableDir <- "/home/fatemeh/TriTrypGenome/tsfm/74sites2/alignments/Logos/tables/"
tables <- list.files(path = TableDir,pattern = "_Table.txt", recursive = FALSE)
# TableDir <- "/home/fatemeh/doc/DavesTable/"
# tables <- list.files(path = TableDir,pattern = "_table.txt", recursive = FALSE)
ldf.74 <- list()
for (i in 1:length(tables)) {
  tablepath <- paste(TableDir,tables[i],sep = "")
  ldf.74[[i]]     <- read.table(tablepath,header=T);
  names(ldf.74)[i] <- paste(gsub("_Table.txt","",tables[i]),".74",sep = "")
}
# MAJOR.74       <- read.table("MAJOR/MAJOR.74_table.txt",header=T);
# INFANTUM.74    <- read.table("INFANTUM/INFANTUM.74_table.txt",header=T);
# MEXICANA.74    <- read.table("MEXICANA/MEXICANA.74_table.txt",header=T);
# VIANNIA.74     <- read.table("VIANNIA/VIANNIA.74_table.txt",header=T);
# ENRIETTII.74   <- read.table("ENRIETTII/ENRIETTII.74_table.txt",header=T);
# LEPTOCRITH.74  <- read.table("LEPTOCRITH/LEPTOCRITH.74_table.txt",header=T);
# AFTRYP.74      <- read.table("AFTRYP/AFTRYP.74_table.txt",header=T);
# AMTRYP.74      <- read.table("AMTRYP/AMTRYP.74_table.txt",header=T);
# ldf.74   <- list(MAJOR.74,INFANTUM.74,MEXICANA.74,VIANNIA.74,AFTRYP.74,AMTRYP.74,LEPTOCRITH.74,ENRIETTII.74);
#classes  <- c('L','I','V','R','C','M','E','Q','Y','W','S','T','P','H','G','D','N','K','F','A');
classes <- levels(factor(ldf.74[[1]]$aa))
df.names <- c("Major (n = 8)","Infantum (n = 3)","Mexicana (n = 2)","Viannia (n = 4)","Af. Tryp. (n = 5)","Am. Tryp. (n = 11)","Lepto/Crith (n = 3)","Enriettii (n = 2)");

gains <- vector()
convs <- vector()
for (class in classes) {
  i <- 1
  filenm <- paste("bubble_",class,".pdf",sep="");
  for (df in ldf.74) {
    gainvals <- (df$gainbits * df$gainfht);
    convvals <- (df$convbits * df$convfht);
    gainvals <- log(gainvals[gainvals>0])
    convvals <- log(convvals[convvals>0])
    gains <- append(gains, values=gainvals);
    convs <- append(convs, values=convvals);
  }
}
bitdatagains <- data.frame(bits = gains,type=rep("Gain",length(gains)))
bitdataconvs <- data.frame(bits = convs,type=rep("Conversion",length(convs)))
bitdata <- rbind(bitdatagains,bitdataconvs)
log(0.5)
log(2)
library(ggplot2)
ggplot(bitdata,
       aes(x = bits,
           fill = type)) + theme_bw() + geom_density(alpha = 0.3) + xlab("log Gain or Conversion Bits") +
  ggtitle("log-letter-height distributions for Gains and Conversions") +
  theme(legend.title = element_blank()) + geom_vline(
    xintercept = log(0.5),
    color = "blue",
    size = 1.5
  )+geom_vline(
    xintercept = log(2),
    color = "red",
    size = 1.5
  )

library("heR.Misc"); ## THIS IS REQUIRED FOR THE BUBBLEPLOT FUNCTION AND MUST BE DOWNLOADED FROM
## http://exposurescience.org/her.html

# THESE ARE THE COORDINATES FOR THE TRNA STRUCTURE BACKBONE IN THE FIGURES

line.x <- c(6.875,6.500,6.125,5.750,5.375,5.000,4.625,4.625,5.000,5.000,2.375
            ,2.750,2.375,2.750,2.500,2.875,3.250,2.875,2.500,2.125,1.750,1.375,1.000
            ,0.625,0.250,0.625,1.000,1.500,1.125,1.500,1.125,1.500,1.125,1.500,1.125
            ,1.500,1.125,0.625,1.000,1.375,1.750,2.125,2.500,2.875,2.375,2.750,2.375
            ,2.750,2.375,4.250,4.250,3.875,4.250,3.875,4.250,4.250,3.875,3.500,3.125
            ,2.750,2.250,1.875,1.500,1.125,1.500,1.875,2.250,2.750,3.125,3.500,3.875
            ,4.250,4.625,5.000,5.375,5.750,6.125,6.500,6.875,7.250,7.625,8.000,8.375);

line.y <- c(8.875,8.500,8.875,8.500,8.875,8.500,8.875,7.375,7.000
            ,3.500,3.500,3.875,4.250,4.625,5.125,5.500,5.875,6.250,6.625
            ,7.000,7.375,7.000,6.625,6.250,5.875,5.500,5.125,4.625,4.250
            ,3.875,3.500,3.125,2.750,2.375,2.000,1.625,1.250,0.750,0.375
            ,0.000,-0.375,0.000,0.375,0.750,1.250,1.625,2.000,2.375,2.750
            ,2.750,3.875,4.250,4.625,5.000,5.375,8.500,8.875,8.500,8.875
            ,8.500,8.000,8.375,8.750,9.125,9.500,9.875,10.250,9.750,10.125
            ,9.750,10.125,9.750,10.125,9.750,10.125,9.750,10.125,9.750,10.125
            ,9.750,10.125,9.750,10.125);



## THESE ARE FOR THE SPRINZL COORD LABELS
coord.labels <- c("1","5","10","14","18","21","25","30","35","40","45","50","55","60","65","70");
coord.labels.x <- c(6.875,5.375,2.375,2.61,1.750,1.00,1.125,1.500,1.750,2.750,3.65,3.870,1.875,2.250,4.250,6.125);
coord.labels.y <- c(8.875,8.875,3.500,5.125,7.375,5.125,3.500,1.625,-0.375,1.625,4.680,8.875,8.375,10.250,9.750,10.125);

xbump <- 0.5;
ybump <- 0.5;

up <- c(13,14,15,16,17);
up.coord.labels <- coord.labels[up];
up.coord.labels.x <- coord.labels.x[up];
up.coord.labels.y <- coord.labels.y[up] + ybump;

dn <- c(1,2,5,9,11,12);
dn.coord.labels <- coord.labels[dn];
dn.coord.labels.x <- coord.labels.x[dn];
dn.coord.labels.y <- coord.labels.y[dn] - ybump;

lt <- c(3,4,7);
lt.coord.labels <- coord.labels[lt];
lt.coord.labels.x <- coord.labels.x[lt] - xbump;
lt.coord.labels.y <- coord.labels.y[lt];

rt <- c(6,8,10);
rt.coord.labels <- coord.labels[rt];
rt.coord.labels.x <- coord.labels.x[rt] + xbump;
rt.coord.labels.y <- coord.labels.y[rt];


alpha=0.5;
fact=0.5;
area=TRUE;
legend=FALSE;  

map2rgb <- function (c) { rgb(t(col2rgb(c))/255,alpha=alpha);}
colormap <- function (g,c) { 
  y <- rep(0,length(g));
  # y[g <  0.48            & c < 0.44]             <- rgb(t(col2rgb("white"))/255,alpha=alpha);
  # y[g >= 0.48 & g < 0.95 & c < 0.44]             <- rgb(t(col2rgb("darkred"))/255,alpha=alpha);
  # y[g >= 0.95            & c < 0.44]             <- rgb(t(col2rgb("red"))/255,alpha=alpha);
  # y[g <  0.48            & c >= 0.44 & c < 0.70] <- rgb(t(col2rgb("darkblue"))/255,alpha=alpha);
  # y[g >= 0.48 & g < 0.95 & c >= 0.44 & c < 0.70] <- rgb(t(col2rgb("darkmagenta"))/255,alpha=alpha);
  # y[g >= 0.95            & c >= 0.44 & c < 0.70] <- map2rgb("deeppink");
  # y[g <  0.48            & c >= 0.70]            <- map2rgb("blue");
  # y[g >= 0.48 & g < 0.95 & c >= 0.70]            <- map2rgb("blueviolet");
  # y[g >= 0.95            & c >= 0.70]            <- map2rgb("magenta");
  lg <- 0.8
  lc <- 0.2
  y[g <  lc             & c < lc]               <- "#FFFFFFC3" #rgb(t(col2rgb("white"))/255,alpha=alpha);
  y[g >= lc & g < lg   & c < lc]               <- "#B87082C3" #rgb(t(col2rgb("red4"))/255,alpha=alpha);
  y[g >= lg             & c < lc]               <- "#A2021DC3" #rgb(t(col2rgb("red"))/255,alpha=alpha);
  y[g <  lc             & c >= lc & c < lg]    <- "#8190AEC3" #rgb(t(col2rgb("royalblue4"))/255,alpha=alpha);
  y[g >= lc & g < lg   & c >= lc & c < lg]    <- "#765C8CC3" #rgb(t(col2rgb("darkorchid4"))/255,alpha=alpha);
  y[g >= lg             & c >= lc & c < lg]    <- "#B700B7C3" #rgb(t(col2rgb("darkorchid1"))/255,alpha=alpha);
  y[g <  lc             & c >= lg]              <- "#083EAEC3" #rgb(t(col2rgb("royalblue1"))/255,alpha=alpha); #map2rgb("blue");
  y[g >= lc & g < lg   & c >= lg]              <- "#190081C3" #rgb(t(col2rgb("purple4"))/255,alpha=alpha);" #map2rgb("blueviolet");
  y[g >= lg             & c >= lg]              <- "#6C008CC3" #rgb(t(col2rgb("purple1"))/255,alpha=alpha); #map2rgb("purple");
  y;
}


widthmap <- function (g,c) {
  y <- rep(1,length(g));
  y[g < 0.48 & c < 0.44] <- 1;
  y[g >= 0.48 | c >= 0.44] <- 2;
  y[g >= 0.95 | c >= 0.70] <- 3;    
  y;
}


for (class in classes) {
  i <- 1
  filenm <- paste("bubble_",class,".pdf",sep="");
  pdf(file=filenm,version="1.4",height=3.5);
  par(mfrow=c(2,4), las=1);
  for (df in ldf.74) {
    gains <- (df$gainbits * df$gainfht);
    convs <- (df$convbits * df$convfht);
    colors <- colormap(gains,convs);
    df$gains <- gains #%%%
    df$convs <- convs #%%%
    
    widths <- widthmap(gains,convs);
    op <- par(mar = rep(0.75, 4))
    bubbleplot(
      df$x[df$aa == class],
      df$y[df$aa == class],
      (df$fbits[df$aa == class] * df$fht[df$aa == class]),
      fact=fact, #0.265165 = sqrt(2*0.375^2)/2
      area=area, 
      fg = rgb(t(col2rgb("black")/255)),
      bg = colors[df$aa == class],
      box=FALSE,
      axes=FALSE,
      lwd=0.5,
      xlim = c(0,8),
      ylim = c(-0.8,10.8),
      main = paste(df.names[i],class)
    );
    par(op);
    lines(line.x,line.y);
    text(labels=up.coord.labels,x=up.coord.labels.x,y=up.coord.labels.y,cex=0.75);
    text(labels=dn.coord.labels,x=dn.coord.labels.x,y=dn.coord.labels.y,cex=0.75);
    text(labels=lt.coord.labels,x=lt.coord.labels.x,y=lt.coord.labels.y,cex=0.75);
    text(labels=rt.coord.labels,x=rt.coord.labels.x,y=rt.coord.labels.y,cex=0.75);
    prime.x <- c(line.x[1],line.x[83]);
    prime.y <- c(line.y[1],line.y[83]);
    text(labels=c("5'","3'"),x=prime.x,y=prime.y,adj=c(-0.9,0),cex=0.75);   
    if (df.names[i] == "Enriettii (n = 2)") {
      legend.x <- rep(df$x[df$aa == "X" &  df$state == "A" & 
                             (df$sprinzl == "69" |  df$sprinzl == "71" |
                                df$sprinzl == "73" )] + c(0.2,0.3,0.4),3);  
      legend.y <- c(rep(df$y[df$aa == "X" &  df$state == "A" & df$sprinzl == "35"]+0.40,3),
                    rep(df$y[df$aa == "X" &  df$state == "A" & df$sprinzl == "31"],3), 
                    rep(df$y[df$aa == "X" &  df$state == "A" & df$sprinzl == "27"]-0.30,3));
      legend.z <- rep(2.2,9);
      legend.c <- colormap(c(0,0,0,0.5,0.5,0.5,3,3,3),c(0,0.5,3,0,0.5,3,0,0.5,3));
      bubbleplot(legend.x, legend.y, legend.z,fact=fact,area=area,bg=legend.c,
                 add=TRUE,
                 box=FALSE,axes=FALSE,lwd=1);
    } 
    i <- i + 1
  }
  dev.off()
}


# all.bubble(MAJOR.72,name="MAJOR.v5/MAJOR.72");
# all.bubble(INFANTUM.72,name="INFANTUM.v5/INFANTUM.72");
# all.bubble(MEXICANA.72,name="MEXICANA.v5/MEXICANA.72");
# all.bubble(VIANNIA.72,name="VIANNIA.v5/VIANNIA.72");
# all.bubble(ENRIETTII.72,name="ENRIETTII.v5/ENRIETTII.72");
# all.bubble(LEPTOCRITH.72,name="LEPTOCRITH.v5/LEPTOCRITH.72");
# all.bubble(AMTRYP.72,name="AMTRYP.v5/AMTRYP.72");
# all.bubble(AFTRYP.72,name="AFTRYP.v5/AFTRYP.72");
# 
# 
# all.bubble(MAJOR.74,name="MAJOR.v5/MAJOR.74");
# all.bubble(INFANTUM.74,name="INFANTUM.v5/INFANTUM.74");
# all.bubble(MEXICANA.74,name="MEXICANA.v5/MEXICANA.74");
# all.bubble(VIANNIA.74,name="VIANNIA.v5/VIANNIA.74");
# all.bubble(ENRIETTII.74,name="ENRIETTII.v5/ENRIETTII.74");
# all.bubble(LEPTOCRITH.74,name="LEPTOCRITH.v5/LEPTOCRITH.74");
# all.bubble(AMTRYP.74,name="AMTRYP.v5/AMTRYP.74");
# all.bubble(AFTRYP.74,name="AFTRYP.v5/AFTRYP.74")