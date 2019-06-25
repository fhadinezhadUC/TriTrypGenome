library(gsubfn)
library(ggplot2)
library(seqinr)
library(gdata)
library(ggplot2)
library(reshape2)
library(plyr)
#source("https://bioconductor.org/biocLite.R")
library(Biostrings)
# annotate.final.geneset.round1() will read the main integrated gene file from integrated_tse_ara.txt
# Keep genes that either their ara score is > 106 or their tse score is > 49
# add column: genefun and fill it in the fallowing order:
                       #1. genes with same tse and ara identity are assign the same idenitity
                       #2. genes with different identity, pseudo|truncated genes, genes with un assigned 
                       #   identity and genes with letter N in their sequence are shown as ??
# add column: note and fill it in the fallowing order:
                       #1. genes with the same ara and tse marked as "T"
                       #2. genes with unassigned identity by both tse and ara are set as "UnAssigned"
                       #3. genes with letter N in their sequence are set as "ContainsN"
                       #4. genes with unmatched identity between ara and tse are set as "Undet" 
                       # these top 4 cases had no overlap!

annotate.final.geneset.round1 <- function(){
  genefilepath <-
    "/home/fatemeh/TriTrypGenome_project/GenePrediction/integrated_tse_ara.txt"
  genefile <-
    read.table(genefilepath, header = TRUE, colClasses = "character")
  #table(genefile$foundby)
  #ara both  tse 
  #750 3597   34 
  # set the ara score of genes not found by aragorn/tse to -1 to convert the column to integer
  genefile[genefile$arascore == "notfound",]$arascore <- -1
  genefile$arascore <- as.integer(genefile$arascore)
  genefile[genefile$tsescore == "notfound",]$tsescore <- -1
  genefile$tsescore <- as.integer(genefile$tsescore)
  # cuttoff score for TSE 50, for ARA 107
  dismiss0 <-
    (genefile$foundby == "both") &
    (genefile$arascore < 107 & genefile$tsescore < 50)
  dismiss1 <-
    (genefile$foundby == "ara") &
    (genefile$arascore < 107) # most of these genes were pseudo or truncated or both
  dismiss2 <- (genefile$foundby == "tse") & (genefile$tsescore < 50)
  dismissedGenes <- genefile[(dismiss0 | dismiss1 | dismiss2),]
  #table(dismissedGenes$foundby)
  #ara both  tse
  #714   18   33(28 of these tse only genes were pseudo| truncated genes)
  
  Anotated_genefile <-
    genefile[(!dismiss0 & !dismiss1 & !dismiss2),]
  #ara both  tse
  #36 3579    1
  
  # genefunc column is added which shows the presented gene identity in table
  # Genes marked as ?? in this column include: pseudo|truncated genes, genes with unmatched identity,
  # genes with unassigned identity|anticodon by any of genefinders, and genes with letter N|n in their sequence.
  
  ambiguty1 <-
    Anotated_genefile$tsefunc == "" &
    Anotated_genefile$arafunc == "" 
  ambiguty2 <-
    Anotated_genefile$tsenote != "notfound" # genes that are not markes as notfound by tse are pseudo|truncated
  ambiguty3 <-
    Anotated_genefile$foundby == "both" &
    Anotated_genefile$tsefunc != Anotated_genefile$arafunc 
  ambiguty4 <- logical(length = nrow(Anotated_genefile))
  for (i in 1:nrow(Anotated_genefile)) {
    if (Anotated_genefile$foundby[i] != "ara")
      ambiguty4[i] <-
        length(grep("n|N", Anotated_genefile$tsegeneseq[i])) == 1
    else
      ambiguty4[i] <-
        length(grep("n|N", Anotated_genefile$arageneseq[i])) == 1
  }
  
  ambiguties <- ambiguty1 | ambiguty2 | ambiguty3 | ambiguty4
  #table(Anotated_genefile[!ambiguties, ]$foundby)
  #ara both
  #36  3535
  Anotated_genefile$genefunc <- ""
  Anotated_genefile[!ambiguties,]$genefunc <-
    Anotated_genefile[!ambiguties,]$arafunc
  Anotated_genefile[ambiguties,]$genefunc <- "#"
  # 45 genes are marked as #from which 1 found by tse and 44 found by both
  # 3571 not marked as #
  # resolving the identity of some of the # genes
  # column note is added for the status of gene's identity
  Anotated_genefile$note <- "T"
  Anotated_genefile[ambiguty1, ]$note <- "UnAssigned" # 2 genes
  Anotated_genefile[ambiguty2, ]$note <-
    Anotated_genefile[ambiguty2, ]$tsenote # 4 truncated 2 pseudo
  Anotated_genefile[ambiguty4, ]$note <- "ContainsN" # 4 genes
  # top 3 lines do not overlap
  Anotated_genefile[ambiguty3, ]$note <-
    "Undet" # 33 genes found by both with unmatch identity
  #table(paste(Anotated_genefile[ambiguty3, ]$tsefunc,Anotated_genefile[ambiguty3, ]$arafunc,sep="|"))
  #E|L G|W I|D ?|L M|L M|O R|S Y|N 
  #1   3   3   3   9   2   1  11 
  
  # based on the results from annotate.final.geneset.round2 we will add 10 Y genes which have same folding by both genefinders
  # with relatively higher tse score than ara score. tse score of 72 (based on density figure it is amongst well predicted genes). 
  # and ara score of 108 which is very close to the aragorns lowscore genes.
  # based on paper: Gene organization and sequence analyses of transfer RNA genes in Trypanosomatid parasites
  # Analysis of the tRNA genes in Tritryps indicated that only the tRNA-Tyr genes contain an intron; which was previously reported 
  # in T. brucei. The intron is 11 bases long in L. major and T. brucei, and 13 bases long in T. cruzi; and as in other organisms, 
  # it is located between bases 37 and 38 (data not shown).
  # based on the TSE annotation the introns of length 13 are annotated between position 37 and 38, however it is annotated at site 35-36
  Ygenes <- Anotated_genefile$genefunc == "#" &
    Anotated_genefile$tsefunc == "Y" &
    Anotated_genefile$arascore > 106 &
    Anotated_genefile$tsescore > 50
  Anotated_genefile[Ygenes, ]$genefunc <- Anotated_genefile[Ygenes, ]$tsefunc
  Anotated_genefile[Ygenes, ]$note <- "resolved"
  Anotated_genefile
}
prepare.tsfm.input <- function(Anotated_genefile){
  library(gsubfn)
  library(ggplot2)
  library(seqinr)
  # table(Anotated_genefile$foundby)
  # ara both  tse 
  # 36 3579    1 
  tsfmInput <- Anotated_genefile[Anotated_genefile$genefunc!="#" & Anotated_genefile$genefunc != "Z",]
  # table(tsfmInput$foundby)
  # ara both 
  # 36 3459
  # remove genes of genomes TrangeliSC58(6 genes) and TcruziCLBrener(11 genes)
  lowqualityGenomes <- tsfmInput$sourceOrg=="TrangeliSC58" | tsfmInput$sourceOrg=="TcruziCLBrener"
  tsfmInput <- tsfmInput[!lowqualityGenomes,]
  # 3478 genes left to be alined
  resultpath <- "/home/fatemeh/TriTrypGenome/GeneAnnotation/"
  write.table(tsfmInput,col.names = TRUE,file = paste(resultpath,"tsfm_input_geneset.txt",sep = ""))
}
create.summary.table <- function() {
  # create a table with columns: Annotation, Intersection, ARAonly, Union
  genefile <- annotate.final.geneset.round1()
  firstcol <- c("#tRNA","#N/#G","Min Gene Length","Max Gene Length","%intron","%G" ,"%C" ,"%T" ,"%A","A",  "C",  "D",  "E",  "F",  "G",  "H",  "I",  "K",  "L",  "M",  "N",  "P",  "Q",  "R",  "S",  "T",  "V",  "W",  "X",  "Y",  "Z", "#")
  resulttable <- data.frame(firstcol,rep(0,length(firstcol)),rep(0,length(firstcol)),rep(0,length(firstcol)))
  names(resulttable) <- c("Annotation", "Intersection", "ARAonly", "Union")
  
  isboth <- genefile$foundby=="both"
  isaraonly <- genefile$foundby=="ara"
  intersectDF <- genefile[isboth,]
  ARAonlyDF <- genefile[isaraonly,]
  genefile$geneseq <- ""
  genefile[isaraonly,]$geneseq <- genefile[isaraonly,]$arageneseq
  genefile[!isaraonly,]$geneseq <- genefile[!isaraonly,]$tsegeneseq
  
  resulttable[1,2:4] <- c(nrow(intersectDF),nrow(ARAonlyDF),nrow(genefile))
  resulttable[2,2:4] <- c(sum(nchar(intersectDF$tsegeneseq))/nrow(intersectDF),sum(nchar(ARAonlyDF$arageneseq))/nrow(ARAonlyDF),sum(nchar(genefile$geneseq))/nrow(genefile)) # sum of gene length / # genes from the first row
  resulttable[3,2:4] <- c(min(nchar(intersectDF$tsegeneseq)),min(nchar(ARAonlyDF$arageneseq)),min(nchar(genefile$geneseq)))
  resulttable[4,2:4] <- c(max(nchar(intersectDF$tsegeneseq)),max(nchar(ARAonlyDF$arageneseq)),max(nchar(genefile$geneseq)))
  resulttable[5, 2:4] <-
    c(
      100 * nrow(bothdf[bothdf$tseintronbegin != 0, ]) / nrow(bothdf),
      100 * nrow(ARAonlyDF[ARAonlyDF$araintronbegin != "nointron", ]) / nrow(ARAonlyDF),
      100 * (nrow(bothdf[bothdf$tseintronbegin != 0, ]) + nrow(ARAonlyDF[ARAonlyDF$araintronbegin !=
                                                                           "nointron", ])) / nrow(genefile)
    )
  resulttable[6,2:4] <- c(Gpercentage(paste(intersectDF$tsegeneseq,collapse = "")),Gpercentage(paste(ARAonlyDF$arageneseq,collapse = "")),Gpercentage(paste(genefile$geneseq,collapse = "")))
  resulttable[7,2:4] <- c(Cpercentage(paste(intersectDF$tsegeneseq,collapse = "")),Cpercentage(paste(ARAonlyDF$arageneseq,collapse = "")),Cpercentage(paste(genefile$geneseq,collapse = "")))
  resulttable[8,2:4] <- c(Tpercentage(paste(intersectDF$tsegeneseq,collapse = "")),Tpercentage(paste(ARAonlyDF$arageneseq,collapse = "")),Tpercentage(paste(genefile$geneseq,collapse = "")))
  resulttable[9,2:4] <- c(Apercentage(paste(intersectDF$tsegeneseq,collapse = "")),Apercentage(paste(ARAonlyDF$arageneseq,collapse = "")),Apercentage(paste(genefile$geneseq,collapse = "")))
  
  # make a table of Class frequencies for each set:
  intersect_ClassFreq <- as.data.frame(table(intersectDF$genefunc))
  names(intersect_ClassFreq) <- c("class","freq")
  intersect_ClassFreq$class <- as.character(intersect_ClassFreq$class)
  AraOnly_ClassFreq <- as.data.frame(table(ARAonlyDF$genefunc))
  names(AraOnly_ClassFreq) <- c("class","freq")
  AraOnly_ClassFreq$class <- as.character(AraOnly_ClassFreq$class)
  Union_ClassFreq <- as.data.frame(table(genefile$genefunc))
  names(Union_ClassFreq) <- c("class","freq")
  Union_ClassFreq$class <- as.character(Union_ClassFreq$class)
  
  for (i in 1:nrow(intersect_ClassFreq)) {
    curr_class <- intersect_ClassFreq$class[i]
    resulttable[resulttable$Annotation==curr_class,][,2] <- intersect_ClassFreq$freq[i]
  }
  for (i in 1:nrow(AraOnly_ClassFreq)) {
    curr_class <- AraOnly_ClassFreq$class[i]
    resulttable[resulttable$Annotation==curr_class,][,3] <- AraOnly_ClassFreq$freq[i]
  }
  for (i in 1:nrow(Union_ClassFreq)) {
    curr_class <- Union_ClassFreq$class[i]
    resulttable[resulttable$Annotation==curr_class,][,4] <- Union_ClassFreq$freq[i]
  }
  
  resulttable$Intersection <- round(resulttable$Intersection,digits = 0)
  resulttable$ARAonly <- round(resulttable$ARAonly,digits = 0)
  resulttable$Union <- round(resulttable$Union,digits = 0)
  resulttable$Annotation <- as.character(resulttable$Annotation)
  Latex.file.prep(resulttable,"/home/fatemeh/TriTrypGenome/Document_Latex/summarytable.txt")
  
}
clustersize.dist.visualize <- function(genefile){
  genefile <- annotate.final.geneset()
  clusterdf <- create.clusterDF(genefile)
  # dataframe c("clusterLen","Tfreq","Lfreq","Ofreq","percOfgenes")
  maxcluslen <- max(nchar(clusterdf$clusterSeq2))
  sizeDF <- data.frame(seq(1,maxcluslen,1),rep(0,maxcluslen))
  names(sizeDF) <- c("ClusterLen","Tfreq")
  sizeDF$Lfreq <- 0
  sizeDF$Ofreq <- 0
  sizeDF$percOfgenes <- 0
  for (i in 1: nrow(sizeDF)) {
    sizeDF$Tfreq[i] <- nrow(clusterdf[clusterdf$class=="T" & nchar(clusterdf$clusterSeq2) == i,])
    sizeDF$Lfreq[i] <- nrow(clusterdf[clusterdf$class=="L" & nchar(clusterdf$clusterSeq2) == i,])
    sizeDF$Ofreq[i] <- nrow(clusterdf[clusterdf$class=="O" & nchar(clusterdf$clusterSeq2) == i,])
    size_i_freq <- sizeDF$Tfreq[i] + sizeDF$Lfreq[i] + sizeDF$Ofreq[i]
    gene_num <- size_i_freq * i 
    sizeDF$percOfgenes[i] <- round((gene_num/nrow(genefile))*100,digits = 0)
  }
  
  sizeDF2 <- melt(sizeDF, id = c("ClusterLen","percOfgenes"))
  # Sort the data by ClusterLen and variable
  # Calculate the cumulative sum of the variable value for each cluster
  # Create the plot
  sizeDF2$ClusterLen <- factor(sizeDF2$ClusterLen)
  sizeDF2$OLTorder <- 0
  sizeDF2[sizeDF2$variable == "Ofreq", ]$OLTorder <- 0
  sizeDF2[sizeDF2$variable == "Lfreq", ]$OLTorder <- 1
  sizeDF2[sizeDF2$variable == "Tfreq", ]$OLTorder <- 2
  sizeDF2_sorted <-
    sizeDF2[order(sizeDF2$ClusterLen, sizeDF2$OLTorder), ]
  df_cumsum <- ddply(sizeDF2_sorted, "ClusterLen",
                     transform, label_ypos = cumsum(value))
  df_cumsum$label_ypos2 <- ""
  df_cumsum$label2 <- ""
  for (i in 1:nrow(df_cumsum)) {
    if (i %% 3 == 0)
    {
      df_cumsum$label_ypos2[i] <- toString(df_cumsum$label_ypos[i])
      df_cumsum$label2[i] <-
        paste(toString(df_cumsum$percOfgenes[i]), "%", sep = "")
    }
    else
    {
      df_cumsum$label_ypos2[i] = ""
      df_cumsum$label2[i] = ""
    }
  }
  df_cumsum$label1 <- ""
  df_cumsum$label1 <- df_cumsum$value
  df_cumsum[df_cumsum$value < 15,]$label1 <- ""
  ggplot(data = df_cumsum, aes(x = ClusterLen, y = value, fill = variable)) +
    geom_bar(stat = "identity") +
    geom_text(
      aes(y = label_ypos, label = label1),
      vjust = 1.3,
      color = "gray20",
      size = 4
    ) +
    scale_fill_brewer(palette = "Paired") + geom_text(aes(y = as.integer(label_ypos2), label =
                                                            label2),
                                                      vjust = -0.3,
                                                      size = 5,color = "darkgreen") + 
    scale_fill_manual(values = c("#00AFBB", "#E7B800", "lightgray"),name = "Genomes",labels = c("Trypanosoma", "Leishmania", "Others")) +
    theme_minimal() + xlab("Cluster Length") + ylab("Frequency")
  
}
clusterOrganization <- function(clusterDF,group1,filen){
  isbig <- nchar(clusterDF$clusterSeq) > 0
  bigclusterDF <- clusterDF[isbig,]
  bigclusterDF$clusterSeq <- paste(bigclusterDF$clusterSeq,bigclusterDF$clusterDir,sep = ",")
  bigclusterDF <- bigclusterDF[order(bigclusterDF$sourceOrg,bigclusterDF$soureSeq,bigclusterDF$clusterBegin),]
  # assign IDs to each cluster of same sourceOrg
  bigclusterDF$LocalClusID <- 0
  bigclusterDF$LocalClusID[1] <- 1
  
  bigclusterDF$seqnum <- 0
  bigclusterDF$seqnum[1] <- 1
  seqnumcount <- 1
  setnumber <- 1
  for (i in 1:(nrow(bigclusterDF)-1)) {
    if(bigclusterDF$sourceOrg[i+1]==bigclusterDF$sourceOrg[i])
    {
      setnumber <- setnumber + 1
      bigclusterDF$LocalClusID[i+1] <- setnumber
    }
    else
    {
      seqnumcount <- 1
      setnumber <- 1
      bigclusterDF$LocalClusID[i+1] <- setnumber
    }
    
    
    if(bigclusterDF$soureSeq[i + 1] != bigclusterDF$soureSeq[i])
    {
      seqnumcount <- seqnumcount + 1
      bigclusterDF$seqnum[i+1] <-  seqnumcount
    }
    else
    {
      bigclusterDF$seqnum[i+1] <- seqnumcount
    }
  }
  
  
  #############
  isgroup1 <- bigclusterDF$sourceOrg %in% group1
  bigclusterDF <- bigclusterDF[isgroup1,]
  ###########
  sources <- nrow(table(bigclusterDF$sourceOrg))
  #bigclusterDF$sourceOrg <- factor(bigclusterDF$sourceOrg)
  bigclusterDF$isbig <- FALSE
  # for (i in 1:nrow(bigclusterDF)) {
  #   bigclusterDF$isbig[i] <- nchar(bigclusterDF$clusterSeq2[i]) > 4
  # }
  # for (i in 1:nrow(bigclusterDF)) {
  #   bigclusterDF$isbig[i] <- all(grepl("[[:upper:]]", strsplit(bigclusterDF$clusterSeq2[i], "")[[1]]))
  # }
  
  bigclusterDF$col <- "0"
  bigclusterDF[(bigclusterDF$seqnum %% 2) == 0,]$col <- "1"
  #bigclusterDF[!bigclusterDF$isbig,]$col <- "2"
  g <- ggplot(
    data = bigclusterDF,
    aes(
      x = bigclusterDF$LocalClusID,
      y = bigclusterDF$sourceOrg,
      label = bigclusterDF$clusterSeq,
      color = bigclusterDF$col
    )
  ) + geom_point() + theme(
    legend.position = "none") + geom_label_repel(
      aes(label = bigclusterDF$clusterSeq),
      box.padding   = 0.1,
      point.padding = 0.1,
      segment.color = 'grey50',force = 1,nudge_x = 0.1
    )
  
  myheight <- sources * 2
  ggsave(file = paste("/home/fatemeh/TriTrypGenome/GeneAnnotation/",filen,".jpg",sep = ""), g,width = 100, height = myheight, units = "cm",limitsize = FALSE)
  
  
}
#_______________________________________________
# annotate.final.geneset.round2() will investigate genes marked as # in the genefun column to resolve the conflicts if possible, 
# based on evidence from clusters, relative genefiders scores, and their secondary sctructure
annotate.final.geneset.round2 <- function(){
  # summary of ambiguty3 as <tseScore|tseFun|araFun|araScore> and the second row as the frequency 
  #43|?|L|113 43|R|S|113 47|M|O|112 48|E|L|113 49|I|D|115 50|M|L|112 54|G|W|113 60|Y|N|101 72|Y|N|108 
  #   3          1          2          1          3          9          3          1         10 
  Anotated_genefile <- annotate.final.geneset.round1()
  # list all the ambiguties 
  ambiguty1 <-
    Anotated_genefile$tsefunc == "" &
    Anotated_genefile$arafunc == "" 
  ambiguty2 <-
    Anotated_genefile$tsenote != "notfound" # genes that are not markes as notfound by tse are pseudo|truncated
  ambiguty3 <-
    Anotated_genefile$foundby == "both" &
    Anotated_genefile$tsefunc != Anotated_genefile$arafunc 
  ambiguty4 <- logical(length = nrow(Anotated_genefile))
  for (i in 1:nrow(Anotated_genefile)) {
    if (Anotated_genefile$foundby[i] != "ara")
      ambiguty4[i] <-
        length(grep("n|N", Anotated_genefile$tsegeneseq[i])) == 1
    else
      ambiguty4[i] <-
        length(grep("n|N", Anotated_genefile$arageneseq[i])) == 1
  }
  Anotated_genefile2 <- Anotated_genefile
  Anotated_genefile2[ambiguty3,]$genefunc <- tolower(Anotated_genefile2[ambiguty3,]$tsefunc)
  Anotated_genefile2 <- Anotated_genefile2[!ambiguty2 & !ambiguty1 & !ambiguty4 ,] 
  # creating a text output to visualize the clusters containing ambiguty3
  clusterDF.textvisualization(Anotated_genefile2)
  group <- c("TcruzimarinkelleiB7","TcruziCLBrenerEsmeraldo-like","TcruziCLBrenerNon-Esmeraldo-like","TcruziEsmeraldo","TcruziSylvioX10-1-2012","TcruziSylvioX10-1","TcruziJRcl4", "TcruzicruziDm28c", "TcruziDm28c")
  # creating a figure of genomes that have will be completed by adding 10 Y genes
  clusterOrganization(clusterDF,group,"Ygenes")
}
# create.clusterDF takes a subset of annotated genefile as input. it creates a dataframe for clusters called clusterDF
# column clusterSeq is the sequence of each cluster made from single aminoacid letter of genes. 
# column clusterDir has a string of + and - which shows the strand of each gene within a cluster, in order.
create.clusterDF <- function(tsfmGeneDF){
  geneset <- assignCluster(tsfmGeneDF,1000)
  m <- matrix(nrow = max(geneset$cluster),ncol = 8)
  clusterDF<- as.data.frame(m)
  names(clusterDF) <- c("sourceOrg","soureSeq","clusterID","clusterSeq", "clusterBegin","clusterEnd","clusterDir","clusterLen")
  for (i in 1:max(geneset$cluster)) {
    clusterDF$sourceOrg[i] <- geneset[geneset$cluster==i,]$sourceOrg[1]
    clusterDF$soureSeq[i] <- geneset[geneset$cluster==i,]$sourceseq[1]
    clusterDF$clusterID[i] <- i
    clusterDF$clusterSeq[i] <- paste(geneset[geneset$cluster==i,]$genefunc,collapse = "")
    clusterDF$clusterBegin[i] <- sort(geneset[geneset$cluster==i,]$begin)[1]
    clusterDF$clusterEnd[i] <- sort(geneset[geneset$cluster==i,]$end)[nrow(geneset[geneset$cluster==i,])]
    clusterDF$clusterDir[i] <- paste(geneset[geneset$cluster==i,]$direction,collapse = "")
    clusterDF$clusterLen[i] <- nrow(geneset[geneset$cluster==i,])
  }
  clusterDF$class <- toupper(substr(clusterDF$sourceOrg,1,1))
  clusterDF[clusterDF$class!="T" & clusterDF$class!="L",]$class <- "O"
  clusterDF
}
assignCluster <- function(geneDF, clusterdistance) {
  
  geneDF$tsebegin <- as.integer(geneDF$tsebegin)
  geneDF$tseend <- as.integer(geneDF$tseend)
  geneDF$arabegin <- as.integer(geneDF$arabegin)
  geneDF$araend <- as.integer(geneDF$araend)
  
  geneDF$begin <- geneDF$tsebegin
  geneDF[geneDF$foundby=="ara",]$begin <- geneDF[geneDF$foundby=="ara",]$arabegin
  geneDF$end <- geneDF$tseend
  geneDF[geneDF$foundby=="ara",]$end <- geneDF[geneDF$foundby=="ara",]$araend
  
  geneDF <-
    geneDF[order(geneDF$sourceOrg, geneDF$sourceseq, geneDF$begin), ]
  geneDF$cluster <- 0
  geneDF$cluster[1] <- 1
  setnumber <- 1
  
  for (i in 1:(nrow(geneDF) - 1)) {
    geneDist <-
      abs(geneDF$begin[i + 1] - geneDF$begin[i])
    if ((geneDist < clusterdistance) &&
        (geneDF$sourceOrg[i] == geneDF$sourceOrg[i + 1]) &&
        (geneDF$sourceseq[i + 1] == geneDF$sourceseq[i]))
      geneDF$cluster[i + 1] <- setnumber
    else
    {
      setnumber <- setnumber + 1
      geneDF$cluster[i + 1] <- setnumber
    }
  }
  geneDF
}
clusterDF.textvisualization <- function(Anotated_genefile2){
  clusterDF <- create.clusterDF(Anotated_genefile2)
  clusterDF2 <- clusterDF[nchar(clusterDF$clusterSeq) > 1,]
  clusters <- names(table(clusterDF2$clusterSeq))
  clusterfreq <- as.integer(table(clusterDF2$clusterSeq))
  clusterfreqdf <- data.frame(clusters,clusterfreq)
  cluster_set <-integer(length = length(clusters))
  clusterfreqdf$clusters <- as.character(clusterfreqdf$clusters)
  clusterfreqdf <- clusterfreqdf[order(nchar(clusterfreqdf$clusters)),]
  setnumber <- 1
  for (i in 1:length(clusters)) {
    for (j in 1:length(clusters)) {
      S1 <-
        substring(tolower(clusterfreqdf$clusters[i]), seq(1, nchar(clusterfreqdf$clusters[i]), 1), seq(1, nchar(clusterfreqdf$clusters[i]), 1))
      S2 <-
        substring(tolower(clusterfreqdf$clusters[j]), seq(1, nchar(clusterfreqdf$clusters[j]), 1), seq(1, nchar(clusterfreqdf$clusters[j]), 1))
      C1 <- setdiff(S1, S2)
      C2 <- setdiff(S2, S1)
      C3 <- intersect(S1,S2)
      if ((length(C1) == 0 | length(C2) == 0) & length(S1) > 2 & length(S2) > 2) {
        if (cluster_set[i] != 0 & cluster_set[j] == 0)
        {
          cluster_set[j] <-  cluster_set[i]
          #setnumber = setnumber + 1
        }
        
        if (cluster_set[j] != 0 & cluster_set[i] == 0)
        {
          cluster_set[i] <-  cluster_set[j]
          #setnumber = setnumber + 1
        }
        if (cluster_set[j] != 0 & cluster_set[i] != 0)
        {
          cluster_set[i] <-  cluster_set[j]
          #setnumber = setnumber + 1
        }
        if (cluster_set[j] == 0 & cluster_set[i] == 0)
        {
          cluster_set[i] = setnumber
          cluster_set[j] = setnumber
          setnumber = setnumber + 1
        }
      }
    }
  }
  clusterfreqdf$cluster_set <- cluster_set
  cluster_setDF <- clusterfreqdf
  cluster_setDF <- cluster_setDF[order(cluster_setDF$cluster_set),]
  # for each cluster add the genomeclasses that have it 
  # and the direction of each classes
  cluster_setDF$Genomeclass <- ""
  cluster_setDF$Tdirs <- ""
  cluster_setDF$Ldirs <- ""
  cluster_setDF$Odirs <- ""
  
  for (i in 1:nrow(cluster_setDF)) {
    currclus <- cluster_setDF$clusters[i]
    cluster_setDF$Genomeclass[i] <- paste(unique(clusterDF2[clusterDF2$clusterSeq==currclus,]$class),collapse = "" )
    cluster_setDF$Tdirs[i] <- paste(unique(clusterDF2[clusterDF2$clusterSeq==currclus & clusterDF2$class=="T",]$clusterDir),collapse = "|" )
    cluster_setDF$Ldirs[i] <- paste(unique(clusterDF2[clusterDF2$clusterSeq==currclus & clusterDF2$class=="L",]$clusterDir),collapse = "|" )
    cluster_setDF$Odirs[i] <- paste(unique(clusterDF2[clusterDF2$clusterSeq==currclus & clusterDF2$class=="O",]$clusterDir),collapse = "|" )
  }
  cluster_setDF <- cluster_setDF[order(cluster_setDF$cluster_set,cluster_setDF$Genomeclass),]
  # for each cluster_set keep those that have at least one lower case letter in the clusters 
  cluster_setDF2 <- cluster_setDF
  cluster_setnumers <- unique(cluster_setDF2$cluster_set)
  for (i in 1:length(cluster_setnumers)) {
    j <- cluster_setnumers[i]
    currclusters <- as.character(cluster_setDF2[cluster_setDF2$cluster_set==j,]$clusters)
    cluStr <- paste(as.character(currclusters),collapse = "")
    if(all(grepl("[[:upper:]]", strsplit(cluStr, "")[[1]])))
    {cluster_setDF2 <- cluster_setDF2[cluster_setDF2$cluster_set!=j,]}
  }
  cluster_setDF2$clusters <- as.character(cluster_setDF2$clusters)
  cluster_setDF2 <- cluster_setDF2[order(cluster_setDF2$cluster_set,cluster_setDF2$clusters),]
  cluster_setDF2
}
gcContent <- function(x) {
  x <- DNAString(toupper(x))
  alf <- alphabetFrequency(x, as.prob = TRUE)
  sum(alf[c("G", "C")]) * 100
}
atContent <- function(x) {
  x <- DNAString(toupper(x))
  alf <- alphabetFrequency(x, as.prob = TRUE)
  sum(alf[c("A", "T")]) * 100
}
Apercentage <- function(x) {
  x <- DNAString(toupper(x))
  alf <- alphabetFrequency(x, as.prob = TRUE)
  sum(alf[c("A")]) * 100
}
Cpercentage <- function(x) {
  x <- DNAString(toupper(x))
  alf <- alphabetFrequency(x, as.prob = TRUE)
  sum(alf[c("C")]) * 100
}
Gpercentage <- function(x) {
  x <- DNAString(toupper(x))
  alf <- alphabetFrequency(x, as.prob = TRUE)
  sum(alf[c("G")]) * 100
}
Tpercentage <- function(x) {
  x <- DNAString(toupper(x))
  alf <- alphabetFrequency(x, as.prob = TRUE)
  sum(alf[c("T")]) * 100
}


Latex.file.prep <- function(df,filename){
  firstline <- ""
  Llines <- paste(firstline, collapse = "&")
  Llines <- paste(Llines,"\\\\",sep="")
  #df$clusters <- as.character(df$clusters)
  for (i in 1:nrow(df)) {
    currline <- paste(df[i,],collapse = "&")
    currline <- paste(currline,"\\\\",sep="")
    Llines <- c(Llines,currline)
  }
  writeLines(Llines, filename,sep = "\n")
  
}
