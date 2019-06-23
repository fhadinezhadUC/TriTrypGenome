library(gsubfn)
library(ggplot2)
library(seqinr)
library(gdata)
library(ggplot2)
library(reshape2)
library(plyr)
#source("https://bioconductor.org/biocLite.R")
library(Biostrings)
annotate.final.geneset <- function(genefilepath){
  # genefilepath is the path to the file integrated_tse_ara.txt
  genefilepath <-
    "/home/fatemeh/TriTrypGenome_project/GenePrediction/integrated_tse_ara.txt"
  genefile <-
    read.table(genefilepath, header = TRUE, colClasses = "character")
  # cuttoff score for TSE 50, for ARA 107
  genefile[genefile$arascore == "notfound", ]$arascore <- -1
  genefile$arascore <- as.integer(genefile$arascore)
  genefile[genefile$tsescore == "notfound", ]$tsescore <- -1
  genefile$tsescore <- as.integer(genefile$tsescore)
  
  dismiss0 <-
    (genefile$foundby == "both") &
    (genefile$arascore < 107 | genefile$tsescore < 50)
  dismiss1 <-
    (genefile$foundby == "ara") &
    (genefile$arascore < 107) # most of these genes were pseudo or truncated or both
  dismiss2 <- (genefile$foundby == "tse") & (genefile$tsescore < 50)
  genefile <- genefile[(!dismiss0 & !dismiss1 & !dismiss2), ]
  
  # genefunc column is added which shows the presented gene identity in table
  # Genes marked as ?? include: pseudo|truncated genes, genes with unmatched identity, genes with unassigned identity|anticodon by any of genefinders.
  
  ambiguty1 <-
    genefile$tsefunc == "" &  genefile$arafunc == "" # 2 genes
  ambiguty2 <-
    genefile$tsenote != "notfound" # 2 genes not the same 2 genes in ambiguty1
  ambiguty3 <-
    genefile$foundby == "both" &
    genefile$tsefunc != genefile$arafunc # 22 genes
  ambiguty4 <- logical(length = nrow(genefile))
  for (i in 1:nrow(genefile)) {
    if (genefile$foundby[i] != "ara")
      ambiguty4[i] <- length(grep("n|N", genefile$tsegeneseq[i])) == 1
    else
      ambiguty4[i] <- length(grep("n|N", genefile$arageneseq[i])) == 1
  }
  
  ambiguties <- ambiguty1 | ambiguty2 | ambiguty3 | ambiguty4
  genefile$genefunc <- ""
  genefile[ambiguties, ]$genefunc <- "??"
  
  #table(genefile[!ambiguties, ]$foundby)
  #ara both 
  #36 3525 
  
  genefile[!ambiguties, ]$genefunc <-
    genefile[!ambiguties, ]$arafunc
  
  # resctoring 10 genes with ambiguties Y|N
  # they were genes with exact same sequence. with relatively high TSE score of 72 and relatively low ARA score of 107.
  # 9 of these genes are from clade Tcruzi and complete 8 cruzi genomes which do not have any Y genes annotated for them in figure.
  genefile[genefile$genefunc=="??"&genefile$tsefunc=="Y",]$genefunc <- "Y"
  genefile
}
prepare.tsfm.input <- function(genefile){
  library(gsubfn)
  library(ggplot2)
  library(seqinr)
  tsfmInput <- genefile[genefile$genefunc!="??" & genefile$genefunc != "Z",]
  resultpath <- "/home/fatemeh/TriTrypGenome_project/GeneAnnotation/"
  write.table(tsfmInput,col.names = TRUE,file = paste(resultpath,"tsfm_input_geneset.txt",sep = ""))
}
create.summary.table <- function(genefile) {
  # create a table with columns: Annotation, Intersection, ARAonly, Union
  prepare.tsfm.input
  firstcol <- c("#tRNA","#N/#G","Min Gene Length","Max Gene Length","%intron","%G" ,"%C" ,"%T" ,"%A","A",  "C",  "D",  "E",  "F",  "G",  "H",  "I",  "K",  "L",  "M",  "N",  "P",  "Q",  "R",  "S",  "T",  "V",  "W",  "X",  "Y",  "Z", "??")
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
  AraOnly_ClassFreq <- as.data.frame(table(ARAonlyDF$genefunc))
  names(AraOnly_ClassFreq) <- c("class","freq")
  Union_ClassFreq <- as.data.frame(table(genefile$genefunc))
  names(Union_ClassFreq) <- c("class","freq")
  intersect_ClassFreq$class <- as.character(intersect_ClassFreq$class)
  AraOnly_ClassFreq$class <- as.character(AraOnly_ClassFreq$class)
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
  Latex.file.prep(resulttable,"resulttable")
  
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




#_______________________________________________

create.clusterDF <- function(tsfmGeneDF){
  geneset <- assignCluster(tsfmGeneDF,1000)
  m <- matrix(nrow = max(geneset$cluster),ncol = 9)
  clusterDF<- as.data.frame(m)
  names(clusterDF) <- c("sourceOrg","soureSeq","clusterID","clusterSeq1","clusterSeq2", "clusterBegin","clusterEnd","clusterDir","clusterLen")
  for (i in 1:max(geneset$cluster)) {
    clusterDF$sourceOrg[i] <- geneset[geneset$cluster==i,]$sourceOrg[1]
    clusterDF$soureSeq[i] <- geneset[geneset$cluster==i,]$sourceseq[1]
    clusterDF$clusterID[i] <- i
    clusterDF$clusterSeq1[i] <- paste(geneset[geneset$cluster==i,]$genecode,collapse = ",")
    clusterDF$clusterSeq2[i] <- paste(geneset[geneset$cluster==i,]$genefunc,collapse = "")
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
  filename <- paste("/home/fatemeh/Leishmania_2019/Leishmania_2019/Results/",clade,".txt",sep = "")
  writeLines(Llines, filename,sep = "\n")
  
}
