# this script has two function split_tRNAgene_with_genomename to split based on each individual genome and
# split_tRNAgene_with_clustername to split them based on clusters
# This script reads in the .fasta file containing aligned tRNA genes of HomoC and TryTryp genomes
# Extract the functional classes from Headers
# Makes a dataframe of "genomename", "sequence", "headers", "funclass"
# Splits the dataframe into a list of data frames based on genomename
# Writes each dataframe as a fasta file in a file with genome's names
# At the end it calls function vidualization passing the list of data frames to visualize the missing functional classes for each genome or cluster of genomes
# Need to run script splitFuncClass.sh to split each file based on functional classes

library(ggplot2)
library(gdata)
split_tRNAgene_with_genomename <- function() {
  resultpath <- "/home/fatemeh/TriTrypGenome_project/tsfm/input1/"
  fastafile <-
    read.table(
      "/home/fatemeh/TriTrypGenome_project/GeneAnnotation/Aligned_TriTryp_Homo.fasta",
      sep = "\n"
    )
  fastafile <- as.character(fastafile$V1)
  sequences <-
    as.character(fastafile[seq(2, length(fastafile), 2)])
  headers <- as.character(fastafile[seq(1, length(fastafile), 2)])
  headers <- gsub("_ara", "_", headers)
  headers <- gsub("\\s", "", headers)
  
  headers2 <-
    lapply(X = headers , function(X)
      gsub("Homo", "HOMO_", X, fixed = TRUE))
  headers2 <- as.character(headers2)
  
  genomeName <-
    lapply(X = headers2 , function(X)
      gsub(">", "", X, fixed = TRUE))
  genomeName <- as.character(genomeName)
  
  genomeName2 <-
    lapply(X = genomeName, function(X)
      unlist(strsplit(X, split = "_"))[1])
  genomeName2 <- as.character(genomeName2)
  
  
  functionalclasses <-
    lapply(X = headers, function(X)
      unlist(strsplit(X, split = "_"))[length(unlist(strsplit(X, split = "_")))])
  
  functionalclasses <- as.character(functionalclasses)
  
  tRNAdf <-
    data.frame(genomeName2, sequences, headers2, functionalclasses)
  names(tRNAdf) <-
    c("genomename", "sequence", "headers", "funclass")
  
  genome_list = split(tRNAdf, f = tRNAdf$genomename)
  
  # make a file for each genome and put the genomes in there
  for (i in 1:length(genome_list)) {
    # write the DF in a file to be processed by the bash
    filepath <-
      paste(resultpath,
            names(genome_list[i]),
            ".fasta",
            sep = "")
    
    write.fwf(
      data.frame(genome_list[[i]]$headers,
                 genome_list[[i]]$sequence),
      filepath,
      sep = "\n",
      colnames = FALSE
    )
    
  }
  vidualization(genome_list, resultpath)
}
split_tRNAgene_with_clustername <- function() {
  resultpath <- "/home/fatemeh/TriTrypGenome_project/tsfm/input2/"
  fastafile <-
    read.table(
      "/home/fatemeh/TriTrypGenome_project/GeneAnnotation/Aligned_TriTryp_Homo.fasta",
      sep = "\n"
    )
  fastafile <- as.character(fastafile$V1)
  sequences <-
    as.character(fastafile[seq(2, length(fastafile), 2)])
  headers <- as.character(fastafile[seq(1, length(fastafile), 2)])
  headers <- gsub("_ara", "_", headers)
  headers <- gsub("\\s", "", headers)
  
  headers2 <-
    lapply(X = headers , function(X)
      gsub("Homo", "HOMO_", X, fixed = TRUE))
  headers2 <- as.character(headers2)
  
  genomeName <-
    lapply(X = headers2 , function(X)
      gsub(">", "", X, fixed = TRUE))
  genomeName <- as.character(genomeName)
  
  genomeName2 <-
    lapply(X = genomeName, function(X)
      unlist(strsplit(X, split = "_"))[1])
  genomeName2 <- as.character(genomeName2)
  
  # we will split genomes based on clusters of genomes
  # add another column "clustername" to tRNAdf to assign the new names to each cluster of genomes:
  clusNames <- character(length = length(headers))
  clusNames <- genomeName2
  for (i in 1:length(clusNames)) {
    clusternames <- clusNames[i]
    if (clusternames == "LspMARLEM2494" |
        clusternames == "LenriettiiLEM3045")
      clusNames[i] <- "LenriettiComplex"
    if (clusternames == "TbruceigambienseDAL972" |
        clusternames == "TbruceiLister427" |
        clusternames == "TbruceiTREU927" |
        clusternames == "TevansiSTIB805" |
        clusternames == "TcongolenseIL3000" |
        clusternames == "TvivaxY486")
      clusNames[i] <- "AfricanTrypanosome"
    #clusternames == "TcruziCLBrener" clusternames == "TrangeliSC58" removed
    if (clusternames == "TgrayiANR4" |
        clusternames == "TcruziCLBrenerEsmeraldo-like" |
        clusternames == "TcruziCLBrenerNon-Esmeraldo-like" |
        clusternames == "TcruzicruziDm28c" |
        clusternames == "TcruziDm28c" |
        clusternames == "TcruziEsmeraldo" |
        clusternames == "TcruziJRcl4" |
        clusternames == "TcruzimarinkelleiB7" |
        clusternames == "TcruziSylvioX10-1" |
        clusternames == "TcruziSylvioX10-1-2012" |
        clusternames == "TcruziTulacl2")
      # clusternames == "TtheileriEdinburgh"
      clusNames[i] <- "AmericanTrypanosome"
    if (clusternames == "CfasciculataCfCl" |
        clusternames == "LseymouriATCC30220" |
        clusternames == "LpyrrhocorisH10")
      clusNames[i] <- "Leishmania1"
    if (clusternames == "LmajorFriedlin" |
        clusternames == "LmajorLV39c5" |
        clusternames == "LmajorSD75" |
        clusternames == "LturanicaLEM423" |
        clusternames == "LarabicaLEM1108" |
        clusternames == "LtropicaL590" |
        clusternames == "LaethiopicaL147" |
        clusternames == "LgerbilliLEM452")
      clusNames[i] <- "Leishmania3"
    if (clusternames == "LdonovaniBHU1220" |
        clusternames == "LdonovaniBPK282A1" |
        clusternames == "LinfantumJPCM5")
      clusNames[i] <- "Leishmania4"
    if (clusternames == "LamazonensisMHOMBR71973M2269" |
        clusternames == "LmexicanaMHOMGT2001U1103")
      clusNames[i] <- "Leishmania2"
    if (clusternames == "LbraziliensisMHOMBR75M2904" |
        clusternames == "LbraziliensisMHOMBR75M2903" |
        clusternames == "LpanamensisMHOMPA94PSC1" |
        clusternames == "LpanamensisMHOMCOL81L13")
      clusNames[i] <- "Leishmania5"
  }
  
  functionalclasses <-
    lapply(X = headers, function(X)
      unlist(strsplit(X, split = "_"))[length(unlist(strsplit(X, split = "_")))])
  
  functionalclasses <- as.character(functionalclasses)
  
  tRNAdf <-
    data.frame(genomeName2,
               sequences,
               headers2,
               functionalclasses,
               clusNames)
  names(tRNAdf) <-
    c("genomename", "sequence", "headers", "funclass", "clusNames")
  
  genome_list = split(tRNAdf, f = tRNAdf$clusNames)
  
  # make a file for each genome and put the genomes in there
  for (i in 1:length(genome_list)) {
    # write the DF in a file to be processed by the bash
    filepath <-
      paste(resultpath,
            names(genome_list[i]),
            ".fasta",
            sep = "")
    
    write.fwf(
      data.frame(genome_list[[i]]$headers,
                 genome_list[[i]]$sequence),
      filepath,
      sep = "\n",
      colnames = FALSE
    )
  }
  vidualization(genome_list, )
  
}
vidualization <- function(genome_list, plotpath) {
  tRNAcountdf <- data.frame(names(genome_list))
  tRNAcountdf$counts <- 0
  names(tRNAcountdf) <- c("genome", "tRNAcounts")
  for (i in 1:length(genome_list)) {
    tRNAcountdf[tolower(tRNAcountdf$genome) == tolower(names(genome_list[i])),]$tRNAcounts <-
      nrow(genome_list[[i]])
  }
  
  tRNAcountdf <- tRNAcountdf[order(tRNAcountdf$tRNAcounts), ]
  tRNAcountdf$genome <-
    factor(tRNAcountdf$genome, levels = tRNAcountdf$genome)
  tRNAcountdf$genome <-
    as.character(tRNAcountdf$genome, levels = tRNAcountdf$genome)
  
  func_classes <-
    c(
      "A",
      "R",
      "N" ,
      "D",
      "C" ,
      "Q",
      "E",
      "G",
      "H",
      "I",
      "L",
      "K",
      "M" ,
      "X",
      "F",
      "P",
      "S",
      "T",
      "W",
      "Y" ,
      "V"
    )
  tRNAcountdf$percent_21 <- 0
  tRNAcountdf$missing <- ''
  for (i in 1:length(genome_list)) {
    print(i)
    if (tolower(names(genome_list[i])) != "homo")
    {
      countsdf <-
        as.data.frame(table(as.character(genome_list[[i]]$funclass)))
      tRNAcountdf[tRNAcountdf$genome == names(genome_list[i]), ]$missing <-
        paste(setdiff(func_classes, as.character(countsdf$Var1)), collapse = " ")
      diffset <-
        setdiff(func_classes, as.character(countsdf$Var1))
      fcount <- 21 - length(diffset)
      tRNAcountdf[tRNAcountdf$genome == names(genome_list[i]), ]$percent_21 <-
        (fcount / 21) * 100
    }
    else
      tRNAcountdf[tolower(tRNAcountdf$genome) == tolower(names(genome_list[i])), ]$percent_21 <-
        100
    
  }
  tRNAcountdf$genome <-
    factor(tRNAcountdf$genome, levels = tRNAcountdf$genome)
  p <- ggplot(data = tRNAcountdf, aes(x = genome, y = percent_21)) +
    geom_bar(stat = "identity", fill = "#56B4E9")  + theme(
      axis.text.x = element_text(angle = 90, hjust = 1),
      plot.title = element_text(hjust = 0.5) ,
      text = element_text(size = 14)
    ) +
    geom_text(
      aes(y = tRNAcountdf$percent_21, label = tRNAcountdf$tRNAcounts),
      vjust = 1.3,
      color = "gray20",
      size = 4
    ) +
    labs(title =
           "Percentage of 22 tRNA functional classes covered by each cluster",
         x =
           "Genome", y = "Percentage of tRNA functional classes(total=22)") + geom_text(
             aes(label =
                   missing),
             vjust = 0.5,
             angle = 90,
             hjust = -0.1,
             color = "blueviolet",
             size = 4
           )
  p
  ggsave(
    paste(plotpath, "funcPerc_clustered2.png", sep = ""),
    width = 14,
    height = 9
  )
}
