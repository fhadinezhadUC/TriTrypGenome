Integrate_Tse_Ara <- function() {
  #http://www.bioconductor.org/packages/release/bioc/vignettes/GenomicRanges/inst/doc/GenomicRangesIntroduction.pdf
  #https://support.bioconductor.org/p/56123/

  library(IRanges)
  library(stringr)
  library(miceadds)
  library(GenomicRanges)
  library(gdata)
  # read ara
  resultpath <- "/home/fatemeh/Leishmania_2019/Leishmania_2019/Results/Integrated_Genes/"
  ara_filenames <-
    list.files(
      "/home/fatemeh/Leishmania_2019/Leishmania_2019/Results/AragornGenes/",
      pattern = "*.ara.out",
      full.names = TRUE
    )
  #just take the ara_filenames and make the tse file for it to read.
  tse_filenames <- ara_filenames
  tse_ss_filenames <- ara_filenames
  for (i in 1:length(ara_filenames)) {
    namearray <- unlist(strsplit(ara_filenames[i], split = "/"))
    nameorg <- namearray[length(namearray)]
    nameorg <- unlist(strsplit(nameorg, split = ".ara.out"))[1]
    
    tse_filenames[i] = paste(
      "/home/fatemeh/Leishmania_2019/Leishmania_2019/Results/tRNAScanGenes/tse.out/",
      nameorg,
      ".tse.out",
      sep = ""
    )
    tse_ss_filenames[i] = paste(
      "/home/fatemeh/Leishmania_2019/Leishmania_2019/Results/tRNAScanGenes/ss.tse.out/",
      nameorg,
      ".SS.tse.out",
      sep = ""
    )
  }
  
  
  aradf <- making_ara_df(ara_filenames[1])
  tsedf <-
    making_tse_df(tse_filenames[1], tse_ss_filenames[1])
  integrated_tse_ara <- integrate(aradf, tsedf)
  
  for (i in 2:length(ara_filenames)) {
    aradf <- making_ara_df(ara_filenames[i])
    tsedf <-
      making_tse_df(tse_filenames[i], tse_ss_filenames[i])
    temp_integrated <- integrate(aradf, tsedf)
    integrated_tse_ara <-
      rbind(integrated_tse_ara, temp_integrated)
  }
  integrated_tse_ara <- addfuncTypes(integrated_tse_ara)
  # I think I did the alignment before this initiator detection step!
  integrated_tse_ara <- initiatorDetecting(integrated_tse_ara)
  
  # decide about the functional class of genes found by both TSE and ARA with different identities. 
  integrated_tse_ara <- fix.Mismatch.Identity(integrated_tse_ara)
  formatoutput(integrated_tse_ara,resultpath)
  # the out put is written in resultpath with name integrated_tse_ara.txt 
  # everytime we need the integrated file we can read is as this:
  # geneDF <- read.table(file = paste(resultpath,"integrated_tse_ara.txt",sep = ""), header = TRUE,colClasses = "character")
  integrated_tse_ara
  
}
# ____________________ formating the final integrated tse ara file with a fixed length format and write it into four file _____________________________

formatoutput <- function(integrated_tse_ara,resultpath) {
  # writing the dataframe in two files :
  # file1 for coordinates: geneid, sourceOrg, sourceseq, sourceSO, direction, tse/aracoordinate, foundby
  # file2for SS:          geneid:\n tse: tsegeneseq \n tseSS/ arageneseq, araSS /
  # file3 for identities:  geneid, tse/araidentity, tse/araac, tseacloc, arascore/tsescore, tsenote, aranote
  # file4 for intron:      geneid, tseintroncoordinate
  
  
  n <- data.frame(
    "GeneId",
    "TseIntronBegin",
    "TseIntronEnd",
    "AraIntronLocStart",
    "AraIntronLocEnd"
  )
  names(n) <- c(
    "GeneId",
    "TseIntronBegin",
    "TseIntronEnd",
    "AraIntronLocStart",
    "AraIntronLocEnd"
  )
  introndf <-  data.frame(
    integrated_tse_ara$geneid,
    integrated_tse_ara$tseintronbegin,
    integrated_tse_ara$tseintronend,
    integrated_tse_ara$araintronbegin,
    integrated_tse_ara$araintronend
  )
  names(introndf) <- names(n)
  
  write.fwf(
    rbind(n, introndf),
    width = c(60, 20, 20, 20, 20),
    colnames = FALSE,
    file = paste(resultpath,"Intron.txt",sep = "")
  )
  #_________________________________________________________________________________
 
  #positive strand
  # isfivematchthreenot <- (Fiveprimend == 0) & (Threeprimend == -1)
  # fivematchthreenot <- both[isfivematchthreenot,]
  # 
  # write.fwf(
  #   data.frame(
  #     paste("ara: ", fivematchthreenot$arageness, fivematchthreenot$arabegin, fivematchthreenot$araend, sep = "-"),
  #     paste("ara: ", fivematchthreenot$arageneseq, sep = "-"),
  #     paste("tse: ", fivematchthreenot$tsegeness, fivematchthreenot$tsebegin, fivematchthreenot$tseend, sep = "-"),
  #     paste("tse: ", fivematchthreenot$tsegeneseq, sep = "-")
  #   ),
  #   file = "/home/fatemeh/Leishmania_Aug2018/fivematchthree-1.txt",
  #   sep = "\n",
  #   colnames = FALSE
  # )
  # 
  # #negative strand
  # isfivematchthreenot <- (Fiveprimend == 1) & (Threeprimend == -1)
  # fivematchthreenot <- both[isfivematchthreenot,]
  # 
  # write.fwf(
  #   data.frame(
  #     paste("ara: ", fivematchthreenot$arageness, fivematchthreenot$arabegin, fivematchthreenot$araend, sep = "-"),
  #     paste("ara: ", fivematchthreenot$arageneseq, sep = "-"),
  #     paste("tse: ", fivematchthreenot$tsegeness, fivematchthreenot$tsebegin, fivematchthreenot$tseend, sep = "-"),
  #     paste("tse: ", fivematchthreenot$tsegeneseq, sep = "-")
  #   ),
  #   file = "/home/fatemeh/Leishmania_Aug2018/five1three-1.txt",
  #   sep = "\n",
  #   colnames = FALSE
  # )
  # 
  # 
  # isfivematchthreenot <- (Fiveprimend == 0) & (Threeprimend == 0)
  # fivematchthreenot <- both[isfivematchthreenot,]
  # 
  # write.fwf(
  #   data.frame(
  #     paste("ara: ", fivematchthreenot$arageness, fivematchthreenot$arabegin, fivematchthreenot$araend, sep = "-"),
  #     paste("ara: ", fivematchthreenot$arageneseq, sep = "-"),
  #     paste("tse: ", fivematchthreenot$tsegeness, fivematchthreenot$tsebegin, fivematchthreenot$tseend, sep = "-"),
  #     paste("tse: ", fivematchthreenot$tsegeneseq, sep = "-")
  #   ),
  #   file = "/home/fatemeh/Leishmania_Aug2018/fivethreematches.txt",
  #   sep = "\n",
  #   colnames = FALSE
  # )
  # table(substring(fivematchthreenot$arageneseq,nchar(fivematchthreenot$arageneseq)-2,nchar(fivematchthreenot$arageneseq)))
  # table(fivematchthreenot$direction)
  # table(fivematchthreenot$tseidentity == fivematchthreenot$araidentity)
  # 
  # isfivematchthreenot <- (Fiveprimend == 1) & (Threeprimend == -2)
  # fivematchthreenot <- both[isfivematchthreenot,]
  # 
  # write.fwf(
  #   data.frame(
  #     paste("ara: ", fivematchthreenot$arageness, fivematchthreenot$arabegin, fivematchthreenot$araend, sep = "-"),
  #     paste("ara: ", fivematchthreenot$arageneseq, sep = "-"),
  #     paste("tse: ", fivematchthreenot$tsegeness, fivematchthreenot$tsebegin, fivematchthreenot$tseend, sep = "-"),
  #     paste("tse: ", fivematchthreenot$tsegeneseq, sep = "-")
  #   ),
  #   file = "/home/fatemeh/Leishmania_Aug2018/five1three-2",
  #   sep = "\n",
  #   colnames = FALSE
  # )
  # 
  # isfivematchthreenot <- (Fiveprimend == 0) & (Threeprimend == -2)
  # fivematchthreenot <- both[isfivematchthreenot,]
  # 
  # write.fwf(
  #   data.frame(
  #     paste("ara: ", fivematchthreenot$arageness, fivematchthreenot$arabegin, fivematchthreenot$araend, sep = "-"),
  #     paste("ara: ", fivematchthreenot$arageneseq, sep = "-"),
  #     paste("tse: ", fivematchthreenot$tsegeness, fivematchthreenot$tsebegin, fivematchthreenot$tseend, sep = "-"),
  #     paste("tse: ", fivematchthreenot$tsegeneseq, sep = "-")
  #   ),
  #   file = "/home/fatemeh/Leishmania_Aug2018/five0three-2",
  #   sep = "\n",
  #   colnames = FALSE
  # )
  # 
  # isfivematchthreenot <- (Fiveprimend == 0) & (Threeprimend == -3)
  # fivematchthreenot <- both[isfivematchthreenot,]
  # 
  # write.fwf(
  #   data.frame(
  #     paste("ara: ", fivematchthreenot$arageness, fivematchthreenot$arabegin, fivematchthreenot$araend, sep = "-"),
  #     paste("ara: ", fivematchthreenot$arageneseq, sep = "-"),
  #     paste("tse: ", fivematchthreenot$tsegeness, fivematchthreenot$tsebegin, fivematchthreenot$tseend, sep = "-"),
  #     paste("tse: ", fivematchthreenot$tsegeneseq, sep = "-")
  #   ),
  #   file = "/home/fatemeh/Leishmania_Aug2018/five0three-3",
  #   sep = "\n",
  #   colnames = FALSE
  # )
  # 
  # isfivematchthreenot <- (Fiveprimend == 1) & (Threeprimend == -3)
  # fivematchthreenot <- both[isfivematchthreenot,]
  # 
  # write.fwf(
  #   data.frame(
  #     paste("ara: ", fivematchthreenot$arageness, fivematchthreenot$arabegin, fivematchthreenot$araend, sep = "-"),
  #     paste("ara: ", fivematchthreenot$arageneseq, sep = "-"),
  #     paste("tse: ", fivematchthreenot$tsegeness, fivematchthreenot$tsebegin, fivematchthreenot$tseend, sep = "-"),
  #     paste("tse: ", fivematchthreenot$tsegeneseq, sep = "-")
  #   ),
  #   file = "/home/fatemeh/Leishmania_Aug2018/five1three-3",
  #   sep = "\n",
  #   colnames = FALSE
  # )
  #__________________________________________________________________________________
  n <- data.frame(
    "GeneId",
    "SourceOrg",
    "SourceSeq",
    "SourceSO",
    "Direction",
    "TseBegin",
    "TseEnd",
    "AraBegin",
    "AraEnd",
    "Foundby"
  )
  
  names(n) <- c(
    "GeneId",
    "SourceOrg",
    "SourceSeq",
    "SourceSO",
    "Direction",
    "TseBegin",
    "TseEnd",
    "AraBegin",
    "AraEnd",
    "Foundby"
  )
  coordinatedf <- data.frame(
    integrated_tse_ara$geneid,
    integrated_tse_ara$sourceOrg,
    integrated_tse_ara$sourceseq,
    integrated_tse_ara$sourceSO,
    integrated_tse_ara$direction,
    integrated_tse_ara$tsebegin,
    integrated_tse_ara$tseend,
    integrated_tse_ara$arabegin,
    integrated_tse_ara$araend,
    integrated_tse_ara$foundby
  )
  names(coordinatedf) <- names(n)
  write.fwf(
    rbind(n, coordinatedf),
    colnames = FALSE,
    width = c(57, 32, 32, 11, 10, 10, 10, 10, 10, 7),
    file = paste(resultpath,"Coordinate.txt",sep = "")
  )
  
  n <-
    data.frame(
      "GeneId",
      "TseIdentity",
      "AraIdentity",
      "TseAc",
      "AraAc",
      "TseAcLoc",
      "AraScore",
      "TseScore",
      "TseNote",
      "AraNote",
      "TseFunc",
      "AraFunc"
    )
  names(n) <-
    c(
      "GeneId",
      "TseIdentity",
      "AraIdentity",
      "TseAc",
      "AraAc",
      "TseAcLoc",
      "AraScore",
      "TseScore",
      "TseNote",
      "AraNote",
      "TseFunc",
      "AraFunc"
    )
  identitydf <- data.frame(
    integrated_tse_ara$geneid,
    integrated_tse_ara$tseidentity,
    integrated_tse_ara$araidentity,
    integrated_tse_ara$tseac,
    integrated_tse_ara$araac,
    integrated_tse_ara$tseacloc,
    integrated_tse_ara$arascore,
    integrated_tse_ara$tsescore,
    integrated_tse_ara$tsenote,
    integrated_tse_ara$aranote,
    integrated_tse_ara$tsefunc,
    integrated_tse_ara$arafunc
  )
  names(identitydf) <- names(n)
  write.fwf(
    rbind(n, identitydf),
    width = c(58, 12, 12, 15, 15, 12, 12, 12, 32, 10,7,7),
    file = paste(resultpath,"Identity.txt",sep = ""),
    colnames = FALSE
  )
  
  write.fwf(
    data.frame(
      paste("geneid: ", integrated_tse_ara$geneid, sep = ""),
      paste("tseseq: ", integrated_tse_ara$tsegeneseq, sep = ""),
      paste("tsess:  ", integrated_tse_ara$tsegeness, sep = ""),
      paste("araseq: ", integrated_tse_ara$arageneseq, sep = ""),
      paste("arass:  ", integrated_tse_ara$arageness, sep = ""),
      character(length = length(integrated_tse_ara$arageness))
    ),
    file = paste(resultpath,"SecondaryS.txt",sep = ""),
    sep = "\n",
    colnames = FALSE
  )
  
  # do we need to remove the tails of the sequences based on their end displacement, before the alignment? 
  genesequences <- character(length = nrow(integrated_tse_ara))
  for (i in 1:nrow(integrated_tse_ara)) {
    genesequences[i] <- integrated_tse_ara$tsegeneseq[i]
    if(integrated_tse_ara$tsegeneseq[i] == "notfound")
    genesequences[i] <- integrated_tse_ara$arageneseq[i]
  }
  
  integrated_tse_ara$genesequence <- genesequences
  write.fwf(
    data.frame(
      paste(">", integrated_tse_ara$geneid, sep = ""),
      integrated_tse_ara$genesequence
    ),
    file = paste(resultpath,"Integrated_Genes.fasta",sep = ""),
    sep = "\n",
    colnames = FALSE
  )
  
  # Write one file to have all the data to read everytime we need the genefile
  write.table(integrated_tse_ara,col.names = TRUE,file = paste(resultpath,"integrated_tse_ara.txt",sep = ""))
}
# ___________________ integrate (union) genes found by tse and ara into integrated_tse_ara dataframe __________________________________________________

integrate <- function(aradf, tsedf) {
  aradf$arabegin <- as.integer(as.character(aradf$arabegin))
  aradf$araend <- as.integer(as.character(aradf$araend))
  tsedf$tsebegin <- as.integer(as.character(tsedf$tsebegin))
  tsedf$tseend <- as.integer(as.character(tsedf$tseend))
  tsedf$tseintronbegin <- (as.character(tsedf$tseintronbegin))
  tsedf$tseintronend <- (as.character(tsedf$tseintronend))
  tsedf$tsesourceseq <- as.character(tsedf$tsesourceseq)
  aradf$arasourceseq <- as.character(aradf$arasourceseq)
  
  # for tse output if it is on reverse strand replace the begin and end
  bigger <- tsedf$tsebegin > tsedf$tseend
  for (i in 1:length(bigger)) {
    if (bigger[i])
    {
      temp <- tsedf$tsebegin[i]
      tsedf$tsebegin[i] <- tsedf$tseend[i]
      tsedf$tseend[i] <- temp
      
      temp <- tsedf$tseintronbegin[i]
      tsedf$tseintronbegin[i] <- tsedf$tseintronend[i]
      tsedf$tseintronend[i] <- temp
    }
  }
  
  aradf$aradirection = as.character(aradf$aradirection)
  tsedf$tsedirection = as.character(tsedf$tsedirection)
  aradf$arasourceseq <- as.character(aradf$arasourceseq)
  tsedf$tsesourceseq <- as.character(tsedf$tsesourceseq)
  
  ara <-
    GRanges(
      seqnames = aradf$arasourceseq ,
      ranges = IRanges(aradf$arabegin, aradf$araend),
      strand = Rle(strand(aradf$aradirection))
    )
  
  tse <-
    GRanges(
      seqnames = tsedf$tsesourceseq,
      ranges = IRanges(tsedf$tsebegin, tsedf$tseend),
      strand = Rle(strand(tsedf$tsedirection))
    )
  overlapps <- as.data.frame(findOverlaps(ara, tse))
  names(overlapps) <- c("ararecord", "tserecord")
  m <- matrix(ncol = ncol(aradf) + ncol(tsedf), nrow = 1)
  overlapdf <- as.data.frame(m)
  names(overlapdf) <- c(names(aradf), names(tsedf))
  for (i in 1:nrow(overlapps)) {
    tempbind <-
      as.data.frame(c(aradf[overlapps$ararecord[i],], tsedf[overlapps$tserecord[i],]))
    overlapdf <- rbind(tempbind, overlapdf)
  }
  integrated_gene_file <- overlapdf
  integrated_gene_file <-
    integrated_gene_file[1:(nrow(integrated_gene_file) - 1), ]
  integrated_gene_file$foundby <- "both"
  
  # dealing with those that do not have overlapps
  diff_ara_tse <-
    setdiff(seq(1, nrow(aradf), 1), overlapps$ararecord)
  diff_tse_ara <-
    setdiff(seq(1, nrow(tsedf), 1), overlapps$tserecord)
  
  if (length(diff_ara_tse) > 0)
    for (i in 1:length(diff_ara_tse)) {
      tempbind <-
        as.data.frame(c(aradf[diff_ara_tse[i],], rep("notfound", ncol(tsedf)), "ara"))
      names(tempbind) <- names(integrated_gene_file)
      tempbind$tsebegin = -1
      tempbind$tseend = -1
      integrated_gene_file <- rbind(tempbind, integrated_gene_file)
    }
  if (length(diff_tse_ara) > 0)
    for (i in 1:length(diff_tse_ara)) {
      tempbind <-
        as.data.frame(c(rep("notfound", ncol(aradf)), tsedf[diff_tse_ara[i],], "tse"))
      names(tempbind) <- names(integrated_gene_file)
      tempbind$arabegin = -1
      tempbind$araend = -1
      integrated_gene_file <- rbind(tempbind, integrated_gene_file)
    }
  
  # now we make integrated_tse_ara with unified geneid for tse and ara
  m <-
    matrix(nrow = nrow(integrated_gene_file),
           ncol = 27)
  integrated_tse_ara <- as.data.frame(m)
  names(integrated_tse_ara) <- c(
    "geneid",
    "sourceOrg",
    "tseidentity",
    "araidentity",
    "direction",
    "tsebegin",
    "tseend",
    "arabegin",
    "araend",
    "tseac",
    "araac",
    "sourceseq",
    "tsescore",
    "arascore",
    "tsegeneseq",
    "arageneseq",
    "tsegeness",
    "arageness",
    "tseintronbegin",
    "tseintronend",
    "araintronbegin",
    ####
    "araintronend",
    #####
    "tseacloc",
    "sourceSO",
    "tsenote",
    ####
    "aranote",
    #####
    "foundby"
  )
  
  
  # make the geneid from sourceOrg+sourceseq+index within that sourceOrg(which is the rownumber of our dataframe)
  
  for (i in 1:ncol(integrated_gene_file)) {
    integrated_gene_file[, i] <- as.character(integrated_gene_file[, i])
  }
  
  settse <-
    (integrated_gene_file$foundby == "tse" |
       integrated_gene_file$foundby == "both")
  integrated_tse_ara$sourceOrg[settse] = integrated_gene_file$tsesourceOrg[settse]
  integrated_tse_ara$sourceseq[settse] = integrated_gene_file$tsesourceseq[settse]
  integrated_tse_ara$direction[settse] = integrated_gene_file$tsedirection[settse]
  setara <- (integrated_gene_file$foundby == "ara")
  integrated_tse_ara$sourceOrg[setara] =  integrated_gene_file$arasourceOrg[setara]
  integrated_tse_ara$sourceseq[setara] = integrated_gene_file$arasourceseq[setara]
  integrated_tse_ara$direction[setara] = integrated_gene_file$aradirection[setara]
  
  integrated_tse_ara$tseidentity = integrated_gene_file$tseidentity
  integrated_tse_ara$araidentity = integrated_gene_file$araidentity
  integrated_tse_ara$tsebegin = integrated_gene_file$tsebegin
  integrated_tse_ara$tseend = integrated_gene_file$tseend
  integrated_tse_ara$arabegin = integrated_gene_file$arabegin
  integrated_tse_ara$araend = integrated_gene_file$araend
  integrated_tse_ara$tseac = integrated_gene_file$tseac
  integrated_tse_ara$araac = integrated_gene_file$araac
  integrated_tse_ara$tsescore = integrated_gene_file$tsescore
  integrated_tse_ara$arascore = integrated_gene_file$arascore
  integrated_tse_ara$tsegeneseq = integrated_gene_file$tsegeneseq
  integrated_tse_ara$arageneseq = integrated_gene_file$arageneseq
  integrated_tse_ara$tsegeness = integrated_gene_file$tsegeness
  integrated_tse_ara$arageness = integrated_gene_file$arageness
  integrated_tse_ara$tseintronbegin = integrated_gene_file$tseintronbegin
  integrated_tse_ara$tseintronend = integrated_gene_file$tseintronend
  integrated_tse_ara$tseacloc = integrated_gene_file$tseacloc
  integrated_tse_ara$sourceSO = integrated_gene_file$arasourceSO
  integrated_tse_ara$foundby = integrated_gene_file$foundby
  integrated_tse_ara$araintronbegin = integrated_gene_file$araintronbegin
  integrated_tse_ara$araintronend = integrated_gene_file$araintronend
  integrated_tse_ara$tsenote = integrated_gene_file$note
  integrated_tse_ara$aranote = integrated_gene_file$aranote
  
  integrated_tse_ara <- integrated_tse_ara[order(
    integrated_tse_ara$sourceOrg,
    integrated_tse_ara$sourceseq,
    integrated_tse_ara$arabegin
  ), ]
  
  for (i in 1:nrow(integrated_tse_ara)) {
    integrated_tse_ara$geneid[i] = paste(integrated_tse_ara$sourceOrg[i],
                                         integrated_tse_ara$sourceseq[i],
                                         i,
                                         sep = "_")
  }
  
  integrated_tse_ara
  
}
#_____________________ reading aragorn's output  (one file) and return its info as aradf dataframe ____________________________________________

making_ara_df <- function(arafilename) {
  namearray <- unlist(strsplit(arafilename, split = "/"))
  nameorg <- namearray[length(namearray)]
  nameorg <- unlist(strsplit(nameorg, split = "\\."))[1]
  
  aragenename = ""
  arasourceOrg = nameorg
  araidentity = ""
  aradirection = ""
  arabegin = ""
  araend = ""
  araintronbegin = ""
  araintronend = ""
  araac = ""
  arasourceSO = ""
  arasourceseq = ""
  arascore = ""
  arageneseq = ""
  arageness = ""
  aranote = ""
  
  geneinfo <-
    data.frame(
      aragenename ,
      arasourceOrg ,
      araidentity ,
      aradirection ,
      arabegin ,
      araend ,
      araintronbegin,
      araintronend,
      araac ,
      arasourceSO ,
      arasourceseq ,
      arascore ,
      arageneseq ,
      arageness,
      aranote
    )
  
  con = file(arafilename, "r")
  tmrinaline = ""
  checktmrna = ""
  while (TRUE) {
    #for (x in 1:54) {
    aranote = ""
    line1 = readLines(con, n = 1)
    
    #read two lines each loop
    if (length(line1) == 0) {
      break
    }
    line2 = readLines(con, n = 1)
    if (length(line2) == 0) {
      break
    }
    # if line1 starts with > and line2 does not start with 0
    temp1 <- grep("^>+", line1, value = TRUE)
    temp2 <- grep("^0+", line2, value = TRUE)
    if (length(temp1) > 0 & length(temp2) == 0)
    {
      # read the first element of line2 to see how many genes are found for this sequence (count)
      # read the next count * 3 lines
      # extract the sourceOrg, sourceseq, sourceSO
      line1array <-
        unlist(strsplit(line1, split = " | ", fixed = TRUE))
      arasourceseq <- substring(line1array[1], 2)
      #arasourceOrg <- substring(line1array[2], 10)
      arasourceSO <- substring(line1array[5], 4)
      
      genecount <-
        as.integer(unlist(strsplit(
          line2, split = " ", fixed = TRUE
        ))[1])
      for (j in 1:genecount) {
        #making the geneID with sourceorg+sourceseq+genenum
        aragenename <-
          paste(arasourceOrg, arasourceseq, j, sep = "_")
        
        trna_line1 <- readLines(con, n = 1)
        trna_line2 <- readLines(con, n = 1)
        trna_line3 <- readLines(con, n = 1)
        
        arageneseq <- trna_line2
        arageness <- trna_line3
        
        trna_line1_array <-
          unlist(strsplit(trna_line1, split = "\\s+"))
        
        #if we had ")i(" then we have the location for intron
        araactemp <- trna_line1_array[6]
        if (length(grep("+[)]i[(]+", araactemp)) != 0)
        {
          araac <- unlist(strsplit(araactemp, split = "i"))[1]
          intr <-
            unlist(strsplit(unlist(
              strsplit(araactemp, split = "i")
            )[2], split = ","))
          araintronbegin <- substring(intr[1], 2)
          araintronend <-
            as.integer(as.character(araintronbegin)) + as.integer(as.character(substring(intr[2], 1, nchar(intr[2]) -
                                                                                           1)))
        }
        else
        {
          araac <- araactemp
          araintronbegin = "nointron"
          araintronend = "nointron"
        }
        araidentity <- substring(trna_line1_array[2], 6)
        if (substring(araidentity, nchar(araidentity)) == "*")
        {
          aranote <- "pseudo"
          araidentity <-
            substring(trna_line1_array[2], 6, nchar(trna_line1_array[2]) - 1)
        }
        else
          aranote = ""
        coordinate <-
          unlist(str_extract_all(trna_line1_array[3], "\\d+"))
        arabegin <- coordinate[1]
        araend <- coordinate[2]
        firstchar = substring(trna_line1_array[3], 1, 1)
        if (firstchar == "c")
        {
          aradirection <- "-"
        }
        else
        {
          aradirection <- "+"
        }
        
        arascore <- trna_line1_array[4]
        
        # make one record in dataframe
        tempdf <- data.frame(
          aragenename ,
          arasourceOrg ,
          araidentity ,
          aradirection ,
          arabegin ,
          araend ,
          araintronbegin,
          araintronend,
          araac,
          arasourceSO ,
          arasourceseq ,
          arascore ,
          arageneseq ,
          arageness,
          aranote
        )
        geneinfo$araintronend <-
          as.character(geneinfo$araintronend)
        geneinfo <- rbind(geneinfo, tempdf)
        
      }
    }
  }
  close(con)
  geneinfo <- geneinfo[2:nrow(geneinfo),]
  geneinfo
}
#_____________________ reading tRNAscan's output  (one file) and return its info as tsedf dataframe ___________________________________________

making_tse_df <- function(tse_filename, tse_ss_filename) {
  # extracting the sourceOrg from filemame
  namearray <- unlist(strsplit(tse_filename, split = "/"))
  nameorg <- namearray[length(namearray)]
  nameorg <- unlist(strsplit(nameorg, split = "\\."))[1]
  
  tsegenename = ""
  tsesourceOrg =  nameorg
  tseidentity = ""
  tsedirection = ""
  tsebegin = ""
  tseend = ""
  tseac = ""
  tsesourceseq = ""
  tsescore = ""
  tsegeneseq = ""
  tsegeness = ""
  tseintronbegin = ""
  tseintronend = ""
  tseacloc = ""
  note = ""
  
  geneinfo <-
    data.frame(
      tsegenename,
      tsesourceOrg,
      tseidentity,
      tsedirection,
      tsebegin,
      tseend,
      tseac,
      tsesourceseq,
      tsescore,
      tsegeneseq,
      tsegeness,
      tseintronbegin,
      tseintronend,
      tseacloc,
      note
    )
  
  con = file(tse_filename, "r")
  con2 = file(tse_ss_filename, "r") # read 6 lines for each tRNA gene
  
  # # # skip the header
  readLines(con, n = 1)
  readLines(con, n = 1)
  readLines(con, n = 1)
  
  while (TRUE) {
    line = readLines(con, n = 1)
    if (length(line) == 0) {
      break
    }
    
    line_array <- unlist(strsplit(line, split = "\\s+"))
    tsesourceseq = line_array[1]
    tsegenename <-
      paste(tsesourceOrg, tsesourceseq, line_array[2], sep = "_")
    tseidentity = line_array[5]
    if (as.integer(line_array[3]) > as.integer(line_array[4]))
    {
      tsedirection = "-"
    }
    else
    {
      tsedirection = "+"
    }
    tsebegin = line_array[3]
    tseend = line_array[4]
    tseac = line_array[6]
    tsescore = line_array[9]
    tseintronbegin = line_array[7]
    tseintronend = line_array[8]
    if (!is.na(line_array[10]))
      note = line_array[10]
    else
      note = "notfound"
    
    # read the first two lines
    ss_line1 = readLines(con2, n = 1)
    ss_line2 = readLines(con2, n = 1)
    tseacloc = unlist(strsplit(ss_line2, split = "\\s+"))[6]
    l = "."
    while (TRUE) {
      l = readLines(con2, n = 1)
      if (l == "")
        break
      if (unlist(strsplit(l, split = "\\s+"))[1] == "Seq:")
        tsegeneseq <- unlist(strsplit(l, split = "\\s+"))[2]
      if (unlist(strsplit(l, split = "\\s+"))[1] == "Str:")
        tsegeness <- unlist(strsplit(l, split = "\\s+"))[2]
      if(unlist(strsplit(l, split = "\\s+"))[1] == "Possible")
        if( length(grep("truncation",unlist(strsplit(l, split = "\\s+"))[2] ))  !=0 )
        {
          if(note=="notfound")
            note="truncated"
          else
            note <- paste("truncated",note, sep = ",") 
        }
        
    }
    
    tempdf <-
      data.frame(
        tsegenename,
        tsesourceOrg,
        tseidentity,
        tsedirection,
        tsebegin,
        tseend,
        tseac,
        tsesourceseq,
        tsescore,
        tsegeneseq,
        tsegeness,
        tseintronbegin,
        tseintronend,
        tseacloc,
        note
      )
    
    geneinfo <- rbind(geneinfo, tempdf)
  }
  close(con)
  close(con2)
  geneinfo <- geneinfo[2:nrow(geneinfo), ]
  geneinfo
}
#_______________________________________________________________________________________________________________________________________________________
addfuncTypes <- function(integrated_tse_ara) {
  integrated_tse_ara$arafunc <- ''
  integrated_tse_ara$tsefunc <- ''
  for (i in 1:nrow(integrated_tse_ara)) {
    id <- tolower(integrated_tse_ara$tseidentity[i])
    if (id == "ala")
      integrated_tse_ara$tsefunc[i] <- "A"
    if (id == "arg")
      integrated_tse_ara$tsefunc[i] <- "R"
    if (id == "asn")
      integrated_tse_ara$tsefunc[i] <- "N"
    if (id == "asp")
      integrated_tse_ara$tsefunc[i] <- "D"
    if (id == "cys")
      integrated_tse_ara$tsefunc[i] <- "C"
    if (id == "gln")
      integrated_tse_ara$tsefunc[i] <- "Q"
    if (id == "glu")
      integrated_tse_ara$tsefunc[i] <- "E"
    if (id == "gly")
      integrated_tse_ara$tsefunc[i] <- "G"
    if (id == "his")
      integrated_tse_ara$tsefunc[i] <- "H"
    if (id == "ile")
      integrated_tse_ara$tsefunc[i] <- "I"
    if (id == "start")
      integrated_tse_ara$tsefunc[i] <- "start"
    if (id == "leu")
      integrated_tse_ara$tsefunc[i] <- "L"
    if (id == "lys")
      integrated_tse_ara$tsefunc[i] <- "K"
    if (id == "met")
      integrated_tse_ara$tsefunc[i] <- "M"
    if (id == "phe")
      integrated_tse_ara$tsefunc[i] <- "F"
    if (id == "pro")
      integrated_tse_ara$tsefunc[i] <- "P"
    if (id == "ser")
      integrated_tse_ara$tsefunc[i] <- "S"
    if (id == "thr")
      integrated_tse_ara$tsefunc[i] <- "T"
    if (id == "trp")
      integrated_tse_ara$tsefunc[i] <- "W"
    if (id == "tyr")
      integrated_tse_ara$tsefunc[i] <-
        "Y"
    if (id == "val")
      integrated_tse_ara$tsefunc[i] <-
        "V"
    if (id == "sec")
      integrated_tse_ara$tsefunc[i] <-
        "Z"
    if (id == "stop")
      integrated_tse_ara$tsefunc[i] <-
        "#"
    if (id == "sup")
      integrated_tse_ara$tsefunc[i] <-
        "?"
    if (id == "pyl")
      integrated_tse_ara$tsefunc[i] <-
        "O"
  }
  for (i in 1:nrow(integrated_tse_ara)) {
    id <- tolower(integrated_tse_ara$araidentity[i])
    if (id == "ala")
      integrated_tse_ara$arafunc[i] <- "A"
    if (id == "arg")
      integrated_tse_ara$arafunc[i] <- "R"
    if (id == "asn")
      integrated_tse_ara$arafunc[i] <- "N"
    if (id == "asp")
      integrated_tse_ara$arafunc[i] <- "D"
    if (id == "cys")
      integrated_tse_ara$arafunc[i] <- "C"
    if (id == "gln")
      integrated_tse_ara$arafunc[i] <- "Q"
    if (id == "glu")
      integrated_tse_ara$arafunc[i] <- "E"
    if (id == "gly")
      integrated_tse_ara$arafunc[i] <- "G"
    if (id == "his")
      integrated_tse_ara$arafunc[i] <- "H"
    if (id == "ile")
      integrated_tse_ara$arafunc[i] <- "I"
    if (id == "start")
      integrated_tse_ara$arafunc[i] <- "start"
    if (id == "leu")
      integrated_tse_ara$arafunc[i] <- "L"
    if (id == "lys")
      integrated_tse_ara$arafunc[i] <- "K"
    if (id == "met")
      integrated_tse_ara$arafunc[i] <- "M"
    if (id == "phe")
      integrated_tse_ara$arafunc[i] <- "F"
    if (id == "pro")
      integrated_tse_ara$arafunc[i] <- "P"
    if (id == "ser")
      integrated_tse_ara$arafunc[i] <- "S"
    if (id == "thr")
      integrated_tse_ara$arafunc[i] <- "T"
    if (id == "trp")
      integrated_tse_ara$arafunc[i] <- "W"
    if (id == "tyr")
      integrated_tse_ara$arafunc[i] <-
        "Y"
    if (id == "val")
      integrated_tse_ara$arafunc[i] <-
        "V"
    if (id == "sec")
      integrated_tse_ara$arafunc[i] <-
        "Z"
    if (id == "stop")
      integrated_tse_ara$arafunc[i] <-
        "#"
    if (id == "sup")
      integrated_tse_ara$arafunc[i] <-
        "?"
    if (id == "pyl")
      integrated_tse_ara$arafunc[i] <-
        "O"
  }
  integrated_tse_ara
    
}
#______________________________________________________________________________________________________________________________________________________
initiatorDetecting <- function(geneDF) {
  library(cluster)    # clustering algorithms
  library(factoextra) # clustering visualization
  library(dendextend) # for comparing two dendrograms
  # GeneDF from the function profile_cm
  # GeneDF$headers <- gsub(">","",GeneDF$headers)
  
  GeneDF$araac <- genefile[genefile$geneid %in% GeneDF$geneid,]$araac
  GeneDF$tseac <- genefile[genefile$geneid %in% GeneDF$geneid,]$tseac
  GeneDF$foundby <- genefile[genefile$geneid %in% GeneDF$geneid,]$foundby
  GeneDF$tseidentity <- genefile[genefile$geneid %in% GeneDF$geneid,]$tseidentity
  geneDF <- GeneDF
  iscat <- geneDF$araac == "(cat)" | geneDF$tseac == "CAT"
  catDF <- geneDF[iscat, ]
  #table(catDF$foundby)
  #ara both 
  #1  184 
  
  geneseqs <- character(length = nrow(catDF))
  for (i in 1:nrow(catDF)) {
    if (catDF$foundby[i] == "both")
      geneseqs[i] <- catDF$sequences[i]
    else if (catDF$foundby[i] == "tse")
      geneseqs[i] <- catDF$sequences[i]
    else if (catDF$foundby[i] == "ara")
      geneseqs[i] <- catDF$sequences[i]
  }
  m <- matrix(nrow = nrow(catDF), ncol = nrow(catDF))
  distanceDF <- as.data.frame(m)
  for (i in 1:length(geneseqs)) {
    for (j in 1:length(geneseqs)) {
      distanceDF[i, j] <- adist(geneseqs[i], geneseqs[j])
    }
  }
  
  # m <- matrix(nrow = nrow(clus1DF), ncol = nrow(clus1DF))
  # distanceDF <- as.data.frame(m)
  # for (i in 1:length(clus1DF$sequences)) {
  #   for (j in 1:length(clus1DF$sequences)) {
  #     distanceDF[i, j] <- adist(clus1DF$sequences[i], clus1DF$sequences[j])
  #   }
  # }
  # m <- matrix(nrow = nrow(clus2DF), ncol = nrow(clus2DF))
  # distanceDF <- as.data.frame(m)
  # for (i in 1:length(clus2DF$sequences)) {
  #   for (j in 1:length(clus2DF$sequences)) {
  #     distanceDF[i, j] <- adist(clus2DF$sequences[i], clus2DF$sequences[j])
  #   }
  # }
  ####################################################### Hierarchical Cluster Analysis ######################################
  # Hierarchical clustering is an alternative approach to k-means clustering for identifying groups in the dataset.
  # we used Agglomerative clustering which works in a bottom-up manner
  # we measure the dissimilarity between two clusters of observations with Ward.D2 which It minimizes the total within-cluster variance.
  # At each step the pair of clusters with minimum between-cluster distance are merged
  
  # libraries
  #library(tidyverse)  # data manipulation
  library(cluster)    # clustering algorithms
  library(factoextra) # clustering visualization
  library(dendextend) # for comparing two dendrograms
  
  # The height of the fusion, provided on the vertical axis, indicates the (dis)similarity between two observations.
  # Note that, conclusions about the proximity of two observations can be drawn only based on the height where branches
  # containing those two observations first are fused. We cannot use the proximity of two observations along the horizontal axis as a criteria of their similarity
  hc <- hclust(as.dist(distanceDF), method = 'ward.D2')
  plot(hc, cex = 0.6, hang = -1)
  
  # The height of the cut to the dendrogram controls the number of clusters obtained. It plays the same role as the k in k-means clustering.
  #
  rect.hclust(hc, k = 3, border = 2:5)
  sub_groups <- cutree(hc, k = 3)
  fviz_cluster(list(data = distanceDF, cluster = sub_groups))
  
  
  # extracting minimum and maximum
  subg1 <- sub_groups == 1
  subg2 <- sub_groups == 2
  subg3 <- sub_groups == 3
  
  cluster1 <-
    paste(min(distanceDF[subg1, subg1]), max(distanceDF[subg1, subg1]), sep = "-")
  cluster2 <-
    paste(min(distanceDF[subg2, subg2]), max(distanceDF[subg2, subg2]), sep = "-")
  cluster3 <-
    paste(min(distanceDF[subg3, subg3]), max(distanceDF[subg3, subg3]), sep = "-")
  
  clusterDF <- data.frame(table(sub_groups))
  clusterDF$editdistanceRange <- c(cluster1, cluster2, cluster3)
  
  clus1DF <- catDF[subg1, ]
  clus2DF <- catDF[subg2, ]
  clus3DF <- catDF[subg3, ]
  
  ## making a table with columns: Clusters, 1-72(A-T), 29-41(GC), 30-40(GC/T),31-39(GF), 32,33,34,35,36,37,38
  ## green green indicates bases conserved in all initiators throughout the three domains
  
  m <- matrix(nrow = 3,
              ncol = 9,
              data = 0)
  summerytable <- as.data.frame(m)
  names(summerytable) <-
    c(
      "Clusters",
      "#tRNAs",
      "11–24(C-G)",
      "54-60(A-A)(T-T)",
      "1-72(A-T)",
      "29-31(GGG)",
      "39-41(CCC/CCT)",
      "#posisInDloop",
      "20A"
    )
  
  summerytable[, 1] <- c("Cluster1", "Cluster2", "Cluster3")
  summerytable[, 2] <- c(nrow(clus1DF), nrow(clus2DF), nrow(clus3DF))
  GC1 = 0
  #>>>>>>>..>>>>........<<<<.>>>>>.......<<<<<....>>>>>.......<<<<<<<<<<<<.
  #
  for (i in 1:nrow(clus1DF)) {
    if ((substring(clus1DF$sequences[i], 11, 11) == "C") &
        (substring(clus1DF$sequences[i], 24, 24) == "G"))
      GC1 <- GC1 + 1
  }
  GC2 = 0
  for (i in 1:nrow(clus2DF)) {
    if ((substring(clus2DF$sequences[i], 11, 11) == "C") &
        (substring(clus2DF$sequences[i], 24, 24) == "G"))
      GC2 <- GC2 + 1
  }
  GC3 = 0
  for (i in 1:nrow(clus3DF)) {
    if ((substring(clus3DF$sequences[i], 11, 11) == "C") &
        (substring(clus3DF$sequences[i], 24, 24) == "G"))
      GC3 <- GC3 + 1
  }
  
  summerytable[, 3] <- c(GC1, GC2, GC3)
  
  AA1 = 0
  for (i in 1:nrow(clus1DF)) {
    if (((substring(clus1DF$sequences[i], 53, 53) == "A") &
         (substring(clus1DF$sequences[i], 59, 59) == "A")) |
        ((substring(clus1DF$sequences[i], 53, 53) == "T") &
         (substring(clus1DF$sequences[i], 59, 59) == "T")))
      AA1 <- AA1 + 1
  }
  AA2 = 0
  for (i in 1:nrow(clus2DF)) {
    if (((substring(clus2DF$sequences[i], 53, 53) == "A") &
         (substring(clus2DF$sequences[i], 59, 59) == "A")) |
        (substring(clus2DF$sequences[i], 53, 53) == "T") &
        (substring(clus2DF$sequences[i], 59, 59) == "T"))
      AA2 <- AA2 + 1
  }
  AA3 = 0
  for (i in 1:nrow(clus3DF)) {
    if (((substring(clus3DF$sequences[i], 53, 53) == "A") &
         (substring(clus3DF$sequences[i], 59, 59) == "A")) |
        (substring(clus3DF$sequences[i], 53, 53) == "T") &
        (substring(clus3DF$sequences[i], 59, 59) == "T"))
      AA3 <- AA3 + 1
  }
  
  summerytable[, 4] <- c(AA1, AA2, AA3)
  
  AT1 = 0
  for (i in 1:nrow(clus1DF)) {
    if ((substring(clus1DF$sequences[i], 1, 1) == "A") &
        (substring(clus1DF$sequences[i], 71, 71) == "T"))
      AT1 <- AT1 + 1
  }
  AT2 = 0
  for (i in 1:nrow(clus2DF)) {
    if ((substring(clus2DF$sequences[i], 1, 1) == "A") &
        (substring(clus2DF$sequences[i], 71, 71) == "T"))
      AT2 <- AT2 + 1
  }
  AT3 = 0
  for (i in 1:nrow(clus3DF)) {
    if ((substring(clus3DF$sequences[i], 1, 1) == "A") &
        (substring(clus3DF$sequences[i], 71, 71) == "T"))
      AT3 <- AT3 + 1
  }
  summerytable[, 5] <- c(AT1, AT2, AT3)
  
  GGG1 = 0
  for (i in 1:nrow(clus1DF)) {
    if ((substring(clus1DF$sequences[i], 29, 31) == "GGG"))
      GGG1 <- GGG1 + 1
  }
  GGG2 = 0
  for (i in 1:nrow(clus2DF)) {
    if ((substring(clus2DF$sequences[i], 29, 31) == "GGG"))
      GGG2 <- GGG2 + 1
  }
  GGG3 = 0
  for (i in 1:nrow(clus3DF)) {
    if ((substring(clus3DF$sequences[i], 29, 31) == "GGG"))
      GGG3 <- GGG3 + 1
  }
  summerytable[, 6] <- c(GGG1, GGG2, GGG3)
  
  CCC1 = 0
  for (i in 1:nrow(clus1DF)) {
    if ((substring(clus1DF$sequences[i], 39, 41) == "CCC") |
        (substring(clus1DF$sequences[i], 39, 41) == "CCT"))
      CCC1 <- CCC1 + 1
  }
  CCC2 = 0
  for (i in 1:nrow(clus2DF)) {
    if ((substring(clus2DF$sequences[i], 39, 41) == "CCC") |
        (substring(clus2DF$sequences[i], 39, 41) == "CCT"))
      CCC2 <- CCC2 + 1
  }
  CCC3 = 0
  for (i in 1:nrow(clus3DF)) {
    if ((substring(clus3DF$sequences[i], 39, 41) == "CCC") |
        (substring(clus3DF$sequences[i], 39, 41) == "CCT"))
      CCC3 <- CCC3 + 1
  }
  summerytable[, 7] <- c(CCC1, CCC2, CCC3)
  
  summerytable[, 8] <- c(7, 8, 1)
  
  
  A120 = 0
  for (i in 1:nrow(clus1DF)) {
    if ((substring(clus1DF$sequences[i], 20, 20) == "A"))
      A120 <- A120 + 1
  }
  
  A220 = 0
  for (i in 1:nrow(clus2DF)) {
    if ((substring(clus2DF$sequences[i], 20, 20) == "A"))
      A220 <- A220 + 1
  }
  A320 = 0
  for (i in 1:nrow(clus3DF)) {
    if ((substring(clus3DF$sequences[i], 20, 20) == "A"))
      A320 <- A320 + 1
  }
  summerytable[, 9] <- c(A120, A220, A320)
  
  #summerytable[3, ] <- c("Cluster3", 2, 2, 2, 0, 0, 0, "8/9", 0)
  for (i in 1:nrow(clus1DF)) {
    geneDF[geneDF$geneid == clus1DF$geneid[i], ]$arafunc <- "X"
    geneDF[geneDF$geneid == clus1DF$geneid[i], ]$tsefunc <- "X"
  }
  #summerytable[3, ] <- c("Cluster3", 2, 2, 2, 0, 0, 0, "8/9", 0)
  
  # add one column to see the min and max distance of sequences in one cluster
  summerytable$distanceRange <- ""
  m <- matrix(nrow = nrow(clus1DF), ncol = nrow(clus1DF))
  distanceDF <- as.data.frame(m)
  for (i in 1:length(clus1DF$sequences)) {
    for (j in 1:length(clus1DF$sequences)) {
      distanceDF[i, j] <- adist(clus1DF$sequences[i], clus1DF$sequences[j])
    }
  }
  summerytable$distanceRange[1] <- paste(min(distanceDF),max(distanceDF),sep = "-")
  m <- matrix(nrow = nrow(clus2DF), ncol = nrow(clus2DF))
  distanceDF <- as.data.frame(m)
  for (i in 1:length(clus2DF$sequences)) {
    for (j in 1:length(clus2DF$sequences)) {
      distanceDF[i, j] <- adist(clus2DF$sequences[i], clus2DF$sequences[j])
    }
  }
  summerytable$distanceRange[2] <- paste(min(distanceDF),max(distanceDF),sep = "-")
  m <- matrix(nrow = nrow(clus1DF), ncol = nrow(clus1DF))
  distanceDF <- as.data.frame(m)
  for (i in 1:length(clus1DF$sequences)) {
    for (j in 1:length(clus1DF$sequences)) {
      distanceDF[i, j] <- adist(clus1DF$sequences[i], clus1DF$sequences[j])
    }
  }
  m <- matrix(nrow = nrow(clus3DF), ncol = nrow(clus3DF))
  distanceDF <- as.data.frame(m)
  for (i in 1:length(clus3DF$sequences)) {
    for (j in 1:length(clus3DF$sequences)) {
      distanceDF[i, j] <- adist(clus3DF$sequences[i], clus3DF$sequences[j])
    }
  }
  summerytable$distanceRange[3] <- paste(min(distanceDF),max(distanceDF),sep = "-")
  geneDF
}
#______________________________________________________________________________________________________________________________________________________
fix.Mismatch.Identity(geneDF) {
  isboth <- geneDF$foundby == "both"
  intersectDF <- geneDF[isboth,]
  mismatchDF <- intersectDF[intersectDF$tsefunc!=intersectDF$arafunc,]
  table(paste(mismatchDF$arafunc,mismatchDF$tsefunc,sep = "|"))
  # D|I L|? L|E L|M N|Y O|M S|R W|G 
  # 5   3   1   9  11   2   1   3 
  D_df <- mismatchDF[mismatchDF$arafunc=="D",]
  araD_df <- geneDF[geneDF$arafunc=="D",]
  tseD_df <- geneDF[geneDF$tsefunc=="D",]
}
