TriTrypAlignment <- function() {
  library(gdata)
  library(readr)
  library(Hmisc)
  library(stringr)
  library(seqinr)
  dirpath <-
    "/home/fatemeh/TriTrypGenome_project/GeneAnnotation/"
  genefilename <-"tsfm_input_geneset.txt"
  geneDF <- read.input(dirpath, genefilename)
  geneDF_Int_noVar <- remove.Vararm(geneDF)
  geneDF_Int_noVar_noIntron <- remove.Intron(geneDF_Int_noVar)
  # read homo gene file, change the headers to match with TryTrip genefile
  HomoDF <- read.homogenes(dirpath,"hg38-tRNAs.fa")
  
  # bind TriTryp genefile with homo genefile and write the result in to a fasta file to be processed by covea
  write.genefile(geneDF_Int_noVar_noIntron,HomoDF,
                 dirpath,
                 "coveainput.fasta")
  system(
    "covea /home/fatemeh/TriTrypGenome_project/GeneAnnotation/TRNA2-euk.cm /home/fatemeh/TriTrypGenome_project/GeneAnnotation/coveainput.fasta > /home/fatemeh/TriTrypGenome_project/GeneAnnotation/Aligned_TriTryp_Homo.covea"
  )
  
  ################################## Covea Editing   ################################## 
  
  covea_filename <- "Aligned_TriTryp_Homo.covea"
  readSeqsIntoDf(dirpath, covea_filename)
  seqDB <-
    read_csv(paste(dirpath, "/trash/coveaDF.txt", sep = ""),
             col_types = cols(.default = col_character()))
  SSDB <-
    read_csv(paste(dirpath, "/trash/coveaDF_SS.txt", sep = ""))
  
  fasta_filename <- "Aligned_TriTryp_Homo.fasta"
  CS_filename <- "Aligned_TriTryp_Homo_structfile.txt"
  
  editAlignment(seqDB, SSDB, dirpath, fasta_filename, CS_filename)
  ################################## Covea Editing   ################################## 
  }

write.genefile <- function(geneDF_Int_noVar_noIntron,HomoDF,resultpath,filename) {
  # this function will write the genefile with removed variable arm and removed intron into a fasta file format called Integrated_Genes_NoVarIntron.fasta
  geneDF_Int_noVar_noIntron$GeneSeq <-
    geneDF_Int_noVar_noIntron$GeneSeqTse
  geneDF_Int_noVar_noIntron[geneDF_Int_noVar_noIntron$foundby == "ara", ]$GeneSeq <-
    geneDF_Int_noVar_noIntron[geneDF_Int_noVar_noIntron$foundby == "ara", ]$GeneSeqAra
  TriTrypDF <- data.frame(
    paste(">", geneDF_Int_noVar_noIntron$GeneID, sep = ""),
    geneDF_Int_noVar_noIntron$GeneSeq
  )
  names(TriTrypDF) <- c("header", "seq")
  names(HomoDF) <- c("header", "seq")
  coveainputDF <- rbind(TriTrypDF, HomoDF)
  write.fwf(
    coveainputDF
    ,
    file = paste(resultpath, filename, sep = ""),
    sep = "\n",
    colnames = FALSE
  )
}
read.input <- function(dirpath, genefilename) {
  # this function will read the prepared gene file as input for from tsfm_input_geneset.txt
  # For genes found by ony ara, it will add _ara<genemodel> at the end of gene id (36 genes)
  # keep all the information needed to make a fasta file for doing the alignment
  
  filepath <- paste(dirpath,genefilename,sep = "")
  genefile <-
    read.table(filepath, header = TRUE, colClasses = "character")
  araonly <- genefile$foundby == "ara"
  genefile$geneid[araonly] <-
    paste(genefile$geneid[araonly], "_ara", genefile$genefunc[araonly], sep = "")
  genefile$geneid[!araonly] <-
    paste(genefile$geneid[!araonly], "_", genefile$genefunc[!araonly], sep = "")
  
  geneseq <- character(length = length(genefile$geneid))
  geness <- character(length = length(genefile$geneid))
  geneDF_fasta <-
    data.frame(
      genefile$geneid,
      genefile$tsegeneseq,
      genefile$arageneseq,
      genefile$tsegeness,
      genefile$arageness,
      genefile$foundby
    )
  names(geneDF_fasta) <-
    c("GeneID",
      "GeneSeqTse",
      "GeneSeqAra",
      "SSTSE",
      "SSARA",
      "foundby")
  
  genefile$arageneseq <- as.character(genefile$arageneseq)
  genefile$tsegeneseq <- as.character(genefile$tsegeneseq)
  geneDF_fasta$GeneID <- as.character(geneDF_fasta$GeneID)
  geneDF_fasta$GeneSeqTse <- as.character(geneDF_fasta$GeneSeqTse)
  geneDF_fasta$GeneSeqAra <- as.character(geneDF_fasta$GeneSeqAra)
  geneDF_fasta$SSARA <- as.character(geneDF_fasta$SSARA)
  geneDF_fasta$SSTSE <- as.character(geneDF_fasta$SSTSE)
  
  geneDF_fasta$GeneSeq <- ""
  geneDF_fasta$GeneSeq <- geneDF_fasta$GeneSeqTse
  geneDF_fasta[geneDF_fasta$foundby=="ara",]$GeneSeq <- toupper(geneDF_fasta[geneDF_fasta$foundby=="ara",]$GeneSeqAra)
  # write.fwf(
  #   data.frame(
  #     paste(">", geneDF_fasta$GeneID, sep = ""),
  #     geneDF_fasta$GeneSeq
  #   ),
  #   file = paste("/home/fatemeh/Leishmania_2019/LeishLatex/", "final_geneset.fasta", sep = ""),
  #   sep = "\n",
  #   colnames = FALSE
  # )
 
  geneDF_fasta
}
read.homogenes <- function(dirpath, genefilename){
  filepath <- paste(dirpath, genefilename, sep = "")
  fastafile <-
    read.fasta(filepath,seqtype="DNA",as.string = TRUE)
  headers <- names(fastafile)
  sequences <- character(length = length(headers))
  funcClasses <- character(length = length(headers))
  for (i in 1:length(headers)) {
    seq <- fastafile[i]
    sequences[i] <- toupper(paste(paste(seq,collapse = ""),collapse = ""))
  }
  identities <- substr(headers,19,21)
  for (i in 1:length(identities)) {
    id <- identities[i]
    id <- tolower(id)
    if (id == "ala")
      funcClasses[i] <- "A"
    if (id == "arg")
      funcClasses[i] <- "R"
    if (id == "asn")
      funcClasses[i] <- "N"
    if (id == "asp")
      funcClasses[i] <- "D"
    if (id == "cys")
      funcClasses[i] <- "C"
    if (id == "gln")
      funcClasses[i] <- "Q"
    if (id == "glu")
      funcClasses[i] <- "E"
    if (id == "gly")
      funcClasses[i] <- "G"
    if (id == "his")
      funcClasses[i] <- "H"
    if (id == "ile")
      funcClasses[i] <- "I"
    if (id == "start")
      funcClasses[i] <- "start"
    if (id == "leu")
      funcClasses[i] <- "L"
    if (id == "lys")
      funcClasses[i] <- "K"
    if (id == "met")
      funcClasses[i] <- "M"
    if (id == "phe")
      funcClasses[i] <- "F"
    if (id == "pro")
      funcClasses[i] <- "P"
    if (id == "ser")
      funcClasses[i] <- "S"
    if (id == "thr")
      funcClasses[i] <- "T"
    if (id == "trp")
      funcClasses[i] <- "W"
    if (id == "tyr")
      funcClasses[i] <-
      "Y"
    if (id == "val")
      funcClasses[i] <-
      "V"
    if (id == "sec")
      funcClasses[i] <-
      "Z"
    if (id == "stop")
      funcClasses[i] <-
      "#"
    if (id == "sup")
      funcClasses[i] <-
      "?"
    if (id == "pyl")
      funcClasses[i] <-
      "O"
    if (id == "ime")
      funcClasses[i] <-
      "X"
  }
  headers2 <- paste(">",headers,"_",funcClasses,sep = "")
  
  geneDF_fasta <-
    data.frame(
      headers2,
      sequences
    )
  names(geneDF_fasta) <-
    c("header",
      "geneseq")
  
  geneDF_fasta$header <- as.character(geneDF_fasta$header)
  geneDF_fasta$geneseq <- as.character(geneDF_fasta$geneseq)
  geneDF_fasta
}
remove.Vararm <- function(geneDF) {
  # This function will remove the variable arms from both secondary structure and sequence
  # For genes found by both genefinders or tseonly: we use the TSE sequences as reference which is always reported up to base 73.
  # We will count number of arms by counting number of occurances of "><" in function countarms. if we had 4 arms, we remove tha third arm starting ">" and ending with "<"
  
  for (i in 1:nrow(geneDF)) {
    # count number of "><"
    
    geneSS <- geneDF[i, ]$SSTSE
    if (geneDF$foundby[i] == "ara")
      geneSS <- geneDF[i, ]$SSARA
    numarms <- countarms(geneSS)
    if (numarms == 4)
    {
      geneseq <- geneDF[i, ]$GeneSeqTse
      geness <- geneDF[i, ]$SSTSE
      
      if (geneDF$foundby[i] == "ara")
      {
        geneseq <- geneDF[i, ]$GeneSeqAra
        geness <- geneDF[i, ]$SSARA
      }
      varcor <-
        removearm(geneseq, geness, as.character(geneDF$foundby[i]))
      firstchunk <-
        substring(geneseq, 1, varcor[1] - 1)
      SS_firstchunk <-
        substring(geness, 1, varcor[1] - 1)
      
      secondchunk <-
        substring(geneseq,
                  varcor[2] + 1,
                  nchar(geneseq))
      SS_secondchunk <-
        substring(geness,
                  varcor[2] + 1,
                  nchar(geness))
      if (geneDF$foundby[i] == "both")
      {
        geneDF$GeneSeqTse[i] <-
          paste(firstchunk, secondchunk, sep = "")
        geneDF$SSTSE[i] <-
          paste(SS_firstchunk, SS_secondchunk, sep = "")
      }
      if (geneDF$foundby[i] == "ara")
      {
        geneDF$GeneSeqAra[i] <-
          paste(firstchunk, secondchunk, sep = "")
        geneDF$SSARA[i] <-
          paste(SS_firstchunk, SS_secondchunk, sep = "")
      }
    }
    
    if (numarms < 3)
    {
      print("We have genes with unusual structure which need to be removed!")
    }
  }
  geneDF
}
remove.Intron <-  function(geneDF_Int_noVar) {
  # This function will remove the introns from both secondary structure and sequence
  # For genes found by both genefinders: we use the TSE sequences as reference which is always reported up to base 73.
  # Based on TSE manuscript, in genes found by tse nucleotides matching the "consensus" tRNA model used in Cove analysis appear in upper case,
  # while introns and other nucleotides in non-conserved positions are printed in lower-case letters. So, we will remove all the lower case letters
  # the result will be saved in a dataframe with only two columns "geneid" and "geneseq"
  # for genes found by aragorn, remove sites that are marked as i in their secondary structure, and reports genes up to position 73
  for (i in 1:nrow(geneDF_Int_noVar)) {
    if (geneDF_Int_noVar$foundby[i] == "both") {
      geneDF_Int_noVar$GeneSeqTse[i] <-
        gsub("[[:lower:]]", "", geneDF_Int_noVar$GeneSeqTse[i])
    }
    if (geneDF_Int_noVar$foundby[i] == "ara") {
      geness <- geneDF_Int_noVar$SSARA[i]
      geness_arr <-
        substring(geness, seq(1, nchar(geness), 1), seq(1, nchar(geness), 1))
      gene <- geneDF_Int_noVar$GeneSeqAra[i]
      p73 <- str_sub(gene, nchar(geness) + 1, nchar(geness) + 1)
      gene <- str_sub(gene, 1, nchar(geness))
      gene_arr <-
        substring(gene, seq(1, nchar(gene), 1), seq(1, nchar(gene), 1))
      isintron <- geness_arr == "i"
      geneDF_Int_noVar$GeneSeqAra[i] <-
        paste(paste(gene_arr[!isintron], collapse = ""), p73, sep = "")
      geneDF_Int_noVar$SSARA[i] <-
        paste(geness_arr[!isintron], collapse = "")
    }
  }
  geneDF_Int_noVar
}
countarms <- function(geneSS) {
  geneSSarr <-
    substring(geneSS, seq(1, nchar(geneSS), 1), seq(1, nchar(geneSS), 1))
  counter = 0
  flag = 0
  vararmend = 0
  for (i in 1:length(geneSSarr)) {
    if ((geneSSarr[i] == "<" | geneSSarr[i] == ")") & flag == 0)
    {
      counter <- counter + 1
      flag = 1
      
    }
    if (geneSSarr[i] == ">" | geneSSarr[i] == "(")
      flag = 0
  }
  
  counter
}
removearm <- function(geneSeq, geneSS,foundby) {
  flag = 0
  counter = 0
  forward = 1
  varbeg = 0
  varend = 0
  armS=">"
  armE="<"
  if(foundby == "ara")
  {
    armS="("
    armE=")"
  }
  geneSSarr <-
    substring(geneSS, seq(1, nchar(geneSS), 1), seq(1, nchar(geneSS), 1))
  for (i in 1:length(geneSSarr)) {
    if ((geneSSarr[i] == armS) & flag == 0)
    {
      counter <- counter + 1
      flag = 1
      if (counter == 3)
      {
        print(i)
        # this is the begining of the variable arm
        # count howmany > you have
        # go forward untill you see that many <
        forward = 0
        varbeg = i
        while (geneSSarr[i] != armE)
        {
          if (geneSSarr[i] == armS)
            forward = forward + 1
          i = i + 1
        }
        while (forward != 0) {
          if (geneSSarr[i] == armE)
            forward = forward - 1
          
          i = i + 1
        }
        varend = i - 1
      }
    }
    if (geneSSarr[i] == armE)
      flag = 0
  }
  # print(varend)
  # print(varbeg)
  varcor <- c(varbeg, varend)
  varcor
}
readSeqsIntoDf <- function(dirpath,covea_filename) {
  # this function will read the sequences from covea file along with their secondary structure into a dataframe
  # will write the result in coveaDF.txt and coveaDF_SS.txt
  coveafilepath <- paste(dirpath, covea_filename, sep = "")
  CS <-
    grep("#=CS +",
         readLines(coveafilepath),
         value = TRUE)
  CStemp <-
    unlist(strsplit(CS, split = "\\s+"))
  CS <- CStemp[seq(2, length(CStemp), 2)]
  CS <- paste(CS, collapse = '')
  CSarr <- substring(CS, seq(1, nchar(CS), 1), seq(1, nchar(CS), 1))
  
  # read the name of the sequences into variable "seqnames"
  SQs <- grep("#=SQ +",
              readLines(coveafilepath),
              value = TRUE)
  seqnames <- character(length = length(SQs))
  for (i in 1:length(SQs)) {
    seqnames[i] <- unlist(strsplit(SQs[i], split = " "))[2]
  }
  
  # define the main data frame as "seqDB" and assign names to it
  m <- matrix(ncol = (length(seqnames) + 1), nrow = length(CSarr))
  seqDB  <- as.data.frame(m)
  names(seqDB) <- c(seqnames, "CS")
  
  for (i in 1:length(seqnames)) {
    pat <- paste("^", seqnames[i], " +", sep = "")
    myseq <-
      grep(pattern = pat,
           readLines(coveafilepath),
           value = TRUE)
    temp <-
      unlist(strsplit(myseq, split = "\\s+"))
    myseq <- temp[seq(2, length(temp), 2)]
    myseq <- paste(myseq, collapse = '')
    myseqarr <-
      substring(myseq, seq(1, nchar(myseq), 1), seq(1, nchar(myseq), 1))
    seqDB[, i] <- myseqarr
  }
  seqDB$CS <- CSarr
  
  #_____________________ reading the #=SS lines into a nother data frame as SSDB ____________
  # define the main data frame as "seqDB" and assign names to it
  
  SSDB = seqDB
  
  SSs <- grep("#=SS +",
              readLines(coveafilepath),
              value = TRUE)
  
  # number of sections is CStemp/2 = 3
  numsec <- length(CStemp) / 2
  end <- length(SSs) / numsec
  for (i in 1:end) {
    if (numsec == 3)
      myseq <-
        paste(unlist(strsplit(SSs[i], split = "\\s+"))[2],
              unlist(strsplit(SSs[end + i], split = "\\s+"))[2],
              unlist(strsplit(SSs[(2 * end) + i], split = "\\s+"))[2],
              sep = "")
    if (numsec == 2)
      myseq <-
        paste(unlist(strsplit(SSs[i], split = "\\s+"))[2], unlist(strsplit(SSs[end +
                                                                                 i], split = "\\s+"))[2], sep = "")
    
    myseqarr <-
      substring(myseq, seq(1, nchar(myseq), 1), seq(1, nchar(myseq), 1))
    SSDB[, i] <- myseqarr
  }
  
  write_csv(seqDB, path  = paste(dirpath, "/trash/coveaDF.txt", sep = ""))
  write_csv(SSDB, path  = paste(dirpath, "/trash/coveaDF_SS.txt", sep = ""))
  
}
writeCovea <- function(SSDB, seqDB, dirpath,cs,fasta_filename,CS_filename) {
  # this function will translates . to - and writes the output as 
  coveaseqs <- character(length = ncol(seqDB) - 1)
  coveass <- character(length = ncol(seqDB) - 1)
  
  seqDB <- seqDB[, -ncol(seqDB)]
  SSDB <- SSDB[, -ncol(SSDB)]
  for (i in 1:(length(coveaseqs))) {
    coveaseqs[i] <-
      as.character(paste((as.data.frame(seqDB[, i])[, 1]), collapse = ''))
    coveass[i] <-
      as.character(paste((as.data.frame(SSDB[, i])[, 1]), collapse = ''))
    
  }
  
  mynames <- names(seqDB)[!names(seqDB) %in% "CS"]
  
  
  # remove "."s from sequences and write them in a fasta file to run with covea
  for (i in 1:length(coveaseqs)) {
    coveaseqs[i] <- gsub("[.]", "-",coveaseqs[i])
  }
  
  write.fwf(
    data.frame(paste(">", mynames, sep = ""),
               coveaseqs),
    file = paste(dirpath, fasta_filename, sep = ""),
    sep = "\n",
    colnames = FALSE
  )
  
  #cs <- paste(seqDB$CS, collapse = '')
  write.fwf(data.frame(cs),
            file = paste(dirpath, CS_filename, sep = ""),
            colnames = FALSE)
  
}
editAlignment <- function(seqDB, SSDB, dirpath,fasta_filename,CS_filename) {
  # this function will remove:
  # 1. sites that have more than 99% gap (removing rows that have more than 99% ".") 
  # 2. sequences with more than 8 gaps
  # at the end it will call the function writeCovea to write 
  delpos = " "
  for (i in 1:nrow(seqDB)) {
    temp <-
      data.frame(table(
        seqDB[i,] == "." |
          seqDB[i,] == "t" |
          seqDB[i,] == "g" | seqDB[i,] == "c" | seqDB[i,] == "a"
      ))
    if (length(temp[temp$Var1 == "FALSE", 2]) != 0)
    {
      if (temp[temp$Var1 == "FALSE", 2] != ncol(seqDB))
      {
        gapperc <-
          (temp[temp$Var1 == "TRUE", 2] / (temp[temp$Var1 == "TRUE", 2] + temp[temp$Var1 ==
                                                                                 "FALSE", 2])) * 100
        if (gapperc > 98)#####################
        {
          delpos = c(delpos, i)
        }
        
      }
    }
    else{
      if (temp[temp$Var1 == "TRUE", 2] == ncol(seqDB))
        delpos = c(delpos, i)
    }
  }
  delpos <- delpos[2:length(delpos)]
  delpos <- unique(delpos)
  delpos = as.integer(delpos)
  seqDB <- seqDB[-delpos,]
  SSDB <- SSDB[-delpos,]
  
  cs <- paste(seqDB$CS, collapse = '')
  
  # remove sequences with more than 8 gaps!
  delpos = " "
  flag = 0
  for (i in 1:ncol(seqDB)) {
    temp <- data.frame(table(seqDB[, i] == "."))
    tempseq <- paste((as.data.frame(seqDB[, i])[, 1]), collapse = '')
    if (length(grep("\\.\\.\\.", tempseq)) == 1)
    {
      if (names(seqDB[, i]) != "CS")
        delpos = c(delpos, i)
    }
    if (length(temp[temp$Var1 == "TRUE", 2]) != 0)
      if (temp[temp$Var1 == "TRUE", 2] > 3)
      {
        if (names(seqDB[, i]) != "CS")
          delpos = c(delpos, i)
      }
    for (x in 1:length(unlist(seqDB[, i]))) {
      if (unlist(seqDB[, i])[x] == "n" | unlist(seqDB[, i])[x] == "N")
        flag = 1
    }
    if (flag == 1)
    {
      print("shit")
      delpos = c(delpos, i)
      flag = 0
    }
  }
  if (length(delpos) > 1)
  {
    delpos <- delpos[2:length(delpos)]
    delpos <- unique(delpos)
    delpos = as.integer(delpos)
    seqDB <- seqDB[,-delpos]
    SSDB <- SSDB[,-delpos]
  }
  
  # remove sites that are lowercase or . in 99% of sequences
  
  # remove sequences with gap in their anticodon
  
  writeCovea(SSDB, seqDB, dirpath,cs,fasta_filename,CS_filename)
  
}
map2sprinzle <- function(dirpath,fasta_filename,CS_filename){
  # this function will read the final aigned sequences with their secondary structure
  # assignes positions to each base according to sprinzle
  seqdf <-
    read.table(paste(dirpath, fasta_filename, sep = ""))
  seq_arr <- as.character(seqdf$V1)
  headers <- seq_arr[seq(1, length(seq_arr), 2)]
  sequences <- seq_arr[seq(2, length(seq_arr), 2)]
  
  con <- file(description=paste(dirpath,CS_filename,sep = ""), open="r")
  cs <- linn <-readLines(con)
  csarr <- unlist(strsplit(cs,split = "\\s"))
  CS <- csarr[length(csarr)]
  CSARR <- substring(CS, seq(1, nchar(CS), 1), seq(1, nchar(CS), 1))
  sprinzlepos <- integer(length = length(CSARR))
  if(substr(CS,1,16)==">>>>>>>..>>>>...")
    for (i in 1:16) 
      sprinzlepos[i] <- i
  # find three most conserved positions from pos 17 to 21

  freqm <- matrix(data = 0,nrow = 4,ncol = 5)
  
  # rows: A,C,T,G
  for (i in 1:length(sequences)) {
    for (j in 1:5) {
      if(substr(sequences[i],j,j)=="A")
      freqm[1,j] <- freqm[1,j]+1
      if(substr(sequences[i],j,j)=="C")
      freqm[2,j] <- freqm[2,j]+1
      if(substr(sequences[i],j,j)=="T")
      freqm[3,j] <- freqm[3,j]+1
      if(substr(sequences[i],j,j)=="G")
      freqm[4,j] <- freqm[4,j]+1
    }
  }
#  [,1] [,2] [,3] [,4] [,5]
#  [1,]   77   65  306  732  745
#  [2,]   90 1796 1704 1455 1088
#  [3,]  392  440  701  649  668
#  [4,] 3017 1275  865  739 1075
  sprinzlepos[17] <- 18
  sprinzlepos[18] <- 19
  sprinzlepos[19] <- 20
  sprinzlepos[20] <- 20
  sprinzlepos[21] <- 21
  if(substr(CS,22,25)=="<<<<")
    for (i in 22:25) 
      sprinzlepos[i] <- i
  sprinzlepos[26] <- 26 # not paired according to sprinzle!
  
  if(substr(CS,27,43)==">>>>>.......<<<<<")
    for (i in 27:43) 
      sprinzlepos[i] <- i
  sprinzlepos[44] <- 44 # not paired according to sprinzle!
  
  # 47 is not present
  sprinzlepos[45] <- 45
  sprinzlepos[46] <- 46
  sprinzlepos[47] <- 48
  
  if(substr(CS,48,64)==">>>>>.......<<<<<")
    for (i in 48:64) 
      sprinzlepos[i] <- i + 1
  if(substr(CS,64,72)=="<<<<<<<<.")
    for (i in 65:72) 
      sprinzlepos[i] <- i + 1
  
}