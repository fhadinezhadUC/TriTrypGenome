Gtest.Genomes <- function(){
  library(FSA)
  library(DescTools)
  genefile <-
    read.table("/home/fatemeh/TriTrypGenome_project/GeneAnnotation/tsfm_input_geneset.txt", header = TRUE, colClasses = "character")
  
  # assign a clade to each gene/row of the geneset
  genefile$clade <- ""
  for (i in 1:nrow(genefile)) {
    clusternames <- genefile$sourceOrg[i]
    if (clusternames == "LspMARLEM2494" |
        clusternames == "LenriettiiLEM3045")
      genefile$clade[i] <- "Lenriettii"
    if (clusternames == "TbruceigambienseDAL972" |
        clusternames == "TbruceiLister427" |
        clusternames == "TbruceiTREU927" |
        clusternames == "TevansiSTIB805" |
        clusternames == "TcongolenseIL3000" |
        clusternames == "TvivaxY486")
      genefile$clade[i] <- "AfricanT"
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
      genefile$clade[i] <- "AmericanT"
    if (clusternames == "CfasciculataCfCl" |
        clusternames == "LseymouriATCC30220" |
        clusternames == "LpyrrhocorisH10")
      genefile$clade[i] <- "Leptospira"
    if (clusternames == "LmajorFriedlin" |
        clusternames == "LmajorLV39c5" |
        clusternames == "LmajorSD75" |
        clusternames == "LturanicaLEM423" |
        clusternames == "LarabicaLEM1108" |
        clusternames == "LtropicaL590" |
        clusternames == "LaethiopicaL147" |
        clusternames == "LgerbilliLEM452")
      genefile$clade[i] <- "Lmajor"
    if (clusternames == "LdonovaniBHU1220" |
        clusternames == "LdonovaniBPK282A1" |
        clusternames == "LinfantumJPCM5")
      genefile$clade[i] <- "Linfantum"
    if (clusternames == "LamazonensisMHOMBR71973M2269" |
        clusternames == "LmexicanaMHOMGT2001U1103")
      genefile$clade[i] <- "LMexicana"
    if (clusternames == "LbraziliensisMHOMBR75M2904" |
        clusternames == "LbraziliensisMHOMBR75M2903" |
        clusternames == "LpanamensisMHOMPA94PSC1" |
        clusternames == "LpanamensisMHOMCOL81L13")
      genefile$clade[i] <- "LViannia"
  }
  
  TryTripclades <- c(
    "Lmajor",
    "Linfantum",
    "LMexicana",
    "LViannia",
    "Lenriettii",
    "Leptospira",
    "AmericanT",
    "AfricanT"
  )
  pvalues <- numeric(length = length(TryTripclades))
  
  # run Gtest for each clade and save the pvalues in pvalues
  GtestResults <- list()
  tempgenefile <- genefile
  for (j in 1:length(TryTripclades)) {
    genefile <- tempgenefile[tempgenefile$clade==TryTripclades[j],]
    comp_df <- data.frame(table(genefile$sourceOrg))
    names(comp_df)<- c("genome","genes")
    comp_df$TotalGp <- 0
    comp_df$TotalCp <- 0
    comp_df$TotalTp <- 0
    comp_df$genesLen <- 0
    for (i in 1:nrow(comp_df)) {
      genes <- paste(genefile[genefile$sourceOrg==comp_df$genome[i],]$genesequence,sep = "",collapse = "")
      comp_df$genesLen[i] <- nchar(genes)
      comp_df$TotalGp[i] <- Gpercentage(genes)
      comp_df$TotalCp[i] <- Cpercentage(genes)
      comp_df$TotalTp[i] <- Tpercentage(genes)
    }
    
    Gfreq <- comp_df$TotalGp * comp_df$genesLen
    Cfreq <- comp_df$TotalCp * comp_df$genesLen
    Tfreq <- comp_df$TotalTp * comp_df$genesLen
    Afreq <- (100 - comp_df$TotalGp - comp_df$TotalCp - comp_df$TotalTp) * comp_df$genesLen
    M <- as.table(cbind(Gfreq, Cfreq, Tfreq,Afreq))
    dimnames(M) <-
      list(
        clades =  clades <-
          comp_df$genome,
        nucl = c("G", "C", "T","A")
      )
    #library(DescTools)
    Xsq <- GTest(M)
    GtestResults[[j]] <- Xsq
    pvalues[j] <- Xsq$p.value
  }
  
  # a dataframe of clades and pvalues
  pvalue_df<- data.frame(TryTripclades,pvalues)
  # adjusted pvalues
  pvalue_df$BH_adjusted_pvalues <- p.adjust(pvalue_df$pvalues, method = "BH")
  
  # results:
  # TryTripclades      pvalues BH_adjusted_pvalues
  # 1        Lmajor 0.0278723370        0.0445957392
  # 2     Linfantum 0.6493118480        0.6493118480
  # 3     LMexicana 0.0001338053        0.0002676105
  # 4      LViannia 0.0001110509        0.0002676105
  # 5    Lenriettii 0.6180969617        0.6493118480
  # 6    Leptospira 0.0721779754        0.0962373006
  # 7     AmericanT 0.0000000000        0.0000000000
  # 8      AfricanT 0.0000000000        0.0000000000
  
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