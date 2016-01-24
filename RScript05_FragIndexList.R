rm(list=ls(all.names=TRUE))
rm(list=objects(all.names=TRUE))
#dev.off()

########################################################################
## This script reads the sequence data from the fasta files, creates
## GCAT ratios for each interval, creates the ggplot objects for the
## sequence plots, and the signals (fluctuation of GC percentages)
########################################################################

## Execute this code as follows:
## nohup R CMD BATCH --no-save '--args Chr="chr7" NCores=11' RScript02_SeqCompTest.R chr7_seqcomp.Rout &

########################################################################
## Run Path definition file                                           ##
########################################################################
RScriptPath <- '~/Project_BAC/RScripts_BAC/'
DataPath <- '~/Project_BAC/Data/'
RDataPath <- '~/Project_BAC/RData/'
RPlotPath <- '~/Project_BAC/Plots/'
Filename.Header <- paste('~/Project_BAC/RScripts_BAC/HeaderFile_BAC.R', sep='')
source(Filename.Header)

source(paste(RScriptPath, 'fn_Library_BAC.R', sep=''))
DataPath.mm52 <- '/z/Proj/newtongroup/snandi/mm52-all7341/intensities_inca34_1pixel/'

Packages_Par <- MyAutoLoads

# source('~/R_Packages/Registration/R/loadPackages.R')
# library(rpart)

# Packages_Par <- c(Packages_Par, 'seqinr')
# ########################################################################
Args <- (commandArgs(TRUE))
for(i in 1:length(Args)){
  eval(parse(text = Args[[i]]))
}

Today <- Sys.Date()
ConversionFactor <- 206
BasePairInterval <- ConversionFactor

Chromosomes <- paste0('chr', c(1:22, "X", "Y"))

Chr <- Chromosomes[10]
for(Chr in Chromosomes){
  
  
  #ChrNum <- 7
  ChrNum <- gsub(pattern = 'chr', replacement = '', x = Chr)
  
  bp.loc <- fn_load_bploc(
    ConversionFactor = ConversionFactor, 
    Filename.bploc = paste0('/ua/snandi/human_nMaps/GC_Content/mm52_all7431.goldOnly.bploc_', Chr)
  )
  
  
  ########################################################################
  ## Load the list of fragements and the number of molecules aligned to them
  ########################################################################
  Filename_fragTable <- paste0('/z/Proj/newtongroup/snandi/mm52-all7341/RData/', Chr, '/', Chr, '_Table.RData')
  load(Filename_fragTable)
  Table <- get(paste0(Chr, '_', 'Table'))
  FragIndices10 <- subset(Table, numMolecules >= 10)[, 'refStartIndex']
  N10 <- length(FragIndices10)
  FragIndices20 <- subset(Table, numMolecules >= 20)[, 'refStartIndex']
  
  Frags10 <- cbind(rep(x = ChrNum, times = N10), FragIndices10)
  Filename <- paste0('/ua/snandi/Project_BAC/RScripts_BAC/SubmitFiles/Frags10_', Chr, '.txt')
  write.table(x = Frags10, file = Filename, quote = FALSE, sep = '\t', row.names = FALSE, col.names = FALSE)
}
