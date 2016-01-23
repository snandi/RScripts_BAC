rm(list=ls(all.names=TRUE))
rm(list=objects(all.names=TRUE))
#dev.off()

## This is the debug version of RScript02

########################################################################
## This script reads the sequence data from the fasta files, creates
## GCAT ratios for each interval, creates the ggplot objects for the
## sequence plots, and the signals (fluctuation of GC percentages)
########################################################################

########################################################################
## Run Path definition file                                           ##
########################################################################
RScriptPath <- '~/Project_BAC/RScripts_BAC/'
DataPath <- '~/Project_BAC/Data/'
RDataPath <- '~/Project_BAC/RData/'
RPlotPath <- '~/Project_BAC/Plots/'
Filename.Header <- paste('~/RScripts/HeaderFile_lmcg.R', sep='')
source(Filename.Header)
source(paste(RScriptPath, 'fn_Library_BAC.R', sep=''))
DataPath.mm52 <- '/z/Proj/newtongroup/snandi/mm52-all7341/intensities_inca34_1pixel/'

source('~/R_Packages/Registration/R/loadPackages.R')

PackagesLoaded <- loadPackages()
Packages <- PackagesLoaded$Packages
Packages_Par <- PackagesLoaded$Packages_Par
Packages_Par <- c(Packages_Par, 'seqinr')
########################################################################
# Args <- (commandArgs(TRUE))
# for(i in 1:length(Args)){
#   eval(parse(text = Args[[i]]))
# }

Today <- Sys.Date()
ConversionFactor <- 206
BasePairInterval <- ConversionFactor
Chr <- 'chr7'
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
#########################################################################

## The text files saved by this function contains the following elements:
## Chr, FragIndex, GC_VarIndex, GC_pct, Length
SeqComp_Signal <- c()

FragIndex <- 109

for(FragIndex in FragIndices10){
  FragmentName <- paste(Chr, '_frag', FragIndex, '_SeqComp', sep='')
  FragmentFilename <- paste(DataPath.mm52, FragmentName, '_GC_Signal.txt', sep='')
  FragmentFilename.RData <- paste(DataPath.mm52, FragmentName, '.RData', sep='')
  
  SeqComp_Signal <- rbind(SeqComp_Signal, read.table(file = FragmentFilename, header = FALSE))
}
colnames(SeqComp_Signal) <- c('Chr', 'FragIndex', 'GC_VarIndex', 'GC_pct', 'Length_Pixels')
str(SeqComp_Signal)
SeqComp_Signal$GC_VarIndex <- round(SeqComp_Signal$GC_VarIndex, 4)
SeqComp_Signal$GC_pct <- round(SeqComp_Signal$GC_pct, 4)

load(FragmentFilename.RData)
SeqComp[['SeqPlot']]

qplot() + geom_histogram(aes(x = GC_VarIndex), data = SeqComp_Signal)

View(subset(SeqComp_Signal, FragIndex > 15200))
