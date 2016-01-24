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
Filename.Header <- paste('~/Project_BAC/RScripts_BAC/HeaderFile_BAC.R', sep='')
source(Filename.Header)

source(paste(RScriptPath, 'fn_Library_BAC.R', sep=''))
DataPath.mm52 <- '/z/Proj/newtongroup/snandi/mm52-all7341/intensities_inca34_1pixel/'

Packages_Par <- MyAutoLoads

# source('~/R_Packages/Registration/R/loadPackages.R')
# library(rpart)
########################################################################
# Arguments: Chr, FragIndex
Args <- (commandArgs(TRUE))
for(i in 1:length(Args)){
  eval(parse(text = Args[[i]]))
}
print(ChrNum)
Chr <- paste0('chr', ChrNum)
print(Chr)
print(FragIndex)

Today <- Sys.Date()
ConversionFactor <- 206
BasePairInterval <- ConversionFactor
#Chr <- 'chr22'
#ChrNum <- 7
#ChrNum <- gsub(pattern = 'chr', replacement = '', x = Chr)


bp.loc <- fn_load_bploc(
  ConversionFactor = ConversionFactor, 
  Filename.bploc = paste0('/ua/snandi/human_nMaps/GC_Content/mm52_all7431.goldOnly.bploc_', Chr)
)

BackbonePixels <- 1
OpticalRes_Factor <- 1

BasePairInterval <- 206*OpticalRes_Factor   ## Length of base pair interval to estimate gcat %

NumBP_Frag <- subset(bp.loc, alignedChr == Chr & alignedFragIndex == FragIndex)[['BasePairLength']] ## Length of frag in BP
#NumSubFrag <- round(NumBP_Frag/BasePairInterval, 0) ## Number of sub fragments
PixelLength_Theo <- subset(bp.loc, alignedChr == Chr & alignedFragIndex == FragIndex)[['PixelLength_Theo']]

SeqComp <- fn_saveSeqComp(
 Chr              = Chr, 
 FragIndex        = FragIndex,
 bp.loc           = bp.loc,
 BasePairInterval = BasePairInterval,
 Save             = TRUE,
 DataPath         = DataPath.mm52
)

