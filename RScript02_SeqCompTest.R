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
#Chr <- 'chr7'
#ChrNum <- 7
ChrNum <- gsub(pattern = 'chr', replacement = '', x = Chr)

bp.loc <- fn_load_bploc(
  ConversionFactor = ConversionFactor, 
  Filename.bploc = paste0('/ua/snandi/human_nMaps/GC_Content/mm52_all7431.goldOnly.bploc_', Chr)
)

########################################################################
## Corresponding Nmaps of the BAC DNA regions exist only for Chr 7
########################################################################

FragIndex <- 5
FragIndices <- c(12437:12447) 
BackbonePixels <- 1
OpticalRes_Factor <- 1

########################################################################
## Load the list of fragements and the number of molecules aligned to them
########################################################################
Filename_fragTable <- paste0('/z/Proj/newtongroup/snandi/mm52-all7341/RData/', Chr, '/', Chr, '_Table.RData')
load(Filename_fragTable)
Table <- get(paste0(Chr, '_', 'Table'))
FragIndices10 <- subset(Table, numMolecules >= 10)[, 'refStartIndex']
FragIndices20 <- subset(Table, numMolecules >= 20)[, 'refStartIndex']
#########################################################################

BasePairInterval <- 206*OpticalRes_Factor   ## Length of base pair interval to estimate gcat %

NumBP_Frag <- subset(bp.loc, alignedChr == Chr & alignedFragIndex == FragIndex)[['BasePairLength']] ## Length of frag in BP
#NumSubFrag <- round(NumBP_Frag/BasePairInterval, 0) ## Number of sub fragments
PixelLength_Theo <- subset(bp.loc, alignedChr == Chr & alignedFragIndex == FragIndex)[['PixelLength_Theo']]

# fn_saveSeqComp(
#  Chr = Chr, 
#  FragIndex = FragIndex,
#  bp.loc = bp.loc,
#  BasePairInterval = BasePairInterval,
#  Save = TRUE,
#  DataPath = DataPath.mm52
# )

#########################################################################
## Parallelized execution of saving sequence compositions list of objects
#########################################################################
## For parallel execution, using doParallel
#cl <- makeCluster(NCores)
#cl <- makePSOCKcluster(NCores)
#doParallel::registerDoParallel(cl)
Time1 <- Sys.time()
# For parallel execution, using doSNOW
cl <- makeCluster(NCores, type = "SOCK")
doSNOW::registerDoSNOW(cl)

foreach(FragIndex = FragIndices10[1:20], .inorder = FALSE, .packages = Packages_Par) %dopar% fn_saveSeqComp(
  Chr = Chr, 
  FragIndex = FragIndex,
  bp.loc = bp.loc,
  BasePairInterval = BasePairInterval,
  Save = TRUE,
  DataPath = DataPath.mm52
)

stopCluster(cl)
print(Sys.time() - Time1)
## The text files saved by this function contains the following elements:
## Chr, FragIndex, GC_VarIndex, GC_pct, Length_kb, Length_Pixels
