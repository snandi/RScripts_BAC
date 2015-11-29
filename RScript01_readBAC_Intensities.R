rm(list=ls(all.names=TRUE))
rm(list=objects(all.names=TRUE))
#dev.off()

########################################################################
## This script registers the intensities of Nmaps aligned to BAC segments
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

########################################################################

Today <- Sys.Date()
ConversionFactor <- 206
BasePairInterval <- ConversionFactor
Chr <- 'Chr7'
ChrNum <- 7

bp.loc <- fn_load_bploc(
  ConversionFactor = ConversionFactor, 
  Filename.bploc = paste0(DataPath, 'bploc_BAC_chr7')
)

########################################################################
## Corresponding Nmaps of the BAC DNA regions exist only for Chr 7
########################################################################
Chr <- 'chr7'
ChrNum <- gsub(pattern = 'chr', replacement = '', x = Chr)

FragIndex <- 1
FragIndices <- c(12437:12447) 
BackbonePixels <- 1
OpticalRes_Factor <- 1

########################################################################
## Load the list of fragements and the number of molecules aligned to them
########################################################################
Filename_fragTable <- paste0('/z/Proj/newtongroup/snandi/mm52-all7341/RData/', Chr, '/', Chr, '_Table.RData')
load(Filename_fragTable)
FragIndices10 <- subset(chr7_Table, numMolecules >= 10)[, 'refStartIndex']
########################################################################

BasePairInterval <- 206*OpticalRes_Factor   ## Length of base pair interval to estimate gcat %

NumBP_Frag <- subset(bp.loc, alignedChr == Chr & alignedFragIndex == FragIndex)[['BasePairLength']] ## Length of frag in BP
#NumSubFrag <- round(NumBP_Frag/BasePairInterval, 0) ## Number of sub fragments
PixelLength_Theo <- subset(bp.loc, alignedChr == Chr & alignedFragIndex == FragIndex)[['PixelLength_Theo']]

# IntensityData_inRange <- fn_saveTruncData(
#   Chr                   = Chr, 
#   FragIndex             = FragIndex, 
#   DataPath.mf           = DataPath.mm52, 
#   Truncate              = FALSE,
#   #  TruncateLength        = 5,
#   StretchPercentAllowed = 50, 
#   Save                  = TRUE, 
#   bp.loc                = bp.loc
# )

#######################################################################
## This for loop just saves the intensities of Nmaps in RData format ##
#######################################################################
# for(FragIndex in FragIndices){
#   print(FragIndex)
#   IntensityData_inRange <- fn_saveTruncData(
#     Chr                   = Chr, 
#     FragIndex             = FragIndex, 
#     DataPath.mf           = DataPath.mm52, 
#     Truncate              = FALSE,
#     #  TruncateLength        = 5,
#     StretchPercentAllowed = 50, 
#     Save                  = TRUE, 
#     bp.loc                = bp.loc
#   )
  fn_smoothByFragment(
    Chr         = Chr, 
    FragIndex   = FragIndex, 
    DataPath.mf = DataPath.mm52, 
    Save        = TRUE
  )
# }

#########################################################################
## Parallelized execution of smoothing the intensity files produced from 
## the above for loop
#########################################################################
## For execution by fragment 
#fn_smoothByFragment(Chr, FragIndex=FragIndex, DataPath.mf, Save=TRUE)

## For parallel execution  
cl <- makeCluster(22)
registerDoParallel(cl)
foreach(FragIndex = FragIndices10[1:44], .inorder=FALSE, .packages=Packages_Par) %dopar% fn_saveTruncData(
      Chr                   = Chr, 
      FragIndex             = FragIndex, 
      DataPath.mf           = DataPath.mm52, 
      Truncate              = FALSE,
      StretchPercentAllowed = 50, 
      Save                  = TRUE, 
      bp.loc                = bp.loc
    )
stopCluster(cl)

cl <- makeCluster(22)
registerDoParallel(cl)
foreach(FragIndex = FragIndices10[1:44], .inorder=FALSE, .packages=Packages_Par) %dopar% fn_smoothByFragment(
  Chr           = Chr, 
  FragIndex     = FragIndex, 
  DataPath.mf   = DataPath.mm52, 
  Save          = TRUE
  )
stopCluster(cl)

