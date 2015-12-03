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
Chr <- 'chr7'

Today <- Sys.Date()
ConversionFactor <- 206
BasePairInterval <- ConversionFactor

bp.loc <- fn_load_bploc(
  ConversionFactor = ConversionFactor, 
  Filename.bploc = paste0('/ua/snandi/human_nMaps/GC_Content/mm52_all7431.goldOnly.bploc_', Chr)
)

########################################################################
## Corresponding Nmaps of the BAC DNA regions exist only for Chr 7
########################################################################
ChrNum <- gsub(pattern = 'chr', replacement = '', x = Chr)

FragIndex <- 12442

FragIndices <- c(12437:12447) 
BackbonePixels <- 1
OpticalRes_Factor <- 1

########################################################################
## Load the list of fragements and the number of molecules aligned to them
########################################################################
Filename_fragTable <- paste0('/z/Proj/newtongroup/snandi/mm52-all7341/RData/', Chr, '/', Chr, '_Table.RData')
load(Filename_fragTable)
Table <- get(paste0(Chr, '_', 'Table'))

bp.loc_BAC <- subset(bp.loc, alignedFragIndex %in% FragIndices)
FragTable <- subset(Table, refStartIndex %in% FragIndices)
########################################################################

BasePairInterval <- 206*OpticalRes_Factor   ## Length of base pair interval to estimate gcat %

NumBP_Frag <- subset(bp.loc, alignedChr == Chr & alignedFragIndex == FragIndex)[['BasePairLength']] ## Length of frag in BP
#NumSubFrag <- round(NumBP_Frag/BasePairInterval, 0) ## Number of sub fragments
PixelLength_Theo <- subset(bp.loc, alignedChr == Chr & alignedFragIndex == FragIndex)[['PixelLength_Theo']]

IntensityData <- fn_loadTruncData(
  Chr = Chr, 
  FragIndex = FragIndex, 
  DataPath.mf = DataPath.mm52, 
  SaveInGlobalEnv = T
)
Regist_PixelFrom <- 46
#Regist_PixelTo <- PixelLength_Theo - 6
Regist_PixelTo <- 85

Data_toRegist <- fn_prepDataForRegist_VarNorm(
  Chr = Chr, 
  FragIndex = FragIndex, 
  DataPath.mf = DataPath.mm52, 
  Regist_PixelFrom = Regist_PixelFrom, 
  Regist_PixelTo = Regist_PixelTo,
  PairwiseSimilarityThreshold = 0.2
)
## Extract all elements of the list
for(Name in names(Data_toRegist)){
  Object <- Data_toRegist[[Name]]
  assign(x=paste('Data', Name, sep='_'), value=Object, envir=.GlobalEnv)
}
PixelPos <- seq(from = Regist_PixelFrom, to = Regist_PixelTo, by = 1)

Data_Intensity_toRegist_D1 <- diff(Data_Intensity_toRegist)
Sim_NoisyCurves <- fn_kma.similarity_mat(
  Mat               = Data_Intensity_toRegist_D1, 
  Xaxis             = PixelPos[-1],
  similarity.method = "d1.pearson"
)
AvgSimilarity <- round(Sim_NoisyCurves[['AvgSimilarity']], 4)

Plot <- fn_plotMultCurves(
  Data = Data_Intensity_toRegist, 
  ColsToPlot = colnames(Data_Intensity_toRegist), 
  XVar = PixelPos, 
  MainTitlePhrase = paste('Avg Similarity:', AvgSimilarity)
)
Plot
