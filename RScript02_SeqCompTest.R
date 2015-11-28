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

FragIndex <- 12440
FragIndices <- c(12437:12447) 
BackbonePixels <- 1
OpticalRes_Factor <- 1

#FragIndices <- scan(file='~/Project_GC_Content/RData/chr13_fragIndexList_Min10.txt')

BasePairInterval <- 206*OpticalRes_Factor   ## Length of base pair interval to estimate gcat %

NumBP_Frag <- subset(bp.loc, alignedChr == Chr & alignedFragIndex == FragIndex)[['BasePairLength']] ## Length of frag in BP
#NumSubFrag <- round(NumBP_Frag/BasePairInterval, 0) ## Number of sub fragments
PixelLength_Theo <- subset(bp.loc, alignedChr == Chr & alignedFragIndex == FragIndex)[['PixelLength_Theo']]

#######################################################################
## This for loop just saves the intensities of Nmaps in RData format ##
#######################################################################
#for(FragIndex in FragIndices){
  print(FragIndex)
  
  ############################ PLOTTING ################################
  FragBP_Start <- bp.loc[which(bp.loc$alignedChr==Chr & 
                                 bp.loc$alignedFragIndex == FragIndex), 
                         'refMapCoordStart'] 
  FragBP_End <- bp.loc[which(bp.loc$alignedChr==Chr & 
                               bp.loc$alignedFragIndex == FragIndex), 
                       'refMapCoordEnd']
  NumBP_Frag <- FragBP_End - FragBP_Start ## Length of frag in BP
  NumSubFrag <- round(NumBP_Frag/BasePairInterval, 0)
  SeqComp <- fn_returnSeqComp(Chr=Chr, FragIndex=FragIndex,
                              Interval=BasePairInterval,
                              numPixels=NumSubFrag,
                              Filename=paste0('~/human_nMaps/SequenceData/', Chr, '.fa'),
                              FragBP_Start=FragBP_Start,
                              FragBP_End=FragBP_End)

 ## PlotGC <- fn_createGCATPlot(SeqComp=SeqComp, xlab='', FragBP_Start, FragBP_End, 
 ##                              FragIndex, Chr)
  GCAT <- SeqComp[['SplitSeq_GCAT']]
  GCAT.Long <- melt(data = GCAT, id.vars = 'Base', measure.vars = c('C', 'G', 'A', 'T'))
colnames(GCAT.Long) <- c('bp', 'base', 'proportion')
str(GCAT.Long)

  #ggplot(aes(x = X1, y = value, fill = X2), data = Data) + geom_bar()
  #qplot(proportion, data = GCAT.Long, geom="bar", fill=factor(base))

SeqPlot <- qplot() + geom_bar(aes(y = proportion, x = bp, fill = base),
                   data = GCAT.Long, stat = 'identity')

GCAT <- as.data.frame(GCAT)
GCAT$Test <- 0
GCAT$Test <- apply(X = GCAT, MARGIN = 1, FUN = function(Row){(Row[1] + Row[2] < 0.35) || (Row[1] + Row[2] > 0.65)})
mean(GCAT$Test)

2*mean(as.matrix(GCAT[,c('C', 'G')]))
