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
RDataPath <- '~/Project_BAC/RData/'
RPlotPath <- '~/Project_BAC/Plots/'
Filename.Header <- paste('~/RScripts/HeaderFile_lmcg.R', sep='')
source(Filename.Header)
source(paste(RScriptPath, 'fn_Library_BAC.R', sep=''))
########################################################################

Today <- Sys.Date()
ConversionFactor <- 206
BasePairInterval <- ConversionFactor
Chr <- 'Chr7'
ChrNum <- 7

bploc.folder <- paste(AspenDataPath, 'alignmentChunks/', sep='')
BPLocs_WholeGenome <- fn_readBPloc(Folder=bploc.folder, File='mm52_all7431.goldOnly.bploc')
BPLocs_LINE_Chr13 <- fn_load_bplocs_LINE(Folder='~/Project_GC_Content/LINE/', 
                                         Chr=Chr, 
                                         RepFamily=RepFamily, 
                                         ConversionFactor=ConversionFactor, 
                                         BPLocs_WholeGenome=BPLocs_WholeGenome)
#View(subset(BPLocs_LINE_Chr13, BP_Length>6000))

#xtable(subset(BPLocs_LINE_Chr13, BP_Length>6000)[15:30,3:8], digits=0)
########################################################################
## Chosing: Chr 13, Frag 2168 - 2184, has 3 LINE segments, of more than
## 6KB length, but only 2168 exists in the dataset. 
## Hence, choosing 3178 through 3190
########################################################################
DownThreshold <- 0.75  ## To eliminate outliers
UpThreshold <- 1.25
Chr <- 'chr13'
FragIndex <- 38
FragIndices <- c(3174:3199) ## These fragments have 2 LINE segments
BackbonePixels <- 1
OpticalRes_Factor <- 1

FragIndices <- scan(file='~/Project_GC_Content/RData/chr13_fragIndexList_Min10.txt')

# OpticalRes_Factor <- 1
# BasePairInterval <- 206*OpticalRes_Factor   ## Length of base pair interval to estimate gcat %
# 
# FragBP_Start <- BPLocs_WholeGenome[which(BPLocs_WholeGenome$alignedChr==Chr & 
#                                            BPLocs_WholeGenome$alignedFragIndex == FragIndex), 
#                                    'refMapCoordStart'] 
# FragBP_End <- BPLocs_WholeGenome[which(BPLocs_WholeGenome$alignedChr==Chr & 
#                                          BPLocs_WholeGenome$alignedFragIndex == FragIndex), 
#                                  'refMapCoordEnd']
# NumBP_Frag <- FragBP_End - FragBP_Start ## Length of frag in BP
# NumSubFrag <- round(NumBP_Frag/BasePairInterval, 0)

#######################################################################
## This for loop just saves the intensities of Nmaps in RData format ##
#######################################################################
for(FragIndex in FragIndices){
  print(FragIndex)
  InputDataPath <- paste(AspenDataPath, 'intensities_inca34_', 
                         BackbonePixels, 'pixel/', sep='')
  IntensityData.Raw <- 
    fn_returnIntervalIntensity(Chr=Chr, FragIndex=FragIndex, 
                               RDataPath=InputDataPath, 
                               TruncateLength=5)
  if(is.null(IntensityData.Raw)){
    next
  } else{
    IntensityData <- 
      fn_eliminateIntensityOutliers(IntensityData=IntensityData.Raw, 
                                    DownThreshold=0.75, UpThreshold=1.25)[['IntensityData']]
    
    NumPixels <- aggregate(IntensityData$PixelNum/OpticalRes_Factor,
                           by=list(IntensityData$MoleculeID), FUN=max)
    
    colnames(NumPixels) <- c('MoleculeID', 'Pixels')
    MoleculesInFrag <- NumPixels
    MoleculesInFrag$Pixels <- 1
    colnames(MoleculesInFrag)[colnames(MoleculesInFrag)=='Pixels'] <- paste('Frag', FragIndex, sep='_')
    if('OverlapMolecules' %notin% ls()){
      OverlapMolecules <- MoleculesInFrag
    } else{
      OverlapMolecules <- merge(OverlapMolecules, MoleculesInFrag, by='MoleculeID', 
                                all=T)    
    }
    
    IntensityData.Aligned.Wide <- reshape(IntensityData[,-1],
                                          timevar='MoleculeID',
                                          idvar='PixelNum',
                                          direction='wide')
    
    Filename <- paste(RDataPath, Chr, '/', Chr, '_frag', FragIndex, '_intensities_BB_', 
                      BackbonePixels, 'pixels.RData', sep='')
    save(IntensityData.Aligned.Wide, file=Filename)
    
    ############################ PLOTTING ################################
    FragBP_Start <- BPLocs_WholeGenome[which(BPLocs_WholeGenome$alignedChr==Chr & 
                                               BPLocs_WholeGenome$alignedFragIndex == FragIndex), 
                                       'refMapCoordStart'] 
    FragBP_End <- BPLocs_WholeGenome[which(BPLocs_WholeGenome$alignedChr==Chr & 
                                             BPLocs_WholeGenome$alignedFragIndex == FragIndex), 
                                     'refMapCoordEnd']
    NumBP_Frag <- FragBP_End - FragBP_Start ## Length of frag in BP
    NumSubFrag <- round(NumBP_Frag/BasePairInterval, 0)
    SeqComp <- fn_returnSeqComp(Chr=Chr, FragIndex=FragIndex,
                                Interval=BasePairInterval,
                                DataPath=DataPath,
                                numPixels=NumSubFrag,
                                paste('/omm/data/sequence/human_wchr-b37/chr3.fa'),
                                FragBP_Start=FragBP_Start,
                                FragBP_End=FragBP_End)
    PlotGC <- fn_createGCATPlot(SeqComp=SeqComp, xlab='', FragBP_Start, FragBP_End, 
                                FragIndex, Chr)
    PixelLen_byFrag <- fn_plotPixelLength_byFrag(NumPixels, FragIndex, Chr)
    Filename.pdf <- paste(RPlotPath, 'human_', Chr, '/', Chr, '_Frag', FragIndex, 
                          '_BB_', BackbonePixels, 'pixels_', Today, '.pdf', sep='')
    try(fn_plotIntensities(IntensityData=IntensityData, PlotGC=PlotGC, 
                       Filename.pdf=Filename.pdf, PixelLen_byFrag))
  }
}


