source('~/RScripts/fn_Library_SN.R')

###############################################################################
## Run fn_Libraries 
###############################################################################
fnLibs <- list.files('~/Project_CurveReg/RScripts_CurveReg/', pattern = '.R')[grep(pattern = 'fn_Library', x = list.files('~/Project_CurveReg/RScripts_CurveReg/', pattern = '.R'))]
for(lib in fnLibs){
  File <- paste0('~/Project_CurveReg/RScripts_CurveReg/', lib)
  source(File)
}

###############################################################################
## Run Registration package libraries
###############################################################################
source('~/R_Packages/Registration/R/fn_Library_Registration.R')
source('~/R_Packages/Registration/R/fn_pairwiseDistance_fdasrvf.R')
source('~/R_Packages/Registration/R/plot_regist_fda.R')
source('~/R_Packages/Registration/R/plot_regist_fdasrvf.R')

###############################################################################
## This function load the bploc file 
###############################################################################
fn_load_bploc <- function(ConversionFactor=206, Filename.bploc){
  bp.loc <- read.table(file=Filename.bploc, header=T, sep=' ')
  bp.loc$BasePairLength <- bp.loc$refMapCoordEnd - bp.loc$refMapCoordStart
  bp.loc$PixelLength_Theo <- round(bp.loc$BasePairLength/ConversionFactor, 0)
  return(bp.loc)
}
###############################################################################

###############################################################################
## This function saves sequence composition related objects
###############################################################################
## The text files saved by this function contains the following elements:
## Chr, FragIndex, GC_VarIndex, GC_pct, Length
fn_saveSeqComp <- function(
  Chr, 
  FragIndex,
  bp.loc,
  BasePairInterval,
  Save = TRUE,
  DataPath
  ){
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

  GCAT <- SeqComp[['SplitSeq_GCAT']]
  GCAT.Long <- melt(data = GCAT, id.vars = 'Base', measure.vars = c('C', 'G', 'A', 'T'))
  colnames(GCAT.Long) <- c('bp', 'base', 'proportion')
  
  GCAT <- as.data.frame(GCAT)
  GCAT$Test <- 0
  GCAT$Test <- apply(X = GCAT, MARGIN = 1, FUN = function(Row){ 
    Val = 0; 
    if(Row[1] + Row[2] < 0.35) {Val = 0.5}; 
    if(Row[1] + Row[2] > 0.65) {Val = 1}; 
    return(Val)}
  )

  GC_VarIndex <- round(mean(GCAT$Test), 2 )

  my.Colors=c('gray25', 'olivedrab1', 'olivedrab4', 'gray56')
  Length_kb <- round((FragBP_End - FragBP_Start)/1000, 4)
  Length_Pixels <- round((FragBP_End - FragBP_Start)/BasePairInterval, 0)
  #if(MainTitle==''){
    MainTitle <- paste(Chr, 'Frag', FragIndex, ',', Length_kb, 'kb Long, Variability:', GC_VarIndex)
  #}

  SeqPlot <- qplot() + geom_bar(aes(y = proportion, x = bp, fill = base),
                   data = GCAT.Long, stat = 'identity') +
  ylab(label = '') + ggtitle(label = MainTitle) + 
  scale_fill_manual(values = my.Colors) +
  theme(legend.position = 'top') 
  
  GC_pct <- round(colMeans(GCAT)[['G']] + colMeans(GCAT)[['C']], 2)
  SeqComp[['SeqPlot']] <- SeqPlot
  SeqComp[['GC_VarIndex']] <- GC_VarIndex
  SeqComp[['GC_pct']] <- GC_pct
  
  FragmentName <- paste(Chr, '_frag', FragIndex, '_SeqComp', sep='')

  FragmentFilename.out <- paste(DataPath, FragmentName, '_GC_Signal.txt', sep='')
  cat(c(Chr, FragIndex, GC_VarIndex, GC_pct, Length_kb, Length_Pixels), file = FragmentFilename.out, append = FALSE)

  if(Save){
    FragmentFilename.out <- paste(DataPath, FragmentName, '.RData', sep='')
    save(SeqComp, file=FragmentFilename.out)
  } else{
    return(SeqComp)
  }

}
###############################################################################

