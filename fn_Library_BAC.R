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



