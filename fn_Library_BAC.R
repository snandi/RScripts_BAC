source('~/RScripts/fn_Library_SN.R')

########################################################
## Run fn_Libraries 
########################################################
fnLibs <- list.files('~/Project_CurveReg/RScripts_CurveReg/', pattern = '.R')[grep(pattern = 'fn_Library', x = list.files('~/Project_CurveReg/RScripts_CurveReg/', pattern = '.R'))]
for(lib in fnLibs){
  File <- paste0('~/Project_CurveReg/RScripts_CurveReg/', lib)
  source(File)
}

########################################################
## Run Registration package libraries
########################################################
source('~/R_Packages/Registration/R/fn_Library_Registration.R')
source('~/R_Packages/Registration/R/fn_pairwiseDistance_fdasrvf.R')
source('~/R_Packages/Registration/R/plot_regist_fda.R')
source('~/R_Packages/Registration/R/plot_regist_fdasrvf.R')




