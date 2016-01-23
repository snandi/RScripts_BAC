## To choose a CRAN Mirror
#chooseCRANmirror(graphics=FALSE)
## Choose 57 for Berkeley, 60 for Iowa, 68 for Revolution R

#install.packages("rneos", repos="http://r-forge.r-project.org/projects/rneos/")

fn_getPackages <- function(){
  Packages <- c(
#     'ade4',                         ## Required by doSNOW            
#     'alr3',
#     'AppliedPredictiveModeling',    ## Used in Coursera machine learning course
#     'animint',                      ## For animated ggplot2 graphics
#     'biganalytics', 
#     'biglm', 
#     'bigmemory', 
    'boot',  
    'car',
#     'caret',                        ## For prediction algorithms like boosting, random forest, etc 
#     'chron', 
#     'cluster', 
# #     'clusterSim',  
#     'clValid', 
#     'DAAG',
#     'devtools',                     ## Only works on R3.2.2 or higher, to install packages from github
# #    'doParallel',
    'doSNOW', 
#     'dtw', 
#     'ElemStatLearn',                ## Used in Coursera machine learning course
# #     'epicalc', 
#     'faraway', 
#     'fAssets', 
#     'fBasics', 
#     'fda', 
#     'fdakma',
#     'fdasrvf',
#     'fda.usc',                      ## For functional data depth
#     'fields', 
#     'fOptions',
    'foreach',
#     'forecast',                     ## for time series modeling  
#     'foreign', 
#     'fpc',                          ## For cluster comparisons ## Wont install on lmcg
# #     'fPortfolio', 
# #     'fRegression', 
# #     'fSeries', 
# #     'fUtilities', 
#     'gam', 
#     'gdata', 
#     'glmnet', 
# #     'gmaps', 
#     'gmodels', 
    'ggplot2',  
    'graphics', 
#     'grid',                         ## Needed by animint
#     'gridExtra', 
#     'gsubfn', 
#     'gtools', 
#     'hexbin',                       ## Needed by animint
#     'HiddenMarkov', 
#     'HMM', 
#     'Hmisc', 
#     'hydroTSM', 
#     'inline', 
#     'iterators', 
#     'jpeg',			                    ## To read/write jpg images
#     'kernlab', 
#     'knitr', 
#     'leaps', 
#     'limma', 
#     'lme4',
#     'lmerTest',            
#     'lmtest',  
# #    'longitudinalData', 
#     'lpSolve', 
#     'lubridate', 
#     'maps', 
#     'maptools', 
#     'mAr', 
#     'matrixStats', 
#     'mclust', 
#     'MEMSS', 
#     'mgcv', 
#     'modeest',
# #    'NCStats',                      ## For chi-square post hoc tests
#     'nlme', 
#     'nortest', 
#     'numDeriv',                     ## To calculate gradient, hessian, jacobian of functions
#     'nws',
#     'orthogonalsplinebasis',
#     'pdfCluster',   	              ## Cluster pdfs from mixtures, non-parametric
#     'PerformanceAnalytics', 
#     'pgmm',                         ## Used in Coursera machine learning course
# #    'phangorn', 
#     'plotrix', 
#     'pls', 
#     'plyr',  
#     'png',                          ## To read and write png images
#     'profr',  
#     'proto',                        ## Needed by animint
#     'psych', 
#     'quantmod', 
#     'quantreg',
#     'rattle',             ## To get better plots of regression tree outputs
#     'R2HTML', 
#     'rbenchmark', 
#     'Rcpp',
#     'readr',  		                  ## Read flat/tabular text files from disk
#     'readstata13',                  ## Read stata 13 objects
    'reshape',  
    'reshape2',  
# #    'RMySQL',   
#     'robustbase', 
#     'robustX',                      ## For multivariate median and other experimental stats
#     'rpart',                        ## Used in Coursera machine learning course
# #    'RJSONIO',                      ## Required by animint
# #    'Rsymphony', 
#     'sandwich', 
#     'scales',                       ## Needed by animint
#     'SenSrivastava',
    'seqinr',                       
#     'servr',                        ## To plot html & js objects produced by animint
#     'sets',                         ## Set operations      
#     'slam',   
#     'snow', 
#     'sos',
#     'speff2trial',
#     'st', 
#     'stats', 
#     'survey', 
#     'svMisc', 
#     'tseries', 
#     'TTR', 
#     'urca',
#     'UsingR',    
#     'vrtest', 
#     'wle', 
# #     'wordcloud', 
# #     'WriteXLS',
#     'xgboost',                      ## Extreme Gradient Boosting
#     'XML', 
#     'xtable', 
#     'xts', 
#     'zipcode', 
     'zoo'
  )
  return(Packages)  
}

Packages <- fn_getPackages()

Packages_Installed <- Packages
## For loop for requiring packages and installing them if something doesnt exist
for(Package in Packages){
  if(require(package=Package, character.only=T) == F){
    print(paste(Package, 'not Installed'))
    Packages_Installed <- Packages_Installed[Packages_Installed != Package]
#     try(install.packages(Package, dependencies = TRUE))
  } else{
    print(paste(Package, 'already exists'))
    require(package=Package, character.only=T)
  }
}

## For parallel processing, when passing the list of packages to load
## in all the cores. Could be different from Packages
MyAutoLoads <- Packages_Installed
Packages_Par <- MyAutoLoads


