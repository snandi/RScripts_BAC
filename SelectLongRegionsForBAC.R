rm( list = ls( all.names = TRUE ) )
rm( list = objects( all.names = TRUE ) )
#dev.off()

########################################################################
## This script loads the subIntervals of different chromosomes to select
## long intervals for ordering BACs
########################################################################

loadPairwiseAlignmentScores <- function( 
  outputFoldername,
  filenamePrefix = 'PairwiseAlignmentScores'
){
  filenamePairs <- paste0( outputFoldername, filenamePrefix, '.RData' )
  load( filenamePairs )
  return( subIntervalsPairs )
}

########################################################################
## Load header files and source functions
########################################################################
library(Registration)
PackagesLoaded <- loadPackagesForRegistration()
Packages <- PackagesLoaded$Packages
Packages_Par <- PackagesLoaded$Packages_Par

Packages_Par <- c(
  'fda',
  'fdakma',
  'ggplot2',
  'gtools',
  'plyr',
  'utils'
)

########################################################################
## Define folder paths and source function libraries
########################################################################
RScriptPathCurveReg <- paste0( '~/Project_CurveReg/RScripts_CurveReg/' )
RScriptPathBAC <- paste0( '~/Project_BAC/RScripts_BAC/' )

# source(paste0(RScriptPath, 'CompareSequencesByAlignmentScore.R'))
# source(paste0(RScriptPath, 'LoadMFlorumSequences.R'))
source(paste(RScriptPathCurveReg, 'fn_Library_mm52.R', sep=''))
source(paste(RScriptPathCurveReg, 'fn_Library_CurveReg.R', sep=''))

source(paste(RScriptPathBAC, 'fn_Library_BAC.R', sep=''))

dataPathMF <- '/z/Proj/newtongroup/snandi/MF_cap348/'
dataPathMM <- '/z/Proj/newtongroup/snandi/mm52-all7341/'

########################################################################
## Command line arguments for some static variables
## nohup R CMD BATCH --no-save '--args FLANK_PIXELS=10 MIN_SUBINTERVAL_LENGTH=50 SLIDING_WINDOW=200 GAP_OPENING=-10 CHROMOSOME_INT=15' MM02_CompareSequencesByAlignmentScore.R &
########################################################################
Args <- (commandArgs(TRUE))
if (length(Args)==0) {
  FLANK_PIXELS            <- 10   ## Num of pixels to chop for punctates
  MIN_SUBINTERVAL_LENGTH  <- 50   ## Pixels
  SLIDING_WINDOW          <- 200  ## Sliding window comparison
  GAP_OPENING             <- -10
  CHROMOSOME_INT          <- 1
} else{
  for(i in 1:length(Args)){
    eval(parse(text = Args[[i]]))
  }
}

########################################################################
## Static public variables
########################################################################
CHROMOSOME_STR          <- paste0('chr', CHROMOSOME_INT) 
PREFIX                  <- 'mm'
CONVERSION_FACTOR       <- 206
MIN_COVERAGE            <- 15
nCores                  <- max( 1, detectCores() - 1 )
#FLANK_PIXELS            <- 10 ## Num of pixels to chop for punctates
#MIN_SUBINTERVAL_LENGTH  <- 50 ## Pixels

#GAP_OPENING             <- -20               # Could be command line
GAP_EXTENSION           <- GAP_OPENING
scoreOnlyBoolean        <- TRUE
#SLIDING_WINDOW          <- 200                # Could be command line
SLIDING_INCREMENT       <- 1                 # Could be command line
MFLORUM_FASTA_FILE      <- '~/newtongroup/MF_cap348/fastaSequence/mesoplasma_florum_1.fasta'

## Define the substitution matrix
sigma <- nucleotideSubstitutionMatrix( match = 1, mismatch = 0, baseOnly = TRUE )
sigma2 <- matrix( data = c( 1,0,0,1,0,1,1,0,0,1,1,0,1,0,0,1 ), nrow = 4 )
attributes( sigma2 ) <- attributes( sigma )                   
sigma2
########################################################################

# outputPath <- paste0( outputPath, CHROMOSOME_STR, '/' )

## Load SubIntervals Dataset
filenameSubInt <- paste0( dataPathMM, 'pairwiseComparison/', PREFIX, 'Subint_', CHROMOSOME_STR, 'Pixels', MIN_SUBINTERVAL_LENGTH, 
                          'Flank', FLANK_PIXELS, 'Coverage', MIN_COVERAGE, '.RData' )
load( filenameSubInt )

subset(bpSubIntervalsCoverage, subIntervals > 3)

# fromInterval <- 580
# toInterval <- 591
# subset( bpSubIntervalsCoverage, fragIndex >= fromInterval & fragIndex <= toInterval )
# Min <- min(subset( bpSubIntervalsCoverage, fragIndex >= fromInterval & fragIndex <= toInterval )$bpLocStart)
# Max <- max(subset( bpSubIntervalsCoverage, fragIndex >= fromInterval & fragIndex <= toInterval )$bpLocEnd)
# Min
# Max
# Max - Min

consecutiveFrags <- diff( unique( bpSubIntervalsCoverage$fragIndex ) )
names( consecutiveFrags ) <- unique( bpSubIntervalsCoverage$fragIndex )[-1]
View(as.data.frame(consecutiveFrags))
