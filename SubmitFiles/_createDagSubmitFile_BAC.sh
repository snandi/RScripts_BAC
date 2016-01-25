#!/bin/bash

## This script creates the dag submit file, for all fragments in the Frag10 file

## Go to /ua/snandi/Project_BAC/RScripts_BAC/SubmitFiles/
## _createDagSubmitFile_BAC.sh Frag10_chr10.txt > dagSubmit_chr10.dag

FILE=$1

LineNum=0
echo "CONFIG dagman_config"
while read line; 
do
    linearray=( $line )
    LineNum=$(( $LineNum + 1 ))
    ChrNum=${linearray[0]}
    FragIndex=${linearray[1]}
    jobID=chr$ChrNum"_"frag$FragIndex
    echo "JOB $jobID /ua/snandi/Project_BAC/RScripts_BAC/SubmitFiles/submitFile_BAC_RS02.txt"
    echo "VARS $jobID ChrNum=\"$ChrNum\""
    echo "VARS $jobID FragIndex=\"$FragIndex\""
    echo " "
done < $FILE



