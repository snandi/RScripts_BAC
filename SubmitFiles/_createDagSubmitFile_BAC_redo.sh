#!/bin/bash

## This script creates the dag submit file, for those fragments that failed in earlier runs

## Go to /ua/snandi/Project_BAC/RScripts_BAC/SubmitFiles/
## _createDagSubmitFile_BAC_redo.sh Frag10_chr10.txt > dagSubmit_chr10_redo.dag

FILE=$1

LineNum=0
echo "CONFIG dagman_config"
while read line; 
do
    linearray=( $line )
    LineNum=$(( $LineNum + 1 ))
    ChrNum=${linearray[0]}
    FragIndex=${linearray[1]}
    DataFile="/z/Proj/newtongroup/snandi/mm52-all7341/intensities_inca34_1pixel/chr$ChrNum"_"frag$FragIndex"_"SeqComp_GC_Signal.txt"
    if [ ! -f "$DataFile" ];
    then
    	jobID=chr$ChrNum"_"frag$FragIndex
	echo "JOB $jobID /ua/snandi/Project_BAC/RScripts_BAC/SubmitFiles/submitFile_BAC_RS02.txt"
	echo "VARS $jobID ChrNum=\"$ChrNum\""
	echo "VARS $jobID FragIndex=\"$FragIndex\""
	echo " "
done < $FILE

