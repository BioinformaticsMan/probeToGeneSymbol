#!/usr/bin/bash

##Description: This script aims to extract the expression data from the matrix file with basic process
## the mean value of a probe will be reserved. 
##Useage: bash run.sh

dir=~/data/GSE131617
matrixFile=$dir/GSE131617-GPL5175_series_matrix.txt
familyFile=$dir/GSE131617_family.soft
name=`basename $matrixFile .txt`
outFile=${name}_anno.txt
traitFile=${name}_trait.txt

python=~/Biosoft/miniconda2/envs/py2/bin/python3
script=probeToGeneSymbol.py

$python $script $matrixFile $familyFile $outFile $traitFile
