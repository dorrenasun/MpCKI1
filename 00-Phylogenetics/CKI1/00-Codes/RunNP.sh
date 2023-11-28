#!/bin/bash
#Directories
PATH_CODE=$PWD
PATH_PARENT=$(cd ../ && pwd)
PATH_IN="$PATH_PARENT/01-Input"
PATH_OUT="$PATH_PARENT/02-Output"

#Extract runID from timestamp
runID=$(awk '{print $2}' ${PATH_OUT}/R2_timestamp.txt)
outDir="${PATH_OUT}/${runID}"
PATH_OUT_COPY="${outDir}/02-Output"

cd $PATH_OUT_COPY

phy_file=$(ls |grep '.phy')
tree_base_np=${PATH_OUT_COPY}/$(basename $phy_file .phy)_np
echo $tree_base_np
echo "Analyzing $phy_file ..."

# iqtree -nt AUTO -ntmax 12 -s $phy_file -m TEST -mset phyml -b 1000 --prefix $tree_base_np

#IPV urls: https://www.ebi.ac.uk:443/interpro//result/InterProScan/iprscan5-R20231115-103500-0336-85892870-p1m/
# https://www.ebi.ac.uk:443/interpro//result/InterProScan/iprscan5-R20231115-103500-0336-85892870-p1m/

#Translate the tree
echo "==============================================="
cd $PATH_CODE
Rscript 05-TranslateTree.R 2
