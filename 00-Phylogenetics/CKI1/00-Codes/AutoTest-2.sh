#!/bin/bash
#Directories
PATH_CODE=$PWD
PATH_IN="$PWD/../01-Input"
PATH_OUT="$PWD/../02-Output"

file_timestamp="../02-Output/R2_timestamp.txt"
runID="R2_$(date +"%Y%m%d_%H%M%S")"
echo "runID, $runID" > $file_timestamp

#Create folder for output
outDir="${PATH_OUT}/${runID}"
mkdir -p $outDir

LOG_FILE="${outDir}/AutoTest-2.log"
{
echo "Copying scripts.."
cd $outDir
PATH_CODE_COPY=$outDir/00-Codes
mkdir -p $PATH_CODE_COPY
cp $PATH_CODE/AutoTest-2.sh $PATH_CODE_COPY

PATH_IN_COPY=$outDir/01-Input
mkdir -p $PATH_IN_COPY

PATH_OUT_COPY=$outDir/02-Output
mkdir -p $PATH_OUT_COPY

#Retrieve sequences from the manual subtree
echo "==============================================="
echo "Retrieving sequences!"

cd $PATH_CODE
Rscript 06-ExtractFromSubtree.R

#Run alignment
echo "==============================================="
echo "Aligning with MAFFT!"
cd $PATH_OUT_COPY

mafft_base="./07-Subtree_LINSI"
mafft_out="${mafft_base}.fasta"
if [[ -f $mafft_out ]]; then
	echo "Previous ${mafft_out} exists! Removed!"
	rm $mafft_out
fi

mafft --version
# mafft --thread 8 --genafpair --maxiterate 1000 05-Subtree.fasta > $mafft_out
mafft --thread 8 --localpair --maxiterate 1000 06-Subtree.fasta > $mafft_out

#Mask alignments again
echo "==============================================="
echo "Masking alignments!"
cd $PATH_CODE
ratio=98

mask_base="${mafft_base}.${ratio}"
if [[ -f "${mask_base}.fasta" ]]; then
	echo "Removing previously masked alignments ..."
	ls |grep "${mask_base}.*"
	rm "${mask_base}.*"
fi

Rscript 04-MaskAlignment.R 2 $ratio


# Test Evolutionary models & Build tree with iqtree
echo "==============================================="
echo "Final tree with fast iqtree2!"
cd $PATH_OUT_COPY
tree_base=${mask_base/"07"/"08"}

if [[ -f "${tree_base}.log" ]]; then
	echo "Removing previous iq-tree outputs ..."
	ls |grep "${tree_base}.*"
	rm "${tree_base}.*"
fi

iqtree -nt AUTO -ntmax 8 -s "${mask_base}.phy" -m TEST -mset phyml -B 1000 -alrt 1000 --prefix $tree_base

#Translate the tree
echo "==============================================="
cd $PATH_CODE
Rscript 05-TranslateTree.R 2

# #Run Self blast
# echo "==============================================="
# echo "Running subtree self blast!"
# cd $PATH_OUT_COPY

# makeblastdb -in 05-Subtree.fasta -parse_seqids -dbtype prot -out 05-Subtree_db
# 	blastp -db 05-Subtree_db -query 05-Subtree.fasta -out 08-Subtree-self_blast.tsv -max_target_seqs 50000 -outfmt "6 std ppos" -num_threads 8
# 	rm 05-Subtree_db.*

# echo "==============================================="
# cd $codeDir
#Rscript 06-BlastMatrix.R


say "I finished!"

echo "Test finished!"
}&>$LOG_FILE
