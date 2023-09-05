#!/bin/bash

#Sources
PATH_CODE=$PWD
#PATH_RESOURCE="$PATH_CODE/../../Resources"
PATH_RESOURCE="$HOME/Library/CloudStorage/OneDrive-Personal/Marchantia/Codes_Analysis/Phylogenetics/Resources"
PATH_GENOME="$PATH_RESOURCE/01-Genomes/Combined/AllGenomes-merged.fasta"
PATH_ONEKP="$PATH_RESOURCE/02-OneKP/Extracted/OneKP-Nonseed-merged.fasta"
PATH_DONG="$PATH_RESOURCE/03-Dong2022/Combined/Dong2022-merged.fasta"
# REPEAT_BLAST=True #True or False
REPEAT_BLAST=False #True or False

#Create folder for output if not exits
PATH_IN="$PWD/../01-Input"
PATH_OUT="$PWD/../02-Output"

mkdir -p $PATH_OUT
file_timestamp=$PATH_OUT/R1_timestamp.txt

echo "Path_OneKP, $PATH_ONEKP" > $file_timestamp
echo "Path_Genome, $PATH_GENOME" >> $file_timestamp
echo "Path_Dong, $PATH_DONG" >> $file_timestamp


echo "==============================================="
echo "REPEAT_BLAST: $REPEAT_BLAST"
BLAST_OUT_ONEKP="01-OneKP-Nonseed-blast.tsv"
BLAST_OUT_GENOME="01-AnnotatedGenomes-blast.tsv"
BLAST_OUT_DONG="01-Dong-blast.tsv"

if [[ $REPEAT_BLAST == True ]] || [[ ! -f "$PATH_OUT/$BLAST_OUT_GENOME" ]] || [[ ! -f "$PATH_OUT/$BLAST_OUT_DONG" ]] || [[ ! -f "$PATH_OUT/$BLAST_OUT_ONEKP" ]]
then
	#Blast against known genomes
	echo "==============================================="
	echo "Running blast against annotated genome!"
	cd $PATH_OUT
	makeblastdb -in $PATH_GENOME -parse_seqids -dbtype prot -out 00-Genomes_db
	blastp -db 00-Genomes_db -query $PATH_IN/Blast_Input.fasta -out $BLAST_OUT_GENOME -evalue '0.001' -max_target_seqs 50000 -outfmt "6 std ppos" -num_threads 8
	rm 00-Genomes_db.*

	#Blast against OneKP
	echo "Running blast against OneKP!"
	makeblastdb -in $PATH_ONEKP -parse_seqids -dbtype prot -out 00-OneKP_db
	blastp -db 00-OneKP_db -query $PATH_IN/Blast_Input.fasta -out $BLAST_OUT_ONEKP -evalue '0.001' -max_target_seqs 50000 -outfmt "6 std ppos" -num_threads 8
	rm 00-OneKP_db.*

	#Blast against Dong42
	echo "Running blast against Dong42!"
	makeblastdb -in $PATH_DONG -parse_seqids -dbtype prot -out 00-Dong_db
	blastp -db 00-Dong_db -query $PATH_IN/Blast_Input.fasta -out $BLAST_OUT_DONG -evalue '0.001' -max_target_seqs 50000 -outfmt "6 std ppos" -num_threads 8
	rm 00-Dong_db.*
	cd $PATH_CODE


else
	echo "REPEAT_BLAST is False; Skipping blast."
fi

#Create a timestamp for a particular run
runID="R1_$(date +"%Y%m%d_%H%M%S")"
echo "runID, $runID" >> $file_timestamp

#Create folder for output
outDir="${PATH_OUT}/${runID}"
mkdir -p $outDir

#Create log file for each run
LOG_FILE="${outDir}/AutoTest-1.log"
{
echo "Copying scripts.."
cd $outDir

PATH_CODE_COPY=$outDir/00-Codes
mkdir -p $PATH_CODE_COPY
cp $PATH_CODE/AutoTest-1.sh $PATH_CODE_COPY

echo "Copying Blast input.."

PATH_IN_COPY=$outDir/01-Input
mkdir -p $PATH_IN_COPY

PATH_OUT_COPY=$outDir/02-Output
mkdir -p $PATH_OUT_COPY

cp $PATH_IN/Blast_Input.fasta $PATH_IN_COPY


echo "==============================================="
echo "Retrieving seqs!"

cd $PATH_CODE
Rscript 01-RetrieveFasta.R


echo "==============================================="
echo "Merging seqs!"
cd $PATH_CODE
Rscript 02-MergeFasta.R


#Scan for HK domains
echo "==============================================="
echo "Checking with hmmsearch!"
cd $PATH_OUT_COPY
cat $PATH_IN/PF00512.hmm $PATH_IN/PF02518.hmm $PATH_IN/PF00072.hmm $PATH_IN/PF03924.hmm > $PATH_IN_COPY/05-Domains.hmm
hmmsearch --domtblout 03-Candidates-all_dom.tsv $PATH_IN_COPY/05-Domains.hmm 02-Candidates-all.fasta > 03-Candidates-all_dom.out


#Filter by HisKA
echo "==============================================="
echo "Filter sequences by HisKA"
cd $PATH_CODE
Rscript 03-FilterbyHmm.R

# align with MAFFT
echo "==============================================="
echo "Making alignment with mafft!"
cd $PATH_OUT_COPY

mafft --version
mafft --retree 2 03-Candidates-filtered.fasta > 04-Candidates-filtered_FFTNS2.fasta

#Mask alignment result
echo "==============================================="
echo "Masking alignments!"
cd $PATH_CODE
ratio=80
Rscript 04-MaskAlignment.R 1 $ratio

# Test Evolutionary models & Build tree with iqtree
echo "==============================================="
echo "Initial tree with iqtree2!"
cd $PATH_OUT_COPY
iqtree -nt AUTO -ntmax 8 -s ./04-Candidates-filtered_FFTNS2.${ratio}.phy -m TEST -mset phyml -B 1000 -alrt 1000 --prefix ./05-Candidates-filtered_FFTNS2_${ratio}

cd $PATH_CODE
Rscript 05-TranslateTree.R 1
say "I finished!"
# echo "Mafft alignment!"
# cd 03-Retrieved



echo "Test finished!"

}&>$LOG_FILE

