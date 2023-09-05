# Source codes for phylogenetic analysis of CKI1 histidine kinase
This directory include scripts used for the phylogenetic analysis of KO (CYP701) homologs.

Working pipeline (All scripts should be run from the 00-Codes folder):

- Step 1: Run "Autotest-1.sh" to generate a initial tree (output file: 04-Candidates-all_FFTNS2_80.treefile.newick)
	1. blastp to search in different local databases. Threshold: e<1e-3. All hits were kept for the initial analysis.
	2. MAFFT alignment using the FFTNS2 algorithm. Columns with >80% gaps were masked from the alignment.
	3. Phylogenetic inference with IQ-TREE 2 and Ultrafast bootstrap.
- Step 2-1: Manually select three subtrees from the initial tree (easily done with Figtree), and save them as a nexus files in the "03-KeyOutput" folder.
	The three subtrees are as follows:
	1. the core clade containing all candidate CKI1 and CKI1L homologs, saved as "06-Subtree-Core.nex";
	2. the larger clade including the core clade and its close neighbors with Marchantia sequences, saved as "06-Subtree-Neighbors.nex"; 
	3. the outgroup clade containing putative CHK sequences, saved as "06-Subtree-Outgroup.nex"
- Step 2-2: Copy the following files from outputs of Step 1 to the "03-KeyOutput" folder:
	- "02-Candidates-all.fasta"
	- "02-Candidates-all.id"
	- "02-Sources.tsv"
	Edit the last file to keep genomes/transcriptomes that will be used in the next step of analysis. 

- Step 3: Run "Autotest-2.sh" to test if the sequence selection works okay.
	1. MAFFT alignment using the LINSI or EINSI algorithm. Columns with >80% gaps were masked from the alignment.
	2. Phylogenetic inference with IQ-TREE 2 and Ultrafast bootstrap
- Step 4: Run standard non-parametric analysis with IQ-TREE 2
```
	iqtree -nt AUTO -ntmax 16 -s ./06-Subtree_EINSI.80.phy -m TEST -mset phyml -b 1000 --prefix ./07-Subtree_EINSI_80_np
```
- Step 5: Adjust the tree branch order manually and reroot the tree. Save the adjusted tree as a newick file "07-Subtree_EINSI_80_np.treefile.ordered.newick"

- Step 6: Run interpro analysis with the sequence files generated in Step 3.

- Step 6: Visualize the tree with "07-Visualization.R". Pay attention to the input files or directory in the script.
- Step 7 (optional): Translate the names of sequence alignments etc. using "08-TranslateAlignment.R". Pay attention to the input files or directory in the script.


The files used in creating the phylogenetic tree of the manuscript is in the directory "ForPublication".
- CKI1-all_seqs.fasta: fasta file for candidate sequences of the final tree
- CKI1-EINSI.fasta: file for sequence alignment, before masking
- CKI1-EINSI-mask80.fasta: file for sequence alignment, masked & used in phylogenetic inference
- CKI1-tree.newick: file for the final phylogenetic tree, in a machine-readable form