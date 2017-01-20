#!/bin/bash



#==================================================================================================
# Created on: 2017-01-20
# Usage: annotate_peaks.sh
# Author: Javier Quilez (GitHub: jaquol)
# Goal: annotating a list of ChIP-seq peaks by finding the closest gene to each peak
#==================================================================================================



#==================================================================================================
# (1) Locating data
#==================================================================================================

# peaks
DATA=$HOME/bioinfo_course/session01/data
ifile1="$DATA/gv_009_02_01_chipseq_peaks.narrowPeak"

# genes
genes="$DATA/gencode.v24.annotation_genes.bed"



#==================================================================================================
# (2) Filtering peaks
#==================================================================================================

# the output of the analyses of this session will be saved in the directory defined with:
ANALYSIS=$HOME/bioinfo_course/session01/analysis/annotate_peaks

# cutoff values
min_enrichment=5
min_qval=10	# expressed as -log10(x)

# filter peaks
tfile1="$ANALYSIS/gv_009_02_01_chipseq_peaks_filtered.narrowPeak"
awk -v min_qval=$min_qval -v min_enrichment=$min_enrichment '($7 > min_enrichment) && ($9 > min_qval)' $ifile1 | grep -v "chrM\|chrY\|chrUn\|luciferase" > $tfile1



#==================================================================================================
# (3) Annotate peaks
#==================================================================================================

# as in the previous analysis we are using bedtools but this time the closest tool, which looks for entries in the lists provided in -a and -b that overlap in terms of genomic coordinates  
ofile="$ANALYSIS/gv_009_02_01_chipseq_peaks_annotated.bed"
bedtools closest -a $tfile1 -b $genes > $ofile



#==================================================================================================
# (4) Remove intermediate files
#==================================================================================================

# rm is another built-in Unix function to delete files and directories --use it with caution!!
rm $tfile1
