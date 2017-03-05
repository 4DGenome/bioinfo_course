#!/bin/bash

##  This script corresponds to the commands in the ../README.md file.
##  It will execute the main routines executed from the cluster

# Create directories and move to working directory

mkdir -p bioinfo_course/session02/
cd bioinfo_course/session02
mkdir -p analysis/{annotation,closest,overlap}

# Copy the data nedded in the session

cp -r /users/GR/mb/evidal/courses/bioinfo_course/session02/* .

# Load bash configuration

source /users/GR/mb/jquilez/.bashrc

# List data

ls -lhr data

# Show some lines of the gene annotation file

head data/gencode.v24.annotation_genes.bed

# Set file paths as variables

genes=data/gencode.v24.annotation_genes.bed
sizes=data/chromosome.sizes

tss=analysis/annotation/gencode.v24.annotation_genes.tss.bed
proms=analysis/annotation/gencode.v24.annotation_genes.promoters.bed

# Get TSS

awk -v OFS="\t" \
	'{if($6 == "+"){tss = $2}else{tss = $3}print $1, tss, tss, $4, $5, $6}' \
	$genes | \
	bedtools sort -i - > \
	$tss

# Compare TSS with original gene annotation

head $tss

head $genes | cut -f 1-6

# Get promoters

bedtools slop \
		 -i $tss \
		 -g $sizes \
		 -b 5000 > \
		 $proms

# Identify promoters with PRBS

bedtools intersect -u -a $proms -b data/gv_009_01_01_chipseq_peaks.bed
bedtools intersect -u -a $proms -b data/gv_009_01_01_chipseq_peaks.bed | wc -l

for i in data/gv*
do
	echo $i
done
for i in data/gv*
do
	basename $i .bed
	bedtools intersect -u -a $proms -b $i | wc -l
done
for i in data/gv*
do
	outfile=analysis/overlap/$(basename $i .bed).tss_overlap.bed
	bedtools intersect -u -a $proms -b $i > $outfile
done

# Number of promoters with peaks

wc -l analysis/overlap/*.tss_overlap.bed

# Shuffle one peak file

bedtools shuffle -i data/gv_092_01_01_chipseq_peaks.bed -g $sizes

# Shuffle and then overlap with promoters

bedtools shuffle -i data/gv_092_01_01_chipseq_peaks.bed -g $sizes | \
		 bedtools intersect -u -a $proms -b - | \
		 wc -l

# Compare with real overlap

wc -l analysis/overlap/gv_092_01_01_chipseq_peaks.tss_overlap.bed

# Do it for all files

for i in data/gv*
do
	name=$(basename $i .bed)
	tssoverlap=analysis/overlap/$name.tss_overlap.bed

	nreal=$(cat $tssoverlap | wc -l)
	nrandom=$(bedtools shuffle -i $i -g $sizes | \
	    bedtools intersect -u -a $proms -b - | wc -l)
					
	echo $name $nreal $nrandom

done > analysis/overlap/promoters_with_PRBD_and_random.txt

# Compute enrichment over expected in random

awk '{print $0, $2 / $3}' \
	analysis/overlap/promoters_with_PRBD_and_random.txt

# Get closest PRBS to TSS

for i in data/gv*bed
do
	name=$(basename $i .bed)
	bedtools closest -d -t first -a $tss -b $i > analysis/closest/$name.closest_tss.bed
done
