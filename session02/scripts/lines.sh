#!/bin/bash

# settings

source /users/GR/mb/jquilez/.bashrc
cd bioinfo_course/session02/

##  prepare annotation

# set some variables pointing to files

genes=data/gencode.v24.annotation_genes.bed
sizes=data/chromosome.sizes

tss=analysis/annotation/gencode.v24.annotation_genes.tss.bed
proms=analysis/annotation/gencode.v24.annotation_genes.promoters.bed

# get TSS

awk -v OFS="\t" \
	'{if($6 == "+"){tss = $2}else{tss = $3}print $1, tss, tss, $4, $5, $6}' \
	$genes | \
	bedtools sort -i - > \
	$tss

# get promoters (+ - 5 Kbp of TSS)

bedtools slop \
		 -i $tss \
		 -g $sizes \
		 -b 5000 > \
		 $proms

##  overlap PRBS with promoters

cat data/chip_PR_T0.txt data/chip_PR_R30.txt | \
	while read line
	do

		bedtools intersect -u -a $proms -b data/$line > analysis/overlap/$line.tss_overlap.bed

	done

wc -l analysis/overlap/*.tss_overlap.bed

# compare with random

cat data/chip_PR_T0.txt data/chip_PR_R30.txt | \
	while read line
	do

		nreal=$(cat analysis/overlap/$line.tss_overlap.bed | wc -l)
		nrandom=$(bedtools shuffle -i data/$line -g $sizes | \
						 bedtools intersect -u -a $proms -b - | wc -l)
		echo $line $nreal $nrandom

	done > analysis/overlap/promoters_with_PRBD_and_random.txt

# get enrichment

awk '{print $0, $2 / $3}' \
	analysis/overlap/promoters_with_PRBD_and_random.txt > \
	analysis/overlap/enrichment_over_random.txt

##  closest promoter to PRBS

cat data/chip_PR_T0.txt data/chip_PR_R30.txt | \
	while read line
	do

		bedtools closest -d -a $tss  -b data/$line > analysis/closest/$line.closest_tss.bed

	done
