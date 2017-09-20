#!/bin/bash
#$ -N average_profiles_binsize10
#$ -q short-sl65
#$ -l virtual_free=10G
#$ -l h_rt=06:00:00
#$ -o /users/project/4DGenome/bioinfo_course/session03/job_out/average_profiles_binsize10_$JOB_ID.out
#$ -e /users/project/4DGenome/bioinfo_course/session03/job_out/average_profiles_binsize10_$JOB_ID.err
#$ -j y
#$ -M javier.quilez@crg.eu
#$ -m abe
#$ -pe smp 10
/software/mb/el7.2/anaconda2/bin/computeMatrix reference-point -S /users/GR/mb/jquilez/data/chipseq/samples/gv_009_02_01_chipseq/profiles/hg38_mmtv/single_end/gv_009_02_01_chipseq.rpm.bw /users/GR/mb/jquilez/data/chipseq/samples/gv_066_01_01_chipseq/profiles/hg38_mmtv/single_end/gv_066_01_01_chipseq.rpm.bw -R /users/project/4DGenome/bioinfo_course/session03/tmp_gv_066_01_01_chipseq.bed -out /users/project/4DGenome/bioinfo_course/session03/data/average_profiles_binsize10.mat --referencePoint=center --binSize=10 --upstream=1000 --downstream=1000 --numberOfProcessors=10
/software/mb/el7.2/anaconda2/bin/plotHeatmap -m /users/project/4DGenome/bioinfo_course/session03/data/average_profiles_binsize10.mat -out /users/project/4DGenome/bioinfo_course/session03/figures/average_profiles_binsize10.pdf --plotFileFormat pdf --heatmapHeight 10 --heatmapWidth 8 --xAxisLabel='' --startLabel='' --endLabel='' --colorMap=Blues --refPointLabel='Peak' --regionsLabel T60-peaks --legendLocation=best --samplesLabel Untreated T60
rm -f  /users/project/4DGenome/bioinfo_course/session03/tmp_gv_066_01_01_chipseq.bed
