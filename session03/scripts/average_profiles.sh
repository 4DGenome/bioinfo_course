#!/bin/bash


#==================================================================================================
# Created on: 2017-09-20
# Usage: ./average_profiles.sh
# Author: Javier Quilez (GitHub: jaquol)
# Goal: generate ChIP-seq tag counts profiles centered the peaks
#==================================================================================================



#==================================================================================================
# CONFIGURATION VARIABLES AND PATHS
#==================================================================================================

# variables
profiles="gv_009_02_01_chipseq gv_066_01_01_chipseq"
profiles_names="Untreated T60"
peaks="gv_066_01_01_chipseq"
peaks_names="T60-peaks"
process=average_profiles
peak_caller=macs2
data_type=chipseq
version=hg38_mmtv
sequencing_type=single_end
peak_calling_mode=with_control
bin_size=10
n_regions=1000
flanking_region=1000

# Paths
SAMPLES=/users/GR/mb/jquilez/data/chipseq/samples
PWD=`pwd`
ANALYSIS=$PWD
JOB_CMD=$ANALYSIS/job_cmd 
JOB_OUT=$ANALYSIS/job_out
mkdir -p $JOB_CMD
mkdir -p $JOB_OUT
mkdir -p $ANALYSIS/data
mkdir -p $ANALYSIS/figures
compute_matrix=`which computeMatrix`
plot_heatmap=`which plotHeatmap`

# Cluster parameters
queue=short-sl65
memory=10G
max_time=06:00:00
slots=10



#==================================================================================================
# COMMANDS
#==================================================================================================

# prepare regions
ibeds=""
for s in $peaks; do
	tbed=$ANALYSIS/tmp_$s.bed
	ibed=$SAMPLES/$s/peaks/$peak_caller/$version/$peak_calling_mode/$sequencing_type/${s}_peaks.narrowPeak
	shuf -n $n_regions $ibed > $tbed 
	ibeds="$ibeds $tbed"
done

# prepare read per million profiles
ibws=""
for s in $profiles; do
	ibw=$SAMPLES/$s/profiles/$version/$sequencing_type/$s.rpm.bw
	ibws="$ibws $ibw"
done

# Build job: parameters
job_name=${process}_binsize${bin_size}
job_file=$JOB_CMD/$job_name.sh
m_out=$JOB_OUT
echo "#!/bin/bash
#$ -N $job_name
#$ -q $queue
#$ -l virtual_free=$memory
#$ -l h_rt=$max_time
#$ -o $m_out/${job_name}_\$JOB_ID.out
#$ -e $m_out/${job_name}_\$JOB_ID.err
#$ -j y
#$ -M javier.quilez@crg.eu
#$ -m abe
#$ -pe smp $slots" > $job_file

# deeptools commands: compute matrix
# -S = bigWig containing the scores (reads per million here) to plot
# -R = BED file containing the regions to plot
omat=$ANALYSIS/data/$job_name.mat
job_cmd="$compute_matrix reference-point \
			-S $ibws \
			-R $ibeds \
			-out $omat \
			--referencePoint=center \
			--binSize=$bin_size \
			--upstream=$flanking_region \
			--downstream=$flanking_region \
			--numberOfProcessors=$slots"
echo $job_cmd >> $job_file

# deeptools commands: plot profile
# -m = matrix from the computeMatrix tool
opdf=$ANALYSIS/figures/${job_name}.pdf
job_cmd="$plot_heatmap \
			-m $omat \
			-out $opdf \
			--plotFileFormat pdf \
			--heatmapHeight 10 \
			--heatmapWidth 8 \
			--xAxisLabel='' \
			--startLabel='' \
			--endLabel='' \
			--colorMap=Blues \
			--refPointLabel='Peak' \
			--regionsLabel $peaks_names \
			--legendLocation=best \
			--samplesLabel $profiles_names"
echo $job_cmd >> $job_file

echo "rm -f $ibeds" >> $job_file

# Submit job
chmod a+x $job_file
#qsub < $job_file
$job_file

