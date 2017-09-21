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
profiles="gv_009_02_01_chipseq gv_066_01_01_chipseq"	# samples for which profiles will be used
profiles_names="Untreated T60"							# names that will be shown in the plot for the profiles
peaks="gv_066_01_01_chipseq"							# samples for which peak regions will be used
peaks_names="T60-peaks"									# names that will be shonw in the plot for the regions
process=average_profiles								# name of this analysis (whatever)
peak_caller=macs2										# variable to tune the selected profile/peak (ignore)
data_type=chipseq 										# variable to tune the selected profile/peak (ignore)
version=hg38_mmtv										# variable to tune the selected profile/peak (ignore)
sequencing_type=single_end								# variable to tune the selected profile/peak (ignore)
peak_calling_mode=with_control							# variable to tune the selected profile/peak (ignore)
bin_size=10												# bin size in bp for computeMatrix
n_regions=1000											# number of random regions selected	
flanking_region=1000									# number of bp up and downstream of the peak region
email=javier.quilez@crg.eu								# email to notify begin, end and abortion of job in the cluster

# Paths (input, output and tools paths are defined)
SAMPLES=/users/GR/mb/jquilez/data/chipseq/samples
ANALYSIS=$HOME/bioinfo_course/session03
JOB_CMD=$ANALYSIS/job_cmd 
JOB_OUT=$ANALYSIS/job_out
mkdir -p $JOB_CMD
mkdir -p $JOB_OUT
mkdir -p $ANALYSIS/data
mkdir -p $ANALYSIS/figures
compute_matrix=`which computeMatrix`
plot_heatmap=`which plotHeatmap`

# Cluster parameters
queue=short-sl7
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
	# comment the line above and uncomment the one below to perform the analysis on the entire dataset of peaks
	#ibeds="$ibeds $ibed"
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
#$ -M $email
#$ -m abe
#$ -pe smp $slots" > $job_file

# deeptools commands: compute matrix
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
			--perGroup \
			--refPointLabel='Peak' \
			--regionsLabel $peaks_names \
			--legendLocation=best \
			--samplesLabel $profiles_names"
echo $job_cmd >> $job_file

# remove intermediate files
echo "rm -f $ibeds" >> $job_file

# Submit job
chmod a+x $job_file
qsub < $job_file
#cat $job_file

