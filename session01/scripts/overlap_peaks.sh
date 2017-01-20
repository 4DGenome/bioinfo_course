#!/bin/bash

# lines starting with the "#" are not executed, which is useful to:
# - document code
# - test your code by running certain blocks only, that is, thos uncommented

# although it is not needed for executing the script, I like including a header as the one below
# at the beginning of my scripts to provide some basic information

#==================================================================================================
# Created on: 2017-01-20
# Usage: overlap_peaks.sh
# Author: Javier Quilez (GitHub: jaquol)
# Goal: calculate overlap between 2 lists of ChIP-seq peaks allowing both q-value and enrichment filtering
#==================================================================================================



#==================================================================================================
# (1) Locating data
#==================================================================================================

# a variable is defined with:
variable_name="variable_content"
# tips for defining variables
# - use only letters and numbers for the variable_name
# - do not use spaces for the variable_name
# - do not use a space between the variable_name and the variable_content

# a variable is referenced by preciding its name with "$". For instance, if we want to print its content
echo $variable_name

# Now we are going to define the paths to the input files, that is, the 2 *.narrowPeak files
ifile1="data/gv_009_02_01_chipseq_peaks.narrowPeak"
ifile2="data/gv_066_01_01_chipseq_peaks.narrowPeak"

# however, using absolute paths instead of relative paths is preferred; otherwise, the paths defined
# above will only be valid as long as this script is executed from the directory above `data`
# the abosulte path will be different for each of you because the path of the home directory is
# different for everybody. We can use the `$HOME` variable, which is a Unix built-in variable
# that stores the user's home directory path
ifile1="$HOME/bioinfo_course/session01/data/gv_009_02_01_chipseq_peaks.narrowPeak"
ifile2="$HOME/bioinfo_course/session01/data/gv_066_01_01_chipseq_peaks.narrowPeak"

# additionally, note that part of the paths above is shared, so we can clean the code a bit with:
DATA=$HOME/bioinfo_course/session01/data
ifile1="$DATA/gv_009_02_01_chipseq_peaks.narrowPeak"
ifile2="$DATA/gv_066_01_01_chipseq_peaks.narrowPeak"

# as a convention I typically use:
# - lower case for variables referencing to text strings or paths to files
# - upper case for variables referencing to directories



#==================================================================================================
# (2) Filtering peaks
#==================================================================================================

# as explained in the accompanying README file, one of the steps of a ChIP-seq pipeline is the identification of regions of signal enrichment (aka peaks or binding sites). For this dataset I called peaks with MACS2 (https://github.com/taoliu/MACS), which provides a list of peaks as narrowPeak files. As you may remember, in 2 important fields are shown in columns 7 and 9 of a narrowPeak:
# - 7th column: enrichment (over the input/control)
# - 9th column: significance of the peak, expressed as the -log10(FDR q-value)
# by default I run MACS so that it provides a list of all the peaks identified that have q-value < 0.05 (no filtering on the enrichment is applied). Such criteria are too permissive for any downstream analysis but I do it this way because I prefer having a comprehensive list to filter from than having to run MACS2 for different thresholds. That is why the next step we will perform is filtering our initial list of peaks based on both enrichment and significance. There are no consensus cutoff values and the empirical way in which I define them is out of the scope of this session. For the time being we will apply very stringent cutoff values:
# - enrichment > 5
# - q-value < 10e-10
# For your analyses I would recommend you applying 2 sets of cutoffs, one more permissive and one more stringent, and see if the results are robust

# cutoff values
min_enrichment=5
min_qval=10	# expressed as -log10(x)

# the output of the analyses of this session will be saved in the directory defined with:
ANALYSIS=$HOME/bioinfo_course/session01/analysis/overlap_peaks
# having a directory specific for the output of the analysis is not required to run the script. However, I think that following a structure as we are doing (e.g. data, scripts, analysis) is a good practice

# we will create 2 files to save the filtered lists of peaks
tfile1="$ANALYSIS/gv_009_02_01_chipseq_peaks_filtered.narrowPeak"
tfile2="$ANALYSIS/gv_066_01_01_chipseq_peaks_filtered.narrowPeak"

# and here comes the code for the filtering...
awk -v min_qval=$min_qval -v min_enrichment=$min_enrichment '($7 > min_enrichment) && ($9 > min_qval)' $ifile1 | grep -v "chrM\|chrY\|chrUn\|luciferase" > $tfile1
# many things are happening here so let's break it down conceptually:
# - open the file with the first list of peaks
# - subset those rows in which both enrichment and -log10(q-value) values are greated than the pre-defined cutoffs
# - exlcude peaks in chromosomes other than autosomes and chrX
# - print the surviving rows into a new file

# the synthax used introduces many novel Unix programming concepts. If you are interested in knowing the details of the synthax the code is explained in the README file.

# and he do the same for the other list of peaks
awk -v min_qval=$min_qval -v min_enrichment=$min_enrichment '($7 > min_enrichment) && ($9 > min_qval)' $ifile2 | grep -v "chrM\|chrY\|chrUn\|luciferase" > $tfile2

# by the way, we can count the number of filtered peaks in each file
wc -l $tfile1
wc -l $tfile2



#==================================================================================================
# (3) Overlap peaks
#==================================================================================================

# here we will compare the genomic coordinates of the peaks in the 2 filtered peaks lists to find what peaks are shared and what peaks are specific of each list/sample. As done before, we will first define a file to save the output of this bit of code
ofile=$ANALYSIS/overlap_peaks.bed
bedtools multiinter -i $tfile1 $tfile2 -header -names sample1 sample2 > $ofile

# the synthax used introduces the BEDtools program. If you are interested in knowing the details of the synthax the code is explained in the README file.



#==================================================================================================
# (4) Remove intermediate files
#==================================================================================================

# rm is another built-in Unix function to delete files and directories --use it with caution!!
rm $tfile1 $tfile2