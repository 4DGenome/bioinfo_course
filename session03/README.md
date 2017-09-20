# Overview

<br>

# Access to the cluster

In [session01](https://github.com/4DGenome/bioinfo_course/tree/master/session01) and [session02](https://github.com/4DGenome/bioinfo_course/tree/master/session02) we worked from the command line in the [CRG cluster](http://www.linux.crg.es/index.php/Main_Page).

As a reminder, the CRG cluster is accessed with:
```bash
ssh -Y USER@ant-login.linux.crg.es
```
Where your USER should be provided by the system administrators of the CRG cluster and has the first letter of your name plus your (first) lastname, all lower-case (e.g. `jquilez` for Javier Quilez).

The command above brings us to our individual home directory in the cluster.

Most of the present session will be executed in [R](https://www.r-project.org/), which, as seen in [session02](https://github.com/4DGenome/bioinfo_course/tree/master/session02), can be run in the command line by typing:
```bash
R
```

[RStudio](https://www.rstudio.com/) offers a more user-friendly interface to R, and the CRG cluster hosts a RStudio server:
```
http://rstudio.linux.crg.es/
```
You will be asked for your USER (see above) and your password. By accessing to the RStudio server we will be working in a node of the CRG cluster from our home directory. You can get your home directory with:
```R
getwd()
```

<br>

# Load data

# Exploratory analysis

## Boxplots

##  PCA

# Differential expression

## P60 vs T0

## E60 vs T0

# Volcano plot

# Venn diagrams

<br>

# \*-seq signal over regions

Very often we want to know if a given protein is binding over a set of genomic regions. We can measure this by looking at the accumulation of sequencing reads of a ChIP-seq experiment targeting that protein over such regions. More generally, this analysis can be extended to other \*-seq experiments.

For this analysis we need at least two files:
1. read counts profiles (preferably normalized by the number of reads so that different profiles from different samples are comparable)
2. set of regions of interest

For (1) we will use the read per million profiles for untreated and treated (60 minutes after R5020; aka. T60) PR ChIP-seq samples:
```
# untreated
ibw1=/users/GR/mb/jquilez/data/chipseq/samples/gv_009_02_01_chipseq/profiles/hg38_mmtv/single_end/gv_009_02_01_chipseq.rpm.bw
# T60
ibw2=/users/GR/mb/jquilez/data/chipseq/samples/gv_066_01_01_chipseq/profiles/hg38_mmtv/single_end/gv_066_01_01_chipseq.rpm.bw
```

For (2) we will use the binding regions identified in the T60 samples:
```
ibed=/users/project/4DGenome/bioinfo_course/session03/tmp_gv_066_01_01_chipseq.bed
```

To generate the plots we will use [deepTools](http://deeptools.readthedocs.io/en/latest/content/tools/computeMatrix.html), a suite of tools for exploring high-throughput sequencing data.

First, we will use [computeMatrix](http://deeptools.readthedocs.io/en/latest/content/tools/computeMatrix.html) to generate a matrix in which:
- each row is a region in (2)
- each region is extended _X_ bp upstream and downstream, and the extended region is partitioned in bins of size _S_ bp, each bin resulting in a column of the matrix
- for each cell of the matrix, the signal from (1) is averaged (mean by default)

The code:
```bash
# output file (it is used later as input file!)
omat=/users/project/4DGenome/bioinfo_course/session03/data/average_profiles_binsize10.mat
# see computeMatrix options for the other parameters
/software/mb/el7.2/anaconda2/bin/computeMatrix reference-point -S $ibw1 $ibw2 -R $ibed -out $omat --referencePoint=center --binSize=10 --upstream=1000 --downstream=1000 --numberOfProcessors=10
```

Second, we plot:
- the resulting matrix as a heatmap (the higher the average number of normalized reads per bin per region, the more intense the color in the used scale)
- an average profile that averages the values for all regions for each bin.

The code:
```bash
# output file
opdf=/users/project/4DGenome/bioinfo_course/session03/figures/average_profiles_binsize10.pdf
/software/mb/el7.2/anaconda2/bin/plotHeatmap -m $omat -out $opdf --plotFileFormat pdf --heatmapHeight 10 --heatmapWidth 8 --xAxisLabel='' --startLabel='' --endLabel='' --colorMap=Blues --perGroup --refPointLabel='Peak' --regionsLabel T60-peaks --legendLocation=best --samplesLabel Untreated T60

# remove intermediate files
rm -f  /users/project/4DGenome/bioinfo_course/session03/tmp_gv_066_01_01_chipseq.bed
```

