# Overview

This session focuses on how to perform some of the analyses using RNA-seq and ChIP-seq data that we are requested more often.
1. Refresh how to access to the cluster and RStudio server
2. Load RNA-seq data into R
3. Exploratory analysis: boxplots and PCA
4. Differential expression analysis (DEA)
5. Interpretation of DEA results: volcano plots and venn diagrams
6. ChIP-seq signal over regions

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

# Load RNA-Seq data
In order to start our practice we need to load the data in R. In particular, we will load the gene counts per sample, the metadata that describes the sample groups, and the gene biotype annotation, in order to work only with protein coding genes.

```R

library(ggplot2)
library(reshape2)

### LOAD DATA

# metadata
design <- read.table("../data/design.txt", header=F, sep="\t", stringsAsFactors=F)
colnames(design) <- c("sample","group")

# counts
count_list <- list()
for(i in 1:nrow(design)){
  sample_file <- list.files("../data", pattern=design$sample[i], full.names = T)
  count_list[[i]] <- read.table(sample_file, header=T, sep="\t", stringsAsFactors=F)
}
counts <- sapply(count_list,"[[",7)
colnames(counts) <- design$sample
rownames(counts) <- sapply(strsplit(count_list[[1]]$Geneid, "\\."), "[[", 1)


### FILTERING NON CODING GENES

# load biotype
biotype <- read.table("../data/gene_biotype_from_biomart.txt", header=T, sep="\t", stringsAsFactors = F)

# get biotype
av_genes_biotype <- biotype[match(rownames(counts),biotype[,1]),2]
table(av_genes_biotype)
table(is.na(av_genes_biotype))

# filter
counts <- counts[(!is.na(av_genes_biotype) & av_genes_biotype=="protein_coding"),]
```

# Exploratory analysis

```R
group_colors <- c("#e41a1c","#377eb8","#4daf4a")
names(group_colors) <- unique(design$group)

sample_colors <- group_colors[design$group]

# melt data for GGPLOT
melt_counts <- melt(counts,varnames = c("gene","sample"))
melt_counts$group <- design[melt_counts$sample,"group"]

```

## Boxplots
Boxplots can summarize in a very simple way the range of our samples. In particular, the main quartiles of sample distributions and outliers are usually plotted. Big differences between replicates are indicative of batch effect biases or normalization problems.

```R
p <- ggplot(melt_counts, aes(x=sample, y=value, fill=group)) + theme(axis.text.x=element_text(angle=-45, hjust=0))
p + geom_boxplot()
p + geom_boxplot() + scale_y_log10()
```

## PCA

```R
mod <- prcomp(counts)
explain_variance <- round((mod$sdev^2/sum(mod$sdev^2))*100,digits=2)

plot(mod$rotation[,1], mod$rotation[,2], xlab=paste("PC1 (",explain_variance[1],"%)") , ylab=paste("PC2 (",explain_variance[2],")"), type="n")
text(mod$rotation[,1], mod$rotation[,2], labels=rownames(mod$rotation), col=sample_colors)
legend("bottomleft",legend=names(group_colors),col=group_colors,lwd=2, cex=0.7)
```

# Differential expression

## P60 vs T0

## E60 vs T0

# Volcano plot

# Venn diagrams

<br>

# ChIP-seq signal over regions

Very often we want to know if a given protein is binding over a set of genomic regions. We can measure this by looking at the accumulation of sequencing reads of a ChIP-seq experiment targeting that protein over such regions. More generally, this analysis can be extended to other \*-seq experiments.

For this analysis we need at least two files:
1. read counts profiles (preferably normalized by the number of reads so that different profiles from different samples are comparable)
2. set of regions of interest

For (1) we will use the read per million profiles for untreated and treated (60 minutes after R5020; aka. T60) PR ChIP-seq samples in T47D cells:
```
# untreated
ibw1=/users/GR/mb/jquilez/data/chipseq/samples/gv_009_02_01_chipseq/profiles/hg38_mmtv/single_end/gv_009_02_01_chipseq.rpm.bw
# T60
ibw2=/users/GR/mb/jquilez/data/chipseq/samples/gv_066_01_01_chipseq/profiles/hg38_mmtv/single_end/gv_066_01_01_chipseq.rpm.bw
```

For (2) we will use the binding regions identified in the T60 samples:
```
ibed=/users/GR/mb/jquilez/data/chipseq/samples/gv_066_01_01_chipseq/peaks/macs2/hg38_mmtv/with_control/single_end/gv_066_01_01_chipseq_peaks.narrowPeak
```

This file contains ~80,000 peaks:
```
wc -l $ibed
```
The computation time of the analysis notably scales up with the number of regions. Therefore, for the sake of speed, here we will reduce the analysis to a subset of 1,000 random peak regions. The random subset of peaks can be generated with:
```bash
# make directory in your home directory to store the output of the analysis
ANALYSIS=$HOME/bioinfo_course/session03
mkdir -p $ANALYSIS/{data,figures}
tbed=$ANALYSIS/data/tmp_gv_066_01_01_chipseq.bed
shuf -n 1000 $ibed --random-source=$ibed > $tbed
```
The `--random-source` option guarantees that the same set of _random_ peaks is selected every time the command is executed and thus ensures the reproducibility of the analysis; indeed the [random selection is not as random as we could think](https://www.random.org/randomness/).

To generate the plots we will use [deepTools](http://deeptools.readthedocs.io/en/latest/index.html), a suite of tools for exploring high-throughput sequencing data.

First, we will use [computeMatrix](http://deeptools.readthedocs.io/en/latest/content/tools/computeMatrix.html) to generate a matrix in which:
- each row is a region in (2)
- each region is extended _X_ bp upstream and downstream, and the extended region is partitioned in bins of size _S_ bp, each bin resulting in a column of the matrix
- for each cell of the matrix, the signal from (1) is averaged (mean by default)

The code:
```bash
# output file (it is used later as input file!)
omat=$ANALYSIS/data/average_profiles_binsize10.mat
# see computeMatrix options for the other parameters
/software/mb/el7.2/anaconda2/bin/computeMatrix reference-point -S $ibw1 $ibw2 -R $tbed -out $omat --referencePoint=center --binSize=10 --upstream=1000 --downstream=1000 --numberOfProcessors=10
```

Second, we use [plotHeatmap](http://deeptools.readthedocs.io/en/latest/content/tools/plotHeatmap.html) to plot:
- the resulting matrix as a heatmap (the higher the average number of normalized reads per bin per region, the more intense the color in the used scale)
- an average profile that averages the values in all regions for each bin.

The code:
```bash
# output file
opdf=$ANALYSIS/figures/average_profiles_binsize10.pdf
/software/mb/el7.2/anaconda2/bin/plotHeatmap -m $omat -out $opdf --plotFileFormat pdf --heatmapHeight 10 --heatmapWidth 8 --xAxisLabel='' --startLabel='' --endLabel='' --colorMap=Blues --perGroup --refPointLabel='Peak' --regionsLabel T60-peaks --legendLocation=best --samplesLabel Untreated T60
```

We remove the intermediate files:
```bash
rm -f $tbed
```

## Bonus: submitting a job to the cluster

Above we run the code in a login node of the CRG cluster. As warned when entering the cluster, we are not supposed to run heavy computation in login nodes (as it would be the case for much larger set of regions or multiple ChIP-seq signal profiles). In such cases, you are encouraged to:
- prepare a job script with the code to be run
- submit the job to a computing node in the cluster

As an example, this can be done with (`cat` prints it but does not execute it):
```
cat /users/project/4DGenome/bioinfo_course/session03/scripts/average_profiles.sh
```
A more user-friendly view of the [same script](https://github.com/4DGenome/bioinfo_course/blob/master/session03/scripts/average_profiles.sh).

Edit the script in the 'variables' section to change things like the email to which notification about the job are sent, number of random regions, etc.

And finally... execute it!
```
/users/project/4DGenome/bioinfo_course/session03/scripts/average_profiles.sh
```

