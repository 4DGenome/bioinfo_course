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
In order to start our practice we need to load the data in R. In particular, we will load the gene counts per sample and the metadata that describes the sample groups. Also, the gene biotype annotation was loaded, in order to work only with protein coding genes.

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
Exploratory analysis is first step in a real analysis. It provides a first sight of the data, and serves to detect some biases that potentially could modify the results of our analysis.

In this case, we will start by defining colors associated to each group, and then we will apply a reshape to our data, ir order to be compatible with [ggplot](http://ggplot2.tidyverse.org/index.html) package, tipically used to plot beatiful charts.

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
PCA (Principal Component Analysis) provides a powerful but intuitive description of our samples. It is defined as an orthogonal rotation in order to match the main variability directions in our data, and how our original variables (genes) contribute to them. In a general case, we will expect a good separation between our groups in the space defined by the first and second principal components, otherwise, higher sources of variability could be present in our data

```R
mod <- prcomp(counts)
explain_variance <- round((mod$sdev^2/sum(mod$sdev^2))*100,digits=2)

plot(mod$rotation[,1], mod$rotation[,2], xlab=paste("PC1 (",explain_variance[1],"%)") , ylab=paste("PC2 (",explain_variance[2],")"), type="n")
text(mod$rotation[,1], mod$rotation[,2], labels=rownames(mod$rotation), col=sample_colors)
legend("bottomleft",legend=names(group_colors),col=group_colors,lwd=2, cex=0.7)
```

# Differential expression

Usually, when are interested in comparing different conditions to see what genes
change their expression. There are several ways to do so. Here we are going to
use
[DESeq2](https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html),
one of the most popular packages to perform RNAseq differential expression.

As our experiment presents three conditions, we can define two different
comparisons: one for each hormone (compared to untreated cells).

## P60 vs T0

First we need to define a design `data.frame`. In this case, we can make use of
the "metadata" already loaded and subset the relevant rows. Remember the meaning
of the verbs `filter` and `mutate` explained in
the
[previous](https://github.com/4DGenome/bioinfo_course/tree/master/session02#basic-interaction-with) [session](https://github.com/4DGenome/bioinfo_course/tree/master/session02#transform-variables).

```R
design_t0_r6 <- filter(design, group %in% c("T47D R6", "T47D T0")) %>%
    mutate(group = relevel(factor(group), "T47D T0"))
           
rownames(design_t0_r6) <- design_t0_r6$sample
```
Some tips to better understand the code:

* `%in%` takes the left hand side vector and tells if each of its elements is
  contained **in** the right hand side one.
* `factor` is function to tell `R` that the variable should be treated as
  categorical. This is important if we want to use the variable to define the
  groups of the comparison.
* `relevel` is a useful function that allows to define the *reference* group
  of a factor. In this case, we want the untreated to be the reference.

It is important that the design object is annotated with the row names matching
the column names of the matrix of counts.

Next thing we need is to up an *"container"* object with the requirements of the
`DESeq2` package. The easiest way to do it is using the `DESeqDataSetFromMatrix` function, that takes
three arguments:

```R
dds_t0_r6 <- DESeqDataSetFromMatrix(countData = counts[,design_t0_r6$sample],
                                    colData = design_t0_r6,
                                    design = ~ group)

dds_t0_r6 <- DESeq(dds_t0_r6)
```
* matrix of counts
* design `data.frame`
* factor defining the comparison (of the design `data.frame`). The `~` symbol is
  required to let `R` know that we are specifying a comparison model.

> Do you know what is doing the `[,]` operator?
>
> It allows us to subset the `[rows, columns]` of a matrix. It takes a vector of
> indexes as arguments (the sample IDs in this case). Whenever one of the `rows`
> or `columns` index vector is empty, it returns all of its elements.

And then feed the `DESeq` function with it (where all the magic happens). The
result is a quite complicated object with a lot of information about the
comparison process. However we can extract the relevant information using the
`results` function. I'm pretty sure you can guess what the rest of this code
chuck does (basically, coerce the output to a `data.frame` and create a new
column from the row names).

```R
res_t0_r6 <- results(dds_t0_r6) %>%
    as.data.frame
res_t0_r6$gene <- rownames(res_t0_r6)
```

Nice, we have the results. Do you want to have a look?

```R
head(res_t0_r6)
```

* **baseMean** is the average of the expression for **all** samples
  (irrespective of their treatment).
* **log2FoldChange** well, it's the log2 Fold Change between the treated and the
  untreated samples.
* **lfcSE** is the standard error of the log2 Fold Change.
* **stat** is the statistical score used in the test (a moderated t-statistic).
* **pvalue** is the raw p-value obtained from the comparison.
* **padj** is the adjusted p-value. We need to modify the raw p-value in order
  to avoid a big false positive rate due to multiple comparisons (we are
  performing as many tests as genes in the data set!)
*  **gene** is the gen ID

In most experiments there is a set of genes that present expression levels so
low that cannot be reliably quantified. Those present a `NA` in the `padj`
column. It could be use full to filter them out.

```R
res_t0_e6 <- filter(res_t0_e6, !is.na(padj))
```
Once we have the numerical results, let's simplify to make them more
understandable. We can use the same approach as in
the
[last session](https://github.com/4DGenome/bioinfo_course/tree/master/session02#transform-variables) to
decide what genes present biological and statistical differences between conditions.

```R
res_t0_e6 <- mutate(res_t0_e6,
                    direction = ifelse(log2FoldChange > 0, "up", "down"),
                    stat_diff = ifelse(padj < .01, "signf", "no-sig"),
                    bio_diff = ifelse(abs(log2FoldChange) > 1, "relevant", "irrelevant"),
                    diff_group = paste(stat_diff, bio_diff),
                    change = ifelse((stat_diff == "signf") & (bio_diff == "relevant"),
                                    "change", "no change"))
```

Remember that one can `group_by` some variables and `summarize` the result. For
instance, counting the number of differentially expressed genes on each
"direction". Or we can store a vector with the gene IDs of the genes going up or down.

```R
group_by(res_t0_e6, change, direction) %>%
    summarize(n = n())

genes_t0_r6_up <- filter(res_t0_r6, change == "change", direction == "up")$gene
genes_t0_r6_down <- filter(res_t0_r6, change == "change", direction == "down")$gene
```

## E60 vs T0

Once we have performed one differential analysis, we can repeat it for the other
hormone. We just need to change the subseting of both the design `data.frame`
and the matrix of counts.

```R
design_t0_e6 <- filter(design, group %in% c("T47D E6", "T47D T0")) %>%
    mutate(group = relevel(factor(group), "T47D T0"))
           
rownames(design_t0_e6) <- design_t0_e6$sample
```
The rest is exactly the same (bu for the object names, because we don't want to lose all the
previous analysis).

```R
dds_t0_e6 <- DESeqDataSetFromMatrix(countData = counts[,design_t0_e6$sample],
                                    colData = design_t0_e6,
                                    design = ~ group)

dds_t0_e6 <- DESeq(dds_t0_e6)

res_t0_e6 <- results(dds_t0_e6) %>%
    as.data.frame
res_t0_e6$gene <- rownames(res_t0_e6)

res_t0_e6 <- filter(res_t0_e6, !is.na(padj))

res_t0_e6 <- mutate(res_t0_e6,
                    direction = ifelse(log2FoldChange > 0, "up", "down"),
                    stat_diff = ifelse(padj < .01, "signf", "no-sig"),
                    bio_diff = ifelse(abs(log2FoldChange) > 1, "relevant", "irrelevant"),
                    diff_group = paste(stat_diff, bio_diff),
                    change = ifelse((stat_diff == "signf") & (bio_diff == "relevant"),
                                    "change", "no change"))
```

How many genes are changing upon estradiol treatment?

```R
group_by(res_t0_e6, change, direction) %>%
    summarize(n = n())

genes_t0_e6_up <- filter(res_t0_e6, change == "change", direction == "up")$gene
genes_t0_e6_down <- filter(res_t0_e6, change == "change", direction == "down")$gene
```

# Volcano plot

A common way to visually summarize a differential expression experiment is via
a volcano plot. It depicts both "dimensions" of the differences: the **magnitude**
of the change (encoded in the Fold Change on the X axis) and its **variability** (encoded in
the p-value on the Y axis). Each point represents one gene.

We can make one of such plots with the following code:

```R
ggplot(res_t0_r6,
       aes(x = log2FoldChange,
           y = -log10(padj))) +
    geom_point(alpha = .1) +
    geom_hline(yintercept = -log10(.01), col = 2, lty = 2) +
    geom_vline(xintercept = c(-1, 1), col = 2, lty = 2) +
    coord_cartesian(ylim = c(0, 50))
```

where the `ggplot` function initializes the plot taking the `data.frame` of the
results and the definition of the variables we want on each axis and then the
`geom_*` functions draw things in the plot area. We have included some
transparency (`alpha`) for the points to get a feeling of the density of the
genes. `_vline` and `_hline` draws vertical and horizontal lines marking the
thresholds selected.

The plot for the other treatment is equivalent:

```R
ggplot(res_t0_e6,
       aes(x = log2FoldChange,
           y = -log10(padj))) +
    geom_point(alpha = .1) +
    geom_hline(yintercept = -log10(.01), col = 2, lty = 2) +
    geom_vline(xintercept = c(-1, 1), col = 2, lty = 2) +
    coord_cartesian(ylim = c(0, 50))
```

> BONUS: Can we put both plots in the same figure, with the same scale?

Sure! But we need some preliminary work to prepare a `data.frame` with the
relevant info (including a new variable `comparison` that labels the
comparison). Please take note of the faceting obtained with the `facet_wrap`
function.


```R
aux_res_t0_r6 <- mutate(res_t0_r6,
                        comparison = "P6 vs T0") %>%
    dplyr::select(gene, log2FoldChange, padj, comparison)

aux_res_t0_e6 <- mutate(res_t0_e6,
                        comparison = "E6 vs T0") %>%
    dplyr::select(gene, log2FoldChange, padj, comparison)

aux_res <- rbind(aux_res_t0_r6, aux_res_t0_e6)

ggplot(aux_res,
       aes(x = log2FoldChange,
           y = -log10(padj))) +
        geom_point(alpha = .1) +
    geom_hline(yintercept = -log10(.01), col = 2, lty = 2) +
    geom_vline(xintercept = c(-1, 1), col = 2, lty = 2) +
    coord_cartesian(ylim = c(0, 50)) +
    facet_wrap(~ comparison, ncol = 1)

```


# Venn diagrams

Volcano plots are fine, but sometimes we want something simpler that summarizes
other aspects of the analysis. If we are interested in what genes are regulated
by both hormones, we can draw a Venn diagram.

There are several packages to do that. Here we are using one of them
(`VennDiagram`), but there are many others.

First we have to define a object via the `venn.diagram` function. The most
important arguments are the first one (one list with genes up regulated in each
condition) and the second one, that have to be set to `NULL` if we want to
display the result in the screen. The rest of the arguments define some
graphical properties of the plot. You can have a look at all of them in the
help page `?venn.diagram`.

```R
venn_up <- venn.diagram(list(R6 = genes_t0_r6_up,
                             E6 = genes_t0_e6_up),
                          NULL,
                          fill = c("blue", "red"),
                          alpha = c(0.5,0.5),
                          cex = 2,
                          main = "Upregualted genes after hormone",
                          main.cex = 2)

grid.newpage()
grid.draw(venn_up)
```
After defining the object, the last two lines draw it in a new grid (i.e. window).

Of course the same can be done for the down regulated genes.

```R
venn_down <- venn.diagram(list(R6 = genes_t0_r6_down,
                               E6 = genes_t0_e6_down),
                          NULL,
                          fill = c("blue", "red"),
                          alpha = c(0.5,0.5),
                          cex = 2,
                          main = "Downregualted genes after hormone",
                          main.cex = 2)

grid.newpage()
grid.draw(venn_down)
```

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

