
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


### FILTERING

# load biotype
biotype <- read.table("../data/gene_biotype_from_biomart.txt", header=T, sep="\t", stringsAsFactors = F)

# get biotype
av_genes_biotype <- biotype[match(rownames(counts),biotype[,1]),2]
table(av_genes_biotype)
table(is.na(av_genes_biotype))

# filter
counts <- counts[(!is.na(av_genes_biotype) & av_genes_biotype=="protein_coding"),]


### DESCRIPTIVE PLOTS

group_colors <- c("#e41a1c","#377eb8","#4daf4a")
names(group_colors) <- unique(design$group)

sample_colors <- group_colors[design$group]

# boxplot
melt_counts <- melt(counts,varnames = c("gene","sample"))
melt_counts$group <- design[melt_counts$sample,"group"]

p <- ggplot(melt_counts, aes(x=sample, y=value, fill=group)) + theme(axis.text.x=element_text(angle=-45, hjust=0))
p + geom_boxplot()
p + geom_boxplot() + scale_y_log10()

# PCA
mod <- prcomp(counts)
explain_variance <- round((mod$sdev^2/sum(mod$sdev^2))*100,digits=2)

plot(mod$rotation[,1], mod$rotation[,2], xlab=paste("PC1 (",explain_variance[1],"%)") , ylab=paste("PC2 (",explain_variance[2],")"), type="n")
text(mod$rotation[,1], mod$rotation[,2], labels=rownames(mod$rotation), col=sample_colors)
legend("bottomleft",legend=names(group_colors),col=group_colors,lwd=2, cex=0.7)
