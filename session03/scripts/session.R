
library(ggplot2)
library(reshape2)
library(dplyr)
library(DESeq2)
library(VennDiagram)

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

# melt data for GGPLOT
melt_counts <- melt(counts,varnames = c("gene","sample"))
melt_counts$group <- design[match(melt_counts$sample,design$sample),"group"]

# boxplot
p <- ggplot(melt_counts, aes(x=sample, y=value, fill=group)) + theme(axis.text.x=element_text(angle=-45, hjust=0))
p + geom_boxplot(outlier.alpha = .1)
p + geom_boxplot(outlier.alpha = .1) + scale_y_log10()

# PCA
mod <- prcomp(counts)
explain_variance <- round((mod$sdev^2/sum(mod$sdev^2))*100,digits=2)

plot(mod$rotation[,1], mod$rotation[,2], xlab=paste("PC1 (",explain_variance[1],"%)") , ylab=paste("PC2 (",explain_variance[2],"%)"), type="n")
text(mod$rotation[,1], mod$rotation[,2], labels=rownames(mod$rotation), col=sample_colors)
legend("bottomleft",legend=names(group_colors),col=group_colors,lwd=2, cex=0.7)

### DIFFERENTIAL ANALYSIS

# R6 vs T0

design_t0_r6 <- filter(design, group %in% c("T47D R6", "T47D T0")) %>%
    mutate(group = relevel(factor(group), "T47D T0"))
           
rownames(design_t0_r6) <- design_t0_r6$sample

dds_t0_r6 <- DESeqDataSetFromMatrix(countData = counts[,design_t0_r6$sample],
                                    colData = design_t0_r6,
                                    design = ~ group)

dds_t0_r6 <- DESeq(dds_t0_r6)

res_t0_r6 <- results(dds_t0_r6) %>%
    as.data.frame
res_t0_r6$gene <- rownames(res_t0_r6)

head(res_t0_r6)

res_t0_r6 <- filter(res_t0_r6, !is.na(padj))

res_t0_r6 <- mutate(res_t0_r6,
                    direction = ifelse(log2FoldChange > 0, "up", "down"),
                    stat_diff = ifelse(padj < .01, "signf", "no-sig"),
                    bio_diff = ifelse(abs(log2FoldChange) > 1, "relevant", "irrelevant"),
                    diff_group = paste(stat_diff, bio_diff),
                    change = ifelse((stat_diff == "signf") & (bio_diff == "relevant"),
                                    "change", "no change"))

group_by(res_t0_r6, change, direction) %>%
    summarize(n = n())

genes_t0_r6_up <- filter(res_t0_r6, change == "change", direction == "up")$gene
genes_t0_r6_down <- filter(res_t0_r6, change == "change", direction == "down")$gene


# E6 vs T0

design_t0_e6 <- filter(design, group %in% c("T47D E6", "T47D T0")) %>%
    mutate(group = relevel(factor(group), "T47D T0"))
           
rownames(design_t0_e6) <- design_t0_e6$sample

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

group_by(res_t0_e6, change, direction) %>%
    summarize(n = n())

genes_t0_e6_up <- filter(res_t0_e6, change == "change", direction == "up")$gene
genes_t0_e6_down <- filter(res_t0_e6, change == "change", direction == "down")$gene


### VOLCANO PLOT

ggplot(res_t0_r6,
       aes(x = log2FoldChange,
           y = -log10(padj))) +
    geom_point(alpha = .1) +
    geom_hline(yintercept = -log10(.01), col = 2, lty = 2) +
    geom_vline(xintercept = c(-1, 1), col = 2, lty = 2) +
    coord_cartesian(ylim = c(0, 50))

ggplot(res_t0_e6,
       aes(x = log2FoldChange,
           y = -log10(padj))) +
    geom_point(alpha = .1) +
    geom_hline(yintercept = -log10(.01), col = 2, lty = 2) +
    geom_vline(xintercept = c(-1, 1), col = 2, lty = 2) +
    coord_cartesian(ylim = c(0, 50))

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

### VENN DIAGRAMS

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



