
# load metadata
design <- read.table("../data/design.txt", header=F, sep="\t", stringsAsFactors=F)
colnames(design) <- c("sample","group")

# load counts
count_list <- list()
for(i in 1:nrow(design)){
  sample_file <- list.files("../data", pattern=design$sample[i], full.names = T)
  count_list[[i]] <- read.table(sample_file, header=T, sep="\t", stringsAsFactors=F)
}
counts <- sapply(count_list,"[[",7)
colnames(counts) <- design$sample
rownames(counts) <- count_list[[1]]$Geneid


