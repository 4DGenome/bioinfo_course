#!/usr/bin/env Rscript

# dependencies

library("tidyverse")

# settings

options(stringsAsFactors = F)

# list files

difexp_file <- list.files("data",
                     pattern = "deseq2",
                     full = TRUE)

tss_file <- list.files("analysis/annotation",
                       pattern = "tss",
                       full = TRUE)

# read data

difexp <- read.delim(difexp_file)

tss <- read.delim(tss_file, header = F)
names(tss) <- c("chr", "start", "end", "id", "score", "strand")

# describe dataset

summary(difexp)

str(difexp)

head(difexp)

str(tss)

# filter rows

difexp <- filter(difexp, !is.na(padj))

# select columns

tss <- select(tss, chr, id)

# merge info

difexp <- inner_join(difexp, tss)

# fist plot (volcano plot)

ggplot(difexp,
       aes(x = log2FoldChange, y = -log10(padj))) +
    geom_point(alpha = .3)

# create / transform columns

difexp <- mutate(difexp,
                 log2FoldChange = - log2FoldChange,
                 stat = - stat,
                 direction = ifelse(log2FoldChange > 0, "up", "down"),
                 stat_diff = ifelse(padj < .01, "signf", "no-sig"),
                 bio_diff = ifelse(abs(log2FoldChange) > 1, "relevant", "irrelevant"),
                 diff_group = paste(stat_diff, bio_diff),
                 change = ifelse((stat_diff == "signf") & (bio_diff == "relevant"),
                                 "change", "no change"))


# improve plot

ggplot(difexp,
       aes(x = log2FoldChange, y = -log10(padj),
           col = diff_group)) +
    geom_point(alpha = .3) +
    ylim(c(0, 50))

# summarize

summarize(difexp,
          n_up = sum(direction == "up"),
          n_tot = n(),
          p_up = n_up / n_tot,
          se_p_up = sqrt(p_up * (1 - p_up) / n()),
          ll = p_up - 2 * se_p_up,
          ul = p_up + 2 * se_p_up)

# filter and summarize

filter(difexp, chr == "chr13") %>%
    summarize(n_up = sum(direction == "up"),
              n_tot = n(),
              p_up = n_up / n_tot,
              se_p_up = sqrt(p_up * (1 - p_up) / n()),
              ll = p_up - 2 * se_p_up,
              ul = p_up + 2 * se_p_up)


# group and summarize

group_by(difexp, chr) %>%
    summarize(n_up = sum(direction == "up"),
              n_tot = n(),
              p_up = n_up / n_tot,
              se_p_up = sqrt(p_up * (1 - p_up) / n()),
              ll = p_up - 2 * se_p_up,
              ul = p_up + 2 * se_p_up)

group_by(difexp, change) %>%
    summarize(n_up = sum(direction == "up"),
              n_tot = n(),
              p_up = n_up / n_tot,
              se_p_up = sqrt(p_up * (1 - p_up) / n()),
              ll = p_up - 2 * se_p_up,
              ul = p_up + 2 * se_p_up)


# add distance to closest peak

closest_files <- list.files("analysis/closest",
                        pattern = ".closest_tss.bed",
                        full = TRUE)

closest_list <- lapply(closest_files, read.delim, header = F)

lapply(closest_list, "[[", "V10") %>% sapply(mean, na.rm = T) %>% rank

closest <- closest_list[[6]]

names(closest) <- c("chr_tss", "start_tss", "end_tss", "id",
                    "score_tss", "strand_tss",
                    "chr_peak", "start_peak", "end_peak", "dis")


# create new variable

closest <- select(closest, id, dis) %>%
    mutate(disKbp = dis / 1e3)

# merge information

dat <- inner_join(difexp, closest)

# explore distances

ggplot(dat,
       aes(x = disKbp)) +
    geom_histogram(bins = 100) +
    xlim(c(0, 500))

# mark those genes with a peak in the promoter and summarize

mutate(dat,
       overlap = ifelse(dis < 5e3, "peak", "no peak")) %>%
    group_by(direction, overlap) %>%
    summarize(n_ch = sum(change == "change"),
              n_tot = n(),
              p_ch = n_ch / n_tot,
              se_p_ch = sqrt(p_ch * (1 - p_ch) / n()),
              ll = p_ch - 2 * se_p_ch,
              ul = p_ch + 2 * se_p_ch) %>%
    na.omit

# auxiliary matrix to get Wald CI 95 % from a lm/glm output

alpha <- .05
wald_mat <- matrix(c(1, 0, 1, qnorm(alpha / 2), 1, qnorm(1 - alpha / 2)), 2)
colnames(wald_mat) <- c("est", "ll", "ul")

# fit a logistic model to check if the probability of a gene to be regulated
# upon hormone treatment is associated with the presence of PRBS in the promoter

mutate(dat,
       overlap = ifelse(dis < 5e3, 1, 0),
       up = ifelse(direction == "up", 1, 0)) %>%
    glm(change == "change" ~ up * overlap, binomial, .) %>%
    summary %$%
    coefficients %>%
    (function(x) x[,1:2] %*% wald_mat) %>%
    exp
