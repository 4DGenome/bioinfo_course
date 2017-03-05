#!/usr/bin/env Rscript

# dependencies

library("tidyverse")
#library("VennDiagram")

# settings

options(stringsAsFactors = F)

# list files

difexp_file <- list.files("data",
                     pattern = "deseq2",
                     full = TRUE)

# read data

difexp <- read.delim(difexp_file)

tss <- read.delim(tss_file, header = F)

names(tss) <- c("chr", "tss", "tss2", "id", "score", "strand")

# subset columns

tss <- select(tss, chr, id)

# describe dataset

summary(difexp)

str(difexp)

head(difexp)

head(tss)

# merge datasets

difexp <- inner_join(difexp, tss)

str(difexp)

# filter rows

difexp <- filter(difexp, !is.na(padj))

# fist plot (volcano plot)

ggplot(difexp,
       aes(x = log2FoldChange, y = -log10(padj))) +
    geom_point(alpha = .3)

# create / transform columns

difexp <- mutate(difexp,
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

lapply(closest_list, "[[", "V17") %>% sapply(mean, na.rm = T) %>% rank

closest <- closest_list[[5]]

names(closest) <- c("chr_tss", "start_tss", "end_tss", "id",
                    "score_tss", "strand_tss",
                    "chr_peak", "start_peak", "end_peak",
                    "id_peak", "score_peak", "strand_peak",
                    "signal_peak", "pval_peak", "qval_peak",
                    "summit_peak", "dis")

closest <- select(closest, id, dis) %>%
    mutate(disMbp = dis / 1e6)

dat <- inner_join(difexp, closest)

# explore distances

ggplot(dat,
       aes(x = dis)) +
    geom_histogram(bins = 100) +
    xlim(c(0, .5e6))

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

ps <- mutate(dat,
       overlap = ifelse(dis < 5e3, "peak", "no peak")) %>%
    group_by(direction, overlap) %>%
    summarize(n_ch = sum(change == "change"),
              n_tot = n(),
              p_ch = n_ch / n_tot,
              se_p_ch = sqrt(p_ch * (1 - p_ch) / n()),
              ll = p_ch - 2 * se_p_ch,
              ul = p_ch + 2 * se_p_ch) %>%
    na.omit %>%
    with(setNames(p_ch, paste(direction, overlap)))

ps / min(ps)

mutate(rel_ch = p_ch / p_ch[(direction == "up") & (overlap == "no peak")])


alpha <- .05
wald_mat <- matrix(c(1, 0, 1, qnorm(alpha / 2), 1, qnorm(1 - alpha / 2)), 2)
colnames(wald_mat) <- c("est", "ll", "ul")

mutate(dat,
       overlap = ifelse(dis < 10e3, 1, 0),
       down = ifelse(direction == "down", 1, 0)) %>%
    glm(change == "change" ~ down * overlap, binomial, .) %>%
    summary %$%
    coefficients %>%
    (function(x) x[,1:2] %*% wald_mat) %>%
    exp
