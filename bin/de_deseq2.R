#! /usr/bin/env Rscript

"Perform DE analysis using DESeq2

Usage:
  de_deseq2.R [--nofilter --output-dir <PATH>] --b1 <PATH> --b2 <PATH> <txi>

Options:
  --nofilter           : Disable filter [defalt: FALSE]
  --output-dir <PATH>  : Output directory [default: .]
  --b1 <PATH>          : Control samples list; same as rMATs
  --b2 <PATH>          : Case samples list; same as rMATs
  <txi>                : Count matrix; tximport object RDS

" -> doc


# %% Prepare
library(DESeq2)
library(tidyverse)

options(stringAsFactors = FALSE)


# %% Subs
readb <- function(x) {
  .content <- suppressWarnings(readLines(x))
  data.frame(path = str_split_1(.content, ","))
}


sampleidfrompath <- function(x) {
  basename(dirname(x))
}


min_replicates <- function(x) {
  x %>%
    table() %>%
    min() %>%
    as.numeric()
}


excludena <- function(x, y) {
  x[, !is.na(y$group)]
}

# %% Main
argv <- docopt::docopt(doc)
nofilter <- argv$`nofilter`
output_dir <- argv$`output_dir`

path_b1 <- argv$`b1`
path_b2 <- argv$`b2`
path_txi <- argv$`txi`

txi <- readRDS(path_txi)

sample_set <- bind_rows(
  readb(path_b1) %>%
    mutate(group = "set1"),
  readb(path_b2) %>%
    mutate(group = "set2")
) %>%
  mutate(sample_id = sampleidfrompath(path)) %>%
  select(-path)

col_data <- data.frame(sample_id = colnames(txi$gene$counts)) %>%
  left_join(sample_set, by = "sample_id") %>%
  column_to_rownames("sample_id")

.rows_na <- col_data %>%
  filter(is.na(group))

if (nrow(.rows_na) > 0) {
  txi$gene$counts <- excludena(txi$gene$counts, col_data)
  txi$gene$abundance <- excludena(txi$gene$abundance, col_data)
  txi$gene$length <- excludena(txi$gene$length, col_data)

  col_data <- col_data %>%
    filter(!is.na(group))

  message(paste0("These ", nrow(.rows_na), " samples were excluded from this comparison:"))
  message(paste(rownames(.rows_na), collapse = ","))
}

col_data$group <- factor(col_data$group)
dds <- DESeq2::DESeqDataSetFromTximport(
  txi$gene,
  col_data,
  ~group
)

CUTOFF_RAW <- 10

if (!nofilter) {
  keep <- rowSums(DESeq2::counts(dds) > CUTOFF_RAW) %>%
    {
      . >= min_replicates(col_data$group)
    }
  dds <- dds[keep, ]
}

dds <- DESeq2::estimateSizeFactors(dds)
dds <- DESeq2::estimateDispersions(dds)
dds <- DESeq2::nbinomWaldTest(dds)

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

DESeq2::results(dds, contrast = c("group", "set2", "set1")) %>%
  data.frame() %>%
  arrange(padj) %>%
  rownames_to_column(var = "feature_id") %>%
  write_tsv(file.path(output_dir, "wald.tsv"))

DESeq2::counts(dds, normalized = TRUE) %>%
  data.frame() %>%
  rownames_to_column(var = "feature_id") %>%
  write_tsv(file.path(output_dir, "normalized_counts.tsv"))
