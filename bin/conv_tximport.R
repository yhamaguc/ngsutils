#! /usr/bin/env Rscript

"Convert any tool results to tximport count matrix

Usage:
  conv_tximport.R --gtf <PATH> --type <TYPE> [--output-dir <PATH>] <input_dirs>...

Options:
  --gtf <PATH>         : GTF file
  --type <TYPE>        : stringtie/kallisto/rsem/salmon
  --output-dir <PATH>  : Output directory [default: .]
  <input_dirs>         : The directories containing count data file;
                         kallisto: abundance.h5, RSEM: quantified.isoforms.results, StringTie: t_data.ctab, Salmon: quant.sf
" -> doc

# %%
# Requires
#
library(tximport)
library(rtracklayer)
library(tidyverse)
library(docopt)


# %%
# Subs
#
load_gtf <- function(path, cols, types = c("transcript")) {
  gtf <- path %>% readGFF(version = 2L, tags = cols, filter = list(type = types))
  gtf <- gtf %>% select(all_of(cols))
  return(gtf)
}


load_data <- function(type, inputs, t2g) {
  sample_names <- inputs %>%
    dirname() %>%
    basename()
  names(inputs) <- sample_names

  # NOTE: For RSEM recommended befor import cut off non-required columns except 1-8
  # cat rsem.isoforms.results | cut -f 1-8
  txi.tx <- inputs %>%
    tximport(
      type = type,
      txIn = TRUE,
      txOut = TRUE
    )

  txi.gene <- inputs %>%
    tximport(
      type = type,
      txIn = TRUE,
      txOut = FALSE,
      tx2gene = t2g,
      countsFromAbundance = "lengthScaledTPM"
    )

  return(
    list(
      transcript = txi.tx,
      gene = txi.gene
    )
  )
}


estimate_ <- function(txi) {
  col_data <- data.frame(row.names = colnames(txi$counts))

  .dds <- DESeq2::DESeqDataSetFromTximport(txi, col_data, ~1)
  .dds <- DESeq2::estimateSizeFactors(.dds)

  .dds
}


counts_ <- function(txi) {
  .counts <- DESeq2::counts(estimate_(txi), normalized = TRUE)

  .counts
}


vst_ <- function(txi) {
  .vst <- DESeq2::vst(estimate_(txi), blind = FALSE)

  .vst
}


# %%
# Main
#

argv <- docopt(doc)

message(argv)

gtf_path <- argv$gtf
output_dir <- argv$`output_dir`
type <- argv$type
input_dirs <- argv$input_dirs

t2p <- list(
  kallisto = "abundance.h5",
  rsem = "quantified.isoforms.results",
  stringtie = "t_data.ctab",
  salmon = "quant.sf"
)

inputs <- list.files(
  input_dirs, pattern = t2p[[type]],
  full.names = TRUE,
  recursive = TRUE
)

t2g <- load_gtf(
  gtf_path,
  cols = c("transcript_id", "gene_id", "gene_name"),
  types = c("exon")
) %>% distinct()

results <- load_data(type, inputs, t2g)

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

saveRDS(results, file = file.path(output_dir, "txi.rds"))

counts_(results$transcript) %>% write.table(
  file = file.path(output_dir, "count_matrix_transcript.tsv"),
  quote = FALSE,
  sep = "\t",
  col.names = NA
)

counts_(results$gene) %>% write.table(
  file = file.path(output_dir, "count_matrix_gene.tsv"),
  quote = FALSE,
  sep = "\t",
  col.names = NA
)

vst_transcript <- vst_(results$transcript) %>%
  SummarizedExperiment::assay()

vst_transcript %>%
  write.table(
    file = file.path(output_dir, "vst_transcript.tsv"),
    quote = FALSE,
    sep = "\t",
    col.names = NA
  )

ggcorrplot::ggcorrplot(cor(vst_transcript, method = "spearman"), lab = TRUE) %>%
  ggsave(filename = file.path(output_dir, "vst_corr_transcript.svg"), plot = .)

vst_gene <- vst_(results$gene) %>%
  SummarizedExperiment::assay()

vst_gene %>%
  write.table(
    file = file.path(output_dir, "vst_gene.tsv"),
    quote = FALSE,
    sep = "\t",
    col.names = NA
  )

ggcorrplot::ggcorrplot(cor(vst_gene, method = "spearman"), lab = TRUE) %>%
  ggsave(filename = file.path(output_dir, "vst_corr_gene.svg"), plot = .)
