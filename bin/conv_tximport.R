#! /usr/bin/env Rscript

"Convert any tool results to tximport count matrix

Usage:
  conv_tximport.R --gtf <PATH> --type <TYPE> [--output-dir <PATH>] <input-dir>...

Arguments:
  <input-dir>  : Directory containing count data files;
                 kallisto: abundance.h5, RSEM: quantified.isoforms.results, StringTie: t_data.ctab, Salmon: quant.sf

Options:
  --gtf <PATH>         : GTF file
  --type <TYPE>        : stringtie/kallisto/rsem/salmon
  --output-dir <PATH>  : Output directory [default: .]
" -> doc

# %%
# Prepare
#
suppressPackageStartupMessages({
  library(tximport)
  library(rtracklayer)
  library(tidyverse)
  library(ggcorrplot)
  library(docopt)
})


# %%
# Constants
#

# NOTE: ggcorrplot's `lab` prints one number per sample pair, so the label count
#   grows with the square of the sample count. A dozen samples is legible; a
#   448-sample cohort would put ~200,000 text elements into a single SVG, which
#   is unreadable and slow to render. Above this many samples the figure is still
#   drawn, only the printed numbers are dropped.
MAX_SAMPLES_FOR_CORR_LABELS <- 30


# %%
# Subs
#
load_gtf <- function(path, cols, types = c("transcript")) {
  gtf <- path |> readGFF(version = 2L, tags = cols, filter = list(type = types))
  gtf <- gtf |> select(all_of(cols))
  return(gtf)
}


load_data <- function(type, inputs, t2g) {
  sample_names <- inputs |>
    dirname() |>
    basename()
  names(inputs) <- sample_names

  # NOTE: For RSEM recommended before import cut off non-required columns except 1-8
  # cat rsem.isoforms.results | cut -f 1-8
  txi.tx <- inputs |>
    tximport(
      type = type,
      txIn = TRUE,
      txOut = TRUE
    )

  txi.gene <- inputs |>
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


draw_corr_ <- function(mat) {
  ggcorrplot(
    cor(mat, method = "spearman"),
    lab = ncol(mat) <= MAX_SAMPLES_FOR_CORR_LABELS
  )
}


# %%
# Main
#
argv <- docopt(doc)

# NOTE: docopt maps <input-dir> to input_dir and --output-dir to output_dir
gtf <- argv$gtf
type <- argv$type
output_dir <- argv$output_dir
input_dir <- argv$input_dir

message("Input: gtf=", gtf, " output_dir=", output_dir)

t2p <- list(
  kallisto = "abundance.h5",
  rsem = "quantified.isoforms.results",
  stringtie = "t_data.ctab",
  salmon = "quant.sf"
)

inputs <- list.files(input_dir,
  pattern = t2p[[type]],
  full.names = TRUE,
  recursive = TRUE
)

t2g <- load_gtf(
  gtf,
  cols = c("transcript_id", "gene_id", "gene_name"),
  types = c("exon")
) |> distinct()

results <- load_data(type, inputs, t2g)

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

saveRDS(results, file = file.path(output_dir, "txi.rds"))

counts_(results$transcript) |> write.table(
  file = file.path(output_dir, "count_matrix_transcript.tsv"),
  quote = FALSE,
  sep = "\t",
  col.names = NA
)

counts_(results$gene) |> write.table(
  file = file.path(output_dir, "count_matrix_gene.tsv"),
  quote = FALSE,
  sep = "\t",
  col.names = NA
)

vst_transcript <- vst_(results$transcript) |>
  SummarizedExperiment::assay()

vst_transcript |>
  write.table(
    file = file.path(output_dir, "vst_transcript.tsv"),
    quote = FALSE,
    sep = "\t",
    col.names = NA
  )

ggsave(
  filename = file.path(output_dir, "vst_corr_transcript.svg"),
  plot = draw_corr_(vst_transcript)
)

vst_gene <- vst_(results$gene) |>
  SummarizedExperiment::assay()

vst_gene |>
  write.table(
    file = file.path(output_dir, "vst_gene.tsv"),
    quote = FALSE,
    sep = "\t",
    col.names = NA
  )

ggsave(
  filename = file.path(output_dir, "vst_corr_gene.svg"),
  plot = draw_corr_(vst_gene)
)
