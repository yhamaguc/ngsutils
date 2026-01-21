#!/bin/bash

# Usage: ./dl_gencode_human.sh <relase>

ver="$1"

BASE_URL="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_${ver}"

echo "Downloading GENCODE v${ver} files..."
prefixes=(
  2wayconspseudos.gff3.gz
  2wayconspseudos.gtf.gz
  annotation.gff3.gz
  annotation.gtf.gz
  basic.annotation.gff3.gz
  basic.annotation.gtf.gz
  chr_patch_hapl_scaff.annotation.gff3.gz
  chr_patch_hapl_scaff.annotation.gtf.gz
  chr_patch_hapl_scaff.basic.annotation.gff3.gz
  chr_patch_hapl_scaff.basic.annotation.gtf.gz
  lncRNA_transcripts.fa.gz
  long_noncoding_RNAs.gff3.gz
  long_noncoding_RNAs.gtf.gz
  metadata.Annotation_remark.gz
  metadata.EntrezGene.gz
  metadata.Exon_supporting_feature.gz
  metadata.Gene_source.gz
  metadata.HGNC.gz
  metadata.PDB.gz
  metadata.PolyA_feature.gz
  metadata.Pubmed_id.gz
  metadata.RefSeq.gz
  metadata.Selenocysteine.gz
  metadata.SwissProt.gz
  metadata.TrEMBL.gz
  metadata.Transcript_source.gz
  metadata.Transcript_supporting_feature.gz
  pc_transcripts.fa.gz
  pc_translations.fa.gz
  polyAs.gff3.gz
  polyAs.gtf.gz
  primary_assembly.annotation.gff3.gz
  primary_assembly.annotation.gtf.gz
  primary_assembly.basic.annotation.gff3.gz
  primary_assembly.basic.annotation.gtf.gz
  promoter_windows.gff3.gz
  tRNAs.gff3.gz
  tRNAs.gtf.gz
  transcript_rankings.txt.gz
  transcripts.fa.gz
)

for p in "${prefixes[@]}"; do
  wget "${BASE_URL}/gencode.v${ver}.${p}" || echo "Warning: Failed to download ${p}"
done

echo "Decompressing files..."
unpigz gencode.v"${ver}".*.gz

indexing_targets=(
  gencode.v"${ver}".basic.annotation.gtf
  gencode.v"${ver}".annotation.gtf
)

for t in "${indexing_targets[@]}"; do
  echo "Converting ${t} to BED and SQLITE format and indexing..."

  if [ ! -f "$t" ]; then
      echo "Skip: $t not found."
      continue
  fi

  gtf2bed.sh "$t"
  gtf2sqlite --gtf "$t"

  bgzip --threads 8 "${t%.gtf}.bed"

  tabix -p bed "${t%.gtf}.bed.gz"
done
  