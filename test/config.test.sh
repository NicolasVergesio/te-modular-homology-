#!/usr/bin/env bash
# =============================================================================
# TEST mínimo: datos SINTÉTICOS (chr1 de 2000 bp, 1 gen, 1 TE). Corre en segundos,
# sin datos grandes ni red. Valida que la instalación y el flujo 01→05 funcionen.
# Uso:  bash run.sh test/config.test.sh
# =============================================================================
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"     # carpeta test/
DATA="${HERE}/data"

export GTF="${DATA}/test.gtf"
export TE_HITS="${DATA}/test.nrph.hits"
export GENOME_FASTA="${DATA}/test.fa"
export DFAM_TAXA="${DATA}/test_dfam_taxa.tsv"
export CHROM_ALIAS=""                       # los hits usan '1' -> se normaliza a chr1 solo
export ENST2UNIPROT=""                      # sin uniprot (queda NA)

export LIFTOVER_CHAIN=""                     # mismo ensamblado -> sin liftOver
export LIFTOVER_TARGET="none"
export CHROMS="1"

export E_VALUE_THRESHOLD="1e-10"
export OVERLAP_BP_THR="30"
export DFAM_ORGANISM="Homo sapiens"

export GENOME="test"
export OUTDIR="${HERE}/resultados_test"

export RSCRIPT="${RSCRIPT:-/home/nvergesio/micromamba/envs/R_analisis_rna_transposones/bin/Rscript}"
export PYTHON="${PYTHON:-python3}"
