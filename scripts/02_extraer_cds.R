#!/usr/bin/env Rscript
# 02_extraer_cds.R — CDS completo por transcripto + QC "sano" + region TE (solape).
#
# Port parametrizado de extraer_regiones.qmd. Ademas AGREGA la columna `sano` a
# hits_filtrados_en_cds.tsv (nace aca y se propaga a los pasos 04/05).
#
# Salidas en OUTDIR: cds_completos.fasta, cds_completos_sanos.fasta, cds_qc.tsv,
#   regiones_te.fasta, regiones_te_sanos.fasta, y hits_filtrados_en_cds.tsv (+col sano).
suppressPackageStartupMessages({
  library(data.table); library(GenomicRanges); library(Biostrings)
  library(Rsamtools); library(GenomicFeatures); library(tidyverse)
})

cfg <- function(k, d = "") { v <- Sys.getenv(k); if (nzchar(v)) v else d }
OUTDIR <- cfg("OUTDIR", "."); GENOME <- cfg("GENOME_FASTA")
stopifnot(nzchar(GENOME), file.exists(GENOME))
METRICS <- file.path(OUTDIR, "metrics.tsv")        # append (el paso 01 lo reinicia)
metric <- function(paso, met, val) cat(sprintf("%s\t%s\t%s\n", paso, met, val),
                                       file = METRICS, append = TRUE)
HITS <- file.path(OUTDIR, "hits_filtrados_en_cds.tsv")       # tabla central (raiz)
ACC  <- file.path(OUTDIR, "01_intersect", "all_cds_coord.tsv")
D02  <- file.path(OUTDIR, "02_extraccion"); dir.create(D02, showWarnings = FALSE, recursive = TRUE)
stopifnot(file.exists(HITS), file.exists(ACC))

df_hits <- read_tsv(HITS, show_col_types = FALSE)
df_cds  <- read_tsv(ACC,  show_col_types = FALSE)

if (!file.exists(paste0(GENOME, ".fai"))) { message("02: indexando genoma (.fai)"); indexFa(GENOME) }
fa <- FaFile(GENOME); gchr <- seqnames(seqinfo(fa))

# --- CDS por transcripto (solo tx con hits) ---
df_cds_filt <- df_cds %>% semi_join(df_hits, by = "transcript_id") %>%
  filter(seqnames %in% gchr)
cds_gr  <- GRanges(df_cds_filt$seqnames,
                   IRanges(df_cds_filt$start, df_cds_filt$end),
                   strand = df_cds_filt$strand)
cds_grl <- split(cds_gr, df_cds_filt$transcript_id)
cds_seq <- extractTranscriptSeqs(fa, cds_grl)
# header del FASTA: agrega uniprot si está disponible (transcript_id|uniprot). Los IDs
# internos siguen siendo el transcript_id; el paso 05 corta el header en el '|'.
uni <- df_hits %>% distinct(transcript_id, uniprot_ids)
uni <- setNames(uni$uniprot_ids, uni$transcript_id)
# header homogéneo: siempre 'tx|uniprot'; si no hay uniprot -> 'tx|NA' (placeholder)
fa_hdr <- function(ids) { u <- uni[ids]; u <- ifelse(is.na(u) | !nzchar(u), "NA", u); paste0(ids, "|", u) }
cds_out <- cds_seq; names(cds_out) <- fa_hdr(names(cds_seq))
writeXStringSet(cds_out, file.path(D02, "cds_completos.fasta"))
message("02: CDS completos: ", length(cds_seq))

# --- QC "transcripto sano": multiplo de 3 + ATG inicial + stop en marco 3' ---
stop_pos <- df_cds_filt %>%
  group_by(transcript_id) %>%
  summarise(seqnames = dplyr::first(seqnames), strand = dplyr::first(strand),
            cds_end_genomic = if (dplyr::first(strand) == "+") max(end) else min(start),
            .groups = "drop") %>%
  mutate(start = if_else(strand == "+", cds_end_genomic + 1L, cds_end_genomic - 3L),
         end   = if_else(strand == "+", cds_end_genomic + 3L, cds_end_genomic - 1L))
stop_gr <- GRanges(stop_pos$seqnames, IRanges(stop_pos$start, stop_pos$end),
                   strand = stop_pos$strand)
stop_codon <- setNames(as.character(getSeq(fa, stop_gr)), stop_pos$transcript_id)

cds_qc <- tibble(
  transcript_id = names(cds_seq), width = width(cds_seq),
  multiplo_3  = width(cds_seq) %% 3 == 0,
  empieza_atg = as.character(subseq(cds_seq, 1, 3)) == "ATG",
  stop_codon  = stop_codon[names(cds_seq)],
  stop_ok     = stop_codon[names(cds_seq)] %in% c("TAA", "TAG", "TGA")
) %>% mutate(sano = multiplo_3 & empieza_atg & stop_ok)
write_tsv(cds_qc, file.path(D02, "cds_qc.tsv"))
cds_sanos <- cds_qc %>% filter(sano) %>% pull(transcript_id)
cds_san <- cds_seq[cds_sanos]; names(cds_san) <- fa_hdr(cds_sanos)
writeXStringSet(cds_san, file.path(D02, "cds_completos_sanos.fasta"))
message("02: transcriptos sanos: ", length(cds_sanos), " de ", nrow(cds_qc))
metric("02", "tx_con_hits", length(cds_seq))     # transcriptos distintos con hits (= nº CDS reconstruidos)
metric("02", "tx_sanos", length(cds_sanos))
metric("02", "tx_no_sanos", nrow(cds_qc) - length(cds_sanos))

# --- AGREGAR `sano` a la tabla de hits (nace aca; se propaga a 04/05) ---
df_hits <- df_hits %>% select(-any_of("sano")) %>%
  left_join(cds_qc %>% select(transcript_id, sano), by = "transcript_id")
write_tsv(df_hits, HITS)
message("02: sano agregado a hits_filtrados_en_cds.tsv (",
        sum(df_hits$sano %in% TRUE), " TRUE / ", sum(df_hits$sano %in% FALSE), " FALSE)")
metric("02", "hits_sanos", sum(df_hits$sano %in% TRUE))
metric("02", "hits_no_sanos", sum(df_hits$sano %in% FALSE))

# --- Region del transposon (solape CDS∩TE), strand del gen ---
te_gr <- GRanges(df_hits$seq_cds,
                 IRanges(df_hits$solapamiento_start, df_hits$solapamiento_end),
                 strand = df_hits$strand_gene)
keep   <- as.character(seqnames(te_gr)) %in% gchr
te_seq <- getSeq(fa, te_gr[keep])
names(te_seq) <- with(df_hits[keep, ],
  sprintf("%s|%s|%s|%s:%d-%d(%s)|cds:%s-%s|aa:%s-%s",
          transcript_id, dfam_id, name_te, seq_cds,
          solapamiento_start, solapamiento_end, strand_gene,
          start_cds_rel, end_cds_rel, start_aa, end_aa))
writeXStringSet(te_seq, file.path(D02, "regiones_te.fasta"))

# QC de correctitud: cada región TE debe ser subcadena EXACTA del CDS de su transcripto.
# Valida toda la cadena de coordenadas (liftOver + intersect + hebra + proyección). Debe dar 100%.
tx_keep <- df_hits$transcript_id[keep]
es_sub <- mapply(function(s, tx) if (tx %in% names(cds_seq))
                   grepl(as.character(s), as.character(cds_seq[[tx]]), fixed = TRUE) else NA,
                 te_seq, tx_keep)
pct_sub <- round(mean(es_sub, na.rm = TRUE) * 100, 2)
metric("02", "regiones_te_subcadena_cds_pct", pct_sub)
message("02: región TE subcadena del CDS de su transcripto: ", pct_sub, "% (chequeo de coordenadas)")
if (pct_sub < 100) warning("02: ", sum(!es_sub, na.rm = TRUE),
                           " regiones TE NO son subcadena del CDS -> revisar coordenadas")

te_sano <- df_hits$transcript_id[keep] %in% cds_sanos
writeXStringSet(te_seq[te_sano], file.path(D02, "regiones_te_sanos.fasta"))
message("02: regiones TE: ", length(te_seq), " (sanas: ", sum(te_sano), ")")
message("02: salidas -> ", OUTDIR, "/hits_filtrados_en_cds.tsv (raiz, +sano); ",
        "02_extraccion/{cds_qc,cds_completos*,regiones_te*}")
