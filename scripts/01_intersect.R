#!/usr/bin/env Rscript
# 01_intersect.R — Intersect TE hits (Dfam nrph) x CDS del GTF -> tabla de hits rica.
#
# Port parametrizado de pan_tro.qmd (front-end de chimp). Reproduce el esquema de
# hits_filtrados_en_cds.tsv (36 columnas) que consumen los pasos 02/04/05.
# Todo se parametriza por env-vars (config.<especie>.sh); nada de rutas fijas.
#
# LiftOver FLEXIBLE (LIFTOVER_TARGET): 'gtf' lifta el CDS del GTF al ensamblado de
# los hits (caso chimp: GTF Ensembl -> panTro6); 'hits' lifta los hits al del GTF;
# 'none' (o LIFTOVER_CHAIN vacio) no lifta (hits y GTF ya en la misma version).
#
# Salidas en OUTDIR: all_cds_coord.tsv, hits_totales_en_cds.tsv, hits_filtrados_en_cds.tsv
suppressPackageStartupMessages({
  library(data.table); library(rtracklayer); library(GenomicRanges)
  library(GenomicFeatures); library(janitor); library(tidyverse)
  library(GenomeInfoDb)
})

cfg <- function(k, d = "") { v <- Sys.getenv(k); if (nzchar(v)) v else d }
GTF        <- cfg("GTF");        stopifnot(nzchar(GTF), file.exists(GTF))
TE_HITS    <- cfg("TE_HITS");    stopifnot(nzchar(TE_HITS), file.exists(TE_HITS))
CHAIN      <- cfg("LIFTOVER_CHAIN")
LO_TARGET  <- cfg("LIFTOVER_TARGET", if (nzchar(CHAIN)) "gtf" else "none")
CHROM_ALIAS<- cfg("CHROM_ALIAS")               # TSV con columnas ucsc/genbank/refseq (opcional)
CHROMS_RAW <- cfg("CHROMS")                    # nombres crudos del GTF a conservar, coma-separados
DFAM_TAXA  <- cfg("DFAM_TAXA")
ENST2UNI   <- cfg("ENST2UNIPROT")              # mapeo precomputado transcript->uniprot (opcional)
EVAL_THR   <- as.numeric(cfg("E_VALUE_THRESHOLD", "1e-10"))
OVL_THR    <- as.numeric(cfg("OVERLAP_BP_THR", "30"))
OUTDIR     <- cfg("OUTDIR", "."); dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)
if (nzchar(CHAIN)) stopifnot(file.exists(CHAIN))
stopifnot(LO_TARGET %in% c("gtf", "hits", "none"))

# --- metrics.tsv: registro largo (paso, metrica, valor). El paso 01 lo REINICIA. ---
METRICS <- file.path(OUTDIR, "metrics.tsv")       # en la RAIZ de resultados (accesible)
cat("paso\tmetrica\tvalor\n", file = METRICS)     # trunca (arranca la corrida)
metric <- function(paso, met, val) cat(sprintf("%s\t%s\t%s\n", paso, met, val),
                                       file = METRICS, append = TRUE)
# subcarpeta de detalle del paso 01 (los archivos intermedios no ensucian la raiz)
D01 <- file.path(OUTDIR, "01_intersect"); dir.create(D01, showWarnings = FALSE, recursive = TRUE)

# ---------------------------------------------------------------------------
# 1) GTF: importar, quedarse con cromosomas de interes, estilo UCSC (chr..., chrM)
# ---------------------------------------------------------------------------
message("01: importando GTF ", basename(GTF))
gtf <- import(GTF)
if (nzchar(CHROMS_RAW)) {
  chroms <- trimws(strsplit(CHROMS_RAW, ",")[[1]])
  gtf <- keepSeqlevels(gtf, chroms, pruning.mode = "coarse")
}
seqlevels(gtf) <- paste0("chr", seqlevels(gtf))
seqlevels(gtf)[seqlevels(gtf) == "chrMT"] <- "chrM"

# CDS de genes/transcriptos protein_coding
cds <- gtf[gtf$gene_biotype == "protein_coding" &
           gtf$transcript_biotype == "protein_coding" &
           as.character(gtf$type) == "CDS"]
message("01: features CDS protein_coding: ", length(cds))
metric("01", "cds_features_protein_coding", length(cds))

# ---------------------------------------------------------------------------
# 2) LiftOver del GTF (si corresponde) -> gr_cds queda en el ensamblado de los hits
# ---------------------------------------------------------------------------
chain <- if (nzchar(CHAIN)) import.chain(CHAIN) else NULL
if (LO_TARGET == "gtf") {
  message("01: liftOver del CDS del GTF (", basename(CHAIN), ")")
  gr_cds <- unlist(liftOver(cds, chain))
  message("01:   CDS mapeados: ", round(length(gr_cds) / length(cds) * 100, 2), "%")
  metric("01", "cds_lifteados_pct", round(length(gr_cds) / length(cds) * 100, 2))
} else {
  gr_cds <- cds
}

# ---------------------------------------------------------------------------
# 3) Hits TE (nrph): fix coords, mapear cromosomas al estilo del GTF (UCSC)
# ---------------------------------------------------------------------------
hits <- fread(TE_HITS) %>% clean_names()
hits <- hits %>% mutate(start = pmin(ali_st, ali_en), end = pmax(ali_st, ali_en))

if (nzchar(CHROM_ALIAS)) {
  # el seq_name de los hits es un accession (genbank); traducir a UCSC (chrN)
  alias <- fread(CHROM_ALIAS) %>% clean_names() %>% as_tibble()
  hits <- hits %>%
    left_join(alias, by = c(number_seq_name = "genbank")) %>%
    dplyr::filter(complete.cases(number_ucsc)) %>%
    dplyr::rename(seqname_ucsc = number_ucsc)
} else {
  # sin alias: normalizar los seqnames de los hits IGUAL que el GTF (chrN, MT->M)
  # para que matcheen aunque vengan como '1' (Ensembl) o ya como 'chr1' (UCSC).
  sn <- as.character(hits$number_seq_name)
  sn <- ifelse(grepl("^chr", sn), sn, paste0("chr", sn))
  sn <- ifelse(sn == "chrMT", "chrM", sn)
  hits <- hits %>% mutate(seqname_ucsc = sn)
}

hits_gr <- makeGRangesFromDataFrame(as.data.frame(hits), keep.extra.columns = TRUE,
                                    seqnames.field = "seqname_ucsc",
                                    start.field = "start", end.field = "end",
                                    strand.field = "strand")
if (LO_TARGET == "hits") {
  message("01: liftOver de los hits (", basename(CHAIN), ")")
  hits_gr <- unlist(liftOver(hits_gr, chain))
}

# ---------------------------------------------------------------------------
# 4) Intersect CDS x hits -> df_intersect (esquema rico)
# ---------------------------------------------------------------------------
ov <- findOverlaps(gr_cds, hits_gr, ignore.strand = TRUE)
cds_hit <- gr_cds[queryHits(ov)]
te_hit  <- hits_gr[subjectHits(ov)]
inter   <- pintersect(cds_hit, te_hit, ignore.strand = TRUE)

df <- data.frame(
  seq_cds = as.character(seqnames(cds_hit)), start_cds = start(cds_hit),
  end_cds = end(cds_hit), strand_gene = as.character(strand(cds_hit)),
  gene_id = cds_hit$gene_id, transcript_id = cds_hit$transcript_id,
  transcritp_v = cds_hit$transcript_version, gene_bio = cds_hit$gene_biotype,
  transcript_bio = cds_hit$transcript_biotype, exon_number = cds_hit$exon_number,
  cds_phase = cds_hit$phase,
  seq_te = as.character(seqnames(te_hit)), start_te = start(te_hit),
  end_te = end(te_hit), strand_te = as.character(strand(te_hit)),
  name_te = te_hit$family_name, dfam_id = te_hit$family_acc,
  e_value = te_hit$e_value, sq_len = te_hit$sq_len, kimura_te = te_hit$kimura_div,
  solapamiento_start = start(inter), solapamiento_end = end(inter),
  coverage = width(inter), stringsAsFactors = FALSE)

# --- anotacion Dfam (name, repClass, repFamily, class_family, kw_raw, species) ---
if (nzchar(DFAM_TAXA) && file.exists(DFAM_TAXA)) {
  df_te_type <- read_tsv(DFAM_TAXA, show_col_types = FALSE)
  df <- df %>% left_join(df_te_type, by = c(dfam_id = "accession"),
                         relationship = "many-to-one")
} else {                                          # sin taxonomia: crear las columnas en NA
  df <- df %>% mutate(name = NA_character_, repClass = NA_character_,
                      repFamily = NA_character_, class_family = NA_character_,
                      kw_raw = NA_character_, species = NA_character_)
}

# --- uniprot desde mapeo precomputado (n_ids, uniprot_ids) ---
if (nzchar(ENST2UNI) && file.exists(ENST2UNI)) {
  m <- read_tsv(ENST2UNI, show_col_types = FALSE)
  mg <- m %>% dplyr::filter(!is.na(uniprot_cons)) %>%
    group_by(ensembl_transcript_id) %>%
    summarise(n_ids = n_distinct(uniprot_cons),
              uniprot_ids = paste0(unique(uniprot_cons), collapse = ";"), .groups = "drop")
  df <- df %>% left_join(mg, by = c(transcript_id = "ensembl_transcript_id"),
                         relationship = "many-to-one")
} else {
  df$n_ids <- NA_integer_; df$uniprot_ids <- NA_character_
}

# ---------------------------------------------------------------------------
# 5) Coordenadas relativas al CDS y aminoacidicas (pmapToTranscripts, frame-aware)
# ---------------------------------------------------------------------------
# Un liftOver correcto conserva, por transcripto, el MISMO nº de exones CDS, en un
# solo cromosoma y una sola hebra. Se excluye (y ANOTA) todo lo demas:
#   roto_chrom_hebra (varios chr/hebra) | multimapeo (out>in, exon duplicado) |
#   perdida_parcial (out<in, exon perdido). Solo 'ok' entra a la proyeccion a aa.
ncds_in <- table(cds$transcript_id)               # nº exones CDS ANTES del liftOver
cds_by_tx_all <- split(gr_cds, gr_cds$transcript_id)
n_strand <- lengths(unique(strand(cds_by_tx_all)))
n_chrom  <- lengths(unique(seqnames(cds_by_tx_all)))
n_out    <- lengths(cds_by_tx_all)
n_in     <- as.integer(ncds_in[names(cds_by_tx_all)])
motivo <- ifelse(n_chrom > 1 | n_strand > 1, "roto_chrom_hebra",
          ifelse(n_out > n_in, "multimapeo",
          ifelse(n_out < n_in, "perdida_parcial", "ok")))
tx_ok  <- names(cds_by_tx_all)[motivo == "ok"]
excl <- data.frame(transcript_id = names(cds_by_tx_all), n_cds_in = n_in,
                   n_cds_out = as.integer(n_out), n_chrom = as.integer(n_chrom),
                   n_strand = as.integer(n_strand), motivo = motivo,
                   stringsAsFactors = FALSE)
excl <- excl[excl$motivo != "ok", ]
write_tsv(excl, file.path(D01, "transcriptos_excluidos_liftover.tsv"))
tab <- table(factor(motivo, levels = c("ok","roto_chrom_hebra","multimapeo","perdida_parcial")))
message("01: liftOver por tx -> ok:", tab["ok"], " roto:", tab["roto_chrom_hebra"],
        " multimapeo:", tab["multimapeo"], " perdida_parcial:", tab["perdida_parcial"],
        " (excluidos anotados en transcriptos_excluidos_liftover.tsv)")
metric("01", "tx_lifteados", length(cds_by_tx_all))
metric("01", "tx_ok", tab["ok"])
metric("01", "tx_excl_roto_chrom_hebra", tab["roto_chrom_hebra"])
metric("01", "tx_excl_multimapeo", tab["multimapeo"])
metric("01", "tx_excl_perdida_parcial", tab["perdida_parcial"])
cds_by_tx <- cds_by_tx_all[tx_ok]

ok <- df$transcript_id %in% tx_ok
overlap_gr <- GRanges(df$seq_cds[ok],
                      IRanges(df$solapamiento_start[ok], df$solapamiento_end[ok]),
                      strand = df$strand_gene[ok])
cds_coords <- pmapToTranscripts(overlap_gr, cds_by_tx[df$transcript_id[ok]])
df$start_cds_rel <- NA_integer_; df$end_cds_rel <- NA_integer_; df$width_cds_rel <- NA_integer_
df$start_cds_rel[ok] <- start(cds_coords)
df$end_cds_rel[ok]   <- end(cds_coords)
df$width_cds_rel[ok] <- width(cds_coords)
df <- df %>% mutate(start_aa = ceiling(start_cds_rel / 3),
                    end_aa   = ceiling(end_cds_rel   / 3))

# ---------------------------------------------------------------------------
# 6) Filtrar y guardar
# ---------------------------------------------------------------------------
df_filt <- df %>% dplyr::filter(e_value <= EVAL_THR, coverage >= OVL_THR,
                                complete.cases(width_cds_rel))

write_tsv(as.data.frame(gr_cds), file.path(D01, "all_cds_coord.tsv"))
write_tsv(df,      file.path(D01, "hits_totales_en_cds.tsv"))
write_tsv(df_filt, file.path(OUTDIR, "hits_filtrados_en_cds.tsv"))   # tabla central en la raiz

# --- metricas del paso 01 ---
metric("01", "solapes_crudos", nrow(df))
metric("01", "solapes_pasan_evalue", sum(df$e_value <= EVAL_THR, na.rm = TRUE))
metric("01", "solapes_pasan_coverage", sum(df$coverage >= OVL_THR, na.rm = TRUE))
metric("01", "hits_filtrados", nrow(df_filt))
# composicion final por clado (informativo; NO se filtra por clado)
comp <- sort(table(df_filt$species), decreasing = TRUE)
for (cl in names(comp)) metric("01", paste0("clado:", cl), comp[[cl]])

message("01: hits totales ", nrow(df), " -> filtrados ", nrow(df_filt),
        " (e_value<=", EVAL_THR, ", coverage>=", OVL_THR, ")")
message("01: salidas -> ", OUTDIR, "/hits_filtrados_en_cds.tsv + metrics.tsv (raiz); ",
        "01_intersect/{all_cds_coord,hits_totales_en_cds,transcriptos_excluidos_liftover}")
