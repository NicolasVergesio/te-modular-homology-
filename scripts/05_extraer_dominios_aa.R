#!/usr/bin/env Rscript
# 05_extraer_dominios_aa.R — dominios proteicos TE-derivados para CD-Hit (dos pasos).
#
# STEP 1  -> 04a_hits_aa.fasta : TODOS los hits, uno por (posicion genomica x isoforma),
#            uniendo los trozos de exon de ese elemento. Crudo y redundante (por las dudas).
# STEP 2  -> 04b_dominios_aa.fasta : dentro de cada LOCUS, dedup por CONTENCION de secuencia
#            (si una AA es subcadena de otra -> queda la mas larga; el header lista todos los
#            TEs, representante primero). Secuencias que solo se solapan parcial se conservan
#            (CD-Hit las junta). Nunca se dedupea entre loci distintos (senal de homologia);
#            los modulos genuinamente distintos ya quedaron en loci distintos por el paso 04.
#
# Uso:  Rscript 05_extraer_dominios_aa.R [clusters_dir] [cds_completos.fasta]
suppressPackageStartupMessages({ library(data.table); library(Biostrings) })

args  <- commandArgs(trailingOnly = TRUE)
CLDIR <- if (length(args) >= 1) args[1] else "clusters_panTro6"
CDS   <- if (length(args) >= 2) args[2] else "cds_completos.fasta"
MAGNA <- file.path(CLDIR, "01_hits_magna.tsv")
FA1   <- file.path(CLDIR, "04a_hits_aa.fasta")
FA2   <- file.path(CLDIR, "04b_dominios_aa.fasta")
TAB2  <- file.path(CLDIR, "04b_dominios_aa.tsv")
TAB3  <- file.path(CLDIR, "04b_productos_tx.tsv")   # tabla larga: 1 fila por hit
stopifnot(file.exists(MAGNA), file.exists(CDS))

# --- traducir CDS (corta el header en '|' -> transcript_id) ---
cds <- readDNAStringSet(CDS)
names(cds) <- sub("[|[:space:]].*$", "", names(cds))
cds <- subseq(cds, 1L, width(cds) %/% 3L * 3L)
prot <- translate(cds, if.fuzzy.codon = "solve", no.init.codon = TRUE)
pv <- setNames(as.character(prot), names(prot))

d <- fread(MAGNA, sep = "\t", na.strings = c("", "NA"))
d[, `:=`(start_aa = as.integer(start_aa), end_aa = as.integer(end_aa))]
if (!"sano" %in% names(d)) d[, sano := NA]
if (!"uniprot_ids" %in% names(d)) d[, uniprot_ids := NA_character_]
if (!"gene_id" %in% names(d)) d[, gene_id := NA_character_]
if (!"rol" %in% names(d)) d[, rol := ""]

# formatea logico/NA -> "TRUE"/"FALSE"/"NA" (para el campo sano de los headers)
fmt_sano <- function(x) ifelse(is.na(x), "NA", ifelse(x, "TRUE", "FALSE"))

# --- convencion de separadores en los campos multi-valor -----------------------
#   ';'  separa TRANSCRIPTOS (un slot por tx, en el mismo orden que la columna tx)
#   ','  separa varios valores DENTRO de un mismo transcripto
#   'NA' ocupa el slot cuando ese transcripto no tiene valor (nunca se omite)
# Asi 'uniprot', 'gene' y 'sano_tx' quedan alineados posicionalmente con 'tx':
# el slot i corresponde siempre al i-esimo transcripto. (La tabla larga
# 04b_productos_tx.tsv trae lo mismo sin codificar, una fila por hit.)
# El paso 01 junta las accessions de un tx con ';', asi que se normalizan a ','.
norm_multi <- function(x) gsub(";", ",", x, fixed = TRUE)

# valores de 'vals' agrupados por transcripto, alineados a utx
por_tx <- function(vals, txs, utx) {
  vapply(utx, function(u) {
    v <- unlist(strsplit(vals[txs == u], "[;,]"))
    v <- unique(v[!is.na(v) & nzchar(v) & v != "NA"])
    if (!length(v)) "NA" else paste(v, collapse = ",")
  }, character(1), USE.NAMES = FALSE)
}

# ============================================================================
# STEP 1 — un hit por (posicion x isoforma): union de sus trozos de exon
# ============================================================================
h <- d[!is.na(start_aa) & !is.na(end_aa),
       .(id_locus = id_locus[1], name_te = name_te[1], rol = rol[1],
         gene_id = gene_id[1], uniprot = uniprot_ids[1], sano = any(sano %in% TRUE),
         aa_start = min(start_aa), aa_end = max(end_aa)),
       by = .(id_pos, transcript_id)]
get_seq <- function(tx, s, e) { p <- pv[[tx]]; if (is.null(p) || is.na(p)) return(NA_character_)
  e <- min(e, nchar(p)); if (s > e) return(""); substr(p, s, e) }
h[, seq := mapply(get_seq, transcript_id, aa_start, aa_end)]
h <- h[!is.na(seq) & nchar(seq) > 0L]
h[, flag_stop := grepl("\\*", seq)]

setorder(h, id_locus, id_pos, transcript_id)
h[, hit_id := paste(id_pos, transcript_id, sep = "__")]
hdr1 <- sprintf(">%s | id_pos=%s | te=%s | locus=%s | gene=%s | tx=%s | uniprot=%s | sano=%s | aa=%d-%d%s",
                h$hit_id, h$id_pos, h$name_te, h$id_locus,
                ifelse(is.na(h$gene_id), "NA", h$gene_id), h$transcript_id,
                ifelse(is.na(h$uniprot), "NA", norm_multi(h$uniprot)), fmt_sano(h$sano),
                h$aa_start, h$aa_end,
                ifelse(h$flag_stop, " STOP_INTERNO", ""))
fa1 <- character(2L * nrow(h)); fa1[c(TRUE, FALSE)] <- hdr1; fa1[c(FALSE, TRUE)] <- h$seq
writeLines(fa1, FA1)

# ============================================================================
# STEP 2 — por locus, dedup por CONTENCION de secuencia
# ============================================================================
mods <- list(); largo <- list(); mid <- 0L
for (loc in unique(h$id_locus)) {
  g <- h[id_locus == loc]
  # mas larga primero; en empate, el representante del paso 04 primero
  g <- g[order(-nchar(seq), rol != "representante")]
  kept <- list()
  for (i in seq_len(nrow(g))) {
    s <- g$seq[i]; absorbed <- FALSE
    for (k in seq_along(kept)) {
      if (grepl(s, kept[[k]]$seq, fixed = TRUE)) {         # s subcadena de un modulo -> absorber
        kept[[k]]$tes  <- c(kept[[k]]$tes,  g$name_te[i])
        kept[[k]]$tx   <- c(kept[[k]]$tx,   g$transcript_id[i])
        kept[[k]]$gene <- c(kept[[k]]$gene, g$gene_id[i])
        kept[[k]]$uni  <- c(kept[[k]]$uni,  g$uniprot[i])
        kept[[k]]$sano <- kept[[k]]$sano || g$sano[i]            # OR -> columna sano (tabla)
        kept[[k]]$sano_vec <- c(kept[[k]]$sano_vec, g$sano[i])   # por-tx -> header 04b encadenado
        kept[[k]]$id_pos_v <- c(kept[[k]]$id_pos_v, g$id_pos[i]) # por-hit -> tabla larga
        kept[[k]]$s_v <- c(kept[[k]]$s_v, g$aa_start[i])
        kept[[k]]$e_v <- c(kept[[k]]$e_v, g$aa_end[i])
        absorbed <- TRUE; break
      }
    }
    if (!absorbed) kept[[length(kept) + 1L]] <- list(
      seq = s, rep = g$name_te[i], tes = g$name_te[i], tx = g$transcript_id[i],
      gene = g$gene_id[i], uni = g$uniprot[i], aa_start = g$aa_start[i], aa_end = g$aa_end[i],
      tx_rep = g$transcript_id[i], sano = g$sano[i], sano_vec = g$sano[i], stop = g$flag_stop[i],
      id_pos_v = g$id_pos[i], s_v = g$aa_start[i], e_v = g$aa_end[i])
  }
  for (m in kept) {
    utx <- unique(m$tx)                          # tx en orden de aparicion
    sano_tx <- paste(fmt_sano(m$sano_vec[match(utx, m$tx)]), collapse = ";")  # sano alineado a utx
    mid <- mid + 1L
    mods[[length(mods) + 1L]] <- data.table(
      mod_id = mid,
      id_locus = loc, gene = paste(por_tx(m$gene, m$tx, utx), collapse = ";"), rep = m$rep,
      tes = paste(unique(m$tes), collapse = ";"),
      n_te = length(unique(m$tes)), tx_rep = m$tx_rep,
      tx = paste(utx, collapse = ";"), n_tx = length(utx),
      uniprot = paste(por_tx(m$uni, m$tx, utx), collapse = ";"),
      aa_start = m$aa_start, aa_end = m$aa_end,
      len_aa = nchar(m$seq), sano = m$sano, flag_stop = m$stop, seq = m$seq,
      sano_tx = sano_tx)
    # tabla larga: una fila por HIT, sin colapsar nada (unico lugar donde TE y tx
    # pueden variar a la vez sin ambiguedad)
    largo[[length(largo) + 1L]] <- data.table(
      mod_id = mid, id_locus = loc, id_pos = m$id_pos_v, te = m$tes, tx = m$tx,
      gene = ifelse(is.na(m$gene), "NA", m$gene),
      uniprot = ifelse(is.na(m$uni), "NA", norm_multi(m$uni)),
      sano = fmt_sano(m$sano_vec), aa_start = m$s_v, aa_end = m$e_v)
  }
}
dist <- rbindlist(mods)
setorder(dist, id_locus, -len_aa)
dist[, producto_id := sprintf("%s_%d", id_locus, seq_len(.N)), by = id_locus]
setcolorder(dist, c("producto_id", "id_locus", "gene", "rep", "tes", "n_te", "tx_rep", "tx",
                    "n_tx", "uniprot", "aa_start", "aa_end", "len_aa", "sano", "flag_stop", "seq",
                    "sano_tx"))

hdr2 <- sprintf(">%s | gene=%s | rep=%s | tes=%s | tx=%s | uniprot=%s | sano=%s | aa=%d-%d | len=%d%s",
                dist$producto_id, dist$gene, dist$rep, dist$tes, dist$tx, dist$uniprot,
                dist$sano_tx, dist$aa_start, dist$aa_end, dist$len_aa,
                ifelse(dist$flag_stop, " STOP_INTERNO", ""))
fa2 <- character(2L * nrow(dist)); fa2[c(TRUE, FALSE)] <- hdr2; fa2[c(FALSE, TRUE)] <- dist$seq
writeLines(fa2, FA2)
setcolorder(dist, c(setdiff(names(dist), c("sano_tx", "seq", "mod_id")), "sano_tx", "seq"))
fwrite(dist[, !"mod_id"], TAB2, sep = "\t")

# --- tabla larga: una fila por hit (sin campos multi-valor) ---
lar <- rbindlist(largo)
lar <- merge(lar, dist[, .(mod_id, producto_id)], by = "mod_id", sort = FALSE)[, !"mod_id"]
setcolorder(lar, c("producto_id", "id_locus", "id_pos", "te", "tx", "gene", "uniprot",
                   "sano", "aa_start", "aa_end"))
setorder(lar, producto_id, tx, id_pos)
fwrite(lar, TAB3, sep = "\t")

# --- metrics + resumen ---
msg <- function(...) cat(..., "\n", file = stderr())
OUTDIR_M <- Sys.getenv("OUTDIR")
if (nzchar(OUTDIR_M)) {
  metric <- function(m, v) cat(sprintf("05\t%s\t%s\n", m, v),
                               file = file.path(OUTDIR_M, "metrics.tsv"), append = TRUE)
  metric("hits_aa_step1", nrow(h))
  metric("dominios_step2", nrow(dist))
  metric("dominios_con_stop_interno", dist[flag_stop == TRUE, .N])
  metric("loci_multidominio", dist[, .N, by = id_locus][N > 1, .N])
  metric("dominios_sanos_sin_stop", dist[sano == TRUE & flag_stop == FALSE, .N])
}
msg("05 STEP1: hits crudos (posicion x isoforma) :", nrow(h), "->", FA1)
msg("05 STEP2: dominios distintos (dedup contencion):", nrow(dist))
msg("05:   loci con >1 dominio :", dist[, .N, by = id_locus][N > 1, .N])
msg("05:   con STOP interno    :", dist[flag_stop == TRUE, .N],
    " | sanos sin stop:", dist[sano == TRUE & flag_stop == FALSE, .N])
msg("05:   longitud AA -> min", min(dist$len_aa), "mediana", as.integer(median(dist$len_aa)),
    "max", max(dist$len_aa))
msg("05: salida ->", FA2, "+", TAB2)
