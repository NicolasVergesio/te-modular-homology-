#!/usr/bin/env Rscript
# 04_clusterizar_loci.R — Reduccion de redundancia de hits TE por locus genomico.
#
# Esquema "anotar, NO borrar" pero con salida NAVEGABLE por niveles. Lee
# hits_filtrados_en_cds.tsv (una fila por hit x CDS/isoforma) y produce una CARPETA
# clusters_<GENOME>/ con tres tablas encadenadas por IDs cruzados + un mini-reporte:
#
#   01_hits_magna.tsv  (todas las filas hit x isoforma)  -> + id_fila, id_pos,
#                       id_locus, rol, flag, locus_id_final. Apunta HACIA ARRIBA.
#   02_posiciones.tsv  (1 fila por POSICION genomica; colapsa la multiplicacion por
#                       isoforma) -> lleva ids_magna (lista) y su id_locus.
#   03_loci.tsv        (1 fila por INSERCION = "hits finales depurados") -> lleva
#                       ids_pos (lista), tipo uni/multihit y columna 'estado' (ok/
#                       revisar) para filtrar el nucleo limpio.
#   REPORTE_clusters.txt  (cascada de conteos: magna -> posiciones -> clusters ->
#                       loci, uni/multi, flags, isoformas/genes).
#
# Los tres niveles:
#   magna (N filas)        redundante: mismo TE fisico contado por cada isoforma/exon.
#   posicion genomica      dedup por coords + cadena + subfamilia (el TE fisico).
#   locus / insercion      posiciones agrupadas por reglas familia-aware (el EVENTO).
# Conteos: nº inserciones = nrow(03_loci); nº isoformas = uniqueN(transcript_id).
#
# Reglas (justificadas en REPORTE_extraccion_regiones.md, secciones "Reduccion de
# redundancia de hits" y "Validacion empirica de las reglas"):
#   Parche #1: misma familia se colapsa solo con solape RECIPROCO >= 50% (no
#              single-linkage transitivo) -> no fusiona L1 en tandem en un solo locus.
#   Parche #2: el SVA es CONTENEDOR de un dominio tipo-Alu (5', antisentido). Un Alu
#              >= 80% contenido en el SVA Y en cadena OPUESTA -> es_subparte_sva, el
#              evento es el SVA; el e-value NO decide. Contenido pero misma cadena o
#              en el tercio 3' -> revisar_sva. Solape < 50% -> modulo Alu independiente.
#              Zona gris 50-80% -> revisar_sva.
#   Parche #3: la raiz de familia sale de repFamily; si esta vacia, cae a repClass y
#              luego a la raiz de name_te (nunca clusteriza por string vacio).
#   Familias distintas anidadas/adyacentes -> se CONSERVAN (homologia modular).
#   Split (raro): misma subfamilia exacta flanqueando otra familia, con hits anidados
#              en medio y span-union <= p90 del ancho de la familia -> un solo locus.
#
# Uso:  Rscript 04_clusterizar_loci.R [entrada.tsv] [GENOME] [outdir]
# Default: hits_filtrados_en_cds.tsv  panTro6  clusters_panTro6/
suppressPackageStartupMessages(library(data.table))

args   <- commandArgs(trailingOnly = TRUE)
IN     <- if (length(args) >= 1) args[1] else "hits_filtrados_en_cds.tsv"
GENOME <- if (length(args) >= 2) args[2] else "panTro6"
OUTDIR <- if (length(args) >= 3) args[3] else sprintf("clusters_%s", GENOME)
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)
MAGNA <- file.path(OUTDIR, "01_hits_magna.tsv")
POS   <- file.path(OUTDIR, "02_posiciones.tsv")
LOCI  <- file.path(OUTDIR, "03_loci.tsv")
REP   <- file.path(OUTDIR, "REPORTE_clusters.txt")

# --- umbrales (respaldados por la distribucion real, ver reporte) ---
SVA_CONT   <- 0.80   # contencion del Alu dentro del SVA para es_subparte_sva
GRAY_LOW   <- 0.50   # por debajo -> modulo Alu independiente; entre GRAY_LOW y SVA_CONT -> revisar
RECIP      <- 0.50   # solape reciproco minimo para colapsar misma familia (Parche #1)
# p90 del ancho de hit por familia -> tope del span-union en un merge/split
FAM_P90 <- c(Alu = 309L, L1 = 1664L, SVA = 464L, MIR = 243L, L2 = 874L,
             `ERVL-MaLR` = 512L)

stopifnot(file.exists(IN))
d <- fread(IN, sep = "\t", na.strings = c("", "NA"))
orig0 <- copy(names(d))                         # columnas originales (copy: names() es por referencia en data.table)
req <- c("seq_te","start_te","end_te","strand_te","name_te","dfam_id","e_value",
         "repClass","repFamily","coverage","transcript_id","gene_id")
miss <- setdiff(req, names(d))
if (length(miss)) stop("Faltan columnas en la entrada: ", paste(miss, collapse=", "))
d[, id_fila := sprintf("m%05d", .I)]            # ID unico por fila de la magna

# --- raiz de familia (Parche #3): repFamily -> repClass -> raiz de name_te ---
fam_root_of <- function(repFamily, repClass, name_te) {
  fr <- ifelse(is.na(repFamily) | !nzchar(repFamily), NA_character_, repFamily)
  fr <- ifelse(is.na(fr) & !is.na(repClass) & nzchar(repClass), repClass, fr)
  root_name <- sub("[_-].*$", "", name_te)
  ifelse(is.na(fr), root_name, fr)
}
d[, fam_root := fam_root_of(repFamily, repClass, name_te)]
d[, is_sva := (!is.na(repFamily) & repFamily == "SVA") | grepl("^SVA", name_te)]
d[, is_alu := !is.na(repFamily) & repFamily == "Alu"]
d[, ev := suppressWarnings(as.numeric(e_value))]
d[, cov_num := suppressWarnings(as.numeric(coverage))]

# --- tabla de HITS genomicos unicos (colapsa la multiplicacion por isoforma) ---
# El hit genomico se define por coords + cadena + subfamilia; e_value es constante por
# hit, coverage (solape con CDS) varia por isoforma -> se toma el maximo (mejor CDS).
hkey_cols <- c("seq_te","start_te","end_te","strand_te","name_te")
d[, hit_key := do.call(paste, c(.SD, sep = "|")), .SDcols = hkey_cols]
h <- d[, .(seq_te = seq_te[1], start_te = start_te[1], end_te = end_te[1],
           strand_te = strand_te[1], name_te = name_te[1], dfam_id = dfam_id[1],
           fam_root = fam_root[1], is_sva = is_sva[1], is_alu = is_alu[1],
           ev = ev[1], cov_num = max(cov_num, na.rm = TRUE)), by = hit_key]

# --- clustering crudo single-linkage por solape genomico (locus_id_inicial) ---
setorder(h, seq_te, start_te, end_te)
h[, prev_maxend := shift(cummax(end_te)), by = seq_te]
h[, nuevo_cl := is.na(prev_maxend) | start_te > prev_maxend]
h[, cl := cumsum(nuevo_cl)]
h[, locus_id_inicial := sprintf("%s:%d-%d", seq_te, min(start_te), max(end_te)), by = cl]

# ------------------------------------------------------------------------------
# Resolver por cluster: devuelve, por hit, un subgrupo (locus final) y un flag.
# ------------------------------------------------------------------------------
resolve_cluster <- function(g) {
  n <- nrow(g)
  if (n == 1L) return(list(sub = 1L, flag = ""))

  parent <- seq_len(n)
  find <- function(x) { while (parent[x] != x) { parent[x] <<- parent[parent[x]]; x <- parent[x] }; x }
  uni  <- function(a, b) { ra <- find(a); rb <- find(b); if (ra != rb) parent[ra] <<- rb }

  # Convencion half-open: features adyacentes comparten coordenada (fin de una =
  # inicio de la siguiente) -> width = end - start, solape sin +1 (asi tocarse = 0).
  s <- g$start_te; e <- g$end_te; w <- e - s
  ov <- function(i, j) max(0L, min(e[i], e[j]) - max(s[i], s[j]))
  flag <- rep("", n)

  sva_idx <- which(g$is_sva)
  alu_idx <- which(g$is_alu)

  # --- Parche #2: SVA contenedor ---
  if (length(sva_idx) && length(alu_idx)) {
    for (a in alu_idx) {
      ovs <- vapply(sva_idx, function(k) ov(a, k), integer(1))
      if (all(ovs == 0L)) next
      k    <- sva_idx[which.max(ovs)]
      cont <- ov(a, k) / w[a]
      opp  <- g$strand_te[a] != g$strand_te[k]
      rel  <- ((s[a] + e[a]) / 2 - s[k]) / w[k]          # 0..1 en coords genomicas
      if (g$strand_te[k] == "-") rel <- 1 - rel          # orienta 5'->3'
      pos3 <- rel > 2/3
      if (cont >= SVA_CONT && opp && !pos3) {
        uni(a, k); flag[a] <- "es_subparte_sva"
      } else if (cont >= SVA_CONT) {
        flag[a] <- "revisar_sva"                          # contenido pero misma cadena / 3'
      } else if (cont < GRAY_LOW) {
        flag[a] <- "modulo_alu_independiente"
      } else {
        flag[a] <- "revisar_sva"                          # zona gris 50-80%
      }
    }
  }

  # --- misma familia ---
  # Dos situaciones distintas del mismo fam_root:
  #  (a) MISMA subfamilia exacta (name_te) -> fragmentacion o split: colapsar a un
  #      locus (span-union), acotado por p90 del ancho de la familia si se conoce.
  #  (b) DISTINTA subfamilia -> categoria B (ambiguedad) solo si solape RECIPROCO
  #      >= 50% (Parche #1); si no, se conservan separadas (L1 en tandem, no uno).
  if (n >= 2L) for (i in 1:(n - 1L)) for (j in (i + 1L):n) {
    if (g$fam_root[i] != g$fam_root[j]) next
    o <- ov(i, j)
    if (g$name_te[i] == g$name_te[j]) {
      lo <- min(s[i], s[j]); hi <- max(e[i], e[j]); spanw <- hi - lo
      p90 <- FAM_P90[g$fam_root[i]]                       # NA con nombre si la familia no esta tabulada
      if (is.na(p90) || spanw <= p90) {
        # hits de OTRA familia estrictamente contenidos entre las dos piezas -> split
        between <- which(seq_len(n) != i & seq_len(n) != j &
                         s >= min(e[i], e[j]) & e <= max(s[i], s[j]) &
                         g$fam_root != g$fam_root[i])
        uni(i, j)
        if (length(between) && o == 0L) {                 # piezas flanqueando un anidado
          nf <- paste(unique(g$fam_root[between]), collapse = "/")
          flag[i] <- paste0("split_por:", nf); flag[j] <- paste0("split_por:", nf)
          for (b in between) if (flag[b] == "")
            flag[b] <- paste0("anidado_en:", g$fam_root[i])
        } else {                                          # fragmentacion simple
          if (flag[i] == "") flag[i] <- "fragmentado"
          if (flag[j] == "") flag[j] <- "fragmentado"
        }
      }
    } else {
      recip <- if (o > 0L) o / min(w[i], w[j]) else 0
      if (recip >= RECIP) {
        uni(i, j)
        if (flag[i] == "") flag[i] <- "ambiguedad_subfamilia"
        if (flag[j] == "") flag[j] <- "ambiguedad_subfamilia"
      }
    }
  }

  # --- familias distintas: anidado vs adyacente (se conservan, solo se etiquetan) ---
  for (i in seq_len(n)) {
    if (flag[i] != "") next
    others <- which(seq_len(n) != i & g$fam_root != g$fam_root[i])
    others <- others[vapply(others, function(k) ov(i, k) > 0L, logical(1))]
    if (!length(others)) next
    cont_by <- others[vapply(others, function(k) s[k] <= s[i] && e[k] >= e[i], logical(1))]
    if (length(cont_by)) flag[i] <- paste0("anidado_en:", g$fam_root[cont_by[1]])
    else                 flag[i] <- "modulo_distinto_adyacente"
  }

  list(sub = vapply(seq_len(n), find, integer(1)), flag = flag)
}

# --- aplicar el resolver a cada cluster inicial ---
setorder(h, cl, start_te, end_te)
h[, `:=`(sub = NA_integer_, flag = "")]
for (c_id in unique(h$cl)) {
  idx <- which(h$cl == c_id)
  r <- resolve_cluster(h[idx])
  set(h, idx, "sub", r$sub)
  set(h, idx, "flag", r$flag)
}
h[, grp := .GRP, by = .(cl, sub)]

# --- locus_id_final: span-UNION del grupo; etiqueta desde el representante ---
# representante = menor e-value; en grupos con SVA el evento ES el SVA (e-value NO decide)
h[, gstart := min(start_te), by = grp]
h[, gend   := max(end_te),   by = grp]
pick_rep <- function(sub) {
  cand <- if (any(sub$is_sva)) which(sub$is_sva) else seq_len(nrow(sub))
  o <- order(sub$ev[cand], -replace(sub$cov_num[cand], is.na(sub$cov_num[cand]), -Inf))
  cand[o[1]]
}
h[, es_representante := FALSE]
h[, rep_name := NA_character_]
for (gg in unique(h$grp)) {
  idx <- which(h$grp == gg)
  r <- pick_rep(h[idx])
  set(h, idx[r], "es_representante", TRUE)
  set(h, idx, "rep_name", h$name_te[idx[r]])
}
h[, locus_id_final := sprintf("%s:%s:%d-%d", seq_te, rep_name, gstart, gend)]

# ==============================================================================
# IDs cortos y estables (ordenados por posicion genomica) + rol legible
# ==============================================================================
setorder(h, seq_te, start_te, end_te, name_te)
h[, id_pos := sprintf("p%05d", .I)]                       # 1 por posicion genomica
loc_order <- h[, .(o = min(.I)), by = locus_id_final][order(o)]
loc_order[, id_locus := sprintf("L%05d", .I)]             # 1 por insercion/locus
h[loc_order, id_locus := i.id_locus, on = "locus_id_final"]
h[, rol := fifelse(es_representante, "representante", "miembro")]

# ==============================================================================
# 01_hits_magna.tsv  — todas las filas hit x isoforma, con punteros hacia arriba
# ==============================================================================
kmap <- h[, .(hit_key, id_pos, id_locus, rol, flag, locus_id_final)]
d2 <- merge(d, kmap, by = "hit_key", all.x = TRUE, sort = FALSE)
magna <- d2[, c(orig0, "id_fila", "id_pos", "id_locus", "rol", "flag",
                "locus_id_final"), with = FALSE]
setorder(magna, id_fila)
fwrite(magna, MAGNA, sep = "\t")

# ==============================================================================
# 02_posiciones.tsv  — 1 fila por posicion genomica, con la lista de ids_magna
# ==============================================================================
mg <- d[, .(ids_magna = paste(id_fila, collapse = ","),
            n_filas = .N, n_iso = uniqueN(transcript_id)), by = hit_key]
pos <- merge(h, mg, by = "hit_key", sort = FALSE)
pos <- pos[order(id_pos), .(id_pos, chr = seq_te, start = start_te, end = end_te,
             strand = strand_te, subfam = name_te, familia = fam_root, evalue = ev,
             id_locus, rol, flag, n_iso, n_filas, ids_magna)]
fwrite(pos, POS, sep = "\t")

# ==============================================================================
# 03_loci.tsv  — 1 fila por insercion (hits finales depurados) + estado ok/revisar
# ==============================================================================
loci <- h[order(id_locus), .(
  chr        = seq_te[1L],
  start      = min(start_te),
  end        = max(end_te),
  tipo       = fifelse(.N > 1L, "multihit", "unihit"),
  n_pos      = .N,
  subfam_rep = rep_name[1L],
  familia    = fam_root[which(es_representante)[1L]],
  evalue_rep = ev[which(es_representante)[1L]],
  id_pos_rep = id_pos[which(es_representante)[1L]],   # posicion REPRESENTANTE del locus
  subfams    = paste(unique(name_te), collapse = ","),
  flags      = { fl <- unique(flag[nzchar(flag)]); if (length(fl)) paste(fl, collapse = ";") else "" },
  estado     = fifelse(any(flag == "revisar_sva"), "revisar", "ok"),
  ids_pos    = paste(id_pos, collapse = ",")          # TODAS las posiciones (rep incluida)
), by = id_locus]
has_sano <- "sano" %in% names(d2)
iso <- d2[order(id_fila), .(n_iso = uniqueN(transcript_id), n_genes = uniqueN(gene_id),
                            n_iso_sano = if (has_sano) uniqueN(transcript_id[sano %in% TRUE]) else NA_integer_,
                            n_filas = .N, ids_magna = paste(id_fila, collapse = ",")),
          by = id_locus]
loci <- merge(loci, iso, by = "id_locus", sort = FALSE)
setcolorder(loci, c("id_locus","chr","start","end","tipo","subfam_rep","familia",
                    "evalue_rep","estado","n_pos","n_iso","n_iso_sano","n_genes","n_filas",
                    "id_pos_rep","ids_pos","subfams","flags","ids_magna"))
setorder(loci, id_locus)
fwrite(loci, LOCI, sep = "\t")

# ==============================================================================
# REPORTE_clusters.txt  — cascada de conteos y desgloses
# ==============================================================================
n_magna    <- nrow(d)
n_pos      <- nrow(h)
n_cl       <- uniqueN(h$cl)
n_cl_multi <- h[, .N, by = cl][N > 1L, .N]
n_cl_uni   <- n_cl - n_cl_multi
n_loci     <- nrow(loci)
n_loci_mul <- loci[tipo == "multihit", .N]
n_loci_uni <- loci[tipo == "unihit",  .N]
n_iso_tot  <- uniqueN(d$transcript_id)
n_gen_tot  <- uniqueN(d$gene_id)
n_ok       <- loci[estado == "ok", .N]
n_rev      <- loci[estado == "revisar", .N]
ft <- h[flag != "", .N, by = .(f = sub(":.*$", "", flag))][order(-N)]

L <- c(
  sprintf("REPORTE DE CLUSTERIZACION DE HITS TE  (genoma: %s)", GENOME),
  sprintf("Generado: %s", format(Sys.time(), "%Y-%m-%d %H:%M")),
  sprintf("Entrada: %s", IN),
  strrep("=", 64),
  "",
  "CASCADA DE REDUCCION DE REDUNDANCIA (cada nivel colapsa el anterior)",
  strrep("-", 64),
  sprintf("  %6d  filas magna        (hit x isoforma/exon; redundante)", n_magna),
  sprintf("  %6d  posiciones genomicas (dedup coords+cadena+subfamilia)", n_pos),
  sprintf("  %6d  clusters iniciales   (solape genomico crudo)", n_cl),
  sprintf("             - unihit : %d", n_cl_uni),
  sprintf("             - multihit: %d", n_cl_multi),
  sprintf("  %6d  LOCI finales / inserciones  (reglas familia-aware)", n_loci),
  sprintf("             - unihit : %d", n_loci_uni),
  sprintf("             - multihit: %d", n_loci_mul),
  "",
  "  Nota: hay mas loci que clusters iniciales porque las reglas VUELVEN",
  "  a separar modulos distintos (familias anidadas/adyacentes, Alu vs SVA)",
  "  que el solape crudo habia juntado.",
  "",
  "TABLA FINAL DEPURADA (03_loci.tsv)",
  strrep("-", 64),
  sprintf("  %6d  inserciones totales", n_loci),
  sprintf("             - estado 'ok'     : %d", n_ok),
  sprintf("             - estado 'revisar': %d  (revision manual)", n_rev),
  sprintf("  %6d  isoformas afectadas (transcript_id distintos)", n_iso_tot),
  sprintf("  %6d  genes afectados     (gene_id distintos)", n_gen_tot),
  "",
  "FLAGS POR POSICION GENOMICA",
  strrep("-", 64)
)
if (nrow(ft)) L <- c(L, sprintf("  %-28s %d", ft$f, ft$N)) else L <- c(L, "  (sin flags)")
L <- c(L, "",
  "SALIDAS",
  strrep("-", 64),
  "  01_hits_magna.tsv  : todas las filas; punteros id_pos / id_locus.",
  "  02_posiciones.tsv  : 1 x posicion; lista ids_magna, apunta a id_locus.",
  "  03_loci.tsv        : 1 x insercion (hits finales); lista ids_pos, estado.",
  "",
  "NAVEGACION: de un locus -> ids_pos -> (en 02) ids_magna -> filas en 01.",
  "            de una fila  -> id_locus / id_pos -> tablas 02 y 03.")
writeLines(L, REP)

# ------------------------------------------------------------------------------
# Resumen breve a stderr
# ------------------------------------------------------------------------------
msg <- function(...) cat(..., "\n", file = stderr())
# metrics.tsv (append; OUTDIR viene del config via el driver)
OUTDIR_M <- Sys.getenv("OUTDIR")
if (nzchar(OUTDIR_M)) {
  metric <- function(m, v) cat(sprintf("04\t%s\t%s\n", m, v),
                               file = file.path(OUTDIR_M, "metrics.tsv"), append = TRUE)
  metric("posiciones_genomicas", n_pos)
  metric("loci_inserciones", n_loci)
  metric("loci_ok", n_ok); metric("loci_revisar", n_rev)
  metric("isoformas_afectadas", n_iso_tot); metric("genes_afectados", n_gen_tot)
}
msg("04: entrada        :", IN, "(", n_magna, "filas magna )")
msg("04: posiciones      :", n_pos)
msg("04: clusters inic.  :", n_cl, "(", n_cl_multi, "multi-hit )")
msg("04: loci finales    :", n_loci, "( ok", n_ok, "/ revisar", n_rev, ")")
msg("04: isoformas/genes :", n_iso_tot, "/", n_gen_tot)
msg("04: salida ->", OUTDIR, "/ {01_hits_magna, 02_posiciones, 03_loci}.tsv + REPORTE_clusters.txt")
