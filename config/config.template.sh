#!/usr/bin/env bash
# =============================================================================
# PLANTILLA de config. Copiar a config.<especie>.sh y completar los <...>.
# Todas las variables se exportan y las leen los scripts R/Python (Sys.getenv /
# os.environ). No hay parser de YAML: a propósito, sin dependencias extra.
# Ver config.panTro6.sh para un ejemplo completo (chimpancé).
# =============================================================================

# --- carpeta con los DATOS de esta especie y donde irán los RESULTADOS ---
DATA="/ruta/a/tus/datos/<especie>"                    # <-- EDITAR (absoluta)

# ============================ INPUTS =========================================
export GTF="${DATA}/anotacion.gtf.gz"                 # GTF Ensembl de la especie
export TE_HITS="${DATA}/hits.nrph.hits.gz"            # hits Dfam (formato nrph)
export GENOME_FASTA="${DATA}/genoma.fa"               # genoma (se indexa .fai solo)
export DFAM_TAXA="/ruta/a/dfam_curated_taxa.tsv"      # GLOBAL de Dfam (helpers/parse_dfam_embl.sh); mismo archivo para todas las especies
export CHROM_ALIAS=""                                 # mapa seqname->UCSC SOLO si los hits usan accesiones GenBank. Si usan chr1/1, dejar vacío.
export ENST2UNIPROT=""                                # mapeo transcripto->uniprot (opcional; helpers/map_uniprot.qmd). Vacío = columna uniprot en NA.

# ====================== LIFTOVER (si hace falta) =============================
# Solo si el GTF y los hits están en ENSAMBLADOS distintos. Si coinciden, dejar vacío/none.
export LIFTOVER_CHAIN=""                              # ruta al .over.chain, o vacío
export LIFTOVER_TARGET="none"                         # gtf (lifta el GTF) | hits (lifta los hits) | none

# --- cromosomas del GTF a conservar (nombres crudos; se les antepone 'chr' y MT->M) ---
export CHROMS="1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,X,Y,MT"   # <-- ajustar (ratón 1-19; grandes simios agregan 2A,2B)

# ============================ FILTROS ========================================
export E_VALUE_THRESHOLD="1e-10"
export OVERLAP_BP_THR="30"

# ===================== VALIDACIÓN DFAM (paso 03, opcional) ===================
export DFAM_ORGANISM="Homo sapiens"                  # organismo para la API de Dfam; algunos (p.ej. 'Pan troglodytes') fallan -> usar un proxy cercano

# ============================ SALIDAS ========================================
export GENOME="<especie>"                            # etiqueta corta -> carpeta clusters_<GENOME>
export OUTDIR="${DATA}/resultados_${GENOME}"          # TODO lo generado va acá (separado por especie)

# ============================ BINARIOS =======================================
export RSCRIPT="${RSCRIPT:-Rscript}"                 # o ruta al Rscript del env con los paquetes
export PYTHON="${PYTHON:-python3}"
