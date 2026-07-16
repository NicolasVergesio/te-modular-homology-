#!/usr/bin/env bash
# =============================================================================
# Config del pipeline TE -> homología modular, para CHIMPANCÉ (panTro6).
# Para otra especie: copiar este archivo, ajustar rutas/cromosomas/liftOver.
# Todas las variables se exportan y las leen los scripts R/Python (Sys.getenv /
# os.environ). No hay parser de YAML: es a propósito, sin dependencias extra.
# =============================================================================
# Ruta ABSOLUTA a los datos de chimp en ESTA máquina (editar para la tuya).
DATA="/home/nvergesio/doctorado/transposon-mediated_modular_homology/pan_tro"
EP="${DATA}/ensembl_pipe"

# ============================ INPUTS =========================================
export GTF="${EP}/Pan_troglodytes.Pan_tro_3.0.116.gtf.gz"   # GTF Ensembl
export TE_HITS="${DATA}/DApanTro2.nrph.hits.gz"             # hits Dfam (nrph)
export GENOME_FASTA="${EP}/panTro6.fa"                      # genoma (para .fai/traducir)
export DFAM_TAXA="${DATA}/dfam_curated_taxa.tsv"           # repClass/repFamily por accession
export CHROM_ALIAS="${EP}/panTro6.chromAlias.txt"          # mapa seqname(genbank)->UCSC (hits)
export ENST2UNIPROT="${EP}/map_tx_to_uniprot.tsv"          # mapeo transcripto->uniprot precomputado

# ====================== LIFTOVER (flexible) ==================================
# gtf  = liftar el CDS del GTF al ensamblado de los hits (caso chimp: Ensembl->panTro6)
# hits = liftar los hits al ensamblado del GTF
# none = hits y GTF ya en la misma versión (o LIFTOVER_CHAIN vacío)
export LIFTOVER_CHAIN="${EP}/panTro5ToPanTro6.over.chain"
export LIFTOVER_TARGET="gtf"

# --- cromosomas del GTF a conservar (nombres crudos; se les antepone 'chr' y MT->M) ---
export CHROMS="1,2A,2B,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X,Y,MT"

# ============================ FILTROS ========================================
export E_VALUE_THRESHOLD="1e-10"
export OVERLAP_BP_THR="30"

# ===================== VALIDACIÓN DFAM (paso 03, opcional) ===================
export DFAM_ORGANISM="Homo sapiens"       # 'Pan troglodytes' da ERROR en Dfam; Alu/L1/... compartidos

# ============================ SALIDAS ========================================
export GENOME="panTro6"                   # etiqueta -> carpeta clusters_<GENOME>
# Carpeta de RESULTADOS dedicada (separada de los datos crudos en ${EP}).
# Adentro: hits_filtrados_en_cds.tsv + metrics.tsv (raiz) + subcarpetas 01_intersect/,
# 02_extraccion/, clusters_<GENOME>/.
export OUTDIR="${EP}/resultados_${GENOME}"

# ============================ BINARIOS =======================================
export RSCRIPT="${RSCRIPT:-/home/nvergesio/micromamba/envs/R_analisis_rna_transposones/bin/Rscript}"
export PYTHON="${PYTHON:-python3}"
