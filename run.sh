#!/usr/bin/env bash
# =============================================================================
# Driver del pipeline TE -> homología modular. Encadena los pasos 01,02,06,07
# (03 validación Dfam y CD-Hit son opcionales, por flag).
#
#   pipeline/run.sh pipeline/config/config.panTro6.sh            # núcleo 01->05
#   pipeline/run.sh pipeline/config/config.panTro6.sh --dfam     # + validación Dfam (lento/red)
#   pipeline/run.sh pipeline/config/config.panTro6.sh --cdhit    # + CD-Hit (requiere cd-hit)
#
# Pasos:
#   01_intersect.R          hits Dfam ∩ CDS (+ liftOver) -> hits_filtrados_en_cds.tsv
#   02_extraer_cds.R        CDS fasta + QC 'sano' (+ col sano en la tabla) + regiones TE
#   03_validar_dfam.py      (opcional) doble validación contra Dfam
#   04_clusterizar_loci.R   redundancia -> clusters_<GENOME>/{01_hits_magna..03_loci}
#   05_extraer_dominios_aa.R dominios AA por producto -> clusters_<GENOME>/04b_dominios_aa.fasta
# =============================================================================
set -euo pipefail
CONFIG="${1:?uso: run.sh <config.sh> [--dfam] [--cdhit]}"; shift || true
# shellcheck source=/dev/null
source "$CONFIG"
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"; S="${HERE}/scripts"

RUN_DFAM=false; RUN_CDHIT=false
for a in "$@"; do case "$a" in
  --dfam)  RUN_DFAM=true;;
  --cdhit) RUN_CDHIT=true;;
  *) echo "flag desconocido: $a" >&2; exit 2;;
esac; done

mkdir -p "$OUTDIR"
echo "==> pipeline TE->homología | config: $CONFIG | OUTDIR: $OUTDIR"

echo "==> 01 intersect (TE ∩ CDS, liftOver=$LIFTOVER_TARGET)"
"$RSCRIPT" "${S}/01_intersect.R"

echo "==> 02 extraer CDS + QC sano + regiones TE"
"$RSCRIPT" "${S}/02_extraer_cds.R"

if $RUN_DFAM; then
  echo "==> 03 validación Dfam (--all)"
  "${PYTHON:-python3}" "${S}/03_validar_dfam.py" --all
else
  echo "==> 03 validación Dfam: OMITIDA (pasar --dfam para correrla)"
fi

echo "==> 04 clustering de redundancia -> clusters_${GENOME}"
"$RSCRIPT" "${S}/04_clusterizar_loci.R" \
  "${OUTDIR}/hits_filtrados_en_cds.tsv" "${GENOME}" "${OUTDIR}/clusters_${GENOME}"

echo "==> 05 dominios AA (STEP1 crudo 04a + STEP2 dedup-contención 04b, input CD-Hit)"
"$RSCRIPT" "${S}/05_extraer_dominios_aa.R" \
  "${OUTDIR}/clusters_${GENOME}" "${OUTDIR}/02_extraccion/cds_completos.fasta"

if $RUN_CDHIT; then
  echo "==> 06 CD-Hit sobre 04b_dominios_aa.fasta"
  bash "${S}/06_cdhit.sh" "${OUTDIR}/clusters_${GENOME}/04b_dominios_aa.fasta" || \
    echo "   (06 CD-Hit falló o no está instalado; se omite)"
fi

echo "==> LISTO. Tablas en ${OUTDIR}; clusters + dominios en ${OUTDIR}/clusters_${GENOME}"
