#!/usr/bin/env bash
# 06_cdhit.sh — CD-Hit sobre los dominios AA TE-derivados (04b_dominios_aa.fasta).
# Agrupa por identidad de secuencia al 60/70/80 % → grupos de homología modular.
#
# Requiere cd-hit en el PATH (o en la variable CDHIT_BIN).
# Uso:  06_cdhit.sh <04b_dominios_aa.fasta> [outdir]
set -euo pipefail
FASTA="${1:?uso: 06_cdhit.sh <04b_dominios_aa.fasta> [outdir]}"
OUT="${2:-$(dirname "$FASTA")/cdhit}"
CDHIT="${CDHIT_BIN:-cd-hit}"
command -v "$CDHIT" >/dev/null 2>&1 || { echo "06: cd-hit no encontrado (instalalo o seteá CDHIT_BIN)"; exit 1; }
mkdir -p "$OUT"

for id in 0.6 0.7 0.8; do
  # word size de cd-hit segun identidad: 4 para 0.6, 5 para >=0.7
  n=4; awk "BEGIN{exit !($id>=0.7)}" && n=5
  "$CDHIT" -i "$FASTA" -o "$OUT/dominios_c${id}.fasta" -c "$id" -n "$n" -d 0 -M 0 -T 0 >/dev/null
  ncl=$(grep -c '^>' "$OUT/dominios_c${id}.fasta.clstr" 2>/dev/null || echo "?")
  echo "06: cd-hit c=${id} -> $OUT/dominios_c${id}.fasta  (${ncl} clusters)"
done
echo "06: CD-Hit listo en $OUT/"
