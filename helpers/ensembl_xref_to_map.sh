#!/usr/bin/env bash
# =============================================================================
# ensembl_xref_to_map.sh — convierte el xref UniProt de Ensembl al mapeo de 2
# columnas que consume el pipeline (variable ENST2UNIPROT del config).
#
#   helpers/ensembl_xref_to_map.sh Especie.Ensamblado.RELEASE.uniprot.tsv.gz \
#       > data/map_tx_to_uniprot.tsv
#
# El input se baja del FTP de Ensembl (tsv/<especie>/), y conviene usar la MISMA
# release que el GTF: así las isoformas coinciden. Alternativa: helpers/map_uniprot.qmd
# (biomaRt), que necesita red; este script no.
#
# REGLA: se conservan TODAS las accessions UniProt de cada transcripto, una por fila.
# NO se prioriza ni se descarta ninguna, a proposito: helpers/validar_vs_uniprot.py
# prueba todas contra el proteoma, se queda con el mejor resultado y deja constancia de
# cual uso en su columna 'uniprot_usado'. Elegir una sola aca seria tirar informacion
# antes de saber cual sirve.
#
# Se aceptan estos db_name (lista blanca, no descarte):
#   Uniprot_isoform    accession de isoforma (P25440-3): la proteina EXACTA de ese
#                      transcripto. Humano tiene 35202; chimp, 13.
#   Uniprot/SWISSPROT  entrada curada.
#   Uniprot/SPTREMBL   entrada predicha por computadora.
# El resto (p.ej. Uniprot_gn_trans_name, que son NOMBRES de gen y no accessions) se ignora.
#
# Dato util al interpretar los resultados: la columna info_type distingue DIRECT (el link
# viene declarado) de SEQUENCE_MATCH (Ensembl lo asigno alineando secuencias). Las
# SEQUENCE_MATCH casi no existen en el FASTA del proteoma de referencia de UniProt
# (chimp 0/24199; humano 57/3709), asi que suelen caer como 'accession ausente del FASTA'
# en la validacion. No es un problema biologico: tienen 100% de identidad, son la MISMA
# proteina con otro accession, de un proteoma redundante que UniProt no distribuye.
# =============================================================================
set -euo pipefail

IN="${1:?uso: ensembl_xref_to_map.sh <Especie...uniprot.tsv.gz> > map_tx_to_uniprot.tsv}"
[ -r "$IN" ] || { echo "no puedo leer: $IN" >&2; exit 1; }

CAT=cat; case "$IN" in *.gz) CAT=zcat;; esac

$CAT "$IN" | awk -F'\t' -v OFS='\t' -v fn="$IN" '
  NR == 1 {
    for (i = 1; i <= NF; i++) col[$i] = i          # ubicar columnas por NOMBRE
    tx = col["transcript_stable_id"]; xr = col["xref"]; db = col["db_name"]
    it = col["info_type"]
    if (!tx || !xr || !db || !it) {
      print "ERROR: falta transcript_stable_id/xref/db_name/info_type en el header." > "/dev/stderr"
      print "       No parece un xref de Ensembl: " fn > "/dev/stderr"
      exit 1
    }
    print "ensembl_transcript_id", "uniprot_cons"
    next
  }
  {
    if ($db != "Uniprot_isoform" && $db != "Uniprot/SWISSPROT" && $db != "Uniprot/SPTREMBL")
      next                                          # lista blanca: el resto se ignora
    t = $tx
    if (!(t in acc)) acc[t] = $xr
    else if (index(";" acc[t] ";", ";" $xr ";") == 0) acc[t] = acc[t] ";" $xr
  }
  END {
    # una FILA por accession (no unidas por ";"): asi el paso 01 cuenta bien la
    # ambiguedad en n_ids. El propio pipeline las junta con ";" en la columna uniprot.
    for (t in acc) { n = split(acc[t], a, ";"); for (i = 1; i <= n; i++) print t, a[i] }
  }
' | { read -r h1 h2; printf "%s\t%s\n" "$h1" "$h2"; sort; }
