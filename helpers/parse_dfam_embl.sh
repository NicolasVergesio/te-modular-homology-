#!/usr/bin/env bash
#
# Parsea el flatfile EMBL curado de Dfam a una tabla TSV.
#
# Extrae, por familia:
#   accession      - AC  (DFxxxxxxx)
#   name           - NM  (ej. MIR, AluY)
#   repClass       - CC RepeatMasker Annotations: Type      (ej. SINE, LINE, LTR, DNA, RC)
#   repFamily      - CC RepeatMasker Annotations: SubType   (ej. Alu, L1, ERVL-MaLR, hAT-Charlie)
#   class_family   - repClass/repFamily  <- mismo formato que repClass/repFamily de UCSC rmsk
#   kw_raw         - linea KW cruda (referencia; en algunas familias es texto libre, no el codigo)
#   species        - CC Species, o el campo OS si falta
#
# La clasificacion se toma del bloque CC (Type/SubType), que es consistente en todas
# las familias. La linea KW se usa solo como fallback porque en algunas familias
# (p.ej. THE1A) contiene una descripcion en texto libre en lugar del codigo Clase/Familia.
#
# NOTA: esta es la clasificacion de 2 niveles estilo RepeatMasker, NO la ruta jerarquica
# completa de Dfam (root;Interspersed_Repeat;Transposable_Element;...). Esa ruta completa
# solo la da la API REST de Dfam, no este flatfile.
#
# Uso:
#   gunzip -k Dfam_curatedonly.embl.gz
#   ./parse_dfam_embl.sh Dfam_curatedonly.embl > dfam_curated_taxa.tsv
# o sin descomprimir a disco:
#   zcat Dfam_curatedonly.embl.gz | ./parse_dfam_embl.sh /dev/stdin > dfam_curated_taxa.tsv

set -euo pipefail

INPUT="${1:-/dev/stdin}"

awk '
BEGIN {
    OFS = "\t"
    print "accession", "name", "repClass", "repFamily", "class_family", "kw_raw", "species"
}
/^AC / { l=$0; sub(/^AC +/,"",l); sub(/;.*/,"",l);             acc=l }
/^NM / { l=$0; sub(/^NM +/,"",l); sub(/[[:space:]]+$/,"",l);   nm=l }
/^KW / { l=$0; sub(/^KW +/,"",l); sub(/\.[[:space:]]*$/,"",l); kw=l }
/^OS / { l=$0; sub(/^OS +/,"",l); sub(/[[:space:]]+$/,"",l);   os=l }
/^CC +Type: /    { l=$0; sub(/^CC +Type: +/,"",l);    sub(/[[:space:]]+$/,"",l); ty=l }
/^CC +SubType: / { l=$0; sub(/^CC +SubType: +/,"",l); sub(/[[:space:]]+$/,"",l); st=l }
/^CC +Species: / { l=$0; sub(/^CC +Species: +/,"",l); sub(/[[:space:]]+$/,"",l); sp=l }
/^\/\// {
    cls=ty; fam=st
    # fallback: sin bloque CC, derivar de KW si tiene formato Clase/Familia
    if (cls=="" && split(kw, a, "/") >= 2) { cls=a[1]; fam=a[2] }
    cf = (fam!="") ? cls"/"fam : cls
    spp = (sp!="") ? sp : os
    print acc, nm, cls, fam, cf, kw, spp
    acc=nm=kw=os=ty=st=sp=""
}
' "$INPUT"
