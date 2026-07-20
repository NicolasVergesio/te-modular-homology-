#!/usr/bin/env python3
"""Valida las secuencias AA inferidas (04a_hits_aa.fasta) contra las proteinas
reales de UniProt, para estimar que fraccion esta 'bien inferida'.

Idea: cada hit del 04a trae uniprot=<acc>, aa=<start>-<end> y la secuencia
inferida (traduccion del CDS en esa ventana). Se busca la proteina real por su
accession en el FASTA de UniProt y se compara:
  - exacto   : inferido == proteina[start-1:end]            (posicion + secuencia)
  - substring: inferido aparece en algun lado de la proteina (tolera coords corridas)
  - stop     : la traduccion tiene STOP interno ('*')       -> traduccion rota
  - otra_iso : aparece en otra isoforma del mismo accession -> referencia era otra
  - mismatch : no aparece en ningun lado                    -> desajuste real

Requiere un FASTA de UniProt CON isoformas (headers >sp|ACC|... / >tr|ACC|...);
para validar hits con accession isoforma-especifica (P12345-3) hace falta que el
FASTA las incluya (proteome 'all isoforms').

El flag 'sano' se lee del propio header del 04a (lo agrega el paso 05), asi que
las metricas se desglosan por sano=TRUE/FALSE sin archivos extra.

Salidas (por defecto, en la misma carpeta del 04a; --no-write las desactiva):
  05_validacion_uniprot.tsv  : UNA FILA POR HIT con su status (todos, no solo mismatch)
  05_validacion_uniprot.txt  : el reporte de porcentajes, igual que lo impreso

Uso:
  validar_vs_uniprot.py UNIPROT.fasta[.gz] 04a_hits_aa.fasta
  validar_vs_uniprot.py ... --out-tsv T.tsv --report R.txt   # rutas explicitas
  validar_vs_uniprot.py ... --no-write                       # solo stdout
"""
import argparse, gzip, os, re, sys
from collections import Counter


def load_uniprot(path):
    """path -> (acc->seq, acc->'sp'/'tr')."""
    op = gzip.open if path.endswith(".gz") else open
    seqs, src, acc, cur, buf = {}, {}, None, None, []
    with op(path, "rt") as fh:
        for line in fh:
            if line.startswith(">"):
                if acc:
                    seqs[acc] = "".join(buf); src[acc] = cur
                m = re.match(r">(sp|tr)\|([^|]+)\|", line)
                cur = m.group(1) if m else "?"
                acc = m.group(2) if m else line[1:].split()[0]
                buf = []
            else:
                buf.append(line.strip())
    if acc:
        seqs[acc] = "".join(buf); src[acc] = cur
    return seqs, src


def parse_04a(path):
    """Itera dicts por hit (campos del header + la secuencia inferida)."""
    with open(path) as fh:
        hdr = None
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                hdr = line
            elif hdr is not None:
                def g(pat, d="NA"):
                    m = re.search(pat, hdr); return m.group(1) if m else d
                yield dict(id=hdr[1:].split(" |")[0].strip(),
                           id_pos=g(r"id_pos=(\S+)"), locus=g(r"locus=(\S+)"),
                           gene=g(r"gene=(\S+)"), tx=g(r"tx=(\S+)"), te=g(r"te=(\S+)"),
                           uniprot=g(r"uniprot=(\S+)"), sano=g(r"sano=(\S+)"),
                           aa_start=int(g(r"aa=(\d+)-\d+", "0")),
                           aa_end=int(g(r"aa=\d+-(\d+)", "0")), seq=line)
                hdr = None


def main():
    ap = argparse.ArgumentParser(description="Valida 04a_hits_aa.fasta vs UniProt.")
    ap.add_argument("uniprot_fasta", help="FASTA UniProt (con isoformas), .gz ok")
    ap.add_argument("hits_04a", help="04a_hits_aa.fasta del pipeline")
    ap.add_argument("--dump-mismatch", metavar="TSV",
                    help="volcar SOLO los mismatch reales a un TSV (subset de --out-tsv)")
    ap.add_argument("--out-tsv", metavar="TSV",
                    help="tabla una-fila-por-hit con su status "
                         "(default: 05_validacion_uniprot.tsv junto al 04a)")
    ap.add_argument("--report", metavar="TXT",
                    help="guardar el reporte de porcentajes "
                         "(default: 05_validacion_uniprot.txt junto al 04a)")
    ap.add_argument("--gene-xref", metavar="TSV",
                    help="xref de Ensembl (gene_stable_id..xref, .gz ok) para el fallback "
                         "por gen; sin esto el fallback usa solo los accessions que ya "
                         "aparecen en el 04a")
    ap.add_argument("--no-write", action="store_true",
                    help="no escribir nada, solo imprimir el reporte")
    args = ap.parse_args()

    if not args.no_write:
        base = os.path.join(os.path.dirname(os.path.abspath(args.hits_04a)),
                            "05_validacion_uniprot")
        args.out_tsv = args.out_tsv or base + ".tsv"
        args.report = args.report or base + ".txt"

    up, src = load_uniprot(args.uniprot_fasta)
    sys.stderr.write(f"uniprot: {len(up)} secuencias\n")
    base_index = {}
    for a, s in up.items():
        base_index.setdefault(re.sub(r"-\d+$", "", a), []).append(s)

    # --- indice gene -> accessions, para el fallback ---
    # Ensembl solo publica el xref tx->UniProt cuando la traduccion coincide con la
    # proteina, asi que muchos transcriptos quedan sin accession aunque su GEN si
    # tenga proteinas conocidas (via otros transcriptos). Se juntan de dos fuentes.
    gene_accs = {}
    def add_gene_acc(gene, acc):
        if gene and gene != "NA" and acc in up:
            gene_accs.setdefault(gene, set()).add(acc)

    if args.gene_xref:
        op = gzip.open if args.gene_xref.endswith(".gz") else open
        with op(args.gene_xref, "rt") as fh:
            for line in fh:
                f = line.rstrip("\n").split("\t")
                if len(f) > 3 and f[0].startswith("ENS"):
                    add_gene_acc(f[0], f[3])
    for h in parse_04a(args.hits_04a):                 # los que ya trae el propio 04a
        for a in re.split(r"[;,]", h["uniprot"]):
            add_gene_acc(h["gene"], a)
    sys.stderr.write(f"fallback por gen: {len(gene_accs)} genes con proteinas UniProt\n")

    def clasificar_una(uni, s, e, seq):
        base = re.sub(r"-\d+$", "", uni)
        ref = up.get(uni) or up.get(base)           # 1ro isoforma exacta, sino canonica
        if ref is None:
            return "no_acc", base
        if seq == ref[s - 1:e]:
            return "exacto", base
        if seq in ref:
            return "substr", base
        if "*" in seq:
            return "stop", base
        if any(seq in bs for bs in base_index.get(base, ())):
            return "otra_iso", base
        return "mismatch", base

    # el header puede traer varias accessions ('A0A2I3TBD6;A0A6D2XCG3'): se prueban
    # todas y gana la mejor (el mapeo tx->UniProt no siempre desambigua a una sola).
    RANK = {"exacto": 0, "substr": 1, "otra_iso": 2, "stop": 3, "mismatch": 4, "no_acc": 5}

    def fuente_de(a):
        """'sp' (SwissProt, curada) / 'tr' (trEMBL, predicha) / None si no esta en el FASTA."""
        base = re.sub(r"-\d+$", "", a)
        return src.get(base) if (a in up or base in up) else None

    def orden(x):
        """Criterio para elegir entre las accessions de un mismo transcripto.
        1o el status (la evidencia manda: si trEMBL da 'exacto' y SwissProt 'mismatch',
        gana trEMBL). Recien ANTE EMPATE se prefiere la curada, y dentro de ella la
        accession de isoforma, que es la proteina exacta del transcripto. Sin este
        desempate el ganador seria el primero de la cadena, o sea el orden del archivo
        de Ensembl: arbitrario, y dejaria la columna 'fuente' al azar en los
        transcriptos que tienen curada y predicha a la vez (7297 en hg38)."""
        (out, base), acc = x
        return (RANK[out],
                0 if src.get(base) == "sp" else 1,
                0 if re.search(r"-\d+$", acc) else 1)

    def clasificar(uni, s, e, seq, gene):
        """-> (status, base_acc, acc_usada, fuentes_disponibles). Si no hay mapeo directo
        tx->UniProt, cae al fallback por gen (proteinas que UniProt tiene para ESE gen, via
        OTROS transcriptos): no valida coordenadas, solo presencia de la secuencia."""
        fuentes = "NA"
        if uni != "NA":
            cands = [a for a in re.split(r"[;,]", uni) if a]
            fs = sorted({f for f in (fuente_de(a) for a in cands) if f})
            fuentes = ",".join(fs) if fs else "NA"
            (out, base), acc = min(((clasificar_una(a, s, e, seq), a) for a in cands),
                                   key=orden)
            if out != "no_acc":                     # hubo mapeo directo utilizable
                return out, base, acc, fuentes
        # --- fallback por gen ---
        accs = gene_accs.get(gene, ())
        for a in accs:
            if seq in up[a]:
                return "otro_tx_del_gen", re.sub(r"-\d+$", "", a), a, fuentes
        if accs:
            return "producto_distinto", None, ";".join(sorted(accs)[:3]), fuentes
        return (("no_acc", None, uni, fuentes) if uni != "NA"
                else ("na", None, "NA", fuentes))

    COLS = ["id", "id_pos", "locus", "gene", "tx", "te", "uniprot", "uniprot_usado",
            "fuente", "fuentes_disponibles", "sano", "aa_start", "aa_end", "len_aa",
            "status", "seq_inferida"]
    tot = Counter(); by_sano = {"TRUE": Counter(), "FALSE": Counter(), "NA": Counter()}
    rows = []
    for h in parse_04a(args.hits_04a):
        out, base, acc, fuentes = clasificar(h["uniprot"], h["aa_start"], h["aa_end"],
                                             h["seq"], h["gene"])
        tot[out] += 1
        by_sano.get(h["sano"], by_sano["NA"])[out] += 1
        h.update(status=out, uniprot_usado=acc, fuente=src.get(base, "NA") if base else "NA",
                 fuentes_disponibles=fuentes, len_aa=len(h["seq"]), seq_inferida=h["seq"])
        rows.append(h)

    linea = []                                   # el reporte se acumula para poder guardarlo
    def w(s=""): linea.append(s); print(s)

    DIRECTO  = ["exacto", "substr", "stop", "otra_iso", "mismatch"]
    FALLBACK = ["otro_tx_del_gen", "producto_distinto", "no_acc", "na"]

    def reporte(c, titulo):
        comp = sum(c[k] for k in DIRECTO)
        fb   = sum(c[k] for k in FALLBACK)
        pc = lambda n: f"{100*n/comp:5.1f}%" if comp else "  n/a"
        pf = lambda n: f"{100*n/fb:5.1f}%" if fb else "  n/a"
        w(f"\n### {titulo}   (hits={sum(c.values())})")
        w(f"  [A] CON mapeo directo tx->UniProt  : {comp}")
        w(f"    exacto en coordenadas            : {c['exacto']:6d} ({pc(c['exacto'])})")
        w(f"    substring (coord corrida)        : {c['substr']:6d} ({pc(c['substr'])})")
        w(f"    STOP interno (traduccion rota)   : {c['stop']:6d} ({pc(c['stop'])})")
        w(f"    rescatado en otra isoforma       : {c['otra_iso']:6d} ({pc(c['otra_iso'])})")
        w(f"    MISMATCH real                    : {c['mismatch']:6d} ({pc(c['mismatch'])})")
        w(f"  [B] SIN mapeo directo (via gen)    : {fb}")
        w(f"    hallada en otra prot. del gen    : {c['otro_tx_del_gen']:6d} ({pf(c['otro_tx_del_gen'])})")
        w(f"    producto distinto al del gen     : {c['producto_distinto']:6d} ({pf(c['producto_distinto'])})")
        w(f"    accession ausente del FASTA      : {c['no_acc']:6d} ({pf(c['no_acc'])})")
        w(f"    sin UniProt para el gen          : {c['na']:6d} ({pf(c['na'])})")

    def glosario():
        w("\n" + "=" * 66)
        w("QUE SIGNIFICA CADA CATEGORIA")
        w("=" * 66)
        w("""
Cada hit es una ventana AA (aa_start-aa_end) de un transcripto, traducida por el
pipeline. 'Validar' = comparar esa secuencia inferida contra la proteina real de
UniProt. Para eso hace falta saber que proteina corresponde al transcripto: ese
puente lo da el xref de Ensembl, que SOLO existe cuando la traduccion ya coincidia
con UniProt. De ahi los dos bloques:

[A] CON mapeo directo — el header traia accession y esta en el FASTA. Es la
    validacion estricta: se compara secuencia Y posicion.
  exacto      la secuencia inferida es identica a proteina[aa_start-1:aa_end].
              Secuencia y coordenadas correctas. Es el resultado buscado.
  substring   la secuencia aparece en la proteina pero en OTRA posicion: la
              traduccion es correcta y las coordenadas estan corridas.
  stop        la traduccion tiene un '*' (STOP) interno -> el CDS anotado esta
              roto o desfasado; no se compara contra nada.
  otra_iso    no coincide con la isoforma referenciada, pero SI aparece en otra
              isoforma de la MISMA entrada UniProt (P12345-1, -2, ...): la
              referencia elegida era la isoforma equivocada.
  MISMATCH    la secuencia inferida NO aparece en ninguna parte de esa proteina
    real      ni de sus isoformas, y no tiene STOP interno. Es decir: tradujimos
              una proteina que UniProt no reconoce para ese transcripto. Causas
              tipicas: el CDS del GTF y la entrada UniProt describen productos
              distintos (frecuente en entradas trEMBL, predichas por computadora
              y no curadas), o la anotacion del CDS esta desplazada. Mirar la
              columna 'fuente': sp=SwissProt (curada, un mismatch pesa) vs
              tr=trEMBL (predicha, un mismatch dice mas de UniProt que nuestro).

[B] SIN mapeo directo — el transcripto no tiene accession utilizable, asi que NO
    se puede validar posicion. Se cae al gen: se prueban las proteinas que UniProt
    tiene para ESE gen (aportadas por otros transcriptos). Es evidencia mas debil.
  hallada en    la secuencia inferida aparece dentro de alguna proteina del gen:
  otra prot.    el modulo TE existe en el producto real aunque este transcripto
  del gen       puntual no este mapeado. Sostiene el hit.
  producto      el gen tiene proteinas conocidas y la secuencia no aparece en
  distinto      ninguna: isoforma con un producto genuinamente distinto (tipico
                de las isoformas nuevas de Ensembl, sin equivalente en UniProt).
  accession     el header traia accession pero no esta en el FASTA descargado
  ausente       (FASTA incompleto o sin isoformas). Problema de datos, no del hit.
  sin UniProt   ni el transcripto ni el gen tienen proteina en UniProt: no hay
  para el gen   absolutamente nada contra que comparar.

DOS COLUMNAS PARA LEER 'fuente' SIN EQUIVOCARSE (en el TSV, no en este resumen):
  uniprot_usado        la accession con la que efectivamente se valido.
  fuente               si ESA es 'sp' (SwissProt, curada) o 'tr' (trEMBL, predicha).
  fuentes_disponibles  que habia para elegir en ese transcripto: 'sp', 'tr' o 'sp,tr'.
Se leen juntas. 'fuente=tr' con 'fuentes_disponibles=tr' significa que trEMBL era lo
unico que habia. 'fuente=tr' con 'fuentes_disponibles=sp,tr' es MUY distinto: habia una
entrada curada y NO fue la que coincidio -> mirar ese caso, la curada dio peor status.
Cuando ambas empatan gana siempre la curada, asi que ese cruce no aparece por azar.

SESGO IMPORTANTE: el bloque [A] existe solo para transcriptos que Ensembl ya habia
verificado identicos a UniProt. Los porcentajes de [A] miden si NUESTRAS coordenadas
y traduccion son correctas; NO son una estimacion insesgada sobre todos los hits.

Sobre sano=FALSE con 'exacto': 'sano' es un QC del CDS (largo multiplo de 3, ATG
inicial, STOP final). Un CDS puede fallar ese QC -por ejemplo empezar sin ATG por
estar 5' incompleto- y aun asi contener intacta la ventana del TE, que suele caer
lejos de los extremos. Por eso hay sano=FALSE que dan exacto: el defecto esta en el
borde del CDS, no en la region que nos interesa. Al reves, los sano=FALSE
concentran los 'stop' y los 'mismatch', que es lo que hace util al flag.""")

    w("=" * 66)
    w("VALIDACION secuencias inferidas (04a) vs UniProt real")
    w("=" * 66)
    w(f"04a     : {args.hits_04a}")
    w(f"uniprot : {args.uniprot_fasta} ({len(up)} secuencias)")
    reporte(tot, "TODOS")
    reporte(by_sano["TRUE"], "solo sano=TRUE")
    reporte(by_sano["FALSE"], "solo sano=FALSE")

    mismatches = [r for r in rows if r["status"] == "mismatch"]
    if mismatches:
        srcs = Counter(m["fuente"] for m in mismatches)
        sanos = Counter(m["sano"] for m in mismatches)
        w("\n" + "-" * 66)
        w(f"MISMATCH reales: {len(mismatches)}  |  fuente: "
          f"sp={srcs['sp']} tr={srcs['tr']}  |  sano: "
          f"TRUE={sanos['TRUE']} FALSE={sanos['FALSE']} NA={sanos['NA']}")
    glosario()

    def dump(path, cols, data):
        with open(path, "w") as fh:
            fh.write("\t".join(cols) + "\n")
            for r in data:
                fh.write("\t".join(str(r[c]) for c in cols) + "\n")
        sys.stderr.write(f"-> {path} ({len(data)} filas)\n")

    if args.out_tsv:
        dump(args.out_tsv, COLS, rows)
    if args.dump_mismatch and mismatches:
        dump(args.dump_mismatch, COLS, mismatches)
    if args.report:
        with open(args.report, "w") as fh:
            fh.write("\n".join(linea) + "\n")
        sys.stderr.write(f"-> {args.report}\n")


if __name__ == "__main__":
    main()
