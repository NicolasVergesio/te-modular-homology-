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

    def clasificar(uni, s, e, seq):
        if uni == "NA":
            return "na", None, "NA"
        cands = [a for a in re.split(r"[;,]", uni) if a]
        best = min(((clasificar_una(a, s, e, seq), a) for a in cands),
                   key=lambda x: RANK[x[0][0]])
        (out, base), acc = best
        return out, base, acc

    COLS = ["id", "id_pos", "locus", "gene", "tx", "te", "uniprot", "uniprot_usado",
            "fuente", "sano", "aa_start", "aa_end", "len_aa", "status", "seq_inferida"]
    tot = Counter(); by_sano = {"TRUE": Counter(), "FALSE": Counter(), "NA": Counter()}
    rows = []
    for h in parse_04a(args.hits_04a):
        out, base, acc = clasificar(h["uniprot"], h["aa_start"], h["aa_end"], h["seq"])
        tot[out] += 1
        by_sano.get(h["sano"], by_sano["NA"])[out] += 1
        h.update(status=out, uniprot_usado=acc, fuente=src.get(base, "NA") if base else "NA",
                 len_aa=len(h["seq"]), seq_inferida=h["seq"])
        rows.append(h)

    linea = []                                   # el reporte se acumula para poder guardarlo
    def w(s=""): linea.append(s); print(s)

    def reporte(c, titulo):
        comp = sum(c[k] for k in ("exacto", "substr", "stop", "otra_iso", "mismatch"))
        p = lambda n: f"{100*n/comp:5.1f}%" if comp else "  n/a"
        w(f"\n### {titulo}   (hits={sum(c.values())})")
        w(f"  sin uniprot en el header (NA)      : {c['na']}")
        w(f"  accession ausente del FASTA        : {c['no_acc']}")
        w(f"  COMPARABLES                        : {comp}")
        w(f"    exacto en coordenadas            : {c['exacto']:6d} ({p(c['exacto'])})")
        w(f"    substring (coord corrida)        : {c['substr']:6d} ({p(c['substr'])})")
        w(f"    STOP interno (traduccion rota)   : {c['stop']:6d} ({p(c['stop'])})")
        w(f"    rescatado en otra isoforma       : {c['otra_iso']:6d} ({p(c['otra_iso'])})")
        w(f"    MISMATCH real                    : {c['mismatch']:6d} ({p(c['mismatch'])})")

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
