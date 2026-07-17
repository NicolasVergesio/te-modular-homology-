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

Uso:
  validar_vs_uniprot.py UNIPROT.fasta[.gz] 04a_hits_aa.fasta [--dump-mismatch out.tsv]
"""
import argparse, gzip, re, sys
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
    """Itera (uniprot, aa_start, aa_end, seq, sano, tx, te) por hit."""
    with open(path) as fh:
        hdr = None
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                hdr = line
            elif hdr is not None:
                def g(pat, d=None):
                    m = re.search(pat, hdr); return m.group(1) if m else d
                yield (g(r"uniprot=(\S+)", "NA"),
                       int(g(r"aa=(\d+)-\d+", "0")), int(g(r"aa=\d+-(\d+)", "0")),
                       line, g(r"sano=(\S+)", "NA"), g(r"tx=(\S+)"), g(r"te=(\S+)"))
                hdr = None


def main():
    ap = argparse.ArgumentParser(description="Valida 04a_hits_aa.fasta vs UniProt.")
    ap.add_argument("uniprot_fasta", help="FASTA UniProt (con isoformas), .gz ok")
    ap.add_argument("hits_04a", help="04a_hits_aa.fasta del pipeline")
    ap.add_argument("--dump-mismatch", metavar="TSV",
                    help="volcar los mismatch reales a un TSV para inspeccion")
    args = ap.parse_args()

    up, src = load_uniprot(args.uniprot_fasta)
    sys.stderr.write(f"uniprot: {len(up)} secuencias\n")
    base_index = {}
    for a, s in up.items():
        base_index.setdefault(re.sub(r"-\d+$", "", a), []).append(s)

    def clasificar(uni, s, e, seq):
        if uni == "NA":
            return "na", None
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

    OUT = ["na", "no_acc", "exacto", "substr", "stop", "otra_iso", "mismatch"]
    tot = Counter(); by_sano = {"TRUE": Counter(), "FALSE": Counter(), "NA": Counter()}
    mismatches = []
    for uni, s, e, seq, sano, tx, te in parse_04a(args.hits_04a):
        out, base = clasificar(uni, s, e, seq)
        tot[out] += 1
        by_sano.get(sano, by_sano["NA"])[out] += 1
        if out == "mismatch":
            mismatches.append((uni, src.get(base, "?"), s, e, sano, tx, te, seq))

    def reporte(c, titulo):
        comp = sum(c[k] for k in ("exacto", "substr", "stop", "otra_iso", "mismatch"))
        p = lambda n: f"{100*n/comp:5.1f}%" if comp else "  n/a"
        print(f"\n### {titulo}   (hits={sum(c.values())})")
        print(f"  no verificable (NA/acc no hallada) : {c['na'] + c['no_acc']}")
        print(f"  COMPARABLES                        : {comp}")
        print(f"    exacto en coordenadas            : {c['exacto']:6d} ({p(c['exacto'])})")
        print(f"    substring (coord corrida)        : {c['substr']:6d} ({p(c['substr'])})")
        print(f"    STOP interno (traduccion rota)   : {c['stop']:6d} ({p(c['stop'])})")
        print(f"    rescatado en otra isoforma       : {c['otra_iso']:6d} ({p(c['otra_iso'])})")
        print(f"    MISMATCH real                    : {c['mismatch']:6d} ({p(c['mismatch'])})")

    print("=" * 66)
    print("VALIDACION secuencias inferidas (04a) vs UniProt real")
    print("=" * 66)
    reporte(tot, "TODOS")
    reporte(by_sano["TRUE"], "solo sano=TRUE")
    reporte(by_sano["FALSE"], "solo sano=FALSE")

    if mismatches:
        srcs = Counter(m[1] for m in mismatches)
        sanos = Counter(m[4] for m in mismatches)
        print("\n" + "-" * 66)
        print(f"MISMATCH reales: {len(mismatches)}  |  fuente: "
              f"sp={srcs['sp']} tr={srcs['tr']}  |  sano: "
              f"TRUE={sanos['TRUE']} FALSE={sanos['FALSE']} NA={sanos['NA']}")

    if args.dump_mismatch and mismatches:
        with open(args.dump_mismatch, "w") as fh:
            fh.write("uniprot\tfuente\taa_start\taa_end\tsano\ttx\tte\tseq_inferida\n")
            for m in mismatches:
                fh.write("\t".join(map(str, m)) + "\n")
        sys.stderr.write(f"mismatches -> {args.dump_mismatch}\n")


if __name__ == "__main__":
    main()
