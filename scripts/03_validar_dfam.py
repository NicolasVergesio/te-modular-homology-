#!/usr/bin/env python3
"""Doble validación de hits TE∩CDS contra el buscador de secuencias de Dfam.

Estrategia de dos niveles (por transcripto/hit):
  Nivel 1  -> se envía el CDS completo del transcripto (coords comparables a
              start_cds_rel/end_cds_rel). Un hit del pipeline se considera
              validado si Dfam devuelve un match de la MISMA repClass cuya
              posición en el CDS solapa con [start_cds_rel, end_cds_rel].
  Nivel 2  -> para los hits NO validados en el nivel 1, se reintenta con una
              ventana GENÓMICA (start_te..end_te ± FLANK) extraída de panTro6.fa;
              da más contexto a nhmmer y rescata TEs cortos/divergentes.

Notas:
  - El organismo "Pan troglodytes" da ERROR en el servidor de Dfam; se usa
    "Homo sapiens" (familias Alu/L1/MIR/L2/… compartidas entre primates).
  - Resumible: cada job se cachea en dfam_cache/ por su clave; re-ejecutar no
    reenvía lo ya resuelto.
  - Cortés con el servidor: envíos seriales con pausa y polling con backoff.

Uso:
  python3 validar_dfam.py            # piloto: 50 transcriptos sanos (semilla fija)
  python3 validar_dfam.py --all      # todos los transcriptos sanos
  python3 validar_dfam.py --n 100    # piloto de N transcriptos
"""
import argparse, csv, json, os, random, sys, time, urllib.parse, urllib.request

BASE   = os.path.dirname(os.path.abspath(__file__))
# Parametrizable por el config del pipeline (env-vars); fallback = dir del script.
DATA   = os.environ.get("OUTDIR", BASE)
TSV    = os.path.join(DATA, "hits_filtrados_en_cds.tsv")          # tabla central (raiz)
QC     = os.path.join(DATA, "02_extraccion", "cds_qc.tsv")
CDSFA  = os.path.join(DATA, "02_extraccion", "cds_completos.fasta")
GENOME = os.environ.get("GENOME_FASTA", os.path.join(BASE, "panTro6.fa"))
FAI    = GENOME + ".fai"
OUTV   = os.path.join(DATA, "03_validacion_dfam"); os.makedirs(OUTV, exist_ok=True)
CACHE  = os.path.join(OUTV, "dfam_cache")
API    = "https://dfam.org/api/searches"
ORG    = os.environ.get("DFAM_ORGANISM", "Homo sapiens")
FLANK  = 200          # nt de flanco a cada lado en la ventana genómica (nivel 2)

# ----------------------------- utilidades FASTA/FAIDX -----------------------
def read_fasta(path):
    d, name, buf = {}, None, []
    for line in open(path):
        line = line.rstrip("\n")
        if line.startswith(">"):
            if name: d[name] = "".join(buf)
            name, buf = line[1:].split()[0].split("|")[0], []   # 'tx|uniprot' -> solo el transcript_id
        else:
            buf.append(line)
    if name: d[name] = "".join(buf)
    return d

def load_fai(path):
    idx = {}
    for line in open(path):
        name, length, offset, linebases, linewidth = line.split()[:5]
        idx[name] = (int(length), int(offset), int(linebases), int(linewidth))
    return idx

_COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")
def revcomp(s): return s.translate(_COMP)[::-1]

def faidx_fetch(fh, fai, chrom, start, end):
    """Devuelve la secuencia genómica [start, end] 1-based inclusive, hebra +."""
    length, offset, lb, lw = fai[chrom]
    start = max(1, start); end = min(length, end)
    seq = []
    # byte de inicio de la base 'start'
    def byte_of(pos):
        return offset + (pos - 1) // lb * lw + (pos - 1) % lb
    fh.seek(byte_of(start))
    need = end - start + 1
    # leemos con margen por los saltos de línea
    raw = fh.read(need + need // lb + 2).decode()
    seq = raw.replace("\n", "")[:need]
    return seq

# ----------------------------- API de Dfam ----------------------------------
def _with_retry(fn, tries=5, base=5):
    """Ejecuta fn() reintentando ante fallos de red (DNS/conexión) con backoff
    exponencial: 5, 10, 20, 40, 80 s. Re-lanza la última excepción si se agotan."""
    last = None
    for k in range(tries):
        try:
            return fn()
        except Exception as e:
            last = e
            if k < tries - 1:
                time.sleep(base * (2 ** k))
    raise last

def _post(url, fields):
    data = urllib.parse.urlencode(fields).encode()
    req = urllib.request.Request(url, data=data, method="POST")
    with urllib.request.urlopen(req, timeout=60) as r:
        return json.load(r)

def _get(url):
    with urllib.request.urlopen(url, timeout=60) as r:
        return json.load(r)

def dfam_search(seq, cutoff="curated", organism=ORG, max_wait=240):
    """Envía una secuencia y espera el resultado. Devuelve el dict final de Dfam
    (con 'results') o {'status':'ERROR'/'TIMEOUT'}."""
    try:
        sub = _with_retry(lambda: _post(API, {"sequence": ">q\n" + seq,
                                              "organism": organism,
                                              "cutoff": cutoff}))
    except Exception as e:
        # Red caída de forma sostenida: no tiene sentido seguir. Salimos limpio
        # sin cachear nada; lo ya resuelto queda en dfam_cache/ y la corrida es
        # reanudable (re-ejecutar salta los jobs cacheados).
        print(f"\n[red] envío falló tras reintentos ({e}). Abortando de forma "
              f"limpia; reanudá cuando vuelva la red.", flush=True)
        sys.exit(1)
    jid = sub.get("id")
    if not jid:
        return {"status": "SUBMIT_FAIL", "raw": sub}
    t0 = time.time()
    while time.time() - t0 < max_wait:
        time.sleep(5)
        try:
            res = _get(f"{API}/{jid}")
        except Exception as e:
            time.sleep(4); continue
        if "results" in res:
            return res
        st = res.get("status")
        if st == "ERROR":
            return {"status": "ERROR", "raw": res}
    return {"status": "TIMEOUT", "id": jid}

def cached_search(key, seq, cutoff="curated", organism=ORG):
    path = os.path.join(CACHE, key + ".json")
    if os.path.exists(path):
        return json.load(open(path))
    res = dfam_search(seq, cutoff=cutoff, organism=organism)
    # Sólo cacheamos resultados exitosos (con 'results'); así un TIMEOUT/ERROR
    # transitorio NO queda congelado y se reintenta en la próxima corrida.
    if "results" in res:
        json.dump(res, open(path, "w"))
    return res

def parse_hits(res):
    """Aplana los hits de Dfam a: [{family,type,acc,evalue,qs,qe,strand}]."""
    out = []
    for r in res.get("results", []):
        for h in r.get("hits", []):
            a, b = int(h["ali_start"]), int(h["ali_end"])
            out.append(dict(family=h["query"], type=h.get("type", ""),
                            acc=h["accession"], evalue=h["e_value"],
                            qs=min(a, b), qe=max(a, b), strand=h["strand"]))
    return out

# ----------------------------- concordancia ---------------------------------
def overlap(a0, a1, b0, b1):
    return max(0, min(a1, b1) - max(a0, b0) + 1)

def norm_class(x):
    return (x or "").strip().lower()

def family_root(name):
    """Prefijo alfabético de un nombre de familia (AluSc8->alu, L1MB->l1, MIRb->mir)."""
    import re
    m = re.match(r"[A-Za-z]+[0-9]*", name or "")
    return (m.group(0) if m else name or "").lower()

def best_match(hits, p_start, p_end, p_class, p_name):
    """Elige el mejor hit de Dfam para un hit del pipeline: prioriza solape de
    posición y misma clase. Devuelve (dict_hit, ov_bp) o (None, 0)."""
    best, best_score = None, -1
    for h in hits:
        ov = overlap(p_start, p_end, h["qs"], h["qe"])
        same_class = norm_class(h["type"]) == norm_class(p_class)
        # puntaje: solape (bp) + bonus por misma clase para desempatar
        score = ov + (0.5 if same_class and ov > 0 else 0)
        if score > best_score:
            best, best_score = (h, ov), score
    if best and best[1] > 0:
        return best
    # sin solape: devolvemos igual el de misma clase si existe (family_only)
    for h in hits:
        if norm_class(h["type"]) == norm_class(p_class):
            return h, 0
    return (hits[0] if hits else None), 0

def concord_label(h, ov, p_class, p_name):
    if h is None:
        return "sin_hit"
    same_class = norm_class(h["type"]) == norm_class(p_class)
    same_fam   = family_root(h["family"]) == family_root(p_name)
    if ov > 0 and (same_class or same_fam):
        return "exact_family" if same_fam else "same_class_pos"
    if ov > 0:
        return "pos_only_diff_class"
    if same_class or same_fam:
        return "family_only_diff_pos"
    return "otro_hit"

# ----------------------------- programa principal ---------------------------
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--all", action="store_true", help="todos los tx sanos")
    ap.add_argument("--n", type=int, default=50, help="tamaño del piloto")
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--out", default=None)
    args = ap.parse_args()

    os.makedirs(CACHE, exist_ok=True)
    rows = list(csv.DictReader(open(TSV), delimiter="\t"))
    sano = {r["transcript_id"] for r in csv.DictReader(open(QC), delimiter="\t")
            if r["sano"] == "TRUE"}
    rows = [r for r in rows if r["transcript_id"] in sano]

    tx_all = sorted(set(r["transcript_id"] for r in rows))
    if args.all:
        tx_sel = set(tx_all)
        tag = "todos"
    else:
        random.seed(args.seed)
        tx_sel = set(random.sample(tx_all, min(args.n, len(tx_all))))
        tag = f"piloto{len(tx_sel)}"
    rows = [r for r in rows if r["transcript_id"] in tx_sel]
    out_path = args.out or os.path.join(OUTV, f"validacion_dfam_{tag}.tsv")

    cds = read_fasta(CDSFA)
    fai = load_fai(FAI)
    gfh = open(GENOME, "rb")

    # ---- Nivel 1: una búsqueda de CDS por transcripto ----
    by_tx = {}
    for r in rows:
        by_tx.setdefault(r["transcript_id"], []).append(r)
    txs = sorted(by_tx)
    print(f"[{tag}] {len(txs)} transcriptos, {len(rows)} hits", flush=True)

    cds_hits = {}     # tx -> lista de hits parseados de Dfam (CDS)
    t0 = time.time()
    for i, tx in enumerate(txs, 1):
        seq = cds.get(tx)
        if not seq:
            cds_hits[tx] = None; continue
        res = cached_search(f"cds_{tx}", seq)
        cds_hits[tx] = parse_hits(res) if "results" in res else None
        st = "ok" if cds_hits[tx] is not None else res.get("status", "?")
        print(f"  N1 {i}/{len(txs)} {tx}: {len(cds_hits[tx]) if cds_hits[tx] is not None else st} hits Dfam", flush=True)

    # ---- Evaluar cada hit del pipeline en nivel 1, y armar cola de nivel 2 ----
    results = []
    for r in rows:
        tx = r["transcript_id"]
        ps, pe = int(r["start_cds_rel"]), int(r["end_cds_rel"])
        pcl, pnm = r["repClass"], r["name_te"]
        hits = cds_hits.get(tx)
        rec = dict(r)  # copiamos campos del pipeline
        if hits is None:
            rec.update(n1_status="job_error")
            h, ov = None, 0
        else:
            h, ov = best_match(hits, ps, pe, pcl, pnm)
            rec.update(n1_status="ok")
        rec.update(n1_family=h["family"] if h else "", n1_type=h["type"] if h else "",
                   n1_acc=h["acc"] if h else "", n1_evalue=h["evalue"] if h else "",
                   n1_pos=f'{h["qs"]}-{h["qe"]}' if h else "", n1_overlap=ov,
                   n1_concord=concord_label(h, ov, pcl, pnm))
        results.append(rec)

    validated_n1 = sum(1 for r in results
                       if r["n1_concord"] in ("exact_family", "same_class_pos"))
    print(f"[N1] validados: {validated_n1}/{len(results)} hits", flush=True)

    # ---- Nivel 2: rescate genómico para hits NO validados en N1 ----
    for r in results:
        r.update(n2_family="", n2_type="", n2_acc="", n2_evalue="",
                 n2_pos="", n2_overlap="", n2_concord="", n2_status="")
        if r["n1_concord"] in ("exact_family", "same_class_pos"):
            continue
        chrom = r["seq_cds"]; s, e = int(r["start_te"]), int(r["end_te"])
        if chrom not in fai:
            r["n2_status"] = "sin_cromosoma"; continue
        ws, we = max(1, s - FLANK), e + FLANK
        win = faidx_fetch(gfh, fai, chrom, ws, we)
        key = f"geno_{chrom}_{ws}_{we}"
        res = cached_search(key, win)
        if "results" not in res:
            r["n2_status"] = res.get("status", "?"); continue
        r["n2_status"] = "ok"
        hits = parse_hits(res)
        # posición del TE del pipeline dentro de la ventana
        tp_s, tp_e = s - ws + 1, e - ws + 1
        h, ov = best_match(hits, tp_s, tp_e, r["repClass"], r["name_te"])
        r.update(n2_family=h["family"] if h else "", n2_type=h["type"] if h else "",
                 n2_acc=h["acc"] if h else "", n2_evalue=h["evalue"] if h else "",
                 n2_pos=f'{h["qs"]}-{h["qe"]}' if h else "", n2_overlap=ov,
                 n2_concord=concord_label(h, ov, r["repClass"], r["name_te"]))

    # ---- estado final por hit ----
    for r in results:
        n1 = r["n1_concord"]; n2 = r["n2_concord"]
        if n1 in ("exact_family", "same_class_pos"):
            r["validacion"] = "validado_cds"
        elif n2 in ("exact_family", "same_class_pos"):
            r["validacion"] = "validado_genomico"
        elif n1 == "family_only_diff_pos" or n2 == "family_only_diff_pos":
            r["validacion"] = "familia_ok_pos_distinta"
        elif r["n1_status"] == "job_error" and r["n2_status"] in ("", "?", "ERROR", "TIMEOUT"):
            r["validacion"] = "error_busqueda"
        else:
            r["validacion"] = "no_encontrado"

    # ---- salida ----
    cols = list(csv.DictReader(open(TSV), delimiter="\t").fieldnames)
    extra = ["n1_status","n1_family","n1_type","n1_acc","n1_evalue","n1_pos","n1_overlap","n1_concord",
             "n2_status","n2_family","n2_type","n2_acc","n2_evalue","n2_pos","n2_overlap","n2_concord",
             "validacion"]
    with open(out_path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols + extra, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for r in results:
            w.writerow(r)

    # ---- resumen ----
    from collections import Counter
    c = Counter(r["validacion"] for r in results)
    dt = time.time() - t0
    print(f"\n=== RESUMEN {tag} ({dt/60:.1f} min) ===")
    print(f"hits: {len(results)}  transcriptos: {len(txs)}")
    for k in ("validado_cds","validado_genomico","familia_ok_pos_distinta","no_encontrado","error_busqueda"):
        print(f"  {k:26s} {c.get(k,0)}")
    print(f"salida -> {out_path}")

    # ---- reporte a archivo (queda junto a la tabla, como el REPORTE_clusters.txt del paso 04) ----
    import datetime
    n_val = c.get("validado_cds", 0) + c.get("validado_genomico", 0)
    n_tot = len(results)
    rep_path = os.path.join(OUTV, "REPORTE_validacion_dfam.txt")
    L = [f"REPORTE DE VALIDACION DFAM  (tag: {tag})",
         f"Generado: {datetime.datetime.now():%Y-%m-%d %H:%M}  ({dt/60:.1f} min)",
         f"Entrada: {TSV}",
         f"Organismo (API Dfam): {ORG}",
         "=" * 60, "",
         "METODO: doble validacion contra el buscador de secuencias de Dfam.",
         "  Nivel 1 (CDS completo): se manda el CDS entero del transcripto.",
         f"  Nivel 2 (rescate genomico): ventana [start_te-{FLANK}, end_te+{FLANK}] del genoma",
         "     (el TE completo + flancos; el solape solo es demasiado corto para re-detectar).",
         "", "RESULTADO", "-" * 60,
         f"  hits validados : {n_tot}    transcriptos: {len(txs)}"]
    for k in ("validado_cds", "validado_genomico", "familia_ok_pos_distinta",
              "no_encontrado", "error_busqueda"):
        L.append(f"    {k:26s} {c.get(k,0)}")
    L.append(f"  -> {n_val}/{n_tot} ({100*n_val/n_tot:.1f}%) validados en algun nivel"
             if n_tot else "  -> sin hits")
    L += ["", "CAVEATS", "-" * 60,
          f"  - Organismo = {ORG} (proxy): 'Pan troglodytes' da ERROR en el servidor de Dfam;",
          "    las familias Alu/L1/MIR/L2/... son compartidas entre primates.",
          "  - Tabla SEPARADA: NO reescribe hits_filtrados_en_cds.tsv. Cruzar por join si hace falta.",
          "  - Paso OPCIONAL (QC): correr con `run.sh <config> --dfam`. Resumible (dfam_cache/)."]
    with open(rep_path, "w") as fh:
        fh.write("\n".join(L) + "\n")
    print(f"reporte -> {rep_path}")

if __name__ == "__main__":
    main()
