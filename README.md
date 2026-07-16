# te-modular-homology

Pipeline reproducible y **parametrizado por especie** que, a partir de **hits de transposones
(Dfam nrph)** + **GTF** + **genoma**, encuentra los TE que aportan secuencia codificante y
extrae los **dominios proteicos TE-derivados** listos para **CD-Hit** (análisis de homología
modular mediada por transposones).

Todo se controla con un único `config.<especie>.sh`; no hay rutas fijas en el código. Se corre
de punta a punta con un solo comando, y escala a otros genomas copiando el config.

---

## Instalación / dependencias

- **R** con: `rtracklayer`, `GenomicRanges`, `GenomicFeatures`, `Biostrings`, `Rsamtools`,
  `data.table`, `tidyverse`, `janitor`.
- **Python 3** (solo stdlib) para la validación Dfam (paso 03, opcional; requiere red).
- *(opcional)* `cd-hit` para el paso 06; `biomaRt` (R) para el helper de UniProt.

No usa `TxDb`/`bedtools`/`yq`: el CDS-por-transcripto se arma directo del GTF; el config son
env-vars de shell (mapea 1:1 a un `config.yaml` de Snakemake si se estandariza después).

## Uso rápido

```bash
git clone <este-repo> && cd te-modular-homology
cp config/config.template.sh config/config.miEspecie.sh   # editar rutas/cromosomas
bash run.sh config/config.miEspecie.sh                     # núcleo 01→05
bash run.sh config/config.miEspecie.sh --dfam --cdhit      # + validación Dfam + CD-Hit
```

## Los pasos

| Paso | Script | Rol | Mundo |
|---|---|---|---|
| 01 | `intersect.R` | Hits Dfam ∩ CDS del GTF (+ liftOver flexible) → tabla de hits rica. | compartido |
| 02 | `extraer_cds.R` | CDS por transcripto + QC `sano` + región TE; agrega col `sano` a la tabla. | compartido |
| 03 | `validar_dfam.py` | *(opcional `--dfam`)* Doble validación contra el buscador de Dfam. | compartido |
| 04 | `clusterizar_loci.R` | Reducción de redundancia → `01_hits_magna, 02_posiciones, 03_loci`. | **genómico** |
| 05 | `extraer_dominios_aa.R` | Dominios AA: `04a_hits_aa` (crudo) + `04b_dominios_aa` (dedup contención). | **proteico** |
| 06 | `cdhit.sh` | *(opcional `--cdhit`)* CD-Hit 60/70/80 % → grupos de homología. | **proteico** |

> **Dos mundos:** el **genómico** (paso 04) cuenta *inserciones* (loci); el **proteico**
> (paso 05) arma *módulos* para homología (dominios). Detalle visual en `docs/`.

## Config (variables principales)

| Variable | Qué es |
|---|---|
| `GTF`, `TE_HITS`, `GENOME_FASTA` | Inputs de la especie (GTF Ensembl, hits Dfam nrph, genoma). |
| `DFAM_TAXA` | `repClass`/`repFamily`/clado por accession — **global de Dfam** (`helpers/parse_dfam_embl.sh`). |
| `CHROM_ALIAS` | *(opcional)* mapa seqname→UCSC, **solo si los hits usan accesiones GenBank**. Si usan `chr1`/`1`, vacío. |
| `ENST2UNIPROT` | *(opcional)* mapeo transcripto→uniprot (`helpers/map_uniprot.qmd`). Vacío = uniprot `NA`. |
| `CHROMS` | Cromosomas del GTF a conservar (crudos; se les antepone `chr`, `MT→M`). |
| `LIFTOVER_CHAIN` + `LIFTOVER_TARGET` | LiftOver flexible: `gtf` \| `hits` \| `none`. |
| `E_VALUE_THRESHOLD`, `OVERLAP_BP_THR` | Filtros del hit (default `1e-10`, `30` nt). |
| `DFAM_ORGANISM` | Organismo para la API de Dfam (paso 03). |
| `GENOME`, `OUTDIR` | Etiqueta corta y carpeta de resultados (`resultados_<GENOME>`, separada por especie). |

## Cómo agregar una especie

1. `cp config/config.template.sh config/config.<especie>.sh` y editar:
   - rutas de `GTF`/`TE_HITS`/`GENOME_FASTA` y la carpeta base (`DATA` → `OUTDIR`).
   - `CHROMS` (ratón `1-19`; grandes simios agregan `2A,2B`).
   - `LIFTOVER_*` (solo si GTF y hits difieren de ensamblado).
   - `CHROM_ALIAS` (solo si los hits usan accesiones GenBank).
2. *(opcional)* Generar el mapeo UniProt con `helpers/map_uniprot.qmd` (biomaRt) y apuntar
   `ENST2UNIPROT`.
3. `bash run.sh config/config.<especie>.sh`.

Cada especie escribe en su propia `resultados_<GENOME>/` — **nunca se mezclan**.

## Estructura de salidas (`<OUTDIR>/`)

```
hits_filtrados_en_cds.tsv   ← tabla central de hits (paso 02 le agrega la col sano)
metrics.tsv                  ← métricas de toda la corrida (paso/metrica/valor)
01_intersect/               ← all_cds_coord, hits_totales, transcriptos_excluidos_liftover
02_extraccion/              ← cds_qc, cds_completos(+sanos), regiones_te(+sanos)
03_validacion_dfam/         ← validacion_dfam_todos.tsv, dfam_cache/   (solo con --dfam)
clusters_<GENOME>/          ← 01_hits_magna, 02_posiciones, 03_loci, 04a_hits_aa, 04b_dominios_aa
```

## Ejemplo: chimpancé (baseline de no-regresión)

`config/config.panTro6.sh` es un ejemplo completo. Con los inputs de chimp (GTF Ensembl
`Pan_tro_3.0`, hits `DApanTro2.nrph`, genoma `panTro6`, chain `panTro5→panTro6`) el pipeline
reproduce exactamente:

```
01 intersect : 3697 hits filtrados (de 11992; excluye 2338 tx con liftOver imperfecto)
02 CDS + QC  : 2456 transcriptos sanos → 2842 hits sanos
04 clustering: 3137 posiciones → 3016 loci (2992 ok / 24 revisar)
05 dominios  : STEP1 3597 hits crudos (04a) → STEP2 3101 dominios distintos (04b)
```

> Los datos grandes (genoma, GTF, hits) **no** están en el repo — son inputs por especie. El
> README indica de dónde bajarlos; los resultados van a `resultados_<GENOME>/` (fuera de git).

## Documentación (`docs/`)

- **`flujo_pipeline.html`** — el flujo completo: los pasos, inputs/outputs, cómo se concatenan,
  y las decisiones de los dos mundos.
- **`redundancia_casos.html`** — los 9 casos de solape (genómico → proteico) + matriz.
- **`pipeline_overview.md`** — el "por qué": los cuatro niveles, las decisiones por mundo,
  Miedo A, el hallazgo de productos alternativos por splicing.
- **`walkthrough_decisiones.md`** — recorrido paso a paso con ejemplos reales y cada decisión.

## Helpers (`helpers/`)

- **`map_uniprot.qmd`** — plantilla Quarto (biomaRt) para generar el `ENST2UNIPROT` de cualquier
  especie (parametrizable por dataset y versión de Ensembl).
- **`parse_dfam_embl.sh`** — genera `dfam_curated_taxa.tsv` desde el flatfile EMBL de Dfam.
