# te-modular-homology

Pipeline reproducible y **parametrizado por especie** que, a partir de **hits de transposones
(Dfam nrph)** + **GTF** + **genoma**, encuentra los TE que aportan secuencia codificante y
extrae los **dominios proteicos TE-derivados** listos para **CD-Hit** (análisis de homología
modular mediada por transposones).

Todo se controla con un único `config.<especie>.sh`; no hay rutas fijas en el código. Se corre
de punta a punta con un solo comando, y escala a otros genomas copiando el config.

---

## Instalación / dependencias

La forma recomendada: crear el entorno conda con todo incluido y verificarlo con el test
sintético (corre en segundos, sin datos grandes ni red).

```bash
mamba env create -f environment.yml      # o conda / micromamba
conda activate te-modular-homology
bash run.sh test/config.test.sh          # debe terminar en "LISTO"
```

El test arma un cromosoma de 2 kb con 1 gen y 1 TE, y debe dar **1 hit → 1 locus → 1 dominio de
51 alaninas**. Si eso sale, la instalación está bien.

Qué instala (por si preferís armarlo a mano):

- **R** con: `rtracklayer`, `GenomicRanges`, `GenomicFeatures`, `Biostrings`, `Rsamtools`,
  `GenomeInfoDb`, `data.table`, `tidyverse`, `janitor`.
- **Python 3** (solo stdlib) para la validación Dfam (paso 03, opcional; requiere red).
- *(opcional)* `cd-hit` para el paso 06; `biomaRt` (R) para el helper de UniProt.

Si ya tenés un R con los paquetes y no querés el env, apuntá `RSCRIPT` en el config al `Rscript`
de ese R.

No usa `TxDb`/`bedtools`/`yq` ni el binario `liftOver` de UCSC (el liftOver es el de
`rtracklayer`): el CDS-por-transcripto se arma directo del GTF; el config son env-vars de shell
(mapea 1:1 a un `config.yaml` de Snakemake si se estandariza después).

> **Quarto** no está en el `environment.yml`: sólo hace falta para renderizar
> `helpers/map_uniprot.qmd`, y conviene instalarlo con el instalador oficial o usar el que ya
> trae RStudio/Positron. Sin Quarto, se pueden correr los chunks del `.qmd` a mano.

## Qué hay que descargar

Sólo **tres archivos son imprescindibles**. El resto se agrega según el caso; el pipeline corre
igual sin ellos, perdiendo la funcionalidad que se indica.

| # | Archivo | ¿Obligatorio? | Variable | Para qué / qué se pierde sin él | Dónde se consigue |
|---|---|---|---|---|---|
| 1 | **GTF** de la especie | **SÍ** | `GTF` | La anotación de CDS: sin esto no hay nada que intersecar. | Ensembl FTP (`Especie.Ensamblado.RELEASE.gtf.gz`) |
| 2 | **Hits Dfam** formato `nrph` | **SÍ** | `TE_HITS` | Las posiciones de los TE en el genoma. | Dfam releases (`*.nrph.hits.gz`) |
| 3 | **Genoma FASTA** | **SÍ** | `GENOME_FASTA` | Extraer las secuencias. El `.fai` **no** se descarga: el paso 02 lo genera solo. | UCSC o Ensembl |
| 4 | `dfam_curated_taxa.tsv` | Recomendado | `DFAM_TAXA` | `repClass`/`repFamily`/clado por accession. **Sin él el pipeline no falla**: esas columnas quedan en `NA` y el paso 04 infiere la familia de la raíz del nombre del TE (`AluSx1` → `AluSx1`), que es peor pero funciona. | Se genera una vez con `helpers/parse_dfam_embl.sh` desde `Dfam_curatedonly.embl.gz`. **Es global de Dfam**: el mismo archivo sirve para todas las especies. |
| 5 | `.over.chain` | **Sólo si** los hits y el GTF están en **ensamblados distintos** | `LIFTOVER_CHAIN` + `LIFTOVER_TARGET` | Llevar unos u otros al ensamblado común. Si ambos ya están en el mismo, no se descarga nada. | UCSC (`liftOver/` del genoma) |
| 6 | `chromAlias.txt` | **Sólo si** los hits usan **accesiones GenBank** (`CM000377.3`) | `CHROM_ALIAS` | Traducir esos nombres a `chr1`. Si los hits ya vienen como `chr1` o `1`, no hace falta: el paso 01 normaliza solo. | UCSC (`<genoma>.chromAlias.txt`) |
| 7 | Mapeo transcripto→UniProt | Opcional | `ENST2UNIPROT` | Llena la columna `uniprot`. Sin él queda `NA` en toda la corrida (no rompe nada aguas abajo). | Ensembl FTP (`tsv/…uniprot.tsv.gz`) o `helpers/map_uniprot.qmd` (biomaRt) |
| 8 | FASTA de UniProt **con isoformas** | Opcional | *(ninguna)* | Sólo para `helpers/validar_vs_uniprot.py`, que se corre aparte. No participa del pipeline. | UniProt (proteoma, "all isoforms") |

> **Cuidado con el par 5–6:** son problemas distintos y se confunden. El `.chain` es para
> **coordenadas** de ensamblados distintos (panTro5 vs panTro6); el `chromAlias` es para
> **nombres** de cromosomas del mismo ensamblado. Podés necesitar uno, otro, los dos o ninguno.

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

## Config: qué hacer con lo que no usás

**Regla: dejá la variable con la ruta vacía (`=""`). Nunca comentes la línea.**

```bash
export LIFTOVER_CHAIN=""        # ✅ correcto: no tengo chain
#export LIFTOVER_CHAIN="..."    # ❌ evitar (ver abajo)
```

Vacío y comentado *parecen* lo mismo, pero no lo son. Las variables que leen R y Python dan igual
en los dos casos (`Sys.getenv()` devuelve `""` tanto si está vacía como si no existe). Pero
`run.sh` corre con `set -u` y usa algunas directamente, así que **comentarlas aborta la corrida**:

```
$ bash run.sh config/config.miEspecie.sh        # con LIFTOVER_TARGET comentada
run.sh: line 33: LIFTOVER_TARGET: unbound variable
```

Con `LIFTOVER_TARGET=""` la misma corrida termina bien. Como la distinción depende de qué variable
es —y no es evidente cuál lee quién—, la regla simple de usar siempre `""` es la que conviene:
funciona para todas. El `config.template.sh` ya viene así.

Casos concretos de "no lo tengo":

| Situación | Qué poner |
|---|---|
| Mismo ensamblado, sin liftOver | `LIFTOVER_CHAIN=""` y `LIFTOVER_TARGET="none"` |
| Los hits ya usan `chr1` o `1` | `CHROM_ALIAS=""` |
| Sin mapeo a UniProt | `ENST2UNIPROT=""` (la columna `uniprot` queda `NA`) |
| Sin taxonomía de Dfam | `DFAM_TAXA=""` (familia inferida del nombre del TE) |
| Todos los cromosomas del GTF | `CHROMS=""` |

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
hits_filtrados_en_cds.tsv    ← tabla central de hits
metrics.tsv                  ← métricas de toda la corrida (paso/metrica/valor)
01_intersect/                ← all_cds_coord, hits_totales_en_cds, transcriptos_excluidos_liftover
02_extraccion/               ← cds_qc, cds_completos(+sanos), regiones_te(+sanos)
03_validacion_dfam/          ← validacion_dfam_todos.tsv, dfam_cache/   (solo con --dfam)
clusters_<GENOME>/           ← 01_hits_magna … 04b, validación, REPORTE
```

### Qué es cada archivo

Son muchos; **los que se usan habitualmente están marcados con ★**. El resto es trazabilidad
(poder auditar de dónde salió cada número) o insumo intermedio.

**Raíz de `<OUTDIR>/`**

| archivo | una fila = | qué es |
|---|---|---|
| ★ `hits_filtrados_en_cds.tsv` | un hit × exón × isoforma | **La tabla central.** Cada TE que cae en un CDS, con sus coordenadas genómicas, las relativas al CDS, las aminoacídicas (`start_aa`/`end_aa`), la anotación Dfam (`repClass`/`repFamily`/`class_family`), UniProt y el flag `sano` que agrega el paso 02. Todo lo demás deriva de acá. |
| `metrics.tsv` | una métrica | `paso · metrica · valor` de toda la corrida. Sirve para comparar corridas entre sí y detectar regresiones de un vistazo. |

**`01_intersect/`** — el cruce crudo, antes de filtrar

| archivo | una fila = | qué es |
|---|---|---|
| `all_cds_coord.tsv` | un exón CDS del GTF | Todos los CDS `protein_coding` con sus coordenadas (ya lifteadas si correspondía). Es el universo contra el que se cruzan los TE; pesado y sólo para auditar. |
| `hits_totales_en_cds.tsv` | un solape TE∩CDS | Todos los solapes **antes** de aplicar `E_VALUE_THRESHOLD` y `OVERLAP_BP_THR`. Comparado con `hits_filtrados_en_cds.tsv` muestra exactamente qué descartó el filtro. |
| `transcriptos_excluidos_liftover.tsv` | un transcripto excluido | Sólo con `LIFTOVER_TARGET="gtf"`: transcriptos que el liftOver rompió, con el motivo (`roto_chrom_hebra`, `multimapeo`, `perdida_parcial`). Vacío = liftOver limpio. |

**`02_extraccion/`** — secuencias y QC

| archivo | una entrada = | qué es |
|---|---|---|
| ★ `cds_qc.tsv` | un transcripto | El QC que define `sano`: `width`, `multiplo_3`, `empieza_atg`, `stop_codon`, `stop_ok`. Acá se ve **por qué** un transcripto quedó marcado como no sano. |
| `cds_completos.fasta` | un transcripto | CDS nucleotídico completo reconstruido uniendo exones. Es el insumo del paso 05 (se traduce acá). |
| `cds_completos_sanos.fasta` | un transcripto | Igual, restringido a los `sano=TRUE`. |
| `regiones_te.fasta` | un hit | Sólo el tramo nucleotídico del TE dentro del CDS (no la proteína). Útil para análisis a nivel ADN. |
| `regiones_te_sanos.fasta` | un hit | Igual, restringido a transcriptos sanos. |

**`clusters_<GENOME>/`** — reducción de redundancia (paso 04) y dominios proteicos (paso 05)

| archivo | una fila = | qué es |
|---|---|---|
| `01_hits_magna.tsv` | un hit × exón × isoforma | La tabla central **más** las columnas de redundancia: `id_pos` (posición genómica), `id_locus` (inserción), `rol` (representante o no), `flag`. No borra nada: anota. |
| `02_posiciones.tsv` | una **posición genómica** de TE | Colapsa las isoformas: varias filas magna del mismo TE en el mismo lugar → una posición. |
| `03_loci.tsv` | una **inserción** (locus) | Colapsa posiciones que son el mismo evento de inserción (SVA que contiene un Alu, L1 fragmentado, etc.). **Es la unidad "genómica": contar inserciones se hace acá.** |
| `REPORTE_clusters.txt` | — | Texto legible con la cascada `filas magna → posiciones → loci` y los casos que quedaron marcados para revisar. Es la explicación de los tres archivos anteriores. |
| `04a_hits_aa.fasta` | un hit (posición × isoforma) | Secuencia AA de **todos** los hits, crudo y redundante. Respaldo y entrada del validador. |
| ★ `04b_dominios_aa.fasta` | un dominio | **La salida principal: el input de CD-Hit.** Dentro de cada locus se deduplica por contención de secuencia. |
| ★ `04b_dominios_aa.tsv` | un dominio | La misma información en tabla, con los campos multi-valor alineados por transcripto (ver *Campos multi-valor*, abajo). |
| ★ `04b_productos_tx.tsv` | un hit (**combinación inserción × transcripto**) | Versión **larga** de la anterior: sin campos multi-valor, una fila por combinación. Es la cómoda para analizar en R y la única donde TE y transcripto pueden variar a la vez sin ambigüedad. |
| `05_validacion_uniprot.tsv` | un hit | Sólo si se corrió `helpers/validar_vs_uniprot.py`: el `status` de cada hit contra UniProt, con qué accession validó y cuáles había. |
| `05_validacion_uniprot.txt` | — | El reporte de porcentajes de lo anterior, con un glosario de cada categoría al final. |

**`cdhit/`** (sólo con `--cdhit`): `dominios_c0.6/0.7/0.8.fasta` y sus `.clstr` — los grupos de
homología a tres umbrales de identidad.

### Las dos unidades que no hay que confundir

- **Inserción / locus** (`03_loci.tsv`) → mundo **genómico**: "cuántos TE aportaron secuencia
  codificante".
- **Dominio** (`04b_dominios_aa`) → mundo **proteico**: "cuántos módulos proteicos distintos hay".

Un locus puede dar más de un dominio (isoformas con productos distintos) y varios hits colapsan
en un dominio. Los números **no** tienen por qué coincidir, y comparar uno contra otro sin aclarar
cuál se está usando es la confusión más fácil de cometer.

### Campos multi-valor (en `04b_dominios_aa`)

Un dominio puede estar sostenido por varios transcriptos, y cada uno tener varias accessions:

| separador | significa |
|---|---|
| `;` | separa **transcriptos**: siempre `n_tx` slots, en el orden de la columna `tx` |
| `,` | separa varios valores **dentro** de un transcripto |
| `NA` | ocupa el slot vacío, nunca se omite |

```
tx      = ENSPTRT00000001312;ENSPTRT00000077350
uniprot = H2PYY3,A0A6D2VWW0;K7BA23,A0A6D2XP50
```

Aplica a `uniprot`, `gene` y `sano`. **`tes` no**: los TEs son otro eje. Si necesitás la
correspondencia exacta sin decodificar nada, usá `04b_productos_tx.tsv`.

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
