# Overview del pipeline TE → homología modular

Documento de referencia: qué hace el pipeline, cómo se corre, y las **decisiones
conceptuales** que hay detrás (el "por qué", que no se ve leyendo el código). Complementa
al `README.md` (que es la referencia técnica de variables/pasos).

---

## Qué hace, en una frase

A partir de **hits de transposones (Dfam nrph) + GTF + genoma**, encuentra los TE que aportan
secuencia codificante, reduce la redundancia a nivel de inserción, y extrae los dominios
proteicos TE-derivados **listos para CD-Hit** (análisis de homología modular).

## Cómo se corre (un comando)

```bash
pipeline/run.sh pipeline/config/config.panTro6.sh            # núcleo 01→05
pipeline/run.sh pipeline/config/config.panTro6.sh --dfam     # + validación Dfam (opcional)
pipeline/run.sh pipeline/config/config.panTro6.sh --cdhit    # + CD-Hit (opcional)
```

## Baseline chimpancé (no-regresión — el pipeline debe reproducir esto)

```
01 intersect : 3697 hits filtrados (de 11992; excluye 2338 tx con liftOver imperfecto)
02 CDS + QC  : 2456 transcriptos sanos → 2842 hits sanos ; 3697 regiones TE
06 clustering: 3137 posiciones genómicas → 3016 loci (2992 ok / 24 revisar)
07 dominios  : STEP1 3597 hits AA crudos (04a) → STEP2 3101 dominios distintos (04b, dedup por contención)
```

> **Decisión (2026-07-14):** el paso 01 excluye los transcriptos cuyo CDS no lifteó de forma
> consistente (nº de exones out≠in, o cae en varios chr/hebra) — motivos `multimapeo` /
> `perdida_parcial` / `roto_chrom_hebra`, anotados en `transcriptos_excluidos_liftover.tsv`.
> Esto bajó los productos con STOP interno de 50 a 9 (esos tx rompían el marco de lectura).
> Antes eran 3915→3156→3284; ahora 3697→3016→3101. Métricas de cada paso en `metrics.tsv`.
> El filtro por clado NO se aplica (todos sobreviven); la composición se reporta en `metrics.tsv`.

---

## El modelo mental: cuatro niveles, cada uno para una pregunta distinta

Todo sale de **una tabla maestra** (`01_hits_magna.tsv`, 3915 filas). Los demás números son
vistas/subconjuntos de ella. Confundir los niveles es la fuente de confusión #1:

```
3697  filas magna (hit × isoforma × trozo de exón)   ← trazabilidad
 3137  posiciones genómicas (dedup coords+cadena+subfamilia)   ← "el TE físico"
  3016  loci / inserciones (reglas familia-aware)   ← ¿CUÁNTAS inserciones? (para contar)
   3101  dominios AA distintos (dedup por contención intra-locus)   ← input de CD-Hit (para homología)
```

- **Contar inserciones** → nivel **locus** (3016).
- **Input de CD-Hit / homología** → nivel **dominio proteico distinto** (3101).
- Cada fila arrastra su `id_locus` y su `sano`, así que se sube y baja entre niveles sin perder el hilo.

## Las decisiones del algoritmo, por mundo

El detalle de qué decide el pipeline en cada mundo. (Visual paso a paso en
`redundancia_casos.html` y `flujo_pipeline.html`.)

**Compartido — pasos 01-02 (calidad y coordenadas):**
- Solo CDS de genes/transcriptos `protein_coding` (fuera UTR, lncRNA, pseudogenes, IG/TR).
- Unificar cromosomas al estilo UCSC + **liftOver flexible** (`gtf|hits|none`); **excluir y anotar** los
  transcriptos con liftOver imperfecto (nº exones out≠in: `multimapeo`/`perdida_parcial`/`roto`).
- Proyección genómico→aminoácido con `pmapToTranscripts` (frame-aware); chequeo de subcadena = 100%.
- Filtros del hit: `e-value ≤ 1e-10` + `coverage ≥ 30 nt`.
- Bandera de calidad `sano` (múltiplo de 3 + ATG + stop en marco) — NO filtra, marca. **Sin filtro por clado.**

**Mundo GENÓMICO — paso 04 (objetivo: contar inserciones):**
- Dedup a **posiciones** por `(coords + cadena + subfamilia)`.
- Cluster inicial por solape (single-linkage), luego reglas **familia-aware**:
  - misma subfamilia fragmentada → 1 locus (span-unión ≤ p90);
  - misma subfam. flanqueando otra (X-Y-X) → piezas = 1 locus (`split`), lo del medio aparte;
  - misma familia distinta subfam. → 1 locus si recíproco ≥50% (`ambiguedad`); si <50% → tándem separado;
  - **SVA contiene Alu** (≥80%, cadena opuesta) → `es_subparte_sva`, 1 locus (SVA);
  - familias distintas anidadas/adyacentes → se **conservan** (loci distintos).
- Representante = menor e-value (en grupos con SVA, gana el SVA). → **loci = inserciones (3016)**.

**Mundo PROTEICO — paso 05 (objetivo: módulos para homología):**
- STEP 1: extraer por **(posición × isoforma)**, uniendo trozos de exón. Crudo → `04a`.
- STEP 2: dentro de cada locus, **dedup por contención de secuencia** (subcadena → queda la larga;
  header lista todos los TEs). Sin mirar familia (el paso 04 ya separó los módulos en loci distintos).
- Conservar secuencias **distintas** (isoformas divergentes, tipo MER41). **Nunca** dedup entre loci
  (eso es la señal de homología). → **dominios distintos (3101)** = input de CD-Hit (`04b`).

## Por qué la reducción de redundancia importa (el "Miedo A")

Una misma inserción física aparece en varias isoformas. Si se mandan las 3915 crudas a CD-Hit,
una sola inserción se cuenta muchas veces e **infla los clusters con falsa diversidad** — justo
el criterio "clusters con ≥2 secuencias distintas" del método del director es vulnerable a eso,
y su filtro de parálogos (a nivel gen) NO lo cubre (son isoformas del mismo gen). Por eso se
deduplica **antes**. Ejemplo real: `LTR105_Mam` cae en 8 filas / 3 UniProt de un mismo gen.

## Qué es "un hit final" y qué NO es homología modular

- **Hit final depurado = una fila de `03_loci.tsv`** (una inserción). El representante da la
  etiqueta+e_value oficial; la extensión es el span-unión de las piezas. `estado=ok` es el núcleo
  limpio (los `revisar` son casos SVA/Alu dudosos).
- **Homología modular = módulos proteicos TE-derivados SIMILARES entre proteínas DISTINTAS**
  (idealmente genes distintos), detectada por secuencia (CD-Hit). Dos productos del **mismo gen**
  (isoformas) **NO** son homología modular entre sí — es splicing alternativo, y es lo que el
  pipeline del director filtra a propósito.

## Hallazgo clave: una inserción puede dar productos proteicos DISTINTOS

Verificado con datos reales traduciendo el CDS de chimp:
- `LTR105_Mam`: 3 isoformas → AA **idénticas** (redundancia real → colapsan a 1 producto).
- `MER41-int` (locus L02594): 2 isoformas de un mismo gen, misma inserción, proteína idéntica
  hasta AA 389 y a partir del TE **dos C-terminales totalmente distintos** (`VELNLDCHCENAKPW…`
  vs `GPPAKFVADQLAGSS`) — el TE aporta **exones alternativos** por splicing. NO es redundancia.

Por eso el paso 05 dedup por **contención de secuencia dentro del mismo locus** (una AA subcadena
de otra → queda la larga; el header lista los TEs): colapsa idénticas y contenidas (LTR105,
Alu-en-SVA, pilas de Alu) y **conserva las distintas** (MER41). Nunca se deduplica entre loci
distintos (eso ES la señal para CD-Hit). Resultado: **82 loci con >1 dominio distinto**.

## `sano` + `flag_stop` (calidad de secuencia)

- `sano` (QC del CDS: múltiplo de 3 + ATG + stop en marco) nace en el paso 02, se propaga a
  hits/loci/productos.
- `flag_stop` (stop interno en el AA extraído) marca traducciones sospechosas.
- Cruzados: **2407 dominios limpios** (sano+sin stop) de 3101; los con STOP interno bajaron a 8
  (la exclusión del liftOver limpió los frameshifts). Se filtran a gusto antes del CD-Hit.

## LiftOver flexible

`LIFTOVER_TARGET = gtf | hits | none`. Chimp usa `gtf` (Ensembl→panTro6). Si otra especie tiene
los hits en otro ensamblado, `hits`; si coinciden, `none`.

## Multi-especie

Copiar `config/config.panTro6.sh`, ajustar `GTF/TE_HITS/GENOME_FASTA/CHROMS/CHROM_ALIAS/
LIFTOVER_*`, generar el mapeo UniProt (bloque BioMart en `pan_tro/_legacy/pan_tro.qmd`, cambiando
el dataset), y correr. Ver README, sección "Cómo agregar una especie".

## Pendientes

- `06_cdhit.sh`: wrapper de CD-Hit 60/70/80 % (hoy `--cdhit` apunta a un script que falta).
- JOIN de la validación Dfam (`validacion_dfam_todos.tsv`) sobre la tabla de clusters, para
  reportar "inserciones validadas" (la validación corrió sobre el subconjunto sano = 2906).
