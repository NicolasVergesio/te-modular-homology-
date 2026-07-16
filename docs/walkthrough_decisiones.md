# Caminata por el pipeline — registro paso a paso + decisiones

Registro de la recorrida "estación por estación" del pipeline (con ejemplos reales) y de las
decisiones que se fueron tomando. Complementa a `pipeline_overview.md` (el resumen) y al
`README.md` (la referencia técnica). Se va actualizando a medida que avanzamos.

---

## Estación 0 — Los insumos

Todo arranca cruzando **dos tablas**:
- **Hits de TE (Dfam nrph)**: una línea = una copia física de un transposón en el genoma.
  Ej: `CM009238.2  DF000000051  AluSz  ... + 2919 3167 ... 9.52`. Cromosoma = accesión GenBank.
  `nrph` = non-redundant per hit (Dfam ya eligió la mejor familia por región).
- **GTF (Ensembl)**: una línea CDS = un exón codificante de un transcripto.
  Ej: `1 ensembl CDS 139599211 139599225 - 0 gene_id "ENSPTRG..." transcript_id "ENSPTRT..."`.
  Cromosoma = número pelado (`1`).

**Pregunta que resuelve el pipeline:** ¿dónde un TE cae dentro de un CDS? = solape de coords.

**Filtro de entrada (paso 01):** solo `type==CDS` + `gene_biotype==protein_coding` +
`transcript_biotype==protein_coding`. Se van UTRs, exones no codificantes, lncRNA, pseudogenes,
y los 256 segmentos IG/TR (inmunoglobulinas / receptores T). Se estudia TE→proteína, así que solo
la parte traducida de transcriptos codificantes reales.

## Estación 1 — El cruce TE ∩ CDS (paso 01)

**1A — Nombres de cromosoma.** El mismo cromosoma tiene 3 nombres: hits `CM009238.2` (GenBank),
GTF `1` (Ensembl), objetivo `chr1` (UCSC). Se unifican a UCSC: al GTF se le antepone `chr` (y
`MT→M`); a los hits se los traduce con la tabla `chromAlias` (GenBank→UCSC). `2A`/`2B` = el chr2
humano está partido en dos en grandes simios. Se descartan scaffolds (solo `CHROMS` principales).

**1B — LiftOver (versiones de ensamblado).** GTF y hits están en builds distintos del genoma;
las coordenadas se corren entre versiones. Un archivo `chain` mapea bloque a bloque; `liftOver`
traduce. Decisión flexible `LIFTOVER_TARGET = gtf|hits|none` (chimp mueve el GTF a panTro6).
Consecuencias: mapea 98,82% de los CDS; y algunos transcriptos quedan mal (ver decisión abajo).
Seguro incorporado: la región TE extraída es subcadena EXACTA del CDS reconstruido (100%).

**1C — Tabla rica, proyección a aa, filtros.** Por cada solape se guardan tres bloques
(gen/CDS, TE completo, solape) — clave: se guarda el **TE completo** (start_te/end_te) y el
**solape** por separado. La posición en aminoácidos se obtiene con `pmapToTranscripts`
(frame-aware: junta exones 5'→3', respeta hebra) y `aa = ceil(nt/3)` — POST-liftOver, sobre el
CDS lifteado. Ej real: un `AluY` = los primeros 18 aa de una proteína. Dos filtros: `e-value ≤
1e-10` (confianza) y `coverage ≥ 30 nt` (10 codones). 11992 solapes → 3697 filtrados.

### Decisiones tomadas en la Estación 1 (2026-07-14)

1. **Exclusión estricta del liftOver.** Un tx solo pasa si `nº exones CDS out == in`, en un solo
   chr y una sola hebra. Se excluyen y ANOTAN (`transcriptos_excluidos_liftover.tsv`) los
   `multimapeo` (out>in, exón duplicado, 1064), `perdida_parcial` (out<in, 1025) y
   `roto_chrom_hebra` (249). Antes solo se sacaban los 249 rotos. **Efecto:** bajó los productos
   con STOP interno de 50 → 9 (esos tx rompían el marco). Baseline 3915→3156→3284 pasó a
   3697→3016→3132.
2. **`metrics.tsv`** por corrida (long: paso/metrica/valor) — cascada completa + composición por
   clado. Los 4 scripts appendan; el paso 01 lo reinicia.
3. **Taxonomía Dfam**: si falta el archivo, ahora crea las columnas en `NA` (antes no las creaba
   y rompía el 06).
4. **Filtro por clado: NO se aplica** — todos sobreviven hasta el final; la composición se
   informa en `metrics.tsv` (Eutheria/Primates/Mammalia/... — 43% son TEs antiguos compartidos).

### Preguntas de la Estación 1 aún abiertas / a decidir
- (ninguna pendiente por ahora)

## Estación 2 — QC del "CDS sano" (paso 02)

Un transcripto es **sano** si su CDS reconstruido está bien formado. Tres chequeos (todos deben
dar bien):
1. **multiplo_3**: largo divisible por 3 (se lee en codones). Ej falla: `ENSPTRT00000005424`
   width 2167 (2167/3 no entero).
2. **empieza_atg**: arranca en `ATG`. Ej falla: `ENSPTRT00000002864`.
3. **stop_ok**: termina en stop en marco (`TAA/TAG/TGA`). Ej falla: `ENSPTRT00000020613` (GTG).

**Detalle del stop:** en Ensembl el codón de stop es una feature APARTE, no está en el CDS. Por
eso no se mira el final del CDS: se leen del genoma los 3 nt inmediatamente 3' del último nt
codificante (respetando la hebra) y se chequea que sean stop.

**Nomenclatura de las métricas (aclarado 2026-07-14):** `tx_*` cuenta TRANSCRIPTOS, `hits_*`
cuenta HITS. Un transcripto puede tener varios hits; al marcarlo sano/no-sano, todos sus hits
heredan el estado.
```
  3697 hits ──(por transcripto)──► 3174 tx distintos  (tx_con_hits)
                                    ├─ 2456 sano   (tx_sanos)    → 2842 hits (hits_sanos)
                                    └─  718 no-sano (tx_no_sanos) →  855 hits (hits_no_sanos)
```
(Se renombró `cds_completos` → `tx_con_hits` para que el prefijo `tx_` sea siempre "cuenta
transcriptos".)

**`sano` es una BANDERA de calidad, no un filtro** — los no-sanos se conservan, anotados. Sirve
para el paso 05: si el marco está roto, la traducción a proteína puede salir con stops internos.
Es el compañero de `flag_stop`. (Fue lo que bajó de 50→9 al excluir los tx mal lifteados.)

### Decisión de organización de salidas (2026-07-14)

Se separaron **datos** de **resultados**. `OUTDIR` pasó a ser una carpeta dedicada
`${EP}/resultados_<GENOME>/` (antes era `${EP}`, que mezclaba genoma/GTF con outputs). Adentro:
```
resultados_<GENOME>/
├── hits_filtrados_en_cds.tsv   ← tabla central (raíz, accesible)
├── metrics.tsv                  ← dashboard (raíz)
├── 01_intersect/               ← all_cds_coord, hits_totales, transcriptos_excluidos
├── 02_extraccion/              ← cds_qc, cds_completos*, regiones_te*
├── 03_validacion_dfam/         ← dfam_cache/, validacion_dfam_todos.tsv (paso 03 opcional)
└── clusters_<GENOME>/          ← 01_hits_magna .. 04_dominios_aa + REPORTE
```
El dir de datos (`ensembl_pipe/`) queda limpio (solo genoma, GTF, chain, alias, uniprot map).
Escala a otras especies: cada una su `resultados_<especie>/`.

## Estación 3 — Extracción de la región del TE (paso 02)

Se extrae el **solape CDS∩TE en nucleótidos, en la hebra del gen** (`regiones_te.fasta`).
Ejemplo real (`L2c_3end`, `chr1:55788638-55788673`, gen en `−`):
- `+` del genoma (guardado en el FASTA): `CTGAGAAAGCAGCATGAAGAAACGCCCAGAGAGAGT`
- extraído (hebra `−` = reverse-complement): `ACTCTCTCTGGGCGTTTCTTCATGCTGCTTTCTCAG`

**Dos hebras en juego (no confundir):** `strand_te` = en qué hebra matcheó Dfam (orientación
del TE); `strand_gene` = en qué hebra se transcribe el gen. La extracción **siempre usa
`strand_gene`** (queremos la secuencia como la lee la proteína). La relación entre ambas dice
si el TE quedó sentido (`==`) o antisentido (`≠`) respecto del gen — dato biológico, pero no
cambia cómo se extrae.

**Header:** `transcript|dfam_id|familia|chr:coords(hebra)|cds:a-b|aa:c-d` — trazabilidad completa.

**Chequeo de subcadena (el "seguro"):** dos caminos independientes que deben concordar — camino A
(coordenadas: liftOver→intersect→pmapToTranscripts dice "posición 28-63") y camino B (secuencia:
se extrae el ADN real). Se verifica que la secuencia B esté en la posición A del CDS. Da 100% →
toda la cadena de coordenadas (incl. hebra) es correcta.

### Decisiones (2026-07-14)
- **Chequeo de subcadena ahora es AUTOMÁTICO por corrida** (antes solo se validó una vez en el
  qmd). El script 02 verifica que cada región TE sea subcadena del CDS de su transcripto y
  reporta el % en `metrics.tsv` (`regiones_te_subcadena_cds_pct`); avisa si baja de 100. Es un
  chequeo de CONTENCIÓN (el pedacito ⊂ el CDS completo), no de igualdad; valida coords+hebra+
  liftOver por dos caminos de extracción independientes (coords del solape vs. coords de todos
  los exones). Chimp: 100%.
- **UniProt en el header de `cds_completos.fasta`**: pasa de `>tx` a `>tx|uniprot`; si no hay
  uniprot usa placeholder `>tx|NA` (header homogéneo). Varios uniprot van separados por `;`
  (Swiss-Prot + TrEMBL). El paso 05 corta el header en el `|`. En chimp 3174/3174 tienen uniprot.

## Estación 4 — Validación Dfam (paso 03, OPCIONAL)

QC independiente: ¿cada hit es de verdad un TE? Se re-busca con el motor de Dfam. Dos niveles:
Nivel 1 (manda el CDS completo del transcripto) y Nivel 2/rescate (ventana genómica
`[start_te-200, end_te+200]` = el TE completo + flancos; el solape solo es muy corto ~64 nt para
re-detectar). Resultado (baseline 2842 sanos): 2001 validado_cds + 827 validado_genomico + 2
familia_ok + 12 no_encontrado = **2828/2842 (99,5%)**.

Caveats: organismo `Homo sapiens` por defecto (`Pan troglodytes` da ERROR en Dfam; familias
compartidas). Escribe una **tabla SEPARADA** (`validacion_dfam_todos.tsv`), NO reescribe la madre.
Resumible (`dfam_cache/`). Lento (red).

### Decisiones (2026-07-14/15)
- Se agregó `REPORTE_validacion_dfam.txt` (desglose + caveats), como el del paso 04.
- **Bug encontrado y corregido:** el cambio del uniprot en el header de `cds_completos.fasta`
  rompió DOS lectores de CDS (el 07 y `validar_dfam.py`): buscaban por `tx` pero el header era
  `tx|uniprot`. Síntoma: `validado_cds=0` (todo caía al Nivel 2). Fix: ambos cortan el header en
  el `|`. Lección: al cambiar un formato de archivo, revisar TODOS sus lectores.

## Estación 5/6 — dominios proteicos (paso 05 rediseñado, 2026-07-16)

Tras estresar todos los casos de solape (ver artefacto `redundancia_casos.html`), el paso 05 quedó
en **dos pasos**:
- **STEP 1 `04a_hits_aa.fasta`**: TODOS los hits, uno por (posición genómica × isoforma), uniendo los
  trozos de exón. Crudo y redundante (respaldo).
- **STEP 2 `04b_dominios_aa.fasta`** (input de CD-Hit): dentro de cada **locus**, dedup por
  **CONTENCIÓN de secuencia** — si una AA es subcadena de otra, queda la más larga y el header lista
  todos los TEs (`rep=...  tes=SVA_B;AluSp`). Header completo: `rep | tes | tx | uniprot(;NA) | aa | len`.

**Por qué contención dentro del locus, sin mirar familia:** validado empíricamente — de 55 loci con
familias mixtas, los 55 son de contención (0 solapes parciales cross-familia). El paso 04 **ya** separó
los módulos genuinamente distintos en loci diferentes (split → L1-Alu-L1 en 2 loci; anidado_en genérico
→ 2 loci; modulo_distinto_adyacente → 2 loci). Así, dentro de un locus la contención siempre = "mismo
evento". Solo el SVA/Alu colapsa entre familias (es_subparte_sva); cualquier otro anidamiento se conserva
(está en loci distintos).

**Casos verificados end-to-end:** L00044 (Alu-en-SVA) → 1 dominio `tes=SVA_B;AluSp`; L01497 (4 Alus) → 1
dominio `tes=FLAM_A;AluJr;AluJr4;AluSx4`; L00417 (L1 fragmentado, 2 piezas no contenidas) → 2 dominios.

**Escalera:** 3597 hits crudos (STEP1) → **3101 dominios distintos** (STEP2); 82 loci multi-dominio;
2407 sanos sin stop interno. Reemplaza el enfoque viejo de "unión por locus" (que metía huecos no-TE).
