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
1-bis. **Cuarto motivo de exclusión: `orden_roto` (2026-07-21).** Los tres motivos de arriba miran
   *cuántos* exones sobrevivieron y *dónde*, pero no **en qué orden**. El liftOver puede correr un
   exón y romper la correspondencia entre `exon_number` y el orden genómico 5'→3'; ahí los pasos 01
   y 02 arman CDS distintos y las coordenadas aminoacídicas quedan corridas sin que ningún QC lo
   note. Es que el paso 02 (`extractTranscriptSeqs`) respeta el orden de la GRangesList
   (= `exon_number`), mientras que `pmapToTranscripts` lo **ignora** y ordena por coordenada
   genómica — verificado pasándole los mismos exones en dos órdenes distintos: devuelve lo mismo.
   Por eso no alcanza con ordenarlos antes de llamarla; hay que excluir el transcripto.

   Que es artefacto y no biología se comprueba mirando el GTF de chimp **antes** del liftOver:
   0 transcriptos desordenados de 49.949; después, 771 de 46.749. Las anotaciones que corren sin
   liftOver también dan 0 (hg38 v116: 0/278.455; v115: 0/211.446; pongo: 0/37.189).

   De esos 771, 741 ya caían en los otros tres motivos: `orden_roto` agrega **30**. Sólo 2 tenían
   hit filtrado, y de esos 2 uno tenía la coordenada mal (`ENSPTRT00000106637`, exón CDS 1 de 18 bp
   a 5,3 Mb de los otros 30 → su hit de SVA_A salía 18 nt = 6 aa corrido, y pasaba igual el QC de
   pertenencia del paso 02, que usa `grepl`). **Efecto en panTro6:** 3697→3695 hits, 3101→3099
   dominios; el resto de la cascada se mueve exactamente −2. hg38 y pongo no se mueven.

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

---

## Estación 7 — Mapeo transcripto → UniProt (2026-07-17/20)

La columna `uniprot` sale de un archivo de dos columnas (`ENST2UNIPROT` del config). **Es opcional
y no participa del filtrado ni del clustering**: si falta, todo corre igual y la columna queda `NA`.

**Ejemplo de la tabla de ENTRADA** — el xref que se baja del FTP de Ensembl (`tsv/<especie>/`):

```
gene_stable_id      transcript_stable_id  protein_stable_id   xref        db_name            info_type       source_identity  xref_identity  linkage_type
ENSG00000189337     ENST00000636564       ENSP00000489835     A0A1B0GTU0  Uniprot/SPTREMBL   DIRECT          100              100            -
ENSG00000204604     ENST00000595646       ENSP00000470381     Q5VIY5      Uniprot/SWISSPROT  DIRECT          -                -              -
ENSG00000204604     ENST00000595646       ENSP00000470381     Q5VIY5-1    Uniprot_isoform    DIRECT          -                -              -
```

**Ejemplo de la tabla de SALIDA** que consume el pipeline:

```
ensembl_transcript_id  uniprot_cons
ENST00000000233        P84085
ENST00000436979        P25440
ENST00000436979        P25440-3
```

Se genera con `helpers/ensembl_xref_to_map.sh entrada.uniprot.tsv.gz > map_tx_to_uniprot.tsv`
(escribe a stdout: el `>` no es opcional). Alternativa: `helpers/map_uniprot.qmd` (biomaRt), que
necesita red — **ojo, ese helper prioriza por FILA y no colapsa por transcripto**, así que puede
dejar varias filas del mismo tx.

### Decisión: se conservan TODAS las accessions de un transcripto

Un transcripto puede mapear a varias entradas de UniProt. En vez de elegir una al generar el map,
**se conservan todas** (una por fila) y decide después quien valida: `helpers/validar_vs_uniprot.py`
las prueba contra el proteoma, se queda con el mejor resultado y deja constancia en su columna
`uniprot_usado`. Elegir antes sería descartar información sin saber cuál sirve.

### Las dos columnas que hay que mirar: `db_name` e `info_type`

- `db_name`: `Uniprot_isoform` (accession de isoforma, p.ej. `P25440-3` → la proteína exacta de ese
  transcripto) · `Uniprot/SWISSPROT` (curada) · `Uniprot/SPTREMBL` (predicha). Hay otros valores
  que **no son accessions** (p.ej. `Uniprot_gn_trans_name`, que son nombres de gen): por eso el
  helper filtra por **lista blanca**, no por descarte.
- `info_type`: `DIRECT` (el vínculo viene declarado) vs `SEQUENCE_MATCH` (Ensembl lo asignó
  alineando secuencias con exonerate, por eso siempre da 100% de identidad).

**Dato práctico sobre `SEQUENCE_MATCH`:** esas accessions casi nunca están en el FASTA del proteoma
de referencia de UniProt, así que suelen terminar como `accession ausente del FASTA` al validar. No
es un problema biológico —tienen 100% de identidad, son la misma proteína con otro identificador,
de un proteoma redundante que UniProt no distribuye—, pero conviene saberlo al leer los resultados.

### Desempate al validar: primero la evidencia, después la curaduría

Cuando un transcripto trae varias accessions, `validar_vs_uniprot.py` elige así:

1. **manda el status**: si una predicha da `exacto` y una curada da `mismatch`, gana la predicha;
2. **sólo ante empate** gana SwissProt sobre trEMBL;
3. dentro de SwissProt, gana la accession de isoforma.

Sin el paso 2 el ganador sería el primero de la cadena, o sea el orden del archivo de Ensembl:
arbitrario. La columna `fuentes_disponibles` del TSV completa la lectura — dice qué había para
elegir (`sp`, `tr`, `sp,tr`), no sólo qué se usó, y así se distingue "validó con trEMBL porque era
lo único" de "validó con trEMBL habiendo una curada que no coincidió".

---

## Estación 8 — Campos multi-valor: la convención de separadores (2026-07-20)

Un dominio de `04b` puede estar sostenido por **varios transcriptos**, y cada uno tener **varias
accessions**. Antes `sano` se escribía con un slot por transcripto pero `uniprot` y `gene` eran
conjuntos aplanados y deduplicados que además **descartaban los `NA`**, así que la cantidad de
valores no coincidía con la de transcriptos y no había forma de saber cuál era de cuál.

**Convención, aplicada a `uniprot`, `gene` y `sano`:**

| separador | significa |
|---|---|
| `;` | separa **transcriptos**: siempre hay `n_tx` slots, en el mismo orden que la columna `tx` |
| `,` | separa varios valores **dentro** de un transcripto |
| `NA` | ocupa el slot vacío — **nunca se omite** |

```
tx      = ENSPTRT00000001312;ENSPTRT00000077350
uniprot = H2PYY3,A0A6D2VWW0;K7BA23,A0A6D2XP50      <- 2 accessions para cada uno de los 2 tx
```

Así `strsplit(campo, ";")` devuelve **siempre** tantos elementos como transcriptos, y la posición
*i* corresponde al *i*-ésimo. Es la misma convención de GFF3 (`;` entre atributos, `,` para varios
valores de un atributo).

**`tes` NO sigue esta convención**, a propósito: los TEs son otro eje. Un dominio puede tener varios
TEs *y* varios transcriptos a la vez, y meterlos en el mismo campo inventaría una correspondencia
que no existe.

### Salida nueva: `04b_productos_tx.tsv` (tabla larga)

Para eso está la tabla larga: **una fila por hit**, sin ningún campo multi-valor.

```
producto_id  id_locus  id_pos  te      tx                  gene                uniprot                sano  aa_start  aa_end
L00001_1     L00001    p00001  AluSc   ENSPTRT00000079696  ENSPTRG00000000028  A0A2I3SLA6             TRUE  119       151
L00002_1     L00002    p00002  7SLRNA  ENSPTRT00000087657  ENSPTRG00000000029  A0A2I3TBD6,A0A6D2XCG3  TRUE  111       155
```

Es el único lugar donde TE y transcripto pueden variar a la vez sin ambigüedad, y es la forma
cómoda de analizar en R. Los campos compactos quedan para el FASTA, donde hace falta una sola línea.

---

## Estación 9 — Validación contra UniProt (`helpers/validar_vs_uniprot.py`)

Compara la secuencia AA **inferida** por el pipeline contra la proteína **real** de UniProt. Corre
aparte del pipeline, después de que existan los `clusters_<GENOME>/`:

```bash
python3 helpers/validar_vs_uniprot.py UNIPROT_iso.fasta.gz clusters_X/04a_hits_aa.fasta \
        [--gene-xref data/Especie.Ensamblado.116.uniprot.tsv.gz]
```

Escribe `05_validacion_uniprot.tsv` (una fila por hit, con su status) y `.txt` (el reporte, que
cierra con un glosario de cada categoría). El reporte separa dos bloques:

- **[A] con mapeo directo tx→UniProt** — validación estricta: se comparan secuencia **y** posición.
  Categorías: `exacto`, `substr` (coordenadas corridas), `stop` (traducción rota), `otra_iso`
  (aparece en otra isoforma de la misma entrada), `mismatch`.
- **[B] sin mapeo directo** — el transcripto no tiene accession utilizable, así que no se puede
  validar posición. Se prueban las proteínas que UniProt tiene para **ese gen** (vía otros
  transcriptos, con `--gene-xref`): `otro_tx_del_gen`, `producto_distinto`, `no_acc`, `na`.

**Limitación estructural que conviene tener presente:** Ensembl sólo publica el xref transcripto→
UniProt cuando la traducción **ya coincide** con la proteína. Por eso el bloque [A] sólo existe
para transcriptos que Ensembl ya había verificado idénticos, y sus porcentajes miden si las
**coordenadas y la traducción del pipeline** son correctas — no son una estimación insesgada sobre
el total de hits. Los transcriptos cuyo producto difiere quedan sin accession y caen en [B].

---

## Pendientes de diseño (abiertos, 2026-07-21)

1. **El liftOver de HITS no tiene validación.** La rama `LIFTOVER_TARGET="gtf"` reporta 6 métricas,
   un mensaje y una tabla de transcriptos excluidos (`roto_chrom_hebra` / `multimapeo` /
   `perdida_parcial` / `orden_roto`). La rama `"hits"` es una sola línea sin nada: `unlist(liftOver(...))`, que
   esconde tres fallos silenciosos — hits no mapeados (desaparecen sin contarse), multimapeo (1 hit
   → N rangos, infla conteos) y fragmentación (un hit que cruza un hueco de la chain se parte, y
   cada fragmento tiene menos bp de solape, pudiendo caer bajo `OVERLAP_BP_THR`). Para arreglarlo
   hay que darle un **ID estable a cada hit antes** del liftOver y decidir qué hacer con los
   fragmentados (¿descartar, unir por min-max, quedarse con el más largo?).
2. **La validación por transcripto corre siempre**, incluso con `LIFTOVER_TARGET="none"`, donde
   `n_out == n_in` por construcción. El log muestra `ok:N roto:0 multimapeo:0` y parece un chequeo
   exitoso, pero es vacío. Convendría silenciarla o marcarla como no aplicable.
3. **Versión general.** Hoy el formato `nrph` de Dfam está asumido en tres puntos (nombres de
   columna, filtro por e-value, cruce con `DFAM_TAXA` por accession). Aceptar un BED/GFF genérico
   de anotaciones propias + variables que indiquen qué columna es cada cosa permitiría cruzar
   cualquier anotación con cualquier GTF.
