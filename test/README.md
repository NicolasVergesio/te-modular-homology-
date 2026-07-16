# Test mínimo

Datos **sintéticos** (chr1 de 2000 bp, 1 gen `protein_coding`, 1 TE que solapa el CDS). Valida
que la instalación y el flujo `01→05` funcionen — corre en **segundos**, sin datos grandes ni red.

## Correr

```bash
bash run.sh test/config.test.sh
```

## Salida esperada

```
01: 1 hit filtrado
02: 1 transcripto sano ; región TE subcadena del CDS = 100%
04: 1 locus
05: STEP1 1 hit → STEP2 1 dominio
```

El dominio (`test/resultados_test/clusters_test/04b_dominios_aa.fasta`):
```
>L00001_1 | rep=TestAlu | tes=TestAlu | tx=TESTT01 | uniprot=NA | aa=34-84 | len=51
AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA   (51 aa)
```

Si ves esos números, la instalación está OK. Los resultados van a `test/resultados_test/`
(ignorado por git).

## Qué contiene `data/`

| Archivo | Qué es |
|---|---|
| `test.fa` | genoma `chr1` (2 kb) con un ORF válido (ATG + 98 codones Ala + stop) |
| `test.gtf` | 1 gen / 1 transcripto `protein_coding`, 1 exón CDS (501-797) |
| `test.nrph.hits` | 1 hit Dfam (`TestAlu`) que solapa el CDS (600-750) |
| `test_dfam_taxa.tsv` | taxonomía del TE del test |

Ejercita: filtro CDS protein_coding, normalización de cromosomas (`1`→`chr1`), intersect,
proyección a aa, QC `sano`, chequeo de subcadena, clustering y extracción de dominios — todo
sin liftOver (mismo ensamblado) ni datos externos.
