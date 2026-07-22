# results

This directory contains selected final benchmark figures that are documented for
public viewing. Large raw data, large intermediate tables, matched-site tables,
correlation intermediates, logs, and temporary outputs are not stored here for
Git publication.

## Main Figure

| File | Dataset | Cell line | Tools | Ground truth | Purpose |
|---|---|---|---|---|---|
| `PRJEB80229_HeLa_total_benchmark.png` | PRJEB80229 RNA004 | HeLa | Dorado, m6Anet, SingleMod, DeepRM | eTAM-seq | Main README figure showing a four-tool HeLa benchmark across read-depth strata. |
| `PRJEB80229_HeLa_total_benchmark.pdf` | PRJEB80229 RNA004 | HeLa | Dorado, m6Anet, SingleMod, DeepRM | eTAM-seq | Full-resolution PDF for the main figure. |

## Supplementary Public Figures

| File | Dataset | Cell line | Tools | Ground truth | Purpose |
|---|---|---|---|---|---|
| `PRJEB80229_HEK293T_total_benchmark.png` | PRJEB80229 RNA004 | HEK293T | Dorado, m6Anet, SingleMod, DeepRM | GLORI | README supplementary figure showing the same benchmark structure for a GLORI comparison. |
| `PRJEB80229_HEK293T_total_benchmark.pdf` | PRJEB80229 RNA004 | HEK293T | Dorado, m6Anet, SingleMod, DeepRM | GLORI | Full-resolution PDF for the PRJEB80229 HEK293T figure. |
| `PRJEB106788_HEK293T_total_benchmark.png` | PRJEB106788 RNA004 | HEK293T | Dorado, m6Anet, SingleMod, DeepRM | GLORI | Additional public HEK293T benchmark figure listed here but not embedded on the README front page. |
| `PRJEB106788_HEK293T_total_benchmark.pdf` | PRJEB106788 RNA004 | HEK293T | Dorado, m6Anet, SingleMod, DeepRM | GLORI | Full-resolution PDF for the PRJEB106788 HEK293T figure. |

## Git Policy

`.gitignore` ignores `results/**` by default and allowlists only this README and
the selected final benchmark PNG/PDF files listed above. Do not add temporary
figures, offset checks, intermediate TSV/CSV files, raw data, logs, or scratch
outputs to Git without a separate review.
