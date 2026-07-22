# Ground-Truth Datasets

The files in this directory are project-curated genomic benchmark tables, not
raw sequencing inputs.

## HEK293T GLORI

`HEK293T_GLORI_rep1_rep2_intersection_mean.tsv` is the intersection of GLORI
replicates GSM6432590 and GSM6432591 from GSE210563 (BioProject PRJNA866163),
with mean depth and mean stoichiometry. Cite Liu et al., *Nature Biotechnology*
41, 355-366 (2023), DOI `10.1038/s41587-022-01487-9`.

## HeLa eTAM-seq

`HeLa_eTAM_FTOm_FTOp_rep1_rep2_rep3_intersection_mean.tsv` is a local table
derived from the HeLa poly(A) WT eTAM-seq results in GSE211303 (BioProject
PRJNA869740), which is part of SuperSeries GSE201064. It represents the
intersection of three FTO- and FTO+ replicate comparisons and retains mean
depth and mean stoichiometry. The source sample groups are GSM6463599-GSM6463601
(FTO-) and GSM6463602-GSM6463604 (FTO+). It is not an unchanged GEO download.
The committed table has 13,147 data rows and MD5
`a5705498d644dbdeb76b5cdeadd452ac`.

Cite Xiao et al., *Nature Biotechnology* 41, 993-1003 (2023), DOI
`10.1038/s41587-022-01587-6`.
