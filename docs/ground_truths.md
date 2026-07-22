# Ground Truths

## HEK293T GLORI

HEK293T predictions are benchmarked against
`metadata/ground_truth_datasets/HEK293T_GLORI_rep1_rep2_intersection_mean.tsv`.
The genomic table is the intersection of GSM6432590 and GSM6432591 from
GSE210563, with mean depth and mean stoichiometry. Validate it with:

```bash
python scripts/02_ground_truth/glori/prepare_glori_sites.py \
    --input metadata/ground_truth_datasets/HEK293T_GLORI_rep1_rep2_intersection_mean.tsv \
    --output /path/to/HEK293T_GLORI_prepared.tsv
```

## HeLa eTAM-seq

HeLa predictions are benchmarked against
`metadata/ground_truth_datasets/HeLa_eTAM_FTOm_FTOp_rep1_rep2_rep3_intersection_mean.tsv`.
This genomic table was derived locally from GSE211303 HeLa poly(A) WT FTO- and
FTO+ rep1-rep3 results. GSE211303 is the eTAM-seq SubSeries within GSE201064;
its BioProject is PRJNA869740. Validate it with:

```bash
python scripts/02_ground_truth/etam_seq/prepare_etam_sites.py \
    --input metadata/ground_truth_datasets/HeLa_eTAM_FTOm_FTOp_rep1_rep2_rep3_intersection_mean.tsv \
    --output /path/to/HeLa_eTAM_prepared.tsv
```

Both curated inputs already use genomic coordinates and do not require
transcript-to-genome conversion. Full method and accession citations are in
`docs/citations.md`.
