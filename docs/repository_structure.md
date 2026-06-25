# Repository Structure

This repository separates configuration, metadata, scripts, scheduler templates,
documentation, and generated outputs.

## configs/

Place shared configuration examples here. Committed examples must use placeholder
paths only. Local path files such as `configs/paths.sh` are ignored by Git.

## metadata/

Place verified metadata tables here. Current templates cover nanopore datasets,
ground-truth datasets, and software versions. Leave fields blank until the
values are known and verified.

## environments/

Place reproducible environment definitions here after tool versions and install
methods are confirmed.

## scripts/

Place future command-line scripts here. Scripts should correspond to documented
pipeline stages and should avoid hidden assumptions about local paths.

Expected future stages:

1. Dorado basecalling.
2. Minimap2 splice-aware alignment.
3. BAM quality control.
4. DeepRM.
5. m6Anet.
6. Dorado modification calling.
7. SingleMod.
8. GLORI preparation for HEK293T.
9. eTAM-seq preparation for HeLa.
10. Site-table standardization.
11. Coordinate-offset matching.
12. Depth-stratified correlation analysis.
13. Visualization.

## slurm/

Place future SLURM job templates here. Templates should call scripts or clearly
document the planned command. Use environment variables or placeholders for
cluster-specific values.

## results/

Place generated outputs here during local or HPC runs. Large result files are
ignored by Git. Small final summary tables or example outputs may be committed
later only if they are intentionally curated and documented.
