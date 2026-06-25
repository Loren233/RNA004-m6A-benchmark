# Repository Rules for Agents

This repository is a scaffold for a reproducible RNA004 direct RNA m6A
benchmark. Keep changes conservative, transparent, and reproducible.

## Do

- Keep raw sequencing data, large intermediates, and generated benchmark results
  out of version control.
- Use placeholder paths in committed examples.
- Put local machine, HPC, and data paths in untracked files derived from
  `configs/paths.example.sh`.
- Clearly label incomplete methods, scripts, and outputs as planned work.
- Record real software versions only after they have been verified.
- Preserve the separation between metadata, configuration, scripts, SLURM job
  templates, documentation, and generated results.

## Do Not

- Do not hard-code personal HPC paths or user-specific directories.
- Do not add POD5, FAST5, FASTQ, BAM, BAI, SAM, BLOW5, eventalign outputs, or
  other large data files.
- Do not invent accession numbers, dataset identifiers, software versions,
  completed metrics, or benchmark results.
- Do not add runnable analysis scripts unless the commands, inputs, outputs,
  dependencies, and assumptions are documented.
- Do not initialize Git, commit, or push unless explicitly requested.

## Planned Pipeline Scope

The planned workflow covers Dorado basecalling, Minimap2 splice-aware alignment,
BAM quality control, DeepRM, m6Anet, Dorado modification calling, SingleMod,
GLORI preparation for HEK293T, eTAM-seq preparation for HeLa, site-table
standardization, coordinate-offset matching, depth-stratified correlation
analysis, and visualization.
