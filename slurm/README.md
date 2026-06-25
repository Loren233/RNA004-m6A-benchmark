# slurm

Future SLURM job templates should be placed here.

Recommended naming pattern:

- `01_dorado_basecalling.sbatch`
- `02_minimap2_alignment.sbatch`
- `03_bam_qc.sbatch`
- `04_deeprm.sbatch`
- `05_m6anet.sbatch`
- `06_dorado_mod_calling.sbatch`
- `07_singlemod.sbatch`
- `08_prepare_glori_hek293t.sbatch`
- `09_prepare_etamseq_hela.sbatch`
- `10_standardize_sites.sbatch`
- `11_match_coordinates.sbatch`
- `12_depth_correlation.sbatch`
- `13_visualization.sbatch`

Committed SLURM templates must use placeholders or environment variables for
project paths, data paths, partitions, accounts, and scratch directories. Do not
commit personal HPC paths, account names, credentials, or cluster-specific
secrets.

Current status: planned work only. No SLURM job scripts are included.

Resource guidance:

- Dorado basecalling is a GPU job.
- DeepRM inference is a GPU job.
- SingleMod neural-network prediction should request a GPU only if the current
  implementation uses CUDA.
- Samtools, Modkit, f5c eventalign, coordinate conversion, table formatting, and
  comparison steps should be CPU jobs unless the program explicitly requires a
  GPU.
- Do not request GPUs for non-GPU tasks.
- Do not substantially increase CPU, memory, or wall-time values without
  evidence from existing scripts or completed runs.
