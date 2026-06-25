# configs

Configuration files for the benchmark will live here.

Use `paths.example.sh` as a template for local path configuration. Copy it to an
untracked file such as `paths.sh` before filling in machine-specific paths.

Planned future files may include dataset-specific YAML or TSV configs for:

- HEK293T nanopore versus GLORI.
- HeLa nanopore versus eTAM-seq.
- Tool parameters for Dorado, Minimap2, DeepRM, m6Anet, and SingleMod.
- Depth bins and coordinate-offset matching rules.
