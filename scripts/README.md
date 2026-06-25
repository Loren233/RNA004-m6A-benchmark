# scripts

Future pipeline and analysis scripts should be placed here.

Planned script groups:

- Basecalling wrappers for Dorado.
- Splice-aware alignment wrappers for Minimap2.
- BAM quality-control summaries.
- Tool-specific preprocessing and output collection for DeepRM, m6Anet, Dorado
  modification calling, and SingleMod.
- GLORI preparation for HEK293T.
- eTAM-seq preparation for HeLa.
- Site-table standardization.
- Coordinate-offset matching.
- Depth-stratified correlation analysis.
- Visualization.

Scripts should be thin, documented command-line entry points. Shared logic
should be factored only after a real implementation exists.

Current status: planned work only. No runnable analysis scripts are included.
