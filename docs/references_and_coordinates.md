# References and Coordinates

Coordinate handling is a core part of this benchmark. Do not hide coordinate
assumptions inside ad hoc column names or permanent offset patches.

## Reference Preparation

- GENCODE v49 transcript FASTA headers should be cleaned to retain transcript
  IDs only.
- Cleaned transcript reference names should match GTF `transcript_id` values.
- Tool outputs aligned to transcript FASTA may use ENST transcript coordinates.
- Tool outputs aligned to GRCh38 genome FASTA may already be genomic.

## Coordinate Systems

- BED start is 0-based.
- BED end is half-open.
- GTF genomic coordinates are 1-based inclusive.
- m6Anet `transcript_position` is treated as 0-based in the current conversion
  script design.

Use explicit column names:

- `start_0based`
- `end`
- `pos_1based`
- `transcript_position_0based`

Avoid ambiguous names such as `position` unless the coordinate system is also
encoded in the name or metadata.

## Strand Handling

Negative-strand transcript conversion must be verified against the GTF exon
structure. Do not use permanent plus-or-minus-one shifts as a substitute for
fixing coordinate definitions.

## Conversion Rules

Apply transcript-to-genome conversion only when the input table is transcript
coordinate based. Do not apply conversion to already genomic Dorado or SingleMod
tables.
