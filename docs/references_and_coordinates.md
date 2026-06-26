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

## Transcript-to-Genome Converter

Dorado, m6Anet, and possibly DeepRM can initially produce transcript-coordinate
site tables when reads are aligned to the cleaned GENCODE transcript FASTA. The
generic converter at `scripts/09_standardization/convert_transcript_to_genome.py`
uses GENCODE v49 exon annotations from `gencode.v49.annotation.gtf` to convert
transcript positions into genomic coordinates before site-table standardization
and GLORI comparison.

Coordinate conventions used by the converter:

- Transcript input positions are normalized to 0-based coordinates.
- Dorado transcript tables use `start_0based`; do not use Dorado `pos_1based`
  as a 0-based transcript coordinate.
- m6Anet `transcript_position` is treated as 0-based.
- GTF exon coordinates are 1-based inclusive.
- Genomic output `pos_1based` is 1-based.
- Genomic output `start_0based` is BED-style 0-based.
- Genomic output `end_0based_exclusive` is BED-style half-open.

The converter does not aggregate duplicate genomic positions. Distinct
transcript isoforms may map to the same `chr`, `pos_1based`, and `strand`;
duplicate-site aggregation remains a separate standardization decision.
