# 03_standardization

This folder contains scripts for standardizing tool-specific RNA modification outputs into a common genomic site-level format.

The scripts should be used in this order.

## 01_standardize_site_tables.py

Purpose:

Convert raw site-level outputs from different tools into a common column structure.

Typical standardized columns include:

```text
chr / ref
pos_1based
start_0based
end
strand
modified_count
depth
stoich
motif
tool

Use this step when a tool output has tool-specific column names or formats.

Examples:

Dorado site-level output
DeepRM site-level output
m6Anet site-level output
SingleMod site-level output
02_convert_transcript_to_genome.py

Purpose:

Convert transcript-level coordinates to genomic coordinates using a GTF annotation.

Use this step when the standardized output is still based on transcript IDs, for example:

ENST00000000233.10    122

After conversion, the output should use genomic coordinates, for example:

chr7    127588532

This step is needed before comparing transcript-level tools with genomic ground truth datasets such as GLORI or eTAM.

This step is not needed if the tool output is already in genomic coordinates.

03_collapse_genomic_sites.py

Purpose:

Collapse duplicated genomic sites after coordinate conversion.

This is needed because multiple transcript isoforms can map to the same genomic site. Without collapsing, the same genomic position may be counted multiple times during benchmarking.

Collapse key:

chr + pos_1based + strand

Collapse rule:

depth: sum
modified_count: sum
stoich: modified_count / depth
motif: unique motifs joined by comma
source_rows: number of input rows collapsed into the genomic site

Use this step after 02_convert_transcript_to_genome.py if duplicated genomic sites are present.

Check duplicates with:

awk -F'\t' 'NR>1{
    k=$1"\t"$2"\t"$5
    c[k]++
}
END{
    dup=0
    for(k in c){
        if(c[k]>1) dup++
    }
    print "duplicated genomic sites:", dup
}' standardized_genomic.tsv

If the duplicate count is 0, this step is optional.

Recommended workflow

For transcript-level tools:

raw tool output
↓
01_standardize_site_tables.py
↓
02_convert_transcript_to_genome.py
↓
03_collapse_genomic_sites.py
↓
04_benchmarking

For tools already in genomic coordinates:

raw tool output
↓
01_standardize_site_tables.py
↓
03_collapse_genomic_sites.py, only if duplicated genomic sites exist
↓
04_benchmarking
