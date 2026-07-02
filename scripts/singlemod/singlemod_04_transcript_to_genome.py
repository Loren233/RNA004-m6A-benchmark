#!/usr/bin/env python3
import argparse
import bisect
import csv
import gzip
import sys
from collections import defaultdict


def open_text(path):
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "r")


def parse_gtf_attributes(attr_text):
    attrs = {}
    for item in attr_text.strip().split(";"):
        item = item.strip()
        if not item:
            continue
        if " " not in item:
            continue
        key, value = item.split(" ", 1)
        attrs[key] = value.strip().strip('"')
    return attrs


def chrom_sort_key(chrom):
    c = chrom
    if c.startswith("chr"):
        c = c[3:]
    order = {
        "M": 25,
        "MT": 25,
        "X": 23,
        "Y": 24,
    }
    if c.isdigit():
        return (0, int(c))
    return (1, order.get(c, 1000), c)


def build_transcript_maps(gtf_path):
    exons = defaultdict(list)

    with open_text(gtf_path) as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue

            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                continue

            chrom, source, feature, start, end, score, strand, frame, attrs_text = fields

            if feature != "exon":
                continue

            attrs = parse_gtf_attributes(attrs_text)
            transcript_id = attrs.get("transcript_id")

            if not transcript_id:
                continue

            exons[transcript_id].append(
                {
                    "chrom": chrom,
                    "start": int(start),
                    "end": int(end),
                    "strand": strand,
                }
            )

    maps = {}
    unversioned_maps = {}
    unversioned_duplicates = set()

    for transcript_id, exon_list in exons.items():
        if not exon_list:
            continue

        chrom = exon_list[0]["chrom"]
        strand = exon_list[0]["strand"]

        if strand == "+":
            exon_list = sorted(exon_list, key=lambda x: (x["start"], x["end"]))
        elif strand == "-":
            exon_list = sorted(exon_list, key=lambda x: (x["start"], x["end"]), reverse=True)
        else:
            continue

        intervals = []
        tx_cursor = 1

        for exon in exon_list:
            exon_len = exon["end"] - exon["start"] + 1
            tx_start = tx_cursor
            tx_end = tx_cursor + exon_len - 1

            intervals.append(
                {
                    "tx_start": tx_start,
                    "tx_end": tx_end,
                    "g_start": exon["start"],
                    "g_end": exon["end"],
                }
            )

            tx_cursor = tx_end + 1

        ends = [x["tx_end"] for x in intervals]

        mapping = {
            "transcript_id": transcript_id,
            "chrom": chrom,
            "strand": strand,
            "intervals": intervals,
            "ends": ends,
            "length": tx_cursor - 1,
        }

        maps[transcript_id] = mapping

        tid_no_version = transcript_id.split(".")[0]
        if tid_no_version in unversioned_maps:
            unversioned_duplicates.add(tid_no_version)
        else:
            unversioned_maps[tid_no_version] = mapping

    for tid in unversioned_duplicates:
        unversioned_maps.pop(tid, None)

    return maps, unversioned_maps


def tx_to_genome(mapping, tx_pos):
    i = bisect.bisect_left(mapping["ends"], tx_pos)
    if i >= len(mapping["intervals"]):
        return None

    interval = mapping["intervals"][i]

    if tx_pos < interval["tx_start"] or tx_pos > interval["tx_end"]:
        return None

    offset = tx_pos - interval["tx_start"]

    if mapping["strand"] == "+":
        g_pos = interval["g_start"] + offset
    else:
        g_pos = interval["g_end"] - offset

    return mapping["chrom"], g_pos, mapping["strand"]


def as_float(x, default=0.0):
    try:
        return float(x)
    except Exception:
        return default


def as_int_float(x, default=0.0):
    try:
        return float(x)
    except Exception:
        return default


def main():
    parser = argparse.ArgumentParser(
        description="Convert SingleMod transcript-level coordinates to genomic coordinates using a GTF annotation."
    )
    parser.add_argument("--singlemod", required=True, help="SingleMod transcript-level TSV")
    parser.add_argument("--gtf", required=True, help="GENCODE GTF or GTF.GZ")
    parser.add_argument("--out-raw", required=True, help="Output raw genomic TSV, one row per transcript-level site")
    parser.add_argument("--out-collapsed", required=True, help="Output collapsed genomic TSV, aggregated by chr, pos, strand")
    parser.add_argument("--unmapped", required=True, help="Output unmapped transcript-level rows")
    args = parser.parse_args()

    print(f"Loading GTF: {args.gtf}", file=sys.stderr)
    maps, unversioned_maps = build_transcript_maps(args.gtf)
    print(f"Transcript maps loaded: {len(maps)}", file=sys.stderr)

    raw_fields = [
        "chr",
        "pos_1based",
        "start_0based",
        "end",
        "site_id",
        "singlemod_stoich",
        "strand",
        "Total",
        "m6A_count",
        "motif",
        "transcript_id",
        "transcript_pos_1based",
        "transcript_site_id",
    ]

    collapsed_fields = [
        "chr",
        "pos_1based",
        "start_0based",
        "end",
        "site_id",
        "singlemod_stoich",
        "strand",
        "Total",
        "m6A_count",
        "motif",
        "transcript_count",
    ]

    required = {"ref", "pos_1based", "singlemod_stoich", "Total", "m6A_count", "motif"}

    collapsed = {}
    total_rows = 0
    mapped_rows = 0
    unmapped_rows = 0

    with open_text(args.singlemod) as sm_in, \
         open(args.out_raw, "w", newline="") as raw_out, \
         open(args.unmapped, "w", newline="") as unmapped_out:

        reader = csv.DictReader(sm_in, delimiter="\t")
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise SystemExit(f"Missing required SingleMod columns: {sorted(missing)}")

        raw_writer = csv.DictWriter(raw_out, delimiter="\t", fieldnames=raw_fields, lineterminator="\n")
        raw_writer.writeheader()

        unmapped_writer = csv.DictWriter(unmapped_out, delimiter="\t", fieldnames=reader.fieldnames, lineterminator="\n")
        unmapped_writer.writeheader()

        for row in reader:
            total_rows += 1

            transcript_id = row["ref"]
            try:
                tx_pos = int(row["pos_1based"])
            except Exception:
                unmapped_writer.writerow(row)
                unmapped_rows += 1
                continue

            mapping = maps.get(transcript_id)
            if mapping is None:
                mapping = unversioned_maps.get(transcript_id.split(".")[0])

            if mapping is None:
                unmapped_writer.writerow(row)
                unmapped_rows += 1
                continue

            converted = tx_to_genome(mapping, tx_pos)
            if converted is None:
                unmapped_writer.writerow(row)
                unmapped_rows += 1
                continue

            chrom, g_pos, strand = converted
            motif = row["motif"]
            total = as_int_float(row["Total"])
            m6a_count = as_int_float(row["m6A_count"])
            stoich = as_float(row["singlemod_stoich"])

            genomic_site_id = f"{chrom}|{g_pos}|{strand}|{motif}"

            out_row = {
                "chr": chrom,
                "pos_1based": g_pos,
                "start_0based": g_pos - 1,
                "end": g_pos,
                "site_id": genomic_site_id,
                "singlemod_stoich": stoich,
                "strand": strand,
                "Total": total,
                "m6A_count": m6a_count,
                "motif": motif,
                "transcript_id": transcript_id,
                "transcript_pos_1based": tx_pos,
                "transcript_site_id": row.get("site_id", ""),
            }

            raw_writer.writerow(out_row)
            mapped_rows += 1

            key = (chrom, g_pos, strand)
            if key not in collapsed:
                collapsed[key] = {
                    "Total": 0.0,
                    "m6A_count": 0.0,
                    "motifs": set(),
                    "transcript_count": 0,
                }

            collapsed[key]["Total"] += total
            collapsed[key]["m6A_count"] += m6a_count
            collapsed[key]["motifs"].add(motif)
            collapsed[key]["transcript_count"] += 1

    print(f"Rows total: {total_rows}", file=sys.stderr)
    print(f"Rows mapped: {mapped_rows}", file=sys.stderr)
    print(f"Rows unmapped: {unmapped_rows}", file=sys.stderr)
    print(f"Collapsed genomic sites: {len(collapsed)}", file=sys.stderr)

    with open(args.out_collapsed, "w", newline="") as collapsed_out:
        writer = csv.DictWriter(collapsed_out, delimiter="\t", fieldnames=collapsed_fields, lineterminator="\n")
        writer.writeheader()

        for chrom, g_pos, strand in sorted(collapsed.keys(), key=lambda x: (chrom_sort_key(x[0]), x[1], x[2])):
            item = collapsed[(chrom, g_pos, strand)]
            total = item["Total"]
            m6a_count = item["m6A_count"]
            stoich = m6a_count / total if total > 0 else 0.0
            motif = ",".join(sorted(item["motifs"]))
            site_id = f"{chrom}|{g_pos}|{strand}|{motif}"

            writer.writerow(
                {
                    "chr": chrom,
                    "pos_1based": g_pos,
                    "start_0based": g_pos - 1,
                    "end": g_pos,
                    "site_id": site_id,
                    "singlemod_stoich": stoich,
                    "strand": strand,
                    "Total": total,
                    "m6A_count": m6a_count,
                    "motif": motif,
                    "transcript_count": item["transcript_count"],
                }
            )


if __name__ == "__main__":
    main()
