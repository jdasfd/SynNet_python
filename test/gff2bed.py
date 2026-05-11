"""
test/gff2bed.py - Step 1: GFF3 to BED converter
output: BED6
chrom\tstart(0-based)\tend\tgene_id\tscore\tstrand

Usage:
    python gff2bed.py -i genome.gff -o genome.bed
    python gff2bed.py -i genome.gff -t gene -k Name --primary-only
    python gff2bed.py -i genome.gff --min-length 100
"""

import sys
import argparse
import re
from pathlib import Path
from typing import Optional


def log_info(msg):
    print(f"[INFO] {msg}", file=sys.stderr)


def log_error(msg):
    print(f"[ERROR] {msg}", file=sys.stderr)


def parse_gff_attributes(attr_str: str) -> dict:
    attrs = {}
    if not attr_str or attr_str == '.':
        return attrs
    for item in attr_str.strip().split(';'):
        item = item.strip()
        if not item:
            continue
        match = re.match(r'([^=:\s]+)[=:](.+)', item)
        if match:
            key, val = match.groups()
            attrs[key.strip()] = val.strip()
    return attrs


def extract_gene_id(attrs: dict, key: str = "ID", parent_key: str = "Parent") -> Optional[str]:
    if key in attrs:
        return attrs[key]
    if parent_key in attrs:
        parents = [p.strip() for p in attrs[parent_key].split(',') if p.strip()]
        if parents:
            return parents[0]
    return attrs.get("ID")


def gff3_to_bed(
        input_file: str,
        output_file: Optional[str] = None,
        *,
        feat_type: str = "mRNA",
        id_key: str = "ID",
        primary_only: bool = False,
        min_length: int = 0,
        verbose: bool = False,
) -> str:
    if output_file is None:
        if input_file == "-":
            output_file = "stdout.bed"
        else:
            output_file = str(Path(input_file).with_suffix(".bed"))

    if input_file == "-":
        fin = sys.stdin
        input_name = "stdin"
    else:
        fin = open(input_file, 'r')
        input_name = input_file

    if output_file == "-":
        fout = sys.stdout
        output_name = "stdout"
    else:
        Path(output_file).parent.mkdir(parents=True, exist_ok=True)
        fout = open(output_file, 'w')
        output_name = output_file

    seen_ids = set() if primary_only else None
    n_lines = 0
    n_written = 0
    n_skipped = 0

    try:
        for line in fin:
            n_lines += 1
            if line.startswith('#') or not line.strip():
                continue

            fields = line.rstrip('\n').split('\t')
            if len(fields) < 9:
                n_skipped += 1
                continue

            seqid, source, feature, start, end, score, strand, phase, attr_str = fields[:9]

            if feature != feat_type:
                continue

            attrs = parse_gff_attributes(attr_str)
            gene_id = extract_gene_id(attrs, key=id_key)

            if not gene_id:
                n_skipped += 1
                continue

            if primary_only:
                if gene_id in seen_ids:
                    continue
                seen_ids.add(gene_id)

            try:
                start_pos = int(start)
                end_pos = int(end)
                length = end_pos - start_pos + 1
                if length < min_length:
                    continue
            except ValueError:
                n_skipped += 1
                continue

            bed_start = start_pos - 1
            bed_end = end_pos
            bed_score = score if score != '.' else '0'

            fout.write(f"{seqid}\t{bed_start}\t{bed_end}\t{gene_id}\t{bed_score}\t{strand}\n")
            n_written += 1

            if verbose and n_written % 1000 == 0:
                log_info(f"Processed {n_written} {feat_type} features...")

    finally:
        if fin is not sys.stdin:
            fin.close()
        if fout is not sys.stdout:
            fout.close()

    if verbose and n_skipped > 0:
        log_info(f"Skipped {n_skipped} malformed/filtered lines")

    log_info(f"Written {n_written} features to {output_name}")
    return output_file


def main():
    parser = argparse.ArgumentParser(
        description="Step 1: Convert GFF3 to BED format",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python gff2bed.py -i genome.gff3 -o genome.bed
  python gff2bed.py -i genome.gff3 -t gene -k Name --primary-only
  python gff2bed.py -i genome.gff3 --min-length 100
        """,
    )

    parser.add_argument("-i", "--input", required=True,
                        help="Input GFF3 file")
    parser.add_argument("-o", "--output",
                        help="Output BED file (default: {input}.bed)")
    parser.add_argument("-t", "--feat-type", default="mRNA",
                        help="Feature type to extract (default: mRNA)")
    parser.add_argument("-k", "--id-key", default="ID",
                        help="Attribute key for gene ID (default: ID)")
    parser.add_argument("--primary-only", action="store_true",
                        help="Keep only one entry per gene ID")
    parser.add_argument("--min-length", type=int, default=0,
                        help="Minimum feature length (default: 0)")
    parser.add_argument("-v", "--verbose", action="store_true",
                        help="Verbose output")

    args = parser.parse_args()

    log_info("Step 1: GFF3 to BED Converter")

    try:
        result = gff3_to_bed(
            args.input,
            args.output,
            feat_type=args.feat_type,
            id_key=args.id_key,
            primary_only=args.primary_only,
            min_length=args.min_length,
            verbose=args.verbose,
        )
        log_info(f"Done! Output: {result}")
    except FileNotFoundError as e:
        log_error(f"File not found: {e}")
        sys.exit(1)
    except Exception as e:
        log_error(f"Failed: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()
