"""
test/gff2bed.py - GFF3 to BED converter
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
from typing import Optional, TextIO

test_dir = Path(__file__).parent
if str(test_dir) not in sys.path:
    sys.path.insert(0, str(test_dir))

from utils.logger import setup_logger, info, warning, error, success, debug


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
                if verbose and n_skipped < 5:
                    print(f"Skipping malformed line {n_lines}: {line[:80]}...", file=sys.stderr)
                n_skipped += 1
                continue

            seqid, source, feature, start, end, score, strand, phase, attr_str = fields[:9]

            if feature != feat_type:
                continue

            attrs = parse_gff_attributes(attr_str)

            gene_id = extract_gene_id(attrs, key=id_key)

            if not gene_id:
                if verbose and n_skipped < 5:
                    print(f"Skipping line {n_lines}: no {id_key} found in '{attr_str[:50]}...'", file=sys.stderr)
                n_skipped += 1
                continue

            # primary_only
            if primary_only:
                if gene_id in seen_ids:
                    continue
                seen_ids.add(gene_id)

            # length_filtering
            try:
                start_pos = int(start)
                end_pos = int(end)
                length = end_pos - start_pos + 1
                if length < min_length:
                    continue
            except ValueError:
                if verbose and n_skipped < 5:
                    print(f"Skipping line {n_lines}: invalid coordinates '{start}-{end}'", file=sys.stderr)
                n_skipped += 1
                continue

            # convert 1-based to 0-based
            bed_start = start_pos - 1
            bed_end = end_pos

            bed_score = score if score != '.' else '0'

            # write into bed 6
            fout.write(f"{seqid}\t{bed_start}\t{bed_end}\t{gene_id}\t{bed_score}\t{strand}\n")
            n_written += 1

            if verbose and n_written % 1000 == 0:
                print(f"✓ Processed {n_written} {feat_type} features...", file=sys.stderr)

    finally:
        if fin is not sys.stdin:
            fin.close()
        if fout is not sys.stdout:
            fout.close()

    if verbose:
        if n_skipped > 0:
            print(f"Skipped {n_skipped} malformed/filtered lines", file=sys.stderr)

    return output_file

def main():
    parser = argparse.ArgumentParser(
        description="GFF3 to BED converter",
        epilog="Examples:\n"
               "  python gff2bed.py -i genome.gff3 -o genome.bed\n"
               "  python gff2bed.py -i genome.gff3 -t gene -k Name --primary-only\n"
               "  cat genome.gff3 | python gff2bed.py -o - > genome.bed",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    parser.add_argument("-i", "--input", default="-",
                        help="Input GFF3 file (default: stdin, use '-' for explicit stdin)")
    parser.add_argument("-o", "--output",
                        help="Output BED file (default: {input}.bed, use '-' for stdout)")
    parser.add_argument("-t", "--feat-type", default="mRNA",
                        help="Feature type to extract (default: mRNA, e.g., gene/exon/CDS)")
    parser.add_argument("-k", "--id-key", default="ID",
                        help="Attribute key for gene ID (default: ID, e.g., Name/gene_id)")
    parser.add_argument("--primary-only", action="store_true",
                        help="Keep only one entry per gene ID (remove duplicates)")
    parser.add_argument("--min-length", type=int, default=0,
                        help="Minimum feature length to keep (default: 0, no filter)")
    parser.add_argument("-v", "--verbose", action="store_true",
                        help="Print progress information to stderr")

    args = parser.parse_args()

    try:
        result = gff3_to_bed(
            args.input,
            args.output,
            feat_type=args.feat_type,
            id_key=args.id_key,
            primary_only=args.primary_only,
            min_length=args.min_length,
            verbose=args.verbose or args.input == "-",
        )
        if args.verbose:
            print(f"Done: {result}", file=sys.stderr)
        return 0
    except FileNotFoundError as e:
        print(f"Error: {e}", file=sys.stderr)
        return 1
    except Exception as e:
        print(f"Unexpected error: {e}", file=sys.stderr)
        if args.verbose:
            import traceback
            traceback.print_exc()
        return 2


if __name__ == "__main__":
    sys.exit(main())
