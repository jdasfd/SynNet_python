"""
test/gff2bed.py - Step 1: GFF3 to BED converter
output: BED6 format
chrom\tstart(0-based)\tend\tgene_id\tscore\tstrand

Usage:
    python gff2bed.py -i seqs/ -s species.lst
    python gff2bed.py -i seqs/ -s species.lst --output-dir beds
    python gff2bed.py -i seqs/Atha.gff -o Atha.bed
"""

import sys
import argparse
import re
from pathlib import Path
from typing import Optional, List


def log_info(msg):
    print(f"[INFO] {msg}", file=sys.stderr)


def log_warn(msg):
    print(f"[WARN] {msg}", file=sys.stderr)


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


def extract_gene_id(attrs: dict, key: str = "ID") -> Optional[str]:
    return attrs.get(key)


def gff3_to_bed(
        input_file: str,
        output_file: Optional[str] = None,
        *,
        feat_type: str = "mRNA",
        id_key: str = "ID",
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

    n_written = 0
    n_skipped = 0

    try:
        for line in fin:
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

            try:
                start_pos = int(start)
                end_pos = int(end)
                length = end_pos - start_pos + 1
                if length < min_length:
                    n_skipped += 1
                    continue
            except ValueError:
                n_skipped += 1
                continue

            bed_start = start_pos - 1
            bed_end = end_pos
            bed_score = score if score != '.' else '0'

            fout.write(f"{seqid}\t{bed_start}\t{bed_end}\t{gene_id}\t{bed_score}\t{strand}\n")
            n_written += 1

            if verbose and n_written % 10000 == 0:
                log_info(f"Written {n_written} features...")

    finally:
        if fin is not sys.stdin:
            fin.close()
        if fout is not sys.stdout:
            fout.close()

    if n_skipped > 0:
        log_info(f"Skipped {n_skipped} lines")
    log_info(f"Written {n_written} features to {output_name}")
    return output_file


def batch_gff2bed(
        input_dir: str,
        species_list_file: str,
        output_dir: Optional[str] = None,
        *,
        feat_type: str = "mRNA",
        id_key: str = "ID",
        min_length: int = 0,
        verbose: bool = False,
) -> List[str]:
    in_dir = Path(input_dir)

    with open(species_list_file, 'r') as f:
        species = [line.strip() for line in f if line.strip() and not line.startswith('#')]

    if len(species) == 0:
        log_error("Species list is empty")
        return []

    if output_dir:
        out_dir = Path(output_dir)
    else:
        out_dir = in_dir

    out_dir.mkdir(parents=True, exist_ok=True)

    results = []
    for sp in species:
        gff_file = in_dir / f"{sp}.gff"
        if not gff_file.exists():
            gff_file = in_dir / f"{sp}.gff3"
        if not gff_file.exists():
            log_warn(f"GFF not found for {sp} in {in_dir}, skipping")
            continue

        bed_file = out_dir / f"{sp}.bed"
        log_info(f"Converting: {gff_file.name} -> {bed_file.name}")

        try:
            result = gff3_to_bed(
                str(gff_file), str(bed_file),
                feat_type=feat_type, id_key=id_key,
                min_length=min_length, verbose=verbose,
            )
            results.append(result)
        except Exception as e:
            log_error(f"Failed to convert {gff_file.name}: {e}")

    return results


def main():
    parser = argparse.ArgumentParser(
        description="Step 1: Convert GFF3 to BED format",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Batch mode (input is directory, requires -s)
  python gff2bed.py -i seqs/ -s species.lst
  python gff2bed.py -i seqs/ -s species.lst --output-dir beds

  # Single file
  python gff2bed.py -i seqs/Atha.gff -o Atha.bed
        """,
    )

    parser.add_argument("-i", "--input", required=True,
                        help="Input GFF3 file or directory (batch mode with -s)")
    parser.add_argument("-o", "--output",
                        help="Output BED file (single file mode only)")
    parser.add_argument("-s", "--species-list",
                        help="Species list file for batch mode (one name per line)")
    parser.add_argument("--output-dir",
                        help="Output directory for BED files (batch mode, default: same as input)")
    parser.add_argument("-t", "--feat-type", default="mRNA",
                        help="Feature type to extract (default: mRNA)")
    parser.add_argument("-k", "--id-key", default="ID",
                        help="Attribute key for gene ID (default: ID)")
    parser.add_argument("--min-length", type=int, default=0,
                        help="Minimum feature length (default: 0)")
    parser.add_argument("-v", "--verbose", action="store_true",
                        help="Verbose output")

    args = parser.parse_args()

    log_info("Step 1: GFF3 to BED Converter")

    try:
        if args.species_list:
            results = batch_gff2bed(
                args.input, args.species_list, args.output_dir,
                feat_type=args.feat_type, id_key=args.id_key,
                min_length=args.min_length, verbose=args.verbose,
            )
            log_info(f"Done! Converted {len(results)} species")
        else:
            input_path = Path(args.input)
            if input_path.is_dir():
                log_error("Input is a directory but no -s species list provided")
                sys.exit(1)

            result = gff3_to_bed(
                args.input, args.output,
                feat_type=args.feat_type, id_key=args.id_key,
                min_length=args.min_length, verbose=args.verbose,
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
