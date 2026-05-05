import re
import sys
from pathlib import Path
from typing import Optional

from synnet.utils.logger import get_logger, info, warning, error
from synnet.utils.io import open_gff

logger = get_logger(__name__)


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
        fin = open_gff(input_file)
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
                    warning(f"Skipping malformed line {n_lines}: {line[:80]}...")
                n_skipped += 1
                continue

            seqid, source, feature, start, end, score, strand, phase, attr_str = fields[:9]

            if feature != feat_type:
                continue

            attrs = parse_gff_attributes(attr_str)

            gene_id = extract_gene_id(attrs, key=id_key)

            if not gene_id:
                if verbose and n_skipped < 5:
                    warning(f"Skipping line {n_lines}: no {id_key} found in '{attr_str[:50]}...'")
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
                if verbose and n_skipped < 5:
                    warning(f"Skipping line {n_lines}: invalid coordinates '{start}-{end}'")
                n_skipped += 1
                continue

            bed_start = start_pos - 1
            bed_end = end_pos

            bed_score = score if score != '.' else '0'

            fout.write(f"{seqid}\t{bed_start}\t{bed_end}\t{gene_id}\t{bed_score}\t{strand}\n")
            n_written += 1

            if verbose and n_written % 1000 == 0:
                info(f"Processed {n_written} {feat_type} features...")

    finally:
        if fin is not sys.stdin:
            fin.close()
        if fout is not sys.stdout:
            fout.close()

    if verbose and n_skipped > 0:
        info(f"Skipped {n_skipped} malformed/filtered lines")

    info(f"Written {n_written} features to {output_name}")

    return output_file
