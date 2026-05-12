import gzip
import sys
from pathlib import Path
from typing import Union, TextIO, Optional, List, Dict, Set, Tuple, Iterator

from synnet.utils.logger import info, warning, error, debug


def open_file(filepath: Union[str, Path], mode: str = 'rt') -> TextIO:
    if str(filepath) == "-":
        return sys.stdin if 'r' in mode else sys.stdout

    path = Path(filepath)

    if path.suffix == ".gz":
        return gzip.open(path, mode=mode)

    return open(path, mode=mode)


def open_gff(filepath: str, mode: str = 'rt') -> TextIO:
    return open_file(filepath, mode)


def check_files(*files: Union[str, Path], required: bool = True) -> bool:
    missing = [str(f) for f in files if not Path(f).exists()]

    if missing:
        msg = f"Missing files: {', '.join(missing)}"
        if required:
            error(msg)
        return False
    return True


def ensure_dir(dirpath: Union[str, Path]) -> Path:
    path = Path(dirpath)
    path.mkdir(parents=True, exist_ok=True)
    return path


def backup_file(filepath: str, suffix: str = ".bak") -> str:
    src = Path(filepath)
    if src.exists():
        dst = src.with_suffix(src.suffix + suffix)
        src.rename(dst)
        return str(dst)
    return str(src)


def read_lines(filepath: Union[str, Path], skip_comments: bool = True, skip_empty: bool = True) -> List[str]:
    with open_file(filepath, 'r') as f:
        lines = []
        for line in f:
            line = line.rstrip('\n')
            if skip_empty and not line:
                continue
            if skip_comments and line.startswith('#'):
                continue
            lines.append(line)
        return lines


def read_tsv(filepath: Union[str, Path], has_header: bool = True, skip_comments: bool = True) -> Tuple[List[str], List[List[str]]]:
    header = []
    rows = []
    with open_file(filepath, 'r') as f:
        for i, line in enumerate(f):
            line = line.rstrip('\n')
            if skip_comments and line.startswith('#'):
                continue
            if not line:
                continue
            parts = line.split('\t')
            if has_header and i == 0 and not header:
                header = parts
            else:
                rows.append(parts)
    return header, rows


def write_tsv(filepath: Union[str, Path], rows: List[List[str]], header: Optional[List[str]] = None) -> None:
    ensure_dir(Path(filepath).parent)
    with open_file(filepath, 'w') as f:
        if header:
            f.write('\t'.join(header) + '\n')
        for row in rows:
            f.write('\t'.join(str(x) for x in row) + '\n')
    info(f"Written: {filepath}")


def read_species_list(list_file: str) -> List[str]:
    species = read_lines(list_file, skip_comments=True, skip_empty=True)
    if len(species) < 2:
        raise ValueError(f"Need >= 2 species, got {len(species)}")
    return species


def read_anchors_file(
        filepath: Union[str, Path],
        min_score: float = 0,
        exclude_lifted: bool = False,
) -> Iterator[Tuple[str, str, float, bool, int]]:
    block_id = 0
    with open_file(filepath, 'r') as f:
        for line in f:
            if line.startswith('###'):
                block_id += 1
                continue

            if line.startswith('#') or not line.strip():
                continue

            parts = line.rstrip('\n').split('\t')
            if len(parts) < 3:
                continue

            gene_a, gene_b = parts[0], parts[1]
            score_str = parts[2]

            is_lifted = score_str.endswith('L')
            try:
                weight = float(score_str.rstrip('L'))
            except ValueError:
                weight = 0.0

            if exclude_lifted and is_lifted:
                continue

            if weight < min_score:
                continue

            yield gene_a, gene_b, weight, is_lifted, block_id


def write_network_tsv(
        filepath: Union[str, Path],
        edges: List[Tuple[str, str, float, bool, str, int]],
) -> None:
    ensure_dir(Path(filepath).parent)
    with open_file(filepath, 'w') as f:
        f.write("source\ttarget\tscore\tis_lifted\tspecies_pair\tblock_id\n")
        for src, tgt, score, is_lifted, sp_pair, block_id in edges:
            f.write(f"{src}\t{tgt}\t{score}\t{is_lifted}\t{sp_pair}\t{block_id}\n")
    info(f"Exported: {filepath}")


def write_stats_file(
        filepath: Union[str, Path],
        stats: Dict[str, any],
        title: str = "Statistics",
) -> None:
    ensure_dir(Path(filepath).parent)
    with open_file(filepath, 'w') as f:
        f.write(f"# {title}\n\n")
        for key, value in stats.items():
            if isinstance(value, dict):
                f.write(f"\n# {key}\n")
                for k, v in value.items():
                    f.write(f"  {k}: {v}\n")
            else:
                f.write(f"{key}: {value}\n")
    info(f"Exported: {filepath}")


def read_synnet_tsv(
        filepath: Union[str, Path],
) -> Tuple[List[Tuple[str, int, str, str]], Set[str], Dict[str, Set[str]]]:
    edges = []
    nodes = set()
    cluster_genes = {}

    with open_file(filepath, 'r') as f:
        header = f.readline()
        for line in f:
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 4:
                continue
            cluster_id, score_str, node1, node2 = parts[0], parts[1], parts[2], parts[3]
            try:
                score = int(score_str)
            except ValueError:
                score = 0
            edges.append((cluster_id, score, node1, node2))
            nodes.add(node1)
            nodes.add(node2)
            if cluster_id not in cluster_genes:
                cluster_genes[cluster_id] = set()
            cluster_genes[cluster_id].add(node1)
            cluster_genes[cluster_id].add(node2)

    return edges, nodes, cluster_genes


def read_bed_file(filepath: Union[str, Path]) -> Dict[str, Dict[str, any]]:
    genes = {}
    with open_file(filepath, 'r') as f:
        for line in f:
            parts = line.rstrip('\n').split('\t')
            if len(parts) >= 4:
                genes[parts[3]] = {
                    'chrom': parts[0],
                    'start': int(parts[1]),
                    'end': int(parts[2]),
                    'strand': parts[5] if len(parts) > 5 else '+',
                }
    return genes


def build_gene_species_map(species_list: List[str], bed_dir: Union[str, Path]) -> Dict[str, str]:
    gene_map = {}
    bed_path = Path(bed_dir)
    for sp in species_list:
        bed_file = bed_path / f"{sp}.bed"
        if not bed_file.exists():
            warning(f"BED file not found: {bed_file}")
            continue
        with open_file(bed_file, 'r') as f:
            for line in f:
                parts = line.rstrip('\n').split('\t')
                if len(parts) >= 4:
                    gene_map[parts[3]] = sp
    return gene_map
