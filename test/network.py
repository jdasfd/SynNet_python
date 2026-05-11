"""
test/network.py - Step 3: Build synteny network from .anchors / .lifted.anchors files

Usage:
    python network.py -s species.lst -d .
    python network.py -s species.lst --use-anchors
    python network.py -s species.lst --min-score 100
    python network.py -s species.lst -o Final_Network --formats tsv,graphml
"""

import sys
import argparse
from pathlib import Path
from typing import List, Dict, Set, Optional, Tuple
from dataclasses import dataclass, field
from collections import defaultdict

try:
    import networkx as nx
    _HAS_NX = True
except ImportError:
    _HAS_NX = False


def log_info(msg):
    print(f"[INFO] {msg}", file=sys.stderr)


def log_warn(msg):
    print(f"[WARN] {msg}", file=sys.stderr)


def log_error(msg):
    print(f"[ERROR] {msg}", file=sys.stderr)


@dataclass
class AnchorEdge:
    source: str
    target: str
    score: float
    is_lifted: bool
    species_pair: str
    block_id: int = 0


@dataclass
class NetworkStats:
    total_nodes: int = 0
    total_edges: int = 0
    lifted_edges: int = 0
    filtered_by_score: int = 0
    species_pair_counts: Dict[str, int] = field(default_factory=dict)


def load_species_list(list_file: str) -> List[str]:
    with open(list_file, 'r') as f:
        species = [line.strip() for line in f if line.strip() and not line.startswith('#')]
    if len(species) < 2:
        raise ValueError(f"Need >= 2 species, got {len(species)}")
    return species


def parse_anchor_line(line: str) -> Optional[Tuple[str, str, float, bool]]:
    if line.startswith('#') or not line.strip():
        return None
    parts = line.rstrip('\n').split('\t')
    if len(parts) < 3:
        return None
    gene_a, gene_b = parts[0], parts[1]
    score_str = parts[2]
    is_lifted = score_str.endswith('L')
    try:
        weight = float(score_str.rstrip('L'))
    except ValueError:
        weight = 0.0
    return gene_a, gene_b, weight, is_lifted


def parse_anchors_file(
        filepath: Path, species_pair: str,
        *, min_score: float = 0, exclude_lifted: bool = False,
) -> Tuple[List[AnchorEdge], int, int]:
    edges = []
    block_id = 0
    n_total = 0
    n_skipped = 0

    with open(filepath, 'r') as f:
        for line in f:
            if line.startswith('###'):
                block_id += 1
                continue
            parsed = parse_anchor_line(line)
            if not parsed:
                continue
            gene_a, gene_b, weight, is_lifted = parsed
            n_total += 1
            if exclude_lifted and is_lifted:
                n_skipped += 1
                continue
            if weight < min_score:
                n_skipped += 1
                continue
            edges.append(AnchorEdge(
                source=gene_a, target=gene_b, score=weight,
                is_lifted=is_lifted, species_pair=species_pair, block_id=block_id,
            ))
    return edges, n_total, n_skipped


def build_network(
        species_list: List[str], work_dir: Path,
        *, use_lifted: bool = True, min_score: float = 0, exclude_lifted: bool = False,
) -> Tuple[List[AnchorEdge], Set[str], NetworkStats]:
    edges = []
    nodes = set()
    stats = NetworkStats()
    suffix = ".lifted.anchors" if use_lifted else ".anchors"

    for i in range(len(species_list) - 1):
        sp_a, sp_b = species_list[i], species_list[i + 1]
        pair_name = f"{sp_a}.{sp_b}"
        anchor_file = work_dir / f"{pair_name}{suffix}"

        if not anchor_file.exists():
            alt_suffix = ".anchors" if use_lifted else ".lifted.anchors"
            alt_file = work_dir / f"{pair_name}{alt_suffix}"
            if alt_file.exists():
                log_warn(f"{anchor_file.name} not found, using {alt_file.name}")
                anchor_file = alt_file
            else:
                log_warn(f"File not found: {pair_name}{suffix}")
                continue

        pair_edges, n_total, n_skipped = parse_anchors_file(
            anchor_file, pair_name, min_score=min_score, exclude_lifted=exclude_lifted,
        )
        stats.filtered_by_score += n_skipped

        for e in pair_edges:
            nodes.add(e.source)
            nodes.add(e.target)
            if e.is_lifted:
                stats.lifted_edges += 1

        edges.extend(pair_edges)
        stats.species_pair_counts[pair_name] = len(pair_edges)
        log_info(f"{anchor_file.name}: {len(pair_edges)} edges (total {n_total}, skipped {n_skipped})")

    stats.total_nodes = len(nodes)
    stats.total_edges = len(edges)
    return edges, nodes, stats


def export_tsv(edges: List[AnchorEdge], output_file: Path):
    with open(output_file, 'w') as f:
        f.write("source\ttarget\tscore\tis_lifted\tspecies_pair\tblock_id\n")
        for e in edges:
            f.write(f"{e.source}\t{e.target}\t{e.score}\t"
                    f"{e.is_lifted}\t{e.species_pair}\t{e.block_id}\n")
    log_info(f"Exported: {output_file}")


def export_graphml(edges: List[AnchorEdge], nodes: Set[str], output_file: Path):
    if not _HAS_NX:
        log_warn("networkx not installed, skip GraphML")
        return
    G = nx.Graph()
    G.add_nodes_from(nodes)
    for e in edges:
        G.add_edge(e.source, e.target, weight=e.score, is_lifted=e.is_lifted,
                    species_pair=e.species_pair, block_id=e.block_id)
    nx.write_graphml(G, output_file)
    log_info(f"Exported: {output_file}")


def export_gexf(edges: List[AnchorEdge], nodes: Set[str], output_file: Path):
    if not _HAS_NX:
        log_warn("networkx not installed, skip GEXF")
        return
    G = nx.Graph()
    G.add_nodes_from(nodes)
    for e in edges:
        G.add_edge(e.source, e.target, weight=e.score, is_lifted=e.is_lifted,
                    species_pair=e.species_pair, block_id=e.block_id)
    nx.write_gexf(G, output_file)
    log_info(f"Exported: {output_file}")


def export_stats(stats: NetworkStats, output_file: Path):
    with open(output_file, 'w') as f:
        f.write("# SynNet Network Statistics\n\n")
        f.write(f"total_nodes: {stats.total_nodes}\n")
        f.write(f"total_edges: {stats.total_edges}\n")
        f.write(f"lifted_edges: {stats.lifted_edges}\n")
        f.write(f"filtered_by_score: {stats.filtered_by_score}\n")
        f.write(f"\n# Per-species-pair edge counts\n")
        for pair, count in stats.species_pair_counts.items():
            f.write(f"  {pair}: {count}\n")
    log_info(f"Exported: {output_file}")


def main():
    parser = argparse.ArgumentParser(
        description="Step 3: Build synteny network from .anchors files",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python network.py -s species.lst -d .                       # default: use .lifted.anchors
  python network.py -s species.lst --use-anchors               # use .anchors instead
  python network.py -s species.lst --min-score 100             # filter by alignment score
  python network.py -s species.lst --exclude-lifted            # exclude lifted anchor edges
  python network.py -s species.lst -o Final_Network --formats tsv,graphml
        """,
    )

    parser.add_argument("-s", "--species-list", required=True,
                        help="Species list file (chain order)")
    parser.add_argument("-d", "--work-dir", default=".",
                        help="Directory containing .anchors files (default: current dir)")
    parser.add_argument("--use-anchors", action="store_true",
                        help="Use .anchors instead of .lifted.anchors")
    parser.add_argument("--exclude-lifted", action="store_true",
                        help="Exclude lifted edges (marked with 'L' suffix)")
    parser.add_argument("--min-score", type=float, default=0,
                        help="Minimum alignment score cutoff (default: 0, no filter)")
    parser.add_argument("-o", "--output-prefix", default="Final_Network",
                        help="Output file prefix (default: Final_Network)")
    parser.add_argument("--formats", type=str, default="tsv",
                        help="Output formats: tsv,graphml,gexf (comma-separated, default: tsv)")
    parser.add_argument("-v", "--verbose", action="store_true",
                        help="Verbose output")

    args = parser.parse_args()

    log_info("Step 3: SynNet Network Builder")

    species = load_species_list(args.species_list)
    log_info(f"Loaded {len(species)} species: {' -> '.join(species)}")

    work_dir = Path(args.work_dir)
    if not work_dir.exists():
        log_error(f"Work directory not found: {work_dir}")
        sys.exit(1)

    use_lifted = not args.use_anchors
    log_info(f"Reading {'lifted ' if use_lifted else ''}anchors from {work_dir}")

    edges, nodes, stats = build_network(
        species, work_dir,
        use_lifted=use_lifted,
        min_score=args.min_score,
        exclude_lifted=args.exclude_lifted,
    )

    if not edges:
        log_error("No edges found. Check input files.")
        sys.exit(1)

    log_info(f"\nNetwork: {stats.total_nodes} nodes, {stats.total_edges} edges "
             f"({stats.lifted_edges} lifted)")
    if stats.filtered_by_score:
        log_info(f"Filtered by score < {args.min_score}: {stats.filtered_by_score} edges")

    formats = [f.strip() for f in args.formats.split(',')]
    prefix = args.output_prefix

    if "tsv" in formats:
        export_tsv(edges, Path(f"{prefix}.tsv"))
    if "graphml" in formats:
        export_graphml(edges, nodes, Path(f"{prefix}.graphml"))
    if "gexf" in formats:
        export_gexf(edges, nodes, Path(f"{prefix}.gexf"))

    export_stats(stats, Path(f"{prefix}.stats.txt"))

    log_info(f"Done! Output: {prefix}.*")


if __name__ == "__main__":
    main()
