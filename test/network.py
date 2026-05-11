"""
test/network.py - Step 3: Build synteny network from anchors files

Input: .anchors or .lifted.anchors files from mcscan output
Output: Final_Network.tsv (node1, node2, score, is_lifted)

Usage:
    python network.py -s species.lst
    python network.py -s species.lst -i jcvi_output -o network_output
    python network.py -s species.lst --no-lifted
    python network.py -s species.lst --min-score 0
"""

import sys
import argparse
from pathlib import Path
from typing import List, Tuple, Set, Dict
from dataclasses import dataclass, field


def log_info(msg):
    print(f"[INFO] {msg}", file=sys.stderr)


def log_warn(msg):
    print(f"[WARN] {msg}", file=sys.stderr)


def log_error(msg):
    print(f"[ERROR] {msg}", file=sys.stderr)


@dataclass
class NetworkStats:
    total_edges: int = 0
    lifted_edges: int = 0
    filtered_edges: int = 0
    total_nodes: int = 0
    pair_counts: Dict[str, int] = field(default_factory=dict)


def parse_anchors_file(
        anchor_file: Path,
        pair_name: str,
        include_lifted: bool = True,
        min_score: int = 0,
) -> Tuple[List[Tuple[str, str, int, bool]], int, int]:
    edges = []
    n_total = 0
    n_filtered = 0

    with open(anchor_file, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue

            parts = line.split('\t')
            if len(parts) < 2:
                continue

            gene1, gene2 = parts[0], parts[1]
            score = 0
            is_lifted = False

            if len(parts) >= 3:
                score_str = parts[2]
                if score_str.endswith('L'):
                    is_lifted = True
                    score_str = score_str[:-1]
                try:
                    score = int(score_str)
                except ValueError:
                    score = 0

            n_total += 1

            if is_lifted and not include_lifted:
                n_filtered += 1
                continue

            if score < min_score:
                n_filtered += 1
                continue

            edges.append((gene1, gene2, score, is_lifted))

    return edges, n_total, n_filtered


def build_network(
        species_list: List[str],
        input_dir: Path,
        include_lifted: bool = True,
        min_score: int = 0,
) -> Tuple[List[Tuple[str, str, int, bool]], Set[str], NetworkStats]:
    edges = []
    nodes = set()
    stats = NetworkStats()

    suffix = ".lifted.anchors"

    for i in range(len(species_list) - 1):
        sp_a, sp_b = species_list[i], species_list[i + 1]
        pair_name = f"{sp_a}.{sp_b}"
        anchor_file = input_dir / f"{pair_name}{suffix}"

        if not anchor_file.exists():
            alt_file = input_dir / f"{pair_name}.anchors"
            if alt_file.exists():
                log_warn(f"{anchor_file.name} not found, using {alt_file.name}")
                anchor_file = alt_file
            else:
                log_warn(f"File not found: {pair_name}{suffix}")
                continue

        pair_edges, n_total, n_filtered = parse_anchors_file(
            anchor_file, pair_name, include_lifted=include_lifted, min_score=min_score
        )

        stats.filtered_edges += n_filtered

        for gene1, gene2, score, is_lifted in pair_edges:
            nodes.add(gene1)
            nodes.add(gene2)
            if is_lifted:
                stats.lifted_edges += 1

        edges.extend(pair_edges)
        stats.pair_counts[pair_name] = len(pair_edges)
        log_info(f"{anchor_file.name}: {len(pair_edges)} edges (total {n_total}, filtered {n_filtered})")

    stats.total_edges = len(edges)
    stats.total_nodes = len(nodes)

    return edges, nodes, stats


def write_network_tsv(
        edges: List[Tuple[str, str, int, bool]],
        output_file: Path,
) -> None:
    with open(output_file, 'w') as f:
        f.write("node1\tnode2\tscore\tis_lifted\n")
        for gene1, gene2, score, is_lifted in edges:
            f.write(f"{gene1}\t{gene2}\t{score}\t{is_lifted}\n")


def write_stats_txt(
        stats: NetworkStats,
        output_file: Path,
) -> None:
    with open(output_file, 'w') as f:
        f.write(f"Total nodes: {stats.total_nodes}\n")
        f.write(f"Total edges: {stats.total_edges}\n")
        f.write(f"Lifted edges: {stats.lifted_edges}\n")
        f.write(f"Filtered edges (by score): {stats.filtered_edges}\n")
        f.write("\nEdges per species pair:\n")
        for pair, count in sorted(stats.pair_counts.items()):
            f.write(f"  {pair}: {count}\n")


def main():
    parser = argparse.ArgumentParser(
        description="Step 3: Build synteny network from anchors files",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python network.py -s species.lst                           # default: jcvi_output -> network_output
  python network.py -s species.lst -i my_input -o my_output
  python network.py -s species.lst --no-lifted
  python network.py -s species.lst --min-score 0.5
        """,
    )

    parser.add_argument("-s", "--species-list", required=True,
                        help="Species list file (determines pairwise order)")
    parser.add_argument("-i", "--input-dir", default="jcvi_output",
                        help="Input directory containing anchors files (default: jcvi_output)")
    parser.add_argument("-o", "--output-dir", default="network_output",
                        help="Output directory for network files (default: network_output)")
    parser.add_argument("--no-lifted", action="store_true",
                        help="Exclude lifted alignments (rows with 'L' suffix in score)")
    parser.add_argument("--min-score", type=int, default=0,
                        help="Minimum score threshold (default: 0)")
    parser.add_argument("-v", "--verbose", action="store_true",
                        help="Verbose output")

    args = parser.parse_args()

    log_info("Step 3: Build Synteny Network")

    with open(args.species_list, 'r') as f:
        species = [line.strip() for line in f if line.strip() and not line.startswith('#')]

    if len(species) < 2:
        log_error("Species list must contain at least 2 species")
        sys.exit(1)

    log_info(f"Species: {' -> '.join(species)}")

    input_dir = Path(args.input_dir)
    if not input_dir.exists():
        log_error(f"Input directory not found: {input_dir}")
        sys.exit(1)

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    include_lifted = not args.no_lifted
    log_info(f"Input dir: {input_dir}")
    log_info(f"Output dir: {output_dir}")
    log_info(f"Include lifted: {include_lifted}")
    log_info(f"Min score: {args.min_score}")

    edges, nodes, stats = build_network(
        species, input_dir,
        include_lifted=include_lifted,
        min_score=args.min_score,
    )

    output_tsv = output_dir / "Final_Network.tsv"
    output_stats = output_dir / "Final_Network.stats.txt"

    write_network_tsv(edges, output_tsv)
    write_stats_txt(stats, output_stats)

    log_info(f"\nNetwork: {stats.total_nodes} nodes, {stats.total_edges} edges ({stats.lifted_edges} lifted)")
    log_info(f"Exported: {output_tsv}")
    log_info(f"Exported: {output_stats}")
    log_info("Done!")


if __name__ == "__main__":
    main()
