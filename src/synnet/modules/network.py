import sys
from pathlib import Path
from typing import List, Dict, Set, Optional, Tuple
from dataclasses import dataclass, field

from synnet.utils.logger import get_logger, info, warning, error, debug
from synnet.utils.io import (
    read_species_list,
    read_anchors_file,
    write_network_tsv,
    write_stats_file,
    ensure_dir,
)

logger = get_logger(__name__)

try:
    import networkx as nx # type: ignore
    _HAS_NX = True
except ImportError:
    _HAS_NX = False


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
                warning(f"{anchor_file.name} not found, using {alt_file.name}")
                anchor_file = alt_file
            else:
                warning(f"File not found: {pair_name}{suffix}")
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
        info(f"{anchor_file.name}: {len(pair_edges)} edges (total {n_total}, filtered {n_filtered})")

    stats.total_edges = len(edges)
    stats.total_nodes = len(nodes)

    return edges, nodes, stats


def write_network_tsv_output(
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


def run_network(
        species_list_file: str,
        *,
        input_dir: str = "jcvi_output",
        output_dir: str = "network_output",
        no_lifted: bool = False,
        min_score: int = 0,
) -> dict:
    info("Build Synteny Network")

    species = read_species_list(species_list_file)
    info(f"Species: {' -> '.join(species)}")

    in_dir = Path(input_dir)
    if not in_dir.exists():
        error(f"Input directory not found: {in_dir}")
        return {"success": False, "error": f"Input directory not found: {in_dir}"}

    out_dir = Path(output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    include_lifted = not no_lifted
    info(f"Input dir: {in_dir}")
    info(f"Output dir: {out_dir}")
    info(f"Include lifted: {include_lifted}")
    info(f"Min score: {min_score}")

    edges, nodes, stats = build_network(
        species, in_dir,
        include_lifted=include_lifted,
        min_score=min_score,
    )

    output_tsv = out_dir / "Final_Network.tsv"
    output_stats = out_dir / "Final_Network.stats.txt"

    write_network_tsv_output(edges, output_tsv)
    write_stats_txt(stats, output_stats)

    info(f"\nNetwork: {stats.total_nodes} nodes, {stats.total_edges} edges ({stats.lifted_edges} lifted)")
    info(f"Exported: {output_tsv}")
    info(f"Exported: {output_stats}")
    info("Done!")

    return {
        "success": True,
        "stats": {
            "nodes": stats.total_nodes,
            "edges": stats.total_edges,
            "lifted_edges": stats.lifted_edges,
        },
    }
