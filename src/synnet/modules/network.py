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
    import networkx as nx
    _HAS_NX = True
except ImportError:
    _HAS_NX = False


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


def parse_anchors_file(
        filepath: Path,
        species_pair: str,
        *,
        min_score: float = 0,
        exclude_lifted: bool = False,
) -> Tuple[List[AnchorEdge], int, int]:
    edges = []
    n_total = 0
    n_skipped = 0

    for gene_a, gene_b, weight, is_lifted, block_id in read_anchors_file(
            filepath, min_score=min_score, exclude_lifted=exclude_lifted):
        n_total += 1
        edges.append(AnchorEdge(
            source=gene_a,
            target=gene_b,
            score=weight,
            is_lifted=is_lifted,
            species_pair=species_pair,
            block_id=block_id,
        ))

    return edges, n_total, n_skipped


def build_network(
        species_list: List[str],
        work_dir: Path,
        *,
        use_lifted: bool = True,
        min_score: float = 0,
        exclude_lifted: bool = False,
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
                warning(f"{anchor_file.name} not found, using {alt_file.name}")
                anchor_file = alt_file
            else:
                warning(f"File not found: {pair_name}{suffix}")
                continue

        debug(f"Reading: {anchor_file}")
        pair_edges, n_total, n_skipped = parse_anchors_file(
            anchor_file, pair_name,
            min_score=min_score,
            exclude_lifted=exclude_lifted,
        )
        stats.filtered_by_score += n_skipped

        for e in pair_edges:
            nodes.add(e.source)
            nodes.add(e.target)
            if e.is_lifted:
                stats.lifted_edges += 1

        edges.extend(pair_edges)
        stats.species_pair_counts[pair_name] = len(pair_edges)
        info(f"{anchor_file.name}: {len(pair_edges)} edges (total {n_total}, skipped {n_skipped})")

    stats.total_nodes = len(nodes)
    stats.total_edges = len(edges)

    return edges, nodes, stats


def export_tsv(edges: List[AnchorEdge], output_file: Path):
    edge_tuples = [
        (e.source, e.target, e.score, e.is_lifted, e.species_pair, e.block_id)
        for e in edges
    ]
    write_network_tsv(output_file, edge_tuples)


def export_graphml(edges: List[AnchorEdge], nodes: Set[str], output_file: Path):
    if not _HAS_NX:
        warning("networkx not installed, skip GraphML")
        return

    G = nx.Graph()
    G.add_nodes_from(nodes)
    for e in edges:
        G.add_edge(e.source, e.target,
                   weight=e.score,
                   is_lifted=e.is_lifted,
                   species_pair=e.species_pair,
                   block_id=e.block_id)
    nx.write_graphml(G, output_file)
    info(f"Exported: {output_file}")


def export_gexf(edges: List[AnchorEdge], nodes: Set[str], output_file: Path):
    if not _HAS_NX:
        warning("networkx not installed, skip GEXF")
        return

    G = nx.Graph()
    G.add_nodes_from(nodes)
    for e in edges:
        G.add_edge(e.source, e.target,
                   weight=e.score,
                   is_lifted=e.is_lifted,
                   species_pair=e.species_pair,
                   block_id=e.block_id)
    nx.write_gexf(G, output_file)
    info(f"Exported: {output_file}")


def export_stats(stats: NetworkStats, output_file: Path):
    stats_dict = {
        "total_nodes": stats.total_nodes,
        "total_edges": stats.total_edges,
        "lifted_edges": stats.lifted_edges,
        "filtered_by_score": stats.filtered_by_score,
        "species_pair_counts": stats.species_pair_counts,
    }
    write_stats_file(output_file, stats_dict, title="SynNet Network Statistics")


def run_network(
        species_list_file: str,
        *,
        work_dir: str = ".",
        use_lifted: bool = True,
        min_score: float = 0,
        exclude_lifted: bool = False,
        output_prefix: str = "Final_Network",
        formats: str = "tsv",
) -> dict:
    info("SynNet Network Builder")

    species = read_species_list(species_list_file)
    info(f"Loaded {len(species)} species: {' -> '.join(species)}")

    wd = Path(work_dir)
    if not wd.exists():
        error(f"Work directory not found: {wd}")
        return {"success": False, "error": f"Work directory not found: {wd}"}

    info(f"Reading {'lifted ' if use_lifted else ''}anchors from {wd}")

    edges, nodes, stats = build_network(
        species, wd,
        use_lifted=use_lifted,
        min_score=min_score,
        exclude_lifted=exclude_lifted,
    )

    if not edges:
        error("No edges found. Check input files.")
        return {"success": False, "error": "No edges found"}

    info(f"\nNetwork: {stats.total_nodes} nodes, {stats.total_edges} edges "
         f"({stats.lifted_edges} lifted)")
    if stats.filtered_by_score:
        info(f"Filtered by score < {min_score}: {stats.filtered_by_score} edges")

    fmt_list = [f.strip() for f in formats.split(',')]
    prefix = output_prefix

    if "tsv" in fmt_list:
        export_tsv(edges, Path(f"{prefix}.tsv"))
    if "graphml" in fmt_list:
        export_graphml(edges, nodes, Path(f"{prefix}.graphml"))
    if "gexf" in fmt_list:
        export_gexf(edges, nodes, Path(f"{prefix}.gexf"))

    export_stats(stats, Path(f"{prefix}.stats.txt"))

    info("Network build completed!")

    return {
        "success": True,
        "stats": {
            "nodes": stats.total_nodes,
            "edges": stats.total_edges,
            "lifted_edges": stats.lifted_edges,
        },
    }
