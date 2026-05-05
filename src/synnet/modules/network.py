import sys
from pathlib import Path
from typing import List, Dict, Set, Optional, Tuple
from dataclasses import dataclass, field
from collections import defaultdict

from synnet.utils.logger import get_logger, info, warning, error, debug

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
    filtered_by_cluster_size: int = 0
    species_pair_counts: Dict[str, int] = field(default_factory=dict)
    n_clusters: int = 0
    largest_cluster: int = 0
    smallest_cluster: int = 0


def load_species_list(list_file: str) -> List[str]:
    with open(list_file, 'r') as f:
        species = [line.strip() for line in f
                   if line.strip() and not line.startswith('#')]
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
        filepath: Path,
        species_pair: str,
        *,
        min_score: float = 0,
        exclude_lifted: bool = False,
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


def cluster_connected_components(edges: List[AnchorEdge], nodes: Set[str]) -> List[Set[str]]:
    if not _HAS_NX:
        adj = defaultdict(set)
        for e in edges:
            adj[e.source].add(e.target)
            adj[e.target].add(e.source)

        visited = set()
        clusters = []
        for node in nodes:
            if node in visited:
                continue
            cluster = set()
            stack = [node]
            while stack:
                n = stack.pop()
                if n in visited:
                    continue
                visited.add(n)
                cluster.add(n)
                stack.extend(adj[n] - visited)
            clusters.append(cluster)
        return sorted(clusters, key=len, reverse=True)

    G = nx.Graph()
    G.add_nodes_from(nodes)
    for e in edges:
        G.add_edge(e.source, e.target)
    return sorted(nx.connected_components(G), key=len, reverse=True)


def cluster_mcl(edges: List[AnchorEdge], nodes: Set[str], inflation: float = 2.0) -> List[Set[str]]:
    if not _HAS_NX:
        error("networkx required for MCL-like clustering")
        return []

    info("Falling back to connected components (install 'mcl' for MCL clustering)")
    return cluster_connected_components(edges, nodes)


def cluster_louvain(edges: List[AnchorEdge], nodes: Set[str]) -> List[Set[str]]:
    if not _HAS_NX:
        error("networkx required for Louvain clustering")
        return []

    try:
        from networkx.algorithms.community import louvain_communities
    except ImportError:
        info("Louvain not available, falling back to connected components")
        return cluster_connected_components(edges, nodes)

    G = nx.Graph()
    G.add_nodes_from(nodes)
    for e in edges:
        if G.has_edge(e.source, e.target):
            G[e.source][e.target]['weight'] += e.score
        else:
            G.add_edge(e.source, e.target, weight=e.score)

    communities = louvain_communities(G, weight='weight')
    return sorted(communities, key=len, reverse=True)


def export_tsv(edges: List[AnchorEdge], output_file: Path):
    with open(output_file, 'w') as f:
        f.write("source\ttarget\tscore\tis_lifted\tspecies_pair\tblock_id\n")
        for e in edges:
            f.write(f"{e.source}\t{e.target}\t{e.score}\t"
                    f"{e.is_lifted}\t{e.species_pair}\t{e.block_id}\n")
    info(f"Exported: {output_file}")


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


def export_clusters(clusters: List[Set[str]], output_file: Path):
    with open(output_file, 'w') as f:
        f.write("gene_id\tcluster_id\tcluster_size\n")
        for cid, cluster in enumerate(clusters, 1):
            size = len(cluster)
            for gene in sorted(cluster):
                f.write(f"{gene}\tCLUSTER_{cid:05d}\t{size}\n")
    info(f"Exported: {output_file} ({len(clusters)} clusters)")


def export_stats(stats: NetworkStats, clusters: List[Set[str]], output_file: Path):
    if clusters:
        sizes = [len(c) for c in clusters]
        stats.n_clusters = len(clusters)
        stats.largest_cluster = max(sizes)
        stats.smallest_cluster = min(sizes)

    with open(output_file, 'w') as f:
        f.write("# SynNet Network Statistics\n\n")
        f.write(f"total_nodes: {stats.total_nodes}\n")
        f.write(f"total_edges: {stats.total_edges}\n")
        f.write(f"lifted_edges: {stats.lifted_edges}\n")
        f.write(f"filtered_by_score: {stats.filtered_by_score}\n")
        f.write(f"filtered_by_cluster_size: {stats.filtered_by_cluster_size}\n")
        f.write(f"total_clusters: {stats.n_clusters}\n")
        f.write(f"largest_cluster: {stats.largest_cluster}\n")
        f.write(f"smallest_cluster: {stats.smallest_cluster}\n")
        f.write(f"\n# Per-species-pair edge counts\n")
        for pair, count in stats.species_pair_counts.items():
            f.write(f"  {pair}: {count}\n")

    info(f"Exported: {output_file}")


def run_network(
        species_list_file: str,
        *,
        work_dir: str = ".",
        use_lifted: bool = True,
        min_score: float = 0,
        exclude_lifted: bool = False,
        output_prefix: str = "Final_Network",
        formats: str = "tsv",
        cluster_method: str = "cc",
        mcl_inflation: float = 2.0,
        min_cluster_size: int = 2,
) -> dict:
    info("SynNet Network Builder")

    species = load_species_list(species_list_file)
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

    info(f"\nClustering method: {cluster_method}")
    if cluster_method == "cc":
        clusters = cluster_connected_components(edges, nodes)
    elif cluster_method == "mcl":
        clusters = cluster_mcl(edges, nodes, inflation=mcl_inflation)
    elif cluster_method == "louvain":
        clusters = cluster_louvain(edges, nodes)
    else:
        clusters = cluster_connected_components(edges, nodes)

    if clusters:
        n_before = len(clusters)
        clusters = [c for c in clusters if len(c) >= min_cluster_size]
        n_filtered = n_before - len(clusters)
        if n_filtered:
            info(f"Filtered {n_filtered} clusters smaller than {min_cluster_size}")
            stats.filtered_by_cluster_size = n_filtered
        export_clusters(clusters, Path(f"{prefix}.clusters.tsv"))
        if clusters:
            sizes = [len(c) for c in clusters]
            info(f"Clusters: {len(clusters)} total, "
                 f"largest={max(sizes)}, smallest={min(sizes)}")
        else:
            warning("No clusters remaining after filtering")

    export_stats(stats, clusters, Path(f"{prefix}.stats.txt"))

    info("Completed!")

    return {
        "success": True,
        "stats": {
            "nodes": stats.total_nodes,
            "edges": stats.total_edges,
            "lifted_edges": stats.lifted_edges,
            "clusters": stats.n_clusters,
        },
    }
