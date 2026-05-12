import sys
from pathlib import Path
from typing import List, Dict, Set, Optional, Tuple
from dataclasses import dataclass, field
from collections import defaultdict

from synnet.utils.logger import get_logger, info, warning, error, debug
from synnet.utils.io import build_gene_species_map as io_build_gene_species_map

logger = get_logger(__name__)

try:
    import networkx as nx
    _HAS_NX = True
except ImportError:
    _HAS_NX = False


@dataclass
class ClusterResult:
    clusters: List[Set[str]]
    method: str
    n_before_filter: int = 0
    filtered_by_size: int = 0
    filtered_by_species: int = 0
    filtered_by_ortholog: int = 0


def cluster_connected_components(edges: List[Tuple[str, str, float]], nodes: Set[str]) -> List[Set[str]]:
    if not _HAS_NX:
        adj = defaultdict(set)
        for src, tgt, _ in edges:
            adj[src].add(tgt)
            adj[tgt].add(src)

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
    for src, tgt, _ in edges:
        G.add_edge(src, tgt)
    return sorted(nx.connected_components(G), key=len, reverse=True)


def cluster_mcl(edges: List[Tuple[str, str, float]], nodes: Set[str], inflation: float = 2.0) -> List[Set[str]]:
    if not _HAS_NX:
        error("networkx required for MCL-like clustering")
        return []

    info("Falling back to connected components (install 'mcl' for MCL clustering)")
    return cluster_connected_components(edges, nodes)


def cluster_louvain(edges: List[Tuple[str, str, float]], nodes: Set[str]) -> List[Set[str]]:
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
    for src, tgt, score in edges:
        if G.has_edge(src, tgt):
            G[src][tgt]['weight'] += score
        else:
            G.add_edge(src, tgt, weight=score)

    communities = louvain_communities(G, weight='weight')
    return sorted(communities, key=len, reverse=True)


def run_clustering(
        edges: List[Tuple[str, str, float]],
        nodes: Set[str],
        method: str = "cc",
        mcl_inflation: float = 2.0,
) -> List[Set[str]]:
    if method == "cc":
        return cluster_connected_components(edges, nodes)
    elif method == "mcl":
        return cluster_mcl(edges, nodes, inflation=mcl_inflation)
    elif method == "louvain":
        return cluster_louvain(edges, nodes)
    else:
        return cluster_connected_components(edges, nodes)


def build_gene_species_map(species_list: List[str], work_dir: str = ".") -> Dict[str, str]:
    return io_build_gene_species_map(species_list, Path(work_dir))


def infer_species_from_map(gene_id: str, gene_species_map: Dict[str, str]) -> Optional[str]:
    if gene_id in gene_species_map:
        return gene_species_map[gene_id]
    base = gene_id.split('.')[0]
    if base in gene_species_map:
        return gene_species_map[base]
    return None


def filter_by_min_size(clusters: List[Set[str]], min_size: int) -> Tuple[List[Set[str]], int]:
    if min_size <= 1:
        return clusters, 0
    filtered = [c for c in clusters if len(c) >= min_size]
    return filtered, len(clusters) - len(filtered)


def filter_by_min_species(
        clusters: List[Set[str]],
        gene_species_map: Dict[str, str],
        min_species: int,
) -> Tuple[List[Set[str]], int]:
    if min_species <= 1:
        return clusters, 0

    filtered = []
    for cluster in clusters:
        species_in_cluster = set()
        for gene in cluster:
            sp = infer_species_from_map(gene, gene_species_map)
            if sp:
                species_in_cluster.add(sp)
        if len(species_in_cluster) >= min_species:
            filtered.append(cluster)

    return filtered, len(clusters) - len(filtered)


def filter_by_ortholog(
        clusters: List[Set[str]],
        gene_species_map: Dict[str, str],
        require_ortholog: bool = True,
) -> Tuple[List[Set[str]], int]:
    if not require_ortholog:
        return clusters, 0

    filtered = []
    for cluster in clusters:
        species_genes = defaultdict(set)
        for gene in cluster:
            sp = infer_species_from_map(gene, gene_species_map)
            if sp:
                species_genes[sp].add(gene)

        is_ortholog = all(len(genes) == 1 for genes in species_genes.values())
        if is_ortholog:
            filtered.append(cluster)

    return filtered, len(clusters) - len(filtered)


def run_cluster(
        edges: List[Tuple[str, str, float]],
        nodes: Set[str],
        species_list: List[str],
        *,
        method: str = "cc",
        mcl_inflation: float = 2.0,
        min_cluster_size: int = 2,
        min_species_count: int = 1,
        ortholog_only: bool = False,
        gene_species_map: Optional[Dict[str, str]] = None,
) -> ClusterResult:
    info(f"Clustering method: {method}")

    clusters = run_clustering(edges, nodes, method=method, mcl_inflation=mcl_inflation)
    n_before = len(clusters)
    info(f"Raw clusters: {n_before}")

    filtered_by_size = 0
    filtered_by_species = 0
    filtered_by_ortholog = 0

    if min_cluster_size > 1:
        clusters, n = filter_by_min_size(clusters, min_cluster_size)
        filtered_by_size = n
        if n:
            info(f"Filtered {n} clusters smaller than {min_cluster_size}")

    if min_species_count > 1 or ortholog_only:
        if gene_species_map is None:
            warning("No gene-species mapping provided, species-based filtering disabled")
            warning("Provide --species-list and --work-dir (with .bed files) to enable it")
        else:
            if min_species_count > 1:
                clusters, n = filter_by_min_species(clusters, gene_species_map, min_species_count)
                filtered_by_species = n
                if n:
                    info(f"Filtered {n} clusters with < {min_species_count} species")

            if ortholog_only:
                clusters, n = filter_by_ortholog(clusters, gene_species_map, require_ortholog=True)
                filtered_by_ortholog = n
                if n:
                    info(f"Filtered {n} non-ortholog clusters (1-to-1 required)")

    if clusters:
        sizes = [len(c) for c in clusters]
        info(f"Final clusters: {len(clusters)}, "
             f"largest={max(sizes)}, smallest={min(sizes)}")
    else:
        warning("No clusters remaining after filtering")

    return ClusterResult(
        clusters=clusters,
        method=method,
        n_before_filter=n_before,
        filtered_by_size=filtered_by_size,
        filtered_by_species=filtered_by_species,
        filtered_by_ortholog=filtered_by_ortholog,
    )


def export_clusters(clusters: List[Set[str]], output_file: Path,
                     gene_species_map: Optional[Dict[str, str]] = None):
    with open(output_file, 'w') as f:
        header = "gene_id\tcluster_id\tcluster_size"
        has_species = gene_species_map is not None
        if has_species:
            header += "\tspecies\tspecies_count"
        f.write(header + "\n")

        for cid, cluster in enumerate(clusters, 1):
            size = len(cluster)
            cluster_label = f"CL{cid}"

            if has_species:
                species_in_cluster = set()
                for gene in cluster:
                    sp = infer_species_from_map(gene, gene_species_map)
                    if sp:
                        species_in_cluster.add(sp)
                sp_str = ",".join(sorted(species_in_cluster))
                sp_count = len(species_in_cluster)
                for gene in sorted(cluster):
                    sp = infer_species_from_map(gene, gene_species_map) or "unknown"
                    f.write(f"{gene}\t{cluster_label}\t{size}\t{sp}\t{sp_count}\n")
            else:
                for gene in sorted(cluster):
                    f.write(f"{gene}\t{cluster_label}\t{size}\n")

    info(f"Exported: {output_file} ({len(clusters)} clusters)")


def export_cluster_summary(clusters: List[Set[str]], output_file: Path,
                            gene_species_map: Optional[Dict[str, str]] = None):
    with open(output_file, 'w') as f:
        header = "cluster_id\tcluster_size"
        has_species = gene_species_map is not None
        if has_species:
            header += "\tspecies\tspecies_count"
        f.write(header + "\n")

        for cid, cluster in enumerate(clusters, 1):
            size = len(cluster)
            cluster_label = f"CL{cid}"

            if has_species:
                species_in_cluster = set()
                for gene in cluster:
                    sp = infer_species_from_map(gene, gene_species_map)
                    if sp:
                        species_in_cluster.add(sp)
                sp_str = ",".join(sorted(species_in_cluster))
                sp_count = len(species_in_cluster)
                f.write(f"{cluster_label}\t{size}\t{sp_str}\t{sp_count}\n")
            else:
                f.write(f"{cluster_label}\t{size}\n")

    info(f"Exported: {output_file} ({len(clusters)} clusters)")
