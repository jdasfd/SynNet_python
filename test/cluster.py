"""
test/cluster.py - Step 4: Cluster and filter synteny network

Independent implementation - no module imports.

Usage:
    python cluster.py -i Final_Network.tsv -s species.lst -d ../seqs
    python cluster.py -i Final_Network.tsv -s species.lst --cluster-num 4
    python cluster.py -i Final_Network.tsv -s species.lst --min-species 2
    python cluster.py -i Final_Network.tsv -s species.lst --ortholog-only
    python cluster.py -i Final_Network.tsv -s species.lst --cluster louvain --cluster-num 3
"""

import sys
import argparse
from pathlib import Path
from typing import List, Dict, Set, Optional, Tuple
from dataclasses import dataclass
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
class ClusterResult:
    clusters: List[Set[str]]
    method: str
    n_before_filter: int = 0
    filtered_by_size: int = 0
    filtered_by_species: int = 0
    filtered_by_ortholog: int = 0


def load_species_list(list_file: str) -> List[str]:
    with open(list_file, 'r') as f:
        species = [line.strip() for line in f if line.strip() and not line.startswith('#')]
    if len(species) < 2:
        raise ValueError(f"Need >= 2 species, got {len(species)}")
    return species


def build_gene_species_map(species_list: List[str], work_dir: str = ".") -> Dict[str, str]:
    gene_map = {}
    wd = Path(work_dir)
    for sp in species_list:
        bed_file = wd / f"{sp}.bed"
        if not bed_file.exists():
            continue
        with open(bed_file, 'r') as f:
            for line in f:
                parts = line.rstrip('\n').split('\t')
                if len(parts) >= 4:
                    gene_map[parts[3]] = sp
    return gene_map


def infer_species(gene_id: str, gene_species_map: Dict[str, str]) -> Optional[str]:
    if gene_id in gene_species_map:
        return gene_species_map[gene_id]
    base = gene_id.split('.')[0]
    if base in gene_species_map:
        return gene_species_map[base]
    return None


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
    log_warn("MCL not available, falling back to connected components (install 'mcl' for MCL)")
    return cluster_connected_components(edges, nodes)


def cluster_louvain(edges: List[Tuple[str, str, float]], nodes: Set[str]) -> List[Set[str]]:
    if not _HAS_NX:
        log_error("networkx required for Louvain clustering")
        return []
    try:
        from networkx.algorithms.community import louvain_communities
    except ImportError:
        log_warn("Louvain not available, falling back to connected components")
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
            sp = infer_species(gene, gene_species_map)
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
            sp = infer_species(gene, gene_species_map)
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
    log_info(f"Clustering method: {method}")

    clusters = run_clustering(edges, nodes, method=method, mcl_inflation=mcl_inflation)
    n_before = len(clusters)
    log_info(f"Raw clusters: {n_before}")

    filtered_by_size = 0
    filtered_by_species = 0
    filtered_by_ortholog = 0

    if min_cluster_size > 1:
        clusters, n = filter_by_min_size(clusters, min_cluster_size)
        filtered_by_size = n
        if n:
            log_info(f"Filtered {n} clusters smaller than {min_cluster_size}")

    if min_species_count > 1 or ortholog_only:
        if gene_species_map is None:
            log_warn("No gene-species mapping, species-based filtering disabled")
            log_warn("Provide -s and -d (with .bed files) to enable it")
        else:
            if min_species_count > 1:
                clusters, n = filter_by_min_species(clusters, gene_species_map, min_species_count)
                filtered_by_species = n
                if n:
                    log_info(f"Filtered {n} clusters with < {min_species_count} species")
            if ortholog_only:
                clusters, n = filter_by_ortholog(clusters, gene_species_map, require_ortholog=True)
                filtered_by_ortholog = n
                if n:
                    log_info(f"Filtered {n} non-ortholog clusters (1-to-1 required)")

    if clusters:
        sizes = [len(c) for c in clusters]
        log_info(f"Final clusters: {len(clusters)}, largest={max(sizes)}, smallest={min(sizes)}")
    else:
        log_warn("No clusters remaining after filtering")

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
                    sp = infer_species(gene, gene_species_map)
                    if sp:
                        species_in_cluster.add(sp)
                sp_count = len(species_in_cluster)
                for gene in sorted(cluster):
                    sp = infer_species(gene, gene_species_map) or "unknown"
                    f.write(f"{gene}\t{cluster_label}\t{size}\t{sp}\t{sp_count}\n")
            else:
                for gene in sorted(cluster):
                    f.write(f"{gene}\t{cluster_label}\t{size}\n")

    log_info(f"Exported: {output_file} ({len(clusters)} clusters)")


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
                    sp = infer_species(gene, gene_species_map)
                    if sp:
                        species_in_cluster.add(sp)
                sp_str = ",".join(sorted(species_in_cluster))
                sp_count = len(species_in_cluster)
                f.write(f"{cluster_label}\t{size}\t{sp_str}\t{sp_count}\n")
            else:
                f.write(f"{cluster_label}\t{size}\n")

    log_info(f"Exported: {output_file} ({len(clusters)} clusters)")


def main():
    parser = argparse.ArgumentParser(
        description="Step 4: Cluster and filter synteny network",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python cluster.py -i Final_Network.tsv -s species.lst -d ../seqs
  python cluster.py -i Final_Network.tsv -s species.lst --cluster-num 4
  python cluster.py -i Final_Network.tsv -s species.lst --min-species 2
  python cluster.py -i Final_Network.tsv -s species.lst --ortholog-only
  python cluster.py -i Final_Network.tsv -s species.lst --cluster louvain
        """,
    )

    parser.add_argument("-i", "--input", required=True,
                        help="Network TSV file (from 'network' command)")
    parser.add_argument("-s", "--species-list", required=True,
                        help="Species list file (for species-based filtering)")
    parser.add_argument("-d", "--bed-dir", default=".",
                        help="Directory containing .bed files (default: current dir)")
    parser.add_argument("-o", "--output-prefix", default="Filtered",
                        help="Output file prefix (default: Filtered)")
    parser.add_argument("--cluster", type=str, default="cc",
                        choices=["cc", "mcl", "louvain"],
                        help="Clustering method: cc, mcl, louvain (default: cc)")
    parser.add_argument("--mcl-inflation", type=float, default=2.0,
                        help="MCL inflation parameter (default: 2.0)")
    parser.add_argument("--cluster-num", type=int, default=2,
                        help="Minimum cluster size to keep (default: 2)")
    parser.add_argument("--min-species", type=int, default=1,
                        help="Minimum number of species in a cluster (default: 1)")
    parser.add_argument("--ortholog-only", action="store_true",
                        help="Keep only 1-to-1 ortholog clusters")
    parser.add_argument("-v", "--verbose", action="store_true",
                        help="Verbose output")

    args = parser.parse_args()

    log_info("Step 4: SynNet Cluster & Filter")

    species = load_species_list(args.species_list)
    log_info(f"Loaded {len(species)} species: {' -> '.join(species)}")

    gene_species_map = build_gene_species_map(species, args.bed_dir)
    if gene_species_map:
        log_info(f"Built gene-species map: {len(gene_species_map)} genes from {len(set(gene_species_map.values()))} species")
    else:
        log_warn(f"No .bed files found in {args.bed_dir}, species-based filtering will be disabled")

    edges = []
    nodes = set()
    with open(args.input, 'r') as f:
        header = f.readline()
        for line in f:
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 3:
                continue
            src, tgt = parts[0], parts[1]
            try:
                score = float(parts[2])
            except ValueError:
                score = 0.0
            edges.append((src, tgt, score))
            nodes.add(src)
            nodes.add(tgt)

    log_info(f"Loaded {len(edges)} edges, {len(nodes)} nodes from {args.input}")

    result = run_cluster(
        edges, nodes, species,
        method=args.cluster,
        mcl_inflation=args.mcl_inflation,
        min_cluster_size=args.cluster_num,
        min_species_count=args.min_species,
        ortholog_only=args.ortholog_only,
        gene_species_map=gene_species_map,
    )

    if result.clusters:
        export_clusters(result.clusters, Path(f"{args.output_prefix}.clusters.tsv"),
                        gene_species_map=gene_species_map)
        export_cluster_summary(result.clusters, Path(f"{args.output_prefix}.clusters.summary.tsv"),
                               gene_species_map=gene_species_map)

    log_info(f"\nFiltering summary:")
    log_info(f"  Raw clusters: {result.n_before_filter}")
    log_info(f"  Filtered by size < {args.cluster_num}: {result.filtered_by_size}")
    log_info(f"  Filtered by species < {args.min_species}: {result.filtered_by_species}")
    log_info(f"  Filtered by ortholog: {result.filtered_by_ortholog}")
    log_info(f"  Final clusters: {len(result.clusters)}")

    log_info(f"\nDone! Output: {args.output_prefix}.*")


if __name__ == "__main__":
    main()
