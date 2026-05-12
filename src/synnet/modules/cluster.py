import sys
from pathlib import Path
from typing import List, Dict, Set, Optional, Tuple
from dataclasses import dataclass, field
from collections import defaultdict

from synnet.utils.logger import get_logger, info, warning, error, debug

logger = get_logger(__name__)

_HAS_COMMUNITY = False
try:
    import community as community_louvain # type: ignore
    _HAS_COMMUNITY = True
except ImportError:
    pass

_HAS_NX = False
try:
    import networkx as nx # type: ignore
    _HAS_NX = True
except ImportError:
    pass


@dataclass
class ClusterStats:
    total_clusters: int = 0
    filtered_by_size: int = 0
    filtered_by_species: int = 0
    final_clusters: int = 0


def load_network(network_file: Path) -> Tuple[List[Tuple[str, str, int]], Set[str], Dict[Tuple[str, str], int]]:
    edges = []
    nodes = set()
    edge_weights = {}

    with open(network_file, 'r') as f:
        header = f.readline()
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) < 2:
                continue
            node1, node2 = parts[0], parts[1]
            try:
                score = int(parts[2]) if len(parts) >= 3 else 1
            except ValueError:
                score = 0
            edges.append((node1, node2, score))
            nodes.add(node1)
            nodes.add(node2)
            edge_weights[(node1, node2)] = score
            edge_weights[(node2, node1)] = score

    return edges, nodes, edge_weights


def cluster_connected_components(
        edges: List[Tuple[str, str, int]],
        nodes: Set[str],
) -> List[Set[str]]:
    parent = {n: n for n in nodes}

    def find(x):
        if parent[x] != x:
            parent[x] = find(parent[x])
        return parent[x]

    def union(x, y):
        px, py = find(x), find(y)
        if px != py:
            parent[px] = py

    for n1, n2, _ in edges:
        union(n1, n2)

    clusters = defaultdict(set)
    for n in nodes:
        clusters[find(n)].add(n)

    return list(clusters.values())


def cluster_louvain(
        edges: List[Tuple[str, str, int]],
        nodes: Set[str],
) -> List[Set[str]]:
    if not _HAS_NX or not _HAS_COMMUNITY:
        raise ImportError("networkx and python-louvain required for louvain method")

    G = nx.Graph()
    G.add_nodes_from(nodes)
    for n1, n2, w in edges:
        G.add_edge(n1, n2, weight=w)

    partition = community_louvain.best_partition(G, weight='weight')

    clusters = defaultdict(set)
    for node, comm_id in partition.items():
        clusters[comm_id].add(node)

    return list(clusters.values())


def cluster_infomap(
        edges: List[Tuple[str, str, int]],
        nodes: Set[str],
) -> List[Set[str]]:
    if not _HAS_NX:
        raise ImportError("networkx required for infomap method")

    try:
        import infomap # type: ignore
    except ImportError:
        raise ImportError("infomap package required for infomap method. Install: pip install infomap")

    G = nx.Graph()
    G.add_nodes_from(nodes)
    for n1, n2, w in edges:
        G.add_edge(n1, n2, weight=w)

    im = infomap.Infomap("--two-level --directed")
    for n1, n2, w in edges:
        im.addLink(n1, n2, w)
    im.run()

    clusters = defaultdict(set)
    for node in im.tree:
        if node.isLeaf:
            clusters[node.moduleIndex()].add(node.physicalId)

    node_to_cluster = {}
    for cluster_id, cluster_nodes in clusters.items():
        for node in cluster_nodes:
            node_to_cluster[node] = cluster_id

    result = []
    for cluster_id, cluster_nodes in clusters.items():
        valid_nodes = {n for n in cluster_nodes if n in nodes}
        if valid_nodes:
            result.append(valid_nodes)

    for node in nodes:
        if node not in node_to_cluster:
            result.append({node})

    return result


def cluster_label_propagation(
        edges: List[Tuple[str, str, int]],
        nodes: Set[str],
) -> List[Set[str]]:
    if not _HAS_NX:
        raise ImportError("networkx required for label_propagation method")

    G = nx.Graph()
    G.add_nodes_from(nodes)
    for n1, n2, _ in edges:
        G.add_edge(n1, n2)

    communities = nx.algorithms.community.label_propagation_communities(G)

    return [set(c) for c in communities]


def cluster_spectral(
        edges: List[Tuple[str, str, int]],
        nodes: Set[str],
        k: int = 10,
) -> List[Set[str]]:
    if not _HAS_NX:
        raise ImportError("networkx required for spectral method")

    try:
        from sklearn.cluster import SpectralClustering # type: ignore
    except ImportError:
        raise ImportError("scikit-learn required for spectral method. Install: pip install scikit-learn")

    G = nx.Graph()
    G.add_nodes_from(nodes)
    for n1, n2, w in edges:
        G.add_edge(n1, n2, weight=w)

    node_list = list(G.nodes())
    n_nodes = len(node_list)

    if n_nodes <= k:
        k = max(2, n_nodes // 2)

    adj_matrix = nx.to_numpy_array(G, nodelist=node_list, weight='weight')

    clustering = SpectralClustering(
        n_clusters=k,
        affinity='precomputed',
        assign_labels='kmeans',
        random_state=42
    )
    labels = clustering.fit_predict(adj_matrix)

    clusters = defaultdict(set)
    for i, label in enumerate(labels):
        clusters[label].add(node_list[i])

    return list(clusters.values())


def load_bed_files(bed_dir: Path, species_list: List[str]) -> Dict[str, str]:
    gene_to_species = {}
    for sp in species_list:
        bed_file = bed_dir / f"{sp}.bed"
        if not bed_file.exists():
            warning(f"BED file not found: {bed_file}")
            continue
        with open(bed_file, 'r') as f:
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) >= 4:
                    gene_id = parts[3]
                    gene_to_species[gene_id] = sp
    return gene_to_species


def count_species_in_cluster(
        genes: Set[str],
        gene_to_species: Dict[str, str],
) -> Tuple[int, Dict[str, int]]:
    species_counts = defaultdict(int)
    unknown_count = 0
    for gene in genes:
        sp = gene_to_species.get(gene, "unknown")
        if sp == "unknown":
            unknown_count += 1
        else:
            species_counts[sp] += 1
    if unknown_count > 0:
        species_counts["unknown"] = unknown_count
    return len(species_counts), dict(species_counts)


def filter_clusters(
        clusters: List[Set[str]],
        gene_to_species: Dict[str, str],
        min_cluster_size: int = 2,
        min_species: int = 1,
) -> Tuple[List[Set[str]], ClusterStats]:
    stats = ClusterStats(total_clusters=len(clusters))

    filtered = []
    for cluster in clusters:
        if len(cluster) < min_cluster_size:
            stats.filtered_by_size += 1
            continue

        n_species, _ = count_species_in_cluster(cluster, gene_to_species)
        if n_species < min_species:
            stats.filtered_by_species += 1
            continue

        filtered.append(cluster)

    stats.final_clusters = len(filtered)
    return filtered, stats


def write_clusters_tsv(
        clusters: List[Set[str]],
        output_file: Path,
) -> Dict[str, str]:
    gene_to_cluster = {}
    with open(output_file, 'w') as f:
        f.write("cluster_id\tgenes\n")
        for i, cluster in enumerate(clusters, 1):
            cluster_id = f"CL{i:07d}"
            genes_str = ",".join(sorted(cluster))
            f.write(f"{cluster_id}\t{genes_str}\n")
            for gene in cluster:
                gene_to_cluster[gene] = cluster_id
    return gene_to_cluster


def write_cluster_summary_tsv(
        clusters: List[Set[str]],
        gene_to_species: Dict[str, str],
        output_file: Path,
) -> None:
    with open(output_file, 'w') as f:
        f.write("cluster_id\tsize\tspecies_count\tspecies_composition\n")
        for i, cluster in enumerate(clusters, 1):
            cluster_id = f"CL{i:07d}"
            n_species, species_counts = count_species_in_cluster(cluster, gene_to_species)
            composition = ",".join(f"{sp}:{c}" for sp, c in sorted(species_counts.items()))
            f.write(f"{cluster_id}\t{len(cluster)}\t{n_species}\t{composition}\n")


def load_gene_list(gene_list_file: Path) -> Set[str]:
    genes = set()
    with open(gene_list_file, 'r') as f:
        for line in f:
            gene = line.strip()
            if gene and not gene.startswith('#'):
                genes.add(gene)
    return genes


def write_target_clusters_tsv(
        gene_to_cluster: Dict[str, str],
        query_genes: Set[str],
        edge_weights: Dict[Tuple[str, str], int],
        output_file: Path,
) -> None:
    target_cluster_ids = set()
    for gene in query_genes:
        if gene in gene_to_cluster:
            target_cluster_ids.add(gene_to_cluster[gene])

    cluster_to_genes = defaultdict(set)
    for gene, cid in gene_to_cluster.items():
        if cid in target_cluster_ids:
            cluster_to_genes[cid].add(gene)

    with open(output_file, 'w') as f:
        f.write("cluster_id\tscore\tnode1\tnode2\n")
        for cluster_id in sorted(target_cluster_ids):
            cluster_genes = cluster_to_genes[cluster_id]
            written_edges = set()
            for gene1 in cluster_genes:
                for gene2 in cluster_genes:
                    if gene1 >= gene2:
                        continue
                    edge_key = (gene1, gene2)
                    if edge_key in edge_weights:
                        score = edge_weights[edge_key]
                        f.write(f"{cluster_id}\t{score}\t{gene1}\t{gene2}\n")
                        written_edges.add(edge_key)

            if not written_edges and len(cluster_genes) >= 2:
                genes_list = sorted(cluster_genes)
                for i in range(len(genes_list) - 1):
                    f.write(f"{cluster_id}\t0\t{genes_list[i]}\t{genes_list[i+1]}\n")


def write_synnet_tsv(
        gene_to_cluster: Dict[str, str],
        query_genes: Set[str],
        edge_weights: Dict[Tuple[str, str], int],
        output_file: Path,
) -> None:
    target_cluster_ids = set()
    for gene in query_genes:
        if gene in gene_to_cluster:
            target_cluster_ids.add(gene_to_cluster[gene])

    cluster_to_genes = defaultdict(set)
    for gene, cid in gene_to_cluster.items():
        if cid in target_cluster_ids:
            cluster_to_genes[cid].add(gene)

    with open(output_file, 'w') as f:
        f.write("cluster_id\tscore\tnode1\tnode2\n")
        for cluster_id in sorted(target_cluster_ids):
            cluster_genes = cluster_to_genes[cluster_id]
            written_edges = set()
            for gene1 in cluster_genes:
                for gene2 in cluster_genes:
                    if gene1 >= gene2:
                        continue
                    edge_key = (gene1, gene2)
                    if edge_key in edge_weights:
                        score = edge_weights[edge_key]
                        f.write(f"{cluster_id}\t{score}\t{gene1}\t{gene2}\n")
                        written_edges.add(edge_key)

            if not written_edges and len(cluster_genes) >= 2:
                genes_list = sorted(cluster_genes)
                for i in range(len(genes_list) - 1):
                    f.write(f"{cluster_id}\t0\t{genes_list[i]}\t{genes_list[i+1]}\n")


def write_synnet_tsv_all(
        clusters: List[Set[str]],
        edge_weights: Dict[Tuple[str, str], int],
        output_file: Path,
) -> None:
    with open(output_file, 'w') as f:
        f.write("cluster_id\tscore\tnode1\tnode2\n")
        for i, cluster in enumerate(clusters, 1):
            cluster_id = f"CL{i:07d}"
            written_edges = set()
            for gene1 in cluster:
                for gene2 in cluster:
                    if gene1 >= gene2:
                        continue
                    edge_key = (gene1, gene2)
                    if edge_key in edge_weights:
                        score = edge_weights[edge_key]
                        f.write(f"{cluster_id}\t{score}\t{gene1}\t{gene2}\n")
                        written_edges.add(edge_key)

            if not written_edges and len(cluster) >= 2:
                genes_list = sorted(cluster)
                for j in range(len(genes_list) - 1):
                    f.write(f"{cluster_id}\t0\t{genes_list[j]}\t{genes_list[j+1]}\n")


def run_cluster(
        input_file: str,
        species_list_file: str,
        bed_dir: str,
        *,
        output_dir: str = "network_output",
        method: str = "cc",
        k: int = 10,
        cluster_size: int = 2,
        min_species: int = 1,
        gene_list: Optional[str] = None,
) -> dict:
    info("Step 4: Cluster Synteny Network")

    species_list_path = Path(species_list_file)
    with open(species_list_path, 'r') as f:
        species = [line.strip() for line in f if line.strip() and not line.startswith('#')]

    info(f"Species: {', '.join(species)}")

    bed_dir_path = Path(bed_dir)
    if not bed_dir_path.exists():
        error(f"BED directory not found: {bed_dir_path}")
        return {"success": False, "error": f"BED directory not found: {bed_dir_path}"}

    info(f"Loading BED files from: {bed_dir_path}")
    gene_to_species = load_bed_files(bed_dir_path, species)
    info(f"Loaded {len(gene_to_species)} gene-species mappings")

    network_file = Path(input_file)
    if not network_file.exists():
        error(f"Network file not found: {network_file}")
        return {"success": False, "error": f"Network file not found: {network_file}"}

    info(f"Loading network: {network_file}")
    edges, nodes, edge_weights = load_network(network_file)
    info(f"Loaded {len(edges)} edges, {len(nodes)} nodes")

    info(f"Clustering (method={method})...")
    try:
        if method == "cc":
            clusters = cluster_connected_components(edges, nodes)
        elif method == "louvain":
            clusters = cluster_louvain(edges, nodes)
        elif method == "infomap":
            clusters = cluster_infomap(edges, nodes)
        elif method == "label_prop":
            clusters = cluster_label_propagation(edges, nodes)
        elif method == "spectral":
            clusters = cluster_spectral(edges, nodes, k=k)
        else:
            raise ValueError(f"Unknown method: {method}")
    except ImportError as e:
        warning(str(e))
        info("Falling back to connected components...")
        clusters = cluster_connected_components(edges, nodes)
    info(f"Found {len(clusters)} raw clusters")

    info(f"Filtering: min_size={cluster_size}, min_species={min_species}")
    filtered_clusters, stats = filter_clusters(
        clusters, gene_to_species,
        min_cluster_size=cluster_size,
        min_species=min_species,
    )

    info(f"  Filtered by size < {cluster_size}: {stats.filtered_by_size}")
    info(f"  Filtered by species < {min_species}: {stats.filtered_by_species}")
    info(f"  Final clusters: {stats.final_clusters}")

    output_dir_path = Path(output_dir)
    output_dir_path.mkdir(parents=True, exist_ok=True)

    output_clusters = output_dir_path / "Filtered.clusters.tsv"
    output_summary = output_dir_path / "Filtered.cluster_summary.tsv"

    gene_to_cluster = write_clusters_tsv(filtered_clusters, output_clusters)
    write_cluster_summary_tsv(filtered_clusters, gene_to_species, output_summary)

    info(f"Exported: {output_clusters}")
    info(f"Exported: {output_summary}")

    output_synnet = output_dir_path / "Clusters.synnet.tsv"

    if gene_list:
        species_dir = species_list_path.parent
        gene_list_file = species_dir / gene_list
        if not gene_list_file.exists():
            error(f"Gene list file not found: {gene_list_file}")
            return {"success": False, "error": f"Gene list file not found: {gene_list_file}"}

        info(f"Loading gene list: {gene_list_file}")
        query_genes = load_gene_list(gene_list_file)
        info(f"Loaded {len(query_genes)} query genes")

        matched = sum(1 for g in query_genes if g in gene_to_cluster)
        info(f"Matched {matched} query genes in clusters")

        output_target = output_dir_path / "Filtered.target.clusters.tsv"
        write_target_clusters_tsv(gene_to_cluster, query_genes, edge_weights, output_target)
        info(f"Exported: {output_target}")

        write_synnet_tsv(gene_to_cluster, query_genes, edge_weights, output_synnet)
        info(f"Exported: {output_synnet}")
    else:
        write_synnet_tsv_all(filtered_clusters, edge_weights, output_synnet)
        info(f"Exported: {output_synnet}")

    info("Done!")
    return {
        "success": True,
        "output_clusters": str(output_clusters),
        "output_summary": str(output_summary),
        "output_synnet": str(output_synnet),
        "stats": {
            "total_clusters": stats.total_clusters,
            "filtered_by_size": stats.filtered_by_size,
            "filtered_by_species": stats.filtered_by_species,
            "final_clusters": stats.final_clusters,
        }
    }
