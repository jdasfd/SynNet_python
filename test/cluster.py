"""
test/cluster.py - Cluster and filter synteny network

Usage:
    python cluster.py -i Final_Network.tsv -s species.lst -d ../seqs
    python cluster.py -i Final_Network.tsv --cluster-num 4
    python cluster.py -i Final_Network.tsv --min-species 2
    python cluster.py -i Final_Network.tsv --ortholog-only
    python cluster.py -i Final_Network.tsv --cluster louvain --cluster-num 3
"""

import sys
import argparse
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent
SRC_DIR = PROJECT_ROOT / "src"
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

from synnet.utils.logger import setup_logger, info, warning, error, success
from synnet.modules.network import load_species_list
from synnet.modules.cluster import (
    run_cluster, export_clusters, export_cluster_summary, build_gene_species_map,
)


def main():
    parser = argparse.ArgumentParser(
        description="Cluster and filter synteny network",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python cluster.py -i Final_Network.tsv -s species.lst -d ../seqs   # cluster with species info
  python cluster.py -i Final_Network.tsv --cluster-num 4             # min cluster size = 4
  python cluster.py -i Final_Network.tsv --min-species 2             # require >= 2 species per cluster
  python cluster.py -i Final_Network.tsv --ortholog-only             # keep only 1-to-1 ortholog clusters
  python cluster.py -i Final_Network.tsv --cluster louvain           # use Louvain clustering
        """,
    )

    parser.add_argument("-i", "--input", required=True,
                        help="Network TSV file (from 'network' command)")
    parser.add_argument("-s", "--species-list", required=True,
                        help="Species list file (for species-based filtering)")
    parser.add_argument("-d", "--work-dir", default=".", required=True,
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
                        help="Keep only 1-to-1 ortholog clusters (each species has exactly 1 gene)")
    parser.add_argument("-v", "--verbose", action="store_true",
                        help="Verbose output")

    args = parser.parse_args()
    setup_logger(level="DEBUG" if args.verbose else "INFO")

    info("SynNet Cluster & Filter")

    species = []
    gene_species_map = None
    if args.species_list:
        species = load_species_list(args.species_list)
        info(f"Loaded {len(species)} species: {' -> '.join(species)}")
        gene_species_map = build_gene_species_map(species, args.work_dir)
        if gene_species_map:
            info(f"Built gene-species map: {len(gene_species_map)} genes from {len(set(gene_species_map.values()))} species")
        else:
            warning(f"No .bed files found in {args.work_dir}, species-based filtering will be disabled")

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

    info(f"Loaded {len(edges)} edges, {len(nodes)} nodes from {args.input}")

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

    info(f"\nFiltering summary:")
    info(f"  Raw clusters: {result.n_before_filter}")
    info(f"  Filtered by size < {args.cluster_num}: {result.filtered_by_size}")
    info(f"  Filtered by species < {args.min_species}: {result.filtered_by_species}")
    info(f"  Filtered by ortholog: {result.filtered_by_ortholog}")
    info(f"  Final clusters: {len(result.clusters)}")

    success(f"\nDone! Output: {args.output_prefix}.*")


if __name__ == "__main__":
    main()
