"""
test/network.py - Build synteny network from .anchors / .lifted.anchors files

Usage:
    python network.py -s species.lst
    python network.py -s species.lst --use-anchors
    python network.py -s species.lst --min-score 100
    python network.py -s species.lst -o network --formats tsv,graphml
    python network.py -s species.lst --cluster mcl --mcl-inflation 2.0
    python network.py -s species.lst --min-species 2
    python network.py -s species.lst --ortholog-only
"""

import sys
import argparse
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent
SRC_DIR = PROJECT_ROOT / "src"
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

from synnet.utils.logger import setup_logger, info, warning, error, success
from synnet.modules.network import (
    load_species_list, build_network,
    export_tsv, export_graphml, export_gexf, export_stats,
)
from synnet.modules.cluster import (
    run_cluster, export_clusters, export_cluster_summary, build_gene_species_map,
)


def main():
    parser = argparse.ArgumentParser(
        description="Build synteny network from .anchors files",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python network.py -s species.lst                       # default: use .lifted.anchors
  python network.py -s species.lst --use-anchors         # use .anchors instead
  python network.py -s species.lst --min-score 100       # filter by alignment score
  python network.py -s species.lst --exclude-lifted      # exclude lifted anchor edges
  python network.py -s species.lst --cluster mcl         # use MCL clustering
  python network.py -s species.lst --cluster louvain     # use Louvain community detection
  python network.py -s species.lst --min-species 2       # keep clusters with >= 2 species
  python network.py -s species.lst --ortholog-only       # keep only 1-to-1 ortholog clusters
        """,
    )

    parser.add_argument("-s", "--species-list", required=True,
                        help="Species list file (chain order)")
    parser.add_argument("-d", "--work-dir", default=".",
                        help="Directory containing .anchors files (default: current dir)")
    parser.add_argument("--bed-dir", default=None,
                        help="Directory containing .bed files (default: same as work-dir)")
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
    parser.add_argument("--cluster", type=str, default="cc",
                        choices=["cc", "mcl", "louvain"],
                        help="Clustering method: cc (connected components), mcl, louvain (default: cc)")
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

    info("SynNet Network Builder")

    species = load_species_list(args.species_list)
    info(f"Loaded {len(species)} species: {' -> '.join(species)}")

    work_dir = Path(args.work_dir)
    if not work_dir.exists():
        error(f"Work directory not found: {work_dir}")
        sys.exit(1)

    use_lifted = not args.use_anchors
    info(f"Reading {'lifted ' if use_lifted else ''}anchors from {work_dir}")

    edges, nodes, stats = build_network(
        species, work_dir,
        use_lifted=use_lifted,
        min_score=args.min_score,
        exclude_lifted=args.exclude_lifted,
    )

    if not edges:
        error("No edges found. Check input files.")
        sys.exit(1)

    info(f"\nNetwork: {stats.total_nodes} nodes, {stats.total_edges} edges "
         f"({stats.lifted_edges} lifted)")
    if stats.filtered_by_score:
        info(f"Filtered by score < {args.min_score}: {stats.filtered_by_score} edges")

    formats = [f.strip() for f in args.formats.split(',')]
    prefix = args.output_prefix

    if "tsv" in formats:
        export_tsv(edges, Path(f"{prefix}.tsv"))
    if "graphml" in formats:
        export_graphml(edges, nodes, Path(f"{prefix}.graphml"))
    if "gexf" in formats:
        export_gexf(edges, nodes, Path(f"{prefix}.gexf"))

    export_stats(stats, Path(f"{prefix}.stats.txt"))

    simple_edges = [(e.source, e.target, e.score) for e in edges]

    gene_species_map = build_gene_species_map(species, args.bed_dir if args.bed_dir else args.work_dir)
    if gene_species_map:
        info(f"Built gene-species map: {len(gene_species_map)} genes from {len(set(gene_species_map.values()))} species")
    else:
        info("No .bed files found, species-based filtering will be disabled")

    result = run_cluster(
        simple_edges, nodes, species,
        method=args.cluster,
        mcl_inflation=args.mcl_inflation,
        min_cluster_size=args.cluster_num,
        min_species_count=args.min_species,
        ortholog_only=args.ortholog_only,
        gene_species_map=gene_species_map if gene_species_map else None,
    )

    if result.clusters:
        export_clusters(result.clusters, Path(f"{prefix}.clusters.tsv"),
                        gene_species_map=gene_species_map if gene_species_map else None)
        export_cluster_summary(result.clusters, Path(f"{prefix}.clusters.summary.tsv"),
                               gene_species_map=gene_species_map if gene_species_map else None)

    info(f"\nFiltering summary:")
    info(f"  Raw clusters: {result.n_before_filter}")
    info(f"  Filtered by size < {args.cluster_num}: {result.filtered_by_size}")
    info(f"  Filtered by species < {args.min_species}: {result.filtered_by_species}")
    info(f"  Filtered by ortholog: {result.filtered_by_ortholog}")
    info(f"  Final clusters: {len(result.clusters)}")

    success(f"\nDone! Output: {prefix}.*")


if __name__ == "__main__":
    main()
