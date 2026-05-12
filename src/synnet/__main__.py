import argparse
import sys
from pathlib import Path

from synnet.utils.logger import setup_logger, info, warning, error
from synnet._version import __version__

COMMANDS = {}


def register(name):
    def decorator(func):
        COMMANDS[name] = func
        return func
    return decorator


@register("gff2bed")
def cmd_gff2bed(args):
    from synnet.modules.gff2bed import gff3_to_bed

    info(f"Converting: {args.input} -> {args.output or 'auto'}")

    try:
        output = gff3_to_bed(
            args.input,
            args.output,
            feat_type=args.feat_type,
            id_key=args.id_key,
            primary_only=args.primary_only,
            min_length=args.min_length,
            verbose=args.verbose,
        )
        info(f"Done: {output}")
        return 0
    except Exception as e:
        error(f"Failed: {e}")
        return 1


@register("mcscan")
def cmd_mcscan(args):
    from synnet.modules.mcscan import run_chain_ortholog

    info(f"MCScan: {args.species_list}")

    try:
        result = run_chain_ortholog(
            args.species_list,
            cscore=args.cscore,
            min_size=args.min_size,
            cpus=args.cpus,
            dry_run=args.dry_run,
            no_intra=args.no_intra,
        )
        if result["success"]:
            info("All comparisons completed successfully")
        else:
            error(f"Some comparisons failed: {result.get('error', 'see above')}")
        return 0 if result["success"] else 1
    except Exception as e:
        error(f"Failed: {e}")
        return 1


@register("network")
def cmd_network(args):
    from synnet.modules.network import run_network

    info(f"Network: {args.species_list}")

    try:
        result = run_network(
            args.species_list,
            work_dir=args.work_dir,
            use_lifted=not args.use_anchors,
            min_score=args.min_score,
            exclude_lifted=args.exclude_lifted,
            output_prefix=args.output_prefix,
            formats=args.formats,
        )
        if result["success"]:
            info("Network built successfully")
        else:
            error(f"Network build failed: {result.get('error', 'see above')}")
        return 0 if result["success"] else 1
    except Exception as e:
        error(f"Failed: {e}")
        return 1


@register("cluster")
def cmd_cluster(args):
    from synnet.modules.cluster import (
        run_cluster, export_clusters, export_cluster_summary, build_gene_species_map,
    )
    from synnet.modules.network import load_species_list

    info(f"Cluster: {args.input}")

    try:
        species = load_species_list(args.species_list)
        info(f"Loaded {len(species)} species: {' -> '.join(species)}")

        gene_species_map = build_gene_species_map(species, args.bed_dir)
        if gene_species_map:
            info(f"Built gene-species map: {len(gene_species_map)} genes")
        else:
            warning("No .bed files found, species-based filtering will be disabled")

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

        info(f"Loaded {len(edges)} edges, {len(nodes)} nodes")

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
            prefix = args.output_prefix
            export_clusters(result.clusters, Path(f"{prefix}.clusters.tsv"),
                            gene_species_map=gene_species_map)
            export_cluster_summary(result.clusters, Path(f"{prefix}.clusters.summary.tsv"),
                                   gene_species_map=gene_species_map)

        info(f"\nFiltering summary:")
        info(f"  Raw clusters: {result.n_before_filter}")
        info(f"  Filtered by size < {args.cluster_num}: {result.filtered_by_size}")
        info(f"  Filtered by species < {args.min_species}: {result.filtered_by_species}")
        info(f"  Filtered by ortholog: {result.filtered_by_ortholog}")
        info(f"  Final clusters: {len(result.clusters)}")

        info("Clustering completed!")
        return 0
    except Exception as e:
        error(f"Failed: {e}")
        return 1


@register("viz")
def cmd_viz(args):
    from synnet.modules.viz import visualize_synnet

    info(f"Viz: {args.input}")

    try:
        result = visualize_synnet(
            synnet_file=args.input,
            species_list_file=args.species_list,
            bed_dir=args.bed_dir,
            output_dir=args.output_dir,
        )
        if result["success"]:
            info("Visualization completed successfully")
            info(f"Output: {result.get('output', '')}")
        else:
            error(f"Visualization failed: {result.get('error', 'see above')}")
        return 0 if result["success"] else 1
    except Exception as e:
        error(f"Failed: {e}")
        return 1


def create_parser():
    parser = argparse.ArgumentParser(
        prog="synnet",
        description=f"SynNet v{__version__}: Synteny Network Builder\n\n"
                    "Pipeline: gff2bed -> mcscan -> network -> cluster -> viz",
        epilog="Use 'synnet <command> --help' for command-specific help.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("-V", "--version", action="version", version=f"%(prog)s {__version__}")

    subparsers = parser.add_subparsers(dest="command", metavar="<command>", help="Available commands")

    p = subparsers.add_parser("gff2bed", help="Step 1: Convert GFF3 to BED format")
    p.add_argument("-i", "--input", required=True, help="Input GFF3 file or directory")
    p.add_argument("-o", "--output", help="Output BED file or directory")
    p.add_argument("-t", "--feat-type", default="mRNA", help="Feature type (default: mRNA)")
    p.add_argument("-k", "--id-key", default="ID", help="Attribute key for gene ID (default: ID)")
    p.add_argument("--primary-only", action="store_true", help="Keep one entry per gene ID")
    p.add_argument("--min-length", type=int, default=0, help="Minimum feature length")
    p.add_argument("-v", "--verbose", action="store_true", help="Verbose output")
    p.set_defaults(func=cmd_gff2bed)

    p = subparsers.add_parser("mcscan", help="Step 2: Chain-wise MCScan (auto-detect sequence type)")
    p.add_argument("-s", "--species-list", required=True,
                   help="Species list file (.lst/.txt), one name per line")
    p.add_argument("--cscore", type=float, default=0.7,
                   help="C-score cutoff (default: 0.7)")
    p.add_argument("--min-size", type=int, default=4,
                   help="Minimum anchors in a cluster (default: 4)")
    p.add_argument("--cpus", type=int, default=4,
                   help="CPU cores (default: 4)")
    p.add_argument("--no-intra", action="store_true",
                   help="Skip intra-species (self) synteny detection")
    p.add_argument("--dry-run", action="store_true",
                   help="Print commands without executing")
    p.add_argument("-v", "--verbose", action="store_true", help="Verbose output")
    p.set_defaults(func=cmd_mcscan)

    p = subparsers.add_parser("network", help="Step 3: Build synteny network from .anchors files")
    p.add_argument("-s", "--species-list", required=True,
                   help="Species list file (chain order)")
    p.add_argument("-d", "--work-dir", default=".",
                   help="Directory containing .anchors files (default: current dir)")
    p.add_argument("--use-anchors", action="store_true",
                   help="Use .anchors instead of .lifted.anchors")
    p.add_argument("--exclude-lifted", action="store_true",
                   help="Exclude lifted edges (marked with 'L' suffix)")
    p.add_argument("--min-score", type=float, default=0,
                   help="Minimum alignment score cutoff (default: 0)")
    p.add_argument("-o", "--output-prefix", default="Final_Network",
                   help="Output file prefix (default: Final_Network)")
    p.add_argument("--formats", type=str, default="tsv",
                   help="Output formats: tsv,graphml,gexf (comma-separated, default: tsv)")
    p.add_argument("-v", "--verbose", action="store_true", help="Verbose output")
    p.set_defaults(func=cmd_network)

    p = subparsers.add_parser("cluster", help="Step 4: Cluster and filter synteny network")
    p.add_argument("-i", "--input", required=True,
                   help="Network TSV file (from 'network' command)")
    p.add_argument("-s", "--species-list", required=True,
                   help="Species list file (for species-based filtering)")
    p.add_argument("-d", "--bed-dir", default=".",
                   help="Directory containing .bed files for species mapping (default: current dir)")
    p.add_argument("-o", "--output-prefix", default="Filtered",
                   help="Output file prefix (default: Filtered)")
    p.add_argument("--cluster", type=str, default="cc",
                   choices=["cc", "mcl", "louvain"],
                   help="Clustering method: cc, mcl, louvain (default: cc)")
    p.add_argument("--mcl-inflation", type=float, default=2.0,
                   help="MCL inflation parameter (default: 2.0)")
    p.add_argument("--cluster-num", type=int, default=2,
                   help="Minimum cluster size to keep (default: 2)")
    p.add_argument("--min-species", type=int, default=1,
                   help="Minimum number of species in a cluster (default: 1)")
    p.add_argument("--ortholog-only", action="store_true",
                   help="Keep only 1-to-1 ortholog clusters")
    p.add_argument("-v", "--verbose", action="store_true", help="Verbose output")
    p.set_defaults(func=cmd_cluster)

    p = subparsers.add_parser("viz", help="Step 5: Visualize synteny network (interactive HTML)")
    p.add_argument("-i", "--input", required=True,
                   help="Clusters.synnet.tsv file (from 'cluster' command)")
    p.add_argument("-s", "--species-list", required=True,
                   help="Species list file (for color coding and species mapping)")
    p.add_argument("-d", "--bed-dir", default=".",
                   help="Directory containing .bed files (default: current dir)")
    p.add_argument("-o", "--output-dir", default=None,
                   help="Output directory (default: same as input file)")
    p.add_argument("-v", "--verbose", action="store_true", help="Verbose output")
    p.set_defaults(func=cmd_viz)

    return parser


def main():
    parser = create_parser()
    args = parser.parse_args()

    setup_logger(level="DEBUG" if args.verbose else "INFO")

    if args.command in COMMANDS:
        return COMMANDS[args.command](args)
    else:
        parser.print_help()
        return 1


if __name__ == "__main__":
    sys.exit(main())
