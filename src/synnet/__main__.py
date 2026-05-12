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
    from synnet.modules.gff2bed import gff3_to_bed, batch_gff2bed

    info("GFF3 to BED Converter")

    try:
        if args.species_list:
            results = batch_gff2bed(
                args.input, args.species_list, args.output_dir,
                feat_type=args.feat_type, id_key=args.id_key,
                min_length=args.min_length, verbose=args.verbose,
            )
            info(f"Done! Converted {len(results)} species")
        else:
            input_path = Path(args.input)
            if input_path.is_dir():
                error("Input is a directory but no -s species list provided")
                return 1

            result = gff3_to_bed(
                args.input, args.output,
                feat_type=args.feat_type, id_key=args.id_key,
                min_length=args.min_length, verbose=args.verbose,
            )
            info(f"Done! Output: {result}")
        return 0
    except Exception as e:
        error(f"Failed: {e}")
        return 1


@register("mcscan")
def cmd_mcscan(args):
    from synnet.modules.mcscan import run_chain_ortholog

    info("SynNet AutoMCScan")

    try:
        result = run_chain_ortholog(
            args.species_list,
            input_dir=args.input_dir,
            output_dir=args.output_dir,
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

    info("Build Synteny Network")

    try:
        result = run_network(
            args.species_list,
            input_dir=args.input_dir,
            output_dir=args.output_dir,
            no_lifted=args.no_lifted,
            min_score=args.min_score,
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
    from synnet.modules.cluster import run_cluster

    info("Cluster Synteny Network")

    try:
        result = run_cluster(
            args.input,
            args.species_list,
            args.bed_dir,
            output_dir=args.output_dir,
            method=args.method,
            k=args.k,
            cluster_size=args.cluster_size,
            min_species=args.min_species,
            gene_list=args.gene_list,
        )
        return 0 if result["success"] else 1
    except Exception as e:
        error(f"Failed: {e}")
        return 1


@register("viz")
def cmd_viz(args):
    from synnet.modules.viz import visualize_synnet

    info("Visualize Synteny Network")

    try:
        result = visualize_synnet(
            args.input,
            args.species_list,
            args.bed_dir,
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

    p = subparsers.add_parser("gff2bed", help="Convert GFF3 to BED format")
    p.add_argument("-i", "--input", required=True,
                   help="Input GFF3 file or directory [required]")
    p.add_argument("-o", "--output",
                   help="Output BED file (single file mode only)")
    p.add_argument("-s", "--species-list",
                   help="Species list file (one name per line)")
    p.add_argument("--output-dir",
                   help="Output directory for BED files (default: same as input)")
    p.add_argument("-t", "--feat-type", default="mRNA",
                   help="Feature type to extract (default: mRNA)")
    p.add_argument("-k", "--id-key", default="ID",
                   help="Attribute key for gene ID (default: ID)")
    p.add_argument("--min-length", type=int, default=0,
                   help="Minimum feature length (default: 0)")
    p.add_argument("-v", "--verbose", action="store_true",
                   help="Verbose output")
    p.set_defaults(func=cmd_gff2bed)

    p = subparsers.add_parser("mcscan", help="Chain-wise MCScan")
    p.add_argument("-s", "--species-list", required=True,
                   help="Species list file [required]")
    p.add_argument("-i", "--input-dir", default=".",
                   help="Directory containing sequence and annotation files (default: current dir)")
    p.add_argument("-o", "--output-dir", default="jcvi_output",
                   help="Output directory for jcvi results (default: jcvi_output)")
    p.add_argument("--cscore", type=float, default=0.7,
                   help="C-score cutoff (default: 0.7)")
    p.add_argument("--min-size", type=int, default=4,
                   help="Minimum anchors in a cluster (default: 4)")
    p.add_argument("--cpus", type=int, default=4,
                   help="CPU cores (default: 4)")
    p.add_argument("--no-intra", action="store_true",
                   help="Skip self synteny detection")
    p.add_argument("--dry-run", action="store_true",
                   help="Print commands without executing")
    p.add_argument("-v", "--verbose", action="store_true",
                   help="Verbose output")
    p.set_defaults(func=cmd_mcscan)

    p = subparsers.add_parser("network", help="Build synteny network")
    p.add_argument("-s", "--species-list", required=True,
                   help="Species list file [required]")
    p.add_argument("-i", "--input-dir", default="jcvi_output",
                   help="Input directory containing anchors files (default: jcvi_output)")
    p.add_argument("-o", "--output-dir", default="network_output",
                   help="Output directory for network files (default: network_output)")
    p.add_argument("--no-lifted", action="store_true",
                   help="Exclude lifted alignments (rows with 'L' suffix in score)")
    p.add_argument("--min-score", type=int, default=0,
                   help="Minimum score threshold (default: 0)")
    p.add_argument("-v", "--verbose", action="store_true",
                   help="Verbose output")
    p.set_defaults(func=cmd_network)

    p = subparsers.add_parser("cluster", help="Cluster synteny network and filter")
    p.add_argument("-i", "--input", required=True,
                   help="Input network TSV file (from network step) [required]")
    p.add_argument("-s", "--species-list", required=True,
                   help="Species list file [required]")
    p.add_argument("-b", "--bed-dir", required=True,
                   help="Directory containing .bed files (for gene-species mapping) [required]")
    p.add_argument("-o", "--output-dir", default="network_output",
                   help="Output directory (default: network_output)")
    p.add_argument("--method",
                   choices=["cc", "louvain", "infomap", "label_prop", "spectral"],
                   default="cc",
                   help="Clustering method (default: cc)")
    p.add_argument("--k", type=int, default=10,
                   help="Number of clusters for spectral method (default: 10)")
    p.add_argument("--cluster-size", type=int, default=2,
                   help="Minimum cluster size (default: 2)")
    p.add_argument("--min-species", type=int, default=1,
                   help="Minimum species count per cluster (default: 1)")
    p.add_argument("--gene-list",
                   help="Gene list file (same directory as species.lst) to filter clusters")
    p.add_argument("-v", "--verbose", action="store_true",
                   help="Verbose output")
    p.set_defaults(func=cmd_cluster)

    p = subparsers.add_parser("viz", help="Visualize synteny network (interactive HTML)")
    p.add_argument("-i", "--input", required=True,
                   help="Clusters.synnet.tsv file (from cluster command) [required]")
    p.add_argument("-s", "--species-list", required=True,
                   help="Species list file [required]")
    p.add_argument("-b", "--bed-dir", default=".",
                   help="Directory containing .bed files (default: current dir)")
    p.add_argument("-o", "--output-dir", default=None,
                   help="Output directory (default: same as input file)")
    p.add_argument("-v", "--verbose", action="store_true",
                   help="Verbose output")
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
