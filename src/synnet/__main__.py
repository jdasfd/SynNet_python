import argparse
import sys
from pathlib import Path

src_dir = Path(__file__).parent.parent
if str(src_dir) not in sys.path:
    sys.path.insert(0, str(src_dir))

from synnet.utils.logger import setup_logger, info, error
from synnet._version import __version__

COMMANDS = {}

def register(name):
    def decorator(func):
        COMMANDS[name] = func
        return func
    return decorator

# ==== command 1: gff2bed (jcvi.annotation.reformat) ====
@register("gff2bed")
def cmd_gff2bed(args):
    from synnet.modules.gff2bed import gff3_to_bed

    info(f"Converting: {args.input} → {args.output or 'auto'}")

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

# ==== command 2: automcscan (chain-wise JCVI ortholog) ====
@register("automcscan")
def cmd_automcscan(args):
    from synnet.modules.mcscan import run_chain_ortholog

    info(f"AutoMCScan: {args.species_list}")

    try:
        result = run_chain_ortholog(
            args.species_list,
            cscore=args.cscore,
            min_size=args.min_size,
            cpus=args.cpus,
            dry_run=args.dry_run,
        )
        if result["success"]:
            info("All comparisons completed successfully")
        else:
            error(f"Some comparisons failed: {result.get('error', 'see above')}")
        return 0 if result["success"] else 1
    except Exception as e:
        error(f"Failed: {e}")
        return 1

# ==== command 3: network (build synteny network) ====
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
            cluster_method=args.cluster,
            mcl_inflation=args.mcl_inflation,
            min_cluster_size=args.cluster_num,
        )
        if result["success"]:
            info("Network built successfully")
        else:
            error(f"Network build failed: {result.get('error', 'see above')}")
        return 0 if result["success"] else 1
    except Exception as e:
        error(f"Failed: {e}")
        return 1

def create_parser():
    parser = argparse.ArgumentParser(
        prog="python -m synnet",
        description=f"SynNet v{__version__}: Synteny Network Builder",
        epilog="Use 'python -m synnet <command> --help' for command-specific help.",
    )
    parser.add_argument("-V", "--version", action="version", version=f"%(prog)s {__version__}")

    subparsers = parser.add_subparsers(dest="command", metavar="<command>", help="Available commands")
    subparsers.required = True

    p = subparsers.add_parser("gff2bed", help="Convert GFF3 to BED format")
    p.add_argument("-i", "--input", required=True, help="Input GFF3 file")
    p.add_argument("-o", "--output", help="Output BED file (default: {input}.bed)")
    p.add_argument("-t", "--feat-type", default="mRNA", help="Feature type (default: mRNA)")
    p.add_argument("-k", "--id-key", default="ID", help="Attribute key for gene ID (default: ID)")
    p.add_argument("--primary-only", action="store_true", help="Keep one entry per gene ID")
    p.add_argument("--min-length", type=int, default=0, help="Minimum feature length")
    p.add_argument("-v", "--verbose", action="store_true", help="Verbose output")
    p.set_defaults(func=cmd_gff2bed)

    p = subparsers.add_parser("automcscan", help="Chain-wise MCScan (auto-detect sequence type)")
    p.add_argument("-s", "--species-list", required=True,
                   help="Species list file (.lst/.txt), one name per line in current dir")
    p.add_argument("--cscore", type=float, default=0.7,
                   help="C-score cutoff (default: 0.7)")
    p.add_argument("--min-size", type=int, default=4,
                   help="Minimum anchors in a cluster (default: 4)")
    p.add_argument("--cpus", type=int, default=4,
                   help="CPU cores (default: 4)")
    p.add_argument("--dry-run", action="store_true",
                   help="Print commands without executing")
    p.add_argument("-v", "--verbose", action="store_true", help="Verbose output")
    p.set_defaults(func=cmd_automcscan)

    p = subparsers.add_parser("network", help="Build synteny network from .anchors files")
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
    p.add_argument("--cluster", type=str, default="cc",
                   choices=["cc", "mcl", "louvain"],
                   help="Clustering method: cc, mcl, louvain (default: cc)")
    p.add_argument("--mcl-inflation", type=float, default=2.0,
                   help="MCL inflation parameter (default: 2.0)")
    p.add_argument("--cluster-num", type=int, default=2,
                   help="Minimum cluster size to keep (default: 2)")
    p.add_argument("-v", "--verbose", action="store_true", help="Verbose output")
    p.set_defaults(func=cmd_network)

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