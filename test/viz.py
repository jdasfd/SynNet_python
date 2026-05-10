"""
test/viz.py - Visualize synteny network

Usage:
    python viz.py -i Final_Network.tsv -s species.lst -d ../seqs
    python viz.py -i Final_Network.tsv --plot-type static --layout circular
    python viz.py -i Final_Network.tsv --plot-type interactive
    python viz.py -i Final_Network.tsv --top-k 100
    python viz.py -i Final_Network.tsv --output-dir viz_output
"""

import sys
import argparse
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent
SRC_DIR = PROJECT_ROOT / "src"
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

from synnet.utils.logger import setup_logger, info, warning, error, success
from synnet.modules.viz import run_viz


def main():
    parser = argparse.ArgumentParser(
        description="Visualize synteny network",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python viz.py -i Final_Network.tsv -s species.lst -d ../seqs  # generate all plots
  python viz.py -i Final_Network.tsv --plot-type static         # static PNG only
  python viz.py -i Final_Network.tsv --plot-type interactive    # interactive HTML only
  python viz.py -i Final_Network.tsv --layout circular          # circular layout
  python viz.py -i Final_Network.tsv --top-k 100                # show top 100 nodes
  python viz.py -i Final_Network.tsv --output-dir viz_output    # output to directory
        """,
    )

    parser.add_argument("-i", "--input", required=True,
                        help="Network TSV file (from 'network' command)")
    parser.add_argument("-s", "--species-list",
                        help="Species list file (for color coding by species)")
    parser.add_argument("-d", "--work-dir", default=".",
                        help="Directory containing .bed files (for species mapping, default: current dir)")
    parser.add_argument("--output-dir", default=".",
                        help="Output directory for plots (default: current dir)")
    parser.add_argument("--plot-type", type=str, default="all",
                        choices=["static", "interactive", "all"],
                        help="Plot type: static (PNG), interactive (HTML), all (default: all)")
    parser.add_argument("--layout", type=str, default="spring",
                        choices=["spring", "circular", "kamada_kawai"],
                        help="Layout algorithm for static plot (default: spring)")
    parser.add_argument("--top-k", type=int, default=None,
                        help="Show only top K nodes by degree (default: show all)")
    parser.add_argument("--interactive", action="store_true",
                        help="Force interactive HTML output (requires pyvis)")
    parser.add_argument("-v", "--verbose", action="store_true",
                        help="Verbose output")

    args = parser.parse_args()
    setup_logger(level="DEBUG" if args.verbose else "INFO")

    info("SynNet Visualization")

    try:
        result = run_viz(
            args.input,
            species_list_file=args.species_list,
            work_dir=args.work_dir,
            output_dir=args.output_dir,
            plot_type=args.plot_type,
            layout=args.layout,
            top_k=args.top_k,
            interactive=args.interactive,
        )
        if result["success"]:
            for fmt, path in result.get("outputs", {}).items():
                info(f"  {fmt}: {path}")
            success("\nDone!")
        else:
            error("Visualization failed")
            sys.exit(1)
    except ImportError as e:
        error(f"Missing dependency: {e}")
        info("Install required packages: pip install matplotlib networkx pyvis")
        sys.exit(1)


if __name__ == "__main__":
    main()
