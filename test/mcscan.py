"""
test/mcscan.py - Step 2: Chain-wise MCScan alignment (auto-detect sequence type)

Usage:
    python mcscan.py -s species.lst
    python mcscan.py -s species.lst --cscore 0.9
    python mcscan.py -s species.lst --min-size 5 --cpus 8
    python mcscan.py -s species.lst --dry-run -v
"""

import sys
import argparse
import subprocess
from pathlib import Path
from typing import List, Optional, Literal
from dataclasses import dataclass


def log_info(msg):
    print(f"[INFO] {msg}", file=sys.stderr)


def log_warn(msg):
    print(f"[WARN] {msg}", file=sys.stderr)


def log_error(msg):
    print(f"[ERROR] {msg}", file=sys.stderr)


CONFIG = {
    "align_soft": "last",
    "dist": 20,
    "no_strip_names": True,
    "no_dotplot": True,
    "prot_suffixes": [".pep", ".prot", ".faa"],
    "cds_suffixes": [".cds"],
}


@dataclass
class SpeciesInfo:
    name: str
    seq_type: Literal["prot", "cds"]
    seq_file: Path
    bed_file: Path


@dataclass
class SpeciesPair:
    species_a: SpeciesInfo
    species_b: SpeciesInfo
    anchors_file: Optional[Path] = None
    n_anchors: int = 0
    status: str = "pending"
    error_msg: str = ""


def detect_seq_type(filepath: Path) -> Optional[Literal["prot", "cds"]]:
    suffix = filepath.suffix.lower()
    if suffix in CONFIG["prot_suffixes"]:
        return "prot"
    elif suffix in CONFIG["cds_suffixes"]:
        return "cds"
    return None


def load_species_from_current_dir(list_file: str) -> List[SpeciesInfo]:
    cwd = Path.cwd()
    species_list = []
    detected_types = set()

    list_path = cwd / list_file
    if not list_path.exists():
        raise FileNotFoundError(f"Species list not found: {list_path}, -s required")

    with open(list_file, 'r') as f:
        names = [line.strip() for line in f if line.strip() and not line.startswith('#')]

    if len(names) < 2:
        raise ValueError(f"Species list must contain >= 2 names, got {len(names)}")

    for name in names:
        seq_file = None
        seq_type = None

        for suffix in CONFIG["prot_suffixes"] + CONFIG["cds_suffixes"]:
            candidate = cwd / f"{name}{suffix}"
            if candidate.exists():
                seq_file = candidate
                seq_type = detect_seq_type(candidate)
                break

        if not seq_file or not seq_type:
            available = [f.name for f in cwd.glob(f"{name}.*") if f.suffix]
            raise FileNotFoundError(
                f"Missing sequence file for '{name}' in current directory.\n"
                f"Expected: {'/'.join(CONFIG['prot_suffixes'])} or {'/'.join(CONFIG['cds_suffixes'])}\n"
                f"Found in {cwd}: {available or 'none'}\n"
            )

        bed_file = cwd / f"{name}.bed"
        if not bed_file.exists():
            raise FileNotFoundError(f"Missing BED file: {bed_file}")

        detected_types.add(seq_type)
        species_list.append(SpeciesInfo(
            name=name,
            seq_type=seq_type,
            seq_file=seq_file,
            bed_file=bed_file,
        ))

    if len(detected_types) > 1:
        raise ValueError(f"Mixed sequence types detected: {detected_types}.")

    log_info(f"Loaded {len(species_list)} species ({list(detected_types)[0]})")
    log_info(f"Working directory: {cwd}")

    return species_list


def generate_chain_pairs(species_list: List[SpeciesInfo]) -> List[SpeciesPair]:
    return [SpeciesPair(species_a=species_list[i], species_b=species_list[i + 1])
            for i in range(len(species_list) - 1)]


def run_jcvi_ortholog(pair: SpeciesPair, *, cscore, min_size, cpus, dry_run) -> SpeciesPair:
    pair.status = "running"
    log_info(f"Running: {pair.species_a.name} vs {pair.species_b.name}")

    cmd = [
        "python", "-m", "jcvi.compara.catalog", "ortholog",
        pair.species_a.name,
        pair.species_b.name,
        "--dbtype", pair.species_a.seq_type,
        "--cpus", str(cpus),
        "--cscore", str(cscore),
        "--min_size", str(min_size),
        "--dist", str(CONFIG["dist"]),
        "--align_soft", CONFIG["align_soft"],
    ]

    if CONFIG["no_strip_names"]:
        cmd.append("--no_strip_names")
    if CONFIG["no_dotplot"]:
        cmd.append("--no_dotplot")

    if dry_run:
        log_info("[DRY-RUN] Skip execution")
        pair.status = "done"
        pair.anchors_file = Path(f"{pair.species_a.name}.{pair.species_b.name}.anchors")
        return pair

    try:
        result = subprocess.run(cmd, capture_output=True, text=True)

        if result.returncode != 0:
            pair.status = "failed"
            pair.error_msg = f"jcvi returned code {result.returncode}"
            log_error(f"JCVI failed: {pair.error_msg}")
            return pair

        prefix = f"{pair.species_a.name}.{pair.species_b.name}"
        anchors = Path(f"{prefix}.anchors")

        if anchors.exists() and anchors.stat().st_size > 0:
            pair.anchors_file = anchors
            pair.n_anchors = sum(1 for line in open(anchors) if not line.startswith('#'))
            pair.status = "done"
            log_info(f"{pair.n_anchors} anchors")
        else:
            pair.status = "failed"
            pair.error_msg = "No anchors generated"
            log_warn(f"{pair.error_msg}")

        return pair

    except Exception as e:
        pair.status = "failed"
        pair.error_msg = str(e)
        log_error(f"Exception: {e}")
        return pair


def run_chain_ortholog(
        species_list_file: str,
        *,
        cscore: float = 0.7,
        min_size: int = 4,
        cpus: int = 4,
        dry_run: bool = False,
) -> dict:
    log_info("Step 2: SynNet AutoMCScan")
    log_info(f"List: {species_list_file}")

    try:
        species = load_species_from_current_dir(species_list_file)
    except (FileNotFoundError, ValueError) as e:
        log_error(f"Failed to load species: {e}")
        return {"success": False, "error": str(e)}

    pairs = generate_chain_pairs(species)

    log_info(f"\nStarting {len(pairs)} comparisons...")

    for i, pair in enumerate(pairs, 1):
        log_info(f"\n[{i}/{len(pairs)}]")
        run_jcvi_ortholog(pair, cscore=cscore, min_size=min_size, cpus=cpus, dry_run=dry_run)

    stats = {
        "total": len(pairs),
        "done": sum(1 for p in pairs if p.status == "done"),
        "failed": sum(1 for p in pairs if p.status == "failed"),
        "anchors": sum(p.n_anchors for p in pairs if p.status == "done"),
    }

    log_info(f"\nResults: {stats['done']}/{stats['total']} done, {stats['anchors']} anchors")

    if stats["done"] > 0:
        log_info(f"\nOutput files:")
        for p in pairs:
            if p.status == "done" and p.anchors_file:
                log_info(f"{p.anchors_file.name} ({p.n_anchors} anchors)")

    log_info("Done!")
    return {"success": stats["failed"] == 0, "stats": stats}


def main():
    parser = argparse.ArgumentParser(
        description="Step 2: Chain-wise MCScan alignment (auto-detect sequence type)",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python mcscan.py -s species.lst                    # default: cscore=0.7, min-size=4
  python mcscan.py -s species.lst --cscore 0.9       # stricter C-score cutoff
  python mcscan.py -s species.lst --min-size 5        # require larger anchor blocks
  python mcscan.py -s species.lst --cpus 8            # use 8 CPU cores
  python mcscan.py -s species.lst --dry-run -v        # preview commands only
        """,
    )

    parser.add_argument("-s", "--species-list", required=True,
                        help="Species list file (.lst/.txt), one name per line")
    parser.add_argument("--cscore", type=float, default=0.7,
                        help="C-score cutoff for filtering anchors (default: 0.7)")
    parser.add_argument("--min-size", type=int, default=4, dest="min_size",
                        help="Minimum anchors in a cluster (default: 4)")
    parser.add_argument("--cpus", type=int, default=4,
                        help="CPU cores for LAST alignment (default: 4)")
    parser.add_argument("--dry-run", action="store_true",
                        help="Print commands without executing")
    parser.add_argument("-v", "--verbose", action="store_true",
                        help="Verbose output")

    args = parser.parse_args()

    result = run_chain_ortholog(
        args.species_list,
        cscore=args.cscore,
        min_size=args.min_size,
        cpus=args.cpus,
        dry_run=args.dry_run,
    )

    sys.exit(0 if result["success"] else 1)


if __name__ == "__main__":
    main()
