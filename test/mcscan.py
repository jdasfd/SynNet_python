"""
test/mcscan.py - Step 2: Chain-wise MCScan alignment (auto-detect sequence type)

Performs both intra-species (self) and inter-species (chain-wise) synteny detection.

Input: species list + .pep/.gff/.bed files in input directory
Output: .anchors/.lifted.anchors files in output directory

Usage:
    python mcscan.py -s species.lst -i seqs
    python mcscan.py -s species.lst -i seqs -o jcvi_output
    python mcscan.py -s species.lst -i seqs --cscore 0.9
    python mcscan.py -s species.lst -i seqs --min-size 5 --cpus 8
    python mcscan.py -s species.lst -i seqs --no-intra    # skip intra-species
    python mcscan.py -s species.lst -i seqs --dry-run -v
"""

import sys
import os
import argparse
import subprocess
import shutil
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


def load_species(list_file: str, input_dir: Path) -> List[SpeciesInfo]:
    species_list = []
    detected_types = set()

    with open(list_file, 'r') as f:
        names = [line.strip() for line in f if line.strip() and not line.startswith('#')]

    if len(names) < 2:
        raise ValueError(f"Species list must contain >= 2 names, got {len(names)}")

    for name in names:
        seq_file = None
        seq_type = None

        for suffix in CONFIG["prot_suffixes"] + CONFIG["cds_suffixes"]:
            candidate = input_dir / f"{name}{suffix}"
            if candidate.exists():
                seq_file = candidate
                seq_type = detect_seq_type(candidate)
                break

        if not seq_file or not seq_type:
            available = [f.name for f in input_dir.glob(f"{name}.*") if f.suffix]
            raise FileNotFoundError(
                f"Missing sequence file for '{name}' in {input_dir}.\n"
                f"Expected: {'/'.join(CONFIG['prot_suffixes'])} or {'/'.join(CONFIG['cds_suffixes'])}\n"
                f"Found: {available or 'none'}\n"
            )

        bed_file = input_dir / f"{name}.bed"
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

    return species_list


def generate_chain_pairs(species_list: List[SpeciesInfo]) -> List[SpeciesPair]:
    return [SpeciesPair(species_a=species_list[i], species_b=species_list[i + 1])
            for i in range(len(species_list) - 1)]


def _symlink_to_dir(files: List[Path], target_dir: Path):
    for f in files:
        link = target_dir / f.name
        if link.exists() or link.is_symlink():
            link.unlink()
        link.symlink_to(f.resolve())


def _copy_anchors_output(output_dir: Path, species_pairs: List[SpeciesPair]):
    anchors_suffixes = [".anchors", ".lifted.anchors", ".last", ".last.filtered"]
    collected = []
    for f in output_dir.iterdir():
        for suf in anchors_suffixes:
            if f.name.endswith(suf):
                collected.append(f)
                break
    return collected


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
        log_info(f"[DRY-RUN] Command: {' '.join(cmd)}")
        pair.status = "done"
        pair.anchors_file = Path(f"{pair.species_a.name}.{pair.species_b.name}.anchors")
        return pair

    try:
        result = subprocess.run(cmd, capture_output=True, text=True)

        if result.returncode != 0:
            pair.status = "failed"
            pair.error_msg = f"jcvi returned code {result.returncode}"
            log_error(f"JCVI failed: {pair.error_msg}")
            if result.stderr:
                log_error(f"stderr: {result.stderr[:500]}")
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


def run_jcvi_self(species: SpeciesInfo, *, cscore, min_size, cpus, dry_run) -> dict:
    log_info(f"Running intra-species: {species.name} vs {species.name}")

    cmd = [
        "python", "-m", "jcvi.compara.catalog", "ortholog",
        species.name,
        species.name,
        "--dbtype", species.seq_type,
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
        log_info(f"[DRY-RUN] Command: {' '.join(cmd)}")
        return {"status": "done", "n_anchors": 0}

    try:
        result = subprocess.run(cmd, capture_output=True, text=True)

        if result.returncode != 0:
            log_error(f"JCVI failed with code {result.returncode}")
            if result.stderr:
                log_error(f"stderr: {result.stderr[:500]}")
            return {"status": "failed", "n_anchors": 0}

        prefix = f"{species.name}.{species.name}"
        anchors = Path(f"{prefix}.anchors")

        if anchors.exists() and anchors.stat().st_size > 0:
            n_anchors = sum(1 for line in open(anchors) if not line.startswith('#'))
            log_info(f"{n_anchors} anchors (intra-species)")
            return {"status": "done", "n_anchors": n_anchors, "file": str(anchors)}
        else:
            log_warn("No anchors generated for intra-species")
            return {"status": "failed", "n_anchors": 0}

    except Exception as e:
        log_error(f"Exception: {e}")
        return {"status": "failed", "n_anchors": 0}


def run_chain_ortholog(
        species_list_file: str,
        *,
        input_dir: str = ".",
        output_dir: str = "jcvi_output",
        cscore: float = 0.7,
        min_size: int = 4,
        cpus: int = 4,
        dry_run: bool = False,
        no_intra: bool = False,
) -> dict:
    log_info("Step 2: SynNet AutoMCScan")
    log_info(f"Species list: {species_list_file}")
    log_info(f"Input dir: {input_dir}")
    log_info(f"Output dir: {output_dir}")

    in_dir = Path(input_dir).resolve()
    out_dir = Path(output_dir).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    try:
        species = load_species(species_list_file, in_dir)
    except (FileNotFoundError, ValueError) as e:
        log_error(f"Failed to load species: {e}")
        return {"success": False, "error": str(e)}

    all_input_files = []
    for sp in species:
        all_input_files.append(sp.seq_file)
        all_input_files.append(sp.bed_file)

    _symlink_to_dir(all_input_files, out_dir)
    log_info(f"Symlinked {len(all_input_files)} input files to {out_dir}")

    original_cwd = Path.cwd()
    os.chdir(out_dir)
    log_info(f"Working directory: {out_dir}")

    intra_results = []
    if not no_intra:
        log_info(f"\n=== Intra-species synteny ({len(species)} species) ===")
        for i, sp in enumerate(species, 1):
            log_info(f"\n[{i}/{len(species)}] Intra: {sp.name}")
            result = run_jcvi_self(sp, cscore=cscore, min_size=min_size, cpus=cpus, dry_run=dry_run)
            intra_results.append({"species": sp.name, **result})

    pairs = generate_chain_pairs(species)

    log_info(f"\n=== Inter-species synteny ({len(pairs)} pairs) ===")

    for i, pair in enumerate(pairs, 1):
        log_info(f"\n[{i}/{len(pairs)}]")
        run_jcvi_ortholog(pair, cscore=cscore, min_size=min_size, cpus=cpus, dry_run=dry_run)

    os.chdir(original_cwd)

    output_files = _copy_anchors_output(out_dir, pairs)

    intra_stats = {
        "total": len(intra_results),
        "done": sum(1 for r in intra_results if r.get("status") == "done"),
        "anchors": sum(r.get("n_anchors", 0) for r in intra_results),
    }

    inter_stats = {
        "total": len(pairs),
        "done": sum(1 for p in pairs if p.status == "done"),
        "failed": sum(1 for p in pairs if p.status == "failed"),
        "anchors": sum(p.n_anchors for p in pairs if p.status == "done"),
    }

    log_info(f"\n=== Results ===")
    if not no_intra:
        log_info(f"Intra-species: {intra_stats['done']}/{intra_stats['total']} done, {intra_stats['anchors']} anchors")
    log_info(f"Inter-species: {inter_stats['done']}/{inter_stats['total']} done, {inter_stats['anchors']} anchors")

    if inter_stats["done"] > 0 or intra_stats["done"] > 0:
        log_info(f"\nOutput files in {out_dir}:")
        if not no_intra:
            for r in intra_results:
                if r.get("status") == "done" and r.get("file"):
                    log_info(f"  {Path(r['file']).name} ({r['n_anchors']} anchors, intra)")
        for p in pairs:
            if p.status == "done" and p.anchors_file:
                log_info(f"  {p.anchors_file.name} ({p.n_anchors} anchors)")

    log_info("Done!")
    return {
        "success": inter_stats["failed"] == 0,
        "intra_stats": intra_stats,
        "inter_stats": inter_stats,
        "output_dir": str(out_dir)
    }


def main():
    parser = argparse.ArgumentParser(
        description="Step 2: Chain-wise MCScan alignment (auto-detect sequence type)",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python mcscan.py -s species.lst -i seqs                   # default output: jcvi_output/
  python mcscan.py -s species.lst -i seqs -o my_output      # custom output directory
  python mcscan.py -s species.lst -i seqs --cscore 0.9      # stricter C-score cutoff
  python mcscan.py -s species.lst -i seqs --min-size 5      # require larger anchor blocks
  python mcscan.py -s species.lst -i seqs --cpus 8          # use 8 CPU cores
  python mcscan.py -s species.lst -i seqs --no-intra        # skip intra-species synteny
  python mcscan.py -s species.lst -i seqs --dry-run -v      # preview commands only
        """,
    )

    parser.add_argument("-s", "--species-list", required=True,
                        help="Species list file (.lst/.txt), one name per line")
    parser.add_argument("-i", "--input-dir", default=".",
                        help="Directory containing .pep/.gff/.bed files (default: current dir)")
    parser.add_argument("-o", "--output-dir", default="jcvi_output",
                        help="Output directory for jcvi results (default: jcvi_output)")
    parser.add_argument("--cscore", type=float, default=0.7,
                        help="C-score cutoff for filtering anchors (default: 0.7)")
    parser.add_argument("--min-size", type=int, default=4, dest="min_size",
                        help="Minimum anchors in a cluster (default: 4)")
    parser.add_argument("--cpus", type=int, default=4,
                        help="CPU cores for LAST alignment (default: 4)")
    parser.add_argument("--no-intra", action="store_true",
                        help="Skip intra-species (self) synteny detection")
    parser.add_argument("--dry-run", action="store_true",
                        help="Print commands without executing")
    parser.add_argument("-v", "--verbose", action="store_true",
                        help="Verbose output")

    args = parser.parse_args()

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

    sys.exit(0 if result["success"] else 1)


if __name__ == "__main__":
    main()
