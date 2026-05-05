import subprocess
from pathlib import Path
from typing import List, Optional, Literal
from dataclasses import dataclass

from synnet.utils.logger import get_logger, info, warning, error, success, debug

logger = get_logger(__name__)

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
        names = [line.strip() for line in f
                 if line.strip() and not line.startswith('#')]

    if len(names) < 2:
        raise ValueError(f"Species list must contain ≥2 names, got {len(names)}")

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

        debug(f"{name}: {seq_file.name} [{seq_type}], {bed_file.name}")

    if len(detected_types) > 1:
        raise ValueError(
            f"Mixed sequence types detected: {detected_types}."
        )
    info(f"Loaded {len(species_list)} species ({list(detected_types)[0]})")
    info(f"Working directory: {cwd}")

    return species_list


def generate_chain_pairs(species_list: List[SpeciesInfo]) -> List[SpeciesPair]:
    pairs = []
    for i in range(len(species_list) - 1):
        pairs.append(SpeciesPair(
            species_a=species_list[i],
            species_b=species_list[i + 1]
        ))
    return pairs


def run_jcvi_ortholog(
        pair: SpeciesPair,
        *,
        cscore: float,
        min_size: int,
        cpus: int,
        dry_run: bool,
) -> SpeciesPair:
    pair.status = "running"
    info(f"Running: {pair.species_a.name} vs {pair.species_b.name}")

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

    debug(f"CMD: {' '.join(cmd)}")

    if dry_run:
        info("[DRY-RUN] Skip execution")
        pair.status = "done"
        pair.anchors_file = Path(f"{pair.species_a.name}.{pair.species_b.name}.anchors")
        return pair

    try:
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
        )

        if result.returncode != 0:
            pair.status = "failed"
            pair.error_msg = f"jcvi returned code {result.returncode}"
            error(f"JCVI failed: {pair.error_msg}")
            return pair

        prefix = f"{pair.species_a.name}.{pair.species_b.name}"
        anchors = Path(f"{prefix}.anchors")

        if anchors.exists() and anchors.stat().st_size > 0:
            pair.anchors_file = anchors
            pair.n_anchors = sum(1 for line in open(anchors) if not line.startswith('#'))
            pair.status = "done"
            success(f"{pair.n_anchors} anchors")
        else:
            pair.status = "failed"
            pair.error_msg = "No anchors generated"
            warning(f"{pair.error_msg}")

        return pair

    except Exception as e:
        pair.status = "failed"
        pair.error_msg = str(e)
        error(f"Exception: {e}")
        return pair


def run_chain_ortholog(
        species_list_file: str,
        *,
        cscore: float = 0.7,
        min_size: int = 4,
        cpus: int = 4,
        dry_run: bool = False,
) -> dict:
    info("SynNet AutoMCScan")
    info(f"List: {species_list_file}")

    try:
        species = load_species_from_current_dir(species_list_file)
    except (FileNotFoundError, ValueError) as e:
        error(f"Failed to load species: {e}")
        return {"success": False, "error": str(e)}

    pairs = generate_chain_pairs(species)

    info(f"\nStarting {len(pairs)} comparisons...")

    for i, pair in enumerate(pairs, 1):
        info(f"\n[{i}/{len(pairs)}]")
        run_jcvi_ortholog(
            pair,
            cscore=cscore,
            min_size=min_size,
            cpus=cpus,
            dry_run=dry_run,
        )

    stats = {
        "total": len(pairs),
        "done": sum(1 for p in pairs if p.status == "done"),
        "failed": sum(1 for p in pairs if p.status == "failed"),
        "anchors": sum(p.n_anchors for p in pairs if p.status == "done"),
    }

    info(f"\nResults: {stats['done']}/{stats['total']} done, {stats['anchors']} anchors")

    if stats["done"] > 0:
        info(f"\nOutput files:")
        for p in pairs:
            if p.status == "done" and p.anchors_file:
                info(f"{p.anchors_file.name} ({p.n_anchors} anchors)")

    success("Completed!")

    return {
        "success": stats["failed"] == 0,
        "stats": stats,
    }
