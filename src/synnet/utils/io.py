import gzip
import sys
from pathlib import Path
from typing import Union, TextIO, Optional, List

def open_gff(filepath: str, mode: str = 'rt') -> TextIO:
    if filepath == "-":
        return sys.stdin

    path = Path(filepath)

    if path.suffix == ".gz":
        return gzip.open(path, mode=mode)

    return open(path, mode=mode)

def check_files(*files: Union[str, Path], required: bool = True) -> bool:
    missing = [str(f) for f in files if not Path(f).exists()]

    if missing:
        msg = f"Missing files: {', '.join(missing)}"
        if required:
            from .logger import error
            error(msg)
        return False
    return True

def ensure_dir(dirpath: Union[str, Path]) -> Path:
    path = Path(dirpath)
    path.mkdir(parents=True, exist_ok=True)
    return path

def backup_file(filepath: str, suffix: str = ".bak") -> str:
    src = Path(filepath)
    if src.exists():
        dst = src.with_suffix(src.suffix + suffix)
        src.rename(dst)
        return str(dst)
    return str(src)

if __name__ == "__main__":
    import tempfile

    print("Testing utils/io.py")

    with tempfile.NamedTemporaryFile(mode='w', suffix='.gff', delete=False) as f:
        f.write("Chr1\t.\tgene\t100\t500\t.\t+\t.\tID=gene001\n")
        tmp_gff = f.name

    with open_gff(tmp_gff) as fin:
        line = fin.readline()
        assert "gene001" in line, "Failed to read GFF"
    print("open_gff: OK")

    assert check_files(tmp_gff) == True, "Existing file should pass"
    assert check_files("nonexistent.txt", required=False) == False, "Missing file should fail"
    print("check_files: OK")

    test_dir = Path(tempfile.mkdtemp()) / "nested" / "dir"
    result = ensure_dir(test_dir)
    assert result.exists(), "Directory should be created"
    print("ensure_dir: OK")

    # clean up
    Path(tmp_gff).unlink()
    print("All io utils tests passed!")