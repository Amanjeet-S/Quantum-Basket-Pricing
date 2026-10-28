#!/usr/bin/env python3
"""Check retained bytes, paper exports, navigation and Python syntax without execution."""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
from pathlib import Path, PurePosixPath
import re
import subprocess
from urllib.parse import unquote, urlsplit


ROOT = Path(__file__).resolve().parents[1]
DIRECTORIES = {"paper", "src", "scripts", "tests", "configs", "data", "docs", "replication", "results", "cpp"}
ROOT_FILES = {"README.md", ".gitignore", "requirements.lock"}
PAPERS = {f"paper/{name}.{suffix}" for name in (
    "generalisation_proofs", "capped_gbm_companion", "constrained_hedging_companion")
    for suffix in ("tex", "pdf")}
GUIDES = ("README.md", "paper/README.md", "scripts/README.md", "replication/README.md",
          "data/README.md", "configs/README.md", "cpp/README.md", "docs/reproducing_paper.md")


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def safe_file(value: object) -> Path:
    if not isinstance(value, str) or not value or "\\" in value:
        raise ValueError(f"Invalid relative path: {value!r}")
    parts = value.split("/")
    if PurePosixPath(value).is_absolute() or any(x in ("", ".", "..") for x in parts):
        raise ValueError(f"Unsafe relative path: {value!r}")
    path = ROOT.joinpath(*parts)
    if path.is_symlink() or not path.resolve(strict=True).is_relative_to(ROOT) or not path.is_file():
        raise ValueError(f"Missing, escaping or symbolic file: {value}")
    return path


def prose(text: str) -> str:
    lines, fence = [], None
    for line in text.splitlines():
        marker = re.match(r"^\s{0,3}(`{3,}|~{3,})", line)
        if marker:
            if fence is None:
                fence = marker.group(1)
            elif marker.group(1)[0] == fence[0] and len(marker.group(1)) >= len(fence):
                fence = None
            continue
        if fence is None:
            lines.append(line)
    return re.sub(r"(`+).*?\1", "", "\n".join(lines), flags=re.S)


def check() -> dict:
    errors, counts = [], {"retained_files": 0, "paper_files": 0, "markdown_files": 0,
                          "local_links": 0, "python_files": 0, "tracked_files": None}
    try:
        spec = importlib.util.spec_from_file_location("qbp_layout", safe_file("replication/reproduce.py"))
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        _, _, files, _ = module.checked_files()
        counts["retained_files"] = len(files)
    except (OSError, ValueError, KeyError, TypeError, ImportError) as error:
        errors.append(f"Preservation: {error}")
    try:
        manifest = json.loads(safe_file("paper/manifest.json").read_text(encoding="utf-8"))
        records = manifest["files"]
        if manifest.get("schema") != "qbp-paper-manifest-v1" or not isinstance(records, list):
            raise ValueError("Unrecognised paper manifest schema")
        names = [record["path"] for record in records]
        if len(names) != len(PAPERS) or set(names) != PAPERS:
            raise ValueError("Paper manifest must contain the three TeX/PDF pairs exactly once")
        for record in records:
            path = safe_file(record["path"])
            if not re.fullmatch(r"[0-9a-f]{64}", record["sha256"]) or digest(path) != record["sha256"]:
                raise ValueError(f"Paper checksum mismatch: {record['path']}")
            if path.suffix == ".pdf" and not path.read_bytes().startswith(b"%PDF-"):
                raise ValueError(f"Invalid PDF header: {record['path']}")
            counts["paper_files"] += 1
    except (OSError, ValueError, KeyError, TypeError) as error:
        errors.append(f"Paper exports: {error}")
    for directory in sorted(DIRECTORIES):
        if not (ROOT / directory).is_dir():
            errors.append(f"Missing directory: {directory}")
    markdown = set(ROOT / name for name in GUIDES) | set((ROOT / "docs").glob("*.md"))
    for path in sorted(markdown):
        try:
            body = prose(path.read_text(encoding="utf-8"))
            counts["markdown_files"] += 1
            for raw in re.findall(r'\]\((<[^>]+>|[^)\s]+)(?:\s+"[^"]*")?\)', body):
                target = raw.removeprefix("<").removesuffix(">")
                if target.startswith("#") or urlsplit(target).scheme or target.startswith("//"):
                    continue
                target = unquote(target.split("#", 1)[0].split("?", 1)[0])
                resolved = (path.parent / target).resolve()
                counts["local_links"] += 1
                if not resolved.is_relative_to(ROOT) or not resolved.exists():
                    errors.append(f"Broken local link: {path.relative_to(ROOT)} -> {target}")
        except (OSError, ValueError) as error:
            errors.append(f"Navigation {path.relative_to(ROOT)}: {error}")
    for directory in ("src", "scripts", "tests", "replication"):
        for path in sorted((ROOT / directory).rglob("*.py")):
            if path.is_relative_to(ROOT / "replication/history"):
                continue
            try:
                compile(path.read_bytes(), str(path), "exec")
                counts["python_files"] += 1
            except (OSError, SyntaxError, ValueError) as error:
                errors.append(f"Python syntax {path.relative_to(ROOT)}: {error}")
    try:
        tracked = subprocess.run(["git", "-C", str(ROOT), "ls-files", "-z"], check=True,
                                 capture_output=True).stdout.decode().split("\0")
        tracked = [name for name in tracked if name]
        counts["tracked_files"] = len(tracked)
        for name in tracked:
            parts = PurePosixPath(name).parts
            if (len(parts) == 1 and name not in ROOT_FILES) or (len(parts) > 1 and parts[0] not in DIRECTORIES):
                errors.append(f"Unexpected tracked root entry: {name}")
    except (OSError, subprocess.CalledProcessError):
        pass  # Source archives without Git still receive all content checks.
    return {"status": "passed" if not errors else "failed", "counts": counts, "errors": errors}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, help="Optional new report file outside this repository")
    args = parser.parse_args()
    if args.output is not None:
        output = args.output.expanduser().resolve()
        if output.is_relative_to(ROOT) or output.exists():
            parser.error("The report must be a new file outside the repository")
    report = check()
    rendered = json.dumps(report, sort_keys=True)
    print(rendered)
    if args.output is not None:
        with output.open("x", encoding="utf-8") as stream:
            stream.write(rendered + "\n")
    return 0 if report["status"] == "passed" else 1


if __name__ == "__main__":
    raise SystemExit(main())
