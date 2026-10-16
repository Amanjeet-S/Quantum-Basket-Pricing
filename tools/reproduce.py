#!/usr/bin/env python3
"""Verify retained files or copy their original layout into an isolated workspace."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path, PurePosixPath
import re
import shutil
import sys


ROOT = Path(__file__).resolve().parents[1]
LAYOUT = ROOT / "tools" / "repository_layout.json"
SUMMARY = "_reproduction_workspace.json"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def relative_path(value: object) -> Path:
    if not isinstance(value, str) or not value or "\\" in value:
        raise ValueError(f"Invalid relative path: {value!r}")
    if PurePosixPath(value).is_absolute() or any(
            part in ("", ".", "..") for part in value.split("/")):
        raise ValueError(f"Unsafe relative path: {value!r}")
    if value == SUMMARY:
        raise ValueError(f"Reserved workspace record path: {value!r}")
    return Path(*PurePosixPath(value).parts)


def checked_files() -> tuple[dict, str, list[tuple[Path, Path, str]], int]:
    layout_hash = sha256(LAYOUT)
    manifest = json.loads(LAYOUT.read_text(encoding="utf-8"))
    if not isinstance(manifest, dict) or manifest.get("schema") != "qbp-repository-layout-v1":
        raise ValueError("Unrecognised repository layout schema.")
    commit = manifest.get("original_commit")
    if not isinstance(commit, str) or not re.fullmatch(r"(?:[0-9a-f]{40}|[0-9a-f]{64})", commit):
        raise ValueError("The original commit must be a hexadecimal Git identifier.")
    records = manifest.get("files")
    if not isinstance(records, list) or not records:
        raise ValueError("The repository layout must contain a non-empty file list.")
    checked, destinations, total_bytes = [], set(), 0
    for record in records:
        if not isinstance(record, dict):
            raise ValueError("Each mapped file must be a JSON object.")
        original = relative_path(record["original_path"])
        current = relative_path(record["current_path"])
        expected = record["sha256"]
        if not isinstance(expected, str) or not re.fullmatch(r"[0-9a-f]{64}", expected):
            raise ValueError(f"Invalid SHA256 for {current}.")
        if original in destinations:
            raise ValueError(f"Duplicate original path: {original}.")
        destinations.add(original)
        source = ROOT / current
        if source.is_symlink() or not source.resolve(strict=True).is_relative_to(ROOT):
            raise ValueError(f"Mapped source is a symlink or leaves the repository: {current}.")
        if not source.is_file():
            raise ValueError(f"Mapped source is not a regular file: {current}.")
        if sha256(source) != expected:
            raise ValueError(f"Source checksum mismatch: {current}.")
        total_bytes += source.stat().st_size
        checked.append((source, original, expected))
    for original in destinations:
        if any(parent in destinations for parent in original.parents):
            raise ValueError(f"Original file paths overlap: {original}.")
    if sha256(LAYOUT) != layout_hash:
        raise ValueError("The repository layout changed during verification.")
    return manifest, layout_hash, checked, total_bytes


def destination(value: Path) -> Path:
    expanded = value.expanduser()
    if expanded.is_symlink():
        raise ValueError("The output directory cannot be a symlink.")
    output = expanded.resolve()
    if output.is_relative_to(ROOT):
        raise ValueError("Choose an output directory outside the source repository.")
    if output.exists() and (not output.is_dir() or any(output.iterdir())):
        raise ValueError("The output directory must be new or empty.")
    return output


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    action = parser.add_mutually_exclusive_group(required=True)
    action.add_argument("--output", type=Path, help="New or empty directory outside this repository")
    action.add_argument("--verify-only", action="store_true", help="Check all retained hashes without copying")
    args = parser.parse_args()
    try:
        output = destination(args.output) if args.output is not None else None
        manifest, layout_hash, checked, total_bytes = checked_files()
        if output is None:
            print(f"Verified {len(checked)} retained files ({total_bytes} bytes); no files copied.")
            return 0
        # No destination is created until every mapped source has been checked.
        output = destination(output)
        output.mkdir(parents=True, exist_ok=True)
        for source, original, expected in checked:
            target = output / original
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(source, target)
            if sha256(target) != expected:
                raise ValueError(f"Copied checksum mismatch: {original}; incomplete workspace retained.")
        summary = {"schema": "qbp-reproduction-workspace-v1",
                   "original_commit": manifest["original_commit"],
                   "layout_sha256": layout_hash, "copied_files": len(checked),
                   "copied_bytes": total_bytes, "source_files_unchanged": True,
                   "scope": "Byte-exact original file layout; no network access or experiment execution."}
        (output / SUMMARY).write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
        print(f"Copied and verified {len(checked)} files at {output}.")
        return 0
    except (OSError, ValueError, KeyError, TypeError) as error:
        print(f"Reproduction failed: {error}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
