"""Portable, lossless access to ordered shards of a compressed research archive.

Shards are consecutive *bytes of the original gzip file*, not independent gzip
streams. Their concatenation therefore preserves the original file exactly.
No uncompressed archive is written while checking or loading an artifact.

``load_json`` returns the whole decoded object, as Python's JSON loader normally
does. ``iter_top_level_array`` instead holds only one array item at a time, apart
from a small read buffer. Use it for large top-level arrays such as ``records``.
"""
from __future__ import annotations

import argparse
from contextlib import contextmanager
import gzip
import hashlib
import io
import json
import os
from pathlib import Path
import tempfile
from typing import Iterator

SCHEMA = "qbp-ordered-compressed-byte-shards-v1"
CHUNK = 1024 * 1024
DEFAULT_PART_BYTES = 32 * 1024 * 1024


def manifest_path(archive: str | Path) -> Path:
    """Return the sidecar name without requiring the canonical archive to exist."""
    p = Path(archive)
    return p if p.name.endswith(".manifest.json") else p.with_name(p.name + ".manifest.json")


def _manifest(archive: str | Path) -> tuple[Path, dict]:
    path = manifest_path(archive)
    data = json.loads(path.read_text(encoding="utf-8"))
    if data.get("schema") != SCHEMA or not data.get("parts"):
        raise ValueError("Unrecognised or empty shard manifest")
    if data["original_bytes"] != sum(part["bytes"] for part in data["parts"]):
        raise ValueError("Manifest part sizes do not sum to the original size")
    if [part["index"] for part in data["parts"]] != list(range(len(data["parts"]))):
        raise ValueError("Manifest parts are not in consecutive byte order")
    return path, data


def _part_path(manifest: Path, part: dict) -> Path:
    relative = Path(part["path"])
    if relative.is_absolute() or ".." in relative.parts:
        raise ValueError("Shard paths must remain below the manifest directory")
    path = (manifest.parent / relative).resolve()
    if not path.is_relative_to(manifest.parent.resolve()):
        raise ValueError("Shard resolves outside the manifest directory")
    return path


def _hash_file(path: Path) -> tuple[str, int]:
    h = hashlib.sha256()
    length = 0
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(CHUNK), b""):
            h.update(block)
            length += len(block)
    return h.hexdigest(), length


def verify_archive(archive: str | Path, *, compare_canonical: bool = True) -> dict:
    """Verify every shard and the SHA-256 of its ordered byte concatenation.

    The canonical file is optional. When present it is checked against the same
    original SHA-256 and byte length. Memory use is bounded by ``CHUNK``.
    """
    path, data = _manifest(archive)
    whole = hashlib.sha256()
    length = 0
    checked = []
    for part in data["parts"]:
        source = _part_path(path, part)
        local = hashlib.sha256()
        size = 0
        with source.open("rb") as stream:
            for block in iter(lambda: stream.read(CHUNK), b""):
                local.update(block)
                whole.update(block)
                size += len(block)
                length += len(block)
        if size != part["bytes"] or local.hexdigest() != part["sha256"]:
            raise ValueError(f"Shard checksum or length mismatch: {part['path']}")
        checked.append({"index": part["index"], "bytes": size, "sha256": local.hexdigest()})
    if length != data["original_bytes"] or whole.hexdigest() != data["original_sha256"]:
        raise ValueError("Ordered shard concatenation differs from the original archive")
    original = path.parent / data["original_filename"]
    canonical_checked = compare_canonical and original.exists()
    if canonical_checked:
        digest, size = _hash_file(original)
        if (digest, size) != (data["original_sha256"], data["original_bytes"]):
            raise ValueError("Canonical archive differs from the manifest")
    return {"status": "passed", "manifest": path.name, "parts": checked,
            "original_bytes": length, "ordered_concatenation_sha256": whole.hexdigest(),
            "canonical_file_checked": bool(canonical_checked)}


def archive_sha256(archive: str | Path) -> str:
    """Return the verified original compressed-byte hash with or without shards."""
    path = Path(archive)
    if manifest_path(path).exists():
        return verify_archive(path)["ordered_concatenation_sha256"]
    return _hash_file(path)[0]


def resolve_archive(root: str | Path, relative: str | Path) -> Path:
    """Return a verified gzip-compatible path, reconstructing only if necessary.

    A missing canonical file is reconstructed in the operating system's
    temporary directory, outside the repository. Its cache directory is named
    by the original compressed SHA-256. Existing cache files are rechecked.
    """
    archive = Path(root) / relative
    sidecar = manifest_path(archive)
    if not sidecar.exists():
        _hash_file(archive)
        return archive
    manifest, data = _manifest(archive)
    verify_archive(archive)
    canonical = manifest.parent / data["original_filename"]
    if canonical.exists():
        return canonical
    cache = Path(tempfile.gettempdir()) / "qbp-verified-archives" / data["original_sha256"]
    return reassemble_archive(archive, cache / data["original_filename"])


class _OrderedShardReader(io.RawIOBase):
    """Non-seeking byte stream across consecutive shard files."""
    def __init__(self, manifest: Path, data: dict):
        super().__init__()
        self._paths = iter(_part_path(manifest, part) for part in data["parts"])
        self._current = None

    def readable(self) -> bool:
        return True

    def readinto(self, target) -> int:
        if self.closed:
            raise ValueError("I/O operation on closed shard stream")
        view = memoryview(target)
        written = 0
        while written < len(view):
            if self._current is None:
                try:
                    self._current = next(self._paths).open("rb")
                except StopIteration:
                    break
            count = self._current.readinto(view[written:])
            if count:
                written += count
            else:
                self._current.close()
                self._current = None
        return written

    def close(self) -> None:
        if self._current is not None:
            self._current.close()
            self._current = None
        super().close()


@contextmanager
def open_archive(archive: str | Path, *, text: bool = True, verify: bool = True):
    """Open decoded gzip content from a canonical archive or its ordered shards.

    ``archive`` is the intended canonical path even when that file is absent.
    Passing the manifest path is also supported. If a manifest exists, checksum
    verification is performed by default before returning decoded content.
    """
    requested = Path(archive)
    sidecar = manifest_path(requested)
    if sidecar.exists():
        manifest, data = _manifest(requested)
        if verify:
            verify_archive(requested)
        canonical = manifest.parent / data["original_filename"]
    else:
        if requested.name.endswith(".manifest.json"):
            raise FileNotFoundError(requested)
        manifest, data, canonical = None, None, requested
    if canonical.exists():
        compressed = canonical.open("rb")
    elif manifest is not None:
        compressed = io.BufferedReader(_OrderedShardReader(manifest, data), buffer_size=CHUNK)
    else:
        raise FileNotFoundError(canonical)
    with compressed:
        with gzip.GzipFile(fileobj=compressed, mode="rb") as decoded:
            if text:
                with io.TextIOWrapper(decoded, encoding="utf-8") as stream:
                    yield stream
            else:
                yield decoded


def load_json(archive: str | Path, *, verify: bool = True):
    """Read a gzip JSON object without writing an expanded file to disk."""
    with open_archive(archive, verify=verify) as stream:
        return json.load(stream)


def iter_jsonl(archive: str | Path, *, verify: bool = True) -> Iterator[object]:
    """Read one JSON Lines record at a time from a gzip archive or its shards."""
    with open_archive(archive, verify=verify) as stream:
        for line in stream:
            if line.strip():
                yield json.loads(line)


class _JSONCursor:
    def __init__(self, stream):
        self.stream, self.buffer, self.position, self.eof = stream, "", 0, False
        self.decoder = json.JSONDecoder()

    def _read(self):
        if self.position:
            self.buffer = self.buffer[self.position:]
            self.position = 0
        block = self.stream.read(64 * 1024)
        self.buffer += block
        self.eof = not block

    def char(self):
        while True:
            while self.position < len(self.buffer) and self.buffer[self.position].isspace():
                self.position += 1
            if self.position < len(self.buffer):
                return self.buffer[self.position]
            if self.eof:
                raise ValueError("Unexpected end of JSON content")
            self._read()

    def expect(self, char):
        actual = self.char()
        if actual != char:
            raise ValueError(f"Expected {char!r}, found {actual!r}")
        self.position += 1

    def value(self):
        self.char()
        while True:
            try:
                result, end = self.decoder.raw_decode(self.buffer, self.position)
                # A numerical token may end at the current chunk boundary even
                # though its remaining digits have not arrived yet.
                if end == len(self.buffer) and not self.eof:
                    self._read()
                    continue
                self.position = end
                return result
            except json.JSONDecodeError:
                if self.eof:
                    raise
                self._read()


def iter_top_level_array(archive: str | Path, key: str, *, verify: bool = True) -> Iterator[object]:
    """Yield one item of a named top-level JSON array at a time.

    The array's items are decoded individually. Fields before it are decoded
    individually too, so choose an early field such as ``records`` for bounded
    memory use. Fields following the selected array are not parsed. The full
    compressed-byte integrity check still runs first when a manifest exists.
    """
    with open_archive(archive, verify=verify) as stream:
        cursor = _JSONCursor(stream)
        cursor.expect("{")
        while cursor.char() != "}":
            name = cursor.value()
            if not isinstance(name, str):
                raise ValueError("A JSON object key must be a string")
            cursor.expect(":")
            if name == key:
                cursor.expect("[")
                if cursor.char() != "]":
                    while True:
                        yield cursor.value()
                        if cursor.char() == "]":
                            break
                        cursor.expect(",")
                cursor.expect("]")
                return
            cursor.value()
            if cursor.char() != "}":
                cursor.expect(",")
        raise KeyError(f"Top-level array {key!r} was not found")


def reassemble_archive(archive: str | Path, output: str | Path) -> Path:
    """Write the original gzip bytes to a new file and check their SHA-256.

    Existing output is accepted only when it already equals the original.
    No decoded JSON or uncompressed research array is written.
    """
    path, data = _manifest(archive)
    verify_archive(archive)
    destination = Path(output)
    if destination.exists():
        if _hash_file(destination) != (data["original_sha256"], data["original_bytes"]):
            raise FileExistsError("Refusing to overwrite a different archive")
        return destination
    temporary = destination.with_name(destination.name + ".partial")
    if temporary.exists():
        raise FileExistsError(temporary)
    destination.parent.mkdir(parents=True, exist_ok=True)
    try:
        with temporary.open("xb") as out:
            for part in data["parts"]:
                with _part_path(path, part).open("rb") as source:
                    for block in iter(lambda: source.read(CHUNK), b""):
                        out.write(block)
        if _hash_file(temporary) != (data["original_sha256"], data["original_bytes"]):
            raise ValueError("Reassembled archive does not equal the original")
        temporary.replace(destination)
    except BaseException:
        temporary.unlink(missing_ok=True)
        raise
    return destination


def package_archive(archive: str | Path, max_part_bytes: int = DEFAULT_PART_BYTES) -> dict:
    """Create immutable ordered compressed-byte shards beside an existing gzip."""
    source = Path(archive)
    if max_part_bytes <= 0 or max_part_bytes > DEFAULT_PART_BYTES:
        raise ValueError("Part size must lie between one byte and 32 MiB")
    manifest = manifest_path(source)
    if manifest.exists():
        report = verify_archive(source)
        if max(part["bytes"] for part in report["parts"]) > max_part_bytes:
            raise ValueError("Existing package exceeds the requested part size")
        return report
    original_stat = source.stat()
    directory = source.with_name(source.name + ".parts")
    directory.mkdir(exist_ok=False)
    parts = []
    whole = hashlib.sha256()
    total = 0
    with source.open("rb") as stream:
        while True:
            first = stream.read(min(CHUNK, max_part_bytes))
            if not first:
                break
            name = directory / f"part-{len(parts):04d}.bin"
            local = hashlib.sha256()
            size = 0
            with name.open("xb") as out:
                block = first
                while block:
                    out.write(block)
                    local.update(block)
                    whole.update(block)
                    size += len(block)
                    total += len(block)
                    if size == max_part_bytes:
                        break
                    block = stream.read(min(CHUNK, max_part_bytes - size))
            parts.append({"index": len(parts), "path": str(name.relative_to(manifest.parent)),
                          "bytes": size, "sha256": local.hexdigest()})
    after = source.stat()
    if (after.st_size, after.st_mtime_ns, after.st_ino) != (original_stat.st_size, original_stat.st_mtime_ns, original_stat.st_ino):
        raise RuntimeError("Original archive changed during packaging")
    data = {"schema": SCHEMA, "representation": "ordered slices of the original compressed gzip bytes",
            "original_filename": source.name, "original_bytes": total,
            "original_sha256": whole.hexdigest(), "maximum_part_bytes": max_part_bytes,
            "parts": parts}
    temporary = manifest.with_name(manifest.name + ".partial")
    temporary.write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    temporary.replace(manifest)
    return verify_archive(source)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("command", choices=("package", "verify", "reassemble", "count"))
    parser.add_argument("archive", type=Path)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--key", default="records")
    args = parser.parse_args()
    if args.command == "package":
        print(json.dumps(package_archive(args.archive), indent=2))
    elif args.command == "verify":
        print(json.dumps(verify_archive(args.archive), indent=2))
    elif args.command == "reassemble":
        if args.output is None:
            parser.error("reassemble requires --output")
        print(reassemble_archive(args.archive, args.output))
    else:
        print(json.dumps({"array": args.key, "items": sum(1 for _ in iter_top_level_array(args.archive, args.key))}))


if __name__ == "__main__":
    main()
