"""Independent portable tests of lossless shard access and reconstruction."""
from pathlib import Path
import argparse
import gzip
import hashlib
import json
import os
import sys
import tempfile


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--out", type=Path)
    parser.add_argument("--stream-production-records", action="store_true")
    args = parser.parse_args()
    root = args.root.resolve()
    sys.path.insert(0, str(root))
    import artifact_archive as loader
    tested = []
    with tempfile.TemporaryDirectory(prefix="qbp-artifact-loader-") as tmp:
        stage = Path(tmp)
        value = {"schema": "loader-test", "large_header": "β" * 70000,
                 "records": [{"index": i, "payload": "é" * (i * 1000),
                              "negative": -10 ** (i + 1), "decimal": i / 10}
                             for i in range(10)], "tail": "preserved"}
        raw = (json.dumps(value, ensure_ascii=False) + "\n").encode("utf-8")
        original = gzip.compress(raw, mtime=0)
        archive = stage / "test.json.gz"
        archive.write_bytes(original)
        report = loader.package_archive(archive, max_part_bytes=128)
        assert report["original_bytes"] == len(original)
        assert report["ordered_concatenation_sha256"] == hashlib.sha256(original).hexdigest()
        assert loader.load_json(archive) == value
        assert list(loader.iter_top_level_array(archive, "records")) == value["records"]
        tested.append("canonical and streaming Unicode, negative integer and decimal decoding")
        archive.unlink()
        assert not archive.exists()
        assert loader.archive_sha256(archive) == hashlib.sha256(original).hexdigest()
        assert loader.load_json(archive) == value
        assert list(loader.iter_top_level_array(archive, "records")) == value["records"]
        with loader.open_archive(archive, text=False) as stream:
            assert stream.read() == raw
        tested.append("canonical-absent JSON and byte decoding from ordered compressed shards")
        output = stage / "restored.json.gz"
        assert loader.reassemble_archive(archive, output) == output
        assert output.read_bytes() == original
        assert loader.reassemble_archive(archive, output) == output
        output.write_bytes(b"different")
        try:
            loader.reassemble_archive(archive, output)
        except FileExistsError:
            pass
        else:
            raise AssertionError("A different destination must not be overwritten")
        tested.append("byteexact compressed reconstruction and refusal to overwrite different bytes")
        manifest = json.loads(loader.manifest_path(archive).read_text())
        first = stage / manifest["parts"][0]["path"]
        expected = first.read_bytes()
        first.write_bytes(bytes([expected[0] ^ 1]) + expected[1:])
        try:
            loader.verify_archive(archive)
        except ValueError:
            pass
        else:
            raise AssertionError("A damaged part must fail verification")
        first.write_bytes(expected)
        tested.append("single-byte shard corruption is detected")
        reversed_manifest = dict(manifest)
        reversed_manifest["parts"] = list(reversed(manifest["parts"]))
        sidecar = loader.manifest_path(archive)
        sidecar.write_text(json.dumps(reversed_manifest))
        try:
            loader.verify_archive(archive)
        except ValueError:
            pass
        else:
            raise AssertionError("Reordered parts must fail verification")
        sidecar.write_text(json.dumps(manifest))
        tested.append("shard reordering is rejected")
        line_archive = stage / "test.jsonl.gz"
        line_archive.write_bytes(gzip.compress(b'{"n":1}\n\n{"n":2}\n', mtime=0))
        loader.package_archive(line_archive, max_part_bytes=8)
        line_archive.unlink()
        assert list(loader.iter_jsonl(line_archive)) == [{"n": 1}, {"n": 2}]
        tested.append("canonical-absent JSON Lines streaming")
        # The production shards are hard-linked into a directory without the
        # canonical archive. This tests repository-style access without making
        # another large copy of either compressed or expanded data.
        production = root / "results/higherdim_precision_comparison.json.gz"
        production_manifest = json.loads(loader.manifest_path(production).read_text())
        isolated = stage / "production"
        isolated.mkdir()
        isolated_archive = isolated / production.name
        loader.manifest_path(isolated_archive).write_text(json.dumps(production_manifest))
        for part in production_manifest["parts"]:
            target = isolated / part["path"]
            target.parent.mkdir(parents=True, exist_ok=True)
            os.link(production.parent / part["path"], target)
        assert not isolated_archive.exists()
        production_report = loader.verify_archive(isolated_archive)
        assert not production_report["canonical_file_checked"]
        with loader.open_archive(isolated_archive) as stream:
            prefix = stream.read(100)
        assert prefix.startswith('{"schema":"qbp-higherdim-complete-precision-comparison-v1"')
        production_records = None
        if args.stream_production_records:
            production_records = sum(1 for _ in loader.iter_top_level_array(isolated_archive, "records"))
        tested.append("production archive loaded from verified hard-linked shards with canonical file absent")
    report = {"status": "passed", "tests": tested,
              "loader_sha256": hashlib.sha256((root / "artifact_archive.py").read_bytes()).hexdigest(),
              "audit_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
              "production_original_bytes": production_report["original_bytes"],
              "production_original_sha256": production_report["ordered_concatenation_sha256"],
              "production_streamed_records": production_records,
              "no_expanded_file_created": True}
    if args.out:
        args.out.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
