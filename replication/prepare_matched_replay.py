"""Restore exact historical bindings for the unchanged matched audit.

The released protocol/input exports contain provenance rebindings. This helper
does not edit them or the checker. It verifies a fixed recovery archive and
the narrowly declared export differences, then prepares an external audit.
"""
from __future__ import annotations

import argparse
import copy
import gzip
import hashlib
import io
import json
from pathlib import Path
import shlex
import sys
import tarfile

MANIFEST_SHA256 = "f46bc0a0efe62c62d30d32646e86b560ff2f6c69294d5ca6083f08e57794e5c8"
ARCHIVE_SHA256 = "ae388304467bc7665655960934b800675e6d6348d42cbb085e52f58a16b22101"
DEFERRED_RESOURCES_SHA256 = "9a562070de16acb4f43fb8b1fe821560b905fd9ad5cf42a6f8b4388c79fcd5af"
RESOURCE_CHECKER_SHA256 = "590e8ddaccc842cc00bc7d52bc0e22318cedd9f312c44a7e2ebf2ef06c05d3db"
PROTOCOL = "results/matched_deferred_chebyshev/matched_chebyshev_protocol.json"
INPUT = "results/deferred_cleanup/matched_B_inputs.json"
RESULT = "results/matched_deferred_chebyshev/matched_chebyshev_resources.json.gz"


class RecoveryError(ValueError):
    """The requested replay does not satisfy the fixed recovery contract."""


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RecoveryError(message)


def sha(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def checked_file(path: Path, expected: str) -> bytes:
    data = path.read_bytes()
    require(sha(data) == expected, f"Checksum mismatch: {path}")
    return data


def scientific_equal(original_protocol, exported_protocol, original_input,
                     exported_input, manifest) -> None:
    """Require equality after exactly the two declared provenance changes."""
    old_key = manifest["original_input_key"]
    new_key = manifest["exported_input_key"]
    require(Path(old_key).is_absolute(), "Original input key must be absolute")
    require(original_protocol["input_path"] == old_key, "Original input key changed")
    require(original_protocol["evidence_sha256"] == {
        old_key: manifest["archive_members"]["matched_B_inputs.json"]
    }, "Original evidence binding changed")
    comparison = copy.deepcopy(original_protocol)
    comparison["input_path"] = new_key
    comparison["evidence_sha256"] = {
        new_key: manifest["exported_sha256"][INPUT]
    }
    require(comparison == exported_protocol,
            "Protocol differs beyond declared input-path/evidence provenance")
    require(original_input["parent_evidence_sha256"]["deferred_cleanup_protocol.json"]
            == manifest["original_deferred_protocol_sha256"],
            "Original parent protocol binding changed")
    comparison = copy.deepcopy(original_input)
    comparison["parent_evidence_sha256"]["deferred_cleanup_protocol.json"] = (
        manifest["exported_deferred_protocol_sha256"])
    require(comparison == exported_input,
            "Input differs beyond declared parent protocol provenance")
    require(original_protocol["source_sha256"] == manifest["scientific_source_sha256"],
            "Scientific source manifest changed")


def archive_members(data: bytes, manifest) -> dict[str, bytes]:
    expected = manifest["archive_members"]
    with tarfile.open(fileobj=io.BytesIO(data), mode="r:gz") as archive:
        members = archive.getmembers()
        require(len(members) == len(expected), "Unexpected recovery archive members")
        require({member.name for member in members} == set(expected),
                "Unexpected recovery archive member names")
        recovered = {}
        for member in members:
            require(member.isfile() and member.size <= 2_000_000,
                    "Recovery archive contains an unsupported member")
            source = archive.extractfile(member)
            require(source is not None, "Recovery archive member is unreadable")
            body = source.read()
            require(sha(body) == expected[member.name],
                    "Recovered original checksum mismatch")
            recovered[member.name] = body
    return recovered


def validate_provenance_chain(repository: Path, manifest) -> int:
    """Invert every declared byte replacement and recover each original hash."""
    for item in manifest["provenance_chain"]:
        data = checked_file(repository / item["exported_repository_path"],
                            item["exported_sha256"])
        for replacement in reversed(item["exact_declared_byte_transform"]):
            source = replacement["from"].encode()
            target = replacement["to"].encode()
            require(source != target and target,
                    "Provenance inverse must be a nonempty distinct replacement")
            require(data.count(target) == replacement["occurrences"],
                    "Provenance inverse occurrence count mismatch")
            data = data.replace(target, source)
        require(sha(data) == item["original_sha256"],
                "Provenance inverse does not recover the original checksum")
    return len(manifest["provenance_chain"])


def prepare(repository: Path, workspace: Path, output: Path, archive: Path,
            manifest_path: Path, deferred_review: Path,
            python: str | None = None) -> dict:
    repository, workspace = repository.resolve(), workspace.resolve()
    output, deferred_review = output.resolve(), deferred_review.resolve()
    require(not output.exists(), "Output must be a new external directory")
    require(not output.is_relative_to(repository) and not output.is_relative_to(workspace),
            "Output must be outside the repository and reconstructed workspace")
    manifest = json.loads(checked_file(manifest_path, MANIFEST_SHA256))
    require(manifest["schema"] == "qbp-matched-original-bindings-recovery-v1",
            "Unsupported recovery schema")
    require(manifest["archive_sha256"] == ARCHIVE_SHA256, "Archive binding changed")
    provenance_links = validate_provenance_chain(repository, manifest)
    recovered = archive_members(checked_file(archive, ARCHIVE_SHA256), manifest)
    exported = {}
    for relative, expected in manifest["exported_sha256"].items():
        canonical = checked_file(repository / relative, expected)
        require(checked_file(workspace / relative, expected) == canonical,
                f"Reconstruction differs: {relative}")
        exported[relative] = canonical
    original_protocol = json.loads(recovered["matched_chebyshev_protocol.json"])
    original_input = json.loads(recovered["matched_B_inputs.json"])
    scientific_equal(original_protocol, json.loads(exported[PROTOCOL]),
                     original_input, json.loads(exported[INPUT]), manifest)
    result = json.loads(gzip.decompress(exported[RESULT]))
    require(result["status"] == "completed", "Authoritative result is incomplete")
    require(result["protocol_sha256"] == sha(recovered["matched_chebyshev_protocol.json"]),
            "Original protocol does not satisfy the result commitment")
    require(result["source_sha256"] == original_protocol["source_sha256"],
            "Result source premises differ")
    layout = json.loads((repository / "replication/source_layout.json").read_bytes())
    mapped = {item["original_path"]: item for item in layout["files"]}
    for relative, expected in manifest["scientific_source_sha256"].items():
        checked_file(workspace / relative, expected)
        require(relative in mapped and mapped[relative]["sha256"] == expected,
                f"Scientific source mapping changed: {relative}")
        checked_file(repository / mapped[relative]["current_path"], expected)
    checker = workspace / "check_matched_deferred_chebyshev.py"
    checked_file(checker, manifest["reconstructed_checker_sha256"])
    checked_file(workspace / "check_deferred_resources_independent.py",
                 RESOURCE_CHECKER_SHA256)
    resource_path = workspace / "results/deferred_cleanup/deferred_cleanup_resources.json.gz"
    checked_file(resource_path, DEFERRED_RESOURCES_SHA256)
    review = json.loads(deferred_review.read_bytes())
    require(review.get("status") == "passed"
            and review.get("resources_sha256") == DEFERRED_RESOURCES_SHA256
            and review.get("independent_checker_sha256") == RESOURCE_CHECKER_SHA256,
            "A matching passed independent resource review is required")
    require([review.get(key) for key in (
        "checked_architecture_profiles", "checked_qpe_plans",
        "checked_complete_variant_charges", "retained_rejections")]
        == [808, 1616, 3232, 12], "Resource review coverage changed")
    command = [python or sys.executable, str(checker), "--repo", str(workspace),
               "--results", str(output), "--deferred-resources", str(resource_path),
               "--deferred-review", str(deferred_review), "--input-evidence",
               str(output / "matched_B_inputs.json"), "--output",
               str(output / "matched_audit.json")]
    output.mkdir(parents=True)
    for name, body in recovered.items():
        (output / name).write_bytes(body)
    (output / "matched_chebyshev_resources.json.gz").write_bytes(exported[RESULT])
    record = {"schema": "qbp-prepared-original-matched-replay-v1",
              "status": "prepared_exact_original_bindings",
              "repository": str(repository), "workspace": str(workspace),
              "archive_sha256": ARCHIVE_SHA256, "manifest_sha256": MANIFEST_SHA256,
              "original_sha256": manifest["archive_members"],
              "exported_sha256": manifest["exported_sha256"],
              "scientific_sources_checked": len(manifest["scientific_source_sha256"]),
              "provenance_links_reconstructed": provenance_links,
              "deferred_review_sha256": sha(deferred_review.read_bytes()),
              "checker_sha256": manifest["reconstructed_checker_sha256"],
              "checker_argv": command, "checker_command": shlex.join(command),
              "scope": "Exact bindings restored; unchanged checker must still be run. "
                       "Preparation is not a successful audit or a new scientific result."}
    (output / "recovery_record.json").write_text(json.dumps(record, indent=2) + "\n")
    return record


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repository", type=Path, default=Path.cwd())
    parser.add_argument("--workspace", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--deferred-review", type=Path, required=True)
    parser.add_argument("--archive", type=Path)
    parser.add_argument("--manifest", type=Path)
    parser.add_argument("--python")
    args = parser.parse_args()
    recovery = args.repository / "replication/recovery"
    archive = args.archive or recovery / "matched_replay_bindings.tar.gz"
    manifest = args.manifest or recovery / "matched_replay_bindings_manifest.json"
    try:
        record = prepare(args.repository, args.workspace, args.output, archive,
                         manifest, args.deferred_review, args.python)
    except (OSError, ValueError, KeyError, tarfile.TarError) as error:
        parser.exit(1, f"Replay preparation failed: {error}\n")
    print(record["checker_command"])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
