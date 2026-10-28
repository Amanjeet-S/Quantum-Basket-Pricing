"""Fail-closed checks for restoring exact matched-audit provenance.

Run from the organised repository root. These test the preparation interface;
the independent matched checker is a separate scientific verification.
"""
from __future__ import annotations

import importlib.util
import copy
import json
import os
from pathlib import Path
import shutil
import tempfile
import unittest

HERE = Path(__file__).resolve().parent
REPOSITORY = Path(os.environ.get("QBP_RECOVERY_REPOSITORY", Path.cwd())).resolve()
HELPER = Path(os.environ.get("QBP_RECOVERY_HELPER", HERE / "prepare_matched_replay.py"))
if not HELPER.exists():
    HELPER = REPOSITORY / "replication/prepare_matched_replay.py"
spec = importlib.util.spec_from_file_location("matched_recovery_under_test", HELPER)
recovery = importlib.util.module_from_spec(spec)
spec.loader.exec_module(recovery)
RECOVERY = HERE if (HERE / "matched_replay_bindings.tar.gz").exists() else (
    REPOSITORY / "replication/recovery")


class MatchedRecoveryTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.fixture_temp = tempfile.TemporaryDirectory(prefix="qbp-recovery-fixture-")
        cls.fixture = Path(cls.fixture_temp.name)
        cls.manifest = json.loads((RECOVERY / "matched_replay_bindings_manifest.json").read_bytes())
        layout = json.loads((REPOSITORY / "replication/source_layout.json").read_bytes())
        mapped = {item["original_path"]: item["current_path"] for item in layout["files"]}
        for relative in cls.manifest["scientific_source_sha256"]:
            target = cls.fixture / relative
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(REPOSITORY / mapped[relative], target)
        for relative in cls.manifest["exported_sha256"]:
            target = cls.fixture / relative
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(REPOSITORY / relative, target)
        for name in ("check_matched_deferred_chebyshev.py", "check_deferred_resources_independent.py"):
            shutil.copyfile(REPOSITORY / mapped[name], cls.fixture / name)
        relative = "results/deferred_cleanup/deferred_cleanup_resources.json.gz"
        shutil.copyfile(REPOSITORY / relative, cls.fixture / relative)

    @classmethod
    def tearDownClass(cls):
        cls.fixture_temp.cleanup()

    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix="qbp-recovery-case-")
        self.root = Path(self.temp.name)
        self.workspace = self.root / "workspace"
        shutil.copytree(self.fixture, self.workspace)
        self.output = self.root / "prepared"
        self.archive = self.root / "bindings.tar.gz"
        self.manifest_path = self.root / "manifest.json"
        shutil.copyfile(RECOVERY / "matched_replay_bindings.tar.gz", self.archive)
        shutil.copyfile(RECOVERY / "matched_replay_bindings_manifest.json", self.manifest_path)
        self.review = self.root / "resource_review.json"
        # Interface fixture only. A separate actual scientific replay supplies
        # the real passed review; this unit test must not relabel an older
        # retained reviewer version as the currently required checker.
        self.review.write_text(json.dumps({
            "unit_test_interface_fixture_only": True,
            "status": "passed",
            "resources_sha256": recovery.DEFERRED_RESOURCES_SHA256,
            "independent_checker_sha256": recovery.RESOURCE_CHECKER_SHA256,
            "checked_architecture_profiles": 808,
            "checked_qpe_plans": 1616,
            "checked_complete_variant_charges": 3232,
            "retained_rejections": 12,
        }))

    def tearDown(self):
        self.temp.cleanup()

    def prepare(self):
        return recovery.prepare(REPOSITORY, self.workspace, self.output, self.archive,
                                self.manifest_path, self.review)

    def reject(self, pattern):
        with self.assertRaisesRegex(recovery.RecoveryError, pattern):
            self.prepare()
        self.assertFalse(self.output.exists())

    def test_exact_originals_and_unchanged_checker_command(self):
        record = self.prepare()
        self.assertEqual(record["scientific_sources_checked"], 9)
        self.assertEqual(record["provenance_links_reconstructed"], 8)
        for name, expected in self.manifest["archive_members"].items():
            self.assertEqual(recovery.sha((self.output / name).read_bytes()), expected)
        command = record["checker_argv"]
        self.assertEqual(command[1], str((self.workspace / "check_matched_deferred_chebyshev.py").resolve()))
        self.assertIn("--input-evidence", command)
        self.assertEqual(len(list(self.output.iterdir())), 4)
        self.assertFalse((self.output / "matched_audit.json").exists())

    def test_changed_exported_protocol_is_rejected(self):
        path = self.workspace / recovery.PROTOCOL
        path.write_bytes(path.read_bytes() + b" ")
        self.reject("Checksum mismatch")

    def test_changed_exported_input_is_rejected(self):
        path = self.workspace / recovery.INPUT
        path.write_bytes(path.read_bytes() + b" ")
        self.reject("Checksum mismatch")

    def test_changed_scientific_source_is_rejected(self):
        path = self.workspace / "qpe_estimator.py"
        path.write_bytes(path.read_bytes() + b"\n# altered fixture\n")
        self.reject("Checksum mismatch")

    def test_changed_archive_is_rejected(self):
        self.archive.write_bytes(self.archive.read_bytes() + b"altered")
        self.reject("Checksum mismatch")

    def test_changed_manifest_is_rejected(self):
        self.manifest_path.write_bytes(self.manifest_path.read_bytes() + b" ")
        self.reject("Checksum mismatch")

    def test_scientific_case_change_is_not_an_allowed_provenance_delta(self):
        files = recovery.archive_members(self.archive.read_bytes(), self.manifest)
        original_protocol = json.loads(files["matched_chebyshev_protocol.json"])
        original_input = json.loads(files["matched_B_inputs.json"])
        original_protocol["cases"][0]["model_id"] = "unregistered-model"
        with self.assertRaisesRegex(recovery.RecoveryError, "beyond declared"):
            recovery.scientific_equal(
                original_protocol, json.loads((self.workspace / recovery.PROTOCOL).read_bytes()),
                original_input, json.loads((self.workspace / recovery.INPUT).read_bytes()),
                self.manifest)

    def test_failed_independent_resource_review_is_rejected(self):
        review = json.loads(self.review.read_bytes())
        review["status"] = "failed"
        self.review.write_text(json.dumps(review))
        self.reject("matching passed independent resource review")

    def chain_fixture(self):
        manifest = copy.deepcopy(self.manifest)
        manifest["provenance_chain"] = manifest["provenance_chain"][:1]
        root = self.root / "chain-repository"
        relative = manifest["provenance_chain"][0]["exported_repository_path"]
        target = root / relative
        target.parent.mkdir(parents=True)
        shutil.copyfile(REPOSITORY / relative, target)
        return root, manifest, target

    def test_changed_provenance_chain_export_is_rejected(self):
        root, manifest, target = self.chain_fixture()
        target.write_bytes(target.read_bytes() + b"\n# altered chain export\n")
        with self.assertRaisesRegex(recovery.RecoveryError, "Checksum mismatch"):
            recovery.validate_provenance_chain(root, manifest)

    def test_provenance_inverse_requires_exact_occurrence_counts(self):
        root, manifest, _ = self.chain_fixture()
        manifest["provenance_chain"][0]["exact_declared_byte_transform"][-1]["occurrences"] += 1
        with self.assertRaisesRegex(recovery.RecoveryError, "occurrence count mismatch"):
            recovery.validate_provenance_chain(root, manifest)

    def test_existing_output_is_not_overwritten(self):
        self.output.mkdir()
        marker = self.output / "keep.txt"
        marker.write_text("preserve")
        with self.assertRaisesRegex(recovery.RecoveryError, "new external directory"):
            self.prepare()
        self.assertEqual(marker.read_text(), "preserve")

    def test_output_within_frozen_workspace_is_rejected(self):
        self.output = self.workspace / "must-not-be-created"
        self.reject("outside the repository and reconstructed workspace")


if __name__ == "__main__":
    unittest.main()
