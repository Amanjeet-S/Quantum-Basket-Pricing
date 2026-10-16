"""Safety and byte-preservation checks for isolated workspace reconstruction."""

import contextlib
import copy
import hashlib
import io
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import reproduce


class ReproductionTests(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory(prefix="qbp-reproduction-test-")
        self.addCleanup(temporary.cleanup)
        self.base = Path(temporary.name).resolve()
        self.root = self.base / "repository"
        (self.root / "tools").mkdir(parents=True)
        (self.root / "src").mkdir()
        self.layout = self.root / "tools/repository_layout.json"
        self.payloads = {"src/model.py": b"print('retained')\n", "src/evidence.bin": b"\x00\xff\n"}
        for name, payload in self.payloads.items():
            (self.root / name).write_bytes(payload)
        self.manifest = {"schema": "qbp-repository-layout-v1", "original_commit": "a" * 40,
                         "files": [{"original_path": name.removeprefix("src/"),
                                    "current_path": name,
                                    "sha256": hashlib.sha256(payload).hexdigest()}
                                   for name, payload in self.payloads.items()]}
        self.original = copy.deepcopy(self.manifest)
        self.output = self.base / "workspace"
        self.write_layout()
        for name, value in (("ROOT", self.root), ("LAYOUT", self.layout)):
            active = patch.object(reproduce, name, value)
            active.start()
            self.addCleanup(active.stop)

    def write_layout(self):
        self.layout.write_text(json.dumps(self.manifest), encoding="utf-8")

    def run_tool(self, *arguments):
        stdout, stderr = io.StringIO(), io.StringIO()
        with patch("sys.argv", ["reproduce.py", *map(str, arguments)]), \
                contextlib.redirect_stdout(stdout), contextlib.redirect_stderr(stderr):
            status = reproduce.main()
        return status, stdout.getvalue(), stderr.getvalue()

    def test_isolated_copy_preserves_bytes_and_records_summary(self):
        self.output.mkdir()  # An existing empty output is also supported.
        status, _, error = self.run_tool("--output", self.output)
        self.assertEqual(status, 0, error)
        for current, payload in self.payloads.items():
            target = self.output / current.removeprefix("src/")
            self.assertEqual(target.read_bytes(), payload)
            self.assertNotEqual(target.stat().st_ino, (self.root / current).stat().st_ino)
            target.write_bytes(b"independent workspace change")
            self.assertEqual((self.root / current).read_bytes(), payload)
        summary = json.loads((self.output / reproduce.SUMMARY).read_text())
        self.assertEqual(summary["original_commit"], self.manifest["original_commit"])
        self.assertEqual(summary["layout_sha256"], reproduce.sha256(self.layout))
        self.assertEqual(summary["copied_files"], 2)
        self.assertEqual(summary["copied_bytes"], sum(map(len, self.payloads.values())))

    def test_verify_only_creates_no_files(self):
        before = {p.relative_to(self.base) for p in self.base.rglob("*")}
        status, message, error = self.run_tool("--verify-only")
        self.assertEqual(status, 0, error)
        self.assertIn("no files copied", message)
        self.assertEqual(before, {p.relative_to(self.base) for p in self.base.rglob("*")})

    def test_bad_later_hash_fails_before_output_creation(self):
        self.manifest["files"][-1]["sha256"] = "0" * 64
        self.write_layout()
        status, _, error = self.run_tool("--output", self.output)
        self.assertEqual(status, 1)
        self.assertIn("checksum mismatch", error)
        self.assertFalse(self.output.exists())

    def test_unsafe_paths_are_rejected(self):
        for field in ("original_path", "current_path"):
            for value in ("../escape", "/absolute", "x/../escape", "x//y", "x\\y"):
                with self.subTest(field=field, value=value):
                    self.manifest = copy.deepcopy(self.original)
                    self.manifest["files"][0][field] = value
                    self.write_layout()
                    self.assertEqual(self.run_tool("--output", self.output)[0], 1)
                    self.assertFalse(self.output.exists())

    def test_duplicate_original_paths_are_rejected(self):
        self.manifest["files"][1]["original_path"] = self.manifest["files"][0]["original_path"]
        self.write_layout()
        self.assertEqual(self.run_tool("--output", self.output)[0], 1)
        self.assertFalse(self.output.exists())

    def test_overlapping_original_file_paths_are_rejected(self):
        self.manifest["files"][1]["original_path"] = "model.py/child"
        self.write_layout()
        self.assertEqual(self.run_tool("--output", self.output)[0], 1)
        self.assertFalse(self.output.exists())

    def test_symlink_source_is_rejected(self):
        link = self.root / "src/link.py"
        link.symlink_to(self.root / "src/model.py")
        self.manifest["files"][0]["current_path"] = "src/link.py"
        self.write_layout()
        self.assertEqual(self.run_tool("--output", self.output)[0], 1)
        self.assertFalse(self.output.exists())

    def test_nonempty_and_source_contained_outputs_are_rejected(self):
        self.output.mkdir()
        sentinel = self.output / "keep.txt"
        sentinel.write_text("preserve")
        for output in (self.output, self.root, self.root / "new-workspace"):
            with self.subTest(output=output):
                self.assertEqual(self.run_tool("--output", output)[0], 1)
        self.assertEqual(sentinel.read_text(), "preserve")
        self.assertFalse((self.root / "new-workspace").exists())


if __name__ == "__main__":
    unittest.main()
