import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

PY_DIR = Path(__file__).resolve().parents[1] / "resources" / "scripts" / "py"
sys.path.insert(0, str(PY_DIR))

import check_outputs  # noqa: E402
import prepare_scatter_inputs as psi  # noqa: E402


class ManifestTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.dir = Path(self.tmp.name)
        for name in ("a.bam", "a.bai", "b.bam", "b.bai"):
            (self.dir / name).write_text("x")

    def tearDown(self):
        self.tmp.cleanup()

    def read(self, text):
        path = self.dir / "m.tsv"
        path.write_text(text, encoding="utf-8")
        return psi.read_manifest(path, self.dir)

    def test_tsv_and_csv(self):
        for sep in ("\t", ","):
            text = sep.join(["sample_id", "bam_path", "index_path"]) + "\n" + sep.join(["S1", "a.bam", "a.bai"]) + "\n"
            samples, crams, _ = self.read(text)
            self.assertEqual(samples, ["S1"])
            self.assertEqual(crams[0], (self.dir / "a.bam").resolve())

    def test_header_aliases_and_blank_rows(self):
        text = "Sample,CRAM,CRAI\nS1,a.bam,a.bai\n,,\nS2,b.bam,b.bai\n"
        self.assertEqual(self.read(text)[0], ["S1", "S2"])

    def test_duplicate_sample(self):
        text = "sample,bam,index\nS1,a.bam,a.bai\nS1,b.bam,b.bai\n"
        with self.assertRaisesRegex(ValueError, "duplicate"):
            self.read(text)

    def test_invalid_sample_id(self):
        with self.assertRaisesRegex(ValueError, "invalid sample ID"):
            self.read("sample,bam,index\nbad id,a.bam,a.bai\n")

    def test_missing_file(self):
        with self.assertRaisesRegex(ValueError, "does not exist"):
            self.read("sample,bam,index\nS1,nope.bam,a.bai\n")

    def test_missing_column(self):
        with self.assertRaisesRegex(ValueError, "exactly one index"):
            self.read("sample,bam,other\nS1,a.bam,a.bai\n")

    def test_no_rows(self):
        with self.assertRaisesRegex(ValueError, "no data rows"):
            self.read("sample,bam,index\n")


class KeyTests(unittest.TestCase):
    def test_wrapper_key(self):
        self.assertEqual(psi.wrapper_key("SwirlPipeline.x"), "SwirlPipeline.x")
        self.assertEqual(psi.wrapper_key("MitochondriaPipeline.x"), "SwirlPipeline.x")
        with self.assertRaises(ValueError):
            psi.wrapper_key("MitochondriaPipelineWrapper.x")
        with self.assertRaises(ValueError):
            psi.wrapper_key("x")

    def test_options_keep_null_and_zero(self):
        out = psi.read_json_options({"SwirlPipeline.a": None, "SwirlPipeline.b": 0.0, "SwirlPipeline.c": False}, "src")
        self.assertEqual(out, {"SwirlPipeline.a": None, "SwirlPipeline.b": 0.0, "SwirlPipeline.c": False})

    def test_conflicting_values(self):
        with self.assertRaisesRegex(ValueError, "Conflicting"):
            psi.merge_inputs({"k": 1}, {"k": 2}, "src")
        psi.merge_inputs({"k": 1}, {"k": 1}, "src")


class FilesJsonTests(unittest.TestCase):
    def test_resolves_relative_and_rejects_missing(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root / "ref.fa").write_text("x")
            good = root / "good.json"
            good.write_text(json.dumps({"MitochondriaPipeline.ref": "ref.fa"}))
            out = psi.read_json_inputs(good, root, resolve_files=True)
            self.assertEqual(out, {"SwirlPipeline.ref": str((root / "ref.fa").resolve())})
            bad = root / "bad.json"
            bad.write_text(json.dumps({"MitochondriaPipeline.ref": "missing.fa"}))
            with self.assertRaisesRegex(ValueError, "does not exist"):
                psi.read_json_inputs(bad, root, resolve_files=True)

    def test_helper_script_keys_skipped(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "s.json"
            path.write_text(json.dumps({"SwirlPipeline.prepare_script": "x.py", "SwirlPipeline.k": "v"}))
            out = psi.read_json_inputs(path, Path(tmp), resolve_files=False)
            self.assertEqual(out, {"SwirlPipeline.k": "v"})


class ReadScatterConfigTests(unittest.TestCase):
    def run_reader(self, root, config):
        cfg = root / "cfg.json"
        cfg.write_text(json.dumps(config))
        return subprocess.run(
            [sys.executable, str(PY_DIR / "read_scatter_config.py"), "--config-json", str(cfg), "--repo-root", str(root)],
            capture_output=True,
            text=True,
        )

    def make_repo(self, root):
        for name in ("wf.wdl", "prep.py", "stage.py", "hc.jar"):
            (root / name).write_text("x")
        (root / "files.json").write_text("{}")
        (root / "scripts.json").write_text(
            json.dumps({"SwirlPipeline.prepare_script": "prep.py", "SwirlPipeline.stage_script": "stage.py"})
        )
        return {
            "workflow": {
                "wdl": "wf.wdl",
                "files": "files.json",
                "scripts": "scripts.json",
                "output_root": "out",
                "dependencies": {"haplocheck_cli": "hc.jar", "tools_path": ""},
                "options": {},
            }
        }

    def test_outputs_eight_lines(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp).resolve()
            result = self.run_reader(root, self.make_repo(root))
            self.assertEqual(result.returncode, 0, result.stderr)
            lines = result.stdout.split("\n")[:-1]
            self.assertEqual(len(lines), 8)
            self.assertEqual(lines[1], str(root / "wf.wdl"))
            self.assertEqual(lines[7], str(root / "hc.jar"))

    def test_missing_options_fails(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp).resolve()
            config = self.make_repo(root)
            del config["workflow"]["options"]
            result = self.run_reader(root, config)
            self.assertEqual(result.returncode, 1)
            self.assertIn("options", result.stderr)


class CheckOutputsTests(unittest.TestCase):
    def populate(self, out, records=2):
        for name in check_outputs.MERGED + ["S1.coverage.tsv", "S1.statistics.tsv"]:
            (out / name).write_text("x")
        (out / "S1.final.vcf").write_text("##h\n#CHROM\n" + "chrM\t1\n" * records)

    def run_check(self, out):
        return subprocess.run(
            [sys.executable, str(PY_DIR / "check_outputs.py"), str(out)], capture_output=True, text=True
        )

    def test_ok(self):
        with tempfile.TemporaryDirectory() as tmp:
            self.populate(Path(tmp))
            result = self.run_check(Path(tmp))
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertIn("S1: 2 VCF records", result.stdout)

    def test_missing_and_empty(self):
        with tempfile.TemporaryDirectory() as tmp:
            out = Path(tmp)
            self.populate(out, records=0)
            (out / "merged_statistics.tsv").unlink()
            result = self.run_check(out)
            self.assertEqual(result.returncode, 1)
            self.assertIn("merged_statistics.tsv", result.stderr)
            self.assertIn("no records", result.stderr)


if __name__ == "__main__":
    unittest.main()
