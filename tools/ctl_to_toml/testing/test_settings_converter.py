"""Regression tests; all file writes are isolated in temporary directories."""

from datetime import date
import importlib.machinery
import importlib.util
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[3]
CLI = Path(__file__).resolve().parents[1] / "settings_converter"
sys.path.insert(0, str(CLI.parent))
from control_document import ControlDocument, parse_legacy_value, tomllib

loader = importlib.machinery.SourceFileLoader("settings_converter", str(CLI))
spec = importlib.util.spec_from_loader(loader.name, loader)
converter = importlib.util.module_from_spec(spec)
loader.exec_module(converter)


class ConversionTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.directory = Path(self.temp.name)

    def run_cli(self, *arguments):
        return subprocess.run(
            [sys.executable, str(CLI), *map(str, arguments)],
            cwd=self.directory, capture_output=True, text=True,
        )

    def test_types_comments_order_and_roundtrip(self):
        source = '''! Section heading

<n> +42 ! integer
<real> 2.0 ! floating
<tiny> 1e-3 ! exponent
<fortran> 1D-3 ! Fortran exponent
<on> T ! enabled
<off> F ! disabled
<stamp> 2026-11-24 12:30:00 ! date
<text> "0012" ! explicit string
<path> /tmp/a#b ! path
! End
'''
        document = ControlDocument.parse(source, ".control")
        output = document.render(".toml")
        values = tomllib.loads(output)
        self.assertEqual(list(values), ["n", "real", "tiny", "fortran", "on", "off", "stamp", "text", "path"])
        self.assertIs(type(values["n"]), int)
        self.assertIs(type(values["real"]), float)
        self.assertEqual(values["tiny"], 0.001)
        self.assertEqual(values["fortran"], 0.001)
        self.assertIs(values["on"], True)
        self.assertIs(values["off"], False)
        self.assertEqual(values["text"], "0012")
        self.assertEqual(values["stamp"], "2026-11-24 12:30:00")
        legacy = ControlDocument.parse(output, ".toml").render(".control")
        self.assertIn("<on> T ! enabled", legacy)
        self.assertIn("<off> F ! disabled", legacy)
        self.assertIn("<text> 0012 ! explicit string", legacy)
        self.assertEqual(
            [(entry.key, entry.comment) for entry in document.entries],
            [(entry.key, entry.comment) for entry in ControlDocument.parse(legacy, ".control").entries],
        )

    def test_toml_string_syntax_and_dates(self):
        source = '''# Heading
path = 'C:\\data\\a#b' # path
text = "a\\\"b" # quote
stamp = 2026-11-24T12:30:00Z # date
single = """abc""" # multiline syntax
'''
        legacy = ControlDocument.parse(source, ".toml").render(".control")
        result = tomllib.loads(ControlDocument.parse(legacy, ".control").render(".toml"))
        self.assertEqual(result["path"], "C:\\data\\a#b")
        self.assertEqual(result["text"], 'a"b')
        self.assertEqual(result["stamp"], "2026-11-24 12:30:00+00:00")
        self.assertEqual(result["single"], "abc")

    def test_quoted_legacy_delimiters(self):
        result = ControlDocument.parse('<text> "a!b#c" ! comment\n', ".control").render(".toml")
        self.assertEqual(tomllib.loads(result)["text"], "a!b#c")
        self.assertIn("# comment", result)

    def test_repository_samples(self):
        for filename, suffix in (("SAMPLE.control", ".control"), ("SAMPLE_toml", ".toml"), ("SAMPLE-coupled_toml", ".toml")):
            with self.subTest(filename=filename):
                document = ControlDocument.parse((ROOT / "route/settings" / filename).read_text(), suffix)
                other = ".toml" if suffix == ".control" else ".control"
                restored = ControlDocument.parse(document.render(other), other)
                self.assertEqual(
                    [(entry.key, entry.value, entry.comment) for entry in document.entries],
                    [(entry.key, entry.value, entry.comment) for entry in restored.entries],
                )

    def test_cli_silence_verbosity_help_and_collisions(self):
        source = self.directory / "settings.control"
        source.write_text("<dt_qsim> 86400 ! seconds\n")
        original = source.read_bytes()
        today = date.today().strftime("%Y%m%d")
        for filename in ("settings.toml", f"settings_{today}.toml", f"settings_{today}_1.toml"):
            result = self.run_cli(source)
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertEqual((result.stdout, result.stderr), ("", ""))
            self.assertTrue((self.directory / filename).exists())
        self.assertEqual(source.read_bytes(), original)
        result = self.run_cli("-v", self.directory / "settings.toml")
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn("Converted 1 settings", result.stderr)
        self.assertTrue((self.directory / f"settings_{today}.control").exists())
        self.assertEqual(source.read_bytes(), original)
        result = self.run_cli("--help")
        self.assertEqual(result.returncode, 0)
        self.assertIn("--verbose", result.stdout)

    def test_invalid_inputs_have_no_output(self):
        cases = [
            ("wrong.control", 'name = "value"'),
            ("wrong.toml", '<name> value ! comment'),
            ("bad.toml", 'name = "unterminated'),
            ("bad.control", '<name> "unterminated ! comment'),
            ("duplicate.control", '<name> 1 !\n<name> 2 !'),
            ("duplicate.toml", 'name = 1\nname = 2'),
            ("empty.control", '! only comments'),
            ("empty.toml", '# only comments'),
            ("table.toml", '[table]\nname = 1'),
            ("array.toml", 'name = [1, 2]'),
            ("bang.toml", 'name = "a!b"'),
            ("newline.toml", 'name = "a\\nb"'),
            ("wrong.txt", '<name> value !'),
        ]
        for filename, content in cases:
            with self.subTest(filename=filename):
                source = self.directory / filename
                source.write_text(content)
                before = set(self.directory.iterdir())
                result = self.run_cli(source)
                self.assertNotEqual(result.returncode, 0)
                self.assertIn("error:", result.stderr)
                self.assertNotIn("Traceback", result.stderr)
                self.assertEqual(set(self.directory.iterdir()), before)
        self.assertNotEqual(self.run_cli(self.directory / "missing.control").returncode, 0)
        self.assertNotEqual(self.run_cli(self.directory).returncode, 0)

    def test_io_errors_and_failed_write_cleanup(self):
        source = self.directory / "input.control"
        with patch.object(Path, "open", side_effect=PermissionError("permission denied")):
            with self.assertRaises(PermissionError):
                converter.write_output(source, "value = 1\n")
        from unittest.mock import mock_open
        handle = mock_open()
        handle.return_value.write.side_effect = OSError("disk full")
        with patch.object(Path, "open", handle), patch.object(Path, "unlink") as unlink:
            with self.assertRaises(OSError):
                converter.write_output(source, "value = 1\n")
            unlink.assert_called_once()

    def test_nonfinite_floating_point_types(self):
        import math
        document = ControlDocument.parse("positive = inf\nnegative = -inf\nmissing = nan\n", ".toml")
        result = tomllib.loads(ControlDocument.parse(document.render(".control"), ".control").render(".toml"))
        self.assertEqual(result["positive"], math.inf)
        self.assertEqual(result["negative"], -math.inf)
        self.assertTrue(math.isnan(result["missing"]))

    def test_existing_logical_spellings(self):
        for text in ("T", "true", ".true.", "TRUE"):
            self.assertIs(parse_legacy_value(text), True)
        for text in ("F", "false", ".false.", "FALSE"):
            self.assertIs(parse_legacy_value(text), False)


if __name__ == "__main__":
    unittest.main()
