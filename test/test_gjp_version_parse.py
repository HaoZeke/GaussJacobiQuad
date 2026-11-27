# BEGIN_HEADER
# -----------------------------------------------------------------------------
# Gauss-Jacobi Quadrature Implementation
# Authors: Rohit Goswami <rgoswami[at]ieee.org>
# Source: GaussJacobiQuad Library
# License: MIT
# GitHub Repository: https://github.com/HaoZeke/GaussJacobiQuad
# -----------------------------------------------------------------------------
# END_HEADER
"""Findings 1-3: pre-release suffixes, a 2-component version, and compact bounds."""

import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "parse_gjp_version.py"


class VersionError(Exception):
    pass


def legacy_parse(version):
    """Parser in e84ad222 meson.build. Mis-handles rc suffixes, 1.2, and bounds."""
    parts = version.split(".")
    major = int(parts[0])
    minor = int(parts[1])
    patch = int(parts[2].split("-")[0].split("+")[0])
    return major, minor, patch


def parse_version(version):
    if SCRIPT.is_file():
        proc = subprocess.run(
            [sys.executable, str(SCRIPT), version],
            capture_output=True,
            text=True,
        )
        if proc.returncode != 0:
            raise VersionError(proc.stderr.strip() or f"exit {proc.returncode}")
        return tuple(int(part) for part in proc.stdout.split())
    return legacy_parse(version)


class VersionParseTests(unittest.TestCase):
    def test_glued_prerelease(self):
        self.assertEqual(parse_version("1.2.0rc1"), (1, 2, 0))

    def test_hyphen_prerelease(self):
        self.assertEqual(parse_version("1.2.0-rc.1"), (1, 2, 0))

    def test_two_components_default_patch(self):
        self.assertEqual(parse_version("1.2"), (1, 2, 0))

    def test_build_metadata(self):
        self.assertEqual(parse_version("1.2.3+gabcdef"), (1, 2, 3))

    def test_dotted_dev_post_rc(self):
        self.assertEqual(parse_version("1.2.0.dev1"), (1, 2, 0))
        self.assertEqual(parse_version("1.2.3.post1"), (1, 2, 3))
        self.assertEqual(parse_version("1.2.0.rc1"), (1, 2, 0))
        self.assertEqual(parse_version("1.2.0rc.1"), (1, 2, 0))

    def test_fourth_component_rejected(self):
        with self.assertRaises(VersionError):
            parse_version("1.2.3.4")

    def test_junk_last_component_rejected(self):
        for version in ("1.2.3&x", "1.2.3/foo"):
            with self.subTest(version=version):
                with self.assertRaises(VersionError):
                    parse_version(version)

    def test_project_version(self):
        self.assertEqual(parse_version("0.2.6"), (0, 2, 6))

    def test_component_above_99_is_rejected(self):
        for version in ("1.100.0", "1.2.100", "100.0.0"):
            with self.subTest(version=version):
                with self.assertRaises(VersionError):
                    parse_version(version)

    def test_in_range_neighbors_do_not_collide(self):
        self.assertNotEqual(10000 * 1 + 100 * 99 + 99, 10000 * 2)
        self.assertEqual(parse_version("1.99.99"), (1, 99, 99))
        self.assertEqual(parse_version("2.0.0"), (2, 0, 0))

    def test_unchecked_packing_collides_past_99(self):
        self.assertEqual(10000 * 1 + 100 * 100 + 0, 10000 * 2)
        self.assertEqual(10000 * 1 + 100 * 2 + 100, 10000 * 1 + 100 * 3)

    def test_apply_keeps_suffix_outside_fpm(self):
        # --apply roots itself at the script's parent, so the copy is the tree.
        for version in ("1.2.0rc1", "1.2.0-rc.1"):
            with self.subTest(version=version):
                with tempfile.TemporaryDirectory() as tmp:
                    root = Path(tmp)
                    _copy_apply_tree(root)
                    proc = _apply(root, version)
                    self.assertEqual(proc.returncode, 0, proc.stderr)
                    self.assertEqual(proc.stdout.strip(), "1 2 0")
                    fpm = (root / "fpm.toml").read_text()
                    meson = (root / "meson.build").read_text()
                    source = (root / "src" / "gjp_version.f90").read_text()
                    self.assertIn('version = "1.2.0"', fpm)
                    self.assertNotIn(version, fpm)
                    town = (root / "towncrier.toml").read_text()
                    self.assertIn('version = "1.2.0"', town)
                    self.assertNotIn(version, town)
                    self.assertIn(f"version: '{version}'", meson)
                    self.assertIn(f'gjp_version_string = "{version}"', source)
                    self.assertIn("gjp_version_major = 1", source)
                    self.assertIn("gjp_version_minor = 2", source)
                    self.assertIn("gjp_version_patch = 0", source)

    def test_apply_rejects_fourth_component_without_writing(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            _copy_apply_tree(root)
            before = {
                rel: (root / rel).read_text()
                for rel in (
                    "fpm.toml",
                    "meson.build",
                    "src/gjp_version.f90",
                )
            }
            proc = _apply(root, "1.2.3.4")
            self.assertNotEqual(proc.returncode, 0)
            for rel, text in before.items():
                self.assertEqual((root / rel).read_text(), text)

    def test_apply_rewrites_towncrier(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            _copy_apply_tree(root)
            before = (root / "towncrier.toml").read_text()
            self.assertIn('version = "0.2.6"', before)
            proc = _apply(root, "9.9.9")
            self.assertEqual(proc.returncode, 0, proc.stderr)
            self.assertEqual(proc.stdout.strip(), "9 9 9")
            after = (root / "towncrier.toml").read_text()
            self.assertIn('version = "9.9.9"', after)
            self.assertNotIn('version = "0.2.6"', after)
            self.assertIn('version = "9.9.9"', (root / "fpm.toml").read_text())


_APPLY_FILES = (
    "fpm.toml",
    "meson.build",
    "apidocs/Doxygen-GaussJacobiQuad.cfg",
    "CITATION.cff",
    "src/gjp_version.f90",
    "scripts/parse_gjp_version.py",
    "towncrier.toml",
)


def _copy_apply_tree(dest: Path) -> None:
    for rel in _APPLY_FILES:
        target = dest / rel
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy(ROOT / rel, target)


def _apply(root: Path, version: str) -> subprocess.CompletedProcess[str]:
    script = root / "scripts" / "parse_gjp_version.py"
    return subprocess.run(
        [sys.executable, str(script), "--apply", version],
        capture_output=True,
        text=True,
    )


if __name__ == "__main__":
    unittest.main()
