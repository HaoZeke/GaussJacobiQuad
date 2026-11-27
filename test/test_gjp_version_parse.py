"""Findings 1-3: pre-release suffixes, a 2-component version, and compact bounds."""

import subprocess
import sys
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


if __name__ == "__main__":
    unittest.main()
