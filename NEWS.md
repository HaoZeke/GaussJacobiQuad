# Changelog

All notable changes to this project are documented in this file.

The format is driven by [Cocogitto](https://docs.cocogitto.io/) from
[Conventional Commits](https://www.conventionalcommits.org/). Older entries
below were produced with towncrier and are kept for historical continuity.

<!-- cog-changelog-marker -->

## Unreleased

Pending fragments migrated from the former `newsfragments/` directory (will be
folded into the next `cog bump` changelog section):

- **Features**: Added a variant of algorithm 665, refactored into `gjp_common`
  to accentuate similarities between 665 and standard GW (#7).
- **Experimental**: Added a single-file exporter for external libraries, which
  concatenates required files into one and updates the headers accordingly (#10).

## 0.1.0 (2023-08-26)

### Features

- Added a CLI for validating against SciPy
- Added an ISO_C_BINDING compatible interface and a C header
- Added an f2py generated wrapper for Python interoperability
- Implemented the Golub-Welsch algorithm with newer LAPACK routines
- Wrote a CLI for obtaining formatted results from SymPy

### Improved Documentation

- Added benchmarks and a discussion on the interfaces (#8)

### Misc

- Initial public release
