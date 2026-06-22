# Release workflow

GaussJacobiQuad releases are driven by [cocogitto (`cog`)](https://docs.cocogitto.io/)
for semver inference, version-file updates, and annotated tags. [towncrier](https://towncrier.readthedocs.io/)
owns the curated user-facing notes in [`NEWS.md`](../NEWS.md) via fragments in
[`newsfragments/`](../newsfragments/). There is **no** separate `tbump` step; that
tool was retired in favour of this single path.

Configuration lives in [`cog.toml`](../cog.toml) and [`towncrier.toml`](../towncrier.toml).
Version strings are patched by [`scripts/set_version.sh`](../scripts/set_version.sh).

## Roles of each tool

| Tool | Owns | Does not own |
|------|------|----------------|
| **cog** | Conventional-commit lint, semver choice, tag `vX.Y.Z`, release commit | `NEWS.md` body (disabled via `disable_changelog = true`) |
| **towncrier** | `NEWS.md` sections from `newsfragments/*` | Git tags / semver |
| **set_version.sh** | `fpm.toml`, `meson.build`, Doxygen `PROJECT_NUMBER`, `CITATION.cff` | Commit message / tag |
| **add_headers.py** | License / version headers in `src/`, `interfaces/`, `scripts/` | Version in build manifests |
| **GitHub Actions** (`release.yml`) | GitHub Release on `v*` tag push (Zenodo DOI minting) | Choosing the version |

## Conventional commits

PR commits are linted by the `commit-lint` workflow (`cog check` via
`cocogitto/cocogitto-action`). Preferred types going forward:

| Type | Semver impact | Example |
|------|---------------|---------|
| `feat` | minor | `feat: add Clenshaw-Curtis nodes` |
| `fix` | patch | `fix: correct off-by-one in imtqlx` |
| `docs`, `chore`, `ci`, `build`, `test`, `refactor`, `perf`, `style` | none (unless `BREAKING CHANGE`) | `chore: tidy pre-commit config` |
| `feat!` / footer `BREAKING CHANGE:` | major | — |

Legacy types (`enh`, `bug`, `maint`, `doc`, `bld`, `tst`, plus Title-Case
`ENH`/`BUG`/`MAINT`/`DOC`/`CI`/`REL`) remain accepted so existing branches can
still pass CI; new work should use the standard lowercase forms above.

Install a local commit-msg hook optionally:

```bash
cog install-hook commit-msg
# or enforce style only via pre-commit (already configured for formatting)
pipx run pre-commit install
```

## Pre-release checklist

Before running `cog bump` on `main`:

1. Every PR intended for the release is merged; CI is green on `main`.
2. User-facing changes have a towncrier fragment under `newsfragments/`, named
   `<issue-or-pr>.<type>.md` where `<type>` is one of `feature`, `bugfix`,
   `doc`, `removal`, `misc`. Example:

   ```bash
   echo "Faster Golub-Welsch path for large n." > newsfragments/42.feature.md
   ```

3. Working tree is clean (`git status` shows nothing to commit).
4. You are on `main` (`cog.toml` sets `branch_whitelist = ["main"]`).
5. `cog` and `towncrier` are available (`cargo install cocogitto` or your distro
   package; `towncrier` via the project conda env in `environment.yml`).

Dry-run first to confirm the inferred version:

```bash
cog bump --auto --dry-run
```

Override with `--patch`, `--minor`, `--major`, or `--version X.Y.Z` when the
commit log would choose the wrong semver (e.g. only `feat(docs)` commits).

## Cutting a release

On an up-to-date `main`:

```bash
# Activate the project env if you use conda/mamba (for towncrier + header scripts)
micromamba activate gaussjacquad   # or: mamba activate gaussjacquad

cog bump --auto --annotated "Release vX.Y.Z

Summary of headline changes for maintainers / Zenodo.
See NEWS.md for the full user-facing list."
```

What happens (from `cog.toml` `pre_bump_hooks`, then cog itself):

1. **Version pick** — semver increment from conventional commits since the latest
   `v*` tag (or from your explicit flag).
2. **`scripts/set_version.sh`** — writes `X.Y.Z` into `fpm.toml`, `meson.build`,
   `apidocs/Doxygen-GaussJacobiQuad.cfg`, and `vX.Y.Z` into `CITATION.cff`.
3. **License headers** — `scripts/add_headers.py` refreshes Fortran/C/Python
   headers (commit / date / citation line).
4. **towncrier** — consumes `newsfragments/*`, prepends a dated section to
   `NEWS.md`, and removes consumed fragments.
5. **Release commit + annotated tag** — `vX.Y.Z` (prefix from `tag_prefix`).
6. **post_bump_hooks** — prints a push reminder.

Inspect, then publish:

```bash
git show HEAD
git push origin main --follow-tags
```

Pushing the tag triggers [`.github/workflows/release.yml`](../.github/workflows/release.yml),
which creates (or skips if already present) a GitHub Release with
`--generate-notes`. That release is what Zenodo watches for DOI minting.

## Install / export sanity after a release

Confirm consumers can still build and install from the tagged tree.

### fpm (library install)

`fpm.toml` sets `[install] library = true` so `fpm install` ships the static
library and `.mod` files, not only executables.

```bash
fpm build
fpm test
fpm install --prefix /tmp/gjq-prefix
# expects libGaussJacobiQuad.a (or equivalent) and modules under the prefix
```

### meson (subproject / system install)

`meson.build` installs the Fortran library (correct name `GaussJacobiQuad`),
executables, the C interop library, the public C header, and exposes
`gaussjacobiquad_dep` via `declare_dependency` for `subproject()` / wrap use.
An optional `gaussjacobiquad.pc` is installed for pkg-config consumers.

```bash
meson setup builddir --prefix=/tmp/gjq-prefix
meson compile -C builddir
meson install -C builddir
pkg-config --modversion gaussjacobiquad   # if PKG_CONFIG_PATH includes prefix/lib/pkgconfig
```

### Single-file export (vendoring)

Unaffected by semver tooling; still:

```bash
python scripts/export_single.py --modulename gjp_algo665
python scripts/export_single.py --modulename gjp_gw --keep-comments
```

Run header refresh (or a full `cog bump`) before vendoring into other projects.

## Recovery / common mistakes

| Situation | Action |
|-----------|--------|
| `cog bump` refused (not on `main`, dirty tree) | Merge/rebase onto `main`, commit or stash, retry |
| Wrong semver inferred | `cog bump --version X.Y.Z` (or `--patch` / `--minor` / `--major`) |
| Forgot a news fragment | Add it before bump; or edit `NEWS.md` in a follow-up `docs:` commit |
| Tag pushed but no GitHub Release | Re-run the `Release` workflow or `gh release create vX.Y.Z --generate-notes` |
| Need only lint commits in a PR | CI runs `cocogitto-action`; locally: `cog check origin/main..HEAD` |

## Deprecated: tbump

The former [`tbump.toml`](../tbump.toml) flow (`tbump` + `REL:` commits + ad-hoc
towncrier invocation) is superseded by `cog bump`. Do not run `tbump` on this
repository; it would diverge version sources and commit conventions.
