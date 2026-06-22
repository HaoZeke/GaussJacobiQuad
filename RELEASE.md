# Release engineering

This project uses **[cocogitto](https://docs.cocogitto.io/)** (`cog`) as the single
release-engineering tool. It enforces [Conventional Commits](https://www.conventionalcommits.org/),
generates [NEWS.md](NEWS.md), and drives semantic version bumps and git tags.

The previous **tbump + towncrier + newsfragments** workflow is retired. Stub
configs remain only as deprecation markers; do not invoke `tbump` or `towncrier`
for new releases.

## Tooling overview

| Concern              | Tool                         | Config / output          |
|----------------------|------------------------------|--------------------------|
| Commit message lint  | `cog verify` / `cog check`   | `cog.toml`               |
| Changelog            | `cog changelog` / `cog bump` | `NEWS.md`                |
| Version bump + tag   | `cog bump`                   | `fpm.toml`, `meson.build`, `CITATION.cff`, Doxygen cfg |
| GitHub Release       | `.github/workflows/release.yml` on `v*` tags | uses changelog notes |
| PR commit lint       | `.github/workflows/commit_lint.yml` | `cog check` on PR range |

Install `cog` (once per machine):

```bash
cargo install cocogitto
# or grab a binary from https://github.com/cocogitto/cocogitto/releases
cog --version
```

Optional local git hook so bad commit subjects never land:

```bash
cog install-hook commit-msg
```

Build/test dependencies stay in [environment.yml](environment.yml) (conda/mamba).
`cog` is intentionally not in that env because it is a standalone Rust binary.

## Commit conventions

Subjects must match Conventional Commits, for example:

```text
feat: add Golub-Welsch CLI benchmark harness
fix: correct off-by-one in IMTQLX rotation
docs: document meson install layout
ci: run cog check on pull requests
chore: tidy export_single script
```

Types that affect semver (see `cog.toml`):

- `feat` / `enh` → **minor**
- `fix` / `bug` / `perf` → **patch**
- any commit with a `BREAKING CHANGE:` footer (or `!` after the type) → **major**
- `docs`, `style`, `refactor`, `test`, `build`, `ci`, `chore`, … → no automatic bump alone

Scopes are optional (`feat(cli): …`). Keep the subject imperative and under ~72
characters. Prefer one logical change per commit so the changelog stays readable.

## Package install layout (fpm and meson)

Both backends must install a **usable library package**, not only executables.

### fpm

[fpm.toml](fpm.toml) sets `[install] library = true`, so `fpm install` ships the
Fortran library and modules for dependents that consume this package via fpm.

```bash
fpm build
fpm install --prefix "$HOME/.local"
```

### meson

[meson.build](meson.build) installs:

- `libGaussJacobiQuad` (core Fortran library; name fixed from the former typo `GaussJacobiQaud`)
- `libgjp_cinterp` (ISO_C_BINDING + C shim)
- CLI executables (`gjp_quad`, `gjp_quad_rec`, `gjp_quad_gw`, `gjp_quad_algo665`, `c_cli_gjpq`)
- C header: `$prefix/include/gaussjacobiquad/GaussJacobiQuadCInterp.h`
- Fortran `.mod` files: `$prefix/include/gaussjacobiquad/*.mod`
- pkg-config file: `gaussjacobiquad.pc` (variable `moddir` points at the module dir)

```bash
meson setup builddir --prefix="$HOME/.local"
meson compile -C builddir
meson install -C builddir

# Consumers (C / mixed)
pkg-config --cflags --libs gaussjacobiquad

# Meson subproject / dependency('gaussjacobiquad') also works in-tree via
# meson.override_dependency in this build file.
```

Verify install quickly:

```bash
test -f "$HOME/.local/include/gaussjacobiquad/GaussJacobiQuadCInterp.h"
pkg-config --exists gaussjacobiquad && pkg-config --modversion gaussjacobiquad
```

## Cutting a release

Perform releases from an up-to-date `main` with a clean worktree.

### 1. Preconditions

```bash
git checkout main
git pull --ff-only origin main
git status   # must be clean
cog check --from-latest-tag --ignore-merge-commits
fpm build && pytest -vvv   # or: meson setup b && meson compile -C b
```

All commits since the latest `v*` tag should pass `cog check`. Fix outliers with
`git rebase -i` / `cog edit` before bumping.

### 2. Preview the next version and changelog

```bash
cog bump --auto --dry-run          # prints suggested semver
cog changelog                      # preview since last tag
# or force a channel:
cog bump --minor --dry-run
cog bump --patch --dry-run
cog bump --major --dry-run
```

If you still need to mention work that predates conventional commits (the old
newsfragments), edit the **Unreleased** section of `NEWS.md` first; `cog bump`
will prepend the generated section above the historical entries.

### 3. Bump, commit, and tag

```bash
cog bump --auto
# equivalents: cog bump --minor | --patch | --major
# or explicit: cog bump --version 0.2.0
```

`cog bump` will:

1. Run `pre_bump_hooks` from [cog.toml](cog.toml) (rewrite version fields in
   `fpm.toml`, `meson.build`, Doxygen config, `CITATION.cff`, and refresh headers).
2. Regenerate the top of [NEWS.md](NEWS.md) from conventional commits.
3. Create a release commit and an annotated-style tag `vX.Y.Z` (prefix from `tag_prefix`).
4. Run `post_bump_hooks` (push branch + tag to `origin`).

If you prefer to push yourself, temporarily comment out the `post_bump_hooks` in
`cog.toml`, run `cog bump`, then:

```bash
git push origin main
git push origin "v$(cog get-version)"
```

### 4. GitHub Release and Zenodo

Pushing the `v*` tag triggers [.github/workflows/release.yml](.github/workflows/release.yml),
which creates a GitHub Release whose notes are taken from the matching section of
`NEWS.md` (falling back to generated notes). With Zenodo–GitHub integration
enabled on the repository, the new GitHub Release mints / updates the DOI.

Confirm:

```bash
gh release view "v$(cog get-version)"
```

### 5. Post-release housekeeping

- Delete any obsolete files under `newsfragments/` that were already migrated.
- Bump `date-released` in `CITATION.cff` if you publish a citable snapshot.
- Announce / update badges if the Zenodo concept DOI badge needs a refresh.

## CI expectations

| Workflow            | Trigger        | What it does                                      |
|---------------------|----------------|---------------------------------------------------|
| `commit_lint.yml`   | pull_request   | Installs `cog`, runs `cog check` on PR commits    |
| `release.yml`       | push `v*` tags | Publishes GitHub Release from `NEWS.md`           |
| `build_test.yml`    | push / PR      | `fpm build` + `pytest` in the conda env           |
| `pre_commit.yml`    | push / PR      | formatting / lint hooks                           |
| `build_docs.yml`    | push / PR      | Doxygen site                                      |

PR authors should use conventional commit subjects so `commit_lint` stays green.
Rebase rather than merge commits when possible (`ignore_merge_commits = true`
softens this, but clean history is easier to changelog).

## Recovering from a bad bump

If `cog bump` fails mid-flight before the tag is pushed:

```bash
git status
# restore versioned files if needed
git checkout -- fpm.toml meson.build CITATION.cff apidocs/Doxygen-GaussJacobiQuad.cfg NEWS.md
git reset --hard HEAD~1   # only if the bump commit was created locally and not pushed
```

If the tag was pushed erroneously, delete it locally and on the remote (coordinate
with anyone who may have already consumed it), fix commits, and re-run the bump
with the correct version:

```bash
git tag -d "vX.Y.Z"
git push origin :refs/tags/vX.Y.Z
# fix history / manifests, then:
cog bump --version X.Y.Z
```

## Why not tbump / towncrier anymore?

- **Duplication**: tbump patched version files; towncrier assembled `NEWS.md`
  from hand-written fragments; neither enforced commit quality on its own.
- **Single source of truth**: conventional commits already describe user-visible
  changes. `cog` turns that stream into both the semver decision and the
  changelog without a parallel fragment tree.
- **CI symmetry**: the same `cog check` gate runs on PRs and before a human
  cuts a release locally.

If you discover a workflow still calling `tbump` or `towncrier`, update it to
`cog` and point at this document.
