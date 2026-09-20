# Releasing levi

This repository mirrors the Bioconductor `devel` branch of levi and adds the
GitHub Actions workflows in `.github/workflows/`. Bioconductor is the source
of truth for the package; releases here are snapshots of it.

## Branches and tags

| Ref | Meaning |
|---|---|
| `master` | current package, kept in step with Bioconductor `devel` |
| `release-1.20`, `v1.20.0` | the version originally submitted to Bioconductor |
| `v2.0.0` | the rewritten package (Bioconductor devel 1.99.0, released as 2.0.0) |
| `gh-pages` | pkgdown site, written by the `pkgdown` workflow; never edit by hand |

## Updating `master` from Bioconductor

1. Push the change to Bioconductor first (`git push origin devel` in the
   Bioconductor clone) and wait for a clean build report.
2. Copy the package into this repository, keeping the GitHub-only files:

   ```bash
   cd ~/Documents/GitHub/levi
   rsync -a --delete \
     --exclude '.git' --exclude '.github' --exclude '.gitattributes' \
     --exclude 'codemeta.json' --exclude 'src/*.o' --exclude 'src/*.so' \
     ~/Documents/Bioconductor/Levi_Home/levi/ ./
   ```

3. Update the `"version"` field of `codemeta.json` to match `DESCRIPTION`.
4. Commit and push `master`. Every push runs `R-CMD-check`, `test-coverage`
   and `pkgdown`.

## Creating a release

The `release` workflow runs when a tag starting with `v` is pushed. It builds
the source tarball, checks it on the Bioconductor devel image and publishes a
GitHub Release whose notes are the matching section of `NEWS`.

Before tagging, make sure the three version numbers agree:

```bash
grep '^Version' DESCRIPTION                # Version: 2.0.1
grep '^Changes in version' NEWS | head -1  # Changes in version 2.0.1
```

Then create an annotated tag and push it:

```bash
git checkout master
git pull
git tag -a v2.0.1 -m "levi 2.0.1: short description of the release"
git push origin v2.0.1
```

Rules:

- The tag must start with `v`, otherwise the workflow does not run.
- The number after `v` must equal the `Version` in `DESCRIPTION` and the
  `Changes in version` heading in `NEWS`. If they differ, the release is still
  created but its notes read "No NEWS section found for version X".
- Bioconductor assigns the even release version itself. While the package is
  in `devel` as `x.99.0`, tag the GitHub release with the version it will
  become (for example `v2.0.0` for `1.99.0`) and say so in the tag message.

## Fixing a wrong tag

If the tag was pushed but the release is wrong, delete the GitHub Release on
the Releases page first, then remove the tag and start again:

```bash
git tag -d v2.0.1                   # local
git push origin :refs/tags/v2.0.1   # GitHub
```

## Continuous integration

| Workflow | Runs on | What it does |
|---|---|---|
| `R-CMD-check` | push and pull request to `master` or `devel` | `R CMD build`, `R CMD check` and `BiocCheck` on `bioconductor/bioconductor_docker:devel` (required); `R CMD check` on macOS and Windows with the Bioconductor release (informative only) |
| `test-coverage` | same | `covr` on the devel image; uploads to Codecov when the `CODECOV_TOKEN` secret exists |
| `pkgdown` | push to `master`, releases | builds the site on the devel image and deploys it to `gh-pages` from a plain runner |
| `release` | tags `v*` | tarball, check and GitHub Release |

Two BiocCheck items are reported but do not fail the job: the bioc-devel
mailing-list subscription and the support-site registration. Both query
external services (the first needs an admin password) and fail at random on a
runner. Any other BiocCheck error fails the build.

The first run after a dependency change is slow (20 to 30 minutes) because the
Suggests are installed from scratch; later runs use the cache.
