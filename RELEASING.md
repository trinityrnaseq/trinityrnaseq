# Releasing Trinity

This is the checklist for putting out a new Trinity release: the GitHub release with tarballs, the Docker image,
and the Bioconda recipe. It is written so that a developer, or an AI agent working for one, can follow it from
start to finish. Do the steps in order. **Nothing public (pushes, tags, releases, Docker Hub, Bioconda) happens
until every check in "Before release" has passed and the maintainer has signed off.**

Lessons from v2.16.0 are marked *(lesson)*. v2.16.0 was released before three items were resolved, which forced
v2.16.1 a few hours later.

---

## 0. Repository layout to know

- Develop on `devel`. `master` holds the latest release and is fast-forwarded from `devel` at release time.
- These submodules must be pushed to their own `master` *before* the superproject that pins them:
  - `Inchworm`, `Chrysalis`, `Butterfly`
  - `trinity-plugins/bamsifter`, `trinity-plugins/seqtk-trinity`
  - `trinity_ext_sample_data` (extended tests)
  - `trinityrnaseq.wiki`, the GitHub wiki repo, which can also be edited on the web
- The version string is set in **three** places, and all three must agree:
  - `Trinity`: `my $VERSION = "Trinity-vX.Y.Z";`
  - `Docker/VERSION.txt`: `X.Y.Z`
  - `Docker/Dockerfile`: `ENV TRINITY_VERSION="X.Y.Z"`
- `Changelog.txt`: the newest section goes at the top, with the header `# Trinity-vX.Y.Z Month D, YYYY`.
- Versioning: use a patch bump (Z) for fixes and a minor bump (Y) for new features or dependency changes.

## 1. Before release: resolve open decisions

Do these before any testing. Re-testing after a late change is cheap; re-releasing is not.

- [ ] **Deprecated or removed options:** remove them outright rather than leaving deprecated no-ops behind. The
      maintainer prefers to break old pipelines with a clear error message than to keep dead options around.
      *(lesson: in v2.16.0 `analyze_diff_expr.pl --include_GOplot` was left as a warning no-op; v2.16.1 removed it.)*
- [ ] **Stale tests:** every directory in `trinity_ext_sample_data/` and every script in
      `sample_data/test_Trinity_Assembly/misc_run_tests/` should either run or be deleted. Tests that reference
      removed scripts must go.
- [ ] **Third-party tool versions:** the Docker image pins exact versions. Bioconda resolves version ranges, so
      users get whatever conda picks.
  - For each tool Trinity drives (salmon, bowtie2, jellyfish, samtools, kallisto, RSEM, …), check the latest
    upstream release.
  - Either *test Trinity against the latest version* (section 2.4) or *pin the range to versions that were tested*,
    in both the Dockerfile and the Bioconda recipe.
  - Never ship a range that includes an untested major version. *(lesson: v2.16.0's recipe allowed salmon 2.x,
    which had only been tested against a fake salmon.)*
- [ ] **Version reporting:** `Trinity --version` must report the version being released. It also queries GitHub's
      `releases/latest` endpoint and compares versions numerically, so a release candidate that is newer than the
      latest published release is reported as newer, not as outdated. Confirm this works (section 2.5).
- [ ] **Changelog:** list every user-visible change, especially removed options and dependency changes.

## 2. Before release: testing

Run everything against the exact commits that will be released, including submodule commits, from a fresh
recursive clone. All of the testing can run without pushing anything.

**Environment notes for the Broad host, where these steps were developed:**
- There is no g++, cmake or java on the host, and host perl lacks threads. Build and test inside Apptainer.
- Use either an Apptainer sandbox built from the candidate `Docker/Dockerfile`, or one built from the published
  image: `apptainer build --sandbox trin docker://trinityrnaseq/trinityrnaseq:<ver>`.
- Use `APPTAINERENV_TMPDIR=/tmp`, `apptainer exec --cleanenv`, `LC_ALL=C`, and `TZ=UTC`, plus
  `TRINITY_HOME=<clone>`, `PICARD_HOME`, `GATK_HOME` and `LD_LIBRARY_PATH=/usr/local/lib` as in the Dockerfile.
- **Never submit to UGER or any grid.** Run everything locally.

### 2.1 Build and unit tests
- [ ] `make && make plugins && make test` all pass.

### 2.2 Full test framework (all tests, not only `make test_all`)
- [ ] `sample_data/test_Trinity_Assembly`: `runMe.sh` plus every `misc_run_tests/*.sh` listed in its `Makefile`.
      Run `./cleanme.pl` between scripts.
- [ ] `trinity_ext_sample_data`: every `test_*` directory, both the ones listed in its top-level `Makefile`
      (`make test_all`) and the extra ones such as Glimma, filter_low_expr, ExN50 and TissueSpecificityGraph.
- [ ] Run each test separately and record name, exit code and runtime, so that one failure doesn't hide the
      others.
- [ ] Also scan the logs for `Error`, `Traceback`, `died` and `Segmentation` lines in tests that still exited 0.

### 2.3 Full-length reconstruction regression benchmark
- [ ] Follow [`benchmarking/FL_regression/README.md`](benchmarking/FL_regression/README.md): mouse and
      S. pombe, 10M read pairs each, assembled with the previous release and with the candidate.
- [ ] Add rows to `benchmarking/FL_regression/results_history.tsv`.
- [ ] A drop of more than about 1% in CDS-based FL genes needs investigation before release.
- Runtime: about 2.3 h for pombe and about 5 h for mouse per assembly, at 6 CPUs.

### 2.4 Third-party tool compatibility (salmon in particular)
salmon is used in two places:
- `util/support_scripts/salmon_runner.pl`: the final expression filtering of every Trinity assembly.
- `util/align_and_estimate_abundance.pl --est_method salmon`. Its outputs feed
  `abundance_estimates_to_matrix.pl`, tximport and DTU.

- [ ] With the salmon version pinned in the Dockerfile, the tests in 2.2 already cover it.
- [ ] With the newest salmon allowed by the Bioconda recipe range, rerun at least:
  - `sample_data/test_Trinity_Assembly/runMe.sh`
  - `misc_run_tests/__run_PE_samples_file.sh` and `misc_run_tests/__runMe_include_long_reads.sh`
  - `trinity_ext_sample_data/test_align_and_estimate_abundance`, `test_tximport`, `test_ExN50`, `test_DTU`

  Official salmon binaries are on the COMBINE-lab/salmon GitHub releases page; verify the `.sha256`. Put the
  binary first on `PATH` inside the sandbox. Compare transcript counts and quant outputs against the run with the
  pinned salmon.
- [ ] Exit codes aren't enough. Also run `salmon_runner.pl` and `filter_transcripts_require_min_cov.pl` with
      each salmon version on the *same* pre-filter assembly (`Trinity --no_salmon`), and compare which transcripts
      are kept.
- [ ] If any of these fail, or the filtering changes materially, cap the version: Trinity's
      `check_required_3rd_party_tool_installations()` check, the Dockerfile, and the Bioconda range (for example
      `salmon >=1.10,<2`).

**Status as of v2.16.1 (2026-11):** Trinity requires salmon 1.x and is tested with 1.10.0.
- salmon 2.x (2.12.0, the Rust rewrite) runs and its exit codes are fine.
- On the sample assembly, though, it assigned zero reads to 11 of 89 transcripts, against 3 for salmon 1.10.
  That included isoforms with 200–490 reads under 1.10, so the final filter kept 74 transcripts instead of 82.
- salmon 2.x also cannot read indexes built by 1.x.
- Re-test before allowing salmon 2.x. Salmon 1.11 and 1.12 Linux binaries threw a locale exception on Ubuntu 24.04,
  which is why the Docker image uses 1.10.0.

### 2.5 Version check
- [ ] `Trinity --version` prints `Trinity version: Trinity-vX.Y.Z`.
- [ ] Before the release is published it should say the version is newer than the latest published release.
      After publishing, it should say it is the latest.

## 3. Release

Each step needs the maintainer's go-ahead. Pushing over https may not have credentials; SSH URLs
(`git@github.com:trinityrnaseq/<repo>.git`) work.

1. **Submodules:** push each changed submodule commit to that repo's `master`. Each push must be a fast-forward:
   if upstream `master` has commits the superproject never pinned, rebase onto them. Then confirm with
   `git ls-remote`.
2. **Wiki:** push it, rebasing first if someone edited it on the web. Then bump the `trinityrnaseq.wiki` pin in
   the superproject.
3. **Changelog:** replace any "(unreleased)" with the release date.
4. **Dockerfile:** set `ENV TRINITY_CO=<sha>` to the commit just *before* the Dockerfile change, so that the
   Dockerfile commit is the release commit. Commit, then push `devel`.
5. **master and tag:**
   - Fast-forward `master` to `devel` and push it.
   - Create the lightweight tag `Trinity-vX.Y.Z` on the same commit and push it.
6. **Tarballs:** build both from a fresh recursive clone of the tag:
   ```bash
   git clone --recursive --branch Trinity-vX.Y.Z https://github.com/trinityrnaseq/trinityrnaseq.git trinityrnaseq-vX.Y.Z
   find trinityrnaseq-vX.Y.Z -name .git -prune -exec rm -rf {} +
   tar czf trinityrnaseq-vX.Y.Z.FULL_with_extendedTestData.tar.gz trinityrnaseq-vX.Y.Z
   rm -rf trinityrnaseq-vX.Y.Z/trinity_ext_sample_data && mkdir trinityrnaseq-vX.Y.Z/trinity_ext_sample_data
   tar czf trinityrnaseq-vX.Y.Z.FULL.tar.gz trinityrnaseq-vX.Y.Z
   ```
   - The `.gitmodules` files and the wiki directory stay in both tarballs.
   - Expected sizes: about 45 MB for FULL and about 350 MB for the version with the extended test data.
   - Extract FULL into an empty directory, then run `make && make plugins` and `sample_data/test_Trinity_Assembly/runMe.sh`.
7. **GitHub release:**
   ```bash
   gh release create Trinity-vX.Y.Z -R trinityrnaseq/trinityrnaseq --title "Trinity Release-vX.Y.Z" \
     --notes-file notes.md trinityrnaseq-vX.Y.Z.FULL.tar.gz trinityrnaseq-vX.Y.Z.FULL_with_extendedTestData.tar.gz
   ```
   - The notes are the Changelog section plus a line naming the tarballs and the Docker tag.
   - Download the FULL asset again and check that its sha256 matches; Bioconda pins it.
8. **Docker:** build on a VM with Docker. The maintainer uses `bhaas-gce3`, reached with
   `gcloud compute ssh --zone us-central1-f bhaas-gce3 --project methods-dev-lab --tunnel-through-iap`.
   - In a fresh clone of the tag, in `Docker/`: run `build_docker.sh`, then `run_test_docker.sh` and
     `run_ext_test_docker.sh`. All must pass.
   - Then run `push_docker.sh`, which pushes both `:X.Y.Z` and `:latest`.
   - Singularity/Apptainer images are no longer built.
9. **Bioconda:** update `recipes/trinity` in bioconda/bioconda-recipes. Open the PR from a branch in that repo;
   the maintainer's account has push access.
   - Set `version`, `sha256` (of the FULL asset) and `build: number: 0`.
   - Re-check that every recipe patch still applies to the new tarball (`git apply --check`), and drop patches
     that upstream made obsolete.
   - Dependency ranges must match what was tested in 2.4.
   - Commits to third-party repos should use the maintainer's GitHub noreply email, not a personal address.

## 4. After release
- [ ] Update `benchmarking/FL_regression/results_history.tsv` if it wasn't already.
- [ ] Comment on, and close, the GitHub issues fixed in the release.
- [ ] Watch the Bioconda PR's CI and answer reviewers.
