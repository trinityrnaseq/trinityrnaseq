# Full-length reconstruction regression benchmark

Before each Trinity release, we assemble the mouse and fission yeast (S. pombe) RNA-Seq data sets from the
original Trinity paper (Grabherr et al., Nat Biotech 2011) with both the previous release and the release
candidate. We then check that the number of full-length (FL) reconstructed reference transcripts is consistent.

The running record of results is in [`results_history.tsv`](results_history.tsv). Add rows for every release.

## Data

| Data set | Source (full ~50M PE reads) | Benchmark subset | Library |
|---|---|---|---|
| mouse | `gs://ctat-public-ftp-resources/trinity_test_rnaseq_datasets/mouse.all.{Left,Right}.fq.gz` | first 10M pairs | 76 bp PE, dUTP strand-specific (`--SS_lib_type RF`) |
| pombe (SP2) | `gs://ctat-public-ftp-resources/trinity_test_rnaseq_datasets/SP2.all.{Left,Right}.fq.bz2` | first 10M pairs | 76 bp PE, dUTP strand-specific (`--SS_lib_type RF`) |

The same bucket also holds `drosoph_50M` (fasta) and `BT474` reads, which are not part of the routine benchmark.

The reference transcript sets are pinned by `prep_data.sh`:

| Reference | File | Entries |
|---|---|---|
| mouse CDS (primary metric) | GENCODE vM39 `pc_transcripts`; the CDS span of each transcript, with identical CDSs collapsed per gene | 185,119 CDS / 22,311 genes |
| mouse cDNA | GENCODE vM39 `transcripts`, transcript_type == protein_coding | 180,384 |
| pombe CDS (primary metric) | PomBase `cds.fa.gz`, protein-coding only | 5,362 |
| pombe cDNA | Ensembl Fungi release-63 ASM294v2 `cdna.all` | 5,175 |

These are the md5 sums of the downloads used for the v2.16.0 benchmark (2026-11). PomBase does not version
`cds.fa.gz`, so check this sum before comparing against older rows.

```
5955d3eb3f579122d9df32b13081b919  gencode.vM39.pc_transcripts.fa.gz
9926a56585ecc29f63b61ad425d6db8f  gencode.vM39.transcripts.fa.gz
bda2d2f4f5d985dc24899afa8de8099a  pombe.cdna.all.fa.gz
b125ac3121f73aaf1ab78d51667c70e7  pombe.cds.fa.gz
```

The CDS-based counts are the metric to track. The cDNA-based counts are low by design: an assembly has to
cover at least 99% of the reference *including annotated UTRs*, and current annotations carry long UTRs and
many long-read isoforms. The cDNA counts are recorded anyway as a secondary signal.

## Procedure

```bash
B=$TRINITY_HOME/benchmarking/FL_regression
D=/path/to/bench_data                      # ~8 GB reads + ~1 GB refs

# 1. data + references (once; reuse across releases)
$B/prep_data.sh $D

# 2. assemble with the previous release and the release candidate (run both under the same conditions)
$B/run_trinity.sh <old Trinity> pombe $D runs/old_pombe
$B/run_trinity.sh <new Trinity> pombe $D runs/new_pombe
$B/run_trinity.sh <old Trinity> mouse $D runs/old_mouse
$B/run_trinity.sh <new Trinity> mouse $D runs/new_mouse

# 3. FL evaluation: use ONE evaluator install (TRINITY_HOME) for all assemblies; blat must be on PATH
for v in old new; do
  $B/fl_eval.sh $TRINITY_HOME runs/${v}_pombe/trinity_out_dir.Trinity.fasta $D/ref/pombe.cds.ref.fa  eval/${v}_pombe_cds  4
  $B/fl_eval.sh $TRINITY_HOME runs/${v}_pombe/trinity_out_dir.Trinity.fasta $D/ref/pombe.cdna.ref.fa eval/${v}_pombe_cdna 4
  $B/fl_eval.sh $TRINITY_HOME runs/${v}_mouse/trinity_out_dir.Trinity.fasta $D/ref/mouse.cds.ref.fa  eval/${v}_mouse_cds  6
  $B/fl_eval.sh $TRINITY_HOME runs/${v}_mouse/trinity_out_dir.Trinity.fasta $D/ref/mouse.cdna.ref.fa eval/${v}_mouse_cdna 6
done
```

Each `fl_eval.sh` prints the `FL.pslx.maps.summary` line, which has these columns: label, FL genes, FL transcripts, FL genes
(multi-isoform genes), FL transcripts (multi-isoform genes). `TrinityStats.txt` has the transcript and gene counts and N50.

Settings used: `--SS_lib_type RF --CPU 6 --max_memory 20G`, with Trinity's default normalization and salmon filtering.
On a shared 28-core host, pombe takes about 2.3 h and mouse about 5 h per assembly. Peak memory is under 20 GB.
Each BLAT chunk against the mouse references uses a few GB.

### Running the versions in containers (no Docker / no UGER)

At the Broad, we ran the comparison locally with unprivileged Apptainer:
- **Previous release:** an Apptainer sandbox of the published image (`apptainer build --sandbox trin_old docker://trinityrnaseq/trinityrnaseq:<old>`), using `/usr/local/bin/Trinity`.
- **Release candidate:** a sandbox built from the candidate `Docker/Dockerfile`, or a `make && make plugins` build of the candidate commit inside an Ubuntu 24.04 sandbox.

Wrap the commands like this:

```bash
export APPTAINERENV_TMPDIR=/tmp
apptainer exec --cleanenv --env LC_ALL=C -B $PWD:$PWD <sandbox> \
    bash -c "cd $PWD && $B/run_trinity.sh /usr/local/bin/Trinity pombe $D runs/old_pombe"
```

`--cleanenv` and `LC_ALL=C` keep the host's locale and TMPDIR settings out of the container.

## Interpreting results

Trinity is not bit-for-bit reproducible across runs: thread scheduling and normalization tie-breaking
cause small differences. For v2.15.2 vs v2.16.0, the CDS-based FL counts differed by 0.2% or less, and
the total transcript counts by a few dozen. We haven't yet measured same-version replicate noise; doing so
(two runs of the same version) would firm up the threshold below.

As a rule of thumb, a change of more than about 1% in CDS-based FL genes on either data set warrants investigation before release. Bisect by
swapping one component (Inchworm, Chrysalis, Butterfly, or the Trinity driver) at a time.

When a reference or setting changes, rerun the previous release under the new conditions. Never compare
counts across rows that used different references or settings.
