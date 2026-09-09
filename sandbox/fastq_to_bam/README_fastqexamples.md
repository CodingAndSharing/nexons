# FASTQ → BAM for nexons — step by step

Turning Oxford Nanopore cDNA reads into the coordinate-sorted, indexed BAM that
`nexons.py` expects, on mouse GRCm39.

Follow the **RECOMMENDED** steps in order. Each one ends in something you can
check, so a failure is caught at the step that caused it rather than an hour
later. Every parameter of the alignment command itself is explained separately
in [`COMMAND_EXPLAINED.md`](COMMAND_EXPLAINED.md).

## What you end up with

```
sandbox/big_data/
├── input/                  ← your FASTQ reads go here
│   ├── PBM75668_pass_L001_barcode01.fastq.gz
│   └── PBM75668_pass_L001_barcode02.fastq.gz
├── reference/              ← genome, annotation, minimap2 index
│   ├── Mus_musculus.GRCm39.dna.primary_assembly.fa.gz
│   ├── Mus_musculus.GRCm39.116.gtf.gz
│   └── Mus_musculus_GRCm39_2026_04.mmi
├── outputs/                ← BAMs land here
│   ├── PBM75668_pass_L001_barcode01.bam
│   └── PBM75668_pass_L001_barcode01.bam.bai
└── tmp/                    ← samtools sort scratch
```

Nothing under `big_data/` is version-controlled — it is excluded by
`sandbox/.gitignore`, which is deliberate: a single barcode here is ~9 GB.

---

## Step 0 (RECOMMENDED) — check it will fit before downloading anything

**Do this first.** minimap2 holds the entire index in RAM for the whole run. If
it does not fit, the kernel kills the process with no message beyond `Killed`,
after you have already spent an hour on the download and the index build.

```bash
free -g              # total and available RAM
nproc                # cores
df -h .              # free disk on the data volume
```

This machine reports **~7.6 GiB total, ~5 GiB available, 16 cores**.

The arithmetic for a whole-genome mouse splice index:

| Quantity | Value | Where it comes from |
|---|---|---|
| Genome size | ~2.7 Gbase | GRCm39 primary assembly |
| Minimizer window (`-x splice`) | `w=5` | Denser than the `w=10` of genomic presets |
| Minimizers stored | ~900 M | ≈ 2 × 2.7e9 / (w+1) |
| Hash table | ~7 GiB | ~8 bytes per minimizer |
| Packed sequence | ~0.7 GiB | 2 bits per base |
| **Resident total** | **~8 GiB** | plus ~1 GiB working memory |

**That does not fit in 5 GiB available**, and the *build* peaks higher still.
So on this laptop as configured, pick one:

- **Raise the memory ceiling** — see [Appendix: more
  memory](#appendix-more-memory). This is WSL2, so the 7.6 GiB is an
  allocation, not the hardware. If the laptop physically has 16 or 32 GB this
  is a one-line config change and the whole-genome path opens up. **Check this
  before accepting the constraint.**
- **Use a subset reference** — index only the chromosomes your GTF covers.
  Covered at [Step 5b](#step-5b-alternative--subset-reference).
- **Align on a bigger machine.** The index cannot be built elsewhere and
  copied back: it has to be resident *during* alignment, so the machine that
  holds it must also do the aligning.

These estimates are arithmetic, not measurements — minimap2 has never been run
in this environment. Measure with `/usr/bin/time -v` on the index build and
trust that over the table above.

---

## Step 1 (RECOMMENDED) — install the tools

```bash
cd sandbox
pixi install
```

Then confirm all four are present:

```bash
pixi run minimap2 --version
pixi run samtools --version | head -1
pixi run seqkit version
pixi run pigz --version
```

If `minimap2` is missing after this, the manifest and the lock file have
drifted apart — check that `pixi.toml` lists it under `[dependencies]` and
re-run `pixi install`.

**Always invoke the pipeline through `pixi run`.** Running the script directly
picks up whatever `minimap2` happens to be on your `PATH`, which may be a
different version than the one pinned here, and version drift in an aligner
changes results.

---

## Step 2 (RECOMMENDED) — create the folder layout

```bash
cd sandbox/big_data
mkdir -p input outputs reference tmp
```

---

## Step 3 (RECOMMENDED) — set the environment variables

The script reads `$REPO/.env` — the repo root, *not* the sandbox. Start from
the template:

```bash
cd /path/to/nexons
cp .env.example .env
```

Then edit `.env` to this:

```bash
# --- data locations ---
FASTQ_DIR=$SANDBOX/big_data/input
BAM_DIR=$SANDBOX/big_data/outputs
TMPDIR=$SANDBOX/big_data/tmp

# --- reference ---
REF_MMI=$SANDBOX/big_data/reference/Mus_musculus_GRCm39_2026_04.mmi
REF_GTF=$SANDBOX/big_data/reference/Mus_musculus.GRCm39.116.gtf.gz
REF_FASTA=

# --- resources ---
THREADS=8
SORT_THREADS=4
SORT_MEM=512M
SECONDARY=no
```

Four things to get right here.

**`$SANDBOX` and `$REPO` are set by the script before it sources `.env`**, so
these paths stay correct no matter which directory you run from. A bare
relative path like `sandbox/big_data/input` is resolved against your *current*
directory and breaks the moment you `cd` anywhere.

**`REF_MMI` is the genome you align TO.** It is never your reads. Reads have no
variable — they are found in `FASTQ_DIR` or passed as arguments. The script
rejects a FASTQ in either reference variable by name, because the bare "does
not exist" error gives no hint as to the cause.

**Set either `REF_MMI` or `REF_FASTA`, not both.** With only `REF_FASTA` the
script builds the index once and reuses it. With `REF_MMI` it uses the index as
given. Leave the other empty.

**`SORT_MEM` is per thread.** `SORT_THREADS=4` with `SORT_MEM=512M` reserves
about 2 GiB in total, not 512 MiB. This is an easy way to request several times
more memory than you intended.

Check it loads:

```bash
cd sandbox
pixi run ./fastq_to_bam/fastq_to_bam.sh
```

It will stop at the first thing that is actually missing. The line
`[env] loaded /path/to/nexons/.env` confirms the file was found.

---

## Step 4 (RECOMMENDED) — download the reference

Mouse GRCm39, Ensembl release 116. Both URLs and sizes below were checked
against the server.

```bash
cd sandbox/big_data/reference

# Genome, 806 MB gzipped, ~2.7 Gbase
curl -L --fail -O \
  https://ftp.ensembl.org/pub/release-116/fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna.primary_assembly.fa.gz

# Annotation for nexons, 108 MB gzipped
curl -L --fail -O \
  https://ftp.ensembl.org/pub/release-116/gtf/mus_musculus/Mus_musculus.GRCm39.116.gtf.gz
```

Use **`primary_assembly`**, not `toplevel`. `toplevel` additionally contains
haplotype and patch scaffolds, which give many reads a second near-identical
place to map: mapping quality collapses toward zero, and what nexons counts
changes quietly rather than failing.

`--fail` matters. Without it curl writes an HTML error page into the `.fa.gz`
and you discover it during the index build.

Verify before going further:

```bash
ls -lh
gzip -t Mus_musculus.GRCm39.dna.primary_assembly.fa.gz && echo "genome OK"
gzip -t Mus_musculus.GRCm39.116.gtf.gz && echo "gtf OK"
```

Two consistency requirements. The **assembly must match** — a GRCm39 GTF
against a GRCm38/mm10 genome puts every coordinate in the wrong place without
erroring. And **chromosome naming must match**: Ensembl uses `1`, `19`, `X`
where UCSC uses `chr1`, `chr19`, `chrX`. Mixing the conventions is not an
error; it produces zero counts everywhere. Both files above are Ensembl, so
they agree with each other.

Do not decompress the genome. minimap2 reads gzip directly.

---

## Step 5 (RECOMMENDED) — build the minimap2 index

Only if you were not handed one. **If a collaborator gave you
`Mus_musculus_GRCm39_2026_04.mmi`, skip the build but check two things**: its
size against Step 0, and which preset built it.

```bash
cd sandbox
pixi run minimap2 -x splice -t 16 \
  -d big_data/reference/Mus_musculus_GRCm39_2026_04.mmi \
     big_data/reference/Mus_musculus.GRCm39.dna.primary_assembly.fa.gz
```

`-x splice` is load-bearing here, not just at alignment time: **`k` and `w` are
baked into the `.mmi`**. If the index was built with a different preset,
minimap2 silently uses the index's values and ignores the `-x splice` you pass
later, with only a warning that is easy to miss. An index built with
`-x map-ont` (`w=10`) gives quietly worse sensitivity on short exons.

Do not add `-I 1G`. That splits the index into batches and then requires
`--split-prefix` at alignment time. At the default `-I 4G` the mouse genome
indexes as a single batch.

Time it and watch the peak, because this is the step most likely to be killed:

```bash
/usr/bin/time -v pixi run minimap2 -x splice -t 16 -d ... 2>&1 | grep -E "Maximum resident|Elapsed"
```

### Step 5b (alternative) — subset reference

If Step 0 said the whole genome will not fit and you cannot raise the ceiling,
index only the chromosomes your GTF actually covers. For a defined gene panel
this costs nothing in accuracy for those genes.

```bash
cd sandbox/big_data/reference

# Example: chromosomes 2 and 11 only. Substitute your own.
for c in 2 11; do
  curl -L --fail -O \
    "https://ftp.ensembl.org/pub/release-116/fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna.chromosome.$c.fa.gz"
done
cat Mus_musculus.GRCm39.dna.chromosome.*.fa.gz > GRCm39_subset.fa.gz

cd ../..
pixi run minimap2 -x splice -t 16 \
  -d big_data/reference/Mus_musculus_GRCm39_2026_04.subset.mmi \
     big_data/reference/GRCm39_subset.fa.gz
```

Concatenating gzip members is valid gzip, so `cat` on `.fa.gz` files works.

Point `REF_MMI` at the `.subset.mmi` and **subset the GTF to the same
chromosomes**, or nexons will look for genes that are not in the BAM:

```bash
zcat Mus_musculus.GRCm39.116.gtf.gz \
  | awk '$1=="2" || $1=="11" || /^#/' \
  | gzip > GRCm39_subset.116.gtf.gz
```

Reads whose true origin is outside the subset will either fail to map or
mis-map onto the included chromosomes. That is acceptable for panel
quantification and **not** acceptable for anything transcriptome-wide.

---

## Step 6 (RECOMMENDED) — put the reads in `input/`

```bash
mv /wherever/PBM75668_pass_L001_barcode01.fastq.gz sandbox/big_data/input/
```

**The script refuses to start if any `*.part`, `*.filepart` or `*.crdownload`
file is present in `FASTQ_DIR`.** This is intentional. A partly-transferred
FASTQ is a valid-looking file that is simply missing its tail; minimap2 aligns
it and exits 0, so the truncation surfaces as missing counts in nexons rather
than as an error.

Verify each file is complete before aligning. A FASTQ record is exactly four
lines, so a complete file has a line count divisible by four:

```bash
cd sandbox/big_data/input
for f in *.fastq*; do
  awk 'END{printf "%s: %d reads, remainder %d -> %s\n", FILENAME, NR/4, NR%4, \
       (NR%4==0 ? "OK" : "TRUNCATED")}' "$f"
done
```

This reads the whole file, so it takes a minute or two per barcode. It is worth
it — it is the difference between a wrong answer and no answer.

### The two files currently in this repo

Both are still named `.part` and sat unchanged for hours, so the transfers are
dead rather than slow. Checked by record count:

| File | Reads | Lines mod 4 | Verdict |
|---|---|---|---|
| `barcode01` | 3,442,720 | 0 | Complete — usable |
| `barcode02` | 3,796,680 | **2** | **Truncated mid-record** — re-fetch |

To use barcode01:

```bash
cd sandbox/big_data
mv PBM75668_pass_L001_barcode01.2AM3fGSx.fastq.gz.part \
   input/PBM75668_pass_L001_barcode01.fastq
```

Note the extension. **Despite the `.gz` in the original name this file is
plain text, not gzip** — naming it `.fastq.gz` would make `zcat`, `gzip -t` and
`seqkit` fail on it. (minimap2 itself would cope; it sniffs the format.) Move
barcode02 out of `input/` until it has been re-fetched.

A record-aligned line count means the file ends on a record boundary. It is
strong evidence, not proof the run had no further reads — compare against the
expected yield from the sequencing report if you have it.

---

## Step 7 (RECOMMENDED) — run the pipeline

```bash
cd sandbox
pixi run ./fastq_to_bam/fastq_to_bam.sh
```

That aligns every FASTQ in `input/`. For one file:

```bash
pixi run ./fastq_to_bam/fastq_to_bam.sh big_data/input/PBM75668_pass_L001_barcode01.fastq
```

Per input, the script runs exactly this — the command from
[`COMMAND_EXPLAINED.md`](COMMAND_EXPLAINED.md), plus the three things that
command omits:

```bash
minimap2 -ax splice --secondary=no -t 8 \
         big_data/reference/Mus_musculus_GRCm39_2026_04.mmi \
         big_data/input/PBM75668_pass_L001_barcode01.fastq \
  | samtools sort -@ 4 -m 512M \
                  -T big_data/tmp/sort.PBM75668_pass_L001_barcode01 \
                  -o big_data/outputs/PBM75668_pass_L001_barcode01.bam -

samtools index -@ 4 big_data/outputs/PBM75668_pass_L001_barcode01.bam
samtools flagstat big_data/outputs/PBM75668_pass_L001_barcode01.bam
```

The additions are the `.bai` index (nexons queries by region and needs it), an
explicit `-T` sort scratch directory, and `-@`/`-m` instead of samtools' rather
conservative one-thread default.

Expect these on stdout:

```
[env]   loaded /path/to/nexons/.env
[ref]   .../Mus_musculus_GRCm39_2026_04.mmi (8.0 GiB)
[cpu]   minimap2 -t 8 | samtools sort -@ 4 -m 512M
[out]   .../big_data/outputs
[in]    1 file(s)
[align] PBM75668_pass_L001_barcode01
[done]  PBM75668_pass_L001_barcode01 in 47m12s -> .../barcode01.bam
```

A `WARNING: this may not fit in memory` line here is Step 0 catching you late;
stop and deal with it rather than hoping.

Re-running is safe. A BAM newer than its input is reported as `[skip]`, so an
interrupted batch resumes rather than redoing finished barcodes. `THREADS=8`
is deliberate on a 16-core box: it leaves cores for the concurrent sort.
`-t16` has both programs contending for every core and typically runs *slower*
end to end.

### Timing

No measured figure exists for this environment — minimap2 has never
successfully run here. For an estimate on your own hardware, use the archived
script's benchmark mode, which aligns *n* reads and extrapolates:

```bash
BENCHMARK=200000 pixi run ./fastq_to_bam/fastq_to_bam.sh.old
```

For scale: barcode01 is 3,442,720 reads / 4.17 Gbase, mean read length
1,211 bp, and its uncompressed FASTQ carries roughly 466 Mbase per GiB. Note
that transfer and alignment overlap — finished barcodes can be aligned while
the rest are still arriving.

---

## Step 8 (RECOMMENDED) — check the outputs

```bash
cd sandbox/big_data/outputs
ls -lh
samtools flagstat PBM75668_pass_L001_barcode01.bam
samtools idxstats PBM75668_pass_L001_barcode01.bam | head
```

What to look for:

| Check | Expected | If it fails |
|---|---|---|
| `.bam` **and** `.bam.bai` present | both | nexons cannot do region queries without the `.bai` |
| Mapped % in `flagstat` | high for cDNA on its own genome | A few percent means wrong assembly, wrong organism, or a truncated input |
| `idxstats` chromosome names | `1`, `2`, `19` … | If these are `chr1`-style, your GTF must match or counts come out zero |
| `secondary` in `flagstat` | 0 | Expected: `--secondary=no`. It means *not measured*, not *no ambiguity* |
| Total reads vs input | comparable | A large shortfall points at a truncated FASTQ |

The mapped-read count is the denominator nexons' totals should be read
against, so record it.

---

## Step 9 (RECOMMENDED) — hand it to nexons

```bash
cd /path/to/nexons
python nexons.py \
  sandbox/big_data/reference/Mus_musculus.GRCm39.116.gtf.gz \
  sandbox/big_data/outputs/*.bam \
  --direction none \
  --outbase sandbox/big_data/outputs/nexons_run1
```

`--direction none` is correct for this library: ONT cDNA is unstranded, so
requiring a strand match would discard roughly half the reads. Use the same
GTF you verified in Step 4 — and if you took the Step 5b subset path, the
*subset* GTF.

---

## Troubleshooting

| Symptom | Cause |
|---|---|
| `Killed`, no other output | Out of memory. Step 0. Usually the index |
| `ERROR: minimap2 is not on PATH` | Not running through `pixi run`, or `pixi install` incomplete |
| `ERROR: REF_MMI looks like sequencing reads` | A FASTQ in a reference variable. Reads go in `FASTQ_DIR` |
| `ERROR: ... does not exist` on a path in `.env` | Relative path resolved against your cwd. Use `$SANDBOX/...` |
| `ERROR: unfinished download(s)` | A `.part` file in `input/`. Finish or move it aside |
| `WARNING: TMPDIR ... is a tmpfs` | Sort spills would consume RAM. Point `TMPDIR` at the data volume |
| Everything maps but counts are zero | Chromosome naming mismatch between GTF and reference |
| `gzip: not in gzip format` | A plain FASTQ under a `.gz` name. Rename it |
| samtools parse error at the pipe | `-a` missing from minimap2, so it emitted PAF not SAM |
| BAM exists but nexons cannot read regions | No `.bai`. Run `samtools index` |

---

## Appendix: more memory

### Raise the WSL allocation

**Check this before accepting the 8 GB ceiling.** This is WSL2, and the
7.6 GiB that `free -h` reports is WSL's *allocation*, not the laptop's RAM — by
default WSL2 takes 50% of Windows' total, historically capped at 8 GB. If the
laptop physically has 16 or 32 GB, most of it is simply not being offered to
Linux, and the whole-genome index that does not fit above may fit after one
config change.

Check the physical total from PowerShell:

```powershell
(Get-CIMInstance Win32_ComputerSystem).TotalPhysicalMemory / 1GB
```

If that is well above 8, create or edit `C:\Users\<you>\.wslconfig`:

```ini
[wsl2]
memory=24GB          # leave Windows 8 GB or so; do not allocate everything
processors=16
swap=16GB            # cheap insurance against an OOM kill mid-index
```

Then, from PowerShell:

```powershell
wsl --shutdown
```

and reopen the shell. Confirm with `free -h`. With ~24 GB the GRCm39 splice
index builds comfortably and the subset path stops being necessary.

Two caveats. WSL's memory is reserved from Windows while in use, so
over-allocating makes Windows itself swap — leave it real headroom. And
`.wslconfig` lives on the Windows side, outside this repo, so it is not
version-controlled with the pipeline; note the setting in your methods if the
analysis depends on it.

### Why `TMPDIR` points away from `/tmp`

`samtools sort` holds `SORT_THREADS × SORT_MEM` in RAM and spills the excess to
`$TMPDIR`. On this machine `/tmp` is a **tmpfs** — RAM, not disk — so spilling
there both runs out of room and steals the memory minimap2 is using for the
index. `.env` therefore points `TMPDIR` at the data volume. A sort that dies on
a full `/tmp` throws away the alignment work that preceded it.

### A note on `JUNC_BED`

minimap2's `--junc-bed` rewards junctions present in an annotation. That
improves accuracy on known transcripts — and biases against unannotated ones. A
cryptic exon is by definition unannotated, so annotated-junction bonuses push
reads toward the canonical splice form and suppress the signal. **Leave it
unset for cryptic exon work.** It is wired up only in the archived script.

---

## Also available

[`fastq_to_bam_examples.ipynb`](fastq_to_bam_examples.ipynb) walks the same
pipeline stage by stage with the preflight checks, a read-length plot and a
throughput measurement, using the pixi environment as its kernel
(**Python (nexons pixi)**). It is cheap and subsampled by default, degrades
rather than fails when tools are missing, and is generated by
`build_fastq_notebook.py` — edit the builder, not the notebook. Clear outputs
before committing.

[`fastq_to_bam.sh.old`](fastq_to_bam.sh.old) is the previous script. It retains
`BENCHMARK=<n>` and `JUNC_BED`, which the current lean script drops.
