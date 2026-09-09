# FASTQ → BAM for nexons

Turning Oxford Nanopore cDNA reads into the coordinate-sorted, indexed BAM that
`nexons.py` expects. Everything here runs inside the sandbox pixi environment.

The pipeline is one script, [`fastq_to_bam.sh`](fastq_to_bam.sh); this file
explains what it does, what to set before running it, and how long to expect it
to take on a 16-core / 8 GB laptop.

## Or run it as a notebook

[`fastq_to_bam_examples.ipynb`](fastq_to_bam_examples.ipynb) walks the same
pipeline stage by stage, with the preflight checks, the read-length and
alignment-QC plots, and a throughput measurement that extrapolates to your full
run. Open it with the **Python (nexons pixi)** kernel:

```bash
cd sandbox
pixi install                      # once
pixi run jupyter lab fastq_to_bam/fastq_to_bam_examples.ipynb
```

It is cheap by default — a 50,000-read subsample, about a minute end to end —
and every long stage is gated behind an explicit flag, so nothing starts a
multi-hour job unless you ask. It also degrades rather than failing if
`minimap2` is not installed yet, so it is safe to run before `pixi install`
finishes.

The notebook is a **build product** of
[`build_fastq_notebook.py`](build_fastq_notebook.py) — prose and code stay in
one reviewable diff, and a notebook broken by a bad merge can be regenerated:

```bash
pixi run python fastq_to_bam/build_fastq_notebook.py
```

Edit the builder, not the `.ipynb`. And **clear the outputs before committing**
(`Kernel → Restart Kernel and Clear All Outputs`) — executed cells store their
output inside the file.

---

## The dataset this was written against

This is a **multiplexed run** — `barcode01` and `barcode02` were both present
in `big_data`, so budget per barcode and multiply by however many arrive.

Measured across the whole of `barcode01` (8.9 GiB, download complete):

| Property | Value |
|---|---|
| Reads | 3,442,720 |
| Total sequence | 4.17 Gbase |
| Mean read length | 1,211 bp |
| Read length range | 13 – 10,460 bp |
| Record size on disk (uncompressed) | ~2,600 bytes/read |
| Yield per GiB of uncompressed FASTQ | ~387,000 reads ≈ 469 Mbase |
| Basecall model (from the `RG` tag) | `dna_r10.4.1_e8.2_400bps_hac@v5.2.0` |
| Library | `SMART-Seq`, one file per barcode |

Two consequences worth knowing before you pick flags:

- **This is cDNA, not direct RNA.** The `dna_` model prefix and the SMART-Seq
  kit both say so. cDNA reads come off both strands, so the aligner must look
  for `GT..AG` on either strand. `-x splice` already sets `-ub` (both strands),
  which is what you want. Do **not** add `-uf` — that is the direct-RNA setting
  and it will lose roughly half your junctions here.
- **The library is unstranded**, so nexons should run with `--direction none`
  (its default). Incidentally, that option's help text advertises `opposing`
  while the code tests for `opposite`; passing `opposing` raises. Not a problem
  for this library, but worth knowing.

### Two things about the file as downloaded

**It is not gzipped**, despite the `.fastq.gz` name — the first bytes are plain
ASCII `@<uuid>`, and `gzip -t` rejects it. Every tool here reads it fine
(minimap2 sniffs the magic bytes rather than trusting the extension), but it is
costing you about 3× the disk it needs to. Once the transfer is finished:

```bash
cd sandbox/big_data
# One file per barcode; give each an honest name, then compress in parallel.
# The run prefix is taken from the filename, so nothing is hardcoded.
for f in *_pass_L001_barcode*.fastq.gz; do
    bc=${f#*_barcode}; bc=barcode${bc%%.*}     # e.g. barcode01
    mv "$f" "${f%%_*}_${bc}.fastq"             # e.g. RUNID_barcode01.fastq
done
pixi run pigz -p 8 *_barcode*.fastq            # → .fastq.gz, ~3x smaller
```

Compressing is optional — minimap2 reads either — but at ~9 GiB per barcode it
pays for itself quickly on a multiplexed run.

**The download was still running** when this was written: `barcode01` had
finished at 8.9 GiB, `barcode02` was still arriving at roughly 9–11 MiB/s, and
both still carried the `.part` suffix that the transfer tool strips on
completion. `fastq_to_bam.sh` refuses to start while any `.part` file is
present in the input directory, because minimap2 will align a truncated FASTQ
without complaining and exit 0, and you will not notice until the counts look
wrong.

That guard is deliberately blunt — it blocks on *any* partial file, not just
the one you asked for. If you want to start aligning finished barcodes while
the rest download, move them into a separate directory and point `FASTQ_DIR`
at that.

---

## Step 0 — install the tools

`minimap2`, `samtools`, `seqkit` and `pigz` are declared in
`sandbox/pixi.toml`. They are **not yet installed** — the lockfile has not been
regenerated, so run this once:

```bash
cd sandbox
pixi install
```

This also picks up the `plotly` / `kaleido` / `streamlit` pins added for
`viz_scripts`. Then prefix commands with `pixi run`, or `pixi shell` once and
drop the prefix.

---

## Step 1 — the `.env` file

The script reads configuration from `.env` at the **repo root** via
`set -a; source .env; set +a`, so every variable it defines becomes available
to the tools without ever being echoed.

I could not read your `.env` — it is masked from my environment — so the
variable names below are what the script *expects*. Reconcile them with what
you already have; anything absent falls back to the default in the right-hand
column.

| Variable | Purpose | Default |
|---|---|---|
| `FASTQ_DIR` | where the reads are | `sandbox/big_data` |
| `BAM_DIR` | where BAMs go | `sandbox/big_data/bam` |
| `REF_FASTA` | reference genome FASTA | *(required)* |
| `REF_MMI` | prebuilt minimap2 index; skips the build | derived from `REF_FASTA` |
| `REF_GTF` | annotation, passed on to nexons | *(unset)* |
| `THREADS` | alignment threads | `nproc` |
| `SORT_THREADS` / `SORT_MEM` | `samtools sort` resources | `4` / `512M` |
| `TMPDIR` | sort spill directory | `sandbox/big_data/tmp` |
| `SECONDARY` | emit secondary alignments (`yes`/`no`) | `no` |
| `JUNC_BED` | annotated junctions to reward | *(unset — see below)* |
| `FTP_HOST` / `FTP_USER` / `FTP_PASS` | download credentials | *(unset)* |

### Why `TMPDIR` points away from `/tmp`

`/tmp` here is a **3.9 GiB tmpfs** — RAM, not disk. `samtools sort` spills
roughly the size of the final BAM, several GiB for one barcode, so sending it
to `/tmp` fails two ways at once: it runs out of room, *and* every byte it does
write occupies memory that minimap2 needs for its index. A sort that dies
part-way throws away the alignment work that preceded it.

So `TMPDIR` defaults to `sandbox/big_data/tmp` on the ext4 volume, which has
**689 GiB free**. **You do not need to enlarge `/tmp` for this pipeline** — it
is never used.

One way that default can be bypassed: `TMPDIR="${TMPDIR:-...}"` honours an
existing value, and many shells and job schedulers export `TMPDIR` already. The
script therefore checks at startup and, before aligning anything:

- refuses to start if `TMPDIR` or `BAM_DIR` has less free space than the run
  needs (budgeting half the total input size on each, since a BAM runs
  0.35–0.5× its uncompressed FASTQ);
- warns if `TMPDIR` resolves onto a tmpfs, naming the RAM cost.

Both print before the first read is aligned, so a space problem costs you
seconds rather than the better part of an hour. Check what it resolved to in
the `[tmp]` line of the startup banner.

### If you still want a bigger `/tmp`

Confirm what you actually have first — the numbers above are what I measured,
and your shell may see something different:

```bash
df -h /tmp && findmnt /tmp
```

If `FSTYPE` is `tmpfs`, resizing costs RAM, which on this machine is the scarce
resource rather than disk. Temporarily, until reboot:

```bash
sudo mount -o remount,size=8G /tmp
```

Persistently, add to `/etc/fstab` (WSL honours it — `mountFsTab` defaults on):

```
tmpfs /tmp tmpfs rw,nosuid,nodev,size=8G 0 0
```

Better on a memory-limited box: make `/tmp` disk-backed instead of larger, so
it costs no RAM at all:

```bash
sudo mkdir -p /var/tmpdisk && sudo chmod 1777 /var/tmpdisk
sudo mount --bind /var/tmpdisk /tmp
```

But for this pipeline all three are unnecessary — setting `TMPDIR` in `.env` is
the same fix without touching system mounts.

### Credentials

Keep them out of the shell and out of git. `.env` is already gitignored at the
repo root. A password passed on a command line is visible to every user on the
machine via `ps`, so prefer a `~/.netrc`:

```
machine <ftp-host>
  login    <username>
  password <password>
```

```bash
chmod 600 ~/.netrc
```

The real host, user and password go in `.env` (`FTP_HOST` / `FTP_USER` /
`FTP_PASS`) or `~/.netrc`, never in a file under version control. This repo is
a fork of a public one, so assume anything committed here is world-readable:
`.env.example` documents the variable names, `.env` holds the values and is
gitignored.

> Two habits worth keeping. A password passed as a command-line argument is
> visible in `ps` to every user on the machine, which is why the commands below
> use `--netrc` rather than `--user`. And if a credential is ever pasted
> somewhere it should not be — a chat window, a ticket, a shared notebook —
> treat it as compromised and have it rotated rather than hoping it went
> unnoticed.

---

## Step 2 — download

Already done in your case; recorded for reproducibility. With `~/.netrc` in
place, no secret appears on the command line:

```bash
mkdir -p sandbox/big_data && cd sandbox/big_data
lftp -e "mirror --verbose --parallel=4 $REMOTE_DIR ; quit" "$FTP_HOST"
# or, per-file:
curl --netrc -O "ftp://$FTP_HOST/$REMOTE_DIR/<run>_pass_L001_barcode01.fastq.gz"
```

`lftp mirror` is worth it over `curl` for a multi-GB transfer: it resumes
(`mirror --continue`) and parallelises across files.

---

## Step 3 — verify the transfer

Do this before spending an hour of CPU on it:

```bash
cd sandbox/big_data
ls -la *.part                     # must return nothing
pixi run seqkit stats -a *.fastq  # read count, N50, mean length
```

If the file really is gzipped, `gzip -t file.fastq.gz` should be silent; a
truncated gzip reports `unexpected end of file`. For a plain FASTQ, a record
count divisible by four is the equivalent check:

```bash
awk 'END{print NR, NR%4}' *_barcode01.fastq   # second number must be 0
```

---

## Step 4 — the reference, and the memory problem

**This is the step that decides whether the job runs on this laptop at all.**

The `splice` preset uses `-k15 -w5`. That `w=5` window is half of `map-ont`'s
`w=10`, so the minimizer table holds roughly twice as many entries, and a
whole-GRCh38 splice index runs to well over ten gigabytes. This machine has
**7.6 GiB of RAM total and ~5.5 GiB free**. A whole-genome `-x splice` index
build will be OOM-killed, usually with no message beyond the shell reporting
`Killed`.

The numbers below are estimates from how minimap2's index scales, not
measurements — I had no minimap2 binary available to measure them. Confirm the
peak yourself with:

```bash
/usr/bin/time -v pixi run minimap2 -x splice -t 8 -d ref.splice.mmi ref.fa \
  2>&1 | grep 'Maximum resident'
```

You have three ways forward.

### Option A — targeted reference (recommended here)

nexons quantifies the genes in a GTF. If that GTF covers one locus or a handful
of genes, indexing 3.1 Gb of genome to align against 60 Mb of it is wasted work
and wasted RAM. Extract the relevant chromosome(s):

```bash
cd sandbox/big_data
# GRCh38 chromosome 19 (UNC13A is at 19p13.11) -- 16.7 MB gzipped, ~59 Mb of sequence
curl -O https://ftp.ensembl.org/pub/current_fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.19.fa.gz
# Ensembl ships plain gzip, not bgzip, so `samtools faidx` cannot read it
# directly -- it needs either an uncompressed or a bgzipped FASTA.
gunzip Homo_sapiens.GRCh38.dna.chromosome.19.fa.gz
pixi run samtools faidx Homo_sapiens.GRCh38.dna.chromosome.19.fa
```

Index build drops to seconds and well under 1 GiB; alignment is several times
faster because there are far fewer seed hits to extend.

The cost is real and you should state it in any methods section: reads whose
true origin is elsewhere in the genome have nowhere else to go, so they either
fail to map or get forced onto the subset. For a paralogous gene family that
inflates apparent expression. Subset only when the analysis is genuinely
targeted, and make sure the GTF you give nexons is subset to match.

### Option B — whole genome, split index

`-I` caps how many reference bases minimap2 loads at once, bounding memory:

```bash
pixi run minimap2 -ax splice -I 1G -t 8 GRCh38.fa reads.fastq > out.sam
```

Memory stays inside the cap, but minimap2 makes one pass over the reads **per
part** — four passes for GRCh38 at `-I 1G` — so wall time multiplies roughly
fourfold. The more serious problem is downstream: each part is aligned
independently, so mapping quality and the primary/secondary flags are only
meaningful *within* a part. A read with hits in two parts emerges with a
primary alignment in each. nexons counts primary alignments, so that
double-counts. If you go this route you must collapse to one alignment per read
first, which is enough extra machinery that Option A or C is usually the better
trade.

### Option C — whole genome, bigger machine

The clean answer if you need genome-wide alignment with trustworthy MAPQ.
Aligning against a prebuilt index still needs the whole index resident, so
building the `.mmi` elsewhere and copying it back does **not** help — the
bigger machine has to do the alignment too. Everything below therefore runs on
that machine, not on the laptop.

Nothing about the pipeline changes. You point `REF_FASTA` at the whole genome
instead of chromosome 19 and run the same script; it builds the index on first
invocation and reuses it afterwards.

**1. Get the reference.** Use `primary_assembly`, *not* `toplevel`:

```bash
cd sandbox/big_data
curl -L --fail -O https://ftp.ensembl.org/pub/current_fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
```

`toplevel` also contains haplotype and patch scaffolds. They inflate the index
and, worse, give reads a second near-identical place to map, which collapses
MAPQ and quietly changes what nexons counts. `primary_assembly` is the one you
want. minimap2 reads gzipped FASTA directly, so leave it compressed.

**2. Build the index** (the memory-hungry step, ~10–20 min):

```bash
pixi run minimap2 -x splice -t 16 \
  -d Homo_sapiens.GRCh38.dna.primary_assembly.splice.mmi \
     Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
```

`-x splice` matters: it sets `k=15, w=5`, and those values are baked into the
`.mmi`. An index built under a different preset is silently used with *its*
k/w, not the one you ask for at alignment time. Do **not** add `-I 1G` here —
that splits the index and forces `--split-prefix` at alignment time, which is
Option B. At the default `-I 4G` the 3.1 Gbase genome indexes as a single
batch, which is what you want.

You can skip this step and let the script do it — it runs exactly this command
when `REF_MMI` is unset. Doing it by hand just separates the one failure you
care about from a four-hour pipeline run.

**3. Run the pipeline:**

```bash
cd sandbox
REF_MMI=big_data/Homo_sapiens.GRCh38.dna.primary_assembly.splice.mmi \
FASTQ_DIR=big_data \
BAM_DIR=big_data/bam \
THREADS=16 SORT_THREADS=4 SORT_MEM=1G \
TMPDIR=big_data/tmp \
  ./fastq_to_bam/fastq_to_bam.sh
```

Set these in `.env` instead if you prefer; the script reads them from there.

**Budget the RAM, don't just scale it.** The index stays resident for the whole
alignment, and `samtools sort` runs *concurrently* with it off the same pipe,
using roughly `SORT_THREADS × SORT_MEM`. So more cores does not mean you should
raise the sort allocation:

| Machine | `THREADS` | `SORT_THREADS` × `SORT_MEM` | Rough peak |
|---|---|---|---|
| 32 GB / 16 core | 16 | 4 × 768M (≈3 GiB) | ~22 GiB |
| 64 GB / 32 core | 32 | 8 × 2G (≈16 GiB) | ~36 GiB |

On a 32 GB machine keep the sort small — the index is the tenant that cannot be
evicted. `TMPDIR` must be on a real disk with ~2× the BAM size free, and it
must not be a tmpfs; the script checks both and aborts up front.

The index sizes and peaks above are **estimates, not measurements** — I had no
minimap2 binary to measure them. minimap2's own figure for a default human
index is ~11 GiB peak; `-x splice` uses `w=5` instead of `w=10`, which roughly
doubles the minimizer count, hence ≥32 GB rather than 16. Confirm with
`/usr/bin/time -v` on the index build before planning a long run around it.

**4. Sanity-check that it was worth it.** The point of the whole genome is
reads that chr19 had nowhere to put, so measure that:

```bash
pixi run samtools idxstats big_data/bam/<sample>.bam \
  | awk '{m+=$3} END {printf "mapped: %d\n", m}'
pixi run samtools view -c -f 4 big_data/bam/<sample>.bam   # still unmapped
```

If the mapped count is within a percent or two of your chr19 run for the genes
in your GTF, Option A was sufficient and you can go back to it.

### Option D — raise the WSL memory limit

**Check this before accepting the 8 GB ceiling.** This is WSL2
(`6.18.33.2-microsoft-standard-WSL2`), and the 7.6 GiB that `free -h` reports
is WSL's *allocation*, not the laptop's RAM — by default WSL2 takes 50% of
Windows' total, historically capped at 8 GB. If the laptop physically has 16 or
32 GB, most of it is simply not being offered to Linux, and the whole-genome
index that "exceeds this machine's RAM" above may fit after one config change.

Check the physical total from PowerShell:

```powershell
(Get-CIMInstance Win32_ComputerSystem).TotalPhysicalMemory / 1GB
```

If that comes back well above 8, create or edit `C:\Users\<you>\.wslconfig`:

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

and reopen the shell. Confirm with `free -h` inside WSL. With ~24 GB the
GRCh38 splice index builds comfortably and Option A stops being a constraint.

Two caveats. WSL's memory is reserved from Windows while in use, so
over-allocating makes Windows itself swap — leave it real headroom. And the
`.wslconfig` file lives on the Windows side, outside this repo, so it is not
version-controlled with the pipeline; note the setting in your methods if the
analysis depends on it.

### A note on `JUNC_BED`

minimap2's `--junc-bed` rewards junctions present in an annotation. That
improves accuracy on known transcripts — and biases against exactly the thing
this project is looking for. A cryptic exon is by definition unannotated, so
annotated-junction bonuses push reads toward the canonical splice form and
suppress the signal. **Leave `JUNC_BED` unset for cryptic exon work.** It is
wired up for the case where you are quantifying known isoforms instead.

---

## Step 5 — run it

```bash
cd sandbox/fastq_to_bam
pixi run ./fastq_to_bam.sh                    # everything in $FASTQ_DIR
pixi run ./fastq_to_bam.sh ../big_data/one.fastq   # a single file
```

The alignment is piped straight into `samtools sort`, so no intermediate SAM is
ever written — that would cost roughly twice the FASTQ size in disk for nothing.
Per input the script writes `$BAM_DIR/<sample>.bam` plus `.bai`, then prints
`flagstat`. Re-running skips any sample whose `.bai` is newer than its FASTQ, so
an interrupted batch resumes.

**Measure before committing to the full run:**

```bash
BENCHMARK=200000 pixi run ./fastq_to_bam.sh
```

That aligns 200k reads, reports Mbase/s and Gbase/hour for your actual
reference and thread count, and converts it to minutes per GiB of FASTQ. Two
minutes of measurement beats any estimate in this file.

---

## Step 6 — check the BAM

```bash
pixi run samtools flagstat sample.bam
pixi run samtools view -c -F 0x904 sample.bam    # primary mapped reads
pixi run samtools view sample.bam | awk '$6 ~ /N/' | head   # spliced (N in CIGAR)
```

Three things to look at:

- **Mapped fraction.** ONT cDNA against a matched whole genome should be well
  above 90%. A low fraction against a subset reference is expected, not a bug.
- **Reads with `N` in the CIGAR.** `N` is the spliced gap. If essentially no
  read has one, the index was built without `-x splice` and nexons will see an
  unspliced pile of reads. This is the most common silent failure.
- **`samtools quickcheck sample.bam`** — exits non-zero on a truncated BAM,
  which is what you get if the sort ran out of disk.

---

## Step 7 — hand it to nexons

```bash
cd ../..
python nexons.py annotation.gtf sandbox/big_data/bam/*.bam \
    --direction none \
    --outbase sandbox/big_data/nexons_out \
    --annotated-bam sandbox/big_data/annotated
```

`--direction none` for this unstranded cDNA library. `--annotated-bam` writes
the per-read classification tags that `sandbox/viz_scripts/nexons_viz.py`
renders — see `viz_scripts/README.md`.

The GTF must use the same contig names as the reference FASTA. Ensembl FASTA
headers are `19`; UCSC are `chr19`. Mixing them yields a BAM that indexes fine
and a nexons run that finds nothing at all.

---

## Time and memory

The input size is measured; the throughputs are not. `-ax splice` is taken at
1.5–4 Mbase/s across 16 threads against a whole genome, and roughly 6–15
Mbase/s against a single chromosome, because there are far fewer seed hits to
extend. I had no minimap2 binary available to measure this, so **treat the
times as order-of-magnitude** and get the real number from
`BENCHMARK=200000` before planning around them.

Per stage, for 16 cores / 8 GB RAM:

| Stage | Targeted reference (chr19) | Whole GRCh38 |
|---|---|---|
| Index build | ~10 s, <1 GiB RAM | **exceeds this machine's RAM** |
| Peak RAM while aligning | ~1–2 GiB | >10 GiB for the index alone |
| Alignment, per Gbase of reads | ~1–3 min | ~4–11 min |
| `samtools sort` + index | ~1–2 min per GiB of BAM | same |

### For one barcode

`barcode01` is 3,442,720 reads / **4.17 Gbase** / 8.9 GiB on disk:

| | chr19 | Whole GRCh38 |
|---|---|---|
| Alignment | 5–12 min | 17–46 min |
| Sort + index | 5–10 min | 5–10 min |
| **Per barcode** | **~12–20 min** | ~25–55 min |

### For the run

This is multiplexed, so multiply. Two barcodes were downloading when this was
written and the total is unknown; the arithmetic is linear in barcode count:

| Barcodes | FASTQ total | chr19 | Whole genome (≥32 GB machine) |
|---|---|---|---|
| 2 | ~18 GiB | 25–40 min | 50 min – 2 h |
| 6 | ~53 GiB | 1.2–2 h | 2.5–5.5 h |
| 12 | ~107 GiB | 2.5–4 h | 5–11 h |

Add the transfer: at the observed 9–11 MiB/s each 8.9 GiB barcode takes about
15 minutes, so on a 12-barcode run the download alone is ~3 hours and is
comparable to the whole alignment cost on chr19. The two overlap, though — you
can start aligning finished barcodes while the rest download, and the script's
`.part` guard plus its skip-if-newer logic make repeated invocation safe.

Disk, so you can plan: **689 GiB free** on the ext4 volume (it was 708 GiB
before the first two barcodes landed — each one costs ~9 GiB), which is ample
even for a 12-barcode run. Expect each BAM at roughly 0.35–0.5× its
uncompressed FASTQ (~3–4 GiB per barcode), sort temporaries about the same
again transiently, and a whole-genome splice `.mmi` over 10 GiB if you build
one. Note that `/tmp`'s 3.9 GiB does not enter into this: it is a RAM-backed
tmpfs and the pipeline does not use it.

**Recommendation: check Option D first, then default to Option A.** If the
laptop physically has more than 8 GB, raising the WSL memory limit is a
one-line config change that removes the constraint entirely and costs nothing
scientifically. Failing that, subset the reference to the chromosomes your GTF
covers (Option A) and the whole thing finishes inside half an hour with RAM to
spare — at the cost, stated in your methods, that reads originating elsewhere
in the genome have nowhere else to go.

---

## Troubleshooting

| Symptom | Cause |
|---|---|
| `Killed` during index build, no other output | OOM. Reference too big for 8 GB — Option A or `-I 1G` |
| Script refuses to start, names a `.part` file | Download unfinished. Wait for it |
| BAM has no `N` in any CIGAR | Index built without `-x splice`. Delete the `.mmi` and rebuild |
| `samtools sort` fails late with a disk error | Sort spilled into a full filesystem. Check the `[tmp]` startup line; set `TMPDIR` onto the big volume in `.env` |
| Script aborts naming free space before aligning | The disk preflight. `TMPDIR` or `BAM_DIR` is too small for ~half the input size |
| `WARNING: ... is on tmpfs` | `TMPDIR` was inherited from the shell onto RAM-backed storage. Unset it or point it at disk |
| `Killed` with plenty of disk free | RAM, not disk. See Option D (WSL memory limit) then Option A |
| nexons reports almost no hits | Contig naming mismatch (`19` vs `chr19`) between GTF and FASTA |
| Mapped fraction far below 90% | Expected against a subset reference; otherwise check the basecall model matches the chemistry |
| nexons' `secondary` bucket is always 0 | `SECONDARY=no` (the default) suppressed them. Set `SECONDARY=yes` |
| `gzip: not in gzip format` | The file is plain FASTQ under a `.gz` name. Harmless; rename it |
