# The alignment command, parameter by parameter

The reference run this pipeline is built around:

```bash
minimap2 -ax splice --secondary=no -t8 \
    Mus_musculus_GRCm39_2026_04.mmi \
    PBM75668_pass_L001_barcode01.fastq.gz \
  | samtools sort -o out.bam -
```

It is two programs joined by a pipe. `minimap2` decides where each read came
from in the mouse genome and writes one text line per alignment; `samtools
sort` reorders those lines by genomic coordinate and packs them into a binary
BAM. Nothing else happens — no filtering, no trimming, no quantification.

`sandbox/fastq_to_bam/fastq_to_bam.sh` runs exactly this command, once per
input FASTQ, and adds the three things it leaves out (see
[What the command omits](#what-the-command-omits)).

---

## The aligner

### `minimap2`

A pairwise aligner built for long, error-prone reads. It works in three
stages: reduce the reference to *minimizers* (short k-mers chosen by a
deterministic rule, so only a fraction of positions are stored), *chain*
minimizer matches shared between a read and a locus into a rough collinear
path, then *align* the gaps between chained anchors at base resolution. The
chaining stage is what lets it tolerate 5–10% error and place a read spanning
tens of kilobases.

### `-a`

Output **SAM** rather than minimap2's default **PAF**.

PAF is one compact line per alignment giving only the coordinate span. It has
no CIGAR string, so it does not record where the insertions, deletions and
introns fall inside the read. `samtools` cannot read it, and nexons needs
exactly the information it omits — the exon/intron structure of each read.

If you drop the `-a`, the pipeline fails immediately with a samtools parse
error rather than producing anything wrong. Cheap mistake to spot.

### `-x splice`

Load the **splice preset**. This is the single most consequential parameter.

`-x` selects a bundle of ~15 tuned settings appropriate to a read type. The
`splice` bundle differs from ordinary long-read presets in one structural way:
it permits very long deletions in an alignment and scores them as **introns**
rather than penalising them as sequencing errors. Without it, a cDNA read
crossing a 10 kb intron cannot be aligned as one read — it gets broken into
separate fragments or soft-clipped, and every splice junction in your data
disappears.

The settings that matter downstream:

| Setting | Value in `splice` | Why |
|---|---|---|
| `-k` | 15 | Minimizer k-mer length |
| `-w` | 5 | Minimizer window — denser than the 10 used for genomic reads, so short exons still get anchors |
| `-G` | 200k | Maximum intron length; alignments needing a longer gap are split |
| non-canonical splice penalty | on | Junctions at `GT..AG` are preferred over arbitrary gaps |

Run `minimap2 --help` for the full expansion of the preset on your build.

Two practical notes. First, `-a` and `-x splice` are conventionally written
together as `-ax splice` — identical to passing them separately. Second,
order matters: a preset assigns values as it is parsed, so any flag intended
to override a preset default must appear **after** `-x splice`, not before.

For an unstranded cDNA library like this one, the preset's default handling of
splice-site orientation is correct. Direct-RNA libraries should add `-uf` to
restrict junction detection to the transcript strand; doing that on cDNA would
discard roughly half your junctions.

### `--secondary=no`

Emit only one alignment per read, suppressing *secondary* alignments.

By default minimap2 reports up to five additional placements for reads that
map about as well in more than one locus — paralogues, repeats, gene families.
Suppressing them has three effects:

- The BAM gets substantially smaller, and sorting gets faster.
- Downstream counting cannot double-count one read against two genes.
- **The multi-mapping rate becomes invisible.** nexons keeps a `secondary`
  bucket in its output; with this flag that bucket reads zero, which looks
  like "no ambiguity" rather than "not measured".

For a defined gene panel this is the right default. If a client asks how
ambiguous the mapping was, re-run one barcode with `--secondary=yes` and read
the bucket, rather than inferring it from the suppressed run.

Note this is separate from *supplementary* alignments (chimeric or
fusion-spanning reads, SAM flag 0x800), which `--secondary=no` does **not**
remove.

### `-t8`

Use 8 threads for the alignment stage. `-t8` and `-t 8` are the same thing.

Three things it does not do:

- It does not parallelise loading the index, which is single-threaded I/O and
  can dominate wall time on a large reference.
- It does not apply to `samtools sort`, which is a separate process on the
  other side of the pipe and takes its own `-@`.
- It does not reduce memory. The index is held once and shared across threads,
  but each thread carries its own alignment buffers, so peak memory rises
  slowly with `-t`.

On a 16-core machine, `-t8` leaves headroom for the sort running concurrently.
Setting `-t16` here would have both programs contending for every core and
typically runs *slower* end to end.

### `Mus_musculus_GRCm39_2026_04.mmi` — first positional argument

The **target**: what you align to. minimap2 always takes the target first and
the query second; swapping them is silently accepted and produces meaningless
output, so this order is worth checking.

The `.mmi` extension means this is a **prebuilt index**, not a FASTA. It was
produced by an earlier `minimap2 -x splice -d ... genome.fa` run. Passing an
index rather than a FASTA skips re-indexing on every invocation, which for a
mammalian genome costs more than the alignment itself.

One trap: **`-k` and `-w` are baked into the `.mmi` at build time.** If this
index was built with a different preset, minimap2 uses the index's values and
ignores the `-k`/`-w` implied by `-x splice` here — with only a warning that
is easy to miss in a long log. An index built with `-x map-ont` (`w=10`) used
for splice alignment gives quietly worse sensitivity on short exons. If you
did not build the index yourself, confirm which preset was used.

`GRCm39` is the mouse assembly, so the annotation you hand nexons afterwards
must also be mouse GRCm39, with matching chromosome names — Ensembl-style
(`1`, `19`, `X`) rather than UCSC-style (`chr1`, `chr19`, `chrX`). A naming
mismatch is not an error; it produces zero counts everywhere.

### `PBM75668_pass_L001_barcode01.fastq.gz` — second positional argument

The **query**: your basecalled reads. `pass` indicates the basecaller's
quality filter has already been applied. minimap2 reads gzip transparently, so
compressed input costs nothing but decompression time, and multiple FASTQs can
be given at once.

---

## The join

### `|`

The pipe is a performance decision, not a formality.

SAM is uncompressed text. One ONT barcode of this size produces tens of
gigabytes of it. Writing that to disk and reading it back would add both the
space and roughly twice the I/O of the entire rest of the job. The pipe hands
minimap2's stdout straight to samtools' stdin, so the SAM exists only in a
64 KiB kernel buffer and the only file written is the final BAM.

The cost is that a failure anywhere in the pipeline loses the whole run — and
that by default **only the last command's exit status is reported**. A
minimap2 that dies half way through leaves you with a valid, silently
truncated BAM and a zero exit code. The script guards this with `set -o
pipefail`; a bare shell invocation of this command does not.

---

## The sorter

### `samtools sort`

minimap2 emits alignments in the order the reads appear in the FASTQ, which is
effectively random with respect to the genome. `sort` reorders them by
coordinate. This is required, not cosmetic: coordinate order is what makes an
index possible, and therefore what makes "fetch the reads over this exon" a
fast operation instead of a full-file scan.

It is also the memory-hungry half of the pipeline. Sorting works in batches
that are held in RAM, written out as temporary files when full, then merged.

### `-o out.bam`

The output path. samtools infers the format from the extension, so `.bam`
gives BAM; `.sam` here would write uncompressed text and `.cram` would need
the reference again at read time.

### `-` (trailing dash)

Read the input from **stdin**.

samtools expects an input filename as its last positional argument. Without
the dash it looks for a file, finds none, and errors — the pipe alone is not
enough. This is the parameter most often left off by accident.

---

## What the command omits

Three gaps, all of which `fastq_to_bam.sh` fills.

**1. No BAM index.** nexons queries the BAM by region, which requires a
companion `.bai`. The command above never creates one:

```bash
samtools index -@ 4 out.bam
```

**2. Sort defaults are conservative and the temp location is unspecified.**
Bare `samtools sort` uses one thread and 768 MiB, and spills to `$TMPDIR` —
which on this machine is a RAM-backed `tmpfs`, so spilled data competes for
the memory minimap2 is using for the index. Being explicit:

```bash
| samtools sort -@ 4 -m 512M -T /path/on/disk/sort.sample -o out.bam -
```

`-m` is **per thread**, so `-@ 4 -m 512M` reserves about 2 GiB total, not
512 MiB. This is a common way to accidentally request far more memory than
intended.

**3. No memory check on the index.** minimap2 holds the whole index resident
for the duration of the run. If it does not fit, the kernel kills the process
with no message beyond `Killed` — indistinguishable at a glance from a crash.

This machine has ~8 GiB total and ~5 GiB available. A splice index for the
whole ~2.7 Gbase mouse genome is several GiB, because the denser `w=5` window
stores roughly twice the minimizers of a genomic preset. Check before running:

```bash
ls -lh Mus_musculus_GRCm39_2026_04.mmi
```

Treat the file size as the resident requirement, and add ~1 GiB of overhead
plus whatever the sort is set to use. If that exceeds available memory, the
options are a larger machine, or a reference subset to the chromosomes your
GTF actually covers. `README_fastqexamples.md` covers both.

---

## The full command as the script runs it

```bash
minimap2 -ax splice --secondary=no -t 8 \
         Mus_musculus_GRCm39_2026_04.mmi \
         PBM75668_pass_L001_barcode01.fastq.gz \
  | samtools sort -@ 4 -m 512M \
                  -T big_data/tmp/sort.barcode01 \
                  -o big_data/bam/barcode01.bam -

samtools index -@ 4 big_data/bam/barcode01.bam
samtools flagstat big_data/bam/barcode01.bam
```

The `flagstat` line is the cheapest evidence the BAM is neither empty nor
truncated, and its mapped-read count is the denominator nexons' totals should
be read against.
