#!/usr/bin/env python
"""
Regenerate sandbox/notebooks/nexons_viz_demo.ipynb.

The notebook is a build product, not something to hand-edit: keeping it in a
script means the prose and the code stay in one reviewable diff, and a broken
notebook can always be rebuilt with

    cd sandbox
    ./.pixi/envs/default/bin/python viz_scripts/build_viz_notebook.py
"""
from pathlib import Path

import nbformat as nbf

HERE = Path(__file__).resolve().parent
OUT = HERE.parent / "notebooks" / "nexons_viz_demo.ipynb"

nb = nbf.v4.new_notebook()
C = nb.cells
md = lambda s: C.append(nbf.v4.new_markdown_cell(s.strip("\n")))
code = lambda s: C.append(nbf.v4.new_code_cell(s.strip("\n")))


# ==========================================================================
md(r"""
# Seeing what nexons decided about every read

`nexons` quantitates long-read RNA-seq against a GTF. For each read it makes a
decision — *this read matches transcript X over its whole length*, *this one is
compatible with the gene but not with any single transcript*, *this one is
ambiguous between two genes* — and then it aggregates those decisions into a
counts table and **throws the per-read detail away**.

That is fine when the counts look sensible. It is unhelpful when they don't,
because a surprising count has at least four different explanations that the
table cannot distinguish:

| The count is low because… | What you'd see per read |
|---|---|
| the transcript really is lowly expressed | few reads over the locus at all |
| reads are truncated and score as `partial` | plenty of reads, most of them short |
| the locus overlaps a paralog and reads score as `multi_gene` | plenty of reads, all ambiguous |
| the library is reverse-stranded and reads hit the opposite strand | plenty of reads, all antisense |

This notebook walks through `nexons_viz`, which puts the per-read decision back
on screen. The workflow is two steps:

1. run `nexons.py` with `--annotated-bam`, which writes a copy of each input
   BAM with nexons' decision stamped onto every read as SAM tags;
2. point `nexons_viz` at those tagged BAMs.

**What you will get:** three figures per gene — a read pileup coloured by
classification, a whole-file composition summary, and a splice-inclusion plot
for a region you nominate — plus the tables behind them.

**How to run this notebook.** Launch JupyterLab from `sandbox/` so the pixi
environment is the kernel:

```bash
cd sandbox
pixi run jupyter lab
```

then open `notebooks/nexons_viz_demo.ipynb` and pick the pixi kernel.
""")

# ==========================================================================
md(r"""
## 0. Paths

Everything is resolved relative to `sandbox/`, found by walking up from
wherever the notebook was opened. Nothing is written outside `sandbox/data/`.

- inputs (GTF + BAMs made by `make_dummy_data.py`) → `sandbox/data/input/`
- tagged BAMs → `sandbox/data/output/annotated/`
- figures and tables → `sandbox/data/output/figures/`
""")

code(r"""
import sys, subprocess
from pathlib import Path

import pandas as pd

here = Path.cwd()
SANDBOX = next((p for p in (here, *here.parents)
                if (p / "data" / "make_dummy_data.py").exists()), None)
assert SANDBOX is not None, (
    f"Could not locate sandbox/ from {here}. Open this notebook from "
    "sandbox/notebooks/ after launching `pixi run jupyter lab` in sandbox/."
)

DATA        = SANDBOX / "data"
DATA_INPUT  = DATA / "input"
DATA_OUTPUT = DATA / "output"
ANNOTATED   = DATA_OUTPUT / "annotated"
FIGURES     = DATA_OUTPUT / "figures"
for d in (ANNOTATED, FIGURES):
    d.mkdir(parents=True, exist_ok=True)

NEXONS = SANDBOX.parent / "nexons.py"
sys.path.insert(0, str(SANDBOX / "viz_scripts"))
sys.path.insert(0, str(DATA))

import nexons_viz as nv

print("sandbox   :", SANDBOX)
print("nexons.py :", NEXONS, "(found)" if NEXONS.exists() else "(MISSING)")
print("viewer    : nexons_viz from", Path(nv.__file__).relative_to(SANDBOX))
""")

# ==========================================================================
md(r"""
## 1. Build the dummy dataset

`make_dummy_data.py` writes three toy datasets under `data/input/`. All of them
use a synthetic gene `UNC13A_TOY` carrying a **cryptic exon** — modelled on the
real UNC13A cryptic exon that appears when TDP-43 is lost in ALS and FTD. Each
dataset has three samples (`control`, `tdp43_depleted`, `mixed`) that differ in
how often reads splice into that cryptic exon.

| dataset | what it is built to show |
|---|---|
| `example1` | the clean case: cryptic inclusion rises from control to depleted |
| `example2` | an overlapping paralog, so many reads become `multi_gene` |
| `example3` | a third of reads on the opposite strand, as with a mis-set library direction |

Skip this cell if `data/input/` is already populated.
""")

code(r"""
import make_dummy_data

if not (DATA_INPUT / "example1" / "control.bam").exists():
    make_dummy_data.main()
else:
    print("data/input already populated - skipping generation")

for ex in sorted(DATA_INPUT.glob("example*")):
    print(ex.name, "->", ", ".join(sorted(p.name for p in ex.iterdir())))
""")

# ==========================================================================
md(r"""
## 2. Run nexons with `--annotated-bam`

This is the step that makes the visualisation possible. `--annotated-bam DIR`
tells nexons to write, alongside its usual counts table, a copy of each input
BAM into `DIR` with five tags added to every read:

| tag | meaning |
|---|---|
| `XC` | the classification: `unique`, `partial`, `gene`, `multi_gene`, `no_hit`, `no_gene`, `secondary`, `unmapped` |
| `XT` | transcript the read was assigned to (empty unless `XC` is `unique` or `partial`) |
| `XG` | gene(s) the read hit, comma-separated when ambiguous |
| `XD` | strand call: `S` same strand as the gene, `O` opposing |
| `XF` | the flex/tolerance actually needed to make the match |

The counts nexons reports are derived from exactly these tags, so a figure
built from the tags cannot disagree with the table — which is the point.
""")

code(r"""
SAMPLES = ["control.bam", "tdp43_depleted.bam", "mixed.bam"]

# Quantitate one example and write tagged BAMs; returns the tagged-BAM dir.
def run_nexons(example: int) -> Path:
    src = DATA_INPUT / f"example{example}"
    out = ANNOTATED / f"example{example}"
    out.mkdir(parents=True, exist_ok=True)
    (src / "results").mkdir(exist_ok=True)
    cmd = [sys.executable, str(NEXONS),
           "--maxtsl", "0",
           "--annotated-bam", str(out),
           "--outbase", "results/nexons_output",
           "unc13a_toy.gtf", *SAMPLES]
    r = subprocess.run(cmd, cwd=src, capture_output=True, text=True)
    if r.returncode:
        raise RuntimeError(f"nexons failed on example{example}:\n{r.stderr[-3000:]}")
    return out

TAGGED = {ex: run_nexons(ex) for ex in (1, 2, 3)}
for ex, d in TAGGED.items():
    print(f"example{ex} ->", ", ".join(sorted(p.name for p in d.glob("*.bam"))))
""")

md(r"""
A quick look at what a tagged read actually carries. If this cell shows tags,
everything downstream will work; if it doesn't, step 2 didn't run.
""")

code(r"""
import pysam

with pysam.AlignmentFile(TAGGED[1] / "control.annotated.bam") as bam:
    for read in bam.head(3):
        print(read.query_name,
              {t: read.get_tag(t) for t in ("XC", "XT", "XG", "XD", "XF")
               if read.has_tag(t)})
""")

# ==========================================================================
md(r"""
## 3. Load the gene model and the reads

Two loaders, both thin:

- `load_gene_models(gtf)` reads the annotation **using nexons' own GTF parser**,
  so the viewer sees exactly the transcripts nexons quantitated against. Note
  `max_tsl=0`: nexons treats `0` as *no transcript-support-level filtering*,
  and the viewer copies that convention so the two never disagree.
- `load_reads(bams, gene, names=...)` returns one tidy row per alignment, with
  its classification, strand call and exon blocks.
""")

code(r"""
GTF = DATA_INPUT / "example1" / "unc13a_toy.gtf"
NAMES = [s.replace(".bam", "") for s in SAMPLES]
# load_reads takes {label: path}, so the panel labels are ours to choose
# rather than being whatever the file happens to be called.
BAMS = {n: TAGGED[1] / f"{n}.annotated.bam" for n in NAMES}

genes = nv.load_gene_models(GTF, max_tsl=0)
gene = genes["UNC13A_TOY"]
print(f"{gene['name']}  {gene['chrom']}:{gene['start']:,}-{gene['end']:,} "
      f"({gene['strand']})  {len(gene['transcripts'])} transcripts")
for tid, tx in gene["transcripts"].items():
    print(f"   {tid:26s} {len(tx['exons']):2d} exons")

reads = nv.load_reads(BAMS, gene)
print(f"\n{len(reads):,} alignments loaded")
reads.head(3)
""")

md(r"""
The cryptic exon is the feature these datasets are built around. Its
coordinates come from the generator rather than being typed in by hand, so the
region highlighted in the figures is guaranteed to be the one the data was
simulated with.
""")

code(r"""
REGION = (make_dummy_data.CHROM, *make_dummy_data.CRYPTIC_EXON)
print("region of interest:  {}:{:,}-{:,}  ({} bp)".format(
    *REGION, REGION[2] - REGION[1]))
""")

# ==========================================================================
md(r"""
## 4. Figure 1 — the read pileup

Every read is one horizontal line: thick blocks are aligned exon segments, thin
connectors are the gaps between them (introns, or deletions). Reads are grouped
by classification and sorted by start position within each group, so a block of
the wrong colour is easy to spot. The gene model sits on top on the same axis.

**Introns are compressed.** The cryptic exon here is ~150 bp inside a ~20 kb
gene — about 0.7% of the span. On a true linear axis it would be roughly one
pixel wide. `gap=300` squashes every intron to 300 bp of screen space while
leaving exons at true scale. Pass `gap=None` for a genuine linear axis when you
need to judge real distances.

Only 150 reads per sample are drawn by default; when reads are dropped, the
panel label says so. Raise `max_reads_per_sample` if you want more, but past a
few hundred reads per panel each read is under a pixel tall and the pileup
becomes a solid block.
""")

code(r"""
fig = nv.plot_read_pileup(reads, gene, region=REGION, max_reads_per_sample=120)
fig.show()
""")

md(r"""
### Reading it

The shaded band is the cryptic exon.

- **Dark blue (`unique`)** — the read matched one transcript along its whole
  length. This is the only class that contributes to a confident
  transcript-level count.
- **Light blue (`partial`)** — the read was compatible with one transcript but
  didn't cover all of it. Usually 5′ truncation: a real feature of long-read
  data, not an error, but it means the read can't discriminate between
  transcripts that differ only in the region it doesn't cover.
- **Green (`gene`)** — the read is compatible with the gene but with no single
  annotated transcript. **This is the interesting class.** A read that splices
  into an unannotated cryptic exon lands here, because the annotation has no
  transcript containing that exon. A green block that consistently covers the
  same unannotated region is the signature of a novel splice event; green reads
  scattered with no common structure are more likely alignment noise.
- **Grey (`secondary`)** — not the primary alignment; nexons does not count
  these. They are drawn so you can see whether a locus is attracting large
  numbers of secondary alignments, which usually means repeat content.

In this dataset the useful comparison is *within* the green and blue groups
across panels: reads that traverse the shaded band versus reads that jump over
it. Zoom into the band in the interactive figure to see it directly — that
comparison is what Figure 3 quantifies.
""")

# ==========================================================================
md(r"""
## 5. Figure 2 — where the reads went

The pileup shows individual reads; this is the whole-file view. Left panel: the
classification breakdown per sample. Right panel: the strand call, computed
over reads nexons actually assigned to the gene (`no call` covers reads with no
strand decision, including secondary alignments and ambiguous assignments).

Read this **before** interpreting any count, as a QC gate.
""")

code(r"""
fig = nv.plot_classification_composition(reads)
fig.show()

counts = nv.classification_table(reads)
counts
""")

md(r"""
### Reading it

What you want to see: bars of similar shape across samples, dominated by
`unique` and `partial`, with the strand panel almost entirely one colour.

What each departure means:

- **A large `multi_gene` (purple) fraction** — reads are ambiguous between
  overlapping genes. Per-gene counts at this locus are not trustworthy; nexons
  cannot attribute those reads and neither can you. Check whether a paralog or
  a readthrough transcript overlaps the gene.
- **A large `no_hit` fraction** — reads are aligning to the locus but their
  splice structure doesn't match the annotation at all. Either the annotation
  is wrong for this sample, or the alignment is.
- **A substantial `opposite strand` fraction** — the library direction is
  probably mis-set. Re-run nexons with the correct `--direction`; the counts
  from this run are wrong, not merely noisy.
- **Bars that differ in shape between samples** — a batch effect in library
  prep or alignment. Any biological comparison across those samples is
  confounded before it starts.

Example 1 is the clean case: identical composition across the three samples,
so the differences Figure 3 finds are differences in splicing rather than in
data quality. Sections 7 and 8 show what the two pathological cases look like.
""")

# ==========================================================================
md(r"""
## 6. Figure 3 — inclusion of the region of interest

This is the quantitative claim: of the reads that could report on the cryptic
exon, how many include it?

**Only reads that span the entire region are counted.** A read that stops
inside the region is evidence of neither inclusion nor skipping, and counting
it as a skip is the easiest way to turn a coverage difference into a splicing
result. Non-spanning reads are reported as their own column so you can see how
much of the data was set aside.

Evidence tier is stacked inside each bar rather than summed away. An inclusion
percentage carried by `unique` reads is a stronger claim than the same
percentage carried by `partial` reads, and the figure keeps that visible.
""")

code(r"""
fig = nv.plot_region_support(reads, region=REGION)
fig.show()

support = nv.region_support_table(reads, region=REGION)
support[["not_spanning", "spanning", "include", "skip", "percent_included"]]
""")

md(r"""
### Reading it

Left panel — absolute spanning-read counts, split into *include* (blue) and
*skip* (orange), each stacked by evidence tier. Right panel — the same thing as
a percentage.

Three checks before believing the percentage:

1. **Is the spanning count large enough?** A 60% inclusion rate from 10
   spanning reads is a rate estimated from 10 observations. The left panel is
   there so the denominator is never hidden.
2. **How much was set aside?** If `not_spanning` is most of the reads, the
   region is at the edge of your coverage and the estimate rests on a small,
   possibly non-random subset.
3. **Which tier carries the signal?** If the difference between samples lives
   entirely in the pale (`other`) blocks, it is being carried by reads nexons
   could not assign confidently.

In this dataset inclusion rises from a few percent in `control` to a majority
in `tdp43_depleted`, with `mixed` in between — the expected direction for a
cryptic exon de-repressed by loss of TDP-43. The signal is carried mostly by
`unique` reads and the spanning denominators are in the hundreds, so all three
checks pass.

Note what the counts table in Section 5 does *not* show: it is identical across
the three samples. The splicing change is invisible at the level of
classification composition and only appears when you ask about the region.
""")

# ==========================================================================
md(r"""
## 7. Example 2 — an overlapping paralog

Same gene, same cryptic exon, but a paralog now overlaps the locus. The
biology hasn't changed; the interpretability has.
""")

code(r"""
genes2 = nv.load_gene_models(DATA_INPUT / "example2" / "unc13a_toy.gtf", max_tsl=0)
reads2 = nv.load_reads({n: TAGGED[2] / f"{n}.annotated.bam" for n in NAMES},
                       genes2["UNC13A_TOY"])
nv.plot_classification_composition(
    reads2, title="Example 2: overlapping paralog").show()
nv.classification_table(reads2)
""")

md(r"""
Roughly half of all reads are now `multi_gene`. Nexons cannot decide whether
they belong to `UNC13A_TOY` or to the paralog, so it declines to attribute
them — which is the correct behaviour, and exactly why the per-gene count
drops even though the same number of reads is present.

The lesson: a count that halves between two annotations is not necessarily a
change in expression. Look at this figure before concluding anything from the
counts table.
""")

# ==========================================================================
md(r"""
## 8. Example 3 — a mis-set library direction
""")

code(r"""
genes3 = nv.load_gene_models(DATA_INPUT / "example3" / "unc13a_toy.gtf", max_tsl=0)
reads3 = nv.load_reads({n: TAGGED[3] / f"{n}.annotated.bam" for n in NAMES},
                       genes3["UNC13A_TOY"])
nv.plot_classification_composition(
    reads3, title="Example 3: antisense contamination").show()
reads3.groupby(["sample", "xd"], observed=True).size().unstack(fill_value=0)
""")

md(r"""
About a third of assigned reads carry `XD=O` — they align to the opposite
strand from the gene. In real data this is almost always a library-direction
setting that doesn't match the protocol, not antisense transcription. The fix
is to re-run nexons with the correct `--direction` rather than to reinterpret
the numbers.

Note that the classification panel looks unremarkable here. Without the strand
panel this failure mode is silent.
""")

# ==========================================================================
md(r"""
## 9. Saving figures

`save_figure` always writes interactive HTML and additionally writes a static
PNG when an image backend is installed. HTML is the more useful of the two: the
pileup is worth zooming into, and the static version can't be.
""")

code(r"""
written = []
written += nv.save_figure(nv.plot_read_pileup(reads, gene, region=REGION,
                                              max_reads_per_sample=120),
                          FIGURES / "notebook_pileup")
written += nv.save_figure(nv.plot_classification_composition(reads),
                          FIGURES / "notebook_composition")
written += nv.save_figure(nv.plot_region_support(reads, region=REGION),
                          FIGURES / "notebook_region_support")
support.to_csv(FIGURES / "notebook_region_support.csv")
counts.to_csv(FIGURES / "notebook_classification_counts.csv")

for p in written:
    print("wrote", Path(p).relative_to(SANDBOX))
""")

# ==========================================================================
md(r"""
## 10. Doing this from the shell, or from a browser

Everything above is also a command-line tool:

```bash
cd sandbox/viz_scripts
python nexons_viz.py \
    --gtf  ../data/input/example1/unc13a_toy.gtf \
    --gene UNC13A_TOY \
    --bam  ../data/output/annotated/example1/*.annotated.bam \
    --name control tdp43_depleted mixed \
    --region UNC13A_TOY_chr:29500-29649 \
    --out  ../data/output/figures/example1
```

and a Streamlit app for exploring without writing code — pick a bundled
example or upload your own tagged BAMs, then change gene and region
interactively:

```bash
cd sandbox
./.pixi/envs/default/bin/streamlit run viz_scripts/app.py
```

### If something goes wrong

| symptom | cause |
|---|---|
| `MissingTagsError` | the BAM came from a plain nexons run; re-run with `--annotated-bam` |
| `no genes in GTF` | transcript-support-level filtering removed everything; use `max_tsl=0` |
| gene not found | pass the gene ID rather than the display name; the error lists what is available |
| pileup is a solid block | too many reads per panel; lower `max_reads_per_sample` |
| no PNG written, HTML fine | the static-image backend isn't installed; the HTML is complete |
""")

nb.metadata["kernelspec"] = {"display_name": "Python 3 (pixi)",
                             "language": "python", "name": "python3"}
OUT.parent.mkdir(parents=True, exist_ok=True)
nbf.write(nb, OUT)
print(f"wrote {OUT} ({len(nb.cells)} cells)")
