# Implementation plan — nexons read-classification viewer

> Status: **awaiting approval**. Nothing in this document has been built yet.
> It is the implementation plan for the viewer described in
> [`README_goal.md`](README_goal.md), and it answers that document's §0 finding
> (nexons does not currently emit annotated BAMs) before building anything on
> top of it.

---

## What is being built

Three things, sharing one library:

| Deliverable | File | Role |
|---|---|---|
| Library + CLI | `sandbox/viz_scripts/nexons_viz.py` | Reads annotated BAMs + GTF, builds the figures. Importable, and runnable from a shell. |
| Notebook | `sandbox/viz_scripts/nexons_viz_demo.ipynb` | Worked example on the sandbox data, plus the interpretation guide for every chart. |
| Streamlit app | `sandbox/viz_scripts/app.py` | Upload-and-explore interface over the same library. |

Plus the prerequisite from `README_goal.md` §0: a `--annotated-bam` option in
`nexons.py` itself, without which none of the above has an input.

## Feasibility and honest caveats

Confidence: **high**. The classification logic already exists in a single
self-contained function (`process_bam_file()` in `nexons.py`), so emitting
per-read tags is a contained addition rather than a reimplementation of the
decision tree. The viewer itself is standard genomic-track rendering, and the
three synthetic datasets in `sandbox/data/input/` were built to exercise
exactly the display paths that matter.

Two limitations to state up front:

- **The Streamlit app will not be click-tested during development.** The build
  environment has no browser. It will be validated headlessly — module imports,
  the render path called directly against the demo fixtures, and
  `streamlit run --headless` confirmed to start and serve — but the first
  person to actually use it in a browser will be you.
- **Read-level rendering has a legibility ceiling** of a few thousand reads per
  panel. Above a configurable cap the library will downsample and state the
  sampling on the figure rather than silently drawing an unreadable panel.

---

## Step 1 — Add `--annotated-bam` tag writer to `nexons.py`

Extend `process_bam_file()` to optionally stream the input BAM back out with
per-read classification tags, and add the `--annotated-bam` / `-a` CLI option.

Tags use the `X?` space the SAM specification reserves for local use:

| Tag | Type | Content |
|---|---|---|
| `XC` | `Z` | `unique` \| `partial` \| `gene` \| `multi_gene` \| `no_hit` \| `no_gene` \| `secondary` \| `unmapped` |
| `XT` | `Z` | assigned transcript ID, or `.` when not transcript-resolved |
| `XG` | `Z` | assigned gene ID (comma-joined when multi-gene) |
| `XD` | `A` | `S` same strand / `O` opposing strand |
| `XF` | `i` | number of boundaries that matched only within `--flex`/`--endflex` tolerance |

Every read in the input is written out exactly once — including unmapped and
secondary reads — so the annotated BAM is a faithful superset of the input, not
a filtered subset.

**Acceptance test:** re-run examples 1–3 and assert that the tag distribution
reproduces the counts in the existing `nexons_output_*.txt` tables *exactly*.
If the tags disagree with the count table, the viewer is lying about the data,
which defeats the entire purpose of the tool.

**Deliverable:** patched `nexons.py`.

## Step 2 — Generate annotated BAMs for the three examples

Run the patched nexons over `sandbox/data/input/example{1,2,3}/`, producing
sorted and indexed tagged BAMs at:

```
sandbox/data/output/annotated/example{1,2,3}/{control,tdp43_depleted,mixed}.annotated.bam
```

Print a per-file tag census (`XC` value counts, `XD` split) and cross-check it
against `nexons_output_gene.txt` / `_unique.txt`.

**Deliverable:** nine annotated BAMs, which become the fixtures for every later
step and for the Streamlit app's demo mode.

## Step 3 — Build the `nexons_viz.py` data layer

The library's input half:

- **GTF reader** — reuses `nexons.read_gtf` so the viewer and the quantifier
  can never disagree about what the annotation says.
- **Annotated-BAM reader** — turns reads into a tidy pandas DataFrame
  (`read_id`, `sample`, exon blocks from the CIGAR via nexons' own
  `get_exons`, `XC`, `XT`, `XG`, `XD`, `XF`, strand).
- **`ShrunkAxis` coordinate transform** — compresses introns to a fixed budget
  so a 128 bp cryptic exon stays visible in a ~90 kb gene. Built from the union
  of annotated exons and observed read blocks; exposes forward and inverse
  mapping plus tick positions in real genomic coordinates.

**Deliverable:** module with inline assertions on the transform (monotonic,
exon widths preserved, round-trip inverse correct).

## Step 4 — Build the read-pileup figure

The headline view. One Plotly panel per sample, stacked on a shared shrunk
x-axis, with the gene model track (all annotated transcripts, exons as blocks)
pinned at the top and reads as rows below — thick segments for exon blocks,
thin connectors for introns, coloured by `XC` classification using a fixed,
colourblind-safe palette shared across every figure in the project.

- Reads grouped by classification, then sorted by start position.
- Hover shows read ID, `XC`, `XT`, `XG`, `XD`, `XF`.
- Optional highlight band for a region of interest (e.g. the cryptic exon).
- One trace per category using NaN-separated coordinate arrays, so a
  3000-read panel stays responsive instead of becoming 3000 shape objects.
- Downsamples above a configurable cap, with the sampling stated in the panel
  title.

**Deliverable:** `plot_read_pileup()` plus `read_pileup_example1.html` / `.png`
for the three-sample `UNC13A_TOY` view.

## Step 5 — Build the composition and region-support figures

Two supporting views that make the pileup quantitative.

- **`plot_classification_composition()`** — per-sample stacked bars of `XC`
  categories, with a second panel for the `XD` strand split. This is the figure
  that makes example2's 50% multi-gene and example3's 33% antisense read at a
  glance.
- **`plot_region_support()`** — for a user-specified region of interest,
  per-sample counts of reads that **include** vs. **skip** it, split by whether
  the supporting read was `unique` or `partial`. This is the
  NMD/truncation-confounder check from the ALS use case in `README_goal.md`
  §2.2: it distinguishes "low cryptic inclusion" from "the evidence was
  degraded before sequencing".

**Deliverable:** both functions plus `composition_all_examples.png` and
`region_support_example1.png`.

## Step 6 — Add a CLI entry point

Make the library runnable standalone:

```bash
python nexons_viz.py \
    --gtf unc13a_toy.gtf \
    --bam control.annotated.bam tdp43_depleted.annotated.bam mixed.annotated.bam \
    --gene UNC13A_TOY \
    --region UNC13A_TOY_chr:29400-29700 \
    --out prefix
```

Writes a self-contained interactive HTML plus static PNGs. argparse help text
names the tag requirement and points at `nexons --annotated-bam`. Handed a BAM
with no `XC` tags, it fails with a clear actionable message rather than a
traceback.

**Deliverable:** working CLI, verified on example1 and example2.

## Step 7 — Write the demo notebook with interpretation guide

`sandbox/viz_scripts/nexons_viz_demo.ipynb`, on the existing
**"Python (nexons pixi)"** kernel and following the sandbox path conventions
already established (`data/input/` for inputs, `data/output/` for products).

Sections:

1. What the viewer needs as input, and how to produce it.
2. Loading example1.
3. Each of the three figure types, rendered inline and interactive.
4. **The interpretation guide** — the substantive part. Per chart: what a
   healthy panel looks like, what each classification colour means
   biologically, and the specific artefacts to look for — partial-read-driven
   signal, multi-gene ambiguity, strand errors, annotation-bounded blind spots
   — each illustrated by switching to whichever of example1/2/3 exhibits it.
5. Reading the same views for the cryptic-exon contrast.

**Deliverable:** notebook, validated by executing every cell.

## Step 8 — Build and validate the Streamlit app

`sandbox/viz_scripts/app.py`:

- Sidebar: either upload a GTF plus one-or-more annotated BAM/BAI pairs, or
  pick one of the bundled sandbox examples.
- Gene selector populated from the GTF; region-of-interest input; sample
  multiselect; classification filter checkboxes.
- The three figures in tabs, plus a download button for the HTML export.
- `st.cache_data` on the BAM parse, so re-rendering on a filter change is
  instant rather than re-reading the BAM.

Adds `plotly`, `streamlit` and `kaleido` to `sandbox/pixi.toml`.

**Validation:** headless only, as noted in the caveats above — imports, direct
render-path call against the demo fixtures, and `streamlit run --headless`
confirmed to start and serve.

**Deliverable:** `app.py` and updated `pixi.toml` / `pixi.lock`.

## Step 9 — Write the viz_scripts README and save artifacts

`sandbox/viz_scripts/README.md` covering the tag contract (what each tag means
and which nexons version emits them), install, the three ways to run the viewer
(library, CLI, Streamlit), the downsampling cap and known limits, and an
explicit statement of what was and was not validated. Cross-referenced from
`README_goal.md` as the implementation of the answers written there.

---

## Dependencies added

`plotly`, `streamlit`, `kaleido` — into `sandbox/pixi.toml`, from the
`conda-forge` channel already configured there.

## Fallback position

If the Streamlit route proves awkward (BAM + BAI pair uploads are the likely
friction point, since Streamlit's uploader hands you file buffers and pysam
wants indexed paths on disk), the plan degrades gracefully: the notebook plus
the interactive HTML export from the CLI already deliver the full
visualization, and the app becomes an optional convenience layer to revisit
rather than a blocker.
