# nexons_viz — read-level visualisation of nexons classifications

`nexons` decides, for every read, which transcript or gene it belongs to, then
aggregates those decisions into a counts table and discards the per-read
detail. `nexons_viz` puts the detail back on screen, so that a surprising count
can be traced to the reads that produced it.

The distinctions the counts table cannot make, and this can:

| A count is lower than expected because… | visible as |
|---|---|
| the transcript is lowly expressed | few reads over the locus |
| reads are truncated | a wall of `partial` calls |
| a paralog overlaps the gene | a wall of `multi_gene` calls |
| the library direction is mis-set | reads on the opposing strand |

---

## The two-step workflow

Visualisation needs per-read decisions, and a plain nexons run does not emit
them. `nexons.py` therefore grew a `--annotated-bam` option (added on this
branch) that writes a tagged copy of each input BAM alongside the counts.

### 1. Quantitate, keeping the per-read calls

```bash
cd sandbox/data/input/example1
python ../../../../nexons.py \
    --maxtsl 0 \
    --annotated-bam ../../output/annotated/example1 \
    --outbase results/nexons_output \
    unc13a_toy.gtf control.bam tdp43_depleted.bam mixed.bam
```

Every read in the output BAMs carries five tags:

| tag | meaning |
|---|---|
| `XC` | classification: `unique`, `partial`, `gene`, `multi_gene`, `no_hit`, `no_gene`, `secondary`, `unmapped` |
| `XT` | assigned transcript (set only for `unique` and `partial`) |
| `XG` | gene(s) hit, comma-separated when ambiguous |
| `XD` | strand call — `S` same strand as the gene, `O` opposing |
| `XF` | the flex actually required for the match |

The counts nexons reports are derived from these same tags, so the figures
cannot disagree with the table.

### 2. Visualise

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

Writes, for each of the three figures, a `.html` (interactive) and a `.png`
(static), plus `_classification_counts.csv` and `_region_support.csv`.

---

## Three ways in

| | best for |
|---|---|
| `nexons_viz.py` as a **library** | scripted or batch analysis; full control over the figures |
| `nexons_viz.py` as a **CLI** | one gene, one command, files on disk |
| `app.py` (**Streamlit**) | exploring interactively — change gene or region without editing code |
| `../notebooks/nexons_viz_demo.ipynb` | learning what the figures mean |

Start with the notebook. It walks the whole workflow on the bundled dummy data
and explains how to read each figure and what each failure mode looks like.

```bash
cd sandbox
pixi run jupyter lab                                  # notebook
pixi run streamlit run viz_scripts/app.py             # browser app
```

---

## The three figures

### 1. Read pileup — `plot_read_pileup`

One horizontal line per read, coloured by classification, grouped by class and
sorted by position within group so a wrongly-coloured block stands out. The
gene model sits on top on the same axis.

Two things worth knowing:

- **Introns are compressed** (`gap=300`). A 150 bp cryptic exon inside a 20 kb
  gene is under 1% of the span and would otherwise be about one pixel wide.
  Pass `gap=None` for a true linear axis when real distances matter.
- **Reads are downsampled** to 150 per sample by default. Beyond a few hundred
  rows per panel each read is under a pixel tall and the pileup becomes a solid
  block. When reads are dropped the panel label says so, with the seed.

### 2. Composition — `plot_classification_composition`

Classification breakdown next to strand call, per sample. This is the QC gate:
read it before interpreting any count. Bars that differ in shape between
samples are a batch effect, and any biological comparison across them is
already confounded.

### 3. Region support — `plot_region_support`

Of the reads that could report on a region, how many include it?

**Only reads spanning the whole region are counted.** A read that stops inside
the region is evidence of neither inclusion nor skipping; scoring it as a skip
is the easiest way to turn a coverage difference into a splicing result.
Non-spanning reads are reported in their own column so the discard is visible.

Evidence tier is stacked inside each bar rather than summed away: the same
inclusion percentage carried by `partial` reads is a weaker claim than one
carried by `unique` reads.

---

## Library use

```python
import nexons_viz as nv

genes = nv.load_gene_models("unc13a_toy.gtf", max_tsl=0)
gene  = genes["UNC13A_TOY"]

reads = nv.load_reads({"control":  "control.annotated.bam",
                       "depleted": "tdp43_depleted.annotated.bam"}, gene)

region = ("UNC13A_TOY_chr", 29500, 29649)
nv.plot_read_pileup(reads, gene, region=region).show()
nv.plot_classification_composition(reads).show()
nv.plot_region_support(reads, region=region).show()

nv.classification_table(reads)
nv.region_support_table(reads, region=region)
```

`load_reads` returns one tidy row per alignment — `sample`, `read`, `chrom`,
`start`, `end`, `strand`, `blocks`, `xc`, `xt`, `xg`, `xd`, `xf` — so anything
the built-in figures don't cover can be done with pandas directly.

Note `max_tsl=0`. Nexons treats `0` as *no transcript-support-level filtering*
and translates it internally to `None`; the viewer copies that convention. A
GTF that yields no genes is nearly always this filter, and the loader says so
in the error.

---

## Dependencies

Declared in `sandbox/pixi.toml`. `plotly` is pinned to 6.x and `kaleido` to
0.2.1 on purpose: plotly 7 exports static images through kaleido ≥1, which
downloads and drives a headless Chrome on first use, and PNG export fails
without it. Plotly 6 with kaleido 0.2.1 renders offline. The HTML output never
depends on this — only the PNGs do.

**The lockfile has not been updated.** The pins above were added to
`pixi.toml` but `pixi` could not be run from the environment these files were
written in, so the packages were installed straight into
`.pixi/envs/default/`. Before committing, run:

```bash
cd sandbox
pixi install
```

to regenerate `pixi.lock` from the manifest.

---

## Files

| file | what it is |
|---|---|
| `nexons_viz.py` | the viewer: data layer, three figures, CLI |
| `app.py` | Streamlit front end |
| `build_viz_notebook.py` | regenerates `../notebooks/nexons_viz_demo.ipynb` |
| `README_goal.md` | the original hackathon brief, with answers |
| `README_plan.md` | the build plan these files were written against |

The notebook is a build product — edit `build_viz_notebook.py` and rebuild
rather than hand-editing the `.ipynb`:

```bash
cd sandbox
./.pixi/envs/default/bin/python viz_scripts/build_viz_notebook.py
```

---

## Troubleshooting

| symptom | cause |
|---|---|
| `MissingTagsError` | BAM came from a plain nexons run; re-run with `--annotated-bam` |
| `no genes found in GTF` | TSL filtering removed everything; use `max_tsl=0` |
| gene not found | pass the gene ID, not the display name; the error lists what's available |
| pileup is a solid block | lower `max_reads_per_sample` |
| HTML written, no PNG | kaleido missing or mismatched; the HTML is complete |
