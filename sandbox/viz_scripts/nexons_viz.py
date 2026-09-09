#!/usr/bin/env python
"""
nexons_viz - read level visualisation of nexons classifications.

Nexons assigns every read in a BAM to a category (unique / partial / gene /
multi_gene / no_hit / no_gene) and then throws that decision away, keeping only
the totals.  Run with ``--annotated-bam`` it instead writes the reads back out
carrying their decision as tags, and this module turns those tagged BAMs into
pictures:

    plot_read_pileup()               every read, drawn, coloured by its call
    plot_classification_composition() what fraction of reads went where
    plot_region_support()            include vs skip over a region of interest

The point of all three is to let a human check whether a number in
``nexons_output_unique.txt`` is supported by reads that look the way they
ought to, before that number becomes a result.

Tags consumed (written by ``nexons.py --annotated-bam``):

    XC Z  classification
    XT Z  assigned transcript, '.' if none
    XG Z  assigned gene(s), comma joined for multi_gene, '.' if none
    XD A  S / O same or opposing strand, '.' where nexons made no call
    XF i  number of exon boundaries that only matched within the flex
          tolerance, -1 where no transcript match was attempted

Usage as a library::

    import nexons_viz as nv
    genes  = nv.load_gene_models("unc13a_toy.gtf")
    reads  = nv.load_reads({"control": "control.annotated.bam"}, gene=genes["UNC13A_TOY"])
    fig    = nv.plot_read_pileup(reads, genes["UNC13A_TOY"])
    fig.show()

Usage from a shell: see ``python nexons_viz.py --help``.
"""

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pysam
import plotly.graph_objects as go
from plotly.subplots import make_subplots

# nexons.py lives two directories up from sandbox/viz_scripts/.  We import it
# rather than reimplementing read_gtf/get_exons so that the viewer and the
# quantifier can never disagree about what the annotation says or where a
# read's exons are.
_REPO_ROOT = Path(__file__).resolve().parents[2]
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

import nexons as _nexons                                          # noqa: E402


class _QuietOptions:
    """nexons' log/warn/debug helpers read a module level `options` global."""
    verbose = False
    quiet = True


if not hasattr(_nexons, "options"):
    _nexons.options = _QuietOptions()


# --------------------------------------------------------------------------
# Presentation constants
# --------------------------------------------------------------------------

#: Classification order, most to least confident.  Used for stacking, legend
#: order and row grouping so every figure in the project reads the same way.
CLASS_ORDER = [
    "unique", "partial", "gene", "multi_gene",
    "no_hit", "no_gene", "secondary", "unmapped",
]

#: Okabe-Ito, which is colourblind safe.  multi_gene is purple to match the
#: existing nexons QC report, where the multi gene slice is already purple.
CLASS_COLOURS = {
    "unique":     "#0072B2",   # blue
    "partial":    "#56B4E9",   # sky blue
    "gene":       "#009E73",   # green
    "multi_gene": "#CC79A7",   # purple
    "no_hit":     "#D55E00",   # vermillion
    "no_gene":    "#E69F00",   # orange
    "secondary":  "#999999",
    "unmapped":   "#CCCCCC",
}

CLASS_MEANING = {
    "unique":     "matched one transcript over its whole length",
    "partial":    "matched one transcript but did not cover all of it",
    "gene":       "compatible with the gene but not with a single transcript",
    "multi_gene": "compatible with more than one gene, so not assigned",
    "no_hit":     "inside a gene but matching no transcript model",
    "no_gene":    "no gene surrounds this alignment",
    "secondary":  "not the primary alignment, not quantitated",
    "unmapped":   "no alignment reported",
}

REQUIRED_TAGS = ("XC", "XT", "XG", "XD", "XF")


class MissingTagsError(RuntimeError):
    """Raised when a BAM has not been through `nexons --annotated-bam`."""


# --------------------------------------------------------------------------
# Inputs
# --------------------------------------------------------------------------

def load_gene_models(gtf_file, max_tsl=None):
    """Read a GTF into nexons' own gene/transcript/exon structure.

    `max_tsl` follows the nexons CLI convention where 0 and None both mean
    "no transcript support level filtering at all".
    """
    if max_tsl == 0:
        max_tsl = None
    genes = _nexons.read_gtf(str(gtf_file), max_tsl)
    if not genes:
        raise ValueError(
            f"No genes read from {gtf_file}. If the GTF has no "
            f"transcript_support_level attributes, use max_tsl=None."
        )
    return genes


def parse_region(text):
    """'chr:1000-2000' -> ('chr', 1000, 2000). Commas and spaces tolerated."""
    if text is None:
        return None
    text = str(text).replace(",", "").replace(" ", "")
    try:
        chrom, span = text.rsplit(":", 1)
        start, end = span.split("-")
        return chrom, int(start), int(end)
    except ValueError:
        raise ValueError(f"Could not parse region '{text}', expected chr:start-end")


def load_reads(bams, gene=None, region=None, require_tags=True):
    """Load tagged reads into a tidy DataFrame, one row per alignment record.

    Parameters
    ----------
    bams : dict {sample_name: path} or list of paths
        Annotated BAMs. A list gets sample names from the file stems.
    gene : gene dict from load_gene_models, optional
        Restrict to reads overlapping this gene's span (uses the BAM index).
    region : (chrom, start, end), optional
        Restrict to this window instead. Ignored if `gene` is given.
    require_tags : bool
        Raise MissingTagsError if the first record carries no XC tag.

    Returns
    -------
    DataFrame with columns sample, read_id, chrom, start, end, reverse, blocks,
    n_blocks, length, xc, xt, xg, xd, xf.  `blocks` holds the read's exon
    blocks as a list of (start, end) pairs, taken from nexons' own get_exons
    so they are identical to the ones the classification was made on.
    """
    if not isinstance(bams, dict):
        bams = {Path(b).name.replace(".annotated.bam", "").replace(".bam", ""): b
                for b in bams}

    fetch_at = None
    if gene is not None:
        fetch_at = (gene["chrom"], gene["start"] - 1, gene["end"])
    elif region is not None:
        fetch_at = (region[0], region[1] - 1, region[2])

    rows = []
    for sample, path in bams.items():
        handle = pysam.AlignmentFile(str(path), "rb")

        if require_tags:
            _assert_tagged(handle, path)

        if fetch_at is not None and handle.has_index():
            try:
                iterator = handle.fetch(*fetch_at)
            except ValueError:
                # Contig not in this BAM's header
                iterator = iter(())
        else:
            iterator = handle.fetch(until_eof=True)

        for read in iterator:
            if read.is_unmapped:
                blocks = []
            else:
                blocks = [tuple(b) for b in _nexons.get_exons(read)]

            rows.append({
                "sample":   sample,
                "read_id":  read.query_name,
                "chrom":    read.reference_name,
                "start":    read.reference_start + 1 if not read.is_unmapped else None,
                "end":      read.reference_end if not read.is_unmapped else None,
                "reverse":  bool(read.is_reverse),
                "blocks":   blocks,
                "n_blocks": len(blocks),
                "length":   read.query_length,
                "xc":       read.get_tag("XC") if read.has_tag("XC") else "untagged",
                "xt":       read.get_tag("XT") if read.has_tag("XT") else ".",
                "xg":       read.get_tag("XG") if read.has_tag("XG") else ".",
                "xd":       read.get_tag("XD") if read.has_tag("XD") else ".",
                "xf":       read.get_tag("XF") if read.has_tag("XF") else -1,
            })
        handle.close()

    reads = pd.DataFrame(rows)
    if reads.empty:
        return reads

    reads["sample"] = pd.Categorical(reads["sample"], categories=list(bams), ordered=True)
    present = [c for c in CLASS_ORDER if c in set(reads["xc"])]
    present += sorted(set(reads["xc"]) - set(present))
    reads["xc"] = pd.Categorical(reads["xc"], categories=present, ordered=True)
    return reads


def _assert_tagged(handle, path):
    """Look at the first few records and complain helpfully if untagged."""
    peek = pysam.AlignmentFile(handle.filename, "rb")
    try:
        for i, read in enumerate(peek.fetch(until_eof=True)):
            if read.has_tag("XC"):
                return
            if i >= 100:
                break
    finally:
        peek.close()
    raise MissingTagsError(
        f"{path} carries no XC tags, so it has not been through nexons.\n"
        f"Produce an annotated BAM first, for example:\n"
        f"    mkdir -p annotated\n"
        f"    python nexons.py --maxtsl 0 --annotated-bam annotated \\\n"
        f"        --outbase annotated/nx your.gtf {Path(path).name}\n"
        f"then point this tool at annotated/{Path(path).stem}.annotated.bam"
    )


# --------------------------------------------------------------------------
# Squashed coordinates
# --------------------------------------------------------------------------

class ShrunkAxis:
    """Piecewise linear genome -> plot coordinate map that shrinks introns.

    A cryptic exon of 128 bp inside a 90 kb gene is roughly one thousandth of
    the axis, which is less than a pixel on a normal screen: on a true linear
    axis the interesting part of these figures is invisible.  This keeps exonic
    regions at their real scale and compresses everything between them to a
    fixed width, which is the standard trick used by isoform viewers.

    Set `gap=None` to disable shrinking and get a plain linear axis, which is
    what you want when the exact intron sizes matter.
    """

    def __init__(self, keep, lo, hi, gap=300, pad=200):
        self.lo = int(lo) - pad
        self.hi = int(hi) + pad
        self.gap = gap

        merged = self._merge([(max(int(s), self.lo), min(int(e), self.hi))
                              for s, e in keep if int(e) >= self.lo and int(s) <= self.hi])

        gen, plot = [self.lo], [0.0]
        cursor, offset = self.lo, 0.0
        for s, e in merged:
            if s > cursor:
                width = (s - cursor) if gap is None else min(gap, s - cursor)
                offset += width
                gen.append(s)
                plot.append(offset)
                cursor = s
            if e > cursor:
                offset += e - cursor
                gen.append(e)
                plot.append(offset)
                cursor = e
        if self.hi > cursor:
            width = (self.hi - cursor) if gap is None else min(gap, self.hi - cursor)
            offset += width
            gen.append(self.hi)
            plot.append(offset)

        self.gen_breaks = np.asarray(gen, dtype=float)
        self.plot_breaks = np.asarray(plot, dtype=float)
        self.kept = merged
        self.width = float(self.plot_breaks[-1])

        # A non-monotonic map would silently scramble the picture.
        assert np.all(np.diff(self.gen_breaks) >= 0), "genomic breakpoints not sorted"
        assert np.all(np.diff(self.plot_breaks) >= 0), "plot breakpoints not sorted"

    @staticmethod
    def _merge(intervals):
        out = []
        for s, e in sorted(intervals):
            if out and s <= out[-1][1] + 1:
                out[-1][1] = max(out[-1][1], e)
            else:
                out.append([s, e])
        return [(s, e) for s, e in out]

    def __call__(self, pos):
        return np.interp(np.asarray(pos, dtype=float), self.gen_breaks, self.plot_breaks)

    def inverse(self, x):
        return np.interp(np.asarray(x, dtype=float), self.plot_breaks, self.gen_breaks)

    def ticks(self, max_ticks=10):
        """Tick positions at kept-region boundaries, labelled genomically."""
        edges = sorted({self.lo, self.hi} | {p for seg in self.kept for p in seg})
        if len(edges) > max_ticks:
            idx = np.linspace(0, len(edges) - 1, max_ticks).round().astype(int)
            edges = [edges[i] for i in sorted(set(idx))]
        return list(self(edges)), [f"{int(e):,}" for e in edges]

    @classmethod
    def for_gene(cls, gene, reads=None, region=None, gap=300, pad=200):
        """Build an axis from a gene's annotated exons plus any observed blocks.

        `region` is always kept at full scale even if nothing is annotated
        there, which is the case that matters: an unannotated cryptic exon is
        exactly the thing that must not be compressed away.
        """
        keep = []
        for tx in gene["transcripts"].values():
            keep.extend((s, e) for s, e in tx["exons"])
        if region is not None:
            keep.append((region[1], region[2]))
        if reads is not None and not reads.empty:
            for blocks in reads["blocks"]:
                keep.extend(blocks)
        lo = min([s for s, _ in keep], default=gene["start"])
        hi = max([e for _, e in keep], default=gene["end"])
        return cls(keep, min(lo, gene["start"]), max(hi, gene["end"]), gap=gap, pad=pad)


# --------------------------------------------------------------------------
# Figure 1: the read pileup
# --------------------------------------------------------------------------

def _hover_for(row):
    return (
        f"<b>{row.read_id}</b><br>"
        f"call: {row.xc}<br>"
        f"transcript: {row.xt}<br>"
        f"gene: {row.xg}<br>"
        f"strand: {'same' if row.xd == 'S' else 'opposite' if row.xd == 'O' else 'no call'}"
        f" ({'-' if row.reverse else '+'} read)<br>"
        f"flexed boundaries: {row.xf if row.xf >= 0 else 'n/a'}<br>"
        f"blocks: {row.n_blocks}, span: {int(row.start):,}-{int(row.end):,}"
        "<extra></extra>"
    )


def plot_read_pileup(reads, gene, axis=None, region=None, max_reads_per_sample=150,
                     gap=300, read_px=4.0, seed=0, title=None):
    """Draw every read, one row each, coloured by its nexons classification.

    The gene's transcript models sit in a track at the top, on the same
    (intron shrunk) x axis, so a read can be compared to the thing it was or
    was not assigned to.  Reads are grouped by classification and sorted by
    start position within each group, which makes a block of one colour easy
    to see and, more usefully, makes a block of the *wrong* colour easy to see.

    `region` (chrom, start, end) draws a highlight band, e.g. over a cryptic
    exon.  Above `max_reads_per_sample` the reads are randomly downsampled and
    the sampling is stated in the panel label rather than done silently.
    """
    if reads.empty:
        raise ValueError("No reads to plot")

    reads = reads[reads["start"].notna()].copy()
    axis = axis or ShrunkAxis.for_gene(gene, reads, region=region, gap=gap)
    samples = [s for s in reads["sample"].cat.categories if (reads["sample"] == s).any()]

    transcripts = list(gene["transcripts"].items())

    panels, labels, panel_px = [], [], []
    for sample in samples:
        sub = reads[reads["sample"] == sample]
        note = ""
        if len(sub) > max_reads_per_sample:
            note = (f"  <sup>{max_reads_per_sample:,} of {len(sub):,} reads shown"
                    f" (random sample, seed {seed})</sup>")
            sub = sub.sample(max_reads_per_sample, random_state=seed)
        sub = sub.sort_values(
            ["xc", "start"], key=lambda c: c.cat.codes if hasattr(c, "cat") else c)
        panels.append(sub)
        labels.append(f"<b>{sample}</b>{note}")
        panel_px.append(max(90.0, read_px * len(sub)))

    # Lay the rows out in pixels rather than fractions: a read has to be a few
    # pixels tall to be a read rather than a smear, and the model track needs a
    # floor or its transcript labels land on top of each other.
    model_px = max(70.0, 26.0 * len(transcripts))
    plot_px = model_px + sum(panel_px)
    row_heights = [model_px / plot_px] + [p / plot_px for p in panel_px]

    fig = make_subplots(
        rows=1 + len(samples), cols=1, shared_xaxes=True,
        row_heights=row_heights, vertical_spacing=0.025,
        subplot_titles=[f"<b>gene model</b>: {gene.get('name', gene['id'])}"] + labels,
    )

    # ---- gene model track -------------------------------------------------
    for i, (tx_id, tx) in enumerate(transcripts):
        y = len(transcripts) - 1 - i
        exons = sorted(tx["exons"])
        fig.add_trace(go.Scatter(
            x=axis([exons[0][0], exons[-1][1]]), y=[y, y],
            mode="lines", line=dict(color="#444444", width=1),
            hoverinfo="skip", showlegend=False), row=1, col=1)
        xs, ys = [], []
        for s, e in exons:
            xs.extend(list(axis([s, e])) + [None])
            ys.extend([y, y, None])
        fig.add_trace(go.Scatter(
            x=xs, y=ys, mode="lines", line=dict(color="#222222", width=11),
            name=tx_id, hovertemplate=f"<b>{tx_id}</b><br>{len(exons)} exons<extra></extra>",
            showlegend=False), row=1, col=1)

    fig.update_yaxes(
        row=1, col=1, tickmode="array",
        tickvals=list(range(len(transcripts))),
        ticktext=[t for t, _ in transcripts][::-1],
        range=[-0.6, len(transcripts) - 0.4], showgrid=False)

    # ---- read panels ------------------------------------------------------
    seen_classes = set()
    for panel_i, (sample, sub) in enumerate(zip(samples, panels), start=2):
        sub = sub.reset_index(drop=True)
        for cls in sub["xc"].cat.categories:
            rows_of_class = sub[sub["xc"] == cls]
            if rows_of_class.empty:
                continue
            colour = CLASS_COLOURS.get(cls, "#777777")

            ix, iy = [], []                       # thin intron connectors
            ex, ey, hover = [], [], []            # thick exon blocks
            for y, row in zip(rows_of_class.index, rows_of_class.itertuples()):
                ix.extend(list(axis([row.start, row.end])) + [None])
                iy.extend([y, y, None])
                text = _hover_for(row)
                for s, e in row.blocks:
                    ex.extend(list(axis([s, e])) + [None])
                    ey.extend([y, y, None])
                    hover.extend([text, text, None])

            show = cls not in seen_classes
            seen_classes.add(cls)
            fig.add_trace(go.Scatter(
                x=ix, y=iy, mode="lines", line=dict(color=colour, width=0.7),
                hoverinfo="skip", legendgroup=cls, showlegend=False),
                row=panel_i, col=1)
            fig.add_trace(go.Scatter(
                x=ex, y=ey, mode="lines", line=dict(color=colour, width=2.0),
                name=f"{cls} ({CLASS_MEANING.get(cls, '')})",
                legendgroup=cls, showlegend=show,
                hovertemplate="%{text}", text=hover),
                row=panel_i, col=1)

        fig.update_yaxes(row=panel_i, col=1, showticklabels=False,
                         showgrid=False, range=[-1, max(len(sub), 10)])

    # ---- region highlight -------------------------------------------------
    n_rows = 1 + len(samples)

    if region is not None:
        x0, x1 = axis([region[1], region[2]])
        # Drawn on every row so the band lines up down the whole figure. It is
        # named in the subtitle rather than annotated per row, which would put
        # the same label over every panel.
        for row in range(1, n_rows + 1):
            fig.add_vrect(x0=x0, x1=x1, fillcolor="#F0E442", opacity=0.35,
                          line_width=0, layer="below", row=row, col=1)

    tickvals, ticktext = axis.ticks()
    fig.update_xaxes(range=[0, axis.width], showgrid=False)
    fig.update_xaxes(row=n_rows, col=1, tickmode="array", tickvals=tickvals,
                     ticktext=ticktext, tickangle=-45,
                     title_text=(f"{gene['chrom']}  (introns compressed to {gap} bp)"
                                 if gap is not None else f"{gene['chrom']}  (linear scale)"))

    heading = title or (f"nexons read classifications: "
                        f"{gene.get('name', gene['id'])}, {len(samples)} sample(s)")
    if region is not None:
        heading += (f"<br><sup>shaded band = region of interest, "
                    f"{region[0]}:{region[1]:,}-{region[2]:,}</sup>")

    top_px, bottom_px = 90, 190
    fig.update_layout(
        title=dict(text=heading, x=0, xanchor="left", y=1, yanchor="top",
                   yref="container", pad=dict(t=18, l=10)),
        height=top_px + bottom_px + model_px + sum(panel_px),
        template="plotly_white", hovermode="closest",
        legend=dict(orientation="h", yanchor="bottom", y=0, x=0,
                    xref="paper", yref="container",
                    title_text="", font=dict(size=10)),
        margin=dict(l=170, r=30, t=top_px, b=bottom_px),
    )
    # Subplot titles: left align them over the plotting area.
    for note in fig.layout.annotations:
        if note.text and note.textangle in (None, 0):
            note.font.size = 11
            note.xanchor = "left"
            note.x = 0
    return fig


# --------------------------------------------------------------------------
# Figure 2: where the reads went
# --------------------------------------------------------------------------

def classification_table(reads, normalise=False):
    """Sample x classification counts, ready to print or save as a CSV."""
    table = (reads.groupby(["sample", "xc"], observed=False)
                  .size().unstack(fill_value=0))
    table = table[[c for c in CLASS_ORDER if c in table.columns]
                  + [c for c in table.columns if c not in CLASS_ORDER]]
    if normalise:
        table = 100 * table.div(table.sum(axis=1).replace(0, np.nan), axis=0)
    return table


def plot_classification_composition(reads, normalise=True, title=None):
    """Stacked bars of where each sample's reads ended up, plus the strand split.

    This is the whole-file health check.  A sample whose bar looks different
    from its neighbours differs in something other than biology most of the
    time: a large purple block means the annotation has overlapping genes, a
    large blue block in the strand panel means the library protocol and the
    --direction setting disagree.
    """
    counts = classification_table(reads, normalise=normalise)
    unit = "% of reads" if normalise else "reads"

    strand = (reads.assign(call=reads["xd"].map(
                  {"S": "same strand", "O": "opposite strand"}).fillna("no call"))
                   .groupby(["sample", "call"], observed=False).size().unstack(fill_value=0))
    for column in ("same strand", "opposite strand", "no call"):
        if column not in strand.columns:
            strand[column] = 0
    if normalise:
        strand = 100 * strand.div(strand.sum(axis=1).replace(0, np.nan), axis=0)

    fig = make_subplots(rows=1, cols=2, horizontal_spacing=0.13,
                        subplot_titles=("Read classification",
                                        "Strand call (nexons-assigned reads only)"))

    for cls in counts.columns:
        fig.add_trace(go.Bar(
            x=list(counts.index.astype(str)), y=counts[cls].values, name=cls,
            marker_color=CLASS_COLOURS.get(cls, "#777777"), legendgroup=cls,
            hovertemplate=f"<b>{cls}</b><br>{CLASS_MEANING.get(cls,'')}"
                          f"<br>%{{x}}: %{{y:.1f}} {unit}<extra></extra>"),
            row=1, col=1)

    strand_colours = {"same strand": "#009E73", "opposite strand": "#0072B2",
                      "no call": "#DDDDDD"}
    for column in ("same strand", "opposite strand", "no call"):
        fig.add_trace(go.Bar(
            x=list(strand.index.astype(str)), y=strand[column].values, name=column,
            marker_color=strand_colours[column], legendgroup=column,
            hovertemplate=f"<b>{column}</b><br>%{{x}}: %{{y:.1f}} {unit}<extra></extra>"),
            row=1, col=2)

    fig.update_layout(
        barmode="stack", template="plotly_white",
        title=title or "Where the reads went",
        height=480, margin=dict(l=95, r=95, t=110, b=95),
        legend=dict(orientation="h", yanchor="bottom", y=-0.35, x=0, font=dict(size=10)),
    )
    fig.update_yaxes(title_text=unit, row=1, col=1)
    fig.update_yaxes(title_text=unit, row=1, col=2)
    fig.update_xaxes(tickangle=0, automargin=True)
    return fig


# --------------------------------------------------------------------------
# Figure 3: does the region of interest get included or skipped
# --------------------------------------------------------------------------

def region_support_table(reads, region, min_frac=0.5):
    """Per-sample include / skip counts over `region`, split by evidence tier.

    Only reads that *span* the region are counted, because a read that stops
    short of it is uninformative about it -- it is neither evidence of
    inclusion nor of skipping.  Counting those reads as "skip" is the single
    easiest way to manufacture a splicing result that is really a coverage
    result, so they are reported separately as `not_spanning`.
    """
    chrom, start, end = region
    span = max(1, end - start + 1)
    out = []
    for sample, sub in reads.groupby("sample", observed=True):
        sub = sub[(sub["chrom"] == chrom) & sub["start"].notna()]
        spanning = sub[(sub["start"] <= start) & (sub["end"] >= end)]
        not_spanning = len(sub) - len(spanning)
        row = {"sample": sample, "not_spanning": not_spanning,
               "spanning": len(spanning)}
        for tier, tier_reads in (("unique", spanning[spanning["xc"] == "unique"]),
                                 ("partial", spanning[spanning["xc"] == "partial"]),
                                 ("other", spanning[~spanning["xc"].isin(["unique", "partial"])])):
            inc = skip = 0
            for blocks in tier_reads["blocks"]:
                covered = sum(max(0, min(e, end) - max(s, start) + 1) for s, e in blocks)
                if covered >= min_frac * span:
                    inc += 1
                elif covered == 0:
                    skip += 1
            row[f"include_{tier}"] = inc
            row[f"skip_{tier}"] = skip
        inc_total = sum(row[f"include_{t}"] for t in ("unique", "partial", "other"))
        skip_total = sum(row[f"skip_{t}"] for t in ("unique", "partial", "other"))
        row["include"] = inc_total
        row["skip"] = skip_total
        row["percent_included"] = (100 * inc_total / (inc_total + skip_total)
                                   if (inc_total + skip_total) else np.nan)
        out.append(row)
    return pd.DataFrame(out).set_index("sample")


def plot_region_support(reads, region, min_frac=0.5, title=None):
    """Include vs skip over a region, with the evidence tier kept visible.

    Two samples can show the same inclusion percentage while one is supported
    by full-length `unique` reads and the other only by `partial` ones.  Those
    are not the same result, and a count table cannot tell them apart, so the
    tier is stacked inside each bar rather than summed away.
    """
    table = region_support_table(reads, region, min_frac=min_frac)
    chrom, start, end = region

    fig = make_subplots(rows=1, cols=2, horizontal_spacing=0.16,
                        column_widths=[0.62, 0.38],
                        subplot_titles=(
                            f"Reads spanning {chrom}:{start:,}-{end:,}",
                            "Percent included"))

    tier_shade = {"unique": 1.0, "partial": 0.55, "other": 0.28}
    for outcome, colour in (("include", "#0072B2"), ("skip", "#D55E00")):
        for tier, opacity in tier_shade.items():
            values = table[f"{outcome}_{tier}"].values
            if not values.any():
                continue
            fig.add_trace(go.Bar(
                x=list(table.index.astype(str)), y=values,
                name=f"{outcome} ({tier})", marker_color=colour, opacity=opacity,
                offsetgroup=outcome, legendgroup=outcome,
                hovertemplate=f"<b>{outcome}, {tier} read</b><br>"
                              f"%{{x}}: %{{y}} reads<extra></extra>"),
                row=1, col=1)

    fig.add_trace(go.Bar(
        x=list(table.index.astype(str)), y=table["percent_included"].values,
        marker_color="#000000", opacity=0.75, showlegend=False,
        text=[f"{v:.1f}%" if pd.notna(v) else "n/a" for v in table["percent_included"]],
        textposition="outside",
        hovertemplate="%{x}: %{y:.1f}% of spanning reads include the region<extra></extra>"),
        row=1, col=2)

    fig.update_layout(
        barmode="stack", template="plotly_white",
        title=title or "Inclusion of the region of interest",
        height=490, margin=dict(l=95, r=95, t=110, b=105),
        legend=dict(orientation="h", yanchor="bottom", y=-0.34, x=0, font=dict(size=10)),
    )
    fig.update_yaxes(title_text="spanning reads", row=1, col=1)
    fig.update_yaxes(title_text="% included", range=[0, 109], row=1, col=2)
    fig.update_xaxes(tickangle=0, automargin=True, tickfont=dict(size=11))
    return fig


# --------------------------------------------------------------------------
# Saving
# --------------------------------------------------------------------------

def save_figure(fig, prefix, png=True, scale=2):
    """Write `prefix`.html always, and `prefix`.png if kaleido is installed."""
    written = []
    html_path = Path(f"{prefix}.html")
    html_path.parent.mkdir(parents=True, exist_ok=True)
    fig.write_html(str(html_path), include_plotlyjs="cdn")
    written.append(html_path)
    if png:
        try:
            fig.write_image(f"{prefix}.png", scale=scale)
            written.append(Path(f"{prefix}.png"))
        except Exception as ex:
            print(f"  (no PNG for {prefix}: {ex})", file=sys.stderr)
    return written


# --------------------------------------------------------------------------
# Command line
# --------------------------------------------------------------------------

def get_options(argv=None):
    parser = argparse.ArgumentParser(
        prog="nexons_viz.py",
        description="Visualise nexons read classifications from annotated BAM files.",
        epilog=(
            "The BAMs must carry nexons' XC/XT/XG/XD/XF tags, which means running\n"
            "the quantitation with --annotated-bam first:\n\n"
            "  mkdir -p annotated\n"
            "  python nexons.py --maxtsl 0 --annotated-bam annotated \\\n"
            "      --outbase annotated/nx genes.gtf sample1.bam sample2.bam\n\n"
            "  python nexons_viz.py --gtf genes.gtf --gene MYGENE \\\n"
            "      --bam annotated/sample1.annotated.bam annotated/sample2.annotated.bam \\\n"
            "      --out figures/mygene\n"
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--gtf", required=True,
                        help="GTF of gene models, the same one used for quantitation")
    parser.add_argument("--bam", required=True, nargs="+",
                        help="One or more annotated BAM files")
    parser.add_argument("--name", nargs="+", default=None,
                        help="Sample labels, in the same order as --bam (default: file names)")
    parser.add_argument("--gene", default=None,
                        help="Gene ID or name to draw (default: the first gene in the GTF)")
    parser.add_argument("--region", default=None,
                        help="Region of interest to highlight and quantify, as chr:start-end")
    parser.add_argument("--out", default="nexons_viz",
                        help="Output prefix for the HTML and PNG files")
    parser.add_argument("--maxtsl", type=int, default=0,
                        help="Maximum transcript support level, 0 for no filtering (default 0)")
    parser.add_argument("--max-reads", type=int, default=200,
                        help="Reads drawn per sample before downsampling (default 200)")
    parser.add_argument("--gap", type=int, default=300,
                        help="Compressed intron width in bases; 0 for a true linear axis")
    parser.add_argument("--no-png", action="store_true",
                        help="Only write the interactive HTML")
    return parser.parse_args(argv)


def _pick_gene(genes, wanted):
    if wanted is None:
        return genes[next(iter(genes))]
    if wanted in genes:
        return genes[wanted]
    for gene in genes.values():
        if gene.get("name") == wanted:
            return gene
    raise SystemExit(
        f"Gene '{wanted}' not found. The GTF contains: "
        + ", ".join(list(genes)[:20]) + (" ..." if len(genes) > 20 else ""))


def main(argv=None):
    options = get_options(argv)

    if options.name and len(options.name) != len(options.bam):
        raise SystemExit(f"Got {len(options.name)} --name values for {len(options.bam)} --bam files")

    genes = load_gene_models(options.gtf, options.maxtsl)
    gene = _pick_gene(genes, options.gene)
    region = parse_region(options.region)

    names = options.name or [Path(b).name.replace(".annotated.bam", "").replace(".bam", "")
                             for b in options.bam]
    try:
        reads = load_reads(dict(zip(names, options.bam)), gene=gene)
    except MissingTagsError as ex:
        raise SystemExit(f"\n{ex}\n")

    if reads.empty:
        raise SystemExit(f"No reads found over {gene['chrom']}:{gene['start']}-{gene['end']}")

    print(f"{len(reads):,} reads over {gene.get('name', gene['id'])} "
          f"from {len(names)} sample(s)")

    gap = None if options.gap == 0 else options.gap
    written = []

    written += save_figure(
        plot_read_pileup(reads, gene, region=region, gap=gap,
                         max_reads_per_sample=options.max_reads),
        f"{options.out}_pileup", png=not options.no_png)

    written += save_figure(
        plot_classification_composition(reads),
        f"{options.out}_composition", png=not options.no_png)

    if region is not None:
        written += save_figure(
            plot_region_support(reads, region),
            f"{options.out}_region_support", png=not options.no_png)
        table = region_support_table(reads, region)
        table.to_csv(f"{options.out}_region_support.csv")
        written.append(Path(f"{options.out}_region_support.csv"))

    counts = classification_table(reads)
    counts.to_csv(f"{options.out}_classification_counts.csv")
    written.append(Path(f"{options.out}_classification_counts.csv"))

    for path in written:
        print(f"  wrote {path}")


if __name__ == "__main__":
    main()
