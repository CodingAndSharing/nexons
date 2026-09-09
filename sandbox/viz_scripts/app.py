#!/usr/bin/env python
"""
Streamlit front end for nexons_viz.

Run it with the sandbox pixi environment:

    cd sandbox
    ./.pixi/envs/default/bin/streamlit run viz_scripts/app.py

    # or, once pixi.toml has been locked:
    pixi run streamlit run viz_scripts/app.py

Then open the URL it prints.  The app either loads one of the bundled sandbox
examples or takes uploaded files.  Uploads need the GTF plus, for every BAM,
its .bai index, because pysam reads indexed files from disk rather than from
a buffer -- the uploaded files are written to a temporary directory first.
"""

import sys
import tempfile
from pathlib import Path

import pandas as pd
import streamlit as st

sys.path.insert(0, str(Path(__file__).resolve().parent))
import nexons_viz as nv                                            # noqa: E402

SANDBOX = Path(__file__).resolve().parents[1]
EXAMPLES = SANDBOX / "data" / "output" / "annotated"
INPUTS = SANDBOX / "data" / "input"

st.set_page_config(page_title="nexons read viewer", layout="wide")


# --------------------------------------------------------------------------
# Cached loaders.  The cache key is the file path plus its size and mtime, so
# a re-uploaded or regenerated file is not served from a stale cache.
# --------------------------------------------------------------------------

@st.cache_data(show_spinner="Reading GTF...")
def cached_gene_models(gtf_path, _stamp):
    genes = nv.load_gene_models(gtf_path, max_tsl=0)
    return {gid: g for gid, g in genes.items()}


@st.cache_data(show_spinner="Reading annotated BAMs...")
def cached_reads(bam_map, gene, _stamp):
    return nv.load_reads(dict(bam_map), gene=gene)


def stamp(paths):
    return tuple((str(p), Path(p).stat().st_size, int(Path(p).stat().st_mtime))
                 for p in paths)


def save_uploads(files, directory):
    written = []
    for upload in files:
        target = Path(directory) / upload.name
        target.write_bytes(upload.getbuffer())
        written.append(target)
    return written


# --------------------------------------------------------------------------
# Sidebar: choose the data
# --------------------------------------------------------------------------

st.sidebar.title("nexons read viewer")
st.sidebar.caption("Reads coloured by the classification nexons gave them")

available = sorted(p.name for p in EXAMPLES.glob("example*")) if EXAMPLES.is_dir() else []
source = st.sidebar.radio(
    "Data source",
    (["Bundled sandbox example"] if available else []) + ["Upload my own files"],
    index=0,
)

gtf_path, bam_map = None, {}

if source == "Bundled sandbox example":
    example = st.sidebar.selectbox("Example", available)
    example_dir = EXAMPLES / example
    gtf_candidates = sorted((INPUTS / example).glob("*.gtf"))
    if not gtf_candidates:
        st.error(f"No GTF found in {INPUTS / example}. "
                 "Run the data generator and the nexons annotation step first.")
        st.stop()
    gtf_path = gtf_candidates[0]
    bams = sorted(example_dir.glob("*.annotated.bam"))
    if not bams:
        st.error(
            f"No annotated BAMs in {example_dir}.\n\n"
            "Generate them with:\n\n"
            "```\nmkdir -p data/output/annotated/EXAMPLE\n"
            "cd data/input/EXAMPLE && python ../../../../nexons.py --maxtsl 0 \\\n"
            "    --annotated-bam ../../output/annotated/EXAMPLE \\\n"
            "    --outbase ../../output/annotated/EXAMPLE/nx *.gtf *.bam\n```")
        st.stop()
    chosen = st.sidebar.multiselect(
        "Samples", [b.name for b in bams], default=[b.name for b in bams])
    bam_map = {b.name.replace(".annotated.bam", ""): b for b in bams if b.name in chosen}

else:
    st.sidebar.markdown("**Upload a GTF and your annotated BAMs**")
    gtf_upload = st.sidebar.file_uploader("GTF", type=["gtf", "gff", "txt"])
    bam_uploads = st.sidebar.file_uploader(
        "Annotated BAMs and their .bai indexes", type=["bam", "bai"],
        accept_multiple_files=True)
    if not gtf_upload or not bam_uploads:
        st.info("Upload a GTF and at least one annotated BAM (plus its .bai) to begin.")
        st.stop()

    workdir = Path(tempfile.mkdtemp(prefix="nexons_viz_"))
    gtf_path = workdir / gtf_upload.name
    gtf_path.write_bytes(gtf_upload.getbuffer())
    saved = save_uploads(bam_uploads, workdir)
    bams = [p for p in saved if p.suffix == ".bam"]
    if not bams:
        st.error("No .bam files among the uploads.")
        st.stop()
    missing = [b.name for b in bams if not (b.with_suffix(".bam.bai").exists()
                                            or b.with_suffix(".bai").exists())]
    if missing:
        st.warning("No index uploaded for: " + ", ".join(missing)
                   + ". These will be read end to end, which is slower.")
    bam_map = {b.name.replace(".annotated.bam", "").replace(".bam", ""): b for b in bams}

if not bam_map:
    st.info("Select at least one sample.")
    st.stop()

# --------------------------------------------------------------------------
# Load
# --------------------------------------------------------------------------

genes = cached_gene_models(str(gtf_path), stamp([gtf_path]))
gene_labels = {f"{g.get('name', gid)}  ({gid})": gid for gid, g in genes.items()}
gene_label = st.sidebar.selectbox("Gene", sorted(gene_labels))
gene = genes[gene_labels[gene_label]]

try:
    reads = cached_reads(tuple((k, str(v)) for k, v in bam_map.items()), gene,
                         stamp(list(bam_map.values())))
except nv.MissingTagsError as ex:
    st.error(str(ex))
    st.stop()

if reads.empty:
    st.warning(f"No reads over {gene['chrom']}:{gene['start']:,}-{gene['end']:,}")
    st.stop()

st.sidebar.markdown("---")
default_region = f"{gene['chrom']}:{gene['start']}-{gene['end']}"
region_text = st.sidebar.text_input(
    "Region of interest (chr:start-end)", value="",
    placeholder=default_region,
    help="A cryptic exon, a junction, any window you want include/skip counts for")
region = None
if region_text.strip():
    try:
        region = nv.parse_region(region_text)
    except ValueError as ex:
        st.sidebar.error(str(ex))

present = [c for c in reads["xc"].cat.categories]
keep_classes = st.sidebar.multiselect("Show classifications", present, default=present)
max_reads = st.sidebar.slider("Max reads drawn per sample", 100, 5000, 1500, step=100)
shrink = st.sidebar.checkbox("Compress introns", value=True,
                             help="Off gives a true linear genomic axis")
gap = 300 if shrink else None

filtered = reads[reads["xc"].isin(keep_classes)]
if filtered.empty:
    st.warning("Every read was filtered out. Re-enable a classification.")
    st.stop()

# --------------------------------------------------------------------------
# Main panel
# --------------------------------------------------------------------------

st.title(gene.get("name", gene["id"]))
top = st.columns(4)
top[0].metric("Reads shown", f"{len(filtered):,}")
top[1].metric("Samples", f"{filtered['sample'].nunique()}")
top[2].metric("Transcripts in model", f"{len(gene['transcripts'])}")
top[3].metric("Unique-assigned", f"{(filtered['xc'] == 'unique').mean() * 100:.1f}%")

tabs = st.tabs(["Read pileup", "Composition", "Region support", "Table"])

with tabs[0]:
    fig = nv.plot_read_pileup(filtered, gene, region=region, gap=gap,
                              max_reads_per_sample=max_reads)
    st.plotly_chart(fig, width="stretch")
    st.download_button(
        "Download this figure as HTML",
        fig.to_html(include_plotlyjs="cdn"),
        file_name=f"{gene.get('name', gene['id'])}_pileup.html",
        mime="text/html")
    st.caption(
        "Each row is one read. Thick blocks are aligned exon blocks, thin lines "
        "are the gaps the aligner called as introns. Colour is the call nexons "
        "made about that read; the models at the top are what it was compared to.")

with tabs[1]:
    normalise = st.checkbox("Show as percentages", value=True)
    st.plotly_chart(nv.plot_classification_composition(filtered, normalise=normalise),
                    width="stretch")
    st.caption(
        "A sample whose bar looks unlike its neighbours usually differs in "
        "something other than biology. Purple means the annotation has "
        "overlapping genes; a large opposite-strand block means the library "
        "protocol and the --direction setting disagree.")

with tabs[2]:
    if region is None:
        st.info("Enter a region of interest in the sidebar to get include vs skip counts.")
    else:
        st.plotly_chart(nv.plot_region_support(filtered, region),
                        width="stretch")
        table = nv.region_support_table(filtered, region)
        st.dataframe(table, width="stretch")
        st.caption(
            "Only reads that span the whole region are counted, because a read "
            "stopping short of it is evidence of neither inclusion nor skipping. "
            "The unique/partial split matters: the same percentage carried by "
            "partial reads is a weaker claim than one carried by full-length reads.")

with tabs[3]:
    counts = nv.classification_table(filtered)
    st.dataframe(counts, width="stretch")
    st.download_button("Download counts as CSV", counts.to_csv(),
                       file_name="classification_counts.csv", mime="text/csv")
    with st.expander("What the classifications mean"):
        st.table(pd.DataFrame(
            [{"classification": c, "meaning": nv.CLASS_MEANING[c]} for c in nv.CLASS_ORDER]))
