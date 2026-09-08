## HAckathon challenge: 

**Goal**

The nexons program in this project allows you to align nanopore RNA-Seq reads to a reference and quantify genes or transcripts with different levels of precision. Later on you can use the count table this produces for differential expression. What we would like to have is a viewer which would take in a set of annotated BAM files from nexons and would provide a comprehensive view of the read alignments for multiple samples over and around that gene, showing the classifications which have been made. This will allow the user to build much more confidence in candidates which were selected from differential expression.


**Questions**

1. Does exist already any tool online that perform this task (github, opensource or commercial)? if so, list them below with pros and cons.

2. what is the advantage of having such visualization? explain it with a use case in ALS area

---

# Answers

## 0. Prerequisite: nexons does not currently produce annotated BAM files

The goal statement says the viewer "would take in a set of annotated BAM files
from nexons". That input does not exist yet, and building it is the first
deliverable — everything in §1 and §2 below assumes it.

Verified against `nexons.py` on this branch:

- The only pysam call in the whole program is
  `pysam.AlignmentFile(bam_file, "rb")` (line 326). nexons opens BAMs
  **read-only** and never writes one.
- Its outputs are three count tables (`_unique.txt`, `_partial.txt`,
  `_gene.txt`), a per-BAM JSON stats file, and a per-BAM QC HTML report.
- The per-read classification **is** computed — every read passes through a
  decision tree that assigns it to one of the `outcomes` counters — but the
  result is immediately folded into an integer and discarded. Nothing
  downstream can recover which read produced which count.

The classification vocabulary already in the code is exactly the payload a
viewer would need:

| Category | Values in `nexons.py` |
|---|---|
| Alignment fate | `Total_Reads`, `No_Alignment`, `Secondary_Alignment`, `Primary_Alignment`, `No_Gene`, `No_Hit`, `Multi_Gene` |
| Assignment level | `Unique`, `Partial`, `Gene` (the three output "slots") |
| Directionality | `Same_Strand_Hit`, `Opposing_Strand_Hit` |

**Proposed first task:** add an `--annotated-bam` / `-a` option that streams the
input BAM back out with per-read tags. Using the `X?`/`Y?`/`Z?` space that the
SAM specification reserves for local use:

| Tag | Type | Content |
|---|---|---|
| `XC` | `Z` | classification — `unique` \| `partial` \| `gene` \| `multi_gene` \| `no_hit` \| `no_gene` |
| `XT` | `Z` | assigned transcript ID, or `.` when not transcript-resolved |
| `XG` | `Z` | assigned gene ID (comma-separated when `multi_gene`) |
| `XD` | `A` | directionality — `S` (same strand) / `O` (opposing) |
| `XF` | `i` | number of exon boundaries that only matched within `--flex`/`--endflex` tolerance |

`XF` is the one that is not in the current counters but is worth adding: it
records how much slack the match needed, which is precisely the thing a human
wants to interrogate when a call looks marginal.

This is a small change (one output handle, one `set_tag` block in the existing
per-read loop) and it is what turns nexons from a counter into something a
viewer — including off-the-shelf ones, see §1 — can render.

---

## 1. Does a tool already exist that does this?

**Short answer: no single tool does the whole job, but the job decomposes into
four capabilities that existing tools each cover part of.** The specific
combination the goal asks for — *per-read classifications from a third-party
quantifier, displayed read-level, across multiple samples, over one gene* — is
not, as far as I can find, provided end to end by anything on GitHub or
commercially. The closest off-the-shelf approximation is IGV with tag-based
grouping, which is worth knowing about because it sets the bar the hackathon
project has to clear.

### A. General genome browsers — read-level, multi-sample, tag-aware

| Tool | Pros | Cons |
|---|---|---|
| **IGV** (desktop) | The closest existing answer. Loads N BAMs as N stacked panels over one locus; **group / colour / sort alignments by BAM tag** is built in, so the moment nexons emits `XC` tags you get most of the requested view for free. Handles spliced long reads, junction track, Sashimi plots. Universally installed; reviewers already trust it. | Knows nothing about nexons' semantics — you re-do the group-by-tag setup by hand every session (session XML helps, but it's manual). No cross-sample summary panel, no link back to the QC report, no per-classification counts. Introns dominate the horizontal space, so a 100 kb gene with a 128 bp cryptic exon is mostly empty. Desktop app; awkward to share a view with a collaborator. |
| **igv.js** | Same rendering model, embeddable in a web page or notebook — the obvious substrate to *build on* rather than compete with. Scriptable, so grouping/colouring can be preset instead of clicked. | A library, not a product: everything above the track (layout, summary panels, sample metadata, classification legend) is yours to write. |
| **JBrowse 2** | Modern plugin architecture, genuinely multi-sample, good long-read support, web-native and shareable by URL. | Heavier to deploy; a custom read-colouring scheme means writing a plugin. |
| **UCSC Genome Browser** | Ubiquitous, trivially shareable. | Not designed for read-level inspection of many BAMs; custom-track upload is clumsy for this. |

### B. Sashimi / junction-level plotters

`ggsashimi`, `rmats2sashimiplot`, `sashimi.py`, `trackplot`, `pyGenomeTracks`,
`Gviz`.

- **Pros:** built for exactly the multi-sample, one-locus layout the goal
  describes; publication-quality static output; junction arcs with read
  support numbers make a splicing change immediately legible.
- **Cons:** they are **junction-centric, not read-centric**. A junction count
  cannot tell you whether the supporting read was a full-length transcript or a
  200 bp truncated fragment that happens to cross that junction — which is the
  entire distinction nexons' `unique` / `partial` split encodes. None of them
  has a concept of a per-read classification tag to overlay. Mostly static
  images, so no drill-down.

### C. Isoform-structure viewers (reads collapsed into isoform rows)

| Tool | Notes |
|---|---|
| **IsoVis** (webserver, NAR 2024, `isomix.org/isovis`) | Visualizes isoform structures + expression, and uniquely overlays ORFs and protein domains so the functional consequence of a structural change is visible. Data stays private to the user; publication-ready export; no programming needed. **Con:** takes isoform models + a count matrix, not BAMs — individual reads are never shown, which is the specific thing the goal asks for. |
| **ScisorWiz** (Bioinformatics 2022) | Renders reads for one gene across any number of cell types, clustered by intron chain / TSS / polyA site, with intron shrinking and SNV marks. Closest existing tool *in visual grammar* to what's wanted. **Con:** built for single-cell long-read pipelines (scisorseqr-shaped input), R command-line, static output, groups by cell type rather than by an arbitrary per-read tag. |
| **Swan** (`swan_vis`) | Graph-based transcriptome representation with DE and differential-isoform-usage built in. **Con:** abstracts away from reads to a splice graph; interpretation is graph-shaped, not alignment-shaped. |
| **isoformant** | Reference-free, single-read resolution, k-mer/UMAP clustering of reads with consensus alignment tracks. Genuinely read-level. **Con:** reference-free clustering answers "what read populations exist", not "how did the annotation-based classifier label each read". |
| **isoespy** (Bioinf. Advances 2026), **IsoTV**, **Iso-Seq Browser**, **MatchAnnot**, **IsoView** | Integrated long-read workflows / older exon-only viewers in the same family. Same structural limitation: transcript models in, not classified reads. |

### D. Cross-sample isoform aggregation (no read view)

`isoSeQL` (Bioinformatics 2026), `SQANTI3`, `gffcompare`, `TAMA merge`,
`cDNA_Cupcake chain_samples.py`.

- **Pros:** `isoSeQL` is purpose-built for the "compare isoform profiles across
  many samples" problem that the older merge tools handle poorly — the merge
  tools are known to be limited to small sample numbers, can produce redundant
  isoform links, and don't track sample metadata. SQANTI3's FSM / ISM / NIC /
  NNC classification is conceptually the closest published analogue to nexons'
  `unique` / `partial` / `gene`.
- **Cons:** all operate on collapsed transcript models and count matrices. They
  are the layer *above* the one the goal is asking for — they will tell you
  the candidate is interesting, not whether you should believe it.

### E. Commercial / vendor

ONT **EPI2ME** (`wf-transcriptomes`), PacBio **SMRT Link**, **Partek Flow**,
**CLC Genomics Workbench**, **Geneious Prime**.

- **Pros:** integrated, supported, non-programmer-friendly, include alignment
  browsers.
- **Cons:** closed pipelines built around their own classifier and their own
  intermediate formats. None of them will ingest a third-party per-read
  classification from nexons; adapting them is not possible from outside.

### Verdict for the hackathon

The gap is real but narrow, and that is good news — it means the project is
tractable in a hackathon and has a defensible reason to exist. Concretely:

1. **Do not rebuild an alignment renderer.** Build on `igv.js` or JBrowse 2.
2. **The differentiator is the classification layer, not the pileup**: reads
   coloured and grouped by nexons' `XC`/`XD` tags, with per-sample stacked
   panels sharing one coordinate axis.
3. **The second differentiator is intron shrinking** (borrow the idea from
   ScisorWiz / the exon-only browsers) — without it a 128 bp cryptic exon in a
   ~90 kb gene is invisible at whole-gene zoom.
4. **The third is closing the loop to the count table**: click a transcript row
   in `nexons_output_unique.txt`, land on the reads that produced that number.
   No existing tool can do this, because no existing tool knows about nexons'
   tables.
5. **Ship the tag writer upstream** (§0). It is independently useful — with
   `XC` tags in the BAM, plain IGV becomes a usable fallback viewer, which is a
   good hedge if the hackathon runs short of time.

---

## 2. What is the advantage of such a visualization? (with an ALS use case)

### 2.1 The general argument

Differential expression gives you a row in a table and a p-value. It does not
tell you *why* that row moved, and the count tables nexons produces are exactly
the kind that can move for reasons that have nothing to do with biology. A
read-level, classification-aware, multi-sample view is the tool that
distinguishes a real isoform change from each of these:

1. **Truncation vs. true isoform switch.** Nanopore reads truncate — degraded
   RNA, non-processive reverse transcription. Truncated reads land in
   nexons' `partial` slot. A candidate whose signal comes from `partial` reads
   in one group and `unique` reads in another is a library-quality artefact
   wearing a biology costume. Only a read-level view colour-coded by
   classification makes that visible at a glance.
2. **Gene-level flat, transcript-level moving (and vice versa).** The demo in
   `sandbox/` shows this directly: gene-level counts identical at 870 across all
   three samples while the transcript split swings from 437/23 to 161/299. If
   your DE was run on the gene table you would have seen nothing.
3. **Multi-gene ambiguity.** Overlapping or antisense genes push reads into
   `Multi_Gene`, where they are counted for neither. A viewer shows you that the
   "low expression" of your candidate is actually 50% of its reads being
   discarded as ambiguous — `example2` in the sandbox reproduces this exactly.
4. **Library directionality errors.** A mis-specified `--direction` silently
   halves or scrambles counts; `example3` reproduces a 33% opposing-strand
   population. Invisible in a count table, obvious in a stranded read view.
5. **Boundary tolerance.** nexons merges exon boundaries within `--flex` /
   `--endflex`. A candidate whose classification flips when those are tightened
   is not a candidate. The proposed `XF` tag makes this inspectable.
6. **Structures the annotation does not contain.** This is the deepest one. A
   count table is bounded by the GTF: a read matching a genuinely novel
   junction has nowhere to go but `gene` or `no_hit`. **The viewer is the only
   place where an unannotated structure can be noticed at all.**

### 2.2 ALS use case: the *UNC13A* cryptic exon

**Background.** TDP-43 is normally nuclear and acts as a repressor of cryptic
exon inclusion; in ALS and FTD it is depleted from the nucleus and aggregates
in the cytoplasm. In *UNC13A* — one of the strongest GWAS hits for ALS/FTD —
nuclear TDP-43 loss de-represses a cryptic exon between canonical exons 20 and
21, introducing a premature termination signal, sending the transcript to
nonsense-mediated decay and reducing UNC13A protein. The top ALS/FTD risk
variants sit *inside the cryptic-exon-harbouring intron* and increase inclusion,
and splice-switching antisense oligonucleotides against the cryptic exon rescue
UNC13A protein and synaptic function in model systems — which is why several
ASO programmes now target it.

Two details from that literature make this the ideal case for the proposed
viewer:

- There is not one cryptic exon but **two**, 128 bp and 178 bp, sharing the same
  3′ end and differing by 50 bp at the 5′ end (hg38 chr19:17642414–17642541 for
  the 128 bp form).
- The cryptic transcript is an **NMD substrate**, so it is depleted before you
  ever sequence it.

**Scenario.** Nanopore direct-RNA or cDNA sequencing of iPSC-derived motor
neurons ± TDP-43 knockdown, or post-mortem motor cortex from ALS cases vs.
controls, n = 4 per group. nexons quantifies against GENCODE; DE on the
transcript table flags a cryptic-exon-containing *UNC13A* transcript as up in
the ALS group. What does the viewer add?

1. **It resolves which cryptic exon.** The 128 bp and 178 bp forms differ only
   at the 5′ acceptor. Any gene-level count merges them; a junction-level
   Sashimi plot distinguishes them only if the annotation happens to contain
   both. A read-level view over the exon 20–21 intron shows, per read, which 5′
   start was used — and whether the split between them tracks the risk
   haplotype carried by each donor. That is a biological result, not a QC step.
2. **It separates "low inclusion" from "NMD ate the evidence."** Because
   cryptic-containing transcripts are degraded, they are both rarer and more
   often truncated. In the count table those two causes are indistinguishable.
   In the viewer, a `unique`-rich cryptic population reads very differently from
   a `partial`-only smear.
3. **It controls for post-mortem RNA quality.** ALS post-mortem tissue varies in
   RIN, and 3′ bias differs per sample. Read-start/end distributions stacked
   across samples answer "is the cryptic call confounded by one sample simply
   having longer reads?" — a question a reviewer *will* ask.
4. **It is the only readout that can catch an off-target ASO product.** If the
   experiment is testing a splice-switching ASO, success is "cryptic exon
   skipped, canonical exon 20–21 junction restored". But an ASO can also induce
   an unintended product — skipping of exon 20 or 21 itself, or use of a third
   cryptic acceptor. Those isoforms are **not in the GTF**, so nexons will
   classify those reads as `gene` or `no_hit` and the count table will report a
   clean success. The viewer is where you see the reads that stopped agreeing
   with the annotation.
5. **It survives contact with a collaborator.** "Cryptic *UNC13A* is up,
   adjusted p = 0.003" is an assertion. A four-panel figure showing every
   supporting read, coloured by classification, over the exon 20–21 intron, is
   evidence — which is precisely the "build much more confidence in candidates"
   the goal statement asks for.

### 2.3 Ready-made development data

The three synthetic datasets already in `sandbox/data/input/` are a usable test
bed for the viewer, and each was built to exercise a different display path:

| Dataset | What it exercises |
|---|---|
| `example1` | The clean case — canonical/cryptic switch across `control`, `tdp43_depleted`, `mixed`; gene-level flat, transcript-level moving |
| `example2` | `Multi_Gene` rendering — an overlapping paralog puts 50% of primary alignments into ambiguous assignment |
| `example3` | Strand rendering — ~33% of reads aligned antisense |

They are small, deterministic, aligner-free and carry no real patient data, so
they can be committed and used in CI for the viewer. The locus is a synthetic
toy (`UNC13A_TOY` on `UNC13A_TOY_chr`) with made-up coordinates modelled on the
real exon-20/21 architecture — it is not real *UNC13A* sequence.

### 2.4 Scope note

Everything above concerns research tooling. Nothing in this repository is
validated for clinical or diagnostic use, and no conclusion about any real
sample's TDP-43 status follows from the synthetic demo data.

---

## Key references

- Ma XR, Prudencio M, Koike Y, *et al.* TDP-43 represses cryptic exon inclusion
  in the FTD–ALS gene *UNC13A*. *Nature* 603, 124–130 (2022).
  doi:10.1038/s41586-022-04424-7
- Koike Y, Pickles S, Estades Ayuso V, *et al.* TDP-43 and other hnRNPs regulate
  cryptic exon inclusion of a key ALS/FTD risk gene, *UNC13A*. *PLoS Biol* 21,
  e3002028 (2023). doi:10.1371/journal.pbio.3002028
- Loss of TDP-43 induces synaptic dysfunction that is rescued by *UNC13A*
  splice-switching ASOs. PMID 38979232 (2024).
- Wan CY, Davis J, Choi J, Clark M. IsoVis — a webserver for visualization and
  annotation of alternative RNA isoforms. *Nucleic Acids Res* 52, W341 (2024).
  doi:10.1093/nar/gkae343 · https://isomix.org/isovis/
- Stein AN, *et al.* ScisorWiz: visualizing differential isoform expression in
  single-cell long-read data. *Bioinformatics* (2022).
  https://github.com/ans4013/ScisorWiz
- Liu C, *et al.* isoSeQL: comparing long-read isoforms across multiple
  datasets. *Bioinformatics* 42, btaf680 (2026).
  https://github.com/christine-liu/isoSeQL
- Reese F, Mortazavi A. Swan: a library for the analysis and visualization of
  long-read transcriptomes. *Bioinformatics* (2021).
- Prjibelski AD, *et al.* Accurate isoform discovery with IsoQuant using long
  reads. *Nat Biotechnol* 41, 915–918 (2023). doi:10.1038/s41587-022-01565-y
