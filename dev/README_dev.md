# Nexons — Developer Setup & Worked Example (ALS / UNC13A cryptic-exon use case)

> Status: **planning document**. This describes the environment and the worked
> example that will be built in `sandbox/`. Nothing here has been executed yet —
> it is the spec to review before the pixi environment, dummy dataset, and
> notebook are actually created.

---

## 1. What Nexons does

Nexxons is a command-line Python application used for the quantification and quality control (QC) of long-read RNA sequencing data. [1] (https://www.youtube.com/watch?v=UGmxsFjjyu0&t=660s)

Input: Starts from BAM files of aligned long-read RNA-Seq.
Provide GTP file (Works best with Ensembl GTF)

![img](/images/nexons_cli.png)

> When Oxford NAnopore devices sequence full-length RNA molecules, they generate raw sequence reads. Scientist align these reads to a reference genome to create a alignment file (a BAM file)*. Nexons takes that BAM file alongside a gene reference file (GTF) to perform 2 critical tasks:
> - **Transcript and Gene Quantification**: IT counts how many times each specific RAN transcript or gene appears in your sample. Ir accounts for the slight per-base inaccuracies that can happen in long-read sequencing using flexible matching for splice sites and transcript ends. And it outputs 3 tables:
>   - Full length transcript matches
>   - Partial but unique transcript matches
>   -  Unique gene-level matches
> In the image below, you have A schematic for a single gene that will be quantitated. and there are   3 different isoforms for that gene that appear in the GTF file at the gene build. (this time is not included any structure of the transcripts but it is included the annotated transcript support level). the first 2 isoforms have TLS1 (very well supported, very well supported) but TLS5 has very little support level (very little evidence this is a real transcript. maybe sketchy isoform that has made it into the datasets).
the lines below in grey and red are reads, matching the data of isoforms. the matches to the reference isoforms are in red.there could be more flexibility at the end, the start or in the middle. you can see a clear match (Full match) or partial match - unique (potentially matching one isoform, not the others). or a read that is compatible with two isoforms, but there is no enough context to distingush between two isoforms, then it will be gene level annotation (because we cant be sure that it comes from this gene isoform or the other).  
> ![img](/images/transcripts_isoforms.png)
> If you apply a 'transcript support level filter (lets include transcripts for levels 1 and 2), then one isoform (TLS5) is excluded, and the final match of the read it will become a UNique partial. 
> ![img](/images/transcripts_isoforms_withexcludedisoform.png)
> - Quality control (QC)**: It automatically generates and interactive HTML QC report that gives researchers key metrics about how well their long-read RNA data aligned and mapped to the genome. 
>   - breakdown of the data in the bam file ![readfate](/images/read_fate.png). the ...
>   - alignment fate ![alignmentfate](/images/alignment_fate.png)
>   - alignment directionality ![alignmentdirectionality](/images/alignment_directionality.png)
>   - read lengths and coverage ![read_length_and_coverage](/images/read_length_and_coverage.png)
>   - exon match flexibility ![exonmatchfelixibility](/images/exon_flex.png)

_* is this using certified software and methods?_
_** alignment fate, why these labels? why no colors (POTENTIAL IMPROVEMENTS), why the expected aligment directionalty is 50/50, or 100/0 or 0/100? why read_lengths picture has short length reads. are those as short as with illumina? wich other datasets have a different distribution? Transcript coverage why is expected a fall at the beginning, at the end. for which percentage in all cells in a normal/standard human genes is this expected? about exon match flexibility what is inner exon flex peak there, what is the stats behind. same for next plot_
_*** how good is this quality control method compare to others in the market?_

**In summary**
Nexons quantifies Oxford Nanopore long-read RNA-seq. 
Given a spliced-aligned BAM and a GTF of gene models, it walks each read's exon/intron structure (from its CIGAR string) and asks: does this read match one annotated transcript completely (`unique`), match one transcript but only partially (`partial`), or only resolve to a gene without pinning a transcript (`gene`)? The three resulting count tables, plus a per-BAM QC HTML report, are its output.

## 2. **Use case**
The property that matters for the use case below: because a nanopore read
typically spans a large fraction of a transcript in one piece, nexons can
assign a read to a *specific splice isoform* rather than just an exon-junction count. That is exactly the resolution needed to separate a normal transcript from one carrying a disease-associated cryptic exon.


### 2.1 Scientific background: TDP-43 loss-of-function and the UNC13A cryptic exon in ALS

> *This section is background for interpreting the demo output — it is not
medical advice, and none of the data used below comes from a real patient.
Any decision about diagnosis or treatment must be made by a qualified
healthcare professional with access to the full clinical picture.*

Amyotrophic lateral sclerosis (ALS) is a fatal neurodegenerative disease of motor neurons. In the great majority of ALS cases (and in the related disorder frontotemporal dementia, FTD), the RNA-binding protein **TDP-43** mislocalizes from the nucleus to the cytoplasm and forms aggregates, depleting the nucleus of functional TDP-43. Because TDP-43 normally binds intronic UG-rich sequences to *suppress* the inclusion of certain non-conserved "cryptic" exons, its nuclear loss lets those exons slip into the mature mRNA of specific target genes.

Two genes have become the leading RNA biomarkers of this mechanism:

- **`STMN2`** — cryptic-exon-containing transcripts truncate the message and reduce functional stathmin-2, a protein needed for axon regeneration.
- **`UNC13A`** — TDP-43 depletion allows inclusion of a cryptic exon between   the (real, hg38) exons 20 and 21 of the transcript, disrupting normal splicing. `UNC13A` is independently an ALS/FTD genetic risk locus, so this   gene sits at the intersection of GWAS risk and the TDP-43 splicing mechanism, which is why it is one of the most studied nexons-style targets in the ALS field.

Because cryptic-exon transcripts are typically low-abundance, alternatively polyadenylated, and NMD-sensitive, they are easy to under-count with short-read/junction-counting pipelines. A long-read tool like nexons that classifies whole reads against whole transcript models is a natural fit for quantifying *how much* of the transcript pool has shifted from the canonical isoform to the cryptic-exon isoform.

#### 2.1.1 The worked example

Because we don't have access to real patient sequencing data, the demo uses a **synthetic toy locus**, `UNC13A_TOY`, that reproduces the structure of the real biology (two transcripts differing by one extra exon, at the same relative position TDP-43 studies report) on made-up coordinates. This keeps the example fast, dependency-free (no genome download), and exactly reproducible, while still being biologically legible.

#### 2.1.2 Input dataset (`sandbox/data/`)

| File | Role |
|---|---|
| `unc13a_toy.gtf` | Gene model with **two transcripts** of gene `UNC13A_TOY`: `UNC13A_TOY-CANONICAL` (20 exons, no cryptic exon) and `UNC13A_TOY-CRYPTIC` (the same 20 exons plus one extra cryptic exon inserted between exon 20 and exon 21, matching the reported real-gene position). |
| `control.bam` | Synthetic nanopore-style reads representing **healthy motor neurons** (TDP-43 nuclear, functional): every read is built to splice out the cryptic exon, i.e. matches `UNC13A_TOY-CANONICAL`. |
| `tdp43_depleted.bam` | Synthetic reads representing **ALS-model motor neurons with nuclear TDP-43 loss**: every read retains the cryptic exon, i.e. matches `UNC13A_TOY-CRYPTIC`. |
| `mixed.bam` *(optional third sample)* | A blend of both isoforms, representing a heterogeneous tissue biopsy, to show fractional cryptic-exon inclusion rather than an all-or-nothing switch. |

Both BAMs will be built with `pysam` directly (no aligner needed, since we
already know the intended splice sites) against a single-contig synthetic
reference (`UNC13A_TOY_chr`), with proper spliced CIGAR strings (`N` operations
at intron positions) so nexons' exon-block parser (`get_exons()`) sees them
exactly as it would see real spliced nanopore alignments.

#### 2.1.3 Command

Run from the pixi environment (see §4), from the `sandbox/` directory:

```bash
pixi run python ../nexons.py \
    data/unc13a_toy.gtf \
    data/control.bam data/tdp43_depleted.bam \
    --maxtsl 0 \
    --outbase data/results/nexons_output
```

(`--maxtsl 0` disables the transcript-support-level filter, which is an
Ensembl/GENCODE annotation-quality field that our toy GTF doesn't populate.)

This will run identically inside the Jupyter notebook (§5), either as a `!`
shell command or via `subprocess.run(...)`, using the pixi-environment kernel.

#### 2.1.4 Expected outputs and how to read them

All files land under `data/results/` with the prefix `nexons_output`:

- **`nexons_output_unique.txt`** — one row per transcript, one column per BAM.
  Columns: `Transcript_ID, Gene_ID, Gene_Name, Chr, Start, End, Strand,
  control.bam, tdp43_depleted.bam`. This is the headline table for the demo.
  Expected pattern:

  | Transcript_ID | ... | control.bam | tdp43_depleted.bam |
  |---|---|---|---|
  | UNC13A_TOY-CANONICAL | ... | high count | ≈0 |
  | UNC13A_TOY-CRYPTIC | ... | ≈0 | high count |

  This is the computational readout of the biology: reads shift from the
  canonical transcript to the cryptic-exon transcript as TDP-43 function is
  lost — the same directional shift reported for `UNC13A` cryptic-exon
  inclusion in TDP-43-depleted neurons and post-mortem ALS tissue.

- **`nexons_output_partial.txt`** — same shape, but also counts reads that
  matched a transcript over only part of its length (e.g. a read that only
  covers the exons flanking the cryptic exon without spanning it fully).
  Superset of `unique.txt`'s counts.

- **`nexons_output_gene.txt`** — collapses both transcripts to the
  `UNC13A_TOY` gene level. Because both isoforms belong to the same gene, this
  table will show comparable *total* signal in both samples — demonstrating
  why gene-level quantitation alone would hide the cryptic-exon switch, and
  why nexons' transcript-level resolution is the point.

- **`nexons_output_control_qc.html` / `nexons_output_tdp43_depleted_qc.html`**
  — per-sample QC reports (read-length histogram, alignment/read-fate
  breakdown, directionality, exon-boundary "flex" observations). Used to
  sanity-check that reads were consumed as intended (e.g. `Unique` outcome
  count equal to the number of synthetic reads generated, `No_Gene`/`No_Hit`
  ≈ 0) rather than to derive the biological conclusion itself.

- **`nexons_output_<sample>_stats.txt`** — the JSON backing each QC HTML
  report, useful for scripted checks in the notebook (e.g. asserting
  `outcomes["Unique"] == n_reads_generated`).

#### 2.1.5. Installing the environment with `pixi`

Everything below runs inside `sandbox/`, so the environment is self-contained
and does not touch anything outside the repo.

1. Install the `pixi` CLI (one-time, machine-level — see
   <https://pixi.sh/latest/#installation>):

   ```bash
   curl -fsSL https://pixi.sh/install.sh | sh
   ```

2. From `sandbox/`, initialize the project and add the dependencies nexons
   and the notebook need:

   ```bash
   cd sandbox
   pixi init .
   pixi add python=3.12 pysam jupyterlab ipykernel pandas matplotlib
   ```

   This writes `sandbox/pixi.toml` (the manifest) and `sandbox/pixi.lock`
   (the resolved, reproducible lock file), and materializes the actual
   environment under `sandbox/.pixi/envs/default`.

3. Register that environment as a Jupyter kernel, installed *inside the pixi
   env's own prefix* (not the user's home directory) so it is auto-discovered
   whenever JupyterLab is launched from this same environment:

   ```bash
   pixi run python -m ipykernel install \
       --prefix .pixi/envs/default \
       --name nexons-sandbox \
       --display-name "Python (nexons pixi)"
   ```

#### 2.1.6. Running the example in JupyterLab

```bash
cd sandbox
pixi run jupyter lab
```

Open `sandbox/nexons_als_demo.ipynb` (to be created in the implementation
step) and pick the **"Python (nexons pixi)"** kernel. The notebook will:

1. Build `data/unc13a_toy.gtf` and the two synthetic BAMs into `data/`.
2. Run the `nexons.py` command from §3.2.
3. Load `nexons_output_unique.txt` / `_partial.txt` / `_gene.txt` with pandas
   and render a bar chart of canonical-vs-cryptic transcript counts per
   sample, next to the gene-level total, to make the "hidden at the gene
   level, visible at the transcript level" point visually explicit.
4. Print the QC outcome summary for each BAM as a sanity check.

#### 2.1.7. Caveats

- All sequence data here is **synthetic** and generated purely to exercise
  the software; it is not derived from any real patient, cell line, or public
  dataset, and no genomic coordinates are real genome coordinates.
- This demo shows the *computational quantitation step only*. Concluding
  anything about a real sample's TDP-43 status requires proper long-read
  library prep, alignment (e.g. minimap2 with spliced settings), and
  statistical comparison across biological replicates — none of which is in
  scope here.
- Nothing in this repository or document is intended for clinical or
  diagnostic use.

---

*Next step once this plan is approved: scaffold `sandbox/pixi.toml`, generate the files in `sandbox/data/`, and write the notebook described in §5.*