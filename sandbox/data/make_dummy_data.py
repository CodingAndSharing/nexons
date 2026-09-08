#!/usr/bin/env python
"""
Generate the synthetic (dummy) long-read datasets used by the nexons ALS /
UNC13A cryptic-exon worked example.

Nothing here is real sequence data. A toy contig `UNC13A_TOY_chr` carries a
simplified 21-exon model of UNC13A plus one cryptic exon inserted into the
intron between exon 20 and exon 21 — the position at which TDP-43 loss of
function de-represses cryptic splicing in ALS/FTD. Reads are written directly
as BAM records with spliced CIGAR strings (N operations for the introns), so no
aligner is involved and the exon structure of every read is known by
construction.

Three example datasets are produced, each with the same biology but different
read populations, so that the nexons QC report looks different in each:

  example1  reference behaviour: most reads assign to a single transcript,
            all reads on the gene's own strand.
  example2  the GTF additionally contains an overlapping paralogous gene and
            most reads are short internal fragments compatible with both genes,
            so they land in "Multi Gene Match" (the purple ring of the
            Alignment Fate plot).
  example3  as example1, but one third of the reads are aligned antisense to
            the gene, so a third of the Alignment Directionality plot becomes
            "Opposite Strand Alignment" (blue).

Each example directory gets three BAMs representing the ALS contrast:
  control.bam          TDP-43 nuclear-functional  -> cryptic exon repressed
  tdp43_depleted.bam   TDP-43 loss of function    -> cryptic exon included
  mixed.bam            intermediate / mixed cell population

Usage
-----
    python make_dummy_data.py                 # build all three examples
    python make_dummy_data.py --example 2     # build one example
    python make_dummy_data.py --outdir data   # change destination root
"""

from __future__ import annotations

import argparse
import random
from pathlib import Path

import pysam

# --------------------------------------------------------------------------
# Toy locus geometry
# --------------------------------------------------------------------------

CHROM = "UNC13A_TOY_chr"
CHROM_LENGTH = 40_000

N_EXONS = 21            # simplified UNC13A: exons 1..21
EXON_LENGTH = 150
EXON_PITCH = 1_000      # exon start-to-start spacing (=> ~850 bp introns)
FIRST_EXON_START = 10_000

GENE_ID = "UNC13A_TOY"
GENE_NAME = "UNC13A_TOY"
TX_CANONICAL = "UNC13A_TOY-CANONICAL"
TX_CRYPTIC = "UNC13A_TOY-CRYPTIC"

PARALOG_GENE_ID = "UNC13A_TOY_PARALOG"
PARALOG_TX = "UNC13A_TOY_PARALOG-201"

# Reads placed here fall outside every annotated gene -> "No Gene" in the QC.
INTERGENIC_REGION = (32_000, 33_000)


def exon(i: int) -> tuple[int, int]:
    """1-based inclusive coordinates of annotated exon `i` (1..N_EXONS)."""
    start = FIRST_EXON_START + (i - 1) * EXON_PITCH
    return (start, start + EXON_LENGTH - 1)


#: The cryptic exon sits inside the intron between exon 20 and exon 21.
CRYPTIC_EXON = (exon(20)[1] + 351, exon(20)[1] + 351 + EXON_LENGTH - 1)

CANONICAL_EXONS = [exon(i) for i in range(1, N_EXONS + 1)]
CRYPTIC_EXONS = (
    [exon(i) for i in range(1, 21)] + [CRYPTIC_EXON] + [exon(21)]
)
#: The paralog shares exons 2..20 exactly but has neither exon 1 nor exon 21,
#: so full-length reads never match it while internal fragments do.
PARALOG_EXONS = [exon(i) for i in range(2, 21)]

GENE_START, GENE_END = CANONICAL_EXONS[0][0], CANONICAL_EXONS[-1][1]


# --------------------------------------------------------------------------
# GTF
# --------------------------------------------------------------------------

def _gtf_attrs(**kw: str) -> str:
    return " ".join(f'{k} "{v}";' for k, v in kw.items())


def _gtf_transcript_lines(gene_id, gene_name, tx_id, exons, strand):
    lines = []
    common = dict(gene_id=gene_id, gene_name=gene_name, gene_biotype="protein_coding")
    lines.append(
        [CHROM, "nexons_toy", "transcript", exons[0][0], exons[-1][1], ".", strand, ".",
         _gtf_attrs(**common, transcript_id=tx_id, transcript_biotype="protein_coding",
                    transcript_support_level="1")]
    )
    for n, (s, e) in enumerate(exons, start=1):
        lines.append(
            [CHROM, "nexons_toy", "exon", s, e, ".", strand, ".",
             _gtf_attrs(**common, transcript_id=tx_id, exon_number=str(n),
                        transcript_biotype="protein_coding", transcript_support_level="1")]
        )
    return lines


def write_gtf(path: Path, include_paralog: bool = False) -> Path:
    """Write the toy annotation. Ensembl-style attributes, TSL=1 throughout."""
    rows = []
    rows.append(
        [CHROM, "nexons_toy", "gene", GENE_START, GENE_END, ".", "+", ".",
         _gtf_attrs(gene_id=GENE_ID, gene_name=GENE_NAME, gene_biotype="protein_coding")]
    )
    rows += _gtf_transcript_lines(GENE_ID, GENE_NAME, TX_CANONICAL, CANONICAL_EXONS, "+")
    rows += _gtf_transcript_lines(GENE_ID, GENE_NAME, TX_CRYPTIC, CRYPTIC_EXONS, "+")

    if include_paralog:
        rows.append(
            [CHROM, "nexons_toy", "gene", PARALOG_EXONS[0][0], PARALOG_EXONS[-1][1],
             ".", "+", ".",
             _gtf_attrs(gene_id=PARALOG_GENE_ID, gene_name=PARALOG_GENE_ID,
                        gene_biotype="protein_coding")]
        )
        rows += _gtf_transcript_lines(PARALOG_GENE_ID, PARALOG_GENE_ID, PARALOG_TX,
                                      PARALOG_EXONS, "+")

    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "wt") as out:
        out.write("#!genome-build nexons_toy_UNC13A\n")
        out.write("#!synthetic annotation - not a real genome\n")
        for r in rows:
            out.write("\t".join(str(x) for x in r) + "\n")
    return path


# --------------------------------------------------------------------------
# Read construction
# --------------------------------------------------------------------------

# nexons defaults: --flex 10 (internal splice sites), --endflex 500 (transcript
# ends). Jitter stays inside those windows so the reads still match, while
# giving the QC flexibility plots a realistic spread around zero.
INNER_JITTER = [-3, -2, -1, 0, 0, 0, 0, 0, 1, 2, 3]
END_JITTER = (-150, 150)
MIN_BLOCK = 40          # never emit a CIGAR match block shorter than this


def _jittered_blocks(rng, exons, jitter_first_start=True, jitter_last_end=True):
    """Apply per-boundary jitter to a list of exon blocks."""
    blocks = [list(e) for e in exons]
    for i, b in enumerate(blocks):
        if i == 0:
            if jitter_first_start:
                b[0] += rng.randint(*END_JITTER)
        else:
            b[0] += rng.choice(INNER_JITTER)
        if i == len(blocks) - 1:
            if jitter_last_end:
                b[1] += rng.randint(*END_JITTER)
        else:
            b[1] += rng.choice(INNER_JITTER)
        # End jitter can be larger than the exon, which would give a
        # zero/negative-length CIGAR block. Keep every block >= MIN_BLOCK bp.
        if b[1] - b[0] + 1 < MIN_BLOCK:
            if i == 0 and len(blocks) > 1:
                b[0] = b[1] - MIN_BLOCK + 1
            else:
                b[1] = b[0] + MIN_BLOCK - 1
    return [(b[0], b[1]) for b in blocks]


def _make_read(rng, name, blocks, reverse=False, secondary=False, softclip=True):
    """Build one spliced BAM record from 1-based inclusive exon blocks."""
    cigar = []
    lead = rng.randint(5, 40) if (softclip and rng.random() < 0.35) else 0
    tail = rng.randint(5, 40) if (softclip and rng.random() < 0.35) else 0
    if lead:
        cigar.append((4, lead))
    for i, (s, e) in enumerate(blocks):
        if i:
            cigar.append((3, s - blocks[i - 1][1] - 1))   # N: intron
        cigar.append((0, e - s + 1))                      # M: exon
    if tail:
        cigar.append((4, tail))

    qlen = lead + tail + sum(e - s + 1 for s, e in blocks)

    a = pysam.AlignedSegment()
    a.query_name = name
    a.query_sequence = "".join(rng.choices("ACGT", k=qlen))
    a.query_qualities = pysam.qualitystring_to_array("I" * qlen)
    a.reference_id = 0
    a.reference_start = blocks[0][0] - 1                   # BAM is 0-based
    a.mapping_quality = 60
    a.cigartuples = cigar
    a.is_unmapped = False
    a.is_paired = False
    a.is_reverse = reverse
    a.is_secondary = secondary
    return a


def _unmapped_read(rng, name):
    a = pysam.AlignedSegment()
    a.query_name = name
    n = rng.randint(150, 1200)
    a.query_sequence = "".join(rng.choices("ACGT", k=n))
    a.query_qualities = pysam.qualitystring_to_array("I" * n)
    a.is_unmapped = True
    a.reference_id = -1
    a.reference_start = -1
    a.mapping_quality = 0
    a.is_paired = False
    return a


# --------------------------------------------------------------------------
# Read populations
# --------------------------------------------------------------------------

def _full_length(rng, cryptic):
    exons = CRYPTIC_EXONS if cryptic else CANONICAL_EXONS
    return _jittered_blocks(rng, exons)


def _partial_3prime(rng, cryptic):
    """Read starting inside an internal exon and running to the annotated 3'
    end. Because it contains (or skips) the cryptic exon it is diagnostic for
    exactly one transcript -> nexons status 'partial'."""
    k = rng.randint(2, 19)
    tail = [exon(i) for i in range(k, 21)]
    if cryptic:
        tail = tail + [CRYPTIC_EXON, exon(21)]
    else:
        tail = tail + [exon(21)]
    blocks = _jittered_blocks(rng, tail, jitter_first_start=False)
    # Start somewhere inside the first exon (a truncated read), which nexons
    # scores as an internal partial match.
    s, e = blocks[0]
    blocks[0] = (s + rng.randint(0, 110), e)
    return blocks


def _internal_fragment(rng):
    """Short internal fragment spanning only exons that the canonical and
    cryptic transcripts share -> ambiguous between transcripts (gene-level),
    and ambiguous between genes when the paralog is annotated."""
    a = rng.randint(3, 15)
    b = min(a + rng.randint(2, 5), 20)
    blocks = _jittered_blocks(rng, [exon(i) for i in range(a, b + 1)],
                              jitter_first_start=False, jitter_last_end=False)
    s, e = blocks[0]
    blocks[0] = (s + rng.randint(0, 110), e)
    s, e = blocks[-1]
    blocks[-1] = (s, e - rng.randint(0, 110))
    return blocks


def _intronic(rng):
    """Unspliced / intron-retaining fragment: matches the gene but no
    transcript -> gene-level count only."""
    i = rng.randint(5, 18)
    lo = exon(i)[1] + 60
    hi = exon(i + 1)[0] - 60
    s = rng.randint(lo, hi - 300)
    return [(s, s + rng.randint(200, 500))]


def _intergenic(rng):
    s = rng.randint(*INTERGENIC_REGION)
    return [(s, s + rng.randint(300, 900))]


#: Read-population recipes per example. Counts are numbers of primary
#: alignments; `cryptic` categories are split by the sample's cryptic fraction.
COMPOSITION = {
    1: dict(full=460, partial=190, internal=180, intronic=40, intergenic=30,
            secondary=60, unmapped=40, reverse_fraction=0.0, paralog=False),
    2: dict(full=350, partial=100, internal=450, intronic=0, intergenic=0,
            secondary=40, unmapped=30, reverse_fraction=0.0, paralog=True),
    3: dict(full=460, partial=190, internal=180, intronic=40, intergenic=30,
            secondary=60, unmapped=40, reverse_fraction=1 / 3, paralog=False),
}

#: Fraction of transcript-assignable reads carrying the cryptic exon.
SAMPLES = {
    "control": 0.05,          # TDP-43 in the nucleus, cryptic exon repressed
    "tdp43_depleted": 0.65,   # TDP-43 loss of function, cryptic exon included
    "mixed": 0.35,            # mixed population
}


def build_bam(path: Path, example: int, cryptic_fraction: float, seed: int) -> Path:
    rng = random.Random(seed)
    comp = COMPOSITION[example]
    rev_p = comp["reverse_fraction"]

    recipes = []
    for n, kind in ((comp["full"], "full"), (comp["partial"], "partial")):
        n_cryptic = round(n * cryptic_fraction)
        recipes += [(kind, True)] * n_cryptic + [(kind, False)] * (n - n_cryptic)
    recipes += [("internal", None)] * comp["internal"]
    recipes += [("intronic", None)] * comp["intronic"]
    recipes += [("intergenic", None)] * comp["intergenic"]
    recipes += [("secondary", None)] * comp["secondary"]
    rng.shuffle(recipes)

    builders = {
        "full": _full_length,
        "partial": _partial_3prime,
    }

    reads = []
    for i, (kind, cryptic) in enumerate(recipes):
        name = f"toyread_{i:05d}_{kind}"
        if kind in builders:
            blocks = builders[kind](rng, cryptic)
        elif kind == "internal":
            blocks = _internal_fragment(rng)
        elif kind == "intronic":
            blocks = _intronic(rng)
        elif kind == "intergenic":
            blocks = _intergenic(rng)
        else:  # secondary alignment of a full-length read
            blocks = _full_length(rng, rng.random() < cryptic_fraction)
        reads.append(
            _make_read(rng, name, blocks,
                       reverse=(rng.random() < rev_p),
                       secondary=(kind == "secondary"))
        )

    for i in range(comp["unmapped"]):
        reads.append(_unmapped_read(rng, f"toyread_unmapped_{i:04d}"))

    header = {"HD": {"VN": "1.6", "SO": "coordinate"},
              "SQ": [{"SN": CHROM, "LN": CHROM_LENGTH}]}

    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(".unsorted.bam")
    with pysam.AlignmentFile(str(tmp), "wb", header=header) as out:
        for r in reads:
            out.write(r)
    pysam.sort("-o", str(path), str(tmp))
    tmp.unlink()
    pysam.index(str(path))
    return path


def build_example(example: int, outroot: Path) -> Path:
    """Create one example directory (GTF + three BAMs + empty results dir)."""
    d = outroot / f"example{example}"
    write_gtf(d / "unc13a_toy.gtf", include_paralog=COMPOSITION[example]["paralog"])
    for i, (sample, cf) in enumerate(SAMPLES.items()):
        build_bam(d / f"{sample}.bam", example, cf, seed=1000 * example + i)
    (d / "results").mkdir(exist_ok=True)
    return d


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--outdir", default=str(Path(__file__).parent / "data/input"),
                   help="destination root (default: sandbox/data)")
    p.add_argument("--example", type=int, choices=sorted(COMPOSITION),
                   action="append", help="build only this example (repeatable)")
    args = p.parse_args(argv)

    wanted = args.example or sorted(COMPOSITION)
    for ex in wanted:
        d = build_example(ex, Path(args.outdir))
        print(f"example{ex}: {d}")


if __name__ == "__main__":
    main()
