#!/usr/bin/env bash
#
# fastq_to_bam.sh -- align ONT cDNA FASTQ to a reference with minimap2 and
# emit a coordinate-sorted, indexed BAM ready for nexons.py.
#
# Usage:
#   ./fastq_to_bam.sh                      # all FASTQs in $FASTQ_DIR
#   ./fastq_to_bam.sh reads.fastq.gz       # one named file
#   BENCHMARK=200000 ./fastq_to_bam.sh     # time a subsample, then stop
#
# Configuration comes from the repo-root .env file. See README_fastqexamples.md.

set -euo pipefail

# --------------------------------------------------------------------------
# 1. Load the environment
# --------------------------------------------------------------------------
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SANDBOX="$(dirname "$HERE")"
REPO="$(dirname "$SANDBOX")"

if [[ -f "$REPO/.env" ]]; then
    # `set -a` exports everything the file defines, so the variables below
    # are visible to the tools we call. Values never get echoed.
    set -a; source "$REPO/.env"; set +a
    echo "[env]   loaded $REPO/.env"
else
    echo "[env]   WARNING: no $REPO/.env found; relying on the shell environment" >&2
fi

# Defaults for anything .env did not set.
FASTQ_DIR="${FASTQ_DIR:-$SANDBOX/big_data}"
BAM_DIR="${BAM_DIR:-$SANDBOX/big_data/bam}"
REF_FASTA="${REF_FASTA:-}"
REF_MMI="${REF_MMI:-}"
REF_GTF="${REF_GTF:-}"
THREADS="${THREADS:-$(nproc)}"
SORT_THREADS="${SORT_THREADS:-4}"
SORT_MEM="${SORT_MEM:-512M}"
# samtools sort spills to $TMPDIR. /tmp here is a 3.9 GiB tmpfs -- RAM, not
# disk -- so spilling there both runs out of room and steals memory from
# minimap2. Default onto the big ext4 volume instead. A spill that fills its
# filesystem kills the run late, after the alignment work is already done.
TMPDIR="${TMPDIR:-$SANDBOX/big_data/tmp}"
MM2_PRESET="${MM2_PRESET:--ax splice}"
BENCHMARK="${BENCHMARK:-}"
# minimap2 emits up to 5 secondary alignments per read by default. nexons has
# its own "secondary" bucket and skips those reads, so suppressing them costs
# no counts and shrinks the BAM -- but it also makes that bucket read zero, so
# the multi-mapping rate stops being visible in the QC. Set SECONDARY=yes to
# keep them.
SECONDARY="${SECONDARY:-no}"

mkdir -p "$BAM_DIR" "$TMPDIR"

# --------------------------------------------------------------------------
# 2. Preflight -- fail now, not forty minutes in
# --------------------------------------------------------------------------
for tool in minimap2 samtools; do
    command -v "$tool" >/dev/null 2>&1 || {
        echo "ERROR: $tool not on PATH." >&2
        echo "  Add it to sandbox/pixi.toml and run 'pixi install', then" >&2
        echo "  invoke this script as: pixi run ./fastq_to_bam.sh" >&2
        exit 1
    }
done

[[ -n "$REF_FASTA" || -n "$REF_MMI" ]] || {
    echo "ERROR: set REF_FASTA (reference genome) or REF_MMI (prebuilt index) in .env" >&2
    exit 1
}
# REF_FASTA is the genome you align TO. Putting reads there is an easy mistake
# to make and the resulting error ("does not exist") does not hint at the cause,
# so name it explicitly.
case "${REF_FASTA,,}" in
    *.fastq|*.fastq.gz|*.fq|*.fq.gz|*.fastq.gz.part)
        echo "ERROR: REF_FASTA looks like sequencing reads, not a reference genome:" >&2
        echo "         $REF_FASTA" >&2
        echo "  REF_FASTA is the genome you align TO -- a .fa/.fasta, e.g." >&2
        echo "         REF_FASTA=$SANDBOX/big_data/Homo_sapiens.GRCh38.dna.chromosome.19.fa" >&2
        echo "  Your reads are not set with a variable. Either pass them as arguments:" >&2
        echo "         pixi run ./fastq_to_bam/fastq_to_bam.sh reads.fastq.gz" >&2
        echo "  or point FASTQ_DIR at the directory holding them (default:" >&2
        echo "         $FASTQ_DIR)." >&2
        exit 1
        ;;
esac
if [[ -n "$REF_FASTA" && ! -f "$REF_FASTA" ]]; then
    echo "ERROR: REF_FASTA=$REF_FASTA does not exist" >&2
    echo "  Paths in .env are not resolved relative to .env -- a bare filename is" >&2
    echo "  looked up in your current directory. Use an absolute path." >&2
    exit 1
fi

# Collect the inputs.
if [[ $# -gt 0 ]]; then
    FASTQS=("$@")
else
    mapfile -t FASTQS < <(find "$FASTQ_DIR" -maxdepth 1 -type f \
        \( -name "*.fastq" -o -name "*.fq" -o -name "*.fastq.gz" -o -name "*.fq.gz" \) \
        | sort)
fi

# An in-progress download is the single most likely way to get a silently
# truncated BAM, because minimap2 will happily align a partial file and exit 0.
shopt -s nullglob
partials=("$FASTQ_DIR"/*.part "$FASTQ_DIR"/*.filepart "$FASTQ_DIR"/*.crdownload)
shopt -u nullglob
if [[ ${#partials[@]} -gt 0 ]]; then
    echo "ERROR: unfinished download(s) in $FASTQ_DIR:" >&2
    printf '  %s\n' "${partials[@]##*/}" >&2
    echo "  Wait for the transfer to finish and lose that suffix before aligning." >&2
    exit 1
fi

[[ ${#FASTQS[@]} -gt 0 ]] || { echo "ERROR: no FASTQ files found in $FASTQ_DIR" >&2; exit 1; }

echo "[in]    ${#FASTQS[@]} FASTQ file(s)"
echo "[out]   $BAM_DIR"
echo "[ref]   ${REF_MMI:-$REF_FASTA}"
echo "[cpu]   $THREADS align / $SORT_THREADS sort x $SORT_MEM"
echo "[tmp]   $TMPDIR ($(df -h "$TMPDIR" | awk 'NR==2{print $4}') free)"

# Disk preflight. `samtools sort` spills roughly the size of the final BAM,
# and the BAM is roughly 0.35-0.5x the uncompressed FASTQ, so budget half the
# input size on each of $TMPDIR and $BAM_DIR. Running out mid-sort discards
# the alignment work that preceded it, so refuse now instead.
in_kb=0
for fq in "${FASTQS[@]}"; do
    in_kb=$(( in_kb + $(du -k "$fq" 2>/dev/null | cut -f1) ))
done
need_kb=$(( in_kb / 2 ))
for pair in "TMPDIR:$TMPDIR" "BAM_DIR:$BAM_DIR"; do
    label=${pair%%:*}; dir=${pair#*:}
    free_kb=$(df -Pk "$dir" | awk 'NR==2{print $4}')
    if [[ "$free_kb" -lt "$need_kb" ]]; then
        echo "ERROR: $label=$dir has $((free_kb/1024)) MiB free;" >&2
        echo "       this run needs about $((need_kb/1024)) MiB there" >&2
        echo "       (~half of $((in_kb/1024)) MiB of input). Point $label at a" >&2
        echo "       bigger filesystem in .env, or free space." >&2
        exit 1
    fi
done

# A tmpfs is RAM, not disk. Spilling a multi-GB sort into it competes with
# minimap2's index for the same 8 GB and lands the machine in swap. This is
# reachable without editing anything: TMPDIR is often already exported by the
# shell, and "${TMPDIR:-default}" honours it.
tmp_fstype=$(findmnt -no FSTYPE --target "$TMPDIR" 2>/dev/null || echo unknown)
if [[ "$tmp_fstype" == "tmpfs" ]]; then
    echo "WARNING: $TMPDIR is on tmpfs (RAM-backed, $(df -h "$TMPDIR" | awk 'NR==2{print $2}'))." >&2
    echo "         Sort spills will consume RAM that minimap2 needs. Unset TMPDIR" >&2
    echo "         or set it to a path on disk." >&2
fi

# --------------------------------------------------------------------------
# 3. Reference index
# --------------------------------------------------------------------------
# Building the index once is worth it whenever more than one FASTQ is aligned:
# minimap2 otherwise re-indexes the reference on every invocation. The index is
# preset-specific -- an index built without -x splice will not give spliced
# alignments no matter what flags you pass at align time.
if [[ -z "$REF_MMI" ]]; then
    REF_MMI="${REF_FASTA%.gz}"; REF_MMI="${REF_MMI%.fa}"; REF_MMI="${REF_MMI%.fasta}"
    REF_MMI="${REF_MMI}.splice.mmi"
fi
if [[ ! -f "$REF_MMI" ]]; then
    echo "[index] building $REF_MMI (this is the memory-hungry step)"
    minimap2 $MM2_PRESET -t "$THREADS" -d "$REF_MMI" "$REF_FASTA" || {
        echo "ERROR: index build failed. If it was killed with no message the" >&2
        echo "  reference is too large for this machine's RAM -- see the" >&2
        echo "  'Memory' section of README_fastqexamples.md." >&2
        exit 1
    }
else
    echo "[index] reusing $REF_MMI"
fi

# --------------------------------------------------------------------------
# 4. Align
# --------------------------------------------------------------------------
# minimap2 writes SAM to stdout and we sort straight from the pipe. Writing an
# intermediate SAM would cost roughly twice the FASTQ size in disk for nothing.
align_one() {
    local fq="$1" sample bam t0 dt
    sample="$(basename "$fq")"
    sample="${sample%.gz}"; sample="${sample%.fastq}"; sample="${sample%.fq}"
    bam="$BAM_DIR/${sample}.bam"

    if [[ -f "$bam.bai" && "$bam" -nt "$fq" ]]; then
        echo "[skip]  $sample -- $bam.bai is newer than the input"
        return 0
    fi

    echo "[align] $sample -> $bam"
    t0=$SECONDS
    minimap2 $MM2_PRESET -t "$THREADS" \
             --secondary="$SECONDARY" \
             ${JUNC_BED:+--junc-bed "$JUNC_BED"} \
             "$REF_MMI" "$fq" \
      | samtools sort -@ "$SORT_THREADS" -m "$SORT_MEM" -T "$TMPDIR/sort.$sample" -o "$bam" -
    samtools index -@ "$SORT_THREADS" "$bam"
    dt=$(( SECONDS - t0 ))
    echo "[done]  $sample in $(( dt/60 ))m$(( dt%60 ))s"

    # A mapped-read count is the cheapest check that the BAM is neither empty
    # nor truncated, and it is the number nexons' totals should be read against.
    samtools flagstat -@ "$SORT_THREADS" "$bam" | head -7 | sed 's/^/        /'
}

# --------------------------------------------------------------------------
# 5. Benchmark mode -- measure, do not guess
# --------------------------------------------------------------------------
if [[ -n "$BENCHMARK" ]]; then
    src="${FASTQS[0]}"
    sub="$TMPDIR/benchmark_${BENCHMARK}reads.fastq"
    echo "[bench] taking $BENCHMARK reads from $(basename "$src")"
    { if [[ "$src" == *.gz ]]; then zcat "$src"; else cat "$src"; fi; } \
        | head -n $(( BENCHMARK * 4 )) > "$sub" || true

    bases=$(awk 'NR%4==2{n+=length($0)} END{print n+0}' "$sub")
    t0=$SECONDS
    minimap2 $MM2_PRESET -t "$THREADS" --secondary="$SECONDARY" "$REF_MMI" "$sub" > /dev/null
    dt=$(( SECONDS - t0 )); (( dt > 0 )) || dt=1

    echo
    awk -v b="$bases" -v d="$dt" -v th="$THREADS" 'BEGIN{
        r = b/d
        printf "[bench] %d bases in %ds\n", b, d
        printf "[bench] throughput: %.2f Mbase/s on %s threads (%.1f Gbase/hour)\n",
               r/1e6, th, r*3600/1e9
        printf "[bench] => %.1f min per GiB of uncompressed FASTQ (~469 Mbase/GiB here)\n",
               (469e6/r)/60
    }'
    echo "[bench] alignment only; add 10-20% for the sort and index."
    rm -f "$sub"
    exit 0
fi

for fq in "${FASTQS[@]}"; do align_one "$fq"; done

echo
echo "All BAMs in $BAM_DIR. Next:"
echo "  nexons.py ${REF_GTF:-<annotation.gtf>} $BAM_DIR/*.bam --direction none"
