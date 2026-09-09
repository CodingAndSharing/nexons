#!/usr/bin/env bash
#
# fastq_to_bam.sh -- align ONT cDNA reads to a genome and emit a sorted,
# indexed BAM for nexons.
#
# The whole pipeline is one command, run once per input FASTQ:
#
#   minimap2 -ax splice --secondary=no -t8 <index.mmi> <reads.fastq.gz> \
#     | samtools sort -o <out.bam> -
#
# Every parameter in it is explained in COMMAND_EXPLAINED.md. This script adds
# only what that command leaves out: a BAM index (nexons needs random access),
# an explicit sort temp directory, a loop over multiple barcodes, and four
# preflight checks that each correspond to a way this has actually failed.
#
# Usage
#   pixi run ./fastq_to_bam.sh                    # every FASTQ in $FASTQ_DIR
#   pixi run ./fastq_to_bam.sh reads1.fastq.gz    # just these
#
# Configuration comes from $REPO/.env or the environment. See .env.example.
#
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SANDBOX="$(dirname "$HERE")"
REPO="$(dirname "$SANDBOX")"

if [[ -f "$REPO/.env" ]]; then
    # `set -a` exports everything the file defines, so minimap2 and samtools
    # inherit it. $SANDBOX and $REPO are already set, so .env can use them.
    set -a; source "$REPO/.env"; set +a
    echo "[env]   loaded $REPO/.env"
else
    echo "[env]   no $REPO/.env; using the environment and defaults" >&2
fi

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
REF_MMI="${REF_MMI:-}"                       # prebuilt minimap2 index
REF_FASTA="${REF_FASTA:-}"                   # or a FASTA, to build one from
FASTQ_DIR="${FASTQ_DIR:-$SANDBOX/big_data/input}"
BAM_DIR="${BAM_DIR:-$SANDBOX/big_data/outputs}"

THREADS="${THREADS:-8}"                      # minimap2 -t
SORT_THREADS="${SORT_THREADS:-4}"            # samtools sort -@
SORT_MEM="${SORT_MEM:-512M}"                 # samtools sort -m, PER THREAD

# samtools sort spills to disk when a batch exceeds SORT_MEM. /tmp on this
# machine is a RAM-backed tmpfs, so spilling there consumes the memory
# minimap2 is holding the index in. Default onto the data volume instead.
TMPDIR="${TMPDIR:-$SANDBOX/big_data/tmp}"

MM2_PRESET="${MM2_PRESET:--ax splice}"       # see COMMAND_EXPLAINED.md
SECONDARY="${SECONDARY:-no}"                 # --secondary=

# ---------------------------------------------------------------------------
# Preflight -- fail now, not forty minutes in
# ---------------------------------------------------------------------------
for tool in minimap2 samtools; do
    command -v "$tool" >/dev/null 2>&1 || {
        echo "ERROR: $tool is not on PATH." >&2
        echo "  Run 'pixi install' in $SANDBOX, then invoke this as:" >&2
        echo "    pixi run ./fastq_to_bam/fastq_to_bam.sh" >&2
        exit 1
    }
done

# A reference is required, as either a prebuilt index or a FASTA to build from.
if [[ -z "$REF_MMI" && -z "$REF_FASTA" ]]; then
    echo "ERROR: set REF_MMI (a prebuilt .mmi) or REF_FASTA (a genome .fa)." >&2
    exit 1
fi

# REF_FASTA/REF_MMI is the genome you align TO. Reads go in $FASTQ_DIR or on
# the command line; they are never a reference. Catching this by name because
# the bare "does not exist" error gives no hint as to the cause.
for var in REF_MMI REF_FASTA; do
    case "${!var,,}" in
        *.fastq|*.fastq.gz|*.fq|*.fq.gz|*.part)
            echo "ERROR: $var looks like sequencing reads, not a reference:" >&2
            echo "         ${!var}" >&2
            echo "  Pass reads as arguments instead:" >&2
            echo "    pixi run ./fastq_to_bam/fastq_to_bam.sh ${!var}" >&2
            exit 1
            ;;
    esac
done

# Build the index if only a FASTA was given. minimap2 would otherwise re-index
# on every invocation, which for a genome costs more than the alignment.
if [[ -z "$REF_MMI" ]]; then
    [[ -f "$REF_FASTA" ]] || {
        echo "ERROR: REF_FASTA=$REF_FASTA does not exist." >&2
        echo "  Paths in .env resolve against your current directory, not" >&2
        echo "  against .env. Use an absolute path or \$SANDBOX/..." >&2
        exit 1
    }
    REF_MMI="${REF_FASTA%.gz}"; REF_MMI="${REF_MMI%.fa}"; REF_MMI="${REF_MMI%.fasta}"
    REF_MMI="${REF_MMI}.splice.mmi"
    if [[ ! -f "$REF_MMI" ]]; then
        echo "[index] building $REF_MMI -- the memory-hungry step"
        minimap2 $MM2_PRESET -t "$THREADS" -d "$REF_MMI" "$REF_FASTA"
    fi
fi
[[ -f "$REF_MMI" ]] || { echo "ERROR: REF_MMI=$REF_MMI does not exist." >&2; exit 1; }

# minimap2 holds the entire index resident for the whole run. If it does not
# fit in available RAM the kernel kills the process with no explanation, so
# compare the two up front. Whole-genome splice indexes run to ~10 GiB.
mmi_mib=$(( $(stat -c %s "$REF_MMI") / 1024 / 1024 ))
avail_mib=$(( $(awk '/MemAvailable/{print $2}' /proc/meminfo) / 1024 ))
# SORT_MEM carries a unit suffix, and it is per thread, not total.
case "${SORT_MEM: -1}" in
    G|g) sort_mib=$(( SORT_THREADS * ${SORT_MEM%[Gg]} * 1024 )) ;;
    M|m) sort_mib=$(( SORT_THREADS * ${SORT_MEM%[Mm]} )) ;;
    K|k) sort_mib=$(( SORT_THREADS * ${SORT_MEM%[Kk]} / 1024 )) ;;
    *)   sort_mib=$(( SORT_THREADS * SORT_MEM / 1024 / 1024 )) ;;
esac
gib() { awk -v m="$1" 'BEGIN{printf "%.1f GiB", m/1024}'; }
if (( mmi_mib + sort_mib + 1024 > avail_mib )); then
    echo "WARNING: this may not fit in memory." >&2
    echo "  index $(gib $mmi_mib) + sort $(gib $sort_mib) + ~1 GiB overhead" >&2
    echo "  vs $(gib $avail_mib) available." >&2
    echo "  If minimap2 is killed with no message, this was why. Options: align" >&2
    echo "  on a bigger machine, subset the reference to the chromosomes your" >&2
    echo "  GTF covers, or lower SORT_THREADS/SORT_MEM. See README_fastqexamples.md." >&2
fi

# Collect inputs: explicit arguments, else everything in $FASTQ_DIR.
if (( $# > 0 )); then
    FASTQS=("$@")
else
    mapfile -t FASTQS < <(find "$FASTQ_DIR" -maxdepth 1 -type f \
        \( -name "*.fastq" -o -name "*.fq" -o -name "*.fastq.gz" -o -name "*.fq.gz" \) \
        | sort)
fi
(( ${#FASTQS[@]} > 0 )) || { echo "ERROR: no FASTQ files in $FASTQ_DIR" >&2; exit 1; }

# An interrupted transfer is a valid-looking FASTQ that is simply missing its
# tail. minimap2 aligns it and exits 0, so the truncation shows up as missing
# counts in nexons rather than as an error. Refuse to start.
shopt -s nullglob
partials=("$FASTQ_DIR"/*.part "$FASTQ_DIR"/*.filepart "$FASTQ_DIR"/*.crdownload)
shopt -u nullglob
if (( ${#partials[@]} > 0 )); then
    echo "ERROR: unfinished download(s) in $FASTQ_DIR:" >&2
    printf '  %s\n' "${partials[@]##*/}" >&2
    echo "  Let the transfer finish, or move these aside before aligning." >&2
    exit 1
fi

mkdir -p "$BAM_DIR" "$TMPDIR"

# A relative TMPDIR -- easily inherited from the calling shell -- makes the
# sort scratch prefix depend on whatever directory the caller was in. Resolve
# it once so the path in the log is the path actually used.
TMPDIR="$(cd "$TMPDIR" && pwd)"

# TMPDIR on a tmpfs means sort spills eat RAM rather than disk.
if [[ "$(stat -f -c %T "$TMPDIR")" == "tmpfs" ]]; then
    echo "WARNING: TMPDIR=$TMPDIR is a tmpfs -- spills consume RAM." >&2
fi

echo "[ref]   $REF_MMI ($(gib $mmi_mib))"
echo "[cpu]   minimap2 -t $THREADS | samtools sort -@ $SORT_THREADS -m $SORT_MEM"
echo "[out]   $BAM_DIR"
echo "[in]    ${#FASTQS[@]} file(s)"

# ---------------------------------------------------------------------------
# Align
# ---------------------------------------------------------------------------
for fq in "${FASTQS[@]}"; do
    [[ -f "$fq" ]] || { echo "ERROR: $fq does not exist" >&2; exit 1; }

    # Strip one or two extensions to get a sample name: foo.fastq.gz -> foo
    sample="$(basename "$fq")"; sample="${sample%.gz}"
    sample="${sample%.fastq}"; sample="${sample%.fq}"
    bam="$BAM_DIR/$sample.bam"

    if [[ -f "$bam" && "$bam" -nt "$fq" ]]; then
        echo "[skip]  $sample -- $bam is newer than its input"
        continue
    fi

    echo "[align] $sample"
    t0=$SECONDS

    # This is the command from COMMAND_EXPLAINED.md. The pipe matters: the
    # uncompressed SAM for one barcode is tens of GiB, and piping it straight
    # into sort means it is never written to disk.
    minimap2 $MM2_PRESET --secondary="$SECONDARY" -t "$THREADS" \
             "$REF_MMI" "$fq" \
      | samtools sort -@ "$SORT_THREADS" -m "$SORT_MEM" \
                      -T "$TMPDIR/sort.$sample" -o "$bam" -

    # nexons reads the BAM by region, which needs the .bai index. The bare
    # piped command does not produce one.
    samtools index -@ "$SORT_THREADS" "$bam"

    dt=$(( SECONDS - t0 ))
    echo "[done]  $sample in $(( dt / 60 ))m$(( dt % 60 ))s -> $bam"

    # A mapped-read count is the cheapest evidence the BAM is neither empty
    # nor truncated, and it is what nexons' totals should be read against.
    samtools flagstat -@ "$SORT_THREADS" "$bam" | head -5 | sed 's/^/        /'
done

echo
echo "Next: quantify with nexons. The GTF must use the same assembly and the"
echo "same chromosome naming as the reference above."
echo "  python $REPO/nexons.py <annotation.gtf> $BAM_DIR/*.bam --direction none"
