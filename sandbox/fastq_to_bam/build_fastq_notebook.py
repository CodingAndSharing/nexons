#!/usr/bin/env python
"""
Regenerate sandbox/fastq_to_bam/fastq_to_bam_examples.ipynb.

The notebook is a build product, not something to hand-edit: keeping it in a
script means the prose and the code stay in one reviewable diff, and a broken
notebook can always be rebuilt with

    cd sandbox
    ./.pixi/envs/default/bin/python fastq_to_bam/build_fastq_notebook.py

Generated code cells use only `#` comments -- never a docstring -- because a
triple-quoted string inside a code cell terminates the raw string wrapping it.
"""
from pathlib import Path

import nbformat as nbf

HERE = Path(__file__).resolve().parent
OUT = HERE / "fastq_to_bam_examples.ipynb"

nb = nbf.v4.new_notebook()
C = nb.cells
md = lambda s: C.append(nbf.v4.new_markdown_cell(s.strip("\n")))
code = lambda s: C.append(nbf.v4.new_code_cell(s.strip("\n")))


# ==========================================================================
md(r"""
# FASTQ → BAM, run step by step

This notebook walks the pipeline in `README_fastqexamples.md` one stage at a
time, so you can see each check pass before committing an hour of CPU to the
full run. It is the same `fastq_to_bam.sh` underneath — the notebook adds the
diagnosis, the plots, and a measured time estimate.

**It is cheap by default.** Every heavy stage is either subsampled or gated
behind an explicit flag in the config cell below. Running the whole notebook
top to bottom on a 50,000-read subsample takes about a minute and produces a
real BAM with real QC plots. Nothing here starts a multi-hour job unless you
set `RUN_FULL = True`.

**Two things it deliberately does not do.** It never prints the contents of
`.env` — see the masked summary in step 1 — and it never needs the FTP
password, because the reads are already on disk.

> **Before you commit this notebook, clear its outputs.** Executed cells store
> their output inside the `.ipynb`, so a printed variable becomes a permanent
> part of the file. `Kernel → Restart Kernel and Clear All Outputs`, or
> `jupyter nbconvert --clear-output --inplace fastq_to_bam_examples.ipynb`.

### Kernel

Use **Python (nexons pixi)**. If it is not in the kernel list, the environment
has not been built yet:

```bash
cd sandbox && pixi install
```

`minimap2`, `samtools`, `seqkit` and `pigz` come from that same environment.
The notebook checks for them in step 2 and skips the stages that need them
rather than failing, so it is safe to run before the install finishes.
""")

# ==========================================================================
md(r"""
## Config

The only cell you should need to edit.
""")

code(r"""
from pathlib import Path
import gzip, os, random, re, shutil, subprocess, sys, time

# --- what to run ---------------------------------------------------------
SUBSAMPLE_READS = 50_000   # reads used for the worked example; None = all
RUN_FULL        = False    # True = align every FASTQ in full (tens of minutes+)
DO_DOWNLOAD_REF = True     # fetch GRCh38 chr19 (16.7 MB) if no reference is set
BENCHMARK_READS = 20_000   # reads used for the throughput measurement

# --- layout --------------------------------------------------------------
SANDBOX  = Path.cwd()
if SANDBOX.name == "fastq_to_bam":
    SANDBOX = SANDBOX.parent
REPO     = SANDBOX.parent
SCRIPT   = SANDBOX / "fastq_to_bam" / "fastq_to_bam.sh"
BIG      = SANDBOX / "big_data"
SUBDIR   = BIG / "subsample"
FIGDIR   = BIG / "qc_figures"

for d in (BIG, SUBDIR, FIGDIR):
    d.mkdir(parents=True, exist_ok=True)

assert SCRIPT.exists(), f"cannot find {SCRIPT} -- run this notebook from sandbox/ or sandbox/fastq_to_bam/"
print(f"repo     {REPO}")
print(f"sandbox  {SANDBOX}")
print(f"script   {SCRIPT.relative_to(REPO)}")
""")

# ==========================================================================
md(r"""
## A helper for long-running commands

`!cmd` buffers, so a twenty-minute alignment looks identical to a hung kernel.
This streams each line as it arrives and returns the exit status.
""")

code(r"""
def sh(cmd, extra_env=None, cwd=None, echo=True, capture=False):
    # Run a shell command, streaming output. Returns (exit_code, text).
    env = {**os.environ, **{k: str(v) for k, v in (extra_env or {}).items()}}
    if echo:
        print(f"$ {cmd}\n")
    p = subprocess.Popen(cmd, shell=True, cwd=str(cwd or SANDBOX), env=env,
                         stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                         text=True, bufsize=1)
    lines = []
    for line in p.stdout:
        if capture:
            lines.append(line)
        print(line, end="")
    return p.wait(), "".join(lines)


def have(tool):
    return shutil.which(tool) is not None
""")

# ==========================================================================
md(r"""
## Step 1 — environment variables, without exposing them

`fastq_to_bam.sh` sources `.env` from the repo root itself. The notebook reads
the same file only so it can show you what is set and apply it to this kernel's
environment.

Values for anything password-shaped are replaced with `••••` in the output
below. That matters more than it looks: a printed value would be written into
the notebook file and pushed with it.
""")

code(r"""
ENV_PATH = REPO / ".env"
SECRETISH = re.compile(r"(pass|pwd|secret|token|key|credential|auth)", re.I)

def load_env(path):
    # Minimal .env parser: KEY=value, ignoring blanks, comments and `export`.
    # Returns (values, note). `is_file()` rather than `exists()` because a
    # credentials file is often deliberately unreadable -- root-owned, mode
    # 600, or replaced by a device node by a sandbox -- and then it exists
    # but cannot be read. Reading it must not abort the notebook: the script
    # sources `.env` itself, so an unreadable copy here is inconvenient, not
    # fatal.
    out = {}
    if not path.is_file():
        return out, "absent"
    try:
        text = path.read_text()
    except OSError as e:
        return out, f"unreadable ({e.strerror})"
    for raw in text.splitlines():
        line = raw.strip()
        if not line or line.startswith("#") or "=" not in line:
            continue
        k, v = line.split("=", 1)
        k = k.strip()
        if k.startswith("export "):
            k = k[len("export "):].strip()
        out[k] = v.strip().strip('"').strip("'")
    return out, None

envvals, note = load_env(ENV_PATH)
if envvals:
    os.environ.update(envvals)
    print(f"{ENV_PATH.name}: {len(envvals)} variable(s) loaded into this kernel\n")
    for k in sorted(envvals):
        v = envvals[k]
        shown = "\u2022" * 8 + " (set)" if (SECRETISH.search(k) and v) else (v or "(empty)")
        print(f"  {k:<14} {shown}")
elif note == "absent":
    print(f"No {ENV_PATH} -- falling back to the defaults in fastq_to_bam.sh.")
    print("Copy the template if you want one:  cp .env.example .env")
else:
    print(f"{ENV_PATH} exists but is {note}.")
    print("Nothing here needs it: the reads are already on disk, and the")
    print("script sources .env in its own shell. Continuing with defaults.")

# Anything not supplied by .env falls back to these, matching the script.
os.environ.setdefault("THREADS", str(os.cpu_count() or 4))
print(f"\nthreads for alignment: {os.environ['THREADS']}")
""")

# ==========================================================================
md(r"""
## Step 2 — preflight

Everything that can be checked in a second, checked in a second. The
interesting rows are usually memory and the temp directory.
""")

code(r"""
def gib(nbytes):
    return nbytes / 1024**3

# --- tools ---
tools = {t: have(t) for t in ("minimap2", "samtools", "seqkit", "pigz")}
print("tools")
for t, ok in tools.items():
    print(f"  {'OK  ' if ok else 'MISSING'}  {t}")
if not all(tools[t] for t in ("minimap2", "samtools")):
    print("\n  -> run `pixi install` in sandbox/; alignment stages will be skipped until then")

# --- memory, and whether this is WSL ---
mem_kb = {}
for line in Path("/proc/meminfo").read_text().splitlines():
    parts = line.split()
    if parts and parts[0].rstrip(":") in ("MemTotal", "MemAvailable", "SwapTotal"):
        mem_kb[parts[0].rstrip(":")] = int(parts[1])
is_wsl = "microsoft" in Path("/proc/version").read_text().lower()

print(f"\nmemory")
print(f"  total      {mem_kb.get('MemTotal', 0)/1024**2:5.1f} GiB")
print(f"  available  {mem_kb.get('MemAvailable', 0)/1024**2:5.1f} GiB")
print(f"  swap       {mem_kb.get('SwapTotal', 0)/1024**2:5.1f} GiB")
if is_wsl:
    print("  WSL2 detected: this total is WSL's allocation, not the machine's RAM.")
    print("  If the host has more, raise `memory=` in C:\\Users\\<you>\\.wslconfig")
    print("  and `wsl --shutdown`. See Option D in the README.")

# --- disk, including the tmpfs trap ---
def fstype(path):
    r = subprocess.run(["findmnt", "-no", "FSTYPE", "--target", str(path)],
                       capture_output=True, text=True)
    return r.stdout.strip() or "unknown"

tmpdir = Path(os.environ.get("TMPDIR") or BIG / "tmp")
tmpdir.mkdir(parents=True, exist_ok=True)
print("\ndisk")
for label, p in (("TMPDIR", tmpdir), ("big_data", BIG)):
    u = shutil.disk_usage(p)
    ft = fstype(p)
    flag = "  <-- RAM-backed! sort spills will eat memory" if ft == "tmpfs" else ""
    print(f"  {label:<9} {gib(u.free):6.1f} GiB free   [{ft}]{flag}")

# --- unfinished downloads ---
partials = sorted(p.name for p in BIG.glob("*.part")) + \
           sorted(p.name for p in BIG.glob("*.filepart"))
print("\ndownloads")
if partials:
    print(f"  {len(partials)} unfinished transfer(s) present:")
    for n in partials:
        print(f"    {n}")
    print("  fastq_to_bam.sh refuses to run while these exist, which is why the")
    print("  worked example below uses its own subsample directory instead.")
else:
    print("  no .part files -- all transfers complete")
""")

# ==========================================================================
md(r"""
## Step 3 — what is actually in `big_data`

Read lengths first, because they set the time estimate and they tell you
whether the library is what you think it is.

Two quirks this handles. The FTP copies arrive **uncompressed despite a
`.fastq.gz` name**, so the reader sniffs the gzip magic bytes rather than
trusting the extension. And sampling only the head of a file biases the mean —
the first 200k reads of `barcode01` average 1,032 bp against 1,211 bp for the
whole file — so for a plain file this seeks to random offsets and resynchronises
to a record boundary instead.
""")

code(r"""
def is_gzip(path):
    with open(path, "rb") as fh:
        return fh.read(2) == b"\x1f\x8b"

def smart_open(path):
    return gzip.open(path, "rt") if is_gzip(path) else open(path, "rt")

def _read_at_boundary(fh):
    # Resynchronise to the next FASTQ record: an '@' line whose +2 line is '+'.
    # A quality line can also start with '@', hence the check.
    while True:
        l1 = fh.readline()
        if not l1:
            return None
        if not l1.startswith("@"):
            continue
        l2 = fh.readline()
        l3 = fh.readline()
        if not l3:
            return None
        if l3.startswith("+"):
            return len(l2.rstrip())

def sample_lengths(path, n=20_000, seed=0):
    # Random-offset sampling for plain files; sequential head for gzipped ones
    # (a gzip stream cannot be seeked into cheaply).
    size = path.stat().st_size
    if is_gzip(path):
        lens = []
        with smart_open(path) as fh:
            for i, line in enumerate(fh):
                if i % 4 == 1:
                    lens.append(len(line.rstrip()))
                if len(lens) >= n:
                    break
        return lens, "head (gzipped: cannot seek)"
    rng = random.Random(seed)
    lens = []
    with open(path, "rt", errors="replace") as fh:
        for _ in range(n):
            fh.seek(rng.randrange(0, max(size - 4096, 1)))
            fh.readline()                     # discard partial line
            L = _read_at_boundary(fh)
            if L:
                lens.append(L)
    return lens, "random offsets"

fastqs = sorted(p for p in BIG.iterdir()
                if p.is_file() and re.search(r"\.(fastq|fq)(\.gz)?(\.part)?$", p.name))
if not fastqs:
    print(f"No FASTQ files in {BIG}")
else:
    print(f"{'file':<52}{'GiB':>7}{'gzip?':>7}")
    for p in fastqs:
        print(f"{p.name[:50]:<52}{gib(p.stat().st_size):7.2f}{str(is_gzip(p)):>7}")

    target = max(fastqs, key=lambda p: p.stat().st_size)
    print(f"\nsampling read lengths from {target.name}")
    t0 = time.time()
    lens, how = sample_lengths(target, n=20_000)
    import statistics as st
    print(f"  method       {how}")
    print(f"  reads seen   {len(lens):,} in {time.time()-t0:.1f}s")
    print(f"  mean         {st.mean(lens):,.0f} bp")
    print(f"  median       {st.median(lens):,.0f} bp")
    print(f"  range        {min(lens):,} - {max(lens):,} bp")
    bytes_per_read = target.stat().st_size / max(len(lens), 1)  # placeholder, refined below
    est_reads = target.stat().st_size / (st.mean(lens) + 2 * 60 + 4)  # seq+qual+headers
    print(f"  est. reads   {est_reads:,.0f}  (~{est_reads*st.mean(lens)/1e9:.2f} Gbase)")
    if target.name.endswith(".part"):
        print("\n  NOTE: this file is still an unfinished transfer. Fine for the")
        print("  worked example below; never use it for real quantification.")
""")

# ==========================================================================
code(r"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import statistics as st

if fastqs and lens:
    fig, ax = plt.subplots(figsize=(7, 4))
    ax.hist(lens, bins=80, color="#4878a8", edgecolor="none")
    ax.axvline(st.mean(lens), color="#c04a3b", lw=1.5,
               label=f"mean {st.mean(lens):,.0f} bp")
    ax.axvline(st.median(lens), color="#e2a03f", lw=1.5, ls="--",
               label=f"median {st.median(lens):,.0f} bp")
    ax.set_yscale("log")
    ax.set_xlabel("read length (bp)")
    ax.set_ylabel("reads per bin (log scale)")
    ax.set_title(f"Read-length distribution — {target.name[:38]}\n"
                 f"{len(lens):,} reads sampled at {how}", fontsize=10)
    ax.legend(frameon=False)
    fig.tight_layout()
    out = FIGDIR / "read_length_distribution.png"
    fig.savefig(out, dpi=140)
    print(f"wrote {out}")
    plt.show()
""")

# ==========================================================================
md(r"""
## Step 4 — the reference

**This is the step that decides whether the run fits on this machine.** The
`splice` preset uses `-w5` against `map-ont`'s `-w10`, so the minimizer table
is roughly twice the size and a whole-GRCh38 splice index runs well past 10
GiB. On 8 GiB it is OOM-killed with no message beyond `Killed`.

The notebook defaults to **chr19** — 59 Mb, indexes in seconds, and where
`UNC13A` lives. State the trade-off in your methods: reads originating
elsewhere in the genome have nowhere else to go, so they either fail to map or
are forced onto the subset. If your GTF spans more chromosomes, add them here
or read Option D in the README about lifting the WSL memory ceiling.
""")

code(r"""
CHR19_URL = ("https://ftp.ensembl.org/pub/current_fasta/homo_sapiens/dna/"
             "Homo_sapiens.GRCh38.dna.chromosome.19.fa.gz")

ref = os.environ.get("REF_FASTA", "").strip()
if ref and Path(ref).exists():
    REF = Path(ref)
    print(f"using REF_FASTA from .env: {REF}")
elif DO_DOWNLOAD_REF:
    gz  = BIG / "Homo_sapiens.GRCh38.dna.chromosome.19.fa.gz"
    REF = gz.with_suffix("")            # .fa
    if REF.exists():
        print(f"already present: {REF.name}")
    else:
        sh(f"curl -L --fail -o {gz} {CHR19_URL}")
        # Ensembl ships plain gzip, not bgzip, so samtools faidx cannot read
        # it directly -- decompress rather than fight it.
        sh(f"gunzip -f {gz}")
    if have("samtools") and not (REF.parent / (REF.name + ".fai")).exists():
        sh(f"samtools faidx {REF}")
    print(f"\nreference: {REF}  ({gib(REF.stat().st_size):.2f} GiB)")
else:
    REF = None
    print("No reference. Set REF_FASTA in .env or set DO_DOWNLOAD_REF = True.")
""")

# ==========================================================================
md(r"""
## Step 5 — a subsample, so the example is cheap

Two reasons to subsample rather than point the script at `big_data` directly.
It keeps this walkthrough to about a minute; and `fastq_to_bam.sh` refuses to
start while any `.part` file sits in its input directory, which is the correct
behaviour but blocks you while later barcodes are still arriving. A separate
input directory sidesteps both.
""")

code(r"""
def subsample(src, dst, n_reads):
    # Copy the first n_reads records. Line-based, so it works on plain or
    # gzipped input and never splits a record.
    kept = 0
    with smart_open(src) as fin, open(dst, "w") as fout:
        for i, line in enumerate(fin):
            fout.write(line)
            if i % 4 == 3:
                kept += 1
                if kept >= n_reads:
                    break
    return kept

SUB = None
if fastqs:
    for stale in SUBDIR.glob("*.fastq"):
        stale.unlink()
    stem = re.sub(r"\.(fastq|fq)(\.gz)?(\.part)?$", "", target.name)
    SUB  = SUBDIR / f"{stem[:40]}_sub{SUBSAMPLE_READS//1000}k.fastq"
    t0 = time.time()
    kept = subsample(target, SUB, SUBSAMPLE_READS)
    print(f"{kept:,} reads -> {SUB.name}")
    print(f"  {gib(SUB.stat().st_size)*1024:.1f} MiB in {time.time()-t0:.1f}s")
""")

# ==========================================================================
md(r"""
## Step 6 — run the pipeline

Straight through `fastq_to_bam.sh`, so the notebook and the command line do
exactly the same thing. Watch the banner: the `[tmp]` line tells you where sort
spills are going, and the preflight refuses to start if there is not enough
room for them.
""")

code(r"""
BAMDIR = SUBDIR / "bam"
bam = None
if SUB and REF and all(tools[t] for t in ("minimap2", "samtools")):
    rc, _ = sh(f"bash {SCRIPT}", extra_env={
        "FASTQ_DIR": SUBDIR, "BAM_DIR": BAMDIR,
        "REF_FASTA": REF,    "REF_MMI": REF.with_suffix(".splice.mmi"),
        "TMPDIR":    BIG / "tmp",
    })
    print(f"\nexit {rc}")
    if rc == 0:
        bams = sorted(BAMDIR.glob("*.bam"))
        bam = bams[0] if bams else None
        print(f"BAM: {bam}")
else:
    missing = [t for t in ("minimap2", "samtools") if not tools[t]]
    print("Skipped." + (f" Not installed: {', '.join(missing)}." if missing else ""))
    print("Run `pixi install` in sandbox/, restart the kernel, and re-run.")
""")

# ==========================================================================
md(r"""
## Step 7 — measure the throughput, then extrapolate

Two minutes of measurement beats any estimate. `BENCHMARK=N` aligns N reads and
reports Mbase/s for *your* reference and thread count; the cell below turns
that into a per-barcode and whole-run figure.
""")

code(r"""
if SUB and REF and all(tools[t] for t in ("minimap2", "samtools")):
    rc, out = sh(f"bash {SCRIPT}", capture=True, extra_env={
        "FASTQ_DIR": SUBDIR, "BAM_DIR": BAMDIR, "REF_FASTA": REF,
        "REF_MMI": REF.with_suffix(".splice.mmi"), "TMPDIR": BIG / "tmp",
        "BENCHMARK": BENCHMARK_READS,
    })
    m = re.search(r"throughput:\s*([\d.]+)\s*Mbase/s", out)
    if m:
        mbs = float(m.group(1))
        gbase_per_barcode = 4.17          # measured for barcode01
        n_barcodes = max(len(fastqs), 1)
        per_bc = gbase_per_barcode * 1000 / mbs / 60      # minutes
        print(f"\n{'':<26}{'align':>9}{'+sort/index':>13}{'total':>9}")
        for label, n in ((f"1 barcode", 1), (f"{n_barcodes} present", n_barcodes),
                         ("12 barcodes", 12)):
            a = per_bc * n
            print(f"{label:<26}{a:8.0f}m{a*0.15:12.0f}m{a*1.15:8.0f}m")
        print(f"\nmeasured {mbs:.2f} Mbase/s against {REF.name}")
        print("Scaled from 4.17 Gbase per barcode (measured on barcode01).")
    else:
        print("\nCould not parse a throughput line from the benchmark output.")
else:
    print("Skipped -- needs minimap2, samtools and a reference.")
""")

# ==========================================================================
md(r"""
## Step 8 — is the BAM any good?

Three questions, in order of how often they catch something:

1. **Do reads have `N` in the CIGAR?** `N` is the spliced gap. If essentially
   none do, the index was built without `-x splice` and nexons will see an
   unspliced pile. This is the most common silent failure.
2. **What fraction mapped?** ONT cDNA against a matched whole genome should be
   well above 90%. Against a single chromosome, a low fraction is expected
   rather than a bug — most reads genuinely come from elsewhere.
3. **Is the file complete?** `samtools quickcheck` catches a BAM truncated by a
   sort that ran out of disk.
""")

code(r"""
qc = {}
if bam and bam.exists() and have("samtools"):
    sh(f"samtools quickcheck -v {bam} && echo 'quickcheck: OK'")
    rc, fs = sh(f"samtools flagstat {bam}", capture=True, echo=False)
    print(fs)
    import pysam
    n_spliced = n_mapped = n_total = 0
    aligned_lens = []
    with pysam.AlignmentFile(bam, "rb") as af:
        for r in af:
            n_total += 1
            if r.is_unmapped:
                continue
            n_mapped += 1
            aligned_lens.append(r.query_alignment_length)
            if r.cigarstring and "N" in r.cigarstring:
                n_spliced += 1
    qc = dict(total=n_total, mapped=n_mapped, spliced=n_spliced,
              aligned_lens=aligned_lens)
    print(f"records        {n_total:,}")
    print(f"mapped         {n_mapped:,} ({100*n_mapped/max(n_total,1):.1f}%)")
    print(f"spliced (N)    {n_spliced:,} ({100*n_spliced/max(n_mapped,1):.1f}% of mapped)")
    if n_mapped and n_spliced / n_mapped < 0.01:
        print("\n  WARNING: almost nothing is spliced. Delete the .mmi and rebuild")
        print("  with -x splice; nexons cannot find junctions that are not there.")
else:
    print("No BAM to check yet.")
""")

# ==========================================================================
code(r"""
if qc:
    fig, axes = plt.subplots(1, 2, figsize=(11, 4))

    a = axes[0]
    cats = ["mapped\nunspliced", "mapped\nspliced", "unmapped"]
    vals = [qc["mapped"] - qc["spliced"], qc["spliced"], qc["total"] - qc["mapped"]]
    bars = a.bar(cats, vals, color=["#4878a8", "#5a9e6f", "#b0b0b0"])
    for b, v in zip(bars, vals):
        a.text(b.get_x() + b.get_width()/2, v, f"{v:,}\n{100*v/max(qc['total'],1):.1f}%",
               ha="center", va="bottom", fontsize=9)
    a.set_ylabel("reads")
    a.set_title("Alignment outcome", fontsize=10)
    a.set_ylim(0, max(vals) * 1.25)
    a.spines[["top", "right"]].set_visible(False)

    b = axes[1]
    b.hist(qc["aligned_lens"], bins=70, color="#4878a8", edgecolor="none")
    b.set_yscale("log")
    b.set_xlabel("aligned length (bp)")
    b.set_ylabel("reads per bin (log scale)")
    b.set_title("Aligned length of mapped reads", fontsize=10)
    b.spines[["top", "right"]].set_visible(False)

    fig.suptitle(f"Alignment QC — {bam.name}", fontsize=11)
    fig.tight_layout()
    out = FIGDIR / "alignment_qc.png"
    fig.savefig(out, dpi=140)
    print(f"wrote {out}")
    plt.show()
""")

# ==========================================================================
md(r"""
## Step 9 — hand the BAM to nexons

`--direction none` because this is unstranded cDNA — the `dna_` basecall model
and the SMART-Seq kit both say so. `--annotated-bam` writes the per-read
classification tags that `sandbox/viz_scripts/nexons_viz.py` renders.

The GTF must use the same contig names as the reference FASTA. Ensembl headers
are `19`, UCSC are `chr19`; mixing them gives you a BAM that indexes cleanly
and a nexons run that finds nothing at all.
""")

code(r"""
gtf = os.environ.get("REF_GTF", "").strip()
cmd = (f"python {REPO/'nexons.py'} {gtf or '<annotation.gtf>'} "
       f"{BAMDIR}/*.bam --direction none "
       f"--outbase {BIG/'nexons_out'} --annotated-bam {BIG/'annotated'}")
print(cmd)

if bam and gtf and Path(gtf).exists():
    # Contig names must match, so check before spending the run.
    import pysam
    with pysam.AlignmentFile(bam, "rb") as af:
        bam_contigs = set(af.references)
    gtf_contigs = {l.split("\t")[0] for l in open(gtf)
                   if l.strip() and not l.startswith("#")}
    shared = bam_contigs & gtf_contigs
    print(f"\ncontigs: {len(bam_contigs)} in BAM, {len(gtf_contigs)} in GTF, "
          f"{len(shared)} shared")
    if not shared:
        print("  STOP: no shared contig names. Compare '19' vs 'chr19' before running.")
    else:
        sh(cmd)
else:
    print("\nSet REF_GTF in .env to run this automatically.")
""")

# ==========================================================================
md(r"""
## Step 10 — the full run

Everything above used a subsample. This runs every FASTQ in `big_data` at full
depth, which is tens of minutes per barcode; see the extrapolation in step 7
for your machine's number. It is gated on `RUN_FULL` in the config cell.

The script skips any sample whose `.bai` is newer than its FASTQ, so an
interrupted run resumes rather than restarting — and it will refuse to start
at all while a `.part` file is present.
""")

code(r"""
if RUN_FULL:
    rc, _ = sh(f"bash {SCRIPT}", extra_env={
        "FASTQ_DIR": BIG, "BAM_DIR": BIG / "bam",
        "REF_FASTA": REF, "REF_MMI": REF.with_suffix(".splice.mmi"),
        "TMPDIR": BIG / "tmp",
    })
    print(f"\nexit {rc}")
else:
    print("RUN_FULL is False. Set it True in the config cell to align everything.")
    print(f"Inputs that would be processed: {len(fastqs)}")
    for p in fastqs:
        flag = "  (unfinished transfer -- would be refused)" if p.name.endswith(".part") else ""
        print(f"  {p.name[:56]:<58}{gib(p.stat().st_size):6.2f} GiB{flag}")
""")

# ==========================================================================
md(r"""
## Troubleshooting

| Symptom | Cause |
|---|---|
| `Killed` during index build | OOM. Subset the reference, or raise the WSL memory limit (README Option D) |
| Script refuses to start, names a `.part` file | Transfer unfinished. Use a separate input directory, as steps 5–6 do |
| `WARNING: ... is on tmpfs` | `TMPDIR` inherited onto RAM-backed storage. Unset it or point it at disk |
| Aborts naming free space | Disk preflight: `TMPDIR` or `BAM_DIR` too small for ~half the input size |
| No `N` in any CIGAR | Index built without `-x splice`. Delete the `.mmi` and rebuild |
| nexons finds nothing | Contig naming mismatch (`19` vs `chr19`). Step 9 checks this |
| Kernel not in the list | `cd sandbox && pixi install`, then reopen |

**Before committing:** clear the outputs. `Kernel → Restart Kernel and Clear
All Outputs`, or

```bash
jupyter nbconvert --clear-output --inplace fastq_to_bam_examples.ipynb
```

`sandbox/fastq_to_bam/check_no_secrets.sh` scans notebook JSON along with
everything else, so a stored output containing a credential will fail the push
— but clearing them is cheaper than finding out that way.
""")

nb.metadata["kernelspec"] = {"display_name": "Python (nexons pixi)",
                             "language": "python", "name": "nexons-sandbox"}
nb.metadata["language_info"] = {"name": "python", "version": "3.12"}
nbf.write(nb, OUT)
print(f"wrote {OUT} ({len(nb.cells)} cells)")
