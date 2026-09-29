"""FASTQ quality control: FastQC-style reports, trimming, adapter removal,
nanopore stats, Dü-novo and BBDTools-style tools."""

from __future__ import annotations

import random
import re
from collections import Counter, defaultdict

from chroma_titan.core import align, bam, io, seq, stats
from chroma_titan.tools._common import *  # noqa: F401,F403
from chroma_titan.tools._common import ALN_SAM, GENES_FA, READS_1, READS_2, READS_SINGLE

QC = "fastq_quality_control"


def _reads(src) -> list:
    return io.parse_fastq(io.as_text(src))


def _qual(x) -> list[int]:
    return [ord(c) - 33 for c in x.qual]


def _write_reads(rs, name="output.fastq", message=""):
    return {"text": io.write_fastq(rs), "filename": name,
            "message": message or f"{len(rs)} reads"}


# ===========================================================================
# FastQC-style reports
# ===========================================================================
@T("fastqc_summary", "FastQC-style summary report", QC, "stats",
   [fq("src"), intin("top_sequences", 5), intin("contaminant_k", 5)],
   ex={"src": READS_SINGLE}, up="iuc/fastqc",
   summary="Per-base quality, GC, N content, duplication and adapter flags.")
def fastqc_summary(src, top_sequences=5, contaminant_k=5):
    """Compact multi-metric QC report for one FASTQ file."""
    rs = _reads(src)
    if not rs:
        return values_message({"reads": 0}, "no reads")
    lens = [len(r.seq) for r in rs]
    qs = [q for r in rs for q in _qual(r)]
    allseq = "".join(r.seq.upper() for r in rs)
    dups = Counter(r.seq.upper() for r in rs)
    st = {"reads": len(rs), "bases": sum(lens), "read_length_min": min(lens),
          "read_length_max": max(lens),
          "read_length_mode": Counter(lens).most_common(1)[0][0],
          "mean_quality": round(stats.mean(qs), 2) if qs else 0,
          "median_quality": stats.median(qs) if qs else 0,
          "mean_gc_percent": round(100 * seq.gc_content(allseq), 2),
          "percent_N": round(100 * allseq.count("N") / max(1, len(allseq)), 4),
          "percent_sequences_with_N": round(100 * sum(1 for r in rs if "N" in r.seq.upper()) / len(rs), 2),
          "percent_duplication": round(100 * sum(v - 1 for v in dups.values() if v > 1) / len(rs), 2),
          "most_duplicated_count": max(dups.values()) if dups else 0,
          "percent_over_Q30": round(100 * sum(1 for q in qs if q >= 30) / max(1, len(qs)), 2),
          "adapters_detected": []}
    kmers = seq.kmer_counts(allseq, k=int(contaminant_k), canonical=True)
    if kmers:
        km, cnt = kmers.most_common(1)[0]
        st["kmer_enrichment_top"] = f"{km}:{cnt}"
        st["kmer_over_represented"] = cnt / max(1, len(rs)) > 0.1
    for name, adapter in (("Illumina TruSeq", "AGATCGGAAGAGC"), ("Nextera", "CTGTCTCTTATACACATCT"),
                          ("small RNA 3'", "TGGAATTCTCGG")):
        n = sum(1 for r in rs if adapter[:12] in r.seq.upper())
        if n:
            st["adapters_detected"].append(f"{name}({n})")
    return values_message(st, f"FastQC-style report for {len(rs)} reads")


@T("fastq_per_base_quality", "Per base sequence quality", QC, "table", [fq("src")],
   ex={"src": READS_SINGLE}, up="iuc/fastqc",
   summary="Mean/median/quartiles of the quality score at every read position.")
def per_base_quality(src):
    """``fastx_quality_stats``-style per-cycle statistics."""
    rs = _reads(src)
    if not rs:
        return table([], "no reads")
    maxlen = max(len(r.qual) for r in rs)
    cols: list[list[int]] = [[] for _ in range(maxlen)]
    for r in rs:
        for i, c in enumerate(r.qual):
            cols[i].append(ord(c) - 33)
    out = []
    for i, vals in enumerate(cols):
        if not vals:
            continue
        out.append({"position": i + 1, "reads": len(vals), "mean": round(stats.mean(vals), 3),
                    "median": stats.median(vals), "q1": round(stats.quantile(vals, 0.25), 2),
                    "q3": round(stats.quantile(vals, 0.75), 2), "min": min(vals), "max": max(vals),
                    "percent_below_Q20": round(100 * sum(1 for v in vals if v < 20) / len(vals), 3)})
    res = table(out, f"quality by position over {maxlen} cycles")
    from chroma_titan.core import plot

    res["figure"] = plot.line({"mean Q": [r["mean"] for r in out], "Q1": [r["q1"] for r in out],
                               "Q3": [r["q3"] for r in out]}, x=[r["position"] for r in out],
                               title="per base quality", xlabel="cycle", ylabel="Phred")
    return res


@T("fastq_quality_hist", "Quality score distribution", QC, "table", [fq("src")],
   ex={"src": READS_SINGLE}, up="fastx_quality_statistics",
   summary="Counts of each Phred score in the whole file.")
def quality_histogram(src):
    """Histogram over all quality bases (``fastx_quality_statistics``)."""
    rs = _reads(src)
    c = Counter(q for r in rs for q in _qual(r))
    if not c:
        return table([], "no qualities")
    tot = sum(c.values())
    out = [{"quality": q, "count": c[q], "ascii": q + 33,
            "percent": round(100 * c[q] / tot, 4)} for q in sorted(c)]
    res = table(out, f"{tot} bases in {len(c)} bins")
    from chroma_titan.core import plot

    res["figure"] = plot.bar([str(r["quality"]) for r in out], [r["count"] for r in out],
                             title="quality score distribution", xlabel="Phred", ylabel="count")
    return res


@T("fastq_n_stats", "N content per read position", QC, "table", [fq("src"), boolean("summary", False)],
   ex={"src": READS_SINGLE}, up="fastx_n_statistics",
   summary="Number and percentage of Ns at each cycle (or a single summary).")
def n_stats(src, summary=False):
    """``fastx_n_statistics`` behaviour."""
    rs = _reads(src)
    if not rs:
        return table([], "no reads")
    if summary:
        total = sum(len(r.seq) for r in rs)
        n = sum(r.seq.upper().count("N") for r in rs)
        return values_message({"reads": len(rs), "bases": total, "n_bases": n,
                              "percent_n": round(100 * n / max(1, total), 5),
                              "reads_with_n": sum(1 for r in rs if "N" in r.seq.upper())},
                             f"{n} N bases")
    maxlen = max(len(r.seq) for r in rs)
    cols = [0] * maxlen
    for r in rs:
        for i, c in enumerate(r.seq.upper()):
            if c == "N":
                cols[i] += 1
    return table([{"position": i + 1, "n_count": v, "percent": round(100 * v / len(rs), 4)}
                  for i, v in enumerate(cols)], "N content per cycle")


@T("fastq_gc_distribution", "Per read GC content distribution", QC, "table",
   [fq("src"), intin("bins", 20)], ex={"src": READS_SINGLE, "bins": 20}, up="fastqgc",
   summary="GC% histogram of reads; a bimodal plot means contamination.")
def fastq_gc_distribution(src, bins=20):
    """Galaxy 'NGS: QC and graphing' GC distribution."""
    rs = _reads(src)
    gs = [100 * seq.gc_content(r.seq) for r in rs if r.seq]
    if not gs:
        return table([], "no reads")
    lo, hi = min(gs), max(gs)
    w = ((hi - lo) / max(1, int(bins))) if hi > lo else 1
    nb = max(1, int(bins))
    c = Counter(min(nb - 1, int((g - lo) / w)) for g in gs)
    out = [{"bin": f"{lo + i * w:.1f}-{lo + (i + 1) * w:.1f}", "reads": c.get(i, 0)} for i in range(nb)]
    res = table(out, f"GC distribution of {len(gs)} reads (mean {stats.mean(gs):.2f}%)")
    from chroma_titan.core import plot

    res["figure"] = plot.bar([r["bin"] for r in out], [r["reads"] for r in out],
                             title="per sequence GC", xlabel="GC %", ylabel="reads")
    return res


@T("fastq_per_sequence_quality", "Per sequence quality scores", QC, "table",
   [fq("src"), intin("bins", 15)], ex={"src": READS_SINGLE}, up="iuc/fastqc",
   summary="Histogram of the mean quality of each read.")
def per_sequence_quality(src, bins=15):
    """Distribution of per-read average quality."""
    rs = _reads(src)
    vals = [stats.mean(_qual(r)) for r in rs if r.qual]
    if not vals:
        return table([], "no qualities")
    lo, hi = min(vals), max(vals)
    nb = max(1, int(bins))
    w = ((hi - lo) / nb) if hi > lo else 1
    c = Counter(min(nb - 1, int((v - lo) / w)) for v in vals)
    return table([{"quality_bin": f"{lo + i * w:.1f}-{lo + (i + 1) * w:.1f}", "reads": c.get(i, 0)}
                  for i in range(nb)], f"mean read quality {stats.mean(vals):.2f}")


@T("fastq_duplication_levels", "Sequence duplication levels", QC, "table",
   [fq("src"), intin("max_level", 10)], ex={"src": READS_SINGLE}, up="iuc/fastqc",
   summary="Percentage of reads seen 1x, 2x, … >N times (PCR duplicates).")
def duplication_levels(src, max_level=10):
    """FastQC 'duplication levels' plot data."""
    rs = _reads(src)
    c = Counter(r.seq.upper() for r in rs)
    freq = Counter(c.values())
    ml = int(max_level)
    out = []
    for level in range(1, ml):
        out.append({"duplication_level": str(level), "reads": freq.get(level, 0) * level,
                    "percent_of_total": round(100 * freq.get(level, 0) * level / max(1, len(rs)), 3)})
    over = sum(v * level for level, v in freq.items() if level >= ml)
    out.append({"duplication_level": f">{ml - 1}", "reads": over,
                "percent_of_total": round(100 * over / max(1, len(rs)), 3)})
    dup_copies = sum(v - 1 for v in c.values() if v > 1)
    out.append({"duplication_level": "unique fraction", "reads": len(rs) - dup_copies,
                "percent_of_total": round(100 * (len(rs) - dup_copies) / max(1, len(rs)), 3)})
    return table(out, f"{len(rs)} reads, {dup_copies} duplicate copies")


@T("fastq_overrepresented", "Overrepresented sequences", QC, "table",
   [fq("src"), number("threshold_percent", 0.1), intin("top", 20)], ex={"src": READS_SINGLE},
   up="iuc/fastqc", summary="Sequences appearing more often than a threshold.")
def overrepresented(src, threshold_percent=0.1, top=20):
    """FastQC 'overrepresented sequences' table."""
    rs = _reads(src)
    c = Counter(r.seq.upper() for r in rs)
    n = max(1, len(rs))
    out = [{"sequence": s[:70], "count": v, "percent": round(100 * v / n, 4),
            "example_id": next((r.id for r in rs if r.seq.upper() == s), "")}
           for s, v in c.most_common(int(top)) if 100 * v / n >= threshold_percent]
    return table(out, f"{len(out)} overrepresented sequences")


@T("fastq_kmer_content", "Overrepresented k-mers", QC, "table",
   [fq("src"), intin("k", 5), intin("top", 10)], ex={"src": READS_SINGLE, "k": 5},
   up="iuc/fastqc", summary="K-mer frequencies at the start of reads (adapter signature).")
def kmer_content(src, k=5, top=10):
    """FastQC 'kmer content' style table."""
    rs = _reads(src)
    kk = int(k)
    c = Counter(r.seq.upper()[:kk] for r in rs if len(r.seq) >= kk)
    n = max(1, sum(c.values()))
    out = [{"kmer": km, "count": v, "percent": round(100 * v / n, 4)} for km, v in c.most_common(int(top))]
    return table(out, f"top {kk}-mers at read start")


@T("fastq_encoding", "Detect the quality encoding", QC, "stats", [fq("src")], ex={"src": READS_SINGLE},
   up="fastq-encoding-detect", summary="Guess Phred+33/+64/Solexa+66 from ASCII ranges.")
def detect_encoding(src):
    """Report the min/max ASCII codes and the inferred encoding."""
    body = io.as_text(src).splitlines()
    qchars = [c for i in range(1, len(body), 4) for c in body[i]]
    if not qchars:
        return values_message({"encoding": "none"}, "no quality lines")
    lo, hi = min(ord(c) for c in qchars), max(ord(c) for c in qchars)
    if lo >= 66 and hi <= 104:
        enc = "Illumina 1.3+ / Phred+64"
    elif lo >= 59:
        enc = "Solexa / older Illumina (Phred+66)"
    else:
        enc = "Phred+33 (Illumina 1.8+, Sanger, nanopore)"
    return values_message({"min_ascii": lo, "max_ascii": hi, "phred_min": lo - 33,
                          "phred_max": hi - 33, "encoding": enc, "quality_lines": len(body) // 4},
                         f"quality encoding: {enc}")


@T("fastq_tile_stats", "Per-lane / per-tile quality summary", QC, "table",
   [fq("src")], ex={"src": READS_SINGLE}, up="fastq_stats",
   summary="Aggregate quality by the flowcell/lane field of the read name when present.")
def tile_stats(src):
    """Group reads by a tile/lane id parsed from the read name."""
    rs = _reads(src)
    groups: dict[str, list[float]] = defaultdict(list)
    for r in rs:
        parts = r.id.replace(":", " ").split()
        tile = parts[1] if len(parts) > 1 else "default"
        groups[tile].append(stats.mean(_qual(r)) if r.qual else 0.0)
    out = [{"tile": k, "reads": len(v), "mean_quality": round(stats.mean(v), 3),
            "min_quality": round(min(v), 2), "max_quality": round(max(v), 2)}
           for k, v in sorted(groups.items(), key=lambda kv: stats.mean(kv[1]))[:200]]
    return table(out, f"{len(out)} tiles/lanes")


@T("fastq_read_length_distribution", "Read length distribution", QC, "table",
   [fq("src"), intin("bins", 0)], ex={"src": READS_SINGLE}, up="fastq_len",
   summary="Histogram of read lengths (with optional binning).")
def read_length_distribution(src, bins=0):
    """FastQC 'sequence length distribution'."""
    rs = _reads(src)
    lens = [len(r.seq) for r in rs]
    if not lens:
        return table([], "no reads")
    c = Counter(lens)
    out = [{"length": k, "reads": v, "percent": round(100 * v / len(rs), 4)} for k, v in sorted(c.items())]
    if bins and len(out) > int(bins):
        w = max(1, len(out) // int(bins))
        merged = []
        for i in range(0, len(out), w):
            grp = out[i:i + w]
            merged.append({"length": (f"{grp[0]['length']}-{grp[-1]['length']}"
                                      if len(grp) > 1 else str(grp[0]["length"])),
                           "reads": sum(g["reads"] for g in grp),
                           "percent": round(sum(g["percent"] for g in grp), 4)})
        out = merged
    from chroma_titan.core import plot

    res = table(out, f"lengths {min(lens)}-{max(lens)} bp")
    res["figure"] = plot.bar([str(r["length"]) for r in out[:60]], [r["reads"] for r in out[:60]],
                             title="read length distribution", xlabel="bp", ylabel="reads")
    return res


# ===========================================================================
# Trimming / filtering
# ===========================================================================
@T("fastq_mott_trimmer", "Trim sequences by quality (Mott's algorithm)", QC, "text",
   [fq("src"), number("threshold", 5.0), number("window", 0.5), intin("min_length", 20),
    choice("which_end", ["five_prime", "three_prime", "both"])],
   ex={"src": READS_SINGLE, "threshold": 1.0, "window": 0.3}, up="dynamictrim",
   summary="Sliding-window quality trimming as used by prinseq/dynamictrim.")
def mott_trimmer(src, threshold=5.0, window=0.5, min_length=20, which_end="both"):
    """Mott's algorithm: drop the end with the worst running quality score."""
    rs = _reads(src)
    out = []
    thr = float(threshold)

    def trim_once(quals, keep_from_end: bool):
        """Return the number of bases to remove from one end."""
        if not quals:
            return 0
        limit = stats.mean(quals) - thr
        score, best, best_i = 0.0, 0.0, 0
        for i, v in enumerate(quals):
            score += limit - v
            if score > best:
                best, best_i = score, i
        return 0 if best <= 0 else best_i + 1

    for r in rs:
        s, q = r.seq.upper(), _qual(r)
        if which_end in ("five_prime", "both"):
            cut = trim_once(q, False)
            s, q = s[cut:], q[cut:]
        if which_end in ("three_prime", "both") and q:
            cut = trim_once(q[::-1], True)
            if cut:
                s, q = s[:-cut], q[:-cut]
        if len(s) >= int(min_length) and s:
            out.append(io.Read(r.id, s, "".join(chr(min(93, v + 33)) for v in q), r.desc))
    return _write_reads(out, "trimmed.fastq", f"{len(out)}/{len(rs)} reads after Mott trimming")


@T("fastq_sliding_trim", "Trim Galore! style quality/adapter trim", QC, "text",
   [fq("src"), intin("quality_cutoff", 20), intin("length_min", 25), boolean("trim_ns", True),
    boolean("retain_untrimmed", False)], ex={"src": READS_SINGLE, "quality_cutoff": 20},
   up="iuc/trim_galore", summary="Trim bases below a Phred cut-off from both ends, drop short reads.")
def sliding_trim(src, quality_cutoff=20, length_min=25, trim_ns=True, retain_untrimmed=False):
    """``trim_galore --quality`` behaviour."""
    rs = _reads(src)
    out, trimmed_bases, reads_trimmed = [], 0, 0
    for r in rs:
        s, q = r.seq, _qual(r)
        before = len(s)
        while q and q[0] < quality_cutoff:
            q, s = q[1:], s[1:]
        while q and q[-1] < quality_cutoff:
            q, s = q[:-1], s[:-1]
        if trim_ns:
            while s and s[0].upper() == "N":
                s, q = s[1:], q[1:]
            while s and s[-1].upper() == "N":
                s, q = s[:-1], q[:-1]
        if len(s) != before:
            reads_trimmed += 1
            trimmed_bases += before - len(s)
        if len(s) >= length_min and s:
            out.append(io.Read(r.id, s, "".join(chr(v + 33) for v in q), r.desc))
        elif retain_untrimmed and len(s) != before:
            out.append(io.Read(r.id, s, "".join(chr(v + 33) for v in q), r.desc))
    return _write_reads(out, "trimmed.fastq",
                       f"{len(out)} reads kept, {trimmed_bases} bp trimmed from {reads_trimmed} reads")


@T("fastq_trim_ns", "Trim Ns from read ends", QC, "text", [fq("src"), intin("min_length", 0)],
   ex={"src": READS_SINGLE}, up="du_novo/prinseq", summary="Remove N bases from 5' and 3' ends.")
def trim_ns_tool(src, min_length=0):
    """Prinseq ``--trim_ns``."""
    rs = _reads(src)
    out = []
    for r in rs:
        s, q = r.seq, r.qual
        while s and s[0].upper() == "N":
            s, q = s[1:], q[1:]
        while s and s[-1].upper() == "N":
            s, q = s[:-1], q[:-1]
        if s and len(s) >= int(min_length):
            out.append(io.Read(r.id, s, q, r.desc))
    return _write_reads(out, "trimmed.fastq", f"{len(out)}/{len(rs)} reads after N trimming")


@T("fastq_trim_by_length", "Filter reads by length", QC, "text",
   [fq("src"), intin("min_len", 0), intin("max_len", 0), choice("mode", ["keep", "discard"])],
   ex={"src": READS_SINGLE, "min_len": 60, "max_len": 160}, up="fastqlenfilter",
   summary="Keep or drop reads outside a length window.")
def filter_by_length(src, min_len=0, max_len=0, mode="keep"):
    """Size selection of reads."""
    rs = _reads(src)
    keep = []
    for r in rs:
        ok = (not min_len or len(r.seq) >= int(min_len)) and (not max_len or len(r.seq) <= int(max_len))
        if ok == (mode == "keep"):
            keep.append(r)
    return _write_reads(keep, "filtered.fastq", f"{len(keep)}/{len(rs)} reads ({mode})")


@T("fastq_cutadapt", "Cutadapt-style adapter trimming", QC, "text",
   [fq("src"), textbox("adapter", "AGATCGGAAGAGC"), intin("min_overlap", 3), number("errors", 0.1),
    boolean("trim_nextera", False), intin("minimum_length", 15)],
   ex={"src": READS_SINGLE, "adapter": "AGATCGGAAGAGC", "min_overlap": 8}, up="iuc/cutadapt",
   summary="Find the best 3' adapter match and clip everything from there.")
def cutadapt_lite(src, adapter="AGATCGGAAGAGC", min_overlap=3, errors=0.1, trim_nextera=False,
                  minimum_length=15):
    """Prefix-match adapter search with mismatches (``cutadapt -a``)."""
    rs = _reads(src)
    adapters = [adapter.upper()] if adapter else []
    if trim_nextera:
        adapters.append("CTGTCTCTTATACACATCT")
    out, n_trimmed, total_bp = [], 0, 0
    for r in rs:
        s = r.seq.upper()
        cut, best_err = None, 1e18
        for a in adapters:
            L = len(a)
            for ov in range(L, max(int(min_overlap), 1) - 1, -1):
                if ov > len(s):
                    continue
                mism = sum(1 for x, y in zip(s[-ov:], a[L - ov:]) if x != y)
                if mism <= errors * ov and mism < best_err:
                    best_err, cut = mism, len(s) - ov
        if cut is None:
            if len(s) >= minimum_length:
                out.append(r)
            continue
        n_trimmed += 1
        total_bp += len(s) - cut
        if cut >= minimum_length:
            out.append(io.Read(r.id, s[:cut], r.qual[:cut], r.desc))
    return _write_reads(out, "trimmed.fastq",
                       f"{n_trimmed} adapter-trimmed reads, {total_bp} bp removed")


@T("fastq_polyx_trim", "Trim poly-A / poly-T tails", QC, "text",
   [fq("src"), choice("tail", ["polyA", "polyT", "both"]), intin("min_run", 5)],
   ex={"src": READS_SINGLE, "tail": "polyA"}, up="iuc/trim_galore",
   summary="Remove homopolymer A/T tails (RNA-seq artefacts).")
def polyx_trim(src, tail="both", min_run=5):
    """Trim homopolymer tails of A (or T) from the 3' end."""
    rs = _reads(src)
    chars = {"polyA": "A", "polyT": "T", "both": "AT"}[tail]
    out, n = [], 0
    for r in rs:
        s = r.seq.upper()
        m = re.search("[%s]{%d,}$" % (chars, max(1, int(min_run))), s)
        if m:
            n += 1
            out.append(io.Read(r.id, s[:m.start()], r.qual[:m.start()], r.desc))
        else:
            out.append(r)
    return _write_reads(out, "trimmed.fastq", f"{n} reads with {chars} tails trimmed")


@T("fastq_rrbs_trim", "Trim for RRBS (MspI site)", QC, "text",
   [fq("src"), textbox("site", "CTAG"), intin("trim_site", 0)], ex={"src": READS_SINGLE, "site": "CTAG"},
   up="iuc/trim_galore", summary="Cut at the MspI/CTAG site (``trim_galore --rrbs``).")
def rrbs_trim(src, site="CTAG", trim_site=0):
    """Remove the RRBS digestion site after adapters."""
    rs = _reads(src)
    out, cut = [], 0
    for r in rs:
        s = r.seq.upper()
        i = s.find(site.upper())
        if i >= 0:
            keep = i + (int(trim_site) if trim_site else len(site))
            out.append(io.Read(r.id, s[:keep], r.qual[:keep], r.desc))
            cut += 1
        else:
            out.append(r)
    return _write_reads(out, "rrbs.fastq", f"{cut}/{len(rs)} reads cut at {site}")


@T("fastq_sortbylength", "Sort reads by length or name", QC, "text",
   [fq("src"), choice("by", ["length", "name", "quality", "sequence"]), boolean("largest_first", True)],
   ex={"src": READS_SINGLE, "by": "length"}, up="iuc/fastx_toolkit",
   summary="Deterministic read ordering (SortByLength/SortBySize).")
def sort_by_length(src, by="length", largest_first=True):
    """``sort_by_length`` for FASTQ."""
    rs = _reads(src)
    key = {"length": lambda r: len(r.seq), "name": lambda r: r.id,
           "quality": lambda r: stats.mean(_qual(r)) if r.qual else 0.0,
           "sequence": lambda r: r.seq.upper()}[by]
    return _write_reads(sorted(rs, key=key, reverse=largest_first), "sorted.fastq",
                       f"sorted {len(rs)} reads by {by}")


@T("fastq_collapser", "Collapse duplicate sequences", QC, "table",
   [fq("src"), intin("min_count", 1)], ex={"src": READS_SINGLE}, up="iuc/fastx_toolkit",
   summary="Unique sequences with abundance and mean quality (``fastx_collapser``).")
def fastq_collapser(src, min_count=1):
    """Dereplication with counts."""
    rs = _reads(src)
    groups: dict[str, list] = defaultdict(list)
    for r in rs:
        groups[r.seq.upper()].append(r)
    out = []
    for s, v in groups.items():
        if len(v) < int(min_count):
            continue
        out.append({"sequence": s[:60], "count": len(v),
                    "abundance_percent": round(100 * len(v) / max(1, len(rs)), 4),
                    "mean_quality": round(stats.mean([stats.mean(_qual(x)) for x in v]), 3),
                    "first_id": v[0].id, "length": len(s)})
    out.sort(key=lambda d: -d["count"])
    return table(out, f"{len(out)} unique sequences from {len(rs)} reads")


@T("fastq_randomize", "Randomly shuffle reads", QC, "text", [fq("src"), intin("seed", 42)],
   ex={"src": READS_SINGLE}, up="iuc/bbtools", summary="Shuffle read order to remove positional bias.")
def randomize_reads(src, seed=42):
    """``randomize`` for FASTQ/FASTA."""
    rs = _reads(src)
    random.Random(seed).shuffle(rs)
    return _write_reads(rs, "shuffled.fastq", f"shuffled {len(rs)} reads")


@T("fastq_subsample", "Subsample reads (bbduk style)", QC, "text",
   [fq("src"), intin("n", 0), number("fraction", 0.1), intin("seed", 1), fq("src2", "")],
   ex={"src": READS_SINGLE, "n": 50, "seed": 3}, up="iuc/bbtools",
   summary="Draw a random subset of reads (mates kept in sync).")
def subsample_reads(src, n=0, fraction=0.1, seed=1, src2=""):
    """``bbduk.sh sample=N`` / ``seqtk sample`` for one or two mate files."""
    rs = _reads(src)
    k = int(n) if n else max(1, int(len(rs) * fraction))
    idx = sorted(random.Random(seed).sample(range(len(rs)), min(k, len(rs))))
    picked = [rs[i] for i in idx]
    res = _write_reads(picked, "sample.fastq", f"sampled {len(picked)}/{len(rs)} reads")
    if src2:
        r2 = _reads(src2)
        res["text2"] = io.write_fastq([r2[i] for i in idx if i < len(r2)])
        res["message"] += f" + {len([i for i in idx if i < len(r2)])} mates"
    return res


@T("fastq_normalise", "Digital normalisation of k-mer coverage", QC, "text",
   [fq("src"), intin("k", 21), number("target_depth", 5.0)],
   ex={"src": READS_SINGLE, "k": 21, "target_depth": 3}, up="iuc/bbnorm",
   summary="Discard reads whose median k-mer coverage exceeds a target (BBNorm-like).")
def digital_normalise(src, k=21, target_depth=5.0):
    """k-mer based normalisation to even out coverage."""
    rs = _reads(src)
    kk = max(1, int(k))
    cov: Counter = Counter()
    for r in rs:
        s = r.seq.upper()
        for i in range(max(0, len(s) - kk + 1)):
            cov[s[i:i + kk]] += 1
    out, depths = [], []
    for r in rs:
        s = r.seq.upper()
        ks = [s[i:i + kk] for i in range(max(1, len(s) - kk + 1))]
        d = stats.median([cov[x] for x in ks]) if ks else 0
        depths.append(d)
        if d <= target_depth:
            out.append(r)
    return _write_reads(out, "normalised.fastq",
                       f"normalised to <= {target_depth}x: {len(out)}/{len(rs)} reads kept, "
                       f"median k-mer depth {stats.median(depths) if depths else 0}")


@T("fastq_error_rate", "Estimated error rate from qualities", QC, "stats", [fq("src")],
   ex={"src": READS_SINGLE}, up="iuc/fastp",
   summary="Convert Phred scores into an expected substitution error rate.")
def error_rate(src):
    """Sum of per-base error probabilities."""
    rs = _reads(src)
    qs = [q for r in rs for q in _qual(r)]
    if not qs:
        return values_message({"bases": 0}, "no qualities")
    exp_errors = sum(10 ** (-q / 10) for q in qs)
    bases = len(qs)
    return values_message({"bases": bases, "reads": len(rs), "mean_quality": round(stats.mean(qs), 3),
                          "expected_substitutions": round(exp_errors, 3),
                          "error_rate_percent": round(100 * exp_errors / bases, 5),
                          "mean_error_probability": round(exp_errors / bases, 8)},
                         f"expected {exp_errors:.1f} sequencing errors")


@T("fastq_prinseq", "Prinseq-like multi-criteria filter", QC, "text",
   [fq("src"), intin("min_len", 50), number("max_ns_frac", 1.0), intin("min_qual_mean", 0),
    number("trim_qual_left", 0.0), number("trim_qual_right", 0.0), intin("derep", 0),
    number("dust_threshold", 0.0), intin("max_homopolymer", 0)],
   ex={"src": READS_SINGLE, "min_len": 50, "min_qual_mean": 20, "max_homopolymer": 12},
   up="du_novo/prinseq", summary="Prinseq++-style combined quality/complexity filters.")
def prinseq_lite(src, min_len=50, max_ns_frac=1.0, min_qual_mean=0, trim_qual_left=0.0,
                 trim_qual_right=0.0, derep=0, dust_threshold=0.0, max_homopolymer=0):
    """Prinseq's filter combination in one tool."""
    rs = _reads(src)
    out, seen = [], Counter()
    for r in rs:
        s, q = r.seq.upper(), _qual(r)
        while trim_qual_left and q and q[0] < trim_qual_left:
            q, s = q[1:], s[1:]
        while trim_qual_right and q and q[-1] < trim_qual_right:
            q, s = q[:-1], s[:-1]
        if len(s) < int(min_len) or not s:
            continue
        if s.count("N") / max(1, len(s)) > max_ns_frac:
            continue
        if min_qual_mean and q and stats.mean(q) < min_qual_mean:
            continue
        if max_homopolymer and seq.longest_homopolymer(s).get("length", 0) > int(max_homopolymer):
            continue
        if dust_threshold and seq.dust_score(s) < dust_threshold:
            continue
        if derep:
            seen[s] += 1
            if seen[s] > int(derep):
                continue
        out.append(io.Read(r.id, s, "".join(chr(min(93, v + 33)) for v in q), r.desc))
    return _write_reads(out, "prinseq.fastq", f"{len(out)}/{len(rs)} reads pass prinseq filters")


@T("fastq_fastp", "fastp-like all-in-one QC", QC, "stats",
   [fq("src"), fq("src2", ""), intin("qualifier", 15), intin("unqualified_percent_limit", 40),
    intin("length_required", 15), textbox("adapter_sequence", ""), boolean("dedup", False),
    intin("cut_window_size", 4), intin("cut_mean_quality", 20)],
   ex={"src": READS_1, "src2": READS_2, "adapter_sequence": "AGATCGGAAGAGC"}, up="iuc/fastp",
   summary="Adapter trim + quality filter + duplication stats in one pass.")
def fastp_lite(src, src2="", qualifier=15, unqualified_percent_limit=40, length_required=15,
               adapter_sequence="", dedup=False, cut_window_size=4, cut_mean_quality=20):
    """Report of a fastp-style run (counts only, no alignment)."""
    rs = _reads(src)
    r2 = _reads(src2) if src2 else []
    adapter = (adapter_sequence or "AGATCGGAAGAGC").upper()
    passed, trimmed_reads, quality_trimmed, lq, dup = 0, 0, 0, 0, 0
    kept = []
    for r in rs:
        s, q = r.seq.upper(), _qual(r)
        i = s.find(adapter[:12])
        if i >= 0:
            s, q = s[:i], q[:i]
            trimmed_reads += 1
        w = max(1, int(cut_window_size))
        while len(q) > w and stats.mean(q[-w:]) < cut_mean_quality:
            q, s = q[:-1], s[:-1]
            quality_trimmed += 1
        if q:
            bad = sum(1 for v in q if v < qualifier)
            if 100 * bad / len(q) > unqualified_percent_limit:
                lq += 1
                continue
        if len(s) < length_required:
            continue
        passed += 1
        kept.append(io.Read(r.id, s, "".join(chr(min(93, v + 33)) for v in q), r.desc))
    if dedup:
        c = Counter(x.seq for x in kept)
        dup = sum(v - 1 for v in c.values() if v > 1)
        kept = list({x.seq: x for x in reversed(kept)}.values())
    out = {"reads_in": len(rs), "reads_in_mate2": len(r2),
           "bases_in": sum(len(x.seq) for x in rs), "reads_passed_filter": passed,
           "reads_failed_low_quality": lq, "adapter_trimmed_reads": trimmed_reads,
           "quality_trimmed_bases": quality_trimmed, "duplicates_removed": dup,
           "reads_written": len(kept), "bases_written": sum(len(x.seq) for x in kept),
           "percent_passed": round(100 * passed / max(1, len(rs)), 2)}
    res = _write_reads(kept, "fastp.out.fastq", f"fastp-style: {passed}/{len(rs)} reads pass")
    res["stats"] = out
    return res


@T("fastq_fastx_toolkit_qc", "FASTX-Toolkit per-cycle statistics", QC, "table",
   [fq("src"), choice("stat", ["quality", "gc", "length", "maxquality", "minquality", "ap", "cp"])],
   ex={"src": READS_SINGLE, "stat": "quality"}, up="iuc/fastx_toolkit",
   summary="``fastx_quality_stats`` per-cycle statistics for one metric.")
def fastx_stats(src, stat="quality"):
    """Per-cycle statistics (quality, GC, A/C percentages)."""
    rs = _reads(src)
    if not rs:
        return table([], "no reads")
    maxlen = max(len(r.seq) for r in rs)
    out = []
    for i in range(maxlen):
        col = [r.seq[i].upper() if i < len(r.seq) else "" for r in rs]
        quals = [_qual(r)[i] if i < len(r.qual) else 0 for r in rs]
        n = sum(1 for c in col if c) or 1
        if stat == "quality":
            val, extra = stats.mean(quals), {"min": min(quals), "max": max(quals)}
        elif stat == "gc":
            val, extra = 100 * sum(1 for c in col if c in "GC") / n, {}
        elif stat == "length":
            val, extra = float(n), {}
        elif stat == "maxquality":
            val, extra = max(quals), {}
        elif stat == "minquality":
            val, extra = min(quals), {}
        else:
            want = {"ap": "A", "cp": "C"}.get(stat, "A")
            val, extra = 100 * sum(1 for c in col if c == want) / n, {}
        out.append({"cycle": i + 1, "reads": n, stat: round(float(val), 4), **extra})
    return table(out, f"fastx '{stat}' statistics over {maxlen} cycles")


@T("fastq_bbsplit_screen", "Read screening against reference sequences", QC, "table",
   [fq("src"), fa("references", GENOME), intin("k", 21), boolean("trim", False),
    intin("mismatches", 3)], ex={"src": READS_SINGLE, "references": GENOME},
   up="iuc/bbsplit", summary="Assign each read to the best reference by k-mer seeding.")
def screen_reads(src, references, k=21, trim=False, mismatches=3):
    """BBMap-style ``bbsplit`` assignment table."""
    ref = io.parse_fasta(io.as_text(references))
    reads = _reads(src)
    mapped = align.map_reads([{"id": r.id, "seq": r.seq} for r in ref],
                             [{"id": r.id, "seq": r.seq} for r in reads],
                             k=int(k), max_mismatches=int(mismatches))
    byref: Counter = Counter()
    ident: dict[str, list[float]] = defaultdict(list)
    unmapped = 0
    for m in mapped:
        if m.get("rname") and m["rname"] != "*":
            byref[m["rname"]] += 1
            ident[m["rname"]].append(m.get("identity", 0.0))
        else:
            unmapped += 1
    out = [{"reference": r.id, "reads": byref.get(r.id, 0),
            "percent": round(100 * byref.get(r.id, 0) / max(1, len(reads)), 3),
            "mean_identity": round(stats.mean(ident[r.id]), 4) if ident.get(r.id) else "",
            "length": len(r.seq)} for r in ref]
    out.append({"reference": "unmapped", "reads": unmapped,
                "percent": round(100 * unmapped / max(1, len(reads)), 3), "mean_identity": "",
                "length": 0})
    return table(out, f"screened {len(reads)} reads against {len(ref)} references")


@T("fastq_qc_duplicate_report", "Duplicate rate from a read set", QC, "table",
   [fq("src", READS_1), fq("src2", ""), intin("kmer_for_key", 0)], ex={"src": READS_SINGLE},
   up="iuc/fastp", summary="Library complexity: unique vs duplicated read pairs.")
def duplicate_report(src, src2="", kmer_for_key=0):
    """Duplicate statistics over reads (or over mate pairs)."""
    rs = _reads(src)
    r2 = _reads(src2) if src2 else []
    k = int(kmer_for_key)

    def key(r, other=None):
        s = r.seq.upper()
        if k:
            s = s[:k]
        return (s, other.seq.upper()[:k] if other and k else other.seq.upper() if other else "")
    c: Counter = Counter()
    for i, r in enumerate(rs):
        c[key(r, r2[i] if i < len(r2) else None)] += 1
    total = sum(c.values())
    dup = sum(v - 1 for v in c.values() if v > 1)
    return table([{"metric": "reads", "value": total},
                  {"metric": "unique", "value": len(c)},
                  {"metric": "duplicates", "value": dup},
                  {"metric": "duplicate_percent", "value": round(100 * dup / max(1, total), 3)},
                  {"metric": "distinct_pairs", "value": len(c)},
                  {"metric": "library_complexity", "value": round(len(c) / max(1, total), 4)}],
                 f"{dup} duplicates among {total} reads")


# ===========================================================================
# BBDTools-style extras
# ===========================================================================
@T("fastq_bbnorm_qc", "Quality-based read binning", "bbtools", "table",
   [fq("src"), intin("bin_size", 5)], ex={"src": READS_SINGLE}, up="iuc/bbtools",
   summary="Group reads by mean quality bin (BBNorm-style diagnostics).")
def quality_bins(src, bin_size=5):
    """Per-quality-bin read counts."""
    rs = _reads(src)
    b = max(1, int(bin_size))
    c: Counter = Counter()
    for r in rs:
        mq = int(stats.mean(_qual(r)) // b) * b if r.qual else 0
        c[mq] += 1
    return table([{"quality_bin": f"{k}-{k + b - 1}", "reads": v,
                   "percent": round(100 * v / max(1, len(rs)), 3)} for k, v in sorted(c.items())],
                 f"{len(rs)} reads in {len(c)} bins of {b}")


@T("fastq_bbmerge_report", "Merge mate pairs and report rates", "bbtools", "table",
   [fq("src", READS_1), fq("src2", READS_2), intin("min_overlap", 11), intin("max_insert", 1000)],
   ex={"src": READS_1, "src2": READS_2}, up="iuc/bbmerge",
   summary="``bbmerge.sh``-style overlap merging statistics.")
def bbmerge_report(src, src2, min_overlap=11, max_insert=1000):
    """Report merging rate and merged length of a pair set."""
    r1, r2 = _reads(src), _reads(src2)
    res = align.merge_pairs([{"id": x.id.split("/")[0], "seq": x.seq, "qual": x.qual} for x in r1],
                            [{"id": x.id.split("/")[0], "seq": x.seq, "qual": x.qual} for x in r2],
                            min_overlap=int(min_overlap))
    merged = res.get("merged", [])
    lens = [len(m["seq"]) for m in merged]
    return table([{"metric": "pairs_in", "value": min(len(r1), len(r2))},
                  {"metric": "merged", "value": len(merged)},
                  {"metric": "unmerged", "value": res.get("n_unmerged", 0)},
                  {"metric": "merge_rate_percent", "value": round(100 * res.get("merge_rate", 0), 2)},
                  {"metric": "mean_merged_length", "value": round(stats.mean(lens), 1) if lens else 0},
                  {"metric": "mean_overlap", "value": round(res.get("mean_overlap", 0), 2)},
                  {"metric": "max_insert_expected", "value": int(max_insert)}],
                 f"bbmerge: {len(merged)} merged pairs")


@T("fastq_bbsim_shard", "Shard reads into N files", "bbtools", "table",
   [fq("src"), intin("shards", 4), intin("seed", 1)], ex={"src": READS_SINGLE, "shards": 3},
   up="iuc/bbtools", summary="``shards.sh``: split a FASTQ into balanced random shards.")
def shard_reads(src, shards=4, seed=1):
    """Randomised sharding report of a read set."""
    rs = _reads(src)
    n = max(1, int(shards))
    random.Random(seed).shuffle(rs)
    buckets = [[] for _ in range(n)]
    for i, r in enumerate(rs):
        buckets[i % n].append(r)
    return table([{"shard": i + 1, "reads": len(b), "bases": sum(len(x.seq) for x in b),
                   "mean_length": round(stats.mean([len(x.seq) for x in b]), 1) if b else 0}
                  for i, b in enumerate(buckets)], f"sharded {len(rs)} reads into {n} files")


@T("fastq_bbinsert_check", "Insert-size and pair sanity check", "bbtools", "stats",
   [sam("src", ALN_SAM), number("sd_cutoff", 3.0)], ex={"src": ALN_SAM}, up="iuc/bbtools",
   summary="Flag abnormal fragment lengths in a paired alignment file.")
def insert_size_check(src, sd_cutoff=3.0):
    """``reapr``/BBTools-style insert size QC."""
    hdr, alns = bam.load(src)
    st = bam.insert_size_stats(alns, sd_cutoff=float(sd_cutoff))
    out = {k: (round(v, 4) if isinstance(v, float) else v) for k, v in st.items()}
    return values_message(out, f"insert size mean {out.get('mean')} sd {out.get('sd')}")


@T("fastq_bblearn_errors", "Error model from qualities", "bbtools", "table",
   [fq("src"), intin("k", 1)], ex={"src": READS_SINGLE}, up="iuc/bbduk",
   summary="Empirical ambiguity rate per quality bin (``errorsfromfastq``).")
def error_model(src, k=1):
    """Bin quality vs. base ambiguity as a proxy error model."""
    rs = _reads(src)
    rows: dict[int, list[int]] = defaultdict(lambda: [0, 0])
    for r in rs:
        s = r.seq.upper()
        for i, c in enumerate(s):
            if i >= len(r.qual):
                break
            q = ord(r.qual[i]) - 33
            rows[q][0] += 1
            if c not in "ACGT":
                rows[q][1] += 1
    out = []
    for q in sorted(rows):
        n, amb = rows[q]
        out.append({"quality": q, "bases": n, "ambiguous": amb,
                    "observed_error_rate": round((amb or 1e-6) / n, 6),
                    "implied_error_rate": round(10 ** (-q / 10), 6)})
    return table(out, f"error model over {len(out)} quality bins")


@T("fastq_bbduk_contaminant", "Scrub contaminant k-mers", "bbtools", "text",
   [fq("src"), bigtext("contaminants", "GGGGGGGGGGGGGGGG\nAAAAAAAAAAAAAAAAAA",
                       label="Contaminant sequences (one per line)"), intin("k", 16),
    intin("mismatches", 0), boolean("trim", True)],
   ex={"src": READS_SINGLE, "contaminants": "AGATCGGAAGAGC", "k": 12}, up="iuc/bbduk",
   summary="``bbduk.sh ktrim=r k=…``: drop or trim reads containing contaminant k-mers.")
def bbduk_contaminant(src, contaminants="", k=16, mismatches=0, trim=True):
    """k-mer based scrubbing of a read set."""
    rs = _reads(src)
    kk = max(4, int(k))
    bad = set()
    for line in io.lines(contaminants):
        s = line.strip().upper()
        for i in range(max(1, len(s) - kk + 1)):
            bad.add(s[i:i + kk])
    out, dropped, trimmed = [], 0, 0
    for r in rs:
        s = r.seq.upper()
        hit = next((i for i in range(max(0, len(s) - kk + 1)) if s[i:i + kk] in bad), None)
        if hit is None:
            out.append(r)
            continue
        if trim and hit > 10:
            out.append(io.Read(r.id, s[:hit], r.qual[:hit], r.desc))
            trimmed += 1
        else:
            dropped += 1
    return _write_reads(out, "scrubbed.fastq",
                       f"{dropped} reads removed, {trimmed} trimmed by {len(bad)} contaminant k-mers")


# ===========================================================================
# Nanopore
# ===========================================================================
@T("nanopore_stats", "Nanopore read statistics", "nanopore", "stats",
   [fq("src", READS_SINGLE), intin("min_length", 200)], ex={"src": READS_SINGLE, "min_length": 50},
   up="iuc/nanoq", summary="Read count, total bases, N50 and quality of long reads.")
def nanopore_stats(src, min_length=200):
    """``nanoq``-like summary of a long-read fastq."""
    rs = [r for r in _reads(src) if len(r.seq) >= int(min_length)]
    if not rs:
        return values_message({"reads": 0, "min_length": int(min_length)}, f"no reads >= {min_length} bp")
    lens = sorted((len(r.seq) for r in rs), reverse=True)
    qs = [stats.mean(_qual(r)) for r in rs if r.qual]
    exp_err = stats.mean([10 ** (-q / 10) for q in qs]) if qs else 0.0
    return values_message({"total_reads": len(rs), "total_bases": sum(lens), "longest_read": max(lens),
                          "shortest_read": min(lens), "mean_read_length": round(stats.mean(lens), 1),
                          "median_read_length": stats.median(lens), "n50": seq.n50(lens),
                          "mean_quality": round(stats.mean(qs), 2) if qs else 0,
                          "estimated_identity_percent": round(100 * (1 - exp_err), 3)},
                         f"{len(rs)} long reads, N50 {seq.n50(lens)}")


@T("nanopore_porechop", "Porechop-style adapter splitting", "nanopore", "table",
   [fq("src", READS_SINGLE), textbox("adapter_5p", "AATGTACTTCGTTCAGTTACGTATTGCT"),
    textbox("adapter_3p", "GCAATACGTAACTGAACGAAGT"), intin("min_subseq", 50)],
   ex={"src": READS_SINGLE}, up="iuc/porechop",
   summary="Detect adapter matches and report the trimmed length of each read.")
def porechop_like(src, adapter_5p="", adapter_3p="", min_subseq=50):
    """Simple adapter detection with per-read cut positions."""
    rs = _reads(src)
    out, n_adapter = [], 0
    for r in rs:
        s = r.seq.upper()
        cut5 = s.find(adapter_5p[:16]) if adapter_5p else -1
        cut3 = s.rfind(adapter_3p[:16]) if adapter_3p else -1
        if cut5 >= 0 or cut3 >= 0:
            n_adapter += 1
        keep = s[max(cut5 + 1, 0):cut3 if cut3 > 0 else len(s)]
        out.append({"read": r.id, "length": len(r.seq), "adapter_5p": cut5, "adapter_3p": cut3,
                    "trimmed_length": len(keep), "kept": len(keep) >= int(min_subseq)})
    return table(out[:500], f"{n_adapter}/{len(rs)} reads contain adapter sequence")


@T("nanopore_quality_by_length", "Quality vs read length", "nanopore", "table",
   [fq("src", READS_SINGLE), intin("bins", 10)], ex={"src": READS_SINGLE}, up="iuc/nanoq",
   summary="Binned mean quality as a function of read length.")
def quality_by_length(src, bins=10):
    """Length-binned quality table for long reads."""
    rs = _reads(src)
    pts = [(len(r.seq), stats.mean(_qual(r))) for r in rs if r.qual]
    if not pts:
        return table([], "no data")
    lo, hi = min(p[0] for p in pts), max(p[0] for p in pts)
    nb = max(1, int(bins))
    w = ((hi - lo) / nb) or 1
    groups: dict[int, list[float]] = defaultdict(list)
    for ln, q in pts:
        groups[int((ln - lo) // w)].append(q)
    return table([{"length_bin": f"{lo + k * w:.0f}-{lo + (k + 1) * w:.0f}", "reads": len(v),
                   "mean_quality": round(stats.mean(v), 3),
                   "median_quality": round(stats.median(v), 2)} for k, v in sorted(groups.items())],
                 f"{len(pts)} reads in {len(groups)} bins")


@T("nanopore_kmer_quals", "Nanopore k-mer quality bias", "nanopore", "table",
   [fq("src", READS_SINGLE), intin("k", 5), intin("top", 15)], ex={"src": READS_SINGLE, "k": 4},
   up="iuc/nanopore", summary="Mean quality per k-mer context (homopolymer/systematic errors).")
def kmer_quality_bias(src, k=5, top=15):
    """Context-dependent error signature table."""
    rs = _reads(src)
    kk = max(1, int(k))
    acc: dict[str, list[float]] = defaultdict(list)
    for r in rs:
        s = r.seq.upper()
        q = _qual(r)
        for i in range(max(0, len(s) - kk + 1)):
            acc[s[i:i + kk]].append(stats.mean(q[i:i + kk]))
    out = [{"kmer": km, "observations": len(v), "mean_quality": round(stats.mean(v), 3),
            "homopolymer": len(set(km)) == 1}
           for km, v in sorted(acc.items(), key=lambda kv: -len(kv[1]))[: int(top)]]
    return table(out, f"quality by {kk}-mer context")


@T("nanopore_filter", "Filter long reads by length and quality", "nanopore", "text",
   [fq("src", READS_SINGLE), intin("min_length", 100), number("min_mean_quality", 7.0),
    intin("max_length", 0)], ex={"src": READS_SINGLE, "min_length": 50}, up="iuc/filtlong",
   summary="``filtlong``-style selection of reads.")
def nanopore_filter(src, min_length=100, min_mean_quality=7.0, max_length=0):
    """Length/quality filter for nanopore data."""
    rs = _reads(src)
    keep = [r for r in rs if len(r.seq) >= int(min_length)
            and (not max_length or len(r.seq) <= int(max_length))
            and (not r.qual or stats.mean(_qual(r)) >= min_mean_quality)]
    return _write_reads(keep, "filtered.fastq",
                       f"{len(keep)}/{len(rs)} reads pass (>= {min_length} bp, Q >= {min_mean_quality})")


@T("nanopore_read_table", "Read-level summary table", "nanopore", "table",
   [fq("src", READS_SINGLE), intin("max_rows", 500)], ex={"src": READS_SINGLE}, up="iuc/poretools",
   summary="One row per read: length, quality, GC, Ns, homopolymers.")
def read_table(src, max_rows=500):
    """Per-read diagnostics table (``poretools``-like)."""
    out = []
    for r in _reads(src)[: int(max_rows)]:
        q = _qual(r)
        out.append({"read": r.id, "length": len(r.seq),
                    "mean_quality": round(stats.mean(q), 3) if q else "",
                    "min_quality": min(q) if q else "",
                    "gc_percent": round(100 * seq.gc_content(r.seq), 3),
                    "n_bases": r.seq.upper().count("N"),
                    "longest_homopolymer": seq.longest_homopolymer(r.seq.upper()).get("length", 0)})
    return table(out, f"{len(out)} reads summarised")


@T("nanopore_fast5_proxy", "Signal-level proxy statistics from FASTQ", "nanopore", "table",
   [fq("src", READS_SINGLE), intin("window", 25)], ex={"src": READS_SINGLE}, up="iuc/nanopack",
   summary="Quality moving average along each read (NanoPack-style diagnostics).")
def signal_proxy(src, window=25):
    """Rolling mean quality per read, as a stand-in for raw-signal plots."""
    rs = _reads(src)
    w = max(2, int(window))
    out = []
    for r in rs[:100]:
        q = _qual(r)
        prof = [round(stats.mean(q[i:i + w]), 2) for i in range(0, max(1, len(q) - w + 1), w)]
        out.append({"read": r.id, "length": len(r.seq), "n_windows": len(prof),
                    "min_window_q": min(prof) if prof else "", "max_window_q": max(prof) if prof else "",
                    "profile": ",".join(map(str, prof[:40]))})
    return table(out, f"quality profile in {w} bp windows")
