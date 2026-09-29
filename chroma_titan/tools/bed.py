"""Genomic interval tools: bedtools-style operations, interval statistics,
lift-over, and conversions in/out of BED."""

from __future__ import annotations

import math
import random
import re
from collections import Counter, defaultdict

from chroma_titan.core import genome as _G
from chroma_titan.core import io, seq, stats, tables
from chroma_titan.core import genome
from chroma_titan.tools._common import *  # noqa: F401,F403
from chroma_titan.tools._common import (ANNOT_GFF, BEDGRAPH, GENES_FA, REGIONS_BED, TARGETS_BED,
                                        VARIANTS_VCF)

BED_SEC = "bed"
OPS = "operate_on_genomic_intervals"


def _ivs(src) -> list:
    return io.parse_bed(io.as_text(src))


def _g(src):
    """Genome dictionary {chrom: sequence} + sizes from a FASTA/text source."""
    seqs, sizes = genome.read_genome(src)
    return seqs, sizes


def _bed(ivs, message=""):
    clean = []
    for x in ivs:
        sc = x.score
        if isinstance(sc, float) and (math.isnan(sc) or math.isinf(sc)):
            sc = 0
        clean.append(io.Interval(x.chrom, int(x.start), int(x.end), x.name or ".", sc,
                                x.strand or ".", list(x.extra or [])))
    return {"text": io.write_bed(clean, bed12=False), "filename": "output.bed",
            "message": message or f"{len(clean)} intervals"}


# ===========================================================================
# core bedtools operations
# ===========================================================================
@T("bedtools_intersect", "Intersect", OPS, "table",
   [bed("a"), bed("b"), choice("wa", ["none", "wa", "wb", "both"]), boolean("strand", False),
    choice("report", ["rows", "unique_a", "count_a", "invert"])],
   ex={"a": REGIONS_BED, "b": TARGETS_BED, "wa": "both"}, up="iuc/bedtools",
   summary="bedtools intersect: overlaps between two interval files.")
def bed_intersect(a, b, wa="none", strand=False, report="rows"):
    """Interval intersection with the usual reporting modes."""
    ra, rb = _ivs(a), _ivs(b)
    hits = genome.intersect(ra, rb, wa=wa in ("wa", "both"), wb=wa in ("wb", "both"), strand=strand)
    if report == "invert":
        keep = genome.subtract(ra, rb, strand=strand)
        out = genome.to_table(keep)
        return {"table": out, "rows": len(out), "text": tables.to_tsv(out),
                "filename": "result.tsv", "message": f"{len(out)} A intervals with no overlap in B"}
    if report == "unique_a":
        uniq = list(dict.fromkeys(h.name or f"{h.chrom}:{h.start}-{h.end}" for h in hits))
        return table([{"chrom": n.split(":")[0], "interval": n} for n in uniq],
                     f"{len(uniq)} distinct A intervals overlap")
    if report == "count_a":
        c = Counter(h.name for h in hits)
        rows = [{"A_interval": iv.name or f"{iv.chrom}:{iv.start + 1}-{iv.end}",
                 "chrom": iv.chrom, "start": iv.start, "end": iv.end,
                 "n_overlaps": c.get(iv.name or "", 0)} for iv in ra]
        return table(rows, f"{sum(r['n_overlaps'] for r in rows)} intersections")
    return table(genome.to_table(hits).to_dict("records") if hits else [],
                 f"{len(hits)} intersections from {len(ra)} x {len(rb)} intervals")


@T("bedtools_intersect_intervals_only", "Intersect intervals (geometry only)", OPS, "bed",
   [bed("a"), bed("b"), boolean("strand", False)], ex={"a": REGIONS_BED, "b": TARGETS_BED},
   up="iuc/bedtools", summary="bedtools intersect -u: A intervals that overlap any B interval.")
def bed_intersect_u(a, b, strand=False):
    """Return only the A records with at least one overlap."""
    ra, rb = _ivs(a), _ivs(b)
    hits = {id(x) for x in genome.intersect(ra, rb, strand=strand)}
    keep = [x for x in ra if id(x) in hits]
    return _bed(keep, f"{len(keep)}/{len(ra)} A intervals overlap B")


@T("bedtools_intersect_v", "Intersect -v (only non-overlapping)", OPS, "bed",
   [bed("a"), bed("b"), boolean("strand", False)], ex={"a": REGIONS_BED, "b": TARGETS_BED},
   up="iuc/bedtools", summary="bedtools intersect -v: A intervals with no overlap in B.")
def bed_intersect_v(a, b, strand=False):
    """Anti-filter of A by B (equivalent to subtract -A)."""
    ra, rb = _ivs(a), _ivs(b)
    keep = genome.subtract(ra, rb, strand=strand)
    return _bed(keep, f"{len(keep)}/{len(ra)} A intervals do not overlap B")


@T("bedtools_union_bed_graphs", "Union of multiple interval files", OPS, "table",
   [multi("files", [REGIONS_BED, TARGETS_BED], [REGIONS_BED, TARGETS_BED])],
   ex={"files": [REGIONS_BED, TARGETS_BED]}, up="iuc/bedtools",
   summary="bedtools unionbedg-style coverage matrix over all intervals.")
def union_bed_graphs(files=None):
    """Split the union of all intervals at every boundary, then score each piece."""
    sets = [_ivs(f) for f in (files or [])]
    by_chrom: dict[str, list[tuple[int, int]]] = defaultdict(list)
    for lst in sets:
        for iv in lst:
            by_chrom[iv.chrom].append((iv.start, iv.end))
    rows = []
    for chrom in sorted(by_chrom):
        bounds = sorted({x for a, b in by_chrom[chrom] for x in (a, b)})
        for k in range(len(bounds) - 1):
            lo, hi = bounds[k], bounds[k + 1]
            if hi <= lo:
                continue
            rec = {"chrom": chrom, "start": lo, "end": hi, "length": hi - lo}
            for idx, lst in enumerate(sets):
                ov = sum(max(0, min(hi, x.end) - max(lo, x.start)) for x in lst if x.chrom == chrom)
                rec[f"file{idx + 1}_coverage"] = round(ov / (hi - lo), 4)
            rows.append(rec)
    return table(rows, f"union of {len(sets)} files -> {len(rows)} pieces")


@T("bedtools_multiinter", "MultiIntersect", OPS, "table",
   [multi("files", [REGIONS_BED, TARGETS_BED], [REGIONS_BED, TARGETS_BED]), boolean("header", True)],
   ex={"files": [REGIONS_BED, TARGETS_BED]}, up="iuc/bedtools",
   summary="bedtools multiinter: which files cover each genomic bin.")
def bed_multiinter(files=None, header=True):
    """Overlapping-set table across several interval files."""
    sets = [_ivs(f) for f in (files or [])]
    rows = genome.multiinter(sets)
    return table(rows, f"{len(rows)} regions from {len(sets)} files")


@T("bedtools_subtract", "Subtract", OPS, "bed", [bed("a"), bed("b"), boolean("whole_A", False)],
   ex={"a": REGIONS_BED, "b": TARGETS_BED}, up="iuc/bedtools",
   summary="bedtools subtract: remove the parts of A covered by B.")
def bed_subtract(a, b, whole_A=False):
    """Interval subtraction (``-A`` keeps whole lines when asked)."""
    ra, rb = _ivs(a), _ivs(b)
    out = ra if whole_A else genome.subtract(ra, rb)
    if whole_A:
        covered = {id(x) for x in genome.intersect(ra, rb)}
        out = [x for x in ra if id(x) not in covered]
    return _bed(out, f"{len(out)} intervals after subtracting B from A")


@T("bedtools_merge", "Merge", BED_SEC, "bed",
   [bed("src"), intin("distance", -1), boolean("strand", False),
    choice("collapse", ["none", "concat", "count", "mean", "min", "max", "sum", "distinct"])],
   ex={"src": REGIONS_BED, "distance": 100, "collapse": "count"}, up="iuc/bedtools",
   summary="bedtools merge with distance and column summaries.")
def bed_merge(src, distance=-1, strand=False, collapse="none"):
    """Merge neighbouring/overlapping intervals."""
    ivs = _ivs(src)
    out = genome.merge(ivs, distance=int(distance), strand=strand,
                       collapse="none" if collapse == "none" else collapse)
    return _bed(out, f"merged {len(ivs)} -> {len(out)} intervals (d={distance})")


@T("bedtools_cluster", "Cluster", BED_SEC, "table", [bed("src"), intin("distance", 1)],
   ex={"src": REGIONS_BED, "distance": 50}, up="iuc/bedtools",
   summary="bedtools cluster: assign a cluster id to intervals close together.")
def bed_cluster(src, distance=1):
    """Cluster overlapping/nearby intervals."""
    rows = genome.cluster(_ivs(src), distance=int(distance))
    return table(rows, f"{len(rows)} intervals in {len({r['cluster'] for r in rows})} clusters")


@T("bedtools_links", "Links", BED_SEC, "table", [bed("src"), intin("distance", 1)],
   ex={"src": REGIONS_BED, "distance": 200}, up="iuc/bedtools",
   summary="bedtools links: pairs of intervals within a distance.")
def bed_links(src, distance=1):
    """All links between clustered intervals."""
    rows = genome.links(_ivs(src), distance=int(distance))
    return table(rows, f"{len(rows)} links")


@T("bedtools_window", "Window", OPS, "table", [bed("a"), bed("b"), intin("window", 500)],
   ex={"a": TARGETS_BED, "b": REGIONS_BED, "window": 1000}, up="iuc/bedtools",
   summary="bedtools window: B intervals within W bp of A intervals.")
def bed_window(a, b, window=500):
    """Fixed-window association between two interval sets."""
    rows = genome.window(_ivs(a), _ivs(b), w=int(window))
    return table(rows, f"{len(rows)} A-B pairs within {window} bp")


@T("bedtools_closest", "Closest", OPS, "table", [bed("a"), bed("b"), boolean("distance_only", False)],
   ex={"a": TARGETS_BED, "b": REGIONS_BED}, up="iuc/bedtools",
   summary="bedtools closest: the nearest B interval for each A interval.")
def bed_closest(a, b, distance_only=False):
    """Neighbour search with signed distance."""
    rows = genome.closest(_ivs(a), _ivs(b))
    if distance_only:
        rows = [{"A": r.get("A_name") or f"{r['A_chrom']}:{r['A_start']}",
                 "B": r.get("B_name") or f"{r['B_chrom']}:{r['B_start']}",
                 "distance": r.get("distance", 0)} for r in rows]
    return table(rows, f"{len(rows)} closest pairs")


@T("bedtools_distance", "Distance", BED_SEC, "table", [bed("a"), bed("b")],
   ex={"a": TARGETS_BED, "b": REGIONS_BED}, up="iuc/bedtools",
   summary="bedtools distance: pairwise genomic distances between two files.")
def bed_distance(a, b):
    """Signed distance matrix (sparse) between A and B intervals."""
    ra, rb = _ivs(a), _ivs(b)
    rows = []
    for x in ra:
        for y in rb:
            rows.append({"A": x.name or f"{x.chrom}:{x.start + 1}",
                         "B": y.name or f"{y.chrom}:{y.start + 1}", "chrom": x.chrom,
                         "distance": genome.distance_between(x, y)})
    rows.sort(key=lambda r: abs(r["distance"]))
    return table(rows, f"{len(rows)} pairwise distances")


@T("bedtools_complement", "Complement", BED_SEC, "bed", [bed("src"), anyfile("genome", GENOME_TXT,
                                                                           label="Genome file (chrom\\tsize)")],
   ex={"src": REGIONS_BED, "genome": GENOME_TXT}, up="iuc/bedtools",
   summary="bedtools complement: everything in the genome that is not covered.")
def bed_complement(src, genome=None):
    """Invert a BED file against a genome definition."""
    sizes = genome_file_sizes(genome)
    out = genome_complement(sizes, _ivs(src))
    return _bed(out, f"{len(out)} uncovered blocks")


def genome_file_sizes(src):
    return genome.parse_genome_file(io.as_text(src) if src else "", mode="sizes")


def genome_complement(sizes, ivs):
    return genome.complement(sizes, ivs)


@T("bedtools_coverage", "Coverage", BED_SEC, "table",
   [bed("a"), bed("b"), choice("mode", ["bases", "count", "fraction"])],
   ex={"a": REGIONS_BED, "b": TARGETS_BED}, up="iuc/bedtools",
   summary="bedtools coverage: overlap bases, counts and fraction per A interval.")
def bed_coverage(a, b, mode="bases"):
    """Per-interval coverage statistics."""
    rows = genome.coverage(_ivs(a), _ivs(b) or None)
    if mode == "count":
        rows = [{"chrom": r["chrom"], "start": r["start"], "end": r["end"], "name": r.get("name", ""),
                 "n_overlapping": r.get("n", 0)} for r in rows]
    elif mode == "fraction":
        rows = [{"chrom": r["chrom"], "start": r["start"], "end": r["end"], "name": r.get("name", ""),
                 "percent_coverage": round(100 * r.get("fraction", 0.0), 4)} for r in rows]
    return table(rows, f"coverage of {len(rows)} intervals")


@T("bedtools_jaccard", "Jaccard and Fisher tests", OPS, "stats", [bed("a"), bed("b")],
   ex={"a": REGIONS_BED, "b": TARGETS_BED}, up="iuc/bedtools",
   summary="bedtools jaccard plus the Fisher/ t-tests on interval overlap.")
def bed_jaccard(a, b):
    """Overlap similarity of two interval files."""
    st = genome.jaccard(_ivs(a), _ivs(b))
    return values_message({k: (round(v, 6) if isinstance(v, float) else v) for k, v in st.items()},
                         "interval Jaccard index")


@T("bedtools_slop", "Slop", BED_SEC, "bed",
   [bed("src"), intin("bp", 100), number("fraction", 0.0), choice("direction", ["both", "left", "right"]),
    boolean("genome_bound", False), anyfile("genome", GENOME_TXT, label="Genome file (optional)")],
   ex={"src": REGIONS_BED, "bp": 100, "direction": "both"}, up="iuc/bedtools",
   summary="bedtools slop: expand intervals, optionally clipped to the genome.")
def bed_slop(src, bp=100, fraction=0.0, direction="both", genome_bound=False, genome=None):
    """Grow intervals by a fixed amount or a fraction of their length."""
    sizes = genome_file_sizes(genome) if genome_bound and genome else None
    out = _G.slop(_ivs(src), sizes, bp=int(bp), direction=direction, fraction=fraction)
    return _bed(out, f"slopped {len(out)} intervals by {bp} bp ({direction})")


@T("bedtools_flank", "Flank", BED_SEC, "bed",
   [bed("src"), intin("left", 100), intin("right", 100), boolean("both", False)],
   ex={"src": REGIONS_BED, "left": 200, "right": 200}, up="iuc/bedtools",
   summary="bedtools flank: upstream/downstream regions relative to the strand.")
def bed_flank(src, left=100, right=100, both=False):
    """Strand-aware flanking intervals."""
    out = genome.flank(_ivs(src), left=int(left), right=int(right), both=both)
    return _bed(out, f"{len(out)} flanking intervals")


@T("bedtools_shift", "Shift", BED_SEC, "bed", [bed("src"), intin("by", 100), boolean("strand_aware", True)],
   ex={"src": REGIONS_BED, "by": 50}, up="iuc/bedtools", summary="bedtools shift: move every interval.")
def bed_shift(src, by=100, strand_aware=True):
    """Translate intervals along the genome."""
    ivs = _ivs(src)
    if not strand_aware:
        out = genome.shift(ivs, by=int(by))
    else:
        out = []
        for iv in ivs:
            d = int(by) * (-1 if iv.strand == "-" else 1)
            out.append(io.Interval(iv.chrom, max(0, iv.start + d), max(0, iv.end + d), iv.name,
                                   iv.score, iv.strand, iv.extra))
    return _bed(out, f"shifted {len(out)} intervals by {by} bp")


@T("bedtools_scale", "Scale", BED_SEC, "bed", [bed("src"), anyfile("genome", GENOME_TXT),
                                             number("factor", 1.5)], ex={"src": REGIONS_BED, "factor": 1.5},
   up="iuc/bedtools", summary="bedtools scale: scale coordinates relative to genome length.")
def bed_scale(src, genome=None, factor=1.5):
    """Scale interval positions (useful for genome-wide normalisation)."""
    sizes = genome_file_sizes(genome)
    ivs = []
    for iv in _ivs(src):
        span = iv.end - iv.start
        c = iv.chrom
        new_span = max(1, int(round(span * factor)))
        mid = (iv.start + iv.end) // 2
        lo = max(0, mid - new_span // 2)
        hi = min(sizes.get(c, lo + new_span), lo + new_span)
        ivs.append(io.Interval(c, lo, hi, iv.name, iv.score, iv.strand, iv.extra))
    return _bed(ivs, f"scaled {len(ivs)} intervals by {factor}")


@T("bedtools_shuffle", "Shuffle", BED_SEC, "bed",
   [bed("src"), anyfile("genome", GENOME_TXT), intin("seed", 42), bed("exclude", ""),
    boolean("include", False)], ex={"src": REGIONS_BED, "genome": GENOME_TXT}, up="iuc/bedtools",
   summary="bedtools shuffle: randomise interval positions inside the genome.")
def bed_shuffle(src, genome=None, seed=42, exclude="", include=False):
    """Random re-positioning, optionally avoiding excluded regions."""
    sizes = genome_file_sizes(genome)
    ivs = _ivs(src)
    excl = _ivs(exclude) if exclude else None
    out = _G.shuffle_in(sizes, ivs, seed=int(seed), excl=excl)
    return _bed(out, f"shuffled {len(out)} intervals (seed={seed})")


@T("bedtools_random", "Random", BED_SEC, "bed",
   [anyfile("genome", GENOME_TXT), intin("n", 10), intin("length", 200), intin("seed", 1),
    textbox("chrom", "")], ex={"genome": GENOME_TXT, "n": 8, "length": 250}, up="iuc/bedtools",
   summary="bedtools random: generate random intervals of a fixed width.")
def bed_random(genome=None, n=10, length=200, seed=1, chrom=""):
    """Synthetic interval set for null models."""
    sizes = genome_file_sizes(genome)
    out = _G.random_intervals(sizes, n=int(n), length=int(length), seed=int(seed),
                              chrom_choice=chrom or None)
    return _bed(out, f"{len(out)} random intervals of {length} bp")


@T("bedtools_binnify", "Binnify genome into equal bins", BED_SEC, "bed",
   [anyfile("genome", GENOME_TXT), intin("binsize", 1000)], ex={"genome": GENOME_TXT, "binsize": 1000},
   up="iuc/bedtools", summary="bedtools binnify: fixed-width bins covering every chromosome.")
def bed_binnify(genome=None, binsize=1000):
    """Tile the genome into equal bins."""
    sizes = genome_file_sizes(genome)
    out = _G.binnify(sizes, binsize=int(binsize))
    return _bed(out, f"{len(out)} bins of {binsize} bp")


@T("bedtools_tidy", "Tidy", BED_SEC, "bed", [bed("src"), intin("max_extra", 6)],
   ex={"src": REGIONS_BED}, up="iuc/bedtools", summary="bedtools tidy: fix coordinates and extra columns.")
def bed_tidy_tool(src, max_extra=6):
    """Repair a BED file (clip starts, drop invalid rows, cap extra columns)."""
    out = genome.tidy(_ivs(src), max_extra=int(max_extra))
    return _bed(out, f"tidied {len(out)} intervals")


@T("bedtools_negativeb", "Negative B (remove B from A, keep strand)", BED_SEC, "bed",
   [bed("a"), bed("b")], ex={"a": REGIONS_BED, "b": TARGETS_BED}, up="iuc/bedtools",
   summary="bedtools negativeb: A intervals minus any that overlap B.")
def bed_negativeb(a, b):
    """Whole-line subtraction of B from A."""
    ra, rb = _ivs(a), _ivs(b)
    covered = {f"{x.chrom}:{x.start}-{x.end}" for x in genome.intersect(ra, rb)}
    out = [x for x in ra if f"{x.chrom}:{x.start}-{x.end}" not in covered]
    return _bed(out, f"{len(out)} A intervals survive")


@T("bedtools_remove_overlaps", "Remove overlapping intervals", BED_SEC, "bed",
   [bed("src"), intin("min_overlap", 1)], ex={"src": REGIONS_BED}, up="iuc/bedtools",
   summary="Drop intervals that overlap another interval by at least N bp.")
def bed_remove_overlaps(src, min_overlap=1):
    """Greedy de-overlapping."""
    out = genome.remove_overlaps(_ivs(src), min_overlap=int(min_overlap))
    return _bed(out, f"{len(out)} non-overlapping intervals")


@T("bedtools_intersect_bp", "Intersect bp", BED_SEC, "stats", [bed("a"), bed("b")],
   ex={"a": REGIONS_BED, "b": TARGETS_BED}, up="iuc/bedtools",
   summary="bedtools intersect -loj style base-pair totals for the two files.")
def bed_intersect_bp(a, b):
    """Base-pair overlap accounting."""
    ra, rb = _ivs(a), _ivs(b)
    ov = sum(genome.intersect_intervals(x, y) for x in ra for y in rb if x.chrom == y.chrom)
    st = {"a_intervals": len(ra), "b_intervals": len(rb),
          "a_bases": sum(x.end - x.start for x in ra), "b_bases": sum(y.end - y.start for y in rb),
          "merged_a_bases": sum(x.end - x.start for x in genome.merge(ra)),
          "merged_b_bases": sum(y.end - y.start for y in genome.merge(rb)),
          "overlap_bases": ov,
          "percent_a_covered": round(100 * ov / max(1, sum(x.end - x.start for x in ra)), 3)}
    return values_message(st, f"{ov} overlapping base pairs")


@T("bedtools_getfasta", "Get Fasta (BED intervals to sequence)", "fetch_sequences_alignments", "fasta",
   [fa("genome", GENOME), bed("intervals", REGIONS_BED), choice("strand", ["same", "plus", "minus"]),
    choice("name_mode", ["bed", "sequence", "chrom_only"]), boolean("split", False)],
   ex={"genome": GENOME, "intervals": REGIONS_BED, "strand": "same"}, up="iuc/bedtools",
   summary="bedtools getfasta: extract the sequence of each interval.")
def bed_getfasta(genome, intervals, strand="same", name_mode="bed", split=False):
    """Interval sequences (BED12 blocks become separate records with split=True)."""
    seqs, sizes = _read_genome(genome)
    ivs = _ivs(intervals)
    out = genome_getfasta(seqs, ivs, strand=strand, name_mode=name_mode)
    recs = [io.Seq(iv.name, (iv.extra or [""])[0], f"{iv.chrom}:{iv.start + 1}-{iv.end}{iv.strand}")
            for iv in out]
    return {"text": io.write_fasta(recs), "filename": "sequences.fasta",
            "message": f"{len(recs)} sequences"}


def _read_genome(src):
    return genome.read_genome(src)


def genome_getfasta(seqs, ivs, strand="same", name_mode="bed"):
    return genome.getfasta(seqs, ivs, strand=strand, name_mode=name_mode)


def genome_read(src):
    return genome.read_genome(src)


@T("bedtools_nuc_content", "Nuc content of intervals", "fetch_sequences_alignments", "table",
   [fa("genome", GENOME), bed("intervals", REGIONS_BED)], ex={"genome": GENOME, "intervals": REGIONS_BED},
   up="iuc/bedtools", summary="Per-interval nucleotide composition (GC%, CpG obs/exp).")
def bed_nuc_content(genome, intervals):
    """Composition of every interval, including CpG observed/expected."""
    seqs, sizes = _read_genome(genome)
    ivs = _ivs(intervals)
    rows = _G.nuc_content(ivs, seqs)
    cut = genome_getfasta(seqs, ivs, strand="same", name_mode="bed")
    for r, iv in zip(rows, cut):
        piece = (iv.extra or [""])[0]
        oe = seq.cpg_obs_exp(piece)
        r["cpg_obs_exp"] = round(oe.get("ratio", 0.0), 3)
        r["gc_percent"] = r.get("GC", 0.0)
        r["longest_homopolymer"] = seq.longest_homopolymer(piece).get("length", 0)
    return table(rows, f"composition of {len(rows)} intervals")


@T("bedtools_bamtobed", "BAM-to-BED", BED_SEC, "bed", [sam("src"), boolean("split", True),
                                                     number("min_mapq", 20.0), intin("flag_filter", 3844)],
   ex={"src": ALN_SAM}, up="iuc/bedtools", summary="bedtools bamtobed: alignments as BED intervals.")
def bamtobed(src, split=True, min_mapq=20.0, flag_filter=3844):
    """Convert a SAM/BAM alignment list into BED."""
    from chroma_titan.core import bam

    hdr, alns = bam.load(src)
    out = []
    for a in alns:
        if a.flag & int(flag_filter) or a.mapq < min_mapq or a.rname in ("*", ""):
            continue
        ref_len = io.ref_length(a.cigar)
        out.append(io.Interval(a.rname, max(0, a.pos - 1), a.pos - 1 + ref_len, a.qname, a.mapq,
                              "-" if a.reverse else "+"))
    return _bed(genome.sort_ivs(out), f"{len(out)} alignment intervals")


@T("bedtools_bed12_blocks", "Explain a BED12 file (block stats)", BED_SEC, "table",
   [bed("src", ANNOT_GFF if False else REGIONS_BED), boolean("as_gtf", False)],
   ex={"src": REGIONS_BED}, up="iuc/bedtools", summary="Per-interval block summary (BED12 aware).")
def bed_block_stats(src, as_gtf=False):  # noqa: ARG001
    """Summarise block structure and score fields of a BED file."""
    ivs = _ivs(src)
    rows = []
    for iv in ivs:
        extra = iv.extra or []
        blocks = 0
        try:
            blocks = int(extra[2]) if len(extra) > 2 else 0
        except (TypeError, ValueError):
            blocks = 0
        rows.append({"chrom": iv.chrom, "start": iv.start, "end": iv.end, "name": iv.name,
                     "score": iv.score, "strand": iv.strand, "length": iv.end - iv.start,
                     "n_extra_columns": len(extra), "n_blocks": blocks})
    return table(rows, f"{len(rows)} BED records")


@T("bedtools_bed_graph_to_bag", "BED graph to genome-wide bins", BED_SEC, "table",
   [anyfile("src", BEDGRAPH, label="BEDGRAPH file"), intin("binsize", 1000),
    choice("stat", ["mean", "sum", "max", "min"])],
   ex={"src": BEDGRAPH, "binsize": 1000, "stat": "mean"}, up="iuc/bedtools",
   summary="Re-bin a bedGraph into fixed-size windows with a summary statistic.")
def bedgraph_binning(src, binsize=1000, stat="mean"):
    """bedGraph → binned table."""
    vals = io.parse_bedgraph(io.as_text(src))
    b = max(1, int(binsize))
    acc: dict[tuple[str, int], list[float]] = defaultdict(list)
    for c, s, e, v in vals:
        for bin_start in range(s // b * b, e, b):
            ov = min(e, bin_start + b) - max(s, bin_start)
            if ov > 0:
                acc[(c, bin_start)].append(v * ov / b)
    rows = []
    for (c, bs), vv in sorted(acc.items()):
        val = {"mean": stats.mean(vv), "sum": sum(vv), "max": max(vv), "min": min(vv)}[stat]
        rows.append({"chrom": c, "start": bs, "end": bs + b, stat: round(val, 5), "n_pieces": len(vv)})
    return table(rows, f"{len(vals)} bedGraph lines binned to {b} bp")


@T("bedtools_summarize", "Interval file summary", BED_SEC, "stats", [bed("src")], ex={"src": REGIONS_BED},
   up="iuc/bedtools", summary="Length statistics of a BED file (total, mean, N50, per-chrom).")
def bed_summary(src):
    """bedtools genomecov-ish summary."""
    st = genome.length_stats(_ivs(src))
    sm = genome.summarize_intervals(_ivs(src))
    out = {**st, **sm}
    return values_message({k: (round(v, 4) if isinstance(v, float) else v) for k, v in out.items()},
                          "interval statistics")


@T("bedtools_window_counts", "Count intervals in sliding windows", BED_SEC, "table",
   [bed("src"), intin("size", 1000), textbox("chrom", "")], ex={"src": REGIONS_BED, "size": 1000},
   up="iuc/bedtools", summary="Number of intervals per fixed window (windowed counts).")
def bed_window_counts(src, size=1000, chrom=""):
    """Windowed interval density."""
    ivs = _ivs(src)
    if chrom:
        ivs = [x for x in ivs if x.chrom == chrom]
    rows = genome.window_counts(ivs, size=int(size))
    return table(rows, f"{len(rows)} windows of {size} bp")


@T("bedtools_blacklist_filter", "Remove intervals overlapping a blacklist", BED_SEC, "bed",
   [bed("src"), bed("blacklist", ""), number("max_frac", 0.5)], ex={"src": REGIONS_BED,
                                                                    "blacklist": TARGETS_BED},
   up="iuc/bedtools", summary="Drop intervals that overlap the blacklist by more than a fraction.")
def bed_blacklist(src, blacklist="", max_frac=0.5):
    """Blacklist filtering as used in ChIP-seq pipelines."""
    bl = _ivs(blacklist) if blacklist else []
    out = genome.blacklist_filter(_ivs(src), bl, max_frac=float(max_frac))
    return _bed(out, f"{len(out)}/{len(_ivs(src))} intervals survive the blacklist")


@T("bedtools_intersect_with_counts_by_name", "Count overlaps per feature name", BED_SEC, "table",
   [bed("a"), bed("b"), boolean("reciprocal", False)], ex={"a": REGIONS_BED, "b": TARGETS_BED},
   up="iuc/bedtools", summary="Count how many B intervals hit each named A feature.")
def bed_count_by_name(a, b, reciprocal=False):
    """Name-keyed overlap counting."""
    ra, rb = _ivs(a), _ivs(b)
    hits = genome.intersect(ra, rb)
    c = Counter()
    for h in hits:
        c[h.name or f"{h.chrom}:{h.start}-{h.end}"] += 1
    rows = [{"feature": (x.name or f"{x.chrom}:{x.start + 1}-{x.end}"), "chrom": x.chrom,
             "start": x.start, "end": x.end, "overlaps": c.get(x.name or "", 0)} for x in ra]
    return table(rows, f"{sum(r['overlaps'] for r in rows)} overlaps counted")


@T("bedtools_gff_to_bed", "GFF to BED (operate on features)", "bed", "bed",
   [gff("src"), textbox("feature", "gene"), boolean("merge_exons", False)],
   ex={"src": ANNOT_GFF, "feature": "gene"}, up="iuc/gtf2bed", summary="Convert GFF features to BED.")
def gff_to_bed(src, feature="gene", merge_exons=False):
    """Annotation → intervals."""
    rows = [r for r in io.parse_gff(io.as_text(src))
            if not isinstance(getattr(r, "score", None), float) or not math.isnan(r.score)]
    if feature:
        rows = [r for r in rows if r.type == feature or
                (r.type in ("gene", "mRNA") and feature in ("gene", "mRNA", "transcript"))]
    ivs = _G.gff_to_intervals(rows)
    if merge_exons:
        ivs = genome.merge(ivs, distance=1)
    return _bed(genome.sort_ivs(ivs), f"{len(ivs)} intervals from {feature or 'all features'}")


@T("bedtools_genes_from_gff", "Get gene intervals from an annotation", "annotation", "bed",
   [gff("src"), textbox("biotype", ""), boolean("tss_only", False), intin("extend", 0)],
   ex={"src": ANNOT_GFF}, up="iuc/gtf_sort", summary="Extract gene bodies (or TSSs) from a GTF/GFF3.")
def genes_from_gff_tool(src, biotype="", tss_only=False, extend=0):
    """GFF gene records as BED, optionally reduced to TSSs."""
    rows = io.parse_gff(io.as_text(src))
    ivs = _G.genes_from_gff(rows)
    if biotype:
        ivs = [x for x in ivs if biotype in " ".join(map(str, x.extra))]
    if tss_only:
        ivs = [io.Interval(x.chrom, x.start if x.strand != "-" else max(x.start, x.end - 1),
                           x.end if x.strand != "-" else min(x.end, x.start + 1),
                           x.name, x.score, x.strand, x.extra) for x in ivs]
    if extend:
        ivs = genome.slop(ivs, None, bp=int(extend), direction="both")
    return _bed(genome.sort_ivs(ivs), f"{len(ivs)} gene regions")


@T("bedtools_tss_list", "TSS list from annotation", "annotation", "bed",
   [gff("src"), intin("upstream", 1000), intin("downstream", 1000)], ex={"src": ANNOT_GFF},
   up="iuc/homer", summary="Promoter windows around every transcript start site.")
def tss_promoters(src, upstream=1000, downstream=1000):
    """Promoter definition by strand-aware window around the TSS."""
    rows = io.parse_gff(io.as_text(src))
    tss = _G.tss_list(rows)
    out = []
    for t in tss:
        lo = t.start - (upstream if t.strand != "-" else downstream)
        hi = t.start + (downstream if t.strand != "-" else upstream)
        out.append(io.Interval(t.chrom, max(0, lo), max(1, hi), f"{t.name}_promoter", 0, t.strand))
    return _bed(genome.sort_ivs(out), f"{len(out)} promoters (+{upstream}/-{downstream})")


@T("bedtools_metagene_profile", "Metagene TSS profile of intervals", "chip_seq", "table",
   [bed("src"), gff("annotation", ANNOT_GFF), intin("bin_size", 100), intin("flank", 2000)],
   ex={"src": TARGETS_BED, "annotation": ANNOT_GFF, "bin_size": 200}, up="iuc/deepTools",
   summary="Density of A intervals around TSSs in fixed bins.")
def metagene_tss(src, annotation, bin_size=100, flank=2000):
    """Aggregate interval density relative to transcription start sites."""
    tss = genome.tss_list(io.parse_gff(io.as_text(annotation)))
    rows = _G.tss_distance_profile(_ivs(src), tss, bin_size=int(bin_size))
    res = table(rows, f"{len(rows)} bins around {len(tss)} TSSs")
    from chroma_titan.core import plot

    if rows:
        res["figure"] = plot.line({"interval density": [r.get("count", 0) for r in rows]},
                                 x=[r.get("bin", i) for i, r in enumerate(rows)],
                                 title="metagene profile", xlabel="distance bin", ylabel="overlaps")
    return res


# ===========================================================================
# liftover
# ===========================================================================
@T("liftover_chain_bed", "Lift Over BED with a chain file", "lift_over", "bed",
   [bed("src"), anyfile("chain", "", label="LiftOver chain file (Galaxy format)"),
    choice("mode", ["dotter", "ucsc_chain", "manual_offset"]), intin("offset", 0),
    boolean("fail_to_original", True)], ex={"src": REGIONS_BED, "mode": "manual_offset", "offset": 100},
   up="iuc/ucsc_liftover", summary="Move intervals between assemblies with a chain or offset.")
def liftover_bed(src, chain="", mode="manual_offset", offset=0, fail_to_original=True):
    """Apply a liftOver chain (chrom/targetStart/offset model) to a BED file."""
    chains = parse_chain(io.as_text(chain)) if chain else []
    ivs = _ivs(src)
    out, failed = [], []
    for iv in ivs:
        hit = match_chain(chains, iv)
        if hit:
            out.append(io.Interval(hit[0], hit[1], hit[2], iv.name, iv.score, iv.strand, iv.extra))
        elif mode == "manual_offset" or not chains:
            d = int(offset)
            out.append(io.Interval(iv.chrom, max(0, iv.start + d), max(1, iv.end + d), iv.name,
                                   iv.score, iv.strand, iv.extra))
        elif fail_to_original:
            out.append(iv)
            failed.append(iv)
    res = _bed(genome.sort_ivs(out), f"lifted {len(out)} intervals ({len(failed)} unchained)")
    return res


def parse_chain(text):
    """Minimal UCSC chain parser: tName tStrand tStart tEnd tSize qName qStrand qStart qEnd qSize."""
    out, block = [], None
    for ln in text.splitlines():
        p = ln.split()
        if ln.startswith("chain") and len(p) >= 13:
            block = {"score": int(p[1]), "t_name": p[2], "t_strand": p[3], "t_start": int(p[4]),
                     "t_end": int(p[5]), "q_name": p[7], "q_strand": p[8], "q_start": int(p[9]),
                     "q_end": int(p[10])}
            out.append(block)
    return out


def match_chain(chains, iv):
    for ch in chains:
        if ch["q_name"] == iv.chrom:
            shift = ch["t_start"] - ch["q_start"]
            return (ch["t_name"], iv.start + shift, iv.end + shift)
    return None


@T("liftover_pdb_chain", "LiftOver positions (VCF style)", "lift_over", "table",
   [vcf("src"), anyfile("chain", ""), intin("offset", 0), number("min_mapping_rate", 0.5)],
   ex={"src": VARIANTS_VCF, "offset": 50}, up="iuc/ucsc_liftover",
   summary="Shift variant coordinates between assemblies and report per-record status.")
def liftover_positions(src, chain="", offset=0, min_mapping_rate=0.5):
    """Variant coordinate conversion with a QC metric for unmapped sites."""
    chains = parse_chain(io.as_text(chain)) if chain else []
    v = io.parse_vcf(io.as_text(src))
    rows, lifted = [], 0
    for rec in v.records[:500]:
        hit = None
        for ch in chains:
            if ch["q_name"] == rec.chrom:
                hit = (ch["t_name"], rec.pos + ch["t_start"] - ch["q_start"])
                break
        if hit:
            rows.append({"chrom": hit[0], "position": hit[1], "ref": rec.ref, "alt": rec.alt,
                         "status": "lifted_chain"})
            lifted += 1
        else:
            rows.append({"chrom": rec.chrom, "position": rec.pos + int(offset), "ref": rec.ref,
                         "alt": rec.alt, "status": "offset_applied" if offset else "unchained"})
    rate = lifted / max(1, len(v.records))
    return table(rows, f"lifted {lifted}/{len(v.records)} records (mapping rate "
                       f"{100 * rate:.1f}%, required >= {100 * min_mapping_rate:.0f}%)")


@T("liftover_chain_summary", "Chain QC: coverage of the query assembly", "lift_over", "table",
   [anyfile("chain", "", label="chain file"), intin("top", 20)], ex={"chain": ""},
   up="iuc/ucsc_liftover", summary="Report chain blocks per chromosome pair.")
def chain_summary(chain="", top=20):
    """Summarise a chain file (no file → empty report)."""
    chains = parse_chain(io.as_text(chain)) if chain else []
    c: Counter = Counter()
    size: Counter = Counter()
    for ch in chains:
        key = f"{ch['q_name']}->{ch['t_name']}"
        c[key] += 1
        size[key] += ch["t_end"] - ch["t_start"]
    rows = [{"pair": k, "n_blocks": v, "covered_bases": size[k],
             "percent_of_target": round(100 * size[k] / max(1, sum(size.values())), 3)}
            for k, v in c.most_common(int(top))]
    return table(rows, f"{len(chains)} chain blocks")


# ===========================================================================
# more BED utilities
# ===========================================================================
@T("bedtools_firstline", "First line per interval (deduplicate coordinates)", BED_SEC, "bed",
   [bed("src")], ex={"src": REGIONS_BED}, up="iuc/bedtools",
   summary="bedtools firstline: drop duplicated coordinate rows.")
def bed_firstline(src):
    """Keep the first occurrence of each chrom:start-end key."""
    seen, out = set(), []
    for iv in _ivs(src):
        k = (iv.chrom, iv.start, iv.end)
        if k not in seen:
            seen.add(k)
            out.append(iv)
    return _bed(out, f"{len(out)} unique coordinates from {len(_ivs(src))} lines")


@T("bedtools_getfasta_split", "Split BED12 blocks into FASTA", "fetch_sequences_alignments", "fasta",
   [fa("genome", GENOME), bed("intervals", REGIONS_BED), intin("padding", 0)],
   ex={"genome": GENOME, "intervals": REGIONS_BED, "padding": 10}, up="iuc/bedtools",
   summary="bedtools getfasta -split style extraction with padding.")
def getfasta_split(genome, intervals, padding=0):
    """Pad each interval then extract sequence."""
    seqs, sizes = _read_genome(genome)
    ivs = [io.Interval(x.chrom, max(0, x.start - int(padding)), x.end + int(padding), x.name,
                       x.score, x.strand, x.extra) for x in _ivs(intervals)]
    out = genome_getfasta(seqs, ivs, strand="same", name_mode="bed")
    recs = [io.Seq(iv.name, (iv.extra or [""])[0]) for iv in out]
    return {"text": io.write_fasta(recs), "filename": "sequences.fasta",
            "message": f"{len(recs)} padded sequences ({padding} bp)"}


@T("bedtools_mask_from_genome", "Mask genome regions (soft/hard)", "bed", "fasta",
   [fa("genome", GENOME), bed("intervals", REGIONS_BED), choice("mode", ["soft", "hard"])],
   ex={"genome": GENOME, "intervals": REGIONS_BED, "mode": "soft"}, up="iuc/bedtools",
   summary="Return the genome with the given intervals lower-cased or N-masked.")
def mask_genome(genome, intervals, mode="soft"):
    """bedtools ``maskfasta`` behaviour."""
    seqs, sizes = genome_read(genome)
    ivs = _ivs(intervals)
    by_chrom: dict[str, list] = defaultdict(list)
    for iv in ivs:
        by_chrom[iv.chrom].append(iv)
    out = []
    for name, s in seqs.items():
        chars = list(s if isinstance(s, str) else "".join(s))
        for iv in by_chrom.get(name, ()):
            for i in range(iv.start, min(iv.end, len(chars))):
                chars[i] = "n" if mode == "soft" else "N"
        out.append(io.Seq(name, "".join(chars)))
    return {"text": io.write_fasta(out), "filename": "masked.fa",
            "message": f"masked {len(ivs)} intervals ({mode})"}


@T("bedtools_genomecov_binned", "Genome coverage in bins", "bed", "table",
   [anyfile("src", BEDGRAPH, label="bedGraph or coverage file"), intin("binsize", 500),
    intin("max_depth", 0)], ex={"src": BEDGRAPH, "binsize": 500}, up="iuc/bedtools",
   summary="bedtools genomecov -bga style depth in fixed bins.")
def genome_cov_bins(src, binsize=500, max_depth=0):
    """Depth profile binned to a fixed resolution."""
    vals = io.parse_bedgraph(io.as_text(src))
    b = max(1, int(binsize))
    acc: dict[tuple[str, int], float] = defaultdict(float)
    for c, s, e, v in vals:
        d = min(float(max_depth), v) if max_depth else v
        for bs in range(s // b * b, e, b):
            ov = min(e, bs + b) - max(s, bs)
            if ov > 0:
                acc[(c, bs)] += d * ov
    rows = [{"chrom": c, "start": bs, "end": bs + b, "mean_depth": round(acc[(c, bs)] / b, 5)}
            for (c, bs) in sorted(acc)]
    return table(rows, f"{len(rows)} bins of {b} bp")


@T("bedtools_zero_removal", "Remove zero-coverage regions", "bed", "bed",
   [anyfile("src", BEDGRAPH), number("threshold", 0.0), boolean("below", True)],
   ex={"src": BEDGRAPH, "threshold": 1.0}, up="iuc/bedtools",
   summary="Keep bedGraph segments above (or below) a value, as BED.")
def zero_removal(src, threshold=0.0, below=True):
    """Keep the low (or high) coverage segments and merge adjacent ones."""
    vals = io.parse_bedgraph(io.as_text(src))
    keep = [c for c in vals if (c[3] < threshold if below else c[3] >= threshold)]
    ivs = [io.Interval(c, s, e, ".", round(v, 4)) for c, s, e, v in keep]
    return _bed(genome.merge(ivs, distance=0),
               f"{len(ivs)} segments {'below' if below else 'above'} {threshold}")


@T("bedtools_merge_counts", "Count intervals per merged bin", BED_SEC, "table",
   [bed("src"), intin("size", 1000), choice("normalisation", ["none", "per_kb", "fraction_of_chrom"]),
    anyfile("genome", GENOME_TXT, label="genome file (for fraction)")],
   ex={"src": REGIONS_BED, "size": 1000}, up="iuc/featurecounts",
   summary="Bin a BED file and count (optionally normalise) intervals per bin.")
def merge_counts(src, size=1000, normalisation="none", genome=None):
    """Coverage-style counting in fixed bins."""
    ivs = _ivs(src)
    sizes = genome_file_sizes(genome) if genome else {}
    bins: dict[tuple[str, int], int] = defaultdict(int)
    for iv in ivs:
        for bs in range(iv.start // size * size, iv.end, size):
            bins[(iv.chrom, bs)] += 1
    rows = []
    for (c, bs), n in sorted(bins.items()):
        rec = {"chrom": c, "start": bs, "end": bs + size, "count": n}
        if normalisation == "per_kb":
            rec["count_per_kb"] = round(1000 * n / size, 4)
        elif normalisation == "fraction_of_chrom":
            rec["fraction"] = round(n / max(1, sizes.get(c, size) / size), 5)
        rows.append(rec)
    return table(rows, f"{len(rows)} bins from {len(ivs)} intervals")


@T("bedtools_restrict_to", "Restrict ranges to a reference list", BED_SEC, "bed",
   [bed("src"), bed("restrict", TARGETS_BED), boolean("invert", False)],
   ex={"src": REGIONS_BED, "restrict": TARGETS_BED}, up="iuc/bedtools",
   summary="Keep only intervals that overlap (or do not overlap) a restriction set.")
def restrict_to(src, restrict="", invert=False):
    """Whole-interval filter by overlap with a restriction file."""
    ivs = _ivs(src)
    hits = {id(x) for x in genome.intersect(ivs, _ivs(restrict))} if restrict else set()
    out = [x for x in ivs if (id(x) in hits) != invert]
    return _bed(out, f"{len(out)} intervals {'outside' if invert else 'inside'} the restriction set")


@T("bedtools_intersect_fractions", "Fraction of each interval overlapped", BED_SEC, "table",
   [bed("a"), bed("b")], ex={"a": REGIONS_BED, "b": TARGETS_BED}, up="iuc/bedtools",
   summary="Per-A overlap fraction and base counts (bedtools coverage -d style).")
def overlap_fractions(a, b):
    """Detailed overlap geometry per A interval."""
    ra, rb = _ivs(a), _ivs(b)
    rows = []
    for x in ra:
        ov = sum(genome.intersect_intervals(x, y) for y in rb if y.chrom == x.chrom)
        n = sum(1 for y in rb if y.chrom == x.chrom and genome.intersect_intervals(x, y) > 0)
        rows.append({"chrom": x.chrom, "start": x.start, "end": x.end, "name": x.name,
                     "interval_length": x.end - x.start, "overlap_bases": ov, "n_hits": n,
                     "overlap_fraction": round(ov / max(1, x.end - x.start), 5)})
    return table(rows, f"overlap fractions for {len(rows)} intervals")
