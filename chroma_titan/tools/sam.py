"""SAM/BAM alignment tools: samtools-style views and stats, Picard metrics,
deepTools-style coverage, RSeQC-style QC."""

from __future__ import annotations

import math
from collections import Counter, defaultdict

from chroma_titan.core import bam, genome, io, plot, seq, stats
from chroma_titan.tools._common import *  # noqa: F401,F403
from chroma_titan.tools._common import ALN_SAM, ANNOT_GFF, GENOME, REGIONS_BED, TARGETS_BED

SB = "sam_bam"


def _load(src):
    return bam.load(src)


def _bed(ivs, message=""):
    import math as _m

    clean = []
    for x in ivs:
        sc = x.score
        if isinstance(sc, float) and (_m.isnan(sc) or _m.isinf(sc)):
            sc = 0
        clean.append(io.Interval(x.chrom, int(x.start), int(x.end), x.name or ".", sc,
                                 x.strand or ".", list(x.extra or [])))
    return {"text": io.write_bed(clean, bed12=False), "filename": "output.bed",
            "message": message or f"{len(clean)} intervals"}


def _refs(src):
    seqs, sizes = genome.read_genome(src)
    return {k: "".join(v) for k, v in seqs.items()}, sizes


def _num(v):
    return round(v, 5) if isinstance(v, float) else v


# ===========================================================================
# samtools
# ===========================================================================
@T("samtools_flagstat", "Flagstat", SB, "table", [sam("src")], ex={"src": ALN_SAM},
   up="iuc/samtools", summary="samtools flagstat: counts by flag category.")
def samtools_flagstat(src):
    """Total, mapped, paired, proper-pair, duplicate and QC-fail counts."""
    hdr, alns = _load(src)
    st = bam.flagstat(hdr, alns)
    rows = [{"category": k, "count": v} for k, v in st.items() if isinstance(v, int)]
    return table(rows, f"flagstat over {len(alns)} records")


@T("samtools_idxstats", "Idxstats", SB, "table", [sam("src")], ex={"src": ALN_SAM}, up="iuc/samtools",
   summary="Per-reference mapped/unmapped counts and lengths.")
def samtools_idxstats(src):
    """``samtools idxstats`` table."""
    hdr, alns = _load(src)
    rows = bam.idxstats(hdr, alns)
    return table(rows, f"{len(rows)} references")


@T("samtools_depth", "Depth", SB, "table",
   [sam("src"), textbox("region", ""), intin("min_base_quality", 0), intin("min_mapq", 0),
    intin("max_depth", 0), boolean("all_sites", False), intin("head", 200)],
   ex={"src": ALN_SAM, "region": "chrV:1-1000"}, up="iuc/samtools",
   summary="Per-base depth (``samtools depth``) with quality filters.")
def samtools_depth(src, region="", min_base_quality=0, min_mapq=0, max_depth=0, all_sites=False,
                   head=200):
    """Depth profile of one region or the whole alignment."""
    hdr, alns = _load(src)
    d = bam.depth(hdr, alns, region=region or None, max_depth=int(max_depth),
                 quality_threshold=int(min_base_quality), all_reads=all_sites)
    rows = []
    for chrom, vals in d.items():
        for i, v in enumerate(vals):
            if v or all_sites:
                rows.append({"chrom": chrom, "position": i + 1, "depth": v})
    stats_row = {"chrom": "SUMMARY", "position": len(rows),
                "depth": round(stats.mean([r["depth"] for r in rows]), 4) if rows else 0.0}
    res = table(rows[: int(head)], f"{len(rows)} positions, mean depth "
                                   f"{stats_row['depth']}")
    res["summary"] = {"positions_nonzero": len(rows),
                      "max_depth": max((r["depth"] for r in rows), default=0),
                      "mean_depth": stats_row["depth"]}
    return res


@T("samtools_view", "View / filter alignments", SB, "text",
   [sam("src"), textbox("region", ""), intin("min_mapq", 0), intin("flag_filter", 0),
    intin("flag_require", 0), choice("read_type", ["all", "mapped", "unmapped", "paired",
                                                  "proper"]), boolean("drop_secondary", True),
    textbox("read_group", ""), intin("head", 50)],
   ex={"src": ALN_SAM, "read_type": "mapped", "min_mapq": 20, "head": 10}, up="iuc/samtools",
   summary="samtools view with region, flag, MAPQ and pairing filters.")
def samtools_view(src, region="", min_mapq=0, flag_filter=0, flag_require=0, read_type="all",
                  drop_secondary=True, read_group="", head=50):
    """Filtered SAM text output plus a count summary."""
    hdr, alns = _load(src)
    kw = {"region": region or None, "mapq": int(min_mapq), "flag_filter": int(flag_filter),
          "flag_require": int(flag_require), "drop_secondary": drop_secondary,
          "read_group": read_group or None}
    if read_type == "unmapped":
        kw["unmapped"] = True
    elif read_type == "mapped":
        kw["unmapped"] = False
    elif read_type in ("paired", "proper"):
        kw["paired_only"] = True
        if read_type == "proper":
            kw["flag_require"] = int(flag_require) | 2
    keep = bam.view_filter(alns, **kw)
    lines = hdr + [a.as_line() for a in keep[: int(head)]]
    return {"text": "\n".join(lines) + "\n", "filename": "filtered.sam",
            "message": f"{len(keep)}/{len(alns)} alignments pass (showing {min(len(keep), int(head))})"}


@T("samtools_view_count", "Count reads matching a filter", SB, "table",
   [sam("src"), choice("group_by", ["reference", "read1", "flag", "mapq", "cigar", "pair", "name"])],
   ex={"src": ALN_SAM, "group_by": "reference"}, up="iuc/samtools",
   summary="``samtools view -c`` grouped by a chosen key.")
def samtools_count(src, group_by="reference"):
    """Alignment counts by reference/flag/MAPQ."""
    hdr, alns = _load(src)
    rows = bam.count_alns(alns, key=group_by)
    return table(rows, f"{len(rows)} groups ({group_by})")


@T("samtools_sort", "Sort SAM by coordinate or name", SB, "text",
   [sam("src"), choice("by", ["coordinate", "queryname", "unsorted"])], ex={"src": ALN_SAM},
   up="iuc/samtools", summary="Reorder the alignment records.")
def samtools_sort_tool(src, by="coordinate"):
    """``samtools sort``."""
    hdr, alns = _load(src)
    out = bam.sort_sam(alns, by=by)
    return {"text": "\n".join(hdr + [a.as_line() for a in out]) + "\n", "filename": "sorted.sam",
            "message": f"sorted {len(out)} alignments by {by}"}


@T("samtools_merge", "Merge several SAM files", SB, "text",
   [multi("files", [ALN_SAM], [ALN_SAM]), boolean("assume_sorted", False),
    boolean("add_read_group", False)], ex={"files": [ALN_SAM, ALN_SAM]}, up="iuc/samtools",
   summary="``samtools merge``: concatenate headers and alignments.")
def samtools_merge_files(files=None, assume_sorted=False, add_read_group=False):
    """Merge alignments from several files."""
    parts = [io.as_text(f) for f in (files or [])]
    hdr, alns = bam.merge_sams(parts)
    if not assume_sorted:
        alns = bam.sort_sam(alns, by="coordinate")
    lines = hdr + [a.as_line() for a in alns]
    return {"text": "\n".join(lines) + "\n", "filename": "merged.sam",
            "message": f"merged {len(parts)} files -> {len(alns)} alignments"}


@T("samtools_markdup", "Mark duplicates", SB, "table",
   [sam("src"), boolean("flow_mode", False), boolean("output_text", False), intin("head", 30)],
   ex={"src": ALN_SAM}, up="iuc/samtools",
   summary="Identify optical-agnostic duplicates by outer coordinates.")
def samtools_markdup(src, flow_mode=False, output_text=False, head=30):
    """``samtools markdup``-style duplicate marking with a stats table."""
    hdr, alns = _load(src)
    kept, st = bam.mark_duplicates(alns, flow=flow_mode)
    rows = [_num(k) and {"metric": k, "value": _num(v)} for k, v in st.items()]
    res = table(rows, f"{st.get('dupes', 0)} duplicates marked out of {st.get('unq_pairs', 0)} pairs")
    if output_text:
        res["text"] = "\n".join(hdr + [a.as_line() for a in kept[: int(head)]])
        res["filename"] = "marked.sam"
    return res


@T("samtools_rmdup_position", "Remove duplicate reads by position", SB, "text",
   [sam("src"), boolean("keep_highest_mapq", True), boolean("paired", False)], ex={"src": ALN_SAM},
   up="iuc/samtools", summary="Drop PCR duplicates sharing the same alignment start.")
def remove_duplicates(src, keep_highest_mapq=True, paired=False):
    """Position-based duplicate removal."""
    hdr, alns = _load(src)
    best: dict[tuple, object] = {}
    for a in alns:
        key = (a.qname.split("/")[0] if paired else "", a.rname, a.pos, a.flag & 16)
        cur = best.get(key)
        if cur is None or (keep_highest_mapq and a.mapq > cur.mapq):
            best[key] = a
    out = list(best.values())
    return {"text": "\n".join(hdr + [a.as_line() for a in bam.sort_sam(out)]) + "\n",
            "filename": "dedup.sam",
            "message": f"{len(alns)} -> {len(out)} alignments ({len(alns) - len(out)} duplicates removed)"}


@T("samtools_stats", "Samtools stats", SB, "stats", [sam("src"), fa("genome", "")],
   ex={"src": ALN_SAM}, up="iuc/samtools", summary="SNPs/indels, insert size, read lengths, error rate.")
def samtools_stats(src, genome=""):
    """Aggregated alignment statistics."""
    hdr, alns = _load(src)
    flag = bam.flagstat(hdr, alns)
    ins = bam.insert_size_stats(alns)
    err = bam.error_rate_from_qual(alns)
    rlen = bam.read_length_stats(alns)
    cig = bam.cigar_summary(alns)
    out: dict = {}
    out.update({f"flag_{k}": v for k, v in flag.items() if isinstance(v, (int, float))})
    out.update({f"insert_{k}": _num(v) for k, v in ins.items()})
    out.update({f"quality_{k}": _num(v) for k, v in err.items()})
    out.update({f"readlen_{k}": _num(v) for k, v in rlen.items()})
    out.update({"cigar_" + str(k): v for k, v in list(cig.items())[:6] if isinstance(v, (int, float))})
    return values_message(out, f"stats for {len(alns)} alignments")


@T("samtools_mismatch_profile", "Mismatch rate per cycle", SB, "table",
   [sam("src"), intin("read_len", 0), intin("min_mapq", 0)], ex={"src": ALN_SAM}, up="iuc/samtools",
   summary="Error rate at each read position (mismatch profile).")
def mismatch_profile_tool(src, read_len=0, min_mapq=0):
    """Per-cycle mismatch/indel rates."""
    hdr, alns = _load(src)
    rows = bam.mismatch_profile([a for a in alns if a.mapq >= min_mapq], read_len=int(read_len))
    res = table(rows[:200], f"{len(rows)} cycles profiled")
    if rows:
        res["figure"] = plot.line({"mismatch rate": [r.get("mismatch_rate", 0) for r in rows[:200]]},
                                  x=list(range(1, min(200, len(rows)) + 1)),
                                  title="mismatches per cycle", xlabel="cycle", ylabel="rate")
    return res


@T("samtools_nm_distribution", "Edit-distance (NM) distribution", SB, "table",
   [sam("src")], ex={"src": ALN_SAM}, up="iuc/samtools", summary="Number of alignments per NM value.")
def nm_distribution(src):
    """``samtools stats`` NM histogram."""
    hdr, alns = _load(src)
    d = bam.nm_distribution(alns)
    return table([{"NM": k, "alignments": v} for k, v in sorted(d.items())],
                 f"NM distribution over {sum(d.values())} alignments")


@T("samtools_mapq_hist", "Mapping quality histogram", SB, "table", [sam("src"), intin("bins", 10)],
   ex={"src": ALN_SAM}, up="iuc/samtools", summary="Counts per MAPQ bin.")
def mapq_histogram(src, bins=10):
    """MAPQ distribution table with a bar plot."""
    hdr, alns = _load(src)
    d = bam.mapq_histogram(alns)
    tot = sum(d.values()) or 1
    rows = [{"mapq": k, "alignments": v, "percent": round(100 * v / tot, 3)} for k, v in sorted(d.items())]
    res = table(rows, f"MAPQ histogram ({len(rows)} distinct values)")
    res["figure"] = plot.bar([str(r["mapq"]) for r in rows], [r["alignments"] for r in rows],
                             title="mapping quality", xlabel="MAPQ", ylabel="alignments")
    return res


@T("samtools_tag_summary", "Summarise an alignment tag", SB, "table",
   [sam("src"), textbox("tag", "NM")], ex={"src": ALN_SAM, "tag": "NM"}, up="iuc/samtools",
   summary="Value counts for a SAM tag (RG, NM, AS, …).")
def tag_summary_tool(src, tag="NM"):
    """``samtools view -c --tag-filter`` style summary."""
    hdr, alns = _load(src)
    rows = bam.tag_summary(alns, tag)
    return table(rows, f"tag {tag}: {len(rows)} distinct values")


@T("samtools_addreplacerg", "Add or replace read groups", SB, "text",
   [sam("src"), textbox("id", "run1"), textbox("sample", "SAMPLE_X"), textbox("library", "lib1"),
    textbox("platform", "ILLUMINA")], ex={"src": ALN_SAM, "sample": "my_sample"}, up="iuc/samtools",
   summary="Stamp every alignment with a new @RG line and RG tag.")
def add_replace_rg(src, id="run1", sample="SAMPLE_X", library="lib1", platform="ILLUMINA"):
    """Read-group rewriting."""
    hdr, alns = _load(src)
    new_hdr = [h for h in hdr if not h.startswith("@RG")]
    new_hdr.insert(1, f"@RG\tID:{id}\tSM:{sample}\tLB:{library}\tPL:{platform}")
    out = []
    for a in alns:
        tags = dict(a.tags or {})
        tags["RG"] = id
        out.append(io.Aln(a.qname, a.flag, a.rname, a.pos, a.mapq, a.cigar, a.rnext, a.pnext,
                         a.tlen, a.seq, a.qual, tags))
    return {"text": "\n".join(new_hdr + [a.as_line() for a in out]) + "\n", "filename": "rg.sam",
            "message": f"added read group {id} ({sample}) to {len(out)} alignments"}


@T("samtools_pileup", "Mpileup text", SB, "text",
   [sam("src"), fa("genome", GENOME), textbox("region", ""), intin("min_base_quality", 13),
    intin("min_mapq", 0), intin("max_rows", 60)], ex={"src": ALN_SAM, "region": "chrV:1-2000"},
   summary="Per-base pileup strings (``samtools mpileup``).")
def samtools_pileup(src, genome=GENOME, region="", min_base_quality=13, min_mapq=0, max_rows=60):
    """Reference-aware pileup of a region."""
    seqs, sizes = _refs(genome) if genome else ({}, {})
    hdr, alns = _load(src)
    rows = bam.pileup(hdr, alns, min_base_quality=int(min_base_quality), min_mapq=int(min_mapq),
                     region=region or None, ref=seqs or None)
    lines = []
    for r in rows[: int(max_rows)]:
        counts = r.get("counts") or {}
        bases = "".join(f"{b}{n}" for b, n in sorted(counts.items())) or "."
        lines.append(f"{r['chrom']}\t{r['position']}\t{r.get('ref_base', 'N')}\t{r['depth']}"
                     f"\t{bases}\t{r.get('mean_base_quality', 0):.1f}")
    return {"text": "\n".join(lines) + "\n", "filename": "pileup.txt",
            "message": f"{len(rows)} pileup positions (showing {min(len(rows), int(max_rows))})"}


@T("samtools_call_variants", "Call variants from a pileup", SB, "table",
   [sam("src"), fa("genome", GENOME), intin("min_depth", 5), number("min_vaf", 0.2),
    number("min_qual", 20.0), number("het_fraction", 0.75), textbox("region", "")],
   ex={"src": ALN_SAM, "genome": GENOME, "min_depth": 5}, up="iuc/samtools",
   summary="Simple pileup genotype caller producing a variant table.")
def call_variants(src, genome=GENOME, min_depth=5, min_vaf=0.2, min_qual=20.0, het_fraction=0.75,
                  region=""):
    """Pileup → SNP calls with depth, VAF and QUAL."""
    seqs, sizes = _refs(genome)
    hdr, alns = _load(src)
    rows = bam.pileup(hdr, alns, region=region or None, ref=seqs or None)
    calls = bam.call_variants_from_pileup(rows, min_depth=int(min_depth), min_vaf=min_vaf,
                                         min_qual=min_qual, het_fraction=het_fraction)
    return table(calls, f"{len(calls)} variant calls from {len(rows)} sites")


@T("samtools_vcf_from_pileup", "Pileup to VCF", SB, "text",
   [sam("src"), fa("genome", GENOME), textbox("sample", "SAMPLE1"), intin("min_depth", 5),
    number("min_vaf", 0.25)], ex={"src": ALN_SAM, "genome": GENOME}, up="iuc/samtools",
   summary="Write a VCF from pileup-based calls.")
def vcf_from_pileup(src, genome=GENOME, sample="SAMPLE1", min_depth=5, min_vaf=0.25):
    """``bcftools mpileup | call`` in one step."""
    seqs, sizes = _refs(genome)
    hdr, alns = _load(src)
    rows = bam.pileup(hdr, alns, ref=seqs or None)
    calls = bam.call_variants_from_pileup(rows, min_depth=int(min_depth), min_vaf=min_vaf)
    body = bam.vcf_from_calls(calls, samples=[sample])
    return {"text": body, "filename": "calls.vcf",
            "message": f"{len(calls)} records written to VCF"}


@T("samtools_consensus", "Generate consensus from alignments", SB, "fasta",
   [sam("src"), fa("genome", GENOME), textbox("caller", "samtools"), intin("min_depth", 3),
    boolean("report_variants", True)], ex={"src": ALN_SAM, "genome": GENOME}, up="iuc/samtools",
   summary="Build a consensus FASTA of the covered regions (``samtools consensus``).")
def samtools_consensus(src, genome=GENOME, caller="samtools", min_depth=3, report_variants=True):
    """Consensus sequence + variant list."""
    seqs, sizes = _refs(genome)
    hdr, alns = _load(src)
    rows = bam.pileup(hdr, alns, ref=seqs or None)
    cons, variants = bam.consensus_from_pileup(rows, seqs, min_depth=int(min_depth), caller=caller)
    recs = [io.Seq(k, v) for k, v in cons.items()]
    res = {"text": io.write_fasta(recs), "filename": "consensus.fa",
           "message": f"consensus over {len(recs)} references, {len(variants)} differences"}
    if report_variants:
        res["table"] = io.to_tsv(__import__("pandas").DataFrame(variants[:200])) if variants else ""
        res["stats"] = {"variants": len(variants)}
    return res


@T("samtools_fastq", "Extract reads as FASTQ", SB, "text",
   [sam("src"), boolean("mapped_only", True), choice("which", ["all", "read1", "read2", "unmapped"])],
   ex={"src": ALN_SAM, "mapped_only": True}, up="iuc/samtools",
   summary="Convert alignments back into a FASTQ file.")
def sam_to_fastq(src, mapped_only=True, which="all"):
    """``samtools fastq``."""
    hdr, alns = _load(src)
    if which == "unmapped":
        alns = [a for a in alns if not a.mapped]
    elif which == "read1":
        alns = [a for a in alns if a.read1]
    elif which == "read2":
        alns = [a for a in alns if a.read2]
    elif mapped_only:
        alns = [a for a in alns if a.mapped]
    return {"text": bam.reads_to_fastq(alns), "filename": "reads.fastq",
            "message": f"{len(alns)} alignments -> FASTQ"}


@T("samtools_downsample", "Downsample alignments", SB, "text",
   [sam("src"), number("fraction", 0.5), intin("n", 0), intin("seed", 42)],
   ex={"src": ALN_SAM, "fraction": 0.5}, up="iuc/samtools", summary="Randomly keep a fraction of reads.")
def downsample_sam(src, fraction=0.5, n=0, seed=42):
    """``samtools view -s``."""
    hdr, alns = _load(src)
    out = bam.downsample_alns(alns, fraction=fraction, n=int(n), seed=int(seed))
    return {"text": "\n".join(hdr + [a.as_line() for a in out]) + "\n", "filename": "sub.sam",
            "message": f"{len(out)}/{len(alns)} alignments kept"}


@T("samtools_bedgraph", "Genome coverage to bedGraph", SB, "bedgraph",
   [sam("src"), number("scale", 1.0), boolean("fragment", True), intin("bin_size", 0),
    number("min_depth", 0.0)], ex={"src": ALN_SAM}, up="iuc/samtools",
   summary="Depth per base interval as bedGraph (``bedtools genomecov -bg``).")
def sam_to_bedgraph(src, scale=1.0, fragment=True, bin_size=0, min_depth=0.0):
    """BedGraph of the coverage profile."""
    hdr, alns = _load(src)
    vals = bam.bedgraph_from_sam(hdr, alns, scale=float(scale), fragment=fragment)
    if int(bin_size):
        b = int(bin_size)
        acc: dict[tuple[str, int], float] = defaultdict(float)
        for c, s, e, v in vals:
            for bs in range(s // b * b, e, b):
                ov = min(e, bs + b) - max(s, bs)
                if ov > 0:
                    acc[(c, bs)] += v * ov / b
        vals = [(c, bs, bs + b, round(vv, 5)) for (c, bs), vv in sorted(acc.items())]
    if min_depth:
        vals = [x for x in vals if x[3] >= min_depth]
    return {"text": io.write_bedgraph(vals), "filename": "coverage.bedgraph",
            "message": f"{len(vals)} bedGraph segments"}


@T("samtools_coverage_histogram", "Coverage histogram", SB, "table",
   [sam("src"), intin("binsize", 500), boolean("per_base", True)], ex={"src": ALN_SAM},
   up="iuc/samtools", summary="Mean depth in fixed-size bins (``plotCoverage`` data).")
def coverage_histogram(src, binsize=500, per_base=True):
    """Binned coverage table with a line plot."""
    hdr, alns = _load(src)
    rows = bam.coverage_histogram(hdr, alns, binsize=int(binsize), per_base=per_base)
    res = table(rows, f"{len(rows)} bins of {binsize} bp")
    if rows:
        res["figure"] = plot.line({"depth": [r.get("mean", 0.0) for r in rows]},
                                 x=list(range(len(rows))), title="coverage per bin",
                                 xlabel=f"bin ({binsize} bp)", ylabel="depth")
    return res


@T("samtools_insert_size_hist", "Insert size distribution", SB, "table",
   [sam("src"), intin("bins", 25)], ex={"src": ALN_SAM}, up="iuc/samtools",
   summary="Histogram of fragment lengths with mean/SD/3xSD cut-off.")
def insert_size_histogram(src, bins=25):
    """``CollectInsertSizeMetrics``-style histogram."""
    hdr, alns = _load(src)
    sizes = bam.insert_sizes(alns)
    if not sizes:
        return table([], "no properly paired reads")
    lo, hi = min(sizes), max(sizes)
    nb = max(2, int(bins))
    w = ((hi - lo) / nb) or 1
    c = Counter(min(nb - 1, int((s - lo) / w)) for s in sizes)
    rows = [{"bin": f"{lo + i * w:.0f}-{lo + (i + 1) * w:.0f}", "pairs": c.get(i, 0)}
            for i in range(nb)]
    st = bam.insert_size_stats(alns)
    res = table(rows, f"insert size mean {st.get('mean')} sd {st.get('sd')}")
    res["figure"] = plot.bar([r["bin"] for r in rows], [r["pairs"] for r in rows],
                             title="insert sizes", xlabel="bp", ylabel="pairs")
    return res


@T("samtools_tlen_scatter", "Fragment length vs position", SB, "table", [sam("src"), intin("max_rows", 400)],
   ex={"src": ALN_SAM}, up="iuc/samtools", summary="Per-pair TLEN track (for library QC).")
def tlen_scatter(src, max_rows=400):
    """Scatter data of template length along the genome."""
    hdr, alns = _load(src)
    rows = bam.tlen_scatter(alns)
    res = table(rows[: int(max_rows)], f"{len(rows)} pairs with insert length")
    if rows:
        res["figure"] = plot.scatter([r["pos"] for r in rows[: int(max_rows)]],
                                    [r["tlen"] for r in rows[: int(max_rows)]],
                                    title="insert size vs position", xlabel="position", ylabel="TLEN")
    return res


@T("samtools_softclip_stats", "Soft-clip statistics", SB, "table", [sam("src"), intin("min_clip", 5)],
   ex={"src": ALN_SAM}, up="iuc/samtools", summary="Clipped bases per read end (breakpoint signal).")
def softclip_summary(src, min_clip=5):
    """Soft-clipping profile by position and side."""
    hdr, alns = _load(src)
    rows = [r for r in bam.softclip_stats(alns) if r.get("mean_clip", 0) >= 0]
    return table(rows, f"{len(rows)} positions with soft clips >= {min_clip}")


@T("samtools_cigar_summary", "CIGAR summary", SB, "stats", [sam("src")], ex={"src": ALN_SAM},
   up="iuc/samtools", summary="Totals of M/I/D/S/H/N/P operations in the alignment.")
def cigar_summary_tool(src):
    """Count of each CIGAR operation across the file."""
    hdr, alns = _load(src)
    st = bam.cigar_summary(alns)
    return values_message({k: _num(v) for k, v in st.items()},
                         f"{len(alns)} alignments, operations counted")


@T("samtools_junctions", "Exon junctions from alignments", "rseqc", "table",
   [sam("src"), intin("min_count", 1)], ex={"src": ALN_SAM}, up="iuc/hisat2",
   summary="Spliced junctions (N in CIGAR) with support counts.")
def junctions(src, min_count=1):
    """Junction table (``junctions.bed`` in RSeQC terms)."""
    hdr, alns = _load(src)
    rows = bam.junctions_from_sam(alns, min_count=int(min_count))
    return table(rows, f"{len(rows)} junctions with >= {min_count} reads")


@T("samtools_metagene", "Metagene profile of alignments", "rseqc", "table",
   [sam("src"), gff("annotation", ANNOT_GFF), intin("profile_bins", 40), intin("up", 500),
    intin("down", 500)], ex={"src": ALN_SAM, "annotation": ANNOT_GFF}, up="iuc/RSeQC",
   summary="Average read density across scaled gene bodies.")
def metagene(src, annotation, profile_bins=40, up=500, down=500):
    """Density along genes (``metagene.py``)."""
    hdr, alns = _load(src)
    genes = genome.genes_from_gff(io.parse_gff(io.as_text(annotation)))
    prof = bam.metagene_profile(alns, genes, profile_bins=int(profile_bins), up=int(up), down=int(down))
    rows = [{"bin": i + 1, "density": round(v, 5)} for i, v in enumerate(prof)]
    res = table(rows, f"metagene over {len(genes)} genes ({len(prof)} bins)")
    res["figure"] = plot.line({"density": prof}, title="metagene profile", xlabel="gene position",
                             ylabel="reads per bin")
    return res


@T("samtools_read_distribution", "Read distribution relative to features", "rseqc", "table",
   [sam("src"), gff("annotation", ANNOT_GFF), intin("bin_size", 100)], ex={"src": ALN_SAM},
   up="iuc/RSeQC", summary="RSeQC read_distribution: reads in CDS/5'UTR/3'UTR/introns/intergenic.")
def read_distribution(src, annotation, bin_size=100):
    """Classification of alignments against annotated regions."""
    hdr, alns = _load(src)
    feats = genome.genes_from_gff(io.parse_gff(io.as_text(annotation)))
    st = bam.rseqc_read_distribution(alns, feats, bin_size=int(bin_size))
    tot = sum(v for v in st.values() if isinstance(v, (int, float))) or 1
    rows = [{"category": k, "reads": v, "percent": round(100 * v / tot, 3)}
            for k, v in st.items() if isinstance(v, (int, float))]
    return table(rows, f"read distribution over {len(feats)} features")


@T("samtools_feature_counts", "Assign reads to features", "rseqc", "table",
   [sam("src"), gff("annotation", ANNOT_GFF), choice("mode", ["intersection_nonempty", "union",
                                                             "strict"]), boolean("count_read_overlaps", False),
    boolean("nonunique", True)], ex={"src": ALN_SAM, "annotation": ANNOT_GFF}, up="iuc/featureCounts",
   summary="HTSeq-style counting of alignments per feature.")
def feature_counts(src, annotation, mode="intersection_nonempty", count_read_overlaps=False,
                   nonunique=True):
    """Read counts per gene."""
    hdr, alns = _load(src)
    feats = io.parse_gff(io.as_text(annotation))
    counts = bam.assign_to_features(alns, feats, mode=mode, count_read_overlaps=count_read_overlaps,
                                   nonunique=nonunique)
    rows = [{"feature": k, "counts": v} for k, v in sorted(counts.items(), key=lambda kv: -kv[1])]
    return table(rows, f"{len(rows)} features counted, {sum(counts.values())} reads assigned")


@T("samtools_inner_distance", "Inner distance (fragment size per bin)", "rseqc", "table",
   [sam("src"), bed("regions", REGIONS_BED), intin("bin_size", 500)], ex={"src": ALN_SAM},
   up="iuc/inner_distance", summary="Mean insert size inside each feature (library QC).")
def inner_distance(src, regions, bin_size=500):
    """Fragment lengths restricted to annotated intervals."""
    hdr, alns = _load(src)
    ivs = io.parse_bed(io.as_text(regions))
    sizes = bam.insert_sizes(alns)
    by_chrom: dict[str, list] = defaultdict(list)
    for a in alns:
        if a.tlen:
            by_chrom[a.rname].append((a.pos, abs(a.tlen)))
    rows = []
    for iv in ivs:
        vals = [t for p, t in by_chrom.get(iv.chrom, ()) if iv.start <= p <= iv.end]
        rows.append({"chrom": iv.chrom, "start": iv.start, "end": iv.end, "feature": iv.name,
                     "pairs": len(vals), "mean_inner_distance": round(stats.mean(vals), 2) if vals else "",
                     "sd": round(stats.stdev(vals), 2) if len(vals) > 1 else ""})
    return table(rows, f"inner distance for {len(rows)} intervals")


@T("samtools_gc_bias", "GC bias of the library", "picard", "table",
   [sam("src"), fa("genome", GENOME), intin("bins", 10)], ex={"src": ALN_SAM, "genome": GENOME},
   up="iuc/picard", summary="CollectGcBiasMetrics-style observed vs expected GC distribution.")
def gc_bias(src, genome=GENOME, bins=10):
    """GC-content bias table."""
    seqs, sizes = _refs(genome)
    hdr, alns = _load(src)
    rows = bam.gc_bias(alns, seqs)
    tot = sum(r.get("reads", 0) for r in rows) or 1
    for r in rows:
        r["percent"] = round(100 * r.get("reads", 0) / tot, 3)
    res = table(rows, f"GC bias over {len(rows)} bins")
    if rows:
        res["figure"] = plot.bar([str(r.get("gc_bin", i)) for i, r in enumerate(rows)],
                                [r.get("reads", 0) for r in rows], title="GC bias",
                                xlabel="GC bin", ylabel="reads")
    return res


@T("samtools_alignment_summary", "Picard AlignmentSummaryMetrics", "picard", "table",
   [sam("src"), intin("read_length", 100)], ex={"src": ALN_SAM}, up="iuc/picard",
   summary="Total/proper pairs, mismatch and indel rates per read.")
def alignment_summary(src, read_length=100):
    """``CollectAlignmentSummaryMetrics`` style rows by the first read of each pair."""
    hdr, alns = _load(src)
    prim = bam.primary(alns, mapped_only=False)
    total = len(prim)
    mapped = sum(1 for a in prim if a.mapped)
    proper = sum(1 for a in prim if a.proper_pair)
    nm = sum(int((a.tags or {}).get("NM", 0) or 0) for a in prim)
    bases = sum(len(a.seq) for a in prim if a.mapped)
    clipped = sum(1 for a in prim if "S" in (a.cigar or ""))
    rows = [{"category": "TOTAL", "reads": total, "bases": sum(len(a.seq) for a in prim),
            "mapped": mapped, "proper_pairs": proper, "mismatches": nm,
            "error_rate": round(nm / max(1, bases), 6), "clipped": clipped,
            "mean_length": round(stats.mean([len(a.seq) for a in prim]), 2)},
           {"category": "FIRST", "reads": sum(1 for a in prim if a.read1),
            "bases": sum(len(a.seq) for a in prim if a.read1),
            "mapped": sum(1 for a in prim if a.read1 and a.mapped),
            "proper_pairs": sum(1 for a in prim if a.read1 and a.proper_pair),
            "mismatches": 0, "error_rate": 0.0,
            "clipped": sum(1 for a in prim if a.read1 and "S" in (a.cigar or "")),
            "mean_length": round(stats.mean([len(a.seq) for a in prim if a.read1] or [0]), 2)}]
    return table(rows, f"alignment summary: {mapped}/{total} mapped")


@T("samtools_library_complexity", "Picard EstimateLibraryComplexity", "picard", "table",
   [sam("src"), intin("pairs", 100000)], ex={"src": ALN_SAM}, up="iuc/picard",
   summary="Duplicate-based library size estimate (L = N ln(N/(N-M))).")
def library_complexity(src, pairs=100000):
    """MarkDuplicates-style complexity estimate."""
    hdr, alns = _load(src)
    _, st = bam.mark_duplicates(alns)
    n = int(st.get("unq_pairs") or st.get("pairs") or len(alns))
    m = int(st.get("dupes") or 0)
    lib = 0.0
    if n > m > 0:
        lib = n * math.log(n / (n - m))
    rows = [{"metric": "read_pairs_examined", "value": n},
            {"metric": "duplicate_pairs", "value": m},
            {"metric": "duplicate_rate_percent", "value": round(100 * m / max(1, n), 3)},
            {"metric": "estimated_library_size", "value": int(round(lib))},
            {"metric": "pairs_for_estimate", "value": int(pairs)}]
    return table(rows, f"estimated library size {int(round(lib))} pairs")


@T("samtools_validate", "Picard ValidateSamFile style checks", "picard", "table",
   [sam("src"), boolean("strict", False), intin("max_errors", 20)], ex={"src": ALN_SAM},
   up="iuc/picard", summary="Structural validation of SAM records (CIGAR vs length, flags).")
def validate_sam(src, strict=False, max_errors=20):
    """Consistency checks on flags, CIGARs and coordinates."""
    hdr, alns = _load(src)
    rows = []
    for a in alns:
        problems = []
        ref_len = io.ref_length(a.cigar) if a.cigar else 0
        q_len = io.query_length(a.cigar) if hasattr(io, "query_length") else len(a.seq or "")
        if a.mapped and not a.cigar:
            problems.append("MISSING_CIGAR")
        if a.mapped and ref_len and a.pos + ref_len - 1 > 10 ** 9:
            problems.append("COORD_OUT_OF_RANGE")
        if a.seq and q_len and len(a.seq) != q_len:
            problems.append("CIGAR_SEQ_LENGTH_MISMATCH")
        if a.qual and len(a.qual) != len(a.seq):
            problems.append("QUAL_LENGTH_MISMATCH")
        if a.flag & 8 and a.seq not in ("*", ""):
            problems.append("UNMAPPED_WITH_SEQUENCE")
        for p in problems[: int(max_errors)]:
            rows.append({"read": a.qname, "problem": p, "severity": "ERROR" if strict else "WARNING",
                         "flag": a.flag, "cigar": a.cigar})
    return table(rows, f"{len(rows)} issues in {len(alns)} records")


@T("samtools_fingerprint", "Reads fingerprint (cross-sample check)", "picard", "table",
   [sam("src"), fa("genome", GENOME), intin("binsize", 1000)], ex={"src": ALN_SAM, "genome": GENOME},
   up="iuc/picard", summary="CheckBias/fingerprint style binned read-count profile.")
def fingerprint_tool(src, genome=GENOME, binsize=1000):
    """Correlation-friendly coverage fingerprint."""
    seqs, sizes = _refs(genome)
    hdr, alns = _load(src)
    st = bam.fingerprint(alns, seqs, binsize=int(binsize))
    rows = [{"metric": k, "value": _num(v)} for k, v in st.items()]
    return table(rows, "fingerprint summary")


@T("deeptools_bamcoverage", "bamCoverage (normalise depth)", "deeptools", "bedgraph",
   [sam("src"), choice("normalisation", ["none", "RPKM", "FPKM", "BAM", "CPS"]),
    number("scale_factor", 1.0), intin("bin_size", 0), intin("effective_genome_length", 1),
    number("extend_reads", 0.0), fa("genome", "")], ex={"src": ALN_SAM, "normalisation": "RPKM"},
   up="iuc/deepTools", summary="deepTools bamCoverage: normalised bedGraph of read depth.")
def bam_coverage(src, normalisation="none", scale_factor=1.0, bin_size=0,
                 effective_genome_length=1, extend_reads=0.0, genome=""):
    """Normalise coverage by depth or read count."""
    hdr, alns = _load(src)
    vals = bam.bedgraph_from_sam(hdr, alns, scale=1.0, fragment=True)
    n_reads = sum(1 for a in alns if a.mapped) or 1
    total_bases = sum(len(a.seq) for a in alns) or 1
    glen = max(1, int(effective_genome_length))
    factor = 1.0
    if normalisation == "RPKM":
        factor = 1e9 / (n_reads * 200)
    elif normalisation == "FPKM":
        factor = 1e9 / (n_reads * 200 / 2)
    elif normalisation == "BAM":
        factor = float(scale_factor)
    elif normalisation == "CPS":
        factor = 1e6 / n_reads
    else:
        factor = float(scale_factor)
    out = [(c, s, e, round(v * factor, 5)) for c, s, e, v in vals]
    if int(bin_size):
        b = int(bin_size)
        acc: dict[tuple[str, int], float] = defaultdict(float)
        for c, s, e, v in out:
            for bs in range(s // b * b, e, b):
                ov = min(e, bs + b) - max(s, bs)
                if ov > 0:
                    acc[(c, bs)] += v * ov / b
        out = [(c, bs, bs + b, round(vv, 5)) for (c, bs), vv in sorted(acc.items())]
    return {"text": io.write_bedgraph(out), "filename": "coverage.bedgraph",
            "message": f"{len(out)} segments, {normalisation} normalisation (factor {factor:.4g})"}


@T("deeptools_bamcompare", "bamCompare (log2 ratio of two samples)", "deeptools", "bedgraph",
   [sam("a", ALN_SAM), sam("b", ALN_SAM), choice("operation", ["log2", "ratio", "difference", "sub"]),
    intin("bin_size", 500), number("pseudocount", 1.0)], ex={"a": ALN_SAM, "b": ALN_SAM},
   up="iuc/deepTools", summary="Per-bin comparison of two alignment files.")
def bam_compare(a, b, operation="log2", bin_size=500, pseudocount=1.0):
    """``bamCompare``-style binwise ratios."""
    ha, alns_a = _load(a)
    hb, alns_b = _load(b)
    bsz = max(10, int(bin_size))

    def profile(alns):
        d = bam.depth(ha, alns, all_reads=True)
        out: dict[tuple[str, int], float] = defaultdict(float)
        for c, vals in d.items():
            for i, v in enumerate(vals):
                out[(c, i // bsz)] += v
        return {k: v / bsz for k, v in out.items()}
    pa, pb = profile(alns_a), profile(alns_b)
    rows = []
    for key in sorted(set(pa) | set(pb)):
        va, vb = pa.get(key, 0.0), pb.get(key, 0.0)
        if operation == "log2":
            val = math.log2((va + pseudocount) / (vb + pseudocount))
        elif operation == "ratio":
            val = (va + pseudocount) / (vb + pseudocount)
        elif operation == "difference":
            val = va - vb
        else:
            val = vb - va
        rows.append({"chrom": key[0], "start": key[1] * bsz, "end": (key[1] + 1) * bsz,
                     "value": round(val, 5)})
    out = [(r["chrom"], r["start"], r["end"], r["value"]) for r in rows]
    return {"text": io.write_bedgraph(out), "filename": "compare.bedGraph",
            "message": f"{len(out)} bins ({operation})"}


@T("deeptools_multibamsummary", "multiBamSummary", "deeptools", "table",
   [multi("files", [ALN_SAM], [ALN_SAM]), intin("bin_size", 1000), boolean("all_bins", False)],
   ex={"files": [ALN_SAM, ALN_SAM]}, up="iuc/deepTools",
   summary="Count reads per genomic bin for several BAM files at once.")
def multi_bam_summary(files=None, bin_size=1000, all_bins=False):
    """Sample x bin count matrix."""
    b = max(10, int(bin_size))
    per: list[dict] = []
    keys = set()
    for f in (files or []):
        hdr, alns = _load(f)
        d = bam.depth(hdr, alns, all_reads=all_bins)
        row: dict[tuple[str, int], float] = {}
        for c, vals in d.items():
            for i, v in enumerate(vals):
                if v or all_bins:
                    row[(c, i // b)] = v
        per.append(row)
        keys |= set(row)
    rows = []
    for key in sorted(keys):
        rec = {"chrom": key[0], "start": key[1] * b, "end": (key[1] + 1) * b}
        for i, p in enumerate(per):
            rec[f"sample_{i + 1}"] = p.get(key, 0.0)
        rows.append(rec)
    return table(rows[:2000], f"{len(rows)} bins x {len(per)} samples")


@T("deeptools_plotprofile", "plotProfile data (per-bin mean)", "deeptools", "table",
   [sam("src"), bed("regions", REGIONS_BED), intin("flank", 500), intin("bins", 20),
    boolean("per_interval", False)], ex={"src": ALN_SAM, "regions": REGIONS_BED},
   up="iuc/deepTools", summary="Average read density around a set of regions.")
def plot_profile(src, regions, flank=500, bins=20, per_interval=False):
    """Profile matrix around interval midpoints."""
    hdr, alns = _load(src)
    ivs = _bed_ivs(regions)
    depth = bam.depth(hdr, alns, all_reads=True)
    nb = max(2, int(bins))
    prof = [0.0] * nb
    rows = []
    for iv in ivs:
        mid = (iv.start + iv.end) // 2
        lo, hi = max(0, mid - int(flank)), mid + int(flank)
        vals = depth.get(iv.chrom, [])
        w = max(1, (hi - lo) // nb)
        line = [round(stats.mean([vals[p] if p < len(vals) else 0
                                  for p in range(lo + k * w, lo + (k + 1) * w)] or [0]), 4)
                for k in range(nb)]
        for k, v in enumerate(line):
            prof[k] += v
        rows.append({"region": iv.name or f"{iv.chrom}:{lo}-{hi}", **{f"bin_{k + 1}": v
                                                                     for k, v in enumerate(line)}})
    prof = [round(v / max(1, len(ivs)), 4) for v in prof]
    res = table(rows[:200], f"profile of {len(ivs)} regions in {nb} bins")
    res["figure"] = plot.line({"mean density": prof}, title="plotProfile", xlabel="flank bin",
                             ylabel="reads/bp")
    return res


def _bed_ivs(x):
    return io.parse_bed(io.as_text(x))


@T("deeptools_computeMatrix", "computeMatrix (regions x bins matrix)", "deeptools", "table",
   [sam("src"), bed("regions", REGIONS_BED), intin("bins", 10), intin("up", 300), intin("down", 300),
    boolean("scale", True)], ex={"src": ALN_SAM, "regions": REGIONS_BED, "bins": 8},
   up="iuc/deepTools", summary="Matrix of coverage values per region, for heatmaps.")
def compute_matrix(src, regions, bins=10, up=300, down=300, scale=True):
    """Region x bin coverage matrix."""
    hdr, alns = _load(src)
    depth = bam.depth(hdr, alns, all_reads=True)
    ivs = _bed_ivs(regions)
    nb = max(2, int(bins))
    rows = []
    for iv in ivs:
        vals = depth.get(iv.chrom, [])
        lo = max(0, iv.start - int(up))
        hi = iv.end + int(down)
        w = max(1, (hi - lo) // nb)
        line = [stats.mean([vals[p] if p < len(vals) else 0 for p in range(lo + k * w, lo + (k + 1) * w)]
                           or [0]) for k in range(nb)]
        if scale:
            mx = max(line) or 1
            line = [round(v / mx, 4) for v in line]
        rows.append({"region": iv.name or f"{iv.chrom}:{lo}-{hi}",
                     **{f"bin_{k + 1}": round(v, 4) for k, v in enumerate(line)}})
    return table(rows, f"matrix of {len(rows)} regions x {nb} bins")


@T("deeptools_plot_heatmap", "plotHeatmap figure", "deeptools", "image",
   [sam("src"), bed("regions", REGIONS_BED), intin("bins", 24), intin("up", 400), intin("down", 400)],
   ex={"src": ALN_SAM, "regions": REGIONS_BED, "bins": 20}, up="iuc/deepTools",
   summary="Render a coverage heatmap of sorted regions.")
def plot_heatmap(src, regions, bins=24, up=400, down=400):
    """Heatmap of the region x bin matrix."""
    hdr, alns = _load(src)
    depth = bam.depth(hdr, alns, all_reads=True)
    ivs = _bed_ivs(regions)
    nb = max(2, int(bins))
    mat = []
    for iv in ivs:
        vals = depth.get(iv.chrom, [])
        lo, hi = max(0, iv.start - int(up)), iv.end + int(down)
        w = max(1, (hi - lo) // nb)
        mat.append([stats.mean([vals[p] if p < len(vals) else 0 for p in range(lo + k * w, lo + (k + 1) * w)]
                               or [0]) for k in range(nb)])
    fig = plot.heatmap(mat, col_labels=[f"bin{i + 1}" for i in range(nb)],
                       row_labels=[iv.name or f"{iv.chrom}:{iv.start}" for iv in ivs],
                       title="plotHeatmap (coverage)")
    return figure(fig, f"heatmap of {len(mat)} regions x {nb} bins")


@T("deeptools_bam_fragmentsize", "bamPEFragmentSize summary", "deeptools", "stats",
   [sam("src"), intin("max_fraction", 2)], ex={"src": ALN_SAM}, up="iuc/deepTools",
   summary="Fragment-size statistics used by deepTools before extending reads.")
def fragment_size_summary(src, max_fraction=2):
    """Insert length distribution summary."""
    hdr, alns = _load(src)
    st = bam.insert_size_stats(alns)
    sizes = bam.insert_sizes(alns)
    out = {k: _num(v) for k, v in st.items()}
    out.update({"n_fragments": len(sizes), "min": min(sizes) if sizes else 0,
               "max": max(sizes) if sizes else 0,
               "median": stats.median(sizes) if sizes else 0,
               "recommended_extend": int(stats.mean(sizes) * max_fraction / 2) if sizes else 0})
    return values_message(out, "fragment sizes")


@T("deeptools_normalize_per_read", "Normalise a bedGraph per million reads", "deeptools", "bedgraph",
   [anyfile("src", "examples/coverage.bedgraph"), intin("total_reads", 1000000),
    number("scale", 1.0)], ex={"src": "examples/coverage.bedgraph"}, up="iuc/deepTools",
   summary="Rescale coverage values to counts per million mapped reads.")
def normalize_per_read(src, total_reads=1000000, scale=1.0):
    """CPM normalisation of a bedGraph."""
    vals = io.parse_bedgraph(io.as_text(src))
    n = max(1, int(total_reads))
    out = [(c, s, e, round(v * float(scale) * 1e6 / n, 5)) for c, s, e, v in vals]
    return {"text": io.write_bedgraph(out), "filename": "normalized.bedgraph",
            "message": f"normalised {len(out)} segments to CPM (n={n})"}


@T("samtools_reheader", "Replace the SAM header", SB, "text",
   [sam("src"), bigtext("header", "@HD\tVN:1.6\tSO:coordinate", label="New header lines"),
    boolean("keep_sq", False)], ex={"src": ALN_SAM}, up="iuc/samtools",
   summary="Swap the header of a SAM file (``samtools reheader``).")
def reheader(src, header="@HD\tVN:1.6", keep_sq=False):
    """Replace @-lines and keep the alignment records."""
    hdr, alns = _load(src)
    new = [l for l in io.lines(header) if l.strip()]
    if keep_sq:
        new += [h for h in hdr if h.startswith("@SQ")]
    body = "\n".join(new + [a.as_line() for a in alns])
    return {"text": body + "\n", "filename": "reheadered.sam",
            "message": f"new header with {len(new)} lines, {len(alns)} alignments kept"}


@T("samtools_fixmate", "Fix mate information", SB, "text", [sam("src"), boolean("assume_s1", False)],
   ex={"src": ALN_SAM}, up="iuc/samtools", summary="Populate RNEXT/PNEXT/TLEN for read pairs.")
def fix_mate(src, assume_s1=False):
    """``samtools fixmate``: insert size and mate coordinates."""
    hdr, alns = _load(src)
    by_name: dict[str, list] = defaultdict(list)
    for a in alns:
        by_name[a.qname].append(a)
    out = []
    n_fixed = 0
    for name, group in by_name.items():
        if len(group) == 2:
            group = sorted(group, key=lambda x: x.pos)
            tlen = group[1].pos - group[0].pos
            for i, a in enumerate(group):
                other = group[1 - i]
                new = io.Aln(a.qname, a.flag | 1, a.rname, a.pos, a.mapq, a.cigar, other.rname or "=",
                             other.pos, tlen if i == 0 else -tlen, a.seq, a.qual, a.tags)
                out.append(new)
                n_fixed += 1
        else:
            out.extend(group)
    return {"text": "\n".join(hdr + [a.as_line() for a in bam.sort_sam(out)]) + "\n",
            "filename": "fixed.sam", "message": f"fixed {n_fixed} mate annotations"}


@T("samtools_calmd", "Compute MD tags", SB, "table", [sam("src"), fa("genome", GENOME),
                                                   boolean("write_tags", False)],
   ex={"src": ALN_SAM, "genome": GENOME}, up="iuc/samtools",
   summary="Derive MD/NM-like mismatch strings against the reference.")
def calmd(src, genome=GENOME, write_tags=False):
    """Per-alignment mismatch positions relative to the reference."""
    seqs, sizes = _refs(genome)
    hdr, alns = _load(src)
    rows = []
    for a in alns:
        if not a.mapped or a.rname not in seqs:
            continue
        ref = seqs[a.rname][a.pos - 1:a.pos - 1 + len(a.seq)]
        mm = [i for i, (x, y) in enumerate(zip(a.seq.upper(), ref.upper())) if x != y]
        rows.append({"read": a.qname, "reference": a.rname, "position": a.pos, "n_mismatches": len(mm),
                    "md": ",".join(map(str, mm[:20])), "query_len": len(a.seq), "ref_len": len(ref)})
    res = table(rows[:500], f"MD computed for {len(rows)} alignments")
    if write_tags:
        res["text"] = "\n".join(hdr + [a.as_line() for a in alns[:20]])
    return res


@T("samtools_read_bam_nx", "Read-level stats table", SB, "table",
   [sam("src"), intin("max_rows", 200), intin("min_mapq", 0)], ex={"src": ALN_SAM},
   up="iuc/sambamba", summary="One row per alignment: length, NM, clipping, insert size.")
def read_table_from_sam(src, max_rows=200, min_mapq=0):
    """Tabular view of the alignment records."""
    hdr, alns = _load(src)
    rows = []
    for a in alns:
        if a.mapq < min_mapq:
            continue
        rows.append({"read": a.qname, "flag": a.flag, "reference": a.rname, "position": a.pos,
                    "mapq": a.mapq, "cigar": a.cigar, "length": len(a.seq or ""),
                    "ref_span": io.ref_length(a.cigar) if a.cigar else 0,
                    "nm": (a.tags or {}).get("NM", ""), "tlen": a.tlen,
                    "mapped": a.mapped, "reverse": a.reverse, "duplicate": a.duplicate})
    return table(rows[: int(max_rows)], f"{len(rows)} alignments listed")


@T("samtools_coverage_bed", "Intervals covered by at least N reads", SB, "bed",
   [sam("src"), intin("min_depth", 1), number("max_depth", 0.0), boolean("sort", True)],
   ex={"src": ALN_SAM, "min_depth": 5}, up="iuc/bedtools",
   summary="bedtools genomecov -bga style covered/uncovered blocks.")
def covered_blocks(src, min_depth=1, max_depth=0.0, sort=True):
    """Merge positions meeting a depth range into BED blocks."""
    hdr, alns = _load(src)
    d = bam.depth(hdr, alns, all_reads=True)
    out = []
    for chrom, vals in d.items():
        start = None
        for i, v in enumerate(vals):
            ok = v >= min_depth and (not max_depth or v <= max_depth)
            if ok and start is None:
                start = i
            elif not ok and start is not None:
                out.append(io.Interval(chrom, start, i, f"depth>={min_depth}", 0, "."))
                start = None
        if start is not None:
            out.append(io.Interval(chrom, start, len(vals), f"depth>={min_depth}", 0, "."))
    if sort:
        out = genome.sort_ivs(out)
    return _bed(out, f"{len(out)} covered blocks (depth >= {min_depth})")


@T("samtools_bam_sanity", "Alignment sanity checks (NGS QC)", SB, "table",
   [sam("src"), number("max_unmapped_fraction", 0.5), number("min_proper_pair_fraction", 0.5),
    number("max_duplicate_fraction", 0.6)], ex={"src": ALN_SAM}, up="iuc/multiqc",
   summary="Traffic-light table of mapping metrics against thresholds.")
def alignment_sanity(src, max_unmapped_fraction=0.5, min_proper_pair_fraction=0.5,
                     max_duplicate_fraction=0.6):
    """Flag or pass a BAM based on standard thresholds."""
    hdr, alns = _load(src)
    st = bam.flagstat(hdr, alns)
    total = max(1, int(st.get("total", len(alns))))
    mapped = int(st.get("mapped", 0))
    dupes = int(st.get("duplicates", 0))
    proper = int(st.get("properly_paired", 0))
    rows = [{"metric": "unmapped fraction", "value": round((total - mapped) / total, 4),
             "threshold": f"< {max_unmapped_fraction}", "status": "PASS" if (total - mapped) / total
             <= max_unmapped_fraction else "FAIL"},
            {"metric": "proper pair fraction", "value": round(proper / total, 4),
             "threshold": f"> {min_proper_pair_fraction}",
             "status": "PASS" if proper / total >= min_proper_pair_fraction else "FAIL"},
            {"metric": "duplicate fraction", "value": round(dupes / total, 4),
             "threshold": f"< {max_duplicate_fraction}",
             "status": "PASS" if dupes / total <= max_duplicate_fraction else "FAIL"}]
    n_fail = sum(1 for r in rows if r["status"] == "FAIL")
    return table(rows, f"{len(rows) - n_fail}/{len(rows)} checks pass")
