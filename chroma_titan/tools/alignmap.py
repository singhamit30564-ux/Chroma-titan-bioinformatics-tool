"""Read mapping, de novo assembly, multiple alignment and BLAST-style search.

Panel sections covered: **mapping**, **assembly**, **multiple_alignments**,
**ncbi_blast_**, **nanopore**, **sequence_contamination_filtering**, **lift_over**.
"""
from __future__ import annotations

import math
from collections import Counter, defaultdict

from chroma_titan.core import align, bam, genome, io, plot, seq, stats, tables
from chroma_titan.tools._common import *  # noqa: F401,F403
from chroma_titan.tools._common import (ALN_SAM, GENES_FA, GENOME, MSA, PROTEINS, READS_1,
                                        READS_2, READS_SINGLE, REF_PROTEINS, TAXMAP)

MAP = "mapping"
ASM = "assembly"
MSA_S = "multiple_alignments"
BLAST = "ncbi_blast_"
NANO = "nanopore"
CONTAM = "sequence_contamination_filtering"
LIFT = "lift_over"


# ---------------------------------------------------------------------------
# local helpers
# ---------------------------------------------------------------------------
def _reads(src) -> list:
    return io.parse_fastq(io.as_text(src))


def _fa(src) -> list:
    return io.parse_fasta(io.as_text(src))


def _rows(src) -> list[dict]:
    """FASTA (or pasted sequence) -> list of ``{id, seq}`` dicts."""
    return [{"id": r.id, "seq": r.seq} for r in _fa(src)]


def _qvals(r) -> list[int]:
    return [ord(c) - 33 for c in (r.qual or "")] if r.qual else []


def _hist(values, bins=10):
    """Simple equal-width histogram -> {label: count}."""
    vals = [float(v) for v in values if isinstance(v, (int, float)) and math.isfinite(float(v))]
    if not vals:
        return {}
    lo, hi = min(vals), max(vals)
    if hi <= lo:
        return {f"{lo:g}": len(vals)}
    w = (hi - lo) / max(1, int(bins))
    out: dict[str, int] = {}
    for v in vals:
        i = min(int((v - lo) / w), int(bins) - 1)
        out[f"{lo + i * w:g}-{lo + (i + 1) * w:g}"] = out.get(f"{lo + i * w:g}-{lo + (i + 1) * w:g}", 0) + 1
    return out


def _aln_rows(x):
    """Alignment rows: accept a list of dicts or aligned FASTA text."""
    if isinstance(x, list):
        return [{"id": r.get("id", f"seq_{i + 1}") if isinstance(r, dict) else r.id,
                 "seq": r.get("seq", "") if isinstance(r, dict) else r.seq}
                for i, r in enumerate(x)]
    return [{"id": r.id, "seq": r.seq} for r in io.parse_fasta(io.as_text(x))]


def _map(genome_src, reads_src, k=21, mism=3, local=False, stride=1, min_mapq=0):
    recs = align.map_reads(_rows(genome_src), [{"id": r.id, "seq": r.seq, "qual": r.qual}
                                               for r in _reads(reads_src)],
                          k=int(k), max_mismatches=int(mism), seed_stride=int(stride),
                          local=bool(local), min_mapq=int(min_mapq))
    return recs


def _sam_hdr(genome_src):
    seqs, sizes = genome.read_genome(genome_src)
    return [(name, sizes.get(name, len("".join(parts)))) for name, parts in seqs.items()]


# ===========================================================================
# short-read mapping
# ===========================================================================
@T("bwa_mem_lite", "Map reads to a reference (BWA-MEM style)", MAP, "file",
   [fq("reads", READS_SINGLE, "Reads (FASTQ)"), fa("ref", GENOME, "Reference FASTA"),
    intin("k", 21, "Seed length k", help="Seeds shorter than k are never extended", min=6),
    intin("max_mismatches", 3, "Maximum mismatches per read", min=0),
    intin("seed_stride", 1, "Seed stride", min=1), intin("min_mapq", 0, "Minimum MAPQ", min=0, max=60),
    boolean("local", False, "Soft-clip ends (local alignment)")],
   ex={"reads": READS_1, "ref": GENOME, "k": 21},
   up="bwa mem", tags=("DNA", "mapping", "SAM"),
   summary="Seed-and-extend mapping of FASTQ reads onto a FASTA reference, SAM output.")
def bwa_mem_lite(reads, ref, k=21, max_mismatches=3, seed_stride=1, min_mapq=0, local=False):
    """``bwa mem``-style mapping producing a SAM record set."""
    recs = _map(ref, reads, k, max_mismatches, local, seed_stride, min_mapq)
    mapped = [r for r in recs if r["rname"] != "*"]
    body = align.map_to_sam(recs, _sam_hdr(ref))
    return {"text": body, "filename": "aln.sam",
            "message": f"{len(mapped)}/{len(recs)} reads mapped "
                       f"(mean identity {round(stats.mean([r['identity'] for r in mapped]), 3) if mapped else 0}%)",
            "stats": {"total_reads": len(recs), "mapped": len(mapped),
                      "unmapped": len(recs) - len(mapped),
                      "mapping_rate_percent": round(100 * len(mapped) / max(1, len(recs)), 3),
                      "mean_mapq": round(stats.mean([r["mapq"] for r in mapped]), 2) if mapped else 0,
                      "mean_mismatches": round(stats.mean([r["tags"].get("NM", 0) for r in mapped]), 4)
                      if mapped else 0}}


@T("bowtie2_end_to_end", "End-to-end mapping with Bowtie2 scoring", MAP, "table",
   [fq("reads", READS_SINGLE, "Reads (FASTQ)"), fa("ref", GENOME, "Reference FASTA"),
    intin("k", 31, "Seed length", min=6), intin("max_mismatches", 4, "Max mismatches", min=0),
    intin("head", 200, "Rows to show", min=1)],
   ex={"reads": READS_1, "ref": GENOME, "k": 25},
   up="bowtie2 --end-to-end", tags=("DNA", "mapping"),
   summary="Full-length read alignment table with edit distance, MAPQ and per-read score.")
def bowtie2_end_to_end(reads, ref, k=31, max_mismatches=4, head=200):
    """End-to-end mapping: every read base must be aligned."""
    recs = _map(ref, reads, k, max_mismatches)
    rows = []
    for r in recs:
        qlen = len(r["seq"]) or 1
        nm = r["tags"].get("NM", 0)
        rows.append({"read": r["qname"], "flag": r["flag"], "reference": r["rname"],
                     "position": r["pos"], "mapq": r["mapq"], "cigar": r["cigar"],
                     "edit_distance": nm, "identity_percent": r["identity"],
                     "score": round(-6 * nm - 2 * (qlen - nm) * 0.0, 3),
                     "end_to_end": r["cigar"].count("S") == 0})
    n_ok = sum(1 for r in rows if r["reference"] != "*")
    return table(rows[:int(head)], f"{n_ok}/{len(rows)} reads aligned end-to-end")


@T("hisat2_spliced", "Spliced / soft-clipped mapping (HISAT2 style)", MAP, "table",
   [fq("reads", READS_1, "Reads (FASTQ)"), fa("ref", GENES_FA, "Transcriptome FASTA"),
    intin("k", 17, "Seed length", min=6), intin("max_mismatches", 5, "Max mismatches", min=0),
    intin("head", 200, "Rows to show", min=1)],
   ex={"reads": READS_1, "ref": GENES_FA, "k": 17},
   up="hisat2", tags=("RNA", "mapping", "spliced"),
   summary="Local mapping with clipped ends; clipped bases are reported as putative junctions.")
def hisat2_spliced(reads, ref, k=17, max_mismatches=5, head=200):
    """Splice-aware stand-in: local alignment with soft clipping."""
    recs = _map(ref, reads, k, max_mismatches, local=True)
    rows = []
    for r in recs:
        clips = [int(n) for n, op in io.cigar_ops(r["cigar"]) if op == "S"]
        rows.append({"read": r["qname"], "reference": r["rname"], "position": r["pos"],
                     "mapq": r["mapq"], "cigar": r["cigar"], "clipped_bases": sum(clips),
                     "n_clips": len(clips), "identity_percent": r["identity"],
                     "spliced": sum(clips) > 10})
    n_sp = sum(1 for r in rows if r["spliced"])
    return table(rows[:int(head)], f"{len(rows)} alignments, {n_sp} with junction-scale clipping")


@T("bwa_index_summary", "Summarise a k-mer index of the reference", MAP, "table",
   [fa("ref", GENOME, "Reference FASTA"), intin("k", 11, "k", min=4), intin("top", 20, "Top k-mers", min=1)],
   ex={"ref": GENOME, "k": 10}, up="bwa index", tags=("index", "k-mer"),
   summary="Distinct k-mers, repeated k-mers and the most frequent seeds of an index.")
def bwa_index_summary(ref, k=11, top=20):
    """Report the k-mer index of a reference."""
    idx = align.build_kmer_index(_rows(ref), int(k))
    counts = Counter({km: len(v) for km, v in idx.items()})
    total = sum(counts.values())
    rows = [{"kmer": km, "occurrences": c, "first_hit": "n/a"} for km, c in
            counts.most_common(int(top))]
    return {"message": f"{len(counts)} distinct {k}-mers over {total} positions", "rows": rows,
            "stats": {"sequences": len(_rows(ref)), "distinct_kmers": len(counts),
                      "total_kmer_positions": total,
                      "unique_kmers": sum(1 for c in counts.values() if c == 1),
                      "repeated_kmers": sum(1 for c in counts.values() if c > 1),
                      "max_occurrences": max(counts.values()) if counts else 0}}


@T("kmer_pseudoalign", "k-mer pseudoalignment to transcripts", MAP, "table",
   [fq("reads", READS_1, "Reads (FASTQ)"), fa("transcripts", GENES_FA, "Transcript FASTA"),
    intin("k", 21, "k-mer size", min=8), number("min_frac", 0.4, "Minimum shared k-mer fraction",
                                                min=0.0, max=1.0), intin("head", 100, "Rows", min=1)],
   ex={"reads": READS_1, "transcripts": GENES_FA, "k": 21},
   up="kallisto / salmon", tags=("pseudoalignment", "quantification"),
   summary="Assign every read to the transcript sharing the most k-mers and tabulate counts.")
def kmer_pseudoalign(reads, transcripts, k=21, min_frac=0.4, head=100):
    """kallisto-flavoured pseudoalignment."""
    tr = _rows(transcripts)
    index = {t["id"]: set(seq.kmer_counts(t["seq"], int(k), canonical=True)) for t in tr}
    counts: Counter = Counter()
    assigned = 0
    for r in _reads(reads):
        ks = set(seq.kmer_counts(r.seq, int(k), canonical=True))
        if not ks:
            continue
        best, best_frac = None, 0.0
        for tid, tset in index.items():
            frac = len(ks & tset) / len(ks)
            if frac > best_frac:
                best, best_frac = tid, frac
        if best is not None and best_frac >= float(min_frac):
            counts[best] += 1
            assigned += 1
    rows = [{"transcript": t["id"], "length": len(t["seq"]), "reads": counts.get(t["id"], 0),
             "eff_len_kb": max(0.1, (len(t["seq"]) - int(k) + 1) / 1000.0),
             "tpm": 0.0} for t in tr]
    denom = sum(r["reads"] / r["eff_len_kb"] for r in rows) or 1.0
    for r in rows:
        r["tpm"] = round(1e6 * (r["reads"] / r["eff_len_kb"]) / denom, 3)
    return table(rows[:int(head)], f"{assigned} reads pseudoaligned to {len(tr)} transcripts")


@T("samtools_stats_lite", "Mapping statistics of a SAM file", MAP, "table",
   [sam("src", ALN_SAM, "SAM file")], ex={"src": ALN_SAM}, up="samtools stats",
   tags=("statistics", "mapping"),
   summary="Read totals, mapping rate, mismatch and MAPQ means for an alignment file.")
def samtools_stats_lite(src):
    """``samtools stats``-like summary."""
    hdr, alns = bam.load(src)
    st = align.align_to_reference_stats(alns)
    rows = [{"metric": k, "value": v} for k, v in st.items()]
    return table(rows, f"{st.get('mapped', 0)}/{st.get('reads', 0)} records mapped")


@T("identity_histogram", "Distribution of alignment identity", MAP, "figure",
   [sam("src", ALN_SAM, "SAM file"), intin("bin_size", 5, "Bin width (identity %)", min=1),
    choice("source", ["sam", "reads"], "sam", "Histogram source")],
   ex={"src": ALN_SAM, "bin_size": 10}, up="samtools stats / plot",
   tags=("histogram", "visualization"),
   summary="Percent-identity histogram of every alignment in the file.")
def identity_histogram(src, bin_size=5, source="sam"):
    """Identity histogram of alignments."""
    hdr, alns = bam.load(src)
    hist = align.similarity_histogram([a for a in alns if a.mapped], bin_size=int(bin_size))
    labels = [str(k) for k in sorted(hist)]
    values = [hist[k] for k in sorted(hist)]
    return plot.bar(labels, values, xlabel=f"identity bin (width {bin_size})", ylabel="alignments",
                   title="Alignment identity distribution")


@T("merge_pairs_lite", "Merge overlapping read pairs", MAP, "table",
   [fq("reads1", READS_1, "R1 FASTQ"), fq("reads2", READS_2, "R2 FASTQ"),
    intin("min_overlap", 11, "Minimum overlap", min=2), intin("max_mismatches", 2, "Max mismatches in overlap", min=0),
    intin("head", 100, "Rows", min=1)],
   ex={"reads1": READS_1, "reads2": READS_2, "min_overlap": 10}, up="fastq-join / vsearch mergepairs",
   tags=("paired-end", "merging"),
   summary="Detect mate overlaps, merge into one read and report the merging rate.")
def merge_pairs_lite(reads1, reads2, min_overlap=11, max_mismatches=2, head=100):
    """Mate merging with mismatch-aware overlap search."""
    res = align.merge_pairs([{"id": r.id, "seq": r.seq, "qual": r.qual} for r in _reads(reads1)],
                           [{"id": r.id, "seq": r.seq, "qual": r.qual} for r in _reads(reads2)],
                           min_overlap=int(min_overlap), max_mismatches=int(max_mismatches))
    merged = res.get("merged", [])
    rows = [{"read": m["id"], "merged_length": len(m["seq"]), "overlap": m["overlap"],
             "mismatches": m["mismatches"], "sequence": m["seq"][:60] + ("..." if len(m["seq"]) > 60 else "")}
            for m in merged]
    return {"message": f"{res['n_merged']} merged / {res['n_unmerged']} unmerged "
                       f"(rate {res['merge_rate']}%)", "rows": rows[:int(head)],
            "stats": {k: v for k, v in res.items() if not isinstance(v, list)}}


@T("insert_size_stats", "Insert size (TLEN) distribution", MAP, "table",
   [sam("src", ALN_SAM, "SAM file"), intin("bins", 10, "Histogram bins", min=2),
    intin("sd_cut", 3.0, "Exclude beyond n SD", min=1.0)],
   ex={"src": ALN_SAM, "bins": 8}, up="picard CollectsInsertSizeMetrics",
   tags=("insert size", "paired-end"),
   summary="Median, mean, standard deviation and outlier counts of the fragment insert sizes.")
def insert_size_stats(src, bins=10, sd_cut=3.0):
    """Insert size distribution from proper pairs."""
    hdr, alns = bam.load(src)
    sizes = sorted(abs(a.tlen) for a in alns if a.mapped and a.tlen)
    if not sizes:
        return table([], "no TLEN values in this SAM (paired data required)")
    mu, sd = stats.mean(sizes), stats.stdev(sizes)
    keep = [s for s in sizes if abs(s - mu) <= float(sd_cut) * sd] if sd else sizes
    rows = [{"metric": "median_insert_size", "value": stats.median(sizes)},
            {"metric": "mean_insert_size", "value": round(mu, 3)},
            {"metric": "std_dev", "value": round(sd, 3)},
            {"metric": "min", "value": min(sizes)}, {"metric": "max", "value": max(sizes)},
            {"metric": "pairs", "value": len(sizes)},
            {"metric": "outliers_removed", "value": len(sizes) - len(keep)},
            {"metric": "within_100bp_of_mean", "value": round(
                100 * sum(1 for s in sizes if abs(s - mu) <= 100) / len(sizes), 3)}]
    res = table(rows, f"median insert size {stats.median(sizes)} over {len(sizes)} pairs")
    res["stats"]["histogram"] = _hist(keep, bins)
    return res


@T("softclip_report", "Soft-clipped bases per alignment", MAP, "table",
   [sam("src", ALN_SAM, "SAM file"), intin("min_clip", 10, "Report clips of at least", min=1),
    intin("head", 100, "Rows", min=1)],
   ex={"src": ALN_SAM, "min_clip": 5}, up="samtools view / clipping analysis",
   tags=("soft clip", "structure"),
   summary="List alignments whose clipped tails suggest chimeric reads or junctions.")
def softclip_report(src, min_clip=10, head=100):
    """Summarise soft clipping across the alignment file."""
    hdr, alns = bam.load(src)
    rows = []
    for a in alns:
        if not a.mapped:
            continue
        ops = io.cigar_ops(a.cigar)
        left = ops[0][0] if ops and ops[0][1] == "S" else 0
        right = ops[-1][0] if len(ops) > 1 and ops[-1][1] == "S" else 0
        if left + right >= int(min_clip):
            rows.append({"read": a.qname, "reference": a.rname, "position": a.pos,
                        "left_clip": left, "right_clip": right, "total_clipped": left + right,
                        "cigar": a.cigar, "mapq": a.mapq,
                        "clipped_sequence": (a.seq[:left] + "..." + a.seq[len(a.seq) - right:]
                                             if right else a.seq[:left])[:80]})
    rows.sort(key=lambda r: -r["total_clipped"])
    return table(rows[:int(head)], f"{len(rows)} alignments with clips >= {min_clip} bp")


@T("mapq_filter", "Filter alignments by MAPQ", MAP, "file",
   [sam("src", ALN_SAM, "SAM file"), intin("min_mapq", 30, "Minimum MAPQ", min=0, max=60),
    boolean("drop_unmapped", False, "Also drop unmapped reads")],
   ex={"src": ALN_SAM, "min_mapq": 20}, up="samtools view -q",
   tags=("filter", "SAM"),
   summary="Write a SAM containing only alignments above a MAPQ threshold.")
def mapq_filter(src, min_mapq=30, drop_unmapped=False):
    """``samtools view -q`` style filtering."""
    hdr, alns = bam.load(src)
    keep = [a for a in alns if a.mapq >= int(min_mapq) and (a.mapped or not drop_unmapped)]
    body = "\n".join(hdr + [a.as_line() for a in keep]) + "\n"
    return {"text": body, "filename": "filtered.sam",
            "message": f"{len(keep)}/{len(alns)} alignments kept at MAPQ >= {min_mapq}"}


@T("unmapped_to_fastq", "Extract unmapped reads", MAP, "file",
   [sam("src", ALN_SAM, "SAM file"), boolean("singletons", True, "Write singletons too")],
   ex={"src": ALN_SAM}, up="samtools fastq -f 4",
   tags=("FASTQ", "extraction"),
   summary="Recover unmapped sequences from an alignment as FASTQ for re-analysis.")
def unmapped_to_fastq(src, singletons=True):
    """Extract unmapped alignments as FASTQ."""
    hdr, alns = bam.load(src)
    out = []
    for a in alns:
        if a.mapped or (not singletons and not a.get("unmapped_mate", False)):
            continue
        if not a.seq or a.seq == "*":
            continue
        out.append(f"@{a.qname}\n{a.seq}\n+\n{a.qual or 'I' * len(a.seq)}")
    body = "\n".join(out) + ("\n" if out else "")
    return {"text": body, "filename": "unmapped.fastq",
            "message": f"{len(out)} unmapped reads extracted"}


# ===========================================================================
# long reads / nanopore
# ===========================================================================
@T("nanopore_run_stats", "Long-read length and quality summary", NANO, "table",
   [fq("src", READS_SINGLE, "FASTQ file"), intin("min_length", 0, "Minimum length", min=0),
    intin("bins", 8, "Length histogram bins", min=2)],
   ex={"src": READS_SINGLE, "bins": 6}, up="nanoplot / poretools",
   tags=("nanopore", "quality control"),
   summary="Total bases, N50, mean quality and Q20 fraction of a long-read run.")
def nanopore_run_stats(src, min_length=0, bins=8):
    """Length/quality metrics of a nanopore run."""
    rs = [r for r in _reads(src) if len(r.seq) >= int(min_length)]
    lens = [len(r.seq) for r in rs]
    quals = [stats.mean(_qvals(r)) for r in rs if r.qual]
    if not rs:
        return table([], "no reads passed the length filter")
    rows = [{"metric": "reads", "value": len(rs)}, {"metric": "total_bases", "value": sum(lens)},
            {"metric": "mean_length", "value": round(stats.mean(lens), 2)},
            {"metric": "median_length", "value": stats.median(lens)},
            {"metric": "longest_read", "value": max(lens)},
            {"metric": "N50", "value": seq.n50(lens)}, {"metric": "L50", "value": seq.l50(lens)},
            {"metric": "mean_quality", "value": round(stats.mean(quals), 3) if quals else 0},
            {"metric": "fraction_q20", "value": round(sum(1 for q in quals if q >= 20) / len(quals), 4)
             if quals else 0},
            {"metric": "mean_homopolymer_run", "value": round(stats.mean(
                [float(seq.longest_homopolymer(r.seq).get("length", 0) or 0) for r in rs]), 3)},
            {"metric": "mean_gc_percent", "value": round(100 * stats.mean(
                [seq.gc_content(r.seq) for r in rs]), 3)}]
    res = table(rows, f"{len(rs)} reads, N50 = {seq.n50(lens)}")
    res["stats"]["length_histogram"] = _hist(lens, bins)
    return res


@T("nanopore_trim_adapters", "Trim adapters and low-quality tails", NANO, "table",
   [fq("src", READS_SINGLE, "FASTQ file"), bigtext("adapters", "AGATCGGAAGAGC", "Adapter sequences (one per line)"),
    intin("quality_cut", 7, "Trim tail below this Q", min=0), intin("min_length", 50, "Minimum read length", min=1)],
   ex={"src": READS_1, "quality_cut": 10, "min_length": 60}, up="porechop / fastp",
   tags=("nanopore", "trimming"),
   summary="Remove adapter occurrences and quality-trim read ends, keeping long enough reads.")
def nanopore_trim_adapters(src, adapters="AGATCGGAAGAGC", quality_cut=7, min_length=50):
    """Adapter and quality trimming for long reads."""
    adapt = [a.strip().upper() for a in (adapters or "").splitlines() if a.strip()]
    rows = []
    kept = 0
    for r in _reads(src):
        s, q = r.seq.upper(), r.qual or ""
        n_ad = 0
        for a in adapt:
            while True:
                i = s.find(a)
                if i < 0:
                    break
                s = s[:i]
                n_ad += 1
        arr = [ord(c) - 33 for c in q][: len(s)]
        while arr and arr[-1] < int(quality_cut):
            arr.pop()
            s = s[:-1]
        dropped = len(r.seq) - len(s)
        keep = len(s) >= int(min_length)
        kept += 1 if keep else 0
        rows.append({"read": r.id, "original_length": len(r.seq), "trimmed_length": len(s),
                    "bases_removed": dropped, "adapters_found": n_ad,
                    "mean_quality": round(stats.mean(arr), 2) if arr else 0,
                    "kept": keep})
    return table(rows, f"{kept}/{len(rows)} reads retained")


@T("nanopore_quality_profile", "Per-cycle quality profile", NANO, "figure",
   [fq("src", READS_SINGLE, "FASTQ file"), intin("min_length", 0, "Minimum read length", min=0)],
   ex={"src": READS_1}, up="porechop / NanoPlot", tags=("nanopore", "quality"),
   summary="Mean base quality at every cycle of the run.")
def nanopore_quality_profile(src, min_length=0):
    """Mean quality per read cycle."""
    rs = [r for r in _reads(src) if len(r.seq) >= int(min_length) and r.qual]
    if not rs:
        return plot.line({"quality": [0.0]}, title="no quality strings found")
    width = min(len(r.qual) for r in rs)
    prof = []
    for i in range(width):
        prof.append(stats.mean([ord(r.qual[i]) - 33 for r in rs]))
    return plot.line({"mean_quality": prof}, x=list(range(1, width + 1)),
                     xlabel="cycle", ylabel="mean quality", title="Per-cycle quality profile")


@T("longread_align", "Map long reads onto a reference", NANO, "table",
   [fq("reads", READS_SINGLE, "Long-read FASTQ"), fa("ref", GENOME, "Reference FASTA"),
    intin("k", 15, "Seed length k", min=6), intin("max_mismatches", 12, "Max mismatches", min=0),
    intin("head", 100, "Rows", min=1)],
   ex={"reads": READS_SINGLE, "ref": GENOME, "k": 13, "max_mismatches": 20},
   up="minimap2", tags=("nanopore", "mapping"),
   summary="Minimap2-style long-read mapping with tolerant seed extension.")
def longread_align(reads, ref, k=15, max_mismatches=12, head=100):
    """Long read mapping with error-tolerant seeds."""
    recs = _map(ref, reads, k, max_mismatches)
    rows = [{"read": r["qname"], "reference": r["rname"], "position": r["pos"], "strand":
             "-" if r["flag"] & 16 else "+", "mapq": r["mapq"], "cigar": r["cigar"],
             "identity_percent": r["identity"], "mismatches": r["tags"].get("NM", 0),
             "read_length": len(r["seq"])} for r in recs]
    n_ok = sum(1 for r in rows if r["reference"] != "*")
    return table(rows[:int(head)], f"{n_ok}/{len(rows)} long reads placed")


@T("mash_screen", "Mash / MinHash distance between sequence sets", CONTAM, "table",
   [fa("query", GENOME, "Query FASTA"), fa("database", GENOME, "Target FASTA"),
    intin("k", 21, "k-mer size", min=4), intin("num", 500, "Hashes", min=50)],
   ex={"query": GENOME, "database": GENOME, "k": 21, "num": 300}, up="mash dist",
   tags=("mash", "distance", "sketching"),
   summary="MinHash sketch of every sequence and its Mash distance to each target.")
def mash_screen(query, database, k=21, num=500):
    """Mash distances between two FASTA sets."""
    q, d = _rows(query), _rows(database)
    sketch = {r["id"]: seq.minhash(r["seq"], int(k), int(num)) for r in q + d}
    rows = []
    for a in q:
        sa = sketch[a["id"]]
        for b in d:
            sb = sketch[b["id"]]
            jac = seq.minhash_jaccard(sa, sb)
            dist = seq.mash_distance(sa, sb, int(k))
            rows.append({"query": a["id"], "target": b["id"], "jaccard": round(jac, 5),
                        "mash_distance": round(float(dist), 6),
                        "shared_hashes": len(set(sa) & set(sb)), "sketch_size": len(sa)})
    rows.sort(key=lambda r: r["mash_distance"])
    return table(rows[:400], f"{len(rows)} query x target comparisons")


# ===========================================================================
# de novo assembly
# ===========================================================================
def _debruijn(seqs, k, min_count=1):
    """Tiny de Bruijn graph: unitigs + graph statistics."""
    count: Counter = Counter()
    nxt: dict[str, Counter] = defaultdict(Counter)
    prev: dict[str, Counter] = defaultdict(Counter)
    for s in seqs:
        s = seq.clean(s)
        for i in range(0, max(0, len(s) - k + 1)):
            a = s[i:i + k]
            if "N" in a:
                continue
            count[a] += 1
            if i + 1 <= len(s) - k:
                b = s[i + 1:i + k + 1]
                if "N" not in b:
                    nxt[a][b] += 1
                    prev[b][a] += 1
    nodes = {a: c for a, c in count.items() if c >= int(min_count)}
    unitigs, used = [], set()
    for start in nodes:
        if start in used:
            continue
        path = [start]
        used.add(start)
        cur = start
        while True:
            outs = [n for n in nxt.get(cur, ()) if n in nodes]
            ins = [p for p in prev.get(cur, ()) if p in nodes]
            if len(outs) != 1 or len(ins) != 1 or outs[0] in used:
                break
            cur = outs[0]
            if len([n for n in nxt.get(cur, ()) if n in nodes]) != 1:
                used.add(cur)
                path.append(cur)
                break
            if len([p for p in prev.get(cur, ()) if p in nodes]) != 1:
                used.add(cur)
                path.append(cur)
                break
            used.add(cur)
            path.append(cur)
        unitigs.append(path[0] + "".join(x[-1] for x in path[1:]))
    stats_d = {"kmers_total": len(count), "solid_kmers": len(nodes),
              "edges": sum(len([n for n in nxt.get(a, ()) if n in nodes]) for a in nodes),
              "branching": sum(1 for a in nodes if len([n for n in nxt.get(a, ()) if n in nodes]) > 1),
              "tips": sum(1 for a in nodes if not [n for n in nxt.get(a, ()) if n in nodes]),
              "mean_kmer_coverage": round(stats.mean([nodes[a] for a in nodes]), 3) if nodes else 0.0}
    return unitigs, stats_d


@T("debruijn_unitigs", "De Bruijn graph unitigs from reads", ASM, "table",
   [fq("reads", READS_SINGLE, "Reads (FASTQ)"), intin("k", 15, "k-mer size", min=5),
    intin("min_coverage", 1, "Minimum k-mer coverage", min=1), intin("min_length", 40, "Minimum contig length", min=1),
    intin("head", 200, "Rows", min=1)],
   ex={"reads": READS_SINGLE, "k": 15, "min_length": 60}, up="spades / velvet",
   tags=("assembly", "de Bruijn"),
   summary="Build a de Bruijn graph of read k-mers and emit maximal non-branching unitigs.")
def debruijn_unitigs(reads, k=15, min_coverage=1, min_length=40, head=200):
    """Unitig construction from a k-mer graph."""
    seqs = [r.seq for r in _reads(reads)]
    unitigs, st = _debruijn(seqs, int(k), int(min_coverage))
    keep = [u for u in unitigs if len(u) >= int(min_length)]
    rows = [{"contig": f"contig_{i + 1}", "length": len(u), "gc_percent": round(100 * seq.gc_content(u), 3),
             "kmers": max(1, len(u) - int(k) + 1)} for i, u in enumerate(keep)]
    lens = [len(u) for u in keep]
    return {"message": f"{len(keep)} contigs >= {min_length} bp (N50 {seq.n50(lens) if lens else 0})",
            "rows": rows[:int(head)],
            "stats": {**st, "contigs": len(keep), "total_length": sum(lens),
                      "N50": seq.n50(lens) if lens else 0, "L50": seq.l50(lens) if lens else 0,
                      "longest_contig": max(lens) if lens else 0}}


@T("assembly_contigs_fasta", "Write the unitig assembly as FASTA", ASM, "file",
   [fq("reads", READS_SINGLE, "Reads (FASTQ)"), intin("k", 15, "k-mer size", min=5),
    intin("min_length", 40, "Minimum contig length", min=1),
    boolean("orientation", False, "Add GC-based orientation to names")],
   ex={"reads": READS_SINGLE, "k": 15, "min_length": 60}, up="spades contigs.fasta",
   tags=("assembly", "FASTA"),
   summary="Export assembled contigs in FASTA with length and coverage in the headers.")
def assembly_contigs_fasta(reads, k=15, min_length=40, orientation=False):
    """Contigs to FASTA text."""
    unitigs, st = _debruijn([r.seq for r in _reads(reads)], int(k))
    lines = []
    n = 0
    for u in unitigs:
        if len(u) < int(min_length):
            continue
        n += 1
        tag = f" gc={round(100 * seq.gc_content(u), 1)}%" if orientation else ""
        lines.append(f">contig_{n} length={len(u)} kmercov={st['mean_kmer_coverage']}{tag}")
        for i in range(0, len(u), 60):
            lines.append(u[i:i + 60])
    return {"text": "\n".join(lines) + "\n", "filename": "contigs.fasta",
            "message": f"{n} contigs, {sum(len(u) for u in unitigs if len(u) >= int(min_length))} bp"}


@T("assembly_metrics", "N50, L50 and contig size distribution", ASM, "table",
   [fa("contigs", GENOME, "Contig FASTA"), intin("bins", 8, "Length histogram bins", min=2),
    intin("min_length", 0, "Ignore contigs shorter than", min=0)],
   ex={"contigs": GENOME, "bins": 5}, up="quast / asmstats", tags=("assembly", "metrics", "N50"),
   summary="Assembly summary: sizes, N50/L50, GC and a contig length histogram.")
def assembly_metrics(contigs, bins=8, min_length=0):
    """Contig statistics."""
    rs = [r for r in _rows(contigs) if len(r["seq"]) >= int(min_length)]
    lens = sorted((len(r["seq"]) for r in rs), reverse=True)
    if not lens:
        return table([], "no contigs passed the length filter")
    total = sum(lens)
    rows = [{"metric": "contigs", "value": len(lens)}, {"metric": "total_length", "value": total},
            {"metric": "max_contig", "value": lens[0]}, {"metric": "min_contig", "value": lens[-1]},
            {"metric": "mean_length", "value": round(stats.mean(lens), 2)},
            {"metric": "median_length", "value": stats.median(lens)},
            {"metric": "N50", "value": seq.n50(lens)}, {"metric": "L50", "value": seq.l50(lens)},
            {"metric": "GC_percent", "value": round(100 * stats.mean([seq.gc_content(r["seq"]) for r in rs]), 3)},
            {"metric": "contigs_per_10kb", "value": round(10000 * len(lens) / total, 3)}]
    res = table(rows, f"N50 = {seq.n50(lens)} over {len(lens)} contigs")
    res["stats"]["length_histogram"] = _hist(lens, bins)
    return res


@T("kmer_spectrum", "Read k-mer frequency spectrum", ASM, "table",
   [fq("reads", READS_SINGLE, "Reads (FASTQ)"), intin("k", 17, "k-mer size", min=4),
    intin("max_count", 20, "Merge counts above", min=2), boolean("canonical", True, "Canonical k-mers")],
   ex={"reads": READS_SINGLE, "k": 15}, up="kmergenie / jellyfish hist",
   tags=("assembly", "k-mer spectrum"),
   summary="Histogram of k-mer multiplicities - the error peak and genomic peak used to pick k.")
def kmer_spectrum(reads, k=17, max_count=20, canonical=True):
    """K-mer count spectrum of a read set."""
    c: Counter = Counter()
    for r in _reads(reads):
        c.update(seq.kmer_counts(r.seq, int(k), canonical=bool(canonical)))
    hist: Counter = Counter()
    for km, n in c.items():
        hist[min(int(max_count), n)] += 1
    rows = [{"count": cc, "distinct_kmers": nn, "label": f">{max_count}" if cc == int(max_count) else str(cc)}
            for cc, nn in sorted(hist.items())]
    peak = max(hist.items(), key=lambda kv: kv[1])[0] if hist else 0
    return {"message": f"{len(c)} distinct {k}-mers, modal multiplicity {peak}", "rows": rows,
            "stats": {"total_kmer_instances": sum(c.values()), "distinct_kmers": len(c),
                      "singleton_kmers": hist.get(1, 0), "modal_count": peak}}


@T("greedy_overlap_assembly", "Greedy overlap-layout-consensus assembly", ASM, "table",
   [fq("reads", READS_SINGLE, "Reads (FASTQ)"), intin("min_overlap", 20, "Minimum overlap", min=6),
    intin("head", 100, "Rows", min=1)],
   ex={"reads": READS_SINGLE, "min_overlap": 25}, up="cap3 / phrap", tags=("assembly", "OLC"),
   summary="Repeatedly join the read pair with the largest exact overlap to build contigs.")
def greedy_overlap_assembly(reads, min_overlap=20, head=100):
    """Greedy OLC assembly."""
    rs = [seq.clean(r.seq) for r in _reads(reads)]
    min_ov = int(min_overlap)
    cands = []
    for i, a in enumerate(rs):
        for j, b in enumerate(rs):
            if i == j:
                continue
            for ov in range(min(len(a), len(b)), min_ov - 1, -1):
                if a[-ov:] == b[:ov]:
                    cands.append((ov, i, j))
                    break
    cands.sort(reverse=True)
    contigs = list(rs)
    used = set()
    joins = 0
    for ov, i, j in cands:
        if i in used or j in used:
            continue
        contigs[i] = contigs[i] + contigs[j][ov:]
        used.add(j)
        joins += 1
    out = [c for i, c in enumerate(contigs) if i not in used and c]
    lens = sorted((len(c) for c in out), reverse=True)
    rows = [{"contig": f"contig_{i + 1}", "length": len(c),
             "gc_percent": round(100 * seq.gc_content(c), 3)} for i, c in enumerate(out)]
    return {"message": f"{joins} joins, {len(out)} contigs (N50 {seq.n50(lens) if lens else 0})",
            "rows": rows[:int(head)],
            "stats": {"contigs": len(out), "N50": seq.n50(lens) if lens else 0,
                      "total_length": sum(lens), "overlap_joins": joins}}


@T("contig_nuc_content", "Nucleotide composition of each contig", ASM, "table",
   [fa("contigs", GENOME, "Contig FASTA"), intin("window", 0, "Window size (0 = whole contig)", min=0)],
   ex={"contigs": GENOME}, up="seqkit stats / emboss geece", tags=("assembly", "composition"),
   summary="Per-contig GC, N content, ambiguity and longest homopolymer.")
def contig_nuc_content(contigs, window=0):
    """Composition per contig."""
    rows = []
    for r in _rows(contigs):
        s = r["seq"].upper()
        w = int(window)
        seg = s[:w] if w > 0 else s
        hp = seq.longest_homopolymer(seg)
        rows.append({"contig": r["id"], "length": len(s),
                     "GC_percent": round(100 * seq.gc_content(s), 3),
                     "N_count": s.count("N"),
                     "lowercase_softmask": sum(1 for c in s if c.islower()),
                     "longest_homopolymer": hp.get("length", 0) if isinstance(hp, dict) else 0,
                     "A_percent": round(100 * s.count("A") / max(1, len(s)), 3),
                     "C_percent": round(100 * s.count("C") / max(1, len(s)), 3)})
    return table(rows, f"{len(rows)} contigs profiled")


@T("scaffold_with_pairs", "Link contigs with paired reads", ASM, "table",
   [fq("reads1", READS_1, "R1 FASTQ"), fq("reads2", READS_2, "R2 FASTQ"),
    fa("contigs", GENOME, "Contig FASTA"), intin("k", 21, "Seed length", min=8),
    intin("min_links", 1, "Minimum links to report", min=1)],
   ex={"reads1": READS_1, "reads2": READS_2, "contigs": GENOME, "k": 21},
   up="scaffolds / OPERA", tags=("assembly", "scaffolding"),
   summary="Count read pairs bridging contig pairs to infer scaffold adjacency.")
def scaffold_with_pairs(reads1, reads2, contigs, k=21, min_links=1):
    """Contig linking by paired-end evidence."""
    ct = _rows(contigs)
    name = {c["id"]: i for i, c in enumerate(ct)}
    links: Counter = Counter()
    pos_on: dict[tuple, list] = defaultdict(list)
    for a, b in zip(_map(contigs, reads1, k, 4), _map(contigs, reads2, k, 4)):
        if a["rname"] == "*" or b["rname"] == "*":
            continue
        key = tuple(sorted([a["rname"], b["rname"]]))
        links[key] += 1
        pos_on[key].append((a["pos"], b["pos"]))
    rows = []
    for (x, y), n in links.most_common():
        if n < int(min_links):
            continue
        ds = [abs(p - q) for p, q in pos_on[(x, y)]]
        rows.append({"contig_1": x, "contig_2": y, "links": n,
                    "same_contig": x == y,
                    "mean_mate_distance": round(stats.mean(ds), 2) if ds else 0,
                    "orientation": "+-" if x != y else "++"})
    return table(rows[:200], f"{len(rows)} linked contig pairs from {sum(links.values())} pairs")


@T("assembly_completeness", "Search conserved genes in an assembly", ASM, "table",
   [fa("contigs", GENOME, "Contig FASTA"), fa("genes", GENES_FA, "Expected gene FASTA"),
    intin("word_size", 11, "Word size", min=6), number("min_identity", 80.0, "Minimum identity %", min=0.0)],
   ex={"contigs": GENOME, "genes": GENES_FA, "word_size": 11}, up="BUSCO-lite",
   tags=("assembly", "completeness", "busco"),
   summary="BLAST each expected gene against the assembly and score the recovered fraction.")
def assembly_completeness(contigs, genes, word_size=11, min_identity=80.0):
    """Gene recovery (BUSCO-like) for an assembly."""
    hits = align.blast_lite(_rows(genes), _rows(contigs), word_size=int(word_size),
                           evalue=1.0, max_hits=20)
    found: dict[str, dict] = {}
    for h in align.best_hits(hits):
        if float(h["pident"]) >= float(min_identity):
            found[h["qaccver"]] = h
    genes_all = _rows(genes)
    rows = []
    for g in genes_all:
        h = found.get(g["id"], {})
        qlen = max(1, len(g["seq"]))
        try:
            cov = 100.0 * (int(h.get("qend", 0)) - int(h.get("qstart", 0)) + 1) / qlen
        except (TypeError, ValueError):
            cov = 0.0
        rows.append({"gene": g["id"], "length": qlen,
                    "status": "complete" if g["id"] in found else "missing",
                    "contig": h.get("saccver", ""), "identity_percent": h.get("pident", 0.0),
                    "coverage_percent": round(cov, 2)})
    n_ok = sum(1 for r in rows if r["status"] == "complete")
    return {"message": f"{n_ok}/{len(rows)} expected genes recovered",
            "rows": rows, "stats": {"complete_percent": round(100 * n_ok / max(1, len(rows)), 2),
                                    "missing": [r["gene"] for r in rows if r["status"] == "missing"]}}


@T("contig_read_coverage", "Map reads back onto contigs for coverage", ASM, "table",
   [fq("reads", READS_SINGLE, "Reads (FASTA/FASTQ)"), fa("contigs", GENOME, "Contig FASTA"),
    intin("k", 21, "Seed length", min=6), intin("bins", 10, "Coverage bins", min=2),
    intin("head", 100, "Rows", min=1)],
   ex={"reads": READS_SINGLE, "contigs": GENOME}, up="bwa + samtools depth",
   tags=("assembly", "coverage"),
   summary="Per-contig mapped-read counts, mean depth and zero-coverage fraction.")
def contig_read_coverage(reads, contigs, k=21, bins=10, head=100):
    """Back-map reads onto contigs."""
    ct = _rows(contigs)
    recs = _map(contigs, reads, k, 4)
    per: dict[str, list] = defaultdict(list)
    for r in recs:
        if r["rname"] != "*":
            per[r["rname"]].append(r)
    rows = []
    for c in ct:
        al = per.get(c["id"], [])
        tot = len(c["seq"]) or 1
        covered = set()
        for a in al:
            covered.update(range(max(0, a["pos"] - 1), min(tot, a["pos"] - 1 + len(a["seq"]))))
        rows.append({"contig": c["id"], "length": tot, "reads": len(al),
                    "mean_depth": round(len(al) * len(al[0]["seq"]) / tot, 3) if al else 0.0,
                    "covered_percent": round(100 * len(covered) / tot, 3),
                    "zero_coverage_bases": tot - len(covered)})
    res = table(sorted(rows, key=lambda r: -r["reads"])[:int(head)],
                f"{sum(r['reads'] for r in rows)} reads placed on {len(ct)} contigs")
    res["stats"]["coverage_histogram"] = _hist([r["mean_depth"] for r in rows], bins)
    return res


# ===========================================================================
# multiple sequence alignment
# ===========================================================================
@T("mafft_align", "Multiple sequence alignment (MAFFT)", MSA_S, "table",
   [anyfile("src", MSA, fmt="fasta", label="unaligned FASTA"),
    choice("algorithm", ["auto", "fft-ns-1", "linsi", "refine"], "auto", "Algorithm"),
    number("gap_open", -1.5, "Gap open penalty", max=0.0), number("gap_extend", -0.5, "Gap extend", max=0.0),
    intin("maxiterate", 0, "Refinement iterations", min=0)],
   ex={"src": MSA, "algorithm": "auto"}, up="mafft", tags=("MSA", "alignment"),
   summary="Center-star progressive MSA with optional profile refinement.")
def mafft_align(src, algorithm="auto", gap_open=-1.5, gap_extend=-0.5, maxiterate=0):
    """MAFFT-style alignment of a FASTA file."""
    rs = _aln_rows(src)
    if len(rs) < 2:
        return table([], "need at least two sequences")
    aln = (align.align_muscle(rs, refinement=int(maxiterate) or 1) if algorithm == "refine"
           else align.align_mafft(rs, algorithm=algorithm, gap_open=float(gap_open),
                                gap_extend=float(gap_extend), maxiterate=int(maxiterate)))
    rows = [{"id": a["id"], "aligned_sequence": a["seq"], "gaps": a["seq"].count("-"),
             "identity_to_reference": a.get("identities_to_reference", 0.0)} for a in aln]
    res = table(rows, f"{len(aln)} sequences, {len(aln[0]['seq'])} columns")
    res["stats"] = align.msa_stats(aln)
    return res


@T("muscle_align", "Multiple sequence alignment (MUSCLE)", MSA_S, "table",
   [anyfile("src", MSA, fmt="fasta", label="unaligned FASTA"), intin("refinement", 1, "Refinement rounds", min=0)],
   ex={"src": MSA, "refinement": 2}, up="muscle", tags=("MSA", "alignment"),
   summary="Progressive alignment with refinement cycles, as muscle -refine does.")
def muscle_align(src, refinement=1):
    """MUSCLE-style alignment."""
    aln = align.align_muscle(_aln_rows(src), refinement=int(refinement))
    rows = [{"id": a["id"], "aligned_sequence": a["seq"], "gaps": a["seq"].count("-")} for a in aln]
    res = table(rows, f"aligned {len(aln)} sequences")
    res["stats"] = align.msa_stats(aln)
    return res


@T("clustal_omega_align", "Multiple sequence alignment (Clustal Omega)", MSA_S, "table",
   [anyfile("src", MSA, fmt="fasta", label="unaligned FASTA"),
    choice("guide_tree", ["upgma", "nj"], "upgma", "Guide tree")],
   ex={"src": MSA}, up="clustalo", tags=("MSA", "alignment"),
   summary="Alignment guided by a distance tree, output in Clustal-style columns.")
def clustal_omega_align(src, guide_tree="upgma"):
    """Clustal Omega style alignment."""
    aln = align.clustal_omega(_aln_rows(src), guide_tree=guide_tree)
    rows = [{"id": a["id"], "aligned_sequence": a["seq"], "gaps": a["seq"].count("-")} for a in aln]
    res = table(rows, f"aligned {len(aln)} sequences")
    res["stats"] = align.msa_stats(aln)
    return res


@T("msa_stats", "Alignment column statistics", MSA_S, "table",
   [anyfile("src", MSA, fmt="fasta", label="aligned FASTA"), intin("gap_threshold", 50, "Gappy column threshold %",
                                                             min=1, max=100)],
   ex={"src": MSA}, up="aliview / msa summary", tags=("MSA", "statistics"),
   summary="Conserved, variable, gappy and gap-only columns plus mean pairwise identity.")
def msa_stats(src, gap_threshold=50):
    """Statistics of an alignment."""
    aln = _aln_rows(src)
    st = align.msa_stats(aln)
    n = st.get("sequences", 0)
    L = st.get("alignment_length", 0)
    thr = n * int(gap_threshold) / 100.0
    gappy = sum(1 for i in range(L) if sum(1 for a in aln if a["seq"][i] in "-.") > thr)
    st["gappy_columns"] = gappy
    st["identity_threshold_percent"] = int(gap_threshold)
    rows = [{"metric": k, "value": round(v, 4) if isinstance(v, float) else v} for k, v in st.items()]
    return table(rows, f"alignment of {n} sequences x {L} columns")


@T("msa_consensus", "Consensus sequence of an alignment", MSA_S, "table",
   [anyfile("src", MSA, fmt="fasta", label="aligned FASTA"),
    choice("mode", ["majority", "IUPAC", "first"], "majority", "Consensus mode"),
    number("threshold", 0.5, "Majority fraction", min=0.5, max=1.0)],
   ex={"src": MSA, "mode": "IUPAC"}, up="emboss cons / samtools consensus",
   tags=("MSA", "consensus"),
   summary="Collapse the alignment to one sequence using majority votes or IUPAC codes.")
def msa_consensus(src, mode="majority", threshold=0.5):
    """Consensus of an MSA."""
    aln = _aln_rows(src)
    cons = align.msa_consensus(aln, mode=mode)
    freqs = align.msa_column_frequencies(aln)
    rows = []
    for i, ch in enumerate(cons):
        d = freqs[i] if i < len(freqs) else {}
        rows.append({"position": i + 1, "consensus": ch,
                    "frequencies": "|".join(f"{k}:{v}" for k, v in sorted(d.items()) if v),
                    "majority_fraction": round(max(d.values()) / max(1, len(aln)), 4) if d else 0.0})
    res = table(rows, f"{mode} consensus, {len(cons)} columns")
    res["stats"] = {"consensus_sequence": cons, "threshold": float(threshold)}
    return res


@T("msa_column_frequencies", "Per-column residue frequencies", MSA_S, "table",
   [anyfile("src", MSA, fmt="fasta", label="aligned FASTA"), intin("start", 1, "First column", min=1),
    intin("end", 0, "Last column (0 = all)", min=0)],
   ex={"src": MSA, "start": 1, "end": 20}, up="seqkit stats on alignment",
   tags=("MSA", "frequencies"),
   summary="Frequency table of every residue at every alignment column.")
def msa_column_frequencies(src, start=1, end=0):
    """Column frequency matrix."""
    aln = _aln_rows(src)
    freqs = align.msa_column_frequencies(aln)
    lo = max(1, int(start))
    hi = int(end) or len(freqs)
    rows = []
    for i in range(lo - 1, min(hi, len(freqs))):
        d = freqs[i]
        rows.append({"column": i + 1,
                    "counts": "|".join(f"{k}:{v}" for k, v in sorted(d.items()) if v),
                    "most_common": max(d, key=lambda k: d[k]) if d else ".",
                    "n_distinct": sum(1 for v in d.values() if v)})
    return table(rows, f"columns {lo}-{min(hi, len(freqs))} of {len(freqs)}")


@T("msa_trim_columns", "Trim gappy alignment columns (trimAl)", MSA_S, "table",
   [anyfile("src", MSA, fmt="fasta", label="aligned FASTA"),
    choice("mode", ["gappyout", "strict", "automated1"], "gappyout", "Trim mode"),
    number("threshold", 0.5, "Gap fraction threshold", min=0.0, max=1.0)],
   ex={"src": MSA, "mode": "gappyout", "threshold": 0.4}, up="trimal",
   tags=("MSA", "trimming", "phylogenetics"),
   summary="Remove poorly aligned columns before building a tree.")
def msa_trim_columns(src, mode="gappyout", threshold=0.5):
    """trimAl-style column trimming."""
    aln = _aln_rows(src)
    trimmed = align.msa_trim(aln, mode=mode, threshold=float(threshold))
    before, after = align.msa_stats(aln), align.msa_stats(trimmed)
    rows = [{"id": a["id"], "aligned_sequence": a["seq"], "gaps": a["seq"].count("-")} for a in trimmed]
    res = table(rows, f"{before['alignment_length']} -> {after['alignment_length']} columns kept")
    res["stats"] = {"before": before, "after": after,
                   "columns_removed": before["alignment_length"] - after["alignment_length"]}
    return res


@T("msa_variable_sites", "Variable columns of an alignment", MSA_S, "table",
   [anyfile("src", MSA, fmt="fasta", label="aligned FASTA"), intin("head", 300, "Rows", min=1),
    boolean("informative_only", False, "Only parsimony-informative sites")],
   ex={"src": MSA}, up="variable sites / seaview", tags=("MSA", "polymorphism"),
   summary="Columns that vary between sequences, with entropy and taxon pattern.")
def msa_variable_sites(src, head=300, informative_only=False):
    """Variable site table."""
    aln = _aln_rows(src)
    vs = align.variable_sites(aln)
    rows = []
    for v in vs:
        pat = str(v.get("pattern", ""))
        states = sorted(set(pat.replace("-", "")))
        if informative_only and len(states) < 2:
            continue
        rows.append({"column": v.get("column"), "entropy": v.get("entropy"),
                    "states": v.get("states", len(states)), "pattern": pat,
                    "consensus": v.get("consensus", "")})
    return table(rows[:int(head)], f"{len(rows)} variable columns")


@T("msa_entropy_profile", "Conservation profile of an alignment", MSA_S, "figure",
   [anyfile("src", MSA, fmt="fasta", label="aligned FASTA"), choice("style", ["line", "bar"], "line", "Plot type"),
    number("threshold", 1.0, "Highlight columns with entropy >=", min=0.0)],
   ex={"src": MSA, "style": "line"}, up="plot of per-column entropy", tags=("MSA", "conservation", "plot"),
   summary="Shannon entropy per alignment column: conserved core versus variable loops.")
def msa_entropy_profile(src, style="line", threshold=1.0):
    """Entropy profile across an alignment."""
    aln = _aln_rows(src)
    ent = align.column_entropy(aln)
    ys = [float(e) for e in ent]
    if not ys:
        return plot.line({"entropy": [0.0]}, title="empty alignment")
    if style == "bar":
        return plot.bar([str(i + 1) for i in range(len(ys))], ys, xlabel="column",
                       ylabel="entropy", title="Per-column entropy")
    return plot.line({"entropy": ys, "threshold": [float(threshold)] * len(ys)},
                     x=list(range(1, len(ys) + 1)), xlabel="alignment column",
                     ylabel="Shannon entropy", title="Conservation profile")


@T("msa_format_convert", "Convert alignment to Clustal / Stockholm / PHYLIP", MSA_S, "file",
   [anyfile("src", MSA, fmt="fasta", label="aligned FASTA"),
    choice("format", ["clustal", "stockholm", "phylip", "fasta"], "clustal", "Output format")],
   ex={"src": MSA, "format": "clustal"}, up="alignment conversion", tags=("MSA", "format", "conversion"),
   summary="Re-serialise an alignment in the format expected by phylogenetics or HMM tools.")
def msa_format_convert(src, format="clustal"):
    """Alignment format conversion."""
    aln = _aln_rows(src)
    ext = {"clustal": "aln", "stockholm": "sto", "phylip": "phy", "fasta": "fasta"}[format]
    body = {"clustal": lambda: align.write_clustal(aln), "stockholm": lambda: align.write_stockholm(aln),
            "phylip": lambda: align.to_phylip(aln), "fasta": lambda: align.alignment_to_fasta(aln)}[format]()
    return {"text": body, "filename": f"alignment.{ext}",
            "message": f"{len(aln)} sequences written as {format}"}


@T("msa_back_translate", "Protein alignment to codon alignment", MSA_S, "table",
   [anyfile("protein_alignment", PROTEINS, fmt="fasta", label="aligned proteins"),
    fa("nucleotides", GENES_FA, "Nucleotide CDS FASTA")],
   ex={"protein_alignment": MSA, "nucleotides": GENES_FA}, up="pal2nal / rev-trans",
   tags=("MSA", "codon", "back-translation"),
   summary="Re-insert codons into a protein alignment to obtain a codon-based alignment.")
def msa_back_translate(protein_alignment, nucleotides):
    """pal2nal-style back translation."""
    aln = align.back_translate_alignment(_aln_rows(protein_alignment), _rows(nucleotides))
    rows = [{"id": a.get("id", f"seq_{i + 1}"), "length": len(a["seq"]),
             "codons": len(a["seq"]) // 3, "n_gaps": a["seq"].count("-"),
             "codon_alignment": a["seq"][:90] + ("..." if len(a["seq"]) > 90 else "")}
            for i, a in enumerate(aln)]
    return table(rows, f"{len(aln)} codon alignments rebuilt")


# ===========================================================================
# pairwise alignment
# ===========================================================================
@T("pairwise_align", "Pairwise alignment (Needleman-Wunsch / Smith-Waterman)", MSA_S, "table",
   [anyfile("seq_a", GENES_FA, fmt="fasta", label="Sequence A"), anyfile("seq_b", GENES_FA, fmt="fasta", label="Sequence B"),
    choice("mode", ["global", "local"], "global", "Alignment mode"),
    choice("matrix", ["", "BLOSUM62", "PAM250"], "", "Substitution matrix (empty = match/mismatch)"),
    number("gap_open", -1.5, "Gap open", max=0.0), number("gap_extend", -0.5, "Gap extend", max=0.0),
    intin("match", 2, "Match score", min=1), intin("mismatch", -1, "Mismatch score", max=0)],
   ex={"seq_a": GENES_FA, "seq_b": GENES_FA, "mode": "local"}, up="needle / water (EMBOSS)",
   tags=("pairwise", "alignment"),
   summary="Dynamic-programming alignment of two sequences with affine gaps and a matrix.")
def pairwise_align(seq_a, seq_b, mode="global", matrix="", gap_open=-1.5, gap_extend=-0.5,
                   match=2, mismatch=-1):
    """DP pairwise alignment."""
    a = _first_seq(seq_a)
    b = _first_seq(seq_b)
    if not a or not b:
        return table([], "could not read two sequences from the inputs")
    res = align.align_pairwise(a, b, mode=mode, matrix=matrix, gap_open=float(gap_open),
                              gap_extend=float(gap_extend), match=int(match), mismatch=int(mismatch))
    rows = [{"metric": k, "value": round(v, 4) if isinstance(v, float) else v}
            for k, v in res.items() if not isinstance(v, str)]
    blocks = "\n".join(f"{k:<16}{v}" for k, v in res.items() if isinstance(v, str))
    out = table(rows, f"score {res.get('score', 0)}, identity {res.get('identity_percent', 0)}%")
    out["stats"] = {"alignment_view": blocks[:4000]}
    return out


def _first_src_text(src) -> str:
    return io.as_text(src)


def _first_seq(src) -> str:
    """Sequence text: first FASTA record, or the raw pasted letters."""
    txt = _first_src_text(src)
    rs = io.parse_fasta(txt)
    if rs:
        return rs[0].seq
    return seq.clean("\n".join(ln for ln in txt.splitlines() if not ln.startswith(">")))


@T("pairwise_identity_matrix", "Pairwise identity / distance matrix", MSA_S, "table",
   [anyfile("src", PROTEINS, fmt="fasta", label="FASTA sequences"),
    choice("model", ["identity", "distance", "blosum"], "identity", "Distance model"),
    boolean("heatmap", True, "Also draw a heatmap")],
   ex={"src": PROTEINS, "model": "identity"}, up="peptide-database / emboss proteity",
   tags=("distance matrix", "pairwise"),
   summary="N x N matrix of pairwise identity or substitution-corrected distance.")
def pairwise_identity_matrix(src, model="identity", heatmap=True):
    """Pairwise matrix over a FASTA file."""
    rs = _aln_rows(src)
    names = [r["id"] for r in rs]
    mat = align.distance_matrix([r["seq"] for r in rs],
                               model="identity" if model == "identity" else "blosum",
                               matrix="BLOSUM62" if model == "blosum" else "")
    rows = [{"sequence": names[i], **{names[j]: round(float(mat[i][j]), 4) for j in range(len(names))}}
            for i in range(len(names))]
    res = table(rows, f"{len(names)} x {len(names)} {model} matrix")
    if heatmap and len(names) <= 40:
        res["stats"]["figure"] = plot.heatmap([[float(v) for v in r] for r in mat], row_labels=names,
                                            col_labels=names, title=f"{model} matrix")
    return res


@T("edit_distance_report", "Edit distance and Hamming distance between sequences", MSA_S, "table",
   [anyfile("src", MSA, fmt="fasta", label="FASTA sequences"),
    intin("head", 200, "Rows", min=1), boolean("upper_only", True, "One row per pair")],
   ex={"src": MSA}, up="emboss distmat / edlib", tags=("distance", "pairwise"),
   summary="Levenshtein and Hamming distances for every pair of sequences in a file.")
def edit_distance_report(src, head=200, upper_only=True):
    """All-pairs edit distances."""
    rs = _aln_rows(src)
    rows = []
    for i in range(len(rs)):
        for j in range(len(rs)):
            if i == j or (upper_only and j < i):
                continue
            ed = align.edit_distance(rs[i]["seq"], rs[j]["seq"])
            hd = align.hamming(rs[i]["seq"], rs[j]["seq"])
            rows.append({"sequence_1": rs[i]["id"], "sequence_2": rs[j]["id"],
                        "edit_distance": ed.get("distance", 0),
                        "matches": ed.get("matches", 0), "insertions": ed.get("insertions", 0),
                        "deletions": ed.get("deletions", 0),
                        "substitutions": ed.get("substitutions", 0),
                        "hamming": hd.get("distance", hd) if isinstance(hd, dict) else hd,
                        "identity_percent": ed.get("identity_percent", 0.0)})
    rows.sort(key=lambda r: r["edit_distance"])
    return table(rows[:int(head)], f"{len(rows)} pairwise distances")


@T("pam_kimura_distances", "Substitution-corrected distances", MSA_S, "table",
   [anyfile("src", MSA, fmt="fasta", label="aligned FASTA"),
    choice("model", ["pam", "kimura", "uncorrected"], "pam", "Substitution model"),
    number("max_distance", 3.0, "Cap distance at", min=0.1)],
   ex={"src": MSA, "model": "kimura"}, up="phylip proml / distmat",
   tags=("evolution", "distance"),
   summary="Pairwise substitution distances corrected for multiple hits.")
def pam_kimura_distances(src, model="pam", max_distance=3.0):
    """Corrected distances between aligned sequences."""
    aln = _aln_rows(src)
    rows = []
    for i in range(len(aln)):
        for j in range(i + 1, len(aln)):
            a, b = aln[i]["seq"], aln[j]["seq"]
            pairs = [(x, y) for x, y in zip(a, b) if x not in "-." and y not in "-."]
            n = len(pairs) or 1
            diffs = sum(1 for x, y in pairs if x != y)
            p = diffs / n
            if model == "pam":
                d = float(align.pam_distance(a, b))
            elif model == "kimura":
                d = -0.75 * math.log(max(1e-9, 1 - 1.333 * p)) if p < 0.75 else float(max_distance)
            else:
                d = p
            rows.append({"sequence_1": aln[i]["id"], "sequence_2": aln[j]["id"], "model": model,
                        "compared_positions": n, "differences": diffs,
                        "p_distance": round(p, 5), "distance": round(min(d, float(max_distance)), 5)})
    return table(rows, f"{len(rows)} pairwise {model} distances")


# ===========================================================================
# BLAST-style search
# ===========================================================================
def _hit_table(hits, head=500):
    rows = []
    for h in hits[:int(head)]:
        rows.append({"qaccver": h.get("qaccver", ""), "saccver": h.get("saccver", ""),
                    "pident": round(float(h.get("pident", 0.0)), 3), "length": h.get("length", 0),
                    "mismatch": h.get("mismatch", 0), "gapopen": h.get("gapopen", 0),
                    "qstart": h.get("qstart", 0), "qend": h.get("qend", 0),
                    "sstart": h.get("sstart", 0), "send": h.get("send", 0),
                    "evalue": float(h.get("evalue", 0.0)), "bitscore": round(float(h.get("bitscore", 0.0)), 2)})
    return rows


@T("blastn_search", "BLASTN: nucleotide query versus nucleotide database", BLAST, "table",
   [fa("query", GENES_FA, "Query FASTA"), fa("database", GENOME, "Database FASTA"),
    intin("word_size", 11, "Word size", min=4), number("evalue", 1e-5, "E-value threshold", min=0.0),
    intin("max_target_seqs", 10, "Max hits per query", min=1), boolean("dust", True, "DUST low-complexity mask"),
    intin("gap_open", -11, "Gap open penalty", max=0)],
   ex={"query": GENES_FA, "database": GENOME, "word_size": 11, "evalue": 1e-5},
   up="blastn", tags=("BLAST", "nucleotide", "search"),
   summary="Word-hit BLASTN search with extension, bit scores and E-values.")
def blastn_search(query, database, word_size=11, evalue=1e-5, max_target_seqs=10, dust=True, gap_open=-11):
    """blastn against a FASTA database."""
    hits = align.blast_lite(_rows(query), _rows(database), word_size=int(word_size),
                           evalue=float(evalue), gap_open=float(gap_open), max_hits=int(max_target_seqs),
                           dust=bool(dust))
    res = table(_hit_table(hits), f"{len(hits)} hits for {len(_rows(query))} queries")
    res["stats"] = {"queries_with_hits": len({h.get('qaccver') for h in hits}),
                    "best_bitscore": round(max((float(h.get("bitscore", 0)) for h in hits), default=0.0), 2)}
    return res


@T("blastp_search", "BLASTP: protein query versus protein database", BLAST, "table",
   [fa("query", PROTEINS, "Query proteins"), fa("database", REF_PROTEINS, "Protein database"),
    choice("matrix", ["BLOSUM62", "BLOSUM45", "PAM250"], "BLOSUM62", "Scoring matrix"),
    number("evalue", 10.0, "E-value threshold", min=0.0), intin("word_size", 3, "Word size", min=2),
    intin("max_target_seqs", 10, "Max hits per query", min=1)],
   ex={"query": PROTEINS, "database": REF_PROTEINS, "word_size": 3, "evalue": 10},
   up="blastp", tags=("BLAST", "protein", "search"),
   summary="Protein-protein search scored with BLOSUM/PAM matrices.")
def blastp_search(query, database, matrix="BLOSUM62", evalue=10.0, word_size=3, max_target_seqs=10):
    """blastp-style protein search."""
    hits = align.blast_lite(_rows(query), _rows(database), word_size=int(word_size),
                           evalue=float(evalue), matrix=matrix, scoring="blastp",
                           gap_open=-11.0, max_hits=int(max_target_seqs))
    return table(_hit_table(hits), f"{len(hits)} protein hits")


@T("blastx_search", "BLASTX: nucleotide query in six frames", BLAST, "table",
   [fa("query", GENES_FA, "Nucleotide query"), fa("database", REF_PROTEINS, "Protein database"),
    intin("word_size", 3, "Word size", min=2), number("evalue", 10.0, "E-value", min=0.0),
    choice("matrix", ["BLOSUM62", "PAM250"], "BLOSUM62", "Matrix")],
   ex={"query": GENES_FA, "database": REF_PROTEINS, "word_size": 3, "evalue": 10},
   up="blastx", tags=("BLAST", "translation", "search"),
   summary="Translate the query in six frames and search a protein database.")
def blastx_search(query, database, word_size=3, evalue=10.0, matrix="BLOSUM62"):
    """Six-frame blastx."""
    frames = []
    for r in _rows(query):
        for fr in range(3):
            fwd = seq.translate(r["seq"][fr:], frame=0)
            rc = seq.reverse_complement(r["seq"])[fr:]
            rows = [{"id": f"{r['id']}|frame{fr + 1}", "seq": fwd},
                    {"id": f"{r['id']}|revframe{fr + 1}", "seq": seq.translate(rc, frame=0)}]
            frames.extend([x for x in rows if x["seq"]])
    hits = align.blast_lite(frames, _rows(database), word_size=int(word_size), evalue=float(evalue),
                           matrix=matrix, scoring="blastp", max_hits=10)
    return table(_hit_table(hits), f"{len(hits)} hits from {len(frames)} translated frames")


@T("tblastn_search", "TBLASTN: protein query versus translated database", BLAST, "table",
   [fa("query", PROTEINS, "Protein query"), fa("database", GENOME, "Nucleotide database"),
    intin("word_size", 3, "Word size", min=2), number("evalue", 10.0, "E-value", min=0.0)],
   ex={"query": PROTEINS, "database": GENOME, "word_size": 3, "evalue": 10},
   up="tblastn", tags=("BLAST", "translation", "search"),
   summary="Search a protein query against the six translated frames of a nucleotide database.")
def tblastn_search(query, database, word_size=3, evalue=10.0):
    """Protein query vs translated subject."""
    subj = []
    for r in _rows(database):
        for fr in range(3):
            for tag, s in ((fr + 1, r["seq"][fr:]), (-fr - 1, seq.reverse_complement(r["seq"])[fr:])):
                prot = seq.translate(s, frame=0)
                if prot:
                    subj.append({"id": f"{r['id']}|frame{tag}", "seq": prot})
    hits = align.blast_lite(_rows(query), subj, word_size=int(word_size), evalue=float(evalue),
                           matrix="BLOSUM62", scoring="blastp", max_hits=10)
    return table(_hit_table(hits), f"{len(hits)} tblastn hits against {len(subj)} translated subjects")


@T("blast_best_hit_per_query", "Keep the best hit of every query", BLAST, "table",
   [fa("query", GENES_FA, "Query FASTA"), fa("database", GENOME, "Database FASTA"),
    choice("by", ["bitscore", "evalue"], "bitscore", "Rank by"),
    intin("word_size", 11, "Word size", min=4), number("evalue", 1e-3, "E-value threshold", min=0.0)],
   ex={"query": GENES_FA, "database": GENOME, "word_size": 11}, up="blast + best-hit filter",
   tags=("BLAST", "filtering"),
   summary="Reduce all-vs-all hits to the single top-scoring alignment per query.")
def blast_best_hit_per_query(query, database, by="bitscore", word_size=11, evalue=1e-3):
    """Best hit per query."""
    hits = align.blast_lite(_rows(query), _rows(database), word_size=int(word_size),
                           evalue=float(evalue), max_hits=50)
    best = align.best_hits(hits)
    rows = _hit_table(best, head=1000)
    if by == "evalue":
        rows.sort(key=lambda r: r["evalue"])
    return table(rows, f"{len(rows)} queries with a hit (of {len(_rows(query))})")


@T("blast_identity_histogram", "Identity distribution of BLAST hits", BLAST, "figure",
   [fa("query", GENES_FA, "Query FASTA"), fa("database", GENOME, "Database FASTA"),
    intin("bin_size", 5, "Identity bin width", min=1), intin("word_size", 11, "Word size", min=4)],
   ex={"query": GENES_FA, "database": GENOME, "bin_size": 5}, up="BLAST histogram",
   tags=("BLAST", "histogram", "plot"),
   summary="Histogram of percent identity across the reported hits.")
def blast_identity_histogram(query, database, bin_size=5, word_size=11):
    """Hit identity histogram."""
    hits = align.blast_lite(_rows(query), _rows(database), word_size=int(word_size), evalue=1.0,
                           max_hits=200)
    hist = align.similarity_histogram(hits, bin_size=int(bin_size))
    return plot.bar([str(k) for k in sorted(hist)], [hist[k] for k in sorted(hist)],
                   xlabel="percent identity bin", ylabel="hits",
                   title=f"BLAST hit identity ({len(hits)} hits)")


@T("blast_make_database", "Prepare a BLAST database", BLAST, "file",
   [fa("sequences", GENOME, "Sequences to index"), intin("word_size", 11, "Word size", min=4),
    boolean("protein", False, "Protein database")],
   ex={"sequences": GENES_FA, "word_size": 11}, up="makeblastdb",
   tags=("BLAST", "database", "index"),
   summary="Write a database descriptor with lengths, alphabets and k-mer index size.")
def blast_make_database(sequences, word_size=11, protein=False):
    """makeblastdb-style metadata."""
    rs = _rows(sequences)
    lens = [len(r["seq"]) for r in rs]
    idx = align.build_kmer_index(rs, int(word_size))
    lines = [">Dbtype", ("Prot" if protein else "Nucl"), ">Db_total_length", str(sum(lens)),
             ">Seq-len", ">Query-blen"]
    for r in rs:
        lines += [str(len(r["seq"])), str(len(r["seq"]))]
    lines += [">Taxid_list"] + ["0"] * len(rs)
    body = "\n".join(lines) + f"\n\n# {len(rs)} sequences, {len(idx)} distinct {word_size}-mer index entries\n"
    return {"text": body, "filename": "blastdb.pnd",
            "message": f"database of {len(rs)} sequences ({sum(lens)} bp, {len(idx)} index keys)"}


@T("blast_tabular_normalise", "Normalise a BLAST tabular output", BLAST, "table",
   [tbl("src", COUNTS, "BLAST output (12 columns, tab-separated)"),
    intin("columns", 12, "Expected columns", min=6), number("min_pident", 0.0, "Minimum identity %", min=0.0),
    number("max_evalue", 1e100, "Maximum E-value", min=0.0)],
   ex={"src": COUNTS, "columns": 12, "min_pident": 0.0}, up="blast -outfmt 6",
   tags=("BLAST", "tabular", "transformation"),
   summary="Rename columns of a -outfmt 6 table and filter by identity and E-value.")
def blast_tabular_normalise(src, columns=12, min_pident=0.0, max_evalue=1e100):
    """Column-name a raw BLAST table and filter it."""
    names = ["qaccver", "saccver", "pident", "length", "mismatch", "gapopen", "qstart", "qend",
             "sstart", "send", "evalue", "bitscore", "slens", "qcovs", "stitle"]
    n = int(columns)
    df = io.read_table(src, sep=None, header=None, names=names[:n])
    rows = df.to_dict("records") if n >= 12 else []
    out = []
    for r in rows:
        try:
            pid = float(r.get("pident") or 0.0)
            ev = float(r.get("evalue") or 0.0)
        except (TypeError, ValueError):
            continue
        if pid < float(min_pident) or ev > float(max_evalue):
            continue
        out.append({k: v for k, v in r.items()})
    return table(out[:500], f"{len(out)}/{len(rows)} alignments retained")


@T("blast_lca_taxonomy", "Assign lowest common ancestor taxonomy", CONTAM, "table",
   [fa("query", GENES_FA, "Query FASTA"), fa("database", GENOME, "Database FASTA"),
    bigtext("taxmap", "", "Taxonomy map: accession<TAB>lineage"), intin("word_size", 11, "Word size", min=4),
    intin("top_hits", 5, "Hits used per query", min=1)],
   ex={"query": GENES_FA, "database": GENOME, "word_size": 11,
       "taxmap": "chrV\tk__root;p__Chlorophyta;c__Ulvophyceae\nchrM\tk__root;p__Magnoliophyta;c__Liliopsida"},
   up="MEGAN / kraken LCA", tags=("taxonomy", "LCA", "BLAST"),
   summary="Classify each query by the lowest common ancestor of its best BLAST hits.")
def blast_lca_taxonomy(query, database, taxmap="", word_size=11, top_hits=5):
    """LCA assignment from BLAST hits."""
    db = _rows(database)
    hits = align.blast_lite(_rows(query), db, word_size=int(word_size), evalue=10.0,
                           max_hits=int(top_hits))
    tax: dict[str, list] = {}
    for ln in (taxmap or "").splitlines():
        if not ln.strip():
            continue
        parts = ln.split("\t") if "\t" in ln else ln.split(None, 1)
        if len(parts) >= 2:
            tax[parts[0].strip()] = [x.strip() for x in parts[1].replace(";", "\t").split("\t") if x.strip()]
    if not tax:
        tax = {d["id"]: [d["id"]] for d in db}
    rows = align.lca_of_hits(hits, tax)
    out = []
    for r in rows:
        out.append({"query": r.get("read", r.get("qaccver", "")), "lca": r.get("lca", ""),
                   "rank_depth": r.get("rank", 0), "n_hits": r.get("n_hits", 0)})
    return table(out, f"{len(out)} queries classified")


@T("contamination_filter", "Remove contaminant / host sequences", CONTAM, "table",
   [fa("reads", GENES_FA, "Sequences to screen"), fa("host", GENOME, "Host or contaminant database"),
    number("min_identity", 90.0, "Identity considered contaminant", min=0.0),
    intin("word_size", 11, "Word size", min=4), choice("keep", ["clean", "hits"], "clean", "Keep")],
   ex={"reads": GENES_FA, "host": GENOME, "min_identity": 80}, up="deconseq / blastn screen",
   tags=("contamination", "filtering"),
   summary="Flag sequences with a high-identity hit in the contaminant database and filter them out.")
def contamination_filter(reads, host, min_identity=90.0, word_size=11, keep="clean"):
    """Best-hit contamination screen."""
    hits = align.blast_lite(_rows(reads), _rows(host), word_size=int(word_size), evalue=1.0,
                           max_hits=5)
    best: dict[str, dict] = {}
    for h in hits:
        if float(h.get("pident", 0)) < float(min_identity):
            continue
        cur = best.get(h["qaccver"])
        if cur is None or float(h["bitscore"]) > float(cur["bitscore"]):
            best[h["qaccver"]] = h
    rows = []
    for r in _rows(reads):
        h = best.get(r["id"])
        flagged = h is not None
        rows.append({"query": r["id"], "length": len(r["seq"]), "classified": flagged,
                    "subject": (h or {}).get("saccver", ""),
                    "identity_percent": (h or {}).get("pident", 0.0),
                    "evalue": (h or {}).get("evalue", 0.0),
                    "action": ("removed" if flagged else "kept") if keep == "clean"
                    else ("kept" if flagged else "removed")})
    n_keep = sum(1 for r in rows if r["action"] == "kept")
    return table(rows, f"{n_keep}/{len(rows)} sequences kept after {keep}-oriented filtering")


@T("dereplicate_sequences", "Collapse identical sequences", CONTAM, "table",
   [fa("src", GENES_FA, "FASTA sequences"), choice("count_by", ["abundance", "first"], "abundance",
                                                  "Representative choice"),
    intin("min_size", 1, "Minimum cluster size", min=1)],
   ex={"src": GENES_FA}, up="vsearch --derep_fulllength", tags=("dereplication", "deduplication"),
   summary="Group sequences by checksum and keep one representative with an abundance count.")
def dereplicate_sequences(src, count_by="abundance", min_size=1):
    """Dereplication by sequence checksum."""
    groups: dict[str, list] = defaultdict(list)
    for r in _rows(src):
        groups[seq.sequence_checksum(r["seq"]) if hasattr(seq, "sequence_checksum")
               else str(hash(r["seq"]) & 0xffffffff)].append(r)
    rows = []
    for key, members in groups.items():
        if len(members) < int(min_size):
            continue
        rep = max(members, key=lambda m: len(m["seq"])) if count_by == "abundance" else members[0]
        rows.append({"sequence_id": rep["id"], "abundance": len(members), "size": len(rep["seq"]),
                    "checksum": key, "members": ",".join(m["id"] for m in members[:10])})
    rows.sort(key=lambda r: -r["abundance"])
    return table(rows, f"{len(groups)} unique sequences from {len(_rows(src))} records")


# ===========================================================================
# coordinate liftover
# ===========================================================================
@T("lift_over_chain", "Lift BED intervals through a chain file", LIFT, "file",
   [bed("src", REGIONS_BED, "BED intervals"), tbl("chain", "", "Chain blocks (tabular)"),
    choice("mode", ["smart", "start", "end"], "smart", "Anchor point"),
    choice("fail", ["drop", "keep"], "keep", "Unmapped intervals")],
   ex={"src": REGIONS_BED, "chain": "chrV\t0\t3000\tchrA\t100\t3100\nchrV\t3000\t6000\tchrA\t3100\t6100\n",
       "mode": "smart"}, up="UCSC liftOver", tags=("lift over", "coordinates", "assembly"),
   summary="Re-map interval coordinates between assemblies using aligned chain blocks.")
def lift_over_chain(src, chain, mode="smart", fail="keep"):
    """liftOver over a simple block chain."""
    ivs = io.parse_bed(io.as_text(src))
    blocks = []
    for ln in io.as_text(chain).splitlines():
        if not ln.strip() or ln.startswith("#"):
            continue
        parts = ln.replace(",", "\t").split()
        if len(parts) < 6:
            continue
        try:
            blocks.append({"s_chrom": parts[0], "s_start": int(float(parts[1])), "s_end": int(float(parts[2])),
                          "t_chrom": parts[3], "t_start": int(float(parts[4])), "t_end": int(float(parts[5]))})
        except ValueError:
            continue
    out: list[io.Interval] = []
    unmapped = 0
    for iv in ivs:
        anchor = iv.end - 1 if mode == "end" else iv.start
        hit = next((b for b in blocks if b["s_chrom"] == iv.chrom and b["s_start"] <= anchor < b["s_end"]),
                   None)
        if hit is None:
            unmapped += 1
            if fail == "keep":
                out.append(io.Interval(iv.chrom, iv.start, iv.end, f"{iv.name or '.'}_unmapped",
                                       iv.score if isinstance(iv.score, (int, float)) and
                                       not math.isnan(iv.score) else 0, iv.strand or "."))
            continue
        shift = iv.start - hit["s_start"]
        if mode == "end":
            shift = iv.end - hit["s_end"]
        ns = max(hit["t_start"], hit["t_start"] + shift)
        ne = min(hit["t_end"], ns + (iv.end - iv.start))
        if ne <= ns:
            ne = min(hit["t_end"], ns + 1)
        out.append(io.Interval(hit["t_chrom"], ns, ne, iv.name or ".",
                               iv.score if isinstance(iv.score, (int, float)) and
                               not math.isnan(iv.score) else 0, iv.strand or "."))
    body = io.write_bed(out, bed12=False)
    return {"text": body, "filename": "lifted.bed",
            "message": f"{len(out) - unmapped}/{len(ivs)} intervals lifted through {len(blocks)} "
                       f"chain blocks ({unmapped} unmapped)"}
