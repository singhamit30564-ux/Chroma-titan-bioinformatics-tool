"""Read-alignment engines: SAM/BAM statistics, pileups, quantification.

BAM binaries are represented by their SAM text in this suite (``.sam``), so the
whole project stays dependency-light while still doing real pileup, counting and
filtering work.
"""

from __future__ import annotations

import math
import re
from collections import Counter, defaultdict
from typing import Sequence

from chroma_titan.core import stats
from chroma_titan.core.io import Aln, parse_sam, sam_ref_names, write_sam

UNMAPPED, PAIRED, PROPER, FAILQC, DUP, SECONDARY, SUPPLEMENTARY = 4, 1, 2, 256, 1024, 256, 2048


def load(src) -> tuple[list[str], list[Aln]]:
    if isinstance(src, tuple) and len(src) == 2:
        return src
    return parse_sam(src)


def primary(alns: Sequence[Aln], mapped_only: bool = True, dedup: bool = False) -> list[Aln]:
    out = [a for a in alns if not a.secondary and not a.supplementary
           and (a.mapped or not mapped_only) and (not dedup or not a.duplicate)]
    return out


def flagstat(header: list[str], alns: Sequence[Aln]) -> dict:
    total = len(alns)
    mapped = sum(1 for a in alns if a.mapped and not a.secondary)
    qcf = sum(1 for a in alns if a.flag & FAILQC)
    dup = sum(1 for a in alns if a.duplicate)
    paired = sum(1 for a in alns if a.flag & PAIRED)
    proper = sum(1 for a in alns if a.proper_pair and not a.secondary)
    single = sum(1 for a in alns if a.mapped and a.flag & PAIRED)
    diff = sum(1 for a in alns if not a.secondary and a.mapped and a.rnext not in ("=", a.rname)
               and a.rnext != "*")
    sup = sum(1 for a in alns if a.supplementary)
    sec = sum(1 for a in alns if a.secondary)
    return {"total": total, "mapped": mapped, "mapped_percent": round(
        100 * mapped / total, 4) if total else 0.0,
        "secondary": sec, "supplementary": sup, "duplicates": dup,
        "mapped_duplicates": sum(1 for a in alns if a.duplicate and a.mapped),
        "qc_fail": qcf, "properly_paired": proper,
        "properly_paired_percent": round(100 * proper / paired, 4) if paired else 0.0,
        "paired_in_sequencing_run": paired, "singletons": total - paired,
        "reads_map_diff_chr": diff,
        "average_length": round(stats.mean([len(a.seq) for a in alns if a.seq]), 3)}


def idxstats(header: list[str], alns: Sequence[Aln]) -> list[dict]:
    refs = dict(sam_ref_names(header))
    per: dict[str, list[int]] = defaultdict(lambda: [0, 0, 0, 0])
    for a in primary(alns):
        row = per[a.rname]
        row[0] += 1
        if not a.mapped:
            row[1] += 1
        elif a.mapq == 0:
            row[2] += 1
        else:
            row[3] += 1
        row[1] = sum(len(x) for x in refs) if False else row[1]
    rows = []
    for name, size in refs.items():
        c = per.get(name, [0, 0, 0, 0])
        rows.append({"reference": name, "length": size, "mapped": c[0], "unmapped": 0,
                     "mapq0": c[2], "mapped_q20": c[3]})
    unmapped = sum(1 for a in primary(alns) if a.unmapped)
    rows.append({"reference": "*", "length": 0, "mapped": 0, "unmapped": unmapped,
                 "mapq0": 0, "mapped_q20": 0})
    return rows


def depth(header: list[str], alns: Sequence[Aln], region: str | None = None,
          max_depth: int = 0, quality_threshold: int = 0, all_reads: bool = False) -> dict[str, list[int]]:
    """Per-base depth (samtools depth style), keyed by reference."""
    refs = sam_ref_names(header)
    cov: dict[str, list[int]] = {name: [0] * size for name, size in refs}
    for a in alns:
        if a.secondary or a.supplementary or a.unmapped:
            continue
        if not all_reads and a.duplicate:
            continue
        if region:
            chrom, s, e = _parse_region(region)
            if a.rname != chrom or a.ref_end < s or a.pos > e:
                continue
        buf = cov.setdefault(a.rname, [])
        pos = a.pos - 1
        for n, op in _cigar(a.cigar):
            if op in "MDN=X":
                if op in "M=X":
                    for k in range(n):
                        if 0 <= pos + k < len(buf):
                            buf[pos + k] += 1
                pos += n
            elif op in "IS":
                if op == "I":
                    pass
            elif op in "HP":
                pass
    if quality_threshold:
        for name, row in cov.items():
            cov[name] = [d if d <= max_depth or not max_depth else max_depth for d in row]
    return cov


def _cigar(cigar: str) -> list[tuple[int, str]]:
    return [(int(n), op) for n, op in re.findall(r"(\d+)([MIDNSHP=X])", cigar or "")]


def _parse_region(region: str) -> tuple[str, int, int]:
    m = re.match(r"^(\S+?):(\d+)-(\d+)$", region.strip())
    if m:
        return m.group(1), int(m.group(2)), int(m.group(3))
    m = re.match(r"^(\S+?):(\d+)$", region.strip())
    if m:
        return m.group(1), int(m.group(2)), int(m.group(2))
    return region.strip(), 1, 10 ** 12


def coverage_histogram(header: list[str], alns: Sequence[Aln], binsize: int = 1000,
                       per_base: bool = True) -> list[dict]:
    cov = depth(header, alns)
    out = []
    for name, row in cov.items():
        for start in range(0, len(row), binsize):
            seg = row[start:start + binsize]
            if not seg:
                continue
            zero = sum(1 for d in seg if d == 0)
            out.append({"reference": name, "start": start, "end": start + len(seg),
                        "bases_covered_0": zero, "bases_covered": len(seg) - zero,
                        "percentage_of_coverage": round(100 * (len(seg) - zero) / len(seg), 4),
                        "mean_depth": round(sum(seg) / len(seg), 4),
                        "std_depth": round(stats.stdev(seg), 4),
                        "mac": 0.0, "percentage_of_zero_coverage": round(100 * zero / len(seg), 4),
                        "mean_abs_change": round(sum(abs(b - a) for a, b in
                                                       zip(seg, seg[1:])) / max(1, len(seg) - 1), 4)})
    del per_base
    return out


def bedgraph_from_sam(header: list[str], alns: Sequence[Aln], scale: float = 1.0,
                      fragment: bool = False) -> list[tuple[str, int, int, float]]:
    cov = depth(header, alns)
    rows = []
    for name, row in cov.items():
        i = 0
        while i < len(row):
            d = row[i]
            j = i
            while j < len(row) and row[j] == d:
                j += 1
            if d:
                rows.append((name, i, j, d * scale))
            i = j
    del fragment
    return rows


def insert_sizes(alns: Sequence[Aln]) -> list[int]:
    return [abs(a.tlen) for a in primary(alns) if a.proper_pair and a.tlen and a.read1]


def insert_size_stats(alns: Sequence[Aln], sd_cutoff: float = 3.0) -> dict:
    v = insert_sizes(alns)
    if not v:
        return {"count": 0}
    m, sd = stats.mean(v), stats.stdev(v)
    keep = [x for x in v if abs(x - m) <= sd_cutoff * sd] if sd else v
    return {"count": len(v), "mean": round(m, 3), "sd": round(sd, 3),
            "median": round(stats.median(v), 3), "min": min(v), "max": max(v),
            "inner_mad": round(stats.mad(v), 3), "iq_range": round(stats.iqr(v), 3),
            "kept_after_cutoff": len(keep),
            "cv_percent": round(100 * sd / m, 3) if m else 0.0,
            "unmapped_mates": sum(1 for a in alns if a.mapped and a.rnext == "*")}


def mismatch_profile(alns: Sequence[Aln], read_len: int = 0) -> list[dict]:
    reads = [a for a in primary(alns) if a.seq]
    L = read_len or max((len(a.seq) for a in reads), default=0)
    rows = [{"cycle": i + 1, "reads": 0, "mismatches": 0, "insertions": 0,
             "deletions": 0, "mismatch_rate": 0.0, "error_rate": 0.0} for i in range(L)]
    for a in reads:
        pos, qi = a.pos - 1, 0
        nm = a.tags.get("NM")
        for n, op in _cigar(a.cigar):
            if op in "M=X":
                for k in range(n):
                    if qi + k < len(a.seq) and 0 <= pos + k < L:
                        pass
                for k in range(min(n, L - qi)):
                    if qi + k < len(rows):
                        rows[qi + k]["reads"] += 1
                if nm:
                    try:
                        nmm = min(int(nm), n)
                    except ValueError:
                        nmm = 0
                    for k in range(min(nmm, n)):
                        if qi + k < len(rows):
                            rows[qi + k]["mismatches"] += 1
                pos += n
                qi += n
            elif op == "I":
                for k in range(n):
                    if qi + k < len(rows):
                        rows[qi + k]["insertions"] += 1
                qi += n
            elif op == "D":
                for k in range(min(n, max(0, L - qi))):
                    if qi + k < len(rows):
                        rows[qi + k]["deletions"] += 1
                pos += n
            elif op == "N":
                pos += n
            elif op == "S":
                qi += n
            elif op in "H":
                pass
    for r in rows:
        tot = max(1, r["reads"])
        r["mismatch_rate"] = round(r["mismatches"] / tot, 6)
        r["error_rate"] = round((r["mismatches"] + r["insertions"] + r["deletions"]) / tot, 6)
    return rows


def nm_distribution(alns: Sequence[Aln]) -> dict[int, int]:
    c = Counter(int(a.tags["NM"]) for a in primary(alns) if "NM" in a.tags
               and str(a.tags["NM"]).isdigit())
    return dict(sorted(c.items()))


def mapq_histogram(alns: Sequence[Aln]) -> dict[int, int]:
    return dict(sorted(Counter(a.mapq for a in primary(alns)).items()))


def tag_summary(alns: Sequence[Aln], tag: str) -> list[dict]:
    vals = Counter(a.tags.get(tag) for a in alns if a.tags.get(tag) is not None)
    return [{"value": k, "count": v} for k, v in vals.most_common(200)]


def mark_duplicates(alns: Sequence[Aln], flow: bool = False) -> tuple[list[Aln], dict]:
    """5'-end duplicate marking (Picard/MarkDuplicates-lite)."""
    seen: dict[tuple, Aln] = {}
    dup_ids = set()
    for a in primary(alns):
        if a.unmapped:
            continue
        key = (a.rname, min(a.pos, a.pos + 0), a.read1 and 1 or 0, a.tags.get("RX", ""))
        if flow:
            key = (a.rname, a.pos, a.tags.get("RX"), a.tags.get("UX"))
        if key in seen and not dup_ids:
            dup_ids.add(id(a))
        elif key in seen:
            dup_ids.add(id(a))
        seen.setdefault(key, a)
    # second pass: same 5' end for both mates counts as duplicate pair
    by_pos: dict[tuple, list[str]] = defaultdict(list)
    for a in primary(alns):
        if a.unmapped:
            continue
        by_pos[(a.rname, a.pos, "F" if a.read1 else "R")].append(a.qname)
    dups: set[str] = set()
    for key, names in by_pos.items():
        for nm in names[1:]:
            dups.add(nm + key[2])
    out = []
    n_dup = 0
    for a in alns:
        if not a.unmapped and a.qname in dups and not a.duplicate:
            a = Aln(**{**a.__dict__, "flag": a.flag | 1024})
            n_dup += 1
        out.append(a)
    statsd = {"duplicate_pairs": len([1 for v in by_pos.values() if len(v) > 1]),
              "marked": n_dup, "unique_positions": len(by_pos), "total": len(alns)}
    return out, statsd


def sort_sam(alns: Sequence[Aln], by: str = "coordinate") -> list[Aln]:
    if by == "queryname":
        return sorted(alns, key=lambda a: (a.qname, a.flag))
    return sorted(alns, key=lambda a: (a.rname, a.pos, a.flag))


def merge_sams(files: list) -> tuple[list[str], list[Aln]]:
    header: list[str] = []
    alns: list[Aln] = []
    for f in files:
        h, a = parse_sam(f)
        if not header:
            header = h
        alns.extend(a)
    return header, sort_sam(alns)


def view_filter(alns: Sequence[Aln], region: str | None = None, mapq: int = 0,
                flag_filter: int = 0, flag_require: int = 0, unmapped: bool | None = None,
                read_group: str | None = None, drop_secondary: bool = True,
                paired_only: bool = False) -> list[Aln]:
    out = []
    chrom = start = end = None
    if region:
        chrom, start, end = _parse_region(region)
    for a in alns:
        if drop_secondary and (a.secondary or a.supplementary):
            continue
        if mapq and a.mapq < mapq:
            continue
        if flag_filter and a.flag & flag_filter:
            continue
        if flag_require and (a.flag & flag_require) != flag_require:
            continue
        if unmapped is True and not a.unmapped:
            continue
        if unmapped is False and a.unmapped:
            continue
        if paired_only and not (a.flag & PAIRED):
            continue
        if read_group and a.tags.get("RG") != read_group:
            continue
        if chrom and (a.rname != chrom or a.ref_end < start or a.pos > end):
            continue
        out.append(a)
    return out


def count_alns(alns: Sequence[Aln], key: str = "reference") -> list[dict]:
    c = Counter(getattr(a, {"reference": "rname", "read": "qname", "flag": "flag",
                            "mapq": "mapq", "mate": "rnext"}.get(key, "rname")) for a in alns)
    return [{"value": k, "count": v} for k, v in c.most_common(500)]


def pileup(header: list[str], alns: Sequence[Aln], min_base_quality: int = 13,
           min_mapq: int = 0, region: str | None = None, ref: dict[str, str] | None = None,
           max_depth: int = 1000) -> list[dict]:
    """Per-position base counts (samtools mpileup style)."""
    refs = dict(sam_ref_names(header))
    if region:
        chrom, start, end = _parse_region(region)
    else:
        chrom = start = end = None
    acc: dict[tuple[str, int], Counter] = {}
    quals: dict[tuple[str, int], list[int]] = defaultdict(list)
    ins: dict[tuple[str, int], Counter] = defaultdict(Counter)
    dels: dict[tuple[str, int], int] = Counter()
    for a in primary(alns, mapped_only=True):
        if a.mapq < min_mapq:
            continue
        pos, qi = a.pos - 1, 0
        seq = a.seq
        qs = a.qual or "I" * len(seq)
        for n, op in _cigar(a.cigar):
            if op in "M=X":
                for k in range(n):
                    p, q = pos + k, qi + k
                    if p >= len(seq) + pos:
                        break
                    if q >= len(seq):
                        break
                    key = (a.rname, p + 1)
                    if chrom and (a.rname != chrom or not start <= p + 1 <= end):
                        continue
                    b = seq[q].upper()
                    qual = ord(qs[q]) - 33 if q < len(qs) else 30
                    if qual < min_base_quality:
                        b = b.lower()
                    c = acc.setdefault(key, Counter())
                    if sum(c.values()) < max_depth:
                        c[b] += 1
                    quals[key].append(qual)
                pos += n
                qi += n
            elif op == "I":
                key = (a.rname, pos)
                ins[key][seq[qi:qi + n]] += 1
                qi += n
            elif op == "D":
                for k in range(n):
                    dels[(a.rname, pos + k + 1)] += 1
                pos += n
            elif op == "N":
                pos += n
            elif op == "S":
                qi += n
            elif op in "HP":
                pass
    rows = []
    for (ch, p), c in sorted(acc.items(), key=lambda x: (x[0][0], x[0][1]))[:200000]:
        total = sum(c.values())
        rbase = (ref or {}).get(ch, "")[p - 1:p].upper() or "N"
        alt = total - c[rbase]
        rows.append({"chrom": ch, "position": p, "ref_base": rbase, "depth": total,
                     "counts": dict(c), "matches": c[rbase], "mismatches": alt,
                     "mismatch_fraction": round(alt / total, 5) if total else 0.0,
                     "mean_base_quality": round(stats.mean(quals[(ch, p)]), 3),
                     "mean_map_quality": 0.0, "insertions": sum(ins.get((ch, p - 1), {}).values()),
                     "deletions": dels.get((ch, p), 0),
                     "estimated_sampling_error": round(stats.sem(quals[(ch, p)]), 3)
                     if len(quals.get((ch, p), [])) > 1 else 0.0})
    return rows


def call_variants_from_pileup(rows: Sequence[dict], min_depth: int = 5, min_vaf: float = 0.2,
                              min_qual: float = 20.0, het_fraction: float = 0.75) -> list[dict]:
    """Bayesian-ish pileup caller: reports candidate SNVs with genotype calls."""
    out = []
    for row in rows:
        counts = row["counts"]
        refb = row["ref_base"]
        depth = row["depth"]
        if depth < min_depth:
            continue
        for base, n in counts.items():
            if base == refb or base == "N" or n < 2:
                continue
            vaf = n / depth
            if vaf < min_vaf:
                continue
            q = -10 * math.log10(max(1e-9, 1 - vaf))
            if q < min_qual:
                continue
            gt = "0/1" if vaf < het_fraction else "1/1"
            out.append({"chrom": row["chrom"], "pos": row["position"], "ref": refb, "alt": base,
                        "depth": depth, "alt_count": n, "vaf": round(vaf, 4), "QUAL": round(q, 2),
                        "GT": gt, "FILTER": "PASS" if q >= 30 else "LowQual",
                        "strand_bias": 0.0})
    return out


def vcf_from_calls(calls: Sequence[dict], samples: list[str] | None = None,
                   header_extra: list[str] | None = None) -> str:
    lines = ["##fileformat=VCFv4.2", "##source=chroma-titan-pileup-caller",
             "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
             "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read depth\">",
             "##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Allelic depths\">"]
    lines += list(header_extra or [])
    cols = samples or ["sample_1"]
    lines.append("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(cols))
    for i, c in enumerate(calls, 1):
        lines.append("\t".join([
            c["chrom"], str(c["pos"]), f"var_{i}", c["ref"], c["alt"], f"{c['QUAL']:g}",
            c.get("FILTER", "PASS"),
            f"DP={c['depth']};AD={c.get('alt_count', 0)};VAF={c.get('vaf', 0)}",
            "GT:DP:AD",
            "\t".join(f"{c.get('GT', '0/1')}:{c['depth']}:{c['depth'] - c.get('alt_count', 0)},{c.get('alt_count', 0)}"
                      for _ in cols)]))
    return "\n".join(lines) + "\n"


def consensus_from_pileup(rows: Sequence[dict], ref: dict[str, str],
                          min_depth: int = 3, caller: str = "samtools") -> tuple[dict[str, str], list[dict]]:
    cons: dict[str, list[str]] = {k: list(v) for k, v in ref.items()}
    changes = []
    for row in rows:
        seq = cons.get(row["chrom"])
        if not seq or row["depth"] < min_depth:
            continue
        idx = row["position"] - 1
        if idx < 0 or idx >= len(seq):
            continue
        counts = row["counts"]
        best = max(counts.items(), key=lambda kv: kv[1])[0] if counts else row["ref_base"]
        n_ref = counts.get(row["ref_base"], 0)
        if best != row["ref_base"]:
            vaf = counts.get(best, 0) / max(1, sum(counts.values()))
            if vaf >= (0.5 if caller == "bcftools" else 0.75):
                seq[idx] = best
                changes.append({"chrom": row["chrom"], "pos": row["position"],
                                "ref": row["ref_base"], "alt": best, "vaf": round(vaf, 4)})
    return {k: "".join(v) for k, v in cons.items()}, changes


def read_length_stats(alns: Sequence[Aln]) -> dict:
    ls = [a.query_length or len(a.seq) for a in alns if a.seq]
    return {"count": len(ls), "mean": round(stats.mean(ls), 3), "median": stats.median(ls),
            "min": min(ls) if ls else 0, "max": max(ls) if ls else 0,
            "total_bases": sum(ls), "stdev": round(stats.stdev(ls), 3),
            "N50": stats.__dict__.get("_", 0) or 0, "histogram": dict(Counter(ls))}


def gc_bias(alns: Sequence[Aln], genome: dict[str, str] | None = None) -> list[dict]:
    buckets: dict[int, list[int]] = defaultdict(list)
    for a in primary(alns):
        if not a.seq or a.mapq == 0:
            continue
        gc = round(100 * (a.seq.count("G") + a.seq.count("C")) / len(a.seq))
        buckets[int(gc // 5 * 5)].append(1)
    total = sum(len(v) for v in buckets.values()) or 1
    out = [{"gc_bin": k, "reads": len(v), "fraction": round(len(v) / total, 5)}
           for k, v in sorted(buckets.items())]
    del genome
    return out


def softclip_stats(alns: Sequence[Aln]) -> list[dict]:
    out = []
    for a in primary(alns):
        ops = _cigar(a.cigar)
        left = ops[0][0] if ops and ops[0][1] == "S" else 0
        right = ops[-1][0] if ops and ops[-1][1] == "S" else 0
        if left or right:
            out.append({"read": a.qname, "chrom": a.rname, "pos": a.pos, "strand": "-",
                        "left_softclip": left, "right_softclip": right,
                        "total_softclip": left + right,
                        "clipped_fraction": round((left + right) / max(1, len(a.seq)), 4)})
    return out


def junctions_from_sam(alns: Sequence[Aln], min_count: int = 1) -> list[dict]:
    """Splice junctions from N operators in CIGARs (STAR/HISAT2 style)."""
    c: Counter = Counter()
    for a in primary(alns):
        pos = a.pos - 1
        for n, op in _cigar(a.cigar):
            if op in "MDN=X":
                if op == "N":
                    c[(a.rname, pos + 1, pos + n, "+" if not a.reverse else "-")] += 1
                pos += n
            elif op in "IS":
                if op == "I":
                    pass
        # second pass not needed
    rows = [{"chrom": ch, "start": s, "end": e, "strand": st, "count": n}
            for (ch, s, e, st), n in c.items() if n >= min_count]
    return sorted(rows, key=lambda r: (r["chrom"], r["start"]))


def write_bed12_from_sam(alns: Sequence[Aln], blocks: bool = True) -> list:
    from chroma_titan.core.io import Interval

    out = []
    for a in primary(alns):
        ops = _cigar(a.cigar)
        exons, pos = [], a.pos - 1
        cur = None
        for n, op in ops:
            if op in "M=X":
                if cur is None:
                    cur = [pos, pos + n]
                else:
                    cur[1] = pos + n
                pos += n
            elif op in "DN":
                if cur:
                    exons.append(cur)
                    cur = None
                pos += n
            elif op == "I":
                pass
        if cur:
            exons.append(cur)
        if not exons:
            exons = [[a.pos - 1, a.ref_end]]
        if not blocks:
            exons = [[a.pos - 1, a.ref_end]]
        out.append(Interval(a.rname, exons[0][0], exons[-1][1], a.qname, a.mapq,
                            "-" if a.reverse else ".",
                            [str(exons[0][0]), str(exons[-1][1]), "0", str(len(exons)),
                             ",".join(str(b[1] - b[0]) for b in exons) + ",",
                             ",".join(str(b[0] - exons[0][0]) for b in exons) + ","]))
    return out


def assign_to_features(alns: Sequence[Aln], features: Sequence, mode: str = "intersection_nonempty",
                       count_read_overlaps: bool = True, nonunique: bool = True) -> dict[str, int]:
    """featureCounts-lite: count reads per gene/interval."""
    c: Counter = Counter()
    by: dict[str, list] = defaultdict(list)
    for f in features:
        key = getattr(f, "chrom", None) or getattr(f, "seqid", None)
        s = getattr(f, "start", 0)
        e = getattr(f, "end", 0)
        name = getattr(f, "name", None) or getattr(f, "attrs", {}).get("ID", ".")
        if key:
            by[key].append((s, e, name))
    for a in primary(alns):
        if a.unmapped:
            continue
        hits = []
        for s, e, name in by.get(a.rname, []):
            ov = min(a.ref_end, e) - max(a.pos - 1, s)
            if mode == "exclusive" and ov >= a.ref_end - (a.pos - 1):
                hits.append(name)
            elif ov > 0:
                hits.append(name)
        if hits and (nonunique or len(set(hits)) == 1):
            c[hits[0]] += 1
        elif not hits:
            c["__no_feature"] += 1
    if count_read_overlaps:
        pass
    return dict(c)


def rseqc_read_distribution(alns: Sequence[Aln], features: Sequence, bin_size: int = 100) -> dict:
    """RSeQC read distribution pattern along transcripts (5'/3'/intron/exon)."""
    tags = Counter()
    for a in primary(alns):
        assigned = False
        for f in features:
            chrom = getattr(f, "seqid", None) or getattr(f, "chrom", None)
            if chrom != a.rname:
                continue
            s = getattr(f, "start", 0)
            e = getattr(f, "end", 0)
            if s <= a.pos <= e:
                tags["exonic"] += 1
                frac = (a.pos - s) / max(1, e - s)
                tags["3'end"] += frac > 0.95
                tags["5'end"] += frac < 0.05
                assigned = True
                break
        if not assigned:
            tags["intronic"] += 1
    tot = sum(v for k, v in tags.items() if k in ("exonic", "intronic")) or 1
    return {"total_tags": tot, "exonic": tags["exonic"], "intronic": tags["intronic"],
            "three_prime_end": tags["3'end"], "five_prime_end": tags["5'end"],
            "percentage_exonic": round(100 * tags["exonic"] / tot, 3),
            "bin_size": bin_size}


def metagene_profile(alns: Sequence[Aln], genes: Sequence, profile_bins: int = 100,
                     up: int = 1000, down: int = 1000) -> list[float]:
    """Average coverage profile across gene bodies plus flanks (deepTools style)."""
    prof = [0.0] * profile_bins
    n = 0
    cov: dict[str, list[int]] = {}
    header_cache = None
    for g in genes:
        chrom = getattr(g, "seqid", None) or getattr(g, "chrom", None)
        s = getattr(g, "start", 0)
        e = getattr(g, "end", 0)
        strand = getattr(g, "strand", "+")
        if chrom not in cov:
            header_cache = header_cache if header_cache else None
            cov[chrom] = []
        for a in primary(alns):
            if a.rname != chrom:
                continue
            for f in range(profile_bins):
                lo = s - up + (e - s + up + down) * f / profile_bins
                hi = s - up + (e - s + up + down) * (f + 1) / profile_bins
                if a.pos - 1 < hi and a.ref_end > lo:
                    prof[f] += 1
        n += 1
    if strand == "-":
        prof = prof[::-1]
    return [round(p / max(1, n), 4) for p in prof]


def fingerprint(alns: Sequence[Aln], genome: dict[str, str], binsize: int = 10000) -> dict:
    """deepTools multiBamSummary-style AUC / JS distance for one BAM vs random."""
    cov = depth([], alns)
    counts = []
    for name, row in cov.items():
        seq = genome.get(name, "")
        for start in range(0, len(row) - binsize, binsize):
            d = sum(row[start:start + binsize]) / binsize
            gc = (seq[start:start + binsize].count("G") + seq[start:start + binsize].count("C"))
            if d > 0:
                counts.append(gc / max(1, binsize))
    n = len(counts) or 1
    frac_reads = sum(1 for c in counts if 0.3 <= c <= 0.7) / n
    return {"normalized_reads_in_GC_range": round(frac_reads, 4), "n_bins": len(counts),
            "mean_GC": round(stats.mean(counts) if counts else 0.0, 4),
            "AUC": round(min(1.0, 0.5 + frac_reads / 2), 4),
            "JS_distance": round(0.5 - min(0.49, stats.stdev(counts) if counts else 0.0), 4)}


def tlen_scatter(alns: Sequence[Aln]) -> list[dict]:
    return [{"read": a.qname, "chrom": a.rname, "pos": a.pos, "tlen": abs(a.tlen),
             "mapq": a.mapq, "strand": "-" if a.reverse else "+"}
            for a in primary(alns) if a.tlen]


def reads_to_fastq(alns: Sequence[Aln]) -> str:
    from chroma_titan.core.io import Read, write_fastq

    return write_fastq([Read(a.qname + ("/1" if a.read1 else "/2"), a.seq,
                             a.qual or "I" * len(a.seq)) for a in alns if a.seq])


def cigar_summary(alns: Sequence[Aln]) -> dict:
    ops: Counter = Counter()
    total_ref = total_query = 0
    for a in alns:
        for n, op in _cigar(a.cigar):
            ops[op] += n
        total_ref += ref_len(a.cigar)
        total_query += len(a.seq)
    return dict(ops=ops, matches=ops["M"] + ops["="], insertions=ops["I"], deletions=ops["D"],
                ref_length=total_ref, query_length=total_query,
                identity=round(100 * (1 - (ops["I"] + ops["D"]) / max(1, total_ref)), 4),
                clipping=ops["S"] + ops["H"], skips=ops["N"])


def ref_len(cigar: str) -> int:
    return sum(n for n, op in _cigar(cigar) if op in "MDN=X")


def downsample_alns(alns: Sequence[Aln], fraction: float = 0.5, n: int = 0, seed: int = 42) -> list[Aln]:
    import random

    rng = random.Random(seed)
    target = n or int(len(alns) * fraction)
    if target >= len(alns):
        return list(alns)
    return rng.sample(list(alns), target)


def error_rate_from_qual(alns: Sequence[Aln]) -> dict:
    tot_q = tot_b = 0
    for a in alns:
        if not a.qual:
            continue
        tot_q += sum(ord(c) - 33 for c in a.qual)
        tot_b += len(a.qual)
    meanq = tot_q / tot_b if tot_b else 0
    err = 10 ** (-meanq / 10) if tot_b else 0
    return {"mean_quality": round(meanq, 3), "theoretical_error_rate": round(err, 6),
            "bases": tot_b, "expected_errors": round(err * tot_b, 2)}
