"""Genomic interval algebra (BEDtools-style) built on :class:`Interval`.

Half-open coordinates (BED style) throughout: ``[start, end)``.
"""

from __future__ import annotations

import math
import random
from collections import Counter, defaultdict
from typing import Iterable, Sequence

from chroma_titan.core.io import Interval, parse_bed, write_bed


def sort_ivs(ivs: Sequence[Interval], by: str = "coordinate") -> list[Interval]:
    out = list(ivs)
    if by == "coordinate":
        out.sort(key=lambda x: (x.chrom, x.start, x.end))
    elif by == "size":
        out.sort(key=lambda x: -x.length)
    elif by == "name":
        out.sort(key=lambda x: (x.chrom, x.name))
    elif by == "score":
        out.sort(key=lambda x: -x.score)
    elif by == "chrom":
        out.sort(key=lambda x: x.chrom)
    return out


def merge(ivs: Sequence[Interval], distance: int = 0, collapse: str = "concat",
          strand: bool = False) -> list[Interval]:
    key = (lambda x: (x.chrom, x.strand)) if strand else (lambda x: x.chrom)
    by: dict = defaultdict(list)
    for iv in ivs:
        by[key(iv)].append(iv)
    out: list[Interval] = []
    for k, rows in by.items():
        rows.sort(key=lambda x: (x.start, x.end))
        cur = None
        for iv in rows:
            if cur and iv.start <= cur.end + distance:
                cur.end = max(cur.end, iv.end)
                cur.score = max(cur.score, iv.score)
                names = cur.extra[0].split("|") if cur.extra and cur.extra[0] else [cur.name]
                if iv.name not in names:
                    cur.extra = ["|".join(names + [iv.name])]
                    cur.name = _collapse_names(names + [iv.name], collapse)
            else:
                if cur:
                    out.append(cur)
                cur = Interval(iv.chrom, iv.start, iv.end, iv.name, iv.score,
                               iv.strand if strand else ".", list(iv.extra))
                if cur.extra and collapse != "keep":
                    cur.extra = []
                cur.name = iv.name
        if cur:
            out.append(out and out[-1] is cur and cur or cur)
            if cur.extra and cur.extra[0]:
                names = cur.extra[0].split("|")
                if len(names) > 1:
                    cur.name = _collapse_names(names, collapse)
    return sort_ivs(out)


def _collapse_names(names: list[str], mode: str) -> str:
    names = [n for n in names if n not in (".", "")]
    if mode == "count":
        return str(len(names))
    if mode == "sum" or mode == "mean":
        vals = []
        for n in names:
            try:
                vals.append(float(n))
            except ValueError:
                pass
        if vals:
            return str(sum(vals) if mode == "sum" else round(sum(vals) / len(vals), 4))
    if mode == "distinct":
        return "|".join(sorted(set(names)))
    return "|".join(names) if names else "."


def intersect(a: Sequence[Interval], b: Sequence[Interval], wa: bool = False,
              wb: bool = False, strand: bool = False) -> list[Interval]:
    idx = _index(b)
    out: list[Interval] = []
    for iv in a:
        for other in idx.get(iv.chrom, []):
            if strand and iv.strand != other.strand != ".":
                if iv.strand != other.strand:
                    continue
            s, e = max(iv.start, other.start), min(iv.end, other.end)
            if e > s:
                extra = ([f"{e - s}"] if wa else []) + ([other.name] if wb else [])
                out.append(Interval(iv.chrom, s, e, iv.name, iv.score,
                                    iv.strand, extra))
    return sort_ivs(out)


def subtract(a: Sequence[Interval], b: Sequence[Interval], strand: bool = False) -> list[Interval]:
    idx = _index(b)
    out: list[Interval] = []
    for iv in a:
        blocks = [iv]
        for other in sorted(idx.get(iv.chrom, []), key=lambda x: x.start):
            nxt = []
            for blk in blocks:
                if other.end <= blk.start or other.start >= blk.end:
                    nxt.append(blk)
                    continue
                if blk.start < other.start:
                    nxt.append(Interval(blk.chrom, blk.start, other.start, blk.name,
                                        blk.score, blk.strand, list(blk.extra)))
                if other.end < blk.end:
                    nxt.append(Interval(blk.chrom, other.end, blk.end, blk.name, blk.score,
                                       blk.strand, list(blk.extra)))
            blocks = nxt
        out.extend(blocks)
    return sort_ivs(out)


def complement(genome: dict[str, int], ivs: Sequence[Interval]) -> list[Interval]:
    by: dict[str, list[Interval]] = defaultdict(list)
    for iv in ivs:
        by[iv.chrom].append(iv)
    out: list[Interval] = []
    for chrom, size in genome.items():
        rows = sorted(by.get(chrom, []), key=lambda x: x.start)
        pos = 0
        for iv in rows:
            if iv.start > pos:
                out.append(Interval(chrom, pos, iv.start, "."))
            pos = max(pos, iv.end)
        if pos < size:
            out.append(Interval(chrom, pos, size, "."))
    return sort_ivs(out)


def closest(a: Sequence[Interval], b: Sequence[Interval]) -> list[dict]:
    idx = _index(b)
    out = []
    for iv in a:
        best = None
        for other in idx.get(iv.chrom, []):
            d = other.start - iv.end if other.start >= iv.end else (
                iv.start - other.end if other.end <= iv.end else 0)
            cand = {"chrom": iv.chrom, "start": iv.start, "end": iv.end, "name": iv.name,
                    "score": iv.score, "strand": iv.strand, "other_chrom": other.chrom,
                    "other_start": other.start, "other_end": other.end,
                    "other_name": other.name, "other_score": other.score,
                    "other_strand": other.strand, "distance": d}
            if best is None or abs(d) < abs(best["distance"]):
                best = cand
        if best is None:
            best = {"chrom": iv.chrom, "start": iv.start, "end": iv.end, "name": iv.name,
                    "score": iv.score, "strand": iv.strand, "other_chrom": "nan",
                    "other_start": -1, "other_end": -1, "other_name": "nan",
                    "other_score": 0, "other_strand": ".", "distance": -1}
        out.append(best)
    return out


def window(a: Sequence[Interval], b: Sequence[Interval], w: int = 100) -> list[dict]:
    idx = _index(b)
    out = []
    for iv in a:
        for other in idx.get(iv.chrom, []):
            if other.end > iv.start - w and other.start < iv.end + w:
                d = other.start - iv.end if other.start >= iv.end else (
                    iv.start - other.end if other.end <= iv.end else 0)
                out.append({"chrom": iv.chrom, "start": iv.start, "end": iv.end,
                            "name": iv.name, "other_name": other.name,
                            "other_start": other.start, "other_end": other.end,
                            "distance": d})
    return out


def slop(ivs: Sequence[Interval], genome: dict[str, int] | None = None,
         bp: int = 100, direction: str = "both", fraction: float = 0.0) -> list[Interval]:
    out = []
    for iv in ivs:
        add = int(iv.length * fraction) if fraction else bp
        a = iv.start - (add if direction in ("both", "left") else 0)
        b = iv.end + (add if direction in ("both", "right") else 0)
        if genome:
            a = max(0, min(a, genome.get(iv.chrom, 10 ** 12)))
            b = max(0, min(b, genome.get(iv.chrom, 10 ** 12)))
        else:
            a = max(0, a)
        out.append(Interval(iv.chrom, a, b, iv.name, iv.score, iv.strand, list(iv.extra)))
    return out


def flank(ivs: Sequence[Interval], left: int = 100, right: int = 100,
          both: bool = False) -> list[Interval]:
    out = []
    for iv in ivs:
        if both:
            out.append(Interval(iv.chrom, max(0, iv.start - left), iv.start,
                               f"{iv.name}_up", iv.score, iv.strand))
            out.append(Interval(iv.chrom, iv.end, iv.end + right, f"{iv.name}_down",
                                iv.score, iv.strand))
        else:
            out.append(Interval(iv.chrom, max(0, iv.start - left), iv.end + right,
                                iv.name, iv.score, iv.strand))
    return out


def shift(ivs: Sequence[Interval], by: int = 0) -> list[Interval]:
    return [Interval(iv.chrom, max(0, iv.start + by), max(1, iv.end + by), iv.name,
                     iv.score, iv.strand, list(iv.extra)) for iv in ivs]


def shuffle_in(genome: dict[str, int], ivs: Sequence[Interval], seed: int = 42,
               excl: Sequence[Interval] | None = None) -> list[Interval]:
    rng = random.Random(seed)
    bad = merge(excl or [], distance=1)
    out = []
    for iv in ivs:
        size = genome.get(iv.chrom, 100000)
        for _ in range(200):
            s = rng.randint(0, max(0, size - iv.length))
            cand = Interval(iv.chrom, s, s + iv.length, iv.name, iv.score, iv.strand)
            if not any(iv2.chrom == cand.chrom and cand.start < iv2.end and cand.end > iv2.start
                       for iv2 in bad):
                out.append(cand)
                break
        else:
            out.append(cand)
    return sort_ivs(out)


def bamsiflike_coverage(b: Sequence[Interval], genome: dict[str, int]) -> list[tuple[str, int, int, float]]:
    rows: list[tuple[str, int, int, float]] = []
    by = defaultdict(list)
    for iv in b:
        by[iv.chrom].append(iv)
    for chrom, size in genome.items():
        depth = [0] * size
        for iv in by.get(chrom, []):
            for p in range(max(0, iv.start), min(size, iv.end)):
                depth[p] += 1
        run, start = None, 0
        for p, d in enumerate(depth + [None]):  # type: ignore[list-item]
            if d != run:
                if run is not None and p > start:
                    rows.append((chrom, start, p, float(run)))
                run, start = d, p
    return rows


def coverage(a: Sequence[Interval], b: Sequence[Interval] | None = None,
             counts: dict[str, float] | None = None, mean: bool = False) -> list[dict]:
    """per-base coverage of A by B (or by given numeric scores)."""
    out = []
    by = defaultdict(list)
    for iv in (b or []):
        by[iv.chrom].append(iv)
    for iv in sort_ivs(a):
        hits = [o for o in by.get(iv.chrom, []) if o.end > iv.start and o.start < iv.end]
        bases = 0
        for o in hits:
            bases += min(iv.end, o.end) - max(iv.start, o.start)
        val = 0.0
        if counts:
            val = sum(float(counts.get(o.name, 0) or 0) for o in hits)
        else:
            val = float(len(hits))
        rec = {"chrom": iv.chrom, "start": iv.start, "end": iv.end, "name": iv.name,
               "hits": len(hits), "covered_bases": bases, "coverage": val}
        if mean:
            rec["bases"] = iv.length
            rec["mean_coverage"] = round(val / iv.length, 4) if iv.length else 0.0
        out.append(rec)
    return out


def jaccard(a: Sequence[Interval], b: Sequence[Interval]) -> dict:
    ma, mb = merge(a), merge(b)
    inter = sum(iv.length for iv in intersect(ma, mb))
    union = sum(iv.length for iv in merge(list(ma) + list(mb)))
    return {"intersection": inter, "union": union,
            "jaccard": round(inter / union, 6) if union else 0.0}


def merge_stat(a: Sequence[Interval], b: Sequence[Interval]) -> dict:
    inter = sum(x.length for x in intersect(merge(a), merge(b)))
    return {"total_a": sum(x.length for x in a), "total_b": sum(x.length for x in b),
            "overlap": inter,
            "percent_a_covered": round(100 * inter / sum(x.length for x in a), 4)
            if a else 0.0}


def multiinter(sets: list[list[Interval]]) -> list[dict]:
    """bedtools multiinter-like: interval, count, membership."""
    pts = set()
    for rows in sets:
        for iv in rows:
            pts.add((iv.chrom, iv.start))
            pts.add((iv.chrom, iv.end))
    by_chrom: dict[str, list] = defaultdict(list)
    for rows in sets:
        for i, iv in enumerate(rows):
            by_chrom[iv.chrom].append((iv, i))
    out = []
    for chrom in sorted({c for c, _ in pts}):
        bounds = sorted({p for c, p in pts if c == chrom})
        for a, b in zip(bounds, bounds[1:]):
            if b <= a:
                continue
            members = sorted({i for iv, i in by_chrom[chrom]
                              if iv.start < b and iv.end > a})
            if members:
                out.append({"chrom": chrom, "start": a, "end": b, "count": len(members),
                            "labels": "|".join("ABCD"[m] if m < 4 else str(m) for m in members),
                            "items": ",".join(str(m + 1) for m in members)})
    return out


def cluster(ivs: Sequence[Interval], distance: int = 1, max_cluster_size: int = 0) -> list[dict]:
    rows = sort_ivs(ivs)
    out, cl, last = [], 1, None
    cur: list[Interval] = []
    for iv in rows:
        if last and iv.start - last <= distance and (not max_cluster_size or len(cur) < max_cluster_size):
            cur.append(iv)
            last = max(last, iv.end)
        else:
            if cur:
                out.append(_emit_cluster(cur, cl, distance))
                cl += 1
            cur, last = [iv], iv.end
    if cur:
        out.append(_emit_cluster(cur, cl, distance))
    return out


def _emit_cluster(rows: list[Interval], cl: int, distance: int) -> dict:
    return {"chrom": rows[0].chrom, "start": min(r.start for r in rows),
            "end": max(r.end for r in rows), "n": len(rows), "cluster": cl,
            "names": "|".join(r.name for r in rows),
            "total_bp": sum(r.length for r in rows),
            "span": max(r.end for r in rows) - min(r.start for r in rows)}


def links(ivs: Sequence[Interval], distance: int = 1) -> list[dict]:
    cl = cluster(ivs, distance)
    out = []
    for c in cl:
        for i in range(c["n"]):
            for j in range(i + 1, c["n"]):
                nm = c["names"].split("|")
                out.append({"chrom": c["chrom"], "name1": nm[i], "start1": c["start"],
                            "end1": c["end"], "name2": nm[j]})
    return out


def getfasta(genome: dict[str, list[str]], ivs: Sequence[Interval], strand: str = "same",
             name_mode: str = "bed") -> list[Interval]:
    """Intervals carrying extracted sequence in ``extra`` (used by FASTA tools)."""
    out = []
    for iv in sort_ivs(ivs):
        seq = "".join(genome.get(iv.chrom, []))[iv.start:iv.end]
        if strand == "reverse" or (strand == "same" and iv.strand == "-"):
            from chroma_titan.core.seq import revcomp

            seq = revcomp(seq)
        label = f"{iv.chrom}:{iv.start}-{iv.end}" if name_mode == "bed" else iv.name
        out.append(Interval(iv.chrom, iv.start, iv.end, label, iv.score, iv.strand, [seq]))
    return out


def nuc_content(ivs: Sequence[Interval], genome: dict[str, list[str]]) -> list[dict]:
    from chroma_titan.core.seq import gc_content

    out = []
    for iv in sort_ivs(ivs):
        seq = "".join(genome.get(iv.chrom, []))[iv.start:iv.end]
        n = len(seq) or 1
        out.append({"chrom": iv.chrom, "start": iv.start, "end": iv.end, "name": iv.name,
                    "length": len(seq), "numA": seq.count("A"), "numC": seq.count("C"),
                    "numG": seq.count("G"), "numT": seq.count("T"), "numN": seq.count("N"),
                    "GC": round(100 * (seq.count("G") + seq.count("C")) / n, 3),
                    "GC_window": round(gc_content(seq), 3),
                    "AT": round(100 * (seq.count("A") + seq.count("T")) / n, 3)})
    return out


def random_intervals(genome: dict[str, int], n: int = 10, length: int = 100,
                     seed: int = 1, chrom_choice: str | None = None) -> list[Interval]:
    rng = random.Random(seed)
    chroms = [chrom_choice] if chrom_choice else list(genome)
    out = []
    for i in range(n):
        c = rng.choice(chroms)
        s = rng.randint(0, max(0, genome[c] - length))
        out.append(Interval(c, s, s + length, f"random_{i + 1}", 0,
                            rng.choice(["+", "-", "."])))
    return sort_ivs(out)


def binnify(genome: dict[str, int], binsize: int = 1000) -> list[Interval]:
    out = []
    for c, size in genome.items():
        for s in range(0, size, binsize):
            out.append(Interval(c, s, min(size, s + binsize), f"{c}_{s // binsize}"))
    return out


def intersect_intervals(a: Interval, b: Interval) -> int:
    return max(0, min(a.end, b.end) - max(a.start, b.start))


def to_table(ivs: Sequence[Interval]) -> "pd.DataFrame":  # noqa: F821
    import pandas as pd

    cols = ["chrom", "start", "end", "name", "score", "strand"]
    rows = [dict(zip(cols, [iv.chrom, iv.start, iv.end, iv.name, iv.score, iv.strand]))
            for iv in ivs]
    return pd.DataFrame(rows, columns=cols)


def bed_to_gff(ivs: Sequence[Interval], source: str = "chroma") -> list:
    from chroma_titan.core.io import GffRow

    return [GffRow(iv.chrom, source, "feature", iv.start + 1, iv.end, iv.score,
                   iv.strand, ".", {"ID": iv.name}) for iv in ivs]


def gff_to_intervals(rows: Iterable) -> list[Interval]:
    return [Interval(r.seqid, r.start - 1, r.end,
                     r.get("Name") or r.get("gene_id") or r.get("ID") or ".",
                     r.score, r.strand, [r.type]) for r in rows]


def tidy(ivs: Sequence[Interval], max_extra: int = 6) -> list[Interval]:
    out = []
    for iv in ivs:
        extra = [x for x in iv.extra[:max_extra] if x not in ("", None)]
        out.append(Interval(iv.chrom, int(iv.start), int(iv.end), iv.name or ".",
                           float(iv.score or 0), iv.strand if iv.strand in "+-." else ".",
                           extra))
    return out


def blacklist_filter(ivs: Sequence[Interval], blacklist: Sequence[Interval],
                     max_frac: float = 0.5) -> list[Interval]:
    keep = []
    for iv in ivs:
        ov = sum(intersect_intervals(iv, b) for b in blacklist if b.chrom == iv.chrom)
        if iv.length and ov / iv.length <= max_frac:
            keep.append(iv)
    return keep


def remove_overlaps(ivs: Sequence[Interval], min_overlap: int = 1) -> list[Interval]:
    """Greedy removal of the shorter interval in each overlapping pair."""
    out, taken = [], []
    by: dict[str, list[tuple[int, int]]] = defaultdict(list)
    for iv in sorted(ivs, key=lambda x: -x.length):
        spans = by.get(iv.chrom, [])
        if any(iv.start < e and iv.end > s for s, e in spans if e - s >= min_overlap):
            continue
        spans.append((iv.start, iv.end))
        out.append(iv)
        taken.append(iv.name)
    return sort_ivs(out)


def scale_bedgraph(vals: list[tuple[str, int, int, float]], factor: float) -> list:
    return [(c, s, e, v * factor) for c, s, e, v in vals]


def window_counts(ivs: Sequence[Interval], size: int = 1000) -> list[dict]:
    acc: dict[tuple[str, int], int] = Counter()
    for iv in ivs:
        for b in range(iv.start // size, iv.end // size + 1):
            acc[(iv.chrom, b)] += 1
    return [{"chrom": c, "start": b * size, "end": (b + 1) * size, "count": n}
            for (c, b), n in sorted(acc.items())]


def tss_distance_profile(ivs: Sequence[Interval], tss: Sequence[Interval],
                         bin_size: int = 100) -> list[dict]:
    prof: Counter = Counter()
    for t in tss:
        for iv in ivs:
            if iv.chrom != t.chrom:
                continue
            mid = iv.start + iv.length // 2
            d = mid - t.start
            if abs(d) <= 10000:
                prof[d // bin_size] += 1
    return [{"bin": k * bin_size, "count": n} for k, n in sorted(prof.items())]


def tss_list(gff_rows) -> list[Interval]:
    out = []
    for r in gff_rows:
        if r.type.lower() in ("gene", "mrna", "transcript"):
            pos = r.start - 1 if r.strand == "+" else r.end
            out.append(Interval(r.seqid, pos, pos + 1,
                                r.get("Name") or r.get("gene") or r.get("ID") or ".",
                                0, r.strand))
    return out


def genes_from_gff(gff_rows) -> list[Interval]:
    rows = [r for r in gff_rows if r.type.lower() in ("gene", "mrna", "transcript")]
    return [Interval(r.seqid, r.start - 1, r.end,
                     r.get("Name") or r.get("gene") or r.get("ID") or ".", r.score, r.strand)
            for r in rows]


def length_stats(ivs: Sequence[Interval]) -> dict:
    ls = [iv.length for iv in ivs]
    if not ls:
        return {"count": 0}
    ls_sorted = sorted(ls)
    n = len(ls_sorted)
    return {"count": n, "total_bp": sum(ls), "min": ls_sorted[0], "max": ls_sorted[-1],
            "mean": round(sum(ls) / n, 3), "median": ls_sorted[n // 2],
            "N50": _n50(ls), "L50": _l50(ls),
            "chroms": len({iv.chrom for iv in ivs}),
            "frac_of_genome": 0.0}


def _n50(ls: list[int]) -> int:
    tot = sum(ls)
    acc = 0
    for L in sorted(ls, reverse=True):
        acc += L
        if acc >= tot / 2:
            return L
    return 0


def _l50(ls: list[int]) -> int:
    tot = sum(ls)
    acc, n = 0, 0
    for L in sorted(ls, reverse=True):
        acc += L
        n += 1
        if acc >= tot / 2:
            return n
    return 0


def _index(ivs: Sequence[Interval]) -> dict[str, list[Interval]]:
    by: dict[str, list[Interval]] = defaultdict(list)
    for iv in ivs:
        by[iv.chrom].append(iv)
    return by


def read_genome(src) -> tuple[dict[str, list[str]], dict[str, int]]:
    """FASTA -> (per-chrom sequence strings, sizes)."""
    from chroma_titan.core.io import as_text, parse_fasta

    recs = parse_fasta(as_text(src))
    seqs = {r.id: r.seq.upper() for r in recs}
    return seqs, {k: len(v) for k, v in seqs.items()}


def parse_genome_file(src, mode: str = "sizes") -> dict:
    """``genome.txt`` (chrom<tab>size) or a FASTA/fai file -> {chrom: size}."""
    from chroma_titan.core.io import as_text

    out: dict[str, int] = {}
    for ln in as_text(src).splitlines():
        if not ln.strip() or ln.startswith(("#", "track")):
            continue
        p = ln.split("\t")
        if len(p) >= 2:
            try:
                out[p[0]] = int(float(p[1]))
            except ValueError:
                continue
        else:
            p2 = ln.split()
            if len(p2) >= 2:
                try:
                    out[p2[0]] = int(float(p2[1]))
                except ValueError:
                    continue
    return out


def read_fasta_genome(src):
    return read_genome(src)


def bedgraph_to_bed(vals, threshold: float = 0.0) -> list[Interval]:
    return [Interval(c, s, e, f"v{_fmt(v)}", v) for c, s, e, v in vals if v >= threshold]


def _fmt(v: float) -> str:
    return str(int(v)) if float(v) == int(v) else f"{v:g}"


def rescale_bedgraph(vals, scale: float = 1.0, offset: float = 0.0):
    return [(c, s, e, v * scale + offset) for c, s, e, v in vals]


def restrict_bedgraph_to_regions(vals, ivs: Sequence[Interval]):
    keep = []
    idx = _index(ivs)
    for c, s, e, v in vals:
        for iv in idx.get(c, []):
            a, b = max(s, iv.start), min(e, iv.end)
            if b > a:
                keep.append((c, a, b, v))
                break
    return keep


def summarize_intervals(ivs: Sequence[Interval]) -> dict:
    st = length_stats(ivs)
    st["mean_per_bp"] = round(st.get("total_bp", 0) / max(1, st.get("chroms", 1)), 2)
    return st


def distance_between(a: Interval, b: Interval) -> int:
    if a.chrom != b.chrom:
        return -1
    if a.end <= b.start:
        return b.start - a.end
    if b.end <= a.start:
        return a.start - b.end
    return 0
