"""Sequence motifs: PWM math, scanning, MEME-style discovery, TOMTOM-style compare."""

from __future__ import annotations

import math
import re
from collections import Counter, defaultdict

from chroma_titan.core.seq import clean, revcomp

DNA = "ACGT"


def pfm_to_pwm(pfm: list[list[float]], pseudocount: float = 1.0,
               bg: list[float] | None = None) -> list[list[float]]:
    """Counts -> log2 odds ratio per position/base."""
    bg = bg or [0.25] * 4
    out = []
    for col in pfm:
        tot = sum(col) + pseudocount * len(col)
        row = []
        for i, c in enumerate(col):
            f = (c + pseudocount) / tot
            row.append(math.log2(f / bg[i % len(bg)]))
        out.append(row)
    return out


def pfm_to_icm(pfm: list[list[float]]) -> list[list[float]]:
    """Counts -> information content per column (bits)."""
    out = []
    for col in pfm:
        tot = sum(col) or 1
        freq = [c / tot for c in col]
        ic = 2.0 + sum(f * math.log2(f) for f in freq if f > 0)
        out.append([round(f * ic, 5) for f in freq])
    return out


def pfm_to_consensus(pfm: list[list[float]], threshold: float = 0.0) -> str:
    """IUPAC consensus of a PFM (threshold = min relative frequency to keep a base)."""
    out = []
    for col in pfm:
        tot = sum(col) or 1
        keep = sorted(DNA[i] for i, c in enumerate(col) if c / tot >= max(threshold, 0.25))
        out.append(_iupac(keep))
    return "".join(out)


def _iupac(bases: list[str]) -> str:
    key = "".join(sorted(set(bases)))
    table = {"A": "A", "C": "C", "G": "G", "T": "T", "AC": "M", "AG": "R", "AT": "W",
             "CG": "S", "CT": "Y", "GT": "K", "ACG": "V", "ACT": "H", "AGT": "D",
             "CGT": "B", "ACGT": "N"}
    return table.get(key, "N")


def count_matrix_from_sites(sites: list[str], alphabet: str = DNA) -> list[list[float]]:
    if not sites:
        return []
    L = max(len(s) for s in sites)
    rows = [[0.0] * len(alphabet) for _ in range(L)]
    for s in sites:
        for i, ch in enumerate(s):
            if ch in alphabet:
                rows[i][alphabet.index(ch)] += 1
    return rows


def read_pfm(text: str, alphabet: str = DNA) -> list[list[float]]:
    """Accepts JASPAR-style (4 rows) or MEME/MotifBase-style (one row per position)."""
    from chroma_titan.core.io import as_text

    rows: list[list[float]] = []
    per_base: dict[str, list[float]] = defaultdict(list)
    for ln in as_text(text).splitlines():
        ln = ln.strip().replace("[", " ").replace("]", " ")
        if not ln or ln[0] in ">#{}" or ln.startswith("MOTIF") or ln.startswith("letter"):
            continue
        parts = ln.replace(",", " ").split()
        head = parts[0].strip(":").upper()
        if head in alphabet and len(parts) > 1 and all(_num(x) for x in parts[1:]):
            per_base[head] = [float(x) for x in parts[1:]]
            continue
        if all(_num(x) for x in parts[:len(alphabet)]):
            rows.append([float(x) for x in parts[:len(alphabet)]])
    if per_base and not rows:
        n = max(len(v) for v in per_base.values())
        rows = [[per_base.get(a, [0] * n)[i] for a in alphabet] for i in range(n)]
    return rows


def _num(x: str) -> bool:
    try:
        float(x)
        return True
    except ValueError:
        return False


def write_pfm(pfm: list[list[float]], name: str = "motif", alphabet: str = DNA,
              style: str = "jaspar") -> str:
    if style == "meme":
        head = ["MOTIF " + name, f"letter-probability matrix: alength= {len(alphabet)} "
                                 f"w= {len(pfm)}"]
        body = ["\t".join(f"{v:.5f}" for v in row) for row in pfm]
        return "\n".join(head + body) + "\n"
    width = max((len(r) for r in pfm), default=0)
    lines = [f">{name}"]
    for i, a in enumerate(alphabet):
        vals = " ".join(f"{int(row[i]) if row[i] == int(row[i]) else row[i]:g}" for row in pfm)
        lines.append(f"[{vals}]" if style == "pfm" else f"{a}  [" + vals + " ]")
    del width
    return "\n".join(lines) + "\n"


def scan_sequence(seq: str, pwm: list[list[float]], threshold: float = 0.0,
                  strand: str = "both", min_hits: int = 0,
                  pfm: list[list[float]] | None = None) -> list[dict]:
    """Score every window with the log-odds matrix; report hits over threshold."""
    from chroma_titan.core.seq import back_transcribe

    s = back_transcribe(clean(seq, "ACGTN"))
    w = len(pwm)
    if w < 1 or len(s) < w:
        return []
    rc_pwm = [row[::-1] for row in pwm][::-1]
    hits = []
    for start in range(len(s) - w + 1):
        win = s[start:start + w]
        if "N" in win:
            continue
        for strand_name, matrix, shown in (("+", pwm, win), ("-", rc_pwm, revcomp(win))):
            if strand == "+" and strand_name == "-":
                continue
            if strand == "-" and strand_name == "+":
                continue
            score = sum(matrix[i][DNA.index(ch)] for i, ch in enumerate(shown) if ch in DNA)
            if score >= threshold:
                hits.append({"start": start, "end": start + w, "strand": strand_name,
                             "score": round(score, 4), "match": win if strand_name == "+" else shown,
                             "position": start + 1})
    hits.sort(key=lambda h: -h["score"])
    if min_hits:
        hits = hits[:min_hits]
    return hits


def pwm_score_max(pwm: list[list[float]]) -> float:
    return sum(max(row) for row in pwm) if pwm else 0.0


def pwm_score_min(pwm: list[list[float]]) -> float:
    return sum(min(row) for row in pwm) if pwm else 0.0


def scan_with_pfm(seq: str, pfm: list[list[float]], percentile: float = 80.0,
                  strand: str = "both") -> list[dict]:
    pwm = pfm_to_pwm(pfm)
    hi = pwm_score_max(pwm)
    lo = pwm_score_min(pwm)
    thr = lo + (hi - lo) * percentile / 100
    return scan_sequence(seq, pwm, thr, strand)


def total_information(pfm: list[list[float]]) -> float:
    return round(sum(r[0] for r in pfm_to_icm(pfm)) * 4, 4) if pfm else 0.0


def bits_per_position(pfm: list[list[float]]) -> list[float]:
    return [round(sum(v), 5) for v in pfm_to_icm(pfm)]


def motif_compare(a: list[list[float]], b: list[list[float]], mode: str = "both",
                  min_overlap: int = 6, correction: bool = True) -> list[dict]:
    """TOMTOM-lite: sliding Pearson correlation of two PFMs (all offsets/strands)."""
    def colnorm(m):
        out = []
        for col in m:
            tot = sum(col) or 1
            mu = tot / len(col)
            v = math.sqrt(sum((x - mu) ** 2 for x in col) / len(col)) or 1e-9
            out.append([(x - mu) / v for x in col])
        return out

    na, nb = colnorm(a), colnorm(b)
    results = []
    for strand, target in (("+", nb), ("-", [row[::-1] for row in nb][::-1])):
        for offset in range(-len(target) + min_overlap, len(na) - min_overlap + 1):
            rows = []
            for i in range(len(na)):
                j = i - offset
                if 0 <= j < len(target):
                    rows.append((na[i], target[j]))
            if len(rows) < min_overlap:
                continue
            num = den1 = den2 = 0.0
            for x, y in rows:
                for k in range(len(x)):
                    num += x[k] * y[k]
                    den1 += x[k] * x[k]
                    den2 += y[k] * y[k]
            r = num / math.sqrt(den1 * den2) if den1 and den2 else 0.0
            n_pos = len(rows) * len(na[0])
            results.append({"mode": mode, "strand": strand, "offset": offset,
                            "overlapping_columns": len(rows),
                            "pearson_r": round(r, 5),
                            "euclidean_distance": round(math.sqrt(sum(
                                (x[k] - y[k]) ** 2 for x, y in rows for k in range(len(x)))), 4),
                            "p_value": _pearson_p(r, n_pos)})
    results.sort(key=lambda d: -abs(d["pearson_r"]))
    return results


def _pearson_p(r: float, n: int) -> float:
    from chroma_titan.core import stats

    if n < 3:
        return float("nan")
    r = max(-0.999999, min(0.999999, r))
    t = r * math.sqrt((n - 2) / (1 - r * r))
    return stats.t_pvalue(t, n - 2)


def meme_discover(seqs: list[str], width: int = 8, sites: str = "zoops",
                  n_iter: int = 30, n_init: int = 8, seed: int = 1,
                  bg: list[float] | None = None) -> list[dict]:
    """MEME-lite: expectation-maximisation motif finder (ZOOPS/OOPS/THRESHOLD)."""
    import random

    rng = random.Random(seed)
    rows = [clean(s, "ACGTN") for s in seqs if len(clean(s, "ACGTN")) >= width]
    if not rows:
        return []
    bg = bg or [0.25] * 4
    best = None
    for init in range(n_init):
        starts = [rng.randrange(len(s) - width + 1) for s in rows]
        counts = [[0.0] * 4 for _ in range(width)]
        for s, st in zip(rows, starts):
            for i, ch in enumerate(s[st:st + width]):
                if ch in DNA:
                    counts[i][DNA.index(ch)] += 1
        model = [[(x + 0.25) / (sum(col) + 1.0) for x in col] for col in counts]
        ll_old = -1e18
        for _ in range(n_iter):
            # E-step: posterior of each position being a site (ZOOPS)
            new_counts = [[0.0] * 4 for _ in range(width)]
            total_ll = 0.0
            npk = 0.0
            for s in rows:
                scores = []
                for st in range(len(s) - width + 1):
                    win = s[st:st + width]
                    sc = 0.0
                    ok = True
                    for i, ch in enumerate(win):
                        if ch not in DNA:
                            ok = False
                            break
                        sc += math.log(model[i][DNA.index(ch)] / bg[DNA.index(ch)])
                    scores.append(sc if ok else -1e9)
                if not scores:
                    continue
                mx = max(scores)
                wts = [math.exp(min(500, sc - mx)) for sc in scores]
                z = sum(wts) + (len(s) - width + 1) * 0.0 if sites == "zoops" else sum(wts)
                if sites == "oops" and len(rows) > 1:
                    z = sum(wts)
                p = [x / z for x in wts] if z else []
                total_ll += mx + (math.log(z) if z else 0.0)
                for st, pi in enumerate(p):
                    if pi < 1e-3:
                        continue
                    npk += pi
                    for i, ch in enumerate(s[st:st + width]):
                        if ch in DNA:
                            new_counts[i][DNA.index(ch)] += pi
            model = [[(x + 0.25) / (sum(col) + 1.0) for x in col] for col in new_counts] \
                if npk > 0.5 else model
            if abs(total_ll - ll_old) < 1e-4:
                break
            ll_old = total_ll
        pfm = [[round(c, 3) for c in col] for col in new_counts]
        score = sum(max(col) for col in pfm_to_pwm(pfm))
        if best is None or score > best["score"]:
            best = {"score": round(score, 4), "log_likelihood": round(ll_old, 3),
                    "width": width, "sites": sites, "pfm": pfm,
                    "consensus": pfm_to_consensus(pfm),
                    "pwm": [[round(v, 4) for v in row] for row in pfm_to_pwm(new_counts)],
                    "ic": total_information(pfm), "n_init": init + 1}
    return [best] if best else []


def sites_for_motif(seq: str, pfm: list[list[float]], threshold: float = 80.0) -> list[str]:
    pwm = pfm_to_pwm(pfm)
    hi = pwm_score_max(pwm)
    thr = hi * threshold / 100
    return [h["match"] for h in scan_sequence(seq, pwm, thr, "+")]


def extract_sites(seqs: list[str], pattern: str, fasta_like: bool = False) -> dict:
    """All occurrences of a consensus/IUPAC pattern, both strands."""
    from chroma_titan.core.seq import iupac_to_regex

    rx = re.compile(iupac_to_regex(pattern))
    sites, rc_sites = [], []
    for s in seqs:
        for m in rx.finditer(clean(s, "ACGTN")):
            sites.append(m.group())
            rc_sites.append(revcomp(m.group()))
    return {"sites": sites, "sites_with_gaps": sites, "consensus": pfm_to_consensus(
        count_matrix_from_sites(sites)) if sites else "",
        "count_matrix": count_matrix_from_sites(sites), "n_sites": len(sites),
        "n_distinct": len(set(sites)), "reverse_complement_sites": rc_sites,
        "fasta_like": fasta_like}


def central_enrichment(seqs: list[str], pattern: str, centre_frac: float = 0.2) -> dict:
    """CEMS-lite: fraction of motif hits falling in the central region."""
    from chroma_titan.core.seq import iupac_to_regex

    rx = re.compile(iupac_to_regex(pattern))
    inside = total = 0
    for s in seqs:
        seq = clean(s, "ACGTN")
        half = len(seq) / 2
        win = len(seq) * centre_frac / 2
        for m in rx.finditer(seq) or []:
            total += 1
            mid = (m.start() + m.end()) / 2
            inside += abs(mid - half) <= win
    return {"n_hits": total, "n_central": inside,
            "central_enrichment": round(inside / total, 4) if total else 0.0,
            "expected_by_chance": round(centre_frac, 4),
            "fold_enrichment": round((inside / total) / centre_frac, 3) if total and centre_frac else 0.0}


def motif_enrichment_count(seqs: list[str], pattern: str) -> dict:
    from chroma_titan.core.seq import iupac_to_regex

    rx = re.compile(iupac_to_regex(pattern))
    hits = sum(len(rx.findall(clean(s, "ACGTN"))) for s in seqs)
    n_with = sum(1 for s in seqs if rx.search(clean(s, "ACGTN")))
    return {"pattern": pattern, "total_hits": hits, "sequences_with_hit": n_with,
            "sequences": len(seqs), "fraction_with_hit": round(n_with / max(1, len(seqs)), 4)}


def jaccard_motifs(a: list[list[float]], b: list[list[float]]) -> float:
    """Similarity by overlapping consensus letters (JASPAR-style)."""
    ca, cb = pfm_to_consensus(a), pfm_to_consensus(b)
    n = min(len(ca), len(cb))
    if not n:
        return 0.0
    best = 0.0
    for off in range(-n + 1, n):
        match = tot = 0
        for i in range(n):
            j = i - off
            if 0 <= j < len(cb):
                tot += 1
                match += ca[i] == cb[j]
        best = max(best, match / tot if tot else 0.0)
    return round(best, 4)
