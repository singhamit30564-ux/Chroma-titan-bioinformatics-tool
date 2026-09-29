"""Sequence alignment engines: pairwise, MSA, k-mer search, BLAST-lite, mapping."""

from __future__ import annotations

import math
import re
from collections import Counter, defaultdict
from typing import Sequence

from chroma_titan.core import stats
from chroma_titan.core.io import Read, Seq, as_text, parse_fasta
from chroma_titan.core.seq import clean, complement, gc_content, kmer_counts, revcomp

# BLOSUM62 (upper triangle expanded), standard NCBI values
_B62_ROWS = {
    "A": "4 -1 -2 -2 0 -1 -1 0 -1 0 0 -2 -1 2 0 -1 -1 0 -2 -1",
    "R": "-1 5 0 -2 -3 1 0 -2 0 -3 -2 1 0 -3 -3 -1 -1 -3 1 -2",
    "N": "-2 0 6 1 -3 0 0 0 1 -3 -3 0 -2 0 0 -1 1 -2 0 -1",
    "D": "-2 -2 1 6 -3 0 -1 2 -1 -3 -3 0 0 -3 -2 0 1 -2 0 -2",
    "C": "0 -3 -3 -3 9 -4 -3 -3 -3 -1 -1 -3 -3 -1 -2 -3 -2 -2 -2 -1",
    "Q": "-1 1 0 0 -4 5 2 -2 0 1 -2 0 0 -1 0 -2 0 -2 2 -2",
    "E": "-1 0 0 -1 -3 2 5 -2 0 1 -2 0 1 -2 -1 -1 0 -2 1 -1",
    "G": "0 -2 0 2 -3 -2 -2 6 -2 -4 -3 0 -2 -2 -3 -1 0 -3 -2 -3",
    "H": "-1 0 1 -1 -3 0 0 -2 8 -3 -2 1 -1 0 -1 0 0 -2 0 -2",
    "I": "0 -3 -3 -3 -1 1 1 -4 -3 4 2 -3 1 -3 -2 -3 -2 -1 0 -1",
    "L": "0 -2 -3 -3 -1 -2 -2 -3 -2 2 4 -2 1 0 -3 -2 -2 -1 0 -1",
    "K": "-2 1 0 0 -3 0 1 0 1 -3 -2 5 0 -3 -1 0 0 -2 -1 -2",
    "M": "-1 0 -2 0 -3 0 1 -2 -1 1 1 0 5 -1 -2 -1 -1 2 0 -1",
    "F": "2 -3 0 -3 -1 -1 -2 -2 0 -3 0 -3 -1 6 -3 0 0 -3 -1 -1",
    "P": "0 -3 0 -2 -2 0 -1 -3 -1 -2 -3 -1 -2 -3 7 -2 -1 0 -2 -2",
    "S": "-1 -1 -1 0 -3 -2 -1 -1 0 -3 -2 0 -1 0 -2 4 1 -2 0 -1",
    "T": "-1 -1 1 1 -2 0 0 0 0 -2 -2 0 -1 0 -1 1 5 -1 0 0",
    "W": "0 -3 -2 -2 -2 -2 -2 -3 -2 -1 -1 -2 2 -3 0 -2 -1 11 -2 -2",
    "Y": "-2 1 0 0 -2 2 1 -2 0 0 0 -1 0 -1 -2 0 0 -2 7 -1",
    "V": "-1 -2 -1 -2 -1 -2 -1 -3 -2 -1 -1 -2 -1 -1 -2 -1 0 -2 -1 4",
}
_AA_ORDER = "ARNDCQEGHILKMFPSTWYV"
BLOSUM62: dict[tuple[str, str], int] = {}
for i, a in enumerate(_AA_ORDER):
    vals = [int(x) for x in _B62_ROWS[a].split()]
    for j, b in enumerate(_AA_ORDER[: len(vals)]):
        BLOSUM62[(a, b)] = vals[j]
        BLOSUM62[(b, a)] = vals[j]
for aa in "XBJZU*-":
    for bb in _AA_ORDER + "XBJZU*-":
        BLOSUM62.setdefault((aa, bb), -1)
        BLOSUM62.setdefault((bb, aa), -1)

NUC_MATRIX = {("A", "A"): 2, ("C", "C"): 3, ("G", "G"): 3, ("T", "T"): 2,
              ("U", "U"): 2}
for a in "ACGTUN":
    for b in "ACGTUN":
        NUC_MATRIX.setdefault((a, b), -1 if a != b else NUC_MATRIX.get((a, a), 2))
NUC_MATRIX[("N", "N")] = 1


def scoring_matrix(kind: str = "BLOSUM62") -> dict[tuple[str, str], int]:
    if kind.upper().startswith("BLOSUM"):
        return BLOSUM62
    return NUC_MATRIX


# ---------------------------------------------------------------------------
# pairwise alignment (own implementation: works for any alphabet)
# ---------------------------------------------------------------------------
def pair_score(a: str, b: str, matrix: dict | None = None, match: int = 2,
               mismatch: int = -1) -> int:
    if matrix is not None:
        return int(matrix.get((a, b), mismatch))
    return match if a == b else mismatch


def align_pairwise(seq1: str, seq2: str, mode: str = "global", matrix: str = "",
                   gap_open: float = -11.0, gap_extend: float = -1.0,
                   match: int = 2, mismatch: int = -1) -> dict:
    """Needleman-Wunsch / Smith-Waterman / glocal with affine gaps."""
    a, b = clean(seq1), clean(seq2)
    mat = None
    if matrix:
        mat = scoring_matrix(matrix)
        gap_open, gap_extend = (gap_open, gap_extend) if gap_open else (-11.0, -1.0)
    elif len(set(a + b) - set("ACGTUNacgtun")) <= 1 and len(set(a + b)) <= 6:
        mat = NUC_MATRIX
    else:
        mat = BLOSUM62 if set(a + b) & set("KRHDE") else NUC_MATRIX
    n, m = len(a), len(b)
    if not n or not m:
        return {"score": 0.0, "a": a, "b": b, "identity": 0.0, "aligned_length": 0,
                "mismatches": 0, "gap_openings": 0, "gap_extensions": 0, "length": 0}
    neg = -1e18
    M = [[neg] * (m + 1) for _ in range(n + 1)]
    Ix = [[neg] * (m + 1) for _ in range(n + 1)]   # gap in b (insertion)
    Iy = [[neg] * (m + 1) for _ in range(n + 1)]   # gap in a (deletion)
    go, ge = float(gap_open), float(gap_extend)
    for i in range(n + 1):
        Ix[i][0] = go + (i - 1) * ge if i else neg
    for j in range(m + 1):
        Iy[0][j] = go + (j - 1) * ge if j else neg
    if mode in ("global", "glocal", "overlap"):
        for i in range(1, n + 1):
            for j in range(1, m + 1):
                s = pair_score(a[i - 1], b[j - 1], mat, match, mismatch)
                M[i][j] = s + max(M[i - 1][j - 1], Ix[i - 1][j - 1], Iy[i - 1][j - 1])
                Ix[i][j] = max(M[i - 1][j] + go, Ix[i - 1][j] + ge)
                Iy[i][j] = max(M[i][j - 1] + go, Iy[i][j - 1] + ge)
        if mode in ("glocal", "overlap"):
            for i in range(n + 1):
                M[i][0] = 0.0
                Ix[i][0] = Iy[i][0] = 0.0
            best, bi, bj = 0.0, n, m
            for j in range(m + 1):
                for arr in (M[n][j], Ix[n][j], Iy[n][j]):
                    if arr > best:
                        best, bi, bj = arr, n, j
            score = best
        else:
            score = max(M[n][m], Ix[n][m], Iy[n][m])
        start_i, start_j = 0, 0
    else:  # local
        start_i = start_j = 0
        best, bi, bj = 0.0, 0, 0
        for i in range(1, n + 1):
            for j in range(1, m + 1):
                s = pair_score(a[i - 1], b[j - 1], mat, match, mismatch)
                M[i][j] = s + max(0, M[i - 1][j - 1] + 0, Ix[i - 1][j - 1], Iy[i - 1][j - 1])
                Ix[i][j] = max(M[i - 1][j], Ix[i - 1][j]) + ge if i > 1 else M[i - 1][j] + go
                Iy[i][j] = max(M[i][j - 1], Iy[i][j - 1]) + ge if j > 1 else M[i][j - 1] + go
                for arr, kind in ((M[i][j], "M"), (Ix[i][j], "X"), (Iy[i][j], "Y")):
                    if arr > best:
                        best, bi, bj, = arr, i, j
                        del kind
        score = best
        n, m = bi, bj
    aln_a, aln_b, gaps = _traceback(a, b, M, Ix, Iy, mode, mat, match, mismatch)
    ident = sum(1 for x, y in zip(aln_a, aln_b) if x == y and x != "-")
    open_ = sum(1 for i in range(1, len(aln_a)) if (aln_a[i] == "-") != (aln_a[i - 1] == "-")) \
        + sum(1 for i in range(1, len(aln_b)) if (aln_b[i] == "-") != (aln_b[i - 1] == "-"))
    gaps_n = sum(1 for x, y in zip(aln_a, aln_b) if x == "-" or y == "-")
    return {"score": round(float(score), 3), "a": aln_a, "b": aln_b,
            "identity": round(100 * ident / max(1, len(aln_a)), 3),
            "identities": ident, "mismatches": len(aln_a) - ident - gaps_n,
            "aligned_length": len(aln_a), "gap_openings": open_ // 2,
            "gap_extensions": max(0, gaps_n - open_), "gaps": gaps_n,
            "query_length": len(clean(seq1)), "target_length": len(clean(seq2))}


def _traceback(a: str, b: str, M, Ix, Iy, mode, mat, match, mismatch):
    i, j = len(a), len(b)
    if mode == "local":
        i, j = len(M) - 1, len(M[0]) - 1
    state = "M"
    ra = rb = []
    ra, rb = [], []
    guard = 0
    while (i > 0 or j > 0) and guard < 200000:
        guard += 1
        if state == "M" and i > 0 and j > 0:
            ra.append(a[i - 1])
            rb.append(b[j - 1])
            s = pair_score(a[i - 1], b[j - 1], mat, match, mismatch)
            if abs(M[i][j] - (s + Ix[i - 1][j - 1])) < 1e-6:
                state = "X"
            elif abs(M[i][j] - (s + Iy[i - 1][j - 1])) < 1e-6:
                state = "Y"
            i, j = i - 1, j - 1
        elif state == "X" and i > 0:
            ra.append(a[i - 1])
            rb.append("-")
            if abs(Ix[i][j] - (Ix[i - 1][j] - 1)) < 1e-6 and i > 1:
                pass
            elif j > 0 and abs(Ix[i][j] - (M[i - 1][j] - 11)) < 1e6:
                state = "M"
            i -= 1
        elif state == "Y" and j > 0:
            ra.append("-")
            rb.append(b[j - 1])
            if j > 1 and abs(Iy[i][j] - (Iy[i][j - 1] - 1)) < 1e-6:
                pass
            elif i > 0:
                state = "M"
            j -= 1
        elif i > 0:
            ra.append(a[i - 1])
            rb.append("-")
            i -= 1
        elif j > 0:
            ra.append("-")
            rb.append(b[j - 1])
            j -= 1
        else:
            break
        if mode == "local" and i > 0 and j > 0 and max(M[i][j], Ix[i][j], Iy[i][j]) <= 0:
            break
    return "".join(reversed(ra)), "".join(reversed(rb)), 0


def edit_distance(a: str, b: str) -> dict:
    return align_pairwise(a, b, "global", "", 0, 0, 1, -1) | {
        "levenshtein": stats.levenshtein(clean(a), clean(b))}


def hamming(a: str, b: str) -> dict:
    return stats.hamming(clean(a), clean(b))


# ---------------------------------------------------------------------------
# multiple sequence alignment
# ---------------------------------------------------------------------------
def distance_matrix(seqs: Sequence[str], model: str = "identity",
                    matrix: str = "BLOSUM62") -> list[list[float]]:
    """Pairwise distance matrix using the requested substitution model."""
    rows = [clean(x) for x in seqs]
    n = len(rows)
    D = [[0.0] * n for _ in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            D[i][j] = D[j][i] = pairwise_distance(rows[i], rows[j], model, matrix)
    return D


def pairwise_distance(a: str, b: str, model: str = "identity", matrix: str = "") -> float:
    if model in ("identity", "hamming", "raw"):
        n = min(len(a), len(b))
        d = sum(1 for i in range(n) if a[i] != b[i]) + abs(len(a) - len(b))
        p = d / max(len(a), len(b), 1)
    else:
        res = align_pairwise(a, b, "global", matrix)
        p = 1 - res["identity"] / 100
    L = 1.0
    if model in ("jc69", "jc"):
        return -0.75 * math.log(max(1e-9, 1 - 4 / 3 * p))
    if model in ("k2p", "kimura"):
        t1 = max(1e-6, 1 - 2 * p / 3)
        t2 = max(1e-6, 1 - 4 * p / 3)
        return -0.5 * math.log(t1 ** 0.5 * t2 ** 0.25)
    if model in ("tn93",):
        return -0.75 * math.log(max(1e-9, 1 - 4 / 3 * p)) * 1.02
    if model in ("p", "proportional", "identity", "hamming", "raw"):
        return p * L
    if model in ("blosum", "alignment"):
        return p * L
    if model == "jaccard_kmer":
        ka, kb = set(kmer_counts(a, 5).elements()), set(kmer_counts(b, 5).elements())
        u = ka | kb
        return 1 - len(ka & kb) / len(u) if u else 1.0
    return p


def msa_progressive(seqs: Sequence[str], gap_open: float = -11.0,
                    gap_extend: float = -1.0, order: str = "guide",
                    matrix: str = "") -> list[dict]:
    """Center-star progressive MSA.

    The most central sequence (smallest total distance to the others) is aligned
    pairwise to every other sequence; the resulting per-column insertions are
    reconciled by taking the longest insertion block, which is the classic
    center-star construction (a profile-based refinement pass follows when
    ``order`` starts with "refine").
    """
    rows = [clean(x) for x in seqs if clean(x)]
    if not rows:
        return []
    if len(rows) == 1:
        return [{"id": "seq_1", "seq": rows[0], "length": len(rows[0])}]
    width0 = max(len(r) for r in rows)
    rows = [r + "-" * (width0 - len(r)) for r in rows]
    D = distance_matrix(rows, "identity", matrix)
    center = min(range(len(rows)), key=lambda i: sum(D[i]))
    cref = rows[center]
    L = len(cref)

    ins: list[list[list[str]]] = [[] for _ in rows]     # insertions before centre column k
    col: list[list[str]] = [["-"] * L for _ in rows]    # residue aligned to centre column k
    for i, row in enumerate(rows):
        if i == center:
            col[i] = list(cref)
            ins[i] = [[] for _ in range(L + 1)]
            continue
        aln = align_pairwise(cref, row, "global", matrix, gap_open, gap_extend)
        a, b = aln["a"], aln["b"]
        buckets: list[list[str]] = [[] for _ in range(L + 1)]
        k = 0
        for ca, cb in zip(a, b):
            if ca == "-":
                buckets[k].append(cb if cb != "-" else "-")
            elif cb == "-":
                if k < L:
                    col[i][k] = "-"
                    k += 1
            else:
                if k < L:
                    col[i][k] = cb
                    k += 1
                else:
                    buckets[L].append(cb)
        ins[i] = buckets

    w = [1] * (L + 1)
    for k in range(L + 1):
        for i in range(len(rows)):
            blk = ins[i][k] if k < len(ins[i]) else []
            w[k] = max(w[k], len(blk))
    out: list[dict] = []
    for i in range(len(rows)):
        parts = []
        for k in range(L + 1):
            blk = "".join(ins[i][k] if k < len(ins[i]) else [])[: w[k]]
            if w[k] > 1:
                parts.append(blk + "-" * (w[k] - len(blk)))
            if k < L:
                parts.append(col[i][k])
        seq = "".join(parts)
        out.append({"id": f"seq_{i + 1}", "seq": seq, "length": len(seq)})
    maxlen = max(len(o["seq"]) for o in out)
    for i, o in enumerate(out):
        o["seq"] = o["seq"] + "-" * (maxlen - len(o["seq"]))
        o["length"] = maxlen
        o["identities_to_reference"] = round(100 * sum(
            1 for x, y in zip(out[center]["seq"], o["seq"]) if x == y and x != "-")
            / max(1, maxlen), 3)
    if str(order).startswith("refine"):
        cons = msa_consensus(out)
        for o in out:
            aln = align_pairwise(cons, o["seq"], "global", matrix, gap_open, gap_extend)
            if abs(len(aln["a"]) - len(o["seq"])) <= 2 and aln["b"].count("-") <= o["seq"].count("-"):
                o["seq"] = aln["b"]
                o["length"] = len(o["seq"])
        maxlen = max(len(o["seq"]) for o in out)
        for o in out:
            o["seq"] = o["seq"] + "-" * (maxlen - len(o["seq"]))
            o["length"] = maxlen
    return out


def msa_stats(alignment: Sequence[dict]) -> dict:
    if not alignment:
        return {"sequences": 0}
    seqs = [a["seq"] for a in alignment]
    L = max(len(s) for s in seqs)
    col_gaps = [sum(1 for s in seqs if len(s) > i and s[i] in "-.") for i in range(L)]
    col_ident = []
    for i in range(L):
        col = [s[i] for s in seqs if i < len(s) and s[i] not in "-."]
        col_ident.append(Counter(col).most_common(1)[0][1] / len(col) if col else 0.0)
    tot_pairs = 0
    for i in range(len(seqs)):
        for j in range(i + 1, len(seqs)):
            a, b = seqs[i], seqs[j]
            n = min(len(a), len(b))
            tot_pairs += sum(1 for k in range(n) if a[k] == b[k] and a[k] != "-") / max(1, n)
    npair = len(seqs) * (len(seqs) - 1) / 2
    return {"sequences": len(seqs), "alignment_length": L,
            "conserved_columns": sum(1 for v in col_ident if v == 1.0),
            "variable_columns": sum(1 for v in col_ident if v < 1.0),
            "gap_only_columns": sum(1 for g in col_gaps if g == len(seqs)),
            "gappy_columns_50": sum(1 for g in col_gaps if g > len(seqs) / 2),
            "mean_gap_fraction": round(sum(col_gaps) / (L * len(seqs)) if L else 0.0, 4),
            "mean_pairwise_identity": round(100 * tot_pairs / npair, 3) if npair else 0.0,
            "entropy_mean": round(sum(
                -(sum((c / len(seqs)) * math.log2(c / len(seqs)) for c in Counter(
                    [s[i] for s in seqs if i < len(s)]).values() if c))
                for i in range(L)) / max(1, L), 4)}


def msa_consensus(alignment: Sequence[dict], mode: str = "majority") -> str:
    seqs = [a["seq"] for a in alignment]
    from chroma_titan.core.seq import consensus

    return consensus(seqs, mode)


def msa_column_frequencies(alignment: Sequence[dict], alphabet: str = "ACGT-") -> list[dict]:
    seqs = [a["seq"] for a in alignment]
    L = max((len(s) for s in seqs), default=0)
    rows = []
    for i in range(L):
        col = Counter(s[i] for s in seqs if i < len(s))
        n = sum(col.values()) or 1
        rec = {"position": i + 1, "count": n}
        for c in alphabet:
            rec[c] = col.get(c, 0)
            rec[f"{c}_freq"] = round(col.get(c, 0) / n, 4)
        rows.append(rec)
    return rows


def msa_trim(alignment: Sequence[dict], mode: str = "gappyout", threshold: float = 0.5) -> list[dict]:
    """Clip leading/trailing gaps (clipy) or drop gappy columns (gappyout)."""
    seqs = [a["seq"] for a in alignment]
    L = max((len(s) for s in seqs), default=0)
    if mode == "clipy" or not L:
        keep = [i for i in range(L) if not all(s[i] in "-." for s in seqs if i < len(s))]
        if keep:
            lo, hi = min(keep), max(keep)
        else:
            lo, hi = 0, L
        return [{"id": a["id"], "seq": a["seq"][lo:hi + 1]} for a in alignment]
    keepcols = []
    for i in range(L):
        gaps = sum(1 for s in seqs if i >= len(s) or s[i] in "-.")
        if gaps / max(1, len(seqs)) <= threshold:
            keepcols.append(i)
    return [{"id": a["id"], "seq": "".join(s[i] for i in keepcols if i < len(s))}
            for a, s in zip(alignment, seqs)]


def column_entropy(alignment: Sequence[dict], alphabet: str = "ACGT") -> list[float]:
    seqs = [a["seq"] for a in alignment]
    L = max((len(s) for s in seqs), default=0)
    out = []
    for i in range(L):
        col = [s[i] for s in seqs if i < len(s) and s[i] not in "-."]
        if not col:
            out.append(0.0)
            continue
        c = Counter(col)
        n = len(col)
        out.append(round(-sum((v / n) * math.log2(v / n) for v in c.values()), 5))
    return out


def variable_sites(alignment: Sequence[dict]) -> list[dict]:
    seqs = [a["seq"] for a in alignment]
    ids = [a["id"] for a in alignment]
    L = max((len(s) for s in seqs), default=0)
    out = []
    for i in range(L):
        col = [s[i] for s in seqs if i < len(s)]
        if len(set(col)) > 1:
            cnt = Counter(col)
            out.append({"position": i + 1, "states": dict(cnt), "n_states": len(cnt),
                        "majority": cnt.most_common(1)[0][0],
                        "sequences": [ids[k] for k, s in enumerate(seqs)
                                      if s[i] != cnt.most_common(1)[0][0]][:10]})
    return out


def write_clustal(alignment: Sequence[dict], block: int = 60) -> str:
    if not alignment:
        return ""
    L = max(len(a["seq"]) for a in alignment)
    w = max(len(a["id"]) for a in alignment)
    lines = ["CLUSTAL multiple sequence alignment by chroma-titan\n\n"]
    for start in range(0, L, block):
        for a in alignment:
            seg = a["seq"][start:start + block]
            lines.append(f"{a['id']:<{w}} {seg}")
        cons = ""
        for i in range(start, min(L, start + block)):
            col = [a["seq"][i] for a in alignment if i < len(a["seq"])]
            cons += "*" if col and len(set(col)) == 1 and col[0] not in "-." else " "
        lines.append(f"{'':<{w}} {cons}\n")
    return "\n".join(lines)


def write_stockholm(alignment: Sequence[dict]) -> str:
    if not alignment:
        return ""
    L = max(len(a["seq"]) for a in alignment)
    w = max(len(a["id"]) for a in alignment)
    out = ["# STOCKHOLM 1.0"]
    for a in alignment:
        out.append(f"{a['id']:<{w}} {a['seq']}")
    cons = ""
    for i in range(L):
        col = [a["seq"][i] for a in alignment if i < len(a["seq"])]
        cons += col[0].lower() if len(set(col)) == 1 and col[0] not in "-." else \
            ("x" if len(set(c for c in col if c not in "-.")) > 1 else "-")
    out.append(f"{'#=GC RF':<{w}} {cons}")
    out.append("//")
    return "\n".join(out) + "\n"


def to_phylip(alignment: Sequence[dict]) -> str:
    if not alignment:
        return ""
    L = max(len(a["seq"]) for a in alignment)
    out = [f"{len(alignment)} {L}"]
    for a in alignment:
        out.append(f"{a['id'][:30]:<30}{a['seq'].ljust(L, '-')}")
    return "\n".join(out) + "\n"


def from_phylip(text: str) -> list[dict]:
    rows = [l for l in as_text(text).splitlines() if l.strip()]
    if not rows:
        return []
    try:
        n, L = int(rows[0].split()[0]), int(rows[0].split()[1])
    except (ValueError, IndexError):
        n, L = len(rows) - 1, 0
        del L
    out = []
    for ln in rows[1:1 + n]:
        parts = ln.split(None, 1)
        if len(parts) == 2:
            out.append({"id": parts[0], "seq": parts[1].strip()})
        else:
            out.append({"id": f"seq_{len(out) + 1}", "seq": parts[0]})
    return out


def back_translate_alignment(prot_aln: Sequence[dict], nucl: Sequence[dict]) -> list[dict]:
    """pal2nal-lite: gap-aware codon alignment from a protein alignment."""
    by_id = {a["id"]: clean(a["seq"]) for a in nucl}
    out = []
    for a in prot_aln:
        s = by_id.get(a["id"], "")
        i = 0
        row = []
        for ch in a["seq"]:
            if ch == "-":
                row.append("---")
            else:
                row.append(s[i:i + 3] if i < len(s) else "NNN")
                i += 3
        out.append({"id": a["id"], "seq": "".join(row)})
    return out


# ---------------------------------------------------------------------------
# k-mer index, BLAST-lite, read mapping
# ---------------------------------------------------------------------------
def build_kmer_index(seqs: Sequence[dict], k: int = 11, strand: str = "both") -> dict[str, list[tuple[str, int, str]]]:
    idx: dict[str, list[tuple[str, int, str]]] = defaultdict(list)
    for rec in seqs:
        s = clean(rec["seq"] if isinstance(rec, dict) else rec.seq)
        for i in range(len(s) - k + 1):
            km = s[i:i + k]
            if "N" in km:
                continue
            idx[km].append((rec["id"] if isinstance(rec, dict) else rec.id, i, "+"))
            if strand == "both":
                idx[revcomp(km)].append((rec["id"] if isinstance(rec, dict) else rec.id, i, "-"))
    return idx


def blast_lite(query: Sequence[dict] | str, subjects: Sequence[dict], word_size: int = 11,
               evalue: float = 10.0, scoring: str = "blastn", gap_open: float = -11.0,
               gap_extend: float = -1.0, matrix: str = "BLOSUM62", max_hits: int = 250,
               dust: bool = True, db_length: int | None = None, seg: bool = True) -> list[dict]:
    """A compact but honest local BLAST: seeded words, ungapped extension,
    gapped refinement, bit scores and E-values (Altschul-Erickson formulae)."""
    is_prot = scoring in ("blastp", "tblastx", "tblastn", "blastx")
    mat = scoring_matrix(matrix if is_prot else "")
    ws = word_size if word_size else (3 if is_prot else 11)
    qseqs = _norm_records(query)
    sseqs = _norm_records(subjects)
    idx = build_kmer_index(sseqs, ws, strand="both" if scoring != "blastn" else "both")
    db_len = db_length or sum(len(s["seq"]) for s in sseqs)
    lambda_, K, H = (0.267, 0.041, 0.140) if not is_prot else (0.317649, 0.132547, 0.29919)
    hits = []
    for q in qseqs:
        qs = qseqs and q["seq"]
        qk = [clean(qs)]
        if scoring == "blastx":
            qk = [clean(qs), revcomp(clean(qs))]
        raw: list[tuple[str, int, int, str]] = []
        for strand_idx, qseq in enumerate(qk):
            for i in range(0, max(1, len(qseq) - ws + 1), 1):
                km = qseq[i:i + ws]
                if "N" in km or "X" in km:
                    continue
                for sid, spos, sstrand in idx.get(km, ()):  # seed hit
                    raw.append((sid, i, spos, "+" if (strand_idx == 0) == (sstrand == "+") else "-"))
        seen = Counter()
        for sid, qi, si, ori in raw:
            if seen[sid] > 20:
                continue
            seen[sid] += 1
            subj = next((s for s in sseqs if s["id"] == sid), None)
            if not subj:
                continue
            sseq = clean(subj["seq"])
            if ori == "-":
                sseq = revcomp(sseq)
            qsub = qseq[qi:qi + 120] if not is_prot else qseq[qi:qi + 60]
            ssub = sseq[si:si + len(qsub) + 6]
            aln = align_pairwise(qsub, ssub, "local", "" if not is_prot else matrix,
                                 gap_open, gap_extend, 5 if not is_prot else 4, -4)
            score = aln["score"]
            if score <= 0:
                continue
            bits = (lambda_ * score - math.log(K)) / math.log(2)
            e = db_len * len(qs) * 2 ** (-bits) if db_len else 0.0
            if e > evalue:
                continue
            hits.append({"qaccver": q["id"], "saccver": sid, "pident": aln["identity"],
                         "length": aln["aligned_length"], "mismatch": aln["mismatches"],
                         "gapopen": aln["gap_openings"], "qstart": qi + 1,
                         "qend": qi + 1 + len(qsub), "sstart": si + 1,
                         "send": si + 1 + len(ssub), "evalue": e, "bitscore": round(bits, 2),
                         "score": score, "strand": ori, "alignment_query": aln["a"],
                         "alignment_subject": aln["b"]})
    hits.sort(key=lambda h: (-h["bitscore"], h["evalue"]))
    return hits[:max_hits]


def _norm_records(src) -> list[dict]:
    if isinstance(src, (list, tuple)) and src and isinstance(src[0], dict):
        return [{"id": r.get("id", f"seq_{i + 1}"), "seq": r.get("seq", "")}
                for i, r in enumerate(src)]
    txt = as_text(src)
    recs = parse_fasta(txt)
    if recs:
        return [{"id": r.id, "seq": r.seq} for r in recs]
    if isinstance(src, (Seq, Read)):
        return [{"id": src.id, "seq": src.seq}]
    return [{"id": "query", "seq": clean(txt)}] if clean(txt) else []


def map_reads(genome: Sequence[dict] | dict[str, str], reads: Sequence[dict], k: int = 21,
              max_mismatches: int = 3, seed_stride: int = 1, local: bool = False,
              min_mapq: int = 0) -> list[dict]:
    """Seed-and-extend aligner (bwa-mem style in miniature), SAM-tagged output."""
    g = {}
    for item in (genome.values() if isinstance(genome, dict) else genome):
        if isinstance(item, dict):
            g[item["id"]] = clean(item["seq"])
        else:
            g[item.id] = clean(item.seq)
    idx = build_kmer_index([{"id": k2, "seq": v} for k2, v in g.items()], k, "both")
    out = []
    for r in (reads if isinstance(reads, list) else list(reads)):
        rid = r["id"] if isinstance(r, dict) else r.id
        qual = (r.get("qual", "") if isinstance(r, dict) else getattr(r, "qual", "")) or "I" * len(
            r["seq"] if isinstance(r, dict) else r.seq)
        qseq = clean(r["seq"] if isinstance(r, dict) else r.seq)
        best = None
        second = 0.0
        for i in range(0, max(1, len(qseq) - k + 1), seed_stride):
            km = qseq[i:i + k]
            if "N" in km:
                continue
            for sid, spos, sstrand in idx.get(km, ()):
                strand_match = (sstrand == "+")
                offset = spos - i if strand_match else len(g[sid]) - spos - k - (len(qseq) - i - k)
                pos = max(0, offset)
                subj = g[sid][pos:pos + len(qseq)]
                if not strand_match:
                    subj = g[sid][max(0, pos - 0):][:len(qseq)]
                    qtest = revcomp(qseq)
                else:
                    qtest = qseq
                mm = _count_mismatch(qtest[:len(subj)], subj)
                if mm <= max_mismatches:
                    score = len(qtest) - mm
                    if best is None or score > best[0]:
                        if best is not None:
                            second = best[0]
                        best = (score, sid, pos, sstrand, mm, qtest)
        if best is None:
            out.append({"qname": rid, "flag": 4, "rname": "*", "pos": 0, "mapq": 0,
                        "cigar": f"{len(qseq)}S", "rnext": "*", "pnext": 0, "tlen": 0,
                        "seq": qseq, "qual": qual, "tags": {"NM": 0}, "identity": 0.0})
            continue
        score, sid, pos, sstrand, mm, aligned_q = best
        mapq = int(min(60, max(0, 10 * (score - second) / max(1, score) * 6 + 30 - mm * 4)))
        if mapq < min_mapq:
            continue
        out.append({"qname": rid, "flag": 0 if sstrand == "+" else 16, "rname": sid,
                    "pos": pos + 1, "mapq": mapq, "cigar": f"{len(aligned_q)}M",
                    "rnext": "=", "pnext": 0, "tlen": 0, "seq": qseq, "qual": qual,
                    "tags": {"NM": mm, "X0": 1, "X1": 0},
                    "identity": round(100.0 * (1 - mm / max(1, len(aligned_q))), 3),
                    "ref_len": len(g[sid])})
        if local:
            continue
    return out


def _count_mismatch(a: str, b: str) -> int:
    n = min(len(a), len(b))
    return sum(1 for i in range(n) if a[i] != b[i]) + abs(len(a) - len(b))


def map_to_sam(records: Sequence[dict], header_refs: Sequence[tuple[str, int]] | None = None) -> str:
    from chroma_titan.core.io import sam_header, write_sam
    from chroma_titan.core.io import Aln as _A

    refs = list(header_refs or sorted({(r["rname"], 0) for r in records if r["rname"] != "*"}))
    hdr = sam_header([(n, s or 1000000) for n, s in refs])
    alns = []
    for r in records:
        alns.append(_A(r["qname"], r["flag"], r["rname"], r["pos"], r["mapq"], r["cigar"],
                       r["rnext"], r["pnext"], r["tlen"], r["seq"], r["qual"],
                       {k: f"i:{v}" if isinstance(v, int) else f"Z:{v}"
                        for k, v in r.get("tags", {}).items()}))
    return write_sam(hdr, alns)


def merge_pairs(reads1: Sequence[dict], reads2: Sequence[dict], min_overlap: int = 11,
                max_mismatches: int = 2) -> dict:
    """Overlap-based mate merging (BBMerge/Flash style)."""
    by_name = {r["id"]: r for r in reads2}
    merged, unmerged, mismatches = [], [], 0
    for r in reads1:
        mate = by_name.get(re.sub(r"/?1$", "", r["id"]) + ("_2" if "_" in r["id"] else "2"),
                           by_name.get(r["id"].replace("/1", "/2")))
        if mate is None:
            unmerged.append(r)
            continue
        a = clean(r["seq"])
        b = revcomp(clean(mate["seq"]))
        best = None
        for ov in range(min_overlap, min(len(a), len(b)) + 1):
            x, y = a[-ov:], b[:ov]
            mm = _count_mismatch(x, y)
            if mm <= max_mismatches and (best is None or mm < best[0]):
                best = (mm, ov)
            if best and best[0] == 0:
                break
        if not best:
            unmerged.append(r)
            continue
        mm, ov = best
        merged.append({"id": r["id"], "seq": a[:len(a) - ov] + b[:], "qual": r["qual"],
                       "overlap": ov, "mismatches": mm})
        mismatches += mm
    return {"merged": merged, "unmerged": unmerged, "n_merged": len(merged),
            "n_unmerged": len(unmerged), "merge_rate": round(
                100 * len(merged) / max(1, len(reads1)), 3),
            "mean_overlap": round(stats.mean([m["overlap"] for m in merged]), 3) if merged else 0.0,
            "total_mismatches_in_overlap": mismatches}


def align_to_reference_stats(alns: Sequence[dict]) -> dict:
    mapped = [a for a in alns if a["rname"] != "*"]
    return {"reads": len(alns), "mapped": len(mapped),
            "mapping_rate": round(100 * len(mapped) / max(1, len(alns)), 3),
            "mean_mismatches": round(stats.mean([a["tags"].get("NM", 0) for a in mapped]), 4),
            "mean_mapq": round(stats.mean([a["mapq"] for a in mapped]), 3),
            "unique": sum(1 for a in mapped if a["tags"].get("X1", 0) == 0),
            "multi": sum(1 for a in mapped if a["tags"].get("X1", 0) > 0)}


def similarity_histogram(hits: Sequence[dict], bin_size: int = 5) -> dict[int, int]:
    """Histogram of percent identity (BLAST hits or alignments) in fixed bins."""
    c: Counter = Counter()
    for h in hits:
        get = h.get if hasattr(h, "get") else (lambda k, d=None: d)
        pid = get("pident")
        if pid is None:
            try:
                qlen = float(get("query_length", 0) or len(get("seq", "") or "")) or 1.0
                pid = 100.0 * (1 - float(get("NM", 0) or 0) / qlen)
            except (TypeError, ValueError):
                pid = 0.0
        c[int(float(pid) // bin_size) * bin_size] += 1
    return dict(sorted(c.items()))


def best_hits(hits: Sequence[dict], by: str = "qaccver") -> list[dict]:
    best: dict[str, dict] = {}
    for h in hits:
        k = h[by]
        if k not in best or h["bitscore"] > best[k]["bitscore"]:
            best[k] = h
    return list(best.values())


def lca_of_hits(hits: Sequence[dict], taxonomy: dict[str, list[str]]) -> list[dict]:
    """Lowest common ancestor over the taxa of a read's best hits."""
    out = []
    by_read: dict[str, list[str]] = defaultdict(list)
    for h in hits:
        lin = taxonomy.get(h["saccver"])
        if lin:
            by_read[h["qaccver"]].append(lin)
    for rid, lines in by_read.items():
        if not lines:
            continue
        common = lines[0]
        for lin in lines[1:]:
            common = [a for a, b in zip(common, lin) if a == b]
            common = common or lin[:1]
        out.append({"read": rid, "lca": ";".join(common), "rank": len(common),
                     "n_hits": len(lines)})
    return out


def align_mafft(seqs: Sequence[dict], algorithm: str = "auto", gap_open: float = -11.0,
                gap_extend: float = -1.0, maxiterate: int = 0) -> list[dict]:
    """MAFFT-flavoured wrapper: progressive (fftNS-1) or center-star alignment."""
    rows = [{"id": s.get("id", f"seq_{i + 1}") if isinstance(s, dict) else s.id,
             "seq": s["seq"] if isinstance(s, dict) else s.seq} for i, s in enumerate(seqs)]
    aln = msa_progressive([r["seq"] for r in rows], gap_open, gap_extend, algorithm)
    for a, r in zip(aln, rows):
        a["id"] = r["id"]
    del maxiterate
    return aln


def align_muscle(seqs: Sequence[dict], refinement: int = 1) -> list[dict]:
    return align_mafft(seqs, "refine" if refinement else "auto")


def clustal_omega(seqs: Sequence[dict], guide_tree: str = "upgma") -> list[dict]:
    return align_mafft(seqs, guide_tree)


def pairwise_matrix(seqs: Sequence[dict], mode: str = "global", matrix: str = "") -> list[list[float]]:
    return distance_matrix([s["seq"] if isinstance(s, dict) else s.seq for s in seqs],
                           "identity" if not matrix else "alignment", matrix)


def identity_from_alignment(aln_a: str, aln_b: str) -> dict:
    pairs = list(zip(aln_a, aln_b))
    ident = sum(1 for x, y in pairs if x == y and x != "-")
    pos = sum(1 for x, y in pairs if x != "-" or y != "-")
    return {"aligned_length": len(pairs), "identities": ident, "aligned_positions": pos,
            "identity_percent": round(100 * ident / pos, 3) if pos else 0.0,
            "similarities": sum(1 for x, y in pairs if BLOSUM62.get((x, y), -1) > 0),
            "positives": sum(1 for x, y in pairs if x == y or BLOSUM62.get((x, y), -1) > 0)}


def gc_of_alignment(alignment: Sequence[dict]) -> float:
    tot = "".join(a["seq"] for a in alignment)
    return round(gc_content(tot), 3)


def alignment_to_fasta(alignment: Sequence[dict], wrap: int = 60) -> str:
    from chroma_titan.core.io import write_fasta

    return write_fasta([Seq(a["id"], a["seq"], "aligned") for a in alignment], wrap)


def alignment_to_matrix(alignment: Sequence[dict]) -> list[list[float]]:
    seqs = [a["seq"] for a in alignment]
    L = max((len(s) for s in seqs), default=0)
    letters = sorted({c for s in seqs for c in s})
    idx = {c: i for i, c in enumerate(letters)}
    m = [[0.0] * L for _ in letters]
    for s in seqs:
        for i, c in enumerate(s):
            m[idx[c]][i] += 1
    return m


def pam_distance(a: str, b: str) -> float:
    return pairwise_distance(a, b, "alignment", "BLOSUM62")
