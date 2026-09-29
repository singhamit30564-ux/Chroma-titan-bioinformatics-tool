"""Nucleotide/protein sequence primitives: composition, ORFs, codons, motifs.

All functions accept upper/lower case mixed strings and ignore whitespace and
FASTA headers, so they can be fed either pasted text or the output of
:func:`chroma_titan.core.io.as_text`.
"""

from __future__ import annotations

import math
import re
from collections import Counter, defaultdict

from chroma_titan.core.io import Seq, as_text, parse_fasta

BASES = "ACGT"
DNA = set("ACGTN")
RNA = set("ACGUN")
AA = "ACDEFGHIKLMNPQRSTVWY"
AA_ALL = "ACDEFGHIKLMNPQRSTVWYBXZ*-"

IUPAC = {
    "A": "A", "C": "C", "G": "G", "T": "T", "U": "U",
    "R": "[AG]", "Y": "[CT]", "S": "[GC]", "W": "[AT]", "K": "[GT]",
    "M": "[AC]", "B": "[CGT]", "D": "[AGT]", "H": "[ACT]", "V": "[ACG]",
    "N": "[ACGT]", "-": "[-.]?", "\.": "[ACGT]",
}

_COMP = {
    "A": "T", "T": "A", "G": "C", "C": "G", "U": "A", "N": "N",
    "R": "Y", "Y": "R", "S": "S", "W": "W", "K": "M", "M": "K",
    "B": "V", "V": "B", "D": "H", "H": "D", "-": "-", ".": ".", "*": "*",
}
_COMP.update({k.lower(): v.lower() for k, v in _COMP.items()})
COMPLEMENT = str.maketrans("".join(_COMP), "".join(_COMP.values()))

# ---------------------------------------------------------------------------
# genetic code
# ---------------------------------------------------------------------------
def _build_codon_table() -> dict[str, str]:
    table = {}
    for line in (
        "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
    ):
        pass
    bases = "TCAG"
    aas = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
    i = 0
    for b1 in bases:
        for b2 in bases:
            for b3 in bases:
                table[b1 + b2 + b3] = aas[i]
                i += 1
    return table


CODON_TABLE = _build_codon_table()
CODONS = list(CODON_TABLE)
START_CODONS = ("ATG", "GTG", "TTG")
STOP_CODONS = ("TAA", "TAG", "TGA")

#: E. coli K12 relative synonymous codon usage (freq per 1000), for back-translation
EOLI_USAGE = {
    "Ala": ["GCT", "GCC", "GCA", "GCG"], "Arg": ["CGT", "CGC", "CGA", "CGG", "AGA", "AGG"],
    "Asn": ["AAT", "AAC"], "Asp": ["GAT", "GAC"], "Cys": ["TGT", "TGC"],
    "Gln": ["CAA", "CAG"], "Glu": ["GAA", "GAG"], "Gly": ["GGT", "GGC", "GGA", "GGG"],
    "His": ["CAT", "CAC"], "Ile": ["ATT", "ATC", "ATA"], "Leu": ["TTA", "TTG", "CTT",
                                                                 "CTC", "CTA", "CTG"],
    "Lys": ["AAA", "AAG"], "Met": ["ATG"], "Phe": ["TTT", "TTC"], "Pro": ["CCT", "CCC",
                                                                          "CCA", "CCG"],
    "Ser": ["TCT", "TCC", "TCA", "TCG", "AGT", "AGC"], "Thr": ["ACT", "ACC", "ACA", "ACG"],
    "Trp": ["TGG"], "Tyr": ["TAT", "TAC"], "Val": ["GTT", "GTC", "GTA", "GTG"],
    "Ter": ["TAA", "TAG", "TGA"],
}
THREE_TO_ONE = {
    "Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C", "Gln": "Q", "Glu": "E",
    "Gly": "G", "His": "H", "Ile": "I", "Leu": "L", "Lys": "K", "Met": "M", "Phe": "F",
    "Pro": "P", "Ser": "S", "Thr": "T", "Trp": "W", "Tyr": "Y", "Val": "V", "Ter": "*",
}
ONE_TO_THREE = {v: k for k, v in THREE_TO_ONE.items()}


# ---------------------------------------------------------------------------
# basic cleaning / conversion
# ---------------------------------------------------------------------------
def clean(s: str, alphabet: str | None = None) -> str:
    """Strip whitespace/headers-ish characters and upper-case a sequence."""
    if s is None:
        return ""
    txt = s if "\n" not in s else "".join(
        ln for ln in as_text(s).splitlines() if not ln.startswith(">"))
    out = re.sub(r"[^A-Za-z*\-\.]", "", txt).upper()
    if alphabet is not None:
        allowed = set(alphabet)
        out = "".join(c for c in out if c in allowed)
    return out


def complement(s: str) -> str:
    return clean(s).translate(COMPLEMENT)


def reverse_complement(s: str) -> str:
    return complement(s)[::-1]


revcomp = reverse_complement


def transcribe(s: str) -> str:
    return clean(s).replace("T", "U")


def back_transcribe(s: str) -> str:
    return clean(s).replace("U", "T")


def is_dna(s: str) -> bool:
    c = clean(s)
    if not c:
        return False
    return sum(1 for x in c if x in "ACGTN") / len(c) > 0.9


def is_rna(s: str) -> bool:
    c = clean(s)
    return bool(c) and sum(1 for x in c if x in "ACGUN") / len(c) > 0.9


def is_protein(s: str) -> bool:
    c = clean(s)
    if not c:
        return False
    return sum(1 for x in c if x in AA) / len(c) > 0.85


def alphabet_of(s: str) -> str:
    if is_rna(s):
        return "ACGUN"
    if is_dna(s):
        return "ACGTN"
    return "".join(sorted(set(clean(s))))


def reverse_complement_fastq(seq: str, qual: str) -> tuple[str, str]:
    return revcomp(seq), qual[::-1]


# ---------------------------------------------------------------------------
# composition
# ---------------------------------------------------------------------------
def composition(s: str, alphabet: str | None = None) -> Counter:
    return Counter(clean(s, alphabet or "".join(sorted(set(clean(s)))) + "*-"))


def counts(s: str) -> Counter:
    return Counter(clean(s))


def gc_content(s: str, ignore_n: bool = False) -> float:
    c = clean(s)
    g = c.count("G") + c.count("C")
    n = len(c) - (c.count("N") if ignore_n else 0)
    return 100.0 * g / n if n else 0.0


def gc_skew(s: str) -> float:
    c = clean(s)
    g, t = c.count("G"), c.count("C")
    return (g - t) / (g + t) if g + t else 0.0


def at_skew(s: str) -> float:
    c = clean(s)
    a, t = c.count("A"), c.count("T")
    return (a - t) / (a + t) if a + t else 0.0


def chargaff(s: str) -> dict[str, float]:
    c = counts(s)
    tot = sum(c[b] for b in "ACGT") or 1
    return {b: c[b] / tot for b in "ACGT"}


def cpg_obs_exp(s: str) -> dict[str, float]:
    c = clean(s)
    n = len(c)
    if n < 2:
        return {"cg": 0.0, "gc": 0.0, "ratio": 0.0, "c_count": 0, "g_count": 0}
    cg = c.count("CG")
    gc = 100.0 * (c.count("G") + c.count("C")) / n
    denom = (c.count("C") * c.count("G")) / n
    return {"cg": cg, "gc": gc, "ratio": (cg * n) / denom if denom else 0.0,
            "c_count": c.count("C"), "g_count": c.count("G")}


def nucleotide_frequencies(s: str) -> dict[str, float]:
    c = counts(s)
    tot = sum(c.values()) or 1
    return {k: v / tot for k, v in sorted(c.items())}


def entropy_of(s: str, k: int = 1, alphabet: str = "") -> float:
    """Shannon entropy (bits per k-mer) of the k-mer profile."""
    if not s:
        return 0.0
    kmers = kmer_list(s, k, alphabet)
    tot = len(kmers) or 1
    h = -sum((v / tot) * math.log2(v / tot) for v in Counter(kmers).values())
    return h


def sliding_windows(s: str, size: int = 100, step: int | None = None,
                    fn=gc_content) -> list[dict]:
    c = clean(s)
    step = step or size
    out = []
    for i in range(0, max(1, len(c) - size + 1), step):
        w = c[i:i + size]
        out.append({"start": i, "end": i + len(w), "value": fn(w)})
    return out


def windowed_stats(s: str, size: int = 100, step: int | None = None) -> list[dict]:
    c = clean(s)
    step = step or size
    out = []
    for i in range(0, max(1, len(c) - size + 1), step):
        w = c[i:i + size]
        cnt = counts(w)
        gc = 100.0 * (cnt["G"] + cnt["C"]) / len(w) if w else 0.0
        out.append({"start": i, "end": i + len(w), "length": len(w), "A": cnt["A"],
                    "C": cnt["C"], "G": cnt["G"], "T": cnt["T"], "N": cnt["N"],
                    "GC_percent": gc, "entropy": entropy_of(w, 1, "ACGT"),
                    "GC_skew": gc_skew(w)})
    return out


# ---------------------------------------------------------------------------
# k-mers
# ---------------------------------------------------------------------------
def kmer_list(s: str, k: int = 3, alphabet: str = "ACGTN") -> list[str]:
    c = clean(s, alphabet) if alphabet else clean(s)
    return [c[i:i + k] for i in range(len(c) - k + 1)]


def kmer_counts(s: str, k: int = 3, alphabet: str = "ACGTN",
                canonical: bool = False, min_count: int = 1) -> Counter:
    km = kmer_list(s, k, alphabet)
    if canonical:
        km = [min(x, revcomp(x)) for x in km if set(x) <= set("ACGT")]
    cnt = Counter(km)
    return Counter({x: v for x, v in cnt.items() if v >= min_count})


def base_counts(s: str) -> Counter:
    return counts(s)


def oligo_frequency(s: str, k: int = 6) -> Counter:
    return Counter(kmer_list(s, k, ""))


def profile_vector(s: str, k: int = 4) -> list[float]:
    c = kmer_counts(s, k, "ACGT")
    tot = sum(c.values()) or 1
    return [c.get(x, 0) / tot for x in sorted(kmer_list("".join(BASES * k), k, "ACGT"))]


def n50(lengths: list[int]) -> int:
    tot = sum(lengths)
    if not tot:
        return 0
    acc = 0
    for L in sorted(lengths, reverse=True):
        acc += L
        if acc >= tot / 2:
            return L
    return 0


def l50(lengths: list[int]) -> int:
    tot = sum(lengths)
    acc, n = 0, 0
    for L in sorted(lengths, reverse=True):
        acc += L
        n += 1
        if acc >= tot / 2:
            return n
    return n


def minhash(s: str, k: int = 21, num: int = 100, seed: int = 42) -> list[int]:
    """Bottom-k sketch (Mash-like) of canonical k-mers."""
    import hashlib

    hs = []
    for km in set(kmer_list(s, k, "ACGT")):
        canon = min(km, revcomp(km))
        d = hashlib.md5(f"{seed}{canon}".encode()).digest()
        hs.append(int.from_bytes(d[:8], "big"))
    return sorted(hs)[:num]


def minhash_jaccard(a: list[int], b: list[int]) -> float:
    sa, sb = set(a), set(b)
    u = sa | sb
    return len(sa & sb) / len(u) if u else 0.0


def mash_distance(a: list[int], b: list[int], k: int = 21) -> float:
    j = minhash_jaccard(a, b)
    if j == 0:
        return float("nan")
    if j == 1:
        return 0.0
    return max(0.0, -math.log(2 * j / (1 + j)) / k)


# ---------------------------------------------------------------------------
# codons / translation
# ---------------------------------------------------------------------------
def codons(s: str, frame: int = 0) -> list[str]:
    c = back_transcribe(clean(s))
    return [c[i:i + 3] for i in range(frame, len(c) - 2, 3)]


def translate(s: str, frame: int = 0, to_stop: bool = False, start_m: str = "M") -> str:
    c = back_transcribe(clean(s))
    out = []
    for cod in codons(c, frame):
        if len(cod) < 3:
            break
        aa = CODON_TABLE.get(cod.replace("N", "N"), None)
        if aa is None or "N" in cod:
            aa = "X"
        if aa == "*":
            if to_stop:
                break
            aa = "*"
        out.append(aa)
    p = "".join(out)
    if p.startswith("*"):
        p = start_m + p[1:]
    return p


def translate_all_frames(s: str, to_stop: bool = True) -> list[dict]:
    c = back_transcribe(clean(s))
    res = []
    for f in range(3):
        res.append({"frame": f + 1, "protein": translate(c, f, to_stop)})
        res.append({"frame": -(f + 1), "protein": translate(revcomp(c), f, to_stop)})
    return res


def find_orfs(s: str, min_len: int = 300, table: str = "Standard",
              required_start: str = "ATG", partial: bool = False) -> list[dict]:
    """Six-frame ORF finder. ``min_len`` is in nucleotides."""
    c = back_transcribe(clean(s))
    orfs: list[dict] = []
    starts = (required_start,) if required_start else ()
    for strand, seq in (("+", c), ("-", revcomp(c))):
        for f in range(3):
            i = f
            while i < len(seq) - 2:
                cod = seq[i:i + 3]
                if cod == "ATG" or (not required_start and cod in START_CODONS):
                    j = i
                    while j < len(seq) - 2:
                        cc = seq[j:j + 3]
                        if cc in STOP_CODONS:
                            orfs.append({
                                "id": f"orf_{len(orfs) + 1}",
                                "strand": strand, "frame": f + 1 if strand == "+" else -(f + 1),
                                "start": i if strand == "+" else len(seq) - j - 3,
                                "end": j + 3 if strand == "+" else len(seq) - i,
                                "length": j + 3 - i, "aa_length": (j + 3 - i) // 3 - 1,
                                "protein": translate(seq[i:j], 0, to_stop=True),
                                "start_codon": cod, "stop_codon": cc,
                            })
                            i = j + 3
                            break
                        j += 3
                    else:
                        if partial and j >= len(seq) - 2:
                            orfs.append({"id": f"orf_partial_{len(orfs) + 1}",
                                         "strand": strand, "frame": f + 1, "start": i,
                                         "end": len(seq) - len(seq) % 3,
                                         "length": len(seq) - i, "aa_length": (len(seq) - i) // 3,
                                         "protein": translate(seq[i:], 0, to_stop=True),
                                         "start_codon": cod, "stop_codon": ""})
                        i = j + 3
                        continue
                else:
                    i += 3
    orfs = [o for o in orfs if o["length"] >= min_len]
    orfs.sort(key=lambda o: -o["length"])
    for n, o in enumerate(orfs, 1):
        o["id"] = f"orf_{n}"
    return orfs


def codon_usage(s: str, frame: int = 0) -> Counter:
    c = back_transcribe(clean(s))
    cnt = Counter()
    for cod in codons(c, frame):
        if len(cod) == 3 and cod not in STOP_CODONS and "N" not in cod:
            cnt[cod] += 1
    return cnt


def rscu_table(s: str) -> list[dict]:
    """Relative synonymous codon usage per amino acid."""
    by_aa: dict[str, Counter] = defaultdict(Counter)
    for cod, n in codon_usage(s).items():
        by_aa[CODON_TABLE.get(cod, "X")][cod] += n
    rows = []
    for aa, cnt in sorted(by_aa.items()):
        tot = sum(cnt.values())
        for cod, n in cnt.most_common():
            rows.append({"amino_acid": ONE_TO_THREE.get(aa, aa), "codon": cod, "count": n,
                         "RSCU": n / (tot / len([c for c in CODON_TABLE
                                                 if CODON_TABLE[c] == aa])) if tot else 0.0})
    return rows


def cai(seq: str, reference_codon_freq: dict[str, float] | None = None) -> float:
    """Codon adaptation index vs a reference frequency table (defaults to E. coli)."""
    w = {}
    ref = reference_codon_freq or _ecoli_w_table()
    c = back_transcribe(clean(seq))
    weights = []
    for cod in codons(c):
        aa = CODON_TABLE.get(cod)
        if not aa or aa == "*":
            continue
        syn = [cc for cc, a in CODON_TABLE.items() if a == aa]
        mx = max(ref.get(cc, 1e-9) for cc in syn) or 1e-9
        w[aa] = ref.get(cod, 1e-9) / mx
        weights.append(max(1e-9, w[aa]))
    if not weights:
        return 0.0
    return math.exp(sum(math.log(x) for x in weights) / len(weights))


def _ecoli_w_table() -> dict[str, float]:
    freq: dict[str, float] = {}
    weights = {"GCT": 0.16, "GCC": 0.23, "GCA": 0.18, "GCG": 0.43, "CGT": 0.38, "CGC": 0.40,
               "CGA": 0.06, "CGG": 0.10, "AGA": 0.04, "AGG": 0.02, "AAT": 0.45, "AAC": 0.55,
               "GAT": 0.63, "GAC": 0.37, "TGT": 0.45, "TGC": 0.55, "CAA": 0.35, "CAG": 0.65,
               "GAA": 0.68, "GAG": 0.32, "GGT": 0.34, "GGC": 0.39, "GGA": 0.11, "GGG": 0.16,
               "CAT": 0.55, "CAC": 0.45, "ATT": 0.51, "ATC": 0.42, "ATA": 0.07, "TTA": 0.13,
               "TTG": 0.13, "CTT": 0.11, "CTC": 0.10, "CTA": 0.04, "CTG": 0.49, "AAA": 0.74,
               "AAG": 0.26, "ATG": 1.0, "TTT": 0.57, "TTC": 0.43, "CCT": 0.16, "CCC": 0.12,
               "CCA": 0.20, "CCG": 0.52, "TCT": 0.15, "TCC": 0.15, "TCA": 0.12, "TCG": 0.15,
               "AGT": 0.15, "AGC": 0.28, "ACT": 0.13, "ACC": 0.30, "ACA": 0.19, "ACG": 0.38,
               "TGG": 1.0, "TAT": 0.55, "TAC": 0.45, "GTT": 0.26, "GTC": 0.22, "GTA": 0.15,
               "GTG": 0.37}
    for cod, aa in CODON_TABLE.items():
        freq[cod] = weights.get(cod, 0.0)
    del aa
    return freq


def effective_number_of_codons(seq: str) -> float:
    """ENc (Wright 1990) computed from within-genome codon bias."""
    cnt = codon_usage(seq)
    if not cnt:
        return float("nan")
    by_aa: dict[str, list[int]] = defaultdict(list)
    for cod, n in cnt.items():
        by_aa[CODON_TABLE.get(cod, "X")].append(n)
    F3 = []
    for aa, ns in by_aa.items():
        n = sum(ns)
        if len(ns) < 2 or n < 2:
            continue
        f = (sum((x / n) ** 2 for x in ns) * n - 1) / (n - 1) if n else 0
        F3.append(max(0.0, min(1.0, f)))
    if not F3:
        return 61.0
    f = sum(F3) / len(F3)
    return max(20.0, min(61.0, 20.0 + 29.0 / (f * (3 - 2 * f) if f * (3 - 2 * f) else 1)))


def gc3(seq: str) -> float:
    c = back_transcribe(clean(seq))
    third = [cod[2] for cod in codons(c) if len(cod) == 3 and cod not in STOP_CODONS]
    if not third:
        return 0.0
    return 100.0 * sum(1 for x in third if x in "GC") / len(third)


def back_translate(protein: str, usage: str = "E. coli", optimise: bool = True) -> str:
    """Reverse translate, picking the most used synonymous codon when optimising."""
    p = clean(protein, AA_ALL)
    table = EOLI_USAGE
    out = []
    for aa1 in p:
        name = THREE_TO_ONE.get(aa1, aa1)
        cods = table.get(name if name in table else ONE_TO_THREE.get(aa1, "Met"), ["ATG"])
        if not cods:
            cods = [c for c, a in CODON_TABLE.items() if a == aa1] or ["NNN"]
        out.append(cods[0] if optimise else cods[-1])
    return "".join(out)


def optimise_codons(seq: str, gc_target: float = 50.0) -> str:
    """Codon optimisation: most-used synonymous codon, nudged toward a GC target."""
    c = back_transcribe(clean(seq))
    out = []
    for i, cod in enumerate(codons(c)):
        aa = CODON_TABLE.get(cod)
        if not aa or "N" in cod:
            out.append(cod if len(cod) == 3 else "NNN")
            continue
        name = THREE_TO_ONE.get(aa, aa)
        cands = EOLI_USAGE.get(name, [cod])
        best = cands[0]
        if gc_target:
            best = min(cands, key=lambda x: abs(100 * (x.count("G") + x.count("C")) / 3 - gc_target))
        out.append(best)
    return "".join(out)


# ---------------------------------------------------------------------------
# motifs / patterns
# ---------------------------------------------------------------------------
def iupac_to_regex(pattern: str) -> str:
    pat = clean(pattern).replace("U", "T")
    out = []
    for c in pat:
        if c == ".":
            out.append(".")
        elif c in "[]{}()":
            out.append(c)
        else:
            out.append(IUPAC.get(c, re.escape(c)))
    return "".join(out)


def find_motifs(s: str, pattern: str, strand: str = "both",
                allow_mismatch: int = 0) -> list[dict]:
    """Locate IUPAC/regex motifs, both strands, with optional Hamming mismatch."""
    seq = back_transcribe(clean(s))
    rc = revcomp(seq)
    hits: list[dict] = []
    plain = None
    if pattern and all(ch in IUPAC and IUPAC[ch] not in ("[ACGT]",) for ch in clean(pattern)) \
            and not any(ch in pattern for ch in "[](){}"):
        plain = "".join(IUPAC.get(ch, ch) for ch in clean(pattern))
    for name, subject in (("+", seq), ("-", rc)):
        if strand == "+" and name == "-":
            continue
        if strand == "-" and name == "+":
            continue
        if plain and allow_mismatch:
            L = len(re.sub(r"[\[\]]", "", plain).replace("[ACGT]", "N"))
            core = "".join(ch for ch in clean(pattern) if ch in "ACGTRYWSKMBDHVN")
            for i in range(len(subject) - len(core) + 1):
                sub = subject[i:i + len(core)]
                mm = 0
                ok = True
                for a, b in zip(core, sub):
                    alts = [b] if b not in IUPAC else [x for x in _expand(b)]
                    if a not in alts:
                        mm += 1
                    if mm > allow_mismatch:
                        ok = False
                        break
                if ok:
                    pos = i if name == "+" else len(seq) - i - len(core)
                    hits.append({"pattern": pattern, "strand": name, "start": pos,
                                 "end": pos + len(core), "match": sub, "mismatches": mm})
        else:
            rx = re.compile(iupac_to_regex(pattern)) if not any(c in pattern for c in "[") \
                else re.compile(pattern)
            for m in rx.finditer(subject):
                pos = m.start() if name == "+" else len(seq) - m.end()
                hits.append({"pattern": pattern, "strand": name, "start": pos,
                             "end": pos + (m.end() - m.start()), "match": m.group(),
                             "mismatches": 0})
    hits.sort(key=lambda h: h["start"])
    return hits


def _expand(letter: str) -> str:
    return {
        "A": "A", "C": "C", "G": "G", "T": "T", "U": "U", "R": "AG", "Y": "CT", "S": "GC",
        "W": "AT", "K": "GT", "M": "AC", "B": "CGT", "D": "AGT", "H": "ACT", "V": "ACG",
        "N": "ACGT",
    }.get(letter, letter)


def expand_iupac(pattern: str) -> list[str]:
    combos = [""]
    for ch in clean(pattern):
        alts = _expand(ch)
        combos = [c + a for c in combos for a in alts]
        if len(combos) > 5000:
            break
    return combos[:5000]


def degenerate_to_regex(pattern: str) -> str:
    return iupac_to_regex(pattern)


def consensus(seqs: list[str], mode: str = "majority") -> str:
    if not seqs:
        return ""
    L = max(len(s) for s in seqs)
    out = []
    for i in range(L):
        col = Counter(s[i] for s in seqs if i < len(s) and s[i] not in "-.")
        if not col:
            out.append("-")
        elif mode == "strict":
            top, n = col.most_common(1)[0]
            out.append(top if n == len(col) else "N")
        elif mode == "ambiguity":
            top = [c for c, n in col.items() if n >= 0.5 * sum(col.values())]
            out.append(top[0] if len(top) == 1 else _collapse("".join(sorted(top))))
        else:
            out.append(col.most_common(1)[0][0])
    return "".join(out)


def _collapse(letters: str) -> str:
    setm = set(letters)
    for k, v in {"ACGT": "N", "ACG": "V", "ACT": "H", "AGT": "D", "CGT": "B", "AG": "R",
                 "CT": "Y", "GC": "S", "AT": "W", "GT": "K", "AC": "M"}.items():
        if setm == set(k):
            return v
    return next(iter(sorted(setm)))


def consensus_iupac(seqs: list[str]) -> str:
    return consensus(seqs, mode="ambiguity")


def find_repeats(s: str, min_len: int = 6) -> list[dict]:
    """Exact tandem-repeat detection via shift comparison (TRF-lite)."""
    c = clean(s)
    out = []
    for unit in range(1, min(30, max(2, len(c) // 2)) + 1):
        i = 0
        while i < len(c) - 2 * unit:
            j = i + unit
            n = 0
            while j + n < len(c) and c[i + n] == c[j + n]:
                n += 1
            reps = n // unit + 1
            if n + unit >= min_len and reps >= 2:
                out.append({"start": i, "end": j + n, "period": unit, "copies": reps,
                            "motif": c[i:i + unit]})
                i = j + n
            else:
                i += 1
    out.sort(key=lambda r: -(r["end"] - r["start"]))
    merged, seen = [], set()
    for r in out:
        key = (r["start"], r["period"])
        if key in seen:
            continue
        seen.add(key)
        merged.append(r)
    return merged[:2000]


def low_complexity(s: str, window: int = 64, level: float = 2.4) -> list[dict]:
    """DUST-style low complexity regions."""
    c = clean(s)
    out, inmask, start = [], False, 0
    step = max(1, window // 4)
    for i in range(0, max(1, len(c) - window + 1), step):
        w = c[i:i + window]
        e = entropy_of(w, 2, "ACGT")
        low = e < level
        if low and not inmask:
            inmask, start = True, i
        elif not low and inmask:
            inmask = False
            out.append({"start": start, "end": i + window, "entropy": round(e, 3)})
    if inmask:
        out.append({"start": start, "end": len(c), "entropy": 0.0})
    return out


def dust_score(s: str, window: int = 64) -> float:
    c = clean(s)
    if len(c) < 3:
        return 0.0
    cnt = Counter(kmer_list(c, 3, "ACGT"))
    n = len(c) - 2
    score = sum(v * (v - 1) / 2 for v in cnt.values())
    return min(1.0, score / (n * (n - 1) / 2)) if n > 1 else 0.0


def mask_low_complexity(s: str, window: int = 64, level: float = 2.4, char: str = "N") -> str:
    c = clean(s)
    keep = list(c)
    for r in low_complexity(c, window, level):
        for i in range(r["start"], min(len(keep), r["end"])):
            keep[i] = char
    return "".join(keep)


def homopolymer_runs(s: str, min_len: int = 3) -> list[dict]:
    c = clean(s)
    out = []
    i = 0
    while i < len(c):
        j = i
        while j < len(c) and c[j] == c[i]:
            j += 1
        if j - i >= min_len:
            out.append({"base": c[i], "start": i, "end": j, "length": j - i})
        i = j
    return out


def longest_homopolymer(s: str) -> dict:
    runs = homopolymer_runs(s, 1)
    if not runs:
        return {"base": "", "start": 0, "end": 0, "length": 0}
    return max(runs, key=lambda r: r["length"])


def n50_of_seqs(seqs: list[str]) -> int:
    return n50([len(clean(s)) for s in seqs])


def max_run_of_n(s: str) -> int:
    best = 0
    for r in homopolymer_runs(s, 1):
        if r["base"] == "N":
            best = max(best, r["length"])
    return best


def chaos_game_points(s: str) -> list[tuple[float, float]]:
    pts, x, y = [], 0.5, 0.5
    corner = {"C": (0.0, 0.0), "A": (1.0, 0.0), "G": (0.0, 1.0), "T": (1.0, 1.0)}
    for ch in clean(s):
        cx, cy = corner.get(ch, (0.5, 0.5))
        x, y = (x + cx) / 2, (y + cy) / 2
        pts.append((x, y))
    return pts


def tetranucleotide_profile(s: str) -> Counter:
    return kmer_counts(s, 4, "ACGT")


def gc_clamp(seq: str) -> int:
    return len(seq) - len(seq.rstrip("GC"))


def sha1_of(seq: str) -> str:
    import hashlib

    return hashlib.sha1(clean(seq).encode()).hexdigest()[:16]


def sequence_md5_id(seq: str) -> str:
    import hashlib

    return hashlib.md5(clean(seq).encode()).hexdigest()


# ---------------------------------------------------------------------------
# thermodynamics / oligo quality
# ---------------------------------------------------------------------------
def melting_temp(seq: str, method: str = "nearest", na_conc: float = 50.0,
                 primer_conc: float = 250e-9) -> float:
    """Tm in C. method: wallace (2/4 rule), gc, breslauer, nearest (SantaLucia 1998)."""
    s = clean(seq, "ACGTN")
    if not s:
        return 0.0
    if method in ("wallace", "2-4"):
        return 2 * (s.count("A") + s.count("T")) + 4 * (s.count("G") + s.count("C"))
    if method == "gc":
        return 64.9 + 41 * (s.count("G") + s.count("C") - 16.4) / max(1, len(s))
    nn = {
        "AA": (-7.9, -22.2), "AT": (-7.2, -20.4), "TA": (-7.2, -21.3), "CA": (-8.5, -22.7),
        "GT": (-8.4, -22.4), "CT": (-7.8, -19.9), "GA": (-8.2, -22.2), "CG": (-10.6, -27.2),
        "GC": (-9.8, -24.4), "GG": (-8.0, -19.9), "AC": (-7.8, -21.0), "TC": (-8.2, -22.2),
        "TG": (-8.5, -22.7), "CC": (-8.4, -22.4), "TT": (-7.9, -22.2), "AG": (-7.8, -21.0),
    }
    dh = ds = 0.0
    for i in range(len(s) - 1):
        pair = s[i:i + 2]
        if pair in nn:
            dh += nn[pair][0]
            ds += nn[pair][1]
        else:
            dh += -8.0
            ds += -20.0
    if s[0] in "AC":
        dh += 0.1
        ds += -2.8
    if s[-1] in "AC":
        dh += 0.1
        ds += -2.8
    if s[0] == "G" or s[-1] == "G":
        pass
    init_dh, init_ds = (0.0, 0.0) if s.startswith("GC") or s.endswith("GC") else (0.0, 0.0)
    dh += init_dh
    ds += init_ds
    r = 1.987
    na = na_conc / 1000.0
    ds_corr = ds + 0.368 * len(s) * math.log(na) if na > 0 else ds
    kelvin = (1000 * dh) / (ds_corr + r * math.log(primer_conc / 4)) if ds_corr else 0.0
    return kelvin - 273.15 + 16.6 * math.log10(na_conc / 50.0) if na_conc else kelvin - 273.15


def gc_clamp_score(seq: str) -> int:
    return gc_clamp(seq)


def wallace_tm(seq: str) -> float:
    return melting_temp(seq, "wallace")


def salt_corrected_tm(seq: str, na_mm: float = 50.0) -> float:
    return melting_temp(seq, "nearest", na_conc=na_mm)


def delta_g_duplex(seq: str) -> float:
    """Crude nearest-neighbour dG at 37C for a self-complementary duplex."""
    s = clean(seq, "ACGT")
    dh = sum(100 * (0.1 if a in "AC" else 0.05) for a in s) - 300 * (len(s) - 1) / 10
    return round(dh / 1000 * 4.184, 2)


def hairpin_score(seq: str, loop_min: int = 3, max_stem: int = 12) -> dict:
    """Stem-loop finder returning a free-energy-ish score (higher = more stable)."""
    s = back_transcribe(clean(seq, "ACGUT"))
    best = {"stem": 0, "loop": 0, "score": 0.0, "start": 0, "delta_g": 0.0}
    for i in range(len(s) - 2 * 4):
        for stem in range(min(max_stem, (len(s) - i) // 3), 1, -1):
            left = s[i:i + stem]
            j = i + stem + loop_min
            if j + stem > len(s):
                continue
            right = s[j:j + stem]
            pairs = sum(1 for a, b in zip(left, right[::-1])
                        if (a, b) in {("A", "T"), ("T", "A"), ("G", "C"), ("C", "G")})
            if pairs >= 3 and pairs >= 0.6 * stem:
                dg = -1.5 * pairs + 0.4 * stem
                if dg < best["delta_g"]:
                    best = {"stem": stem, "loop": j - (i + stem), "score": round(-dg, 2),
                            "start": i, "delta_g": round(dg, 2)}
                break
    return best


def hairpin_delta_g(seq: str) -> float:
    return hairpin_score(seq)["delta_g"]


def self_dimer_score(seq1: str, seq2: str | None = None) -> dict:
    """3' end complementarity score using a small ungapped alignment."""
    a = back_transcribe(clean(seq1, "ACGUT"))[-8:]
    b = back_transcribe(clean(seq2 or seq1, "ACGUT"))[-8:]
    rc = revcomp(b)
    m = sum(1 for x, y in zip(a, rc) if x == y)
    return {"dimer_score": m, "delta_g": round(-1.4 * m, 2),
            "three_prime_overlap": m, "seq1_tail": a, "seq2_tail": rc}


def dna_shape_proxy(seq: str) -> list[dict]:
    """Very small DNA-shape surrogate: propeller/twist/minor-width heuristics."""
    rows = []
    for i, dinuc in enumerate(kmer_list(seq, 2, "ACGT")):
        rows.append({
            "position": i, "dinucleotide": dinuc,
            "roll": round(2.0 + 3.0 * (dinuc.count("G") + dinuc.count("C")) / 2, 2),
            "minor_width": round(3.0 + 1.5 * (dinuc.count("A") + dinuc.count("T")), 2),
            "propeller": round(6.0 + 2.0 * (dinuc in ("GA", "AG", "TC", "CT")), 2),
        })
    return rows


def nucleosome_position_proxy(seq: str) -> list[float]:
    return [round(0.5 + 0.5 * math.sin(i * 0.6) * (gc_content(seq[i:i + 10]) / 100), 3)
            for i in range(0, max(1, len(seq) - 10), 5)][:400]


# ---------------------------------------------------------------------------
# restriction enzymes / in-silico cloning
# ---------------------------------------------------------------------------
RESTRICTION_ENZYMES: dict[str, tuple[str, int, int]] = {
    # name: (recognition (IUPAC), cut offset after 5' of top strand, cut on bottom)
    "EcoRI": ("GAATTC", 1, 5), "BamHI": ("GGATCC", 1, 5), "HindIII": ("AAGCTT", 1, 5),
    "NotI": ("GCGGCCGC", 2, 6), "EcoRV": ("GATATC", 3, 3), "PstI": ("CTGCAG", 5, 1),
    "SalI": ("GTCGAC", 1, 5), "XhoI": ("CTCGAG", 1, 5), "SpeI": ("ACTAGT", 1, 5),
    "NcoI": ("CCATGG", 1, 5), "NheI": ("GCTAGC", 1, 5), "KpnI": ("GGTACC", 5, 1),
    "SmaI": ("CCCGGG", 3, 3), "XbaI": ("TCTAGA", 1, 5), "SacI": ("GAGCTC", 5, 1),
    "ClaI": ("ATCGAT", 2, 4), "EcoNI": ("CCTNNAGG", 4, 5), "BglII": ("AGATCT", 1, 5),
    "AvrII": ("CCTAGG", 1, 5), "PvuII": ("CGATCG", 3, 3), "AatII": ("GACGTC", 1, 5),
    "ApaI": ("GGGCCC", 5, 1), "AscI": ("GGCGCGCC", 1, 7), "BsaI": ("GGTCTC", 7, 11),
    "BbsI": ("GAAGAC", 8, 12), "Esp3I": ("CCTCTTC", 5, 6), "PacI": ("TTAATTAA", 5, 3),
    "SwaI": ("ATTTAAAT", 4, 4), "MluI": ("ACGCGT", 1, 5), "SbfI": ("CCTGCAGG", 2, 6),
    "ZraI": ("GACGTC", 3, 3), "HinfI": ("GANTC", 1, 4), "HaeIII": ("GGCC", 2, 2),
    "AluI": ("AGCT", 2, 2), "TaqI": ("TCGA", 1, 3), "MseI": ("TTAA", 1, 3),
    "DpnI": ("GATC", 3, 1), "DpmI": ("GATGAC", 5, 8), "FokI": ("GGATG", 9, 13),
    "BsmBI": ("CGTCTC", 6, 7), "LguI": ("ACGGT", 5, 5), "AsiSI": ("GCGATCGC", 3, 5),
}


def enzyme_table() -> list[dict]:
    rows = []
    for name, (site, c1, c2) in RESTRICTION_ENZYMES.items():
        rows.append({"name": name, "recognition": site, "cut_top": c1, "cut_bottom": c2,
                     "is_6cut": len(site.replace("N", "")) >= 6,
                     "frequency_kb": round(4 ** len(site.replace("N", ".")) / 1000.0, 3)})
    return rows


def cut_positions(s: str, enzymes: list[str] | None = None, min_cut: int = 0) -> list[dict]:
    """All cut sites of the given (or all) enzymes."""
    seq = back_transcribe(clean(s))
    out = []
    for name, (site, c1, c2) in (RESTRICTION_ENZYMES.items()):
        if enzymes and name not in enzymes:
            continue
        rx = re.compile(iupac_to_regex(site))
        for m in rx.finditer(seq):
            pos = m.start() + c1
            out.append({"enzyme": name, "site": m.group(), "position": pos,
                        "start": m.start(), "end": m.end(), "strand": "+",
                        "overhang": m.group()[c1:c2] if c2 > c1 else ""})
    out.sort(key=lambda r: (r["enzyme"], r["position"]))
    if min_cut:
        keep = {}
        for r in out:
            keep[r["enzyme"]] = keep.get(r["enzyme"], 0) + 1
        out = [r for r in out if keep[r["enzyme"]] >= min_cut]
    return out


def cut_frequencies(s: str, enzymes: list[str] | None = None) -> list[dict]:
    sites = cut_positions(s, enzymes)
    cnt = Counter(x["enzyme"] for x in sites)
    L = len(back_transcribe(clean(s))) or 1
    rows = []
    for name in (enzymes or RESTRICTION_ENZYMES):
        n = cnt.get(name, 0)
        rows.append({"enzyme": name, "cuts": n, "per_kb": round(1000 * n / L, 3),
                     "blunt": RESTRICTION_ENZYMES[name][1] == RESTRICTION_ENZYMES[name][2]
                     if name in RESTRICTION_ENZYMES else False})
    return rows


def fragment_sizes(s: str, enzyme: str = "EcoRI", circular: bool = False) -> list[int]:
    seq = back_transcribe(clean(s))
    site, c1, _ = RESTRICTION_ENZYMES.get(enzyme, (None, 0, 0))
    if not site:
        return []
    cuts = [m.start() + c1 for m in re.finditer(iupac_to_regex(site), seq)]
    if not cuts:
        return [len(seq)]
    if circular:
        cuts = sorted(set(cuts + [0]))
        frags = [cuts[i + 1] - cuts[i] for i in range(len(cuts) - 1)]
        frags.append(len(seq) - cuts[-1] + cuts[0])
        return [f for f in frags if f > 0]
    bounds = [0] + cuts + [len(seq)]
    return [bounds[i + 1] - bounds[i] for i in range(len(bounds) - 1) if bounds[i + 1] > bounds[i]]


def digest(s: str, enzyme: str = "EcoRI") -> list[Seq]:
    seq = back_transcribe(clean(s))
    site, c1, c2 = RESTRICTION_ENZYMES.get(enzyme, (None, 0, 0))
    if not site:
        return [Seq("fragment_1", seq, f"unknown enzyme {enzyme}")]
    cuts = [m.start() + c1 for m in re.finditer(iupac_to_regex(site), seq)]
    bounds = [0] + cuts + [len(seq)]
    out = []
    for i in range(len(bounds) - 1):
        a, b = bounds[i], bounds[i + 1]
        if b > a:
            out.append(Seq(f"{enzyme}_{i + 1}_fragment", seq[a:b], f"length={b - a}"))
    return out


def ligation_score(left: str, right: str) -> float:
    a = clean(left)[-6:]
    b = clean(right)[:6]
    return sum(1 for x, y in zip(a, revcomp(b)) if x == y)


def unique_cutters(s: str, enzymes: list[str] | None = None) -> list[str]:
    rows = cut_frequencies(s, enzymes)
    return [r["enzyme"] for r in rows if r["cuts"] == 1]


# ---------------------------------------------------------------------------
# FASTA-level helpers used by many tools
# ---------------------------------------------------------------------------
def read_seq_records(src) -> list[Seq]:
    txt = as_text(src)
    recs = parse_fasta(txt)
    if recs:
        return recs
    s = clean(txt)
    return [Seq("sequence_1", s, "inline sequence")] if s else []


def stats_of_records(recs: list[Seq]) -> list[dict]:
    rows = []
    for r in recs:
        c = counts(r.seq)
        rows.append({
            "id": r.id, "length": len(r.seq), "GC_percent": round(gc_content(r.seq), 3),
            "A": c["A"], "C": c["C"], "G": c["G"], "T": c["T"], "N": c["N"],
            "ambiguity": round(100 * sum(v for k, v in c.items() if k not in "ACGT") /
                               (len(r.seq) or 1), 3),
            "longestORF_nt": max((o["length"] for o in find_orfs(r.seq, 1)), default=0),
            "entropy": round(entropy_of(r.seq, 1, "ACGT"), 4),
            "n_protein": len(translate(r.seq)),
        })
    return rows


def longest_orf_protein(src) -> str:
    recs = read_seq_records(src)
    best = ""
    for r in recs:
        for o in find_orfs(r.seq, 1):
            if len(o["protein"]) > len(best):
                best = o["protein"]
    return best


def sequence_checksum(s: str) -> str:
    return sha1_of(s)


def gc_profile_json(src, size: int = 100, step: int | None = None) -> list[dict]:
    recs = read_seq_records(src)
    if not recs:
        return []
    return windowed_stats(recs[0].seq, size, step)


def find_n_blocks(s: str, min_len: int = 1) -> list[dict]:
    return [r for r in homopolymer_runs(s, min_len) if r["base"] == "N"]


def degenerate_primer_count(pattern: str) -> int:
    mult = {"R": 2, "Y": 2, "S": 2, "W": 2, "K": 2, "M": 2, "B": 3, "D": 3, "H": 3,
            "V": 3, "N": 4}
    n = 1
    for c in clean(pattern):
        n *= mult.get(c, 1)
    return n


def primer3_lite(seq: str, target_len: int = 21, tm_min: float = 58.0, tm_max: float = 62.0,
                 max_prod: int = 400) -> list[dict]:
    """Tiny primer picker: forward/reverse pairs around a target with Tm/product checks."""
    s = back_transcribe(clean(seq))
    out = []
    for i in range(0, max(1, len(s) - target_len)):
        fwd = s[i:i + target_len]
        if len(fwd) < target_len or fwd.count("N") or "ACGT" and not set(fwd) <= set("ACGT"):
            continue
        tf = melting_temp(fwd)
        if not (tm_min <= tf <= tm_max):
            continue
        for span in range(120, min(max_prod, len(s) - i) + 1, 20):
            end = i + span
            rv_region = s[end:end + target_len]
            if len(rv_region) < target_len:
                break
            rev = revcomp(rv_region)
            tr = melting_temp(rev)
            if not (tm_min <= tr <= tm_max):
                continue
            hp = max(abs(hairpin_score(fwd)["delta_g"]), abs(hairpin_score(rev)["delta_g"]))
            dim = self_dimer_score(fwd, rev)["dimer_score"]
            if hp > -6 or dim > 4:
                continue
            out.append({"start": i, "end": end + target_len, "product_size": span,
                        "forward": fwd, "reverse": rev, "tm_forward": round(tf, 2),
                        "tm_reverse": round(tr, 2), "gc_forward": round(gc_content(fwd), 1),
                        "gc_reverse": round(gc_content(rev), 1),
                        "hairpin_dg": round(hp, 2), "dimer": dim})
            break
        if len(out) >= 25:
            break
    return out


def design_primers(seq: str, target_len: int = 22, tm_min: float = 58.0,
                   tm_max: float = 62.0) -> list[dict]:
    return primer3_lite(seq, target_len, tm_min, tm_max)


def restriction_map(src, enzymes: list[str] | None = None) -> list[dict]:
    recs = read_seq_records(src)
    rows = []
    for r in recs:
        for h in cut_positions(r.seq, enzymes):
            rows.append({"sequence": r.id, "length": len(r.seq), **h})
    return rows


def subsequence(seq: str, start: int, end: int, strand: str = "+") -> str:
    c = clean(seq)
    s = c[max(0, start):max(0, end)]
    return revcomp(s) if strand == "-" else s


def extract_region(recs: list[Seq], region: str) -> list[Seq]:
    m = re.match(r"^([^:]+):(\d+)-(\d+)([+-])?$", region.strip())
    if not m:
        m = re.match(r"^([^:]+):(\d+)([+-])?$", region.strip())
        if m:
            region = f"{m.group(1)}:{m.group(2)}-{int(m.group(2)) + 1}"
            return extract_region(recs, region)
        return []
    chrom, a, b, strand = m.group(1), int(m.group(2)), int(m.group(3)), m.group(4) or "+"
    out = []
    for r in recs:
        if r.id != chrom and r.id.split(".")[0] != chrom:
            continue
        sub = r.seq[max(0, a - 1):b]
        if strand == "-":
            sub = revcomp(sub)
        out.append(Seq(f"{chrom}:{a}-{b}{strand if strand != '+' else ''}", sub,
                       f"extracted from {r.id}"))
    return out
