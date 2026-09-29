"""Protein/biochemistry primitives: scales, mass, pI, digestion, motifs."""

from __future__ import annotations

import math
import re
from collections import Counter

from chroma_titan.core.seq import clean

# average residue masses (monoisotopic residues used for peptide mass)
AA_MASS = {
    "A": 71.0788, "R": 156.1011, "N": 114.1038, "D": 115.0886, "C": 103.0092,
    "E": 129.1155, "Q": 128.1307, "G": 57.0519, "H": 137.1411, "I": 113.1594,
    "L": 113.1594, "K": 128.1741, "M": 131.1986, "F": 147.1766, "P": 97.1167,
    "S": 87.0782, "T": 101.1051, "W": 186.2132, "Y": 163.1760, "V": 99.1326,
}
AA_MONO = dict(AA_MASS)
RESIDUE_WATER = 18.0153
H_PROTON = 1.00728
AVG_AA_MASS = {
    "A": 71.08, "R": 156.19, "N": 114.19, "D": 115.09, "C": 103.14, "E": 129.12,
    "Q": 128.13, "G": 57.05, "H": 137.14, "I": 113.16, "L": 113.16, "K": 128.17,
    "M": 131.20, "F": 147.18, "P": 97.12, "S": 87.08, "T": 101.11, "W": 186.21,
    "Y": 163.18, "V": 99.13, "X": 110.0, "B": 113.0, "Z": 128.0, "*": 0.0, "-": 0.0,
}

HYDROPHOBIVITY = {  # Kyte-Doolittle
    "A": 1.8, "R": -4.5, "N": -3.5, "D": -3.5, "C": 2.5, "Q": -3.5, "E": -3.5, "G": -0.4,
    "H": -3.2, "I": 4.5, "L": 3.8, "K": -3.9, "M": 1.9, "F": 2.8, "P": -1.6, "S": -0.8,
    "T": -0.7, "W": -0.9, "Y": -1.3, "V": 4.2,
}
HOPP_WOODS = {k: v / 3.0 for k, v in HYDROPHOBIVITY.items()}  # boundary hydrophobicity
CHARGE = {  # net charge at pH 7 (approx)
    "D": -1, "E": -1, "K": 1, "R": 1, "H": 0.1, "C": 0.04, "Y": -0.01,
}
PK = {  # pKa values (EMBOSS-like)
    "Nterm": 8.6, "Cterm": 3.6, "D": 3.9, "E": 4.07, "H": 6.0, "C": 8.18, "Y": 10.1,
    "K": 10.54, "R": 12.48,
}
VOLUME = {"A": 88.6, "R": 173.4, "N": 114.1, "D": 111.1, "C": 108.5, "Q": 143.8,
          "E": 138.4, "G": 60.1, "H": 153.2, "I": 166.7, "L": 166.7, "K": 168.6,
          "M": 162.9, "F": 189.9, "P": 112.7, "S": 89.0, "T": 116.1, "W": 227.8,
          "Y": 193.6, "V": 140.0}
CHOU_FASMAN = {  # (helix, strand, turn) propensities
    "A": (1.42, 0.66, 0.87), "L": (1.34, 1.37, 0.60), "E": (1.51, 0.37, 0.74),
    "A2": (1.0, 1.0, 1.0), "M": (1.45, 1.05, 0.60), "Q": (1.11, 1.10, 1.10),
    "K": (1.07, 0.74, 1.19), "H": (0.87, 0.62, 1.00), "F": (1.13, 1.38, 0.60),
    "R": (0.98, 0.93, 0.95), "P": (0.12, 0.39, 1.37), "Y": (0.69, 1.47, 1.14),
    "V": (1.33, 1.70, 0.62), "I": (1.08, 1.62, 0.47), "D": (1.01, 0.54, 1.46),
    "T": (0.59, 1.36, 0.77), "W": (0.67, 1.32, 1.12), "N": (0.67, 0.89, 1.56),
    "S": (0.77, 0.75, 1.43), "G": (0.39, 0.93, 1.42), "C": (0.70, 1.19, 1.19),
}
FLEXIBILITY = {"G": 0.054, "P": 0.106, "S": 0.063, "A": 0.060, "T": 0.066, "V": 0.061,
               "L": 0.064, "I": 0.063, "N": 0.081, "D": 0.079, "C": 0.070, "E": 0.070,
               "Q": 0.077, "K": 0.073, "M": 0.058, "H": 0.067, "F": 0.061, "R": 0.080,
               "Y": 0.064, "W": 0.066}
INSTABILITY = {"V": 0, "V2": 0, "H": 3, "D": 2, "L": 1, "R": 0, "M": 2, "F": 2,
               "K": -1, "I": -1, "Y": 1, "C": 2, "N": 0, "P": 1, "G": 0, "S": -1,
               "A": 0, "T": 0, "E": 0, "Q": 0, "W": 2}
POLARITY = {"A": 8.1, "R": 10.5, "N": 11.6, "D": 13.0, "C": 5.5, "Q": 10.5, "E": 12.3,
            "G": 9.0, "H": 10.4, "I": 5.2, "L": 4.9, "K": 11.3, "M": 5.7, "F": 5.2,
            "P": 8.0, "S": 9.2, "T": 8.6, "W": 5.4, "Y": 6.2, "V": 5.9}
WOLF = {"A": 1.05, "R": 1.20, "N": 1.05, "D": 1.05, "C": 1.05, "Q": 1.05, "E": 1.05,
        "G": 1.00, "H": 1.00, "I": 1.00, "L": 1.00, "K": 1.10, "M": 1.05, "F": 1.00,
        "P": 1.10, "S": 1.00, "T": 1.05, "W": 1.00, "Y": 1.05, "V": 1.00}
TRANSMEMBRANE_OCTANOL = HYDROPHOBIVITY

PROSITE_MOTIFS = {
    "N_glycosylation": r"N[^P][ST]",
    "CK2_phosphorylation": r"[ST]..[ED]",
    "PKC_phosphorylation": r"[ST].{0,1}[KR]",
    "myristoylation": r"^G",
    "amideribosylation": r"G{x}",
    "SH3_binding": r"P..[KR]",
    "NLS_classical": r"K[KR]..[KR]",
    "tyrosine_kinase": r"[FYWL]..G.[FYWHM][KR]",
    "EGF_like": r"C.{3,40}C.{20,30}[EQDG][FLIV][RG][CT]C",
    "Zn_finger_C2H2": r"C..C.{12}H.{3,5}H",
    "EF_hand": r"L..[ILV][GSTA].{11}[IVL].{9}[IVL]",
    "thioredoxin": r"W.C.G.P.T.C",
    "Ras_family": r"[GST].[GST].[KR][LIVM]",
    "ATP_binding_ploop": r"G.{0,2}G.K[ST]",
}

PROTEIN_FAMS = {
    "kinase": r"[LIVM].[LIVMFY]{2}[GSTA]",
    "homeobox": r"[KR][KR][VIL].{20}[QE][^C][KR][VIL][FY][KR][KR]",
}


def protein_mass(seq: str, mode: str = "average") -> float:
    s = clean(seq, "ACDEFGHIKLMNPQRSTVWYXBZ*-")
    table = AA_MASS if mode == "mono" else AVG_AA_MASS
    return sum(table.get(c, 110.0) for c in s) + (RESIDUE_WATER if s else 0.0)


def molecular_weight(seq: str, type: str = "protein", mode: str = "average") -> float:
    if type == "dna":
        return dna_mass(seq, mode)
    return protein_mass(seq, mode)


def dna_mass(seq: str, mode: str = "average", stranded: str = "double") -> float:
    water = {"A": 313.21, "C": 289.18, "G": 329.21, "T": 304.2, "U": 306.17, "N": 307.0}
    s = clean(seq, "ACGTUN")
    per = sum(water.get(c, 307.0) + (H2O if False else 18.015) - 18.015 for c in s)
    del per
    total = sum(water.get(c, 307.0) for c in s) + 79.98 * max(0, len(s) - 1)
    return total if stranded == "single" else total * 2 - 2.016


H2O = 18.0153


def peptide_mz(seq: str, charge: int = 1, mode: str = "mono") -> float:
    table = AA_MASS if mode == "mono" else AVG_AA_MASS
    m = sum(table.get(c, 110.0) for c in clean(seq, "ACDEFGHIKLMNPQRSTVWYX")) + RESIDUE_WATER
    return (m + charge * H_PROTON) / charge


def isoelectric_point(seq: str) -> float:
    s = clean(seq, "ACDEFGHIKLMNPQRSTVWY")
    lo, hi = 0.0, 14.0
    for _ in range(80):
        mid = (lo + hi) / 2
        if net_charge(s, mid) > 0:
            lo = mid
        else:
            hi = mid
    return round((lo + hi) / 2, 3)


def net_charge(seq: str, pH: float = 7.0) -> float:
    s = clean(seq, "ACDEFGHIKLMNPQRSTVWYCH")
    q = 0.0
    q += 1 / (1 + 10 ** (pH - PK["Nterm"])) * s.count("K") / max(1, s.count("K")) if False else 0
    for aa, pka in (("K", PK["K"]), ("R", PK["R"]), ("H", PK["H"])):
        q += s.count(aa) / (1 + 10 ** (pH - pka))
    for aa, pka in (("D", PK["D"]), ("E", PK["E"]), ("C", PK["C"]), ("Y", PK["Y"])):
        q -= s.count(aa) / (1 + 10 ** (pka - pH))
    q += 1 / (1 + 10 ** (pH - PK["Nterm"]))
    q -= 1 / (1 + 10 ** (PK["Cterm"] - pH))
    return round(q, 4)


def charge_at_pH(seq: str, pH: float = 7.0) -> float:
    return net_charge(seq, pH)


def gravy(seq: str, window: int = 1, scale: dict | None = None) -> float:
    s = clean(seq, "ACDEFGHIKLMNPQRSTVWY")
    sc = scale or HYDROPHOBIVITY
    if not s:
        return 0.0
    prof = hydrophobicity_profile(s, window, sc)
    return round(sum(prof) / len(prof), 4)


def hydrophobicity_profile(seq: str, window: int = 9, scale: dict | None = None) -> list[float]:
    s = clean(seq, "ACDEFGHIKLMNPQRSTVWYX")
    sc = scale or HYDROPHOBIVITY
    if window % 2 == 0:
        window += 1
    half = window // 2
    vals = [sc.get(c, 0.0) for c in s]
    if len(vals) <= window:
        return [round(sum(vals) / len(vals), 3)] if vals else []
    out = []
    for i in range(len(vals)):
        a, b = max(0, i - half), min(len(vals), i + half + 1)
        out.append(round(sum(vals[a:b]) / (b - a), 4))
    return out


def hydrophobicity_moment(seq: str, window: int = 18, helix: float = 100.0) -> list[float]:
    """Hydrophobic moment of a putative amphipathic helix (Eisenberg)."""
    d = {"A": 0.62, "R": -0.99, "N": -0.78, "D": -0.90, "C": 0.29, "Q": -0.70, "E": -0.74,
         "G": 0.16, "H": -0.40, "I": 1.38, "L": 1.06, "K": -0.99, "M": 0.64, "F": 1.19,
         "P": 0.12, "S": -0.18, "T": -0.05, "W": 0.81, "Y": 0.26, "V": 1.08}
    s = clean(seq, "ACDEFGHIKLMNPQRSTVWY")
    rad = math.radians(helix)
    out = []
    for i in range(0, max(1, len(s) - window + 1), 1):
        seg = s[i:i + window]
        if not seg:
            continue
        h = sum(d.get(c, 0.0) for c in seg) / len(seg)
        mu_r = complex(sum(d.get(c, 0.0) * math.cos(k * rad) for k, c in enumerate(seg)),
                       sum(d.get(c, 0.0) * math.sin(k * rad) for k, c in enumerate(seg)))
        out.append(round(abs(mu_r) / len(seg), 4))
        del h
    return out


def instability_index(seq: str) -> float:
    s = clean(seq, "ACDEFGHIKLMNPQRSTVWY")
    if len(s) < 2:
        return 0.0
    tot = 0.0
    for i in range(len(s) - 1):
        di = INSTABILITY.get(s[i])
        if di is not None:
            tot += di
    return round(10.0 * tot / max(1, len(s) - 1), 3)


def aliphatic_index(seq: str) -> float:
    s = clean(seq, "ACDEFGHIKLMNPQRSTVWY")
    n = len(s) or 1
    a, v, i_, l = (s.count(x) * 100 / n for x in "AVIL")
    bturn = (i_ + l) / (i_ + l + v) if (i_ + l + v) else 0
    return round(a + 2.9 * v + 3.98 * bturn, 3)


def aromaticity(seq: str) -> float:
    s = clean(seq, "ACDEFGHIKLMNPQRSTVWY")
    return round(100 * sum(s.count(c) for c in "FWY") / (len(s) or 1), 3)


def molecular_formula(seq: str) -> str:
    comp = {
        "A": "C3H5NO", "R": "C6H14N4O", "N": "C4H6N2O", "D": "C4H5NO3", "C": "C3H5NOS",
        "E": "C5H7NO3", "Q": "C5H8N2O2", "G": "C2H3NO", "H": "C6H7N3O", "I": "C6H11NO",
        "L": "C6H11NO", "K": "C6H12N2O", "M": "C5H9NOS", "F": "C9H9NO", "P": "C5H7NO",
        "S": "C3H5NO2", "T": "C4H7NO2", "W": "C11H10N2O", "Y": "C9H9NO2", "V": "C5H9NO",
    }
    tot: Counter = Counter()
    for c in clean(seq, "ACDEFGHIKLMNPQRSTVWY"):
        for tok in re.findall(r"([A-Z])(\d*)", comp.get(c, "")):
            tot[tok[0]] += int(tok[1] or 1)
    tot["H"] += 1
    tot["O"] += 1
    order = ["C", "H", "N", "O", "S"]
    return "".join(f"{e}{tot[e] if tot[e] != 1 else ''}" for e in order if tot[e])


def aa_composition(seq: str) -> dict[str, float]:
    s = clean(seq, "ACDEFGHIKLMNPQRSTVWYXBZ*")
    n = len(s) or 1
    c = Counter(s)
    return {aa: round(100 * c[aa] / n, 3) for aa in sorted(c)}


def protein_length(seq: str) -> int:
    return len(clean(seq, "ACDEFGHIKLMNPQRSTVWYXBZ*-"))


def secondary_structure(seq: str) -> str:
    """Chou-Fasman style prediction returning H/E/C string."""
    s = clean(seq, "ACDEFGHIKLMNPQRSTVWY")
    if not s:
        return ""
    hh = [0.0] * len(s)
    ee = [0.0] * len(s)
    for i, c in enumerate(s):
        p = CHOU_FASMAN.get(c, (1, 1, 1))
        for j in range(max(0, i - 3), min(len(s), i + 4)):
            hh[j] += p[0]
            ee[j] += p[1]
    out = []
    for i, c in enumerate(s):
        if hh[i] / 7 > 1.03 and hh[i] > ee[i]:
            out.append("H")
        elif ee[i] / 7 > 1.0 and ee[i] >= hh[i]:
            out.append("E")
        else:
            out.append("C")
    return "".join(out)


def disorder_propensity(seq: str) -> list[float]:
    """IUPred-lite: low complexity + charge/hydropathy balance per window."""
    s = clean(seq, "ACDEFGHIKLMNPQRSTVWY")
    if not s:
        return []
    win = 21
    out = []
    for i in range(len(s)):
        a, b = max(0, i - win // 2), min(len(s), i + win // 2 + 1)
        seg = s[a:b]
        n = len(seg) or 1
        charged = sum(seg.count(c) for c in "DEKR") / n
        order = sum(HYDROPHOBIVITY.get(c, 0.0) ** 2 for c in seg) / n
        out.append(round(min(1.0, max(0.0, 0.9 - charged * 0.6 - math.sqrt(order) * 0.15
                                    + 0.25 * (seg.count("P") + seg.count("G")) / n)), 3))
    return out


def transmembrane_helices(seq: str, window: int = 19, cutoff: float = 1.6) -> list[dict]:
    prof = hydrophobicity_profile(seq, window, TRANSMEMBRANE_OCTANOL)
    s = clean(seq, "ACDEFGHIKLMNPQRSTVWY")
    out, inside, start = [], False, 0
    for i, v in enumerate(prof):
        if v > cutoff and not inside:
            inside, start = True, i
        elif v <= cutoff and inside:
            inside = False
            out.append({"start": max(0, start - window // 2), "end": min(len(s), i + window // 2),
                        "mean_hydrophobicity": round(sum(prof[start:i]) / max(1, i - start), 3)})
    if inside:
        out.append({"start": start, "end": len(prof), "mean_hydrophobicity": 0.0})
    return out


def signal_peptide_score(seq: str) -> dict:
    """GVH-hydrophobicity test in the first 25 residues (von Heijne-lite)."""
    s = clean(seq, "ACDEFGHIKLMNPQRSTVWY")
    n_region = s[:25]
    if len(n_region) < 10:
        return {"score": 0.0, "predicted": False, "window": n_region}
    gvh = (sum(TRANSMEMBRANE_OCTANOL.get(c, 0.0) for c in n_region) +
           2 * n_region.count("G") + 2 * n_region.count("V") + n_region.count("H")) / len(n_region)
    has_clv = bool(re.search(r"[AVST][AVST]X[AVST]", s[:35])) and n_region[:2].count("K") == 0
    score = round(min(1.0, max(0.0, (gvh - 0.4) / 1.6)) + (0.25 if has_clv else 0.0), 3)
    return {"score": score, "predicted": score >= 0.55, "window": n_region,
            "cleavage_site_possible": has_clv}


def disulfides(seq: str) -> list[dict]:
    s = clean(seq, "ACDEFGHIKLMNPQRSTVWY")
    pos = [i for i, c in enumerate(s) if c == "C"]
    return [{"cysteine_1": a + 1, "cysteine_2": b + 1, "separation": b - a}
            for a, b in zip(pos, pos[1:]) if 3 <= b - a <= 25]


def n_glyc_sites(seq: str) -> list[dict]:
    return [{"position": m.start() + 1, "motif": m.group()}
            for m in re.finditer(r"N[^P][ST]", clean(seq, "ACDEFGHIKLMNPQRSTVWY"))]


def prosite_sites(seq: str, families: list[str] | None = None) -> list[dict]:
    s = clean(seq, "ACDEFGHIKLMNPQRSTVWY")
    out = []
    for name, pat in PROSITE_MOTIFS.items():
        if families and name not in families:
            continue
        for m in re.finditer(pat, s):
            out.append({"motif": name, "pattern": pat, "start": m.start() + 1,
                        "end": m.end(), "match": m.group()})
    return out


PROTEASES = {
    "trypsin": (r"(?<=[KR])(?!P)", "C-term"),
    "chymotrypsin": (r"(?<=[FYWHM])(?!P)", "C-term"),
    "arg-c": (r"(?<=R)(?!P)", "C-term"),
    "lys-c": (r"(?<=K)(?!P)", "C-term"),
    "gluc": (r"(?<=D)", "C-term"),
    "asp-n": (r"(?=[DE])", "N-term"),
    "pepsin_ph1.3": (r"(?<=[FL])(?=.)", "C-term"),
    "proteinase_k": (r"(?<=[ACT])\b", "C-term"),
    "thermolysin": (r"(?=[FLIVMY])(?!P)", "N-term"),
    "clostripain": (r"(?<=R)", "C-term"),
    "staphylococcal_peptidase": (r"(?<=[DE])(?!P)", "C-term"),
    "bnps_hydroxylase": (r"(?<=[FYWH])", "N-term"),
}


def digest(seq: str, enzyme: str = "trypsin", min_len: int = 6, max_len: int = 4000,
            missed: int = 0) -> list[str]:
    s = clean(seq, "ACDEFGHIKLMNPQRSTVWYXBZ")
    pat, side = PROTEASES.get(enzyme, PROTEASES["trypsin"])
    parts = [p for p in re.split(pat, s) if p]
    if missed and side == "C-term":
        for _ in range(missed):
            if len(parts) > 1:
                merged = []
                for i in range(0, len(parts) - 1, 2):
                    merged.append(parts[i] + (parts[i + 1] if i + 1 < len(parts) else ""))
                if len(parts) % 2:
                    merged.append(parts[-1])
                parts = merged
    out = [p for p in parts if min_len <= len(p) <= max_len]
    return out


def peptide_masses(peptides: list[str], mode: str = "mono") -> list[float]:
    return [round(peptide_mz(p, 1, mode), 4) for p in peptides]


def mass_to_charge(mass: float, charge: int = 1) -> float:
    return round((mass + charge * H_PROTON) / max(1, charge), 4)


def peptide_from_mass(target: float, seq: str, tolerance: float = 0.5,
                      enzyme: str = "trypsin") -> list[dict]:
    out = []
    for p in digest(seq, enzyme, min_len=4):
        m = peptide_mz(p, 1, "mono")
        if abs(m - target) <= tolerance:
            out.append({"peptide": p, "mass": round(m, 4), "delta": round(m - target, 4),
                        "length": len(p), "score": 100 / (1 + abs(m - target))})
    return sorted(out, key=lambda x: abs(x["delta"]))


def pseudo_neutral_mass(seq: str) -> float:
    return round(protein_mass(seq, "mono") + sum(
        1 for c in clean(seq) if c in "KR") * 0.0, 4)


def amino_acid_property_table() -> list[dict]:
    rows = []
    one = {"Alanine": "A", "Arginine": "R", "Asparagine": "N", "Aspartate": "D",
           "Cysteine": "C", "Glutamine": "Q", "Glutamate": "E", "Glycine": "G",
           "Histidine": "H", "Isoleucine": "I", "Leucine": "L", "Lysine": "K",
           "Methionine": "M", "Phenylalanine": "F", "Proline": "P", "Serine": "S",
           "Threonine": "T", "Tryptophan": "W", "Tyrosine": "Y", "Valine": "V"}
    three = {v: k for k, v in one.items()}
    for name, code in sorted(one.items()):
        rows.append({"amino_acid": name, "one_letter": code, "three_letter": three[code],
                     "monoisotopic_mass": AA_MASS[code], "average_mass": AVG_AA_MASS[code],
                     "hydrophobicity": HYDROPHOBIVITY[code], "charge_pH7": CHARGE.get(code, 0.0),
                     "helix_propensity": CHOU_FASMAN[code][0], "strand_propensity": CHOU_FASMAN[code][1],
                     "turn_propensity": CHOU_FASMAN[code][2], "residue_volume": VOLUME[code],
                     "flexibility": FLEXIBILITY[code], "polarity": POLARITY[code]})
    return rows


def charge_histogram(seq: str, pH_range: tuple[float, float] = (0.0, 14.0), steps: int = 28) -> list[dict]:
    lo, hi = pH_range
    return [{"pH": round(lo + (hi - lo) * i / (steps - 1), 3),
             "net_charge": net_charge(seq, lo + (hi - lo) * i / (steps - 1))}
            for i in range(steps)]


def buffer_capacity_estimate(seq: str) -> float:
    s = clean(seq, "ACDEFGHIKLMNPQRSTVWY")
    return float(sum(s.count(c) for c in "DEKRH"))


def protein_alignment_identity(a: str, b: str) -> float:
    n = min(len(a), len(b))
    if not n:
        return 0.0
    return round(100 * sum(1 for x, y in zip(a, b) if x == y) / n, 3)


def cys_count(seq: str) -> dict:
    s = clean(seq, "ACDEFGHIKLMNPQRSTVWY")
    return {"free_cysteines": s.count("C"), "count_C": s.count("C"),
            "count_M": s.count("M"), "count_W": s.count("W")}


def extinction_coefficient(seq: str, pH: float = 8.0) -> dict:
    s = clean(seq, "ACDEFGHIKLMNPQRSTVWYCH")
    trp, tyr, cyst = s.count("W"), s.count("Y"), s.count("C")
    ext = 5500 * trp + 1490 * tyr
    red = ext + 125 * cyst
    mw = protein_mass(s, "average") or 1.0
    return {"extinction_reduced": red, "extinction_oxidised": ext,
            "A1M_1cm_reduced": round(red / 1000, 4),
            "E1percent_1cm_reduced": round(red / mw * 10, 4),
            "estimated_concentration_uM_for_A1": round(1e6 / red, 3) if red else 0.0,
            "aromatic_residues": trp + tyr, "half_life_at_pH": 0.0}


def charge_density(seq: str) -> float:
    s = clean(seq, "ACDEFGHIKLMNPQRSTVWY")
    return round(sum(CHARGE.get(c, 0.0) for c in s) / (len(s) or 1), 5)


def net_charge_profile(seq: str, window: int = 21) -> list[float]:
    s = clean(seq, "ACDEFGHIKLMNPQRSTVWY")
    half = window // 2
    out = []
    for i in range(len(s)):
        seg = s[max(0, i - half):min(len(s), i + half + 1)]
        out.append(round(sum(CHARGE.get(c, 0.0) for c in seg) / len(seg), 4))
    return out
