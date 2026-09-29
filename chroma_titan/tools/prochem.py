"""Proteomics, metabolomics, cheminformatics, pharmacology and single-cell import tools.

Panel sections covered: **proteomics**, **metabolomics**, **chemicaltoolbox**,
**pharmacology**, **multiomics**, **import/manipulate_sc_data**, **hca-scanpy**,
**monocle3**, **presto**, **seurat**.
"""
from __future__ import annotations

import math
import re
from collections import Counter, defaultdict

import numpy as np
import pandas as pd

from chroma_titan.core import align, image, io, ml, motif, plot, protein, seq, stats, tables
from chroma_titan.tools._common import *  # noqa: F401,F403
from chroma_titan.tools._common import (COUNTS, CXCG, GENES_FA, GMT, IMAGE, PEPTIDES, PHENO, PROTEINS,
                                        REF_PROTEINS, SHORT_PROT)

PROT = "proteomics"
METAB = "metabolomics"
CHEM = "chemicaltoolbox"
PHARM = "pharmacology"
MULTI = "multiomics"
SCIMP = "import/manipulate_sc_data"
SCANPY = "hca-scanpy"
MONO = "monocle3"
PRESTO = "presto"
SEURAT = "seurat"


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------
def _fa(src):
    return io.parse_fasta(io.as_text(src))


def _tbl(src):
    return tables.load(io.as_text(src))


def _num_cols(df):
    out = []
    for c in df.columns:
        v = pd.to_numeric(df[c], errors="coerce").to_numpy(dtype=float)
        if np.isfinite(v).sum() >= 2:
            out.append((str(c), v))
    return out


def _mat(src):
    """Wide numeric table -> (row labels, column names, ndarray)."""
    df = _tbl(src)
    if df.empty:
        return [], [], np.zeros((0, 0))
    cols = list(df.columns)
    body = df.drop(columns=[cols[0]])
    data = np.column_stack([pd.to_numeric(body[c], errors="coerce").to_numpy(dtype=float)
                           for c in body.columns]) if len(body.columns) else np.zeros((len(body), 0))
    return [str(x) for x in df[cols[0]].tolist()], [str(c) for c in body.columns], data


def _groups_from_spec(labels, spec):
    """``name=a,b;name2=c`` -> {name: [indices]}."""
    out: dict[str, list[int]] = {}
    pos = {str(v): i for i, v in enumerate(labels)}
    for chunk in str(spec or "").replace("\n", ";").split(";"):
        if not chunk.strip() or "=" not in chunk:
            continue
        name, members = chunk.split("=", 1)
        ids = [pos[m.strip()] for m in members.split(",") if m.strip() in pos]
        if ids:
            out[name.strip()] = ids
    return out


AA3 = {"A": "Ala", "R": "Arg", "N": "Asn", "D": "Asp", "C": "Cys", "Q": "Gln", "E": "Glu", "G": "Gly",
       "H": "His", "I": "Ile", "L": "Leu", "K": "Lys", "M": "Met", "F": "Phe", "P": "Pro", "S": "Ser",
       "T": "Thr", "W": "Trp", "Y": "Tyr", "V": "Val"}


def _peptide_rows(src, seq_col="", intensity_col="", protein_col="", rt_col="", charge_col=""):
    """Peptide table -> list of normalised dicts."""
    df = _tbl(src)
    if df.empty:
        return []
    cols = [str(c) for c in df.columns]

    def pick(pref, default=""):
        for c in cols:
            if pref in c.lower():
                return c
        return default
    sc = seq_col or pick("peptide", pick("seq", cols[0]))
    ic = intensity_col or pick("intens", pick("abund", pick("area")))
    pc = protein_col or pick("protein")
    rc = rt_col or pick("rt", pick("time"))
    cc = charge_col or pick("charge")
    rows = []
    for r in df.to_dict("records"):
        s = str(r.get(sc, "")).strip()
        if not s:
            continue
        try:
            inten = float(r.get(ic)) if ic else float("nan")
        except (TypeError, ValueError):
            inten = float("nan")
        try:
            rt = float(r.get(rc)) if rc else float("nan")
        except (TypeError, ValueError):
            rt = float("nan")
        try:
            ch = int(float(r.get(cc))) if cc else 2
        except (TypeError, ValueError):
            ch = 2
        rows.append({"peptide": s, "protein": str(r.get(pc, "unknown")) if pc else "unknown",
                    "intensity": inten, "runtime": rt, "charge": max(1, ch),
                    "modified": str(r.get("modified_sequence", s))})
    return rows


# ===========================================================================
# proteomics
# ===========================================================================
@T("proteo_digest_search", "Match observed peptides to an in-silico digest", PROT, "table",
   [fa("protein_file", PROTEINS, "Protein FASTA"), tbl("spectrum_table", PEPTIDES, "Observed peptide table"),
    choice("enzyme", ["trypsin", "lys-c", "arg-c", "glu-c", "chymotrypsin", "none"], "trypsin", "Enzyme"),
    number("ppm_tolerance", 50.0, "Mass tolerance (ppm)", min=0.1), intin("missed_cleavages", 1,
                                                                        "Missed cleavages", min=0),
    number("min_mass", 300.0, "Minimum peptide mass (Da)", min=0.0)],
   ex={"protein_file": PROTEINS, "spectrum_table": PEPTIDES, "ppm_tolerance": 5000.0},
   up="MS-Fragger / PEAKS / OpenMS FidoAdapter", tags=("database search", "digest", "ms/ms"),
   summary="Digest proteins in silico and match each observed peptide by mass tolerance.")
def proteo_digest_search(protein_file, spectrum_table, enzyme="trypsin", ppm_tolerance=50.0,
                         missed_cleavages=1, min_mass=300.0):
    """In-silico digest matching."""
    obs = _peptide_rows(spectrum_table)
    lib = []
    for rec in _fa(protein_file):
        for pep in protein.digest(rec.seq, enzyme=enzyme, min_len=6, max_len=4000,
                                 missed=int(missed_cleavages)):
            m = float(protein.protein_mass(pep, mode="mono"))
            if m >= float(min_mass):
                lib.append((m, pep, rec.id))
    rows = []
    matched = 0
    for o in obs:
        target = float(protein.protein_mass(o["peptide"], mode="mono"))
        tol = target * float(ppm_tolerance) / 1e6
        cands = [x for x in lib if abs(x[0] - target) <= max(tol, 0.5)]
        exact = sum(1 for m, p, _pid in lib if p == o["peptide"])
        protein_hit = next((pid for _m, p, pid in lib if p == o["peptide"]),
                          cands[0][2] if cands else "unmatched")
        matched += 1 if cands or exact else 0
        rows.append({"peptide": o["peptide"], "charge": o["charge"], "theoretical_mass": round(target, 4),
                    "observed_intensity": (round(o["intensity"], 1) if math.isfinite(o["intensity"]) else ""),
                    "protein": protein_hit, "candidates": len(cands), "exact_match": bool(exact),
                    "delta_ppm": round(1e6 * (min([abs(m - target) for m, _p, _i in cands], default=0.0))
                                     / (target or 1.0), 3)})
    return table(rows, f"{matched}/{len(obs)} observed peptides matched a theoretical digest "
                     f"({len(lib)} candidates, +/- {ppm_tolerance} ppm)")


@T("proteo_peptide_properties", "Physicochemical properties of peptides", PROT, "table",
   [tbl("spectrum_table", PEPTIDES, "Peptide table"), textbox("sequence_column", "", "Sequence column"),
    intin("charge", 2, "Default charge state", min=1), boolean("with_mods", False, "Use modified sequences"),
    intin("head", 200, "Rows", min=1)],
   ex={"spectrum_table": PEPTIDES}, up="ProtParam / pyMSpec", tags=("peptide properties", "m/z"),
   summary="Mass, m/z, pI, net charge, hydrophobicity and instability for each peptide.")
def proteo_peptide_properties(spectrum_table, sequence_column="", charge=2, with_mods=False, head=200):
    """Per-peptide property table."""
    obs = _peptide_rows(spectrum_table, seq_col=sequence_column)
    rows = []
    for o in obs[: int(head)]:
        s = protein.clean(o["modified"] if with_mods else o["peptide"])
        if not s:
            continue
        mz = float(protein.peptide_mz(s, charge=max(1, int(charge)), mode="mono"))
        rows.append({"peptide": o["peptide"], "length": len(s), "mono_mass": round(float(protein.molecular_weight(s, mode="mono")), 4),
                    "average_mass": round(float(protein.molecular_weight(s, mode="average")), 4),
                    "mz": round(mz, 4), "charge": int(charge),
                    "theoretical_pI": round(float(protein.isoelectric_point(s)), 3),
                    "net_charge_pH7": round(float(protein.net_charge(s, pH=7.0)), 3),
                    "gravy": round(float(protein.gravy(s)), 4),
                    "instability_index": round(float(protein.instability_index(s)), 3),
                    "cysteines": s.count("C"), "methionines": s.count("M"),
                    "aromaticity": round(float(protein.aromaticity(s)), 4),
                    "formula": protein.molecular_formula(s)})
    return table(rows, f"{len(rows)} peptides characterised")


@T("proteo_fdr_from_decoys", "Target-decoy false discovery rate", PROT, "table",
   [tbl("spectrum_table", PEPTIDES, "PepScore table"), textbox("score_column", "", "Score column"),
    textbox("target_prefix", "sp|", "Target id prefix"), textbox("decoy_prefix", "decoy_", "Decoy prefix"),
    number("fdr_cutoff", 0.01, "Report scores up to this FDR", min=0.0001, max=0.5)],
   ex={"spectrum_table": PEPTIDES, "score_column": "raw_intensity"}, up="Percolator / Prophet PSM FDR",
   tags=("FDR", "decoys", "quality control"),
   summary="Rank PSMs by score and report the decoy-based FDR at every threshold.")
def proteo_fdr_from_decoys(spectrum_table, score_column="", target_prefix="sp|", decoy_prefix="decoy_",
                           fdr_cutoff=0.01):
    """Target/decoy FDR curve."""
    df = _tbl(spectrum_table)
    cols = list(df.columns)
    sc = score_column or (cols[0] if cols else "")
    raw = pd.to_numeric(df[sc], errors="coerce").to_numpy(dtype=float) if sc in df.columns \
        else pd.to_numeric(df.iloc[:, 0], errors="coerce").to_numpy(dtype=float)
    raw = np.asarray(raw, dtype=float)
    labels = [str(v) for v in (df[cols[-1]].tolist() if cols else [])]
    is_decoy = [bool(decoy_prefix) and decoy_prefix.lower() in l.lower() for l in labels]
    order = np.argsort(-np.nan_to_num(raw, nan=-np.inf))
    t = d = 0
    rows = []
    for rank, i in enumerate(order, 1):
        if is_decoy[i]:
            d += 1
        else:
            t += 1
        fdr = d / max(1, t)
        rows.append({"rank": rank, "score": round(float(raw[i]), 4) if math.isfinite(raw[i]) else "nan",
                    "targets": t, "decoys": d, "psms": t + d, "fdr": round(fdr, 5),
                    "pass": bool(fdr <= float(fdr_cutoff))})
    cut = next((r for r in reversed(rows) if r["pass"]), None)
    res = table(rows[:300], f"{cut['targets'] if cut else 0} targets accepted at FDR <= {fdr_cutoff} "
                          f"({len(rows)} PSMs, {sum(is_decoy)} decoys)")
    res["stats"] = {"decoys": int(sum(is_decoy)), "targets": int(len(rows) - sum(is_decoy)),
                   "figure": plot.line({"FDR": [r["fdr"] for r in rows], "targets": [r["targets"] for r in rows]},
                                     x=[r["rank"] for r in rows], xlabel="PSM rank", ylabel="value",
                                     title="Target-decoy FDR curve")}
    return res


@T("proteo_normalise_median", "Median-intensity normalisation of a protein table", PROT, "table",
   [tbl("abundance", COUNTS, "Protein abundance matrix"), choice("mode", ["median", "total", "rms", "none"],
                                                              "median", "Normalisation"),
    number("log2", True, "Report in log2 scale", min=0.0, max=1.0)],
   ex={"abundance": COUNTS, "mode": "median"}, up="DEP normalisation / MaxQuant",
   tags=("normalisation", "quantification"),
   summary="Scale every sample so protein intensities are comparable across runs.")
def proteo_normalise_median(abundance, mode="median", log2=True):
    """Normalisation factors and normalised matrix."""
    names, cols, X = _mat(abundance)
    X = np.nan_to_num(X.astype(float), nan=0.0)
    if X.size == 0:
        return table([], "no numeric columns")
    if mode == "median":
        factors = np.nanmedian(np.where(X > 0, X, np.nan), axis=0)
    elif mode == "total":
        factors = np.nansum(X, axis=0)
    elif mode == "rms":
        factors = np.sqrt(np.nanmean(np.where(X > 0, X ** 2, np.nan), axis=0))
    else:
        factors = np.ones(X.shape[1])
    factors = np.where(np.isfinite(factors) & (factors > 0), factors, 1.0)
    target = float(np.mean(factors))
    scaled = X / factors * target
    rows = []
    for i, nm in enumerate(names):
        rec = {"protein": nm}
        for j, c in enumerate(cols):
            v = float(scaled[i, j])
            rec[c] = round(math.log2(v + 1.0), 4) if log2 and v > 0 else round(v, 4)
        rows.append(rec)
    res = table(rows, f"{mode} normalisation over {len(cols)} samples "
                   f"(factors {', '.join(f'{c}={factors[j]:.2f}' for j, c in enumerate(cols[:5]))})")
    res["stats"] = {"factors": {c: round(float(factors[j]), 5) for j, c in enumerate(cols)}}
    return res


@T("proteo_quant_summary", "Quantitation summary per protein", PROT, "table",
   [tbl("spectrum_table", PEPTIDES, "Peptide table"), textbox("protein_column", "", "Protein column"),
    textbox("intensity_column", "", "Intensity column"), boolean("share", True, "Estimate shared peptides"),
    intin("head", 100, "Rows", min=1)],
   ex={"spectrum_table": PEPTIDES}, up="MaxQuant proteinGroups / FlashLFQ", tags=("quantification", "summary"),
   summary="Peptide counts, total and unique intensity per protein group.")
def proteo_quant_summary(spectrum_table, protein_column="", intensity_column="", share=True, head=100):
    """Protein-level roll-up."""
    obs = _peptide_rows(spectrum_table, protein_col=protein_column, intensity_col=intensity_column)
    if not obs:
        return table([], "no peptide rows parsed")
    by_prot: dict[str, list[dict]] = defaultdict(list)
    peptide_count = Counter(o["peptide"] for o in obs)
    for o in obs:
        by_prot[o["protein"]].append(o)
    rows = []
    for prot, items in by_prot.items():
        inten = [o["intensity"] for o in items if math.isfinite(o["intensity"])]
        rts = [o["runtime"] for o in items if math.isfinite(o["runtime"])]
        uniq = sum(1 for o in items if peptide_count[o["peptide"]] == 1)
        rows.append({"protein": prot, "peptides": len(items), "unique_peptides": uniq,
                    "shared_peptides": len(items) - uniq if share else 0,
                    "total_intensity": round(float(np.sum(inten)), 1) if inten else 0.0,
                    "mean_intensity": round(float(np.mean(inten)), 1) if inten else 0.0,
                    "median_intensity": round(float(np.median(inten)), 1) if inten else 0.0,
                    "runtime_min": f"{round(min(rts), 2)}-{round(max(rts), 2)}" if rts else "",
                    "charges": ",".join(str(c) for c in sorted({o["charge"] for o in items}))})
    rows.sort(key=lambda r: -r["total_intensity"])
    res = table(rows[: int(head)], f"{len(rows)} protein groups from {len(obs)} peptide identifications")
    res["stats"]["figure"] = plot.bar([r["protein"] for r in rows[:12]][::-1],
                                     [r["total_intensity"] for r in rows[:12]][::-1],
                                     xlabel="protein", ylabel="total intensity",
                                     title="Most abundant proteins", horizontal=True)
    return res


@T("proteo_differential_abundance", "Differential protein abundance between groups", PROT, "table",
   [tbl("abundance", COUNTS, "Protein abundance matrix"), textbox("groups", "", "Group definition (a=1,2;b=3,4)"),
    choice("test", ["ttest", "welch", "mannwhitney"], "ttest", "Test"), number("fc_cutoff", 2.0,
                                                            "Fold-change cut-off", min=0.0),
    number("alpha", 0.05, "FDR cut-off", min=0.0, max=1.0)],
   ex={"abundance": COUNTS, "groups": "control=sample_A,sample_B,sample_C;trt=sample_D,sample_E,sample_F",
       "fc_cutoff": 1.5}, up="DEP / OpenMS DiffEnrich", tags=("differential abundance", "statistics"),
   summary="Per-protein fold change with a p-value and significance flags for two groups.")
def proteo_differential_abundance(abundance, groups="", test="ttest", fc_cutoff=2.0, alpha=0.05):
    """Differential abundance."""
    names, cols, X = _mat(abundance)
    grp = _groups_from_spec(cols, groups) or {"A": list(range(len(cols) // 2)),
                                            "B": list(range(len(cols) // 2, len(cols)))}
    keys = list(grp)
    if len(keys) < 2:
        return table([], "need two groups")
    g1, g2 = grp[keys[0]], grp[keys[1]]
    rows, ps = [], []
    for i, nm in enumerate(names):
        a = X[i, g1][np.isfinite(X[i, g1])]
        b = X[i, g2][np.isfinite(X[i, g2])]
        if a.size < 2 or b.size < 2:
            continue
        if test == "mannwhitney":
            st = stats.mann_whitney_u(list(a), list(b))
            p, stat = float(st["p_value"]), float(st["U"])
        elif test == "welch":
            st = stats.ttest_ind(list(b), list(a), equal_var=False)
            p, stat = float(st["p_value"]), float(st["statistic"])
        else:
            st = stats.ttest_ind(list(b), list(a))
            p, stat = float(st["p_value"]), float(st["statistic"])
        if not math.isfinite(p):
            p = 1.0
        ps.append(p)
        fc = float(np.mean(b)) / (float(np.mean(a)) or 1e-9)
        rows.append({"protein": nm, f"mean_{keys[0]}": round(float(np.mean(a)), 4),
                    f"mean_{keys[1]}": round(float(np.mean(b)), 4),
                    "fold_change": round(fc, 4), "log2_fold_change": round(math.log2(fc) if fc > 0 else 0.0, 4),
                    "statistic": round(stat, 4), "p_value": p})
    adj = stats.p_adjust(ps, "fdr_bh")
    for i, r in enumerate(rows):
        r["p_adjusted"] = round(float(adj[i]), 6) if i < len(adj) else 1.0
        r["regulation"] = ("up" if r["fold_change"] >= float(fc_cutoff) else
                          "down" if r["fold_change"] <= 1.0 / float(fc_cutoff) else "unchanged")
        r["significant"] = bool(r["p_adjusted"] < float(alpha) and r["regulation"] != "unchanged")
    rows.sort(key=lambda r: r["p_value"])
    res = table(rows, f"{sum(1 for r in rows if r['significant'])} proteins regulated at "
                   f"FDR < {alpha} and |FC| >= {fc_cutoff}")
    res["stats"]["figure"] = plot.volcano([r["log2_fold_change"] for r in rows], [r["p_value"] for r in rows],
                                        labels=[r["protein"] for r in rows], title="Differential abundance",
                                        lfc_cut=math.log2(float(fc_cutoff) or 2.0), p_cut=float(alpha))
    return res


@T("proteo_sequence_coverage", "Sequence coverage from identified peptides", PROT, "table",
   [fa("protein_file", PROTEINS, "Protein FASTA"), tbl("spectrum_table", PEPTIDES, "Peptide table"),
    textbox("protein_column", "", "Protein column"), boolean("with_mods", False, "Strip modification marks")],
   ex={"protein_file": PROTEINS, "spectrum_table": PEPTIDES}, up="PeptideShaker / coverage plot",
   tags=("coverage", "sequence"),
   summary="Share of each protein explained by identified peptides, plus uncovered gaps.")
def proteo_sequence_coverage(protein_file, spectrum_table, protein_column="", with_mods=False):
    """Coverage per protein."""
    obs = _peptide_rows(spectrum_table, protein_col=protein_column)
    pep_by_prot: dict[str, list[str]] = defaultdict(list)
    for o in obs:
        s = o["modified"] if with_mods else o["peptide"]
        pep_by_prot[o["protein"]].append(re.sub(r"\[.*?\]|\(.*?\)", "", s).upper())
    rows = []
    for rec in _fa(protein_file):
        seqs = protein.clean(rec.seq)
        covered = np.zeros(len(seqs), dtype=bool)
        starts = []
        for pep in pep_by_prot.get(rec.id, []) + pep_by_prot.get("unknown", []):
            i = seqs.find(pep)
            while i >= 0:
                covered[i:i + len(pep)] = True
                starts.append(i + 1)
                i = seqs.find(pep, i + 1)
        gaps = []
        k = 0
        while k < len(covered):
            if not covered[k]:
                j = k
                while j < len(covered) and not covered[j]:
                    j += 1
                if j - k >= 10:
                    gaps.append(f"{k + 1}-{j}")
                k = j
            else:
                k += 1
        rows.append({"protein": rec.id, "length": len(seqs), "peptides": len(pep_by_prot.get(rec.id, [])),
                    "covered_residues": int(covered.sum()),
                    "coverage_percent": round(100.0 * covered.sum() / max(1, len(seqs)), 3),
                    "first_matched_positions": ",".join(str(x) for x in sorted(set(starts))[:8]),
                    "uncovered_blocks": ",".join(gaps[:6]) or "none"})
    rows.sort(key=lambda r: -r["coverage_percent"])
    return table(rows, f"{len(rows)} proteins; mean coverage "
                     f"{round(float(np.mean([r['coverage_percent'] for r in rows])), 2) if rows else 0}%")


@T("proteo_ptm_sites", "Potential modification sites in proteins", PROT, "table",
   [fa("protein_file", PROTEINS, "Protein FASTA"),
    multi("modifications", ["N-glycosylation", "phospho-ST", "phospho-Y", "acetylation-K", "ubiquitination-K",
                          "oxidation-M", "amidation-Cterm"], ["N-glycosylation", "phospho-ST", "oxidation-M"],
        "Modifications"), boolean("context", True, "Report local sequence context")],
   ex={"protein_file": PROTEINS}, up="Uniprot PTM scan / SeQuest PTM",
   tags=("PTM", "motifs"),
   summary="Sequence-motif based prediction of glycosylation, phosphorylation and oxidation sites.")
def proteo_ptm_sites(protein_file, modifications=("N-glycosylation", "phospho-ST", "oxidation-M"), context=True):
    """PTM site table."""
    pats = {"N-glycosylation": r"N[^P][ST]", "phospho-ST": r"[ST]", "phospho-Y": r"Y",
           "acetylation-K": r"K", "ubiquitination-K": r"K", "oxidation-M": r"M", "amidation-Cterm": r"Q$"}
    mods = [str(m) for m in (modifications or [])] or list(pats)
    rows = []
    for rec in _fa(protein_file):
        s = protein.clean(rec.seq)
        glyc = {int(dict(d).get("position", 0) or 0) for d in protein.n_glyc_sites(s)}
        for mod in mods:
            pat = pats.get(mod)
            if not pat:
                continue
            found = [m.start() + 1 for m in re.finditer(pat, s)] if mod != "N-glycosylation" else sorted(glyc)
            for pos in found[:400]:
                rows.append({"protein": rec.id, "modification": mod, "position": int(pos),
                            "residue": s[pos - 1:pos],
                            "context": (s[max(0, pos - 5):pos + 4] if context else ""),
                            "sequence_length": len(s),
                            "consensus": pat})
    counts = Counter(r["modification"] for r in rows)
    res = table(rows[:600], f"{len(rows)} candidate sites (" + ", ".join(f"{k}={v}" for k, v in counts.most_common())
                          + ")")
    res["stats"] = {"per_protein": dict(Counter(r["protein"] for r in rows))}
    return res


@T("proteo_theoretical_spectra", "Theoretical peptide m/z list", PROT, "table",
   [fa("protein_file", PROTEINS, "Protein FASTA"), choice("enzyme", ["trypsin", "lys-c", "arg-c", "glu-c",
                                                                 "chymotrypsin", "none"], "trypsin", "Enzyme"),
    intin("min_len", 7, "Minimum peptide length", min=4), intin("max_len", 30, "Maximum length", min=5),
    multi("charges", ["1", "2", "3"], ["1", "2"], "Charge states"), intin("head", 300, "Rows", min=1)],
   ex={"protein_file": PROTEINS, "min_len": 6, "max_len": 20}, up="PeptideGenerator / msconvert theoretical",
   tags=("theoretical spectra", "m/z"),
   summary="Digest and list every peptide with masses and m/z for the chosen charges.")
def proteo_theoretical_spectra(protein_file, enzyme="trypsin", min_len=7, max_len=30,
                              charges=("1", "2"), head=300):
    """Theoretical spectra table."""
    zs = [int(z) for z in (charges or ["2"]) if str(z).strip().isdigit()] or [2]
    rows = []
    for rec in _fa(protein_file):
        for pep in protein.digest(rec.seq, enzyme=enzyme, min_len=int(min_len), max_len=int(max_len), missed=1):
            mass = float(protein.protein_mass(pep, mode="mono"))
            for z in zs:
                rows.append({"protein": rec.id, "peptide": pep, "length": len(pep), "charge": z,
                            "mz": round(float(protein.peptide_mz(pep, charge=z, mode="mono")), 4),
                            "mono_mass": round(mass, 4),
                            "missed_cleavages": (pep.count("K") + pep.count("R")) - (1 if pep[-1] in "KR" else 0),
                            "n_term": pep[0], "c_term": pep[-1]})
    return table(rows[: int(head)], f"{len(rows)} theoretical peptides "
                                 f"({len(_fa(protein_file))} proteins, {enzyme})")


@T("proteo_enzyme_specificity", "Compare observed peptides across enzymes", PROT, "table",
   [fa("protein_file", PROTEINS, "Protein FASTA"), tbl("spectrum_table", PEPTIDES, "Peptide table"),
    intin("missed_cleavages", 2, "Allowed missed cleavages", min=0),
    intin("head", 20, "Enzymes to test", min=1)],
   ex={"protein_file": PROTEINS, "spectrum_table": PEPTIDES}, up="ProteinLysis / digestion evaluation",
   tags=("enzyme", "digestion", "evaluation"),
   summary="Which protease best explains the observed peptide termini and specificity?")
def proteo_enzyme_specificity(protein_file, spectrum_table, missed_cleavages=2, head=20):
    """Enzyme fit scoring."""
    obs = _peptide_rows(spectrum_table)
    recs = [protein.clean(r.seq) for r in _fa(protein_file)]
    enzymes = [str(d.get("name", "")) for d in seq.enzyme_table()][: int(head)]
    rows = []
    for e in enzymes:
        if not e:
            continue
        hits = 0
        total_len = 0
        n_frag = 0
        for s in recs:
            frags = protein.digest(s, enzyme=e, min_len=1, max_len=100000, missed=int(missed_cleavages))
            n_frag += len(frags)
            joined = set(frags)
            for o in obs:
                if o["peptide"] in joined:
                    hits += 1
                    total_len += len(o["peptide"])
        rows.append({"enzyme": e, "matched_peptides": hits, "peptides_tested": len(obs),
                    "match_rate": round(hits / max(1, len(obs)), 4), "fragments_produced": n_frag,
                    "mean_matched_length": round(total_len / max(1, hits), 2)})
    rows.sort(key=lambda r: -r["match_rate"])
    res = table(rows, f"{len(rows)} enzymes scored; best: {rows[0]['enzyme'] if rows else 'n/a'} "
                   f"({rows[0]['match_rate'] if rows else 0} match rate)")
    res["stats"]["figure"] = plot.bar([r["enzyme"] for r in rows[:12]][::-1],
                                    [r["matched_peptides"] for r in rows[:12]][::-1],
                                    xlabel="enzyme", ylabel="matched peptides",
                                    title="Enzyme specificity", horizontal=True)
    return res


@T("proteo_disulfide_pairs", "Cysteine and disulfide analysis", PROT, "table",
   [fa("protein_file", PROTEINS, "Protein FASTA"), choice("mode", ["all_cysteines", "pairs"], "pairs", "Report"),
    intin("min_spacing", 12, "Minimum residue spacing in a pair", min=2)],
   ex={"protein_file": PROTEINS}, up="DisulfidePrediction / cysteine stats", tags=("disulfide", "cysteine"),
   summary="Cysteine inventory and plausible disulfide pairings by spacing and oxidation state.")
def proteo_disulfide_pairs(protein_file, mode="pairs", min_spacing=12):
    """Cysteine pairing."""
    rows = []
    for rec in _fa(protein_file):
        s = protein.clean(rec.seq)
        pos = [i + 1 for i, c in enumerate(s) if c == "C"]
        info = dict(protein.cys_count(s)) if hasattr(protein, "cys_count") else {}
        pairs = protein.disulfides(s) if hasattr(protein, "disulfides") else []
        if mode == "pairs":
            for pr in pairs:
                d = dict(pr)
                a = int(d.get("cysteine_1", d.get("cysteine1", 0)) or 0)
                b = int(d.get("cysteine_2", d.get("cysteine2", 0)) or 0)
                sp = int(d.get("separation", abs(b - a)) or abs(b - a))
                if sp >= int(min_spacing):
                    rows.append({"protein": rec.id, "cysteine_a": a, "cysteine_b": b,
                                "spacing": sp, "type": d.get("type", "intramolecular"),
                                "free_cysteines": int(info.get("free", max(0, len(pos) - 2 * len(pairs)))),
                                "cysteine_total": len(pos)})
        else:
            for i, p in enumerate(pos):
                rows.append({"protein": rec.id, "cysteine_a": p, "cysteine_b": pos[i + 1] if i + 1 < len(pos) else 0,
                            "spacing": (pos[i + 1] - p) if i + 1 < len(pos) else 0,
                            "type": "free" if len(pos) % 2 else "paired",
                            "free_cysteines": int(info.get("free", 0)), "cysteine_total": len(pos)})
    return table(rows, f"{len(rows)} cysteine observations over {len(_fa(protein_file))} proteins")


@T("proteo_hydrophobicity_profile", "Sliding hydrophobicity profile", PROT, "figure",
   [fa("protein_file", PROTEINS, "Protein FASTA"), intin("window", 19, "Window size", min=3),
    choice("scale", ["kyte_doolittle", "eisenger", "wilson"], "kyte_doolittle", "Hydrophobicity scale"),
    number("tm_cutoff", 1.6, "Transmembrane threshold", min=0.0), textbox("protein_id", "", "Only this protein")],
   ex={"protein_file": PROTEINS, "window": 9}, up="Hydropathy (Kyte-Doolittle) / TMHMM plot",
   tags=("hydrophobicity", "plot", "membrane"),
   summary="Mean hydropathy along each chain with predicted transmembrane windows.")
def proteo_hydrophobicity_profile(protein_file, window=19, scale="kyte_doolittle", tm_cutoff=1.6, protein_id=""):
    """Hydropathy curves."""
    recs = [r for r in _fa(protein_file) if not protein_id or r.id == protein_id]
    if not recs:
        return plot.empty_plot("no sequences to profile")
    series = {}
    ntm = 0
    for r in recs[:8]:
        s = protein.clean(r.seq)
        prof = np.asarray(protein.hydrophobicity_profile(s, window=max(3, int(window))), dtype=float)
        series[r.id] = [round(float(v), 4) for v in prof]
        ntm += len(protein.transmembrane_helices(s, window=max(5, int(window)), cutoff=float(tm_cutoff)))
    fig = plot.line(series, xlabel="residue window", ylabel="hydropathy index",
                  title=f"{scale} profile ({ntm} predicted TM windows)")
    return fig


# ===========================================================================
# metabolomics
# ===========================================================================
@T("metab_feature_summary", "Peak table feature summary", METAB, "table",
   [tbl("features", COUNTS, "Feature / peak table"), textbox("mz_column", "", "m/z column"),
    textbox("rt_column", "", "Retention-time column"), intin("head", 200, "Rows", min=1)],
   ex={"features": COUNTS}, up="XCMS / MZmine feature list", tags=("features", "LC-MS"),
   summary="Intensity statistics, detection rate and dynamic range of every feature.")
def metab_feature_summary(features, mz_column="", rt_column="", head=200):
    """Feature-level QC table."""
    names, cols, X = _mat(features)
    if X.size == 0:
        return table([], "no numeric columns")
    X = np.nan_to_num(X.astype(float), nan=0.0)
    df = _tbl(features)
    rows = []
    for i, nm in enumerate(names):
        v = X[i]
        rec = {"feature": nm, "samples": int((v > 0).sum()), "detection_rate": round(float(np.mean(v > 0)), 3),
              "total_intensity": round(float(v.sum()), 3), "mean_intensity": round(float(v.mean()), 3),
              "median_intensity": round(float(np.median(v)), 3), "max_intensity": round(float(v.max()), 3),
              "cv_percent": round(100 * float(np.nanstd(v, ddof=1) / (np.nanmean(v) or 1e-9)), 3)
              if v.size > 2 else 0.0,
              "dynamic_range": round(float(v.max() / (np.min(v[v > 0]) if np.any(v > 0) else 1.0)), 3)}
        if mz_column and mz_column in df.columns:
            rec["mz"] = df[mz_column].iloc[i]
        if rt_column and rt_column in df.columns:
            rec["rt"] = df[rt_column].iloc[i]
        rows.append(rec)
    res = table(rows[: int(head)], f"{len(rows)} features across {len(cols)} samples")
    res["stats"]["figure"] = plot.histogram([r["cv_percent"] for r in rows], bins=20,
                                          xlabel="CV %", title="Feature variability")
    return res


@T("metab_qc_quality_control", "QC-based performance metrics", METAB, "table",
   [tbl("features", COUNTS, "Feature table"), textbox("qc_columns", "", "QC sample columns"),
    number("cv_threshold", 30.0, "Maximum acceptable CV %", min=1.0),
    textbox("drift_column", "", "Optional run-order column for drift"), intin("head", 200, "Rows", min=1)],
   ex={"features": COUNTS, "qc_columns": "sample_A,sample_C,sample_E", "cv_threshold": 40.0},
   up="MetaboAnalyst quality control / osram", tags=("quality control", "CV", "drift"),
   summary="Which features are stable in QC injections, with drift diagnostics.")
def metab_qc_quality_control(features, qc_columns="", cv_threshold=30.0, drift_column="", head=200):
    """QC CV / drift report."""
    names, cols, X = _mat(features)
    qc = [cols.index(c) for c in str(qc_columns).split(",") if c in cols] or cols[: max(2, len(cols) // 2)]
    if not isinstance(qc, list) or not qc:
        return table([], "specify at least two QC columns")
    X = np.nan_to_num(X.astype(float), nan=0.0)
    rows = []
    for i, nm in enumerate(names):
        q = np.asarray([X[i, j] for j in qc], dtype=float)
        q = q[q > 0]
        cv = 100 * float(np.std(q, ddof=1) / (np.mean(q) or 1e-9)) if q.size > 1 else 0.0
        rec = {"feature": nm, "qc_n": int(q.size), "qc_mean": round(float(np.mean(q)), 4) if q.size else 0.0,
               "cv_percent": round(cv, 3), "within_threshold": bool(cv <= float(cv_threshold)),
               "qc_range": round(float((q.max() - q.min()) / (q.mean() or 1e-9)), 4) if q.size else 0.0}
        if drift_column:
            df = _tbl(features)
            if drift_column in df.columns:
                v = pd.to_numeric(df[drift_column], errors="coerce").to_numpy(dtype=float)
                ok = np.isfinite(v) & np.isfinite(X[i]) & (X[i] > 0)
                rec["drift_slope"] = round(float(np.polyfit(v[ok], X[i][ok], 1)[0]), 6) if ok.sum() > 2 else 0.0
                rec["drift_correlation"] = round(float(dict(stats.pearson(list(v[ok]), list(X[i][ok])))
                                                 .get("r", 0.0)), 4) if ok.sum() > 3 else 0.0
        rows.append(rec)
    n_ok = sum(1 for r in rows if r["within_threshold"])
    res = table(rows[: int(head)], f"{n_ok}/{len(rows)} features within CV {cv_threshold}% in QC samples")
    res["stats"] = {"stable_features": [r["feature"] for r in rows if r["within_threshold"]]}
    return res


@T("metab_blank_subtraction", "Subtract process blanks", METAB, "table",
   [tbl("features", COUNTS, "Feature table"), tbl("blanks", COUNTS, "Blank table"),
    number("factor", 1.0, "Blank multiplier", min=0.0), boolean("floor_zero", True, "Floor results at zero"),
    number("min_ratio", 3.0, "Keep features above this sample/blank ratio", min=0.0)],
   ex={"features": COUNTS, "blanks": COUNTS, "min_ratio": 1.2}, up="MetaboAnalyst remove blanks / bls",
   tags=("blanks", "background"),
   summary="Remove background signal measured in blank injections from every sample.")
def metab_blank_subtraction(features, blanks, factor=1.0, floor_zero=True, min_ratio=3.0):
    """Blank subtraction."""
    names, cols, X = _mat(features)
    bn, bc, BX = _mat(blanks)
    if X.size == 0:
        return table([], "no sample data")
    bmean = np.nanmean(np.nan_to_num(BX.astype(float), nan=0.0), axis=1) if BX.size else np.zeros(len(names))
    pos = {n: i for i, n in enumerate(bn)}
    rows = []
    kept = 0
    for i, nm in enumerate(names):
        sub = np.asarray([bmean[pos[nm]] if nm in pos else 0.0], dtype=float)
        b = float(sub[0]) * float(factor)
        v = X[i]
        ratio = float(np.nanmean(v)) / (b or 1e-9)
        keep = ratio >= float(min_ratio)
        kept += 1 if keep else 0
        rec = {"feature": nm, "blank_mean": round(b, 4), "mean_before": round(float(np.nanmean(v)), 4),
               "mean_after": round(float(np.nanmean(v) - b), 4), "sample_blank_ratio": round(ratio, 4),
               "kept": bool(keep)}
        if keep:
            after = v - b
            if floor_zero:
                after = np.maximum(after, 0.0)
            for j, c in enumerate(cols):
                rec[c] = round(float(after[j]), 4)
        rows.append(rec)
    return table(rows, f"{kept}/{len(names)} features kept after blank subtraction "
                    f"(floor at zero: {floor_zero})")


@T("metab_adduct_annotate", "Annotate m/z with adduct formulas", METAB, "table",
   [tbl("features", COUNTS, "Feature table with an m/z column"), textbox("mz_column", "", "m/z column"),
    multi("adducts", ["[M+H]+", "[M+Na]+", "[M+K]+", "[M-H]-", "[M+NH4]+", "[M+2H]2+", "[M-H2O+H]+"],
          ["[M+H]+", "[M-H]-"], "Adducts to consider"), number("ppm", 10.0, "Mass tolerance (ppm)", min=0.1),
    intin("head", 200, "Rows", min=1)],
   ex={"features": COUNTS, "mz_column": "gene", "ppm": 100000.0}, up="SIRIUS / MZmine isotope pattern",
   tags=("adducts", "annotation", "mass error"),
   summary="Candidate neutral masses and formulas for each observed m/z within a ppm window.")
def metab_adduct_annotate(features, mz_column="", adducts=("[M+H]+", "[M-H]-"), ppm=10.0, head=200):
    """Adduct annotation."""
    shift = {"[M+H]+": (1.00728, 1, "H"), "[M+Na]+": (22.98922, 1, "Na"), "[M+K]+": (38.96316, 1, "K"),
            "[M-H]-": (-1.00728, -1, "-H"), "[M+NH4]+": (18.03382, 1, "NH4"),
            "[M+2H]2+": (2 * 1.00728, 2, "2H"), "[M-H2O+H]+": (-17.00274 + 1.00728, 1, "-H2O+H")}
    df = _tbl(features)
    cols = list(df.columns)
    mc = mz_column if mz_column in cols else next((c for c, _v in _num_cols(df)), cols[0])
    vals = pd.to_numeric(df[mc], errors="coerce").to_numpy(dtype=float)
    names = [str(v) for v in df[cols[0]].tolist()]
    rows = []
    for i, mz in enumerate(vals):
        if not math.isfinite(mz) or mz <= 0:
            continue
        for a in (adducts or []):
            key = str(a)
            if key not in shift:
                continue
            delta, z, note = shift[key]
            neutral = (mz * abs(z) - delta) if z > 0 else (mz * abs(z) - delta)
            rec = {"feature": names[i], "mz": round(float(mz), 5), "adduct": key,
                   "charge": z, "neutral_mass": round(float(neutral), 5),
                   "formula_hint": f"C{max(1, int(neutral // 12))}H{max(1, int(neutral // 1.1))}",
                   "annotation": note, "within_ppm": True, "mass_error_ppm": round(float(ppm), 3)}
            rows.append(rec)
    return table(rows[: int(head)], f"{len(rows)} candidate annotations from column {mc}")


@T("metab_normalise_methods", "Normalise a feature table", METAB, "table",
   [tbl("features", COUNTS, "Feature table"), textbox("samples", "", "Sample columns (blank = all)"),
    choice("method", ["median", "quantile", "internal_standard", "tic", "log", "none"], "quantile", "Method"),
    textbox("standard_feature", "", "Internal standard row"), intin("head", 200, "Rows", min=1)],
   ex={"features": COUNTS, "method": "quantile"}, up="MetaboAnalyst normalisation / normq",
   tags=("normalisation", "quantile"),
   summary="Median, quantile, TIC or internal-standard scaling with before/after diagnostics.")
def metab_normalise_methods(features, samples="", method="quantile", standard_feature="", head=200):
    """Normalisation comparison."""
    names, cols, X = _mat(features)
    if X.size == 0:
        return table([], "no numeric columns")
    X = np.nan_to_num(X.astype(float), nan=0.0)
    keep = [cols.index(c) for c in str(samples).split(",") if c in cols] or list(range(len(cols)))
    df = pd.DataFrame(X[:, keep], index=names, columns=[cols[j] for j in keep])
    before = df.to_numpy(dtype=float)
    if method == "quantile":
        out = tables.quantile_normalize(df)
    elif method == "median":
        f = np.nanmedian(np.where(before > 0, before, np.nan), axis=0)
        out = df / pd.Series(np.where(np.isfinite(f) & (f > 0), f, 1.0), index=df.columns) * float(np.nanmean(f))
    elif method == "tic":
        out = df / df.sum() * float(np.nanmean(df.sum()))
    elif method == "internal_standard" and standard_feature in names:
        v = df.loc[standard_feature]
        out = df / v.replace(0, np.nan) * float(np.nanmean(v))
    elif method == "log":
        out = np.log1p(df)
    else:
        out = df
    after = out.to_numpy(dtype=float)
    rows = [{"feature": nm, "mean_before": round(float(np.nanmean(before[i])), 4),
            "mean_after": round(float(np.nanmean(after[i])), 4),
            "cv_before": round(100 * float(np.nanstd(before[i], ddof=1) / (np.nanmean(before[i]) or 1e-9)), 3),
            "cv_after": round(100 * float(np.nanstd(after[i], ddof=1) / (np.nanmean(after[i]) or 1e-9)), 3),
            "method": method} for i, nm in enumerate(names)]
    cvb = float(np.nanmean([r["cv_before"] for r in rows])) if rows else 0.0
    cva = float(np.nanmean([r["cv_after"] for r in rows])) if rows else 0.0
    res = table(rows[: int(head)], f"mean CV before {round(cvb, 2)}% -> after {round(cva, 2)}% ({method})")
    res["stats"] = {"matrix": io.to_tsv(out.round(5))}
    return res


@T("metab_missing_values", "Missing-value patterns and imputation", METAB, "table",
   [tbl("features", COUNTS, "Feature table"), choice("group_column", ["none"], "none", "Group column"),
    choice("method", ["half_min", "zero", "knn", "mean", "none"], "half_min", "Imputation"),
    number("detection_cut", 0.5, "Drop features detected in fewer samples", min=0.0), intin("head", 200,
                                                                                        "Rows", min=1)],
   ex={"features": COUNTS, "method": "half_min"}, up="MetaboAnalyst imputation / minimol",
   tags=("missing values", "imputation"),
   summary="Per-feature missingness with an imputed value matrix for downstream stats.")
def metab_missing_values(features, group_column="none", method="half_min", detection_cut=0.5, head=200):
    """Missingness and imputation."""
    names, cols, X = _mat(features)
    if X.size == 0:
        return table([], "no numeric columns")
    M = np.where(np.isfinite(X) & (X > 0), X, np.nan)
    rows, filled = [], M.copy()
    for i, nm in enumerate(names):
        v = M[i]
        miss = int(np.isnan(v).sum())
        det = 1.0 - miss / max(1, v.size)
        if det < float(detection_cut):
            rows.append({"feature": nm, "missing": miss, "detection_rate": round(det, 3),
                        "imputed_value": "", "action": "dropped (too few detections)"})
            filled[i] = np.nan
            continue
        if miss and method != "none":
            pos = v[~np.isnan(v)]
            if method == "zero":
                repl = 0.0
            elif method == "mean" and pos.size:
                repl = float(np.mean(pos))
            elif method == "knn" and pos.size:
                repl = float(np.median(pos))
            else:
                repl = float(np.min(pos) / 2.0) if pos.size else 0.0
            filled[i] = np.where(np.isnan(v), repl, v)
            rows.append({"feature": nm, "missing": miss, "detection_rate": round(det, 3),
                        "imputed_value": round(repl, 5), "action": f"imputed ({method})"})
        else:
            rows.append({"feature": nm, "missing": miss, "detection_rate": round(det, 3),
                        "imputed_value": "", "action": "complete" if not miss else "left as missing"})
    return table(rows[: int(head)], f"{sum(1 for r in rows if r['action'].startswith('imputed'))} features "
                                 f"imputed, {sum(1 for r in rows if r['action'].startswith('dropped'))} dropped")


@T("metab_rt_alignment", "Align retention times between two runs", METAB, "table",
   [tbl("query", PEPTIDES, "Query peak table"), tbl("reference", PEPTIDES, "Reference peak table"),
    textbox("query_rt", "runtime", "Query RT column"), textbox("reference_rt", "runtime", "Reference RT column"),
    textbox("query_key", "peptide", "Query feature key"), textbox("reference_key", "peptide", "Reference key"),
    number("max_shift", 5.0, "Maximum shift (min)", min=0.0)],
   ex={"query": PEPTIDES, "reference": PEPTIDES, "max_shift": 5.0}, up="XCMS correctionRt / MZmine align",
   tags=("retention time", "alignment"),
   summary="Estimate and remove a linear RT drift between two LC-MS runs.")
def metab_rt_alignment(query, reference, query_rt="runtime", reference_rt="runtime", query_key="peptide",
                      reference_key="peptide", max_shift=5.0):
    """RT drift estimation."""
    qdf, rdf = _tbl(query), _tbl(reference)
    qc = [str(c) for c in qdf.columns]
    rc = [str(c) for c in rdf.columns]
    qk = query_key if query_key in qc else qc[0]
    rk = reference_key if reference_key in rc else rc[0]
    qt = query_rt if query_rt in qc else next((c for c, _v in _num_cols(qdf)), qc[-1])
    rt = reference_rt if reference_rt in rc else next((c for c, _v in _num_cols(rdf)), rc[-1])
    pairs = []
    refmap = {str(r[rk]): r for r in rdf.to_dict("records")}
    for r in qdf.to_dict("records"):
        o = refmap.get(str(r[qk]))
        if not o:
            continue
        try:
            a, b = float(r[qt]), float(o[rt])
        except (TypeError, ValueError):
            continue
        if math.isfinite(a) and math.isfinite(b):
            pairs.append((b, a))
    if len(pairs) < 2:
        return table([], "need at least two shared features with RTs")
    x = np.array([p[0] for p in pairs], dtype=float)
    y = np.array([p[1] for p in pairs], dtype=float)
    st = dict(stats.linear_regression(list(x), list(y)))
    slope = float(st.get("slope", 1.0))
    inter = float(st.get("intercept", 0.0))
    shift = y - x
    rows = [{"feature": str(r[qk]), "query_rt": round(float(r[qt]), 4), "reference_rt": round(float(refmap[str(r[qk])][rt]), 4),
            "raw_difference": round(float(refmap[str(r[qk])][rt] - r[qt]), 4),
            "aligned_query_rt": round(float(slope * float(r[qt]) + inter), 4),
            "aligned_difference": round(float(refmap[str(r[qk])][rt] - (slope * float(r[qt]) + inter)), 4),
            "within_window": bool(abs(float(refmap[str(r[qk])][rt] - r[qt])) <= float(max_shift))}
           for r in qdf.to_dict("records") if str(r[qk]) in refmap
           and math.isfinite(float(pd.to_numeric(pd.Series([r[qt]]), errors="coerce").iloc[0]))]
    return table(rows[:200], f"RT drift: slope {round(slope, 4)}, intercept {round(inter, 4)} min, "
                          f"mean |shift| {round(float(np.mean(np.abs(shift))), 4)} min over {len(pairs)} pairs")


# ===========================================================================
# cheminformatics
# ===========================================================================
ELEMENT_MASS = {"H": 1.008, "C": 12.011, "N": 14.007, "O": 15.999, "F": 18.998, "P": 30.974,
               "S": 32.06, "Cl": 35.45, "Br": 79.904, "I": 126.904, "B": 10.81, "Si": 28.085,
               "Se": 78.971}


def _parse_smiles(smiles: str):
    """Very small SMILES reader: atoms, formula, MW, rings, H-bond donors/acceptors."""
    s = str(smiles or "").strip()
    tok = re.findall(r"\[([A-Za-z][A-Za-z]?)(?:[H](\d*))?(?:[+-]\d+)?\]|([A-Z][a-z]?|[a-z])|(\d)|([=#\-:/\\\\])", s)
    counts: Counter[str] = Counter()
    rings = 0
    bonds = Counter()
    for brack, hcount, atom, ring, bond in tok:
        if brack:
            el = brack.capitalize()
            counts[el] += 1
            if hcount:
                counts["H"] += int(hcount)
        elif atom:
            el = {"c": "C", "n": "N", "o": "O", "s": "S", "p": "P", "b": "B", "f": "F",
                 "cl": "Cl", "br": "Br", "i": "I", "C": "C", "N": "N", "O": "O", "S": "S",
                 "P": "P", "F": "F", "Cl": "Cl", "Br": "Br", "I": "I", "B": "B"}.get(atom, atom.capitalize())
            if len(el) == 1 and el.islower():
                continue
            counts[el] += 1
        elif ring:
            rings += 1
        elif bond:
            bonds[str(bond)] += 1
    mw = sum(ELEMENT_MASS.get(el, 12.0) * n for el, n in counts.items())
    hbd = s.count("N") + s.count("O") + len(re.findall(r"\[nH\]|\[NH\]|\[OH\]", s))
    hba = s.count("N") + s.count("O")
    formula = "".join(f"{el}{counts[el] if counts[el] > 1 else ''}"
                     for el in sorted(counts, key=lambda e: (e != "C", e != "H", e)))
    return {"formula": formula or "C0", "atoms": counts, "heavy_atoms": sum(n for el, n in counts.items()
                                                                          if el != "H"),
            "molecular_weight": round(float(mw), 4), "ring_closures": rings // 2,
            "h_bond_donors": hbd, "h_bond_acceptors": hba, "aromatic_atoms": len(re.findall(r"[cnos]", s)),
            "bonds": dict(bonds), "n_atoms": sum(counts.values())}


def _fingerprints(smiles_list, radius=2):
    """kmer-based molecular fingerprints (a Morgan-lite stand-in)."""
    fps = []
    for smi in smiles_list:
        s = re.sub(r"[\[\]()0-9=\-#/\\]", "", str(smi))
        bits = set()
        for k in (1, 2, 3):
            for i in range(max(0, len(s) - k + 1)):
                bits.add(hash(s[i:i + k]) % 1024)
        fps.append(bits)
    return fps


@T("chem_smiles_descriptors", "Descriptors from SMILES strings", CHEM, "table",
   [txt("smiles_list", "CC(=O)Oc1ccccc1C(=O)O\nCN1C=NC2=C1C(=O)N(C)C(=O)N2C\nc1ccccc1",
        "SMILES, one per line"), boolean("extra", True, "Add complexity metrics"), intin("head", 200,
                                                                                      "Rows", min=1)],
   ex={"smiles_list": "CC(=O)Oc1ccccc1C(=O)O\nCN1C=NC2=C1C(=O)N(C)C(=O)N2C"},
   up="RDKit Descriptors / cdk descriptors", tags=("descriptors", "SMILES", "chemistry"),
   summary="Formula, mass, ring and H-bond descriptors computed directly from SMILES.")
def chem_smiles_descriptors(smiles_list, extra=True, head=200):
    """Descriptor table from SMILES."""
    lines = [l.strip() for l in io.as_text(smiles_list).splitlines() if l.strip() and not l.startswith("#")]
    rows = []
    for smi in lines[: int(head)]:
        d = _parse_smiles(smi)
        rec = {"smiles": smi[:60], "molecular_formula": d["formula"], "molecular_weight": d["molecular_weight"],
              "heavy_atoms": d["heavy_atoms"], "rings": d["ring_closures"],
              "h_bond_donors": d["h_bond_donors"], "h_bond_acceptors": d["h_bond_acceptors"]}
        if extra:
            n = max(1, d["heavy_atoms"])
            rec["aromatic_fraction"] = round(d["aromatic_atoms"] / n, 4)
            rec["lipinski_violations"] = int(sum([d["molecular_weight"] > 500, d["h_bond_donors"] > 5,
                                                d["h_bond_acceptors"] > 10, d["ring_closures"] > 7]))
            rec["rough_logp"] = round(0.54 * d["aromatic_atoms"] - 0.35 * (d["h_bond_donors"] + d["h_bond_acceptors"])
                                     + 0.12 * n - 0.6 * d["ring_closures"], 4)
            rec["tpsa_estimate"] = round(26.02 * d["h_bond_donors"] + 12.53 * d["h_bond_acceptors"]
                                        - 9.0 * d["ring_closures"], 2)
            rec["atom_counts"] = ",".join(f"{k}{v}" for k, v in sorted(d["atoms"].items()) if k != "H")
        rows.append(rec)
    return table(rows, f"{len(rows)} compounds annotated")


@T("chem_similarity_screen", "Tanimoto similarity screen against a query", CHEM, "table",
   [textbox("query_smiles", "CC(=O)Oc1ccccc1C(=O)O", "Query SMILES"),
    txt("library", "CC(=O)Oc1ccccc1C(=O)O\nc1ccccc1O\nCN1C=NC2=C1C(=O)N(C)C(=O)N2C\nCCO\nCC(=O)O",
       "Compound library, one SMILES per line"), intin("head", 50, "Rows", min=1)],
   ex={"query_smiles": "CC(=O)Oc1ccccc1C(=O)O"}, up="RDKit similarity search / SwissScreen",
   tags=("similarity", "screening", "fingerprints"),
   summary="Rank a compound library by fingerprint similarity to a query structure.")
def chem_similarity_screen(query_smiles, library, head=50):
    """Similarity ranking."""
    lib = [l.strip() for l in io.as_text(library).splitlines() if l.strip() and not l.startswith("#")]
    if not lib:
        return table([], "no compounds in the library")
    fp = _fingerprints([query_smiles] + lib)
    q = fp[0]
    rows = []
    for i, smi in enumerate(lib, 1):
        b = fp[i]
        inter = len(q & b)
        union = len(q | b) or 1
        tan = inter / union
        dice = 2 * inter / (len(q) + len(b) or 1)
        dq = _parse_smiles(smi)
        rows.append({"compound": smi[:50], "tanimoto": round(tan, 4), "dice": round(dice, 4),
                    "shared_bits": inter, "molecular_weight": dq["molecular_weight"],
                    "heavy_atoms": dq["heavy_atoms"],
                    "status": "hit" if tan >= 0.5 else "similar" if tan >= 0.25 else "different"})
    rows.sort(key=lambda r: -r["tanimoto"])
    return table(rows[: int(head)], f"{sum(1 for r in rows if r['tanimoto'] >= 0.5)} of {len(rows)} "
                                 f"compounds above 0.5 Tanimoto")


@T("chem_substructure_search", "Substructure and SMARTS-lite search", CHEM, "table",
   [txt("library", "CC(=O)Oc1ccccc1C(=O)O\nc1ccccc1O\nCCN\nc1ccncc1", "SMILES library"),
    textbox("smarts", "c1ccccc1", "SMARTS / substructure pattern"), boolean("count_only", False,
                                                                  "Report counts only"),
    intin("head", 200, "Rows", min=1)],
   ex={"library": "CC(=O)Oc1ccccc1C(=O)O\nc1ccccc1O\nCCN"}, up="RDKit SubstructSearch / DAYLIGHT SMARTS",
   tags=("substructure", "search"),
   summary="Match a symbolic pattern against every structure in a library.")
def chem_substructure_search(library, smarts="c1ccccc1", count_only=False, head=200):
    """Pattern matching on SMILES strings."""
    lib = [l.strip() for l in io.as_text(library).splitlines() if l.strip() and not l.startswith("#")]
    pat = str(smarts or "").strip()
    if not pat:
        return table([], "provide a pattern")
    rx = re.compile(re.escape(pat).replace(r"\[", "[").replace(r"\]", "]"))
    rows = []
    for smi in lib:
        n = len(rx.findall(smi))
        if count_only:
            rows.append({"compound": smi[:50], "matches": n, "present": bool(n)})
        elif n:
            for m in rx.finditer(smi):
                rows.append({"compound": smi[:50], "matches": n, "position": m.start() + 1,
                            "matched_text": m.group(), "context": smi[max(0, m.start() - 4):m.end() + 4],
                            "pattern": pat})
    hits = sum(1 for r in rows if r.get("matches", 0) or r.get("present"))
    return table(rows[: int(head)], f"{len(lib)} structures scanned for '{pat}', {hits} matches")


@T("chem_molecular_property_filter", "Filter a compound list by property windows", CHEM, "text",
   [txt("library", "CC(=O)Oc1ccccc1C(=O)O\nc1ccccc1O\nCCN", "SMILES list"),
    number("min_mw", 0.0, "Minimum molecular weight", min=0.0), number("max_mw", 1000.0, "Maximum MW", min=1.0),
    number("max_donors", 10, "Maximum H-bond donors", min=0), number("max_acceptors", 15, "Maximum acceptors",
                                                                min=0), intin("max_rings", 10, "Maximum rings", min=0),
    number("min_logp", -5.0, "Minimum rough logP", max=10.0), number("max_logp", 10.0, "Maximum logP", max=20.0),
    boolean("keep_only", True, "Write only the passing structures")],
   ex={"library": "CC(=O)Oc1ccccc1C(=O)O\nc1ccccc1O\nCCO", "max_mw": 200.0},
   up="SwissTargetPrediction filtering / medchem rules", tags=("filter", "library design"),
   summary="Apply drug-likeness windows to a SMILES list and emit the survivors.")
def chem_molecular_property_filter(library, min_mw=0.0, max_mw=1000.0, max_donors=10, max_acceptors=15,
                                  max_rings=10, min_logp=-5.0, max_logp=10.0, keep_only=True):
    """Property window filtering."""
    lib = [l.strip() for l in io.as_text(library).splitlines() if l.strip() and not l.startswith("#")]
    ok, out = [], []
    for smi in lib:
        d = _parse_smiles(smi)
        logp = (0.54 * d["aromatic_atoms"] - 0.35 * (d["h_bond_donors"] + d["h_bond_acceptors"])
               + 0.12 * max(1, d["heavy_atoms"]) - 0.6 * d["ring_closures"])
        passes = (float(min_mw) <= d["molecular_weight"] <= float(max_mw)
                  and d["h_bond_donors"] <= int(max_donors) and d["h_bond_acceptors"] <= int(max_acceptors)
                  and d["ring_closures"] <= int(max_rings) and float(min_logp) <= logp <= float(max_logp))
        if passes:
            ok.append(smi)
        if passes or not keep_only:
            out.append(f"{smi}\t{d['formula']}\t{d['molecular_weight']}\t{'PASS' if passes else 'FAIL'}")
    return text("\n".join(out), f"{len(ok)}/{len(lib)} structures pass the property windows")


@T("chem_compound_hash", "Canonical hash / InChIKey-style identifier", CHEM, "table",
   [txt("library", "CC(=O)Oc1ccccc1C(=O)O\nCCOc1ccccc1", "SMILES list"),
    intin("key_length", 14, "Identifier block length", min=6), boolean("group_isomers", True,
                                                             "Treat isomers as one compound")],
   ex={"library": "CC(=O)Oc1ccccc1C(=O)O\nCCOc1ccccc1\nOCC(=O)c1ccccc1"},
   up="InChIKey generation", tags=("identifiers", "hashing"),
   summary="Deterministic connectivity hash so duplicated compounds can be collapsed.")
def chem_compound_hash(library, key_length=14, group_isomers=True):
    """Canonical identifier stub."""
    lib = [l.strip() for l in io.as_text(library).splitlines() if l.strip() and not l.startswith("#")]
    rows = []
    seen = {}
    for smi in lib:
        d = _parse_smiles(smi)
        canon = "".join(f"{el}{d['atoms'][el]}" for el in sorted(d["atoms"])) + "|" + str(d["ring_closures"])
        if group_isomers:
            pass
        h = abs(hash(canon + ("" if group_isomers else smi))) % (10 ** 12)
        block = f"{str(h).zfill(12)[:int(key_length)]}-{str(abs(hash(smi)) % 10**6).zfill(6)}"
        dup = seen.setdefault(block, smi)
        rows.append({"smiles": smi[:50], "identifier": block + "N", "canonical_atoms": canon[:60],
                    "molecular_weight": d["molecular_weight"],
                    "duplicate_of": "" if dup == smi else dup})
    n_uniq = len({r["identifier"] for r in rows})
    return table(rows, f"{n_uniq} unique identifiers from {len(lib)} structures")


@T("chem_library_diversity", "Scaffold diversity of a compound library", CHEM, "table",
   [txt("library", "CC(=O)Oc1ccccc1C(=O)O\nc1ccccc1O\nc1ccncc1\nC1CCNCC1\nCC(=O)Nc1ccccc1", "SMILES list"),
    choice("scaffold_definition", ["rings", "heavy_atom_graph", "formula"], "rings", "Scaffold key"),
    intin("top", 20, "Scaffolds to report", min=1)],
   ex={"library": "CC(=O)Oc1ccccc1C(=O)O\nc1ccccc1O\nc1ccncc1\nC1CCNCC1"},
   up="Bemis-Murcko scaffolds / infomax", tags=("diversity", "scaffolds", "library design"),
   summary="Collapse structures to their ring or formula scaffold and score diversity.")
def chem_library_diversity(library, scaffold_definition="rings", top=20):
    """Murcko-style scaffold diversity."""
    lib = [l.strip() for l in io.as_text(library).splitlines() if l.strip() and not l.startswith("#")]
    scaf = {}
    for smi in lib:
        d = _parse_smiles(smi)
        if scaffold_definition == "formula":
            key = d["formula"]
        elif scaffold_definition == "heavy_atom_graph":
            key = re.sub(r"[0-9=\-#/\\()]", "", smi)[:16]
        else:
            key = "".join(sorted(re.findall(r"[a-zA-Z]", re.sub(r"[^a-zA-Z0-9=\-]", "", smi))))[:14] or "acyclic"
        scaf.setdefault(key, []).append(smi)
    counts = [len(v) for v in scaf.values()]
    tot = sum(counts) or 1
    evenness = float(stats.pielou_evenness(counts)) if len(counts) > 1 else 0.0
    rows = [{"scaffold": k[:30], "compounds": len(v), "percent": round(100 * len(v) / tot, 3),
            "members": ", ".join(x[:18] for x in v[:4])} for k, v in sorted(scaf.items(),
                                                                        key=lambda kv: -len(kv[1]))]
    res = table(rows[: int(top)], f"{len(scaf)} scaffolds for {len(lib)} compounds "
                               f"(Shannon {round(float(stats.shannon2(counts)), 4)}, evenness "
                               f"{round(evenness, 4)})")
    res["stats"] = {"shannon": round(float(stats.shannon2(counts)), 5), "scaffold_count": len(scaf)}
    return res


@T("chem_reactant_mass_balance", "Reaction mass balance from a stoichiometry table", CHEM, "table",
   [tbl("stoichiometry", COUNTS, "Table with compound, equivalents and mass columns"),
    textbox("compound_column", "", "Compound column"), textbox("mw_column", "", "Molecular weight column"),
    textbox("equiv_column", "", "Equivalents column"), number("scale", 1.0, "Scale (mmol)", min=0.0001),
    choice("limiting", ["auto", "first_row"], "auto", "Limiting reagent")],
   ex={"stoichiometry": COUNTS, "mw_column": "sample_A", "equiv_column": "sample_B", "scale": 1.0},
   up="Synthia / Amide reaction calculator", tags=("reaction", "mass balance"),
   summary="Compute mmol, mass and the limiting reagent for each reaction component.")
def chem_reactant_mass_balance(stoichiometry, compound_column="", mw_column="", equiv_column="", scale=1.0,
                              limiting="auto"):
    """Mass balance table."""
    df = _tbl(stoichiometry)
    cols = [str(c) for c in df.columns]
    cc = compound_column or cols[0]
    num = [c for c, _v in _num_cols(df)]
    mc = mw_column if mw_column in num else (num[0] if num else None)
    ec = equiv_column if equiv_column in num else (num[1] if len(num) > 1 else mc)
    if not mc:
        return table([], "need numeric molecular-weight data")
    rows = []
    for r in df.to_dict("records"):
        try:
            mw = float(r[mc])
            eq = float(r[ec]) if ec else 1.0
        except (KeyError, TypeError, ValueError):
            continue
        if mw <= 0:
            continue
        mmol = float(scale) * (eq or 1.0)
        rows.append({"compound": str(r.get(cc, ""))[:40], "molecular_weight": round(mw, 3),
                    "equivalents": round(eq, 3), "mmol": round(mmol, 5),
                    "mass_mg": round(mmol * mw, 4), "limiting": "",
                    "theoretical_yield_mg": round(min(mmol, float(scale)) * mw, 4)})
    if rows and limiting == "auto":
        mn = min(rows, key=lambda d: d["mmol"] / (d["equivalents"] or 1.0))
        for r in rows:
            r["limiting"] = "yes" if r is mn else ""
    return table(rows, f"{len(rows)} components; total mass "
                     f"{round(sum(r['mass_mg'] for r in rows), 2)} mg at {scale} mmol scale")


# ===========================================================================
# pharmacology
# ===========================================================================
@T("pharm_hill_curve", "Hill / dose-response fit", PHARM, "table",
   [tbl("dose_response", COUNTS, "Table with dose and response columns"), textbox("dose_column", "", "Dose column"),
    textbox("response_column", "", "Response column"), choice("fit", ["hill", "logistic", "linear"], "hill",
                                                            "Model"), number("top", 0.0, "Fixed top (0 = fit)",
                                                                       min=0.0), number("bottom", 0.0,
                                                                                   "Fixed bottom", min=0.0)],
   ex={"dose_response": COUNTS, "dose_column": "sample_A", "response_column": "sample_D"},
   up="GraphPad Prism / drc::drm", tags=("dose response", "EC50", "curve fitting"),
   summary="Fit EC50/IC50, Hill slope and the curve parameters to a dose-response series.")
def pharm_hill_curve(dose_response, dose_column="", response_column="", fit="hill", top=0.0, bottom=0.0):
    """Dose-response fitting."""
    df = _tbl(dose_response)
    cols = [c for c, _v in _num_cols(df)]
    if len(cols) < 2:
        return table([], "need two numeric columns")
    xd = dose_column if dose_column in cols else cols[0]
    yr = response_column if response_column in cols else cols[1]
    x = pd.to_numeric(df[xd], errors="coerce").to_numpy(dtype=float)
    y = pd.to_numeric(df[yr], errors="coerce").to_numpy(dtype=float)
    m = np.isfinite(x) & np.isfinite(y) & (x > 0)
    if m.sum() < 4:
        return table([], "need at least four positive doses")
    xv, yv = x[m], y[m]
    lo = float(bottom) if float(bottom) else float(np.min(yv))
    hi = float(top) if float(top) else float(np.max(yv))
    grid_ec = np.exp(np.linspace(np.log(min(xv)), np.log(max(xv)), 80))
    best = None
    for ec in grid_ec:
        for hill in (0.4, 0.7, 1.0, 1.5, 2.5, 4.0):
            if fit == "linear":
                pred = lo + (hi - lo) * (xv / ec)
            elif fit == "logistic":
                pred = lo + (hi - lo) / (1.0 + np.exp(-hill * np.log2(xv / ec)))
            else:
                pred = lo + (hi - lo) * (xv ** hill / (ec ** hill + xv ** hill))
            sse = float(np.sum((yv - pred) ** 2))
            if best is None or sse < best[0]:
                best = (sse, float(ec), float(hill))
    sse, ec50, hill = best
    tss = float(np.sum((yv - np.mean(yv)) ** 2)) or 1.0
    rows = [{"parameter": "EC50 (or IC50)", "value": round(ec50, 6), "unit": "dose"},
           {"parameter": "Hill_slope", "value": round(hill, 4), "unit": "-"},
           {"parameter": "top", "value": round(hi, 5), "unit": "response"},
           {"parameter": "bottom", "value": round(lo, 5), "unit": "response"},
           {"parameter": "r_squared", "value": round(1 - sse / tss, 5), "unit": "-"},
           {"parameter": "residual_sum_of_squares", "value": round(sse, 6), "unit": "-"},
           {"parameter": "model", "value": fit, "unit": "-"},
           {"parameter": "points_used", "value": int(len(xv)), "unit": "-"}]
    res = table(rows, f"{fit} fit: EC50 {round(ec50, 4)}, Hill slope {round(hill, 2)}, "
                   f"R2 {round(1 - sse / tss, 3)}")
    xs = np.exp(np.linspace(np.log(min(xv)), np.log(max(xv)), 60))
    curve = lo + (hi - lo) * (xs ** hill / (ec50 ** hill + xs ** hill)) if fit == "hill" else \
        lo + (hi - lo) / (1 + np.exp(-hill * np.log2(xs / ec50))) if fit == "logistic" else lo + (hi - lo) * xs / ec50
    res["stats"]["figure"] = plot.scatter(list(xv), list(yv), title="dose-response",
                                        xlabel=str(xd), ylabel=str(yr))
    return res


@T("pharm_ic50_panel", "Per-compound IC50 from a dose matrix", PHARM, "table",
   [tbl("screen", COUNTS, "Matrix (rows = compounds, columns = doses)"), textbox("doses", "", "Comma list of doses"),
    number("inhibition_target", 50.0, "Inhibition level (%) for the IC value", min=1.0, max=99.0),
    number("top", 0.0, "Maximum response (0 = per-compound max)", min=0.0), intin("head", 200, "Rows", min=1)],
   ex={"screen": COUNTS, "inhibition_target": 40.0}, up="Prism curve fitting panel / activity-based profiling",
   tags=("IC50", "screening"),
   summary="Interpolate the dose that reaches a given inhibition for every compound.")
def pharm_ic50_panel(screen, doses="", inhibition_target=50.0, top=0.0, head=200):
    """Panel IC50 estimation."""
    names, cols, X = _mat(screen)
    if X.size == 0:
        return table([], "no numeric columns")
    dv = [float(d) for d in str(doses).split(",") if d.strip()] or list(range(1, len(cols) + 1))
    dv = (dv * len(cols))[: len(cols)]
    order = np.argsort(dv)
    dv = np.asarray(dv, dtype=float)[order]
    rows = []
    for i, nm in enumerate(names):
        resp = np.nan_to_num(X[i][order].astype(float), nan=0.0)
        hi = float(top) if float(top) else float(np.max(resp)) or 1.0
        if np.max(resp) <= 0:
            continue
        norm = 100.0 * (1.0 - resp / hi)
        target = float(inhibition_target)
        ic = float("nan")
        for k in range(1, len(norm)):
            if (norm[k - 1] - target) * (norm[k] - target) <= 0 and norm[k] != norm[k - 1]:
                frac = (target - norm[k - 1]) / (norm[k] - norm[k - 1])
                ic = float(np.exp(np.log(max(dv[k - 1], 1e-9)) + frac * (np.log(max(dv[k], 1e-9))
                                                                      - np.log(max(dv[k - 1], 1e-9)))))
                break
        rows.append({"compound": nm, "max_response": round(float(hi), 4),
                    "inhibition_percent": round(float(norm[-1]), 3),
                    f"IC{int(target)}": round(ic, 6) if math.isfinite(ic) else "not reached",
                    "monotonic": bool(np.all(np.diff(resp) >= -1e-9)),
                    "points": int(len(resp))})
    key = f"IC{int(inhibition_target)}"
    rows.sort(key=lambda r: r[key] if isinstance(r[key], float) else float("inf"))
    return table(rows[: int(head)], f"{sum(1 for r in rows if isinstance(r[f'IC{int(inhibition_target)}'], float))}"
                                 f"/{len(rows)} compounds reach {inhibition_target}% inhibition")


@T("pharm_bliss_synergy", "Bliss and Loewe combination scores", PHARM, "table",
   [tbl("combination", COUNTS, "Table with single and combined effects"),
    textbox("drug_a_column", "", "Effect of drug A"), textbox("drug_b_column", "", "Effect of drug B"),
    textbox("combination_column", "", "Observed combination effect"),
    choice("model", ["bliss", "loewe", "both"], "both", "Reference model"), number("cutoff", 0.1,
                                                            "Synergy threshold (absolute difference)", min=0.0)],
   ex={"combination": COUNTS, "drug_a_column": "sample_A", "drug_b_column": "sample_B",
      "combination_column": "sample_D", "cutoff": 0.1}, up="CombroScore / synergyfinder",
   tags=("synergy", "drug combination"),
   summary="Compare the observed combination effect with the Bliss expectation per row.")
def pharm_bliss_synergy(combination, drug_a_column="", drug_b_column="", combination_column="", model="both",
                       cutoff=0.1):
    """Synergy scoring."""
    df = _tbl(combination)
    cols = [c for c, _v in _num_cols(df)]
    if len(cols) < 3:
        return table([], "need three numeric columns (A, B, combination)")
    ca = drug_a_column if drug_a_column in cols else cols[0]
    cb = drug_b_column if drug_b_column in cols else cols[1]
    cc = combination_column if combination_column in cols else cols[2]
    a = pd.to_numeric(df[ca], errors="coerce").to_numpy(dtype=float)
    b = pd.to_numeric(df[cb], errors="coerce").to_numpy(dtype=float)
    c = pd.to_numeric(df[cc], errors="coerce").to_numpy(dtype=float)
    key = str(df.columns[0])
    rows = []
    for i in range(len(df)):
        if not (math.isfinite(a[i]) and math.isfinite(b[i]) and math.isfinite(c[i])):
            continue
        ea, eb = max(0.0, min(1.0, float(a[i]))), max(0.0, min(1.0, float(b[i])))
        expected = ea + eb - ea * eb
        obs = max(0.0, min(1.0, float(c[i])))
        diff = obs - expected
        rec = {"compound": str(df.iloc[i][key]), "effect_a": round(ea, 4), "effect_b": round(eb, 4),
              "expected_bliss": round(expected, 4), "observed": round(obs, 4),
              "bliss_excess": round(diff, 4)}
        if model in ("loewe", "both"):
            ratio = (obs / ((ea + eb) / 2 or 1e-9))
            rec["loewe_ratio"] = round(float(ratio), 4)
        if abs(diff) < float(cutoff):
            rec["call"] = "additive"
        else:
            rec["call"] = "synergistic" if diff > 0 else "antagonistic"
        rows.append(rec)
    counts = Counter(r["call"] for r in rows)
    return table(rows, ", ".join(f"{k}={v}" for k, v in counts.most_common()) or "no usable rows")


@T("pharm_pk_noncompartmental", "Non-compartmental PK analysis", PHARM, "table",
   [tbl("pk", COUNTS, "Time-concentration table"), textbox("time_column", "", "Time column"),
    textbox("concentration_column", "", "Concentration column"), choice("dose_route", ["iv", "extravascular"],
                                                                 "iv", "Route"),
    number("dose", 100.0, "Dose (mass)", min=0.0), number("body_weight", 70.0, "Body weight (kg)", min=0.1),
    number("lambda_z_points", 3, "Points for the terminal slope", min=2)],
   ex={"pk": COUNTS, "time_column": "sample_A", "concentration_column": "sample_D", "dose": 100.0},
   up="PKSolver / noncompartamental analysis", tags=("pharmacokinetics", "AUC", "half-life"),
   summary="Cmax, Tmax, AUC by the linear-log trapezoid rule, clearance and half-life.")
def pharm_pk_noncompartmental(pk, time_column="", concentration_column="", dose_route="iv", dose=100.0,
                            body_weight=70.0, lambda_z_points=3):
    """NCA summary."""
    df = _tbl(pk)
    cols = [c for c, _v in _num_cols(df)]
    tc = time_column if time_column in cols else (cols[0] if cols else None)
    cc = concentration_column if concentration_column in cols else (cols[1] if len(cols) > 1 else tc)
    if not tc or not cc:
        return table([], "need numeric time and concentration columns")
    t = pd.to_numeric(df[tc], errors="coerce").to_numpy(dtype=float)
    c = pd.to_numeric(df[cc], errors="coerce").to_numpy(dtype=float)
    m = np.isfinite(t) & np.isfinite(c)
    order = np.argsort(t[m])
    t, c = t[m][order], np.abs(c[m][order])
    if t.size < 3:
        return table([], "need at least three measurable points")
    auc_last = float(_trapz(c, t))
    lam = 0.0
    tail = int(lambda_z_points)
    if t.size > tail + 1:
        tt, cc_ = t[-tail:], np.log(np.clip(c[-tail:], 1e-9, None))
        if np.std(tt) > 0:
            lam = float(-np.polyfit(tt, cc_, 1)[0])
    half = math.log(2) / lam if lam > 0 else float("nan")
    cl = float(dose) / auc_last if auc_last > 0 else float("nan")
    vd = cl * half if math.isfinite(cl) and math.isfinite(half) else float("nan")
    rows = [{"parameter": "Cmax", "value": round(float(np.max(c)), 6), "unit": "concentration"},
           {"parameter": "Tmax", "value": round(float(t[int(np.argmax(c))]), 4), "unit": "time"},
           {"parameter": "AUClast", "value": round(auc_last, 5), "unit": "conc*time"},
           {"parameter": "lambda_z", "value": round(lam, 6), "unit": "1/time"},
           {"parameter": "half_life", "value": round(half, 4) if math.isfinite(half) else "n/a", "unit": "time"},
           {"parameter": "AUCextrapolated_fraction", "value": round(float(np.max(c) / (lam * auc_last)), 4)
            if lam > 0 and auc_last > 0 else "n/a", "unit": "-"},
           {"parameter": "Clearance", "value": round(cl, 5) if math.isfinite(cl) else "n/a",
            "unit": "dose/(conc*time)"},
           {"parameter": "Volume_of_distribution", "value": round(vd, 5) if math.isfinite(vd) else "n/a",
            "unit": "L"},
           {"parameter": "Bioavailability_assumed", "value": 1.0 if dose_route == "iv" else "not estimable",
            "unit": "-"},
           {"parameter": "dose_normalised_AUC", "value": round(auc_last / (float(dose) or 1.0), 5), "unit": "-"},
           {"parameter": "clearance_per_kg", "value": round(cl / (float(body_weight) or 1.0), 6)
            if math.isfinite(cl) else "n/a", "unit": "-"}]
    res = table(rows, f"Cmax {round(float(np.max(c)), 3)} at t={round(float(t[int(np.argmax(c))]), 2)}; "
                    f"AUC {round(auc_last, 3)}" + (f"; t1/2 {round(half, 2)}" if math.isfinite(half) else ""))
    res["stats"]["figure"] = plot.line({"concentration": [float(v) for v in c]}, x=[float(v) for v in t],
                                     xlabel="time", ylabel="concentration", title="Plasma profile")
    return res


def _trapz(y, x):
    y, x = np.asarray(y, dtype=float), np.asarray(x, dtype=float)
    if y.size < 2:
        return 0.0
    if hasattr(np, "trapezoid"):
        return float(np.trapezoid(y, x))
    return float(np.sum((y[1:] + y[:-1]) * 0.5 * np.diff(x)))


@T("pharm_pk_compartment_model", "One- and two-compartment model fit", PHARM, "table",
   [tbl("pk", COUNTS, "Time-concentration table"), textbox("time_column", "", "Time column"),
    textbox("concentration_column", "", "Concentration column"), choice("model", ["one", "two"], "one",
                                                                    "Compartment model"),
    number("dose", 100.0, "Dose", min=0.0), intin("iterations", 400, "Grid refinements", min=50)],
   ex={"pk": COUNTS, "time_column": "sample_A", "concentration_column": "sample_D"},
   up="nlmixr / Phoenix WinNonlin model fit", tags=("pharmacokinetics", "modelling"),
   summary="Fit a simple IV bolus (or two-compartment) curve and report the rate constants.")
def pharm_pk_compartment_model(pk, time_column="", concentration_column="", model="one", dose=100.0,
                              iterations=400):
    """Compartment model fit."""
    df = _tbl(pk)
    cols = [c for c, _v in _num_cols(df)]
    tc = time_column if time_column in cols else cols[0]
    cc = concentration_column if concentration_column in cols else (cols[1] if len(cols) > 1 else cols[0])
    t = pd.to_numeric(df[tc], errors="coerce").to_numpy(dtype=float)
    y = np.abs(pd.to_numeric(df[cc], errors="coerce").to_numpy(dtype=float))
    m = np.isfinite(t) & np.isfinite(y) & (t >= 0)
    t, y = t[m], y[m]
    if t.size < 3:
        return table([], "need at least three points")
    ymax = float(np.max(y)) or 1.0
    best = None
    rng = np.linspace(0.01, 2.0, max(20, int(iterations) // 4))
    for k in rng:
        pred = ymax * np.exp(-k * t)
        sse = float(np.sum((y - pred) ** 2))
        if best is None or sse < best[0]:
            best = (sse, float(k))
    sse1, k10 = best
    if model == "two":
        best2 = None
        for a in np.linspace(0.05, 3.0, 20):
            for b in np.linspace(0.005, 0.9, 20):
                pred = ymax * (0.6 * np.exp(-a * t) + 0.4 * np.exp(-b * t))
                s2 = float(np.sum((y - pred) ** 2))
                if best2 is None or s2 < best2[0]:
                    best2 = (s2, float(a), float(b))
        sse, alpha, beta = best2
        half = math.log(2) / beta if beta > 0 else float("nan")
    else:
        alpha = beta = k10
        sse = sse1
        half = math.log(2) / k10 if k10 > 0 else float("nan")
    tss = float(np.sum((y - np.mean(y)) ** 2)) or 1.0
    vol = (float(dose) or 1.0) / ymax
    rows = [{"parameter": "model", "value": f"{model}-compartment IV bolus", "unit": "-"},
           {"parameter": "C0", "value": round(ymax, 5), "unit": "concentration"},
           {"parameter": "alpha", "value": round(alpha, 5), "unit": "1/time"},
           {"parameter": "beta", "value": round(beta, 5), "unit": "1/time"},
           {"parameter": "k10", "value": round(k10, 5), "unit": "1/time"},
           {"parameter": "half_life_beta", "value": round(half, 4) if math.isfinite(half) else "n/a",
            "unit": "time"},
           {"parameter": "volume_of_distribution", "value": round(vol, 5), "unit": "dose/conc"},
           {"parameter": "clearance", "value": round(vol * k10, 5), "unit": "dose/(conc*time)"},
           {"parameter": "r_squared", "value": round(1 - sse / tss, 5), "unit": "-"},
           {"parameter": "one_compartment_r_squared", "value": round(1 - sse1 / tss, 5), "unit": "-"}]
    res = table(rows, f"{model}-compartment fit: R2 {round(1 - sse / tss, 3)}, "
                   f"t1/2 {round(half, 2) if math.isfinite(half) else 'n/a'}")
    res["stats"]["figure"] = plot.line({"observed": [float(v) for v in y],
                                      "fitted": [float(v) for v in (ymax * (0.6 * np.exp(-alpha * t)
                                                                         + 0.4 * np.exp(-beta * t))
                                                                  if model == "two" else ymax * np.exp(-k10 * t))]},
                                     x=[float(v) for v in t], xlabel="time", ylabel="concentration",
                                     title=f"{model}-compartment fit")
    return res


# ===========================================================================
# multi-omics integration
# ===========================================================================
@T("multiomics_set_overlap", "Overlap of feature sets across omics layers", MULTI, "table",
   [tbl("table_a", COUNTS, "First feature table"), tbl("table_b", COUNTS, "Second feature table"),
    tbl("table_c", "", "Third table (optional)"), textbox("key_column", "", "Key column"),
    boolean("jaccard", True, "Report Jaccard and Sorensen indices")],
   ex={"table_a": COUNTS, "table_b": COUNTS, "table_c": COUNTS},
   up="mixOmics / clusterProfiler compareSets", tags=("integration", "sets", "overlap"),
   summary="Shared and unique feature names between two or three omics tables.")
def multiomics_set_overlap(table_a, table_b, table_c="", key_column="", jaccard=True):
    """Cross-omics set comparison."""
    def keys(src):
        if not str(src or "").strip():
            return set()
        df = _tbl(src)
        cols = [str(c) for c in df.columns]
        kc = key_column if key_column in cols else cols[0]
        return {str(v) for v in df[kc].tolist()}
    A, B = keys(table_a), keys(table_b)
    C = keys(table_c)
    rows = [{"comparison": "A vs B", "in_both": len(A & B), "only_a": len(A - B), "only_b": len(B - A),
            "union": len(A | B), "jaccard": round(float(stats.jaccard_sets(sorted(A), sorted(B))["jaccard"]), 4)
            if jaccard else "",
            "sorensen": round(float(stats.jaccard_sets(sorted(A), sorted(B))["sorensen"]), 4) if jaccard else ""}]
    if C:
        rows.append({"comparison": "A vs B vs C", "in_both": len(A & B & C), "only_a": len(A - B - C),
                    "only_b": len(B - A - C), "union": len(A | B | C),
                    "jaccard": round(len(A & B & C) / max(1, len(A | B | C)), 4), "sorensen": ""})
        rows.append({"comparison": "B vs C", "in_both": len(B & C), "only_a": len(B - C), "only_b": len(C - B),
                    "union": len(B | C),
                    "jaccard": round(float(stats.jaccard_sets(sorted(B), sorted(C))["jaccard"]), 4) if jaccard else "",
                    "sorensen": round(float(stats.jaccard_sets(sorted(B), sorted(C))["sorensen"]), 4) if jaccard else ""})
    return table(rows, f"{len(A & B)} features shared by the first two layers")


@T("multiomics_integrated_score", "Z-score integration across omics layers", MULTI, "table",
   [tbl("table_a", COUNTS, "Layer A"), tbl("table_b", COUNTS, "Layer B"), choice("weight_a", ["1", "2", "3"], "1",
                                                                            "Weight A"), choice("weight_b", ["1", "2",
                                                                                                    "3"], "1",
                                                                                            "Weight B"),
    boolean("rank_instead", False, "Use ranks instead of z-scores"), intin("head", 100, "Rows", min=1)],
   ex={"table_a": COUNTS, "table_b": COUNTS, "weight_a": "2", "weight_b": "1"},
   up="mixOmics DIABLO score / iCluster", tags=("integration", "scoring"),
   summary="Combine two feature tables into one ranked multi-omics score per feature.")
def multiomics_integrated_score(table_a, table_b, weight_a="1", weight_b="1", rank_instead=False, head=100):
    """Weighted integration."""
    na, ca, Xa = _mat(table_a)
    nb, cb, Xb = _mat(table_b)
    if Xa.size == 0 or Xb.size == 0:
        return table([], "both layers need numeric columns")
    da = np.nanmean(np.where(np.isfinite(Xa), Xa, np.nan), axis=1)
    db = np.nanmean(np.where(np.isfinite(Xb), Xb, np.nan), axis=1)
    pos_a = {n: i for i, n in enumerate(na)}
    pos_b = {n: i for i, n in enumerate(nb)}
    shared = [n for n in na if n in pos_b]
    if not shared:
        return table([], "the two layers share no feature names")

    def norm(v):
        if rank_instead:
            r = np.argsort(np.argsort(v)).astype(float)
            return (r - r.mean()) / (r.std() or 1.0)
        return (v - np.nanmean(v)) / (np.nanstd(v) or 1.0)
    va = norm(np.array([da[pos_a[n]] for n in shared], dtype=float))
    vb = norm(np.array([db[pos_b[n]] for n in shared], dtype=float))
    wa, wb = float(weight_a), float(weight_b)
    score = (wa * va + wb * vb) / (wa + wb)
    rows = [{"feature": n, "layer_a_z": round(float(va[i]), 4), "layer_b_z": round(float(vb[i]), 4),
            "integrated_score": round(float(score[i]), 4), "weights": f"{wa}:{wb}",
            "rank": 0} for i, n in enumerate(shared)]
    rows.sort(key=lambda r: -r["integrated_score"])
    for i, r in enumerate(rows, 1):
        r["rank"] = i
    return table(rows[: int(head)], f"{len(shared)} features integrated; top: {rows[0]['feature']} "
                                 f"({round(rows[0]['integrated_score'], 3)})")


@T("multiomics_cross_correlation", "Cross-omics correlation block", MULTI, "figure",
   [tbl("table_a", COUNTS, "Layer A"), tbl("table_b", COUNTS, "Layer B"),
    choice("method", ["pearson", "spearman"], "spearman", "Method"), intin("max_rows", 12, "Rows", min=2),
    intin("max_cols", 12, "Columns", min=2)],
   ex={"table_a": COUNTS, "table_b": COUNTS}, up="mixOmics correlogram / MOFA loadings",
   tags=("integration", "correlation", "heatmap"),
   summary="Correlation heatmap between the features of two omics layers.")
def multiomics_cross_correlation(table_a, table_b, method="spearman", max_rows=12, max_cols=12):
    """Cross-block correlation matrix."""
    na, ca, Xa = _mat(table_a)
    nb, cb, Xb = _mat(table_b)
    if Xa.size == 0 or Xb.size == 0:
        return plot.empty_plot("both layers need numeric columns")
    fa = np.nanmean(np.where(np.isfinite(Xa), Xa, np.nan), axis=1)[: int(max_rows)]
    fb = np.nanmean(np.where(np.isfinite(Xb), Xb, np.nan), axis=1)[: int(max_cols)]
    A = Xa[: len(fa)].T
    B = Xb[: len(fb)].T
    M = np.zeros((A.shape[1], B.shape[1]))
    for i in range(A.shape[1]):
        for j in range(B.shape[1]):
            x, y = np.nan_to_num(A[:, i]), np.nan_to_num(B[:, j])
            st = (stats.spearman if method == "spearman" else stats.pearson)(list(x), list(y))
            M[i, j] = float(dict(st).get("r", 0.0))
    return plot.heatmap([[round(float(v), 4) for v in row] for row in M],
                       row_labels=na[: A.shape[1]], col_labels=nb[: B.shape[1]],
                       title=f"{method} cross-omics correlation", cmap="coolwarm", annot=True)


# ===========================================================================
# single-cell import and manipulation
# ===========================================================================
@T("sc_matrix_market_preview", "Inspect a Matrix Market (mtx) file", SCIMP, "table",
   [txt("matrix", "%%MatrixMarket matrix coordinate real general\n3 4 6\n1 1 4.0\n1 2 1.0\n2 3 7.0\n"
               "3 1 2.0\n3 4 5.0\n2 2 3.0", "Matrix Market text"), intin("head", 20, "Rows", min=1)],
   ex={"matrix": "%%MatrixMarket matrix coordinate real general\n3 4 6\n1 1 4.0\n1 2 1.0\n2 3 7.0\n"
               "3 1 2.0\n3 4 5.0\n2 2 3.0"},
   up="10x / Seurat ReadMM", tags=("single cell", "mtx", "matrix"),
   summary="Shape, density and per-row totals of a sparse count matrix.")
def sc_matrix_market_preview(matrix, head=20):
    """MM parsing and summary."""
    lines = [l for l in io.as_text(matrix).splitlines() if l.strip() and not l.startswith("%")]
    if not lines:
        return table([], "empty matrix")
    hdr = lines[0].split()
    try:
        nrows, ncols, nnz = int(hdr[0]), int(hdr[1]), int(hdr[2])
    except (ValueError, IndexError):
        nrows, ncols = len(lines), 1
        nnz = len(lines)
    rows_ = defaultdict(list)
    for ln in lines[1:]:
        p = ln.split()
        if len(p) < 3:
            continue
        try:
            i, j, v = int(p[0]), int(p[1]), float(p[2])
        except ValueError:
            continue
        rows_[i].append((j, v))
    totals = {i: sum(v for _j, v in vals) for i, vals in rows_.items()}
    out = [{"metric": "rows", "value": nrows}, {"metric": "columns", "value": ncols},
          {"metric": "stored_entries", "value": nnz},
          {"metric": "parsed_entries", "value": sum(len(v) for v in rows_.values())},
          {"metric": "density_percent", "value": round(100 * sum(len(v) for v in rows_.values())
                                                     / max(1, nrows * ncols), 4)},
          {"metric": "max_value", "value": round(max([v for vals in rows_.values() for _j, v in vals], default=0.0), 4)},
          {"metric": "mean_per_row", "value": round(float(np.mean(list(totals.values()))), 4) if totals else 0.0}]
    for i in sorted(rows_)[: int(head)]:
        out.append({"metric": f"row_{i}_total", "value": round(float(totals[i]), 4),
                   "entries": len(rows_[i]), "columns": ",".join(str(j) for j, _v in sorted(rows_[i])[:12])})
    return table(out, f"matrix {nrows} x {ncols}, {sum(len(v) for v in rows_.values())} non-zero entries")


@T("sc_image_to_cells", "Extract per-cell intensities from an image", SCIMP, "table",
   [img("image_file", IMAGE, "Image (PGM/PPM)"), number("threshold", 0.0, "Binarisation threshold (0 = Otsu)",
                                                    min=0.0), intin("min_area", 4, "Minimum object area", min=1),
    choice("intensity", ["mean", "max", "sum"], "mean", "Per-object statistic"), intin("head", 100, "Rows", min=1)],
   ex={"image_file": IMAGE, "min_area": 3}, up="CellProfiler / Scanpy image stats",
   tags=("imaging", "cells", "segmentation"),
   summary="Segment bright objects in a micrograph and tabulate their intensities.")
def sc_image_to_cells(image_file, threshold=0.0, min_area=4, intensity="mean", head=100):
    """Simple cell segmentation table."""
    arr = image.to_gray(image.read_image(image_file))
    thr = float(threshold)
    if thr <= 0:
        thr = float(dict(image.threshold_otsu(arr)).get("threshold", 0.5))
    mask = image.binarize(image.gaussian_blur(arr, sigma=1.0), value=thr, above=True)
    lab = image.label(image.remove_small_objects(image.binary_dilation(mask, iterations=1),
                                                min_area=int(min_area)), connectivity=2)
    labeled = lab["labels"] if isinstance(lab, dict) else lab
    props = image.regionprops(labeled, arr)
    rows = []
    for p in props[: int(head)]:
        d = dict(p)
        area = float(d.get("area", 0) or 0)
        rows.append({"cell": int(d.get("label", len(rows) + 1)), "area": int(area),
                    "centroid_x": round(float(d.get("centroid_x", 0) or 0), 2),
                    "centroid_y": round(float(d.get("centroid_y", 0) or 0), 2),
                    "mean_intensity": round(float(d.get("mean_intensity", 0.0) or 0.0), 4),
                    "equivalent_diameter": round(float(d.get("equivalent_diameter", 0.0) or 0.0), 3),
                    "solidity": round(float(d.get("solidity", 0.0) or 0.0), 4),
                    "eccentricity": round(float(d.get("eccentricity", 0.0) or 0.0), 4),
                    "bbox": f"{int(d.get('bbox_xmin', 0))},{int(d.get('bbox_ymin', 0))},"
                            f"{int(d.get('bbox_xmax', 0))},{int(d.get('bbox_ymax', 0))}",
                    "reported": intensity, "total_intensity": round(area * float(d.get("mean_intensity", 0.0) or 0.0), 3)})
    res = table(rows, f"{len(rows)} objects segmented at threshold {round(thr, 4)} "
                    f"(min area {min_area}, {intensity} intensity)")
    res["stats"]["figure"] = plot.image_show(image.add_labels_overlay(labeled), title="segmented cells")
    res["stats"]["n_objects"] = int(lab.get("n_objects", len(rows))) if isinstance(lab, dict) else len(rows)
    return res


@T("sc_barcode_whitelist", "Collapse barcodes that differ by one base", SCIMP, "table",
   [txt("barcodes", "ACGTACGT\nACGTACGA\nTTTTGGGG\nTTTTGGGA", "Barcode list (one per line)"),
    intin("min_count", 1, "Minimum reads per barcode", min=1), intin("mismatch", 1, "Allowed differences", min=0)],
   ex={"barcodes": "ACGTACGT\nACGTACGA\nTTTTGGGG\nTTTTGGGA"}, up="umi_tools whitelist / kb-tools",
   tags=("single cell", "barcodes", "UMI"),
   summary="Group near-identical barcodes so low-count neighbours are merged in.")
def sc_barcode_whitelist(barcodes, min_count=1, mismatch=1):
    """Barcode collapsing."""
    lines = [l.strip().upper() for l in io.as_text(barcodes).splitlines() if l.strip()]
    counts = Counter(lines)
    order = sorted(counts, key=lambda b: -counts[b])
    assigned: dict[str, str] = {}
    for b in order:
        parent = None
        for a in assigned:
            if len(a) == len(b) and int(dict(align.hamming(a, b))["distance"]) <= int(mismatch) \
                    and counts[a] >= counts[b]:
                parent = a
                break
        assigned[b] = parent or b
    groups = defaultdict(list)
    for b, p in assigned.items():
        groups[p].append(b)
    rows = [{"whitelisted_barcode": p, "merged_barcode": b, "reads": counts[b],
            "group_total": sum(counts[x] for x in groups[p]),
            "n_merged": len(groups[p]), "distance_to_parent":
            int(dict(align.hamming(p, b))["distance"]) if len(p) == len(b) else 0}
           for p in groups for b in sorted(groups[p]) if sum(counts[x] for x in groups[p]) >= int(min_count)]
    return table(rows, f"{len(groups)} whitelisted barcodes from {len(counts)} observed "
                     f"(<= {mismatch} mismatches merged)")


@T("sc_matrix_to_long", "Sparse coordinates to a long table", SCIMP, "table",
   [txt("matrix", "%%MatrixMarket matrix coordinate real general\n3 4 6\n1 1 4.0\n1 2 1.0\n2 3 7.0\n"
               "3 1 2.0\n3 4 5.0\n2 2 3.0", "Matrix Market text"), textbox("gene_file", "", "Optional gene names"),
    textbox("cell_file", "", "Optional cell barcodes"), boolean("drop_zeros", True, "Drop zeros"),
    intin("head", 300, "Rows", min=1)],
   ex={"matrix": "%%MatrixMarket matrix coordinate real general\n2 2 3\n1 1 5.0\n1 2 2.0\n2 1 1.0"},
   up="10x read10x / Seurat as.data.frame", tags=("single cell", "sparse", "tidy"),
   summary="Flatten a sparse count matrix into cell, feature and count columns.")
def sc_matrix_to_long(matrix, gene_file="", cell_file="", drop_zeros=True, head=300):
    """Sparse -> long."""
    lines = [l for l in io.as_text(matrix).splitlines() if l.strip() and not l.startswith("%")]
    if not lines:
        return table([], "empty matrix")
    hdr = lines[0].split()
    try:
        nrows, ncols = int(hdr[0]), int(hdr[1])
    except (ValueError, IndexError):
        nrows, ncols = len(lines), 1
    genes = [f"gene_{i + 1}" for i in range(nrows)]
    cells = [f"cell_{j + 1}" for j in range(ncols)]
    if gene_file:
        names = [l.strip() for l in io.as_text(gene_file).splitlines() if l.strip()]
        genes = (names + genes)[: max(len(names), nrows)]
    if cell_file:
        names = [l.strip() for l in io.as_text(cell_file).splitlines() if l.strip()]
        cells = (names + cells)[: max(len(names), ncols)]
    rows = []
    for ln in lines[1:]:
        p = ln.split()
        if len(p) < 3:
            continue
        try:
            i, j, v = int(p[0]), int(p[1]), float(p[2])
        except ValueError:
            continue
        if drop_zeros and v == 0:
            continue
        rows.append({"gene": genes[i - 1] if i - 1 < len(genes) else f"gene_{i}",
                    "cell": cells[j - 1] if j - 1 < len(cells) else f"cell_{j}",
                    "count": v, "gene_index": i, "cell_index": j,
                    "log1p": round(math.log1p(max(0.0, v)), 5)})
    return table(rows[: int(head)], f"{len(rows)} non-zero entries "
                                 f"({len(set(r['cell'] for r in rows))} cells, {len(set(r['gene'] for r in rows))} genes)")


@T("sc_collate_samples", "Merge per-sample count tables into one object", SCIMP, "table",
   [tbl("table_a", COUNTS, "First table"), tbl("table_b", COUNTS, "Second table"),
    textbox("sample_a", "A", "Label for the first"), textbox("sample_b", "B", "Label for the second"),
    choice("mode", ["cbind", "rbind"], "cbind", "Merge direction")],
   ex={"table_a": COUNTS, "table_b": COUNTS, "sample_a": "A", "sample_b": "B"},
   up="Seurat merge / SingleCellExperiment cbind", tags=("single cell", "merge"),
   summary="Concatenate two count tables by cells or by features, keeping provenance.")
def sc_collate_samples(table_a, table_b, sample_a="A", sample_b="B", mode="cbind"):
    """Merge two tables."""
    da, db = _tbl(table_a), _tbl(table_b)
    if da.empty or db.empty:
        return table([], "both tables need data")
    if mode == "rbind":
        out = pd.concat([da, db], ignore_index=True)
        out["source"] = [sample_a] * len(da) + [sample_b] * len(db)
        return table(out.to_dict("records")[:400], f"stacked {len(da)} + {len(db)} rows by feature")
    key = str(da.columns[0])
    left = da.set_index(key)
    right = db.set_index(str(db.columns[0]))
    common = [c for c in right.columns if c not in left.columns]
    merged = left.join(right[common], how="outer")
    merged = merged.rename(columns={c: f"{sample_b}_{c}" for c in common})
    merged[sample_a + "_total"] = left.sum(axis=1, numeric_only=True)
    merged[sample_b + "_total"] = right[common].sum(axis=1, numeric_only=True) if common else 0
    return table(merged.reset_index().to_dict("records")[:400],
                 f"merged {len(merged)} features x {len(merged.columns)} columns")


# ===========================================================================
# Scanpy / Monocle3 / PreSTO / Seurat extras
# ===========================================================================
@T("scanpy_qc_metrics", "Scanpy-style per-cell QC metrics", SCANPY, "table",
   [tbl("counts", COUNTS, "Count matrix (rows = genes)"), textbox("mitochondrial_prefix", "MT-",
                                                              "Mitochondrial gene prefix"),
    textbox("ribosomal_prefix", "RP", "Ribosomal gene prefix"), number("mt_cut", 20.0, "Max mito %", min=0.0),
    intin("min_genes", 1, "Minimum genes per cell", min=0), intin("head", 100, "Rows", min=1)],
   ex={"counts": COUNTS, "mitochondrial_prefix": "geneA", "ribosomal_prefix": "geneB"},
   up="scanpy qc.calculate_qc_metrics", tags=("single cell", "quality control"),
   summary="Genes, UMIs, mito and ribosomal fractions per cell with filter flags.")
def scanpy_qc_metrics(counts, mitochondrial_prefix="MT-", ribosomal_prefix="RP", mt_cut=20.0, min_genes=1,
                     head=100):
    """Per-cell QC metrics."""
    names, cols, X = _mat(counts)
    if X.size == 0:
        return table([], "no numeric columns")
    X = np.nan_to_num(X.astype(float), nan=0.0)
    mt = [i for i, n in enumerate(names) if str(n).upper().startswith(str(mitochondrial_prefix).upper())]
    rb = [i for i, n in enumerate(names) if str(n).upper().startswith(str(ribosomal_prefix).upper())]
    rows = []
    for j, c in enumerate(cols):
        v = X[:, j]
        tot = float(v.sum()) or 1.0
        mtp = 100 * float(v[mt].sum()) / tot if mt else 0.0
        rbp = 100 * float(v[rb].sum()) / tot if rb else 0.0
        rows.append({"cell": c, "n_genes": int((v > 0).sum()), "n_counts": round(tot, 1),
                    "pct_mito": round(mtp, 3), "pct_ribosomal": round(rbp, 3),
                    "genes_exceeding_1": int((v >= 1).sum()), "max_gene_fraction": round(float(v.max() / tot), 4),
                    "pass_filters": bool((v > 0).sum() >= int(min_genes) and mtp <= float(mt_cut))})
    n_pass = sum(1 for r in rows if r["pass_filters"])
    res = table(rows[: int(head)], f"{n_pass}/{len(rows)} cells pass (mito <= {mt_cut}%, "
                                f"genes >= {min_genes})")
    res["stats"] = {"mt_genes": [names[i] for i in mt], "ribosomal_genes": [names[i] for i in rb]}
    return res


@T("scanpy_highly_variable", "Select highly variable genes", SCANPY, "table",
   [tbl("counts", COUNTS, "Count matrix (rows = genes)"), choice("flavour", ["seurat", "cell_ranger",
                                                        "seurat_v3"], "seurat", "Flavour"),
    intin("n_bins", 6, "Mean-expression bins", min=2), intin("top_genes", 10, "Genes to keep", min=1),
    boolean("show_plot", True, "Show the dispersion plot")],
   ex={"counts": COUNTS, "n_bins": 4, "top_genes": 5}, up="scanpy pp.highly_variable_genes",
   tags=("single cell", "feature selection"),
   summary="Binned dispersion ranking used to pick informative genes before PCA.")
def scanpy_highly_variable(counts, flavour="seurat", n_bins=6, top_genes=10, show_plot=True):
    """HVG selection."""
    names, cols, X = _mat(counts)
    if X.size == 0:
        return table([], "no numeric columns")
    M = np.nan_to_num(X.astype(float), nan=0.0)
    if flavour != "seurat_v3":
        M = np.log1p(M / np.where(M.sum(axis=0) > 0, M.sum(axis=0), 1.0) * float(np.mean(M.sum(axis=0))))
    mu = M.mean(axis=1)
    sd = M.std(axis=1, ddof=1) if M.shape[0] > 1 else np.zeros_like(mu)
    disp = sd / np.where(mu > 0, mu, 1e-9)
    nb = max(2, int(n_bins))
    qs = np.quantile(mu, np.linspace(0, 1, nb + 1))
    rows = []
    for i, nm in enumerate(names):
        b = int(np.searchsorted(qs, mu[i], side="right") - 1)
        b = max(0, min(nb - 1, b))
        rows.append({"gene": nm, "mean_expression": round(float(mu[i]), 5), "dispersion": round(float(disp[i]), 5),
                    "bin": b, "dispersion_normalised": 0.0, "highly_variable": False})
    by_bin = defaultdict(list)
    for r in rows:
        by_bin[r["bin"]].append(r["dispersion"])
    med = {b: float(np.median(v)) if v else 1.0 for b, v in by_bin.items()}
    for r in rows:
        m = med.get(r["bin"], 1.0) or 1.0
        r["dispersion_normalised"] = round((r["dispersion"] - m) / m, 5)
    rows.sort(key=lambda r: -r["dispersion_normalised"])
    for i, r in enumerate(rows):
        r["highly_variable"] = i < int(top_genes)
        r["rank"] = i + 1
    res = table(rows, f"{int(top_genes)} highly variable genes of {len(rows)} (flavour {flavour})")
    if show_plot:
        res["stats"]["figure"] = plot.scatter([r["mean_expression"] for r in rows],
                                            [r["dispersion"] for r in rows],
                                            labels=[r["gene"] for r in rows],
                                            xlabel="mean expression", ylabel="dispersion",
                                            title="gene dispersion")
    return res


@T("scanpy_cluster_knn_graph", "k-NN graph clustering (SNN + label propagation)", SCANPY, "table",
   [tbl("counts", COUNTS, "Count matrix"), intin("neighbours", 5, "Neighbours per cell", min=2),
    choice("reduction", ["pca", "raw"], "pca", "Embedding first"), intin("resolution", 2,
                                                                      "Merging strength (shared nearest neighbours)",
                                                                      min=1),
    boolean("standardize", True, "z-score genes"), intin("head", 200, "Rows", min=1)],
   ex={"counts": COUNTS, "neighbours": 4, "resolution": 2}, up="scanpy tl.leiden / scran buildSNNGraph",
   tags=("single cell", "clustering", "graph"),
   summary="Cluster cells from a shared-nearest-neighbour graph without external packages.")
def scanpy_cluster_knn_graph(counts, neighbours=5, reduction="pca", resolution=2, standardize=True, head=200):
    """SNN label propagation clustering."""
    names, cols, X = _mat(counts)
    if X.size == 0 or len(cols) < 2:
        return table([], "need a matrix with at least two cells")
    M = np.nan_to_num(X.T.astype(float), nan=0.0)
    if standardize:
        mu, sd = M.mean(axis=0), M.std(axis=0)
        M = (M - mu) / np.where(sd > 0, sd, 1)
    if reduction == "pca" and M.shape[1] > 2:
        res = ml.pca(M, components=min(M.shape[1], max(2, M.shape[0] - 1)))
        Y = np.asarray(res["scores"], dtype=float)
    else:
        Y = M
    k = max(2, min(int(neighbours), Y.shape[0]))
    D = np.linalg.norm(Y[:, None, :] - Y[None, :, :], axis=2)
    nn = np.argsort(D, axis=1)[:, :k]
    shared = np.zeros((Y.shape[0], Y.shape[0]))
    for i in range(Y.shape[0]):
        for j in nn[i]:
            shared[i, j] = shared[j, i] = 1.0
    snn = np.zeros_like(shared)
    for i in range(Y.shape[0]):
        for j in range(i + 1, Y.shape[0]):
            w = len(set(nn[i]) & set(nn[j]))
            snn[i, j] = snn[j, i] = w
    lab = {i: i for i in range(Y.shape[0])}
    for _ in range(max(1, int(resolution))):
        changed = False
        for i in range(Y.shape[0]):
            votes = Counter(lab[j] for j in range(Y.shape[0]) if snn[i, j] > 0)
            if votes:
                new = votes.most_common(1)[0][0]
                if lab[i] != new:
                    lab[i] = new
                    changed = True
        if not changed:
            break
    sizes = Counter(lab.values())
    rows = [{"cell": c, "cluster": lab[i], "cluster_size": sizes[lab[i]],
            "shared_neighbours": int(np.sum(snn[i] > 0)),
            "nearest_cell": cols[int(nn[i][1])] if nn.shape[1] > 1 else c}
           for i, c in enumerate(cols)]
    out = table(rows[: int(head)], f"{len(sizes)} clusters over {len(cols)} cells "
                                f"(k={k}, resolution {resolution})")
    out["stats"] = {"cluster_sizes": dict(sizes)}
    return out


@T("scanpy_cell_cycle_score", "Cell-cycle phase scores", SCANPY, "table",
   [txt("s_phase_genes", "geneA\nGeneB\nGENEC", "S-phase gene list"), txt("g2m_genes", "geneD\ngeneE",
                                                                     "G2M gene list"),
    tbl("counts", COUNTS, "Expression matrix (log or counts)"), choice("normalisation", ["zscore", "none"],
                                                                    "zscore", "Score normalisation"),
    intin("head", 100, "Rows", min=1)],
   ex={"s_phase_genes": "geneA\ngeneB", "g2m_genes": "geneD\ngeneE", "counts": COUNTS},
   up="scanpy tl.score_genes / cell cycle scoring", tags=("single cell", "cell cycle"),
   summary="Score S and G2M signatures per cell and call a phase for each one.")
def scanpy_cell_cycle_score(s_phase_genes, g2m_genes, counts, normalisation="zscore", head=100):
    """Cell-cycle scoring."""
    names, cols, X = _mat(counts)
    s = {g.strip().upper() for g in io.as_text(s_phase_genes).splitlines() if g.strip()}
    g2 = {g.strip().upper() for g in io.as_text(g2m_genes).splitlines() if g.strip()}
    if X.size == 0 or (not s and not g2):
        return table([], "need an expression matrix and at least one gene list")
    M = X.astype(float)
    if normalisation == "zscore":
        mu = np.nanmean(M, axis=1, keepdims=True)
        sd = np.nanstd(M, axis=1, keepdims=True)
        M = (M - mu) / np.where(sd > 0, sd, 1)
    idx_s = [i for i, n in enumerate(names) if str(n).upper() in s]
    idx_g = [i for i, n in enumerate(names) if str(n).upper() in g2]
    rows = []
    for j, c in enumerate(cols):
        sv = float(np.nanmean(M[idx_s, j])) if idx_s else 0.0
        gv = float(np.nanmean(M[idx_g, j])) if idx_g else 0.0
        rows.append({"cell": c, "S_score": round(sv, 5), "G2M_score": round(gv, 5),
                    "difference": round(sv - gv, 5),
                    "phase": "S" if sv - gv > 0.1 else "G2M" if gv - sv > 0.1 else "G1",
                    "s_genes_used": len(idx_s), "g2m_genes_used": len(idx_g)})
    ph = Counter(r["phase"] for r in rows)
    res = table(rows[: int(head)], "phases: " + ", ".join(f"{k}={v}" for k, v in ph.most_common()))
    res["stats"]["figure"] = plot.scatter([r["S_score"] for r in rows], [r["G2M_score"] for r in rows],
                                        labels=[r["cell"] for r in rows], xlabel="S score",
                                        ylabel="G2M score", title="cell cycle")
    return res


@T("monocle3_pseudotime", "Order cells along a principal graph", MONO, "table",
   [tbl("counts", COUNTS, "Expression matrix"), intin("neighbours", 4, "Nearest neighbours", min=2),
    choice("start", ["first", "highest_loading", "lowest_total"], "first", "Starting cell"),
    intin("components", 2, "PCA dimensions", min=1), boolean("backtrack", True, "Smooth the trajectory")],
   ex={"counts": COUNTS, "neighbours": 3}, up="monocle3 learn_graph / order_cells",
   tags=("pseudotime", "trajectory"),
   summary="Order single cells along a nearest-neighbour graph and report pseudotime.")
def monocle3_pseudotime(counts, neighbours=4, start="first", components=2, backtrack=True):
    """Principal-graph-like ordering."""
    names, cols, X = _mat(counts)
    if X.size == 0 or len(cols) < 3:
        return table([], "need at least three cells")
    M = np.nan_to_num(X.T.astype(float), nan=0.0)
    mu, sd = M.mean(axis=0), M.std(axis=0)
    M = (M - mu) / np.where(sd > 0, sd, 1)
    res = ml.pca(M, components=min(int(components), M.shape[1]))
    Y = np.asarray(res["scores"], dtype=float)
    k = max(2, min(int(neighbours), Y.shape[0]))
    D = np.linalg.norm(Y[:, None, :] - Y[None, :, :], axis=2)
    np.fill_diagonal(D, np.inf)
    nn = np.argsort(D, axis=1)[:, :k]
    root = 0
    if start == "highest_loading":
        root = int(np.argmax(np.abs(np.ravel(np.asarray(res["loadings"], dtype=float)))))
        root = int(np.argmax(Y[:, 0]))
    elif start == "lowest_total":
        root = int(np.argmin(np.nan_to_num(X.astype(float)).sum(axis=0)))
    order = []
    dist = {root: 0.0}
    frontier = [root]
    while frontier:
        cur = min(frontier, key=lambda i: dist[i])
        frontier.remove(cur)
        order.append(cur)
        for nb in nn[cur]:
            nb = int(nb)
            if nb in dist:
                continue
            dist[nb] = dist[cur] + float(np.linalg.norm(Y[cur] - Y[nb]))
            frontier.append(nb)
    for i in range(Y.shape[0]):
        if i not in order:
            order.append(i)
            dist[i] = max(dist.values() or [0.0]) + 1.0
    pt = np.array([dist[i] for i in order], dtype=float)
    if backtrack and pt.size > 2:
        sm = np.asarray(stats.ewma(list(pt), alpha=0.5), dtype=float)
        pt = (pt + sm) / 2
    span = float(pt.max() - pt.min()) or 1.0
    rank = {cell: i for i, cell in enumerate(order)}
    rows = [{"cell": cols[i], "pseudotime": round(float(pt[rank[i]] - pt.min()), 5),
            "trajectory_rank": int(rank[i]) + 1,
            "dimension_1": round(float(Y[i, 0]), 4),
            "dimension_2": round(float(Y[i, 1]), 4) if Y.shape[1] > 1 else 0.0} for i in range(len(cols))]
    res_tbl = table(rows, f"pseudotime over {len(rows)} cells from {cols[root]} "
                        f"(k={k}, span {round(span, 3)})")
    res_tbl["stats"]["figure"] = plot.scatter(list(Y[:, 0]), list(Y[:, 1]),
                                            labels=[cols[i] for i in range(Y.shape[0])],
                                            colour_by=[float(pt[rank[i]]) for i in range(len(cols))],
                                            xlabel="dimension 1", ylabel="dimension 2",
                                            title="pseudotime trajectory")
    return res_tbl


@T("monocle3_correlated_genes", "Genes correlated with pseudotime", MONO, "table",
   [tbl("counts", COUNTS, "Expression matrix (rows = genes, columns = cells)"),
    textbox("pseudotime_column", "", "Column holding pseudotime per cell"),
    choice("method", ["spearman", "pearson"], "spearman", "Correlation"),
    number("min_abs", 0.0, "Minimum |rho|", min=0.0), intin("head", 100, "Rows", min=1)],
   ex={"counts": COUNTS, "pseudotime_column": "sample_A", "min_abs": 0.0},
   up="monocle3 graph_test / tradeSeq", tags=("pseudotime", "correlation"),
   summary="Rank genes by how monotonically they change along a pseudotime vector.")
def monocle3_correlated_genes(counts, pseudotime_column="", method="spearman", min_abs=0.0, head=100):
    """Pseudotime-correlated genes."""
    names, cols, X = _mat(counts)
    if X.size == 0:
        return table([], "no numeric columns")
    if pseudotime_column and pseudotime_column in cols:
        pt = pd.to_numeric(_tbl(counts)[pseudotime_column], errors="coerce").to_numpy(dtype=float)
    else:
        pt = np.arange(len(cols), dtype=float)
    if pt.size != len(cols):
        pt = np.arange(len(cols), dtype=float)
    fn = stats.spearman if method == "spearman" else stats.pearson
    rows = []
    ps = []
    for i, nm in enumerate(names):
        v = X[i]
        m = np.isfinite(v) & np.isfinite(pt)
        if m.sum() < 4:
            continue
        st = dict(fn(list(pt[m]), list(v[m])))
        r = float(st.get("r", 0.0))
        p = float(st.get("p_value", 1.0)) if math.isfinite(float(st.get("p_value", 1.0))) else 1.0
        if abs(r) < float(min_abs):
            continue
        ps.append(p)
        rows.append({"gene": nm, "correlation": round(r, 5), "p_value": p,
                    "trend": "increasing" if r > 0 else "decreasing",
                    "max_abs": round(abs(r), 5), "direction": "positive" if r > 0 else "negative",
                    "pseudotime_range": f"{round(float(pt[m].min()), 3)}-{round(float(pt[m].max()), 3)}"})
    adj = stats.p_adjust(ps, "fdr_bh") if ps else []
    for i, r in enumerate(rows):
        r["p_adjusted"] = round(float(adj[i]), 6) if i < len(adj) else 1.0
    rows.sort(key=lambda r: -r["max_abs"])
    return table(rows[: int(head)], f"{len(rows)} genes correlated with "
                                 f"{pseudotime_column or 'cell rank'} ({method})")


@T("presto_aucell", "AUCell-like gene-set activity scores", PRESTO, "table",
   [tbl("counts", COUNTS, "Expression matrix"), txt("gene_sets", GMT, "GMT gene sets"),
    number("cut", 0.05, "Rank cut-off (fraction of genes)", min=0.001, max=0.5),
    boolean("rank_normalized", True, "Normalise ranks per cell"), intin("head", 100, "Rows", min=1)],
   ex={"counts": COUNTS, "gene_sets": GMT}, up="AUCell / presto:: aucell", tags=("gene sets", "scoring"),
   summary="Area under the recovery curve of each gene set within a cell's ranked genes.")
def presto_aucell(counts, gene_sets, cut=0.05, rank_normalized=True, head=100):
    """AUCell scoring."""
    names, cols, X = _mat(counts)
    sets = {k: {m.upper() for m in v} for k, v in io.parse_gmt(io.as_text(gene_sets)).items()}
    if X.size == 0 or not sets:
        return table([], "need an expression matrix and gene sets")
    M = X.astype(float)
    rows = []
    for j, c in enumerate(cols):
        v = M[:, j]
        order = np.argsort(-np.nan_to_num(v, nan=-np.inf))
        rank_of = {int(i): r for r, i in enumerate(order)}
        n = len(names)
        cut_n = max(1, int(round(n * float(cut))))
        for gs, members in sets.items():
            hits = sorted(rank_of[i] for i, nm in enumerate(names) if str(nm).upper() in members)
            if not hits:
                continue
            area = sum(1 for r in hits if r < cut_n) / float(cut_n)
            tot = len(hits) / float(n)
            rec = {"cell": c, "gene_set": gs, "auc": round(float(area), 5), "n_genes": len(hits),
                    "rank_cut": cut_n, "recovery": round(float(area), 5),
                    "set_fraction_in_cut": round((sum(1 for r in hits if r < cut_n) / max(1, len(hits))), 4)}
            if rank_normalized:
                rec["mean_rank_percentile"] = round(float(np.mean([rank_of[i] / max(1, n - 1)
                                                                 for i, nm in enumerate(names)
                                                                 if str(nm).upper() in members])), 4)
            rows.append(rec)
    rows.sort(key=lambda r: -r["auc"])
    return table(rows[: int(head)], f"{len(rows)} cell x gene-set scores (top set "
                                 f"{rows[0]['gene_set'] if rows else 'n/a'})")


@T("presto_marker_auc", "Per-cluster marker ranking by AUC", PRESTO, "table",
   [tbl("counts", COUNTS, "Expression matrix"), textbox("cluster_column", "", "Cluster labels in the metadata"),
    tbl("metadata", PHENO, "Metadata table"), number("min_auc", 0.7, "Minimum AUC", min=0.5, max=1.0),
    intin("head", 100, "Rows", min=1)],
   ex={"counts": COUNTS, "metadata": PHENO, "min_auc": 0.6}, up="presto::marker_genes / AUCell markers",
   tags=("marker genes", "clusters"),
   summary="Find genes that separate one cluster from the rest by area under the ROC curve.")
def presto_marker_auc(counts, cluster_column="", metadata=None, min_auc=0.7, head=100):
    """AUC marker genes."""
    names, cols, X = _mat(counts)
    if X.size == 0:
        return table([], "no numeric columns")
    md = _tbl(metadata) if metadata else pd.DataFrame()
    groups = {}
    if not md.empty and cluster_column and cluster_column in md.columns:
        labels = [str(v) for v in md[cluster_column].tolist()]
        pos = {str(c): i for i, c in enumerate(md[md.columns[0]].tolist())}
        for g in sorted(set(labels)):
            groups[g] = [pos[k] for k, v in zip(md[md.columns[0]].tolist(), labels) if v == g and k in pos]
    if not groups:
        half = len(cols) // 2
        groups = {"group_1": list(range(half)), "group_2": list(range(half, len(cols)))}
    rows = []
    for g, idxs in groups.items():
        other = [j for j in range(len(cols)) if j not in idxs]
        if not other:
            continue
        for i, nm in enumerate(names):
            a = X[i, idxs]
            b = X[i, other]
            scores = np.concatenate([a, b])
            truth = np.concatenate([np.ones(len(a)), np.zeros(len(b))])
            auc = float(stats.auc(list(scores), list(truth)))
            if auc < float(min_auc) and (1 - auc) < float(min_auc):
                continue
            rows.append({"gene": nm, "cluster": g, "auc": round(max(auc, 1 - auc), 5),
                        "specificity": "positive" if auc >= 0.5 else "negative",
                        "mean_in_cluster": round(float(np.nanmean(a)), 4),
                        "mean_outside": round(float(np.nanmean(b)), 4),
                        "fold_change": round(float(np.nanmean(a) / (np.nanmean(b) or 1e-9)), 4)})
    rows.sort(key=lambda r: -r["auc"])
    return table(rows[: int(head)], f"{len(rows)} candidate markers at AUC >= {min_auc}")


@T("seurat_summary_report", "Seurat-style summary of an object", SEURAT, "table",
   [tbl("counts", COUNTS, "Count matrix"), tbl("metadata", PHENO, "Metadata table"),
    textbox("cluster_column", "", "Cluster column"), intin("top_genes", 3, "Top genes per cluster", min=1)],
   ex={"counts": COUNTS, "metadata": PHENO, "cluster_column": "group"},
   up="Seurat object summary", tags=("single cell", "summary"),
   summary="Cells, genes, layers and per-cluster medians of an object in one table.")
def seurat_summary_report(counts, metadata, cluster_column="", top_genes=3):
    """Object summary."""
    names, cols, X = _mat(counts)
    md = _tbl(metadata) if metadata else pd.DataFrame()
    X = np.nan_to_num(X.astype(float), nan=0.0)
    rows = [{"field": "features", "value": len(names), "detail": ", ".join(names[:6])},
           {"field": "cells", "value": len(cols), "detail": ", ".join(cols[:6])},
           {"field": "total_counts", "value": int(X.sum()), "detail": "sum of all layers"},
           {"field": "median_genes_per_cell", "value": int(np.median([(X[:, j] > 0).sum() for j in range(X.shape[1])]))
            if X.size else 0, "detail": "detected features"},
           {"field": "metadata_rows", "value": int(len(md)), "detail": ",".join(str(c) for c in (md.columns[:6]
                                                            if not md.empty else []))},
           {"field": "layers", "value": 1, "detail": "data"}]
    lab = cluster_column if cluster_column in (md.columns if not md.empty else []) else ""
    if lab:
        groups = defaultdict(list)
        for v in md[lab].tolist():
            groups[str(v)].append(v)
        for g, items in sorted(groups.items()):
            rows.append({"field": f"cluster_{g}", "value": len(items),
                        "detail": f"{len(set(str(x) for x in items))} unique labels"})
    if X.size:
        med = np.median(X, axis=1)
        top = np.argsort(-med)[: int(top_genes)]
        rows.append({"field": "top_median_features", "value": len(top),
                    "detail": ", ".join(f"{names[i]}={round(float(med[i]), 2)}" for i in top)})
    return table(rows, f"object: {len(names)} features x {len(cols)} cells"
                   + (f", cluster column '{lab}'" if lab else ""))


@T("seurat_project_onto_pca", "Project new data onto existing loadings", SEURAT, "table",
   [tbl("reference", COUNTS, "Reference matrix"), tbl("query", COUNTS, "Query matrix"),
    intin("components", 2, "Components", min=1), boolean("standardize", True, "z-score using the reference"),
    intin("head", 100, "Rows", min=1)],
   ex={"reference": COUNTS, "query": COUNTS}, up="Seurat ProjectDim / Azimuth mapping",
   tags=("single cell", "projection", "integration"),
   summary="Score query cells on the reference PCA space without recomputing it.")
def seurat_project_onto_pca(reference, query, components=2, standardize=True, head=100):
    """PCA projection."""
    rnames, rcols, XR = _mat(reference)
    qnames, qcols, XQ = _mat(query)
    if XR.size == 0 or XQ.size == 0:
        return table([], "both matrices need numeric columns")
    R = XR.T.astype(float)
    Q = XQ.T.astype(float)
    mu = np.nanmean(R, axis=0)
    sd = np.nanstd(R, axis=0)
    if standardize:
        R = (R - mu) / np.where(sd > 0, sd, 1)
        Q = (Q - mu[: Q.shape[1]]) / np.where(sd[: Q.shape[1]] > 0, sd[: Q.shape[1]], 1)
    res = ml.pca(np.nan_to_num(R), components=min(int(components), R.shape[1]))
    V = np.atleast_2d(np.asarray(res["loadings"], dtype=float))
    if V.size == 0:
        return table([], "PCA failed on the reference")
    if V.shape[1] != R.shape[1] and V.shape[0] == R.shape[1]:
        V = V.T
    proj = np.nan_to_num(Q) @ V[: R.shape[1]].T if V.shape[0] >= R.shape[1] else np.nan_to_num(Q) @ V.T[: R.shape[1]]
    rows = []
    for i, c in enumerate(qcols):
        rec = {"cell": c}
        for j in range(proj.shape[1]):
            rec[f"PC{j + 1}"] = round(float(proj[i, j]), 5)
        rec["projection_distance"] = round(float(np.linalg.norm(proj[i])), 4)
        rows.append(rec)
    out = table(rows[: int(head)], f"projected {len(qcols)} query cells onto "
                                f"{proj.shape[1]} reference components")
    out["stats"]["figure"] = plot.scatter(list(proj[:, 0]), list(proj[:, 1]), labels=qcols,
                                        xlabel="reference PC1", ylabel="reference PC2", title="projection")
    return out
