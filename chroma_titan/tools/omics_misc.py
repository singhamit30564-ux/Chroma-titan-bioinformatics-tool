"""Metagenomics, microbiome, epigenetics, proteomics, metabolomics and cheminformatics.

Panel sections covered: **metagenomic_analysis**, **mothur**, **qiime2**,
**biodiversity_data_exploration**, **planttribes**, **virology**, **epigenetics**,
**str_fm__microsatellite_analysis**, **sccaf**, **du_novo**,
**transposon_insertion_sequencing**, **mimodd**, **gemini**, **proteomics**,
**metabolomics**, **chemicaltoolbox**, **pharmacology**, **multiomics**,
**import/manipulate_sc_data**, **hca-scanpy**, **monocle3**, **presto**, **seurat**,
**emboss**, **bbtools**.
"""
from __future__ import annotations

import math
import re
from collections import Counter, defaultdict

import numpy as np
import pandas as pd

from chroma_titan.core import align, genome, io, ml, motif, plot, protein, seq, stats, tables
from chroma_titan.tools._common import *  # noqa: F401,F403
from chroma_titan.tools._common import (ALN_SAM, ANNOT_GFF, BEDGRAPH, COUNTS, CXCG, GENES_FA, GENOME,
                                        GENOME_TXT, GMT, IMAGE, JASPAR, MSA, PAIRS, PEPTIDES,
                                        PHENO, PROTEINS, PFM, READS_1, READS_2, READS_SINGLE,
                                        REGIONS_BED, SHORT_DNA, SHORT_PROT, STR, TAXMAP, TREE,
                                        VARIANTS_VCF)

MO = "metagenomic_analysis"
MOTH = "mothur"
Q2 = "qiime2"
BDE = "biodiversity_data_exploration"
PT = "planttribes"
VIRO = "virology"
EPI = "epigenetics"
STRSEC = "str_fm__microsatellite_analysis"
SCCAF = "sccaf"
DUNOVO = "du_novo"
TN = "transposon_insertion_sequencing"
MIM = "mimodd"
GEM = "gemini"
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
EMB = "emboss"
BB = "bbtools"

RANKS = ["kingdom", "phylum", "class", "order", "family", "genus", "species"]


# ---------------------------------------------------------------------------
# shared helpers
# ---------------------------------------------------------------------------
def _fa(src):
    return io.parse_fasta(io.as_text(src))


def _bed(src):
    return io.parse_bed(io.as_text(src))


def _gff(src):
    return io.parse_gff(io.as_text(src))


def _vcf(src):
    return io.parse_vcf(io.as_text(src))


def _tbl(src):
    return tables.load(io.as_text(src))


def _mat(src):
    """Wide numeric table -> (row labels, column names, ndarray)."""
    df = _tbl(src)
    if df.empty:
        return [], [], np.zeros((0, 0))
    cols = list(df.columns)
    idc = cols[0]
    body = df.drop(columns=[idc])
    data = np.column_stack([pd.to_numeric(body[c], errors="coerce").to_numpy(dtype=float)
                           for c in body.columns]) if len(body.columns) else np.zeros((len(body), 0))
    return [str(x) for x in df[idc].tolist()], [str(c) for c in body.columns], data


def _num_cols(df):
    out = []
    for c in df.columns:
        v = pd.to_numeric(df[c], errors="coerce").to_numpy(dtype=float)
        if np.isfinite(v).sum() >= 2:
            out.append((str(c), v))
    return out


def _lineage_table(src, key_col="", lin_col=""):
    """Two-column (id, lineage) table -> {id: 'k__Bacteria;p__...'}.

    Falls back to the first and last columns, so it works on Mothur/QIIME
    taxonomy files as well as on ``examples/taxmap.tsv``.
    """
    df = _tbl(src)
    if df.empty:
        return {}
    cols = list(df.columns)
    kc = key_col if key_col in cols else cols[0]
    lc = lin_col if lin_col in cols else (cols[1] if len(cols) > 1 else cols[0])
    return {str(r[kc]): str(r[lc]) for r in df.to_dict("records") if r.get(kc) not in (None, "")}


def _rank_of(lineage: str, rank: str) -> str:
    """Return the label of ``rank`` in a ';'-separated lineage string."""
    if not lineage:
        return "unclassified"
    parts = [p.strip() for p in re.split(r"[;,]", lineage) if p.strip()]
    key = (rank or "genus").lower()
    prefix = {"kingdom": "k__", "phylum": "p__", "class": "c__", "order": "o__",
              "family": "f__", "genus": "g__", "species": "s__"}.get(key)
    if prefix:
        for p in parts:
            if p.lower().startswith(prefix):
                return p.split("__", 1)[-1].strip() or "unclassified"
    labels = [p.split("__", 1)[-1].strip() for p in parts if "__" in p] or parts
    if key.isdigit() and 0 <= int(key) < len(labels):
        return labels[int(key)]
    return labels[0] if labels else "unclassified"


def _rank_label(rank: str) -> str:
    key = (rank or "genus").lower()
    if key.isdigit():
        return RANKS[min(int(key), len(RANKS) - 1)]
    return key if key in RANKS else "genus"


def _cxcg(src):
    """Bismark / MethylDackel CXCG text -> list of per-cytosine dicts."""
    out = []
    for ln in io.as_text(src).splitlines():
        if not ln.strip() or ln.startswith("#"):
            continue
        p = ln.split("\t") if "\t" in ln else ln.split()
        if len(p) < 5:
            continue
        try:
            pos = int(float(p[1]))
        except ValueError:
            continue
        nums = []
        for tok in p[4:]:
            try:
                nums.append(float(tok))
            except ValueError:
                nums.append(float("nan"))
        frac = next((x for x in nums if math.isfinite(x) and 0.0 <= x <= 1.0), float("nan"))
        big = [x for x in nums[1:] if math.isfinite(x) and x > 1.0]
        methylated = unmethylated = float("nan")
        if big:
            if len(big) >= 2:
                methylated, total = big[-2], big[-1]
                if total >= methylated:
                    unmethylated = total - methylated
                else:
                    methylated, unmethylated = big[0], big[1]
            else:
                methylated = unmethylated = float("nan")
        if math.isnan(methylated) and math.isfinite(frac):
            methylated, unmethylated = round(frac * 100), round((1 - frac) * 100)
        ctx = next((t.upper() for t in p[2:5] if re.fullmatch(r"[CGH]{2,3}", t)), "CpG")
        if ctx == "CG":
            ctx = "CpG"
        out.append({"chrom": p[0], "pos": pos, "strand": p[2] if len(p) > 2 else ".",
                    "context": ctx, "fraction": frac, "methylated": methylated,
                    "unmethylated": unmethylated})
    return out


def _reads(src):
    return io.parse_fastq(io.as_text(src))


# ===========================================================================
# metagenomics / amplicon analysis
# ===========================================================================
@T("meta_count_table_summary", "OTU / feature table summary", MO, "table",
   [tbl("src", COUNTS, "Feature table (rows = features)"), boolean("per_sample", True, "One row per sample")],
   ex={"src": COUNTS}, up="qiime summarize table / phyloseq smry",
   tags=("abundance table", "summary"),
   summary="Library sizes, feature counts and singleton fractions of an abundance table.")
def meta_count_table_summary(src, per_sample=True):
    """Summarise an abundance table."""
    names, cols, X = _mat(src)
    if X.size == 0:
        return table([], "empty table")
    X = np.nan_to_num(X, nan=0.0)
    rows = []
    if per_sample:
        for j, c in enumerate(cols):
            v = X[:, j]
            rows.append({"sample": c, "total_reads": int(v.sum()), "features": int((v > 0).sum()),
                        "singletons": int((v == 1).sum()), "max": int(v.max()) if v.size else 0,
                        "mean_per_feature": round(float(v.mean()), 3),
                        "dominant_fraction": round(float(v.max() / v.sum()), 4) if v.sum() else 0.0})
    else:
        rows.append({"features": len(names), "samples": len(cols), "total_reads": int(X.sum()),
                    "mean_depth": round(float(X.sum() / max(1, X.shape[1])), 2),
                    "zero_cells": int((X == 0).sum())})
    return table(rows, f"{len(names)} features x {len(cols)} samples, {int(X.sum())} reads total")


@T("meta_alpha_diversity", "Within-sample diversity indices", MO, "table",
   [tbl("src", COUNTS, "Abundance table"), choice("index", ["shannon", "simpson", "inverse_simpson",
                                       "chao1", "gini", "pielou"], "shannon", "Diversity index"),
    intin("head", 50, "Subsample depth (0 = full)", min=0)],
   ex={"src": COUNTS, "index": "shannon"}, up="vegan::diversity / phyloseq estimate richness",
   tags=("alpha diversity", "ecology"),
   summary="Richness and evenness per sample, optionally rarefied to a common depth.")
def meta_alpha_diversity(src, index="shannon", head=50):
    """Alpha diversity per sample."""
    names, cols, X = _mat(src)
    if X.size == 0:
        return table([], "empty table")
    X = np.nan_to_num(X, nan=0.0)
    depth = int(head)
    rows = []
    for j, c in enumerate(cols):
        v = X[:, j]
        if depth and v.sum() > depth:
            v = np.asarray(stats.rarefied_counts(list(v), depth, seed=j), dtype=float)
        counts = [float(x) for x in v if x > 0]
        val = {"shannon": lambda: float(stats.shannon2(counts)),
               "simpson": lambda: float(stats.simpson(counts)),
               "inverse_simpson": lambda: float(stats.inverse_simpson(counts)),
               "chao1": lambda: float(stats.chao1(counts)),
               "gini": lambda: float(stats.gini(counts)),
               "pielou": lambda: float(stats.pielou_evenness(counts))}[index]()
        rows.append({"sample": c, "index": index, "value": round(val, 5),
                    "observed": len(counts), "total": int(v.sum()),
                    "shannon": round(float(stats.shannon2(counts)), 5),
                    "simpson": round(float(stats.simpson(counts)), 5),
                    "chao1": round(float(stats.chao1(counts)), 3),
                    "pielou_evenness": round(float(stats.pielou_evenness(counts)), 5)})
    rows.sort(key=lambda r: -r["value"])
    res = table(rows, f"{index} diversity for {len(rows)} samples")
    res["stats"]["figure"] = plot.bar([r["sample"] for r in rows], [r["value"] for r in rows],
                                     xlabel="sample", ylabel=index, title=f"{index} diversity", horizontal=True)
    return res


@T("meta_rarefaction", "Rarefaction curves of observed features", MO, "figure",
   [tbl("src", COUNTS, "Abundance table"), intin("steps", 12, "Curve points", min=3)],
   ex={"src": COUNTS, "steps": 8}, up="vegan::specpool / iNEXT",
   tags=("rarefaction", "diversity", "plot"),
   summary="Observed richness as a function of sequencing depth for every sample.")
def meta_rarefaction(src, steps=12):
    """Rarefaction curves."""
    _names, cols, X = _mat(src)
    if X.size == 0:
        return plot.empty_plot("empty table")
    X = np.nan_to_num(X, nan=0.0)
    series = {}
    for j, c in enumerate(cols[:12]):
        curve = stats.rarefaction(list(X[:, j]), steps=max(3, int(steps)), seed=j)
        pts = [dict(d) for d in curve]
        if pts:
            series[c] = [float(d.get("richness", 0.0)) for d in pts]
    if not series:
        return plot.empty_plot("no sequences to subsample")
    xs = [float(d.get("size", i)) for i, d in enumerate([dict(x) for x in stats.rarefaction(
        list(X[:, 0]), steps=max(3, int(steps)), seed=0)])]
    return plot.line(series, x=xs, xlabel="subsampled reads", ylabel="observed features",
                    title="Rarefaction curves")


@T("meta_beta_brays", "Bray-Curtis dissimilarity between samples", MO, "figure",
   [tbl("src", COUNTS, "Abundance table"), boolean("relativise", True, "Convert to relative abundance"),
    choice("metric", ["bray_curtis", "jaccard", "euclidean"], "bray_curtis", "Metric")],
   ex={"src": COUNTS, "metric": "bray_curtis"}, up="vegan::vegdist / phyloseq distance",
   tags=("beta diversity", "distance", "heatmap"),
   summary="Sample-by-sample community dissimilarity matrix as a heatmap.")
def meta_beta_brays(src, relativise=True, metric="bray_curtis"):
    """Community distance heatmap."""
    _names, cols, X = _mat(src)
    if X.shape[1] < 2:
        return plot.empty_plot("need at least two samples")
    X = np.nan_to_num(X.astype(float), nan=0.0)
    if relativise:
        tot = X.sum(axis=0)
        X = X / np.where(tot > 0, tot, 1.0)
    n = X.shape[1]
    D = np.zeros((n, n))
    for i in range(n):
        for j in range(i + 1, n):
            a, b = X[:, i], X[:, j]
            if metric == "jaccard":
                st = stats.jaccard_sets([k for k in range(X.shape[0]) if a[k] > 0],
                                       [k for k in range(X.shape[0]) if b[k] > 0])
                d = 1.0 - float(dict(st).get("jaccard", 0.0))
            elif metric == "euclidean":
                d = float(np.linalg.norm(a - b))
            else:
                d = float(np.sum(np.abs(a - b)) / max(1e-12, np.sum(a + b)))
            D[i, j] = D[j, i] = d
    return plot.heatmap([[round(float(v), 4) for v in row] for row in D], row_labels=cols, col_labels=cols,
                       title=f"{metric} dissimilarity", cmap="viridis", annot=True)


@T("meta_ordination_mds", "MDS / PCoA ordination of a distance matrix", MO, "figure",
   [tbl("src", COUNTS, "Abundance table"), choice("metric", ["bray_curtis", "euclidean"], "bray_curtis",
                                                 "Metric"), textbox("colour_column", "", "Colour by a table column")],
   ex={"src": COUNTS}, up="vegan::metaMDS / cmdscale", tags=("ordination", "PCoA", "plot"),
   summary="Two-dimensional ordination of community samples with a stress value.")
def meta_ordination_mds(src, metric="bray_curtis", colour_column=""):
    """PCoA-like ordination of samples."""
    _names, cols, X = _mat(src)
    if X.shape[1] < 3:
        return plot.empty_plot("need at least three samples")
    X = np.nan_to_num(X.astype(float), nan=0.0)
    tot = X.sum(axis=0)
    X = X / np.where(tot > 0, tot, 1.0) * 100.0
    n = X.shape[1]
    D = [[0.0] * n for _ in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            a, b = X[:, i], X[:, j]
            d = (float(np.sum(np.abs(a - b)) / max(1e-12, np.sum(a + b))) if metric == "bray_curtis"
                 else float(np.linalg.norm(a - b)))
            D[i][j] = D[j][i] = round(d, 6)
    res = ml.mds(D, components=2)
    E = np.asarray(res.get("embedding", []), dtype=float)
    if E.size == 0:
        return plot.empty_plot("ordination failed")
    colours = None
    if colour_column:
        df = _tbl(src)
        if colour_column in df.columns:
            cats = sorted({str(v) for v in df[colour_column].tolist()})
            mp = {c: i for i, c in enumerate(cats)}
            colours = [float(mp[str(v)]) for v in df[colour_column].tolist()][: E.shape[0]]
    return plot.scatter(list(E[:, 0]), list(E[:, 1]), labels=cols[: E.shape[0]], colour_by=colours,
                       xlabel=f"{metric} axis 1", ylabel=f"{metric} axis 2", title=f"PCoA ({metric})")


@T("meta_taxa_summary_by_rank", "Taxon counts at a lineage rank", MO, "table",
   [tbl("taxonomy", TAXMAP, "Taxonomy table (id, lineage)"),
    choice("rank", RANKS + ["0", "1", "2", "3", "4", "5", "6"], "genus", "Rank"),
    intin("head", 50, "Rows", min=1), boolean("bar", True, "Show a bar plot")],
   ex={"taxonomy": TAXMAP, "rank": "genus"}, up="kraken-report / phylosep tax_table",
   tags=("taxonomy", "abundance", "plot"),
   summary="How many sequences fall into each taxon at the chosen rank.")
def meta_taxa_summary_by_rank(taxonomy, rank="genus", head=50, bar=True):
    """Counts per taxon at one rank."""
    tax = _lineage_table(taxonomy)
    lab = _rank_label(rank)
    counts = Counter(_rank_of(v, lab) for v in tax.values())
    tot = sum(counts.values()) or 1
    rows = [{"taxon": t, "rank": lab, "n_sequences": c, "percent": round(100 * c / tot, 3)}
            for t, c in counts.most_common()]
    out = table(rows[: int(head)], f"{len(rows)} {lab} taxa over {tot} sequences")
    if bar and rows:
        top = rows[:12][::-1]
        out["stats"]["figure"] = plot.bar([r["taxon"] for r in top], [r["n_sequences"] for r in top],
                                        xlabel=lab, ylabel="sequences", title=f"Top {lab} taxa",
                                        horizontal=True)
    return out


@T("meta_filter_by_lineage", "Keep or drop sequences by taxonomic lineage", MO, "text",
   [fa("records", GENES_FA, "FASTA sequences"), tbl("taxonomy", TAXMAP, "Taxonomy table"),
    textbox("pattern", "Proteobacteria", "Lineage substring"), choice("mode", ["keep", "remove"], "keep",
                                                                     "Mode")],
   ex={"records": GENES_FA, "taxonomy": TAXMAP, "pattern": "Proteobacteria"},
   up="qiime taxa filter-table / phyloseq prune", tags=("taxonomy", "filter", "fasta"),
   summary="Subset a FASTA by matching lineage text, reporting what was kept.")
def meta_filter_by_lineage(records, taxonomy, pattern="Proteobacteria", mode="keep"):
    """Filter sequences on their lineage string."""
    tax = _lineage_table(taxonomy)
    recs = _fa(records)
    keep = []
    dropped = 0
    for r in recs:
        lin = tax.get(r.id, "")
        hit = bool(pattern) and pattern.lower() in lin.lower()
        if (hit and mode == "keep") or (not hit and mode == "remove"):
            keep.append(r)
        else:
            dropped += 1
    body = io.write_fasta(keep) if keep else ""
    return text(body, f"{len(keep)} of {len(recs)} records kept ({mode} on '{pattern}', {dropped} dropped)")


@T("meta_abundance_at_rank", "Abundance table aggregated to a taxonomic rank", MO, "table",
   [tbl("src", COUNTS, "Feature table"), tbl("taxonomy", TAXMAP, "Taxonomy table"),
    choice("rank", RANKS, "genus", "Rank"), boolean("relative", True, "Relative abundance")],
   ex={"src": COUNTS, "taxonomy": TAXMAP, "rank": "genus"},
   up="phylosef tax_glom / qiime agglomerate", tags=("taxonomy", "aggregation"),
   summary="Sum sample abundances per taxon at a chosen rank.")
def meta_abundance_at_rank(src, taxonomy, rank="genus", relative=True):
    """Glance features up to a rank."""
    names, cols, X = _mat(src)
    tax = _lineage_table(taxonomy)
    lab = _rank_label(rank)
    groups: dict[str, list[int]] = defaultdict(list)
    for i, nm in enumerate(names):
        groups[_rank_of(tax.get(nm, ""), lab) if tax else nm].append(i)
    rows = []
    for g, idxs in groups.items():
        sub = X[idxs].sum(axis=0) if idxs else np.zeros(len(cols))
        tot = float(X.sum()) or 1.0
        rec = {"taxon": g, "features": len(idxs), "total": int(np.nansum(sub))}
        for j, c in enumerate(cols):
            v = float(sub[j]) if j < len(sub) and math.isfinite(sub[j]) else 0.0
            rec[c] = round(v / tot, 6) if relative else int(v)
        rows.append(rec)
    rows.sort(key=lambda r: -r["total"])
    return table(rows, f"{len(rows)} {lab} taxa after agglomeration")


@T("meta_lca_of_hits", "Lowest common ancestor from BLAST-like hits", MO, "table",
   [fa("query_records", GENES_FA, "Query FASTA"), fa("subject_records", PROTEINS, "Subject FASTA"),
    tbl("taxonomy", TAXMAP, "Taxonomy table"), intin("max_hits", 5, "Hits per query", min=1),
    choice("scoring", ["blastn", "blastp"], "blastn", "Scoring")],
   ex={"query_records": GENES_FA, "subject_records": GENES_FA, "taxonomy": TAXMAP, "max_hits": 3},
   up="MEGAN LCA / kraken2 classify", tags=("taxonomy", "LCA", "classification"),
   summary="Assign each query to the lowest common ancestor of its best hits.")
def meta_lca_of_hits(query_records, subject_records, taxonomy, max_hits=5, scoring="blastn"):
    """LCA classification of BLAST-like hits."""
    q = _fa(query_records)
    subs = [{"id": r.id, "seq": r.seq} for r in _fa(subject_records)]
    hits = align.blast_lite([{"id": r.id, "seq": r.seq} for r in q], subs, scoring=scoring,
                           max_hits=int(max_hits) * max(1, len(subs)))
    tax = _lineage_table(taxonomy)
    lin = {k: [p for p in re.split(r"[;,]", v) if p.strip()] for k, v in tax.items()}
    lca = align.lca_of_hits(hits, lin) if lin else []
    rows = []
    for h in lca:
        d = dict(h)
        rows.append({"query": d.get("q", d.get("query", "")), "taxon": d.get("lca", d.get("taxon", "")),
                    "n_hits": d.get("n_hits", d.get("count", "")),
                    "identity": round(float(d.get("pident", d.get("identity", 0.0)) or 0.0), 3)})
    if not rows:
        by_q = defaultdict(list)
        for h in hits:
            d = dict(h)
            by_q[d.get("q", d.get("query", ""))].append(d)
        for qname, hs in by_q.items():
            hsq = sorted(hs, key=lambda d: -float(d.get("pident", 0.0)))[: int(max_hits)]
            lines = [lin.get(d.get("s", d.get("subject", "")), []) for d in hsq]
            lines = [x for x in lines if x]
            common = []
            for parts in zip(*lines) if lines else []:
                if len(set(parts)) == 1:
                    common.append(parts[0])
                else:
                    break
            rows.append({"query": qname, "taxon": ";".join(x.split("__", 1)[-1] for x in common) or "unclassified",
                        "n_hits": len(hsq),
                        "identity": round(float(hsq[0].get("pident", 0.0)), 3) if hsq else 0.0})
    return table(rows[:200], f"{len(rows)} queries classified by LCA")


@T("meta_denoise_asvs", "Denoise sequences by merging near-identical variants", MO, "text",
   [fa("records", GENES_FA, "Amplicon FASTA"), intin("mismatches", 1, "Allowed differences", min=0),
    intin("min_size", 1, "Minimum cluster size", min=1)],
   ex={"records": GENES_FA, "mismatches": 1}, up="DADA2 / deblur / mothur pre.cluster",
   tags=("denoising", "ASV", "clustering"),
   summary="Collapse sequences that differ by a few substitutions into one ASV.")
def meta_denoise_asvs(records, mismatches=1, min_size=1):
    """Single-linkage denoising on Hamming distance."""
    recs = [r for r in _fa(records) if r.seq]
    used = [False] * len(recs)
    order = sorted(range(len(recs)), key=lambda i: -len(recs[i].seq))
    out, clusters = [], []
    for i in order:
        if used[i]:
            continue
        used[i] = True
        members = [i]
        for j in order:
            if used[j] or len(recs[j].seq) != len(recs[i].seq):
                continue
            hd = align.hamming(recs[i].seq, recs[j].seq)
            if int(dict(hd).get("distance", 99)) <= int(mismatches):
                used[j] = True
                members.append(j)
        if len(members) >= int(min_size):
            clusters.append((i, members))
    for i, members in clusters:
        for k in members:
            out.append(io.Seq(id=recs[i].id if k == i else f"{recs[i].id}_v{k}",
                             seq=recs[k].seq, desc=f"asv={i};size={len(members)}"))
    body = io.write_fasta(out)
    return text(body, f"{len(clusters)} ASVs from {len(recs)} sequences "
                     f"(largest cluster {max((len(m) for _i, m in clusters), default=0)})")


@T("meta_chimera_screen", "Flag chimeric sequences", MO, "table",
   [fa("records", GENES_FA, "Sequences (first = parent set)"), intin("min_gain", 2, "Minimum split advantage",
                                             min=1)],
   ex={"records": GENES_FA}, up="vsearch --uchime-denovo / mothur chimera",
   tags=("chimera", "quality control"),
   summary="Detect sequences whose halves match two different parents better than themselves.")
def meta_chimera_screen(records, min_gain=2):
    """Simple chimera detection."""
    recs = _fa(records)
    if len(recs) < 2:
        return table([], "need at least two sequences")
    rows = []
    for i, r in enumerate(recs):
        s = r.seq.upper()
        if len(s) < 20:
            continue
        mid = len(s) // 2
        left, right = s[:mid], s[mid:]
        best_l = best_r = (-1, 0)
        for j, p in enumerate(recs):
            if j == i:
                continue
            ps = p.seq.upper()
            dl = int(dict(align.hamming(left, ps[:len(left)] if len(ps) >= len(left) else
                                       ps[-len(left):])).get("distance", 99)) if len(ps) >= len(left) else 99
            dr = int(dict(align.hamming(right, ps[:len(right)] if len(ps) >= len(right) else
                                       ps[-len(right):])).get("distance", 99)) if len(ps) >= len(right) else 99
            if dl > best_l[1]:
                best_l = (j, dl)
            if dr > best_r[1]:
                best_r = (j, dr)
        own = 0
        chimera = (best_l[0] != best_r[0] and best_l[0] >= 0 and best_r[0] >= 0
                   and best_l[1] + best_r[1] <= max(2, int(min_gain)))
        rows.append({"sequence": r.id, "length": len(s), "left_parent": recs[best_l[0]].id if best_l[0] >= 0 else "",
                    "left_mismatches": best_l[1], "right_parent": recs[best_r[0]].id if best_r[0] >= 0 else "",
                    "right_mismatches": best_r[1], "chimera_flag": bool(chimera)})
    n = sum(1 for r in rows if r["chimera_flag"])
    return table(rows[:200], f"{n} of {len(rows)} sequences flagged as chimeric")


@T("meta_prevalence_filter", "Remove features seen only in negative controls", MO, "table",
   [tbl("src", COUNTS, "Feature table"), textbox("blank_columns", "", "Blank / control columns"),
    number("min_fold", 2.0, "Minimum enrichment in real samples", min=0.0),
    choice("action", ["report", "keep"], "report", "Action")],
   ex={"src": COUNTS, "blank_columns": "sample_A,sample_B"}, up="decontam prevalence",
   tags=("contamination", "filter"),
   summary="Compare feature prevalence between control and real samples with Fisher's test.")
def meta_prevalence_filter(src, blank_columns="", min_fold=2.0, action="report"):
    """Prevalence-based decontamination."""
    names, cols, X = _mat(src)
    blanks = [cols.index(c) for c in str(blank_columns).split(",") if c in cols] or cols[: max(1, len(cols) // 3)]
    real = [j for j in range(len(cols)) if j not in blanks]
    rows = []
    keep_idx = []
    for i, nm in enumerate(names):
        bp = int(sum(1 for j in blanks if X[i, j] > 0))
        rp = int(sum(1 for j in real if X[i, j] > 0))
        bmean = float(np.nanmean([X[i, j] for j in blanks])) if blanks else 0.0
        rmean = float(np.nanmean([X[i, j] for j in real])) if real else 0.0
        fold = rmean / bmean if bmean > 0 else (math.inf if rmean > 0 else 1.0)
        orr, p = stats.fisher_exact(rp, len(real) - rp, bp, max(0, len(blanks) - bp))
        suspect = bool(p < 0.05 and fold < float(min_fold))
        rows.append({"feature": nm, "prevalence_blank": round(bp / max(1, len(blanks)), 3),
                    "prevalence_real": round(rp / max(1, len(real)), 3),
                    "mean_blank": round(bmean, 3), "mean_real": round(rmean, 3),
                    "fold": round(fold, 4) if math.isfinite(fold) else "inf",
                    "p_value": round(float(p), 6), "contaminant_flag": suspect})
        if action == "keep" and not suspect:
            keep_idx.append(i)
    out = table(rows[:300], f"{sum(1 for r in rows if r['contaminant_flag'])} of {len(rows)} features "
                          f"look like contaminants")
    if action == "keep":
        out["stats"]["kept_features"] = [names[i] for i in keep_idx]
    return out


@T("meta_kraken_report", "Kraken-style hierarchical report", MO, "table",
   [tbl("taxonomy", TAXMAP, "Taxonomy table"), intin("min_percent", 0.0, "Minimum percent", min=0.0),
    boolean("cumulative", True, "Cumulative counts down the lineage")],
   ex={"taxonomy": TAXMAP, "min_percent": 0.0}, up="kraken-report / pavian",
   tags=("taxonomy", "report", "metagenomics"),
   summary="Per-taxon percentages at every rank, like a kraken report.", version="1.0.0")
def meta_kraken_report(taxonomy, min_percent=0.0, cumulative=True):
    """Layered taxonomic report."""
    tax = _lineage_table(taxonomy)
    tot = max(1, len(tax))
    rows = []
    for depth, rk in enumerate(RANKS):
        counts = Counter(_rank_of(v, rk) for v in tax.values())
        for t, c in counts.most_common():
            pct = 100.0 * c / tot
            if pct < float(min_percent):
                continue
            parent = _rank_of(next(v for v in tax.values() if t in v), RANKS[max(0, depth - 1)]) if depth else ""
            rows.append({"depth": depth + 1, "rank": rk, "taxon": t, "n_sequences": c,
                        "percent": round(pct, 3), "parent": parent,
                        "cumulative": int(c) if cumulative else int(c)})
    rows.sort(key=lambda r: (r["depth"], -r["n_sequences"]))
    return table(rows[:400], f"{len(rows)} report lines over {tot} sequences")


# ===========================================================================
# Mothur-style 16S tools
# ===========================================================================
@T("mothur_get_current", "Unique sequences with abundances", MOTH, "table",
   [fa("records", GENES_FA, "FASTA file"), boolean("rename", True, "Rename to a/b/c ids"),
    intin("head", 200, "Rows", min=1)],
   ex={"records": GENES_FA}, up="mothur get.current / count.seqs", tags=("dereplication", "16S"),
   summary="Dereplicate a FASTA and report how often each unique sequence occurs.")
def mothur_get_current(records, rename=True, head=200):
    """Dereplication table."""
    seen: dict[str, list[str]] = defaultdict(list)
    for r in _fa(records):
        seen[r.seq.upper()].append(r.id)
    rows = [{"sequence": (chr(ord("a") + i) if rename and i < 26 else f"uniq_{i + 1}"),
            "abundance": len(ids), "length": len(s), "members": ",".join(ids[:6]),
            "percent": round(100 * len(ids) / max(1, sum(len(v) for v in seen.values())), 3)}
           for i, (s, ids) in enumerate(sorted(seen.items(), key=lambda kv: (-len(kv[1]), kv[0])))]
    return table(rows[: int(head)], f"{len(rows)} unique sequences from "
                                  f"{sum(len(v) for v in seen.values())} input records")


@T("mothur_count_seqs", "Sequence length and composition summary", MOTH, "table",
   [fa("records", GENES_FA, "FASTA file"), intin("bin", 100, "Length bin size", min=1)],
   ex={"records": GENES_FA, "bin": 200}, up="mothur count.seqs", tags=("summary", "length"),
   summary="Length distribution, GC and ambiguity statistics of a nucleotide library.")
def mothur_count_seqs(records, bin=100):
    """Library summary."""
    recs = _fa(records)
    if not recs:
        return table([], "no sequences")
    L = [len(r.seq) for r in recs]
    gc = [round(100 * seq.gc_content(r.seq), 3) for r in recs]
    nbins = max(1, int((max(L) - min(L)) // max(1, int(bin))) + 1)
    hist = {}
    for x in L:
        hist[f"{(x // int(bin)) * int(bin)}"] = hist.get(f"{(x // int(bin)) * int(bin)}", 0) + 1
    amb = sum(1 for r in recs if re.search(r"[^ACGTacgt]", r.seq))
    res = table([{"metric": "sequences", "value": len(recs)},
                {"metric": "bases", "value": int(sum(L))},
                {"metric": "min_length", "value": int(min(L))}, {"metric": "max_length", "value": int(max(L))},
                {"metric": "mean_length", "value": round(float(np.mean(L)), 2)},
                {"metric": "median_length", "value": float(np.median(L))},
                {"metric": "mean_gc_percent", "value": round(float(np.mean(gc)), 3)},
                {"metric": "sequences_with_ambiguity", "value": amb},
                {"metric": "uniformity", "value": round(float(np.std(L) / (np.mean(L) or 1)), 4)}],
               f"{len(recs)} sequences, mean length {round(float(np.mean(L)), 1)} bp")
    res["stats"]["length_histogram"] = hist
    return res


@T("mothur_classify_seqs", "Classify sequences against a reference taxonomy", MOTH, "table",
   [fa("records", GENES_FA, "Query FASTA"), fa("reference", GENES_FA, "Reference FASTA"),
    tbl("taxonomy", TAXMAP, "Reference taxonomy"), number("cutoff", 60.0, "Minimum identity %", min=0.0),
    intin("head", 100, "Rows", min=1)],
   ex={"records": GENES_FA, "reference": GENES_FA, "taxonomy": TAXMAP}, up="mothur classify.seqs",
   tags=("classification", "taxonomy"),
   summary="Best reference match per query with identity and the assigned lineage.")
def mothur_classify_seqs(records, reference, taxonomy, cutoff=60.0, head=100):
    """Nearest-neighbour classification."""
    refs = _fa(reference)
    tax = _lineage_table(taxonomy)
    rows = []
    for r in _fa(records):
        best, bid = "", -1.0
        for p in refs:
            if p.id == r.id or not p.seq:
                continue
            hd = dict(align.hamming(r.seq[: len(p.seq)], p.seq[: len(r.seq)]))
            alen = max(1, min(len(r.seq), len(p.seq)))
            ident = 100.0 * (alen - int(hd.get("distance", alen))) / alen
            if ident > bid:
                bid, best = ident, p.id
        lin = tax.get(best, "")
        rows.append({"query": r.id, "reference": best, "identity_percent": round(bid, 3),
                    "lineage": lin or "unclassified",
                    "confidence": round(bid / 100.0, 4),
                    "classified": bool(bid >= float(cutoff))})
    n = sum(1 for r in rows if r["classified"])
    return table(rows[: int(head)], f"{n}/{len(rows)} queries classified above {cutoff}% identity")


@T("mothur_dist_seqs", "Pairwise sequence distances", MOTH, "table",
   [fa("records", GENES_FA, "FASTA file"), choice("model", ["identity", "p", "k2p", "raw"], "identity",
                                                 "Distance model"), intin("head", 60, "Pair rows", min=1)],
   ex={"records": GENES_FA, "model": "p"}, up="mothur dist.seqs / phylip dist",
   tags=("distance matrix", "pairwise"),
   summary="Every pair of sequences scored by an evolutionary or identity distance.")
def mothur_dist_seqs(records, model="identity", head=60):
    """Pairwise distances."""
    recs = _fa(records)
    if len(recs) < 2:
        return table([], "need at least two sequences")
    D = align.distance_matrix([r.seq for r in recs], model={"identity": "identity", "p": "p",
                                                           "k2p": "k2p", "raw": "identity"}[model])
    rows = []
    for i in range(len(recs)):
        for j in range(i + 1, len(recs)):
            rows.append({"sequence_a": recs[i].id, "sequence_b": recs[j].id,
                        "distance": round(float(D[i][j]), 6), "model": model})
    rows.sort(key=lambda r: r["distance"])
    d = [r["distance"] for r in rows]
    res = table(rows[: int(head)], f"{len(rows)} pairs, mean distance {round(float(np.mean(d)), 5)}, "
                                 f"range {round(min(d), 5)}-{round(max(d), 5)}")
    res["stats"]["matrix"] = [[round(float(v), 5) for v in row] for row in D]
    res["stats"]["labels"] = [r.id for r in recs]
    return res


@T("mothur_trim_seqs", "Trim sequences to a start and end position", MOTH, "text",
   [fa("records", GENES_FA, "FASTA file"), intin("start", 1, "Start position", min=1),
    intin("end", 300, "End position", min=1), boolean("left_justify", False, "Trim from the left")],
   ex={"records": GENES_FA, "start": 1, "end": 200}, up="mothur trim.seqs", tags=("trimming", "16S"),
   summary="Cut every sequence to the same alignment window for comparable regions.")
def mothur_trim_seqs(records, start=1, end=300, left_justify=False):
    """Positional trimming."""
    s0, e0 = int(start), int(end)
    out = []
    for r in _fa(records):
        sq = r.seq
        if left_justify:
            i = 0
            while i < len(sq) and sq[i] in "-.":
                i += 1
            sq = sq[i:]
        piece = sq[max(0, s0 - 1): max(0, e0)]
        if piece:
            out.append(io.Seq(id=r.id, seq=piece, desc=r.desc))
    return text(io.write_fasta(out), f"trimmed {len(out)} sequences to {s0}-{e0} "
                                    f"(lengths {min((len(r.seq) for r in out), default=0)}-"
                                    f"{max((len(r.seq) for r in out), default=0)})")


@T("mothur_screen_seqs", "Screen sequences by length and ambiguity", MOTH, "text",
   [fa("records", GENES_FA, "FASTA file"), intin("min_length", 10, "Minimum length", min=1),
    intin("max_length", 100000, "Maximum length", min=1), intin("max_ambiguous", 0, "Maximum N bases", min=0),
    textbox("patterns", "", "Sequence motifs to exclude (comma separated)")],
   ex={"records": GENES_FA, "min_length": 50, "max_ambiguous": 0}, up="mothur screen.seqs",
   tags=("quality control", "filter"),
   summary="Drop sequences that are too short, too long, ambiguous or contain bad motifs.")
def mothur_screen_seqs(records, min_length=10, max_length=100000, max_ambiguous=0, patterns=""):
    """Quality screening."""
    pats = [p.strip() for p in str(patterns).split(",") if p.strip()]
    keep, drop = [], Counter()
    for r in _fa(records):
        if len(r.seq) < int(min_length):
            drop["too_short"] += 1
            continue
        if len(r.seq) > int(max_length):
            drop["too_long"] += 1
            continue
        if r.seq.upper().count("N") > int(max_ambiguous):
            drop["ambiguous"] += 1
            continue
        if any(p.upper() in r.seq.upper() for p in pats):
            drop["contains_pattern"] += 1
            continue
        keep.append(r)
    return text(io.write_fasta(keep), f"{len(keep)} sequences kept, {sum(drop.values())} removed "
                                     f"({dict(drop)})")


@T("mothur_pre_cluster", "Cluster sequences at a distance cutoff", MOTH, "table",
   [fa("records", GENES_FA, "FASTA file"), number("distance", 0.03, "Maximum distance", min=0.0),
    choice("method", ["average", "maximum", "nearest"], "average", "Linkage"), intin("head", 200, "Rows", min=1)],
   ex={"records": GENES_FA, "distance": 0.1}, up="mothur pre.cluster / vsearch cluster",
   tags=("clustering", "OTU"),
   summary="Greedy agglomerative clustering of sequences into OTU-like groups.")
def mothur_pre_cluster(records, distance=0.03, method="average", head=200):
    """Cluster at a distance threshold."""
    recs = _fa(records)
    if not recs:
        return table([], "no sequences")
    D = align.distance_matrix([r.seq for r in recs], model="identity")
    labels = [-1] * len(recs)
    nxt = 0
    for i in range(len(recs)):
        if labels[i] >= 0:
            continue
        labels[i] = nxt
        for j in range(i + 1, len(recs)):
            if labels[j] >= 0:
                continue
            if float(D[i][j]) <= float(distance):
                labels[j] = nxt
        nxt += 1
    sizes = Counter(labels)
    reps = {}
    rows = [{"cluster": chr(ord("A") + k) if k < 26 else f"OBU{k + 1}", "size": v,
            "representative": reps.setdefault(k, next(recs[i].id for i in range(len(recs)) if labels[i] == k)),
            "max_distance_within": round(max([float(D[i][j]) for i in range(len(recs)) for j in range(i + 1, len(recs))
                                            if labels[i] == labels[j] == k] or [0.0]), 5),
            "members": ",".join(recs[i].id for i in range(len(recs)) if labels[i] == k)[:80]}
           for k, v in sorted(sizes.items())]
    return table(rows[: int(head)], f"{len(rows)} clusters from {len(recs)} sequences at d <= {distance}")


# ===========================================================================
# QIIME 2 style workflows
# ===========================================================================
@T("qiime_dereplicate", "Dereplicate and count feature sequences", Q2, "table",
   [fa("records", GENES_FA, "Feature FASTA"), boolean("relabelling", True, "Renumber feature ids")],
   ex={"records": GENES_FA}, up="qiime quality-control dereplicate", tags=("dereplication", "features"),
   summary="Collapse identical sequences into features with a size column.")
def qiime_dereplicate(records, relabelling=True):
    """Dereplicate features."""
    groups: dict[str, list[str]] = defaultdict(list)
    for r in _fa(records):
        groups[r.seq.upper()].append(r.id)
    rows = [{"feature_id": f"feature_{i + 1}" if relabelling else ids[0], "size": len(ids),
            "length": len(s), "original_ids": ",".join(ids[:5])}
           for i, (s, ids) in enumerate(sorted(groups.items(), key=lambda kv: (-len(kv[1]), kv[0])))]
    return table(rows, f"{len(rows)} dereplicated features from {sum(len(v) for v in groups.values())} sequences")


@T("qiime_fastq_stats", "Per-length and per-quality FASTQ summary", Q2, "table",
   [fq("reads", READS_SINGLE, "FASTQ reads"), intin("quality_window", 4, "Report first n bases individually",
                                              min=1)],
   ex={"reads": READS_1}, up="qiime demux summarize / fastp", tags=("quality control", "reads"),
   summary="Read counts, length profile and mean quality of an amplicon library.")
def qiime_fastq_stats(reads, quality_window=4):
    """Demultiplex-style QC summary."""
    recs = _reads(reads)
    if not recs:
        return table([], "no reads parsed")
    L = np.array([len(r.seq) for r in recs], dtype=float)

    def qual(r):
        q = [max(0, ord(c) - 33) for c in (r.qual or "")]
        return float(np.mean(q)) if q else 0.0
    Q = np.array([qual(r) for r in recs], dtype=float)
    rows = [{"metric": "reads", "value": len(recs)},
            {"metric": "bases", "value": int(L.sum())},
            {"metric": "mean_length", "value": round(float(L.mean()), 2)},
            {"metric": "mean_quality", "value": round(float(Q.mean()), 3)},
            {"metric": "reads_below_q20", "value": int(np.sum(Q < 20))},
            {"metric": "duplicate_rate", "value": round(1 - len({r.seq.upper() for r in recs}) / len(recs), 4)},
            {"metric": "min_length", "value": int(L.min())}, {"metric": "max_length", "value": int(L.max())}]
    for i in range(min(int(quality_window), int(L.min()))):
        q = np.mean([ord(r.qual[i]) - 33 for r in recs if len(r.qual) > i])
        n = np.mean([r.seq[i].upper() == "N" for r in recs if len(r.seq) > i])
        rows.append({"metric": f"cycle_{i + 1}_mean_quality", "value": round(float(q), 3)})
        rows.append({"metric": f"cycle_{i + 1}_n_rate", "value": round(float(n), 5)})
    return table(rows, f"{len(recs)} reads, mean quality {round(float(Q.mean()), 2)}")


@T("qiime_trim_primers", "Cut primer sequences from reads", Q2, "table",
   [fq("reads", READS_1, "FASTQ reads"), textbox("forward_primer", "", "Forward primer"),
    textbox("reverse_primer", "", "Reverse primer"), boolean("remove_primers", True, "Cut primers off"),
    choice("failure", ["keep", "discard"], "discard", "If the primer is absent")],
   ex={"reads": READS_1, "forward_primer": "ACG", "reverse_primer": "CGT"}, up="qiime cutadapt trim-pairs",
   tags=("primers", "trimming", "reads"),
   summary="Locate and remove amplicon primers, reporting how many reads were affected.")
def qiime_trim_primers(reads, forward_primer="", reverse_primer="", remove_primers=True, failure="discard"):
    """Primer trimming."""
    recs = _reads(reads)
    rows = []
    kept = 0
    for r in recs:
        s = r.seq.upper()
        i = s.find(forward_primer.upper()) if forward_primer else 0
        j = s.rfind(seq.revcomp(reverse_primer).upper()) if reverse_primer else len(s)
        found = bool(forward_primer) and i >= 0
        if (forward_primer and i < 0) or (reverse_primer and j < 0):
            rows.append({"read": r.id, "status": "primer_not_found", "start": -1, "end": -1,
                        "trimmed_length": 0})
            if failure == "discard":
                continue
        start = i + len(forward_primer) if (remove_primers and forward_primer and i >= 0) else max(0, i)
        end = j if (remove_primers and reverse_primer and j >= 0) else (j + len(reverse_primer)
                                                                       if j >= 0 else len(s))
        end = max(start, min(end, len(s)))
        kept += 1
        rows.append({"read": r.id, "status": "trimmed" if found else "kept",
                    "start": start, "end": end, "trimmed_length": end - start})
    return table(rows[:200], f"{kept}/{len(recs)} reads usable after primer removal")


@T("qiime_filter_features", "Filter a feature table by prevalence and depth", Q2, "table",
   [tbl("src", COUNTS, "Feature table"), number("min_prevalence", 0.5, "Minimum fraction of samples", min=0.0),
    intin("min_total", 10, "Minimum total count", min=0), intin("min_per_sample", 0, "Minimum count per sample",
                                                              min=0)],
   ex={"src": COUNTS, "min_prevalence": 0.5, "min_total": 100}, up="qiime feature-table filter-features",
   tags=("filter", "table"),
   summary="Drop rare or low-abundance features and report what survives.")
def qiime_filter_features(src, min_prevalence=0.5, min_total=10, min_per_sample=0):
    """Feature filtering."""
    names, cols, X = _mat(src)
    X = np.nan_to_num(X.astype(float), nan=0.0)
    rows, keep = [], []
    for i, nm in enumerate(names):
        v = X[i]
        prev = float(np.mean(v > 0))
        tot = float(np.nansum(v))
        ok = prev >= float(min_prevalence) and tot >= float(min_total) and (not min_per_sample
                                                                           or float(np.min(v)) >= int(min_per_sample))
        rows.append({"feature": nm, "prevalence": round(prev, 3), "total": int(tot),
                    "min_in_sample": int(np.min(v)) if v.size else 0, "kept": bool(ok)})
        if ok:
            keep.append(i)
    res = table(rows, f"{len(keep)}/{len(names)} features pass (prevalence >= {min_prevalence}, "
                    f"total >= {min_total})")
    res["stats"]["kept"] = [names[i] for i in keep]
    return res


@T("qiime_taxa_barplot", "Stacked bar plot of taxon composition", Q2, "table",
   [tbl("src", COUNTS, "Feature table"), tbl("taxonomy", TAXMAP, "Taxonomy table"),
    choice("rank", RANKS, "phylum", "Rank"), intin("top", 5, "Number of taxa to show", min=2)],
   ex={"src": COUNTS, "taxonomy": TAXMAP, "rank": "phylum", "top": 4}, up="qiime taxa barplot",
   tags=("barplot", "composition", "plot"),
   summary="Relative composition of every sample, grouped at a taxonomic rank.")
def qiime_taxa_barplot(src, taxonomy, rank="phylum", top=5):
    """Composition table + grouped bars."""
    names, cols, X = _mat(src)
    tax = _lineage_table(taxonomy)
    lab = _rank_label(rank)
    X = np.nan_to_num(X.astype(float), nan=0.0)
    per_taxon: dict[str, np.ndarray] = {}
    for i, nm in enumerate(names):
        t = _rank_of(tax.get(nm, ""), lab) if tax else nm
        per_taxon[t] = per_taxon.get(t, np.zeros(len(cols))) + X[i]
    order = sorted(per_taxon, key=lambda t: -float(per_taxon[t].sum()))
    show = order[: int(top)]
    other = np.zeros(len(cols))
    for t in order[int(top):]:
        other += per_taxon[t]
    series = {}
    rows = []
    for t in show + (["other"] if other.sum() else []):
        v = other if t == "other" else per_taxon[t]
        for j, c in enumerate(cols):
            tot = float(X[:, j].sum()) or 1.0
            rec = {"sample": c, "taxon": t, "relative_abundance": round(float(v[j]) / tot, 5),
                   "counts": int(v[j])}
            rows.append(rec)
        series[t] = [float(v[j]) / (float(X[:, j].sum()) or 1.0) for j in range(len(cols))]
    res = table(rows, f"composition of {len(cols)} samples across {len(show)} {lab} taxa")
    res["stats"]["figure"] = plot.grouped_bar(list(cols), series, ylabel="relative abundance",
                                             title=f"{lab} composition")
    return res


@T("qiime_merge_metadata", "Join a feature table with sample metadata", Q2, "table",
   [tbl("src", PHENO, "Metadata table"), tbl("extra", COUNTS, "Table to join"),
    textbox("left_key", "", "Key column in metadata"), textbox("right_key", "", "Key column in the other table"),
    choice("how", ["left", "inner"], "left", "Join type")],
   ex={"src": PHENO, "extra": PHENO, "left_key": "sample", "right_key": "sample"},
   up="qiime metadata join / phyloseq merge", tags=("metadata", "join"),
   summary="Attach one tabular annotation to another by a shared key column.")
def qiime_merge_metadata(src, extra, left_key="", right_key="", how="left"):
    """Metadata join."""
    dl, dr = _tbl(src), _tbl(extra)
    if dl.empty or dr.empty:
        return table([], "one of the tables is empty")
    kl = left_key or list(dl.columns)[0]
    kr = right_key or list(dr.columns)[0]
    right = {str(r[kr]): r for r in dr.to_dict("records")} if kr in dr.columns else {}
    rows = []
    matched = 0
    for r in dl.to_dict("records"):
        key = str(r.get(kl))
        other = right.get(key)
        if other:
            matched += 1
            rows.append({**{k: v for k, v in r.items()},
                        **{f"r_{k}": v for k, v in other.items() if k != kr}})
        elif how == "left":
            rows.append({**{k: v for k, v in r.items()}, "matched": False})
    return table(rows[:300], f"{matched} of {len(dl)} rows matched on {kl} = {kr}")


# ===========================================================================
# biodiversity and marker genes
# ===========================================================================
@T("biodiv_checklist", "Species checklist from a taxonomy table", BDE, "table",
   [tbl("taxonomy", TAXMAP, "Taxonomy table"), choice("rank", RANKS, "species", "Rank"),
    boolean("counts", True, "Include sequence counts")],
   ex={"taxonomy": TAXMAP, "rank": "species"}, up="GBIF checklist / spocc", tags=("checklist", "biodiversity"),
   summary="Unique names at a rank with how often each was observed.")
def biodiv_checklist(taxonomy, rank="species", counts=True):
    """Checklist generation."""
    tax = _lineage_table(taxonomy)
    lab = _rank_label(rank)
    agg: dict[str, list[str]] = defaultdict(list)
    for k, v in tax.items():
        agg[_rank_of(v, lab)].append(k)
    rows = [{"taxon": t, "rank": lab, "n_records": len(ids),
            "records": ",".join(ids[:6]),
            "distinct_higher": len({_rank_of(tax[i], "family") for i in ids if i in tax})}
           for t, ids in sorted(agg.items())]
    if not counts:
        rows = [{"taxon": r["taxon"], "rank": r["rank"]} for r in rows]
    return table(rows, f"{len(rows)} {lab} names from {len(tax)} records")


@T("biodiv_marker_match", "Match barcodes to a reference library", BDE, "table",
   [fa("records", GENES_FA, "Query barcodes"), fa("reference", GENES_FA, "Reference barcodes"),
    number("min_identity", 97.0, "Minimum identity %", min=0.0), intin("max_results", 5, "Hits per query", min=1)],
   ex={"records": GENES_FA, "reference": GENES_FA, "min_identity": 90.0}, up="BLAST / IDTaxa",
   tags=("barcode", "identification"),
   summary="Rank reference matches per query and flag confident identifications.")
def biodiv_marker_match(records, reference, min_identity=97.0, max_results=5):
    """Best-hit matching."""
    refs = _fa(reference)
    rows = []
    for r in _fa(records):
        cand = []
        for p in refs:
            if p.id == r.id:
                continue
            alen = max(1, min(len(r.seq), len(p.seq)))
            hd = dict(align.hamming(r.seq[:alen], p.seq[:alen]))
            ident = 100.0 * (alen - int(hd.get("distance", alen))) / alen
            cand.append((ident, p.id))
        cand.sort(reverse=True)
        top = cand[: int(max_results)]
        for rank, (ident, rid) in enumerate(top, 1):
            rows.append({"query": r.id, "rank": rank, "reference": rid, "identity_percent": round(ident, 3),
                        "match": bool(ident >= float(min_identity)) if rank == 1 else ""})
    n_ok = sum(1 for r in rows if r["match"] is True)
    return table(rows[:300], f"{n_ok} queries have a confident match at >= {min_identity}% identity")


@T("biodiv_barcode_gap", "Distance-based species delimitation", BDE, "table",
   [fa("records", GENES_FA, "Barcode alignment"), number("threshold", 0.03, "Delimitation distance", min=0.0),
    intin("head", 60, "Pair rows", min=1)],
   ex={"records": GENES_FA, "threshold": 0.2}, up="ABGD / spades", tags=("delimitation", "distance"),
   summary="Compare within- and between-group barcode distances to find a barcoding gap.")
def biodiv_barcode_gap(records, threshold=0.03, head=60):
    """Barcoding-gap statistics."""
    recs = _fa(records)
    if len(recs) < 3:
        return table([], "need at least three sequences")
    D = np.asarray(align.distance_matrix([r.seq for r in recs], model="identity"), dtype=float)
    pairs = [(float(D[i][j]), i, j) for i in range(len(recs)) for j in range(i + 1, len(recs))]
    within = sorted(d for d, _i, _j in pairs if d <= float(threshold))
    between = sorted(d for d, _i, _j in pairs if d > float(threshold))
    rows = [{"metric": "pairs", "value": len(pairs)},
            {"metric": "below_threshold", "value": len(within)},
            {"metric": "mean_within", "value": round(float(np.mean(within)), 5) if within else 0.0},
            {"metric": "max_within", "value": round(float(max(within)), 5) if within else 0.0},
            {"metric": "mean_between", "value": round(float(np.mean(between)), 5) if between else 0.0},
            {"metric": "min_between", "value": round(float(min(between)), 5) if between else 0.0},
            {"metric": "gap_present", "value": bool(within and between and min(between) > max(within))}]
    res = table(rows, f"barcoding gap {'detected' if rows[-1]['value'] else 'not detected'} "
                    f"at d = {threshold}")
    res["stats"]["figure"] = plot.histogram([p[0] for p in pairs], bins=20, xlabel="distance",
                                           title="Pairwise distance distribution")
    return res


@T("biodiv_endemism_score", "Range-size proxy from taxon record counts", BDE, "table",
   [tbl("taxonomy", TAXMAP, "Taxonomy table"), choice("rank", RANKS, "genus", "Rank"),
    intin("head", 60, "Rows", min=1)],
   ex={"taxonomy": TAXMAP, "rank": "genus"}, up="spocc / endemicity index",
   tags=("biodiversity", "endemism"),
   summary="Score each taxon by how narrowly it is represented in the sample set.")
def biodiv_endemism_score(taxonomy, rank="genus", head=60):
    """Endemism-like index."""
    tax = _lineage_table(taxonomy)
    lab = _rank_label(rank)
    per_taxon: dict[str, set] = defaultdict(set)
    fam_of: dict[str, set] = defaultdict(set)
    for k, v in tax.items():
        t = _rank_of(v, lab)
        per_taxon[t].add(k)
        fam_of[t].add(_rank_of(v, "family"))
    n = max(1, len(tax))
    rows = []
    for t, ids in per_taxon.items():
        c = len(ids)
        rows.append({"taxon": t, "rank": lab, "records": c, "families": len(fam_of[t]),
                    "range_fraction": round(c / n, 5),
                    "endemism_score": round(1.0 / (1.0 + math.log2(max(1, c * max(1, len(fam_of[t]))))), 5)})
    rows.sort(key=lambda r: -r["endemism_score"])
    return table(rows[: int(head)], f"{len(rows)} {lab} taxa scored")


# ===========================================================================
# plant / RAD-seq tools
# ===========================================================================
@T("planttribes_enzyme_cuts", "Restriction enzyme cut frequencies", PT, "table",
   [fa("records", GENES_FA, "FASTA sequences"), multi("enzymes", ["EcoRI", "HindIII", "BamHI", "SbfI", "EcoT22I",
                                                                "PstI", "NsiI"], ["EcoRI", "HindIII", "SbfI"],
                                                   "Enzymes"), intin("min_cut", 0, "Minimum cuts per read", min=0)],
   ex={"records": GENES_FA, "enzymes": ["EcoRI", "HindIII"]}, up="TASSEL / ipyrad --enzymes",
   tags=("RAD-seq", "restriction enzyme"),
   summary="Count how often each enzyme cuts every sequence in a RAD library.")
def planttribes_enzyme_cuts(records, enzymes=("EcoRI", "HindIII"), min_cut=0):
    """Enzyme cut table."""
    names = [str(e) for e in (enzymes or [])]
    recs = _fa(records)
    rows = []
    for r in recs:
        counts = seq.cut_frequencies(r.seq, enzymes=names)
        rec = {"sequence": r.id, "length": len(r.seq)}
        keep = True
        for c in counts:
            d = dict(c)
            name = str(d.get("enzyme", d.get("name", "")))
            n = int(d.get("cuts", d.get("count", 0)) or 0)
            rec[name] = n
            rec[f"{name}_per_kb"] = round(float(d.get("per_kb", 0.0) or 0.0), 4)
            if int(min_cut) and n < int(min_cut):
                keep = False
        if keep:
            rows.append(rec)
    return table(rows, f"{len(rows)} sequences x {len(names)} enzymes")


@T("planttribes_align_contigs", "Digest contigs into RAD loci", PT, "text",
   [fa("records", GENES_FA, "Contig FASTA"), choice("enzyme", ["EcoRI", "HindIII", "BamHI", "PstI"], "EcoRI",
                                                   "Enzyme"), intin("min_len", 30, "Minimum locus length", min=1),
    boolean("sticky_end", True, "Require a compatible overhang")],
   ex={"records": GENES_FA, "enzyme": "EcoRI", "min_len": 50}, up="ipyrad / UNEAF1 demultiplex",
   tags=("RAD-seq", "digestion"),
   summary="Simulated enzyme digestion of contigs into plausibly-sized loci.")
def planttribes_align_contigs(records, enzyme="EcoRI", min_len=30, sticky_end=True):
    """Digest contigs into loci."""
    out = []
    for r in _fa(records):
        frags = seq.fragment_sizes(r.seq, enzyme=enzyme)
        pieces = seq.digest(r.seq, enzyme=enzyme)
        for i, f in enumerate(pieces):
            s = f.seq if isinstance(f, io.Seq) else str(f)
            if len(s) < int(min_len):
                continue
            tag = "R" if (sticky_end and s[:1].upper() in "AG") else "C"
            out.append(io.Seq(id=f"{r.id}_{i + 1}{tag}", seq=s, desc=f"enzyme={enzyme};loci={len(frags)}"))
    return text(io.write_fasta(out), f"{len(out)} loci from {len(_fa(records))} contigs "
                                    f"({enzyme}, >= {min_len} bp)")


# ===========================================================================
# virology
# ===========================================================================
@T("viro_coverage_by_target", "Read coverage per genomic target", VIRO, "table",
   [sam("alignment", ALN_SAM, "SAM alignment"), bed("targets", REGIONS_BED, "BED targets"),
    number("min_depth", 1.0, "Minimum depth to call covered", min=0.0)],
   ex={"alignment": ALN_SAM, "targets": REGIONS_BED, "min_depth": 1.0}, up="ivar trim / nCoV-2019 depth",
   tags=("coverage", "virus", "amplicon"),
   summary="Mean and minimum depth over each target region of a viral panel.")
def viro_coverage_by_target(alignment, targets, min_depth=1.0):
    """Per-target coverage of an alignment."""
    _header, alns = io.parse_sam(io.as_text(alignment))
    ivs = _bed(targets)
    hits = Counter()
    for a in alns:
        if a.unmapped or not a.rname:
            continue
        for k, iv in enumerate(ivs):
            if iv.chrom == a.rname and a.pos + 1 <= iv.end and a.ref_end >= iv.start + 1:
                hits[k] += 1
    rows = []
    for k, iv in enumerate(ivs):
        dep = float(hits[k])
        rows.append({"chrom": iv.chrom, "start": iv.start, "end": iv.end, "name": iv.name or f"target_{k + 1}",
                    "reads": int(dep), "breadth": round(float(dep) / max(1, len(alns)), 4),
                    "covered": bool(dep >= float(min_depth))})
    cov = sum(1 for r in rows if r["covered"])
    res = table(rows, f"{cov}/{len(rows)} targets above depth {min_depth} "
                     f"({len(alns)} alignments)")
    res["stats"]["figure"] = plot.bar([str(r["name"]) for r in rows], [r["reads"] for r in rows],
                                     xlabel="target", ylabel="reads", title="Reads per target", horizontal=True)
    return res


@T("viro_consensus_from_variants", "Build a consensus genome from a VCF", VIRO, "text",
   [anyfile("reference", GENOME, fmt="", label="Reference FASTA"), vcf("variants", VARIANTS_VCF, "VCF file"),
    choice("ambiguity", ["majority", "iupac", "reference"], "majority", "Heterozygous sites"),
    number("min_quality", 0.0, "Minimum variant quality", min=0.0)],
   ex={"reference": GENOME, "variants": VARIANTS_VCF}, up="bcftools consensus / iVar consensus",
   tags=("consensus", "variants", "virus"),
   summary="Apply called variants to a reference to obtain a consensus sequence.")
def viro_consensus_from_variants(reference, variants, ambiguity="majority", min_quality=0.0):
    """Consensus builder."""
    gsrc, _sizes = genome.read_genome(reference)
    vcfobj = _vcf(variants)
    recs = list(vcfobj.records) if hasattr(vcfobj, "records") else list(vcfobj)
    applied = 0
    out = []
    for chrom, bases in gsrc.items():
        seq_list = list(bases)
        for r in recs:
            d = r if isinstance(r, dict) else {k: getattr(r, k, None) for k in ("chrom", "pos", "ref", "alt",
                                                                                "qual")}
            if str(d.get("chrom", d.get("chrom", ""))) != chrom:
                continue
            try:
                q = float(d.get("qual") or 0)
            except (TypeError, ValueError):
                q = 0.0
            if q < float(min_quality):
                continue
            p = int(d.get("pos") or 0) - 1
            alt = str(d.get("alt") or "")
            ref = str(d.get("ref") or "")
            if p < 0 or p >= len(seq_list) or not alt or alt in (".",):
                continue
            if ambiguity == "iupac" and "," in alt:
                iupac = {"AG": "R", "CT": "Y", "CG": "S", "AT": "W", "GT": "K", "AC": "M",
                        "GAT": "D", "ACT": "H", "AGT": "V", "CGT": "B"}
                seq_list[p] = iupac.get("".join(sorted(set(alt.replace(",", "")))), alt[0])
            elif len(alt) > len(ref) and alt.startswith(ref):
                seq_list[p] = alt[0]
            else:
                seq_list[p] = alt.split(",")[0][0]
            applied += 1
        out.append(io.Seq(id=chrom, seq="".join(seq_list), desc=f"consensus;variants_applied={applied}"))
    return text(io.write_fasta(out), f"consensus written; {applied} variant positions applied "
                                    f"({len(recs)} records considered)")


@T("viro_drift_distances", "Antigenic drift matrix between isolates", VIRO, "table",
   [fa("proteins", PROTEINS, "Protein FASTA"), choice("model", ["identity", "p", "k2p"], "identity", "Model"),
    number("cluster_cutoff", 0.05, "Cluster distance", min=0.0)],
   ex={"proteins": PROTEINS, "model": "identity"}, up="Nextclade / phydrift",
   tags=("evolution", "distance", "drift"),
   summary="Pairwise divergence between isolates plus a simple drift cluster label.")
def viro_drift_distances(proteins, model="identity", cluster_cutoff=0.05):
    """Drift matrix."""
    recs = _fa(proteins)
    if len(recs) < 2:
        return table([], "need at least two sequences")
    D = np.asarray(align.distance_matrix([r.seq for r in recs], model=model), dtype=float)
    lab = [0] * len(recs)
    nxt = 0
    for i in range(len(recs)):
        if lab[i]:
            continue
        nxt += 1
        lab[i] = nxt
        for j in range(i + 1, len(recs)):
            if D[i][j] <= float(cluster_cutoff):
                lab[j] = nxt
    rows = [{"isolate_a": recs[i].id, "isolate_b": recs[j].id, "distance": round(float(D[i][j]), 6),
            "percent_difference": round(100 * float(D[i][j]), 3), "model": model}
           for i in range(len(recs)) for j in range(i + 1, len(recs))]
    res = table(rows[:300], f"{len(rows)} pairs, mean distance {round(float(D[np.triu_indices(len(recs), 1)].mean()), 5)}")
    res["stats"]["clusters"] = {recs[i].id: int(lab[i]) for i in range(len(recs))}
    res["stats"]["figure"] = plot.heatmap([[round(float(v), 4) for v in row] for row in D],
                                        row_labels=[r.id for r in recs], col_labels=[r.id for r in recs],
                                        title="drift distances", annot=True)
    return res


@T("viro_variant_load", "Allele fractions and variant load", VIRO, "table",
   [vcf("variants", VARIANTS_VCF, "VCF file"), intin("min_depth", 0, "Minimum depth", min=0),
    number("min_fraction", 0.0, "Minor-allele fraction cut-off", min=0.0)],
   ex={"variants": VARIANTS_VCF, "min_fraction": 0.2}, up="lofreq report / ivar variants",
   tags=("variants", "quasispecies"),
   summary="Per-site allele balance and the intra-host variant load of a sample.")
def viro_variant_load(variants, min_depth=0, min_fraction=0.0):
    """Variant load summary."""
    vcfobj = _vcf(variants)
    keys = list(vcfobj.samples) if vcfobj.samples else ["sample"]
    rows = []
    for rec in vcfobj.records:
        fmt = dict((rec.samples[0] if rec.samples else {}) or {})
        info = dict(rec.info or {})
        dep = _to_int(fmt.get("DP", info.get("DP")), 0)
        ad = [int(x) for x in re.findall(r"\d+", str(fmt.get("AD", info.get("AD", ""))))]
        af = _to_float(fmt.get("AF", info.get("AF")), None)
        vaf = (ad[1] / dep) if dep and len(ad) > 1 else (af if af is not None else (1.0 if dep == 0 else 0.0))
        if dep < int(min_depth) or vaf < float(min_fraction):
            continue
        ref_base = str(rec.ref)[:1].upper()
        alt_base = str(rec.alt).split(",")[0][:1].upper()
        ti = {"AG", "CT"}
        rows.append({"chrom": rec.chrom, "pos": rec.pos, "id": rec.id, "ref": rec.ref,
                    "alt": str(rec.alt).replace(",", "/"), "depth": dep,
                    "allele_fraction": round(float(vaf), 4), "quality": rec.qual, "filter": rec.filter,
                    "counts": "/".join(str(x) for x in ad),
                    "substitution_type": "transition" if (ref_base + alt_base) in ti and len(ref_base) == 1
                    else "transversion"})
    ts = sum(1 for r in rows if r["substitution_type"] == "transition")
    tv = len(rows) - ts
    res = table(rows, f"{len(rows)} variants pass (Ti/Tv = {round(ts / max(1, tv), 3)})")
    res["stats"]["figure"] = (plot.histogram([r["allele_fraction"] for r in rows], bins=15,
                                            xlabel="allele fraction", title="Variant load") if rows
                             else plot.empty_plot("no variants"))
    return res


def _to_int(v, default=0):
    try:
        return int(float(str(v).split(",")[0]))
    except (TypeError, ValueError):
        return default


def _to_float(v, default=None):
    try:
        return float(str(v).split(",")[0])
    except (TypeError, ValueError):
        return default


# ===========================================================================
# epigenetics / methylation
# ===========================================================================
@T("epi_cpg_islands", "Find CpG islands in a genome", EPI, "text",
   [anyfile("reference", GENOME, fmt="", label="Genome FASTA"), intin("min_length", 200, "Minimum length", min=20),
    number("gc_fraction", 0.5, "Minimum GC fraction", min=0.0), number("obs_exp", 0.6, "Minimum CpG obs/exp",
                                                                    min=0.0), intin("window", 100, "Scan window",
                                                                                 min=20)],
   ex={"reference": GENOME, "min_length": 200, "gc_fraction": 0.4, "obs_exp": 0.4},
   up="cpgislandfinder / newcpgreport (EMBOSS)", tags=("CpG island", "promoter", "annotation"),
   summary="Sliding-window detection of CpG-rich, GC-rich regions reported as BED.")
def epi_cpg_islands(reference, min_length=200, gc_fraction=0.5, obs_exp=0.6, window=100):
    """CpG island detection."""
    gsrc, _sizes = genome.read_genome(reference)
    w = max(20, int(window))
    out = []
    for chrom, bases in gsrc.items():
        s = "".join(bases).upper()
        cur = None
        for i in range(0, max(1, len(s) - w + 1), w):
            sub = s[i:i + w]
            st = dict(seq.cpg_obs_exp(sub))
            gc = float(st.get("gc", 100.0 * seq.gc_content(sub)))
            n_cg = float(st.get("cg", 0.0))
            n_c = float(st.get("c_count", 0.0))
            n_g = float(st.get("g_count", 0.0))
            exp = n_c * n_g / max(1, len(sub))
            oe = n_cg / exp if exp > 0 else 0.0
            ok = gc / 100.0 >= float(gc_fraction) and oe >= float(obs_exp)
            if ok:
                if cur and i == cur[2]:
                    cur = (cur[0], cur[1], i + w)
                else:
                    if cur:
                        out.append(cur)
                    cur = (chrom, i, i + w)
            elif cur:
                out.append(cur)
                cur = None
        if cur:
            out.append(cur)
    ivs = [io.Interval(c, a, b, f"CpG_island_{k + 1}", 0, ".")
          for k, (c, a, b) in enumerate(out) if b - a >= int(min_length)]
    body = io.write_bed(ivs, bed12=False) if ivs else ""
    return text(body, f"{len(ivs)} CpG islands >= {min_length} bp (GC >= {100 * gc_fraction}%, "
                     f"obs/exp >= {obs_exp})")


@T("epi_methylation_context_summary", "Methylation by cytosine context", EPI, "table",
   [anyfile("methyl", CXCG, fmt="", label="Bismark / MethylDackel output"),
    choice("context", ["all", "CpG", "CHG", "CHH"], "all", "Context"), number("min_coverage", 0.0,
                                                        "Minimum per-cytosine coverage", min=0.0)],
   ex={"methyl": CXCG}, up="bismark2bedGraph / MethylDackel",
   tags=("DNA methylation", "context"),
   summary="Mean methylation and coverage of cytosines in each sequence context.")
def epi_methylation_context_summary(methyl, context="all", min_coverage=0.0):
    """Context-level methylation stats."""
    rows_in = _cxcg(methyl)
    if not rows_in:
        return table([], "no methylation rows parsed")
    by_ctx: dict[str, list[dict]] = defaultdict(list)
    for r in rows_in:
        if context != "all" and r["context"] != context.upper():
            continue
        by_ctx[r["context"]].append(r)
    rows = []
    for ctx, items in sorted(by_ctx.items()):
        f = [r["fraction"] for r in items if math.isfinite(r["fraction"])]
        cov = [r["methylated"] + r["unmethylated"] for r in items
               if math.isfinite(r["methylated"]) and math.isfinite(r["unmethylated"])]
        cov = [c for c in cov if c >= float(min_coverage)]
        rows.append({"context": ctx, "cytosines": len(items),
                    "mean_methylation": round(float(np.mean(f)), 5) if f else 0.0,
                    "median_methylation": round(float(np.median(f)), 5) if f else 0.0,
                    "fully_methylated": int(sum(1 for x in f if x >= 0.8)),
                    "fully_unmethylated": int(sum(1 for x in f if x <= 0.2)),
                    "mean_coverage": round(float(np.mean(cov)), 2) if cov else 0.0,
                    "chromosomes": len({r["chrom"] for r in items})})
    res = table(rows, f"{len(rows_in)} cytosines in {len(rows)} contexts")
    res["stats"]["figure"] = plot.bar([r["context"] for r in rows], [r["mean_methylation"] for r in rows],
                                     xlabel="context", ylabel="methylation", title="Mean methylation")
    return res


@T("epi_methylation_per_region", "Methylation over annotated regions", EPI, "table",
   [anyfile("methyl", CXCG, fmt="", label="Methylation table"),
    bed("regions", REGIONS_BED, "BED regions"), number("min_fraction", 0.0, "Report only regions above", min=0.0),
    intin("head", 200, "Rows", min=1)],
   ex={"methyl": CXCG, "regions": REGIONS_BED}, up="methylprep / DSS merge",
   tags=("DNA methylation", "regions"),
   summary="Aggregate cytosine-level methylation into per-region betas.")
def epi_methylation_per_region(methyl, regions, min_fraction=0.0, head=200):
    """Per-region methylation."""
    rows_in = _cxcg(methyl)
    ivs = _bed(regions)
    by_chrom: dict[str, list[dict]] = defaultdict(list)
    for r in rows_in:
        by_chrom[r["chrom"]].append(r)
    out = []
    for iv in ivs:
        items = [r for r in by_chrom.get(iv.chrom, []) if iv.start < r["pos"] <= iv.end]
        if not items:
            out.append({"chrom": iv.chrom, "start": iv.start, "end": iv.end, "name": iv.name or ".",
                       "cytosines": 0, "mean_methylation": 0.0, "coverage": 0})
            continue
        f = [r["fraction"] for r in items if math.isfinite(r["fraction"])]
        cov = int(sum((r["methylated"] + r["unmethylated"]) for r in items
                    if math.isfinite(r["methylated"]) and math.isfinite(r["unmethylated"])))
        mean_f = float(np.mean(f)) if f else 0.0
        if mean_f < float(min_fraction):
            continue
        out.append({"chrom": iv.chrom, "start": iv.start, "end": iv.end, "name": iv.name or ".",
                   "cytosines": len(items), "mean_methylation": round(mean_f, 5),
                   "methylated_calls": round(sum(r["fraction"] for r in items if math.isfinite(r["fraction"])), 2),
                   "coverage": cov, "high_methylation": bool(mean_f >= 0.5)})
    return table(out[: int(head)], f"{len(out)} regions summarised from {len(rows_in)} cytosines")


@T("epi_dmr_scan", "Sliding-window DMR detection", EPI, "text",
   [anyfile("methyl", CXCG, fmt="", label="Methylation table"),
    intin("window", 50, "Window size (cytosines)", min=3), intin("step", 25, "Step", min=1),
    number("min_diff", 0.2, "Minimum mean methylation difference", min=0.0),
    number("min_fraction", 0.5, "Hyper-methylation threshold", min=0.0), boolean("bed12", False, "BED12 output")],
   ex={"methyl": CXCG, "window": 20, "step": 10, "min_diff": 0.1},
   up="DSS / methylKit dmr", tags=("DMR", "differential methylation"),
   summary="Call contiguous stretches of cytosines that cross a methylation threshold.")
def epi_dmr_scan(methyl, window=50, step=25, min_diff=0.2, min_fraction=0.5, bed12=False):
    """Simple DMR caller."""
    rows_in = _cxcg(methyl)
    if not rows_in:
        return text("", "no methylation rows parsed")
    by_chrom: dict[str, list[dict]] = defaultdict(list)
    for r in rows_in:
        by_chrom[r["chrom"]].append(r)
    ivs = []
    for chrom, items in sorted(by_chrom.items()):
        items.sort(key=lambda r: r["pos"])
        w, st = max(3, int(window)), max(1, int(step))
        for i in range(0, max(1, len(items) - w + 1), st):
            win = items[i:i + w]
            if len(win) < w:
                break
            f = [r["fraction"] for r in win if math.isfinite(r["fraction"])]
            if len(f) < w:
                continue
            mean_f = float(np.mean(f))
            spread = float(np.max(f) - np.min(f))
            if mean_f >= float(min_fraction) and spread >= float(min_diff):
                ivs.append(io.Interval(chrom, win[0]["pos"] - 1, win[-1]["pos"],
                                      f"DMR_{len(ivs) + 1}", int(round(100 * mean_f)), "."))
    body = io.write_bed(genome.merge(ivs, distance=1), bed12=bool(bed12)) if ivs else ""
    merged = genome.merge(ivs, distance=1) if ivs else []
    return text(body, f"{len(ivs)} windows called, {len(merged)} DMRs after merging "
                    f"(>= {w} cytosines, mean >= {min_fraction}, spread >= {min_diff})")


@T("epi_bisulfite_extractor", "Per-cytosine extraction table", EPI, "table",
   [anyfile("methyl", CXCG, fmt="", label="Bismark CXCG output"),
    choice("context", ["all", "CpG", "CHG", "CHH"], "all", "Context"),
    number("min_fraction", 0.0, "Minimum methylation", min=0.0),
    number("max_fraction", 1.0, "Maximum methylation", min=0.0), intin("head", 300, "Rows", min=1)],
   ex={"methyl": "examples/methylation.cxcg", "context": "CHG"}, up="bismark --methylation_extractor",
   tags=("bisulfite", "cytosine", "table"),
   summary="Clean per-cytosine table with context, fraction and call for each site.")
def epi_bisulfite_extractor(methyl, context="all", min_fraction=0.0, max_fraction=1.0, head=300):
    """Methylation extractor."""
    items = _cxcg(methyl)
    rows = []
    for r in items:
        if context != "all" and r["context"] != context.upper():
            continue
        f = r["fraction"]
        if not math.isfinite(f) or f < float(min_fraction) or f > float(max_fraction):
            continue
        call = "methylated" if f >= 0.5 else "unmethylated"
        rows.append({"chrom": r["chrom"], "start": r["pos"] - 1, "end": r["pos"], "strand": r["strand"],
                    "context": r["context"], "methylation_fraction": round(f, 4), "call": call,
                    "methylated_reads": r["methylated"], "unmethylated_reads": r["unmethylated"]})
    return table(rows[: int(head)], f"{len(rows)} cytosines reported ({context}, "
                                 f"{min_fraction}-{max_fraction})")


@T("epi_coverage_track_stats", "Coverage track summary per region", EPI, "table",
   [anyfile("graph", BEDGRAPH, fmt="", label="bedGraph / bigWig-like track"),
    bed("regions", REGIONS_BED, "BED regions"), choice("mode", ["mean", "median", "max", "sum"], "mean",
                                                     "Statistic")],
   ex={"graph": BEDGRAPH, "regions": REGIONS_BED}, up="computeMatrix / bigWigStatistics",
   tags=("coverage", "chromatin"),
   summary="Reduce a signal track to one value per region for downstream modelling.")
def epi_coverage_track_stats(graph, regions, mode="mean"):
    """Region-wise track statistics."""
    vals = io.parse_bedgraph(io.as_text(graph))
    by_chrom: dict[str, list[tuple[int, int, float]]] = defaultdict(list)
    for c, s, e, v in vals:
        by_chrom[c].append((int(s), int(e), float(v)))
    fn = {"mean": np.mean, "median": np.median, "max": np.max, "sum": np.sum}[mode]
    rows = []
    for iv in _bed(regions):
        seg = [v for s, e, v in by_chrom.get(iv.chrom, []) if e > iv.start and s < iv.end]
        rows.append({"chrom": iv.chrom, "start": iv.start, "end": iv.end, "name": iv.name or ".",
                    "segments": len(seg), mode: round(float(fn(seg)), 5) if seg else 0.0,
                    "span_bp": iv.end - iv.start})
    return table(rows, f"{len(rows)} regions scored by {mode}")


# ===========================================================================
# STR / microsatellites
# ===========================================================================
@T("str_allele_table", "Long-format STR allele table", STRSEC, "table",
   [tbl("profile", STR, "STR profile (locus rows, sample columns)"), boolean("expand", True, "One row per allele")],
   ex={"profile": STR}, up="GeneMapper / STRAF", tags=("STR", "microsatellite"),
   summary="Tidy an ampFLPSTR-style profile into one row per allele call.")
def str_allele_table(profile, expand=True):
    """Tidy the STR profile."""
    df = _tbl(profile)
    if df.empty:
        return table([], "empty profile")
    cols = [str(c) for c in df.columns]
    meta = [c for c in cols if c.lower() in ("locus", "marker", "motif", "chr", "chrom", "start")]
    samples = [c for c in cols if c not in meta]
    rows = []
    for r in df.to_dict("records"):
        locus = str(r.get(meta[0] if meta else cols[0]))
        motif = str(r.get("motif", r.get("Motif", "")))
        chrom = str(r.get("chr", r.get("Chr", r.get("chrom", ""))))
        for s in samples:
            cell = str(r.get(s, "")).strip()
            alleles = [a for a in re.split(r"[/,;|]", cell) if a.strip()] if expand else [cell]
            for a in alleles:
                rows.append({"locus": locus, "motif": motif, "chrom": chrom, "sample": s,
                            "allele": a.strip(), "repeat_size": len(a.strip()),
                            "numeric_value": float(a) if re.fullmatch(r"-?\d+(\.\d+)?", a.strip()) else float("nan")})
    return table(rows[:600], f"{len(rows)} allele calls across {len({r['sample'] for r in rows})} samples "
                          f"and {len({r['locus'] for r in rows})} loci")


@T("str_locus_statistics", "Allelic diversity per STR locus", STRSEC, "table",
   [tbl("profile", STR, "STR profile"), boolean("expected_heterozygosity", True, "Report Ho vs He")],
   ex={"profile": STR}, up="STRAF / poppr", tags=("STR", "population genetics"),
   summary="Allele counts, observed and expected heterozygosity per microsatellite.")
def str_locus_statistics(profile, expected_heterozygosity=True):
    """Per-locus diversity."""
    df = _tbl(profile)
    cols = [str(c) for c in df.columns]
    meta = [c for c in cols if c.lower() in ("locus", "motif", "chr", "start", "chrom")]
    samples = [c for c in cols if c not in meta]
    rows = []
    for r in df.to_dict("records"):
        locus = str(r.get(meta[0] if meta else cols[0]))
        calls = []
        for s in samples:
            cell = str(r.get(s, "")).strip()
            alleles = [a.strip() for a in re.split(r"[/,;|]", cell) if a.strip()]
            if alleles:
                calls.append(alleles)
        pool = [a for c in calls for a in c]
        n = len(pool) or 1
        freqs = Counter(pool)
        he = 1.0 - sum((c / n) ** 2 for c in freqs.values())
        ho = float(np.mean([1 if len(set(c)) > 1 else 0 for c in calls])) if calls else 0.0
        rows.append({"locus": locus, "samples": len(calls), "alleles": len(freqs),
                    "mean_allele_size": round(float(np.mean([len(a) for a in pool])), 3) if pool else 0.0,
                    "observed_heterozygosity": round(ho, 4),
                    "expected_heterozygosity": round(he, 4) if expected_heterozygosity else "",
                    "polymorphic_info_content": round(he, 4),
                    "most_common_allele": freqs.most_common(1)[0][0] if freqs else "",
                    "shannon": round(float(stats.shannon2(list(freqs.values()))), 4)})
    return table(rows, f"{len(rows)} loci typed")


@T("str_profile_match", "Compare two STR profiles", STRSEC, "table",
   [tbl("profile", STR, "STR profile"), textbox("sample_a", "", "First sample"), textbox("sample_b", "",
                                                                                         "Second sample"),
    intin("tolerance", 1, "Allele size tolerance", min=0)],
   ex={"profile": STR, "sample_a": "sample_A", "sample_b": "sample_B"}, up="match.py / CODIS comparison",
   tags=("STR", "matching", "forensics"),
   summary="Count matching, off-by-one and discordant loci between two profiles.")
def str_profile_match(profile, sample_a="", sample_b="", tolerance=1):
    """Profile comparison."""
    df = _tbl(profile)
    cols = [str(c) for c in df.columns]
    samples = [c for c in cols if c.lower() not in ("locus", "motif", "chr", "start", "chrom")]
    a = sample_a if sample_a in samples else (samples[0] if samples else cols[0])
    b = sample_b if sample_b in samples else (samples[1] if len(samples) > 1 else cols[-1])
    lcol = cols[0]
    rows, match = [], 0
    for r in df.to_dict("records"):
        la = {x.strip() for x in re.split(r"[/,;|]", str(r.get(a, ""))) if x.strip()}
        lb = {x.strip() for x in re.split(r"[/,;|]", str(r.get(b, ""))) if x.strip()}

        def nums(ss):
            return [float(x) for x in ss if re.fullmatch(r"-?\d+(\.\d+)?", x)]
        na, nb = nums(la), nums(lb)
        shared = len(la & lb)
        near = sum(1 for x in na for y in nb if 0 < abs(x - y) <= int(tolerance))
        status = "identical" if la and la == lb else ("partial" if shared or near else "discordant")
        match += 1 if shared else 0
        rows.append({"locus": str(r.get(lcol)), "profile_a": "/".join(sorted(la)), "profile_b": "/".join(sorted(lb)),
                    "shared_alleles": shared, "near_matches": near, "status": status})
    res = table(rows, f"{match}/{len(rows)} loci share at least one allele between {a} and {b}")
    res["stats"]["figure"] = plot.bar([r["locus"] for r in rows], [r["shared_alleles"] for r in rows],
                                     xlabel="locus", ylabel="shared alleles", title="Profile comparison")
    return res


@T("str_relatedness", "Pairwise relatedness from shared alleles", STRSEC, "table",
   [tbl("profile", STR, "STR profile"), number("min_loci", 1, "Minimum informative loci", min=1)],
   ex={"profile": STR}, up="King relatedness / LRmix", tags=("relatedness", "kinship"),
   summary="Estimate IBS sharing between samples to flag duplicates and relatives.")
def str_relatedness(profile, min_loci=1):
    """IBS / LR-style relatedness."""
    df = _tbl(profile)
    cols = [str(c) for c in df.columns]
    samples = [c for c in cols if c.lower() not in ("locus", "motif", "chr", "start", "chrom")]
    calls: dict[str, list[set]] = {s: [] for s in samples}
    for r in df.to_dict("records"):
        for s in samples:
            calls[s].append({x.strip() for x in re.split(r"[/,;|]", str(r.get(s, ""))) if x.strip()})
    rows = []
    for i in range(len(samples)):
        for j in range(i + 1, len(samples)):
            a, b = samples[i], samples[j]
            ibs0 = ibs1 = ibs2 = 0
            for ca, cb in zip(calls[a], calls[b]):
                sh = len(ca & cb)
                if sh >= 2:
                    ibs2 += 1
                elif sh == 1:
                    ibs1 += 1
                else:
                    ibs0 += 1
            tot = ibs0 + ibs1 + ibs2
            if tot < int(min_loci):
                continue
            phi = (ibs2 + 0.5 * ibs1) / tot if tot else 0.0
            rows.append({"sample_a": a, "sample_b": b, "loci": tot, "ibs2": ibs2, "ibs1": ibs1, "ibs0": ibs0,
                        "sharing_fraction": round(phi, 4),
                        "relationship": "identical" if phi > 0.95 else ("parent-offspring/sib" if phi > 0.4
                                                                      else ("unrelated" if phi < 0.15 else "distant"))})
    rows.sort(key=lambda r: -r["sharing_fraction"])
    return table(rows, f"{len(rows)} sample pairs scored")


@T("str_expansion_test", "Test for repeat expansions at a locus", STRSEC, "table",
   [tbl("profile", STR, "STR profile"), textbox("locus", "", "Locus name"), number("normal_mean", 0.0,
                                                            "Reference mean size (0 = cohort)", min=0.0),
    number("n_sds", 3.0, "Number of standard deviations", min=1.0)],
   ex={"profile": STR, "locus": "VSTR1"}, up="hipSTRyda / ExpansionHunter", tags=("repeat expansion", "STR"),
   summary="Flag samples whose allele sizes are outliers for the locus repeat structure.")
def str_expansion_test(profile, locus="", normal_mean=0.0, n_sds=3.0):
    """Expansion screening."""
    df = _tbl(profile)
    cols = [str(c) for c in df.columns]
    lcol = cols[0]
    samples = [c for c in cols if c.lower() not in (lcol.lower(), "motif", "chr", "start", "chrom")]
    rows_wanted = [r for r in df.to_dict("records")
                  if not locus or str(r.get(lcol)) == str(locus)]
    if not rows_wanted:
        return table([], f"locus '{locus}' not found (available: {', '.join(str(r.get(lcol)) for r in df.to_dict('records')[:6])})")
    out = []
    for r in rows_wanted:
        vals: dict[str, float] = {}
        for s in samples:
            nums = [float(x) for x in re.findall(r"\d+(?:\.\d+)?", str(r.get(s, "")))]
            if nums:
                vals[s] = float(np.mean(nums))
        if not vals:
            continue
        sizes = np.array(list(vals.values()), dtype=float)
        mu = float(normal_mean) if float(normal_mean) > 0 else float(np.median(sizes))
        sd = float(np.std(sizes, ddof=1)) if sizes.size > 1 else 0.0
        for s, v in vals.items():
            z = (v - mu) / sd if sd > 0 else 0.0
            out.append({"locus": str(r.get(lcol)), "sample": s, "mean_allele_size": round(v, 3),
                       "cohort_median": round(mu, 3), "cohort_sd": round(sd, 3),
                       "deviation_sd": round(z, 4), "expansion_called": bool(sd > 0 and z >= float(n_sds)),
                       "motif": str(r.get("motif", r.get("Motif", "")))})
    return table(out, f"{sum(1 for r in out if r['expansion_called'])} expansions at "
                   f">= {n_sds} SD over {len(out)} calls")


# ===========================================================================
# small specialty toolkits
# ===========================================================================
@T("sccaf_split_contigs", "Split contigs into fragment files", SCCAF, "text",
   [fa("contigs", GENOME, "Contig FASTA"), intin("max_bp", 500, "Fragment size", min=100),
    boolean("keep_order", True, "Number fragments in order")],
   ex={"contigs": GENES_FA, "max_bp": 300}, up="sccaf split_contigs", tags=("scCaF", "FASTA", "split"),
   summary="Break long sequences into fixed-size fragments with sequential ids.")
def sccaf_split_contigs(contigs, max_bp=500, keep_order=True):
    """Split FASTA into chunks."""
    size = max(10, int(max_bp))
    out = []
    for r in _fa(contigs):
        for k in range(0, len(r.seq), size):
            piece = r.seq[k:k + size]
            idx = k // size + 1
            out.append(io.Seq(id=(f"{r.id}_{idx:03d}" if keep_order else f"frag_{len(out) + 1}"),
                             seq=piece, desc=f"parent={r.id};start={k + 1}"))
    return text(io.write_fasta(out), f"{len(out)} fragments of <= {size} bp from {len(_fa(contigs))} sequences")


@T("dunovo_duplex_consensus", "Duplex consensus from paired reads", DUNOVO, "table",
   [fq("reads1", READS_1, "R1 FASTQ"), fq("reads2", READS_2, "R2 FASTQ"), intin("min_overlap", 8,
                                                            "Minimum overlap", min=2),
    number("min_quality", 20.0, "Base quality cut", min=0.0)],
   ex={"reads1": READS_1, "reads2": READS_2, "min_overlap": 6}, up="Duplex Consensus (du novo)",
   tags=("consensus", "duplex", "reads"),
   summary="Merge overlapping read pairs into one consensus with per-read quality stats.")
def dunovo_duplex_consensus(reads1, reads2, min_overlap=8, min_quality=20.0):
    """Pair overlap consensus."""
    a = _reads(reads1)
    b = _reads(reads2)
    res = dict(align.merge_pairs([{"id": r.id, "seq": r.seq, "qual": r.qual} for r in a],
                                [{"id": r.id, "seq": r.seq, "qual": r.qual} for r in b],
                                min_overlap=int(min_overlap), max_mismatches=2))
    rows = []
    for i, m in enumerate(res.get("merged") or []):
        d = dict(m)
        s_, q_ = str(d.get("seq", "")), str(d.get("qual", ""))
        qv = np.array([max(0, ord(c) - 33) for c in q_], dtype=float) if q_ else np.zeros(1)
        rows.append({"read": d.get("id", f"duplex_{i + 1}"), "length": len(s_),
                    "overlap": res.get("mean_overlap", 0), "mean_quality": round(float(qv.mean()), 3),
                    "low_quality_bases": int(np.sum(qv < float(min_quality))),
                    "both_strands": bool(q_.strip())})
    out = table(rows, f"{res.get('n_merged', 0)} duplex consensus reads from {len(a)} + {len(b)} reads "
                    f"(merge rate {res.get('merge_rate', 0)}%)")
    out["stats"] = {"n_unmerged": res.get("n_unmerged", 0), "mean_overlap": res.get("mean_overlap", 0.0),
                   "overlap_mismatches": res.get("total_mismatches_in_overlap", 0),
                   "fasta": io.write_fasta([io.Seq(id=str(dict(m).get('id', f'duplex_{i + 1}')),
                                                  seq=str(dict(m).get('seq', '')))
                                           for i, m in enumerate(res.get("merged") or [])])}
    return out


# ===========================================================================
# epigenetics / methylation
# ===========================================================================
@T("tnseq_hgc_flux", "Gene-level insertion counts (HGC flux)", TN, "table",
   [sam("alignment", ALN_SAM, "Transposon insertion sites"), gff("annotation", ANNOT_GFF, "GFF3 genes"),
    intin("min_insertions", 0, "Report genes with at least n insertions", min=0)],
   ex={"alignment": ALN_SAM, "annotation": ANNOT_GFF}, up="TRANSIT-PROC / ESME flux",
   tags=("TnSeq", "fitness"),
   summary="Normalise insertion counts per gene to obtain a simple flux estimate.")
def tnseq_hgc_flux(alignment, annotation, min_insertions=0):
    """HGC-style flux."""
    _h, alns = io.parse_sam(io.as_text(alignment))
    sites = Counter()
    for a in alns:
        if a.rname:
            sites[(a.rname, a.pos // 1000)] += 1
    genes = [g for g in _gff(annotation) if g.type == "gene"]
    counts = {}
    for g in genes:
        n = sum(v for (c, kb), v in sites.items() if c == g.seqid and g.start // 1000 <= kb <= g.end // 1000)
        counts[g.attrs.get("ID", g.attrs.get("gene_id", f"{g.seqid}:{g.start}")) or f"gene_{g.start}"] = n
    hgc = {k: max(1, int(v / 3) + 1) for k, v in counts.items()}
    tot_reads = max(1, sum(counts.values()))
    rows = []
    for k, n in counts.items():
        exp = hgc[k] / max(1, sum(hgc.values())) * tot_reads
        rows.append({"gene": k, "insertions": n, "hgc": hgc[k], "expected": round(exp, 3),
                    "log2_ratio": round(math.log2((n + 0.5) / (exp + 0.5)), 4),
                    "fitness": round(float(np.log2((n + 0.5) / (exp + 0.5))), 4)})
    rows = [r for r in rows if r["insertions"] >= int(min_insertions)]
    rows.sort(key=lambda r: r["log2_ratio"])
    return table(rows, f"{len(rows)} genes, {tot_reads} insertions placed")


@T("tnseq_log2fc", "Log2 fold change between insertion libraries", TN, "table",
   [tbl("counts_a", COUNTS, "Condition A gene counts"), tbl("counts_b", COUNTS, "Condition B gene counts"),
    number("pseudocount", 1.0, "Pseudocount", min=0.0), number("alpha", 0.05, "Significance threshold", min=0.0)],
   ex={"counts_a": COUNTS, "counts_b": COUNTS}, up="ESME / TRANSIT", tags=("TnSeq", "fold change"),
   summary="Normalised log2 fold changes with per-gene significance for screens.")
def tnseq_log2fc(counts_a, counts_b, pseudocount=1.0, alpha=0.05):
    """Screen log2FC."""
    na, ca, Xa = _mat(counts_a)
    nb, cb, Xb = _mat(counts_b)
    if Xa.size == 0 or Xb.size == 0:
        return table([], "both tables need data")
    va = np.nansum(Xa, axis=1).astype(float)
    vb = np.nansum(Xb, axis=1).astype(float)
    n = min(len(va), len(vb))
    va, vb = va[:n], vb[:n]
    fa = va / max(1e-9, va.sum())
    fb = vb / max(1e-9, vb.sum())
    pc = float(pseudocount)
    lfc = np.log2((fb + pc) / (va / max(1e-9, va.sum()) + pc))
    ps = []
    rows = []
    for i in range(n):
        orr, p = stats.fisher_exact(int(round(vb[i])), int(round(vb.sum() - vb[i])),
                                   int(round(va[i])), int(round(va.sum() - va[i])))
        ps.append(p)
        rows.append({"gene": na[i] if i < len(na) else f"gene_{i + 1}", "count_a": int(va[i]),
                    "count_b": int(vb[i]), "log2_fold_change": round(float(lfc[i]), 4),
                    "p_value": round(float(p), 6), "odds_ratio": round(float(orr), 4)})
    adj = stats.p_adjust(ps, "fdr_bh")
    for i, r in enumerate(rows):
        r["p_adjusted"] = round(float(adj[i]), 6) if i < len(adj) else 1.0
        r["phenotype"] = ("attenuated" if r["log2_fold_change"] < -1 else
                         "hypertrophic" if r["log2_fold_change"] > 1 else "neutral")
        r["significant"] = bool(r["p_adjusted"] < float(alpha))
    rows.sort(key=lambda r: r["log2_fold_change"])
    return table(rows, f"{sum(1 for r in rows if r['significant'])} significant genes of {len(rows)}")


@T("tnseq_ess", "Essentiality score per gene", TN, "table",
   [tbl("counts_a", COUNTS, "T0 insertion counts"), tbl("counts_b", COUNTS, "End-point counts"),
    number("cut_low", -1.0, "Essential cut-off", min=-10.0), number("cut_high", 1.0, "Fitness cut-off", max=10.0)],
   ex={"counts_a": COUNTS, "counts_b": COUNTS}, up="TRANSIT-PROC / BaDIeS", tags=("essentiality", "TnSeq"),
   summary="Combine fold change and insertion count into an essentiality classification.")
def tnseq_ess(counts_a, counts_b, cut_low=-1.0, cut_high=1.0):
    """Essentiality calls."""
    na, _ca, Xa = _mat(counts_a)
    _nb, _cb, Xb = _mat(counts_b)
    if Xa.size == 0:
        return table([], "no counts")
    va = np.nansum(Xa, axis=1).astype(float)
    vb = np.nansum(Xb, axis=1).astype(float) if Xb.size else va
    n = min(len(va), len(vb))
    va, vb = va[:n], vb[:n]
    rows = []
    for i in range(n):
        lfc = math.log2((vb[i] + 1) / (va[i] + 1))
        z = float((lfc - np.mean([math.log2((vb[j] + 1) / (va[j] + 1)) for j in range(n)]))
                  / (np.std([math.log2((vb[j] + 1) / (va[j] + 1)) for j in range(n)]) or 1.0))
        rows.append({"gene": na[i] if i < len(na) else f"gene_{i + 1}", "log2_fold_change": round(lfc, 4),
                    "ess_score": round(z, 4), "t0_counts": int(va[i]), "endpoint_counts": int(vb[i]),
                    "classification": ("essential" if lfc <= float(cut_low) else
                                      "advantageous" if lfc >= float(cut_high) else "neutral")})
    rows.sort(key=lambda r: r["ess_score"])
    return table(rows, f"{sum(1 for r in rows if r['classification'] == 'essential')} essential genes")


@T("mimodd_mirna_hairpins", "miRNA hairpin detection", MIM, "table",
   [fa("records", GENES_FA, "Precursor sequences"), number("min_mfe", -3.0, "Minimum folding energy", max=0.0),
    intin("min_stem", 8, "Minimum stem length", min=4), number("loop_min", 3, "Minimum loop", min=1)],
   ex={"records": GENES_FA, "min_mfe": -1.0}, up="MiModD / RNAfold", tags=("miRNA", "hairpin", "RNA"),
   summary="Stem-loop prediction with free energy and mature-arm extraction.")
def mimodd_mirna_hairpins(records, min_mfe=-3.0, min_stem=8, loop_min=3):
    """Hairpin screening."""
    rows = []
    for r in _fa(records):
        s = seq.clean(r.seq.upper()).replace("U", "T")
        for w in range(int(min_stem) * 2 + int(loop_min), min(len(s), int(min_stem) * 2 + 40) + 1):
            for i in range(0, max(1, len(s) - w + 1)):
                win = s[i:i + w]
                d = dict(seq.hairpin_score(win, loop_min=int(loop_min), max_stem=int(min_stem) + 4))
                dg = float(d.get("delta_g", 0.0) or 0.0)
                stem = int(d.get("stem", 0) or 0)
                if dg <= float(min_mfe) and stem >= int(min_stem):
                    start = i + int(d.get("start", 0) or 0)
                    rows.append({"sequence": r.id, "start": start + 1, "end": start + stem * 2 + 1,
                                "length": w, "stem": stem, "loop": int(d.get("loop", 0) or 0),
                                "delta_g": round(dg, 3), "hairpin_score": round(float(d.get("score", 0.0) or 0.0), 3),
                                "five_p_arm": win[:stem], "three_p_arm": seq.revcomp(win[len(win) - stem:])})
                    break
            if rows and rows[-1]["sequence"] == r.id:
                break
    return table(rows, f"{len(rows)} hairpins below {min_mfe} kcal/mol")


@T("gemini_annotate_regions", "Annotate variants falling in features", GEM, "table",
   [vcf("variants", VARIANTS_VCF, "VCF file"), gff("annotation", ANNOT_GFF, "GFF3 annotation"),
    choice("max_impact", ["high", "medium", "low", "modifier"], "modifier", "Report impacts up to"),
    boolean("gene_summary", False, "Aggregate per gene")],
   ex={"variants": VARIANTS_VCF, "annotation": ANNOT_GFF}, up="gemini annotate / VEP",
   tags=("variant annotation", "GEMINI"),
   summary="Intersect a VCF with genes and summarise the coding consequences.")
def gemini_annotate_regions(variants, annotation, max_impact="modifier", gene_summary=False):
    """Gemini-like annotation."""
    vcfobj = _vcf(variants)
    genes = [g for g in _gff(annotation) if g.type == "gene"]
    order = {"HIGH": 0, "MODERATE": 1, "LOW": 2, "MODIFIER": 3}
    limit = order.get(str(max_impact).upper(), 3)
    rows = []
    for rec in vcfobj.records:
        near = None
        inside = None
        for g in genes:
            if g.seqid != rec.chrom:
                continue
            if g.start <= rec.pos <= g.end:
                inside = g
                break
            if abs((g.start + g.end) // 2 - rec.pos) < 500:
                near = near or g
        g = inside or near
        name = (g.attrs.get("ID") or g.attrs.get("gene_id") or f"{g.seqid}:{g.start}") if g else ""
        impact = "HIGH" if inside and rec.pos == inside.start else (
            "MODERATE" if inside else "LOW" if near else "MODIFIER")
        rows.append({"chrom": rec.chrom, "pos": rec.pos, "id": rec.id, "ref": rec.ref,
                    "alt": str(rec.alt).replace(",", "/"), "gene": name,
                    "effect": "coding" if inside else ("upstream" if near else "intergenic"),
                    "impact": impact, "quality": rec.qual, "filter": rec.filter,
                    "is_snv": bool(getattr(rec, "is_snv", False))})
    rows = [r for r in rows if order.get(r["impact"], 3) <= limit]
    if gene_summary:
        agg = Counter(r["gene"] for r in rows if r["gene"])
        return table([{"gene": k, "variants": v} for k, v in agg.most_common()],
                     f"{len(rows)} variants across {len(agg)} genes")
    return table(rows, f"{len(rows)} variants annotated against {len(genes)} genes")


@T("gemini_chrom_state", "Coverage-based chromatin state per interval", GEM, "table",
   [anyfile("graph", BEDGRAPH, fmt="", label="Signal track"), bed("regions", REGIONS_BED, "Intervals"),
    number("low", 0.2, "Low threshold", min=0.0), number("high", 0.8, "High threshold", min=0.0)],
   ex={"graph": BEDGRAPH, "regions": REGIONS_BED}, up="chromHMM-lite / gemini chrom_state",
   tags=("chromatin", "signal"),
   summary="Bucket regions into repressed, poised and active by their signal level.")
def gemini_chrom_state(graph, regions, low=0.2, high=0.8):
    """Discrete chromatin states."""
    vals = io.parse_bedgraph(io.as_text(graph))
    by_chrom: dict[str, list[tuple[int, int, float]]] = defaultdict(list)
    for c, s, e, v in vals:
        by_chrom[c].append((int(s), int(e), float(v)))
    rows = []
    for iv in _bed(regions):
        seg = [v for s, e, v in by_chrom.get(iv.chrom, []) if e > iv.start and s < iv.end]
        m = float(np.mean(seg)) if seg else 0.0
        rows.append({"chrom": iv.chrom, "start": iv.start, "end": iv.end, "name": iv.name or ".",
                    "signal": round(m, 4),
                    "state": ("A_quiescent" if m <= float(low) else "C_active" if m >= float(high)
                             else "B_poised"),
                    "segments": len(seg)})
    counts = Counter(r["state"] for r in rows)
    return table(rows, "states: " + ", ".join(f"{k}={v}" for k, v in sorted(counts.items())))


@T("emboss_geece", "Coding potential of each DNA strand", EMB, "table",
   [fa("records", GENES_FA, "Nucleotide FASTA"), intin("min_len", 100, "Minimum ORF length", min=10),
    choice("genetic_code", ["Standard", "Vertebrate Mitochondrial"], "Standard", "Genetic code")],
   ex={"records": GENES_FA, "min_len": 60}, up="EMBOSS geece", tags=("coding density", "ORF"),
   summary="Fraction of each sequence covered by ORFs on the forward and reverse strands.")
def emboss_geece(records, min_len=100, genetic_code="Standard"):
    """G+T content / coding density."""
    rows = []
    for r in _fa(records):
        fwd = sum(int(dict(o).get("length", 0) or 0) for o in
                 seq.find_orfs(r.seq, min_len=max(10, int(min_len) // 3), table=genetic_code))
        rev = sum(int(dict(o).get("length", 0) or 0) for o in
                 seq.find_orfs(seq.revcomp(r.seq), min_len=max(10, int(min_len) // 3), table=genetic_code))
        L = max(1, len(r.seq))
        rows.append({"sequence": r.id, "length": L,
                    "coding_density_forward": round(fwd / L, 4), "coding_density_reverse": round(rev / L, 4),
                    "coding_density_combined": round(min(L, fwd + rev) / L, 4),
                    "gt_content": round(float((r.seq.upper().count("G") + r.seq.upper().count("T")) / L), 4),
                    "gc_content": round(100 * float(seq.gc_content(r.seq)), 3)})
    return table(rows, f"{len(rows)} sequences scanned for coding density")


@T("emboss_newcleft", "Find restriction enzymes that leave no cut", EMB, "table",
   [fa("records", GENES_FA, "Nucleotide FASTA"), intin("min_size", 200, "Minimum fragment size", min=20),
    intin("max_enzymes", 25, "Enzymes to test", min=1)],
   ex={"records": GENES_FA, "min_size": 100}, up="EMBOSS newcleft", tags=("restriction enzyme", "cloning"),
   summary="Report enzymes whose cuts leave at least one very large fragment.")
def emboss_newcleft(records, min_size=200, max_enzymes=25):
    """Cleft site search."""
    enz = [str(d.get("name", "")) for d in seq.enzyme_table()[: int(max_enzymes)]]
    rows = []
    for r in _fa(records):
        for e in enz:
            if not e:
                continue
            sizes = seq.fragment_sizes(r.seq, enzyme=e)
            if not sizes:
                continue
            rows.append({"sequence": r.id, "enzyme": e, "n_fragments": len(sizes),
                        "longest_fragment": int(max(sizes)), "mean_fragment": round(float(np.mean(sizes)), 1),
                        "leaves_large_fragment": bool(max(sizes) >= int(min_size))})
    rows = [x for x in rows if x["leaves_large_fragment"]]
    rows.sort(key=lambda d: -d["longest_fragment"])
    return table(rows[:300], f"{len(rows)} enzyme/sequence combinations leave a fragment >= {min_size} bp")


@T("emboss_transeq", "Translate nucleotides in six frames", EMB, "text",
   [fa("records", GENES_FA, "Nucleotide FASTA"), choice("frame", ["1", "2", "3", "-1", "-2", "-3", "all"], "all",
                                                      "Frame"), boolean("to_stop", False, "Stop at first stop codon"),
    intin("min_len", 30, "Minimum translation length", min=1)],
   ex={"records": GENES_FA, "frame": "1", "min_len": 20}, up="EMBOSS transeq", tags=("translation", "ORF"),
   summary="Six-frame translation with a length filter, as a FASTA result.")
def emboss_transeq(records, frame="all", to_stop=False, min_len=30):
    """Six-frame translation."""
    out = []
    for r in _fa(records):
        frames = [0, 1, 2, 3, 4, 5] if frame == "all" else [[0, 1, 2, 3, 4, 5][["1", "2", "3", "-1", "-2",
                                                                              "-3"].index(frame)]]
        for f in frames:
            real = f % 3
            src = r.seq if f < 3 else seq.revcomp(r.seq)
            prot = seq.translate(src[real:], frame=0, to_stop=bool(to_stop))
            prot = str(prot).replace("*", "")
            if len(prot) >= int(min_len):
                out.append(io.Seq(id=f"{r.id}_{f + 1}", seq=prot, desc=f"frame={f + 1};length={len(prot)}"))
    return text(io.write_fasta(out), f"{len(out)} translations from {len(_fa(records))} sequences")


@T("emboss_einverted", "Find inverted (palindromic) repeats", EMB, "table",
   [fa("records", GENES_FA, "Nucleotide FASTA"), intin("threshold", 8, "Minimum repeat length", min=3),
    intin("period", 100000, "Maximum period", min=2), intin("head", 200, "Rows", min=1)],
   ex={"records": GENES_FA, "threshold": 6}, up="EMBOSS einverted", tags=("palindrome", "repeats"),
   summary="Direct and inverted repeats within a sequence, with their positions.")
def emboss_einverted(records, threshold=8, period=100000, head=200):
    """Inverted repeat finder."""
    rows = []
    for r in _fa(records):
        s = r.seq.upper()
        found = seq.find_repeats(s, min_len=int(threshold)) if hasattr(seq, "find_repeats") else []
        for d in found:
            dd = dict(d)
            a, b = int(dd.get("start", 0)), int(dd.get("end", dd.get("stop", 0)))
            if b <= a or b - a < int(threshold) or (b - a) > int(period):
                continue
            left = s[a:b]
            right = s[max(0, b): b + (b - a)]
            rows.append({"sequence": r.id, "start": a + 1, "end": b, "length": b - a,
                        "type": dd.get("type", "inverted" if right and left == seq.revcomp(right) else "direct"),
                        "period": b - a, "identity": round(100.0 * sum(1 for x, y in zip(left, right)
                                                                     if x == y) / max(1, len(left)), 2),
                        "match": left[:40]})
    return table(rows[: int(head)], f"{len(rows)} repeats of at least {threshold} bp")


@T("emboss_getorf", "Extract open reading frames", EMB, "text",
   [fa("records", GENES_FA, "Nucleotide FASTA"), intin("min_len", 90, "Minimum nucleotide length", min=9),
    choice("frame", ["all", "1", "2", "3"], "all", "Frame"),
    textbox("required_start", "ATG", "Required start codon"), boolean("partial", True, "Include partials")],
   ex={"records": GENES_FA, "min_len": 60, "required_start": ""}, up="EMBOSS getorf", tags=("ORF", "genes"),
   summary="All ORFs above a size threshold, written out as a nucleotide FASTA.")
def emboss_getorf(records, min_len=90, frame="all", required_start="ATG", partial=True):
    """ORF extraction."""
    out = []
    for r in _fa(records):
        frames = [1, 2, 3] if frame == "all" else [int(frame)]
        for f in frames:
            src = r.seq[f - 1:]
            want = max(1, int(min_len) // 3)
            for o in seq.find_orfs(src, min_len=want, required_start=required_start or "ATG",
                                 partial=bool(partial)):
                d = dict(o)
                if int(d.get("frame", f) or f) != f:
                    continue
                if str(d.get("strand", "+")) != "+":
                    continue
                nt = int(d.get("length", 0) or 0)
                if nt < int(min_len):
                    continue
                out.append(io.Seq(id=f"{r.id}_{d.get('id', len(out) + 1)}", seq=str(d.get("sequence",
                             r.seq[int(d.get("start", 0)):int(d.get("end", 0))])),
                             desc=f"frame={f};length={nt};aa_length={d.get('aa_length', '')};"
                                  f"protein={str(d.get('protein', ''))[:60]}"))
    return text(io.write_fasta(out), f"{len(out)} ORFs of at least {min_len} nt from "
                                    f"{len(_fa(records))} sequences")


@T("emboss_pepstats", "Protein physicochemical statistics", EMB, "table",
   [fa("proteins", PROTEINS, "Protein FASTA"), number("pH", 7.0, "pH for charge", min=0.0),
    choice("mode", ["average", "mono"], "average", "Mass mode")],
   ex={"proteins": PROTEINS, "pH": 7.0}, up="EMBOSS pepstats / ProtParam", tags=("protein properties", "pI"),
   summary="Molecular weight, pI, charge, instability and hydrophobicity per chain.")
def emboss_pepstats(proteins, pH=7.0, mode="average"):
    """pepstats-style table."""
    rows = []
    for r in _fa(proteins):
        s = protein.clean(r.seq)
        if not s:
            continue
        comp = protein.aa_composition(s)
        rows.append({"sequence": r.id, "length": len(s), "molecular_weight": round(float(protein.molecular_weight(s, mode=mode)), 2),
                    "isoelectric_point": round(float(protein.isoelectric_point(s)), 3),
                    "net_charge_at_pH": round(float(protein.net_charge(s, pH=float(pH))), 3),
                    "gravy": round(float(protein.gravy(s)), 4),
                    "instability_index": round(float(protein.instability_index(s)), 3),
                    "aromaticity": round(float(protein.aromaticity(s)), 4),
                    "cysteines": s.count("C"), "extinction_coefficient_reduced":
                    round(float(dict(protein.extinction_coefficient(s)).get("extinction_reduced", 0.0)), 2),
                    "alanine_percent": round(comp.get("A", 0.0), 3)})
    return table(rows, f"{len(rows)} protein chains characterised")


@T("emboss_water", "Local pairwise alignment of two sequences", EMB, "table",
   [textbox("sequence_a", SHORT_DNA, "Sequence A"), textbox("sequence_b", SHORT_DNA, "Sequence B"),
    choice("mode", ["local", "global"], "local", "Algorithm"), number("gap_open", -10.0, "Gap open", max=0.0),
    number("gap_extend", -0.5, "Gap extend", max=0.0), choice("matrix", ["", "BLOSUM62", "PAM250"], "",
                                                             "Scoring matrix")],
   ex={"sequence_a": "ACGTACGTAGCTAGCTAGCTAGCA", "sequence_b": "ACGTACGTAGCTAGCTAGCTAGCAT"},
   up="EMBOSS water / needle", tags=("pairwise alignment", " Smith-Waterman"),
   summary="Optimal alignment with identity, gaps and the raw score.")
def emboss_water(sequence_a, sequence_b, mode="local", gap_open=-10.0, gap_extend=-0.5, matrix=""):
    """Pairwise alignment."""
    a = str(sequence_a or "").strip()
    b = str(sequence_b or "").strip()
    if not a or not b:
        return table([], "two sequences are required")
    st = align.align_pairwise(a, b, mode=mode, matrix=matrix or "", gap_open=float(gap_open),
                            gap_extend=float(gap_extend))
    d = dict(st)
    rows = [{"metric": k, "value": (round(float(v), 4) if isinstance(v, float) else str(v)[:400])}
            for k, v in d.items() if k not in ("aligned_a", "aligned_b", "aln_a", "aln_b", "cigar")]
    ident = float(d.get("identity", d.get("percent_identity", 0.0)) or 0.0)
    res = table(rows, f"{mode} alignment of {len(a)} x {len(b)}: score {d.get('score', 'n/a')}, "
                    f"identity {round(ident, 2)}%")
    res["stats"]["alignment"] = {k: str(v)[:200] for k, v in d.items() if "aln" in k or "align" in k}
    return res


@T("bbduk_trim", "Quality and adapter trimming (kmer filter)", BB, "text",
   [fq("reads", READS_SINGLE, "FASTQ reads"), number("qual_cutoff", 20.0, "Quality cut-off", min=0.0),
    intin("min_len", 25, "Minimum read length", min=1), textbox("adapters", "AGATCGGAAGAGC",
                                                              "Adapter sequences (comma separated)"),
    boolean("trim_qlow", True, "Trim low-quality tails")],
   ex={"reads": READS_1, "qual_cutoff": 20.0, "min_len": 30}, up="BBMap bbduk.sh", tags=("trimming", "FASTQ"),
   summary="Trim adapter matches and poor-quality tails, dropping short leftovers.")
def bbduk_trim(reads, qual_cutoff=20.0, min_len=25, adapters="", trim_qlow=True):
    """bbduk-like trimming."""
    recs = _reads(reads)
    pats = [a.strip().upper() for a in str(adapters).split(",") if a.strip()]
    cut = float(qual_cutoff)
    out = []
    dropped = Counter()
    adapter_hits = 0
    for r in recs:
        s, q = r.seq.upper(), r.qual
        trimmed_adapter = 0
        for p in pats:
            i = s.find(p)
            if i >= 0:
                s, q = s[:i], q[:i]
                trimmed_adapter += 1
                adapter_hits += 1
        if trim_qlow and q:
            k = len(q)
            while k > 0 and ord(q[k - 1]) - 33 < cut:
                k -= 1
            s, q = s[:k], q[:k]
        if len(s) < int(min_len):
            dropped["too_short" if not trimmed_adapter else "adapter_and_short"] += 1
            continue
        out.append(io.Read(id=r.id, seq=s, qual=q or "I" * len(s), desc=r.desc))
    return text(io.write_fastq(out), f"{len(out)} reads kept of {len(recs)} "
                                    f"({dict(dropped)} dropped, {adapter_hits} adapter trims)")


@T("bbrep_repearness", "Repearness of each read", BB, "table",
   [fq("reads", READS_SINGLE, "FASTQ reads"), intin("kmer", 6, "k-mer size", min=2),
    number("threshold", 0.3, "Repeat fraction cut-off", min=0.0), intin("head", 200, "Rows", min=1)],
   ex={"reads": READS_1, "kmer": 5}, up="BBMap bbrep", tags=("repeats", "complexity"),
   summary="Fraction of a read made of repeated k-mers, for filtering low-complexity data.")
def bbrep_repearness(reads, kmer=6, threshold=0.3, head=200):
    """Repearness scoring."""
    recs = _reads(reads)
    k = max(2, int(kmer))
    rows = []
    for r in recs:
        cnt = seq.kmer_counts(r.seq.upper(), k=k, alphabet="ACGT", min_count=1)
        tot = sum(cnt.values()) or 1
        repeated = sum(v for v in cnt.values() if v > 1)
        rows.append({"read": r.id, "length": len(r.seq), "distinct_kmers": len(cnt),
                    "total_kmers": tot, "repeat_fraction": round(repeated / tot, 4),
                    "repetitive": bool(repeated / tot >= float(threshold)),
                    "dust_score": round(float(seq.dust_score(r.seq.upper())), 3)})
    n = sum(1 for r in rows if r["repetitive"])
    res = table(rows[: int(head)], f"{n}/{len(rows)} reads are >= {threshold} repetitive (k={k})")
    res["stats"]["figure"] = plot.histogram([r["repeat_fraction"] for r in rows], bins=20,
                                          xlabel="repeat fraction", title="Read repearness")
    return res
