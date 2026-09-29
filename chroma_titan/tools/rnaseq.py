"""RNA-seq quantification and differential expression, single-cell and spatial omics,
Hi-C and multi-omics integration.

Panel sections covered: **rna_seq**, **expression_tools**, **single_cell**,
**spatial**, **hicexplorer**, **rseqc**, **iwtomics**, **multiomics**,
**graph_display_data** (a couple), **presto**.
"""
from __future__ import annotations

import math
import re
from collections import Counter, defaultdict

import numpy as np
import pandas as pd

from chroma_titan.core import bam, genome, io, ml, motif, plot, seq, stats, tables
from chroma_titan.tools._common import *  # noqa: F401,F403
from chroma_titan.tools._common import (ALN_SAM, ANNOT_GFF, BEDGRAPH, COUNTS, GENES_FA, GENOME,
                                         GMT, PAIRS, PHENO, REGIONS_BED)

RNA = "rna_seq"
EXPR = "expression_tools"
SC = "single_cell"
SPATIAL = "spatial"
HIC = "hicexplorer"
RSEQC = "rseqc"
IWT = "iwtomics"
MULTI = "multiomics"
GDD = "graph_display_data"


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------
def _matrix(src):
    """Counts table -> (numeric DataFrame indexed by feature, sample columns)."""
    df = tables.load(io.as_text(src))
    if df.empty:
        return df, []
    cols = list(df.columns)
    idc = cols[0] if str(cols[0]).lower() in ("gene", "id", "feature", "gene_id", "cell",
                                             "barcode", "name") else cols[0]
    body = df.drop(columns=[idc]) if idc in df.columns else df
    num = {}
    for c in body.columns:
        vals = []
        for v in body[c].tolist():
            try:
                vals.append(float(v))
            except (TypeError, ValueError):
                vals.append(float("nan"))
        num[str(c)] = vals
    out = pd.DataFrame(num, index=[str(x) for x in df[idc].tolist()])
    return out, list(out.columns)


def _bed(src):
    return io.parse_bed(io.as_text(src)) if str(src or "").strip() else []


def _groups(df, spec: str) -> dict[str, list[str]]:
    """Parse ``A=s1,s2;B=s3,s4``; empty spec splits the columns in half."""
    cols = list(df.columns)
    out: dict[str, list[str]] = {}
    for part in str(spec or "").replace(" ", "").split(";"):
        if not part:
            continue
        if "=" in part:
            name, ids = part.split("=", 1)
            out[name or "group"] = [x for x in ids.split(",") if x in cols]
        else:
            out[part] = [x for x in cols if x == part]
    if not out:
        half = max(1, len(cols) // 2)
        out = {"control": cols[:half], "treatment": cols[half:]}
    return {k: v for k, v in out.items() if v}


def _load_pheno(src, id_col: str = "") -> dict[str, str]:
    """phenotype table -> {sample: group label}."""
    df = tables.load(io.as_text(src))
    if df.empty:
        return {}
    cols = list(df.columns)
    ids = id_col or cols[0]
    lab = next((c for c in cols if re.search(r"group|label|condition|treatment|class|phenotype",
                                            str(c), re.I)), cols[-1])
    out = {}
    for r in df.to_dict("records"):
        out[str(r.get(ids))] = str(r.get(lab))
    return out


# ===========================================================================
# expression matrices
# ===========================================================================
@T("counts_summary", "Summary of a counts matrix", EXPR, "table",
   [tbl("src", COUNTS, "Counts matrix"), intin("top", 12, "Rows", min=1)],
   ex={"src": COUNTS}, up="HTSeq-count / featureCounts output QC",
   tags=("expression", "quality control"),
   summary="Library sizes, genes detected, zero counts and the most abundant features.")
def counts_summary(src, top=20):
    """Per-sample library statistics."""
    df, samples = _matrix(src)
    if df.empty:
        return table([], "empty counts matrix")
    rows = []
    for s in samples:
        v = df[s].to_numpy(dtype=float)
        rows.append({"sample": s, "total_counts": int(np.nansum(v)), "features": int(len(v)),
                    "detected": int(np.sum(v > 0)), "zeros": int(np.sum(v == 0)),
                    "mean": round(float(np.nanmean(v)), 2), "median": round(float(np.nanmedian(v)), 2),
                    "max": int(np.nanmax(v)) if np.isfinite(v).any() else 0,
                    "top_feature": df.index[int(np.nanargmax(v))] if np.isfinite(v).any() else ""})
    return table(rows, f"{len(df)} features x {len(samples)} samples")


@T("counts_cpm", "Normalise counts to CPM", EXPR, "table",
   [tbl("src", COUNTS, "Counts matrix"), boolean("log", False, "log2(1+CPM)"),
    intin("digits", 4, "Decimals", min=1)],
   ex={"src": COUNTS, "log": True}, up="edgeR cpm / DESeq2 vst",
   tags=("normalisation", "expression"),
   summary="Library-size normalisation to counts per million, optionally log transformed.")
def counts_cpm(src, log=False, digits=4):
    """CPM normalisation."""
    df, samples = _matrix(src)
    out = tables.cpm(df)
    if log:
        out = tables.log_transform(out, base=2, pseudocount=1)
    rows = [{"feature": i, **{c: round(float(out.loc[i, c]), int(digits)) for c in samples}}
            for i in out.index]
    return {"message": f"CPM matrix ({'log2 ' if log else ''}{len(rows)} x {len(samples)})",
            "rows": rows[:500]}


@T("counts_tpm", "Normalise to TPM using transcript lengths", EXPR, "table",
   [tbl("src", COUNTS, "Counts matrix"), fa("lengths", GENES_FA, "Transcript FASTA (lengths)"),
    intin("digits", 2, "Decimals", min=1)],
   ex={"src": COUNTS, "lengths": GENES_FA}, up="salmon / featureCounts + TPM",
   tags=("normalisation", "TPM", "expression"),
   summary="Length-corrected transcripts per million with a FASTA of transcript sequences.")
def counts_tpm(src, lengths, digits=2):
    """TPM normalisation from a FASTA of transcripts."""
    import pandas as pd

    df, samples = _matrix(src)
    fa = {r.id: len(r.seq) for r in io.parse_fasta(io.as_text(lengths))}
    if not fa:
        return table([], "no transcript lengths parsed")
    lser = pd.Series({i: float(fa.get(i, max(1, int(np.nanmean(list(fa.values())))))) for i in df.index})
    eff = df.astype(float).T / lser
    out = (eff / eff.sum() * 1e6).T
    rows = [{"feature": i, **{c: round(float(out.loc[i, c]), int(digits)) for c in samples}} for i in out.index]
    return {"message": f"TPM for {len(rows)} features ({len(fa)} lengths used)", "rows": rows[:500]}


@T("counts_quantile_normalize", "Quantile normalisation of a matrix", EXPR, "table",
   [tbl("src", COUNTS, "Counts matrix"), intin("digits", 3, "Decimals", min=1)],
   ex={"src": COUNTS}, up="limma normalizeQuantiles", tags=("normalisation", "array"),
   summary="Force every sample onto the same distribution before comparing groups.")
def counts_quantile_normalize(src, digits=3):
    """Quantile normalisation."""
    df, samples = _matrix(src)
    out = tables.quantile_normalize(df)
    rows = [{"feature": i, **{c: round(float(out.loc[i, c]), int(digits)) for c in samples}}
            for i in out.index]
    return {"message": f"quantile-normalised {len(rows)} x {len(samples)}", "rows": rows[:500]}


@T("counts_zscore", "Row z-scores of a matrix", EXPR, "table",
   [tbl("src", COUNTS, "Counts matrix"), intin("axis", 1, "Axis (1 = per sample, 0 = per feature)", min=0, max=1),
    intin("digits", 3, "Decimals", min=1)],
   ex={"src": COUNTS, "axis": 0}, up="pheatmap / row z-score", tags=("transformation", "expression"),
   summary="Centre and scale each row (or column) so patterns are visible in a heatmap.")
def counts_zscore(src, axis=1, digits=3):
    """z-score transform."""
    df, samples = _matrix(src)
    out = tables.zscore(df, axis=int(axis))
    rows = [{"feature": i, **{c: round(float(out.loc[i, c]), int(digits)) for c in samples}}
            for i in out.index]
    return {"message": f"z-scores (axis {axis})", "rows": rows[:500]}


@T("counts_log_transform", "Log transform a matrix", EXPR, "table",
   [tbl("src", COUNTS, "Counts matrix"), choice("base", ["2", "10", "e"], "2", "Log base"),
    intin("pseudocount", 1, "Pseudocount", min=0)],
   ex={"src": COUNTS, "base": "2"}, up="log2 transform (limma)", tags=("transformation", "expression"),
   summary="Variance-stabilising log transform with a configurable pseudocount.")
def counts_log_transform(src, base="2", pseudocount=1):
    """Log transform."""
    df, samples = _matrix(src)
    b = {"2": 2.0, "10": 10.0, "e": math.e}[base]
    out = tables.log_transform(df, base=b, pseudocount=int(pseudocount))
    rows = [{"feature": i, **{c: round(float(out.loc[i, c]), 4) for c in samples}} for i in out.index]
    return {"message": f"log{base} of {len(rows)} features", "rows": rows[:500]}


@T("counts_fold_change", "Fold change between two columns", EXPR, "table",
   [tbl("src", COUNTS, "Counts matrix"), textbox("column_a", "", "Column A (or group spec A=s1,s2)"),
    textbox("column_b", "", "Column B"), number("pseudocount", 1.0, "Pseudocount", min=0.0),
    number("min_abs_log2fc", 0.0, "Filter |log2FC| above", min=0.0)],
   ex={"src": COUNTS, "column_a": "sample_A", "column_b": "sample_D"},
   up="DESeq2 results / edgeR", tags=("fold change", "expression"),
   summary="Per-feature log2 fold change between two samples, columns or sample groups.")
def counts_fold_change(src, column_a="", column_b="", pseudocount=1.0, min_abs_log2fc=0.0):
    """Fold change between columns or groups."""
    df, samples = _matrix(src)
    gspec = _groups(df, f"{column_a}" if "=" in str(column_a) else "")
    if column_a and column_b and column_a in samples and column_b in samples:
        a, b = df[column_a].to_numpy(dtype=float), df[column_b].to_numpy(dtype=float)
        names = [f"{column_a} vs {column_b}"]
    else:
        ga = [s for s in samples if s in str(column_a).split(",")] or samples[: max(1, len(samples) // 2)]
        gb = [s for s in samples if s in str(column_b).split(",")] or samples[len(samples) // 2:]
        a = df[ga].mean(axis=1).to_numpy(dtype=float)
        b = df[gb].mean(axis=1).to_numpy(dtype=float)
        names = [f"{'+'.join(ga)} vs {'+'.join(gb)}"]
    rows = []
    pc = float(pseudocount)
    for i, (x, y) in enumerate(zip(a, b)):
        fc = math.log2((x + pc) / (y + pc)) if (y + pc) > 0 and (x + pc) > 0 else 0.0
        if abs(fc) < float(min_abs_log2fc):
            continue
        rows.append({"feature": df.index[i], "value_a": round(float(x), 3), "value_b": round(float(y), 3),
                    "log2_fold_change": round(fc, 5), "fold_change": round((x + pc) / (y + pc), 4)
                    if (y + pc) > 0 else "", "contrast": names[0]})
    rows.sort(key=lambda r: -abs(r["log2_fold_change"]))
    return table(rows[:500], f"{len(rows)} features (contrast {names[0]})")


@T("counts_sample_correlation", "Sample-to-sample correlation matrix", EXPR, "figure",
   [tbl("src", COUNTS, "Counts matrix"), choice("method", ["pearson", "spearman"], "pearson", "Method"),
    boolean("log", True, "Log2 transform first")],
   ex={"src": COUNTS, "method": "spearman"}, up="DESeq2 sampleDistance / corrplot",
   tags=("correlation", "samples", "heatmap"),
   summary="Heatmap of pairwise correlations between libraries, the batch-effect smoke test.")
def counts_sample_correlation(src, method="pearson", log=True):
    """Sample correlation heatmap."""
    df, samples = _matrix(src)
    mat = tables.log_transform(df, base=2, pseudocount=1) if log else df
    X = mat.to_numpy(dtype=float)
    if X.shape[1] < 2:
        return plot.empty_plot("need at least two samples")
    if method == "spearman":
        ranks = np.apply_along_axis(lambda c: np.argsort(np.argsort(c)).astype(float), 0, X)
        C = np.corrcoef(ranks.T)
    else:
        C = np.corrcoef(X.T)
    return plot.heatmap([[round(float(v), 4) for v in row] for row in C], row_labels=samples,
                       col_labels=samples, title=f"{method} sample correlation", annot=True)


@T("counts_pca_samples", "PCA of samples", EXPR, "figure",
   [tbl("src", COUNTS, "Counts matrix"), intin("components", 2, "Components", min=2, max=4),
    choice("scale", ["none", "zscore", "log"], "log", "Pre-processing")],
   ex={"src": COUNTS, "scale": "zscore"}, up="DESeq2 plotPCA / prcomp",
   tags=("PCA", "dimensionality reduction", "samples"),
   summary="Project libraries onto principal components to reveal replicate structure.")
def counts_pca_samples(src, components=2, scale="log"):
    """Sample PCA scatter."""
    df, samples = _matrix(src)
    mat = df
    if scale == "log":
        mat = tables.log_transform(mat, base=2, pseudocount=1)
    elif scale == "zscore":
        mat = tables.zscore(tables.log_transform(mat, base=2, pseudocount=1), axis=0)
    res = ml.pca(mat.to_numpy(dtype=float), components=int(components))
    sc = np.asarray(res["scores"], dtype=float)
    ratio = np.asarray(res.get("explained_variance_ratio") or [], dtype=float)
    return plot.scatter(list(sc[:, 0]), list(sc[:, 1]), labels=samples,
                       xlabel=f"PC1 ({round(100 * (ratio[0] if len(ratio) else 0), 1)}%)",
                       ylabel=f"PC2 ({round(100 * (ratio[1] if len(ratio) > 1 else 0), 1)}%)",
                       title="Sample PCA")


@T("counts_clustering_heatmap", "Clustered heatmap of features", EXPR, "figure",
   [tbl("src", COUNTS, "Counts matrix"), intin("top_features", 20, "Most variable features", min=3),
    choice("method", ["average", "complete", "single", "ward"], "average", "Linkage"),
    boolean("zscore_rows", True, "Row z-score")],
   ex={"src": COUNTS, "top_features": 10}, up="pheatmap / ComplexHeatmap",
   tags=("clustering", "heatmap", "expression"),
   summary="Order rows by hierarchical clustering of the most variable genes and draw the matrix.")
def counts_clustering_heatmap(src, top_features=20, method="average", zscore_rows=True):
    """Clustered expression heatmap."""
    df, samples = _matrix(src)
    mat = tables.zscore(tables.log_transform(df, base=2, pseudocount=1), axis=0)
    var = mat.var(axis=1).sort_values(ascending=False)
    keep = list(var.index[: int(top_features)])
    sub = mat.loc[keep]
    X = sub.to_numpy(dtype=float)
    Z = ml.linkage(X, method=method, labels=keep)
    order = _leaf_order(Z)
    m = [[round(float(v), 3) for v in sub.loc[keep[i]].tolist()] for i in order] if Z else \
        [[round(float(v), 3) for v in row] for row in X]
    labels = [keep[i] for i in order] if Z else keep
    if not zscore_rows:
        m = [[round(float(v), 3) for v in df.loc[l].tolist()] for l in labels]
    return plot.heatmap(m, row_labels=labels, col_labels=samples,
                       title=f"top {len(labels)} variable features ({method} linkage)")


def _leaf_order(Z):
    """Leaf ordering of a scipy-style linkage matrix without scipy."""
    Z = np.asarray(Z, dtype=float)
    n = len(Z) + 1
    children = {i: [i] for i in range(n)}
    for k, row in enumerate(Z):
        a, b = int(row[0]), int(row[1])
        children[n + k] = children.get(a, [a]) + children.get(b, [b])
    return children.get(n + len(Z) - 1, list(range(n)))


@T("counts_dedupe_merge", "Merge replicate columns", EXPR, "table",
   [tbl("src", COUNTS, "Counts matrix"), textbox("suffix", "_rep", "Replicate suffix to merge"),
    choice("agg", ["mean", "sum", "median", "max"], "mean", "Aggregation"),
    intin("digits", 3, "Decimals", min=1)],
   ex={"src": COUNTS, "suffix": "_rep"}, up="collate / sum columns",
   tags=("aggregation", "expression"),
   summary="Collapse sample columns sharing a prefix (replicates) into one value per feature.")
def counts_dedupe_merge(src, suffix="_rep", agg="mean", digits=3):
    """Merge replicate columns by name prefix."""
    df, samples = _matrix(src)
    groups: dict[str, list[str]] = defaultdict(list)
    for s in samples:
        base = re.sub(r"[_-]?\d+$", "", str(s)) if not suffix else re.sub(re.escape(suffix) + r"\d*$", "", str(s))
        groups[base or s].append(s)
    fn = {"mean": np.nanmean, "sum": np.nansum, "median": np.nanmedian, "max": np.nanmax}[agg]
    rows = []
    for name, cols in groups.items():
        vals = df[cols].to_numpy(dtype=float)
        col = fn(vals, axis=1)
        rows.append({"feature": "", "column": name, "merged_from": ",".join(cols),
                    **{str(i): round(float(v), int(digits)) for i, v in zip(df.index, col)}})
    rows[0]["feature"] = f"{len(groups)} merged columns via {agg}"
    return table(rows, f"{len(samples)} columns merged into {len(groups)}")


@T("counts_filter_low_expression", "Filter low-count features", EXPR, "table",
   [tbl("src", COUNTS, "Counts matrix"), intin("min_counts", 10, "Minimum counts per sample", min=0),
    intin("min_samples", 1, "Minimum number of samples", min=1),
    number("min_cpm", 0.0, "Minimum CPM", min=0.0)],
   ex={"src": COUNTS, "min_counts": 1500, "min_samples": 2}, up="DESeq2 row filtering",
   tags=("filtering", "expression"),
   summary="Drop features that are not expressed well enough for statistical testing.")
def counts_filter_low_expression(src, min_counts=10, min_samples=1, min_cpm=0.0):
    """Low-count filtering."""
    df, samples = _matrix(src)
    keep, dropped = [], []
    lib = df.to_numpy(dtype=float).sum(axis=0)
    for i, row in df.iterrows():
        v = row.to_numpy(dtype=float)
        n_ok = int(np.sum(v >= float(min_counts)))
        cpm_ok = int(np.sum(1e6 * v / np.where(lib > 0, lib, 1) >= float(min_cpm)))
        (keep if n_ok >= int(min_samples) and cpm_ok >= int(min_samples) else dropped).append(i)
    sub = df.loc[keep]
    rows = [{"feature": i, **{c: float(sub.loc[i, c]) for c in samples}} for i in keep]
    return table(rows, f"{len(keep)} features kept, {len(dropped)} filtered out")


@T("counts_tidy_long", "Long format of an expression matrix", EXPR, "table",
   [tbl("src", COUNTS, "Counts matrix"), textbox("id_name", "feature", "Name of the id column"),
    intin("head", 500, "Rows", min=1)],
   ex={"src": COUNTS}, up="tidyr pivot_longer", tags=("tidy", "transformation"),
   summary="Unroll the matrix into feature / sample / value rows for plotting tools.")
def counts_tidy_long(src, id_name="feature", head=500):
    """Melt to long format."""
    df, samples = _matrix(src)
    rows = []
    for i in df.index:
        for c in samples:
            rows.append({id_name: i, "sample": c, "value": float(df.loc[i, c])})
    return table(rows[:int(head)], f"{len(rows)} observations")


@T("counts_add_calc_column", "Add a computed column to a table", EXPR, "table",
   [tbl("src", COUNTS, "Table"), textbox("name", "ratio", "New column name"),
    textbox("expression", "(sample_A + sample_B) / sample_D", "pandas-style expression"),
    intin("digits", 4, "Decimals", min=1)],
   ex={"src": COUNTS, "name": "ratio", "expression": "sample_A / sample_D"},
   up="Galaxy Add/redo calculation", tags=("transformation", "columns"),
   summary="Evaluate an arithmetic expression over the numeric columns of a table.")
def counts_add_calc_column(src, name="value", expression="", digits=4):
    """Column arithmetic on a table."""
    df = tables.load(io.as_text(src))
    out = tables.add_column(df, name, expression) if expression else df
    rows = out.to_dict("records")
    for r in rows:
        if isinstance(r.get(name), float):
            r[name] = round(r[name], int(digits))
    return table(rows[:500], f"column '{name}' added to {len(rows)} rows")


@T("counts_feature_stats", "Descriptive statistics per feature", EXPR, "table",
   [tbl("src", COUNTS, "Counts matrix"), intin("head", 200, "Rows", min=1),
    boolean("cv_sort", True, "Sort by coefficient of variation")],
   ex={"src": COUNTS}, up="row stats (base R)", tags=("statistics", "expression"),
   summary="Mean, median, spread, CV and detection rate for every gene.")
def counts_feature_stats(src, head=200, cv_sort=True):
    """Per-row summary."""
    df, samples = _matrix(src)
    rows = []
    for i in df.index:
        v = df.loc[i].to_numpy(dtype=float)
        mu = float(np.nanmean(v))
        sd = float(np.nanstd(v, ddof=1)) if len(v) > 1 else 0.0
        rows.append({"feature": i, "mean": round(mu, 3), "median": round(float(np.nanmedian(v)), 3),
                    "min": float(np.nanmin(v)), "max": float(np.nanmax(v)),
                    "stdev": round(sd, 3), "cv": round(sd / mu, 5) if mu > 0 else 0.0,
                    "n_samples": int(len(v)), "detected": int(np.sum(v > 0))})
    if cv_sort:
        rows.sort(key=lambda r: -r["cv"])
    return table(rows[:int(head)], f"{len(rows)} features summarised")


# ===========================================================================
# differential expression
# ===========================================================================
@T("deseq_lite", "Differential expression (negative-binomial-free DESeq2 style)", RNA, "table",
   [tbl("src", COUNTS, "Counts matrix"), textbox("groups", "", "Group spec, e.g. control=sample_A,sample_B;trt=sample_D,sample_E"),
    txt("phenotype_file", "", "Phenotype table (optional)"),
    choice("test", ["wald", "welch", "exact"], "wald", "Test"),
    choice("adjust", ["fdr_bh", "bonferroni", "none"], "fdr_bh", "p adjustment"),
    number("alpha", 0.05, "Significance threshold", min=0.0, max=1.0)],
   ex={"src": COUNTS, "groups": "control=sample_A,sample_B,sample_C;trt=sample_D,sample_E,sample_F"},
   up="DESeq2 / edgeR", tags=("differential expression", "statistics", "RNA-seq"),
   summary="Per-gene log2 fold change, test statistic, raw and adjusted p-values.")
def deseq_lite(src, groups="", phenotype_file="", test="wald", adjust="fdr_bh", alpha=0.05):
    """Differential expression on a counts table."""
    df, samples = _matrix(src)
    if phenotype_file:
        lab = _load_pheno(phenotype_file)
        gspec = ";".join(f"{g}={','.join([s for s in samples if lab.get(s) == g])}"
                        for g in sorted(set(lab.values())))
        grp = _groups(df, gspec)
    else:
        grp = _groups(df, groups)
    names = list(grp)
    if len(names) < 2:
        return table([], "need two groups (e.g. control=...;treatment=...)")
    g1, g2 = grp[names[0]], grp[names[1]]
    lib1 = df[g1].to_numpy(dtype=float).sum(axis=0)
    lib2 = df[g2].to_numpy(dtype=float).sum(axis=0)
    norm1 = df[g1].to_numpy(dtype=float) / np.where(lib1 > 0, lib1, 1) * np.mean(lib1)
    norm2 = df[g2].to_numpy(dtype=float) / np.where(lib2 > 0, lib2, 1) * np.mean(lib1)
    rows, ps = [], []
    for k, gene in enumerate(df.index):
        a, b = norm1[k], norm2[k]
        pc = 1.0
        lfc = math.log2((np.nanmean(b) + pc) / (np.nanmean(a) + pc))
        if test == "exact":
            kk = int(round(np.nansum(b)))
            n_tot = int(round(np.nansum(a) + np.nansum(b)))
            p, stat = float(stats.binom_two_sided(kk, max(1, n_tot), 0.5)), 0.0
        else:
            st = stats.ttest_ind(list(b), list(a), equal_var=(test != "welch"))
            p, stat = float(st["p_value"]), float(st["statistic"])
        if not math.isfinite(p):
            p = 1.0
        ps.append(p)
        rows.append({"feature": gene, "base_mean": round(float(np.nansum(df.loc[gene].to_numpy(dtype=float)))
                                                     / max(1, len(samples)), 3),
                    "mean_control": round(float(np.nanmean(a)), 3), "mean_treatment": round(float(np.nanmean(b)), 3),
                    "log2_fold_change": round(lfc, 5), "statistic": round(stat, 4),
                    "p_value": p, "contrast": f"{names[1]} vs {names[0]}"})
    adj = stats.p_adjust(ps, adjust) if adjust != "none" else list(ps)
    for i, r in enumerate(rows):
        q = float(adj[i]) if i < len(adj) else 1.0
        r["p_adjusted"] = round(q, 6)
        r["significant"] = bool(q < float(alpha) and abs(r["log2_fold_change"]) > 0.5)
    rows.sort(key=lambda r: r["p_value"])
    n_sig = sum(1 for r in rows if r["significant"])
    res = table(rows, f"{n_sig}/{len(rows)} features significant at FDR < {alpha}")
    res["stats"] = {"test": test, "groups": {k: v for k, v in grp.items()},
                   "n_significant": n_sig, "alpha": float(alpha)}
    return res


@T("limma_moderated_t", "limma-style moderated t-test", RNA, "table",
   [tbl("src", COUNTS, "Log-expression matrix"), textbox("groups", "", "Group spec"),
    number("prior_df", 3.0, "Empirical Bayes prior df", min=0.5),
    number("alpha", 0.05, "FDR threshold", min=0.0, max=1.0)],
   ex={"src": COUNTS, "groups": "control=sample_A,sample_B,sample_C;trt=sample_D,sample_E,sample_F"},
   up="limma eBayes", tags=("differential expression", "moderated statistics"),
   summary="Shrink the gene-wise variance toward a pooled estimate and recompute t statistics.")
def limma_moderated_t(src, groups="", prior_df=3.0, alpha=0.05):
    """Moderated t-test."""
    df, samples = _matrix(src)
    grp = _groups(df, groups)
    names = list(grp)
    if len(names) < 2:
        return table([], "need two groups")
    g1, g2 = grp[names[0]], grp[names[1]]
    M = tables.log_transform(df, base=2, pseudocount=1).to_numpy(dtype=float)
    d0 = float(prior_df)
    rows, ps = [], []
    for k, gene in enumerate(df.index):
        a, b = M[k][: len(g1)] if False else np.nan_to_num(M[k, [samples.index(s) for s in g1]]), \
            np.nan_to_num(M[k, [samples.index(s) for s in g2]])
        if len(a) < 2 or len(b) < 2:
            continue
        va, vb = float(np.var(a, ddof=1)), float(np.var(b, ddof=1))
        sp = ((len(a) - 1) * va + (len(b) - 1) * vb) / max(1, len(a) + len(b) - 2)
        s0 = sp * d0 / max(1e-9, (d0 + 1))  # pooled prior variance
        se = math.sqrt(max(1e-9, (sp + s0) * (1 / len(a) + 1 / len(b))))
        tstat = (float(np.mean(b)) - float(np.mean(a))) / se
        dfree = d0 + len(a) + len(b) - 2
        p = stats.t_pvalue(abs(tstat), dfree)
        ps.append(p)
        rows.append({"feature": gene, "log2_fold_change": round(float(np.mean(b) - np.mean(a)), 4),
                    "t": round(tstat, 4), "df": round(dfree, 2),
                    "p_value": p if math.isfinite(p) else 1.0,
                    "prior_df": d0, "residual_sd": round(math.sqrt(sp), 4)})
    adj = stats.p_adjust([q if math.isfinite(q) else 1.0 for q in ps], "fdr_bh")
    for i, r in enumerate(rows):
        q = float(adj[i]) if i < len(adj) else 1.0
        r["p_adjusted"] = round(q, 6)
        r["significant"] = bool(q < alpha)
    rows.sort(key=lambda r: r["p_value"])
    return table(rows, f"{sum(1 for r in rows if r['significant'])} genes at FDR < {alpha}")


@T("volcano_plot", "Volcano plot of differential results", RNA, "figure",
   [tbl("results", COUNTS, "Results table"), textbox("log2fc_column", "", "log2FC column (empty = auto)"),
    textbox("p_column", "", "p-value column (empty = auto)"),
    number("lfc_cut", 1.0, "log2FC cut", min=0.0), number("p_cut", 0.05, "p cut", min=0.0),
    textbox("label_column", "", "Column with point labels")],
   ex={"results": COUNTS, "log2fc_column": "sample_A", "p_column": "sample_D", "lfc_cut": 0.2, "p_cut": 0.6},
   up="EnhancedVolcano / DESeq2 results plot", tags=("plot", "volcano", "differential expression"),
   summary="Plot log2 fold change against -log10 p-value with significance thresholds.")
def volcano_plot(results, log2fc_column="", p_column="", lfc_cut=1.0, p_cut=0.05, label_column=""):
    """Volcano plot from a results table."""
    df = tables.load(io.as_text(results))
    cols = [c for c in df.columns]
    fc_col = log2fc_column or next((c for c in cols if "log2" in c.lower() or "fold" in c.lower()), cols[1] if len(cols) > 1 else cols[0])
    p_col = p_column or next((c for c in cols if c.lower() in ("p_value", "pvalue", "p.adj", "padj", "q_value")),
                            cols[-1])
    lab = label_column or next((c for c in cols if c.lower() in ("feature", "gene", "id")), "")
    x, y, names = [], [], []
    for r in df.to_dict("records"):
        try:
            f = float(r[fc_col])
            p = float(r[p_col])
        except (KeyError, TypeError, ValueError):
            continue
        if not (math.isfinite(f) and math.isfinite(p)):
            continue
        x.append(f)
        y.append(-math.log10(max(p, 1e-300)))
        names.append(str(r.get(lab, "")) if lab else "")
    if not x:
        return plot.empty_plot("no numeric fold change / p-value columns found")
    return plot.volcano(x, y, labels=names, title="Volcano plot", lfc_cut=float(lfc_cut),
                       p_cut=float(p_cut))


@T("ma_plot", "MA plot (fold change versus mean expression)", RNA, "figure",
   [tbl("results", COUNTS, "Results table"), textbox("mean_column", "", "Mean expression column"),
    textbox("log2fc_column", "", "log2FC column")],
   ex={"results": COUNTS, "mean_column": "sample_A", "log2fc_column": "sample_D"},
   up="edgeR plotMD / limma", tags=("plot", "MA", "differential expression"),
   summary="Detect intensity-dependent bias by plotting fold change against abundance.")
def ma_plot(results, mean_column="", log2fc_column=""):
    """MA plot from a results table."""
    df = tables.load(io.as_text(results))
    cols = list(df.columns)
    m_col = mean_column or next((c for c in cols if "mean" in c.lower() or "base" in c.lower()), cols[1] if len(cols) > 1 else cols[0])
    a_col = log2fc_column or next((c for c in cols if "log2" in c.lower() or "fold" in c.lower()), cols[-1])
    xs, ys = [], []
    for r in df.to_dict("records"):
        try:
            xs.append(math.log10(max(1e-9, float(r[m_col]))))
            ys.append(float(r[a_col]))
        except (KeyError, TypeError, ValueError):
            continue
    if not xs:
        return plot.empty_plot("no usable columns")
    return plot.ma_plot(xs, ys, title="MA plot")


@T("gene_set_counts", "Aggregate expression over GMT gene sets", RNA, "table",
   [tbl("src", COUNTS, "Counts matrix"), txt("sets", GMT, "GMT gene sets"),
    choice("agg", ["mean", "sum", "median", "max"], "mean", "Aggregation"),
    boolean("zscore", True, "z-score the profile across samples")],
   ex={"src": COUNTS, "sets": GMT}, up="GSVA / roasters", tags=("pathway", "gene sets", "GMT"),
   summary="Score each pathway by summarising the expression of its member genes.")
def gene_set_counts(src, sets, agg="mean", zscore=True):
    """Gene-set activity scores."""
    df, samples = _matrix(src)
    gsets = io.parse_gmt(io.as_text(sets))
    fn = {"mean": np.nanmean, "sum": np.nansum, "median": np.nanmedian, "max": np.nanmax}[agg]
    base = tables.log_transform(df, base=2, pseudocount=1)
    rows = []
    for name, members in gsets.items():
        present = [m for m in dict.fromkeys(members) if m in base.index]
        if not present:
            rows.append({"gene_set": name, "genes_found": 0, "genes_total": len(set(members))})
            continue
        prof = fn(base.loc[present].to_numpy(dtype=float), axis=0)
        rec = {"gene_set": name, "genes_found": len(present), "genes_total": len(set(members)),
              "coverage": round(len(present) / max(1, len(set(members))), 3)}
        for s, v in zip(samples, prof):
            rec[s] = round(float(v), 4)
        if zscore and np.nanstd(prof) > 0:
            zs = (prof - np.nanmean(prof)) / np.nanstd(prof)
            rec["profile_z"] = ",".join(f"{v:.3f}" for v in zs)
        rows.append(rec)
    return table(rows, f"{len(rows)} gene sets scored")


@T("gsea_lite", "Rank-based GSEA of a gene set", RNA, "table",
   [tbl("results", COUNTS, "Ranked results table"), textbox("rank_column", "", "Ranking column"),
    txt("sets", GMT, "GMT gene sets"), intin("permutations", 200, "Permutations", min=10),
    number("min_size", 2, "Minimum set size", min=1)],
   ex={"results": COUNTS, "rank_column": "sample_A", "permutations": 100},
   up="GSEA / fgsea", tags=("pathway", "GSEA", "enrichment"),
   summary="Running-sum enrichment score with a permutation null for every gene set.")
def gsea_lite(results, rank_column="", sets=GMT, permutations=200, min_size=2):
    """GSEA-lite on a ranked list."""
    df = tables.load(io.as_text(results))
    cols = list(df.columns)
    rc = rank_column or next((c for c in cols if "log2" in c.lower() or "stat" in c.lower() or "fold" in c.lower()),
                            cols[-1])
    lab = next((c for c in cols if c.lower() in ("feature", "gene", "id")), cols[0])
    pairs = []
    for r in df.to_dict("records"):
        try:
            pairs.append((str(r[lab]), float(r[rc])))
        except (KeyError, TypeError, ValueError):
            continue
    if not pairs:
        return table([], "no numeric ranking column found")
    pairs.sort(key=lambda x: -x[1])
    names = [n for n, _v in pairs]
    values = np.array([v for _n, v in pairs], dtype=float)
    gsets = io.parse_gmt(io.as_text(sets))
    rng = np.random.default_rng(1)
    rows = []
    for gname, members in gsets.items():
        mem = [m for m in dict.fromkeys(members) if m in names]
        if len(mem) < int(min_size):
            continue
        idx = np.array([names.index(m) for m in mem])
        es, nes = _running_score(idx, len(names), np.abs(values))
        null = []
        for _ in range(int(permutations)):
            p = rng.permutation(len(names))[: len(idx)]
            null.append(_running_score(np.sort(p), len(names), np.abs(values))[0])
        sd = float(np.std(null)) or 1e-9
        rows.append({"gene_set": gname, "size": len(mem), "enrichment_score": round(es, 5),
                    "normalised_score": round(nes, 5), "mean_null": round(float(np.mean(null)), 5),
                    "p_value": round(float((1 + sum(1 for v in null if abs(v) >= abs(es)))
                                          / (1 + len(null))), 5),
                    "leading_edge": ",".join(names[i] for i in idx[:6])})
    rows.sort(key=lambda r: -abs(r["normalised_score"]))
    return table(rows, f"{len(rows)} sets tested with {permutations} permutations")


def _running_score(idx, n, weights):
    """Weighted running-sum enrichment score for the hit positions in idx."""
    hit = np.zeros(n, dtype=float)
    ii = np.asarray(idx, dtype=int)
    hit[ii] = 1.0
    w = np.asarray(weights, dtype=float)[:n]
    p_plus = np.cumsum(hit * w) / max(1e-9, float(np.sum(hit * w)))
    p_zero = np.cumsum(1.0 - hit) / max(1, n - len(ii))
    diff = p_plus - p_zero
    max_pos, min_neg = float(np.max(diff)), float(np.min(diff))
    es = max_pos if abs(max_pos) >= abs(min_neg) else min_neg
    return es, es / max(1e-9, math.sqrt(max(1, len(ii))))


@T("differential_fisher_on_off", "Fisher test of presence/absence per feature", RNA, "table",
   [tbl("src", COUNTS, "Counts matrix"), textbox("groups", "", "Group spec"),
    intin("cutoff", 1, "Count above which a sample counts as on", min=0)],
   ex={"src": COUNTS, "groups": "control=sample_A,sample_B,sample_C;trt=sample_D,sample_E,sample_F",
       "cutoff": 2000}, up="Fisher exact on/off (edgeR-ancestral)",
   tags=("differential expression", "Fisher"),
   summary="Test on/off patterns between groups when only detection matters.")
def differential_fisher_on_off(src, groups="", cutoff=1):
    """Presence/absence Fisher test."""
    df, samples = _matrix(src)
    grp = _groups(df, groups)
    names = list(grp)
    if len(names) < 2:
        return table([], "need two groups")
    g1, g2 = grp[names[0]], grp[names[1]]
    rows = []
    for i in df.index:
        a = int(sum(1 for s in g1 if float(df.loc[i, s]) > float(cutoff)))
        b = int(sum(1 for s in g2 if float(df.loc[i, s]) > float(cutoff)))
        orr, p = stats.fisher_exact(a, len(g1) - a, b, len(g2) - b)
        rows.append({"feature": i, f"on_{names[0]}": a, f"on_{names[1]}": b,
                    "odds_ratio": round(float(orr), 4) if math.isfinite(orr) else "inf",
                    "p_value": p})
    rows.sort(key=lambda r: r["p_value"])
    return table(rows, f"{len(rows)} features tested (count > {cutoff} = on)")


@T("permutation_group_test", "Permutation test between two sample groups", RNA, "table",
   [tbl("src", COUNTS, "Counts matrix"), textbox("groups", "", "Group spec"),
    intin("permutations", 500, "Permutations", min=20), choice("stat", ["mean_difference", "t", "sum"],
                                                              "mean_difference", "Statistic"),
    intin("seed", 1, "Seed", min=0)],
   ex={"src": COUNTS, "groups": "control=sample_A,sample_B;trt=sample_D,sample_E", "permutations": 200},
   up="permutation test / signal", tags=("statistics", "permutation"),
   summary="Distribution-based p-values for group contrasts when sample numbers are tiny.")
def permutation_group_test(src, groups="", permutations=500, stat="mean_difference", seed=1):
    """Global permutation test on the matrix."""
    df, samples = _matrix(src)
    grp = _groups(df, groups)
    names = list(grp)
    if len(names) < 2:
        return table([], "need two groups")
    A, B = grp[names[0]], grp[names[1]]
    ia = [samples.index(s) for s in A]
    ib = [samples.index(s) for s in B]
    X = df.to_numpy(dtype=float)

    def score(a_idx, b_idx):
        if stat == "t":
            out = []
            for k in range(X.shape[0]):
                st = stats.ttest_ind(list(X[k, b_idx]), list(X[k, a_idx]))
                out.append(abs(float(st["statistic"])))
            return float(np.sum(out))
        if stat == "sum":
            return float(abs(np.sum(X[:, b_idx]) - np.sum(X[:, a_idx])))
        return float(np.sum(np.abs(np.nanmean(X[:, b_idx], axis=1) - np.nanmean(X[:, a_idx], axis=1))))

    obs = score(ia, ib)
    all_idx = list(range(len(samples)))
    rng = np.random.default_rng(int(seed))
    null = []
    for _ in range(int(permutations)):
        perm = rng.permutation(all_idx)
        null.append(score(perm[: len(ia)], perm[len(ia):len(ia) + len(ib)]))
    p = (1 + sum(1 for v in null if v >= obs)) / (1 + len(null))
    rows = [{"metric": "observed_statistic", "value": round(obs, 5)},
            {"metric": "null_mean", "value": round(float(np.mean(null)), 5)},
            {"metric": "null_sd", "value": round(float(np.std(null)), 5)},
            {"metric": "permutations", "value": int(permutations)},
            {"metric": "p_value", "value": round(p, 5)},
            {"metric": "groups", "value": f"{names[0]} (n={len(A)}) vs {names[1]} (n={len(B)})"}]
    res = table(rows, f"permutation p = {round(p, 4)} ({stat})")
    res["stats"] = {"null_histogram": dict(Counter(round(v / 100) * 100 for v in null))}
    return res


@T("iwtonomics_interval_test", "IWTomics-style interval test on profiles", IWT, "table",
   [tbl("src", COUNTS, "Matrix of profiles (features x samples)"), textbox("groups", "", "Group spec"),
    choice("stat", ["mean", "max", "sum"], "mean", "Aggregation"), intin("neighbourhood", 2,
                                                                       "Smoothing window (features)", min=1),
    intin("permutations", 200, "Permutations", min=10)],
   ex={"src": COUNTS, "groups": "control=sample_A,sample_B,sample_C;trt=sample_D,sample_E,sample_F"},
   up="IWTomics (iwtomics)", tags=("wavelet", "interval testing"),
   summary="Smooth the profile along the feature order and find intervals differing between groups.")
def iwtonomics_interval_test(src, groups="", stat="mean", neighbourhood=2, permutations=200):
    """Interval-level group testing along ordered features."""
    df, samples = _matrix(src)
    grp = _groups(df, groups)
    names = list(grp)
    if len(names) < 2:
        return table([], "need two groups")
    A, B = [samples.index(s) for s in grp[names[0]]], [samples.index(s) for s in grp[names[1]]]
    X = df.to_numpy(dtype=float)
    agg = {"mean": np.nanmean, "max": np.nanmax, "sum": np.nansum}[stat]
    a = agg(X[:, A], axis=1)
    b = agg(X[:, B], axis=1)
    d = b - a
    w = int(neighbourhood)
    sm = np.array([float(np.mean(d[max(0, i - w):i + w + 1])) for i in range(len(d))])
    rng = np.random.default_rng(3)
    null = []
    for _ in range(int(permutations)):
        p = rng.permutation(len(samples))
        aa, bb = p[: len(A)], p[len(A):len(A) + len(B)]
        dd = agg(X[:, bb], axis=1) - agg(X[:, aa], axis=1)
        s2 = np.array([float(np.mean(dd[max(0, i - w):i + w + 1])) for i in range(len(dd))])
        null.append(float(np.max(np.abs(s2))))
    rows = []
    for i in range(len(sm)):
        p = (1 + sum(1 for v in null if abs(sm[i]) <= v)) / (1 + len(null))
        if abs(sm[i]) < 1e-6:
            continue
        rows.append({"interval_start_feature": i + 1, "interval_end_feature": i + 1,
                    "window": f"{max(1, i - w + 1)}-{min(len(sm), i + w + 1)}",
                    "mean_difference": round(float(sm[i]), 5),
                    "raw_difference": round(float(d[i]), 5), "p_value": round(float(p), 5),
                    "direction": "up" if sm[i] > 0 else "down"})
    rows.sort(key=lambda r: r["p_value"])
    return table(rows[:200], f"{len(rows)} intervals, {sum(1 for r in rows if r['p_value'] < 0.05)} with p<0.05")


@T("rseqc_gene_body_coverage", "Gene body coverage profile", RSEQC, "figure",
   [sam("alignments", ALN_SAM, "Alignments (SAM)"), gff("annotation", ANNOT_GFF, "GFF3 annotation"),
    intin("bins", 20, "Body bins", min=5)],
   ex={"alignments": ALN_SAM, "annotation": ANNOT_GFF}, up="RSeQC gene_body_coverage",
   tags=("RNA-seq", "quality control", "coverage"),
   summary="Read density along the normalised gene body: 3' bias shows up immediately.")
def rseqc_gene_body_coverage(alignments, annotation, bins=20):
    """Metagene coverage along transcripts."""
    hdr, alns = bam.load(alignments)
    genes = genome.genes_from_gff(io.parse_gff(io.as_text(annotation)))
    nb = max(5, int(bins))
    acc = np.zeros(nb)
    ngenes = 0
    for g in genes:
        hits = [a for a in alns if a.rname == g.chrom and a.mapped and g.start <= a.pos - 1 < g.end]
        if not hits:
            continue
        ngenes += 1
        L = max(1, g.end - g.start)
        for a in hits:
            frac = min(0.9999, max(0.0, (a.pos - 1 - g.start) / L))
            acc[min(nb - 1, int(frac * nb))] += 1
    if ngenes == 0:
        return plot.empty_plot("no reads inside annotated genes")
    return plot.line({"coverage": [float(v) for v in acc]}, x=list(range(1, nb + 1)),
                     xlabel=f"gene body ({nb} equal bins)", ylabel="reads",
                     title=f"Gene body coverage over {ngenes} genes")


@T("rseqc_junction_annotation", "Splice junctions from alignments and annotation", RSEQC, "table",
   [sam("alignments", ALN_SAM, "Alignments (SAM)"), gff("annotation", ANNOT_GFF, "GFF3 annotation"),
    intin("min_intron", 20, "Minimum intron length", min=1), intin("head", 200, "Rows", min=1)],
   ex={"alignments": ALN_SAM, "annotation": ANNOT_GFF}, up="RSeQC junction_annotation",
   tags=("RNA-seq", "junctions", "splicing"),
   summary="Infer N-gaps in CIGARs (or introns between blocks) and check them against the annotation.")
def rseqc_junction_annotation(alignments, annotation, min_intron=20, head=200):
    """Junction detection and annotation."""
    hdr, alns = bam.load(alignments)
    rows_gff = io.parse_gff(io.as_text(annotation))
    known = set()
    exons: dict[str, list] = defaultdict(list)
    for r in rows_gff:
        if r.type in ("exon", "CDS"):
            exons[str(r.seqid)].append((r.start, r.end))
    for chrom, lst in exons.items():
        lst.sort()
        for i in range(len(lst) - 1):
            known.add((chrom, lst[i][1] + 1, lst[i + 1][0] - 1))
    found = Counter()
    for a in alns:
        if not a.mapped:
            continue
        pos = a.pos - 1
        for n, op in io.cigar_ops(a.cigar):
            if op in "MDN=X":
                if op == "N" and n >= int(min_intron):
                    found[(a.rname, pos + 1, pos + n)] += 1
                pos += n
    rows = []
    for (chrom, s, e), n in found.most_common():
        rows.append({"junction": f"{chrom}:{s}-{e}", "chrom": chrom, "start": s, "end": e,
                    "intron_length": e - s + 1, "reads": n,
                    "known": bool((chrom, s, e) in known),
                    "canonical": ""})
    n_known = sum(1 for r in rows if r["known"])
    return table(rows[:int(head)], f"{len(rows)} junctions ({n_known} matching annotated introns)")


@T("rseqc_read_distribution", "Where reads fall relative to features", RSEQC, "table",
   [sam("alignments", ALN_SAM, "Alignments (SAM)"), gff("annotation", ANNOT_GFF, "GFF3 annotation"),
    bed("transposons", REGIONS_BED, "Other feature BED (optional)"), intin("head", 100, "Rows", min=1)],
   ex={"alignments": ALN_SAM, "annotation": ANNOT_GFF}, up="RSeQC read_distribution",
   tags=("RNA-seq", "quality control", "features"),
   summary="Split alignments into CDS, exon-intron, upstream and downstream categories.")
def rseqc_read_distribution(alignments, annotation, transposons="", head=100):
    """Read distribution over feature classes."""
    hdr, alns = bam.load(alignments)
    rows_gff = io.parse_gff(io.as_text(annotation))
    cds = [r for r in rows_gff if r.type == "CDS"]
    ex = [r for r in rows_gff if r.type == "exon"]
    tx = _bed(transposons) if transposons else []
    cats = Counter()
    per_read = []
    for a in alns:
        if not a.mapped:
            cats["unmapped"] += 1
            continue
        s, e = a.pos - 1, a.ref_end
        in_cds = any(r.seqid == a.rname and not (e <= r.start - 1 or s >= r.end) for r in cds)
        in_exon = any(r.seqid == a.rname and not (e <= r.start - 1 or s >= r.end) for r in ex)
        in_other = any(iv.chrom == a.rname and not (e <= iv.start or s >= iv.end) for iv in tx)
        if in_cds:
            cats["CDS_exons"] += 1
        elif in_exon:
            cats["exons_non_coding"] += 1
        elif in_other:
            cats["other_features"] += 1
        else:
            cats["introns_or_intergenic"] += 1
        per_read.append({"read": a.qname, "reference": a.rname, "position": a.pos,
                        "in_CDS": in_cds, "in_exon": in_exon, "other_feature": in_other})
    total = max(1, sum(cats.values()))
    rows = [{"category": k, "reads": v, "percent": round(100 * v / total, 3)}
            for k, v in sorted(cats.items(), key=lambda kv: -kv[1])]
    res = table(rows, f"{len(alns)} records classified")
    res["stats"]["example_reads"] = per_read[: int(head)]
    return res


@T("rseqc_inner_distance", "Inner distance / insert size distribution", RSEQC, "figure",
   [sam("alignments", ALN_SAM, "Alignments (SAM)"), intin("bins", 25, "Bins", min=5),
    intin("max_size", 500, "Truncate at", min=10)],
   ex={"alignments": ALN_SAM}, up="RSeQC inner_distance", tags=("RNA-seq", "insert size", "plot"),
   summary="Fragment inner-distance histogram - the Fritsch-Hill criterion for strandedness.")
def rseqc_inner_distance(alignments, bins=25, max_size=500):
    """Inner distance histogram."""
    hdr, alns = bam.load(alignments)
    vals = [abs(int(a.tlen)) for a in alns if a.mapped and a.tlen][:5000]
    vals = [v for v in vals if 0 < v <= int(max_size)]
    if not vals:
        return plot.empty_plot("no proper-pair TLEN values")
    return plot.histogram([float(v) for v in vals], bins=int(bins), xlabel="inner distance (bp)",
                         title=f"Insert size (n={len(vals)})")


@T("rna_seq_mapping_qc", "RNA-seq alignment QC metrics", RSEQC, "table",
   [sam("alignments", ALN_SAM, "Alignments (SAM)"), intin("head", 100, "Rows", min=1)],
   ex={"alignments": ALN_SAM}, up="RSeQC / hisat2 reporting", tags=("RNA-seq", "quality control"),
   summary="Mapped, duplicate, low-MAPQ and clipping statistics for an RNA-seq library.")
def rna_seq_mapping_qc(alignments, head=100):
    """Alignment QC summary."""
    hdr, alns = bam.load(alignments)
    mapped = [a for a in alns if a.mapped]
    names = Counter(a.qname for a in mapped)
    dup = sum(v - 1 for v in names.values() if v > 1)
    clip = sum(1 for a in mapped if "S" in a.cigar or "H" in a.cigar)
    rows = [{"metric": "total_records", "value": len(alns)},
            {"metric": "mapped", "value": len(mapped)},
            {"metric": "mapping_rate_percent", "value": round(100 * len(mapped) / max(1, len(alns)), 3)},
            {"metric": "duplicates", "value": dup},
            {"metric": "clipped", "value": clip},
            {"metric": "mapq_below_20", "value": sum(1 for a in mapped if a.mapq < 20)},
            {"metric": "reverse_strand_percent", "value": round(100 * sum(1 for a in mapped if a.reverse)
                                                             / max(1, len(mapped)), 3)},
            {"metric": "mean_mapq", "value": round(stats.mean([a.mapq for a in mapped]), 3) if mapped else 0}]
    return table(rows, f"{len(mapped)}/{len(alns)} mapped, {dup} duplicates")


# ===========================================================================
# single cell
# ===========================================================================
def _cell_matrix(src):
    """Cells x genes table -> (labels, matrix)."""
    df = tables.load(io.as_text(src))
    cols = list(df.columns)
    idc = cols[0] if cols[0].lower() in ("gene", "id", "cell", "barcode") else None
    if idc:
        df = df.set_index(idc)
    X = df.to_numpy(dtype=float)
    return list(df.index), list(df.columns), X


@T("seurat_normalise", "SCTransform-lite normalisation of a cell matrix", SC, "table",
   [tbl("src", COUNTS, "Cell x gene matrix"), intin("scale_factor", 10000, "Scale factor", min=100),
    intin("clip", 0, "Clip residuals at (0 = off)", min=0), intin("head", 200, "Rows", min=1)],
   ex={"src": COUNTS, "scale_factor": 5000}, up="Seurat NormalizeData / SCTransform",
   tags=("single cell", "normalisation"),
   summary="Normalise each cell to a total, log1p transform and return the expression matrix.")
def seurat_normalise(src, scale_factor=10000, clip=0, head=200):
    """Log-normalisation of a cell-by-gene matrix."""
    df = tables.load(io.as_text(src))
    cols = list(df.columns)
    idc = cols[0] if cols and cols[0].lower() in ("gene", "id", "cell", "barcode") else None
    mat = df.set_index(idc) if idc else df
    X = mat.to_numpy(dtype=float)
    tot = np.nansum(X, axis=1)
    norm = X / np.where(tot > 0, tot, 1.0)[:, None] * float(scale_factor)
    lg = np.log1p(norm)
    if clip:
        lg = np.clip(lg, -float(clip), float(clip))
    out = np.round(lg, 4)
    rows = [{"cell": str(mat.index[i]), "total_counts": int(tot[i]),
            "genes_detected": int(np.sum(X[i] > 0)),
            **{str(c): float(out[i, j]) for j, c in enumerate(mat.columns[:12])}}
           for i in range(min(len(mat.index), int(head)))]
    return table(rows, f"normalised {X.shape[0]} cells x {X.shape[1]} genes")


@T("seurat_qc_metrics", "Per-cell QC: counts, genes and mitochondrial share", SC, "table",
   [tbl("src", COUNTS, "Cell x gene matrix"), textbox("mito_prefix", "MT-", "Mitochondrial gene prefix"),
    number("min_genes", 1.0, "Minimum genes detected", min=0.0), intin("head", 200, "Rows", min=1)],
   ex={"src": COUNTS, "mito_prefix": "gene"}, up="Seurat dataQC / scDblFinder",
   tags=("single cell", "quality control"),
   summary="QC table per cell with total counts, detected genes and mito fraction plus a filter.")
def seurat_qc_metrics(src, mito_prefix="MT-", min_genes=1.0, head=200):
    """Cell QC metrics."""
    df = tables.load(io.as_text(src))
    cols = list(df.columns)
    idc = cols[0] if cols and cols[0].lower() in ("gene", "id", "cell", "barcode") else None
    mat = df.set_index(idc) if idc else df
    X = mat.to_numpy(dtype=float)
    is_mito = np.array([str(str(c)).upper().startswith(str(mito_prefix).upper()) for c in mat.columns]) \
        if mito_prefix else np.zeros(len(mat.columns), dtype=bool)
    rows = []
    for i, name in enumerate(mat.index):
        tot = float(np.nansum(X[i]))
        det = int(np.sum(X[i] > 0))
        mt = float(np.nansum(X[i][is_mito])) if is_mito.any() else 0.0
        rows.append({"cell": str(name), "n_genes": det, "n_counts": int(tot),
                    "mito_counts": int(mt), "mito_percent": round(100 * mt / tot, 3) if tot else 0.0,
                    "pass_qc": bool(det >= float(min_genes))})
    n_pass = sum(1 for r in rows if r["pass_qc"])
    return table(rows[:int(head)], f"{n_pass}/{len(rows)} cells pass QC")


@T("seurat_find_clusters", "Cluster cells on a k-NN graph", SC, "table",
   [tbl("src", COUNTS, "Cell x gene matrix"), intin("pcs", 2, "PCs to use", min=1),
    intin("resolution", 3, "Cluster count (k)", min=2), intin("neighbours", 5, "k-NN neighbours", min=2),
    intin("seed", 1, "Seed", min=0)],
   ex={"src": COUNTS, "pcs": 2, "resolution": 3}, up="Seurat FindNeighbors/FindClusters",
   tags=("single cell", "clustering", "PCA"),
   summary="PCA, mutual k-NN graph and label propagation into clusters (Louvain-lite).")
def seurat_find_clusters(src, pcs=2, resolution=3, neighbours=5, seed=1):
    """Cluster cells after PCA + kNN graph."""
    df = tables.load(io.as_text(src))
    cols = list(df.columns)
    idc = cols[0] if cols and cols[0].lower() in ("gene", "id", "cell", "barcode") else None
    mat = df.set_index(idc) if idc else df
    X = mat.to_numpy(dtype=float)
    if X.shape[0] < 3:
        return table([], "need at least three cells")
    res = ml.pca(X, components=min(int(pcs), min(X.shape)))
    Y = np.asarray(res["scores"], dtype=float)
    n = Y.shape[0]
    adj = defaultdict(set)
    for i in range(n):
        d = np.linalg.norm(Y - Y[i], axis=1)
        d[i] = 1e18
        for j in np.argsort(d)[: int(neighbours)]:
            adj[i].add(int(j))
            adj[int(j)].add(i)
    labels = list(range(n))
    for _ in range(12):  # label propagation
        for i in range(n):
            nb = [labels[j] for j in adj[i]]
            if nb:
                labels[i] = Counter(nb).most_common(1)[0][0]
    # merge to the requested number of clusters
    uniq = sorted(set(labels))
    if len(uniq) > int(resolution):
        groups = defaultdict(list)
        for i, l in enumerate(labels):
            groups[l].append(i)
        cents = {l: np.mean(Y[v], axis=0) for l, v in groups.items()}
        kk = ml.kmeans(Y, k=max(2, int(resolution)), seed=int(seed))
        labels = [int(x) for x in np.asarray(kk["labels"])]
    rows = [{"cell": str(mat.index[i]), "cluster": int(labels[i]),
            "PC1": round(float(Y[i, 0]), 4), "PC2": round(float(Y[i, 1]), 4) if Y.shape[1] > 1 else 0.0,
            "knn_degree": len(adj[i])} for i in range(n)]
    sizes = Counter(labels)
    res_tbl = table(rows, f"{len(sizes)} clusters over {n} cells "
                          f"({', '.join(f'{k}:{v}' for k, v in sorted(sizes.items()))})")
    res_tbl["stats"] = {"explained_variance_ratio": [round(float(v), 4) for v in
                                                np.asarray(res.get("explained_variance_ratio", []))]}
    return res_tbl


@T("seurat_marker_genes", "Marker genes per cluster", SC, "table",
   [tbl("src", COUNTS, "Cell x gene matrix"), tbl("clusters", "", "Cell / cluster table"),
    number("min_auc", 0.0, "Minimum AUC", min=0.0), intin("head", 300, "Rows", min=1),
    number("min_logfc", 0.0, "Minimum |log2FC|", min=0.0)],
   ex={"src": COUNTS, "clusters": COUNTS}, up="Seurat FindAllMarkers / presto",
   tags=("single cell", "markers", "statistics"),
   summary="Rank-sum test of each gene inside versus outside a cluster with fold changes.")
def seurat_marker_genes(src, clusters, min_auc=0.0, head=300, min_logfc=0.0):
    """Marker detection by Wilcoxon rank-sum (presto-lite)."""
    df = tables.load(io.as_text(src))
    cols = list(df.columns)
    idc = cols[0] if cols and cols[0].lower() in ("gene", "id", "cell", "barcode") else None
    mat = df.set_index(idc) if idc else df
    cl = tables.load(io.as_text(clusters)) if clusters else None
    assign: dict[str, str] = {}
    if cl is not None and not cl.empty:
        cc = list(cl.columns)
        for r in cl.to_dict("records"):
            assign[str(r[cc[0]])] = str(r[cc[1]]) if len(cc) > 1 else ""
    groups: dict[str, list[int]] = defaultdict(list)
    for i, name in enumerate(mat.index):
        groups[str(assign.get(str(name), "all"))].append(i)
    X = np.log1p(np.nan_to_num(mat.to_numpy(dtype=float)))
    rows = []
    for gid, idxs in groups.items():
        if len(idxs) == 0 or len(idxs) == len(X):
            continue
        other = [i for i in range(len(X)) if i not in set(idxs)]
        for j, gene in enumerate(mat.columns):
            a, b = X[idxs, j], X[other, j]
            st = stats.mann_whitney_u(list(a), list(b))
            auc = float(st["U1"]) / max(1.0, len(a) * len(b))
            lfc = float(np.mean(a) - np.mean(b))
            if float(st["p_value"]) > 0.5 or abs(lfc) < float(min_logfc) or auc < float(min_auc):
                continue
            rows.append({"cluster": gid, "gene": str(gene), "avg_log2FC": round(lfc, 4),
                        "p_value": float(st["p_value"]), "pct_in": round(100 * float(np.mean(a > 0)), 2),
                        "pct_out": round(100 * float(np.mean(b > 0)), 2), "AUC": round(auc, 4)})
    rows.sort(key=lambda r: r["p_value"])
    return table(rows[:int(head)], f"{len(rows)} candidate markers")


@T("seurat_umap", "UMAP / t-SNE embedding of cells", SC, "figure",
   [tbl("src", COUNTS, "Cell x gene matrix"), choice("method", ["umap", "tsne", "pca"], "umap", "Method"),
    intin("neighbours", 8, "Neighbour count / perplexity", min=2), intin("components", 2, "Dimensions", min=2)],
   ex={"src": COUNTS, "method": "umap"}, up="Seurat RunUMAP / scanpy",
   tags=("single cell", "embedding", "plot"),
   summary="Two-dimensional embedding of cells coloured by cluster index.")
def seurat_umap(src, method="umap", neighbours=8, components=2):
    """Cell embedding plot."""
    df = tables.load(io.as_text(src))
    cols = list(df.columns)
    idc = cols[0] if cols and cols[0].lower() in ("gene", "id", "cell", "barcode") else None
    mat = df.set_index(idc) if idc else df
    X = mat.to_numpy(dtype=float)
    X = np.log1p(X)
    if X.shape[0] < 3:
        return plot.empty_plot("need at least three cells")
    if method == "tsne":
        res = ml.tsne(X, components=int(components), perplexity=max(1.0, min(float(neighbours), X.shape[0] - 1)))
    elif method == "pca":
        res = ml.pca(X, components=int(components))
    else:
        res = ml.umap_like(X, components=int(components), neighbours=int(neighbours))
    E = np.asarray(res["embedding"] if "embedding" in res else res["scores"], dtype=float)
    labels = [str(i) for i in mat.index]
    grp = [i % 3 for i in range(len(labels))]
    return plot.scatter(list(E[:, 0]), list(E[:, 1]), labels=labels, colour_by=grp,
                       xlabel="dim 1", ylabel="dim 2", title=f"{method} embedding")


@T("scanpy_score_genes", "Score a gene set per cell", SC, "table",
   [tbl("src", COUNTS, "Cell x gene matrix"), txt("sets", GMT, "Gene sets"), intin("ctrl", 50,
                                                           "Background genes", min=1),
    boolean("scale", True, "z-score before averaging"), intin("head", 200, "Rows", min=1)],
   ex={"src": COUNTS, "sets": GMT}, up="scanpy score_genes", tags=("single cell", "scoring"),
   summary="Mean expression of a signature relative to random control genes per cell.")
def scanpy_score_genes(src, sets, ctrl=50, scale=True, head=200):
    """Gene-set scores per cell."""
    df = tables.load(io.as_text(src))
    cols = list(df.columns)
    idc = cols[0] if cols and cols[0].lower() in ("gene", "id", "cell", "barcode") else None
    mat = df.set_index(idc) if idc else df
    X = np.log1p(np.nan_to_num(mat.to_numpy(dtype=float)))
    gsets = io.parse_gmt(io.as_text(sets))
    if scale and X.shape[1] > 1:
        Xs = (X - X.mean(axis=1, keepdims=True)) / np.where(X.std(axis=1, keepdims=True) > 0,
                                                            X.std(axis=1, keepdims=True), 1)
    else:
        Xs = X
    gene_pos = {str(g): j for j, g in enumerate(mat.columns)}
    rows = []
    for name, members in gsets.items():
        idx = [gene_pos[m] for m in dict.fromkeys(members) if m in gene_pos]
        if not idx:
            continue
        bg = [j for j in range(Xs.shape[1]) if j not in set(idx)][: int(ctrl)]
        for i, cell in enumerate(mat.index):
            score = float(np.mean(Xs[i, idx]) - (np.mean(Xs[i, bg]) if bg else 0.0))
            rows.append({"cell": str(cell), "gene_set": name, "score": round(score, 5),
                        "n_genes": len(idx)})
    return table(rows[:int(head)], f"{len(rows)} cell scores for {len(gsets)} sets")


@T("monocle_trajectory", "Order cells along a pseudotime trajectory", SC, "table",
   [tbl("src", COUNTS, "Cell x gene matrix"), intin("components", 2, "PCs", min=1),
    textbox("root_cell", "", "Root cell id (empty = most extreme cell)"), intin("head", 200, "Rows", min=1)],
   ex={"src": COUNTS, "components": 2}, up="Monocle3 order_cells", tags=("single cell", "trajectory"),
   summary="Pseudotime from geodesic distance along the principal curve of the embedding.")
def monocle_trajectory(src, components=2, root_cell="", head=200):
    """Pseudotime ordering."""
    df = tables.load(io.as_text(src))
    cols = list(df.columns)
    idc = cols[0] if cols and cols[0].lower() in ("gene", "id", "cell", "barcode") else None
    mat = df.set_index(idc) if idc else df
    X = np.log1p(np.nan_to_num(mat.to_numpy(dtype=float)))
    res = ml.pca(X, components=min(int(components), min(X.shape)))
    Y = np.asarray(res["scores"], dtype=float)
    if Y.shape[0] < 2:
        return table([], "need at least two cells")
    start = 0
    if root_cell:
        names = [str(i) for i in mat.index]
        start = names.index(root_cell) if root_cell in names else 0
    else:
        start = int(np.argmin(Y[:, 0]))
    order = [start]
    remaining = set(range(len(Y))) - {start}
    cur = start
    while remaining:
        d = np.linalg.norm(Y[list(remaining)] - Y[cur], axis=1)
        nxt = list(remaining)[int(np.argmin(d))]
        order.append(nxt)
        remaining.discard(nxt)
        cur = nxt
    rows = []
    for rank, i in enumerate(order, start=1):
        rows.append({"cell": str(mat.index[i]), "pseudotime": rank,
                    "PC1": round(float(Y[i, 0]), 4),
                    "PC2": round(float(Y[i, 1]), 4) if Y.shape[1] > 1 else 0.0,
                    "distance_from_root": round(float(np.sum(np.linalg.norm(np.diff(Y[order[: rank + 1]], axis=0), axis=1))), 5)
                    if rank > 1 else 0.0})
    return table(rows[:int(head)], f"{len(rows)} cells ordered from {rows[0]['cell']}")


@T("presto_wilcoxon_markers", "Presto-style Wilcoxon marker ranking", SC, "table",
   [tbl("src", COUNTS, "Matrix"), textbox("group_column", "", "Column holding group labels"),
    intin("head", 300, "Rows", min=1)],
   ex={"src": PHENO, "group_column": "label"}, up="presto wilcoxauc",
   tags=("single cell", "rank test", "markers"),
   summary="Fast rank-sum ranking of features between two groups defined by a label column.")
def presto_wilcoxon_markers(src, group_column="", head=300):
    """Wilcoxon marker ranking on a wide table."""
    df = tables.load(io.as_text(src))
    cols = list(df.columns)
    gcol = group_column or next((c for c in cols if "group" in c.lower() or "label" in c.lower()
                                or "type" in c.lower() or "condition" in c.lower()), "")
    if not gcol or gcol not in cols:
        return table([], "no group label column found - set group_column")
    labs = [str(v) for v in df[gcol].tolist()]
    uniq = sorted(set(labs))
    if len(uniq) < 2:
        return table([], "the label column needs at least two groups")
    g1, g2 = uniq[0], uniq[1]
    rows = []
    for c in cols:
        if c == gcol:
            continue
        try:
            v = [float(x) for x in df[c].tolist()]
        except (TypeError, ValueError):
            continue
        a = [v[i] for i in range(len(v)) if labs[i] == g1]
        b = [v[i] for i in range(len(v)) if labs[i] == g2]
        if len(a) < 2 or len(b) < 2:
            continue
        st = stats.mann_whitney_u(a, b)
        rows.append({"feature": c, "group_1": g1, "group_2": g2, "U": float(st["U"]),
                    "z": round(float(st["z"]), 4), "p_value": float(st["p_value"]),
                    "median_1": round(float(np.median(a)), 4), "median_2": round(float(np.median(b)), 4),
                    "log2fc": round(math.log2((abs(np.mean(b)) + 1e-6) / (abs(np.mean(a)) + 1e-6)), 4)})
    rows.sort(key=lambda r: r["p_value"])
    return table(rows[:int(head)], f"{len(rows)} features tested ({g1} vs {g2})")


# ===========================================================================
# spatial, Hi-C and multi-omics
# ===========================================================================
@T("spatial_spot_clusters", "Cluster spatial spots into domains", SPATIAL, "table",
   [tbl("coordinates", COUNTS, "Spot table (x, y, optionally expression)"), intin("k", 3, "Domains", min=2),
    intin("neighbours", 4, "Neighbours to smooth over", min=0), intin("head", 500, "Rows", min=1)],
   ex={"coordinates": COUNTS, "k": 3, "neighbours": 2}, up="Scanpy spatial / Seurat SPOTlight",
   tags=("spatial", "clustering", "domains"),
   summary="Smooth spot expression over its spatial neighbours then cluster into domains.")
def spatial_spot_clusters(coordinates, k=3, neighbours=4, head=500):
    """Spatial domain calling."""
    df = tables.load(io.as_text(coordinates))
    cols = list(df.columns)
    xs = np.arange(len(df), dtype=float)
    ys = np.zeros(len(df))
    for j, c in enumerate(cols[:2]):
        try:
            vals = np.array([float(v) for v in df[c].tolist()], dtype=float)
            if np.isfinite(vals).all():
                if j == 0:
                    xs = vals
                else:
                    ys = vals
        except (TypeError, ValueError):
            pass
    expr_cols = [c for c in cols if c not in cols[:2]]
    M = df[expr_cols].to_numpy(dtype=float) if expr_cols else np.column_stack([xs, ys])
    M = np.nan_to_num(M)
    if neighbours and M.shape[0] > 2:
        Sm = np.zeros_like(M)
        for i in range(M.shape[0]):
            d = np.linalg.norm(np.column_stack([xs, ys]) - np.column_stack([xs[i], ys[i]]), axis=1)
            d[i] = 1e18
            nb = list(np.argsort(d)[: max(1, int(neighbours))]) + [i]
            Sm[i] = np.mean(M[nb], axis=0)
        M = Sm
    kk = ml.kmeans(M, k=max(2, int(k)), seed=1)
    labels = np.asarray(kk["labels"], dtype=int)
    rows = [{"spot": str(r) if r else i, "x": round(float(xs[i]), 3), "y": round(float(ys[i]), 3),
            "domain": int(labels[i])} for i, r in enumerate(df.index if not expr_cols else df[df.columns[0]].tolist())]
    for i, r in enumerate(rows):
        for j, c in enumerate(expr_cols[:4]):
            r[c] = round(float(M[i, j]), 4) if j < M.shape[1] else 0.0
    res = table(rows[:int(head)], f"{len(set(labels))} spatial domains over {len(df)} spots")
    res["stats"] = {"domain_sizes": dict(Counter(int(x) for x in labels))}
    return res


@T("spatial_expression_overlay", "Spot plot coloured by a marker", SPATIAL, "figure",
   [tbl("coordinates", COUNTS, "Spot table"), textbox("marker", "", "Marker / expression column"),
    textbox("x_column", "", "x column"), textbox("y_column", "", "y column")],
   ex={"coordinates": COUNTS, "marker": "sample_D"}, up="scanpy plot_spatial",
   tags=("spatial", "plot", "expression"),
   summary="Scatter the tissue coordinates with colour encoding marker expression.")
def spatial_expression_overlay(coordinates, marker="", x_column="", y_column=""):
    """Spatial expression scatter."""
    df = tables.load(io.as_text(coordinates))
    cols = list(df.columns)
    xs_col = x_column or cols[0]
    ys_col = y_column or (cols[1] if len(cols) > 1 else cols[0])
    def numeric(c):
        try:
            return np.array([float(v) for v in df[c].tolist()], dtype=float)
        except (KeyError, TypeError, ValueError):
            return np.arange(len(df), dtype=float)
    xs, ys = numeric(xs_col), numeric(ys_col)
    mk = marker or next((c for c in cols if c not in (xs_col, ys_col)), None)
    vals = numeric(mk) if mk else None
    return plot.scatter(list(xs), list(ys), colour_by=list(vals) if vals is not None else None,
                       xlabel=str(xs_col), ylabel=str(ys_col),
                       title=f"Spatial overlay ({mk or 'no marker'})")


@T("hic_contact_matrix", "Contact matrix from a pairs file", HIC, "table",
   [txt("pairs", PAIRS, "Hi-C pairs file"), intin("bins", 1000, "Bin size", min=100),
    textbox("chrom", "", "Chromosome filter"), boolean("normalise", True, "Vanilla coverage normalisation"),
    intin("head", 50, "Rows", min=1)],
   ex={"pairs": PAIRS, "bins": 1000}, up="cooler dump / HiCExplorer", tags=("Hi-C", "contacts", "matrix"),
   summary="Bin ligation events into a matrix and optionally apply vanilla-ICE normalisation.")
def hic_contact_matrix(pairs, bins=1000, chrom="", normalise=True, head=50):
    """Binned contact matrix."""
    size = int(bins)
    recs = io.parse_pairs(io.as_text(pairs)) if hasattr(io, "parse_pairs") else []
    mat: dict[tuple, int] = Counter()
    cov: Counter = Counter()
    chrs = set()
    for r in recs:
        c1 = str(r.get("chrom1", r.get("chrom", "")))
        c2 = str(r.get("chrom2", r.get("chrom", "")))
        if chrom and c1 != chrom:
            continue
        p1, p2 = int(float(r.get("pos1", r.get("start", 0)) or 0)), int(float(r.get("pos2", r.get("end", 0)) or 0))
        b1, b2 = p1 // size, p2 // size
        key = tuple(sorted([f"{c1}:{b1}", f"{c2}:{b2}"]))
        mat[key] += 1
        cov[f"{c1}:{b1}"] += 1
        cov[f"{c2}:{b2}"] += 1
        chrs.update([c1, c2])
    bins_list = sorted({k for pair in mat for k in pair})
    if not bins_list:
        return table([], "no ligation events parsed")
    grid = {b: i for i, b in enumerate(bins_list)}
    M = np.zeros((len(bins_list), len(bins_list)))
    for (a, b), n in mat.items():
        M[grid[a], grid[b]] = n
        if a != b:
            M[grid[b], grid[a]] = n
    if normalise:
        v = np.array([max(1, cov[b]) for b in bins_list], dtype=float)
        outer = np.sqrt(np.outer(v, v))
        M = M / np.where(outer > 0, outer, 1) * 1000.0
    rows = []
    for i, a in enumerate(bins_list[: int(head)]):
        rows.append({"bin": a, **{b: round(float(M[i, grid[b]]), 4) for b in bins_list[: min(12, len(bins_list))]}})
    res = table(rows, f"{len(bins_list)} x {len(bins_list)} contact matrix from {int(sum(mat.values()))} pairs")
    res["stats"]["figure"] = plot.heatmap([[round(float(v), 3) for v in row] for row in M[: min(40, len(M))]],
                                        row_labels=bins_list[: min(40, len(M))],
                                        col_labels=bins_list[: min(40, len(M))], title="Hi-C contacts")
    return res


@T("hic_insulation_profile", "Insulation score along the genome", HIC, "figure",
   [txt("pairs", PAIRS, "Hi-C pairs"), intin("bins", 1000, "Bin size", min=100),
    intin("window", 3, "Insulation windows", min=1)],
   ex={"pairs": PAIRS, "bins": 1000, "window": 2}, up="insulation4phipson / HiCExplorer insulation",
   tags=("Hi-C", "TAD", "insulation"),
   summary="Boundary detection from the fraction of contacts crossing each bin.")
def hic_insulation_profile(pairs, bins=1000, window=3):
    """Insulation profile."""
    size = int(bins)
    recs = io.parse_pairs(io.as_text(pairs))
    mat: Counter = Counter()
    for r in recs:
        c1 = str(r.get("chrom1", r.get("chrom", "")))
        p1 = int(float(r.get("pos1", r.get("start", 0)) or 0)) // size
        p2 = int(float(r.get("pos2", r.get("end", 0)) or 0)) // size
        c2 = str(r.get("chrom2", r.get("chrom", "")))
        mat[(c1, p1, c2, p2)] += 1
    by_chrom: dict[str, list] = defaultdict(set)
    for c1, b1, c2, b2 in mat:
        by_chrom[c1].add(b1)
        by_chrom[c2].add(b2)
    series = {}
    for c, bs in by_chrom.items():
        bl = sorted(bs)
        idx = {b: i for i, b in enumerate(bl)}
        n = len(bl)
        M = np.zeros((n, n))
        for (c1, b1, c2, b2), v in mat.items():
            if c1 in idx and c2 in idx:
                M[idx[b1], idx[b2]] += v
        w = max(1, int(window))
        iso = []
        for i in range(n):
            lo, hi = max(0, i - w), min(n, i + w + 1)
            inside = float(np.sum(M[lo:hi, lo:hi]))
            crossing = float(sum(M[i, j] for j in range(n) if j < lo or j >= hi))
            iso.append(round(crossing / max(1e-9, inside + crossing), 5))
        series[c] = iso
    if not series:
        return plot.empty_plot("no contacts parsed")
    return plot.line(series, x=list(range(1, max(len(v) for v in series.values()) + 1)),
                    xlabel="bin index", ylabel="insulation", title="Insulation profile")


@T("hic_compartments", "A/B compartment eigenvector", HIC, "table",
   [txt("pairs", PAIRS, "Hi-C pairs"), intin("bins", 1000, "Bin size", min=100),
    choice("chrom", ["all"], "all", "Chromosome"), intin("head", 100, "Rows", min=1)],
   ex={"pairs": PAIRS, "bins": 1000}, up="eigenvector_c / HiCExplorer",
   tags=("Hi-C", "compartments"),
   summary="First eigenvector of the correlation matrix - the classic A/B compartment signal.")
def hic_compartments(pairs, bins=1000, chrom="all", head=100):
    """Compartment eigenvector."""
    size = int(bins)
    recs = io.parse_pairs(io.as_text(pairs))
    mat: Counter = Counter()
    for r in recs:
        c1 = str(r.get("chrom1", r.get("chrom", "")))
        c2 = str(r.get("chrom2", r.get("chrom", "")))
        if chrom != "all" and c1 != chrom:
            continue
        b1 = int(float(r.get("pos1", r.get("start", 0)) or 0)) // size
        b2 = int(float(r.get("pos2", r.get("end", 0)) or 0)) // size
        key = tuple(sorted([b1, b2]))
        mat[(c1, key)] += 1
    if not mat:
        return table([], "no contacts")
    c = next(iter(mat))[0]
    bl = sorted({k[1][0] for k in mat if k[0] == c} | {k[1][1] for k in mat if k[0] == c})
    idx = {b: i for i, b in enumerate(bl)}
    M = np.zeros((len(bl), len(bl)))
    for (cc, (b1, b2)), v in mat.items():
        if cc == c:
            M[idx[b1], idx[b2]] = v
            M[idx[b2], idx[b1]] = v
    C = np.corrcoef(M + np.eye(len(M)) * 1e-6) if len(M) > 1 else M
    vals, vecs = np.linalg.eigh(C)
    ev = vecs[:, -1] if len(vals) else np.zeros(len(bl))
    rows = [{"chromosome": c, "bin": b, "start": b * size, "end": (b + 1) * size,
            "eigenvector": round(float(ev[idx[b]]), 5),
            "compartment": "A" if ev[idx[b]] > 0 else "B",
            "contacts": int(sum(v for (cc, (x, y)), v in mat.items() if cc == c and b in (x, y)))}
           for b in bl]
    return table(rows[:int(head)], f"{len(rows)} bins, {sum(1 for r in rows if r['compartment'] == 'A')} A-like")


@T("multiomics_merge_tables", "Join two omics tables by feature", MULTI, "table",
   [tbl("left", COUNTS, "Left table"), tbl("right", PHENO, "Right table"),
    textbox("key_left", "", "Left key column"), textbox("key_right", "", "Right key column"),
    choice("how", ["inner", "left", "right", "outer"], "inner", "Join type"),
    number("min_abs_corr", 0.0, "Report |correlation| above", min=0.0)],
   ex={"left": COUNTS, "right": COUNTS, "key_left": "gene", "key_right": "gene"},
   up="multi-omics integration (MOFA-lite)", tags=("multiomics", "join", "integration"),
   summary="Merge expression with another data layer and correlate the shared columns.")
def multiomics_merge_tables(left, right, key_left="", key_right="", how="inner", min_abs_corr=0.0):
    """Merge two tables and cross-correlate."""
    import pandas as pd

    dl, dr = tables.load(io.as_text(left)), tables.load(io.as_text(right))
    cl, cr = list(dl.columns), list(dr.columns)
    kl = key_left or cl[0]
    kr = key_right or cr[0]
    merged = dl.merge(dr, how=how, left_on=kl, right_on=kr, suffixes=("", "_r"))
    num_l = [c for c in cl if c != kl]
    num_r = [c for c in cr if c != kr]
    rows = []
    for a in num_l:
        for b in num_r:
            if a not in merged.columns or b not in merged.columns:
                continue
            x = pd.to_numeric(merged[a], errors="coerce").to_numpy(dtype=float)
            y = pd.to_numeric(merged[b], errors="coerce").to_numpy(dtype=float)
            m = np.isfinite(x) & np.isfinite(y)
            if m.sum() < 3:
                continue
            st = stats.pearson(list(x[m]), list(y[m]))
            if abs(float(st["r"])) < float(min_abs_corr):
                continue
            rows.append({"left_column": a, "right_column": b, "pearson_r": round(float(st["r"]), 4),
                        "p_value": float(st.get("p_value", 1.0)), "n": int(m.sum())})
    rows.sort(key=lambda r: -abs(r["pearson_r"]))
    return table(rows[:400], f"{len(merged)} merged rows; {len(rows)} cross-omics correlations")


@T("multiomics_rank_product", "Rank product consensus across datasets", MULTI, "table",
   [tbl("left", COUNTS, "Results table 1"), tbl("right", COUNTS, "Results table 2"),
    textbox("score_left", "", "Ranking column in table 1"), textbox("score_right", "", "Ranking column in table 2"),
    textbox("key_left", "", "Feature column (left)"), textbox("key_right", "", "Feature column (right)"),
    intin("head", 200, "Rows", min=1)],
   ex={"left": COUNTS, "right": COUNTS, "score_left": "sample_A", "score_right": "sample_D",
       "key_left": "gene", "key_right": "gene"},
   up="RankProd / tie within replicates", tags=("multiomics", "meta-analysis", "ranks"),
   summary="Combine two ranked lists with the rank product to find consistently top features.")
def multiomics_rank_product(left, right, score_left="", score_right="", key_left="", key_right="", head=200):
    """Rank product across two result tables."""
    def rank_map(src, key, score):
        df = tables.load(io.as_text(src))
        cols = list(df.columns)
        k = key or cols[0]
        s = score or cols[-1]
        recs = []
        for r in df.to_dict("records"):
            try:
                recs.append((str(r[k]), abs(float(r[s]))))
            except (KeyError, TypeError, ValueError):
                continue
        recs.sort(key=lambda x: -x[1])
        n = max(1, len(recs))
        return {name: i + 1 for i, (name, _v) in enumerate(recs)}, n

    ra, na = rank_map(left, key_left, score_left)
    rb, nb = rank_map(right, key_right, score_right)
    shared = [f for f in ra if f in rb]
    rows = [{"feature": f, "rank_1": ra[f], "rank_2": rb[f],
            "rank_product": ra[f] * rb[f],
            "geometric_mean_rank": round(math.sqrt(ra[f] * rb[f]), 3)} for f in shared]
    rows.sort(key=lambda r: r["rank_product"])
    return table(rows[:int(head)], f"rank product over {len(shared)} shared features "
                                  f"({na} x {nb} inputs)")
