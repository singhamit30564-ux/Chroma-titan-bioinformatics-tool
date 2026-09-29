"""Statistics, machine learning, network/graph displays and interactive views.

Panel sections covered: **statistics**, **machine_learning**,
**graph_display_data**, **interactive_tools**.
"""
from __future__ import annotations

import math
from collections import Counter, defaultdict

import numpy as np
import pandas as pd

from chroma_titan.core import io, ml, plot, stats, tables
from chroma_titan.tools._common import *  # noqa: F401,F403
from chroma_titan.tools._common import COUNTS, GMT, PHENO, REGIONS_BED, BEDGRAPH, PAIRS

STAT = "statistics"
MLSEC = "machine_learning"
GRAPH = "graph_display_data"
INTERACTIVE = "interactive_tools"


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------
def _frame(src):
    """Table -> DataFrame (first column kept as a normal column)."""
    df = tables.load(io.as_text(src))
    return df


def _num_cols(df):
    out = []
    for c in df.columns:
        v = pd.to_numeric(df[c], errors="coerce").to_numpy(dtype=float)
        if np.isfinite(v).sum() >= 2:
            out.append((str(c), v))
    return out


def _mat(src, index_col=True):
    """Wide numeric matrix -> (index labels, column names, ndarray)."""
    df = _frame(src)
    if df.empty:
        return [], [], np.zeros((0, 0))
    cols = list(df.columns)
    idc = cols[0] if index_col else None
    body = df.drop(columns=[idc]) if idc and idc in df.columns else df
    rows = [str(x) for x in (df[idc].tolist() if idc else range(len(body)))]
    data = np.column_stack([pd.to_numeric(body[c], errors="coerce").to_numpy(dtype=float)
                           for c in body.columns]) if len(body.columns) else np.zeros((len(body), 0))
    return rows, [str(c) for c in body.columns], data


def _groups_from_labels(labels: list[str]) -> dict[str, list[int]]:
    out: dict[str, list[int]] = defaultdict(list)
    for i, v in enumerate(labels):
        out[str(v)].append(i)
    return dict(sorted(out.items()))


def _label_columns(df):
    return [c for c in df.columns
            if not np.isfinite(pd.to_numeric(df[c], errors="coerce").to_numpy(dtype=float)).all()]


# ===========================================================================
# descriptive statistics
# ===========================================================================
@T("stats_describe_columns", "Descriptive statistics of every numeric column", STAT, "table",
   [tbl("src", COUNTS, "Data table"), boolean("include_index", False, "Treat first column as data")],
   ex={"src": COUNTS}, up="datamash describe / summary", tags=("descriptive statistics",),
   summary="n, mean, median, spread, quantiles and skewness for each numeric column.")
def stats_describe_columns(src, include_index=False):
    """Per-column descriptive statistics."""
    df = _frame(src)
    rows = []
    for name, v in _num_cols(df):
        f = v[np.isfinite(v)]
        rows.append({"column": name, "n": int(f.size), "mean": round(float(np.mean(f)), 5),
                    "median": round(float(np.median(f)), 5), "stdev": round(float(stats.stdev(f)), 5),
                    "min": round(float(np.min(f)), 5), "max": round(float(np.max(f)), 5),
                    "q1": round(float(stats.quantile(f, 0.25)), 5), "q3": round(float(stats.quantile(f, 0.75)), 5),
                    "iqr": round(float(stats.iqr(f)), 5) if f.size > 2 else 0.0,
                    "cv": round(float(stats.cv(f)), 5) if np.mean(f) else 0.0,
                    "skewness": round(float(stats.skewness(f)), 4) if f.size > 2 else 0.0,
                    "kurtosis": round(float(stats.kurtosis(f)), 4) if f.size > 3 else 0.0})
    return table(rows, f"{len(rows)} numeric columns described")


@T("stats_group_summary", "Summarise a value column by group", STAT, "table",
   [tbl("src", PHENO, "Data table"), textbox("value_column", "", "Value column"),
    textbox("group_column", "", "Group column"), choice("agg", ["mean", "median", "sum", "count", "stdev",
                                                                "min", "max"], "mean", "Statistic"),
    boolean("sem", True, "Also report SEM")],
   ex={"src": PHENO, "value_column": "yield", "group_column": "group"},
   up="datamash groupby / aggregate", tags=("aggregation", "groups"),
   summary="Group-wise aggregation with sample size, spread and standard error.")
def stats_group_summary(src, value_column="", group_column="", agg="mean", sem=True):
    """Grouped summary table."""
    df = _frame(src)
    cols = list(df.columns)
    vcol = value_column or next((c for c in cols if c != (group_column or cols[0])), cols[-1])
    gcol = group_column or cols[0]
    groups: dict[str, list[float]] = defaultdict(list)
    for r in df.to_dict("records"):
        try:
            groups[str(r[gcol])].append(float(r[vcol]))
        except (KeyError, TypeError, ValueError):
            continue
    fn = {"mean": np.mean, "median": np.median, "sum": np.sum, "count": lambda x: float(len(x)),
          "stdev": (lambda x: float(np.std(x, ddof=1)) if len(x) > 1 else 0.0),
          "min": np.min, "max": np.max}[agg]
    rows = []
    for g, vals in sorted(groups.items()):
        v = np.array([x for x in vals if math.isfinite(x)], dtype=float)
        if not v.size:
            continue
        rec = {"group": g, "n": int(v.size), f"{agg}_{vcol}": round(float(fn(v)), 5),
              "mean": round(float(np.mean(v)), 5), "stdev": round(float(np.std(v, ddof=1)), 5)
              if v.size > 1 else 0.0}
        if sem:
            rec["sem"] = round(float(stats.sem(v)), 5) if v.size > 1 else 0.0
        rows.append(rec)
    return table(rows, f"{len(rows)} groups over column {vcol}")


@T("stats_ttest", "Two-group t-test (and Welch)", STAT, "table",
   [tbl("src", COUNTS, "Matrix (rows = features)"), textbox("group_a", "", "Group A samples"),
    textbox("group_b", "", "Group B samples"), choice("variant", ["student", "welch", "paired"], "student",
                                                      "Test variant")],
   ex={"src": COUNTS, "group_a": "sample_A,sample_B,sample_C", "group_b": "sample_D,sample_E,sample_F"},
   up="t.test / statsmodels", tags=("t-test", "statistics"),
   summary="Student, Welch or paired t-test per feature with the mean difference and effect size.")
def stats_ttest(src, group_a="", group_b="", variant="student"):
    """t-tests across a matrix."""
    rows_, cols_, X = _mat(src)
    ia = [cols_.index(c) for c in str(group_a).split(",") if c in cols_] or list(range(len(cols_) // 2))
    ib = [cols_.index(c) for c in str(group_b).split(",") if c in cols_] or list(range(len(cols_) // 2, len(cols_)))
    out = []
    ps = []
    for i, name in enumerate(rows_):
        a = X[i, ia][np.isfinite(X[i, ia])]
        b = X[i, ib][np.isfinite(X[i, ib])]
        if a.size < 2 or b.size < 2:
            continue
        if variant == "paired" and a.size == b.size:
            st = stats.paired_ttest(list(a), list(b))
        elif variant == "welch":
            st = stats.ttest_ind(list(b), list(a), equal_var=False)
        else:
            st = stats.ttest_ind(list(b), list(a))
        p = float(st["p_value"])
        ps.append(p)
        out.append({"feature": name, "mean_a": round(float(np.mean(a)), 4), "mean_b": round(float(np.mean(b)), 4),
                   "difference": round(float(np.mean(b) - np.mean(a)), 4),
                   "t": round(float(st["statistic"]), 4), "df": st.get("df", ""), "p_value": p,
                   "cohens_d": round(float(st.get("cohens_d", 0.0) or 0.0), 4)})
    adj = stats.p_adjust(ps, "fdr_bh")
    for k, r in enumerate(out):
        r["p_adjusted"] = round(float(adj[k]), 6) if k < len(adj) else 1.0
    return table(out, f"{sum(1 for r in out if r['p_adjusted'] < 0.05)}/{len(out)} features at FDR < 0.05")


@T("stats_rank_tests", "Mann-Whitney and Wilcoxon tests", STAT, "table",
   [tbl("src", COUNTS, "Matrix"), textbox("group_a", "", "Group A columns"),
    textbox("group_b", "", "Group B columns"), choice("test", ["mann_whitney", "wilcoxon"], "mann_whitney",
                                                       "Test")],
   ex={"src": COUNTS, "group_a": "sample_A,sample_B,sample_C", "group_b": "sample_D,sample_E,sample_F"},
   up="wilcox.test / exactRankTests", tags=("nonparametric", "statistics"),
   summary="Distribution-free group comparison for small or non-normal samples.")
def stats_rank_tests(src, group_a="", group_b="", test="mann_whitney"):
    """Non-parametric tests."""
    rows_, cols_, X = _mat(src)
    ia = [cols_.index(c) for c in str(group_a).split(",") if c in cols_] or list(range(len(cols_) // 2))
    ib = [cols_.index(c) for c in str(group_b).split(",") if c in cols_] or list(range(len(cols_) // 2, len(cols_)))
    out, ps = [], []
    for i, name in enumerate(rows_):
        a, b = X[i, ia], X[i, ib]
        if test == "wilcoxon":
            if a.size != b.size or a.size < 3:
                continue
            st = stats.wilcoxon_signed_rank(list(b - a))
            p, stat = float(st["p_value"]), float(st.get("V", st.get("statistic", 0)))
            extra = {}
        else:
            if a.size < 2 or b.size < 2:
                continue
            st = stats.mann_whitney_u(list(a), list(b))
            p, stat = float(st["p_value"]), float(st["U"])
            extra = {"median_a": round(float(np.median(a)), 4), "median_b": round(float(np.median(b)), 4)}
        ps.append(p)
        out.append({"feature": name, "statistic": round(stat, 4), "p_value": p, **extra})
    adj = stats.p_adjust(ps, "fdr_bh")
    for k, r in enumerate(out):
        r["p_adjusted"] = round(float(adj[k]), 6) if k < len(adj) else 1.0
    return table(out, f"{len(out)} features tested with {test}")


@T("stats_correlation_matrix", "Correlation matrix of columns", STAT, "figure",
   [tbl("src", COUNTS, "Numeric table"), choice("method", ["pearson", "spearman", "kendall"], "pearson",
                                                 "Method"), boolean("annot", True, "Show values"),
    intin("top", 0, "Limit to first n columns (0 = all)", min=0)],
   ex={"src": COUNTS, "method": "spearman"}, up="corrplot / Hmisc rcorr",
   tags=("correlation", "heatmap"),
   summary="Pairwise correlation coefficients for every column pair as a heatmap.")
def stats_correlation_matrix(src, method="pearson", annot=True, top=0):
    """Correlation heatmap."""
    _rows, cols, X = _mat(src, index_col=False)
    if X.shape[1] < 2:
        return plot.empty_plot("need at least two numeric columns")
    if top:
        X, cols = X[:, : int(top)], cols[: int(top)]
    n = len(cols)
    M = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            if i == j:
                M[i, j] = 1.0
                continue
            a, b = X[:, i], X[:, j]
            m = np.isfinite(a) & np.isfinite(b)
            if m.sum() < 3:
                continue
            fn = {"pearson": stats.pearson, "spearman": stats.spearman, "kendall": stats.kendall}[method]
            st = fn(list(a[m]), list(b[m]))
            M[i, j] = float(st.get("r", st.get("tau", 0.0)))
    return plot.heatmap([[round(float(v), 4) for v in row] for row in M], row_labels=cols, col_labels=cols,
                       title=f"{method} correlation", annot=bool(annot))


@T("stats_correlation_pair", "Correlate two columns with a p-value", STAT, "table",
   [tbl("src", COUNTS, "Table"), textbox("column_x", "", "Column X"), textbox("column_y", "", "Column Y"),
    choice("method", ["pearson", "spearman", "kendall"], "pearson", "Method")],
   ex={"src": COUNTS, "column_x": "sample_A", "column_y": "sample_B"},
   up="cor.test", tags=("correlation", "statistics"),
   summary="Correlation coefficient, test statistic and p-value for one pair of columns.")
def stats_correlation_pair(src, column_x="", column_y="", method="pearson"):
    """Single pairwise correlation."""
    df = _frame(src)
    cols = [c for c, _v in _num_cols(df)]
    cx = column_x or (cols[0] if cols else None)
    cy = column_y or (cols[1] if len(cols) > 1 else None)
    if not cx or not cy:
        return table([], "need two numeric columns")
    a = pd.to_numeric(df[cx], errors="coerce").to_numpy(dtype=float)
    b = pd.to_numeric(df[cy], errors="coerce").to_numpy(dtype=float)
    m = np.isfinite(a) & np.isfinite(b)
    st = {"pearson": stats.pearson, "spearman": stats.spearman, "kendall": stats.kendall}[method](list(a[m]), list(b[m]))
    rows = [{"metric": k, "value": round(float(v), 6) if isinstance(v, float) else v}
            for k, v in st.items()]
    return table(rows, f"{method}({cx}, {cy}) = {round(float(st.get('r', st.get('tau', 0.0))), 4)}, "
                       f"p = {round(float(st.get('p_value', 1.0)), 5)}")


@T("stats_regression", "Linear regression of y on x", STAT, "table",
   [tbl("src", COUNTS, "Table"), textbox("x_column", "", "Predictor"), textbox("y_column", "", "Response"),
    intin("order", 1, "Polynomial order", min=1, max=4), boolean("loess", False, "Also fit LOESS")],
   ex={"src": COUNTS, "x_column": "sample_A", "y_column": "sample_D"},
   up="lm / statsmodels OLS", tags=("regression", "statistics"),
   summary="Slope, intercept, R^2, residual standard error and the slope p-value.")
def stats_regression(src, x_column="", y_column="", order=1, loess=False):
    """Linear / polynomial regression."""
    df = _frame(src)
    cols = [c for c, _v in _num_cols(df)]
    cx = x_column or cols[0]
    cy = y_column or (cols[1] if len(cols) > 1 else cols[0])
    x = pd.to_numeric(df[cx], errors="coerce").to_numpy(dtype=float)
    y = pd.to_numeric(df[cy], errors="coerce").to_numpy(dtype=float)
    m = np.isfinite(x) & np.isfinite(y)
    if m.sum() < 3:
        return table([], "not enough numeric pairs")
    fit = stats.linear_regression(list(x[m]), list(y[m]), order=int(order))
    rows = [{"metric": k, "value": round(float(v), 6) if isinstance(v, float) else v}
            for k, v in fit.items() if k != "coefficients"]
    if str(fit.get("coefficients")):
        rows.append({"metric": "coefficients", "value": str(fit["coefficients"])[:120]})
    res = table(rows, f"R2 = {round(float(fit.get('r_squared', 0.0)), 4)}, "
                      f"slope = {round(float(fit.get('slope', 0.0)), 5)}")
    if loess:
        try:
            sm = stats.loess_fit(list(x[m]), list(y[m]))
            res["stats"] = {"loess_points": len(sm) if hasattr(sm, "__len__") else 0}
        except Exception:  # noqa: BLE001
            pass
    return res


@T("stats_contingency", "Chi-square and Fisher tests on a contingency table", STAT, "table",
   [tbl("src", COUNTS, "Counts table (first column = row labels)"),
    choice("test", ["chi_square", "fisher", "g_test"], "chi_square", "Test"),
    textbox("row_a", "", "Row A"), textbox("row_b", "", "Row B"),
    textbox("col_a", "", "Column A"), textbox("col_b", "", "Column B")],
   ex={"src": COUNTS, "test": "chi_square"}, up="chisq.test / fisher.test",
   tags=("contingency table", "statistics"),
   summary="Test association in a 2x2 slice or the whole table, with expected counts.")
def stats_contingency(src, test="chi_square", row_a="", row_b="", col_a="", col_b=""):
    """Contingency table testing."""
    df = _frame(src)
    num = [(c, pd.to_numeric(df[c], errors="coerce").fillna(0).to_numpy(dtype=float)) for c in df.columns
           if np.isfinite(pd.to_numeric(df[c], errors="coerce")).sum() >= 2]
    if len(num) < 2 or len(df) < 2:
        return table([], "need at least two numeric columns and two rows")
    rowlab = [str(v) for v in (df.iloc[:, 0].tolist() if df.iloc[:, 0].dtype == object
                              else range(len(df)))]
    idx = {n: i for i, n in enumerate(rowlab)}
    i1 = idx.get(str(row_a), 0) if row_a else 0
    i2 = idx.get(str(row_b), 1) if row_b else 1
    cnames = [c for c, _v in num]
    j1 = cnames.index(col_a) if col_a in cnames else 0
    j2 = cnames.index(col_b) if col_b in cnames else 1
    a = int(num[j1][1][i1])
    b = int(num[j1][1][i2])
    c = int(num[j2][1][i1])
    d = int(num[j2][1][i2])
    if test == "fisher":
        orr, p = stats.fisher_exact(a, b, c, d)
        return table([{"metric": "table", "value": f"[[{a},{b}],[{c},{d}]]"},
                     {"metric": "odds_ratio", "value": round(float(orr), 5)},
                     {"metric": "p_value", "value": round(float(p), 8)}],
                    f"Fisher exact p = {round(float(p), 6)}")
    obs = [a, b, c, d]
    exp = [a + c, b + d, a + b, c + d]
    chi = stats.chi_square(obs, exp)
    g_stat = 2.0 * sum(o * math.log(o / e) for o, e in zip(obs, exp) if o > 0 and e > 0)
    g_p = float(stats.chi2_pvalue(g_stat, 1)) if hasattr(stats, "chi2_pvalue") else 1.0
    g = {"statistic": g_stat, "p_value": g_p, "df": 1}
    rows = [{"test": "chi_square", "statistic": round(float(dict(chi).get("statistic", 0.0)), 5),
            "p_value": round(float(dict(chi).get("p_value", 1.0)), 6),
            "df": dict(chi).get("df", 1), "table": f"[[{a},{b}],[{c},{d}]]"},
           {"test": "g_test", "statistic": round(float(dict(g).get("statistic", 0.0)), 5),
            "p_value": round(float(dict(g).get("p_value", 1.0)), 6),
            "df": dict(g).get("df", 1), "table": f"[[{a},{b}],[{c},{d}]]"}]
    return table(rows, f"2x2 test on rows {rowlab[i1]}/{rowlab[i2]}, columns {cnames[j1]}/{cnames[j2]}")

@T("stats_anova", "One-way ANOVA and Kruskal-Wallis by group", STAT, "table",
   [tbl("src", PHENO, "Data table"), textbox("value_column", "", "Value column"),
    textbox("group_column", "", "Group column"), choice("test", ["anova", "kruskal", "levene", "bartlett"],
                                                        "anova", "Test")],
   ex={"src": PHENO, "value_column": "yield", "group_column": "group"},
   up="aov / kruskal.test / car::leveneTest", tags=("ANOVA", "variance", "statistics"),
   summary="Multi-group comparison with group means, F or H statistic and p-value.")
def stats_anova(src, value_column="", group_column="", test="anova"):
    """Multi-group tests."""
    df = _frame(src)
    cols = list(df.columns)
    vcol = value_column or next((c for c in cols if c != (group_column or cols[0])), cols[-1])
    gcol = group_column or cols[0]
    groups: dict[str, list[float]] = defaultdict(list)
    for r in df.to_dict("records"):
        try:
            groups[str(r[gcol])].append(float(r[vcol]))
        except (KeyError, TypeError, ValueError):
            continue
    if len(groups) < 2:
        return table([], "need at least two groups")
    data = [v for _k, v in sorted(groups.items())]
    fn = {"anova": stats.anova, "kruskal": stats.kruskal, "levene": stats.levene,
          "bartlett": stats.bartlett}[test]
    st = fn(*data) if test in ("anova", "kruskal") else fn(data)
    rows = [{"group": k, "n": len(v), "mean": round(float(np.mean(v)), 5),
            "stdev": round(float(np.std(v, ddof=1)), 5) if len(v) > 1 else 0.0}
            for k, v in sorted(groups.items())]
    res = table(rows, f"{test}: " + ", ".join(f"{k}={round(float(v), 5) if isinstance(v, float) else v}"
                                             for k, v in dict(st).items()))
    res["stats"] = {k: (round(float(v), 6) if isinstance(v, float) else v) for k, v in dict(st).items()}
    return res


@T("stats_p_adjust_table", "Multiple-testing correction of a p-value column", STAT, "table",
   [tbl("src", COUNTS, "Results table"), textbox("p_column", "", "p-value column"),
    choice("method", ["fdr_bh", "bonferroni", "holm", "holm_bonferroni"], "fdr_bh", "Method"),
    number("alpha", 0.05, "Significance threshold", min=0.0, max=1.0), intin("head", 300, "Rows", min=1)],
   ex={"src": COUNTS, "p_column": "sample_A", "method": "fdr_bh"},
   up="p.adjust / BH", tags=("multiple testing", "FDR"),
   summary="Adjust a column of p-values and flag the features that stay significant.")
def stats_p_adjust_table(src, p_column="", method="fdr_bh", alpha=0.05, head=300):
    """p-value adjustment."""
    df = _frame(src)
    cols = list(df.columns)
    pc = p_column or next((c for c in cols if "p" == str(c).lower()[:1] or "value" in str(c).lower()), cols[-1])
    raw = pd.to_numeric(df[pc], errors="coerce").to_numpy(dtype=float)
    finite = np.isfinite(raw)
    vals = [float(v) if math.isfinite(v) else 1.0 for v in raw]
    adj = stats.p_adjust(vals, method)
    rows = []
    for i, r in enumerate(df.to_dict("records")):
        rec = {k: v for k, v in r.items()}
        rec["p_value"] = round(float(raw[i]), 6) if finite[i] else float("nan")
        rec["p_adjusted"] = round(float(adj[i]), 6) if i < len(adj) else 1.0
        rec["significant"] = bool(finite[i] and adj[i] < float(alpha))
        rows.append(rec)
    return table(rows[:int(head)], f"{sum(1 for r in rows if r['significant'])}/{len(rows)} significant "
                                  f"after {method} at alpha={alpha}")


@T("stats_normality", "Normality and variance homogeneity tests", STAT, "table",
   [tbl("src", COUNTS, "Table"), intin("head", 30, "Columns to test", min=1),
    boolean("variance_tests", True, "Include Levene / Bartlett")],
   ex={"src": COUNTS}, up="shapiro / bartlett.test", tags=("normality", "assumption checks"),
   summary="Skew, kurtosis, Jarque-Bera p-value and equality-of-variance diagnostics per column.")
def stats_normality(src, head=30, variance_tests=True):
    """Normality diagnostics."""
    df = _frame(src)
    cols = _num_cols(df)[: int(head)]
    rows = []
    for name, v in cols:
        f = [float(x) for x in v[np.isfinite(v)]]
        if len(f) < 4:
            continue
        jb = stats.jarque_bera_p(f)
        rows.append({"column": name, "n": len(f), "skewness": round(float(stats.skewness(f)), 4),
                    "kurtosis": round(float(stats.kurtosis(f)), 4),
                    "jarque_bera_p": round(float(jb), 6) if math.isfinite(float(jb)) else "n/a",
                    "normal_at_0.05": bool(math.isfinite(float(jb)) and float(jb) > 0.05),
                    "mad": round(float(stats.mad(f)), 4)})
    res = table(rows, f"{len(rows)} columns checked")
    if variance_tests and len(cols) >= 2:
        data = [[float(x) for x in v[np.isfinite(v)]] for _c, v in cols]
        try:
            res["stats"]["levene"] = {k: (round(float(x), 5) if isinstance(x, float) else x)
                                     for k, x in dict(stats.levene(*data)).items()}
            res["stats"]["bartlett"] = {k: (round(float(x), 5) if isinstance(x, float) else x)
                                       for k, x in dict(stats.bartlett(*data)).items()}
        except (ValueError, TypeError, ZeroDivisionError):
            pass
    return res

@T("stats_distribution_plot", "Histogram with density overlay", STAT, "figure",
   [tbl("src", COUNTS, "Table"), textbox("column", "", "Column to plot"), intin("bins", 25, "Bins", min=3),
    boolean("log_y", False, "Log y axis"), boolean("multi", False, "Overlay all numeric columns")],
   ex={"src": COUNTS, "column": "sample_A", "bins": 12}, up="hist / ggplot2 density",
   tags=("distribution", "histogram", "plot"),
   summary="Empirical distribution of a value, or of every column for comparison.")
def stats_distribution_plot(src, column="", bins=25, log_y=False, multi=False):
    """Distribution plots."""
    df = _frame(src)
    cols = _num_cols(df)
    if not cols:
        return plot.empty_plot("no numeric columns")
    if multi:
        series = {n: [float(v) for v in v[np.isfinite(v)]] for n, v in cols[:8]}
        return plot.multi_histogram(series, bins=int(bins), xlabel="value", title="Distributions")
    name, v = next(((n, x) for n, x in cols if n == column), cols[0])
    f = v[np.isfinite(v)]
    return plot.histogram([float(x) for x in f], bins=int(bins), xlabel=name,
                         title=f"Distribution of {name}", log_y=bool(log_y))


@T("stats_qqplot", "Quantile-quantile plot against the normal distribution", STAT, "figure",
   [tbl("src", COUNTS, "Table"), textbox("column", "", "Column"), choice("dist", ["normal", "uniform",
                                                             "exponential"], "normal", "Reference distribution")],
   ex={"src": COUNTS, "column": "sample_A"}, up="qqnorm / qqplot", tags=("Q-Q plot", "diagnostics"),
   summary="Check normality of a metric visually - heavy tails bend away from the line.")
def stats_qqplot(src, column="", dist="normal"):
    """Q-Q plot."""
    df = _frame(src)
    cols = _num_cols(df)
    v = next((x for n, x in cols if n == column), cols[0][1] if cols else np.array([]))
    f = v[np.isfinite(v)]
    if f.size < 3:
        return plot.empty_plot("not enough values")
    return plot.qq_plot([float(x) for x in f], dist=dist, title=f"Q-Q plot vs {dist}")


@T("stats_two_sample_ks", "Kolmogorov-Smirnov test between columns", STAT, "table",
   [tbl("src", COUNTS, "Table"), textbox("column_a", "", "Column A"), textbox("column_b", "", "Column B")],
   ex={"src": COUNTS, "column_a": "sample_A", "column_b": "sample_D"}, up="ks.test",
   tags=("KS test", "distribution"),
   summary="Maximum ECDF difference and its p-value between two samples.")
def stats_two_sample_ks(src, column_a="", column_b=""):
    """Two-sample KS test."""
    df = _frame(src)
    cols = [c for c, _v in _num_cols(df)]
    ca = column_a or cols[0]
    cb = column_b or (cols[1] if len(cols) > 1 else cols[0])
    a = pd.to_numeric(df[ca], errors="coerce").to_numpy(dtype=float)
    b = pd.to_numeric(df[cb], errors="coerce").to_numpy(dtype=float)
    st = stats.kstest_2samp(list(a[np.isfinite(a)]), list(b[np.isfinite(b)]))
    rows = [{"metric": k, "value": round(float(v), 6) if isinstance(v, float) else v}
            for k, v in dict(st).items()]
    return table(rows, f"KS({ca}, {cb}) p = {round(float(dict(st).get('p_value', 1.0)), 5)}")


@T("stats_outlier_detection", "IQR, MAD and z-score outliers", STAT, "table",
   [tbl("src", COUNTS, "Table"), textbox("column", "", "Column"),
    number("iqr_factor", 1.5, "IQR factor", min=0.1), number("mad_factor", 3.5, "Modified z threshold", min=0.5),
    number("z_threshold", 3.0, "z-score threshold", min=1.0),
    choice("method", ["iqr", "mad", "zscore"], "iqr", "Method")],
   ex={"src": COUNTS, "column": "sample_A", "method": "mad"}, up="outlier detection (base R)",
   tags=("outliers", "quality control"),
   summary="Flag rows whose value deviates from the robust centre of the distribution.")
def stats_outlier_detection(src, column="", iqr_factor=1.5, mad_factor=3.5, z_threshold=3.0, method="iqr"):
    """Outlier flags."""
    df = _frame(src)
    cols = _num_cols(df)
    name, v = next(((n, x) for n, x in cols if n == column), cols[0] if cols else ("", np.array([])))
    if v.size == 0:
        return table([], "no numeric column found")
    f = v[np.isfinite(v)]
    if method == "iqr":
        flags = stats.outliers_iqr(v, k=float(iqr_factor))
    elif method == "mad":
        flags = stats.outliers_mad(v, k=float(mad_factor))
    else:
        z = (v - np.mean(f)) / (np.std(f) or 1.0)
        flags = {"low": [i for i in range(len(v)) if math.isfinite(v[i]) and v[i] < np.mean(f) - z_threshold * np.std(f)],
                "high": [i for i in range(len(v)) if math.isfinite(v[i]) and v[i] > np.mean(f) + z_threshold * np.std(f)]}
    bad = set(list(flags.get("low", [])) + list(flags.get("high", []))) if isinstance(flags, dict) else set(flags)
    lab = df.columns[0]
    rows = [{"row": str(df.iloc[i][lab]), "column": name, "value": round(float(v[i]), 5)
             if math.isfinite(v[i]) else "nan",
             "flag": "high" if (isinstance(flags, dict) and i in set(flags.get("high", []))) else "low",
             "method": method} for i in sorted(bad)]
    return table(rows, f"{len(rows)} outliers in {name} ({method}); median "
                       f"{round(float(np.median(f)), 4)}, MAD {round(float(stats.mad(f)), 4)}")


@T("stats_bootstrap_ci", "Bootstrap confidence interval of a statistic", STAT, "table",
   [tbl("src", COUNTS, "Table"), textbox("column", "", "Column"), choice("stat", ["mean", "median",
                                        "difference", "ratio"], "mean", "Statistic"),
    textbox("second_column", "", "Optional second column"), intin("resamples", 1000, "Bootstrap replicates", min=50),
    number("level", 0.95, "Confidence level", min=0.5, max=0.999), intin("seed", 1, "Seed", min=0)],
   ex={"src": COUNTS, "column": "sample_A", "resamples": 400},
   up="boot / simpleboot", tags=("bootstrap", "confidence interval"),
   summary="Resampling confidence interval and standard error for a column statistic.")
def stats_bootstrap_ci(src, column="", stat="mean", second_column="", resamples=1000, level=0.95, seed=1):
    """Bootstrap CI."""
    df = _frame(src)
    cols = _num_cols(df)
    name, v = next(((n, x) for n, x in cols if n == column), cols[0] if cols else ("", np.array([])))
    v = v[np.isfinite(v)]
    if v.size < 3:
        return table([], "need at least three values")
    if stat in ("difference", "ratio") and second_column:
        w = dict(cols).get(second_column, np.array([]))
        w = w[np.isfinite(w)]
        n = min(len(v), len(w))
        pairs = (v[:n], w[:n])
        fn = (lambda a, b: float(np.mean(a) - np.mean(b))) if stat == "difference" else \
            (lambda a, b: float(np.mean(a) / (np.mean(b) or 1e-9)))
        obs = fn(pairs[0], pairs[1])
        rng = np.random.default_rng(int(seed))
        boots = []
        for _ in range(int(resamples)):
            idx = rng.integers(0, n, n)
            boots.append(fn(pairs[0][idx], pairs[1][idx]))
    else:
        fn = np.mean if stat == "mean" else np.median
        obs = float(fn(v))
        rng = np.random.default_rng(int(seed))
        boots = [float(fn(v[rng.integers(0, v.size, v.size)])) for _ in range(int(resamples))]
    lo = float(np.quantile(boots, (1 - float(level)) / 2))
    hi = float(np.quantile(boots, 1 - (1 - float(level)) / 2))
    rows = [{"metric": "statistic", "value": stat}, {"metric": "estimate", "value": round(obs, 6)},
            {"metric": f"{int(level * 100)}%_ci_low", "value": round(lo, 6)},
            {"metric": f"{int(level * 100)}%_ci_high", "value": round(hi, 6)},
            {"metric": "bootstrap_se", "value": round(float(np.std(boots)), 6)},
            {"metric": "replicates", "value": int(resamples)}]
    return table(rows, f"bootstrap {stat} = {round(obs, 4)} [{round(lo, 4)}, {round(hi, 4)}]")


@T("stats_power_analysis", "Power / sample size calculation", STAT, "table",
   [number("effect_size", 0.5, "Cohen's d", min=0.01), number("alpha", 0.05, "Significance level", min=0.0001, max=0.5),
    number("power", 0.8, "Target power", min=0.5, max=0.999),
    intin("n", 0, "Fixed per-group size (0 = solve for n)", min=0)],
   ex={"effect_size": 0.8, "alpha": 0.05, "power": 0.8}, up="pwr.t.test / power.t.test",
   tags=("power analysis", "design"),
   summary="Required group size for a given effect, or the power of a planned experiment.")
def stats_power_analysis(effect_size=0.5, alpha=0.05, power=0.8, n=0):
    """Power analysis."""
    d = float(effect_size)
    st = stats.power_analysis(effect=d, alpha=float(alpha), power=float(power),
                            n=int(n) if int(n) > 0 else None)
    rows = [{"metric": k, "value": round(float(v), 6) if isinstance(v, float) else str(v)[:80]}
            for k, v in dict(st).items()]
    za, zb = stats.norm_ppf(1 - float(alpha) / 2), stats.norm_ppf(float(power))
    n_per = int(math.ceil(2 * ((za + zb) / (d or 1e-9)) ** 2))
    rows.append({"metric": "n_per_group_closed_form", "value": n_per})
    if int(n) > 0:
        z = abs(d) / math.sqrt(2.0 / float(n))
        achieved = float(stats.norm_cdf(z - za)) if hasattr(stats, "norm_cdf") else 0.0
        rows.append({"metric": "power_at_given_n", "value": round(max(0.0, min(1.0, achieved)), 5)})
    return table(rows, f"d = {d}, alpha = {alpha}, power = {power} -> about {n_per} per group")

@T("stats_auc_roc", "ROC curve and AUC for a classifier score", STAT, "table",
   [tbl("src", COUNTS, "Table"), textbox("score_column", "", "Score column"),
    textbox("truth_column", "", "Truth column (label)"),
    textbox("positive", "", "Positive class label")],
   ex={"src": COUNTS, "score_column": "sample_A", "truth_column": "sample_D"},
   up="pROC / roc()", tags=("ROC", "AUC", "classification"),
   summary="Sensitivity/specificity trade-off curve with the area under it.")
def stats_auc_roc(src, score_column="", truth_column="", positive=""):
    """ROC/AUC from a table."""
    df = _frame(src)
    cols = list(df.columns)
    sc = score_column or next((c for c, _v in _num_cols(df)), cols[0])
    tc = truth_column or next((c for c in cols if c != sc), sc)
    y = df[tc].tolist()
    uniq = sorted({str(v) for v in y})
    pos = positive or (uniq[-1] if len(uniq) > 1 else "1")
    lab = np.array([1.0 if str(v) == str(pos) else 0.0 for v in y])
    s = pd.to_numeric(df[sc], errors="coerce").to_numpy(dtype=float)
    m = np.isfinite(s)
    if lab[m].sum() == 0 or (1 - lab[m]).sum() == 0:
        return table([], "the truth column needs both classes")
    curve = stats.roc_curve(list(s[m]), list(lab[m]))
    auc = stats.auc(list(s[m]), list(lab[m])) if hasattr(stats, "auc") else float("nan")
    rows = [{"metric": "AUC", "value": round(float(auc), 5)},
            {"metric": "positives", "value": int(lab[m].sum())},
            {"metric": "negatives", "value": int((1 - lab[m]).sum())},
            {"metric": "points_on_curve", "value": len(curve) if hasattr(curve, "__len__") else 0},
            {"metric": "best_threshold", "value": str(_best_cut(curve))[:60]}]
    return table(rows, f"AUC = {round(float(auc), 4)} for {sc} vs {tc}")


def _best_cut(curve):
    """Threshold maximising Youden's J from an ROC curve table."""
    best, val = None, -1.0
    try:
        for r in curve:
            d = dict(r)
            j = float(d.get("tpr", d.get("sensitivity", 0.0))) - float(d.get("fpr", d.get("1-specificity", 0.0)))
            if j > val:
                val, best = j, d.get("threshold", d.get("cutoff"))
    except (TypeError, ValueError):
        return None
    return best


@T("stats_mutual_information", "Mutual information feature scores", STAT, "table",
   [tbl("src", COUNTS, "Feature matrix"), textbox("target_column", "", "Target column"),
    intin("bins", 10, "Histogram bins", min=2), intin("head", 100, "Rows", min=1)],
   ex={"src": COUNTS, "target_column": "sample_A", "bins": 4},
   up="sklearn mutual_info_regression", tags=("feature selection", "information theory"),
   summary="Information each column shares with a target, for nonlinear feature ranking.")
def stats_mutual_information(src, target_column="", bins=10, head=100):
    """Mutual information scores."""
    df = _frame(src)
    cols = _num_cols(df)
    names = [n for n, _v in cols]
    tname = target_column if target_column in names else (names[-1] if names else None)
    if not tname:
        return table([], "no numeric columns")
    t = dict(cols)[tname]
    rows = []
    for n, v in cols:
        if n == tname:
            continue
        m = np.isfinite(v) & np.isfinite(t)
        if m.sum() < 4:
            continue
        mi = float(stats.mutual_information(list(v[m]), list(t[m]), bins=int(bins)))
        hv = float(stats.shannon([float(x) for x in np.histogram(v[m], bins=int(bins))[0]]))
        rows.append({"feature": n, "target": tname, "mutual_information": round(mi, 6),
                    "feature_entropy": round(hv, 5),
                    "normalized": round(mi / hv, 5) if hv > 1e-12 else 0.0})
    rows.sort(key=lambda r: -r["mutual_information"])
    return table(rows[:int(head)], f"mutual information against {tname} ({len(rows)} features)")

@T("stats_time_series_trend", "Trend, smoothing and autocorrelation", STAT, "table",
   [tbl("src", COUNTS, "Table with a value column"), textbox("value_column", "", "Value column"),
    number("alpha", 0.3, "EWMA alpha", min=0.01, max=1.0), intin("lag", 1, "Autocorrelation lag", min=1),
    boolean("cusum", True, "Also run CUSUM change detection")],
   ex={"src": COUNTS, "value_column": "sample_A"}, up="trend / cusum (statistics)",
   tags=("time series", "smoothing", "trend"),
   summary="EWMA smoothing, Durbin-Watson residual statistic and cumulative-sum change points.")
def stats_time_series_trend(src, value_column="", alpha=0.3, lag=1, cusum=True):
    """Trend diagnostics."""
    df = _frame(src)
    cols = [c for c, _v in _num_cols(df)]
    vc = value_column if value_column in cols else (cols[0] if cols else None)
    if not vc:
        return table([], "no numeric column")
    v = np.nan_to_num(pd.to_numeric(df[vc], errors="coerce").to_numpy(dtype=float), nan=0.0,
                     posinf=0.0, neginf=0.0)
    sm = np.asarray(stats.ewma(list(v), alpha=float(alpha)), dtype=float)
    resid = v - sm
    slope = float(np.polyfit(np.arange(v.size, dtype=float), v, 1)[0]) if v.size > 2 else 0.0
    dw = float(stats.durbin_watson(list(resid))) if v.size > 2 else float("nan")
    cc = stats.cross_correlation(list(v), list(v), max_lag=max(1, int(lag)))
    lag_val = 0.0
    for k, rec in enumerate(cc or []):
        if int(dict(rec).get("lag", k)) == int(lag):
            lag_val = float(dict(rec).get("ccf", dict(rec).get("correlation", 0.0)) or 0.0)
            break
    rows = [{"metric": "column", "value": vc}, {"metric": "n", "value": int(v.size)},
            {"metric": "trend_slope", "value": round(slope, 6)},
            {"metric": "first_value", "value": round(float(v[0]), 5)},
            {"metric": "last_value", "value": round(float(v[-1]), 5)},
            {"metric": "ewma_last", "value": round(float(sm[-1]), 5)},
            {"metric": f"autocorrelation_lag{int(lag)}", "value": round(lag_val, 5) if math.isfinite(lag_val) else 0.0},
            {"metric": "durbin_watson", "value": round(dw, 5) if math.isfinite(dw) else "n/a"}]
    if cusum:
        cs = np.asarray(stats.cusum(list(v)), dtype=float)
        rows.append({"metric": "cusum_max_deviation", "value": round(float(np.max(np.abs(cs))) if cs.size else 0.0, 5)})
    return table(rows, f"trend of {vc}: slope {round(slope, 5)}")

@T("stats_summary_of_bedgraph", "Distribution of a coverage track", STAT, "table",
   [anyfile("graph", BEDGRAPH, fmt="", label="bedGraph / WIG track"), intin("bins", 20, "Histogram bins", min=3),
    choice("weight", ["none", "length"], "length", "Weight by interval length")],
   ex={"graph": BEDGRAPH, "bins": 12}, up="deepTools plotCoverage / summary",
   tags=("coverage", "distribution"),
   summary="Signal statistics of a genomic track: mean, percentiles and zero fraction.")
def stats_summary_of_bedgraph(graph, bins=20, weight="length"):
    """Coverage track summary."""
    vals = io.parse_bedgraph(io.as_text(graph))
    if not vals:
        return table([], "no bedGraph values parsed")
    v = np.array([x[3] for x in vals], dtype=float)
    w = np.array([x[2] - x[1] for x in vals], dtype=float) if weight == "length" else np.ones(len(v))
    wmean = float(np.sum(v * w) / max(1e-9, np.sum(w)))
    rows = [{"metric": "intervals", "value": len(v)},
            {"metric": "total_bp", "value": int(np.sum(w))},
            {"metric": "chromosomes", "value": len({x[0] for x in vals})},
            {"metric": "weighted_mean", "value": round(wmean, 5)},
            {"metric": "mean", "value": round(float(np.mean(v)), 5)},
            {"metric": "median", "value": round(float(np.median(v)), 5)},
            {"metric": "max", "value": round(float(np.max(v)), 5)},
            {"metric": "zero_fraction", "value": round(float(np.mean(v == 0)), 4)},
            {"metric": "q95", "value": round(float(np.quantile(v, 0.95)), 5)}]
    res = table(rows, f"{len(v)} intervals, weighted mean {round(wmean, 4)}")
    hist, edges = np.histogram(v, bins=int(bins))
    res["stats"]["histogram"] = {f"{round(float(edges[i]), 3)}-{round(float(edges[i + 1]), 3)}": int(hist[i])
                                for i in range(len(hist))}
    return res


# ===========================================================================
# machine learning
# ===========================================================================
def _xy(src, xcols, ycol):
    rows_, cols_, X = _mat(src)
    if xcols:
        keep = [cols_.index(c) for c in str(xcols).split(",") if c in cols_]
        if keep:
            X = X[:, keep]
    df = _frame(src)
    ycol_list = None
    if ycol:
        y_raw = df[ycol].tolist() if ycol in df.columns else rows_
        cats = sorted({str(v) for v in y_raw})
        mp = {c: i for i, c in enumerate(cats)}
        ycol_list = np.array([mp[str(v)] for v in y_raw], dtype=float)
    else:
        ycol_list = np.arange(X.shape[0], dtype=float)
    X = np.nan_to_num(X, nan=0.0, posinf=0.0, neginf=0.0)
    return rows_, X, ycol_list, cats if ycol else []


@T("ml_kmeans", "k-means clustering of samples or features", MLSEC, "table",
   [tbl("src", COUNTS, "Matrix"), intin("k", 3, "Clusters", min=2), intin("iters", 100, "Max iterations", min=5),
    choice("transpose", ["rows", "columns"], "rows", "Cluster"), intin("seed", 1, "Seed", min=0),
    boolean("standardize", True, "z-score features")],
   ex={"src": COUNTS, "k": 3, "transpose": "columns"}, up="kmeans / scikit-learn KMeans",
   tags=("clustering", "k-means"),
   summary="Lloyd's algorithm with inertia, silhouette score and cluster sizes.")
def ml_kmeans(src, k=3, iters=100, transpose="rows", seed=1, standardize=True):
    """k-means clustering."""
    rows_, cols_, X = _mat(src)
    labels_names = rows_ if transpose == "rows" else cols_
    M = X if transpose == "rows" else X.T
    if M.shape[0] < 2:
        return table([], "need at least two observations")
    if standardize:
        mu, sd = np.nanmean(M, axis=0), np.nanstd(M, axis=0)
        M = (M - mu) / np.where(sd > 0, sd, 1)
    res = ml.kmeans(np.nan_to_num(M), k=int(k), iters=int(iters), seed=int(seed))
    lab = np.asarray(res["labels"], dtype=int)
    per = res.get("silhouette")
    pers = np.atleast_1d(np.asarray(per, dtype=float)) if per is not None else np.array([])
    rows = [{"observation": str(n), "cluster": int(l),
             "centroid_distance": round(float(np.linalg.norm(M[i] - np.asarray(res["centroids"])[l])), 4),
             "silhouette": round(float(pers[i]), 4) if i < pers.size else ""}
            for i, (n, l) in enumerate(zip(labels_names, lab))]
    mean_sil = round(float(np.mean(pers)), 4) if pers.size else "n/a"
    out = table(rows, f"{len(set(lab.tolist()))} clusters, mean silhouette {mean_sil}")
    out["stats"] = {"inertia": res.get("inertia"), "cluster_sizes": res.get("cluster_sizes"),
                   "k": int(k)}
    return out


@T("ml_dbscan", "DBSCAN density clustering with eps curve", MLSEC, "table",
   [tbl("src", COUNTS, "Matrix"), number("eps", 1.0, "Neighbourhood radius", min=0.0),
    intin("min_samples", 3, "Minimum samples", min=1), boolean("standardize", True, "z-score features")],
   ex={"src": COUNTS, "eps": 2.0, "min_samples": 2}, up="cluster::dbscan / sklearn DBSCAN",
   tags=("clustering", "density"),
   summary="Density-based clustering that leaves noise unassigned, with an eps diagnostic curve.")
def ml_dbscan(src, eps=1.0, min_samples=3, standardize=True):
    """DBSCAN clustering."""
    rows_, cols_, X = _mat(src)
    M = np.nan_to_num(X.astype(float))
    if standardize:
        mu, sd = np.nanmean(M, axis=0), np.nanstd(M, axis=0)
        M = (M - mu) / np.where(sd > 0, sd, 1)
    res = ml.dbscan(M, eps=float(eps), min_samples=int(min_samples))
    lab = np.asarray(res["labels"], dtype=int)
    rows = [{"observation": n, "cluster": int(l), "noise": bool(l < 0)} for n, l in zip(rows_, lab)]
    curve = ml.dbscan_eps_curve(M, k=max(2, int(min_samples))) if hasattr(ml, "dbscan_eps_curve") else []
    out = table(rows, f"{res.get('n_clusters', len(set(lab)))} clusters, {res.get('noise', 0)} noise points")
    out["stats"] = {"sizes": res.get("sizes"), "eps_curve_points": len(curve)}
    return out


@T("ml_hierarchical", "Hierarchical clustering with cophenetic correlation", MLSEC, "figure",
   [tbl("src", COUNTS, "Matrix"), choice("method", ["average", "complete", "single", "ward"], "average",
                                        "Linkage"), intin("clusters", 3, "Cut into k clusters", min=2),
    boolean("rows", True, "Cluster rows (else columns)")],
   ex={"src": COUNTS, "method": "ward", "clusters": 3}, up="hclust / scipy linkage",
   tags=("clustering", "dendrogram"),
   summary="Linkage dendrogram of the observations with the cut height reported.")
def ml_hierarchical(src, method="average", clusters=3, rows=True):
    """Dendrogram + cut."""
    rnames, cols, X = _mat(src)
    M = X if rows else X.T
    names = rnames if rows else cols
    if M.shape[0] < 3:
        return plot.empty_plot("need at least three observations")
    Z = ml.linkage(np.nan_to_num(M.astype(float)), method=method, labels=[str(n) for n in names])
    lab = ml.fcluster(Z, t=max(2, int(clusters)), criterion="maxclust")
    return plot.dendrogram(Z, labels=[str(n) for n in names],
                          title=f"{method} linkage ({len(set(lab))} clusters at maxclust={clusters})")


@T("ml_gmm_bic", "Gaussian mixture model selection by BIC", MLSEC, "table",
   [tbl("src", COUNTS, "Matrix"), intin("max_components", 4, "Maximum components", min=2),
    boolean("standardize", True, "z-score features")],
   ex={"src": COUNTS, "max_components": 3}, up="mclust / GaussianMixture",
   tags=("mixture model", "BIC", "clustering"),
   summary="Fit mixtures of increasing complexity and pick the BIC optimum.")
def ml_gmm_bic(src, max_components=4, standardize=True):
    """GMM with BIC selection."""
    rows_, cols_, X = _mat(src)
    M = np.nan_to_num(X.astype(float))
    if standardize:
        mu, sd = np.nanmean(M, axis=0), np.nanstd(M, axis=0)
        M = (M - mu) / np.where(sd > 0, sd, 1)
    curves = ml.gaussian_mixture_bic(M, max_components=int(max_components))
    rows = []
    for c in curves:
        rows.append({**{k: (round(float(v), 5) if isinstance(v, (int, float, np.floating)) else v)
                       for k, v in dict(c).items()}})
    best = min((c for c in curves if "bic" in dict(c)), key=lambda c: float(dict(c)["bic"]), default=None)
    out = table(rows, f"best model: {dict(best).get('components') if best else '?'} components")
    if best is not None and int(dict(best).get("components", 0)) >= 2:
        fit = ml.gmm_em(M, components=int(dict(best)["components"]))
        out["stats"] = {"labels": [int(x) for x in np.asarray(fit.get("labels", []), dtype=int)]}
    return out


@T("ml_pca", "Principal component analysis", MLSEC, "table",
   [tbl("src", COUNTS, "Matrix"), intin("components", 2, "Components", min=1),
    boolean("standardize", True, "z-score features"), intin("head", 200, "Rows", min=1)],
   ex={"src": COUNTS, "components": 3}, up="prcomp / sklearn PCA", tags=("PCA", "dim reduction"),
   summary="Scores, loadings and explained variance of a principal component analysis.")
def ml_pca(src, components=2, standardize=True, head=200):
    """PCA of a matrix."""
    rows_, cols_, X = _mat(src)
    M = np.nan_to_num(X.astype(float))
    if standardize:
        mu, sd = np.nanmean(M, axis=0), np.nanstd(M, axis=0)
        M = (M - mu) / np.where(sd > 0, sd, 1)
    res = ml.pca(M, components=int(components))
    sc = np.asarray(res["scores"], dtype=float)
    ratio = np.asarray(res.get("explained_variance_ratio", []), dtype=float)
    rows = []
    for i, name in enumerate(rows_):
        rec = {"observation": name}
        for j in range(sc.shape[1]):
            rec[f"PC{j + 1}"] = round(float(sc[i, j]), 5)
        rows.append(rec)
    load = np.atleast_2d(np.asarray(res.get("loadings", []), dtype=float))
    out = table(rows[:int(head)], "PC1 explains "
                                 f"{round(100 * (ratio[0] if ratio.size else 0), 1)}% of the variance")
    if load.size:
        nfeat = len(cols_)
        if load.shape[1] != nfeat and load.shape[0] == nfeat:
            load = load.T
        ld = {f"PC{j + 1}": {cols_[k]: round(float(load[j, k]), 4) for k in range(min(nfeat, load.shape[1]))}
             for j in range(load.shape[0])}
    else:
        ld = {}
    out["stats"] = {"explained_variance_ratio": [round(float(v), 5) for v in ratio], "loadings": ld}
    return out


@T("ml_embedding_plot", "t-SNE / UMAP / MDS embedding", MLSEC, "figure",
   [tbl("src", COUNTS, "Matrix"), choice("method", ["tsne", "umap", "mds"], "umap", "Method"),
    intin("components", 2, "Dimensions", min=2, max=3), intin("neighbours", 5, "Perplexity / neighbours", min=2),
    textbox("label_column", "", "Colour by a column")],
   ex={"src": COUNTS, "method": "tsne", "neighbours": 3}, up="Rtsne / uwot / cmdscale",
   tags=("t-SNE", "UMAP", "embedding"),
   summary="Non-linear 2-D projection of a feature matrix, optionally coloured by group.")
def ml_embedding_plot(src, method="umap", components=2, neighbours=5, label_column=""):
    """2-D embedding scatter."""
    df = _frame(src)
    rows_, cols_, X = _mat(src)
    M = np.nan_to_num(X.astype(float))
    if M.shape[0] < 4:
        return plot.empty_plot("need at least four observations")
    if method == "tsne":
        res = ml.tsne(M, components=int(components),
                     perplexity=max(1.0, min(float(neighbours), (M.shape[0] - 1) / 3.0)))
    elif method == "mds":
        D = np.zeros((M.shape[0], M.shape[0]))
        for i in range(M.shape[0]):
            D[i] = np.linalg.norm(M - M[i], axis=1)
        res = ml.mds(D, components=int(components))
    else:
        res = ml.umap_like(M, components=int(components), neighbours=int(neighbours))
    E = np.asarray(res.get("embedding", []), dtype=float)
    if E.size == 0 or E.shape[0] != M.shape[0]:
        return plot.empty_plot("embedding did not converge")
    colours = None
    if label_column and label_column in df.columns:
        cats = sorted({str(v) for v in df[label_column].tolist()})
        mp = {c: i for i, c in enumerate(cats)}
        colours = [float(mp[str(v)]) for v in df[label_column].tolist()]
    return plot.scatter(list(E[:, 0]), list(E[:, 1]), labels=[str(n) for n in rows_], colour_by=colours,
                       xlabel="dimension 1", ylabel="dimension 2", title=f"{method} embedding")

@T("ml_classifier", "Train a classifier (tree, random forest, logistic, kNN, Naive Bayes)", MLSEC,
   "table",
   [tbl("src", COUNTS, "Matrix"), textbox("label_column", "", "Column with class labels"),
    choice("model", ["decision_tree", "random_forest", "logistic", "naive_bayes", "knn", "svm"],
           "random_forest", "Model"), intin("test_percent", 25, "Test split %", min=0, max=60),
    intin("depth", 3, "Max depth / k / neighbours", min=1), intin("seed", 1, "Seed", min=0)],
   ex={"src": COUNTS, "label_column": "sample_A", "model": "decision_tree", "depth": 2},
   up="caret train / sklearn", tags=("classification", "model evaluation"),
   summary="Train a model on a feature matrix and report accuracy, confusion counts and importances.")
def ml_classifier(src, label_column="", model="random_forest", test_percent=25, depth=3, seed=1):
    """Classifier training and evaluation."""
    df = _frame(src)
    rows_, cols_, X = _mat(src)
    M = np.nan_to_num(X.astype(float))
    if M.shape[1] == 0 or M.shape[0] < 4:
        return table([], "need at least four rows and one numeric column")
    if label_column and label_column in df.columns:
        raw = [str(v) for v in df[label_column].tolist()]
    else:
        med = np.nanmedian(M[:, 0])
        raw = ["high" if v > med else "low" for v in M[:, 0]]
    cats = sorted(set(raw))
    if len(cats) < 2:
        return table([], "the label column must contain at least two classes")
    mp = {c: i for i, c in enumerate(cats)}
    y = np.array([float(mp[v]) for v in raw])
    split = max(2, int(round(len(M) * (1 - float(test_percent) / 100.0))))
    tr, te = slice(0, split), slice(split, len(M))
    held_out = te.start < len(M)
    if model == "decision_tree":
        res = ml.decision_tree(M[tr], y[tr], max_depth=int(depth))
    elif model == "random_forest":
        res = ml.random_forest(M[tr], y[tr], trees=10, max_depth=int(depth), seed=int(seed))
    elif model == "logistic":
        res = ml.logistic_regression(M[tr], y[tr])
    elif model == "naive_bayes":
        res = ml.naive_bayes(M[tr], y[tr], M[te] if held_out else M[tr])
    elif model == "knn":
        res = ml.knn_classify(M[tr], y[tr], M[te] if held_out else M[tr], k=max(1, int(depth)))
    else:
        res = ml.predict_tree(ml.svm_sgd(M[tr], y[tr]), M[te] if held_out else M[tr]) \
            if hasattr(ml, "predict_tree") and False else ml.svm_sgd(M[tr], y[tr])
    pred = np.asarray(dict(res).get("predictions", []), dtype=float).ravel()
    truth = y[te] if (held_out and len(pred) == len(y[te])) else y[tr][: len(pred)]
    acc = float(np.mean([1.0 if round(float(p)) == float(t) else 0.0 for p, t in zip(pred, truth)])) \
        if len(pred) else 0.0
    conf = stats.confusion_matrix([float(v) for v in truth], [round(float(p)) for p in pred],
                                 labels=[float(mp[c]) for c in cats])
    rows = [{"metric": "model", "value": model}, {"metric": "n_train", "value": int(split)},
            {"metric": "n_test", "value": int(len(M) - split) if held_out else 0},
            {"metric": "classes", "value": ",".join(cats)},
            {"metric": "accuracy_evaluated_on", "value": "test" if held_out and len(pred) == len(y[te]) else "train"},
            {"metric": "accuracy", "value": round(acc, 5)},
            {"metric": "training_accuracy",
             "value": round(float(dict(res).get("accuracy", dict(res).get("training_accuracy", 0.0)) or 0.0), 5)},
            {"metric": "oob_accuracy", "value": round(float(dict(res).get("oob_accuracy", 0.0) or 0.0), 5)},
            {"metric": "auc", "value": round(float(dict(res).get("auc", 0.0) or 0.0), 5)}]
    imp = np.ravel(np.asarray(dict(res).get("feature_importances", dict(res).get("coefficients", [])), dtype=float))
    out = table(rows, f"{model}: accuracy {round(100 * acc, 1)}% on {len(pred)} scored rows")
    out["stats"] = {"confusion": conf,
                   "importances": {cols_[i]: round(float(v), 5) for i, v in enumerate(imp) if i < len(cols_)}}
    return out

@T("ml_regression", "Regression with cross-validation", MLSEC, "table",
   [tbl("src", COUNTS, "Matrix"), textbox("target_column", "", "Target column"),
    choice("model", ["linear", "ridge", "lasso", "gradient_boosting", "regression_tree"], "ridge", "Model"),
    intin("folds", 4, "Cross-validation folds", min=2), number("alpha", 1.0, "Regularisation", min=0.0),
    intin("trees", 40, "Trees / iterations", min=2)],
   ex={"src": COUNTS, "target_column": "sample_A", "model": "linear"},
   up="glmnet / xgboost / caret train", tags=("regression", "cross-validation"),
   summary="Fit a continuous target and compare in-sample with cross-validated error.")
def ml_regression(src, target_column="", model="ridge", folds=4, alpha=1.0, trees=40):
    """Regression with CV."""
    rows_, cols_, X = _mat(src)
    df = _frame(src)
    M = np.nan_to_num(X.astype(float))
    y = np.arange(M.shape[0], dtype=float)
    feat_names = list(cols_)
    if target_column and target_column in cols_:
        j = cols_.index(target_column)
        y = M[:, j].copy()
        M = np.delete(M, j, axis=1)
        feat_names = [c for i, c in enumerate(cols_) if i != j]
    if M.shape[1] == 0 or M.shape[0] < 3:
        return table([], "need at least three rows and one predictor")
    fit = _regression_fit(model, M, y, float(alpha), int(trees))
    pred = np.asarray(fit.get("predictions", np.zeros(len(y))), dtype=float).ravel()
    resid = y - (pred[: len(y)] if len(pred) == len(y) else np.zeros(len(y)))
    ss = float(np.sum(resid ** 2))
    tss = float(np.sum((y - np.mean(y)) ** 2)) or 1.0
    coefs = np.ravel(np.asarray(fit.get("coefficients", []), dtype=float)) if "coefficients" in fit else np.array([])
    rows = []
    for i, n in enumerate(feat_names):
        rows.append({"feature": n, "coefficient": round(float(coefs[i + 1]), 6) if i + 1 < len(coefs) else 0.0})
    if not rows and "tree" in str(fit):
        rows = [{"feature": "tree_nodes", "coefficient": float(_count_nodes(fit.get("tree", fit)))}]
    rows += [{"feature": "_intercept", "coefficient": round(float(coefs[0]), 6) if len(coefs) else 0.0},
            {"feature": "_R_squared_in_sample", "coefficient": round(1 - ss / tss, 5)},
            {"feature": "_RMSE_in_sample", "coefficient": round(float(np.sqrt(ss / max(1, len(y)))), 5)}]
    cv = {}
    if model in ("linear", "ridge", "lasso", "regression_tree"):
        def fit_predict(Xtr, ytr, Xte, yte):
            f = _regression_fit(model, Xtr, np.asarray(ytr, dtype=float), float(alpha), int(trees))
            p = _regress_predict(model, Xtr, np.asarray(ytr, dtype=float), Xte, f)
            return float(np.mean((p - np.asarray(yte, dtype=float)) ** 2))
        try:
            cv = ml.cross_validate(fit_predict, M, y, folds=int(folds))
        except (ValueError, TypeError, ZeroDivisionError):
            cv = {}
    if cv:
        rows.append({"feature": "_cv_mean_squared_error",
                    "coefficient": round(float(dict(cv).get("mean", 0.0)), 5)})
    out = table(rows, f"{model}: in-sample R2 {round(1 - ss / tss, 4)}"
                      + (f", CV MSE {round(float(dict(cv).get('mean', 0.0)), 4)}" if cv else ""))
    out["stats"] = {**({"cross_validation": {k: (round(float(v), 5) if isinstance(v, float) else v)
                                            for k, v in dict(cv).items()}} if cv else {}),
                   "n": int(len(y)), "features": feat_names}
    return out


def _regression_fit(model, X, y, alpha, trees):
    if model == "linear":
        return ml.linear_regression_fit(X, y)
    if model == "lasso":
        return ml.lasso_regression(X, y, alpha=alpha, iters=max(2, int(trees)))
    if model == "gradient_boosting":
        return ml.gradient_boosting(X, y, trees=max(2, int(trees)))
    if model == "regression_tree":
        return {"tree": ml.regression_tree(X, y, max_depth=3, min_samples=2)}
    return ml.ridge_regression(X, y, alpha=alpha)


def _regress_predict(model, Xtr, ytr, Xte, fit=None):
    Xte = np.atleast_2d(np.asarray(Xte, dtype=float))
    if model == "regression_tree":
        tree = (fit or {}).get("tree") or ml.regression_tree(Xtr, ytr, max_depth=3, min_samples=2)
        return np.array([_tree_walk(tree, row) for row in Xte], dtype=float)
    if fit is None:
        fit = _regression_fit(model, Xtr, ytr, 1.0, 10)
    coef = np.ravel(np.asarray(fit.get("coefficients", []), dtype=float))
    if coef.size == 0:
        return np.zeros(len(Xte))
    if coef.size == Xte.shape[1] + 1:
        return coef[0] + Xte @ coef[1:]
    return Xte @ coef[: Xte.shape[1]]


def _tree_walk(node, row):
    n = dict(node or {})
    guard = 0
    while "leaf" not in n and guard < 200:
        guard += 1
        f = int(n.get("feature", 0))
        v = float(row[f]) if f < len(row) else 0.0
        n = dict(n["left"] if v <= float(n.get("threshold", 0.0)) else n["right"])
    return float(n.get("leaf", 0.0))


def _count_nodes(node):
    n = dict(node or {})
    if "leaf" in n:
        return 1
    return 1 + _count_nodes(n.get("left")) + _count_nodes(n.get("right"))

@T("ml_cluster_validation", "Silhouette and gap-style cluster validation", MLSEC, "table",
   [tbl("src", COUNTS, "Matrix"), intin("k_min", 2, "Smallest k", min=2), intin("k_max", 6, "Largest k", min=2),
    choice("metric", ["silhouette", "inertia", "both"], "both", "Metric"), intin("seed", 1, "Seed", min=0)],
   ex={"src": COUNTS, "k_min": 2, "k_max": 4}, up="NbClust / clusterCrit",
   tags=("clustering", "validation"),
   summary="Sweep the number of clusters and report silhouette and inertia for each k.")
def ml_cluster_validation(src, k_min=2, k_max=6, metric="both", seed=1):
    """Cluster number sweep."""
    rows_, cols_, X = _mat(src)
    M = np.nan_to_num(X.astype(float))
    mu, sd = np.nanmean(M, axis=0), np.nanstd(M, axis=0)
    M = (M - mu) / np.where(sd > 0, sd, 1)
    rows = []
    for k in range(int(k_min), min(int(k_max), max(2, M.shape[0] - 1)) + 1):
        res = ml.kmeans(M, k=k, seed=int(seed))
        lab = np.asarray(res["labels"], dtype=int)
        sil = res.get("silhouette")
        sil_val = float(np.mean(sil)) if sil is not None and np.ndim(sil) else float(sil or 0.0)
        rows.append({"k": k, "silhouette": round(sil_val, 5), "inertia": round(float(res.get("inertia", 0.0)), 5),
                    "cluster_sizes": ",".join(str(v) for v in np.bincount(lab)),
                    "empty_clusters": int(k - len(set(lab.tolist())))})
    best = max(rows, key=lambda r: r["silhouette"]) if rows else {}
    out = table(rows, f"best k by silhouette: {best.get('k', 'n/a')} ({best.get('silhouette', 0)})")
    out["stats"] = {"metric": metric, "scores": {r["k"]: [r["silhouette"], r["inertia"]] for r in rows}}
    return out


@T("ml_anomaly_detection", "z-score anomaly detection per feature", MLSEC, "table",
   [tbl("src", COUNTS, "Matrix"), number("threshold", 3.0, "z threshold", min=1.0),
    choice("robust", ["median", "mean"], "median", "Centre"), intin("head", 200, "Rows", min=1)],
   ex={"src": COUNTS, "threshold": 1.5}, up="anomaly detection (statistics)",
   tags=("anomaly", "outliers"),
   summary="Flag observations that are unusual in any feature using robust z-scores.")
def ml_anomaly_detection(src, threshold=3.0, robust="median", head=200):
    """Anomaly flags over a matrix."""
    rows_, cols_, X = _mat(src)
    M = np.nan_to_num(X.astype(float))
    if M.size == 0:
        return table([], "empty matrix")
    zs = np.zeros_like(M)
    for j in range(M.shape[1]):
        col = M[:, j]
        if robust == "median":
            centre = float(np.nanmedian(col))
            spread = 1.4826 * float(stats.mad([float(v) for v in col])) or (float(np.nanstd(col)) or 1.0)
        else:
            centre, spread = float(np.nanmean(col)), (float(np.nanstd(col)) or 1.0)
        zs[:, j] = (col - centre) / (spread or 1.0)
    rows = []
    for i, name in enumerate(rows_):
        j = int(np.nanargmax(np.abs(zs[i])))
        score = float(np.nanmax(np.abs(zs[i])))
        if score >= float(threshold):
            rows.append({"observation": name, "max_abs_z": round(score, 4), "feature": cols_[j],
                        "n_flags": int(np.sum(np.abs(zs[i]) >= float(threshold))),
                        "mean_abs_z": round(float(np.nanmean(np.abs(zs[i]))), 4)})
    rows.sort(key=lambda r: -r["max_abs_z"])
    out = table(rows[:int(head)], f"{len(rows)} anomalous observations (|z| >= {threshold})")
    out["stats"] = {"scores": {str(n): round(float(np.nanmax(np.abs(zs[i]))), 4)
                              for i, n in enumerate(rows_)}}
    return out

@T("ml_split_and_encode", "Train/test split and one-hot preview", MLSEC, "table",
   [tbl("src", PHENO, "Data table"), textbox("label_column", "", "Label column"),
    number("test_fraction", 0.25, "Test fraction", min=0.05, max=0.5), intin("seed", 7, "Seed", min=0),
    intin("head", 200, "Rows", min=1)],
   ex={"src": PHENO, "label_column": "group"}, up="caret createDataPartition / fastDummies",
   tags=("preprocessing", "split"),
   summary="Split a table into training and test parts and encode its categorical column.")
def ml_split_and_encode(src, label_column="", test_fraction=0.25, seed=7, head=200):
    """Split + encoding preview."""
    df = _frame(src)
    cols = list(df.columns)
    lc = label_column if label_column in cols else cols[0]
    labs = [str(v) for v in df[lc].tolist()]
    encoded = ml.one_hot(labs)
    n = len(df)
    order = list(range(n))
    rng = np.random.default_rng(int(seed))
    rng.shuffle(order)
    cut = int(round(n * (1 - float(test_fraction))))
    rows = [{"row": int(i), "split": "train" if k < cut else "test", "label": labs[i],
            "numeric_columns": len([c for c, _v in _num_cols(df)])}
           for k, i in enumerate(order)]
    out = table(rows[:int(head)], f"{cut} training / {n - cut} test rows; "
                                f"{len(encoded[0]) if encoded else 0} one-hot columns for {lc}")
    out["stats"] = {"one_hot_preview": encoded[:10], "columns": cols}
    return out

@T("ml_feature_selection_cv", "Permutation importance for a fitted model", MLSEC, "table",
   [tbl("src", COUNTS, "Matrix"), textbox("target_column", "", "Target column"),
    choice("model", ["linear", "ridge", "tree"], "linear", "Base model"),
    intin("repeats", 10, "Permutations per feature", min=2), intin("seed", 2, "Seed", min=0)],
   ex={"src": COUNTS, "target_column": "sample_A", "model": "linear"},
   up="caret varImp / permutation importance", tags=("feature selection", "importance"),
   summary="How much each feature matters, measured by the error increase when shuffled.")
def ml_feature_selection_cv(src, target_column="", model="linear", repeats=10, seed=2):
    """Permutation importance."""
    rows_, cols_, X = _mat(src)
    M = np.nan_to_num(X.astype(float))
    y = np.arange(M.shape[0], dtype=float)
    feat_names = list(cols_)
    if target_column and target_column in cols_:
        j = cols_.index(target_column)
        y = M[:, j].copy()
        M = np.delete(M, j, axis=1)
        feat_names = [c for i, c in enumerate(cols_) if i != j]
    if M.shape[1] == 0 or M.shape[0] < 3:
        return table([], "need at least three rows and one feature")

    def score(Xa, ya):
        p = _regress_predict(model, Xa, np.asarray(ya, dtype=float), Xa)
        return float(np.mean((p - np.asarray(ya, dtype=float)) ** 2))

    res = ml.permutation_importance(score, M, y, n_repeats=int(repeats), seed=int(seed))
    imp = res.get("importances", []) if isinstance(res, dict) else list(res)
    rows = [{"feature": feat_names[i] if i < len(feat_names) else f"feature_{i}",
            "importance": round(float(dict(d).get("mean_importance", 0.0)), 6),
            "std": round(float(dict(d).get("std", 0.0)), 6),
            "n_repeats": int(dict(d).get("n_repeats", repeats))} for i, d in enumerate(imp)]
    rows.sort(key=lambda r: -r["importance"])
    out = table(rows, f"baseline MSE {res.get('baseline_score') if isinstance(res, dict) else 'n/a'}; "
                    f"{len(rows)} features permuted")
    out["stats"] = {"baseline_score": res.get("baseline_score") if isinstance(res, dict) else None}
    return out

def _edge_list(src):
    """Table -> list of (source, target, weight) edges."""
    df = _frame(src)
    edges = []
    for r in df.to_dict("records"):
        vals = [v for v in r.values() if v not in (None, "")]
        if len(vals) < 2:
            continue
        a, b = str(vals[0]), str(vals[1])
        w = 1.0
        for v in vals[2:]:
            try:
                w = float(v)
                break
            except (TypeError, ValueError):
                continue
        if a and b and a != b:
            edges.append((a, b, w))
    return edges


@T("network_from_edge_list", "Graph metrics from an edge list", GRAPH, "table",
   [tbl("edges", COUNTS, "Edge list (two or three columns)"), boolean("directed", False, "Directed"),
    intin("head", 200, "Rows", min=1)],
   ex={"edges": COUNTS}, up="igraph / networkx", tags=("network", "graph"),
   summary="Degree distribution, components, density and self-loop counts of a graph.")
def network_from_edge_list(edges, directed=False, head=200):
    """Basic graph statistics."""
    el = _edge_list(edges)
    if not el:
        return table([], "no edges found")
    nodes = sorted({x for a, b, _w in el for x in (a, b)})
    deg_out, deg_in, weight = Counter(), Counter(), defaultdict(float)
    adj = defaultdict(set)
    for a, b, w in el:
        deg_out[a] += 1
        deg_in[b] += 1
        weight[a] += w
        weight[b] += w
        adj[a].add(b)
        if not directed:
            adj[b].add(a)
            deg_out[b] += 0
    seen, comps = set(), 0
    sizes = []
    for n in nodes:
        if n in seen:
            continue
        comps += 1
        stack, size = [n], 0
        while stack:
            cur = stack.pop()
            if cur in seen:
                continue
            seen.add(cur)
            size += 1
            stack.extend(adj[cur])
        sizes.append(size)
    m = len(el)
    n = len(nodes)
    dens = m / (n * (n - 1) / 2) if not directed and n > 1 else m / (n * (n - 1)) if n > 1 else 0.0
    rows = [{"node": x, "degree": deg_out[x] + (0 if directed else deg_in[x]),
            "in_degree": deg_in[x], "out_degree": deg_out[x], "strength": round(float(weight[x]), 4),
            "neighbours": len(adj[x])} for x in nodes]
    rows.sort(key=lambda r: -r["degree"])
    out = table(rows[:int(head)], f"{n} nodes, {m} edges, {comps} components, density {round(dens, 5)}")
    out["stats"] = {"components": comps, "component_sizes": sizes, "edges": m, "nodes": n,
                   "isolated": sum(1 for x in nodes if not adj[x])}
    return out


@T("network_communities", "Label-propagation communities", GRAPH, "table",
   [tbl("edges", COUNTS, "Edge list"), intin("iterations", 20, "Sweeps", min=1), intin("head", 300, "Rows", min=1)],
   ex={"edges": COUNTS, "iterations": 15}, up="igraph cluster_label_propagation",
   tags=("network", "communities"),
   summary="Partition a graph into communities by iterative neighbour voting.")
def network_communities(edges, iterations=20, head=300):
    """Community detection by label propagation."""
    el = _edge_list(edges)
    if not el:
        return table([], "no edges found")
    adj = defaultdict(set)
    for a, b, _w in el:
        adj[a].add(b)
        adj[b].add(a)
    nodes = sorted(adj)
    lab = {n: i for i, n in enumerate(nodes)}
    for _ in range(int(iterations)):
        changed = False
        for n in nodes:
            votes = Counter(lab[m] for m in adj[n])
            if votes:
                new = votes.most_common(1)[0][0]
                if new != lab[n]:
                    lab[n] = new
                    changed = True
        if not changed:
            break
    sizes = Counter(lab.values())
    rows = [{"node": n, "community": lab[n], "members": sizes[lab[n]], "degree": len(adj[n])} for n in nodes]
    out = table(rows[:int(head)], f"{len(sizes)} communities over {len(nodes)} nodes "
                                f"(largest {max(sizes.values())})")
    out["stats"] = {"community_sizes": dict(sizes)}
    return out


@T("network_betweenness", "Betweenness centrality of nodes", GRAPH, "table",
   [tbl("edges", COUNTS, "Edge list"), intin("samples", 200, "Random walks", min=10), intin("head", 100, "Rows", min=1)],
   ex={"edges": COUNTS, "samples": 100}, up="igraph betweenness / networkx",
   tags=("network", "centrality"),
   summary="Approximate betweenness by sampling shortest paths between node pairs.")
def network_betweenness(edges, samples=200, head=100):
    """Sampling-based betweenness."""
    el = _edge_list(edges)
    if not el:
        return table([], "no edges found")
    adj = defaultdict(set)
    for a, b, _w in el:
        adj[a].add(b)
        adj[b].add(a)
    nodes = sorted(adj)
    between = Counter()
    rng = np.random.default_rng(1)

    def bfs_path(s, t):
        prev = {s: None}
        q = [s]
        while q:
            nxt = []
            for u in q:
                if u == t:
                    path, cur = [], u
                    while cur is not None:
                        path.append(cur)
                        cur = prev[cur]
                    return list(reversed(path))
                for v in adj[u]:
                    if v not in prev:
                        prev[v] = u
                        nxt.append(v)
            q = nxt
        return []

    n = len(nodes)
    for _ in range(min(int(samples), n * (n - 1))):
        s, t = nodes[int(rng.integers(n))], nodes[int(rng.integers(n))]
        if s == t:
            continue
        p = bfs_path(s, t)
        for mid in p[1:-1]:
            between[mid] += 1
    rows = [{"node": x, "betweenness": between[x], "degree": len(adj[x]),
            "normalised": round(between[x] / max(1, int(samples)), 5)} for x in nodes]
    rows.sort(key=lambda r: -r["betweenness"])
    return table(rows[:int(head)], f"top node: {rows[0]['node']}" if rows else "no paths found")


@T("network_upset_plot", "Intersection sizes of gene or feature sets", GRAPH, "table",
   [txt("sets", GMT, "GMT gene sets"), intin("min_size", 1, "Minimum intersection", min=1),
    boolean("show_plot", True, "Draw a bar plot")],
   ex={"sets": GMT, "min_size": 1}, up="UpSetR / ggupset", tags=("sets", "intersections", "plot"),
   summary="Enumerate every combination of sets with its intersection size (UpSet table).")
def network_upset_plot(sets, min_size=1, show_plot=True):
    """UpSet-style intersections."""
    gsets = {k: set(v) for k, v in io.parse_gmt(io.as_text(sets)).items()}
    names = list(gsets)
    rows = []
    for mask in range(1, 1 << len(names)):
        combo = [names[i] for i in range(len(names)) if mask & (1 << i)]
        inter = set.intersection(*(gsets[c] for c in combo)) if combo else set()
        if len(inter) < int(min_size):
            continue
        rows.append({"sets": "|".join(combo), "n_sets": len(combo), "size": len(inter),
                    "elements": ",".join(sorted(inter)[:8]),
                    "exclusive": len(inter - set.union(*(gsets[c] for c in names if c not in combo)))
                    if len(combo) < len(names) else len(inter)})
    rows.sort(key=lambda r: -r["size"])
    out = table(rows[:200], f"{len(rows)} non-empty intersections over {len(names)} sets")
    if show_plot and rows:
        out["stats"]["figure"] = plot.bar([r["sets"] for r in rows[:12]][::-1],
                                         [r["size"] for r in rows[:12]][::-1],
                                         xlabel="set intersection", ylabel="elements",
                                         title="UpSet intersections", horizontal=True)
    return out


@T("network_venn", "Two-set Venn comparison", GRAPH, "table",
   [tbl("left", COUNTS, "Left table"), tbl("right", COUNTS, "Right table"),
    textbox("key_left", "", "Key column (left)"), textbox("key_right", "", "Key column (right)")],
   ex={"left": COUNTS, "right": COUNTS, "key_left": "gene", "key_right": "gene"},
   up="venneuler / ggVennDiagram", tags=("sets", "venn", "comparison"),
   summary="Shared and unique elements of two key columns with a Jaccard index.")
def network_venn(left, right, key_left="", key_right=""):
    """Venn comparison of two tables."""
    dl, dr = _frame(left), _frame(right)
    kl = key_left or list(dl.columns)[0]
    kr = key_right or list(dr.columns)[0]
    a = {str(v) for v in dl[kl].tolist()} if kl in dl.columns else set()
    b = {str(v) for v in dr[kr].tolist()} if kr in dr.columns else set()
    st = tables.jaccard_index(sorted(a), sorted(b))
    rows = [{"comparison": "only in left", "count": len(a - b), "elements": ",".join(sorted(a - b)[:10])},
            {"comparison": "shared", "count": len(a & b), "elements": ",".join(sorted(a & b)[:10])},
            {"comparison": "only in right", "count": len(b - a), "elements": ",".join(sorted(b - a)[:10])}]
    res = table(rows, f"|A|={len(a)}, |B|={len(b)}, intersection={len(a & b)}, "
                     f"Jaccard={round(float(dict(st).get('jaccard', 0.0)), 4)}")
    res["stats"] = {**{k: (round(float(v), 5) if isinstance(v, float) else v) for k, v in dict(st).items()},
                   "figure": plot.pie([r["comparison"] for r in rows], [r["count"] for r in rows],
                                    title="Set composition")}
    return res


@T("network_flow_table", "Sankey-style flow table between stages", GRAPH, "table",
   [tbl("src", COUNTS, "Table with two categorical columns"), textbox("from_column", "", "Source column"),
    textbox("to_column", "", "Target column"), intin("head", 200, "Rows", min=1)],
   ex={"src": PHENO, "from_column": "group", "to_column": "treatment"},
   up="networkD3 sankey", tags=("network", "flows", "sankey"),
   summary="Aggregate transitions between two categorical columns into weighted flows.")
def network_flow_table(src, from_column="", to_column="", head=200):
    """Flow aggregation."""
    df = _frame(src)
    cols = list(df.columns)
    fc = from_column or cols[0]
    tc = to_column or (cols[1] if len(cols) > 1 else cols[0])
    flows = Counter()
    for r in df.to_dict("records"):
        flows[(str(r.get(fc)), str(r.get(tc)))] += 1
    tot = sum(flows.values()) or 1
    rows = [{"source": a, "target": b, "value": n, "percent": round(100 * n / tot, 3)}
            for (a, b), n in flows.most_common()]
    return table(rows[:int(head)], f"{len(rows)} distinct flows over {tot} records")


@T("interactive_scatter_explorer", "Scatter view with per-point table", INTERACTIVE, "table",
   [tbl("src", COUNTS, "Matrix"), textbox("x_column", "", "x column"), textbox("y_column", "", "y column"),
    textbox("colour_column", "", "Colour column"), textbox("label_column", "", "Label column"),
    intin("head", 200, "Rows", min=1)],
   ex={"src": COUNTS, "x_column": "sample_A", "y_column": "sample_D", "colour_column": "sample_B"},
   up="plotly / interactive scatter", tags=("interactive", "scatter", "tooltip"),
   summary="Scatter of two columns coloured by a third, plus the underlying point table.")
def interactive_scatter_explorer(src, x_column="", y_column="", colour_column="", label_column="", head=200):
    """Scatter with a data table."""
    df = _frame(src)
    numc = [c for c, _v in _num_cols(df)]
    if len(numc) < 2:
        return table([], "need two numeric columns")
    cols = list(df.columns)
    cx = x_column if x_column in numc else numc[0]
    cy = y_column if y_column in numc else numc[1]
    cc = colour_column if colour_column in numc else None
    namec = label_column if label_column in cols else cols[0]
    x = pd.to_numeric(df[cx], errors="coerce").to_numpy(dtype=float)
    y = pd.to_numeric(df[cy], errors="coerce").to_numpy(dtype=float)
    cv = pd.to_numeric(df[cc], errors="coerce").to_numpy(dtype=float) if cc else None
    m = np.isfinite(x) & np.isfinite(y)
    fig = plot.scatter(list(x[m]), list(y[m]), colour_by=list(cv[m]) if cv is not None else None,
                      labels=[str(v) for v in df[namec].tolist()][0:int(m.sum())],
                      xlabel=str(cx), ylabel=str(cy), title=f"{cy} vs {cx}")
    keep = np.where(m)[0][: int(head)]
    rows = [{"point": str(df.iloc[i][namec]), str(cx): round(float(x[i]), 4), str(cy): round(float(y[i]), 4),
            **({str(cc): round(float(cv[i]), 4) if math.isfinite(cv[i]) else "nan"} if cv is not None else {})}
           for i in keep]
    out = table(rows, f"{int(m.sum())} points plotted, {len(rows)} listed")
    out["stats"]["figure"] = fig
    return out

@T("interactive_heatmap_viewer", "Annotated heatmap of a matrix", INTERACTIVE, "figure",
   [tbl("src", COUNTS, "Matrix"), boolean("zscore", True, "Row z-score"),
    choice("cmap", ["viridis", "magma", "coolwarm", "Greys"], "coolwarm", "Colour map"),
    intin("max_rows", 40, "Maximum rows", min=3), intin("max_cols", 40, "Maximum columns", min=2)],
   ex={"src": COUNTS, "zscore": True, "cmap": "viridis"}, up="ComplexHeatmap / plotly heatmap",
   tags=("interactive", "heatmap", "matrix"),
   summary="Display a matrix with value labels and row ordering by mean expression.")
def interactive_heatmap_viewer(src, zscore=True, cmap="coolwarm", max_rows=40, max_cols=40):
    """Matrix viewer with annotation."""
    rows_, cols_, X = _mat(src)
    if X.size == 0:
        return plot.empty_plot("empty matrix")
    M = np.nan_to_num(X.astype(float))
    order = list(np.argsort(-np.nanmean(M, axis=1)))[: int(max_rows)]
    M = M[order][:, : int(max_cols)]
    if zscore:
        mu = np.nanmean(M, axis=1, keepdims=True)
        sd = np.nanstd(M, axis=1, keepdims=True)
        M = (M - mu) / np.where(sd > 0, sd, 1)
    return plot.heatmap([[round(float(v), 3) for v in row] for row in M],
                       row_labels=[rows_[i] for i in order][: max_rows],
                       col_labels=cols_[: max_cols], title=f"matrix ({len(order)} x {M.shape[1]})",
                       cmap=cmap, annot=True)


@T("interactive_track_view", "Coverage track over annotated features", INTERACTIVE, "figure",
   [anyfile("graph", BEDGRAPH, fmt="", label="bedGraph / WIG"), gff("annotation", ANNOT_GFF, "GFF3"),
    intin("bins", 40, "Bins per feature", min=5), choice("mode", ["mean", "max"], "mean", "Aggregation")],
   ex={"graph": BEDGRAPH, "annotation": ANNOT_GFF, "bins": 20}, up="IGV / deepTools plotProfile",
   tags=("interactive", "genome browser", "coverage"),
   summary="Metaprofile of a signal track resampled across every annotated feature.")
def interactive_track_view(graph, annotation, bins=40, mode="mean"):
    """Per-feature resampled profile."""
    vals = io.parse_bedgraph(io.as_text(graph))
    tracks: dict[str, list] = defaultdict(list)
    for chrom, s, e, v in vals:
        tracks[chrom].append((int(s), int(e), float(v)))
    nb = max(5, int(bins))
    matrix = []
    labels = []
    for r in io.parse_gff(io.as_text(annotation)):
        if r.type not in ("gene", "CDS", "mRNA"):
            continue
        prof = np.zeros(nb)
        for s, e, v in tracks.get(r.seqid, []):
            if e <= r.start - 1 or s >= r.end:
                continue
            L = max(1, r.end - r.start)
            for i in range(nb):
                lo, hi = r.start + i * L // nb, r.start + (i + 1) * L // nb
                ov = max(0, min(e, hi) - max(s, lo))
                prof[i] += v * ov / max(1, hi - lo)
        if prof.sum():
            matrix.append(list(prof))
            labels.append(f"{r.seqid}:{r.start}-{r.end}")
    if not matrix:
        return plot.empty_plot("no signal inside the annotated features")
    return plot.profile_matrix(np.array(matrix[:20], dtype=float), labels=labels[:20],
                             title=f"metaprofile ({mode}, {nb} bins)",
                             xlabel="relative position within feature")


@T("interactive_pivot_table", "Pivot a table into a matrix", INTERACTIVE, "table",
   [tbl("src", COUNTS, "Table"), textbox("index_column", "", "Row variable"),
    textbox("value_column", "", "Column to keep as an extra id"),
    choice("agg", ["mean", "sum", "count", "max"], "mean", "Aggregation")],
   ex={"src": COUNTS, "index_column": "gene"}, up="tidyr pivot_wider / reshape2",
   tags=("interactive", "pivot", "table"),
   summary="Melt a wide table and pivot it back with an aggregation function.")
def interactive_pivot_table(src, index_column="", value_column="", agg="mean"):
    """Melt then pivot."""
    df = _frame(src)
    cols = list(df.columns)
    if len(cols) < 2:
        return table([], "need at least two columns")
    idx = index_column if index_column in cols else cols[0]
    numeric = [c for c, _v in _num_cols(df)]
    extra = value_column if value_column in cols and value_column != idx else ""
    keep = [c for c in numeric if c != extra] or numeric
    long = tables.melt(df, id_vars=[idx] + ([extra] if extra else []), value_vars=keep,
                      var_name="variable", value_name="value")
    if long.empty:
        return table([], "nothing left to pivot")
    piv = tables.pivot(long, index=idx, columns="variable", values="value", agg=agg)
    head = [str(c) for c in piv.columns]
    rows = []
    for arr in piv.to_numpy():
        rec = {head[0]: str(arr[0])}
        for k, v in zip(head[1:], arr[1:]):
            ok = isinstance(v, (int, float, np.floating)) and not isinstance(v, bool) and math.isfinite(float(v))
            rec[k] = round(float(v), 5) if ok else v
        rows.append(rec)
    return table(rows, f"pivoted {len(rows)} rows x {len(head) - 1} value columns using {agg}")
