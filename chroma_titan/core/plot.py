"""Matplotlib figure builders (Agg backend, no display needed).

Tools return the figure object; the web app / CLI turn it into PNG.  Keeping the
plot code here means a tool is just "compute + call a plot builder".
"""

from __future__ import annotations

import math
from typing import Any, Sequence

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

PALETTE = ["#2a9d8f", "#e76f51", "#264653", "#f4a261", "#8ab17d", "#e9c46a",
           "#b56576", "#6d597a", "#355070", "#b5838d"]


def _fig(width: float = 8.5, height: float = 4.2) -> tuple[Any, Any]:
    fig, ax = plt.subplots(figsize=(width, height), dpi=110)
    fig.patch.set_facecolor("white")
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.grid(alpha=0.25, lw=0.6)
    return fig, ax


def _arr(x: Sequence[float]) -> np.ndarray:
    return np.asarray(list(x), dtype=float)


def scatter(x, y, xlabel: str = "x", ylabel: str = "y", title: str = "",
            labels: Sequence[str] | None = None, colour_by: Sequence[float] | None = None,
            size: float = 18) -> Any:
    fig, ax = _fig()
    x, y = _arr(x), _arr(y)
    if colour_by is not None:
        sc = ax.scatter(x, y, c=_arr(colour_by), s=size, cmap="viridis")
        fig.colorbar(sc, ax=ax, shrink=0.8)
    else:
        ax.scatter(x, y, s=size, color=PALETTE[0], alpha=0.8, edgecolor="none")
    if labels:
        for xi, yi, li in list(zip(x, y, labels))[:60]:
            ax.annotate(str(li), (xi, yi), fontsize=6, alpha=0.8)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title or f"{ylabel} vs {xlabel}", fontsize=10)
    return fig


def line(series: dict[str, Sequence[float]], x: Sequence[float] | None = None,
         xlabel: str = "position", ylabel: str = "value", title: str = "",
         legend: bool = True) -> Any:
    fig, ax = _fig()
    xs = _arr(x) if x is not None else np.arange(max(len(v) for v in series.values()))
    for i, (name, ys) in enumerate(series.items()):
        ys = _arr(ys)
        ax.plot(xs[:len(ys)], ys, label=name, lw=1.4, color=PALETTE[i % len(PALETTE)])
    if legend and len(series) > 1:
        ax.legend(fontsize=7, frameon=False)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title, fontsize=10)
    return fig


def bar(categories: Sequence[Any], values: Sequence[float], xlabel: str = "category",
        ylabel: str = "value", title: str = "", horizontal: bool = False,
        colour: str | None = None) -> Any:
    fig, ax = _fig(9, 4.4)
    cats = [str(c) for c in categories]
    vals = _arr(values)
    idx = np.arange(len(cats))
    c = colour or PALETTE[0]
    if horizontal:
        ax.barh(idx, vals, color=c)
        ax.set_yticks(idx)
        ax.set_yticklabels(cats, fontsize=6)
        ax.invert_yaxis()
        ax.set_xlabel(ylabel)
    else:
        ax.bar(idx, vals, color=c)
        ax.set_xticks(idx)
        ax.set_xticklabels(cats, rotation=90 if len(cats) > 8 else 0, fontsize=6)
        ax.set_ylabel(ylabel)
    ax.set_title(title, fontsize=10)
    return fig


def grouped_bar(groups: Sequence[str], series: dict[str, Sequence[float]],
                ylabel: str = "value", title: str = "") -> Any:
    fig, ax = _fig(9, 4.2)
    n, w = len(groups), 0.8 / max(1, len(series))
    for i, (name, vals) in enumerate(series.items()):
        ax.bar(np.arange(n) + i * w, _arr(vals), width=w, label=name,
               color=PALETTE[i % len(PALETTE)])
    ax.set_xticks(np.arange(n) + w * (len(series) - 1) / 2)
    ax.set_xticklabels([str(g) for g in groups], rotation=90, fontsize=6)
    ax.set_ylabel(ylabel)
    ax.legend(fontsize=7, frameon=False)
    ax.set_title(title, fontsize=10)
    return fig


def boxplot(data: dict[str, Sequence[float]], ylabel: str = "value", title: str = "") -> Any:
    fig, ax = _fig(7, 4.2)
    keys = list(data)
    ax.boxplot([_arr(data[k]) for k in keys], labels=[str(k)[:20] for k in keys],
               patch_artist=True,
               boxprops=dict(facecolor="#e9c46a", alpha=0.7), flierprops=dict(ms=2))
    ax.set_ylabel(ylabel)
    ax.set_title(title, fontsize=10)
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right", fontsize=7)
    return fig


def histogram(values: Sequence[float], bins: int = 40, xlabel: str = "value",
              title: str = "", log_y: bool = False, colour: str | None = None) -> Any:
    fig, ax = _fig()
    ax.hist(_arr(values), bins=int(bins), color=colour or PALETTE[0], alpha=0.85)
    if log_y:
        ax.set_yscale("log")
    ax.set_xlabel(xlabel)
    ax.set_ylabel("count")
    ax.set_title(title, fontsize=10)
    return fig


def multi_histogram(series: dict[str, Sequence[float]], bins: int = 40,
                    xlabel: str = "value", title: str = "") -> Any:
    fig, ax = _fig()
    for i, (name, vals) in enumerate(series.items()):
        ax.hist(_arr(vals), bins=bins, alpha=0.55, label=name,
                color=PALETTE[i % len(PALETTE)])
    ax.set_xlabel(xlabel)
    ax.set_ylabel("count")
    ax.legend(fontsize=7, frameon=False)
    ax.set_title(title, fontsize=10)
    return fig


def density(series: dict[str, Sequence[float]], xlabel: str = "value", title: str = "",
            bw: float = 0.25) -> Any:
    fig, ax = _fig()
    for i, (name, vals) in enumerate(series.items()):
        v = _arr(vals)
        if len(v) < 2:
            continue
        lo, hi = v.min() - 3 * bw, v.max() + 3 * bw
        grid = np.linspace(lo, hi, 400)
        kde = np.exp(-0.5 * ((grid[:, None] - v[None, :]) / bw) ** 2).sum(axis=1) / (
            len(v) * bw * math.sqrt(2 * math.pi))
        ax.plot(grid, kde, label=name, lw=1.4, color=PALETTE[i % len(PALETTE)])
    ax.set_xlabel(xlabel)
    ax.set_ylabel("density")
    ax.legend(fontsize=7, frameon=False)
    ax.set_title(title, fontsize=10)
    return fig


def heatmap(matrix, row_labels=None, col_labels=None, title: str = "",
            cmap: str = "viridis", annot: bool = False, figsize=(8, 6)) -> Any:
    fig, ax = _fig(*figsize)
    m = np.asarray(matrix, dtype=float)
    im = ax.imshow(m, aspect="auto", cmap=cmap)
    fig.colorbar(im, ax=ax, shrink=0.8)
    if row_labels is not None:
        n = len(row_labels)
        ax.set_yticks(np.linspace(0, max(1, m.shape[0] - 1), min(n, 25)).astype(int))
        ax.set_yticklabels([str(x) for x in list(row_labels)[:25]], fontsize=6)
    if col_labels is not None:
        n = len(col_labels)
        ax.set_xticks(np.linspace(0, max(1, m.shape[1] - 1), min(n, 25)).astype(int))
        ax.set_xticklabels([str(x) for x in list(col_labels)[:25]], rotation=90, fontsize=6)
    ax.set_title(title, fontsize=10)
    return fig


def volcano(logfc, pvals, labels: Sequence[str] | None = None, title: str = "Volcano plot",
            lfc_cut: float = 1.0, p_cut: float = 0.05) -> Any:
    fig, ax = _fig()
    lfc, p = _arr(logfc), np.clip(_arr(pvals), 1e-300, 1)
    neg = -np.log10(p)
    sig = (np.abs(lfc) > lfc_cut) & (p < p_cut)
    ax.scatter(lfc[~sig], neg[~sig], s=8, color="#adb5bd", alpha=0.7, label="ns")
    up, dn = sig & (lfc > 0), sig & (lfc < 0)
    ax.scatter(lfc[up], neg[up], s=12, color=PALETTE[1], label=f"up ({int(up.sum())})")
    ax.scatter(lfc[dn], neg[dn], s=12, color=PALETTE[2], label=f"down ({int(dn.sum())})")
    ax.axvline(lfc_cut, ls=":", lw=0.8, color="k")
    ax.axvline(-lfc_cut, ls=":", lw=0.8, color="k")
    ax.axhline(-math.log10(p_cut), ls=":", lw=0.8, color="k")
    if labels:
        order = np.argsort(-neg)[:20]
        for i in order:
            ax.annotate(str(labels[i]), (lfc[i], neg[i]), fontsize=6)
    ax.set_xlabel("log2 fold change")
    ax.set_ylabel("-log10 p")
    ax.legend(fontsize=7, frameon=False)
    ax.set_title(title, fontsize=10)
    return fig


def ma_plot(mean_expr, logfc, pvals=None, title: str = "MA plot") -> Any:
    fig, ax = _fig()
    m, l = _arr(mean_expr), _arr(logfc)
    sig = np.ones(len(m), bool) if pvals is None else _arr(pvals) < 0.05
    ax.scatter(m[~sig], l[~sig], s=7, color="#adb5bd", alpha=0.6)
    ax.scatter(m[sig], l[sig], s=10, color=PALETTE[1], alpha=0.9)
    ax.axhline(0, lw=0.8, color="k", ls=":")
    ax.set_xlabel("average expression")
    ax.set_ylabel("log2 fold change")
    ax.set_title(title, fontsize=10)
    return fig


def qq_plot(values, dist: str = "normal", title: str = "Q-Q plot") -> Any:
    fig, ax = _fig(5.5, 4.4)
    v = np.sort(_arr(values))
    n = len(v)
    if n < 2:
        return fig
    p = (np.arange(1, n + 1) - 0.5) / n
    if dist == "uniform":
        theo = p
    else:
        from chroma_titan.core.stats import norm_ppf

        theo = np.array([norm_ppf(float(x)) for x in p])
        if dist == "lognormal":
            v = np.log(v + 1e-12)
    ax.scatter(theo, v, s=9, color=PALETTE[0])
    lo, hi = min(v.min(), theo.min()), max(v.max(), theo.max())
    ax.plot([lo, hi], [lo, hi], ls="--", lw=0.9, color="k")
    ax.set_xlabel("theoretical quantiles")
    ax.set_ylabel("sample quantiles")
    ax.set_title(title, fontsize=10)
    return fig


def manhattan(chrom: Sequence[str], pos: Sequence[float], pvals: Sequence[float],
              title: str = "Manhattan plot") -> Any:
    fig, ax = _fig(11, 4.0)
    ch = [str(c) for c in chrom]
    uniq = sorted(set(ch), key=lambda s: (len(s), s))
    colour = {c: PALETTE[i % 2] for i, c in enumerate(uniq)}
    offsets, acc = {}, 0
    for c in uniq:
        n = sum(1 for x in ch if x == c)
        offsets[c] = acc
        acc += n
    xs = np.array([offsets[c] + i for i, c in enumerate(ch)])
    ys = -np.log10(np.clip(_arr(pvals), 1e-300, 1))
    for i, c in enumerate(ch):
        ax.scatter(xs[i], ys[i], s=7, color=colour[c], alpha=0.85)
    ax.axhline(-math.log10(5e-8), ls="--", lw=0.8, color="crimson")
    ticks = [offsets[c] + sum(1 for x in ch if x == c) / 2 for c in uniq]
    ax.set_xticks(ticks)
    ax.set_xticklabels([str(c)[:8] for c in uniq], fontsize=6, rotation=90)
    ax.set_ylabel("-log10 p")
    ax.set_xlabel("genomic position")
    ax.set_title(title, fontsize=10)
    return fig


def circos(chrom: Sequence[str], start: Sequence[int], end: Sequence[int],
           values: Sequence[float], title: str = "Genomic overview") -> Any:
    fig = plt.figure(figsize=(6.5, 6.5), dpi=110)
    ax = fig.add_subplot(111, projection="polar")
    total = sum(max(1, e - s) for s, e in zip(_arr(start), _arr(end)))
    acc, out = 0.0, []
    for s, e in zip(_arr(start), _arr(end)):
        out.append(acc + (e - s) / 2)
        acc += (e - s)
    theta = np.array(out) / max(1.0, total) * 2 * math.pi
    ax.scatter(theta, np.clip(_arr(values), 0, None), s=6, c=PALETTE[0], alpha=0.7)
    ax.set_ylim(0, max(1.0, float(np.max(values)) if len(values) else 1.0))
    ax.set_title(title, fontsize=10, pad=14)
    return fig


def dendrogram(linkage: Any, labels: Sequence[str], title: str = "Dendrogram",
               orientation: str = "top") -> Any:
    fig, ax = _fig(9, 4)
    _manual_dendrogram(ax, linkage, list(labels), orientation)
    ax.set_title(title, fontsize=10)
    ax.set_ylabel("distance")
    return fig


def _manual_dendrogram(ax, Z, labels, orientation: str = "top") -> None:
    """Draw a dendrogram from a linkage matrix [[a, b, dist, count], ...]."""
    n = len(labels)
    xs: dict[int, float] = {i: float(i) for i in range(n)}
    ys: dict[int, float] = {i: 0.0 for i in range(n)}
    for k, row in enumerate(Z):
        a, b, d = int(row[0]), int(row[1]), float(row[2])
        xa, xb, ya, yb = xs[a], xs[b], ys.get(a, 0.0), ys.get(b, 0.0)
        ax.plot([xa, xa], [ya, d], lw=0.9, color="#264653")
        ax.plot([xb, xb], [yb, d], lw=0.9, color="#264653")
        ax.plot([xa, xb], [d, d], lw=0.9, color="#264653")
        xs[n + k] = (xa + xb) / 2
        ys[n + k] = d
    ax.set_xticks([xs[i] for i in range(n)])
    ax.set_xticklabels([str(l)[:14] for l in labels], rotation=90, fontsize=6)
    ax.grid(alpha=0.2)


def profile_matrix(matrix, labels: Sequence[str] | None = None, title: str = "Metaprofile",
                   xlabel: str = "position relative to feature") -> Any:
    fig, ax = _fig()
    m = np.asarray(matrix, dtype=float)
    if m.ndim == 1:
        m = m[None, :]
    x = np.arange(m.shape[1])
    for i, row in enumerate(m[:40]):
        ax.plot(x, row, lw=1.0, color=PALETTE[i % len(PALETTE)],
                label=(labels[i] if labels and i < len(labels) else None))
    if m.shape[0] > 1:
        ax.plot(x, np.nanmean(m, axis=0), lw=2.2, color="k", label="mean")
    ax.set_xlabel(xlabel)
    ax.set_ylabel("signal")
    if labels:
        ax.legend(fontsize=6, frameon=False, ncol=2)
    ax.set_title(title, fontsize=10)
    return fig


def sequence_logo(counts_matrix, letters: str = "ACGT", title: str = "Sequence logo") -> Any:
    """Classic ggplot2-style sequence logo from a position x letter count matrix."""
    m = np.asarray(counts_matrix, dtype=float)
    fig, ax = _fig(9, 2.8)
    freq = m / np.clip(m.sum(axis=1, keepdims=True), 1e-9, None)
    ic = 2.0 + (freq * np.log2(np.clip(freq, 1e-9, None))).sum(axis=1)
    order = {c: i for i, c in enumerate(letters)}
    for pos in range(m.shape[0]):
        y = 0.0
        stack = sorted(range(m.shape[1]), key=lambda j: freq[pos, j])
        for j in stack:
            h = freq[pos, j] * ic[pos] if pos < len(ic) else 0.0
            if h <= 0:
                continue
            ax.add_patch(plt.Rectangle((pos, y), 0.85, h, color=PALETTE[j % len(PALETTE)]))
            if h > 0.12:
                ax.text(pos + 0.42, y + h / 2, letters[j] if j < len(letters) else "?",
                        ha="center", va="center", fontsize=7, color="white",
                        fontweight="bold")
            y += h
    ax.set_xlim(-0.2, m.shape[0] - 0.8)
    ax.set_ylim(0, max(0.2, float(np.max(ic)) if len(ic) else 1.0))
    ax.set_xticks(np.arange(m.shape[0]))
    ax.set_xticklabels([str(i + 1) for i in range(m.shape[0])], fontsize=6)
    ax.set_ylabel("bits")
    ax.set_title(title, fontsize=10)
    ax.grid(alpha=0.15)
    return fig


def dotplot(seq1: str, seq2: str | None = None, k: int = 5, title: str = "Dot plot") -> Any:
    from chroma_titan.core.seq import clean, revcomp

    a = clean(seq1)
    b = clean(seq2) if seq2 else a
    ka = [a[i:i + k] for i in range(len(a) - k + 1)]
    kb = [b[i:i + k] for i in range(len(b) - k + 1)]
    index: dict[str, list[int]] = {}
    for i, km in enumerate(kb):
        index.setdefault(km, []).append(i)
    xs, ys = [], []
    for i, km in enumerate(ka):
        for j in index.get(km, ()):
            xs.append(j)
            ys.append(i)
    fig, ax = _fig(6.5, 5.5)
    ax.scatter(xs, ys, s=1.5, color=PALETTE[0])
    if not seq2:
        ax.set_xlabel("self (reverse complement)")
        ax.set_ylabel("forward")
    else:
        ax.set_xlabel("sequence 2")
        ax.set_ylabel("sequence 1")
    ax.invert_yaxis()
    ax.set_title(title, fontsize=10)
    del revcomp
    return fig


def coverage_tracks(regions: dict[str, Sequence[float]], bins: Sequence[tuple[str, int, int]] | None = None,
                    title: str = "Coverage") -> Any:
    fig, ax = _fig(10, 3.4)
    for i, (name, vals) in enumerate(regions.items()):
        v = _arr(vals)
        ax.fill_between(np.arange(len(v)), v, alpha=0.55, label=name,
                        color=PALETTE[i % len(PALETTE)])
    ax.set_xlabel("bin")
    ax.set_ylabel("coverage")
    ax.legend(fontsize=7, frameon=False)
    ax.set_title(title, fontsize=10)
    return fig


def scatter_matrix(df, columns: Sequence[str], title: str = "Scatter matrix") -> Any:
    cols = list(columns)[:5]
    n = len(cols)
    fig, axes = plt.subplots(n, n, figsize=(2.3 * n, 2.3 * n), squeeze=False)
    for i, ci in enumerate(cols):
        for j, cj in enumerate(cols):
            ax = axes[i][j]
            x = pd_num(df[cj])
            y = pd_num(df[ci])
            if i == j:
                ax.hist(x.dropna(), bins=15, color=PALETTE[0])
            else:
                ax.scatter(x, y, s=5, color=PALETTE[2], alpha=0.7)
            ax.set_xticks([])
            ax.set_yticks([])
            if i == n - 1:
                ax.set_xlabel(str(cj)[:10], fontsize=6, rotation=90)
            if j == 0:
                ax.set_ylabel(str(ci)[:10], fontsize=6)
    fig.suptitle(title, fontsize=10)
    fig.tight_layout()
    return fig


def pd_num(s) -> np.ndarray:
    import pandas as pd

    return pd.to_numeric(pd.Series(s), errors="coerce")


def pie(labels: Sequence[Any], values: Sequence[float], title: str = "") -> Any:
    fig, ax = _fig(5.2, 4.4)
    vals = _arr(values)
    keep = vals > 0
    ax.pie(vals[keep], labels=[str(l) for l, k in zip(labels, keep) if k],
           autopct="%1.0f%%", textprops=dict(fontsize=6),
           colors=[PALETTE[i % len(PALETTE)] for i in range(int(keep.sum()))])
    ax.set_title(title, fontsize=10)
    return fig


def violin(groups: dict[str, Sequence[float]], title: str = "Distribution by group") -> Any:
    fig, ax = _fig(7.5, 4.4)
    keys = list(groups)
    data = [_arr(groups[k]) for k in keys]
    parts = ax.violinplot(data, showmeans=True, showextrema=False)
    for i, pc in enumerate(parts["bodies"]):
        pc.set_facecolor(PALETTE[i % len(PALETTE)])
        pc.set_alpha(0.6)
    ax.set_xticks(np.arange(1, len(keys) + 1))
    ax.set_xticklabels([str(k)[:14] for k in keys], rotation=45, ha="right", fontsize=7)
    ax.set_title(title, fontsize=10)
    return fig


def heatmap_with_tree(matrix, row_labels=None, col_labels=None, title: str = "") -> Any:
    return heatmap(matrix, row_labels, col_labels, title=title, cmap="RdYlBu_r")


def image_show(arr, title: str = "", cmap: str = "gray") -> Any:
    fig, ax = _fig(5.5, 5.0)
    im = ax.imshow(np.asarray(arr), cmap=cmap)
    fig.colorbar(im, ax=ax, shrink=0.85)
    ax.set_title(title, fontsize=10)
    ax.set_xticks([])
    ax.set_yticks([])
    return fig


def bar_error(categories, values, errors, ylabel: str = "value", title: str = "") -> Any:
    fig, ax = _fig()
    ax.bar(np.arange(len(list(categories))), _arr(values), yerr=_arr(errors),
           capsize=3, color=PALETTE[0], alpha=0.85)
    ax.set_xticks(np.arange(len(list(categories))))
    ax.set_xticklabels([str(c) for c in categories], rotation=90, fontsize=6)
    ax.set_ylabel(ylabel)
    ax.set_title(title, fontsize=10)
    return fig


def scatter_3d(x, y, z, labels: Sequence[str] | None = None, title: str = "") -> Any:
    fig = plt.figure(figsize=(6.5, 5), dpi=110)
    ax = fig.add_subplot(111, projection="3d")
    ax.scatter(_arr(x), _arr(y), _arr(z), s=12, c=PALETTE[0], alpha=0.8)
    ax.set_xlabel("PC1")
    ax.set_ylabel("PC2")
    ax.set_zlabel("PC3")
    ax.set_title(title, fontsize=10)
    return fig


def empty_plot(title: str = "no data") -> Any:
    fig, ax = _fig()
    ax.text(0.5, 0.5, title, ha="center", va="center", fontsize=12, color="#888")
    ax.set_xticks([])
    ax.set_yticks([])
    return fig
