"""Table helpers shared by the text, filter/sort and datamash tool groups."""

from __future__ import annotations

import math
import re
from typing import Any, Sequence

import pandas as pd

from chroma_titan.core.io import as_text, read_table


def load(src: Any, sep: str | None = None, header: Any = "infer",
         names: list[str] | None = None) -> pd.DataFrame:
    return read_table(src, sep=sep, header=header, names=names)


def to_df(rows: Sequence[dict] | Sequence[Sequence] | dict | pd.DataFrame,
          columns: list[str] | None = None) -> pd.DataFrame:
    if isinstance(rows, pd.DataFrame):
        return rows
    if isinstance(rows, dict):
        return pd.DataFrame(rows)
    if rows and isinstance(rows[0], dict):
        return pd.DataFrame(list(rows))
    frame = pd.DataFrame(list(rows))
    if columns:
        if frame.shape[1] == len(columns):
            frame.columns = columns
        else:
            frame = pd.DataFrame(columns=columns)
    return frame


def numeric_only(df: pd.DataFrame) -> pd.DataFrame:
    return df.select_dtypes("number")


def col(df: pd.DataFrame, key: Any, default=None) -> pd.Series:
    """Resolve a column by name or integer index."""
    if isinstance(key, str):
        if key in df.columns:
            return df[key]
        for c in df.columns:
            if str(c).lower() == key.lower():
                return df[c]
        m = re.match(r"^(?:col|c)[_-]?(\d+)$", key, re.I)
        if m:
            i = int(m.group(1)) - 1
            if 0 <= i < df.shape[1]:
                return df.iloc[:, i]
        if default is not None:
            return default
        raise KeyError(f"no column {key!r}; have {list(df.columns)}")
    return df.iloc[:, int(key)]


def values(df: pd.DataFrame, key: Any) -> list[float]:
    s = col(df, key)
    return [float(x) for x in pd.to_numeric(s, errors="coerce").dropna()]


def split_column(df: pd.DataFrame, column: str, sep: str = ",", into: list[str] | None = None,
                 maxsplit: int = -1) -> pd.DataFrame:
    out = df.copy()
    parts = out[column].astype(str).str.split(sep, n=maxsplit, expand=True)
    names = into or [f"{column}_{i + 1}" for i in range(parts.shape[1])]
    for i, nm in enumerate(names):
        if i < parts.shape[1]:
            out[nm] = parts[i]
    return out


def merge_columns(df: pd.DataFrame, columns: list[str], into: str = "merged",
                  sep: str = "_") -> pd.DataFrame:
    out = df.copy()
    out[into] = out[columns].astype(str).agg(sep.join, axis=1)
    return out


def add_column(df: pd.DataFrame, name: str, expression: str) -> pd.DataFrame:
    """Add a column computed from a safe pandas expression (see ``df.eval``)."""
    out = df.copy()
    try:
        out[name] = out.eval(expression, engine="python")
    except Exception as exc:  # noqa: BLE001
        raise ValueError(f"cannot evaluate {expression!r}: {exc}") from None
    return out


def rename_columns(df: pd.DataFrame, mapping: dict[str, str]) -> pd.DataFrame:
    return df.rename(columns=mapping)


def select_columns(df: pd.DataFrame, columns: list[str]) -> pd.DataFrame:
    missing = [c for c in columns if c not in df.columns]
    if missing:
        raise KeyError(f"missing columns: {missing}")
    return df[columns]


def filter_rows(df: pd.DataFrame, expression: str) -> pd.DataFrame:
    try:
        return df.query(expression, engine="python")
    except Exception as exc:  # noqa: BLE001
        raise ValueError(f"bad filter expression {expression!r}: {exc}") from None


def sort_rows(df: pd.DataFrame, by: list[str] | str, ascending: bool | list[bool] = True) -> pd.DataFrame:
    return df.sort_values(by=by, ascending=ascending, kind="stable")


def dedupe(df: pd.DataFrame, subset: list[str] | None = None, keep: str = "first") -> pd.DataFrame:
    return df.drop_duplicates(subset=subset, keep=keep)


def value_counts(df: pd.DataFrame, column: str, normalize: bool = False,
                 top: int = 0) -> pd.DataFrame:
    vc = df[column].value_counts(normalize=normalize)
    if top:
        vc = vc.head(top)
    return vc.reset_index() if hasattr(vc, "reset_index") else pd.DataFrame(
        {"value": vc.index, "count": vc.values})


def groupby_agg(df: pd.DataFrame, by: list[str], funcs: dict[str, str | list[str]]) -> pd.DataFrame:
    agg = df.groupby(by, dropna=False).agg(funcs)
    agg.columns = ["_".join(str(p) for p in c if str(p)) for c in agg.columns.to_flat_index()] \
        if isinstance(agg.columns, pd.MultiIndex) else list(agg.columns)
    return agg.reset_index()


def pivot(df: pd.DataFrame, index: str, columns: str, values: str,
          agg: str = "mean") -> pd.DataFrame:
    return df.pivot_table(index=index, columns=columns, values=values, aggfunc=agg,
                          fill_value=0).reset_index()


def melt(df: pd.DataFrame, id_vars: list[str], value_vars: list[str] | None = None,
         var_name: str = "sample", value_name: str = "value") -> pd.DataFrame:
    return df.melt(id_vars=id_vars, value_vars=value_vars, var_name=var_name,
                   value_name=value_name)


def cpm(mat: pd.DataFrame) -> pd.DataFrame:
    """Counts per million."""
    tot = mat.sum(axis=0).replace(0, math.nan)
    return (mat / tot * 1e6).fillna(0)


def tpm(mat: pd.DataFrame, lengths: pd.Series) -> pd.DataFrame:
    """Transcripts per minute from a counts matrix + effective lengths (bp)."""
    rate = mat.div(lengths.reindex(mat.index).clip(lower=1), axis=0)
    return (rate / rate.sum(axis=0) * 1e6).fillna(0)


def quantile_normalize(mat: pd.DataFrame) -> pd.DataFrame:
    rank_mean = mat.stack().groupby(mat.rank(method="first", axis=0).stack().astype(int)).mean()
    return mat.rank(method="first", axis=0).apply(
        lambda col: col.map(rank_mean)).reindex(mat.index)


def zscore(mat: pd.DataFrame, axis: int = 1) -> pd.DataFrame:
    m = mat.mean(axis=axis).replace(0, math.nan)
    s = mat.std(axis=axis, ddof=1)
    if axis == 1:
        return ((mat.sub(m, axis=0)).div(s.replace(0, math.nan), axis=0)).fillna(0)
    return ((mat.sub(m, axis=1)).div(s.replace(0, math.nan), axis=1)).fillna(0)


def log_transform(mat: pd.DataFrame, base: float = 2, pseudocount: float = 1) -> pd.DataFrame:
    import numpy as np

    return np.log2(mat + pseudocount) if base == 2 else np.log(mat + pseudocount) / math.log(base)


def fold_change(a: pd.Series, b: pd.Series, pseudocount: float = 1.0) -> pd.DataFrame:
    import numpy as np

    fc = (a + pseudocount) / (b + pseudocount)
    return pd.DataFrame({"log2FoldChange": np.log2(fc), "meanAB": (a + b) / 2,
                         "groupA": a.values, "groupB": b.values})


def set_operations(df_a: pd.DataFrame, df_b: pd.DataFrame, column: str,
                   operation: str = "subtract") -> pd.DataFrame:
    sa, sb = set(df_a[column].astype(str)), set(df_b[column].astype(str))
    keep = {"subtract": sa - sb, "intersect": sa & sb, "union": sa | sb,
            "complement": sb - sa}.get(operation, sa - sb)
    return df_a[df_a[column].astype(str).isin(keep)] if operation != "union" else \
        pd.DataFrame({column: sorted(keep)})


def jaccard_index(a: Sequence[str], b: Sequence[str]) -> dict:
    sa, sb = set(a), set(b)
    u = sa | sb
    return {"intersection": len(sa & sb), "union": len(u),
            "jaccard": round(len(sa & sb) / len(u), 6) if u else 0.0}


def window_matrix(df: pd.DataFrame, value_col: str, size: int = 100, step: int | None = None,
                  label_col: str | None = None, fn: str = "mean") -> pd.DataFrame:
    step = step or size
    out = []
    for start in range(0, max(1, len(df) - size + 1), step):
        chunk = df.iloc[start:start + size]
        vals = pd.to_numeric(chunk[value_col], errors="coerce").dropna()
        rec = {"start": start, "end": min(start + size, len(df)),
               "n": int(len(chunk)),
               "value": float(getattr(vals, fn)()) if len(vals) else float("nan")}
        if label_col and label_col in chunk.columns:
            rec["labels"] = "|".join(map(str, chunk[label_col].head(20)))
        out.append(rec)
    return pd.DataFrame(out)


def tidy_long(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    out.columns = [re.sub(r"\s+", "_", str(c).strip()) for c in out.columns]
    for c in out.select_dtypes("number").columns:
        out[c] = pd.to_numeric(out[c], errors="coerce")
    return out.dropna(how="all")
