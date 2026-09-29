"""Statistics: descriptive tests and p-values without a SciPy dependency.

Implemented with standard series/continued-fraction algorithms so the suite runs
on numpy + the standard library only.  Every function returns plain python
numbers so tool results stay JSON friendly.
"""

from __future__ import annotations

import math
import random
from collections import Counter
from typing import Sequence

Number = float


def _as_list(x) -> list[float]:
    if x is None:
        return []
    if hasattr(x, "tolist"):
        x = x.tolist()
    if isinstance(x, (str, bytes)):
        x = [float(v) for v in str(x).replace(",", " ").split() if v]
    out = []
    for v in list(x):
        try:
            f = float(v)
        except (TypeError, ValueError):
            continue
        if not math.isnan(f):
            out.append(f)
    return out


# ---------------------------------------------------------------------------
# incomplete gamma / beta (Numerical Recipes style)
# ---------------------------------------------------------------------------
def _gser(a: float, x: float) -> float:
    """Regularised lower incomplete gamma P(a,x)."""
    ap, total, term = a, 1.0 / a, 1.0 / a
    for _ in range(1000):
        ap += 1
        term *= x / ap
        total += term
        if abs(term) < abs(total) * 1e-15:
            break
    return total * math.exp(-x + a * math.log(x) - math.lgamma(a))


def _gcf(a: float, x: float) -> float:
    """Continued fraction for the upper incomplete gamma Q(a,x)."""
    tiny = 1e-300
    b = x + 1.0 - a
    c = 1.0 / tiny
    d = h = 1.0 / b
    for i in range(1, 1000):
        an = -i * (i - a)
        b += 2.0
        d = an * d + b
        if abs(d) < tiny:
            d = tiny
        c = b + an / c
        if abs(c) < tiny:
            c = tiny
        d = 1.0 / d
        delta = d * c
        h *= delta
        if abs(delta - 1.0) < 1e-15:
            break
    return math.exp(-x + a * math.log(x) - math.lgamma(a)) * h


def _gammq(a: float, x: float) -> float:
    """Upper incomplete gamma Q(a,x) = 1 - P(a,x)."""
    if x < 0 or a <= 0:
        return float("nan")
    if x == 0:
        return 1.0
    val = 1.0 - _gser(a, x) if x < a + 1 else _gcf(a, x)
    return max(0.0, min(1.0, val))


def gamma_p(a: float, x: float) -> float:
    return max(0.0, min(1.0, _gser(a, x)))


def _betacf(a: float, b: float, x: float) -> float:
    qab, qap, qam = a + b, a + 1.0, a - 1.0
    c, d = 1.0, 1.0 - qab * x / qap
    if abs(d) < 1e-300:
        d = 1e-300
    d = 1.0 / d
    h = d
    for m in range(1, 300):
        m2 = 2 * m
        aa = m * (b - m) * x / ((qam + m2) * (a + m2))
        d = 1.0 + aa * d
        if abs(d) < 1e-300:
            d = 1e-300
        c = 1.0 + aa / c
        if abs(c) < 1e-300:
            c = 1e-300
        d = 1.0 / d
        h *= d * c
        aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2))
        d = 1.0 + aa * d
        if abs(d) < 1e-300:
            d = 1e-300
        c = 1.0 + aa / c
        if abs(c) < 1e-300:
            c = 1e-300
        d = 1.0 / d
        de = d * c
        h *= de
        if abs(de - 1.0) < 1e-14:
            break
    return h


def betainc(a: float, b: float, x: float) -> float:
    """Regularised incomplete beta I_x(a,b)."""
    if x <= 0:
        return 0.0
    if x >= 1:
        return 1.0
    lbeta = math.lgamma(a + b) - math.lgamma(a) - math.lgamma(b) + a * math.log(x) + b * math.log(1 - x)
    bt = math.exp(lbeta)
    if x < (a + 1) / (a + b + 2):
        return max(0.0, min(1.0, bt * _betacf(a, b, x) / a))
    return max(0.0, min(1.0, 1 - bt * _betacf(b, a, 1 - x) / b))


# ---------------------------------------------------------------------------
# distributions
# ---------------------------------------------------------------------------
def norm_cdf(z: float) -> float:
    return 0.5 * (1 + math.erf(z / math.sqrt(2)))


def norm_sf(z: float) -> float:
    return 1 - norm_cdf(z)


def norm_ppf(p: float) -> float:
    """Inverse normal CDF (Acklam's rational approximation)."""
    if not 0 < p < 1:
        return float("nan")
    a = [-3.969683028665376e+01, 2.209460984245205e+02, -2.759285104469687e+02,
         1.383577518672690e+02, -3.066479806614716e+01, 2.506628277459239e+00]
    b = [-5.447609879822406e+01, 1.615858368580409e+02, -1.556989798598866e+02,
         6.680131188771972e+01, -1.328068155043206e+01]
    c = [-7.784894002430293e-03, -3.223964580411365e-01, -2.400758277161838e+00,
         -2.549732539343734e+00, 4.374664141464968e+00, 2.938163982698783e+00]
    d = [7.784695709041462e-03, 3.224671290700398e-01, 2.445134137142996e+00,
         3.754408661907416e+00]
    plow, phigh = 0.02425, 1 - 0.02425
    if p < plow:
        q = math.sqrt(-2 * math.log(p))
        return (((((c[0] * q + c[1]) * q + c[2]) * q + c[3]) * q + c[4]) * q + c[5]) / \
               ((((d[0] * q + d[1]) * q + d[2]) * q + d[3]) * q + 1)
    if p > phigh:
        q = math.sqrt(-2 * math.log(1 - p))
        return -(((((c[0] * q + c[1]) * q + c[2]) * q + c[3]) * q + c[4]) * q + c[5]) / \
               ((((d[0] * q + d[1]) * q + d[2]) * q + d[3]) * q + 1)
    q = p - 0.5
    r = q * q
    return (((((a[0] * r + a[1]) * r + a[2]) * r + a[3]) * r + a[4]) * r + a[5]) * q / \
           (((((b[0] * r + b[1]) * r + b[2]) * r + b[3]) * r + b[4]) * r + 1)


def t_pvalue(t: float, df: float) -> float:
    if df <= 0 or math.isnan(t):
        return float("nan")
    return betainc(df / 2, 0.5, df / (df + t * t))


def chi2_pvalue(x: float, dof: float) -> float:
    if dof <= 0 or x < 0:
        return float("nan")
    return _gammq(dof / 2, x / 2)


def f_pvalue(f: float, df1: float, df2: float) -> float:
    if f <= 0:
        return 1.0
    return betainc(df2 / 2, df1 / 2, df2 / (df2 + df1 * f))


def poisson_pmf(k: int, lam: float) -> float:
    if lam <= 0:
        return 1.0 if k == 0 else 0.0
    return math.exp(-lam + k * math.log(lam) - math.lgamma(k + 1))


def poisson_sf(k: int, lam: float) -> float:
    """P(X >= k) for X~Poisson(lam)."""
    if k <= 0:
        return 1.0
    return min(1.0, sum(poisson_pmf(i, lam) for i in range(k, min(k + 5000,
                                                                  int(lam + 40 * math.sqrt(lam + 1)) + 20))))


def poisson_cdf(k: int, lam: float) -> float:
    if k < 0:
        return 0.0
    return min(1.0, sum(poisson_pmf(i, lam) for i in range(0, min(k + 1, int(lam + 40 * math.sqrt(lam + 1)) + 20))))


def binom_pmf(k: int, n: int, p: float) -> float:
    if not 0 <= k <= n:
        return 0.0
    p = min(max(p, 1e-300), 1 - 1e-300)
    return math.exp(math.lgamma(n + 1) - math.lgamma(k + 1) - math.lgamma(n - k + 1)
                    + k * math.log(p) + (n - k) * math.log(1 - p))


def binom_sf(k: int, n: int, p: float) -> float:
    if k <= 0:
        return 1.0
    return min(1.0, sum(binom_pmf(i, n, p) for i in range(k, n + 1)))


def binom_two_sided(k: int, n: int, p: float = 0.5) -> float:
    pk = binom_pmf(k, n, p)
    return min(1.0, sum(binom_pmf(i, n, p) for i in range(n + 1) if binom_pmf(i, n, p) <= pk * 1.0001))


def fisher_exact(a: int, b: int, c: int, d: int) -> tuple[float, float]:
    """Odds ratio + two-sided p (sum of tables as/less likely)."""
    n = a + b + c + d
    r1, r2, c1 = a + b, c + d, a + c

    def p_at(x: int) -> float:
        # hypergeometric pmf: C(r1,x) C(r2, c1-x) / C(n, c1)
        return math.exp(math.lgamma(r1 + 1) + math.lgamma(r2 + 1) + math.lgamma(c1 + 1)
                        + math.lgamma(n - c1 + 1) - math.lgamma(x + 1)
                        - math.lgamma(r1 - x + 1) - math.lgamma(c1 - x + 1)
                        - math.lgamma(r2 - c1 + x + 1) - math.lgamma(n + 1))

    lo, hi = max(0, c1 - r2), min(r1, c1)
    if n <= 0 or hi < lo:
        return float("nan"), 1.0
    a = max(lo, min(hi, a))
    p0 = p_at(a)
    total = sum(p_at(x) for x in range(lo, hi + 1) if p_at(x) <= p0 * 1.0000001)
    orr = (a * d) / (b * c) if b and c else float("inf") if a and d else 0.0
    return orr, min(1.0, total)


# ---------------------------------------------------------------------------
# descriptive
# ---------------------------------------------------------------------------
def mean(x) -> float:
    v = _as_list(x)
    return sum(v) / len(v) if v else float("nan")


def median(x) -> float:
    v = sorted(_as_list(x))
    n = len(v)
    if not n:
        return float("nan")
    return v[n // 2] if n % 2 else (v[n // 2 - 1] + v[n // 2]) / 2


def geomean(x) -> float:
    v = [a for a in _as_list(x) if a > 0]
    return math.exp(sum(math.log(a) for a in v) / len(v)) if v else 0.0


def harmonic_mean(x) -> float:
    v = [a for a in _as_list(x) if a > 0]
    return len(v) / sum(1 / a for a in v) if v else 0.0


def variance(x, sample: bool = True) -> float:
    v = _as_list(x)
    if len(v) < 2:
        return 0.0
    m = mean(v)
    return sum((a - m) ** 2 for a in v) / (len(v) - (1 if sample else 0))


def stdev(x, sample: bool = True) -> float:
    return math.sqrt(variance(x, sample))


def sem(x) -> float:
    v = _as_list(x)
    return stdev(v) / math.sqrt(len(v)) if len(v) > 1 else 0.0


def cv(x) -> float:
    m = mean(x)
    return 100 * stdev(x) / m if m else float("nan")


def mad(x) -> float:
    v = _as_list(x)
    if not v:
        return 0.0
    return median([abs(a - median(v)) for a in v])


def quantile(x, q: float) -> float:
    v = sorted(_as_list(x))
    if not v:
        return float("nan")
    idx = q * (len(v) - 1)
    lo, hi = int(math.floor(idx)), int(math.ceil(idx))
    return v[lo] + (v[hi] - v[lo]) * (idx - lo)


def iqr(x) -> float:
    return quantile(x, 0.75) - quantile(x, 0.25)


def percentile_of(x, value: float) -> float:
    v = sorted(_as_list(x))
    if not v:
        return float("nan")
    return 100 * sum(1 for a in v if a <= value) / len(v)


def skewness(x) -> float:
    v = _as_list(x)
    if len(v) < 3:
        return 0.0
    m, s = mean(v), stdev(v)
    return sum(((a - m) / s) ** 3 for a in v) / len(v) if s else 0.0


def kurtosis(x, excess: bool = True) -> float:
    v = _as_list(x)
    if len(v) < 4:
        return 0.0
    m, s = mean(v), stdev(v)
    k = sum(((a - m) / s) ** 4 for a in v) / len(v) if s else 0.0
    return k - 3 if excess else k


def summary(x, label: str = "") -> dict:
    v = _as_list(x)
    n = len(v)
    if not n:
        return {"label": label, "count": 0}
    return {"label": label, "count": n, "mean": mean(v), "median": median(v),
            "min": min(v), "max": max(v), "std": stdev(v), "variance": variance(v),
            "sem": sem(v), "cv_percent": cv(v), "IQR": iqr(v), "MAD": mad(v),
            "q25": quantile(v, 0.25), "q75": quantile(v, 0.75), "sum": sum(v),
            "geomean": geomean(v), "skew": skewness(v), "kurtosis": kurtosis(v),
            "unique": len(set(v)), "missing": 0,
            "shapiro_jb_p": jarque_bera_p(v)}


def describe(values: dict[str, Sequence]) -> list[dict]:
    return [summary(v, k) for k, v in values.items()]


# ---------------------------------------------------------------------------
# tests
# ---------------------------------------------------------------------------
def ttest_1samp(x, mu: float = 0.0) -> dict:
    v = _as_list(x)
    if len(v) < 2:
        return {"statistic": float("nan"), "p_value": float("nan"), "df": 0}
    m, s = mean(v), stdev(v)
    t = (m - mu) / (s / math.sqrt(len(v))) if s else float("inf")
    return {"statistic": t, "p_value": t_pvalue(t, len(v) - 1), "df": len(v) - 1,
            "mean": m, "difference": m - mu}


def ttest_ind(a, b, equal_var: bool = True) -> dict:
    x, y = _as_list(a), _as_list(b)
    if len(x) < 2 or len(y) < 2:
        return {"statistic": float("nan"), "p_value": float("nan")}
    n1, n2 = len(x), len(y)
    m1, m2 = mean(x), mean(y)
    if equal_var:
        v1, v2 = variance(x), variance(y)
        sp = math.sqrt(((n1 - 1) * v1 + (n2 - 1) * v2) / (n1 + n2 - 2))
        t = (m1 - m2) / (sp * math.sqrt(1 / n1 + 1 / n2)) if sp else float("nan")
        df = n1 + n2 - 2
    else:
        v1, v2 = variance(x) / n1, variance(y) / n2
        den = math.sqrt(v1 + v2)
        t = (m1 - m2) / den if den else float("nan")
        df = (v1 + v2) ** 2 / (v1 ** 2 / (n1 - 1) + v2 ** 2 / (n2 - 1)) if v1 and v2 else n1 + n2 - 2
    return {"statistic": t, "p_value": t_pvalue(t, df), "df": df,
            "mean_difference": m1 - m2, "n1": n1, "n2": n2,
            "cohens_d": (m1 - m2) / math.sqrt((variance(x) + variance(y)) / 2)
            if (variance(x) or variance(y)) else 0.0}


def paired_ttest(a, b) -> dict:
    x, y = _as_list(a), _as_list(b)
    d = [p - q for p, q in zip(x, y)]
    res = ttest_1samp(d, 0.0)
    res["n_pairs"] = len(d)
    return res


def mann_whitney_u(a, b) -> dict:
    x, y = _as_list(a), _as_list(b)
    if not x or not y:
        return {"U": 0.0, "p_value": float("nan")}
    allv = sorted(x + y)
    ranks = _ranks(allv)
    rank_of = {}
    for val, r in zip(allv, ranks):
        rank_of.setdefault(val, r)
    # average ranks handled via _ranks on combined list
    combined = x + y
    rk = _ranks(combined)
    rx = sum(rk[:len(x)])
    n1, n2 = len(x), len(y)
    u1 = rx - n1 * (n1 + 1) / 2
    u2 = n1 * n2 - u1
    u = min(u1, u2)
    mu = n1 * n2 / 2
    tie = sum(t ** 3 - t for t in Counter(combined).values())
    sigma = math.sqrt((n1 * n2 / 12) * ((n1 + n1 + n2 + 1) - 0)) if False else \
        math.sqrt(((n1 * n2) / 12) * ((n1 + n2 + 1) - tie / (n1 + n2) / (n1 + n2 - 1))) \
        if n1 + n2 > 1 else 1.0
    z = (u - mu) / sigma if sigma else 0.0
    return {"U": u, "U1": u1, "U2": u2, "z": z,
            "p_value": max(0.0, min(1.0, 2 * norm_sf(abs(z)))),
            "n1": n1, "n2": n2, "median1": median(x), "median2": median(y)}


def _ranks(vals: Sequence[float]) -> list[float]:
    order = sorted(range(len(vals)), key=lambda i: vals[i])
    rk = [0.0] * len(vals)
    i = 0
    while i < len(order):
        j = i
        while j + 1 < len(order) and vals[order[j + 1]] == vals[order[i]]:
            j += 1
        avg = (i + j) / 2 + 1
        for k in range(i, j + 1):
            rk[order[k]] = avg
        i = j + 1
    return rk


def wilcoxon_signed_rank(a, b=None) -> dict:
    x = _as_list(a)
    y = _as_list(b) if b is not None else [0.0] * len(x)
    d = [p - q for p, q in zip(x, y) if p != q]
    if not d:
        return {"W": 0.0, "p_value": 1.0, "n": 0}
    rk = _ranks([abs(v) for v in d])
    pos = sum(r for r, v in zip(rk, d) if v > 0)
    neg = sum(r for r, v in zip(rk, d) if v < 0)
    w = min(pos, neg)
    n = len(d)
    mu = n * (n + 1) / 4
    sigma = math.sqrt(n * (n + 1) * (2 * n + 1) / 24)
    z = (w - mu) / sigma if sigma else 0.0
    return {"W": w, "z": z, "p_value": max(0.0, min(1.0, 2 * norm_sf(abs(z)))), "n": n}


def kstest_1samp(x, cdf=norm_cdf, mu: float = 0.0, sd: float = 1.0) -> dict:
    v = sorted(_as_list(x))
    if not v:
        return {"D": 0.0, "p_value": 1.0}
    n = len(v)
    d = max(max(abs(norm_cdf((x0 - mu) / sd) - (i + 1) / n),
                abs(norm_cdf((x0 - mu) / sd) - i / n)) for i, x0 in enumerate(v))
    return {"D": d, "p_value": kolmogorov_p(d, n), "n": n}


def kstest_2samp(a, b) -> dict:
    x, y = sorted(_as_list(a)), sorted(_as_list(b))
    if not x or not y:
        return {"D": 0.0, "p_value": 1.0}
    allv = sorted(set(x + y))
    d = 0.0
    for t in allv:
        ex = sum(1 for v in x if v <= t) / len(x)
        ey = sum(1 for v in y if v <= t) / len(y)
        d = max(d, abs(ex - ey))
    ne = len(x) * len(y) / (len(x) + len(y))
    return {"D": d, "p_value": kolmogorov_p(d, ne), "n1": len(x), "n2": len(y)}


def kolmogorov_p(d: float, n: float) -> float:
    if n <= 0 or d <= 0:
        return 1.0
    lam = (math.sqrt(n) + 0.12 + 0.11 / math.sqrt(n)) * d
    s = 0.0
    for j in range(1, 101):
        term = 2 * (-1) ** (j - 1) * math.exp(-2 * j * j * lam * lam)
        s += term
        if abs(term) < 1e-12:
            break
    return max(0.0, min(1.0, s))


def chi_square(observed, expected=None) -> dict:
    o = _as_list(observed)
    if expected is None:
        e = [sum(o) / len(o)] * len(o)
    else:
        e = _as_list(expected)
    stat = sum((a - b) ** 2 / b for a, b in zip(o, e) if b)
    dof = max(1, len(o) - 1)
    return {"statistic": stat, "p_value": chi2_pvalue(stat, dof), "df": dof}


def contingency_chi2(table: list[list[int]]) -> dict:
    rows = len(table)
    cols = len(table[0])
    n = sum(sum(r) for r in table)
    rt = [sum(r) for r in table]
    ct = [sum(table[i][j] for i in range(rows)) for j in range(cols)]
    stat = 0.0
    for i in range(rows):
        for j in range(cols):
            e = rt[i] * ct[j] / n if n else 0
            if e:
                stat += (table[i][j] - e) ** 2 / e
    dof = max(1, (rows - 1) * (cols - 1))
    exp = [[round(rt[i] * ct[j] / n, 3) if n else 0 for j in range(cols)] for i in range(rows)]
    return {"statistic": stat, "p_value": chi2_pvalue(stat, dof), "df": dof,
            "expected": exp, "cramers_v": math.sqrt(stat / (n * min(rows - 1, cols - 1)))
            if n and min(rows, cols) > 1 else 0.0}


def anova(*groups) -> dict:
    gs = [_as_list(g) for g in groups if _as_list(g)]
    if len(gs) < 2:
        return {"F": float("nan"), "p_value": float("nan")}
    allv = [v for g in gs for v in g]
    gm = mean(allv)
    k = len(gs)
    n = len(allv)
    ssb = sum(len(g) * (mean(g) - gm) ** 2 for g in gs)
    ssw = sum(sum((v - mean(g)) ** 2 for v in g) for g in gs)
    df1, df2 = k - 1, max(1, n - k)
    f = (ssb / df1) / (ssw / df2) if ssw else float("inf")
    return {"F": f, "p_value": f_pvalue(f, df1, df2), "df_between": df1, "df_within": df2,
            "eta_squared": ssb / (ssb + ssw) if (ssb + ssw) else 0.0, "groups": k, "n": n}


def kruskal(*groups) -> dict:
    gs = [_as_list(g) for g in groups if _as_list(g)]
    allv = [v for g in gs for v in g]
    if len(gs) < 2 or len(allv) < 3:
        return {"H": 0.0, "p_value": 1.0}
    rk = _ranks(allv)
    i = 0
    h = 0.0
    for g in gs:
        s = sum(rk[i:i + len(g)])
        h += s * s / len(g)
        i += len(g)
    n = len(allv)
    h = 12 / (n * (n + 1)) * h - 3 * (n + 1)
    tie = sum(t ** 3 - t for t in Counter(allv).values())
    if tie:
        h /= 1 - tie / (n ** 3 - n)
    dof = max(1, len(gs) - 1)
    return {"H": h, "p_value": chi2_pvalue(h, dof), "df": dof, "groups": len(gs), "n": n}


def pearson(a, b) -> dict:
    x, y = _as_list(a), _as_list(b)
    n = min(len(x), len(y))
    x, y = x[:n], y[:n]
    if n < 3:
        return {"r": float("nan"), "p_value": float("nan"), "n": n}
    mx, my = mean(x), mean(y)
    sxy = sum((p - mx) * (q - my) for p, q in zip(x, y))
    sx = math.sqrt(sum((p - mx) ** 2 for p in x))
    sy = math.sqrt(sum((q - my) ** 2 for q in y))
    r = sxy / (sx * sy) if sx and sy else 0.0
    r = max(-0.999999, min(0.999999, r))
    t = r * math.sqrt((n - 2) / (1 - r * r))
    return {"r": r, "r_squared": r * r, "t": t, "p_value": t_pvalue(t, n - 2), "n": n,
            "covariance": sxy / (n - 1)}


def spearman(a, b) -> dict:
    x, y = _as_list(a), _as_list(b)
    n = min(len(x), len(y))
    rx, ry = _ranks(x[:n]), _ranks(y[:n])
    res = pearson(rx, ry)
    res["method"] = "spearman"
    return res


def kendall(a, b) -> dict:
    x, y = _as_list(a), _as_list(b)
    n = min(len(x), len(y))
    conc = disc = tiesx = tiesy = 0
    for i in range(n):
        for j in range(i + 1, n):
            dx = (x[j] - x[i])
            dy = (y[j] - y[i])
            if dx == 0 and dy == 0:
                continue
            if dx == 0:
                tiesx += 1
                continue
            if dy == 0:
                tiesy += 1
                continue
            if dx * dy > 0:
                conc += 1
            else:
                disc += 1
    n0 = n * (n - 1) / 2
    denom = math.sqrt((n0 - tiesx) * (n0 - tiesy))
    tau = (conc - disc) / denom if denom else 0.0
    var = 2 * (2 * n + 5) / (9 * n * (n - 1)) if n > 1 else 1
    z = tau / math.sqrt(var) if var else 0.0
    return {"tau": tau, "z": z, "p_value": max(0.0, min(1.0, 2 * norm_sf(abs(z)))),
            "concordant": conc, "discordant": disc, "n": n}


def correlation_matrix(data: dict[str, Sequence], method: str = "pearson") -> list[dict]:
    keys = list(data)
    rows = []
    for i, ka in enumerate(keys):
        for kb in keys[i:]:
            res = spearman(data[ka], data[kb]) if method == "spearman" else pearson(
                data[ka], data[kb])
            rows.append({"variable_1": ka, "variable_2": kb, "coefficient": round(
                res.get("r", res.get("tau", 0.0)), 6), "p_value": res.get("p_value"),
                "n": res.get("n", 0)})
    return rows


def linear_regression(x, y, order: int = 1) -> dict:
    xs, ys = _as_list(x), _as_list(y)
    n = min(len(xs), len(ys))
    if n < order + 2:
        return {"slope": float("nan"), "intercept": float("nan"), "r_squared": 0.0}
    xs, ys = xs[:n], ys[:n]
    import numpy as np

    X = np.vander(np.array(xs, dtype=float), order + 1)
    coef, res, rank, sv = np.linalg.lstsq(X, np.array(ys, dtype=float), rcond=None)
    pred = X @ coef
    ss_res = float(((np.array(ys) - pred) ** 2).sum())
    ss_tot = float(((np.array(ys) - mean(ys)) ** 2).sum())
    r2 = 1 - ss_res / ss_tot if ss_tot else 1.0
    se = math.sqrt(ss_res / max(1, n - order - 1))
    # np.vander puts the highest power first: coef == [a_order, ..., a1, a0]
    slope = float(coef[-2]) if order >= 1 else float(coef[0])
    try:
        xtx_inv = np.linalg.pinv(X.T @ X)
        se_slope = se * math.sqrt(max(1e-300, float(xtx_inv[len(coef) - 2, len(coef) - 2])))
    except Exception:  # noqa: BLE001
        se_slope = se or 1e-12
    tstat = slope / (se_slope or 1e-12)
    return {"coefficients": [float(c) for c in coef[::-1]], "intercept": float(coef[-1]),
            "slope": slope, "r_squared": r2, "residual_stderr": se,
            "p_value_slope": t_pvalue(tstat, n - order - 1), "n": n,
            "pearson_r": pearson(xs, ys)["r"] if order == 1 else float("nan")}


def loess_fit(x, y, span: float = 0.5) -> list[dict]:
    xs, ys = _as_list(x), _as_list(y)
    n = min(len(xs), len(ys))
    if n < 3:
        return []
    xs, ys = xs[:n], ys[:n]
    order = max(2, int(span * n))
    out = []
    for i in range(n):
        d = sorted(range(n), key=lambda j: abs(xs[j] - xs[i]))[:order]
        wmax = max(1e-9, abs(xs[d[-1]] - xs[i]))
        w = [(1 - (abs(xs[j] - xs[i]) / wmax) ** 3) ** 3 for j in d]
        sw = sum(w) or 1e-9
        mx = sum(wi * xs[j] for wi, j in zip(w, d)) / sw
        my = sum(wi * ys[j] for wi, j in zip(w, d)) / sw
        num = den = 0.0
        for wi, j in zip(w, d):
            num += wi * (xs[j] - mx) * (ys[j] - my)
            den += wi * (xs[j] - mx) ** 2
        slope = num / den if den else 0.0
        out.append({"x": xs[i], "y": my + slope * (xs[i] - mx), "raw": ys[i]})
    return out


def logistic_fit(doses, response) -> dict:
    """4-parameter Hill fit by gradient descent -> EC50/IC50 and hill slope."""
    x, y = _as_list(doses), _as_list(response)
    if len(x) < 5:
        return {"ec50": float("nan"), "hill_slope": float("nan"), "r_squared": 0.0}
    lo, hi = min(y), max(y)
    p = [lo, hi, float(sorted(x)[len(x) // 2]), 1.0]

    def model(xv, p):
        bot, top, ec, h = p
        return bot + (top - bot) / (1 + (ec / max(1e-12, xv)) ** h if xv else (top - bot) + bot)

    lr, best = 1e-2, None
    for _ in range(4000):
        grads = [0.0] * 4
        err = 0.0
        for xi, yi in zip(x, y):
            m = model(max(1e-9, xi), p)
            d = m - yi
            err += d * d
            top, bot = p[1], p[0]
            e = max(1e-9, p[2])
            h = p[3]
            ratio = (e / max(1e-9, xi)) ** h
            denom = (1 + ratio) ** 2
            grads[0] += d * (1 / denom)
            grads[1] += d * (1 - 1 / denom)
            grads[2] += d * (-top + bot) * h * ratio / (e * denom)
            grads[3] += d * (-(top - bot) * ratio * math.log(max(1e-9, e / max(1e-9, xi))) / denom)
        grads = [g / len(x) for g in grads]
        step = lr / (1 + _ / 1500)
        p = [p[i] - step * grads[i] for i in range(4)]
        if max(abs(g) for g in grads) < 1e-9:
            break
    ss_tot = sum((yi - mean(y)) ** 2 for yi in y)
    best = p
    return {"bottom": best[0], "top": best[1], "ec50": best[2], "hill_slope": best[3],
            "r_squared": 1 - err / ss_tot if ss_tot else 1.0, "n": len(x)}


def power_analysis(effect: float, alpha: float = 0.05, power: float = 0.8,
                   n: int | None = None) -> dict:
    za, zb = norm_ppf(1 - alpha / 2), norm_ppf(power)
    if n is None:
        n = math.ceil(2 * ((za + zb) / effect) ** 2) if effect else float("inf")
        return {"n_per_group": n, "effect_size": effect, "alpha": alpha, "power": power}
    z = abs(effect) * math.sqrt(n / 2)
    return {"power": max(0.0, min(1.0, norm_cdf(z - za))), "n_per_group": n,
            "effect_size": effect, "alpha": alpha}


def boot_ci(x, stat=mean, resamples: int = 1000, ci: float = 95, seed: int = 1) -> dict:
    v = _as_list(x)
    if len(v) < 2:
        return {"estimate": float("nan")}
    rng = random.Random(seed)
    vals = []
    for _ in range(resamples):
        sample = [rng.choice(v) for _ in range(len(v))]
        try:
            vals.append(stat(sample))
        except Exception:
            continue
    if not vals:
        return {"estimate": float("nan")}
    lo = quantile(vals, (100 - ci) / 200)
    hi = quantile(vals, 1 - (100 - ci) / 200)
    return {"estimate": stat(v), "ci_low": lo, "ci_high": hi, "ci": ci,
            "std_error": stdev(vals), "resamples": resamples}


def permutation_test(a, b, statistic=None, n: int = 2000, seed: int = 7) -> dict:
    x, y = _as_list(a), _as_list(b)
    stat = statistic or (lambda p, q: mean(p) - mean(q))
    if not x or not y:
        return {"observed": float("nan"), "p_value": float("nan")}
    obs = stat(x, y)
    pool = x + y
    rng = random.Random(seed)
    hits = 0
    for _ in range(n):
        rng.shuffle(pool)
        if abs(stat(pool[:len(x)], pool[len(x):])) >= abs(obs) - 1e-12:
            hits += 1
    return {"observed": obs, "p_value": (hits + 1) / (n + 1), "permutations": n,
            "effect_direction": "greater" if obs > 0 else "less"}


def p_adjust(pvals, method: str = "fdr_bh") -> list[float]:
    p = _as_list(pvals)
    n = len(p)
    if not n:
        return []
    order = sorted(range(n), key=lambda i: p[i])
    adj = [0.0] * n
    if method in ("fdr_bh", "fdr", "bh"):
        prev = 1.0
        for k in range(n - 1, -1, -1):
            i = order[k]
            val = p[i] * n / (k + 1)
            prev = min(prev, val)
            adj[i] = max(0.0, min(1.0, prev))
    elif method == "bonferroni":
        for i in order:
            adj[i] = min(1.0, p[i] * n)
    elif method == "holm":
        prev = 0.0
        for k, i in enumerate(order):
            prev = max(prev, min(1.0, p[i] * (n - k)))
            adj[i] = prev
    elif method == "fdr_by":
        c = sum(1.0 / i for i in range(1, n + 1))
        prev = 1.0
        for k in range(n - 1, -1, -1):
            i = order[k]
            prev = min(prev, p[i] * n * c / (k + 1))
            adj[i] = max(0.0, min(1.0, prev))
    else:
        adj = list(p)
    return [round(a, 10) for a in adj]


def shannon(counts) -> float:
    v = [c for c in _as_list(counts) if c > 0]
    tot = sum(v)
    return -sum((c / tot) * math.log(c / tot) for c in v) if tot else 0.0


def shannon2(counts) -> float:  # log base 2 variant
    v = [c for c in _as_list(counts) if c > 0]
    tot = sum(v)
    return -sum((c / tot) * math.log2(c / tot) for c in v) if tot else 0.0


def simpson(counts) -> float:
    v = _as_list(counts)
    tot = sum(v)
    return 1 - sum((c / tot) ** 2 for c in v) if tot > 1 else 0.0


def simpson_dominance(counts) -> float:
    v = _as_list(counts)
    tot = sum(v)
    return sum((c / tot) ** 2 for c in v) if tot else 0.0


def inverse_simpson(counts) -> float:
    d = simpson_dominance(counts)
    return 1 / d if d else 0.0


def pielou_evenness(counts) -> float:
    v = [c for c in _as_list(counts) if c > 0]
    h = shannon(v)
    return h / math.log(len(v)) if len(v) > 1 and h else 0.0


def berger_parker(counts) -> float:
    v = _as_list(counts)
    return max(v) / sum(v) if v and sum(v) else 0.0


def chao1(counts) -> float:
    v = _as_list(counts)
    s = len([c for c in v if c > 0])
    f1 = sum(1 for c in v if c == 1)
    f2 = sum(1 for c in v if c == 2)
    return s + f1 * (f1 - 1) / (2 * (f2 + 1))


def ace(counts) -> float:
    rare = [c for c in counts if 0 < c < 10]
    s_rare, n_rare = len(rare), sum(rare)
    if not rare:
        return float(len([c for c in counts if c > 0]))
    gamma2 = max(0.0, (n_rare * (sum((c / n_rare) ** 2 * s_rare for c in rare))) - 1) if n_rare else 0
    return sum(1 for c in counts if c > 0) + s_rare / max(1e-9, 1 - gamma2) if n_rare else s_rare


def fisher_alpha(counts) -> float:
    n = sum(counts)
    s = len([c for c in counts if c > 0])
    if n <= 0 or s <= 0:
        return 0.0
    a = 1.0
    for _ in range(100):
        lhs = a * math.log((a + n) / a) if a > 0 else 0
        if abs(lhs - s) < 1e-9:
            break
        a = a * (s / lhs) if lhs else a * 1.1
    return round(a, 6)


def rank_abundance(counts) -> list[dict]:
    v = sorted([c for c in _as_list(counts) if c > 0], reverse=True)
    tot = sum(v) or 1
    return [{"rank": i + 1, "abundance": c, "relative": c / tot,
             "cumulative": sum(v[:i + 1]) / tot} for i, c in enumerate(v)]


def rarefaction(counts, steps: int = 20, seed: int = 3) -> list[dict]:
    """Expected-species rarefaction (Mao tau style) computed by subsampling."""
    pool = []
    for i, c in enumerate(_as_list(counts)):
        pool.extend([i] * int(c))
    rng = random.Random(seed)
    n = len(pool)
    if n == 0:
        return []
    out = []
    for k in range(1, n + 1, max(1, n // steps)):
        sample = rng.sample(pool, k)
        out.append({"size": k, "richness": len(set(sample)),
                    "shannon": shannon2(list(Counter(sample).values())),
                    "total_samples": n})
    return out


def rarefied_counts(counts, depth: int, seed: int = 0) -> list[int]:
    rng = random.Random(seed)
    pool = []
    for i, c in enumerate(_as_list(counts)):
        pool.extend([i] * int(c))
    rng.shuffle(pool)
    take = Counter(pool[: min(depth, len(pool))])
    return [take.get(i, 0) for i in range(len(counts))]


def morans_i(values: list[float], neighbours: list[list[int]],
             lattice: int = 1) -> dict:
    n = len(values)
    if n < 2:
        return {"I": float("nan")}
    m = mean(values)
    v0 = sum((x - m) ** 2 for x in values)
    num = 0.0
    w = 0.0
    for i, nb in enumerate(neighbours):
        for j in nb:
            if 0 <= j < n:
                num += (values[i] - m) * (values[j] - m)
                w += 1
    if not v0 or not w:
        return {"I": 0.0}
    return {"I": (n / w) * num / v0, "expectation": -1 / (n - 1)}


def jaro_similarity(a: str, b: str) -> float:
    from difflib import SequenceMatcher

    return SequenceMatcher(None, a, b).ratio()


def hamming(a: str, b: str) -> dict:
    n = min(len(a), len(b))
    d = sum(1 for i in range(n) if a[i] != b[i]) + abs(len(a) - len(b))
    return {"distance": d, "identity": 1 - d / max(len(a), len(b), 1),
            "length": max(len(a), len(b)), "compared": n}


def levenshtein(a: str, b: str) -> int:
    if a == b:
        return 0
    if not a:
        return len(b)
    if not b:
        return len(a)
    prev = list(range(len(b) + 1))
    for i, ca in enumerate(a, 1):
        cur = [i]
        for j, cb in enumerate(b, 1):
            cur.append(min(prev[j] + 1, cur[j - 1] + 1, prev[j - 1] + (ca != cb)))
        prev = cur
    return prev[-1]


def jaccard_sets(a, b) -> dict:
    sa, sb = set(a), set(b)
    u = sa | sb
    return {"jaccard": len(sa & sb) / len(u) if u else 0.0, "intersection": len(sa & sb),
            "union": len(u), "n_a": len(sa), "n_b": len(sb),
            "sorensen": 2 * len(sa & sb) / (len(sa) + len(sb)) if (sa or sb) else 0.0}


def spearman_from_ranks(a, b) -> float:
    return spearman(a, b)["r"]


def zscores(x) -> list[float]:
    v = _as_list(x)
    s = stdev(v)
    if not s:
        return [0.0] * len(v)
    m = mean(v)
    return [(a - m) / s for a in v]


def outliers_iqr(x, k: float = 1.5) -> dict:
    v = _as_list(x)
    q1, q3 = quantile(v, 0.25), quantile(v, 0.75)
    iqr = q3 - q1
    lo, hi = q1 - k * iqr, q3 + k * iqr
    idx = [i for i, a in enumerate(v) if a < lo or a > hi]
    return {"lower": lo, "upper": hi, "n_outliers": len(idx), "indices": idx[:500],
            "values": [v[i] for i in idx[:500]]}


def outliers_mad(x, k: float = 3.5) -> dict:
    v = _as_list(x)
    m, md = median(v), mad(v)
    if not md:
        return {"n_outliers": 0, "indices": [], "values": [], "threshold": 0.0}
    z = [0.6745 * (a - m) / md for a in v]
    idx = [i for i, a in enumerate(z) if abs(a) > k]
    return {"n_outliers": len(idx), "indices": idx[:500], "values": [v[i] for i in idx[:500]],
            "threshold": k, "max_z": max(abs(a) for a in z) if z else 0.0}


def durbin_watson(residuals) -> float:
    r = _as_list(residuals)
    if len(r) < 2:
        return float("nan")
    num = sum((r[i] - r[i - 1]) ** 2 for i in range(1, len(r)))
    den = sum(x * x for x in r)
    return num / den if den else float("nan")


def vif(matrix: list[list[float]]) -> list[float]:
    import numpy as np

    X = np.array(matrix, dtype=float)
    out = []
    for j in range(X.shape[1]):
        y = X[:, j]
        A = np.delete(X, j, axis=1)
        if A.shape[1] == 0:
            out.append(1.0)
            continue
        A1 = np.column_stack([A, np.ones(len(A))])
        coef = np.linalg.lstsq(A1, y, rcond=None)[0]
        pred = A1 @ coef
        ss_res = float(((y - pred) ** 2).sum())
        ss_tot = float(((y - y.mean()) ** 2).sum())
        r2 = 1 - ss_res / ss_tot if ss_tot else 0.0
        out.append(1 / (1 - min(0.999999, r2)))
    return out


def normality_tests(x) -> dict:
    v = _as_list(x)
    jb = jarque_bera_p(v)
    s, k = skewness(v), kurtosis(v)
    n = len(v)
    ks = kstest_1samp(v)
    return {"jarque_bera": jb, "skew": s, "kurtosis": k, "n": n,
            "shapiro_like_p": jb, "ks_stat": ks["D"], "ks_p": ks["p_value"],
            "normal": jb > 0.05 if n >= 8 else None}


def jarque_bera_p(v: list[float]) -> float:
    n = len(v)
    if n < 8:
        return float("nan")
    jb = n / 6 * (skewness(v) ** 2 + (kurtosis(v) ** 2) / 4)
    return chi2_pvalue(jb, 2)


def levene(*groups) -> dict:
    gs = [_as_list(g) for g in groups if _as_list(g)]
    if len(gs) < 2:
        return {"statistic": float("nan"), "p_value": float("nan")}
    z = [abs(v - median(g)) for g in gs for v in g]
    zt = mean(z)
    zb = [mean([abs(v - median(g)) for v in g]) for g in gs]
    n = len(z)
    k = len(gs)
    num = (n - k) * sum(len(g) * (mean([abs(v - median(g)) for v in g]) - zt) ** 2
                         for g in gs)
    den = sum((abs(v - median(g)) - zb[i]) ** 2 for i, g in enumerate(gs) for v in g)
    stat = num / den if den else float("nan")
    return {"statistic": stat, "p_value": f_pvalue(stat, k - 1, n - k), "df1": k - 1,
            "df2": n - k, "variances": [round(variance(g), 6) for g in gs]}


def bartlett(*groups) -> dict:
    gs = [_as_list(g) for g in groups if len(_as_list(g)) > 1]
    k = len(gs)
    if k < 2:
        return {"statistic": float("nan"), "p_value": float("nan")}
    n = sum(len(g) for g in gs)
    sp = sum((len(g) - 1) * variance(g) for g in gs) / (n - k)
    num = (n - k) * math.log(sp) - sum((len(g) - 1) * math.log(max(1e-300, variance(g)))
                                       for g in gs)
    c = sum(1 / (len(g) - 1) for g in gs) - 1 / (n - k)
    stat = num / (1 + c / (3 * (k - 1))) if c else num
    return {"statistic": stat, "p_value": chi2_pvalue(stat, k - 1), "df": k - 1,
            "pooled_variance": sp}


def spearman_rank_correlation(a, b) -> float:
    return spearman(a, b)["r"]


def auc(scores, labels) -> float:
    """ROC AUC via Mann-Whitney equivalence."""
    pos = [s for s, l in zip(scores, labels) if l]
    neg = [s for s, l in zip(scores, labels) if not l]
    if not pos or not neg:
        return float("nan")
    gt = sum(1 for p in pos for q in neg if p > q) + 0.5 * sum(1 for p in pos for q in neg if p == q)
    return gt / (len(pos) * len(neg))


def roc_curve(scores, labels) -> list[dict]:
    pairs = sorted(zip(scores, labels), key=lambda x: -x[0])
    P = sum(1 for _, l in pairs if l)
    N = len(pairs) - P
    tp = fp = 0
    out = [{"threshold": float("inf"), "fpr": 0.0, "tpr": 0.0, "tp": 0, "fp": 0}]
    for s, l in pairs:
        if l:
            tp += 1
        else:
            fp += 1
        out.append({"threshold": s, "fpr": fp / N if N else 0.0, "tpr": tp / P if P else 0.0,
                    "tp": tp, "fp": fp})
    return out


def confusion_matrix(y_true, y_pred, labels=None) -> list[dict]:
    lab = list(labels or sorted(set(y_true) | set(y_pred)))
    idx = {v: i for i, v in enumerate(lab)}
    M = [[0] * len(lab) for _ in lab]
    for t, p in zip(y_true, y_pred):
        if t in idx and p in idx:
            M[idx[t]][idx[p]] += 1
    rows = []
    for i, l in enumerate(lab):
        tp = M[i][i]
        fp = sum(M[j][i] for j in range(len(lab))) - tp
        fn = sum(M[i]) - tp
        tn = sum(sum(r) for r in M) - tp - fp - fn
        rows.append({"label": l, "tp": tp, "fp": fp, "fn": fn, "tn": tn,
                     "precision": tp / (tp + fp) if tp + fp else 0.0,
                     "recall": tp / (tp + fn) if tp + fn else 0.0,
                     "specificity": tn / (tn + fp) if tn + fp else 0.0,
                     "f1": 2 * tp / (2 * tp + fp + fn) if (2 * tp + fp + fn) else 0.0,
                     "support": sum(M[i])})
    return rows


def mutual_information(x, y, bins: int = 10) -> float:
    xs, ys = _as_list(x), _as_list(y)
    n = min(len(xs), len(ys))
    xs, ys = xs[:n], ys[:n]
    if not n:
        return 0.0

    def dig(vals):
        lo, hi = min(vals), max(vals)
        rng = (hi - lo) or 1
        return [int((v - lo) / rng * bins) % bins for v in vals]

    dx, dy = dig(xs), dig(ys)
    cx, cy, cxy = Counter(dx), Counter(dy), Counter(zip(dx, dy))
    mi = 0.0
    for (a, b), k in cxy.items():
        pxy = k / n
        mi += pxy * math.log(pxy / ((cx[a] / n) * (cy[b] / n)))
    return mi


def information_content(counts_per_position: list[list[int]]) -> list[float]:
    out = []
    for col in counts_per_position:
        tot = sum(col) or 1
        e = -(sum((c / tot) * math.log2(c / tot) for c in col if c) or 0.0)
        out.append(2.0 - e)
    return out


def gini(x) -> float:
    v = sorted(_as_list(x))
    n = len(v)
    if n < 2:
        return 0.0
    cum = 0.0
    for i, val in enumerate(v, 1):
        cum += i * val
    tot = sum(v)
    return (2 * cum) / (n * tot) - (n + 1) / n if tot else 0.0


def ewma(x, alpha: float = 0.3) -> list[float]:
    v = _as_list(x)
    out = []
    cur = None
    for a in v:
        cur = a if cur is None else alpha * a + (1 - alpha) * cur
        out.append(cur)
    return out


def cusum(x) -> list[float]:
    v = _as_list(x)
    m = mean(v)
    out, s = [], 0.0
    for a in v:
        s += a - m
        out.append(s)
    return out


def cross_correlation(a, b, max_lag: int = 10) -> list[dict]:
    x, y = _as_list(a), _as_list(b)
    n = min(len(x), len(y))
    x, y = x[:n], y[:n]
    out = []
    for lag in range(-max_lag, max_lag + 1):
        if lag >= 0:
            xs, ys = x[lag:], y[:n - lag]
        else:
            xs, ys = x[:n + lag], y[-lag:]
        out.append({"lag": lag, "correlation": pearson(xs, ys)["r"] if len(xs) > 2 else 0.0})
    return out
