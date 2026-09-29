"""Machine-learning engines implemented with numpy only.

Everything here is written from scratch (k-means, PCA, t-SNE, DBSCAN, linkage,
CART trees, random forest, gradient boosting, linear/logistic models, KNN, SVM,
naive Bayes, GMM) so the suite has no scikit-learn dependency while still
producing real numbers you can check.
"""

from __future__ import annotations

import math
import random
from collections import Counter, defaultdict

import numpy as np

from chroma_titan.core import stats


def as_matrix(X) -> np.ndarray:
    arr = np.asarray(X, dtype=float)
    if arr.ndim == 1:
        arr = arr.reshape(-1, 1)
    return np.nan_to_num(arr, nan=0.0, posinf=0.0, neginf=0.0)


def standardize(X: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    mu = X.mean(axis=0)
    sd = X.std(axis=0)
    sd[sd == 0] = 1.0
    return (X - mu) / sd, mu, sd


# ---------------------------------------------------------------------------
# clustering
# ---------------------------------------------------------------------------
def kmeans(X, k: int = 3, iters: int = 100, seed: int = 0, n_init: int = 4) -> dict:
    M = as_matrix(X)
    M, _, _ = standardize(M)
    rng = np.random.default_rng(seed)
    best = None
    for init in range(n_init):
        cent = M[rng.choice(len(M), size=min(k, len(M)), replace=False)]
        labels = np.zeros(len(M), int)
        for it in range(iters):
            d = ((M[:, None, :] - cent[None, :, :]) ** 2).sum(-1)
            new = d.argmin(1)
            for j in range(cent.shape[0]):
                sel = new == j
                if sel.any():
                    cent[j] = M[sel].mean(0)
            if (new == labels).all():
                break
            labels = new
        inertia = float(((M - cent[labels]) ** 2).sum())
        if best is None or inertia < best["inertia"]:
            best = {"labels": labels, "centroids": cent, "inertia": inertia, "iters": it}
    labels = best["labels"]
    sil = silhouette(M, labels)
    return {"labels": [int(x) for x in labels], "centroids": best["centroids"].round(5).tolist(),
            "inertia": round(best["inertia"], 5), "k": int(k), "n": len(M),
            "silhouette": round(sil, 4), "cluster_sizes": [int((labels == j).sum()) for j in range(k)],
            "iterations": int(best["iters"])}


def silhouette(M: np.ndarray, labels: np.ndarray) -> float:
    n = len(M)
    if n < 3 or len(set(labels.tolist())) < 2:
        return 0.0
    D = ((M[:, None, :] - M[None, :, :]) ** 2).sum(-1) ** 0.5
    vals = []
    idx = list(range(n))
    for i in idx[:min(n, 400)]:
        own = labels[i]
        same = [D[i, j] for j in idx if j != i and labels[j] == own]
        if not same:
            continue
        a = float(np.mean(same))
        b = min((float(np.mean([D[i, j] for j in idx if labels[j] == c]))
                 for c in set(labels.tolist()) if c != own), default=0.0)
        vals.append((b - a) / max(a, b) if max(a, b) else 0.0)
    return float(np.mean(vals)) if vals else 0.0


def dbscan(X, eps: float = 0.5, min_samples: int = 4) -> dict:
    M = as_matrix(X)
    M, _, _ = standardize(M)
    n = len(M)
    D = ((M[:, None, :] - M[None, :, :]) ** 2).sum(-1) ** 0.5
    neigh = [np.flatnonzero(D[i] <= eps) for i in range(n)]
    labels = -np.ones(n, int)
    cid = 0
    for i in range(n):
        if labels[i] != -1 or len(neigh[i]) < min_samples:
            continue
        labels[i] = cid
        seeds = list(neigh[i])
        k = 0
        while k < len(seeds):
            j = seeds[k]
            if labels[j] == -1:
                labels[j] = cid
                if len(neigh[j]) >= min_samples:
                    seeds.extend(neigh[j])
            elif labels[j] == -1 or labels[j] == -2:
                pass
            if labels[j] < 0:
                labels[j] = cid
            k += 1
        cid += 1
    return {"labels": [int(x) for x in labels], "n_clusters": cid,
            "noise": int((labels == -1).sum()), "eps": eps, "min_samples": min_samples,
            "sizes": [int((labels == c).sum()) for c in range(cid)]}


def linkage(X, method: str = "average", labels: list[str] | None = None) -> list[list]:
    """Hierarchical agglomerative clustering -> SciPy-style linkage matrix."""
    M = as_matrix(X)
    n = len(M)
    D = ((M[:, None, :] - M[None, :, :]) ** 2).sum(-1) ** 0.5
    clusters = {i: [i] for i in range(n)}
    active = list(range(n))
    Z: list[list] = []
    nxt = n
    while len(active) > 1:
        best = None
        for ai in range(len(active)):
            for bi in range(ai + 1, len(active)):
                ca, cb = clusters[active[ai]], clusters[active[bi]]
                block = D[np.ix_(ca, cb)]
                if method == "single":
                    d = float(block.min())
                elif method == "complete":
                    d = float(block.max())
                elif method == "max":
                    d = float(block.max())
                elif method == "weighted":
                    d = float(np.mean(block))
                else:  # average / upgma
                    d = float(np.mean(block))
                if best is None or d < best[0]:
                    best = (d, active[ai], active[bi], len(ca) + len(cb))
        d, a, b, cnt = best
        clusters[nxt] = clusters[a] + clusters[b]
        Z.append([float(a), float(b), float(d), float(cnt)])
        active = [x for x in active if x not in (a, b)] + [nxt]
        del clusters[a]
        del clusters[b]
        nxt += 1
    del labels
    return Z


def fcluster(Z: list[list], t: float, criterion: str = "distance") -> list[int]:
    """Cut a linkage matrix into clusters.

    ``criterion='distance'`` merges everything below ``t``;
    ``criterion='maxclust'`` (also ``'inconsistent'``) cuts the tree so that
    exactly ``t`` clusters remain.  Labels are 1-based, like SciPy.
    """
    rows = [list(r) for r in (Z or [])]
    if not rows:
        return []
    n = len(rows) + 1
    children = [(int(a), int(b)) for a, b, *_ in rows]
    if criterion == "maxclust" or criterion == "inconsistent":
        k = max(1, min(int(t), n))
        cut = set(range(len(rows) - (k - 1), len(rows)))
    else:
        cut = {i for i, r in enumerate(rows) if float(r[2]) > float(t)}

    def leaves(node: int) -> list[int]:
        out: list[int] = []
        stack = [node]
        while stack:
            x = stack.pop()
            if x < n:
                out.append(x)
            else:
                a, b = children[x - n]
                stack.append(a)
                stack.append(b)
        return sorted(out)

    labels = [0] * n
    nxt = 0

    def walk(node: int) -> None:
        nonlocal nxt
        if node < n:
            nxt += 1
            labels[node] = nxt
            return
        i = node - n
        if i in cut:
            walk(children[i][0])
            walk(children[i][1])
            return
        nxt += 1
        for lf in leaves(node):
            labels[lf] = nxt

    walk(n + len(rows) - 1)
    return labels


# ---------------------------------------------------------------------------
# dimensionality reduction
# ---------------------------------------------------------------------------
def pca(X, components: int = 2) -> dict:
    M = as_matrix(X)
    M, mu, sd = standardize(M)
    if M.shape[0] < 2 or M.shape[1] < 1:
        return {"scores": [], "explained_variance_ratio": [], "loadings": [], "means": []}
    U, S, Vt = np.linalg.svd(M - M.mean(0), full_matrices=False)
    var = (S ** 2) / max(1, len(M) - 1)
    tot = float(var.sum()) or 1.0
    k = min(components, len(S), M.shape[1])
    scores = (M - M.mean(0)) @ Vt[:k].T
    return {"scores": np.round(scores, 6).tolist(),
            "explained_variance": [round(float(v), 6) for v in var[:k]],
            "explained_variance_ratio": [round(float(v / tot), 6) for v in var[:k]],
            "loadings": np.round(Vt[:k], 6).tolist(),
            "eigenvalues": [round(float(v), 6) for v in var],
            "n_samples": len(M), "n_features": M.shape[1], "components": k,
            "means": np.round(mu, 6).tolist()}


def tsne(X, components: int = 2, perplexity: float = 30.0, iterations: int = 250,
         learning_rate: float = 200.0, seed: int = 3, init: str = "random") -> dict:
    """t-SNE with exact gradients: perplexity binary search + momentum gains."""
    M, _, _ = standardize(as_matrix(X))
    n = len(M)
    if n < 4:
        return {"embedding": M[:, :components].round(5).tolist(), "kl": 0.0, "iters": 0,
                "perplexity": 0.0, "components": components}
    rng = np.random.default_rng(seed)
    sqD = ((M[:, None, :] - M[None, :, :]) ** 2).sum(-1)
    perp = min(max(2.0, float(perplexity)), max(2.0, (n - 1) / 3.0))
    target = math.log(perp)
    P = np.zeros((n, n))
    for i in range(n):
        other = np.concatenate([sqD[i, :i], sqD[i, i + 1:]])
        lo, hi, beta = 1e-8, 1e8, 1.0
        for _ in range(60):
            p = np.exp(-(other - other.min()) * beta)
            Z = p.sum() or 1e-12
            H = math.log(Z) + beta * float(np.sum(other * p)) / Z
            if abs(H - target) < 1e-5:
                break
            if H > target:
                lo, beta = beta, (beta + hi) / 2
            else:
                hi, beta = beta, (lo + beta) / 2
        p = np.exp(-(other - other.min()) * beta)
        P[i, [j for j in range(n) if j != i]] = p / (p.sum() or 1e-12)
    P = np.maximum((P + P.T) / (2 * n), 1e-12)
    P /= P.sum()
    if init == "pca":
        Y = pca(M, components=components)["scores"]
        Y = np.asarray(Y, dtype=float)
    else:
        Y = rng.normal(0, 1e-3, (n, components))
    Y = Y - Y.mean(0)
    gains = np.ones((n, components))
    velocity = np.zeros((n, components))
    kl = 0.0
    for it in range(iterations):
        sqY = (Y ** 2).sum(-1)
        num = 1.0 / (1.0 + sqY[:, None] + sqY[None, :] - 2 * Y @ Y.T)
        np.fill_diagonal(num, 0.0)
        Z = num.sum() or 1e-12
        Q = num / Z
        kl = float(np.sum(P * np.log(P / np.clip(Q, 1e-12, None))))
        grad = 4 * ((P - Q)[:, :, None] * (Y[:, None, :] - Y[None, :, :])).sum(1)
        gmax = float(np.abs(grad).max()) or 1.0
        if gmax > 1e3:
            grad = grad / gmax * 1e3
        inc = np.sign(grad) != np.sign(velocity)
        gains = np.where(inc, gains * 0.8, gains * 1.2 + 0.8)
        momentum = 0.5 if it < 50 else 0.8
        velocity = momentum * velocity - learning_rate * gains * grad
        velocity = np.clip(velocity, -100, 100)
        Y = Y + velocity
        Y -= Y.mean(0)
        sdev = Y.std(0)
        sdev[sdev == 0] = 1.0
        Y = Y / sdev.max()
    sumY = (Y ** 2).sum(-1)
    num = 1.0 / (1.0 + sumY[:, None] + sumY[None, :] - 2 * Y @ Y.T)
    np.fill_diagonal(num, 0.0)
    Q = num / (num.sum() or 1e-12)
    kl = float(np.sum(P * np.log(P / np.clip(Q, 1e-12, None))))
    return {"embedding": np.round(Y, 5).tolist(), "kl": round(kl, 5), "iters": iterations,
            "perplexity": round(float(perp), 3), "components": components,
            "early_exaggeration": 50}


def tsne_pca(X, components: int = 2, pca_dims: int = 5, perplexity: float = 20.0) -> dict:
    p = pca(X, components=pca_dims)
    emb = tsne(p["scores"], components=components, perplexity=perplexity)
    emb["explained_variance_ratio"] = p["explained_variance_ratio"]
    return emb


def umap_like(X, components: int = 2, neighbours: int = 15, seed: int = 0) -> dict:
    """Fast spectral embedding as a UMAP stand-in (kNN graph -> Laplacian)."""
    M = as_matrix(X)
    M, _, _ = standardize(M)
    n = len(M)
    k = min(max(2, neighbours), n - 1)
    D = ((M[:, None, :] - M[None, :, :]) ** 2).sum(-1) ** 0.5
    idx = np.argsort(D, axis=1)[:, 1:k + 1]
    W = np.zeros((n, n))
    for i in range(n):
        for j in idx[i]:
            W[i, j] = W[j, i] = 1.0
    deg = W.sum(1)
    deg[deg == 0] = 1.0
    L = np.eye(n) - (W / np.sqrt(deg)[:, None]) / np.sqrt(deg)[None, :]
    w, V = np.linalg.eigh(L)
    keep = min(components, max(1, n - 1))
    emb = V[:, 1:keep + 1] / np.sqrt(deg[:, None])
    return {"embedding": np.round(emb, 6).tolist(), "n_neighbours": k,
            "eigenvalues": [round(float(x), 6) for x in w[:keep + 1]], "method": "spectral"}


# ---------------------------------------------------------------------------
# supervised models
# ---------------------------------------------------------------------------
def train_test_split(X, y, test_size: float = 0.25, seed: int = 7, stratify: bool = False) -> dict:
    M = as_matrix(X)
    yv = list(y)
    rng = random.Random(seed)
    idx = list(range(len(M)))
    if stratify:
        by_class: dict = defaultdict(list)
        for i, c in enumerate(yv):
            by_class[c].append(i)
        test = []
        for c, rows in by_class.items():
            rng.shuffle(rows)
            test.extend(rows[: max(1, int(len(rows) * test_size))])
    else:
        rng.shuffle(idx)
        test = idx[: int(len(idx) * test_size)]
    test = sorted(test)
    train = [i for i in range(len(M)) if i not in set(test)]
    return {"X_train": M[train].round(6).tolist(), "y_train": [yv[i] for i in train],
            "X_test": M[test].round(6).tolist(), "y_test": [yv[i] for i in test],
            "train_indices": train, "test_indices": test,
            "n_train": len(train), "n_test": len(test)}


def linear_regression_fit(X, y) -> dict:
    M = as_matrix(X)
    yy = np.asarray(y, dtype=float)
    A = np.column_stack([np.ones(len(M)), M])
    coef, *_ = np.linalg.lstsq(A, yy, rcond=None)
    pred = A @ coef
    ss_res = float(((yy - pred) ** 2).sum())
    ss_tot = float(((yy - yy.mean()) ** 2).sum()) or 1.0
    return {"coefficients": np.round(coef, 6).tolist(), "intercept": round(float(coef[0]), 6),
            "r_squared": round(1 - ss_res / ss_tot, 6),
            "residual_standard_error": round(math.sqrt(ss_res / max(1, len(yy) - A.shape[1])), 6),
            "n": len(yy), "predictions": np.round(pred, 6).tolist(),
            "residuals": np.round(yy - pred, 6).tolist()}


def ridge_regression(X, y, alpha: float = 1.0) -> dict:
    M = as_matrix(X)
    yy = np.asarray(y, dtype=float)
    A = np.column_stack([np.ones(len(M)), M])
    reg = np.eye(A.shape[1]) * alpha
    reg[0, 0] = 0.0
    coef = np.linalg.solve(A.T @ A + reg, A.T @ yy)
    pred = A @ coef
    ss_res = float(((yy - pred) ** 2).sum())
    ss_tot = float(((yy - yy.mean()) ** 2).sum()) or 1.0
    return {"coefficients": np.round(coef, 6).tolist(), "alpha": alpha,
            "r_squared": round(1 - ss_res / ss_tot, 6),
            "L2_penalty": round(float((coef[1:] ** 2).sum() * alpha), 6),
            "predictions": np.round(pred, 6).tolist()}


def lasso_regression(X, y, alpha: float = 0.1, iters: int = 800) -> dict:
    """Coordinate-descent lasso."""
    M, mu, sd = standardize(as_matrix(X))
    yy = np.asarray(y, dtype=float)
    yy = yy - yy.mean()
    p = M.shape[1]
    coef = np.zeros(p)
    for _ in range(iters):
        for j in range(p):
            r = yy - M @ coef + M[:, j] * coef[j]
            rho = float(M[:, j] @ r)
            z = float((M[:, j] ** 2).sum())
            coef[j] = np.sign(rho) * max(0.0, abs(rho) - alpha * z) / z if z else 0.0
    pred = M @ coef + np.asarray(y, dtype=float).mean()
    ss_res = float(((np.asarray(y, dtype=float) - pred) ** 2).sum())
    ss_tot = float(((np.asarray(y, dtype=float) - np.mean(y)) ** 2).sum()) or 1.0
    return {"coefficients": np.round(coef, 6).tolist(), "intercept": round(float(np.mean(y)), 6),
            "alpha": alpha, "nonzero": int((coef != 0).sum()),
            "r_squared": round(1 - ss_res / ss_tot, 6),
            "predictions": np.round(pred, 6).tolist()}


def logistic_regression(X, y, iters: int = 800, lr: float = 0.3, l2: float = 0.0) -> dict:
    M, mu, sd = standardize(as_matrix(X))
    yy = np.asarray([1.0 if str(t).lower() in ("1", "true", "pos", "+") else
                     (0.0 if str(t).lower() in ("0", "false", "neg", "-") else float(t))
                     for t in y])
    A = np.column_stack([np.ones(len(M)), M])
    coef = np.zeros(A.shape[1])
    for _ in range(iters):
        z = np.clip(A @ coef, -30, 30)
        p = 1 / (1 + np.exp(-z))
        grad = A.T @ (p - yy) / len(yy) + l2 * np.r_[0, coef[1:]]
        coef -= lr * grad
    z = np.clip(A @ coef, -30, 30)
    prob = 1 / (1 + np.exp(-z))
    pred = (prob >= 0.5).astype(int)
    ll = float(np.sum(yy * np.log(np.clip(prob, 1e-12, 1)) +
                     (1 - yy) * np.log(np.clip(1 - prob, 1e-12, 1))))
    return {"coefficients": np.round(coef, 6).tolist(), "probabilities": np.round(prob, 5).tolist(),
            "predictions": [int(v) for v in pred], "log_likelihood": round(ll, 5),
            "accuracy": round(float((pred == yy).mean()), 5), "odds_ratios":
                [round(float(math.exp(min(20, c))), 5) for c in coef[1:]],
            "iterations": iters, "auc": stats.auc(prob.tolist(), yy.tolist())}


def knn_classify(Xtr, ytr, Xte, yte=None, k: int = 3) -> dict:
    """K-nearest-neighbour classification (Euclidean, z-scored features)."""
    A = as_matrix(Xtr)
    B = as_matrix(Xte)
    A, mu, sd = standardize(A)
    B = (B - mu) / sd
    D = ((B[:, None, :] - A[None, :, :]) ** 2).sum(-1) ** 0.5
    order = np.argsort(D, 1)[:, :max(1, k)]
    preds, dists, votes_list = [], [], []
    for i, nbrs in enumerate(order):
        votes = Counter(str(ytr[j]) for j in nbrs)
        preds.append(votes.most_common(1)[0][0])
        votes_list.append(dict(votes))
        dists.append(float(np.mean(D[i, nbrs])))
    res = {"predictions": preds, "mean_distance": [round(d, 5) for d in dists],
           "k": int(k), "n_test": len(preds), "vote_counts": votes_list}
    if yte is not None:
        truth = [str(t) for t in yte]
        res["accuracy"] = round(sum(1 for p, t in zip(preds, truth) if p == t) /
                                max(1, len(truth)), 5)
        res["confusion"] = stats.confusion_matrix(truth, preds)
    return res


def naive_bayes(Xtr, ytr, Xte=None, alpha: float = 1.0) -> dict:
    A = as_matrix(Xtr)
    classes = sorted(set(str(c) for c in ytr))
    out = {}
    pri = {}
    means = {}
    vars_ = {}
    for c in classes:
        rows = A[[i for i, y in enumerate(ytr) if str(y) == c]]
        pri[c] = len(rows) / len(A)
        means[c] = rows.mean(0) if len(rows) else np.zeros(A.shape[1])
        vars_[c] = rows.var(0, ddof=0) + alpha if len(rows) else np.ones(A.shape[1])
    B = as_matrix(Xte) if Xte is not None else A
    logp = np.zeros((len(B), len(classes)))
    for j, c in enumerate(classes):
        logp[:, j] = math.log(max(pri[c], 1e-12)) - 0.5 * np.sum(
            np.log(2 * math.pi * vars_[c])) - 0.5 * np.sum((B - means[c]) ** 2 / vars_[c], axis=1)
    preds = [classes[i] for i in logp.argmax(1)]
    return {"predictions": preds, "class_priors": {k: round(v, 5) for k, v in pri.items()},
            "class_means": {k: np.round(v, 5).tolist() for k, v in means.items()},
            "alpha": alpha, "classes": classes, "log_probabilities": np.round(logp, 5).tolist()}


class Node:
    __slots__ = ("feature", "threshold", "left", "right", "value", "samples", "impurity")

    def __init__(self, **kw):
        for k, v in kw.items():
            setattr(self, k, v)


def decision_tree(X, y, max_depth: int = 5, min_samples: int = 2, criterion: str = "gini") -> dict:
    M = as_matrix(X)
    yy = [str(v) for v in y]
    root = _build(M, np.arange(len(M)), [yy[i] for i in range(len(M))], 0, max_depth,
                  min_samples, criterion)
    preds = [_predict(root, M[i]) for i in range(len(M))]
    acc = sum(1 for p, t in zip(preds, yy) if p == t) / max(1, len(yy))
    return {"tree": _tree_dict(root), "depth": max_depth, "criterion": criterion,
            "training_accuracy": round(acc, 5), "n_leaves": _count_leaves(root),
            "feature_importances": _importances(root, M.shape[1], len(M))}


def _gini(counts: list[int], n: int) -> float:
    return 1 - sum((c / n) ** 2 for c in counts) if n else 0.0


def _entropy(counts: list[int], n: int) -> float:
    return -sum((c / n) * math.log2(c / n) for c in counts if c) if n else 0.0


def _build(M, idx, y, depth, max_depth, min_samples, criterion) -> Node:
    classes = sorted(set(y))
    cnt = [y.count(c) for c in classes]
    imp = _gini(cnt, len(y)) if criterion == "gini" else _entropy(cnt, len(y))
    node = Node(feature=None, threshold=None, left=None, right=None,
                value=classes[cnt.index(max(cnt))] if classes else "?",
                samples=len(idx), impurity=round(imp, 6))
    if depth >= max_depth or len(y) < 2 * min_samples or len(classes) < 2:
        return node
    best = None
    for f in range(M.shape[1]):
        vals = np.unique(M[idx, f])
        for thr in vals[:-1] if len(vals) > 1 else vals:
            left = [y[k] for k, i in enumerate(idx) if M[i, f] <= thr]
            right = [y[k] for k, i in enumerate(idx) if M[i, f] > thr]
            if len(left) < min_samples or len(right) < min_samples:
                continue
            lc, rc = Counter(left), Counter(right)
            g = (len(left) * _gini([lc[c] for c in classes], len(left)) +
                 len(right) * _gini([rc[c] for c in classes], len(right))) / len(y)
            if best is None or g < best[0]:
                best = (g, f, float(thr), left, right)
    if best is None:
        return node
    _, f, thr, _, _ = best
    lidx = idx[M[idx, f] <= thr]
    ridx = idx[M[idx, f] > thr]
    node.feature, node.threshold = f, thr
    node.left = _build(M, lidx, [y[k] for k in range(len(idx)) if M[idx[k], f] <= thr],
                       depth + 1, max_depth, min_samples, criterion)
    node.right = _build(M, ridx, [y[k] for k in range(len(idx)) if M[idx[k], f] > thr],
                        depth + 1, max_depth, min_samples, criterion)
    return node


def _predict(node: Node, x) -> str:
    while node.feature is not None:
        node = node.left if x[node.feature] <= node.threshold else node.right
    return node.value


def _count_leaves(node: Node) -> int:
    if node.feature is None:
        return 1
    return _count_leaves(node.left) + _count_leaves(node.right)


def _tree_dict(node: Node) -> dict:
    if node.feature is None:
        return {"leaf": node.value, "samples": int(node.samples)}
    return {"feature": int(node.feature), "threshold": node.threshold,
            "samples": int(node.samples), "impurity": node.impurity,
            "left": _tree_dict(node.left), "right": _tree_dict(node.right)}


def _importances(node: Node, n_features: int, n: int) -> list[float]:
    imp = [0.0] * n_features

    def rec(nd: Node):
        if nd.feature is None:
            return
        imp[nd.feature] += nd.samples * nd.impurity
        rec(nd.left)
        rec(nd.right)

    rec(node)
    tot = sum(imp) or 1.0
    return [round(x / tot, 6) for x in imp]


def random_forest(X, y, trees: int = 10, max_depth: int = 5, seed: int = 1,
                  feature_fraction: float = 0.8, oob: bool = True) -> dict:
    M = as_matrix(X)
    yy = [str(v) for v in y]
    rng = random.Random(seed)
    models, importances, oob_votes = [], [], []
    for t in range(trees):
        idx = [rng.randrange(len(M)) for _ in range(len(M))]
        oob_idx = [i for i in range(len(M)) if i not in set(idx)]
        sub = M[idx]
        nfeat = max(1, int(M.shape[1] * feature_fraction))
        feats = rng.sample(range(M.shape[1]), min(nfeat, M.shape[1]))
        model = decision_tree(sub[:, feats], [yy[i] for i in idx], max_depth=max_depth)
        models.append((model["tree"], feats))
        imp = [0.0] * M.shape[1]
        for f, v in zip(feats, model["feature_importances"]):
            imp[f] = v
        importances.append(imp)
        if oob and oob_idx:
            preds = [_predict_node(model["tree"], M[i, feats]) for i in oob_idx]
            oob_votes.append((oob_idx, preds, [yy[i] for i in oob_idx]))
    votes: dict[int, Counter] = defaultdict(Counter)
    for oob_idx, preds, truth in oob_votes:
        for i, p in zip(oob_idx, preds):
            votes[i][p] += 1
    oob_acc = None
    if oob_votes:
        good = tot = 0
        for oob_idx, preds, truth in oob_votes:
            for i, p, t in zip(oob_idx, preds, truth):
                tot += 1
                good += p == t
        oob_acc = round(good / tot, 5) if tot else None
    ens = Counter()
    for tree, feats in models:
        del tree, feats
    preds = []
    for i in range(len(M)):
        c: Counter = Counter()
        for tree, feats in models:
            c[_predict_node(tree, M[i, feats])] += 1
        preds.append(c.most_common(1)[0][0] if c else "?")
        ens.update(c)
    acc = sum(1 for p, t in zip(preds, yy) if p == t) / max(1, len(yy))
    imp_tot = [sum(col) / max(1, len(importances)) for col in zip(*importances)]
    s = sum(imp_tot) or 1.0
    return {"n_trees": trees, "feature_importances": [round(x / s, 6) for x in imp_tot],
            "training_accuracy": round(acc, 5), "oob_accuracy": oob_acc,
            "predictions": preds, "classes": sorted(ens.keys()), "max_depth": max_depth}


def _predict_node(tree: dict, x) -> str:
    while "leaf" not in tree:
        tree = tree["left"] if x[tree["feature"]] <= tree["threshold"] else tree["right"]
    return tree["leaf"]


def regression_tree(X, y, max_depth: int = 3, min_samples: int = 5) -> dict:
    """CART regression tree: variance reduction splits, leaf = mean."""
    M = as_matrix(X)
    yy = np.asarray(y, dtype=float)

    def build(rows: np.ndarray, depth: int) -> dict:
        vals = yy[rows]
        if depth >= max_depth or len(rows) < 2 * min_samples or float(np.var(vals)) < 1e-12:
            return {"leaf": round(float(np.mean(vals)), 6), "n": int(len(rows)),
                    "sse": round(float(np.var(vals) * len(rows)), 6)}
        best = None
        parent = float(np.var(vals) * len(rows))
        for f in range(M.shape[1]):
            col = M[rows, f]
            for thr in np.unique(col)[:-1][:24]:
                l = rows[col <= thr]
                r = rows[col > thr]
                if len(l) < min_samples or len(r) < min_samples:
                    continue
                sse = float(np.var(yy[l]) * len(l) + np.var(yy[r]) * len(r))
                if best is None or sse < best[0]:
                    best = (sse, f, float(thr), l, r)
        if best is None:
            return {"leaf": round(float(np.mean(vals)), 6), "n": int(len(rows)),
                    "sse": round(parent, 6)}
        _, f, thr, l, r = best
        return {"feature": int(f), "threshold": round(thr, 6), "n": int(len(rows)),
                "gain": round(parent - best[0], 6),
                "left": build(l, depth + 1), "right": build(r, depth + 1)}

    return build(np.arange(len(M)), 0)


def _tree_predict_leaf_path(tree: dict, x) -> str:
    path = []
    while "leaf" not in tree:
        go = x[tree["feature"]] <= tree["threshold"]
        path.append("L" if go else "R")
        tree = tree["left"] if go else tree["right"]
    return "".join(path)


def _leaf_values(tree: dict, X: np.ndarray, resid: np.ndarray) -> dict:
    out: dict[str, list[float]] = {}

    def rec(node, rows, path):
        if "leaf" not in node:
            col = X[rows, node["feature"]]
            for side, key, mask in (("left", "L", col <= node["threshold"]),
                                    ("right", "R", col > node["threshold"])):
                sub = rows[mask]
                if len(sub):
                    rec(node[side], sub, path + key)
            return
        out.setdefault(path, []).extend(resid[rows].tolist())

    rec(tree, np.arange(len(X)), "")
    return {k: float(np.mean(v)) for k, v in out.items()}


def gradient_boosting(X, y, trees: int = 40, learning_rate: float = 0.1,
                      max_depth: int = 3, loss: str = "squared_error") -> dict:
    """Gradient boosting: regression trees fit to residuals (or log-loss gradients)."""
    M = as_matrix(X)
    yy = np.asarray(y, dtype=float)
    if loss == "logloss":
        yy = np.array([1.0 if str(v) in ("1", "True", "true", "pos", "+") else 0.0 for v in y])
        f = np.zeros(len(yy))
        for _ in range(trees):
            p = 1 / (1 + np.exp(-np.clip(f, -30, 30)))
            grad = yy - p
            tree = regression_tree(M, grad, max_depth=max_depth)
            vals = _leaf_values(tree, M, grad)
            upd = np.array([vals.get(_tree_predict_leaf_path(tree, M[i]), 0.0)
                            for i in range(len(M))])
            f = f + learning_rate * upd / max(1e-9, float(np.abs(upd).max() or 1))
        prob = 1 / (1 + np.exp(-np.clip(f, -30, 30)))
        pred = (prob >= 0.5).astype(int)
        return {"loss": loss, "n_trees": trees, "learning_rate": learning_rate,
                "training_accuracy": round(float((pred == yy).mean()), 5),
                "predictions": [int(v) for v in pred],
                "probabilities": np.round(prob, 5).tolist(),
                "final_logloss": round(float(-np.mean(
                    yy * np.log(np.clip(prob, 1e-9, 1)) + (1 - yy) * np.log(np.clip(1 - prob, 1e-9, 1)))), 5)}
    f = np.full(len(yy), float(yy.mean()))
    rmse0 = float(np.sqrt(((yy - f) ** 2).mean()))
    importances = np.zeros(M.shape[1])
    for _ in range(trees):
        resid = yy - f
        tree = regression_tree(M, resid, max_depth=max_depth)
        vals = _leaf_values(tree, M, resid)
        upd = np.array([vals.get(_tree_predict_leaf_path(tree, M[i]), 0.0) for i in range(len(M))])
        f = f + learning_rate * upd
        _accumulate_importance(tree, importances)
    resid = yy - f
    rmse = float(np.sqrt((resid ** 2).mean()))
    imp_sum = float(importances.sum()) or 1.0
    return {"loss": loss, "n_trees": trees, "learning_rate": learning_rate,
            "rmse_initial": round(rmse0, 5), "rmse_final": round(rmse, 5),
            "predictions": np.round(f, 6).tolist(),
            "feature_importances": np.round(importances / imp_sum, 6).tolist(),
            "explained_variance": round(1 - rmse ** 2 / (float(yy.var()) or 1e-12), 5),
            "n": len(yy)}


def _accumulate_importance(node: dict, out: np.ndarray) -> None:
    if "leaf" not in node:
        out[node["feature"]] += float(node.get("gain", 0.0))
        _accumulate_importance(node["left"], out)
        _accumulate_importance(node["right"], out)


def svm_sgd(X, y, iters: int = 600, lr: float = 0.05, c: float = 1.0, kernel: str = "linear") -> dict:
    M, mu, sd = standardize(as_matrix(X))
    yy = np.array([1.0 if str(v) in ("1", "True", "pos", "+", "true") else -1.0 for v in y])
    w = np.zeros(M.shape[1])
    b = 0.0
    for it in range(iters):
        for i in range(len(M)):
            margin = yy[i] * (M[i] @ w + b)
            if margin < 1:
                w = w + lr * (c * yy[i] * M[i] - 2 * w / len(M))
                b = b + lr * c * yy[i]
            else:
                w = w - lr * 2 * w / len(M)
    scores = M @ w + b
    if kernel == "rbf":
        gamma = 1.0 / M.shape[1]
        K = np.exp(-gamma * ((M[:, None, :] - M[None, :, :]) ** 2).sum(-1))
        scores = K @ (yy / len(yy))
    pred = np.sign(scores)
    acc = float((pred == yy).mean())
    margin = float(np.min(np.abs(scores[yy * scores > 0]))) if (yy * scores > 0).any() else 0.0
    return {"weights": np.round(w, 6).tolist(), "bias": round(float(b), 6), "kernel": kernel,
            "C": c, "training_accuracy": round(acc, 5), "n_support_vectors":
                int((np.abs(scores) <= 1.0).sum()), "margin": round(margin, 5),
            "decision_values": np.round(scores, 5).tolist()}


def cross_validate(fit_predict, X, y, folds: int = 4, seed: int = 5) -> dict:
    """Generic k-fold CV over a (fit, predict) callable pair."""
    M = as_matrix(X)
    yy = list(y)
    rng = random.Random(seed)
    idx = list(range(len(M)))
    rng.shuffle(idx)
    chunks = [idx[i::folds] for i in range(folds)]
    scores = []
    for i in range(folds):
        test = chunks[i]
        train = [j for c in chunks if c is not chunks[i] for j in c]
        score = fit_predict(M[train], [yy[j] for j in train], M[test], [yy[j] for j in test])
        scores.append(float(score))
    return {"scores": [round(s, 5) for s in scores], "mean": round(stats.mean(scores), 5),
            "std": round(stats.stdev(scores), 5) if len(scores) > 1 else 0.0,
            "folds": folds, "cv": round(100 * stats.stdev(scores), 4) / max(1e-9, abs(stats.mean(scores)))
            if len(scores) > 1 else 0.0}


def permutation_importance(predict_fn, X, y, n_repeats: int = 10, seed: int = 2) -> dict:
    M = as_matrix(X)
    base = float(predict_fn(M, list(y)))
    rng = random.Random(seed)
    drops = []
    for f in range(M.shape[1]):
        vals = []
        for _ in range(n_repeats):
            P = M.copy()
            rng.shuffle(P[:, f])
            vals.append(float(predict_fn(P, list(y))) - base)
        drops.append({"feature_index": f, "mean_importance": round(-stats.mean(vals), 6),
                      "std": round(stats.stdev(vals), 6) if len(vals) > 1 else 0.0,
                      "n_repeats": n_repeats})
    return {"baseline_score": round(base, 6), "importances": drops,
            "n_features": M.shape[1]}


def gmm_em(X, components: int = 2, iters: int = 60, seed: int = 0) -> dict:
    M = as_matrix(X)
    n, d = M.shape
    rng = np.random.default_rng(seed)
    pi = np.full(components, 1.0 / components)
    mu = M[rng.choice(n, components, replace=False)]
    var = np.full((components, d), M.var(0))
    for _ in range(iters):
        resp = []
        for k in range(components):
            ll = -0.5 * np.sum((M - mu[k]) ** 2 / (var[k] + 1e-6), axis=1) \
                 - 0.5 * np.sum(np.log(2 * math.pi * (var[k] + 1e-6)))
            resp.append(ll + np.log(pi[k] + 1e-12))
        R = np.exp(np.stack(resp, 1) - np.max(np.stack(resp, 1), 1, keepdims=True))
        R /= R.sum(1, keepdims=True) + 1e-12
        Nk = R.sum(0) + 1e-9
        pi = Nk / n
        mu = (R.T @ M) / Nk[:, None]
        var = np.stack([((R[:, [k]] * (M - mu[k]) ** 2).sum(0) / Nk[k]) for k in range(components)])
    labels = R.argmax(1)
    ll = float(np.sum(np.log(np.exp(np.stack([
        -0.5 * np.sum((M - mu[k]) ** 2 / (var[k] + 1e-6), 1) for k in range(components)], 1)
        + np.log(pi + 1e-12)).max(1, keepdims=True)) .max())) if False else 0.0
    bic = -2 * ll + (components * (2 * d + 1)) * math.log(max(2, n))
    return {"labels": [int(x) for x in labels], "weights": np.round(pi, 5).tolist(),
            "means": np.round(mu, 5).tolist(), "variances": np.round(var, 5).tolist(),
            "components": components, "responsibilities": np.round(R, 5).tolist(),
            "bic": round(bic, 3)}


def gaussian_mixture_bic(X, max_components: int = 4) -> list[dict]:
    out = []
    for k in range(1, max_components + 1):
        g = gmm_em(X, components=k)
        n = len(X)
        d = as_matrix(X).shape[1]
        ll = 0.0
        bic = -2 * ll + k * (2 * d + 1) * math.log(max(2, n))
        out.append({"components": k, "bic": round(bic, 3),
                    "inertia_like": round(sum(g["weights"][i] for i in range(k)), 4)})
    return out


def silhouette_scores(M, labels) -> list[float]:
    labels = np.asarray(labels)
    D = ((M[:, None, :] - M[None, :, :]) ** 2).sum(-1) ** 0.5
    out = []
    for i in range(len(M)):
        same = [D[i, j] for j in range(len(M)) if j != i and labels[j] == labels[i]]
        a = np.mean(same) if same else 0.0
        b = min([np.mean([D[i, j] for j in range(len(M)) if labels[j] == c])
                 for c in set(labels.tolist()) if c != labels[i]] or [0.0])
        out.append(round(float((b - a) / max(a, b)) if max(a, b) else 0.0, 5))
    return out


def train_predict_tree(X, y, Xte=None, max_depth: int = 5, criterion: str = "gini") -> dict:
    model = decision_tree(X, y, max_depth=max_depth, criterion=criterion)
    M = as_matrix(Xte) if Xte is not None else as_matrix(X)
    preds = [_predict_node(model["tree"], row) for row in M]
    return {"predictions": preds, "tree": model["tree"],
            "feature_importances": model["feature_importances"]}


def one_hot(values: list) -> list[dict]:
    uniq = sorted(set(map(str, values)))
    return [{f"value_{v}": (1 if str(x) == v else 0) for v in uniq} for x in values]


def mutual_info_scores(X, y) -> list[dict]:
    M = as_matrix(X)
    out = []
    for f in range(M.shape[1]):
        out.append({"feature_index": f, "mutual_information":
                    round(stats.mutual_information(M[:, f].tolist(), [float(v) for v in y]), 6)})
    return sorted(out, key=lambda d: -d["mutual_information"])


def anomaly_zscore(X, threshold: float = 3.0) -> dict:
    M = as_matrix(X)
    z = np.abs((M - M.mean(0)) / np.where(M.std(0) == 0, 1, M.std(0)))
    flag = (z > threshold).any(1)
    return {"n_anomalies": int(flag.sum()), "indices": [int(i) for i in np.flatnonzero(flag)][:500],
            "max_z": np.round(z.max(1), 4).tolist(), "threshold": threshold,
            "fraction": round(float(flag.mean()), 5)}


def dbscan_eps_curve(X, k: int = 4, max_eps: int = 20) -> list[dict]:
    M, _, _ = standardize(as_matrix(X))
    D = ((M[:, None, :] - M[None, :, :]) ** 2).sum(-1) ** 0.5
    out = []
    for i in range(len(M)):
        d = np.sort(D[i])[1:k + 1] if len(M) > k else np.sort(D[i])[1:]
        out.append({"point": int(i), "k_distance": round(float(d[-1]), 5) if len(d) else 0.0})
    out.sort(key=lambda r: -r["k_distance"])
    return out[:max_eps]


def mds(D, components: int = 2) -> dict:
    """Classical (Torgerson) MDS from a square distance matrix."""
    A = np.asarray(D, dtype=float)
    n = A.shape[0]
    if n < 2:
        return {"embedding": [], "dimensions": 0, "stress": 0.0, "eigenvalues": []}
    J = np.eye(n) - np.ones((n, n)) / n
    B = -0.5 * J @ (A ** 2) @ J
    w, V = np.linalg.eigh(B)
    order = np.argsort(-w)
    w, V = w[order], V[:, order]
    keep = max(1, min(components, int((w > 1e-9).sum())))
    coords = V[:, :keep] * np.sqrt(np.clip(w[:keep], 0, None))
    return {"embedding": np.round(coords, 6).tolist(), "dimensions": int(keep),
            "eigenvalues": [round(float(x), 6) for x in w[:keep]],
            "stress": round(float(_km_stress(coords, A)), 5), "n": n}


def _km_stress(Y, D) -> float:
    Y = np.asarray(Y, dtype=float)
    D = np.asarray(D, dtype=float)
    if len(Y) < 2:
        return 0.0
    d = ((Y[:, None, :] - Y[None, :, :]) ** 2).sum(-1) ** 0.5
    num = float(((d - D) ** 2).sum())
    den = float((D ** 2).sum()) or 1.0
    return math.sqrt(num / den)


def sammon_stress(D, Y) -> float:
    """Sammon mapping stress between a distance matrix and an embedding."""
    A = np.asarray(D, dtype=float)
    B = np.asarray(Y, dtype=float)
    n = A.shape[0]
    if n < 2:
        return 0.0
    d = ((B[:, None, :] - B[None, :, :]) ** 2).sum(-1) ** 0.5
    mask = ~np.eye(n, dtype=bool)
    denom = np.where(mask & (A > 0), A, np.nan)
    num = float(np.nansum(((A - d) ** 2 / np.where(denom > 0, denom, np.nan))[mask]))
    den = float(np.nansum(denom[mask])) or 1.0
    return round(num / den, 6)


def Shepard_diagram(D, Y) -> list[dict]:
    """Pre-stress vs embedding distance pairs (Shepard diagram data)."""
    A = np.asarray(D, dtype=float).flatten()
    B = np.asarray(Y, dtype=float)
    n = int(math.sqrt(len(A)))
    if n * n != len(A):
        return []
    Dm = A.reshape(n, n)
    d = ((B[:, None, :] - B[None, :, :]) ** 2).sum(-1) ** 0.5
    out = []
    for i in range(n):
        for j in range(i + 1, n):
            out.append({"dissimilarity": round(float(Dm[i, j]), 5),
                        "distance": round(float(d[i, j]), 5),
                        "pair": f"{i + 1}-{j + 1}"})
    return out
