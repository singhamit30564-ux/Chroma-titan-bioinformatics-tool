"""Phylogenetics: distance methods (NJ/UPGMA), tree maths and Newick surgery."""

from __future__ import annotations

import copy
import itertools
import math
import random
from typing import Sequence

from chroma_titan.core.io import as_text, newick_leaves, parse_newick, write_newick
from chroma_titan.core.seq import clean


# ---------------------------------------------------------------------------
# tree building
# ---------------------------------------------------------------------------
def neighbor_joining(names: Sequence[str], D: Sequence[Sequence[float]]) -> dict:
    """Saitou-Nei neighbour joining on a distance matrix -> nested-dict tree."""
    dist = [[float(x) for x in row] for row in D]
    nodes = list(names)
    while len(nodes) > 2:
        n = len(nodes)
        r = [sum(dist[i][j] for j in range(n) if j != i) / max(1, n - 2) for i in range(n)]
        best = None
        for i in range(n):
            for j in range(i + 1, n):
                q = (n - 2) * dist[i][j] - r[i] - r[j]
                if best is None or q < best[0]:
                    best = (q, i, j)
        _, i, j = best
        di = 0.5 * dist[i][j] + (r[i] - r[j]) / (2 * (n - 2)) if n > 2 else dist[i][j] / 2
        dj = dist[i][j] - di
        new_name = f"node{len(names) + len([1 for x in nodes if str(x).startswith('node')])}"
        merged = {"name": new_name, "length": None, "children": [
            _as_node(nodes[i], max(0.0, di)), _as_node(nodes[j], max(0.0, dj))]}
        newdist = []
        for k in range(n):
            if k in (i, j):
                continue
            newdist.append(0.5 * (dist[i][k] + dist[j][k] - dist[i][j]))
        nodes = [nodes[k] for k in range(n) if k not in (i, j)] + [merged]
        m = len(nodes)
        nd = [[0.0] * m for _ in range(m)]
        old_idx = [k for k in range(n) if k not in (i, j)]
        for a, ka in enumerate(old_idx):
            for b, kb in enumerate(old_idx):
                nd[a][b] = dist[ka][kb]
        for a in range(len(old_idx)):
            nd[a][m - 1] = nd[m - 1][a] = newdist[a]
        dist = nd
    if len(nodes) == 2:
        d = dist[0][1] if dist else 0.0
        return {"name": "root", "length": None,
                "children": [_as_node(nodes[0], d / 2), _as_node(nodes[1], d / 2)]}
    return _as_node(nodes[0], 0.0)


def _as_node(x, length) -> dict:
    if isinstance(x, dict):
        out = dict(x)
        out["length"] = length
        return out
    return {"name": str(x), "length": float(length), "children": []}


def upgma(names: Sequence[str], D: Sequence[Sequence[float]], weighted: bool = False) -> dict:
    """UPGMA / WPGMA average linkage ultrametric tree."""
    dist = [[float(x) for x in row] for row in D]
    nodes = [{"name": str(n), "length": 0.0, "children": [], "_h": 0.0, "_n": 1} for n in names]
    while len(nodes) > 1:
        n = len(nodes)
        best = None
        for i in range(n):
            for j in range(i + 1, n):
                if best is None or dist[i][j] < best[0]:
                    best = (dist[i][j], i, j)
        d, i, j = best
        hi, hj = nodes[i]["_h"], nodes[j]["_h"]
        ni, nj = nodes[i]["_n"], nodes[j]["_n"]
        h = d / 2
        parent = {"name": f"node{len(names) + 1}", "length": None, "children": [
            {k: v for k, v in nodes[i].items() if not k.startswith("_")},
            {k: v for k, v in nodes[j].items() if not k.startswith("_")}]}
        parent["children"][0]["length"] = round(max(0.0, h - hi), 8)
        parent["children"][1]["length"] = round(max(0.0, h - hj), 8)
        parent["_h"] = h
        parent["_n"] = ni + nj
        newdist = []
        for k in range(n):
            if k in (i, j):
                continue
            if weighted:
                newdist.append(0.5 * (dist[i][k] + dist[j][k]))
            else:
                newdist.append((ni * dist[i][k] + nj * dist[j][k]) / (ni + nj))
        keep = [nodes[k] for k in range(n) if k not in (i, j)] + [parent]
        m = len(keep)
        nd = [[0.0] * m for _ in range(m)]
        old = [k for k in range(n) if k not in (i, j)]
        for a, ka in enumerate(old):
            for b, kb in enumerate(old):
                nd[a][b] = dist[ka][kb]
        for a in range(len(old)):
            nd[a][m - 1] = nd[m - 1][a] = newdist[a]
        nodes = keep
        dist = nd
    root = nodes[0]
    return {k: v for k, v in root.items() if not k.startswith("_")}


def nj_tree_from_alignment(alignment: list[dict], model: str = "p") -> dict:
    from chroma_titan.core.align import distance_matrix

    names = [a["id"] for a in alignment]
    seqs = [a["seq"] for a in alignment]
    return neighbor_joining(names, distance_matrix(seqs, model))


def bootstrap_alignment(aln: list[dict], n: int = 100, seed: int = 1,
                        method: str = "nj", model: str = "p") -> list[dict]:
    from chroma_titan.core.align import distance_matrix

    rng = random.Random(seed)
    L = max((len(a["seq"]) for a in aln), default=0)
    names = [a["id"] for a in aln]
    seqs = [a["seq"] for a in aln]
    trees = []
    for _ in range(n):
        cols = [rng.randrange(L) for _ in range(L)]
        resampled = ["".join(s[c] for c in cols if c < len(s)) for s in seqs]
        D = distance_matrix(resampled, model)
        trees.append(neighbor_joining(names, D) if method == "nj" else upgma(names, D))
    return trees


def clade_frequencies(trees: list[dict]) -> dict[frozenset, int]:
    freq: dict[frozenset, int] = {}
    for t in trees:
        for cl in all_splits(t):
            freq[cl] = freq.get(cl, 0) + 1
    return dict(sorted(freq.items(), key=lambda kv: -kv[1]))


def all_splits(node: dict, universe: frozenset | None = None) -> list[frozenset]:
    if universe is None:
        universe = frozenset(newick_leaves(node))
    out: list[frozenset] = []

    def rec(n: dict) -> frozenset:
        if not n.get("children"):
            s = frozenset([n.get("name") or ""])
        else:
            s = frozenset().union(*[rec(c) for c in n["children"]])
        if s and s != universe:
            out.add(s) if isinstance(out, set) else out.append(s)
        return s

    rec(node)
    return [x for x in out if 1 < len(x) < len(universe)]


def majority_rule_consensus(trees: list[dict], cutoff: float = 0.5) -> dict:
    """Consensus tree containing clades supported by >= cutoff of the trees."""
    if not trees:
        return {"name": "consensus", "children": [], "length": None}
    leaves = newick_leaves(trees[0])
    freq = clade_frequencies(trees)
    n = len(trees)
    kept = [cl for cl, c in freq.items() if c / n >= cutoff]
    kept.sort(key=len)
    root: dict = {"name": "consensus", "length": None, "children": []}
    placed: set[str] = set()
    for cl in kept:
        node = {"name": "|".join(sorted(cl)), "length": None, "children": []}
        for child in list(root["children"]):
            if set(child["name"].split("|")) <= cl:
                node["children"].append(child)
                root["children"].remove(child)
        root["children"].append(node)
        placed |= cl
    for leaf in leaves:
        if leaf not in placed:
            root["children"].append({"name": leaf, "length": 1.0, "children": []})
    return root


def robinson_foulds(t1: dict, t2: dict) -> dict:
    a = {frozenset(s) for s in all_splits(t1)}
    b = {frozenset(s) for s in all_splits(t2)}
    norm = {min_clade(s) for s in a}, {min_clade(s) for s in b}
    inter = len(norm[0] & norm[1])
    total = max(1, len(norm[0] | norm[1]))
    return {"splits_shared": inter, "splits_total": total,
            "RF_distance": round(1 - inter / total, 5),
            "RF_normalized": round((len(norm[0] - norm[1]) + len(norm[1] - norm[0])) / total, 5)}


def min_clade(s: frozenset) -> frozenset:
    all_leaves = None
    del all_leaves
    return s


def tree_stats(node: dict, name: str = "") -> dict:
    leaves = newick_leaves(node)
    heights = {}

    def depth(n: dict, d: int) -> int:
        if not n.get("children"):
            heights[n.get("name", "")] = d
            return d
        return 1 + max(depth(c, d + 1) for c in n["children"])

    maxdepth = depth(node, 0)
    blens = _branch_lengths(node)
    coph = cophenetic_distances(node)
    offs = [len(r) for r in coph]
    return {"tree_name": name or node.get("name", "tree"), "n_tips": len(leaves),
            "n_internal": len(leaves) - 2 if len(leaves) > 1 else 0,
            "tree_length": round(sum(v for v in blens if v is not None), 6),
            "max_depth": maxdepth, "min_leaf_depth": min(heights.values()) if heights else 0,
            "ultrametric": _is_ultrametric(node, leaves),
            "n_cherries": _cherries(node), "ladderized_length_var": round(
                stats_var([len(x) for x in leaves]) or 0.0, 4),
            "colless_like": _colless(node), "mean_root_to_tip": round(
                stats_mean(heights.values()) if heights else 0.0, 4)}


def _branch_lengths(node: dict) -> list[float | None]:
    out = [node.get("length")]
    for c in node.get("children", []) or []:
        out.extend(_branch_lengths(c))
    return out


def _is_ultrametric(node: dict, leaves: list[str]) -> bool:
    ds = []

    def rec(n, acc):
        acc = acc + (n.get("length") or 0.0)
        if not n.get("children"):
            ds.append(acc)
            return
        for c in n["children"]:
            rec(c, acc)

    rec(node, 0.0)
    return bool(ds) and (max(ds) - min(ds)) < 1e-6


def _cherries(node: dict) -> int:
    n = 0
    for c in node.get("children", []) or []:
        if len(node.get("children", [])) >= 2:
            pass
    for c in node.get("children", []) or []:
        n += _cherries(c)
    kids = node.get("children") or []
    if len(kids) == 2 and all(not k.get("children") for k in kids):
        n += 1
    return n


def _colless(node: dict) -> int:
    def count_leaves(n):
        return 1 if not n.get("children") else sum(count_leaves(c) for c in n["children"])

    s = 0
    kids = node.get("children") or []
    if len(kids) == 2:
        a, b = kids
        s += abs(count_leaves(a) - count_leaves(b))
    for k in kids:
        s += _colless(k)
    return s


def stats_mean(vals) -> float:
    vals = [v for v in vals if v is not None]
    return sum(vals) / len(vals) if vals else 0.0


def stats_var(vals) -> float:
    vals = [v for v in vals if v is not None]
    if len(vals) < 2:
        return 0.0
    m = stats_mean(vals)
    return sum((x - m) ** 2 for x in vals) / (len(vals) - 1)


def _tip_pair_matrix(node: dict) -> tuple[list[list[float]], dict[str, float]]:
    """Tip-to-tip distances through the MRCA, plus root-to-tip depths."""
    leaves = newick_leaves(node)
    idx = {nm: i for i, nm in enumerate(leaves)}
    n = len(leaves)
    D = [[0.0] * n for _ in range(n)]
    depths: dict[str, float] = {}

    def rec(nd: dict, acc: float) -> list[str]:
        e = acc + float(nd.get("length") or 0.0)
        if not nd.get("children"):
            nm = str(nd.get("name") or "")
            depths[nm] = e
            return [nm]
        below: list[str] = []
        for c in nd["children"]:
            tips_c = rec(c, e)
            for t1 in tips_c:
                d1 = depths.get(t1, e)
                for t2 in below:
                    a, b = idx.get(t1, -1), idx.get(t2, -1)
                    if a >= 0 and b >= 0 and a != b:
                        d = max(0.0, d1 + depths.get(t2, e) - 2 * e)
                        D[a][b] = D[b][a] = round(d, 8)
            below.extend(tips_c)
        return below

    rec(node, 0.0)
    return D, depths


def cophenetic_distances(node: dict) -> list[list[float]]:
    """Square cophenetic (MRCA-depth) distance matrix over the tips."""
    return _tip_pair_matrix(node)[0]


def patristic_matrix(node: dict) -> list[list[float]]:
    """Tip-to-tip path lengths: on a rooted tree this is the MRCA path."""
    return _tip_pair_matrix(node)[0]


def root_to_tip_depths(node: dict) -> dict[str, float]:
    """Distance from the root to every tip."""
    return _tip_pair_matrix(node)[1]


def _adjacency(node: dict):
    """Undirected view of a nested tree: names, adjacency list and tip indices."""
    names: list[str] = []
    adj: dict[int, list[tuple[int, float]]] = {}

    def conv(nd: dict) -> int:
        me = len(names)
        names.append(str(nd.get("name") or ""))
        adj.setdefault(me, [])
        for c in nd.get("children") or []:
            cid = conv(c)
            l = float(c.get("length") or 0.0)
            adj[me].append((cid, l))
            adj.setdefault(cid, []).append((me, l))
        return me

    conv(node)
    tips = set(newick_leaves(node))
    tip_idx = {i for i, nm in enumerate(names) if nm in tips}
    return names, adj, tip_idx


def _distances_from(start: int, adj: dict[int, list[tuple[int, float]]]) -> dict[int, float]:
    out = {start: 0.0}
    stack = [start]
    while stack:
        u = stack.pop()
        for v, l in adj.get(u, []):
            if v not in out:
                out[v] = out[u] + l
                stack.append(v)
    return out


def _tree_diameter(node: dict) -> tuple[float, str, str]:
    D, _ = _tip_pair_matrix(node)
    leaves = newick_leaves(node)
    best: tuple[float, str, str] = (0.0, "", "")
    for i in range(len(leaves)):
        for j in range(i + 1, len(leaves)):
            if D[i][j] > best[0]:
                best = (D[i][j], leaves[i], leaves[j])
    return best


def _build_from(target: int, names: list[str], adj: dict[int, list[tuple[int, float]]],
                tip_idx: set[int]) -> dict:
    def build(u: int, parent: int | None) -> dict:
        kids = []
        for v, l in adj.get(u, []):
            if v == parent:
                continue
            child = build(v, u)
            child["length"] = round(l, 8)
            kids.append(child)
        out: dict = {"name": names[u] if (u in tip_idx or not kids) else "", "length": None}
        if kids:
            out["children"] = kids
        return out

    return build(target, None)


def reroot_at(node: dict, name: str) -> dict:
    """Rebuild the nested tree with the node called ``name`` as its root.

    Re-rooting *on a tip* keeps that tip as a leaf of a fresh root node, so the
    leaf count is unchanged (the branch that led to it is kept on the other side).
    """
    names, adj, tip_idx = _adjacency(node)
    target = next((i for i, nm in enumerate(names) if nm == name), None)
    if target is None:
        return node
    built = _build_from(target, names, adj, tip_idx)
    if target in tip_idx and built.get("children"):
        rest = built["children"][0]
        return {"name": "", "length": None,
                "children": [{"name": names[target], "length": 0.0}, rest]}
    return built

def ladderize(node: dict, reverse: bool = False) -> dict:
    out = copy.deepcopy(node)

    def rec(n: dict):
        if n.get("children"):
            for c in n["children"]:
                rec(c)
            n["children"].sort(key=lambda c: _leaf_count(c), reverse=reverse)

    rec(out)
    return out


def _leaf_count(n: dict) -> int:
    return 1 if not n.get("children") else sum(_leaf_count(c) for c in n["children"])


def prune_tree(node: dict, keep: Sequence[str]) -> dict:
    keep = set(keep)

    def rec(n: dict):
        if not n.get("children"):
            return dict(n) if n.get("name") in keep else None
        kids = [c for c in (rec(x) for x in n["children"]) if c]
        if not kids:
            return None
        if len(kids) == 1:
            return kids[0]
        return {**{k: v for k, v in n.items() if k != "children"}, "children": kids}

    return rec(node) or {"name": "empty", "children": []}


def subtree_by_clade(node: dict, clade_name: str) -> dict | None:
    if node.get("name") == clade_name:
        return node
    for c in node.get("children", []) or []:
        got = subtree_by_clade(c, clade_name)
        if got:
            return got
    return None


def root_at_midpoint(node: dict) -> dict:
    """Root at the node nearest the midpoint of the longest tip-to-tip path."""
    total, a, b = _tree_diameter(node)
    if total <= 0 or not a or not b:
        return node
    names, adj, tip_idx = _adjacency(node)
    ia = next((i for i, nm in enumerate(names) if nm == a), None)
    ib = next((i for i, nm in enumerate(names) if nm == b), None)
    if ia is None or ib is None:
        return node
    best, best_u = None, None
    for u in range(len(names)):
        da = _distances_from(u, adj)
        if ia not in da or ib not in da:
            continue
        imbalance = abs(da[ia] - da[ib])
        if best is None or imbalance < best:
            best, best_u = imbalance, u
    if best_u is None:
        return node
    return _build_from(best_u, names, adj, tip_idx)


def reroot_at_leaf(node: dict, leaf: str) -> dict:
    """Re-root the tree on ``leaf`` (topology preserved, branch lengths kept)."""
    return reroot_at(node, leaf)



def _find(node: dict, name: str) -> dict | None:
    if node.get("name") == name:
        return node
    for c in node.get("children", []) or []:
        got = _find(c, name)
        if got:
            return got
    return None


def add_bootstrap_values(tree: dict, freqs: dict[frozenset, int], n_trees: int) -> dict:
    out = copy.deepcopy(tree)

    def rec(n: dict):
        if n.get("children"):
            for c in n["children"]:
                rec(c)
            cl = frozenset(newick_leaves(n))
            pct = round(100 * freqs.get(cl, 0) / max(1, n_trees), 1)
            base = (n.get("name") or "").split(":")[0]
            n["name"] = f"{base}:{pct}" if base else f"clade:{pct}"
            n["support"] = pct

    rec(out)
    return out


def tree_to_nexus(trees: list[dict]) -> str:
    lines = ["#NEXUS", "begin taxa;", f"\tdimensions ntax={len(newick_leaves(trees[0]))};",
             "\ttaxlabels", "\t\t" + " ".join(newick_leaves(trees[0])) + ";", "end;",
             "begin trees;", "\ttree tree1 = " + write_newick(trees[0])]
    for i, t in enumerate(trees[1:], 2):
        lines.append(f"\ttree tree{i} = " + write_newick(t))
    lines.append("end;")
    return "\n".join(lines) + "\n"


def read_nexus(text: str) -> list[dict]:
    import re

    out = []
    for m in re.finditer(r"tree\s+\S+\s*=\s*(\[[XRN]*\])?\s*([^;\n]+);", as_text(text)):
        nw = m.group(2).strip()
        nw = re.sub(r"^\[[^\]]*\]", "", nw).replace("[and1]", "")
        try:
            out.append(parse_newick(nw))
        except Exception:  # noqa: BLE001
            continue
    return out


def distance_from_tree(tree: dict, model: str = "patristic") -> list[list[float]]:
    """Pairwise tip distances; ``cophenetic`` and ``patristic`` agree on rooted trees."""
    return _tip_pair_matrix(tree)[0]


def tree_length(tree: dict) -> float:
    return round(sum(v for v in _branch_lengths(tree) if v is not None), 6)


def tip_labels(tree: dict) -> list[str]:
    return newick_leaves(tree)


def random_tree(names: Sequence[str], seed: int = 1, mode: str = "yule") -> dict:
    """Yule/constant-birth-death birth process tree (for simulations/tests)."""
    rng = random.Random(seed)
    nodes = [{"name": n, "length": round(rng.uniform(0.01, 0.5), 4), "children": []}
             for n in names]
    while len(nodes) > 1:
        a = nodes.pop(rng.randrange(len(nodes)))
        if nodes:
            b = nodes.pop(rng.randrange(len(nodes)))
            nodes.append({"name": f"node{len(nodes) + 1}", "length": round(rng.uniform(0.01, 0.4), 4),
                          "children": [a, b]})
        else:
            nodes.append(a)
    root = nodes[0]
    root["name"] = "root"
    return root


def _binomial(n: int, p: float, rng) -> int:
    """Small-seed-friendly binomial draw (numpy Generator APIs vary)."""
    if n <= 0 or p <= 0:
        return 0
    if p >= 1:
        return n
    return sum(1 for _ in range(int(n)) if rng.random() < p)


def simulate_sequences(tree: dict, length: int = 300, rate: float = 1.0,
                       seed: int = 4, alphabet: str = "ACGT") -> dict[str, str]:
    """Simulate aligned sequences along a tree under a Jukes-Cantor process."""
    rng = random.Random(seed)
    out: dict[str, str] = {}

    def rec(node: dict, parent_seq: str) -> str:
        blen = (node.get("length") or 0.0) * rate
        seq = list(parent_seq)
        n_sub = _binomial(len(seq), min(0.95, 1 - math.exp(-0.75 * blen)), rng) if blen > 0 else 0
        for _ in range(n_sub):
            i = rng.randrange(len(seq))
            choices = [c for c in alphabet if c != seq[i]]
            seq[i] = rng.choice(choices)
        cur = "".join(seq)
        if not node.get("children"):
            out[node.get("name") or "tip"] = cur
            return cur
        return "".join(rng.choice([cur, rec(c, cur)]) for c in node["children"][:1]) + \
            "".join(rec(c, cur)[1:] for c in node["children"][1:]) if len(node["children"]) > 1 \
            else rec(node["children"][0], cur)

    root_seq = "".join(rng.choice(alphabet) for _ in range(length))
    for c in tree.get("children") or [tree]:
        rec(c, root_seq)
    return out


def jukes_cantor(d: float) -> float:
    return -0.75 * math.log(max(1e-9, 1 - 4 / 3 * d))


def kimura_two_parameter(transitions: int, transversions: int, sites: int) -> dict:
    p = transitions / sites if sites else 0.0
    q = transversions / sites if sites else 0.0
    a = 0.5 * math.log(max(1e-9, 1 - 2 * p - q))
    b = 0.25 * math.log(max(1e-9, 1 - 2 * q))
    return {"P": round(p, 6), "Q": round(q, 6), "K2P_distance": round(-(a + b), 6),
            "transition_rate": round(-0.5 * (a - b), 6), "transversion_rate": round(-b, 6),
            "kappa": round((2 * p) / q if q else float("inf"), 4)}
