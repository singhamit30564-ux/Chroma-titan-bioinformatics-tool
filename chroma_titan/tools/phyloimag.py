"""Phylogenetics, protein structure, imaging, Hi-C and geospatial / climate tools.

Panel sections covered: **phylogenetics**, **protein_modeling**, **imaging**,
**chromosome_conformation**, **gis_data_handling**, **climate_analysis**,
**compute_indicators_for_satellite_remote_sensing**.
"""
from __future__ import annotations

import math
import re
from collections import Counter, defaultdict

import numpy as np
import pandas as pd

from chroma_titan.core import align, genome, image, io, ml, phylo, plot, protein, seq, stats, tables
from chroma_titan.tools._common import *  # noqa: F401,F403
from chroma_titan.tools._common import (COUNTS, GENES_FA, IMAGE, MSA, PAIRS, PHENO, PROTEINS,
                                        REF_PROTEINS, REGIONS_BED, SHORT_PROT, TREE)

PHY = "phylogenetics"
STRUCT = "protein_modeling"
IMG = "imaging"
HIC = "chromosome_conformation"
GIS = "gis_data_handling"
CLIM = "climate_analysis"
RS = "compute_indicators_for_satellite_remote_sensing"


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------
def _fa(src):
    return io.parse_fasta(io.as_text(src))


def _aln(src):
    txt = io.as_text(src)
    recs = io.parse_fasta(txt)
    if recs:
        return [{"id": r.id, "name": r.id, "seq": r.seq, "sequence": r.seq} for r in recs]
    return align.from_phylip(txt) if hasattr(align, "from_phylip") else []


def _tbl(src):
    return tables.load(io.as_text(src))


def _mat(src):
    df = _tbl(src)
    if df.empty:
        return [], [], np.zeros((0, 0))
    cols = list(df.columns)
    body = df.drop(columns=[cols[0]])
    data = np.column_stack([pd.to_numeric(body[c], errors="coerce").to_numpy(dtype=float)
                           for c in body.columns]) if len(body.columns) else np.zeros((len(body), 0))
    return [str(x) for x in df[cols[0]].tolist()], [str(c) for c in body.columns], data


def _num_cols(df):
    out = []
    for c in df.columns:
        v = pd.to_numeric(df[c], errors="coerce").to_numpy(dtype=float)
        if np.isfinite(v).sum() >= 2:
            out.append((str(c), v))
    return out


def _read_newick(source):
    txt = io.as_text(source).strip()
    if not txt:
        return None
    return phylo.parse_newick(txt)


def _img(src):
    return image.to_gray(image.read_image(src))


# ===========================================================================
# phylogenetics
# ===========================================================================
@T("phylo_nj_from_alignment", "Neighbour-joining tree from an alignment", PHY, "text",
   [anyfile("alignment_file", MSA, fmt="", label="Aligned sequences (FASTA/Phylip)"),
    choice("model", ["p", "identity", "k2p", "jc69", "tn93", "jaccard_kmer"], "p", "Substitution model"),
    boolean("root_midpoint", True, "Midpoint root"), intin("min_length", 4, "Minimum sequence length", min=2)],
   ex={"alignment_file": MSA, "model": "p"}, up="FastTree / MEGA NJ / ape NJ",
   tags=("neighbour joining", "tree"),
   summary="Distance matrix plus neighbour-joining tree, written as Newick.")
def phylo_nj_from_alignment(alignment_file, model="p", root_midpoint=True, min_length=4):
    """NJ tree from an alignment."""
    aln = [d for d in _aln(alignment_file) if len(str(d.get("sequence", ""))) >= int(min_length)]
    if len(aln) < 3:
        return text("", "need at least three sequences")
    D = align.distance_matrix([str(d["sequence"]) for d in aln], model=model)
    tree = phylo.neighbor_joining([str(d.get("id", f"seq_{i + 1}")) for i, d in enumerate(aln)], D)
    if root_midpoint:
        try:
            tree = phylo.root_at_midpoint(tree)
        except Exception:  # noqa: BLE001
            pass
    return text(io.write_newick(tree), f"NJ tree of {len(aln)} taxa, length "
                                       f"{round(float(phylo.tree_length(tree)), 5)} ({model})")


@T("phylo_upgma", "UPGMA / WPGMA distance clustering tree", PHY, "text",
   [anyfile("alignment_file", MSA, fmt="", label="Alignment"), boolean("weighted", False, "Use WPGMA")],
   ex={"alignment_file": MSA}, up="ape::upgma / MEGA UPGMA", tags=("UPGMA", "tree"),
   summary="Ultrametric phenetic tree from an alignment distance matrix.")
def phylo_upgma(alignment_file, weighted=False):
    """UPGMA tree."""
    aln = _aln(alignment_file)
    if len(aln) < 3:
        return text("", "need at least three sequences")
    D = align.distance_matrix([str(d["sequence"]) for d in aln], model="identity")
    tree = phylo.upgma([str(d.get("id", f"seq_{i + 1}")) for i, d in enumerate(aln)], D,
                       weighted=bool(weighted))
    return text(io.write_newick(tree), f"UPGMA ({'weighted' if weighted else 'unweighted'}) tree of "
                                       f"{len(aln)} taxa")


@T("phylo_bootstrap", "Non-parametric bootstrap support for a tree", PHY, "text",
   [anyfile("alignment_file", MSA, fmt="", label="Alignment"), intin("replicates", 100, "Bootstrap replicates",
                                                                 min=5), choice("method", ["nj", "upgma"], "nj",
                                                                          "Tree method"),
    intin("seed", 1, "Seed", min=0), choice("model", ["p", "identity", "k2p"], "p", "Distance model")],
   ex={"alignment_file": MSA, "replicates": 25}, up="phangorn bootstrap / seqboot",
   tags=("bootstrap", "support", "tree"),
   summary="Resample alignment columns and annotate the consensus tree with support values.")
def phylo_bootstrap(alignment_file, replicates=100, method="nj", seed=1, model="p"):
    """Bootstrap support."""
    aln = _aln(alignment_file)
    if len(aln) < 4:
        return text("", "need at least four sequences for a meaningful bootstrap")
    trees = phylo.bootstrap_alignment(aln, n=int(replicates), seed=int(seed), method=method, model=model)
    freqs = phylo.clade_frequencies(trees) if trees else {}
    base = phylo.neighbor_joining([str(d.get("id", f"t{i + 1}")) for i, d in enumerate(aln)],
                                align.distance_matrix([str(d["sequence"]) for d in aln], model=model))
    scored = phylo.add_bootstrap_values(base, freqs, max(1, int(replicates))) if freqs else base
    consensus = phylo.majority_rule_consensus(trees, cutoff=0.5) if trees else None
    out = io.write_newick(scored)
    if consensus is not None:
        out += "\n" + io.write_newick(consensus)
    return text(out, f"{len(trees)} bootstrap trees; {len(freqs)} distinct clades, "
                    f"consensus appended as the second line")


@T("phylo_consensus_tree", "Majority-rule consensus of a tree set", PHY, "text",
   [txt("trees", "", "Newick trees, one per line"), number("cutoff", 0.5, "Minimum clade frequency", min=0.0),
    boolean("support_values", True, "Label clades with their frequency")],
   ex={"trees": "((a:0.1,b:0.2)c:0.3,d:0.4);\n((a:0.1,b:0.2)c:0.3,d:0.4);\n((a:0.1,d:0.4)c:0.3,b:0.2);"},
   up="consense / APE consensus", tags=("consensus", "trees"),
   summary="Collapse a set of trees into one topology with clade frequencies.")
def phylo_consensus_tree(trees, cutoff=0.5, support_values=True):
    """Consensus of multiple trees."""
    lines = [l.strip() for l in io.as_text(trees).splitlines() if l.strip() and ";" in l]
    parsed = []
    for l in lines:
        try:
            parsed.append(phylo.parse_newick(l))
        except Exception:  # noqa: BLE001
            continue
    if not parsed:
        return text("", "no Newick trees could be parsed")
    cons = phylo.majority_rule_consensus(parsed, cutoff=float(cutoff))
    freqs = phylo.clade_frequencies(parsed)
    body = io.write_newick(cons) if cons is not None else ""
    note = f"{len(parsed)} trees -> consensus with cutoff {cutoff}; {len(freqs)} clades observed"
    if not support_values:
        note += " (support values suppressed)"
    return text(body, note)


@T("phylo_robinson_foulds", "Robinson-Foulds distance between two trees", PHY, "table",
   [textbox("tree_a", "((seq_1,seq_2),(seq_3,seq_4));", "Newick tree A"),
    textbox("tree_b", "((seq_1,seq_3),(seq_2,seq_4));", "Newick tree B"),
    boolean("normalised", True, "Divide by the maximum")],
   ex={"tree_a": "((seq_1:0.12,seq_2:0.08)AB:0.21,(seq_3:0.17,seq_4:0.15)CD:0.09)root;",
      "tree_b": "((seq_1:0.12,seq_3:0.17)AC:0.2,(seq_2:0.08,seq_4:0.15)BD:0.1)root;"},
   up="ape::RF / treedist", tags=("tree comparison", "RF distance"),
   summary="Count the splits that differ between two topologies.")
def phylo_robinson_foulds(tree_a, tree_b, normalised=True):
    """RF distance."""
    t1, t2 = _read_newick(tree_a), _read_newick(tree_b)
    if t1 is None or t2 is None:
        return table([], "provide two Newick trees")
    d = dict(phylo.robinson_foulds(t1, t2))
    rows = [{"metric": k, "value": round(float(v), 6) if isinstance(v, float) else str(v)[:60]}
            for k, v in d.items()]
    rf = float(d.get("RF_distance", d.get("rf", 0.0)) or 0.0)
    if normalised:
        rows.append({"metric": "expected_max_pairs", "value": int(len(phylo.tip_labels(t1)) *
                                                                (len(phylo.tip_labels(t1)) - 3) // 2)})
    return table(rows, f"Robinson-Foulds distance {rf} between trees with "
                     f"{len(phylo.tip_labels(t1))} and {len(phylo.tip_labels(t2))} tips")


@T("phylo_tree_operations", "Root, ladderize, prune or subset a tree", PHY, "text",
   [anyfile("newick", TREE, fmt="", label="Newick tree"),
    choice("operation", ["midpoint_root", "reroot_leaf", "ladderize", "prune", "subtree", "strip_lengths"],
           "midpoint_root", "Operation"), textbox("leaf", "", "Rerooting leaf / clade name"),
    textbox("keep", "", "Comma list of taxa to keep (prune)"),
    boolean("report_stats", True, "Append tree statistics")],
   ex={"newick": TREE, "operation": "ladderize", "leaf": "seq_3"},
   up="ape root/collapse / ete3 prune", tags=("tree manipulation", "Newick"),
   summary="Apply a single topological operation and return the modified tree.")
def phylo_tree_operations(newick, operation="midpoint_root", leaf="", keep="", report_stats=True):
    """Tree surgery."""
    t = _read_newick(newick)
    if t is None:
        return text("", "no tree provided")
    tips = [str(x) for x in phylo.tip_labels(t)]
    if operation == "midpoint_root":
        out = phylo.root_at_midpoint(t)  # nearest node to the diameter midpoint
    elif operation == "reroot_leaf":
        out = phylo.reroot_at_leaf(t, leaf or tips[0])
    elif operation == "ladderize":
        out = phylo.ladderize(t, reverse=False)
    elif operation == "prune":
        wanted = [x.strip() for x in str(keep).split(",") if x.strip()] or tips[: max(2, len(tips) - 1)]
        out = phylo.prune_tree(t, wanted)
    elif operation == "subtree":
        out = phylo.subtree_by_clade(t, leaf or "") or t
    else:
        out = _strip_lengths(dict(t))
    body = io.write_newick(out)
    if report_stats:
        st = dict(phylo.tree_stats(out, name=operation))
        body += "\n# " + ", ".join(f"{k}={round(float(v), 5) if isinstance(v, float) else v}"
                                 for k, v in st.items())
    return text(body, f"{operation} applied ({len(tips)} input tips)")


def _strip_lengths(node):
    n = dict(node)
    n["length"] = None
    kids = n.get("children") or []
    n["children"] = [_strip_lengths(c) for c in kids]
    return n


@T("phylo_distance_from_tree", "Patristic and cophenetic distance matrices", PHY, "table",
   [anyfile("newick", TREE, fmt="", label="Newick tree"), choice("measure", ["patristic", "cophenetic"],
                                                             "patristic", "Measure"),
    boolean("heatmap", True, "Draw the matrix"), number("ultrametric_cutoff", 0.05, "Ultrametricity check",
                                                     min=0.0)],
   ex={"newick": TREE, "measure": "patristic"}, up="ape::cophenetic.phylo / diagonalize",
   tags=("distance matrix", "tree", "heatmap"),
   summary="Pairwise tip distances implied by a tree, with an ultrametricity test.")
def phylo_distance_from_tree(newick, measure="patristic", heatmap=True, ultrametric_cutoff=0.05):
    """Tree-implied distances."""
    t = _read_newick(newick)
    if t is None:
        return table([], "no tree provided")
    tips = [str(x) for x in phylo.tip_labels(t)]
    depths = dict(phylo.root_to_tip_depths(t))
    D = np.asarray(phylo.distance_from_tree(t, model=measure), dtype=float)
    if D.size == 0 or D.shape[0] != len(tips):
        D = np.asarray(phylo.patristic_matrix(t), dtype=float)
    rows = []
    for i in range(len(tips)):
        for j in range(i + 1, len(tips)):
            rows.append({"tip_a": tips[i], "tip_b": tips[j], "distance": round(float(D[i, j]), 6),
                        "measure": measure, "depth_a": round(float(depths.get(tips[i], 0.0)), 6),
                        "depth_b": round(float(depths.get(tips[j], 0.0)), 6)})
    spread = float(np.ptp(list(depths.values()))) if depths else 0.0
    res = table(rows, f"{len(tips)} tips, mean {measure} distance "
                   f"{round(float(np.mean([r['distance'] for r in rows])), 5) if rows else 0.0}; "
                   f"depth spread {round(spread, 5)}")
    res["stats"] = {"ultrametric": bool(spread <= float(ultrametric_cutoff)),
                   "labels": tips, "matrix": [[round(float(v), 5) for v in row] for row in D]}
    if heatmap:
        res["stats"]["figure"] = plot.heatmap([[round(float(v), 4) for v in row] for row in D],
                                            row_labels=tips, col_labels=tips,
                                            title=f"{measure} distances", cmap="viridis", annot=True)
    return res


@T("phylo_simulate_evolution", "Simulate sequences along a tree", PHY, "text",
   [anyfile("newick", TREE, fmt="", label="Newick tree (or leave empty for a random tree)"),
    intin("length", 300, "Alignment length", min=10), number("rate", 1.0, "Substitutions per site", min=0.0),
    intin("seed", 4, "Seed", min=0), choice("random_tree_mode", ["yule", "pda"], "yule", "Random tree model"),
    intin("ntaxa", 4, "Taxa for a random tree", min=3)],
   ex={"newick": TREE, "length": 120, "rate": 0.4}, up="evolver / Seq-Gen",
   tags=("simulation", "alignment", "evolution"),
   summary="Generate an alignment by evolving a sequence down a topology.")
def phylo_simulate_evolution(newick, length=300, rate=1.0, seed=4, random_tree_mode="yule", ntaxa=4):
    """Simulated alignment."""
    t = _read_newick(newick)
    tips = [str(x) for x in phylo.tip_labels(t)] if t else []
    if not t or len(tips) < 2:
        names = [f"taxon_{i + 1}" for i in range(max(3, int(ntaxa)))]
        t = phylo.random_tree(names, seed=int(seed), mode=random_tree_mode)
    seqs = phylo.simulate_sequences(t, length=int(length), rate=float(rate), seed=int(seed))
    recs = [io.Seq(id=k, seq=v) for k, v in sorted(dict(seqs).items())]
    body = io.write_fasta(recs)
    L = max((len(r.seq) for r in recs), default=0)
    diffs = 0
    for a in range(len(recs)):
        for b in range(a + 1, len(recs)):
            diffs += sum(1 for x, y in zip(recs[a].seq, recs[b].seq) if x != y)
    pairs = max(1, len(recs) * (len(recs) - 1) // 2)
    return text(body, f"{len(recs)} simulated sequences of {L} bp; mean pairwise differences "
                    f"{round(diffs / pairs, 2)}")


@T("phylo_tree_statistics", "Topology statistics of one or more trees", PHY, "table",
   [txt("trees", TREE, "Newick tree(s)"), boolean("splits", True, "Report split spectrum"),
    boolean("colless", True, "Report balance indices")],
   ex={"trees": TREE}, up="ape::NNI stats / treebalance", tags=("tree statistics", "balance"),
   summary="Tips, clades, tree length and a balance measure for each input tree.")
def phylo_tree_statistics(trees, splits=True, colless=True):
    """Tree stats table."""
    lines = [l.strip() for l in io.as_text(trees).splitlines() if l.strip() and ";" in l]
    rows = []
    for k, l in enumerate(lines, 1):
        try:
            t = phylo.parse_newick(l)
        except Exception:  # noqa: BLE001
            continue
        st = dict(phylo.tree_stats(t, name=f"tree_{k}"))
        tips = [str(x) for x in phylo.tip_labels(t)]
        rec = {"tree": f"tree_{k}", "tips": len(tips),
              "internal_nodes": st.get("n_internal", 0),
              "tree_length": round(float(st.get("tree_length", phylo.tree_length(t)) or 0.0), 6),
              "ultrametric": st.get("ultrametric", ""),
              "cherries": st.get("n_cherries", ""),
              "mean_root_to_tip": st.get("mean_root_to_tip", ""),
              "name": st.get("tree_name", "")}
        if splits:
            sp = phylo.all_splits(t)
            rec["splits"] = len(sp)
            rec["trivial_splits_removed"] = len(sp)
        if colless:
            rec["colless_like"] = st.get("colless_like", _balance(t))
        rows.append(rec)
    return table(rows, f"{len(rows)} tree(s) described")


def _balance(node):
    n = dict(node or {})
    kids = n.get("children") or []
    if len(kids) < 2:
        return 0
    a = len(_tips_of(kids[0]))
    b = len(_tips_of(kids[1]))
    return abs(a - b) + _balance(kids[0]) + _balance(kids[1])


def _tips_of(node):
    n = dict(node or {})
    kids = n.get("children") or []
    if not kids:
        return [str(n.get("name", ""))]
    out = []
    for c in kids:
        out.extend(_tips_of(c))
    return out


@T("phylo_pairwise_distances", "Pairwise distances between sequences", PHY, "table",
   [anyfile("alignment_file", MSA, fmt="", label="Alignment"),
    choice("model", ["p", "identity", "k2p", "jc69", "tn93", "jaccard_kmer"], "p", "Model"), boolean("matrix", True,
                                                                         "Emit the full matrix"),
    intin("head", 60, "Pair rows", min=1)],
   ex={"alignment_file": MSA, "model": "p"}, up="ape::dist.dna / phylip protdist",
   tags=("distance matrix", "models"),
   summary="Uncorrected and model-corrected pairwise distances for every pair.")
def phylo_pairwise_distances(alignment_file, model="p", matrix=True, head=60):
    """Distance table."""
    aln = _aln(alignment_file)
    if len(aln) < 2:
        return table([], "need at least two sequences")
    seqs = [str(d["sequence"]) for d in aln]
    ids = [str(d.get("id", f"seq_{i + 1}")) for i, d in enumerate(aln)]
    D = np.asarray(align.distance_matrix(seqs, model=model), dtype=float)
    rows = []
    for i in range(len(ids)):
        for j in range(i + 1, len(ids)):
            rows.append({"sequence_a": ids[i], "sequence_b": ids[j], "distance": round(float(D[i, j]), 6),
                        "model": model, "sites_compared": min(len(seqs[i]), len(seqs[j])),
                        "percent_difference": round(100.0 * float(D[i, j]), 4)})
    rows.sort(key=lambda r: r["distance"])
    res = table(rows[: int(head)], f"{len(rows)} pairs, mean distance "
                                f"{round(float(np.mean([r['distance'] for r in rows])), 5) if rows else 0.0}"
                                f" with the {model} model")
    if matrix:
        res["stats"] = {"labels": ids, "matrix": [[round(float(v), 6) for v in row] for row in D],
                       "figure": plot.heatmap([[round(float(v), 4) for v in row] for row in D],
                                            row_labels=ids, col_labels=ids, title=f"{model} distances",
                                            cmap="viridis", annot=True)}
    return res


@T("phylo_read_nexus", "Read a multi-tree NEXUS file", PHY, "table",
   [txt("nexus", "#NEXUS\nbegin trees;\ntree t1 = ((a:0.1,b:0.2)c:0.3,d:0.4);\ntree t2 = ((a:0.1,b:0.3)c:0.2,d:0.4);\nend;",
        "NEXUS text"), boolean("consensus", True, "Also build a consensus")],
   ex={"nexus": "#NEXUS\nbegin trees;\ntree t1 = ((seq_1:0.12,seq_2:0.08)AB:0.21,(seq_3:0.17,seq_4:0.15)CD:0.09)root;\n"
                "tree t2 = ((seq_1:0.12,seq_2:0.09)AB:0.2,(seq_3:0.16,seq_4:0.15)CD:0.1)root;\nend;"},
   up="phytools read.nexus / DendroPy", tags=("NEXUS", "trees"),
   summary="Parse the trees stored in a NEXUS block and summarise them.")
def phylo_read_nexus(nexus, consensus=True):
    """NEXUS tree reader."""
    trees = phylo.read_nexus(io.as_text(nexus))
    if not trees:
        return table([], "no trees found in the NEXUS text")
    rows = []
    for i, t in enumerate(trees, 1):
        st = dict(phylo.tree_stats(t, name=f"tree_{i}"))
        rows.append({"tree": i, "tips": len(phylo.tip_labels(t)),
                    "length": round(float(phylo.tree_length(t)), 6),
                    "name": st.get("name", ""), "nodes": st.get("nodes", ""),
                    "rooted": st.get("rooted", "")})
    res = table(rows, f"{len(trees)} trees parsed")
    if consensus and len(trees) > 1:
        cons = phylo.majority_rule_consensus(list(trees), cutoff=0.5)
        res["stats"]["consensus_newick"] = io.write_newick(cons) if cons is not None else ""
    res["stats"]["nexus_round_trip"] = len(phylo.tree_to_nexus(list(trees)).splitlines())
    return res


@T("phylo_dating_check", "Root-to-tip distances against sampling order", PHY, "table",
   [anyfile("newick", TREE, fmt="", label="Newick tree (ultrametric or not)"),
    tbl("dates", PHENO, "Table with taxon and date columns"), textbox("taxon_column", "", "Taxon column"),
    textbox("date_column", "", "Date / order column"), choice("regression", ["pearson", "spearman"], "spearman",
                                                             "Correlation")],
   ex={"newick": TREE, "dates": PHENO, "taxon_column": "sample", "date_column": "yield"},
   up="TempEst / liTo root-to-tip regression", tags=("molecular clock", "dating"),
   summary="Test whether root-to-tip divergence correlates with sampling time.")
def phylo_dating_check(newick, dates, taxon_column="", date_column="", regression="spearman"):
    """Root-to-tip regression."""
    t = _read_newick(newick)
    if t is None:
        return table([], "no tree provided")
    depths = dict(phylo.root_to_tip_depths(t))
    df = _tbl(dates)
    cols = list(df.columns)
    tc = taxon_column if taxon_column in cols else cols[0]
    dc = date_column if date_column in cols else next((c for c, _v in _num_cols(df)), cols[-1])
    pairs = []
    for r in df.to_dict("records"):
        nm = str(r[tc])
        key = next((k for k in depths if k == nm or k.endswith(nm) or nm.endswith(str(k))), None)
        try:
            dv = float(r[dc])
        except (TypeError, ValueError):
            continue
        if key is not None and math.isfinite(dv):
            pairs.append((dv, depths[key], key, nm))
    if len(pairs) < 3:
        return table([{"tip": k, "root_to_tip": round(float(v), 6), "date": "unmatched"}
                     for k, v in sorted(depths.items())],
                    "fewer than three tips match the date table - tip depths only")
    x = np.array([p[0] for p in pairs], dtype=float)
    y = np.array([p[1] for p in pairs], dtype=float)
    st = dict((stats.spearman if regression == "spearman" else stats.pearson)(list(x), list(y)))
    fit = dict(stats.linear_regression(list(x), list(y)))
    rows = [{"tip": nm, "date": round(d, 4), "root_to_tip": round(v, 6),
            "residual": round(float(v - (float(fit.get("slope", 0.0)) * d + float(fit.get("intercept", 0.0)))), 6),
            "clock_like": True} for d, v, _k, nm in pairs]
    res = table(rows, f"root-to-tip vs date: r = {round(float(st.get('r', 0.0)), 4)}, "
                    f"p = {round(float(st.get('p_value', 1.0)), 5)}, slope {round(float(fit.get('slope', 0.0)), 6)}")
    res["stats"] = {"correlation": {k: (round(float(v), 6) if isinstance(v, float) else v) for k, v in st.items()},
                   "figure": plot.scatter(list(x), list(y), labels=[p[3] for p in pairs],
                                        xlabel="date / order", ylabel="root-to-tip distance",
                                        title="clock test")}
    return res


# ===========================================================================
# protein structure and properties
# ===========================================================================
@T("protein_domain_sites", "Prosite-style motif and domain sites", STRUCT, "table",
   [fa("protein_file", PROTEINS, "Protein FASTA"), intin("head", 400, "Rows", min=1)],
   ex={"protein_file": PROTEINS}, up="InterProScan / PaSE", tags=("domains", "motifs"),
   summary="Conserved pattern matches along each chain with their positions.")
def protein_domain_sites(protein_file, head=400):
    """Motif sites."""
    rows = []
    for rec in _fa(protein_file):
        s = protein.clean(rec.seq)
        for d in protein.prosite_sites(s):
            dd = dict(d)
            rows.append({"protein": rec.id, "motif": dd.get("motif", ""), "pattern": dd.get("pattern", ""),
                        "start": dd.get("start", 0), "end": dd.get("end", 0),
                        "match": str(dd.get("match", ""))[:30], "length": len(s),
                        "relative_position": round(float(dd.get("start", 0)) / max(1, len(s)), 4)})
    return table(rows[: int(head)], f"{len(rows)} pattern hits in {len(_fa(protein_file))} chains")


@T("protein_secondary_structure", "Predicted secondary structure", STRUCT, "table",
   [fa("protein_file", PROTEINS, "Protein FASTA"), intin("window", 9, "Hydropathy window", min=3),
    number("helix_cutoff", 0.6, "Helix propensity cut-off", min=0.0),
    boolean("profile", True, "Attach a hydrophobicity profile figure")],
   ex={"protein_file": PROTEINS}, up="PSIPRED / Jpred", tags=("secondary structure", "prediction"),
   summary="Helix, strand and coil assignment from propensities and hydropathy.")
def protein_secondary_structure(protein_file, window=9, helix_cutoff=0.6, profile=True):
    """Secondary structure prediction."""
    rows = []
    for rec in _fa(protein_file):
        s = protein.clean(rec.seq)
        if not s:
            continue
        ss = protein.secondary_structure(s)
        ss = str(ss) if not isinstance(ss, str) else ss
        counts = Counter(ss)
        prof = np.asarray(protein.hydrophobicity_profile(s, window=max(3, int(window))), dtype=float)
        rows.append({"protein": rec.id, "length": len(s), "helices": ss.count("H"), "strands": ss.count("E"),
                    "coils": ss.count("C"), "helix_fraction": round(ss.count("H") / max(1, len(ss)), 4),
                    "strand_fraction": round(ss.count("E") / max(1, len(ss)), 4),
                    "mean_hydropathy": round(float(np.mean(prof)), 4),
                    "helix_propensity_over_cutoff": int(np.sum(prof > float(helix_cutoff))),
                    "assignment": ss[:60] + ("..." if len(ss) > 60 else "")})
    res = table(rows, f"{len(rows)} chains folded by propensity")
    if profile and rows:
        first = _fa(protein_file)[0]
        prof = protein.hydrophobicity_profile(protein.clean(first.seq), window=max(3, int(window)))
        res["stats"]["figure"] = plot.line({first.id: [round(float(v), 4) for v in prof]},
                                          xlabel="residue", ylabel="hydropathy",
                                          title="hydropathy profile")
    return res


@T("protein_transmembrane", "Transmembrane helix prediction", STRUCT, "table",
   [fa("protein_file", PROTEINS, "Protein FASTA"), intin("window", 19, "Sliding window", min=7),
    number("cutoff", 1.6, "Hydropathy cut-off", min=0.0), boolean("topology", True, "Report inside/outside")],
   ex={"protein_file": PROTEINS, "cutoff": 1.0}, up="TMHMM / Phobius", tags=("membrane", "topology"),
   summary="Predicted TM segments with length, centred hydropathy and orientation.")
def protein_transmembrane(protein_file, window=19, cutoff=1.6, topology=True):
    """TM helices."""
    rows = []
    for rec in _fa(protein_file):
        s = protein.clean(rec.seq)
        hel = protein.transmembrane_helices(s, window=max(5, int(window)), cutoff=float(cutoff)) if s else []
        for i, h in enumerate(hel):
            d = dict(h)
            a = int(d.get("start", 0)) + 1
            b = int(d.get("end", 0))
            rows.append({"protein": rec.id, "helix": i + 1, "start": a, "end": b, "length": max(0, b - a + 1),
                        "mean_hydropathy": round(float(d.get("mean_hydrophobicity",
                                                          d.get("score", 0.0)) or 0.0), 4),
                        "sequence": str(d.get("sequence", s[a - 1:b]))[:24],
                        "orientation": (d.get("orientation") or ("outside" if i % 2 == 0 else "inside"))
                        if topology else ""})
        if not hel and topology:
            rows.append({"protein": rec.id, "helix": 0, "start": 0, "end": 0, "length": 0,
                        "mean_hydropathy": round(float(protein.gravy(s)), 4) if s else 0.0,
                        "sequence": "", "orientation": "soluble"})
    n = len({r["protein"] for r in rows if r["helix"]})
    return table(rows, f"{n} membrane proteins with {sum(1 for r in rows if r['helix'])} predicted helices")


@T("protein_disorder", "Disorder propensity profile", STRUCT, "figure",
   [fa("protein_file", PROTEINS, "Protein FASTA"), number("cutoff", 0.5, "Disorder cut-off", min=0.0, max=1.0),
    intin("smoothing", 5, "Smoothing window", min=1), textbox("protein_id", "", "Only this protein")],
   ex={"protein_file": PROTEINS, "cutoff": 0.45}, up="IUPred / Disopred", tags=("disorder", "prediction"),
   summary="Per-residue disorder probability with disordered region calls.")
def protein_disorder(protein_file, cutoff=0.5, smoothing=5, protein_id=""):
    """Disorder profile."""
    recs = [r for r in _fa(protein_file) if not protein_id or r.id == protein_id]
    if not recs:
        return plot.empty_plot("no sequences")
    series = {}
    for r in recs[:6]:
        s = protein.clean(r.seq)
        d = np.asarray(protein.disorder_propensity(s), dtype=float)
        if not d.size:
            continue
        w = max(1, int(smoothing))
        pad = w // 2
        sm = np.convolve(np.pad(d, pad, mode="edge"), np.ones(w) / w, mode="valid")[: d.size]
        series[r.id] = [round(float(v), 4) for v in sm]
    if not series:
        return plot.empty_plot("no residues to score")
    return plot.line(series, xlabel="residue", ylabel=f"disorder probability (cutoff {cutoff})",
                    title="disorder propensity")


@T("protein_charge_properties", "Charge, pI and buffer behaviour", STRUCT, "table",
   [fa("protein_file", PROTEINS, "Protein FASTA"), number("pH", 7.0, "Working pH", min=0.0),
    number("ionic_strength", 0.15, "Ionic strength (M)", min=0.0), boolean("profile", True,
                                                                    "Charge vs pH curve")],
   ex={"protein_file": PROTEINS, "pH": 6.5}, up="IPC / EMBOSS charge", tags=("charge", "pI", "purification"),
   summary="Net charge and pI per chain, plus the pH window where it buffers best.")
def protein_charge_properties(protein_file, pH=7.0, ionic_strength=0.15, profile=True):
    """Charge table."""
    rows = []
    curves = {}
    phs: list[float] = []
    for rec in _fa(protein_file):
        s = protein.clean(rec.seq)
        if not s:
            continue
        hist = protein.charge_histogram(s, pH_range=(max(0.5, float(pH) - 3.0), min(13.5, float(pH) + 3.0)),
                                      steps=13)
        curves[rec.id] = [round(float(dict(d).get("net_charge", 0.0)), 3) for d in hist]
        phs = [float(dict(d).get("pH", 0.0)) for d in hist]
        rows.append({"protein": rec.id, "length": len(s),
                    "net_charge": round(float(protein.net_charge(s, pH=float(pH))), 3),
                    "charge_density": round(float(protein.charge_density(s)), 5),
                    "isoelectric_point": round(float(protein.isoelectric_point(s)), 3),
                    "buffer_capacity": round(float(protein.buffer_capacity_estimate(s)), 4),
                    "acidic": s.count("D") + s.count("E"), "basic": s.count("K") + s.count("R") + s.count("H"),
                    "at_pH": float(pH), "ionic_strength_M": float(ionic_strength)})
    res = table(rows, f"{len(rows)} chains at pH {pH}")
    if profile and curves and phs:
        res["stats"]["figure"] = plot.line(curves, x=[round(float(v), 2) for v in phs],
                                         xlabel="pH", ylabel="net charge", title="charge vs pH")
    return res


@T("protein_composition_table", "Amino-acid composition and physicochemistry", STRUCT, "table",
   [fa("protein_file", PROTEINS, "Protein FASTA"), choice("normalisation", ["percent", "count", "per_1000"],
                                                     "percent", "Normalisation"),
    boolean("properties", True, "Add instability / GRAVY / aromaticity")],
   ex={"protein_file": PROTEINS}, up="ProtComp / pepstats composition", tags=("composition", "statistics"),
   summary="Composition of every residue class with global physicochemical indices.")
def protein_composition_table(protein_file, normalisation="percent", properties=True):
    """Composition matrix."""
    recs = _fa(protein_file)
    aas = "ACDEFGHIKLMNPQRSTVWY"
    rows = []
    for rec in recs:
        s = protein.clean(rec.seq)
        comp = dict(protein.aa_composition(s)) if s else {}
        n = max(1, len(s))
        rec_out = {"protein": rec.id, "length": len(s)}
        for a in aas:
            v = float(comp.get(a, 0.0))
            rec_out[a] = round(v, 4) if normalisation == "percent" else (
                round(1000.0 * v / 100.0, 2) if normalisation == "per_1000" else round(n * v / 100.0, 0))
        if properties:
            rec_out["instability"] = round(float(protein.instability_index(s)), 3)
            rec_out["gravy"] = round(float(protein.gravy(s)), 4)
            rec_out["aromaticity"] = round(float(protein.aromaticity(s)), 4)
            rec_out["aliphatic_index"] = round(float(protein.aliphatic_index(s)), 3)
            rec_out["cysteines"] = s.count("C")
        rows.append(rec_out)
    tot = Counter()
    for rec in recs:
        tot.update(protein.clean(rec.seq))
    res = table(rows, f"{len(rows)} proteins; most common residue: "
                   f"{tot.most_common(1)[0][0] if tot else 'n/a'}")
    res["stats"] = {"overall_counts": dict(tot), "normalisation": normalisation,
                   "figure": plot.bar(aas, [round(100.0 * tot[a] / max(1, sum(tot.values())), 2) for a in aas],
                                    xlabel="residue", ylabel="percent", title="overall composition")}
    return res


@T("protein_structure_summary_figure", "Structure-feature summary figure", STRUCT, "figure",
   [fa("protein_file", PROTEINS, "Protein FASTA"), intin("window", 15, "Window", min=3),
    number("tm_cutoff", 1.4, "TM hydropathy cut-off", min=0.0),
    choice("track", ["hydropathy", "charge", "disorder"], "hydropathy", "Track to plot"),
    textbox("protein_id", "", "Restrict to one protein")],
   ex={"protein_file": PROTEINS, "window": 9, "track": "charge"},
   up="UCSF ChimeraX / PyMOL feature view", tags=("structure", "plot", "features"),
   summary="Windowed feature track along a chain, for a quick structural overview.")
def protein_structure_summary_figure(protein_file, window=15, tm_cutoff=1.4, track="hydropathy", protein_id=""):
    """Feature track figure."""
    recs = [r for r in _fa(protein_file) if not protein_id or r.id == protein_id]
    if not recs:
        return plot.empty_plot("no sequences")
    series = {}
    for r in recs[:6]:
        s = protein.clean(r.seq)
        w = max(3, int(window))
        if track == "charge":
            vals = np.asarray(protein.net_charge_profile(s, window=w), dtype=float)
        elif track == "disorder":
            vals = np.asarray(protein.disorder_propensity(s), dtype=float)
        else:
            vals = np.asarray(protein.hydrophobicity_profile(s, window=w), dtype=float)
        if vals.size:
            series[r.id] = [round(float(v), 4) for v in vals]
    if not series:
        return plot.empty_plot("no residues to plot")
    return plot.line(series, xlabel="residue", ylabel=track,
                    title=f"{track} track (window {window}, TM cut-off {tm_cutoff})")


@T("protein_signal_peptide", "Signal peptide and targeting prediction", STRUCT, "table",
   [fa("protein_file", PROTEINS, "Protein FASTA"), number("cutoff", 0.5, "Decision cut-off", min=0.0, max=1.0),
    intin("n_window", 15, "N-region window", min=5)],
   ex={"protein_file": PROTEINS}, up="SignalP / Phobius", tags=("signal peptide", "localisation"),
   summary="N-region hydrophobicity and cleavage-site score for secretion prediction.")
def protein_signal_peptide(protein_file, cutoff=0.5, n_window=15):
    """Signal peptide scan."""
    rows = []
    for rec in _fa(protein_file):
        s = protein.clean(rec.seq)
        d = dict(protein.signal_peptide_score(s)) if s else {}
        score = float(d.get("score", d.get("probability", 0.0)) or 0.0)
        rows.append({"protein": rec.id, "score": round(score, 4),
                    "predicted": bool(score >= float(cutoff)),
                    "cleavage_position": d.get("cleavage_site", d.get("cut", "")),
                    "n_region": s[: max(5, int(n_window))],
                    "n_region_hydropathy": round(float(protein.gravy(s[: max(5, int(n_window))])), 4) if s else 0.0,
                    "length": len(s)})
    n = sum(1 for r in rows if r["predicted"])
    return table(rows, f"{n}/{len(rows)} proteins carry a predicted signal peptide at cutoff {cutoff}")


@T("protein_mass_and_formula", "Mass, formula and extinction coefficients", STRUCT, "table",
   [fa("protein_file", PROTEINS, "Protein FASTA"), choice("mode", ["average", "mono"], "average", "Mass mode"),
    number("pH", 8.0, "pH for extinction", min=0.0), boolean("reduced", True, "Assume reduced cysteines")],
   ex={"protein_file": PROTEINS}, up="peptide mass / expasy compute pI", tags=("mass", "formula"),
   summary="Molecular weight, empirical formula, extinction and molar absorptivity.")
def protein_mass_and_formula(protein_file, mode="average", pH=8.0, reduced=True):
    """Mass table."""
    rows = []
    for rec in _fa(protein_file):
        s = protein.clean(rec.seq)
        if not s:
            continue
        ec = dict(protein.extinction_coefficient(s, pH=float(pH)))
        mass = float(protein.protein_mass(s, mode=mode))
        rows.append({"protein": rec.id, "length": len(s), "formula": protein.molecular_formula(s),
                    "molecular_weight": round(mass, 4), "mode": mode,
                    "extinction_reduced": round(float(ec.get("extinction_reduced", 0.0)), 2),
                    "extinction_oxidised": round(float(ec.get("extinction_oxidised", 0.0)), 2),
                    "E1_percent_1cm": round(float(ec.get("E1percent_1cm_reduced", 0.0)), 4),
                    "A1M_1cm": round(float(ec.get("A1M_1cm_reduced", 0.0)), 2),
                    "aromatic_residues": ec.get("aromatic_residues", s.count("W") + s.count("Y") + s.count("C")),
                    "assumed_state": "reduced" if reduced else "oxidised",
                    "moles_per_l_for_A1": round(float(ec.get("estimated_concentration_uM_for_A1", 0.0)), 3)})
    return table(rows, f"{len(rows)} chains in {mode} mass mode")


@T("protein_align_identities", "Pairwise protein identity matrix", STRUCT, "table",
   [fa("protein_file", PROTEINS, "Query proteins"), fa("reference_file", REF_PROTEINS, "Reference proteins"),
    choice("alignment", ["pairwise", "hamming"], "pairwise", "Method"), number("min_identity", 0.0,
                                                                       "Report pairs above", min=0.0),
    intin("head", 100, "Rows", min=1)],
   ex={"protein_file": PROTEINS, "reference_file": PROTEINS, "min_identity": 20.0},
   up="mmseqs easy-search / blastp identity", tags=("identity", "homology"),
   summary="Identity of every query/reference pair using scoring or k-mer overlap.")
def protein_align_identities(protein_file, reference_file, alignment="pairwise", min_identity=0.0, head=100):
    """Identity matrix."""
    queries = _fa(protein_file)
    refs = _fa(reference_file)
    rows = []
    for q in queries:
        qs = protein.clean(q.seq)
        for r in refs:
            rs = protein.clean(r.seq)
            if not qs or not rs:
                continue
            if alignment == "hamming":
                n = min(len(qs), len(rs))
                ident = 100.0 * sum(1 for a, b in zip(qs[:n], rs[:n]) if a == b) / max(1, n)
            else:
                ident = float(protein.protein_alignment_identity(qs, rs)) if hasattr(protein,
                                                                              "protein_alignment_identity") else 0.0
            if ident < float(min_identity):
                continue
            rows.append({"query": q.id, "reference": r.id, "identity_percent": round(float(ident), 3),
                        "query_length": len(qs), "reference_length": len(rs),
                        "method": alignment, "ortholog_call": bool(ident >= 70.0)})
    rows.sort(key=lambda r: -r["identity_percent"])
    return table(rows[: int(head)], f"{len(rows)} pairs above {min_identity}% identity")


@T("protein_hydrophobicity_moment", "Hydrophobic moment of helical segments", STRUCT, "table",
   [fa("protein_file", PROTEINS, "Protein FASTA"), intin("helix_period", 3.6, "Residues per turn", min=2.0),
    number("window", 18, "Window size", min=6), intin("head", 200, "Rows", min=1)],
   ex={"protein_file": PROTEINS, "window": 12}, up="WEB-based helical wheel / EMBOSS heel",
   tags=("hydrophobic moment", "helix"),
   summary="Amphipathicity of each window, the classic test for facial helices.")
def protein_hydrophobicity_moment(protein_file, helix_period=3.6, window=18, head=200):
    """Hydrophobic moment table."""
    rows = []
    for rec in _fa(protein_file):
        s = protein.clean(rec.seq)
        if not s:
            continue
        w = max(6, int(window))
        mom = np.asarray(protein.hydrophobicity_moment(s, window=w, helix=float(helix_period)), dtype=float)
        for i, v in enumerate(mom):
            rows.append({"protein": rec.id, "window_start": i + 1, "window_end": min(len(s), i + w),
                        "moment": round(float(v), 4), "mean_hydropathy":
                        round(float(protein.gravy(s[i:i + w])), 4),
                        "amphipathic": bool(float(v) > 0.45), "helical_turns": round(w / float(helix_period), 2)})
    rows.sort(key=lambda r: -r["moment"])
    res = table(rows[: int(head)], f"{len(rows)} windows; most amphipathic: "
                                f"{rows[0]['protein']} @{rows[0]['window_start']}" if rows else "no windows")
    if rows:
        res["stats"]["figure"] = plot.histogram([r["moment"] for r in rows], bins=25,
                                              xlabel="hydrophobic moment", title="moment distribution")
    return res


# ===========================================================================
# imaging
# ===========================================================================
@T("image_quality_summary", "Basic image statistics and histogram", IMG, "table",
   [img("image_file", IMAGE, "Image (PGM/PPM)"), intin("bins", 16, "Histogram bins", min=2),
    boolean("per_channel", True, "Split colour channels")],
   ex={"image_file": IMAGE, "bins": 12}, up="Fiji / ImageJ summary stats", tags=("statistics", "histogram"),
   summary="Mean, dynamic range, saturation and the intensity histogram of an image.")
def image_quality_summary(image_file, bins=16, per_channel=True):
    """Image summary."""
    arr = np.asarray(image.read_image(image_file), dtype=float)
    st = dict(image.image_stats(_img(image_file)))
    gray = image.to_gray(arr)
    hist = image.histogram(gray, bins=int(bins))
    rows = [{"metric": k, "value": round(float(v), 5) if isinstance(v, float) else str(v)[:40]}
            for k, v in st.items()]
    rows.append({"metric": "shape", "value": "x".join(str(x) for x in arr.shape)})
    rows.append({"metric": "saturated_fraction", "value": round(float(np.mean(gray >= np.max(gray))
                                                                if gray.size else 0.0), 5)})
    res = table(rows, f"image {arr.shape}: mean {round(float(np.nanmean(gray)), 2)}")
    res["stats"]["histogram"] = [{"bin": f"{d['bin_min']}-{d['bin_max']}", "count": d["count"]} for d in hist]
    res["stats"]["figure"] = plot.histogram([float(v) for v in gray.flatten()][:4000], bins=int(bins),
                                          xlabel="intensity", title="pixel intensities")
    if per_channel and arr.ndim == 3:
        res["stats"]["channels"] = [{"channel": i, "mean": round(float(np.nanmean(arr[..., i])), 4)}
                                  for i in range(arr.shape[-1])]
    return res


@T("image_threshold_cells", "Binarise and count objects", IMG, "table",
   [img("image_file", IMAGE, "Image"), choice("method", ["otsu", "manual", "background"], "otsu", "Threshold method"),
    number("value", 0.5, "Manual threshold (0-1 fraction)", min=0.0), boolean("invert", False, "Invert"),
    intin("min_area", 5, "Minimum object area", min=1), boolean("show", True, "Show the mask")],
   ex={"image_file": IMAGE, "min_area": 6}, up="Cellpose / ilastik all-pixels",
   tags=("segmentation", "threshold", "counting"),
   summary="Global threshold, morphological cleanup and an object count.")
def image_threshold_cells(image_file, method="otsu", value=0.5, invert=False, min_area=5, show=True):
    """Threshold + count."""
    arr = _img(image_file)
    if invert:
        arr = image.invert(arr)
    thr = float(value) if method == "manual" else float(dict(image.threshold_otsu(arr)).get("threshold", 0.5))
    if method == "background":
        thr = float(np.mean(arr) + np.std(arr))
    mask = image.binarize(image.gaussian_blur(arr, sigma=1.0), value=thr, above=True)
    cleaned = image.remove_small_objects(image.binary_erosion(image.binary_dilation(mask, iterations=1),
                                                            iterations=1), min_area=int(min_area))
    lab = image.label(cleaned, connectivity=2)
    labeled = lab["labels"]
    props = image.regionprops(labeled, arr)
    rows = [{"object": int(dict(p).get("label", i + 1)), "area": int(dict(p).get("area", 0)),
            "centroid_x": round(float(dict(p).get("centroid_x", 0)), 2),
            "centroid_y": round(float(dict(p).get("centroid_y", 0)), 2),
            "mean_intensity": round(float(dict(p).get("mean_intensity", 0.0)), 4),
            "solidity": round(float(dict(p).get("solidity", 0.0)), 4),
            "diameter": round(float(dict(p).get("equivalent_diameter", 0.0)), 3)}
           for i, p in enumerate(props)]
    res = table(rows, f"{len(rows)} objects at threshold {round(thr, 2)} "
                   f"(background fraction {round(float(np.mean(np.asarray(cleaned) <= 0)), 4)})")
    if show:
        res["stats"]["figure"] = plot.image_show(image.add_labels_overlay(labeled), title="labelled objects")
    res["stats"]["threshold"] = round(thr, 4)
    return res


@T("image_measure_objects", "Morphometry of segmented objects", IMG, "table",
   [img("image_file", IMAGE, "Image"), number("threshold", 0.0, "Threshold (0 = Otsu)", min=0.0),
    intin("min_area", 8, "Minimum area", min=1), choice("shape_metric", ["solidity", "eccentricity",
                                                              "extent"], "solidity", "Shape metric"),
    intin("head", 100, "Rows", min=1)],
   ex={"image_file": IMAGE, "min_area": 10, "shape_metric": "eccentricity"},
   up="CellProfiler MeasureObjectSizeShape", tags=("morphometry", "shape"),
   summary="Area, perimeter, eccentricity and solidity for each detected object.")
def image_measure_objects(image_file, threshold=0.0, min_area=8, shape_metric="solidity", head=100):
    """Morphometric measurements."""
    arr = _img(image_file)
    thr = float(threshold) if float(threshold) > 0 else float(dict(image.threshold_otsu(arr)).get("threshold", 0.5))
    mask = image.binarize(arr, value=thr, above=True)
    lab = image.label(image.remove_small_objects(mask, min_area=int(min_area)), connectivity=2)
    props = image.regionprops(lab["labels"], arr)
    rows = []
    for p in props:
        d = dict(p)
        rows.append({"object": int(d.get("label", len(rows) + 1)), "area": int(d.get("area", 0)),
                    "perimeter": round(float(d.get("perimeter_estimate", 0.0)), 3),
                    "major_axis": round(float(d.get("major_axis_length", 0.0)), 3),
                    "minor_axis": round(float(d.get("minor_axis_length", 0.0)), 3),
                    shape_metric: round(float(d.get(shape_metric, 0.0)), 4),
                    "mean_intensity": round(float(d.get("mean_intensity", 0.0)), 4),
                    "bbox": f"{int(d.get('bbox_xmin', 0))},{int(d.get('bbox_ymin', 0))}-"
                            f"{int(d.get('bbox_xmax', 0))},{int(d.get('bbox_ymax', 0))}",
                    "diameter": round(float(d.get("equivalent_diameter", 0.0)), 3)})
    res = table(rows[: int(head)], f"{len(rows)} objects measured (min area {min_area})")
    if rows:
        res["stats"]["figure"] = plot.histogram([float(r["area"]) for r in rows], bins=15,
                                              xlabel="area (px)", title="object areas")
    return res


@T("image_texture_features", "Texture descriptors (GLCM and coarseness)", IMG, "table",
   [img("image_file", IMAGE, "Image"), intin("levels", 8, "Grey levels", min=2),
    intin("distance", 1, "GLCM offset", min=1), choice("angle", ["0.0", "0.7854", "1.5708"], "0.0",
                                                     "Direction (radians)"),
    boolean("multi_scale", True, "Multi-scale statistics")],
   ex={"image_file": IMAGE, "levels": 4}, up="Ilastik texture / skimage features",
   tags=("texture", "GLCM"),
   summary="Contrast, homogeneity, energy and entropy at one or several scales.")
def image_texture_features(image_file, levels=8, distance=1, angle="0.0", multi_scale=True):
    """Texture table."""
    arr = _img(image_file)
    gl = dict(image.glcm(arr, distance=int(distance), angle=float(angle), levels=int(levels))[0])
    rows = [{"metric": "glcm_ASM", "value": round(float(gl.get("asm", 0.0)), 5), "scale": 1},
           {"metric": "glcm_contrast", "value": round(float(gl.get("contrast", 0.0)), 5), "scale": 1},
           {"metric": "glcm_correlation", "value": round(float(gl.get("correlation", 0.0)), 5), "scale": 1},
           {"metric": "glcm_homogeneity", "value": round(float(gl.get("homogeneity", 0.0)), 5), "scale": 1},
           {"metric": "glcm_energy", "value": round(float(gl.get("energy", 0.0)), 5), "scale": 1}]
    if multi_scale:
        for lv in (2, 4, 8, 16):
            for d in image.texture_stats(arr, levels=int(lv)):
                dd = dict(d)
                rows.append({"metric": f"{dd.get('level', lv)}_entropy",
                            "value": round(float(dd.get("entropy", 0.0)), 5), "scale": int(lv),
                            "mean": round(float(dd.get("mean", 0.0)), 4),
                            "std": round(float(dd.get("std", 0.0)), 4), "size": dd.get("size", "")})
    res = table(rows, f"texture of {arr.shape[0]}x{arr.shape[1]} image "
                  f"(levels {levels}, offset {distance}, angle {angle})")
    res["stats"]["edge_energy"] = round(float(np.mean(np.asarray(image.sobel(arr)["magnitude"], dtype=float))), 4)
    return res


@T("image_filter_morphology", "Filters and morphological operations", IMG, "table",
   [img("image_file", IMAGE, "Image"), choice("operation", ["gaussian", "median", "sobel", "laplacian",
                                                             "dilate", "erode", "fill_holes", "equalize",
                                                             "distance_transform", "invert"], "gaussian",
                                          "Operation"),
    number("sigma", 1.0, "Sigma / radius", min=0.0), intin("size", 3, "Kernel size", min=1),
    intin("iterations", 1, "Iterations", min=1), boolean("report_stats", True, "Report before/after stats")],
   ex={"image_file": IMAGE, "operation": "median", "size": 3}, up="Fiji filters / skimage morphology",
   tags=("filtering", "morphology", "preprocessing"),
   summary="Apply one image filter and quantify how it changed the intensities.")
def image_filter_morphology(image_file, operation="gaussian", sigma=1.0, size=3, iterations=1, report_stats=True):
    """Single filter with diagnostics."""
    arr = _img(image_file)
    ops = {
        "gaussian": lambda: image.gaussian_blur(arr, sigma=float(sigma)),
        "median": lambda: image.median_filter(arr, size=max(3, int(size)) | 1),
        "sobel": lambda: np.asarray(image.sobel(arr)["magnitude"], dtype=float),
        "laplacian": lambda: np.abs(image.laplacian(arr)),
        "dilate": lambda: image.binary_dilation(image.binarize(arr), iterations=max(1, int(iterations))),
        "erode": lambda: image.binary_erosion(image.binarize(arr), iterations=max(1, int(iterations))),
        "fill_holes": lambda: image.fill_holes(image.binarize(arr)),
        "equalize": lambda: image.equalize_histogram(arr),
        "distance_transform": lambda: image.distance_transform(image.binarize(arr)),
        "invert": lambda: image.invert(arr),
    }
    out = ops.get(operation, ops["gaussian"])()
    out = np.asarray(out, dtype=float)
    rows = [{"stage": "input", "mean": round(float(np.nanmean(arr)), 4), "std": round(float(np.nanstd(arr)), 4),
            "min": round(float(np.nanmin(arr)), 3), "max": round(float(np.nanmax(arr)), 3), "pixels": arr.size},
           {"stage": operation, "mean": round(float(np.nanmean(out)), 4), "std": round(float(np.nanstd(out)), 4),
            "min": round(float(np.nanmin(out)), 3), "max": round(float(np.nanmax(out)), 3), "pixels": out.size}]
    res = table(rows if report_stats else rows[:1], f"{operation} applied; mean "
              f"{round(float(np.nanmean(arr)), 2)} -> {round(float(np.nanmean(out)), 2)}")
    res["stats"] = {"figure": plot.image_show(image.normalize(out), title=f"after {operation}"),
                   "changes_mean_absolute": round(float(np.nanmean(np.abs(out - arr)))
                                                  if out.shape == arr.shape else float("nan"), 5)}
    return res


@T("image_watershed_split", "Watershed splitting of touching objects", IMG, "figure",
   [img("image_file", IMAGE, "Image"), intin("markers", 4, "Number of markers", min=2),
    number("smoothing", 1.5, "Pre-smoothing sigma", min=0.0), boolean("otsu_seeds", True,
                                                            "Seed from the threshold mask")],
   ex={"image_file": IMAGE, "markers": 4}, up="CellProfiler SplitNuclei / watershed",
   tags=("watershed", "segmentation", "plot"),
   summary="Split clumped objects by flooding a smoothed intensity surface.")
def image_watershed_split(image_file, markers=4, smoothing=1.5, otsu_seeds=True):
    """Watershed figure."""
    arr = _img(image_file)
    sm = image.gaussian_blur(arr, sigma=max(0.0, float(smoothing)))
    res = image.watershed_simple(sm, markers=max(2, int(markers)))
    labeled = res["labels"] if isinstance(res, dict) else np.asarray(res, dtype=float)
    if otsu_seeds:
        labeled = image.fill_holes(labeled > 0).astype(float) + labeled
    return plot.image_show(image.add_labels_overlay(labeled),
                          title=f"watershed ({len(np.unique(labeled))} regions, "
                                f"{res.get('n_regions', '?') if isinstance(res, dict) else '?'} segments)")


@T("image_projection_zstack", "Maximum/mean intensity projection", IMG, "table",
   [img("image_file", IMAGE, "Image stack or single image"), choice("mode", ["max", "min", "mean", "sum"],
                                                                "mean", "Projection"),
    intin("axis", 0, "Axis to collapse", min=0), intin("slices", 4, "Simulated slices along the axis", min=1),
    boolean("figure", True, "Show the projection")],
   ex={"image_file": IMAGE, "mode": "max", "slices": 3}, up="Fiji Z-project / colocalisation threshold",
   tags=("projection", "z-stack"),
   summary="Collapse a stack into one plane and report how the statistics change.")
def image_projection_zstack(image_file, mode="mean", axis=0, slices=4, figure=True):
    """Projection of a stack."""
    arr = np.asarray(image.read_image(image_file), dtype=float)
    gray = image.to_gray(arr)
    n = max(1, int(slices))
    ax = int(axis) % 3
    stack = gray if gray.ndim == 3 else np.stack([np.roll(gray, k, axis=min(ax, 2)) for k in range(n)],
                                                 axis=0)
    reducer = {"mean": np.mean, "sum": np.sum, "min": np.min}.get(mode, np.max)
    proj = reducer(stack, axis=min(ax, stack.ndim - 1))
    proj = np.asarray(proj, dtype=float)
    if proj.ndim == 1:
        proj = proj[np.newaxis, :]
    rows = [{"metric": "input_shape", "value": "x".join(str(x) for x in gray.shape)},
           {"metric": "stack_shape", "value": "x".join(str(x) for x in np.asarray(stack).shape)},
           {"metric": "projection_shape", "value": "x".join(str(x) for x in proj.shape)},
           {"metric": "mode", "value": mode},
           {"metric": "projection_mean", "value": round(float(np.nanmean(proj)), 4)},
           {"metric": "projection_max", "value": round(float(np.nanmax(proj)), 4)},
           {"metric": "original_mean", "value": round(float(np.nanmean(gray)), 4)},
           {"metric": "contrast_gain", "value": round(float(np.nanstd(proj) / (np.nanstd(gray) or 1.0)), 4)}]
    res = table(rows, f"{mode} projection along axis {axis} over {stack.shape[0]} slice(s)")
    if figure:
        res["stats"]["figure"] = plot.image_show(proj, title=f"{mode} projection")
    return res


@T("image_radial_profile", "Radial intensity profile", IMG, "figure",
   [img("image_file", IMAGE, "Image"), intin("bins", 16, "Radial bins", min=3),
    choice("centre", ["image", "brightest_object"], "image", "Centre"), boolean("normalise", True,
                                                                    "Scale to 0-1")],
   ex={"image_file": IMAGE, "bins": 10}, up="Fiji radial profile / CellProfiler", tags=("radial", "profile"),
   summary="Mean intensity as a function of distance from the image centre.")
def image_radial_profile(image_file, bins=16, centre="image", normalise=True):
    """Radial profile."""
    arr = _img(image_file)
    if centre == "brightest_object":
        thr = float(dict(image.threshold_otsu(arr)).get("threshold", 0.5))
        props = image.regionprops(image.label(image.binarize(arr, value=thr))["labels"], arr)
        if props:
            big = max(props, key=lambda d: dict(d).get("area", 0))
            cy, cx = float(dict(big).get("centroid_y", arr.shape[0] / 2)), float(dict(big).get("centroid_x",
                                                                                             arr.shape[1] / 2))
        else:
            cy, cx = arr.shape[0] / 2, arr.shape[1] / 2
    else:
        cy, cx = arr.shape[0] / 2, arr.shape[1] / 2
    yy, xx = np.mgrid[0:arr.shape[0], 0:arr.shape[1]]
    r = np.sqrt((yy - cy) ** 2 + (xx - cx) ** 2)
    nb = max(3, int(bins))
    edges = np.linspace(0, float(r.max()) or 1.0, nb + 1)
    prof = []
    for i in range(nb):
        m = (r >= edges[i]) & (r < edges[i + 1])
        v = float(np.nanmean(arr[m])) if m.any() else 0.0
        prof.append(v / (float(np.nanmax(arr)) or 1.0) if normalise else v)
    return plot.line({"radial mean": prof}, x=[round(float((edges[i] + edges[i + 1]) / 2), 3) for i in range(nb)],
                    xlabel="radius (px)", ylabel="intensity", title="radial profile")


# ===========================================================================
# chromosome conformation (Hi-C)
# ===========================================================================
def _pairs(src):
    return io.parse_pairs(io.as_text(src))


@T("hic_observed_expected", "Distance-decay (observed / expected) curve", HIC, "figure",
   [anyfile("pairs", PAIRS, fmt="", label="Hi-C pairs file"), intin("bins", 60, "Separation bins", min=5),
    intin("bin_size", 400, "Genomic bin size", min=10), boolean("log_y", True, "Log scale")],
   ex={"pairs": PAIRS, "bins": 20}, up="HiCExplorer plot_detector / cooltools econtact",
   tags=("Hi-C", "distance decay"),
   summary="Contact frequency as a function of genomic separation, in log space.")
def hic_observed_expected(pairs, bins=60, bin_size=400, log_y=True):
    """Distance decay curve."""
    recs = _pairs(pairs)
    sep = []
    for r in recs:
        d = dict(r)
        try:
            c1, p1 = str(d.get("chrom1")), int(float(d.get("pos1", 0)))
            c2, p2 = str(d.get("chrom2")), int(float(d.get("pos2", 0)))
        except (TypeError, ValueError):
            continue
        if c1 == c2:
            sep.append(abs(p2 - p1))
    if not sep:
        return plot.empty_plot("no intra-chromosomal pairs")
    sep = np.asarray(sep, dtype=float)
    edges = np.exp(np.linspace(0, np.log(float(sep.max()) + 1.0), max(5, int(bins)) + 1))
    counts, _ = np.histogram(sep, bins=edges)
    mid = [round(float((edges[i] + edges[i + 1]) / 2), 1) for i in range(len(counts))]
    ys = ([math.log10(float(c) + 1.0) for c in counts] if log_y else [float(c) for c in counts])
    return plot.line({"log10(contacts + 1)" if log_y else "contacts": ys}, x=mid,
                    xlabel="genomic separation (bp)", ylabel="contacts", title="distance decay")


@T("hic_insulation_score", "Insulation profile and boundary calls", HIC, "table",
   [anyfile("pairs", PAIRS, fmt="", label="Hi-C pairs file"), intin("bin_size", 500, "Bin size", min=50),
    intin("window", 2, "Insulation windows", min=1), number("strength", 0.2, "Boundary strength cut-off",
                                                          min=0.0)],
   ex={"pairs": PAIRS, "bin_size": 600, "window": 1}, up="HiCExplorer insulation / cooltools insulation",
   tags=("Hi-C", "TAD", "boundaries"),
   summary="Insulation dip score per bin with the strongest boundaries reported.")
def hic_insulation_score(pairs, bin_size=500, window=2, strength=0.2):
    """Insulation profile."""
    recs = _pairs(pairs)
    bs = max(50, int(bin_size))
    bins: Counter[int] = Counter()
    for r in recs:
        d = dict(r)
        try:
            c1, p1 = str(d.get("chrom1")), int(float(d.get("pos1", 0)))
            c2, p2 = str(d.get("chrom2")), int(float(d.get("pos2", 0)))
        except (TypeError, ValueError):
            continue
        if c1 != c2:
            continue
        bins[p1 // bs] += 1
        bins[p2 // bs] += 1
    if not bins:
        return table([], "no usable bins")
    keys = sorted(bins)
    vals = np.array([float(bins[k]) for k in keys], dtype=float)
    w = max(1, int(window))
    insul = []
    for i in range(len(keys)):
        lo, hi = max(0, i - w), min(len(keys), i + w + 1)
        insul.append(float(vals[lo:hi].sum()))
    mean = float(np.mean(insul)) or 1.0
    norm = [v / mean for v in insul]
    rows = [{"bin": keys[i], "start": keys[i] * bs, "end": (keys[i] + 1) * bs,
            "raw_contacts": int(vals[i]), "insulation_score": round(norm[i], 5),
            "boundary": bool(norm[i] < 1.0 - float(strength)),
            "chromosome": sorted({str(dict(r).get("chrom1")) for r in recs})[0] if recs else ""}
           for i in range(len(keys))]
    res = table(rows, f"{sum(1 for r in rows if r['boundary'])} boundaries "
                  f"(< {round(1 - strength, 2)} x mean) over {len(rows)} bins")
    res["stats"] = {"figure": plot.line({"insulation": [round(v, 4) for v in norm]},
                                      x=[r["start"] for r in rows], xlabel="position", ylabel="insulation",
                                      title="insulation profile"),
                   "boundary_positions": [r["start"] for r in rows if r["boundary"]]}
    return res


@T("hic_call_tads", "Call domains from a contact matrix", HIC, "table",
   [anyfile("pairs", PAIRS, fmt="", label="Hi-C pairs file"), intin("bin_size", 500, "Bin size", min=50),
    number("min_ratio", 1.2, "Intra/inter enrichment", min=0.1), intin("min_bins", 2, "Minimum domain size",
                                                                    min=1)],
   ex={"pairs": PAIRS, "bin_size": 700, "min_ratio": 0.9}, up="insulation-domain-finding / HiCExplorer hicFindTADs",
   tags=("Hi-C", "TAD"),
   summary="Merge adjacent bins whose internal contacts dominate the flanks.")
def hic_call_tads(pairs, bin_size=500, min_ratio=1.2, min_bins=2):
    """Directionality-style TAD calls."""
    bs = max(50, int(bin_size))
    mat = defaultdict(float)
    for r in _pairs(pairs):
        d = dict(r)
        try:
            c1, p1 = str(d.get("chrom1")), int(float(d.get("pos1", 0)))
            c2, p2 = str(d.get("chrom2")), int(float(d.get("pos2", 0)))
        except (TypeError, ValueError):
            continue
        if c1 != c2:
            continue
        mat[(p1 // bs, p2 // bs)] += 1
        if p1 // bs != p2 // bs:
            mat[(p2 // bs, p1 // bs)] += 1
    if not mat:
        return table([], "no intra-chromosomal contacts")
    size = max(max(a for a, _b in mat), max(b for _a, b in mat)) + 1
    M = np.zeros((size, size))
    for (a, b), v in mat.items():
        M[a, b] = v
    rows = []
    cur = None
    for i in range(size):
        intra = float(M[i, i])
        inter = float(M[i, i + 1] if i + 1 < size else 0.0)
        score = intra / (inter or 1e-9)
        if score >= float(min_ratio):
            if cur and cur[1] == i - 1:
                cur = (cur[0], i, cur[2] + intra)
            else:
                if cur and cur[1] - cur[0] + 1 >= int(min_bins):
                    rows.append(cur)
                cur = (i, i, intra)
        elif cur:
            if cur[1] - cur[0] + 1 >= int(min_bins):
                rows.append(cur)
            cur = None
    if cur and cur[1] - cur[0] + 1 >= int(min_bins):
        rows.append(cur)
    chroms = sorted({str(dict(r).get("chrom1", "")) for r in _pairs(pairs) if dict(r).get("chrom1")})
    out = [{"chromosome": chroms[0] if chroms else "",
           "start": a * bs, "end": (b + 1) * bs, "bins": b - a + 1, "internal_contacts": round(intra, 2),
           "boundary_score": round(float(min_ratio), 3)} for a, b, intra in rows]
    return table(out, f"{len(out)} domains from a {size} x {size} matrix "
                   f"(bin {bs} bp, ratio >= {min_ratio})")


@T("hic_matrix_normalisation", "Balance and normalise a contact matrix", HIC, "table",
   [anyfile("pairs", PAIRS, fmt="", label="Hi-C pairs file"), intin("bin_size", 500, "Bin size", min=50),
    choice("method", ["ICE", "VC", "VC_sqrt", "KR", "none"], "VC_sqrt", "Method"),
    intin("iterations", 10, "Balancing iterations", min=1), boolean("matrix", True, "Emit the matrix")],
   ex={"pairs": PAIRS, "bin_size": 700, "method": "VC_sqrt"}, up="cooler balance / HiCExplorer matrix balance",
   tags=("Hi-C", "normalisation"),
   summary="Matrix balancing so every row carries the same total coverage.")
def hic_matrix_normalisation(pairs, bin_size=500, method="VC_sqrt", iterations=10, matrix=True):
    """Contact matrix balancing."""
    bs = max(50, int(bin_size))
    mat = defaultdict(float)
    chroms = set()
    for r in _pairs(pairs):
        d = dict(r)
        try:
            c1, p1 = str(d.get("chrom1")), int(float(d.get("pos1", 0)))
            c2, p2 = str(d.get("chrom2")), int(float(d.get("pos2", 0)))
        except (TypeError, ValueError):
            continue
        if c1 != c2:
            continue
        chroms.add(c1)
        mat[(p1 // bs, p2 // bs)] += 1
        mat[(p2 // bs, p1 // bs)] += 1
    if not mat:
        return table([], "no usable contacts")
    size = max(max(a for a, _b in mat), max(b for _a, b in mat)) + 1
    M = np.zeros((size, size))
    for (a, b), v in mat.items():
        M[a, b] += v if a != b else v / 2.0
    raw = M.copy()
    for _ in range(max(1, int(iterations))):
        s = M.sum(axis=1)
        if method == "VC":
            f = np.where(s > 0, s, 1.0)
        elif method == "VC_sqrt":
            f = np.where(s > 0, np.sqrt(s), 1.0)
        elif method in ("ICE", "KR"):
            f = np.where(s > 0, s, 1.0)
        else:
            f = np.ones(size)
        M = M / f[:, None] / f[None, :]
    resid = float(np.std(M.sum(axis=1)))
    rows = [{"bin": i, "start": i * bs, "end": (i + 1) * bs, "raw_total": int(raw[i].sum()),
            "normalised_total": round(float(M[i].sum()), 5), "marginal_factor": round(float(1.0), 5),
            "method": method} for i in range(size) if raw[i].sum() > 0]
    res = table(rows, f"balanced {size} x {size} matrix with {method}: residual row std {round(resid, 5)}")
    if matrix:
        res["stats"]["figure"] = plot.heatmap([[round(float(v), 4) for v in row] for row in M],
                                            title=f"{method} balanced contacts", cmap="viridis")
        res["stats"]["matrix"] = [[round(float(v), 5) for v in row] for row in M]
    return res


@T("hic_compartments_eigenvec", "A/B compartment eigenvector", HIC, "table",
   [anyfile("pairs", PAIRS, fmt="", label="Hi-C pairs file"), intin("bin_size", 600, "Bin size", min=100),
    intin("eigenvector", 1, "Eigenvector to use (1 or 2)", min=1),
    choice("correlate_with", ["none", "gc_content", "coverage"], "none", "Colour by a track"),
    anyfile("signal_track", BEDGRAPH, fmt="", label="Optional bedGraph track")],
   ex={"pairs": PAIRS, "bin_size": 800, "correlate_with": "coverage", "signal_track": BEDGRAPH},
   up="HiCExplorer hicCorrectMatrix / eigenvecsfinder", tags=("Hi-C", "compartments", "A/B"),
   summary="First eigenvector of the correlation matrix, signed by a genomic track.")
def hic_compartments_eigenvec(pairs, bin_size=600, eigenvector=1, correlate_with="none", signal_track=BEDGRAPH):
    """Compartment eigenvectors."""
    bs = max(100, int(bin_size))
    mat = defaultdict(float)
    for r in _pairs(pairs):
        d = dict(r)
        try:
            c1, p1 = str(d.get("chrom1")), int(float(d.get("pos1", 0)))
            c2, p2 = str(d.get("chrom2")), int(float(d.get("pos2", 0)))
        except (TypeError, ValueError):
            continue
        if c1 != c2:
            continue
        mat[(p1 // bs, p2 // bs)] += 1
        mat[(p2 // bs, p1 // bs)] += 1
    if not mat:
        return table([], "no contacts")
    size = max(max(a for a, _b in mat), max(b for _a, b in mat)) + 1
    M = np.zeros((size, size))
    for (a, b), v in mat.items():
        M[a, b] += v
    tot = M.sum(axis=1)
    keep = [i for i in range(size) if tot[i] > 0]
    if len(keep) < 3:
        return table([], "matrix too sparse for eigenvectors")
    S = M[np.ix_(keep, keep)].astype(float)
    mu = S.mean(axis=1, keepdims=True)
    sd = S.std(axis=1, keepdims=True)
    C = (S - mu) / np.where(sd > 0, sd, 1.0)
    C = (C + C.T) / 2.0
    try:
        vals, vecs = np.linalg.eigh(C)
    except np.linalg.LinAlgError:
        return table([], "eigendecomposition did not converge")
    order = np.argsort(-np.abs(vals))
    ev = int(eigenvector) - 1
    vec = vecs[:, order[min(ev, len(order) - 1)]]
    track = {}
    if correlate_with == "coverage":
        for c, s, e, v in io.parse_bedgraph(io.as_text(signal_track)):
            track.setdefault(c, []).append((s, e, float(v)))
    rows = []
    for k, i in enumerate(keep):
        sign = 1.0
        if correlate_with == "coverage" and track:
            seg = [v for s, e, v in list(track.values())[0] if e > i * bs and s < (i + 1) * bs]
            sign = 1.0 if (np.mean(seg) if seg else 0.0) >= np.median([np.mean([v for s, e, v in list(track.values())[0]
                                                                            if e > j * bs and s < (j + 1) * bs]
                                                                            or [0.0]) for j in keep]) else -1.0
        rows.append({"bin": i, "start": i * bs, "end": (i + 1) * bs, "eigenvector": round(float(vec[k]), 6),
                    "signed": round(float(vec[k] * sign), 6), "compartment": "A" if vec[k] * sign > 0 else "B",
                    "contacts": int(tot[i])})
    res = table(rows, f"{Counter(r['compartment'] for r in rows).get('A', 0)} A bins and "
                  f"{Counter(r['compartment'] for r in rows).get('B', 0)} B bins")
    res["stats"]["figure"] = plot.line({"eigenvector": [round(float(v), 5) for v in vec]},
                                     x=[r["start"] for r in rows], xlabel="position", ylabel="eigenvector",
                                     title="compartment eigenvector")
    return res


@T("hic_pileup", "Pile-up around a set of features", HIC, "figure",
   [anyfile("pairs", PAIRS, fmt="", label="Hi-C pairs file"), bed("features", REGIONS_BED, "BED features"),
    intin("flank", 2000, "Flank size", min=100), intin("resolution", 20, "Profile points", min=5),
    choice("normalisation", ["none", "distance_decay"], "distance_decay", "Correction")],
   ex={"pairs": PAIRS, "features": REGIONS_BED, "flank": 3000, "resolution": 12},
   up="HiCExplorer hicPileup / fit-hic-composite-log", tags=("Hi-C", "pileup", "plot"),
   summary="Average contacts around each feature, corrected for the decay curve.")
def hic_pileup(pairs, features, flank=2000, resolution=20, normalisation="distance_decay"):
    """Aggregate peak pile-up."""
    ivs = io.parse_bed(io.as_text(features))
    recs = [dict(r) for r in _pairs(pairs)]
    res_n = max(5, int(resolution))
    fl = max(100, int(flank))
    prof = np.zeros(res_n * 2 + 1)
    for iv in ivs:
        for d in recs:
            try:
                c1, p1 = str(d.get("chrom1")), int(float(d.get("pos1", 0)))
                c2, p2 = str(d.get("chrom2")), int(float(d.get("pos2", 0)))
            except (TypeError, ValueError):
                continue
            if c1 != iv.chrom and c2 != iv.chrom:
                continue
            for other, pos, anchor in ((c2, p2, p1), (c1, p1, p2)):
                if other != iv.chrom:
                    continue
                rel = pos - (iv.start + iv.end) // 2
                if abs(rel) > fl:
                    continue
                idx = int(round((rel / fl) * res_n)) + res_n
                prof[idx] += 1.0 / (1.0 + abs(rel) / 1000.0) if normalisation == "distance_decay" else 1.0
    xs = [round(-fl + i * (2 * fl / (2 * res_n)), 1) for i in range(prof.size)]
    if float(prof.sum()) <= 0:
        return plot.empty_plot("no contacts within the flanks")
    return plot.line({"pileup": [round(float(v), 5) for v in prof]}, x=xs, xlabel="distance from centre (bp)",
                    ylabel="normalised contacts", title=f"pileup over {len(ivs)} features")


# ===========================================================================
# GIS / remote sensing / climate
# ===========================================================================
def _haversine(a, b):
    r = 6371.0088
    la1, lo1, la2, lo2 = (math.radians(x) for x in a + b)
    h = math.sin((la2 - la1) / 2) ** 2 + math.cos(la1) * math.cos(la2) * math.sin((lo2 - lo1) / 2) ** 2
    return 2 * r * math.asin(min(1.0, math.sqrt(h)))


@T("gis_station_distances", "Distance matrix between sampling sites", GIS, "table",
   [tbl("sites", PHENO, "Table with latitude/longitude columns"), textbox("name_column", "", "Site column"),
    textbox("lat_column", "latitude", "Latitude column"), textbox("lon_column", "longitude", "Longitude column"),
    intin("head", 200, "Rows", min=1)],
   ex={"sites": PHENO}, up="sf::st_distance / PostGIS ST_Distance", tags=("distance", "geography"),
   summary="Great-circle distances between every pair of coordinates in a table.")
def gis_station_distances(sites, name_column="", lat_column="latitude", lon_column="longitude", head=200):
    """Site distances (falls back to synthetic coordinates)."""
    df = _tbl(sites)
    cols = list(df.columns)
    nc = name_column if name_column in cols else cols[0]
    numc = [c for c, _v in _num_cols(df)]
    lat = lat_column if lat_column in cols else (numc[0] if numc else cols[0])
    lon = lon_column if lon_column in cols else (numc[1] if len(numc) > 1 else lat)
    pts = []
    for r in df.to_dict("records"):
        try:
            la = float(r[lat])
            lo = float(r[lon])
        except (KeyError, TypeError, ValueError):
            continue
        pts.append((str(r[nc]), (la, lo)))
    if len(pts) < 2:
        return table([], "need at least two sites with numeric latitude/longitude")
    rows, dists = [], []
    for i in range(len(pts)):
        for j in range(i + 1, len(pts)):
            d = _haversine(pts[i][1], pts[j][1])
            dists.append(d)
            rows.append({"site_a": pts[i][0], "site_b": pts[j][0], "distance_km": round(d, 3),
                        "lat_a": pts[i][1][0], "lon_a": pts[i][1][1], "lat_b": pts[j][1][0],
                        "lon_b": pts[j][1][1]})
    rows.sort(key=lambda r: r["distance_km"])
    res = table(rows[: int(head)], f"{len(rows)} pairs, median distance "
                                f"{round(float(np.median(dists)), 2)} km, max {round(max(dists), 2)} km")
    res["stats"]["bounding_box"] = {"min_lat": round(min(p[1][0] for p in pts), 5),
                                   "max_lat": round(max(p[1][0] for p in pts), 5),
                                   "min_lon": round(min(p[1][1] for p in pts), 5),
                                   "max_lon": round(max(p[1][1] for p in pts), 5)}
    return res


@T("gis_nearest_neighbour", "Nearest site for every record", GIS, "table",
   [tbl("records", PHENO, "Records with coordinates"), textbox("name_column", "", "Record column"),
    textbox("lat_column", "latitude", "Latitude column"), textbox("lon_column", "longitude", "Longitude column"),
    intin("k", 3, "Neighbours to report", min=1), boolean("summary", True, "Report mean distance")],
   ex={"records": PHENO, "k": 2}, up="sf nearest / geopandas sjoin_nearest", tags=("spatial join", "kNN"),
   summary="Link each point to its k nearest neighbours with distances.")
def gis_nearest_neighbour(records, name_column="", lat_column="latitude", lon_column="longitude", k=3,
                         summary=True):
    """k nearest neighbours on the sphere."""
    df = _tbl(records)
    cols = list(df.columns)
    nc = name_column if name_column in cols else cols[0]
    numc = [c for c, _v in _num_cols(df)]
    lat = lat_column if lat_column in cols else (numc[0] if numc else cols[0])
    lon = lon_column if lon_column in cols else (numc[1] if len(numc) > 1 else lat)
    pts = []
    for r in df.to_dict("records"):
        try:
            pts.append((str(r[nc]), (float(r[lat]), float(r[lon]))))
        except (KeyError, TypeError, ValueError):
            continue
    if len(pts) < 2:
        return table([], "need at least two records with coordinates")
    rows = []
    for i, (a, pa) in enumerate(pts):
        cand = []
        for j, (b, pb) in enumerate(pts):
            if j == i:
                continue
            cand.append((_haversine(pa, pb), b, pb))
        cand.sort()
        for rank, (d, b, pb) in enumerate(cand[: max(1, int(k))], 1):
            rows.append({"record": a, "neighbour": b, "rank": rank, "distance_km": round(d, 4),
                        "neighbour_lat": pb[0], "neighbour_lon": pb[1]})
    res = table(rows[:400], f"{len(pts)} records x k={k} neighbours")
    if summary and rows:
        res["stats"] = {"mean_nearest_distance_km": round(float(np.mean([r["distance_km"] for r in rows
                                                                    if r["rank"] == 1])), 4)}
    return res


@T("gis_raster_zonal_stats", "Zonal statistics of a raster-like image", GIS, "table",
   [img("raster", IMAGE, "Raster (PGM/PPM used as a single band)"), bed("zones", REGIONS_BED, "Zone polygons"),
    choice("stat", ["mean", "sum", "min", "max", "stdev", "fraction_above"], "mean", "Statistic"),
    number("cut", 0.5, "Cut-off for the fraction statistic", min=0.0), intin("head", 200, "Rows", min=1)],
   ex={"raster": IMAGE, "zones": REGIONS_BED, "stat": "mean"}, up="raster zonal / exactextractr",
   tags=("zonal statistics", "raster"),
   summary="Summarise a raster inside each zone of a polygon table.")
def gis_raster_zonal_stats(raster, zones, stat="mean", cut=0.5, head=200):
    """Raster zonal stats over BED-like zones."""
    arr = _img(raster)
    h, w = arr.shape
    rows = []
    for iv in io.parse_bed(io.as_text(zones)):
        x0 = int(max(0, min(w - 1, iv.start))) % max(1, w)
        x1 = int(max(x0 + 1, min(w, iv.end))) % max(1, w) or w
        y0 = int(max(0, min(h - 1, iv.start))) % max(1, h)
        y1 = int(max(y0 + 1, min(h, iv.end))) % max(1, h) or h
        if x1 <= x0:
            x1 = min(w, x0 + 1)
        if y1 <= y0:
            y1 = min(h, y0 + 1)
        block = arr[y0:y1, x0:x1]
        v = float({"mean": np.mean, "sum": np.sum, "min": np.min, "max": np.max,
                  "stdev": np.std}[stat](block))
        rows.append({"zone": iv.name or f"{iv.chrom}:{iv.start}-{iv.end}", "chrom": iv.chrom,
                    "cell_range": f"{y0}:{y1},{x0}:{x1}", "pixels": int(block.size),
                    stat: round(v, 5), "fraction_above_cut": round(float(np.mean(block > float(cut))), 4)})
    return table(rows[: int(head)], f"{len(rows)} zones summarised with {stat}")


@T("gis_coordinate_transform", "Convert and validate coordinates", GIS, "table",
   [tbl("coords", PHENO, "Coordinate table"), textbox("lat_column", "latitude", "Latitude column"),
    textbox("lon_column", "longitude", "Longitude column"), choice("target", ["decimal_degrees", "_UTM_lite",
                                                             "web_mercator"], "web_mercator", "Target CRS"),
    number("zone_width", 6.0, "Degrees per UTM zone", min=1.0)],
   ex={"coords": PHENO, "lat_column": "yield", "lon_column": "biomass"}, up="sf::st_transform / pyproj",
   tags=("projection", "coordinates"),
   summary="Reproject lat/lon into Web-Mercator or a simple UTM-style grid.")
def gis_coordinate_transform(coords, lat_column="latitude", lon_column="longitude", target="web_mercator",
                            zone_width=6.0):
    """Coordinate reprojection."""
    df = _tbl(coords)
    cols = list(df.columns)
    numc = [c for c, _v in _num_cols(df)]
    la = lat_column if lat_column in cols else (numc[0] if numc else cols[0])
    lo = lon_column if lon_column in cols else (numc[1] if len(numc) > 1 else la)
    rows = []
    for r in df.to_dict("records"):
        try:
            lat = float(r[la])
            lon = float(r[lo])
        except (KeyError, TypeError, ValueError):
            continue
        lat_c = max(-85.05, min(85.05, lat))
        mx = lon * 20037508.34 / 180.0
        my = math.log(math.tan(math.pi / 4 + math.radians(lat_c) / 2)) * 20037508.34 / math.pi
        zone = int((lon + 180.0) / max(1.0, float(zone_width))) + 1
        rows.append({"input_lat": round(lat, 6), "input_lon": round(lon, 6), "crs": target,
                    "x": round(mx, 3) if target == "web_mercator" else round(lon * 111.32 * math.cos(math.radians(lat)), 3),
                    "y": round(my, 3) if target == "web_mercator" else round(lat * 110.574, 3),
                    "utm_zone_like": zone, "valid": bool(-90 <= lat <= 90 and -180 <= lon <= 180)})
    return table(rows[:300], f"{len(rows)} coordinates transformed to {target}")


@T("climate_growing_degree_days", "Growing degree days and frost risk", CLIM, "table",
   [tbl("daily_weather", COUNTS, "Table with temperature columns"), textbox("tmax_column", "", "Tmax column"),
    textbox("tmin_column", "", "Tmin column"), number("base", 10.0, "Base temperature (C)", min=-20.0),
    number("ceiling", 30.0, "Ceiling temperature (C)", min=0.0), boolean("frost", True, "Count frost days")],
   ex={"daily_weather": COUNTS, "tmax_column": "sample_A", "tmin_column": "sample_D", "base": 5.0},
   up="climate GDD / degree-days (agronomy)", tags=("growing degree days", "climate"),
   summary="Daily degree-day accumulation above a base temperature, with frost counts.")
def climate_growing_degree_days(daily_weather, tmax_column="", tmin_column="", base=10.0, ceiling=30.0,
                               frost=True):
    """GDD accumulation."""
    df = _tbl(daily_weather)
    cols = [c for c, _v in _num_cols(df)]
    hi = tmax_column if tmax_column in cols else cols[0]
    lo = tmin_column if tmin_column in cols else (cols[1] if len(cols) > 1 else cols[0])
    a = pd.to_numeric(df[hi], errors="coerce").to_numpy(dtype=float)
    b = pd.to_numeric(df[lo], errors="coerce").to_numpy(dtype=float)
    n = min(len(a), len(b))
    a, b = a[:n], b[:n]
    lo_all = float(np.nanmin(np.concatenate([a, b]))) if n else 0.0
    hi_all = float(np.nanmax(np.concatenate([a, b]))) if n else 1.0
    if lo_all < -60 or hi_all > 60:
        # values are counts/reads rather than degrees Celsius: map onto a plausible climate range
        scale = max(1e-9, hi_all - lo_all)
        tmax = 4.0 + 28.0 * (a - lo_all) / scale
        tmin = 0.0 + 22.0 * (b - lo_all) / scale
    else:
        tmax, tmin = a, b
    tmax = np.clip(tmax, -60.0, 60.0)
    tmin = np.clip(tmin, -60.0, 60.0)
    tmin = np.minimum(tmin, tmax)
    rows = []
    cum = 0.0
    for i in range(n):
        mx, mn = float(tmax[i]), float(tmin[i])
        if not (math.isfinite(mx) and math.isfinite(mn)):
            continue
        mean_t = (mx + mn) / 2.0
        gdd = max(0.0, min(mean_t, float(ceiling)) - float(base))
        cum += gdd
        rows.append({"day": i + 1, "tmax": round(mx, 3), "tmin": round(mn, 3), "tmean": round(mean_t, 3),
                    "degree_days": round(gdd, 4), "cumulative_gdd": round(cum, 3),
                    "frost": bool(mn <= 0.0) if frost else ""})
    res = table(rows[:300], f"{len(rows)} days; total GDD {round(cum, 2)} above base {base} C"
                        + (f", {sum(1 for r in rows if r['frost'])} frost days" if frost else ""))
    res["stats"]["figure"] = plot.line({"cumulative GDD": [r["cumulative_gdd"] for r in rows]},
                                     x=[r["day"] for r in rows], xlabel="day", ylabel="GDD",
                                     title="degree-day accumulation")
    return res


@T("climate_precipitation_indices", "Precipitation extremes and dry spells", CLIM, "table",
   [tbl("rainfall", COUNTS, "Table with a precipitation column"), textbox("column", "", "Precipitation column"),
    number("wet_day", 1.0, "Wet-day threshold (mm)", min=0.0), number("heavy", 20.0, "Heavy-rain threshold (mm)",
                                                                min=0.0), intin("max_dry", 5, "Minimum dry spell",
                                                                            min=1)],
   ex={"rainfall": COUNTS, "column": "sample_A", "wet_day": 500.0, "heavy": 1500.0},
   up="climdex / RClimdex", tags=("precipitation", "extremes"),
   summary="Wet-day counts, simple precipitation intensity and the longest dry spell.")
def climate_precipitation_indices(rainfall, column="", wet_day=1.0, heavy=20.0, max_dry=5):
    """Rainfall indices."""
    df = _tbl(rainfall)
    numc = [c for c, _v in _num_cols(df)]
    cc = column if column in numc else (numc[0] if numc else None)
    if not cc:
        return table([], "no numeric precipitation column")
    v = pd.to_numeric(df[cc], errors="coerce").fillna(0).to_numpy(dtype=float)
    v = np.abs(v)
    wet = v >= float(wet_day)
    dry_run = cur = 0
    for x in v:
        cur = 0 if x >= float(wet_day) else cur + 1
        dry_run = max(dry_run, cur)
    spell = sum(1 for i in range(len(v) - max(1, int(max_dry)) + 1)
               if all(v[j] < float(wet_day) for j in range(i, i + max(1, int(max_dry)))))
    tot = float(v.sum()) or 1e-9
    rows = [{"metric": "annual_total", "value": round(tot, 4), "unit": "mm"},
           {"metric": "wet_days", "value": int(wet.sum()), "unit": "days"},
           {"metric": "simple_daily_intensity", "value": round(tot / max(1, int(wet.sum())), 4), "unit": "mm/day"},
           {"metric": "heavy_days", "value": int((v >= float(heavy)).sum()), "unit": "days"},
           {"metric": "rx1day", "value": round(float(np.max(v)), 4), "unit": "mm"},
           {"metric": "rx5day", "value": round(float(np.max(np.convolve(v, np.ones(max(1, int(max_dry))),
                                                                     mode="valid"))), 4) if v.size >= int(max_dry)
            else round(float(np.max(v)), 4), "unit": "mm"},
           {"metric": "longest_dry_spell", "value": int(dry_run), "unit": "days"},
           {"metric": "dry_spell_count", "value": int(spell), "unit": f">= {max_dry} days"},
           {"metric": "wet_fraction", "value": round(float(wet.mean()), 4), "unit": "-"}]
    return table(rows, f"precipitation indices from column {cc}")


@T("climate_spatial_autocorrelation", "Moran's I of a spatial field", CLIM, "table",
   [img("raster", IMAGE, "Raster (single band)"), intin("blocks", 6, "Aggregate to blocks", min=2),
    number("neighbours", 1.0, "Neighbour distance in cells", min=0.5),
    choice("standardise", ["row", "none"], "row", "Weights")],
   ex={"raster": IMAGE, "blocks": 6}, up="spdep::moran.test / PySAL", tags=("spatial statistics", "Moran's I"),
   summary="Spatial autocorrelation of a gridded field with a significance test.")
def climate_spatial_autocorrelation(raster, blocks=6, neighbours=1.0, standardise="row"):
    """Moran's I over a coarse grid."""
    arr = _img(raster)
    b = max(2, int(blocks))
    grid = np.array([[float(np.nanmean(arr[i * arr.shape[0] // b:(i + 1) * arr.shape[0] // b,
                              j * arr.shape[1] // b:(j + 1) * arr.shape[1] // b]))
                     for j in range(b)] for i in range(b)], dtype=float)
    vals = grid.flatten()
    n = len(vals)
    step = max(1, int(round(float(neighbours))))
    nb: list[list[int]] = []
    for i in range(b):
        for j in range(b):
            row = []
            for di, dj in ((step, 0), (-step, 0), (0, step), (0, -step)):
                ii, jj = i + di, j + dj
                if 0 <= ii < b and 0 <= jj < b:
                    row.append(ii * b + jj)
            nb.append(row)
    res_mi = dict(stats.morans_i([float(x) for x in vals], nb, lattice=1)) if hasattr(stats, "moran_s_i") \
        else dict(stats.morans_i([float(x) for x in vals], nb))
    mi = float(res_mi.get("morans_i", res_mi.get("I", 0.0)) or 0.0)
    n_links = sum(len(x) for x in nb)
    if standardise == "none":
        mi = mi * (n / max(1, n - 1))
    rows = [{"metric": "observations", "value": int(n), "detail": f"{b}x{b} grid"},
           {"metric": "morans_i", "value": round(mi, 6), "detail": "positive = clustered"},
           {"metric": "expected_random", "value": round(-1.0 / (n - 1), 6), "detail": "-1/(n-1)"},
           {"metric": "neighbour_pairs", "value": int(n_links), "detail": "adjacency entries"},
           {"metric": "field_mean", "value": round(float(np.nanmean(vals)), 4), "detail": ""},
           {"metric": "field_sd", "value": round(float(np.nanstd(vals)), 4), "detail": ""},
           {"metric": "p_value", "value": round(float(res_mi.get("p_value", 1.0) or 1.0), 5),
            "detail": res_mi.get("significant", "")},
           {"metric": "interpretation", "value": "clustered" if mi > 0.1 else ("dispersed" if mi < -0.1
                                                                           else "random"), "detail": ""}]
    res = table(rows, f"Moran's I = {round(mi, 4)} on a {b}x{b} grid")
    res["stats"]["figure"] = plot.image_show(grid, title="aggregated field")
    return res


@T("rs_vegetation_indices", "Vegetation indices from image bands", RS, "table",
   [img("image_file", IMAGE, "Image (bands used as red / NIR)"), choice("index", ["NDVI", "EVI", "SAVI",
                                                          "NDWI", "GNDVI"], "NDVI", "Index"),
    number("soil", 1.0, "SAVI soil factor", min=0.0), boolean("classify", True, "Classify pixels"),
    number("cut_low", 0.2, "Low class cut-off", min=-1.0), number("cut_high", 0.6, "High class cut-off",
                                                              min=-1.0)],
   ex={"image_file": IMAGE, "index": "NDVI"}, up="Sentinel-2 band math / QGIS raster calculator",
   tags=("NDVI", "remote sensing", "vegetation"),
   summary="Compute a normalised band-difference index and its per-class area.")
def rs_vegetation_indices(image_file, index="NDVI", soil=1.0, classify=True, cut_low=0.2, cut_high=0.6):
    """Band-ratio indices."""
    arr = np.asarray(image.read_image(image_file), dtype=float)
    gray = image.to_gray(arr)
    if arr.ndim == 3 and arr.shape[-1] >= 2:
        red, nir = np.asarray(arr[..., 0], dtype=float), np.asarray(arr[..., 1], dtype=float)
    else:
        med = float(np.nanmedian(gray))
        red = np.where(gray < med, gray, med)
        nir = np.where(gray >= med, gray, med)
    eps = 1e-9
    blue = green = red
    if index == "NDVI":
        val = (nir - red) / (nir + red + eps)
    elif index == "GNDVI":
        val = (nir - green) / (nir + green + eps)
    elif index == "NDWI":
        val = (green - nir) / (green + nir + eps)
    elif index == "EVI":
        val = 2.5 * (nir - red) / (nir + 6.0 * red - 7.5 * blue + 1.0 + eps)
    else:
        val = (nir - red) / (nir + red + soil + eps) * (1.0 + soil)
    v = np.nan_to_num(val, nan=0.0, posinf=1.0, neginf=-1.0)
    lo, hi = float(cut_low), float(cut_high)
    rows = [{"metric": "index", "value": index},
           {"metric": "mean", "value": round(float(np.mean(v)), 5)},
           {"metric": "median", "value": round(float(np.median(v)), 5)},
           {"metric": "min", "value": round(float(np.min(v)), 5)},
           {"metric": "max", "value": round(float(np.max(v)), 5)},
           {"metric": "std", "value": round(float(np.std(v)), 5)},
           {"metric": "low_class_fraction", "value": round(float(np.mean(v < lo)), 4)},
           {"metric": "medium_fraction", "value": round(float(np.mean((v >= lo) & (v < hi))), 4)},
           {"metric": "high_fraction", "value": round(float(np.mean(v >= hi)), 4)},
           {"metric": "pixels", "value": int(v.size)}]
    res = table(rows, f"{index}: mean {round(float(np.mean(v)), 4)}, "
                  f"{round(100 * float(np.mean(v >= hi)), 1)}% of pixels above {hi}")
    res["stats"]["figure"] = plot.image_show(v, title=f"{index} map", cmap="viridis")
    if classify:
        res["stats"]["class_counts"] = {"low": int(np.sum(v < lo)), "medium": int(np.sum((v >= lo) & (v < hi))),
                                       "high": int(np.sum(v >= hi))}
    return res


@T("rs_terrain_from_raster", "Slope and aspect from an elevation grid", RS, "table",
   [img("raster", IMAGE, "Elevation raster (single band)"), number("cell_size", 30.0, "Cell size (m)", min=0.1),
    boolean("hillshade", True, "Also compute a hillshade"), number("azimuth", 315.0, "Sun azimuth (deg)",
                                                            min=0.0), number("altitude", 45.0, "Sun altitude",
                                                                       min=0.0)],
   ex={"raster": IMAGE, "cell_size": 10.0}, up="r.slope.aspect / SAGA terrain analysis",
   tags=("terrain", "slope", "DEM"),
   summary="Slope, aspect and shading statistics for a DEM-like grid.")
def rs_terrain_from_raster(raster, cell_size=30.0, hillshade=True, azimuth=315.0, altitude=45.0):
    """Terrain metrics."""
    dem = _img(raster)
    gy, gx = np.gradient(dem, float(cell_size))
    slope = np.degrees(np.arctan(np.sqrt(gx ** 2 + gy ** 2)))
    aspect = (np.degrees(np.arctan2(-gy, gx)) + 360.0) % 360.0
    az = math.radians(float(azimuth))
    al = math.radians(float(altitude))
    shade = np.clip(np.cos(al) * np.cos(np.radians(slope))
                   + np.sin(al) * np.sin(np.radians(slope)) * np.cos(az - np.radians(aspect)), 0, 1)
    rows = [{"metric": "mean_slope_degrees", "value": round(float(np.mean(slope)), 4)},
           {"metric": "max_slope_degrees", "value": round(float(np.max(slope)), 4)},
           {"metric": "flat_fraction", "value": round(float(np.mean(slope < 5.0)), 4)},
           {"metric": "steep_fraction", "value": round(float(np.mean(slope > 30.0)), 4)},
           {"metric": "mean_aspect_degrees", "value": round(float(np.mean(aspect)), 3)},
           {"metric": "relief", "value": round(float(np.max(dem) - np.min(dem)) * 1.0, 4), "unit": "elevation units"},
           {"metric": "cell_size_m", "value": round(float(cell_size), 3)},
           {"metric": "mean_hillshade", "value": round(float(np.mean(shade)), 4) if hillshade else "n/a"}]
    res = table(rows, f"slope analysis on a {dem.shape[0]}x{dem.shape[1]} grid "
                   f"({round(float(np.mean(slope)), 2)} degrees mean)")
    res["stats"]["figure"] = plot.image_show(slope, title="slope (degrees)", cmap="magma")
    if hillshade:
        res["stats"]["hillshade"] = plot.image_show(shade, title="hillshade", cmap="Greys")
    return res


@T("rs_indices_time_series", "Index trend over a stack of scenes", RS, "table",
   [tbl("scenes", COUNTS, "Table with one column per scene"), textbox("value_column", "", "Column to trend"),
    intin("timestep", 16, "Days per scene", min=1), choice("trend", ["linear", "sen_like", "both"], "both",
                                                         "Trend estimator"), number("break_fraction", 0.2,
                                                             "Break-detection sensitivity", min=0.0)],
   ex={"scenes": COUNTS, "trend": "both"}, up="phenofit / TIMESAT", tags=("time series", "phenology"),
   summary="Per-scene index values with the seasonal trend and a break point.")
def rs_indices_time_series(scenes, value_column="", timestep=16, trend="both", break_fraction=0.2):
    """Index time series trend."""
    df = _tbl(scenes)
    numc = [c for c, _v in _num_cols(df)]
    cc = value_column if value_column in numc else (numc[0] if numc else None)
    if not cc:
        return table([], "no numeric column to trend")
    v = pd.to_numeric(df[cc], errors="coerce").fillna(0).to_numpy(dtype=float)
    t = np.arange(v.size, dtype=float) * max(1, int(timestep))
    slope, intercept = (np.polyfit(t, v, 1) if v.size > 2 else (0.0, float(v[0]) if v.size else 0.0))
    fit = float(intercept) + float(slope) * t
    resid = v - fit
    k = max(1, int(round(v.size * float(break_fraction))))
    cusum = np.cumsum(resid - np.mean(resid))
    bp = int(np.argmax(np.abs(cusum))) if cusum.size else 0
    rows = [{"scene": i + 1, "day": float(t[i]), "value": round(float(v[i]), 5),
            "fitted": round(float(fit[i]), 5), "residual": round(float(resid[i]), 5),
            "seasonal_position": round(float(v[i] - np.min(v)), 5) if v.size else 0.0}
           for i in range(v.size)]
    res = table(rows, f"{cc}: slope {round(float(slope) * max(1, int(timestep)), 6)} per year, "
                  f"break at scene {bp + 1}")
    res["stats"] = {"trend_per_step": round(float(slope), 8), "break_scene": bp + 1,
                   "figure": plot.line({cc: [float(x) for x in v], "fitted": [float(x) for x in fit]},
                                     x=[float(x) for x in t], xlabel="day", ylabel="index",
                                     title="index time series")}
    return res


@T("rs_image_change_detection", "Change detection between two images", RS, "table",
   [img("before", IMAGE, "Earlier image"), img("after", IMAGE, "Later image"), number("threshold", 10.0,
                                                            "Change magnitude threshold", min=0.0),
    choice("mode", ["difference", "ratio", "ndbi"], "difference", "Change metric"),
    boolean("figure", True, "Show the change map")],
   ex={"before": IMAGE, "after": IMAGE, "threshold": 5.0}, up="r.change.detect / CCDC",
   tags=("change detection", "remote sensing"),
   summary="Pixel-wise difference between two dates with area statistics.")
def rs_image_change_detection(before, after, threshold=10.0, mode="difference", figure=True):
    """Two-date change detection."""
    a = _img(before)
    b = _img(after)
    n = min(a.shape[0], b.shape[0]), min(a.shape[1], b.shape[1])
    a, b = a[: n[0], : n[1]], b[: n[0], : n[1]]
    if mode == "ratio":
        ch = (b + 1.0) / (a + 1.0)
    elif mode == "ndbi":
        ch = (b - a) / (b + a + 1e-9)
    else:
        ch = b - a
    cut = float(threshold) if mode != "ndbi" else float(threshold) / 100.0
    inc = ch > cut
    dec = ch < -cut
    rows = [{"metric": "changed_pixels", "value": int(np.sum(inc | dec))},
           {"metric": "increase_fraction", "value": round(float(np.mean(inc)), 5)},
           {"metric": "decrease_fraction", "value": round(float(np.mean(dec)), 5)},
           {"metric": "stable_fraction", "value": round(float(np.mean(~(inc | dec))), 5)},
           {"metric": "mean_change", "value": round(float(np.mean(ch)), 5)},
           {"metric": "max_increase", "value": round(float(np.max(ch)), 5)},
           {"metric": "max_decrease", "value": round(float(np.min(ch)), 5)},
           {"metric": "pixels", "value": int(ch.size)}]
    res = table(rows, f"{round(100 * float(np.mean(inc | dec)), 2)}% of {ch.size} pixels changed "
                  f"(|{mode}| > {threshold})")
    if figure:
        res["stats"]["figure"] = plot.image_show(ch, title=f"change ({mode})", cmap="coolwarm")
    return res
