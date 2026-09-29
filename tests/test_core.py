"""Round-trip tests for the engines in :mod:`chroma_titan.core`.

These are the primitives every tool is built on, so they are tested
independently of the registry.
"""

from __future__ import annotations

import math
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from chroma_titan.core import (align, genome, image, io, ml, motif, phylo, plot,  # noqa: E402
                               protein, seq, stats, tables)


# --------------------------------------------------------------------------- io
def test_fasta_round_trip():
    recs = io.parse_fasta(io.as_text("examples/genes.fa"))
    assert recs and all(r.seq and r.id for r in recs)
    again = io.parse_fasta(io.write_fasta(recs))
    assert [r.seq for r in again] == [r.seq for r in recs]


def test_fastq_parsing_and_quality():
    reads = io.parse_fastq(io.as_text("examples/reads_single.fastq"))
    assert reads and len(reads[0].seq) == len(reads[0].qual)
    q = [ord(c) - 33 for c in reads[0].qual]
    assert min(q) >= 0 and max(q) < 100


def test_bed_gff_vcf_and_sam_parsers():
    ivs = io.parse_bed(io.as_text("examples/regions.bed"))
    assert ivs and all(i.end >= i.start for i in ivs)
    gff = io.parse_gff(io.as_text("examples/annotation.gff"))
    assert gff and all(g.end >= g.start for g in gff)
    vcf = io.parse_vcf(io.as_text("examples/variants.vcf"))
    assert vcf.records and all(v.ref for v in vcf.records)
    assert "CHROM" in vcf.as_text()
    headers, records = io.parse_sam(io.as_text("examples/alignments.sam"))
    assert headers and str(headers[0]).startswith("@")
    assert records and isinstance(getattr(records[0], "flag", 0), int)


def test_resolve_accepts_example_names():
    assert Path(io.resolve("regions.bed")).name == "regions.bed"
    assert Path(io.resolve("@examples/regions.bed")).exists()


# --------------------------------------------------------------------- sequence
def test_gc_content_and_composition():
    assert seq.gc_content("GGCC") == 100.0
    assert seq.gc_content("ATAT") == 0.0
    assert sum(seq.composition("AATTCG").values()) == 6
    assert seq.composition("AATT")["A"] == 2


def test_reverse_complement_and_orfs():
    assert seq.reverse_complement("ATCG") == "CGAT"
    orfs = seq.find_orfs("ATGAAATTTGGGTAA", min_len=6)
    assert orfs and orfs[0]["start"] >= 0
    assert orfs[0]["protein"].startswith("M")


def test_translation_follows_the_genetic_code():
    assert seq.translate("ATGAAATTTGGATAA", to_stop=True) == "MKFG"
    frames = seq.translate_all_frames("ATGAAATTTGGATAA")
    assert len(frames) == 6


def test_digest_and_melting_temperature():
    pieces = seq.digest("GAATTCGGGAATTCAAA" * 2, "EcoRI")
    assert isinstance(pieces, list) and pieces
    assert math.isfinite(float(seq.melting_temp("ATGCATGCATGC")))


def test_codon_and_nucleotide_helpers():
    repeats = seq.find_repeats("ACGTACGTACGT")
    assert repeats and repeats[0]["period"] == 4 and repeats[0]["copies"] == 3
    obs_exp = seq.cpg_obs_exp("CGCGCGCGATATATAT")
    assert set(obs_exp) >= {"cg", "gc", "ratio"}


# ----------------------------------------------------------------------- genome
def test_interval_operations():
    a = io.parse_bed(io.as_text("examples/regions.bed"))
    b = io.parse_bed(io.as_text("examples/targets.bed"))
    assert isinstance(genome.intersect(a, b), list)
    assert all(m.end >= m.start for m in genome.merge(a))
    assert 0.0 <= float(genome.jaccard(a, b)["jaccard"]) <= 1.0
    assert isinstance(genome.closest(a, b), list)


def test_coverage_and_windows():
    ivs = io.parse_bed(io.as_text("examples/regions.bed"))
    cov = genome.coverage(ivs)
    assert cov and {"chrom", "start", "end", "coverage"} <= set(cov[0])
    wins = genome.window_counts(ivs, size=1000)
    assert wins and all(w["count"] >= 0 for w in wins)


def test_length_stats_and_summary_are_dicts():
    ivs = io.parse_bed(io.as_text("examples/regions.bed"))
    stats_row = genome.length_stats(ivs)
    assert isinstance(stats_row, dict)
    assert stats_row["count"] == len(ivs)
    assert stats_row["total_bp"] >= stats_row["max"]
    assert "N50" in genome.summarize_intervals(ivs)


def test_binnify_and_shuffle():
    _seqs, sizes = genome.read_genome("examples/genome.fa")
    assert sizes and all(v > 0 for v in sizes.values())
    bins = genome.binnify(sizes, 1000)
    assert bins and all(b.end > b.start for b in bins)
    ivs = io.parse_bed(io.as_text("examples/regions.bed"))
    shuffled = genome.shuffle_in({"chrV": 60000}, ivs, seed=7)
    assert len(shuffled) == len(ivs)
    assert all(s.end > s.start for s in shuffled)


# ----------------------------------------------------------------------- tables
def test_table_reading_and_text_round_trip():
    df = tables.load(io.as_text("examples/counts.tsv"))
    assert df.shape[0] >= 6
    again = tables.load(tables.as_text(df))
    assert list(again.columns) == list(df.columns)


def test_log_transform_and_normalisations():
    df = tables.load(io.as_text("examples/counts.tsv"))
    mat = df.drop(columns=[df.columns[0]])
    lg = tables.log_transform(mat)
    assert lg.shape == mat.shape and float(lg.max().max()) <= float(mat.max().max())
    per_library = tables.cpm(mat).sum(axis=0)  # each sample column is scaled to one million
    assert np.allclose(per_library, 1e6)
    z = tables.zscore(mat, axis=1)
    assert abs(float(np.mean(np.asarray(z.iloc[0], dtype=float)))) < 1e-6


def test_groupby_and_tidy_helpers():
    df = tables.load(io.as_text("examples/phenotypes.tsv"))
    out = tables.groupby_agg(df, by=["group"], funcs={"yield": "mean"})
    assert not out.empty and len(out) == df["group"].nunique()
    long = tables.melt(df, id_vars=["group"], value_vars=["yield"])
    assert {"group", "value"} <= set(long.columns)


# ------------------------------------------------------------------------ stats
def test_basic_statistics_are_sane():
    st = stats.ttest_ind([1, 2, 3, 4], [5, 6, 7, 8])
    assert st["p_value"] < 0.05 and st["statistic"] < 0
    assert stats.pearson([1, 2, 3], [2, 4, 6])["r"] == pytest.approx(1.0, abs=1e-4)
    assert 0 <= stats.mann_whitney_u([1, 2, 3], [4, 5, 6])["p_value"] <= 1
    fit = stats.linear_regression([1, 2, 3], [2.0, 4.0, 6.0])
    assert fit["slope"] == pytest.approx(2.0, abs=1e-6)


def test_p_adjust_and_fisher():
    adj = stats.p_adjust([0.01, 0.04, 0.5], "fdr_bh")
    assert len(adj) == 3 and adj[-1] >= adj[0]
    odds, p = stats.fisher_exact(5, 1, 1, 5)
    assert odds > 1 and 0 <= p <= 1


def test_diversity_indices():
    counts = [10, 8, 6, 4, 2, 1, 1]
    assert stats.shannon2(counts) > 0
    assert 0 <= stats.simpson(counts) <= 1
    assert stats.chao1(counts) >= sum(1 for c in counts if c > 0)
    rare = stats.rarefaction(counts, steps=4, seed=1)
    assert rare and {"size", "richness"} <= set(rare[0])


def test_correlation_matrix_and_describe():
    df = tables.load(io.as_text("examples/counts.tsv")).drop(columns=["gene"])
    cols = {c: df[c].to_numpy(dtype=float).tolist() for c in df.columns[:3]}
    corr = stats.correlation_matrix(cols)
    assert corr and len(corr) >= 2
    desc = stats.describe({"a": [1.0, 2.0, 3.0]})
    assert desc and isinstance(desc[0], dict)


# ----------------------------------------------------------------- align / phylo
def test_pairwise_alignment_and_identity():
    res = align.align_pairwise("ACGTACGT", "ACGTACGA")
    assert res["identity"] > 40 and res["aligned_length"] >= 8
    assert align.hamming("ACGT", "ACGA")["distance"] == 1


def test_distance_matrix_is_symmetric():
    D = np.asarray(align.distance_matrix(["ACGT", "ACGA", "AGGT"], model="p"), dtype=float)
    assert D.shape == (3, 3)
    assert np.allclose(D, D.T)
    assert list(np.diag(D)) == [0.0, 0.0, 0.0]


def test_newick_round_trip_and_distances():
    text = "((a:0.1,b:0.2)c:0.3,(d:0.4,e:0.2)f:0.1)root;"
    tree = phylo.parse_newick(text)
    assert io.write_newick(tree) == text
    tips = phylo.newick_leaves(tree)
    assert tips == ["a", "b", "d", "e"]
    D = phylo.patristic_matrix(tree)
    assert D[tips.index("a")][tips.index("b")] == pytest.approx(0.3)
    assert D[tips.index("a")][tips.index("d")] == pytest.approx(0.9)
    assert phylo.root_to_tip_depths(tree)["d"] == pytest.approx(0.5)


def test_neighbour_joining_and_upgma_produce_trees():
    names = ["a", "b", "c", "d"]
    D = [[0, .2, .4, .6], [.2, 0, .4, .6], [.4, .4, 0, .6], [.6, .6, .6, 0]]
    for tree in (phylo.neighbor_joining(names, D), phylo.upgma(names, D)):
        assert len(phylo.newick_leaves(tree)) == 4
        assert phylo.tree_length(tree) > 0


def test_tree_operations_are_safe():
    tree = phylo.parse_newick("((a:0.1,b:0.2)c:0.3,(d:0.4,e:0.2)f:0.1)root;")
    assert sorted(phylo.newick_leaves(phylo.root_at_midpoint(tree))) == ["a", "b", "d", "e"]
    assert sorted(phylo.newick_leaves(phylo.reroot_at_leaf(tree, "a"))) == ["a", "b", "d", "e"]
    assert phylo.ladderize(tree) is not None
    assert sorted(phylo.newick_leaves(phylo.prune_tree(tree, ["a", "b", "d"]))) == ["a", "b", "d"]
    assert phylo.robinson_foulds(tree, phylo.parse_newick("((a,b),(d,e));"))["RF_distance"] == 0.0


def test_bootstrap_and_consensus_run():
    aln = [{"id": n, "seq": s} for n, s in
           [("a", "ACGTT"), ("b", "ACATT"), ("c", "AGGTT"), ("d", "AGGAC"), ("e", "AGGAT")]]
    trees = phylo.bootstrap_alignment(aln, n=9, seed=1, method="nj", model="p")
    assert len(trees) == 9
    freqs = phylo.clade_frequencies(trees)
    assert freqs and max(freqs.values()) <= 9
    consensus = phylo.majority_rule_consensus(trees, cutoff=0.5)
    assert consensus is not None


# ----------------------------------------------------------------------- protein
def test_protein_properties():
    chain = "MKVLMFVLAALAGVSSAAADAEIKKLRVLIKQEEQKLVESRELEKLKLENRL"
    assert protein.clean(chain) == chain
    assert 4000 < protein.protein_mass(chain) < 20000
    assert 0 < protein.isoelectric_point(chain) < 14
    assert abs(protein.net_charge(chain, pH=1.0)) > abs(protein.net_charge(chain, pH=13.0))
    assert protein.molecular_formula(chain).startswith("C")
    assert 0 <= protein.gravy(chain) <= 5 or protein.gravy(chain) < 0


def test_protein_annotations():
    chain = "MKVLMFVLAALAGVSSAAADAEIKKLRVLIKQEEQKLVESRELEKLKLENRL"
    assert 0.0 <= protein.signal_peptide_score(chain)["score"] <= 1.0
    assert protein.transmembrane_helices(chain, window=9, cutoff=1.0)
    assert len(protein.disorder_propensity(chain)) == len(chain)
    assert isinstance(protein.prosite_sites(chain), list)
    hist = protein.charge_histogram(chain, pH_range=(4.0, 10.0), steps=5)
    assert len(hist) == 5 and "pH" in hist[0]


# ------------------------------------------------------------------------- image
def test_image_pipeline():
    arr = image.read_image("examples/cells.pgm")
    gray = image.to_gray(arr)
    assert gray.ndim == 2 and gray.size == 48 * 48
    thr = image.threshold_otsu(gray)["threshold"]
    mask = image.binarize(gray, value=float(thr), above=True)
    labelled = image.label(mask)["labels"]
    assert labelled.shape == gray.shape
    props = image.regionprops(labelled, gray)
    if props:
        assert "area" in props[0] and props[0]["area"] > 0
    assert image.gaussian_blur(gray, sigma=1.0).shape == gray.shape
    assert np.isfinite(np.asarray(image.sobel(gray)["magnitude"], dtype=float)).all()
    assert image.resize_nearest(gray, scale=0.5).shape == (24, 24)


# ------------------------------------------------------------ motif / plot / ml
def test_motif_matrix_pipeline():
    pfm = motif.read_pfm(io.as_text("examples/motif.pfm"))
    assert pfm and len(pfm) >= 4 and len(pfm[0]) == 4  # one row per position, four bases
    pwm = motif.pfm_to_pwm(pfm)
    assert len(pwm) == len(pfm)
    consensus = motif.pfm_to_consensus(pfm)
    assert consensus and set(consensus) <= set("ACGTRYKMSWBN")
    assert motif.total_information(pfm) >= 0


def test_plots_render_to_png_bytes():
    data = io.figure_bytes(plot.line({"a": [1, 2, 3]}, title="t"))
    assert data[:4] == b"\x89PNG"
    heat = io.figure_bytes(plot.heatmap([[1.0, 2.0], [3.0, 4.0]], row_labels=["x", "y"],
                                      col_labels=["c1", "c2"]))
    assert heat[:4] == b"\x89PNG"


def test_ml_helpers():
    X = [[0.0, 0.0], [0.1, 0.1], [5.0, 5.0], [5.1, 4.9]]
    km = ml.kmeans(X, k=2, seed=1)
    assert len(km["labels"]) == 4 and len(set(km["labels"])) <= 2
    pca = ml.pca(X, components=2)
    assert np.asarray(pca["scores"]).shape[0] == 4
    Z = ml.linkage(X, method="average")
    assert ml.fcluster(Z, t=1.0, criterion="distance").count(1) >= 1
