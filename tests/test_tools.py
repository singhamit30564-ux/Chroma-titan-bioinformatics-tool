"""Registry, panel and tool-execution tests for Chroma Titan.

``pytest tests/test_tools.py`` runs the fast structural checks plus a sample of
tools; ``pytest -m slow`` additionally executes the bundled example of *every*
registered tool (the same thing ``python -m chroma_titan check`` does).
"""

from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from chroma_titan import datasets, panel, registry  # noqa: E402
from chroma_titan.core import io  # noqa: E402

ALL_TOOLS = registry.all_tools()
IDS = sorted(t.id for t in ALL_TOOLS)


# ---------------------------------------------------------------------------
# the panel itself
# ---------------------------------------------------------------------------
def test_at_least_300_tools_are_registered():
    assert len(ALL_TOOLS) >= 300, f"only {len(ALL_TOOLS)} tools registered"


def test_tool_ids_are_unique_and_bound_to_a_function():
    import re as _re

    seen = set()
    for t in ALL_TOOLS:
        assert t.id not in seen, f"duplicate tool id {t.id}"
        seen.add(t.id)
        assert _re.fullmatch(r"[A-Za-z0-9_]+", t.id), f"tool id {t.id!r} must be alphanumeric"
        assert callable(t.fn), f"{t.id} has no implementation"
        assert t.fn.__module__.startswith("chroma_titan.tools"), f"{t.id}: wrong module"


def test_every_section_is_a_real_galaxy_section():
    for t in ALL_TOOLS:
        assert t.section in panel.SECTION_LABELS, f"{t.id} uses unknown section {t.section!r}"


def test_every_panel_section_has_at_least_one_tool():
    used = {t.section for t in ALL_TOOLS}
    missing = [s for s in panel.SECTION_LABELS if s not in used]
    assert not missing, f"sections without tools: {missing}"


def test_group_of_is_consistent_with_the_all_tools_layout():
    for label, secs in panel.ALL_TOOLS_GROUPS:
        for sec in secs:
            assert panel.group_of(sec) == label
    assert panel.group_of("not_a_section") == "Other Tools"


def test_panel_views_only_reference_known_sections():
    for view, spec in panel.PANEL_VIEWS.items():
        groups = spec["groups"] or {}
        for label, secs in groups.items():
            assert label, f"{view} has an empty group name"
            for sec in secs:
                assert sec in panel.SECTION_LABELS, f"{view}/{label}: unknown section {sec}"


def test_panel_layout_deduplicates_tools():
    layout = panel.panel_layout("all_tools")
    flat = [tid for _g, rows in layout for _s, ids in rows for tid in ids]
    assert len(flat) == len(set(flat)) == registry.count()


def test_search_finds_tools_by_keyword():
    assert registry.search("coverage")
    assert registry.search("bed intersect")
    hits = registry.search("motif")
    assert all("motif" in (t.id + t.name + t.summary + " ".join(t.tags)).lower() for t in hits)
    assert registry.search("definitely-not-a-tool") == []


# ---------------------------------------------------------------------------
# tool declarations
# ---------------------------------------------------------------------------
def test_tool_declarations_are_complete():
    for t in ALL_TOOLS:
        assert t.name and t.summary, f"{t.id} needs a name and a summary"
        assert t.output in registry.OUTPUTS, f"{t.id}: bad output kind {t.output!r}"
        assert t.example, f"{t.id} must ship a runnable example"
        assert t.version and isinstance(t.tags, tuple)


def test_examples_reference_bundled_data():
    for t in ALL_TOOLS:
        for key, value in t.example.items():
            if isinstance(value, str) and value.startswith("examples/"):
                target = Path(io.EXAMPLE_DIR) / value.split("examples/")[-1]
                assert target.exists(), f"{t.id}.{key}: example file {value} is missing"


def test_inputs_have_labels_and_valid_kinds():
    for t in ALL_TOOLS:
        for i in t.inputs:
            assert i.kind in registry.KINDS, f"{t.id}.{i.name}: unknown kind {i.kind}"
            assert i.label, f"{t.id}.{i.name}: every input needs a label"
            assert i.name != "fn"


def test_widget_defaults_are_renderable():
    """The Streamlit form is built from these values, so they must be coherent."""
    pool = {e["path"] for e in datasets.entries() if e["present"]}
    for t in ALL_TOOLS:
        for i in t.inputs:
            if i.kind in ("choice", "multi"):
                if not i.options:
                    continue  # the app falls back to a free-text field
                if i.kind == "choice" and i.default is not None:
                    assert i.default in i.options, f"{t.id}.{i.name}: default outside options"
                if i.kind == "multi" and isinstance(i.default, (list, tuple)):
                    stray = [d for d in i.default if d not in i.options]
                    assert not stray, f"{t.id}.{i.name}: defaults {stray[:2]} outside options"
            if i.kind in ("int", "number", "slider") and i.default is not None:
                value = float(i.default)
                if i.min is not None:
                    assert value >= float(i.min), f"{t.id}.{i.name}: default below min"
                if i.max is not None:
                    assert value <= float(i.max), f"{t.id}.{i.name}: default above max"
            if i.kind == "file" and str(i.default or "").startswith("examples/"):
                assert i.default in pool, f"{t.id}.{i.name}: {i.default} is not bundled"


def test_to_dict_is_json_serialisable():
    for t in ALL_TOOLS[:40]:
        json.dumps(t.to_dict(), default=str)


# ---------------------------------------------------------------------------
# execution
# ---------------------------------------------------------------------------
SAMPLE = [
    "fasta_gc_content", "bedtools_intersect", "bcftools_view", "deseq_lite", "volcano_plot",
    "presto_aucell", "phylo_nj_from_alignment", "image_threshold_cells", "proteo_peptide_properties",
    "chem_smiles_descriptors", "pharm_ic50_panel", "multiomics_integrated_score", "datamash_column_stats",
    "mothur_count_seqs", "str_locus_statistics", "hic_insulation_score", "counts_feature_stats",
    "samtools_feature_counts", "gff3_feature_counts", "climate_growing_degree_days",
]


def _result_shape_ok(res) -> bool:
    if hasattr(res, "savefig"):  # a bare matplotlib figure
        return True
    if isinstance(res, pd.DataFrame):
        return True
    if isinstance(res, dict):
        return any(k in res for k in ("table", "text", "image", "stats", "message"))
    return bool(res)


@pytest.mark.parametrize("tool_id", SAMPLE)
def test_sample_tools_produce_results(tool_id):
    t = registry.get_tool(tool_id)
    assert tool_id in registry.REGISTRY, f"{tool_id} is not registered"
    res = registry.run(tool_id, **t.example)
    assert res, f"{tool_id} returned nothing"
    assert _result_shape_ok(res), f"{tool_id}: unexpected result shape {type(res)}"
    if isinstance(res, dict) and "table" in res:
        assert isinstance(res["table"], pd.DataFrame)


def test_running_a_tool_by_id_coerces_loose_values():
    res = registry.run("fasta_gc_content", src="examples/genes.fa")
    assert "table" in res
    assert len(res["table"]) == len(io.parse_fasta(io.as_text("examples/genes.fa")))


@pytest.mark.slow
def test_cli_check_runs_the_whole_registry():
    out = _cli("check")
    tail = [ln for ln in out.splitlines() if "ok," in ln]
    assert tail, out[-800:]
    assert " 0 failed" in tail[-1], out[-3000:]
    assert int(tail[-1].split()[0]) >= 300


def test_run_rejects_unknown_tool():
    with pytest.raises(KeyError):
        registry.get_tool("no_such_tool_here")


def test_registry_rejects_unknown_input_names():
    from chroma_titan.registry import Input, tool

    with pytest.raises(ValueError, match="no Input declaration"):

        @tool(id="bad_tool", name="Bad", section="bed",
              inputs=[Input("other", kind="text", default="x")])
        def bad_tool(required_arg):  # no Input for required_arg
            return required_arg


@pytest.mark.slow
def test_every_example_runs():
    import matplotlib.pyplot as plt

    failures = []
    for tid in IDS:
        t = registry.get_tool(tid)
        try:
            res = registry.run(tid, **t.example)
            if res is None or (isinstance(res, dict) and not res):
                failures.append(f"{tid}: empty result")
        except Exception as exc:  # noqa: BLE001
            failures.append(f"{tid}: {type(exc).__name__}: {exc}")
        finally:
            plt.close("all")
    assert not failures, "failed tools:\n" + "\n".join(failures)


# ---------------------------------------------------------------------------
# datasets catalogue
# ---------------------------------------------------------------------------
def test_datasets_catalogue_matches_disk():
    assert datasets.count() == len(datasets.CATALOGUE)
    on_disk = set(datasets.on_disk())
    assert set(datasets.names()) == on_disk, "data/examples and datasets.CATALOGUE are out of sync"


def test_datasets_load_returns_the_right_shape():
    df = datasets.load("counts.tsv")
    assert isinstance(df, pd.DataFrame) and df.shape[0] >= 10 and "sample_A" in df.columns
    recs = datasets.load("genes.fa")
    assert {"id", "sequence", "length"} <= set(recs.columns)
    assert datasets.load("cells.pgm").size > 0
    assert datasets.preview("regions.bed").count("\n") > 0


def test_datasets_tools_using():
    users = datasets.tools_using("counts.tsv")
    assert users and all(u in registry.REGISTRY for u in users)


# ---------------------------------------------------------------------------
# command line
# ---------------------------------------------------------------------------
def _cli(*args: str) -> str:
    out = subprocess.run([sys.executable, "-m", "chroma_titan", *args], cwd=ROOT,
                        capture_output=True, text=True, timeout=900)
    assert out.returncode == 0, out.stderr[-2000:]
    return out.stdout


def test_cli_count_and_sections():
    assert int(_cli("count").strip()) == len(ALL_TOOLS)
    sections = _cli("sections")
    assert "total tools in" in sections
    assert "rna_seq" in sections


def test_cli_info_and_run():
    info = _cli("info", "fasta_gc_content")
    assert "fasta_gc_content" in info and "inputs" in info
    payload = json.loads(_cli("run", "fasta_gc_content", "--src=examples/genes.fa", "--json"))
    assert payload["table"] or payload["message"]


def test_cli_search_and_examples():
    assert "bedtools" in _cli("search", "intersect")
    listing = _cli("examples")
    for name in ("genome.fa", "counts.tsv", "cells.pgm", "variants.vcf"):
        assert name in listing


@pytest.mark.slow
def test_cli_check_runs_the_whole_registry():
    out = _cli("check")
    tail = [ln for ln in out.splitlines() if "ok," in ln]
    assert tail, out[-800:]
    assert " 0 failed" in tail[-1], out[-3000:]
    assert int(tail[-1].split()[0]) >= 300
