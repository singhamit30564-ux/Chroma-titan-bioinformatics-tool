"""Smoke tests for the Streamlit panel app (skipped when Streamlit is absent)."""

from __future__ import annotations

import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

st = pytest.importorskip("streamlit")
from streamlit.testing.v1 import AppTest  # noqa: E402

APP = str(ROOT / "app.py")
TOOLS = ["fasta_gc_content", "bedtools_intersect", "volcano_plot", "image_threshold_cells",
         "phylo_nj_from_alignment", "chem_smiles_descriptors", "bcftools_query", "deseq_lite"]


@pytest.fixture(scope="module")
def app():
    return AppTest.from_file(APP, default_timeout=240).run()


def _nav(app):
    """The sidebar page picker (other radios belong to tool file inputs)."""
    return [r for r in app.radio if r.label == "Browse"][0]


def test_app_starts_on_the_tool_panel(app):
    assert not app.exception
    assert app.session_state["page"] == "Tool panel"
    assert _nav(app).value == "Tool panel"
    assert "Group" in [s.label for s in app.selectbox]


@pytest.mark.parametrize("tool_id", TOOLS)
def test_tool_opens_runs_and_shows_output(app, tool_id):
    app.session_state["jump"] = tool_id
    app.run()
    assert not app.exception, f"{tool_id}: {app.exception}"
    assert app.title[0].value  # a tool page has a title
    buttons = [b for b in app.button if b.label == "Run tool"]
    assert buttons, f"{tool_id}: no run button"
    app = buttons[0].click().run()
    assert not app.error, f"{tool_id}: {app.error[0].value if app.error else ''}"
    result = (app.session_state.get("results") or {}).get(tool_id)
    assert result, f"{tool_id}: nothing was stored for this run"


def test_all_tools_and_datasets_pages(app):
    for page in ("All tools", "Datasets", "About this panel"):
        _nav(app).set_value(page).run()
        assert not app.exception, page
        assert not app.error, f"{page}: {app.error[0].value if app.error else ''}"
    _nav(app).set_value("Datasets").run()
    assert app.dataframe, "the datasets page should list the catalogue"
    assert app.code, "the datasets page should preview a file"
