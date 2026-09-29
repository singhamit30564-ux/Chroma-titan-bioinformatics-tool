"""Chroma Titan -- a Galaxy-style bioinformatics tool panel in the browser.

Run it with::

    streamlit run app.py

The sidebar reproduces the tool panel of https://usegalaxy.org/ (the same groups
and sections, read from the usegalaxy-playbook panel configuration), the main
panel renders a tool form from the tool's declared inputs, executes the tool
function and shows the datasets it returns -- tables, figures, text and stats --
with download buttons, exactly like a Galaxy job output.

The tools themselves live in :mod:`chroma_titan.tools` and are plain Python
functions, so the very same tool can be run from the command line::

    python -m chroma_titan run fasta_gc_content --src examples/genome.fa
    python -m chroma_titan check          # run every tool's bundled example
"""

from __future__ import annotations

import json
import traceback
from typing import Any

import pandas as pd
import streamlit as st

from chroma_titan import datasets, panel, registry
from chroma_titan.core import io as cio

# ---------------------------------------------------------------------------
# page setup
# ---------------------------------------------------------------------------
st.set_page_config(page_title="Chroma Titan - Galaxy-style tool panel",
                   page_icon="🧬", layout="wide", initial_sidebar_state="expanded")

PAGES = ["Tool panel", "All tools", "Datasets", "About this panel"]
EXT_BY_FMT = {
    "fasta": [".fa", ".fasta", ".faa", ".fna", ".txt"],
    "fastq": [".fastq", ".fq", ".txt"],
    "bed": [".bed", ".narrowPeak", ".txt"],
    "gff": [".gff", ".gff3", ".gtf", ".txt"],
    "vcf": [".vcf", ".bcf", ".txt"],
    "sam": [".sam", ".bam", ".txt"],
    "tsv": [".tsv", ".csv", ".txt", ".tab"],
    "pgm": [".pgm", ".ppm", ".png", ".jpg", ".txt"],
    "txt": [".txt", ".tsv", ".csv", ".fa", ".fasta", ".bed", ".gff", ".vcf", ".sam",
            ".fastq", ".nwk", ".newick", ".pgm", ".gmt", ".obo", ".pairs", ".cxcg",
            ".bedgraph", ".pfm", ".json", ".md", ".fasta"],
}
EXAMPLES = [e["path"] for e in datasets.entries() if e["present"]]


def _init_state() -> None:
    ss = st.session_state
    ss.setdefault("page", "Tool panel")
    ss.setdefault("tool_id", None)
    ss.setdefault("history", [])
    ss.setdefault("results", {})
    ss.setdefault("errors", {})
    ss.setdefault("query", "")
    ss.setdefault("view", "all_tools")


def _remember(tool_id: str) -> None:
    hist = [h for h in st.session_state["history"] if h != tool_id]
    st.session_state["history"] = ([tool_id] + hist)[:12]


def _fig_to_bytes(obj: Any) -> bytes | None:
    """Render a matplotlib figure (or accept PNG bytes) for ``st.image``."""
    if isinstance(obj, (bytes, bytearray)):
        return bytes(obj)
    if hasattr(obj, "savefig"):
        try:
            data = cio.figure_bytes(obj)
            cio.close(obj)
            return data
        except Exception:  # noqa: BLE001 - a broken figure must not break the page
            return None
    return None


# ---------------------------------------------------------------------------
# sidebar: the tool panel
# ---------------------------------------------------------------------------
def sidebar() -> None:
    ss = st.session_state
    with st.sidebar:
        st.markdown("### 🧬 Chroma Titan")
        st.caption(f"{registry.count()} tools · {len(panel.SECTION_LABELS)} panel sections · "
                   "layout mirrors usegalaxy.org")
        page = st.radio("Browse", PAGES, index=PAGES.index(ss["page"] or "Tool panel"),
                        key="_nav")
        if page != ss["page"]:
            ss["page"] = page
        if page != "Tool panel":
            st.caption("Pick a page above; the tool panel itself is on the first page.")
            return

        st.divider()
        view_names = {v: panel.PANEL_VIEWS[v]["name"] for v in panel.view_ids()}
        jump = ss.pop("jump", None)
        if jump and jump in registry.REGISTRY and registry.get_tool(jump).group not in \
                [g for g, _s in panel.panel_layout(ss["view"])]:
            ss["view"] = "all_tools"
        chosen = st.selectbox("Panel view", list(view_names),
                              index=list(view_names).index(ss["view"]),
                              format_func=lambda k: view_names[k], key="_view")
        ss["view"] = chosen
        query = st.text_input("Search tools", value=ss["query"], key="_query",
                              placeholder="coverage, bed, motif, NDVI…",
                              help="Matches tool id, name, summary, tags and section")
        ss["query"] = query

        layout = panel.panel_layout(chosen)
        if jump and jump in registry.REGISTRY:
            _jump_to(jump, layout)
        if query.strip():
            hits = registry.search(query, limit=200)
            st.caption(f"{len(hits)} tool(s) match “{query}”")
            if not hits:
                st.info("No tool matches that search.")
            for t in hits[:60]:
                if st.button(t.name, key=f"q_{t.id}", help=f"{t.section_label} · {t.id}",
                             width="stretch"):
                    st.session_state["jump"] = t.id
                    ss["query"] = ""
                    st.rerun()
            st.divider()
        if not layout:
            st.info("This panel view has no sections yet.")
            return
        group = st.selectbox("Group", [g for g, _s in layout], key="_group")
        secs = dict(layout)[group]
        labels = [f"{sl}  ({len(ids)})" for sl, ids in secs]
        pick = st.selectbox("Section", labels, key="_section", format_func=lambda x: x)
        if pick not in labels:
            pick = labels[0]
        seclabel, ids = secs[labels.index(pick)]
        options = [f"{_tool_name(i)}  ·  {i}" for i in ids]
        sel = st.selectbox("Tool", options, key="_tool", format_func=lambda x: x)
        if sel not in options:
            sel = options[0]
        chosen_id = ids[options.index(sel)]
        if chosen_id != ss["tool_id"]:
            ss["tool_id"] = chosen_id
            ss["results"].pop(chosen_id, None)
            ss["errors"].pop(chosen_id, None)
        st.caption(f"Section “{seclabel}”: {len(ids)} tools")
        if ss["history"]:
            st.divider()
            st.caption("Recently used")
            for h in ss["history"]:
                if h in registry.REGISTRY and st.button(_tool_name(h), key=f"hist_{h}", help=h,
                                                        width="stretch"):
                    ss["tool_id"] = h
                    st.rerun()


def _jump_to(tool_id: str, layout: list[tuple[str, list[tuple[str, list[str]]]]]) -> None:
    """Point the panel selects at ``tool_id`` (used by search and 'open' buttons)."""
    t = registry.get_tool(tool_id)
    ss = st.session_state
    for label, secs in layout:
        for seclabel, ids in secs:
            if tool_id in ids:
                ss["_group"], ss["_section"], ss["_tool"] = label, f"{seclabel}  ({len(ids)})", \
                    f"{t.name}  ·  {tool_id}"
                ss["tool_id"] = tool_id
                return
    ss["tool_id"] = tool_id


def _tool_name(tool_id: str) -> str:
    try:
        return registry.get_tool(tool_id).name
    except KeyError:  # pragma: no cover
        return tool_id


def _group_matches(group: tuple[str, list[tuple[str, list[str]]]], query: str) -> bool:
    if not query.strip():
        return True
    q = query.lower()
    _label, secs = group
    for _sl, ids in secs:
        for tid in ids:
            t = registry.REGISTRY.get(tid)
            if t and all(w in " ".join([t.id, t.name, t.summary, " ".join(t.tags)]).lower()
                         for w in q.split()):
                return True
    return True  # keep the tree navigable; the search box lists the hits directly


# ---------------------------------------------------------------------------
# main: the tool form (a Galaxy "tool page")
# ---------------------------------------------------------------------------
def tool_page() -> None:
    ss = st.session_state
    tid = ss["tool_id"]
    if not tid or tid not in registry.REGISTRY:
        st.title("🧬 Chroma Titan")
        st.markdown(
            f"""**{registry.count()} tools** organised exactly like the
            [usegalaxy.org](https://usegalaxy.org/) tool panel.

            Pick a *group → section → tool* in the sidebar to open a tool form, or search for
            a keyword such as `coverage`, `motif`, `VCF`, `UMAP`, `NDVI`.

            Every tool is a plain Python function registered with
            :func:`chroma_titan.registry.tool`; bundled example data in ``data/examples``
            means most tools run with one click.""")
        _quick_links()
        return
    t = registry.get_tool(tid)
    _remember(tid)
    st.markdown(f"`{t.group}` › `{t.section_label}`")
    st.title(t.name)
    cols = st.columns([3, 1])
    with cols[0]:
        st.markdown(t.summary or "")
        if t.upstream:
            st.caption(f"Galaxy-equivalent wrapper: **{t.upstream}**")
    with cols[1]:
        st.metric("tool id", t.id, help="python -m chroma_titan run " + t.id)
    tags = list(t.tags) + [t.output, t.version]
    st.caption(" · ".join(f"`{x}`" for x in tags if x))

    left, right = st.columns([2.0, 1.25], gap="large")
    with left:
        values = _render_form(t)
    with right:
        _render_output(t)

    with st.expander("Tool declaration (what the form is generated from)"):
        st.dataframe(_inputs_frame(t), width="stretch", hide_index=True)
        st.code(json.dumps(t.to_dict(), indent=2, default=str), language="json")
    with st.expander("Run from the command line"):
        st.code(_cli_command(t, values), language="bash")


def _quick_links() -> None:
    st.subheader("Jump to a tool")
    samples = [tid for tid in ("fasta_gc_content", "bedtools_intersect", "bcftools_view", "scanpy_umap",
                              "phylo_nj_from_alignment", "image_threshold_cells", "proteo_peptide_properties",
                              "volcano_plot") if tid in registry.REGISTRY]
    found = [s for s in samples if s in registry.REGISTRY]
    cols = st.columns(min(4, max(1, len(found))))
    for i, tid in enumerate(found):
        with cols[i % len(cols)]:
            if st.button(f"{_tool_name(tid)}", key=f"quick_{tid}", width="stretch"):
                st.session_state["tool_id"] = tid
                st.session_state["jump"] = tid
                st.rerun()
    st.divider()
    counts = pd.DataFrame(_section_counts(), columns=["group", "section", "tools"])
    st.dataframe(counts, width="stretch", hide_index=True, height=280)


def _section_counts() -> list[list[Any]]:
    rows: dict[tuple[str, str], int] = {}
    for t in registry.all_tools():
        key = (t.group, t.section_label)
        rows[key] = rows.get(key, 0) + 1
    return [[g, s, n] for (g, s), n in sorted(rows.items())]


def _inputs_frame(t: Any) -> pd.DataFrame:
    return pd.DataFrame([{
        "id": i.name, "kind": i.kind, "label": i.label or i.name.replace("_", " ").title(),
        "default": _short(i.default), "options": ", ".join(map(str, i.options)) if i.options else "",
        "range": f"{i.min}…{i.max}" if i.min is not None or i.max is not None else "",
        "help": i.help,
    } for i in t.inputs])


def _short(v: Any, width: int = 46) -> str:
    if v is None:
        return ""
    s = str(v).replace("\n", "\\n")
    return s if len(s) <= width else s[: width - 1] + "…"


def _render_form(t: Any) -> dict:
    """Widgets for every declared input; returns the collected values."""
    st.subheader("Inputs")
    values: dict[str, Any] = {}
    for inp in t.inputs:
        values[inp.name] = _widget(t, inp)
    col_a, col_b, _rest = st.columns([1, 1, 3])
    run = col_a.button("Run tool", type="primary", icon="🚀", width="stretch")
    reset = col_b.button("Reset to example", icon="🔄", width="stretch",
                        disabled=not bool(t.example))
    if reset:
        _reset_to_example(t)
        st.rerun()
    if run:
        _run(t, values)
    return values


def _reset_to_example(t: Any) -> None:
    """Throw away edited widget values and put the bundled example back."""
    ss = st.session_state
    for inp in t.inputs:
        base = f"{t.id}__{inp.name}"
        example = t.example.get(inp.name)
        if inp.kind == "file":
            example_text = (isinstance(example, str) and example not in EXAMPLES
                            and bool(example))
            ss[base + "_mode"] = "Paste text" if example_text else "Bundled example"
            ss[base + "_up"] = None
            if example_text:
                ss[base + "_paste"] = str(example)
            if isinstance(example, str) and example in EXAMPLES:
                ss[base + "_sel"] = example
            else:
                ss.pop(base + "_sel", None)
            continue
        if example is None:
            ss.pop(base, None)
            if inp.kind == "multifile":
                ss.pop(base + "_multi", None)
            continue
        ss[base] = _cast_example(inp, example)
    st.toast("Inputs reset to the bundled example")


def _cast_example(inp: Any, value: Any) -> Any:
    try:
        if inp.kind == "bool":
            return bool(value)
        if inp.kind == "int":
            return int(float(value))
        if inp.kind in ("number", "slider"):
            return float(value)
        if inp.kind == "multi":
            opts = list(inp.options)
            vals = value if isinstance(value, (list, tuple)) else [value]
            return [v for v in vals if v in opts] or opts[:1]
        if inp.kind == "choice":
            return value if value in list(inp.options) else inp.options[0]
        if inp.kind in ("multifile",):
            return []
        return str(value)
    except (TypeError, ValueError, IndexError):
        return inp.default


def _widget(t: Any, inp: Any) -> Any:
    """One Streamlit widget per declared Galaxy-style input."""
    key = f"{t.id}__{inp.name}"
    label = inp.label or inp.name.replace("_", " ").title()
    help_txt = inp.help or (f"default: {_short(inp.default, 120)}" if inp.default not in (None, "") else "")
    if inp.kind == "file":
        return _file_widget(label, inp, key, help_txt)
    if inp.kind == "multifile":
        opts = EXAMPLES
        return st.multiselect(label, opts, default=_default_list(inp, opts),
                              key=key + "_multi", help=help_txt or "Pick one or more bundled files")
    if inp.kind == "code":
        return st.text_area(label, value=str(inp.default or ""), height=170, key=key, help=help_txt)
    if inp.kind == "text":
        return st.text_input(label, value=str(inp.default or ""), key=key, help=help_txt)
    if inp.kind == "bool":
        return st.checkbox(label, value=bool(inp.default), key=key, help=help_txt)
    if inp.kind in ("choice",):
        opts = list(inp.options)
        if not opts:
            return st.text_input(label, value=str(inp.default or ""), key=key, help=help_txt)
        default = inp.default if inp.default in opts else opts[0]
        return st.selectbox(label, opts, index=opts.index(default), key=key, help=help_txt)
    if inp.kind == "multi":
        opts = list(inp.options)
        if not opts:
            # the tool reads the column/key names from the data, so offer free text
            default = inp.default if isinstance(inp.default, (list, tuple)) else []
            return st.text_input(label, value=", ".join(str(d) for d in default), key=key,
                                help=(help_txt + "  " if help_txt else "") + "comma separated")
        return st.multiselect(label, opts, default=_default_list(inp, opts), key=key, help=help_txt)
    if inp.kind == "slider":
        lo = float(inp.min if inp.min is not None else 0.0)
        hi = float(inp.max if inp.max is not None else max(lo + 1.0, lo))
        val = float(inp.default if inp.default is not None else lo)
        val = min(max(val, lo), hi)
        step = float(inp.step or ((hi - lo) / 100.0 or 1.0))
        return st.slider(label, min_value=lo, max_value=hi, value=val, step=step, key=key, help=help_txt)
    if inp.kind == "int":
        val = int(inp.default if inp.default is not None else 0)
        lo = int(inp.min) if inp.min is not None else None
        hi = int(inp.max) if inp.max is not None else None
        if lo is not None:
            val = max(val, lo)
        if hi is not None:
            val = min(val, hi)
        return st.number_input(label, value=val, min_value=lo, max_value=hi, step=1, key=key,
                              help=help_txt, format="%d")
    # number
    val = float(inp.default if inp.default is not None else 0.0)
    lo = float(inp.min) if inp.min is not None else None
    hi = float(inp.max) if inp.max is not None else None
    if lo is not None:
        val = max(val, lo)
    if hi is not None:
        val = min(val, hi)
    return st.number_input(label, value=float(val), min_value=lo, max_value=hi,
                          step=float(inp.step or 0.1), key=key, help=help_txt, format="%.6g")


def _default_list(inp: Any, opts: list) -> list:
    default = inp.default if isinstance(inp.default, (list, tuple)) else []
    keep = [d for d in default if d in opts]
    return keep or ([opts[0]] if opts and inp.kind == "multi" and default else keep)


def _file_widget(label: str, inp: Any, key: str, help_txt: str) -> Any:
    """Bundled example, uploaded file or pasted text -- Galaxy's dropdown of histories."""
    default = str(inp.default or "")
    pool = EXAMPLES
    pasted_default = bool(default) and default not in pool and len(default) < 4000
    idx = pool.index(default) if default in pool else 0
    modes = ["Bundled example", "Upload a file", "Paste text"]
    mode = st.radio(label, modes, index=2 if pasted_default and not default.startswith("examples/") else 0,
                    horizontal=True, key=key + "_mode", label_visibility="visible")
    if mode == "Upload a file":
        up = st.file_uploader(f"{label} - pick a file", type=EXT_BY_FMT.get(inp.fmt) or None,
                             key=key + "_up")
        if up is not None:
            return up
        st.caption("No file chosen yet, so the bundled example is used.")
        return pool[idx]
    if mode == "Paste text":
        body = st.text_area(f"{label} (pasted)", value="" if default.startswith("examples/") else default,
                           height=170, key=key + "_paste",
                           help="Pasted content is used exactly as if it were the input file")
        return body if body.strip() else pool[idx]
    return st.selectbox(label, pool, index=idx, key=key + "_sel",
                       help=help_txt or "Files bundled in data/examples")


def _run(t: Any, values: dict) -> None:
    ss = st.session_state
    try:
        kw = registry.coerce(t, {k: v for k, v in values.items() if v is not None})
        res = t.fn(**kw)
        ss["results"][t.id] = _normalise(res)
        ss["errors"].pop(t.id, None)
    except Exception as exc:  # noqa: BLE001 - surface the failure like a Galaxy job error
        ss["errors"][t.id] = traceback.format_exc()
        st.session_state.setdefault("results", {}).pop(t.id, None)
        st.error(f"{type(exc).__name__}: {exc}")
        with st.expander("Traceback"):
            st.code(str(ss["errors"][t.id]))


def _normalise(res: Any) -> dict:
    """Turn a tool return value into renderable parts."""
    if res is None:
        return {"text": "", "message": "tool returned nothing"}
    if isinstance(res, str):
        return {"text": res, "message": ""}
    if isinstance(res, pd.DataFrame):
        return {"table": res, "message": "", "filename": "result.tsv"}
    if isinstance(res, dict):
        out = dict(res)
        stats = out.get("stats")
        figures: list[bytes] = []
        if isinstance(stats, dict):
            clean: dict[str, Any] = {}
            for k, v in stats.items():
                img = _fig_to_bytes(v) if (hasattr(v, "savefig") or isinstance(v, (bytes, bytearray))) else None
                if k in ("figure", "plot", "image") or img:
                    if img:
                        figures.append(img)
                    continue
                clean[k] = v
            out["stats"] = clean
            out["_figures"] = figures
        else:
            out["_figures"] = figures
        return out
    # plots returned directly (tools that return a matplotlib figure)
    img = _fig_to_bytes(res)
    if img:
        return {"image": img, "message": "", "filename": "figure.png", "_figures": [img]}
    return {"text": str(res), "message": ""}


def _render_output(t: Any) -> None:
    ss = st.session_state
    st.subheader("Output")
    err = ss["errors"].get(t.id)
    res = ss["results"].get(t.id)
    if err and not res:
        st.error("Job exited with an error")
        with st.expander("Traceback", expanded=True):
            st.code(err)
        return
    if not res:
        st.info("Fill the inputs and press **Run tool**. The bundled example values are "
                "already valid, so one click is enough.")
        if t.example:
            st.caption("Tool defaults: " + ", ".join(f"`{k}={_short(v, 30)}`"
                                                  for k, v in t.example.items()))
        return
    if res.get("message"):
        st.success(str(res["message"]))
    for img in res.get("_figures", []) or []:
        st.image(img, width="stretch")
    if res.get("image") and not res.get("_figures"):
        st.image(res["image"], width="stretch")
        st.download_button("Download figure (PNG)", res["image"],
                          file_name=res.get("filename", "figure.png"), mime="image/png")
    if res.get("table") is not None:
        df = res["table"]
        if isinstance(df, pd.DataFrame) and not df.empty:
            st.dataframe(df, width="stretch", hide_index=True,
                         height=min(430, 44 + 31 * len(df)))
            st.download_button(f"Download table ({len(df)} rows)", df.to_csv(index=False).encode(),
                              file_name=res.get("filename", "result.tsv"), mime="text/tab-separated-values")
        elif isinstance(df, pd.DataFrame):
            st.info("The tool returned an empty table.")
    if res.get("text"):
        body = str(res["text"])
        st.text_input("Preview", value=_short(body, 200), disabled=True)
        st.code(body if len(body) < 20000 else body[:20000] + "\n… truncated",
               language="text" if not body.lstrip().startswith(("{", "[")) else "json")
        st.download_button("Download text output", body.encode(),
                          file_name=res.get("filename", "result.txt"), mime="text/plain")
    stats = res.get("stats")
    if isinstance(stats, dict) and stats:
        with st.expander(f"Statistics and extra datasets ({len(stats)})", expanded=False):
            for k, v in stats.items():
                st.markdown(f"**{k}**")
                _render_stat(v)
    with st.expander("Raw result (JSON)"):
        st.code(json.dumps(_jsonable({k: v for k, v in res.items() if not k.startswith("_")})),
               language="json")


def _render_stat(v: Any) -> None:
    if isinstance(v, pd.DataFrame):
        st.dataframe(v, width="stretch", hide_index=True)
    elif isinstance(v, list) and v and all(isinstance(x, dict) for x in v):
        st.dataframe(pd.DataFrame(v), width="stretch", hide_index=True,
                    height=min(360, 44 + 31 * len(v)))
    elif isinstance(v, dict):
        st.json(v, expanded=2)
    elif isinstance(v, (int, float, str, bool)) or v is None:
        st.write(v if not isinstance(v, str) or len(v) < 700 else v[:700] + "…")
    elif isinstance(v, (list, tuple)):
        st.write(list(v)[:80])
    else:
        st.write(str(v)[:500])


def _jsonable(obj: Any) -> Any:
    if isinstance(obj, pd.DataFrame):
        return obj.head(50).to_dict("records")
    if isinstance(obj, dict):
        return {k: _jsonable(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [_jsonable(x) for x in list(obj)[:50]]
    if isinstance(obj, (bytes, bytearray)):
        return f"<{len(obj)} bytes>"
    return obj


def _cli_command(t: Any, values: dict) -> str:
    parts = [f"python -m chroma_titan run {t.id}"]
    for k, v in (values or {}).items():
        if isinstance(v, (str, int, float, bool)) or v is None:
            parts.append(f"--{k}={json.dumps(v)}" if v is not None else f"--{k}=true")
        elif isinstance(v, (list, tuple)):
            parts.append("--" + k + "=" + json.dumps(list(v)))
    return " \\\n  ".join(parts)


# ---------------------------------------------------------------------------
# other pages
# ---------------------------------------------------------------------------
def all_tools_page() -> None:
    st.title("All tools")
    rows = [t.to_dict(with_example=True) for t in registry.all_tools()]
    for r in rows:
        r["inputs"] = len(r.get("inputs", []))
        r["tags"] = ", ".join(r.get("tags") or [])
        r["example"] = ", ".join(f"{k}={_short(v, 18)}" for k, v in (r.get("example") or {}).items())
    df = pd.DataFrame(rows)[["id", "name", "group", "section_label", "output", "inputs", "tags",
                            "upstream", "version", "summary", "example"]]
    g = st.selectbox("Filter by panel group", ["(all)"] + sorted(df["group"].unique()))
    view = df if g == "(all)" else df[df["group"] == g]
    q = st.text_input("Filter by keyword", "")
    if q:
        low = q.lower()
        mask = view.astype(str).apply(lambda c: c.str.lower().str.contains(low, na=False)).any(axis=1)
        view = view[mask]
    st.caption(f"{len(view)} of {len(df)} tools")
    st.dataframe(view.sort_values(["group", "section_label", "id"]), width="stretch",
                hide_index=True, height=600)
    st.download_button("Download the full tool table (CSV)", df.to_csv(index=False).encode(),
                      file_name="chroma_titan_tools.csv", mime="text/csv")
    open_id = st.selectbox("Open a tool in the panel", [""] + sorted(df["id"]))
    if open_id:
        st.session_state["tool_id"] = open_id
        st.session_state["jump"] = open_id
        st.session_state["page"] = "Tool panel"
        st.rerun()


def datasets_page() -> None:
    st.title("Bundled example datasets")
    st.markdown("These files live in ``data/examples`` and are the default input of every tool, "
                "so any tool can be run without uploading anything. They were generated by "
                "``scripts/make_examples.py``.")
    rows = datasets.table_rows()
    df = pd.DataFrame(rows)
    st.dataframe(df[["name", "format", "section", "lines", "bytes", "description"]],
                width="stretch", hide_index=True, height=460)
    name = st.selectbox("Preview a dataset", list(df["name"]))
    st.code(datasets.preview(name, rows=25), language="text")
    p = datasets.path(name)
    try:
        with open(p, "rb") as fh:
            st.download_button(f"Download {name}", fh.read(), file_name=name)
    except OSError as exc:  # pragma: no cover
        st.warning(f"could not read {name}: {exc}")
    used = datasets.tools_using(name, limit=25)
    st.caption(f"{len(used)} tools use this file, e.g. " + ", ".join(f"`{u}`" for u in used[:10]))


def about_page() -> None:
    st.title("About this panel")
    st.markdown(f"""
Chroma Titan implements **{registry.count()} tools** in the taxonomy of
[usegalaxy.org](https://usegalaxy.org/): {len(panel.SECTION_LABELS)} panel sections inside
{len(panel.ALL_TOOLS_GROUPS)} groups, plus {len(panel.PANEL_VIEWS) - 1} specialised panel views
(Microbiology, Single Cell, Imaging, Biodiversity).

* each tool is a small Python function in ``chroma_titan/tools/*.py``, registered with the
  ``@tool`` decorator, declaring its inputs, output kind, tags, upstream wrapper and an example
* ``chroma_titan/core/*.py`` holds the engines (sequence, table, interval, alignment, statistics,
  plotting, imaging, phylogenetics), so tools run in-process without external binaries
* the panel layout in ``chroma_titan/panel.py`` follows ``panel_views/all_tools.yml`` from the
  ``galaxyproject/usegalaxy-playbook`` deployment repository

### Using it from the command line

```bash
python -m chroma_titan count                       # number of registered tools
python -m chroma_titan search "coverage"           # tool search
python -m chroma_titan info bedtools_coverage      # one tool's declaration
python -m chroma_titan run fasta_gc_content --src=examples/genome.fa
python -m chroma_titan sections                    # per-section tool counts
python -m chroma_titan check                       # execute every tool's example
python -m chroma_titan docs                        # write docs/TOOLS.md
```

### Adding a tool

```python
from chroma_titan.tools._common import T, fa, io, genome, text

@T("my_tool", "My new tool", "fasta_fastq", "text",
   [fa("src", "examples/genome.fa", "FASTA file")],
   ex={{"src": "examples/genome.fa"}}, up="some galaxy wrapper",
   tags=("demo",), summary="Short panel description.")
def my_tool(src):
    "One line docstring used as the fallback summary."
    recs = io.parse_fasta(io.as_text(src))
    return text(f"{{len(recs)}} sequences", f"read {{len(recs)}} records")
```

Then ``python -m chroma_titan check my_tool`` (or the full ``check``) proves it runs, and it
appears in the panel above automatically.
""")
    st.divider()
    st.markdown("Licensed under the Apache License 2.0 - see ``LICENSE``.")


# ---------------------------------------------------------------------------
def main() -> None:
    _init_state()
    registry.load()
    sidebar()
    page = st.session_state["page"]
    if page == "All tools":
        all_tools_page()
    elif page == "Datasets":
        datasets_page()
    elif page == "About this panel":
        about_page()
    else:
        tool_page()


try:
    main()
except SystemExit:  # pragma: no cover
    raise
except Exception as exc:  # noqa: BLE001 - Streamlit should always render something
    st.error(f"The app hit an error: {type(exc).__name__}: {exc}")
    with st.expander("Traceback"):
        st.code(traceback.format_exc())
