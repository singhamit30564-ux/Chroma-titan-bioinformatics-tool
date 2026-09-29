"""Command line interface: ``python -m chroma_titan <command>``.

Commands
--------
list, search, sections, info, run, examples, count, check, docs, app
"""

from __future__ import annotations

import argparse
import ast
import json
import sys
import traceback
from pathlib import Path

from chroma_titan import panel, registry


def _print_tools(tools, width: int = 118) -> None:
    cur = None
    for t in sorted(tools, key=lambda x: (x.group, x.section, x.id)):
        key = (t.group, t.section_label)
        if key != cur:
            cur = key
            print(f"\n[{key[0]}] {key[1]}")
        print(f"  {t.id:<34} {t.summary[:width - 46]}")


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(prog="chroma_titan",
                                 description="Chroma Titan -- Galaxy-style bioinformatics tool suite")
    sub = ap.add_subparsers(dest="cmd")
    p = sub.add_parser("list", help="list registered tools")
    p.add_argument("--section", default="")
    p.add_argument("--group", default="")
    p = sub.add_parser("search", help="substring search over tools")
    p.add_argument("query")
    p.add_argument("--section", default="")
    sub.add_parser("sections", help="panel sections with tool counts")
    sub.add_parser("count", help="number of registered tools")
    p = sub.add_parser("info", help="show one tool's declaration")
    p.add_argument("tool_id")
    p.add_argument("--json", action="store_true")
    p = sub.add_parser("run", help="run a tool: run <id> [--key value ...]")
    p.add_argument("tool_id")
    p.add_argument("rest", nargs="*")
    p.add_argument("--json", action="store_true")
    sub.add_parser("examples", help="bundled example datasets")
    p = sub.add_parser("check", help="run every tool's example (self test)")
    p.add_argument("--section", default="")
    p.add_argument("--verbose", "-v", action="store_true")
    p.add_argument("--failfast", action="store_true")
    p = sub.add_parser("docs", help="generate the tool reference (docs/TOOLS.md)")
    p.add_argument("--out", default="docs/TOOLS.md")
    p = sub.add_parser("app", help="print the streamlit command")
    args, extra = ap.parse_known_args(argv)
    if args.cmd == "run" and extra:
        # tool parameters are not argparse options: ``--key value`` / ``--key=value``
        args.rest = list(args.rest) + list(extra)
    elif extra:
        ap.error(f"unrecognized arguments: {' '.join(extra)}")
    registry.load()

    if args.cmd == "count":
        print(registry.count())
        return 0
    if args.cmd == "list":
        _print_tools(registry.search("", args.section))
        print(f"\n{len(registry.search('', args.section))} tools")
        return 0
    if args.cmd == "search":
        _print_tools(registry.search(args.query, args.section))
        return 0
    if args.cmd == "sections":
        rows: dict[str, int] = {}
        for t in registry.all_tools():
            rows[t.section] = rows.get(t.section, 0) + 1
        for sec in panel.ordered_sections():
            if sec in rows:
                print(f"{rows[sec]:4d}  {sec:<44} {panel.section_label(sec)}")
        print(f"{'':4}  {'-' * 60}\n{sum(rows.values()):4d}  total tools in {len(rows)} sections")
        return 0
    if args.cmd == "info":
        t = registry.get_tool(args.tool_id)
        d = t.to_dict()
        d["example"] = t.example
        if args.json:
            print(json.dumps(d, indent=2, default=str))
            return 0
        print(f"{t.name}  ({t.id})")
        print(f"group   : {t.group}")
        print(f"section : {t.section}  [{t.section_label}]")
        print(f"upstream: {t.upstream or '-'}   version: {t.version}   output: {t.output}")
        print(f"summary : {t.summary}")
        print("inputs  :")
        for i in t.inputs:
            print(f"  - {i.name:<18} {i.kind:<10} default={i.default!r:<40} "
                  f"{'options=' + str(i.options) if i.options else ''}")
        if t.example:
            print("example :")
            print("  " + " ".join(f"--{k} {_short(v)}" for k, v in t.example.items()))
        return 0
    if args.cmd == "run":
        kw = {}
        toks = list(args.rest)
        i = 0
        while i < len(toks):
            tok = str(toks[i])
            if tok.startswith("--"):
                k, _, v = tok[2:].partition("=")
                if not v and i + 1 < len(toks) and not str(toks[i + 1]).startswith("--"):
                    v = str(toks[i + 1])
                    i += 1
                kw[k] = _parse(v or "true")
            elif "=" in tok:
                k, _, v = tok.partition("=")
                kw[k] = _parse(v)
            i += 1
        try:
            res = registry.run(args.tool_id, **kw)
        except Exception as exc:  # noqa: BLE001
            print(f"ERROR {type(exc).__name__}: {exc}", file=sys.stderr)
            traceback.print_exc()
            return 2
        print(json.dumps(_jsonable(res), indent=2, default=str) if args.json else _render(res))
        return 0
    if args.cmd == "examples":
        from chroma_titan.core.io import EXAMPLE_DIR

        for p in sorted(EXAMPLE_DIR.glob("*")):
            print(f"{p.name:<26} {p.stat().st_size:>9,} bytes")
        return 0
    if args.cmd == "check":
        return _check(args.section, args.verbose, args.failfast)
    if args.cmd == "docs":
        return _docs(Path(args.out))
    if args.cmd == "app":
        print("streamlit run app.py")
        return 0
    ap.print_help()
    return 1


def _short(v) -> str:
    s = str(v).replace("\n", " ")
    return (s[:28] + "…") if len(s) > 30 else s


def _parse(v: str):
    try:
        return ast.literal_eval(v)
    except (ValueError, SyntaxError):
        return v


def _jsonable(obj):
    if isinstance(obj, dict):
        return {k: _jsonable(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [_jsonable(v) for v in obj][:200]
    if hasattr(obj, "to_dict"):
        try:
            return _jsonable(json.loads(obj.to_json()))
        except Exception:  # noqa: BLE001
            return str(obj)[:500]
    if isinstance(obj, bytes):
        return f"<{len(obj)} bytes>"
    return obj


def _render(res) -> str:
    if isinstance(res, dict):
        parts = []
        if "message" in res:
            parts.append(str(res["message"]))
        if "stats" in res and isinstance(res["stats"], dict):
            parts.append("\n".join(f"  {k}: {v}" for k, v in res["stats"].items()))
        if "table" in res:
            df = res["table"]
            parts.append(df.head(40).to_string(index=False) + f"\n  ... {len(df)} rows")
        if "text" in res and "table" not in res:
            body = str(res["text"])
            parts.append(body if len(body) < 4000 else body[:4000] + "\n... [truncated]")
        if "image" in res:
            parts.append(f"  <image {len(res['image'])} bytes>")
        return "\n".join(parts) or json.dumps(_jsonable(res), indent=2, default=str)
    return str(res)


def _check(section: str = "", verbose: bool = False, failfast: bool = False) -> int:
    tools = registry.search("", section)
    ok = fail = skipped = 0
    errors: list[tuple[str, str]] = []
    for t in tools:
        if not t.example:
            skipped += 1
            if verbose:
                print(f"SKIP {t.id}: no example")
            continue
        try:
            res = registry.run(t.id, **t.example)
            bad = res is None or (isinstance(res, dict) and not res) or \
                (isinstance(res, str) and not res.strip())
            if bad:
                fail += 1
                errors.append((t.id, "returned empty result"))
            else:
                ok += 1
                if verbose:
                    print(f"OK   {t.id}")
        except Exception as exc:  # noqa: BLE001
            fail += 1
            msg = f"{type(exc).__name__}: {exc}"
            errors.append((t.id, msg))
            if failfast:
                print(f"FAIL {t.id}: {msg}")
                traceback.print_exc()
                return 1
        finally:
            try:  # keep matplotlib from accumulating open figures over a long run
                import matplotlib.pyplot as _plt

                _plt.close("all")
            except Exception:  # noqa: BLE001
                pass
    print(f"\n{ok} ok, {fail} failed, {skipped} without example "
          f"({len(tools)} tools in section {section or 'ALL'})")
    if errors:
        print("\nfailures:")
        by: dict[str, list[str]] = {}
        for tid, msg in errors:
            kind = msg.split(":")[0]
            by.setdefault(kind, []).append(tid)
        for kind, ids in sorted(by.items(), key=lambda kv: -len(kv[1])):
            print(f"  {kind:<28} {len(ids):4d}  e.g. {', '.join(ids[:6])}")
        if verbose:
            for tid, msg in errors:
                print(f"  ! {tid}: {msg[:160]}")
    return 1 if errors else 0


def _docs(out: Path) -> int:
    import collections

    tools = registry.all_tools()
    by_sec: dict[str, list] = collections.defaultdict(list)
    for t in tools:
        by_sec[t.section].append(t)
    lines = ["# Chroma Titan tool reference",
             "",
             f"{len(tools)} tools, grouped like the tool panel of [usegalaxy.org](https://usegalaxy.org/).",
             "This file is generated: `python -m chroma_titan docs`.", "", "| Section | Group | Tools |",
             "|---|---|---|"]
    for sec in panel.ordered_sections():
        if sec in by_sec:
            lines.append(f"| [`{sec}`](#{sec}) | {panel.group_of(sec)} | {len(by_sec[sec])} |")
    lines.append(f"| **total** |  | **{len(tools)}** |")
    lines.append("")
    for sec in panel.ordered_sections():
        if sec not in by_sec:
            continue
        lines.append(f"## {sec}")
        lines.append("")
        lines.append(f"*{panel.section_label(sec)}* — in Galaxy group *{panel.group_of(sec)}*.")
        lines.append("")
        lines.append("| Tool | What it does | Upstream Galaxy tool |")
        lines.append("|---|---|---|")
        for t in sorted(by_sec[sec], key=lambda x: x.id):
            up = f"`{t.upstream}`" if t.upstream else "—"
            lines.append(f"| [`{t.id}`](#{t.id}) | {t.summary[:150]} | {up} |")
        lines.append("")
        for t in sorted(by_sec[sec], key=lambda x: x.id):
            lines.append(f"### {t.id}")
            lines.append("")
            lines.append(f"**{t.name}** — {t.summary}")
            lines.append("")
            lines.append("| parameter | kind | default | options |")
            lines.append("|---|---|---|---|")
            for i in t.inputs:
                lines.append(f"| `{i.name}` | {i.kind} | `{i.default}` | "
                             f"{', '.join(map(str, i.options)) or '—'} |")
            if t.upstream:
                lines.append("")
                lines.append(f"Galaxy Tool Shed: `{t.upstream}`")
            lines.append("")
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text("\n".join(lines) + "\n")
    print(f"wrote {out} ({len(lines)} lines, {len(tools)} tools)")
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
