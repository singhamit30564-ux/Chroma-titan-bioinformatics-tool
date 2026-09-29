"""Tool registry: declarations, discovery, search and dispatch.

A *tool* here is the same idea as a Galaxy tool: a short description of the
inputs it takes, the output it produces and the work it performs.  Instead of
an XML wrapper plus a command line, the work is a Python function, so tools can
be executed by the web UI, the CLI and the test-suite without a job queue.
"""

from __future__ import annotations

import dataclasses
import importlib
import inspect
import pkgutil
from typing import Any, Callable, Iterable

# Input kinds understood by the registry and by the Streamlit renderer.
KINDS = (
    "text",  # single line of text
    "code",  # multi line text / raw sequence
    "file",  # path (bundled example or uploaded file)
    "multifile",  # several paths
    "int",
    "number",
    "slider",
    "choice",
    "multi",  # multi select
    "bool",
)

#: output kind a tool declares (drives how the panel/app renders the result)
OUTPUTS = (
    "text",
    "table",
    "stats",
    "figure",
    "image",
    "text_table",
    "sequences",
    "fasta",
    "bed",
    "gff",
    "vcf",
    "bedgraph",
    "file",
    "json",
)


@dataclasses.dataclass(frozen=True)
class Input:
    """One tool parameter, mirroring a ``<param>`` element of a Galaxy tool."""

    name: str
    kind: str = "text"
    label: str | None = None
    help: str = ""
    default: Any = None
    options: tuple = ()
    min: float | None = None
    max: float | None = None
    step: float | None = None
    fmt: str = ""  # accepted file format/suffix for kind == "file"

    def __post_init__(self) -> None:
        if self.kind not in KINDS:
            raise ValueError(f"unknown input kind {self.kind!r} for {self.name!r}")


@dataclasses.dataclass
class Tool:
    id: str
    name: str
    section: str
    summary: str
    inputs: tuple[Input, ...]
    output: str = "text"
    example: dict = dataclasses.field(default_factory=dict)
    upstream: str = ""
    version: str = "1.0.0"
    tags: tuple[str, ...] = ()
    fn: Callable[..., Any] | None = None

    # -- metadata helpers -------------------------------------------------
    @property
    def group(self) -> str:
        from chroma_titan.panel import group_of

        return group_of(self.section)

    @property
    def section_label(self) -> str:
        from chroma_titan.panel import section_label

        return section_label(self.section)

    @property
    def categories(self) -> tuple[str, str]:
        return (self.group, self.section_label)

    def to_dict(self, with_example: bool = True) -> dict:
        d = {
            "id": self.id,
            "name": self.name,
            "group": self.group,
            "section": self.section,
            "section_label": self.section_label,
            "summary": self.summary,
            "output": self.output,
            "version": self.version,
            "upstream": self.upstream,
            "tags": list(self.tags),
            "inputs": [dataclasses.asdict(i) for i in self.inputs],
        }
        if with_example:
            d["example"] = {k: _short(v) for k, v in self.example.items()}
        return d

    def __call__(self, **kwargs: Any) -> Any:
        assert self.fn is not None
        return self.fn(**kwargs)


def _short(value: Any, width: int = 90) -> Any:
    if isinstance(value, str) and len(value) > width:
        return value[: width - 1] + "…"
    return value


REGISTRY: dict[str, Tool] = {}
_LOADED = False


def tool(
    id: str,
    name: str,
    section: str,
    *,
    output: str = "text",
    inputs: Iterable[Input] = (),
    example: dict | None = None,
    upstream: str = "",
    version: str = "1.0.0",
    tags: tuple[str, ...] = (),
    summary: str = "",
):
    """Register ``fn`` as a tool.  ``summary`` defaults to the docstring."""

    inputs = tuple(inputs)
    for inp in inputs:
        if inp.name == "fn":
            raise ValueError(f"{id}: 'fn' is a reserved input name")

    def deco(fn: Callable[..., Any]):
        doc = (fn.__doc__ or "").strip()
        first = doc.split("\n\n")[0].strip().replace("\n", " ") if doc else ""
        summary_txt = (summary or first or name).strip()
        sig = inspect.signature(fn)
        declared = {i.name for i in inputs}
        extra: list[Input] = []
        for pname, param in sig.parameters.items():
            if pname in ("self", "kwargs"):
                continue
            if pname not in declared:
                if param.default is inspect.Parameter.empty:
                    raise ValueError(
                        f"tool {id}: parameter {pname!r} has no Input declaration"
                    )
                extra.append(Input(pname, **infer_spec(param.default)))
                declared.add(pname)
        all_inputs = tuple(inputs) + tuple(extra)
        t = Tool(
            id=id,
            name=name,
            section=section,
            summary=summary_txt,
            inputs=all_inputs,
            output=output,
            example=dict(example or {}),
            upstream=upstream,
            version=version,
            tags=tuple(tags),
            fn=fn,
        )
        if id in REGISTRY:
            raise ValueError(f"duplicate tool id: {id}")
        REGISTRY[id] = t
        fn.__tool__ = t  # type: ignore[attr-defined]
        return fn

    return deco


# -- input convenience constructors ------------------------------------
def In(name: str, **kw: Any) -> Input:
    return Input(name, **kw)


def code(name: str, default: str = "", **kw: Any) -> Input:
    return Input(name, kind="code", default=default, **kw)


def f(name: str, fmt: str = "txt", default: str = "", **kw: Any) -> Input:
    return Input(name, kind="file", fmt=fmt, default=default, **kw)


def num(name: str, default: float = 0.0, **kw: Any) -> Input:
    kw.setdefault("step", 0.1)
    return Input(name, kind="number", default=default, **kw)


def integer(name: str, default: int = 0, **kw: Any) -> Input:
    kw.setdefault("step", 1)
    return Input(name, kind="int", default=default, **kw)


def flag(name: str, default: bool = False, **kw: Any) -> Input:
    return Input(name, kind="bool", default=default, **kw)


def pick(name: str, options: list | tuple, default: Any = None, **kw: Any) -> Input:
    if not options:
        raise ValueError("choice inputs need options")
    return Input(
        name,
        kind="choice",
        options=tuple(options),
        default=default if default is not None else options[0],
        **kw,
    )


def some(name: str, options: list | tuple, default: list | None = None, **kw: Any) -> Input:
    return Input(
        name, kind="multi", options=tuple(options), default=list(default or options[:1]), **kw
    )


def slide(name: str, lo: float, hi: float, default: float | None = None, **kw: Any) -> Input:
    return Input(name, kind="slider", min=lo, max=hi, default=lo if default is None else default, **kw)


def infer_spec(default: Any) -> dict:
    """Guess a widget spec for an undeclared parameter from its default value."""
    if isinstance(default, bool):
        return {"kind": "bool"}
    if isinstance(default, int):
        return {"kind": "int"}
    if isinstance(default, float):
        return {"kind": "number"}
    if isinstance(default, (list, tuple)):
        return {"kind": "multi", "options": tuple(str(x) for x in default),
                "default": list(default)}
    if isinstance(default, str):
        low = default.lower()
        if low.startswith("examples/") or low.endswith((".fa", ".fasta", ".fastq", ".bed",
                                                        ".gff", ".vcf", ".sam", ".tsv", ".csv",
                                                        ".txt", ".pgm", ".nwk", ".gmt", ".faa")):
            return {"kind": "file", "fmt": low.rsplit(".", 1)[-1]}
    return {"kind": "text"}


# -- discovery -----------------------------------------------------------
def load() -> None:
    """Import every module in ``chroma_titan.tools`` (idempotent)."""
    global _LOADED
    if _LOADED:
        return
    import chroma_titan.tools as pkg

    for mod in pkgutil.iter_modules(pkg.__path__):
        importlib.import_module(f"chroma_titan.tools.{mod.name}")
    _LOADED = True


def all_tools() -> list[Tool]:
    load()
    return list(REGISTRY.values())


def get_tool(tool_id: str) -> Tool:
    load()
    try:
        return REGISTRY[tool_id]
    except KeyError:  # pragma: no cover - defensive
        raise KeyError(f"no tool {tool_id!r}; try one of {', '.join(sorted(REGISTRY))[:200]}") from None


def tool_ids() -> list[str]:
    load()
    return sorted(REGISTRY)


def count() -> int:
    load()
    return len(REGISTRY)


def _haystack(t: Tool) -> str:
    return " ".join(
        [t.id, t.name, t.summary, t.section, t.section_label, t.group, " ".join(t.tags)]
    ).lower()


def search(query: str = "", section: str = "", limit: int = 0) -> list[Tool]:
    """Case-insensitive substring search over id/name/summary/category."""
    load()
    q = (query or "").strip().lower()
    out: list[Tool] = []
    for t in REGISTRY.values():
        if section and t.section != section:
            continue
        if q:
            words = q.split()
            hay = _haystack(t)
            if not all(w in hay for w in words):
                continue
        out.append(t)
    out.sort(key=lambda t: (t.section, t.id))
    return out[:limit] if limit else out


# -- value coercion ------------------------------------------------------
def coerce(tool_obj: Tool, values: dict) -> dict:
    """Turn raw UI/CLI values into what the tool function expects."""
    from chroma_titan.core import io as _io

    kw: dict[str, Any] = {}
    declared = {i.name: i for i in tool_obj.inputs}
    for key, raw in (values or {}).items():
        spec = declared.get(key)
        if spec is None:
            kw[key] = raw
            continue
        kw[key] = _coerce_one(spec, raw, _io)
    for name, spec in declared.items():
        if name not in kw and spec.default is not None:
            kw[name] = spec.default
    return kw


def _coerce_one(spec: Input, raw: Any, _io) -> Any:
    if spec.kind in ("file",):
        if raw is None or raw == "":
            return None
        if not _io.looks_like_path(raw):
            return raw  # pasted text payload
        return _io.resolve(raw)
    if spec.kind == "multifile":
        if raw is None or raw == "":
            return []
        if isinstance(raw, (list, tuple)):
            return [_io.resolve(r) if _io.looks_like_path(r) else r for r in raw]
        return [_io.resolve(p) if _io.looks_like_path(p) else p
                for p in str(raw).replace(",", " ").split()]
    if raw is None or raw == "":
        return None
    if spec.kind in ("int", "slider") and spec.step == int(spec.step or 1):
        return int(float(raw))
    if spec.kind in ("number", "slider"):
        return float(raw)
    if spec.kind == "bool":
        if isinstance(raw, str):
            return raw.strip().lower() in ("1", "true", "yes", "y", "on")
        return bool(raw)
    if spec.kind == "multi":
        if isinstance(raw, str):
            return [x for x in raw.replace(",", " ").split() if x]
        return list(raw)
    return raw


def run(tool_id: str, **values: Any) -> Any:
    """Execute a tool by id with loosely typed values."""
    t = get_tool(tool_id)
    assert t.fn is not None
    return t.fn(**coerce(t, values))


# Backwards friendly alias used by app/CLI.
run_tool = run
