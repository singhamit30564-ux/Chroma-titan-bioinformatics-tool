"""General text tools: Text Manipulation, Filter and Sort, Join/Subtract/Group,
Datamash aggregation, Expression tools, Collection operations, Get/Send data.

These mirror the small but heavily used Galaxy tools that operate on delimited
text tables (`cut`, `paste`, `join`, `head`, `tail`, `grep`, `datamash`, ...).
"""

from __future__ import annotations

import math
import random
import re
from collections import Counter, OrderedDict, defaultdict

from chroma_titan.core import io, stats, tables
from chroma_titan.tools._common import *  # noqa: F401,F403

# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------


def _lines(src) -> list[str]:
    return [ln for ln in io.as_text(src).splitlines() if ln.strip()]


def _fields(src, sep: str | None = None) -> list[list[str]]:
    out = []
    for ln in _lines(src):
        out.append(ln.split(sep) if sep else ln.split("\t") if "\t" in ln else ln.split())
    return out


def _unparse(rows: list[list[str]], sep: str = "\t") -> str:
    return "\n".join(sep.join(r) for r in rows) + "\n"


def _head_of(rows, header=True):
    """Header row of a parsed table, or None when it is numeric data."""
    if not rows or header is False:
        return None
    return rows[0] if _is_header(rows, "auto") else None


def _ci(spec, head, ncol: int) -> int:
    """Resolve a column reference given by name, 1-based index or 'last'."""
    s = str(spec).strip()
    if head:
        for i, h in enumerate(head):
            if str(h).strip() == s or str(h).strip().lower() == s.lower():
                return i
    if s.lower() in ("last", "l"):
        return ncol - 1
    if s.lower() in ("first", "f"):
        return 0
    idx = _col_index(s, ncol)
    return idx[0] if idx else 0


def _cis(spec, head, ncol: int) -> list[int]:
    """Resolve a multi-column reference (names or numbers)."""
    s = str(spec).strip()
    if s in ("", "all", "*", "none"):
        return list(range(ncol))
    parts = [p for p in re.split(r"[,;\s]+", s) if p]
    if head and parts and all(p in {str(h).strip() for h in head} for p in parts):
        return [i for i, h in enumerate(head) if str(h).strip() in parts]
    return _col_index(s, ncol) or list(range(ncol))


def _is_header(rows: list[list[str]], header: str) -> bool:
    if header == "yes":
        return True
    if header == "no":
        return False
    if not rows:
        return False
    return not all(_num(c) for c in rows[0] if c)


def _num(x) -> bool:
    try:
        float(x)
        return True
    except (TypeError, ValueError):
        return False


def _col_index(spec: str, ncols: int) -> list[int]:
    """Parse Galaxy-style column selection: 1,3,5-7,c3."""
    out: list[int] = []
    for part in str(spec).replace(",", " ").split():
        m = re.match(r"^c?(\d+)(?:-c?(\d+))?$", part.strip(), re.I)
        if not m:
            continue
        a = int(m.group(1)) - 1
        b = int(m.group(2)) - 1 if m.group(2) else a
        out.extend(range(a, b + 1))
    return [i for i in out if 0 <= i < ncols]


# ===========================================================================
# Text Manipulation
# ===========================================================================
@T("text_cut_columns", "Cut columns from a table", "text_manipulation",
   inputs=[textbox("columns", "1,3", help="e.g. 1,3,5-7 (1-based)"), tbl("src")],
   ex={"src": PHENO, "columns": "1,3"}, up="text_processing")
def cut_columns(src, columns="1,3"):
    """Keep the selected columns of a delimited table (POSIX ``cut``)."""
    rows = _fields(src)
    if not rows:
        return text("", "empty input")
    idx = _cis(columns, _head_of(rows), max(len(r) for r in rows))
    return text(_unparse([[r[i] for i in idx if i < len(r)] for r in rows]),
                f"kept {len(idx)} of {len(rows)} lines")


@T("text_paste_columns", "Paste tool (join files side by side)", "text_manipulation",
   inputs=[tbl("src"), tbl("other"), choice("delimiter", ["tab", "comma", "space", "newline"]),
           boolean("whole_line", True)],
   ex={"src": PHENO, "other": REGIONS_BED, "delimiter": "tab", "whole_line": True},
   up="text_processing")
def paste_columns(src, other, delimiter="tab", whole_line=True):
    """Concatenate two files line by line (``paste``)."""
    d = {"tab": "\t", "comma": ",", "space": " ", "newline": "\n"}[delimiter]
    a, b = _lines(src), _lines(other)
    n = max(len(a), len(b))
    out = []
    for i in range(n):
        left = a[i] if i < len(a) else ""
        right = b[i] if i < len(b) else ""
        if not whole_line:
            right = right.split("\t")[0]
        out.append(left + d + right)
    return text("\n".join(out) + "\n", f"pasted {len(a)} + {len(b)} lines")


@T("text_join_files", "Join two delimited files on a column", "text_manipulation",
   inputs=[tbl("left"), tbl("right"), intin("left_key", 0), intin("right_key", 0),
           choice("mode", ["inner", "outer", "left outer"]), boolean("header", True)],
   ex={"left": PHENO, "right": REGIONS_BED, "left_key": 0, "right_key": 3},
   up="join")
def join_files(left, right, left_key=0, right_key=0, mode="inner", header=True):
    """Equi-join on one key column (``join``/SQL inner-outer)."""
    lrows, rrows = _fields(left), _fields(right)
    if not lrows or not rrows:
        return text("", "empty input")
    lhdr = _is_header(lrows, "auto") and header
    rhdr = _is_header(rrows, "auto") and header
    ldata, rdata = lrows[1:] if lhdr else lrows, rrows[1:] if rhdr else rrows
    rindex: dict[str, list[list[str]]] = defaultdict(list)
    for r in rdata:
        if right_key < len(r):
            rindex[r[right_key].split("|")[0]].append(r)
    lkeys = lrows[0] if lhdr else [f"col{i + 1}" for i in range(max(len(x) for x in ldata))]
    rkeys = rrows[0] if rhdr else [f"col{i + 1}" for i in range(max(len(x) for x in rdata))]
    used_r: set[str] = set()
    out = []
    for row in ldata:
        if left_key >= len(row):
            continue
        k = row[left_key].split("|")[0]
        matches = rindex.get(k, [])
        if matches:
            used_r.add(k)
            for m in matches:
                out.append(row + [c for i, c in enumerate(m) if i != right_key])
        elif mode in ("outer", "left outer"):
            out.append(row + ["." for i in range(len(rkeys) - 1)])
    if mode == "outer":
        lset = {r[left_key].split("|")[0] for r in ldata if left_key < len(r)}
        for k, rows in rindex.items():
            if k not in lset:
                for r in rows:
                    out.append(["." for _ in lkeys] + r)
    hdr = [lkeys[i] for i in range(min(len(lkeys), left_key + 1))] if header else []
    head = (list(OrderedDict.fromkeys(lkeys + [c for i, c in enumerate(rkeys) if i != right_key]))
            if header else [])
    body = _unparse(out)
    return text(("\t".join(head) + "\n" if head else "") + body,
                f"{len(out)} joined rows ({mode})")


@T("text_join_with_loose_cutoff", "Relaxed join on a prefix", "text_manipulation",
   inputs=[tbl("left"), tbl("right"), intin("prefix_length", 8), intin("cutoff", 2),
           boolean("ignore_case", True)],
   ex={"left": PHENO, "right": REGIONS_BED, "prefix_length": 8, "cutoff": 1},
   up="loose_join")
def relaxed_join(left, right, prefix_length=8, cutoff=2, ignore_case=True):
    """Join two tables on their first column allowing up to ``cutoff`` mismatches."""
    lrows, rrows = _fields(left), _fields(right)
    if not lrows or not rrows:
        return text("", "empty input")
    rkeys = [r[0] for r in rrows if r]
    out = []
    for lr in lrows:
        if not lr:
            continue
        lk = lr[0].lower() if ignore_case else lr[0]
        best, bestd = None, cutoff + 1
        for rk in rkeys:
            key = rk.lower()[:prefix_length]
            a = lk[:prefix_length]
            d = sum(1 for x, y in zip(a, key) if x != y) + abs(len(a) - len(key))
            if d < bestd:
                best, bestd = rk, d
        if best is not None and bestd <= cutoff:
            rr = next(r for r in rrows if r[0] == best)
            out.append(lr + rr)
    return text(_unparse(out), f"{len(out)} relaxed-join rows")


@T("text_replace", "Find and replace (regular expressions)", "text_manipulation",
   inputs=[textbox("pattern"), textbox("replacement"), txt("src"), boolean("regex", True),
           boolean("case_insensitive", False)],
   ex={"src": PHENO, "pattern": "sample_", "replacement": "S", "regex": False},
   up="text_processing")
def replace_text(src, pattern="", replacement="", regex=True, case_insensitive=False):
    """``sed``-style substitution over every line of a text dataset."""
    flags = re.IGNORECASE if case_insensitive else 0
    rx = re.compile(pattern if regex else re.escape(pattern), flags)
    n = 0
    out = []
    for ln in io.as_text(src).splitlines():
        new, k = rx.subn(replacement, ln)
        n += k
        out.append(new)
    return text("\n".join(out) + "\n", f"{n} replacements")


@T("text_change_case", "Change case of a text dataset", "text_manipulation",
   inputs=[txt("src"), choice("mode", ["upper", "lower", "title", "sentence", "swap", "inverse"])],
   ex={"src": PHENO, "mode": "upper"}, up="change_case")
def change_case(src, mode="upper"):
    """Convert text to upper/lower/title case (Galaxy 'Change Case')."""
    lines = io.as_text(src).splitlines()
    fn = {"upper": str.upper, "lower": str.lower, "title": str.title, "swap": str.swapcase,
          "sentence": lambda s: s[:1].upper() + s[1:].lower(),
          "inverse": lambda s: "".join(c.upper() if c.islower() else c.lower() for c in s)}[mode]
    return text("\n".join(fn(l) for l in lines) + "\n", f"{len(lines)} lines changed to {mode}")


@T("text_concatenate", "Concatenate datasets end to end", "text_manipulation",
   inputs=[txt("first"), txt("second"), boolean("blank_line", False)],
   ex={"first": REGIONS_BED, "second": TARGETS_BED}, up="concat")
def concatenate(first, second, blank_line=False):
    """Stack two text datasets vertically."""
    a, b = io.as_text(first).rstrip("\n"), io.as_text(second).rstrip("\n")
    sep = "\n\n" if blank_line else "\n"
    body = a + sep + b + "\n"
    return text(body, f"{len(a.splitlines())} + {len(b.splitlines())} lines")


@T("text_head", "Select first N lines", "text_manipulation",
   inputs=[txt("src"), intin("n", 10)], ex={"src": PHENO, "n": 5}, up="head")
def head_tool(src, n=10):
    """``head -n``: the first n lines of a text file."""
    lines = _lines(src)
    return text("\n".join(lines[:int(n)]) + "\n", f"{min(n, len(lines))}/{len(lines)} lines")


@T("text_tail", "Select last N lines", "text_manipulation",
   inputs=[txt("src"), intin("n", 10), boolean("skip_header", False)],
   ex={"src": PHENO, "n": 5}, up="tail")
def tail_tool(src, n=10, skip_header=False):
    """``tail -n``: the last n lines, optionally keeping the header."""
    lines = _lines(src)
    head, body = ([lines[0]], lines[1:]) if skip_header else ([], lines)
    return text("\n".join(head + body[-int(n):]) + "\n", f"{min(n, len(lines))}/{len(lines)} lines")


@T("text_insert_lines", "Insert lines before/after every Nth line", "text_manipulation",
   inputs=[txt("src"), intin("interval", 2), textbox("text_line", "## separator"),
           choice("position", ["before", "after"])],
   ex={"src": PHENO, "interval": 3, "text_line": "## chunk", "position": "before"},
   up="insert_lines")
def insert_lines(src, interval=2, text_line="## separator", position="before"):
    """Insert a fixed line every ``interval`` lines."""
    lines = _lines(src)
    out = []
    for i, ln in enumerate(lines, 1):
        if i % interval == 1 and position == "before":
            out.append(text_line)
        out.append(ln)
        if i % interval == 0 and position == "after":
            out.append(text_line)
    return text("\n".join(out) + "\n", f"inserted {len(out) - len(lines)} lines")


def _text_wrap(src, width=80, indent=0, newlines=False):
    """Fold long lines to a fixed width."""
    out = []
    for ln in _lines(src):
        for start in range(0, max(1, len(ln)), width):
            out.append(" " * indent + ln[start:start + width])
        if newlines:
            out.append("")
    return text("\n".join(out) + "\n", f"wrapped to {width} columns")


@T("text_wrap", "Wrap text to fixed width", "text_manipulation",
   inputs=[txt("src"), intin("width", 60), intin("indent", 0), boolean("newlines", False)],
   ex={"src": PHENO, "width": 40}, up="textwrapper")
def wrap_tool(src, width=60, indent=0, newlines=False):
    """``fold``-style line wrapping."""
    return _text_wrap(src, width, indent, newlines)


@T("text_word_counter", "Word/character/line counter", "text_manipulation",
   inputs=[txt("src")], ex={"src": PHENO}, up="word_counter")
def word_counter(src):
    """``wc``: lines, words and characters, plus the top words."""
    body = io.as_text(src)
    words = re.findall(r"[A-Za-z0-9_\-\.]+", body)
    counts = Counter(w.lower() for w in words)
    rows = [{"metric": "lines", "value": len(body.splitlines())},
            {"metric": "words", "value": len(words)},
            {"metric": "characters", "value": len(body)},
            {"metric": "unique_words", "value": len(counts)},
            {"metric": "blank_lines", "value": sum(1 for l in body.splitlines() if not l.strip())}]
    res = table(rows, f"{len(body.splitlines())} lines / {len(words)} words")
    res["top_words"] = dict(counts.most_common(10))
    return res


@T("text_grep", "Filter lines matching a pattern", "text_manipulation",
   inputs=[txt("src"), textbox("pattern", "gene"), boolean("invert", False),
           boolean("ignore_case", True), boolean("whole_line", False), intin("max_results", 0)],
   ex={"src": PHENO, "pattern": "treated", "invert": False}, up="text_processing")
def grep_tool(src, pattern="", invert=False, ignore_case=True, whole_line=False, max_results=0):
    """``grep``: keep (or drop with -v) lines matching a regular expression."""
    rx = re.compile(pattern if not whole_line else re.escape(pattern),
                    re.IGNORECASE if ignore_case else 0)
    keep = []
    for ln in io.as_text(src).splitlines():
        hit = bool(rx.search(ln))
        if hit != invert:
            keep.append(ln)
    if max_results:
        keep = keep[:int(max_results)]
    return text("\n".join(keep) + "\n", f"{len(keep)} matching lines")


@T("text_tr", "Translate or delete characters", "text_manipulation",
   inputs=[txt("src"), textbox("set1", "ACGT"), textbox("set2", "TGCA"),
           choice("mode", ["translate", "delete", "squeeze", "complement"])],
   ex={"src": PHENO, "set1": "aeiou", "set2": "", "mode": "delete"}, up="text_processing")
def tr_tool(src, set1="", set2="", mode="translate"):
    """``tr``: translate set1→set2, delete characters, or squeeze repeats."""
    body = io.as_text(src)
    if mode == "delete":
        out = body.translate(str.maketrans("", "", set1))
    elif mode == "squeeze":
        out = re.sub("[" + re.escape(set1) + "]+", set1[:1] or "+", body)
    elif mode == "complement":
        out = body.translate(str.maketrans("ACGTacgt", "TGCAtgca"))
    else:
        tbl_ = {}
        for i, a in enumerate(set1):
            tbl_[ord(a)] = set2[i] if i < len(set2) else None
        out = body.translate(tbl_)
    return text(out, f"tr {mode} applied")


@T("text_sort", "Sort a text file line by line", "text_manipulation",
   inputs=[txt("src"), intin("column", -1), boolean("numeric", False), boolean("reverse", False),
           boolean("unique", False), boolean("header", False)],
   ex={"src": PHENO, "numeric": False, "column": 1}, up="sort")
def sort_lines(src, column=-1, numeric=False, reverse=False, unique=False, header=False):
    """``sort -k``: lexicographic or numeric sort of lines / a column."""
    lines = _lines(src)
    body = lines[1:] if header and lines else lines
    key = (lambda l: _field(l, column, numeric)) if column != -1 else \
        ((lambda l: float(l)) if numeric else (lambda l: l))
    try:
        body.sort(key=key, reverse=reverse)
    except ValueError:
        body.sort(key=str, reverse=reverse)
    if unique:
        body = list(OrderedDict.fromkeys(body))
    out = ([lines[0]] if header and lines else []) + body
    return text("\n".join(out) + "\n", f"sorted {len(body)} lines")


def _field(line: str, column: int, numeric: bool):
    parts = line.split("\t") if "\t" in line else line.split()
    i = int(column) - 1 if column > 0 else (len(parts) + int(column) if column < 0 else 0)
    val = parts[i] if 0 <= i < len(parts) else ""
    if numeric:
        try:
            return float(val)
        except ValueError:
            return float("nan")
    return val


@T("text_uniq", "Report unique/adjacent duplicate lines", "text_manipulation",
   inputs=[txt("src"), boolean("count", True), boolean("ignore_case", False),
           boolean("adjacent_only", False)],
   ex={"src": PHENO, "count": True}, up="uniq")
def uniq_tool(src, count=True, ignore_case=False, adjacent_only=False):
    """``uniq -c``: collapse repeated lines and count them."""
    lines = _lines(src)
    out: "OrderedDict[str,int]" = OrderedDict()
    prev = object()
    for ln in lines:
        k = ln.lower() if ignore_case else ln
        if adjacent_only:
            if k == prev:
                out[k] = out[k] + 1
            else:
                out[k] = out.get(k, 0) + 1
        else:
            out[k] = out.get(k, 0) + 1
        prev = k
    rows = [{"line": k, "count": v} for k, v in out.items()] if count else \
        [{"line": k} for k in out]
    return table(rows, f"{len(out)} unique lines from {len(lines)}")


@T("text_random_lines", "Shuffle or randomly sample lines", "text_manipulation",
   inputs=[txt("src"), intin("n", 0), number("fraction", 1.0), intin("seed", 1),
           boolean("header_kept", False)],
   ex={"src": PHENO, "n": 4, "seed": 7}, up="shuffle")
def random_lines(src, n=0, fraction=1.0, seed=1, header_kept=False):
    """``shuf``: shuffle lines, or draw a random subset."""
    lines = _lines(src)
    head = [lines[0]] if header_kept and lines else []
    body = lines[len(head):]
    rng = random.Random(seed)
    rng.shuffle(body)
    if n:
        body = body[:int(n)]
    elif fraction and fraction < 1:
        body = body[: max(1, int(len(body) * fraction))]
    return text("\n".join(head + body) + "\n", f"{len(body)} random lines")


@T("text_subtract", "Subtract one text file from another", "text_manipulation",
   inputs=[txt("main"), txt("sub"), boolean("whole_line", True)],
   ex={"main": REGIONS_BED, "sub": TARGETS_BED}, up="subtract")
def subtract_text(main, sub, whole_line=True):
    """Remove every line of ``main`` that also occurs in ``sub``."""
    drop = set()
    for ln in _lines(sub):
        drop.add(ln.strip() if whole_line else ln.split("\t")[0])
    keep = [ln for ln in _lines(main)
            if (ln.strip() if whole_line else ln.split("\t")[0]) not in drop]
    return text("\n".join(keep) + "\n", f"{len(keep)} lines remain ({len(_lines(sub))} subtracted)")


@T("text_compare_two_files", "Compare two files line by line", "text_manipulation",
   inputs=[txt("first"), txt("second")], ex={"first": REGIONS_BED, "second": TARGETS_BED},
   up="compare")
def compare_files(first, second):
    """``cmp``-style diff: identical lines, only-in-A and only-in-B counts."""
    a, b = [l.strip() for l in _lines(first)], [l.strip() for l in _lines(second)]
    sa, sb = set(a), set(b)
    rows = [{"metric": "lines_in_first", "value": len(a)},
            {"metric": "lines_in_second", "value": len(b)},
            {"metric": "identical", "value": sum(1 for x, y in zip(a, b) if x == y)},
            {"metric": "shared", "value": len(sa & sb)},
            {"metric": "only_in_first", "value": len(sa - sb)},
            {"metric": "only_in_second", "value": len(sb - sa)},
            {"metric": "jaccard", "value": round(len(sa & sb) / len(sa | sb), 4) if sa | sb else 0.0}]
    return table(rows, "file comparison")


@T("text_split_on_column", "Split a table into chunks by column value", "text_manipulation",
   inputs=[tbl("src"), intin("column", 1), intin("max_groups", 20)],
   ex={"src": PHENO, "column": 2}, up="split")
def split_on_column(src, column=1, max_groups=20):
    """Group lines by the value of one column (like ``split`` per category)."""
    groups: dict[str, list[str]] = defaultdict(list)
    for ln in _lines(src):
        f = ln.split("\t") if "\t" in ln else ln.split()
        key = f[column - 1] if 0 < column <= len(f) else "NA"
        groups[key].append(ln)
    rows = [{"group": k, "n_lines": len(v), "first_line": v[0][:60]}
            for k, v in list(groups.items())[:int(max_groups)]]
    return table(rows, f"{len(groups)} groups")


@T("text_column_formatter", "Reformat a table to aligned fixed width", "text_manipulation",
   inputs=[tbl("src"), intin("pad", 12), boolean("header", True)], ex={"src": PHENO, "pad": 14},
   up="table_formatter")
def column_formatter(src, pad=12, header=True):
    """Pretty-print a delimited table with aligned columns."""
    rows = _fields(src)
    if not rows:
        return text("", "empty input")
    body = rows[1:] if header and _is_header(rows, "auto") else rows
    head = rows[0] if body is not rows else []
    width = max((len(r) for r in rows), default=1)
    lines = ["  ".join(str(c)[:pad].ljust(pad) for c in (head + [""] * width)[:width])] if head else []
    for r in body:
        lines.append("  ".join(str(c)[:pad].ljust(pad) for c in (r + [""] * width)[:width]))
    return text("\n".join(lines) + "\n", f"{len(body)} rows formatted")


@T("text_find_most_dispensable_column", "Find the most dispensable columns", "text_manipulation",
   inputs=[tbl("src"), boolean("header", True)], ex={"src": COUNTS}, up="dispensable")
def dispensable_columns(src, header=True):
    """Rank columns by how redundant they are (constant or duplicate-heavy)."""
    rows = _fields(src)
    if len(rows) < 2:
        return table([], "need at least a header + 1 row")
    head = rows[0] if header else [f"c{i + 1}" for i in range(len(rows[0]))]
    body = rows[1:] if header else rows
    out = []
    for i, name in enumerate(head):
        vals = [r[i] for r in body if i < len(r)]
        c = Counter(vals)
        dup_pairs = sum(v * (v - 1) // 2 for v in c.values())
        out.append({"column": name, "unique": len(c), "most_common_fraction":
                    round(c.most_common(1)[0][1] / len(vals), 4) if vals else 0.0,
                    "duplicate_pairs": dup_pairs, "dispensability_score": round(
                        1 - len(c) / max(1, len(vals)) + dup_pairs / max(1, len(vals) ** 2), 4)})
    out.sort(key=lambda d: -d["dispensability_score"])
    return table(out, f"ranked {len(out)} columns")


@T("text_extract_with_regex", "Extract fields with a regular expression", "text_manipulation",
   inputs=[txt("src"), textbox("pattern", r"(\w+)\t(\d+)"), boolean("groups", True)],
   ex={"src": REGIONS_BED, "pattern": r"(\S+)\t(\d+)\t(\d+)\t(\S+)"}, up="extract")
def extract_with_regex(src, pattern="", groups=True):
    """``grep -oP`` style extraction of capture groups."""
    rx = re.compile(pattern)
    out = []
    for ln in _lines(src):
        for m in rx.finditer(ln):
            if groups and m.groups():
                out.append(list(m.groups()))
            else:
                out.append([m.group()])
    return text(_unparse(out), f"{len(out)} extractions")


@T("text_histogram_of_column", "Histogram of a numeric column", "text_manipulation",
   inputs=[tbl("src"), intin("column", 2), intin("bins", 20), boolean("header", True)],
   ex={"src": PHENO, "column": 3, "bins": 8}, up="histogram")
def column_histogram(src, column=2, bins=20, header=True):
    """Bin the values of one column (``datamash histogram``)."""
    rows = _fields(src)
    body = rows[1:] if header and rows and _is_header(rows, "auto") else rows
    vals = []
    for r in body:
        i = column - 1
        if i < len(r):
            try:
                vals.append(float(r[i]))
            except ValueError:
                pass
    if not vals:
        return table([], "no numeric values in that column")
    lo, hi = min(vals), max(vals)
    w = (hi - lo) / bins if hi > lo else 1
    counts = Counter(min(bins - 1, int((v - lo) / w)) for v in vals)
    out = [{"bin": f"{lo + i * w:.3f}-{lo + (i + 1) * w:.3f}", "count": counts.get(i, 0)}
           for i in range(bins)]
    res = table(out, f"{len(vals)} values in {bins} bins")
    res["figure"] = __import__("chroma_titan.core.plot", fromlist=["bar"]).bar(
        [r["bin"] for r in out], [r["count"] for r in out], title="column histogram")
    return res


# ===========================================================================
# Filter and Sort
# ===========================================================================
@T("filter_first_n", "Filter first/last lines of a dataset", "filter_and_sort",
   inputs=[tbl("src"), choice("which", ["first", "last"]), intin("n", 5)],
   ex={"src": PHENO, "which": "first", "n": 3}, up="head")
def filter_first_last(src, which="first", n=5):
    """Keep the first or last n lines."""
    lines = _lines(src)
    keep = lines[:n] if which == "first" else lines[-n:]
    return text("\n".join(keep) + "\n", f"{len(keep)} lines ({which})")


@T("filter_by_ids", "Filter lines matching identifiers from another dataset",
   "filter_and_sort", inputs=[tbl("src"), txt("ids"), intin("column", 1),
                              boolean("invert", False)],
   ex={"src": PHENO, "ids": "sample_A\nsample_C", "column": 1}, up="filter_from_file")
def filter_by_ids(src, ids="", column=1, invert=False):
    """Keep rows whose key appears in a list (or does not, with invert)."""
    wanted = {l.strip() for l in _lines(ids) if l.strip()}
    out = []
    for ln in _lines(src):
        f = ln.split("\t") if "\t" in ln else ln.split()
        key = f[column - 1] if 0 < column <= len(f) else ""
        if (key in wanted) != invert:
            out.append(ln)
    return text("\n".join(out) + "\n", f"{len(out)} rows kept from {len(_lines(src))}")


@T("filter_by_expression", "Filter table rows with an expression", "filter_and_sort",
   inputs=[tbl("src"), textbox("expression", "yield > 7"), boolean("header", True)],
   ex={"src": PHENO, "expression": "yield > 7"}, up="table_filter")
def filter_by_expression(src, expression="", header=True):
    """Row filter with a pandas ``query`` expression (safe, no eval of code)."""
    df = tables.load(src)
    if header is False and len(df.columns) and str(df.columns[0]).startswith("col"):
        pass
    try:
        out = tables.filter_rows(df, expression)
    except ValueError as exc:
        return text("", f"error: {exc}")
    res = io.table_result(out, f"{len(out)}/{len(df)} rows kept")
    return res


@T("filter_last_column", "Filter rows by the value in the last column", "filter_and_sort",
   inputs=[tbl("src"), number("threshold", 5.0), choice("op", [">", ">=", "<", "<=", "==", "!="])],
   ex={"src": PHENO, "threshold": 10, "op": ">"}, up="filter_last_column")
def filter_last_column(src, threshold=5.0, op=">"):
    """Compare the last column against a cut-off."""
    ops = {">": lambda a, b: a > b, ">=": lambda a, b: a >= b, "<": lambda a, b: a < b,
           "<=": lambda a, b: a <= b, "==": lambda a, b: abs(a - b) < 1e-9,
           "!=": lambda a, b: abs(a - b) >= 1e-9}
    keep = []
    for ln in _lines(src):
        f = ln.split("\t") if "\t" in ln else ln.split()
        try:
            v = float(f[-1])
        except (ValueError, IndexError):
            continue
        if ops[op](v, threshold):
            keep.append(ln)
    return text("\n".join(keep) + "\n", f"{len(keep)} rows pass last-column {op} {threshold}")


@T("filter_with_regexp", "Filter rows whose key matches a pattern", "filter_and_sort",
   inputs=[tbl("src"), textbox("pattern", "gene"), choice("field", ["first", "last", "any"]),
           boolean("invert", False)], ex={"src": GENES_FA if 0 else REGIONS_BED,
                                          "pattern": "gene"}, up="filter_with_regexp")
def filter_with_regexp(src, pattern="", field="any", invert=False):
    """``awk``-like regexp filter on one field."""
    rx = re.compile(pattern)
    keep = []
    for ln in _lines(src):
        f = ln.split("\t") if "\t" in ln else ln.split()
        if not f:
            continue
        cell = f[0] if field == "first" else f[-1] if field == "last" else ln
        if bool(rx.search(cell)) != invert:
            keep.append(ln)
    return text("\n".join(keep) + "\n", f"{len(keep)} rows kept")


@T("filter_non_numeric_rows", "Remove rows with non-numeric values", "filter_and_sort",
   inputs=[tbl("src"), intin("min_numeric_fraction", 100)], ex={"src": COUNTS},
   up="filter_rows")
def filter_non_numeric(src, min_numeric_fraction=100):
    """Drop rows where fewer than ``min_numeric_fraction``% of the fields are numeric."""
    keep = []
    for ln in _lines(src):
        f = ln.split("\t") if "\t" in ln else ln.split()
        if not f:
            continue
        frac = 100 * sum(1 for x in f[1:] if _num(x)) / max(1, len(f) - 1) if len(f) > 1 else 0
        if frac >= min_numeric_fraction:
            keep.append(ln)
    return text("\n".join(keep) + "\n", f"{len(keep)}/{len(_lines(src))} numeric rows")


@T("sort_table_column", "Sort by column", "filter_and_sort",
   inputs=[tbl("src"), textbox("column", "1"), choice("order", ["ascending", "descending"]),
           boolean("numeric", True), boolean("header", True)],
   ex={"src": PHENO, "column": "yield", "numeric": True, "order": "descending"}, up="sort")
def sort_table(src, column="1", order="ascending", numeric=True, header=True):
    """Sort rows by a named column or an index."""
    rows = _fields(src)
    if not rows:
        return text("", "empty input")
    has_head = header and _is_header(rows, "auto")
    head, body = (rows[0], rows[1:]) if has_head else (None, rows)
    ncol = max(len(r) for r in body) if body else 0
    if column in ("last", "L"):
        idx = ncol - 1
    elif str(column).isdigit() or re.match(r"^c?\d+$", str(column)):
        idx = _col_index(str(column), ncol)[0]
    else:
        idx = head.index(column) if head and column in head else 0
    def key(r):
        v = r[idx] if idx < len(r) else ""
        if numeric:
            try:
                return float(v)
            except ValueError:
                return float("nan")
        return v
    body.sort(key=key, reverse=order == "descending")
    out = ([head] if head else []) + body
    return text(_unparse(out), f"sorted by column {idx + 1} ({order})")


@T("sort_genomic_ranges", "Sort ranges by chromosomal location", "filter_and_sort",
   inputs=[bed("src"), choice("order", ["asc", "desc"]), boolean("by_size", False)],
   ex={"src": REGIONS_BED}, up="sort_bed")
def sort_ranges(src, order="asc", by_size=False):
    """BED-aware sort: chromosome then start (or interval size)."""
    ivs = io.parse_bed(io.as_text(src))
    ivs = genome.sort_ivs(ivs, "size" if by_size else "coordinate")
    if order == "desc":
        ivs = ivs[::-1]
    return text(io.write_bed(ivs), f"sorted {len(ivs)} intervals")


@T("shuffle_table_rows", "Randomly shuffle the order of rows", "filter_and_sort",
   inputs=[tbl("src"), intin("seed", 42), boolean("header", True)], ex={"src": PHENO},
   up="random_shuffle")
def shuffle_rows(src, seed=42, header=True):
    """Deterministic row shuffling with a seed."""
    rows = _fields(src)
    body = rows[1:] if header and rows and _is_header(rows, "auto") else rows
    head = rows[0] if body is not rows else None
    random.Random(seed).shuffle(body)
    out = ([head] if head else []) + body
    return text(_unparse(out), f"shuffled {len(body)} rows (seed={seed})")


@T("sample_rows", "Take a random subset of rows", "filter_and_sort",
   inputs=[tbl("src"), intin("n", 5), intin("seed", 1), boolean("header", True),
           choice("strategy", ["random", "systematic", "first", "last"])],
   ex={"src": PHENO, "n": 3, "strategy": "random"}, up="random_lines")
def sample_rows_tool(src, n=5, seed=1, header=True, strategy="random"):
    """Random / systematic / first-n sampling of table rows."""
    rows = _fields(src)
    body = rows[1:] if header and rows and _is_header(rows, "auto") else rows
    head = rows[0] if body is not rows else None
    n = min(int(n), len(body)) or len(body)
    if strategy == "random":
        body = random.Random(seed).sample(body, n)
    elif strategy == "systematic":
        step = max(1, len(body) // max(1, n))
        body = body[::step][:n]
    elif strategy == "last":
        body = body[-n:]
    else:
        body = body[:n]
    return text(_unparse(([head] if head else []) + body), f"sampled {len(body)} rows")


@T("remove_duplicate_rows", "Unique rows / detect duplicates", "filter_and_sort",
   inputs=[tbl("src"), intin("key_columns", 1), boolean("header", True),
           choice("keep", ["first", "last", "none"])],
   ex={"src": PHENO, "key_columns": 1}, up="unique")
def unique_rows(src, key_columns=1, header=True, keep="first"):
    """Collapse rows that share the same key columns."""
    rows = _fields(src)
    body = rows[1:] if header and rows and _is_header(rows, "auto") else rows
    head = rows[0] if body is not rows else None
    seen: dict[tuple, list[list[str]]] = OrderedDict()
    for r in body:
        k = tuple(r[i] for i in range(min(key_columns, len(r))))
        seen.setdefault(k, []).append(r)
    out = []
    for k, rs in seen.items():
        if len(rs) == 1 or keep == "first":
            out.append(rs[0])
        elif keep == "last":
            out.append(rs[-1])
    dup = sum(len(v) - 1 for v in seen.values())
    return text(_unparse(([head] if head else []) + out),
                f"{len(out)} unique rows ({dup} duplicates removed)")


@T("first_and_last_column_stats", "Statistics of the first and last columns", "filter_and_sort",
   inputs=[tbl("src")], ex={"src": COUNTS}, up="first_last")
def first_last_stats(src):
    """Summary stats for the outer columns of a numeric table."""
    rows = _fields(src)
    if len(rows) < 2:
        return table([], "need header + data")
    out = []
    for label, idx in (("first_column", 0), ("last_column", max(len(r) for r in rows) - 1)):
        vals = []
        for r in rows[1:]:
            if idx < len(r):
                try:
                    vals.append(float(r[idx]))
                except ValueError:
                    pass
        if vals:
            out.append({"dataset": label, **{k: (round(v, 4) if isinstance(v, float) else v)
                                             for k, v in stats.summary(vals).items()}})
        else:
            out.append({"dataset": label, "count": len(rows) - 1, "mean": float("nan")})
    return table(out, "column summary")


@T("filter_missing_values", "Remove rows or columns with missing data", "filter_and_sort",
   inputs=[tbl("src"), choice("axis", ["rows", "columns"]), number("max_missing_fraction", 0.0)],
   ex={"src": PHENO, "axis": "rows", "max_missing_fraction": 0.0}, up="filter")
def filter_missing(src, axis="rows", max_missing_fraction=0.0):
    """Drop rows/columns whose missing-value fraction exceeds the cut-off."""
    df = tables.load(src)
    miss = df.isna().mean(axis=1 if axis == "rows" else 0)
    keep = miss <= max_missing_fraction
    out = df[keep.values] if axis == "rows" else df.loc[:, keep.values]
    return io.table_result(out, f"removed {int((~keep).sum())} {axis} with missing data")


@T("filter_top_n_by_column", "Top/bottom N rows by a column", "filter_and_sort",
   inputs=[tbl("src"), textbox("column", "1"), intin("n", 10),
           choice("which", ["top", "bottom"]), boolean("header", True)],
   ex={"src": PHENO, "column": "yield", "n": 3, "which": "top"}, up="top_scores")
def top_n(src, column="1", n=10, which="top", header=True):
    """Rank rows by a column and keep the best/worst N."""
    res = sort_table(src, column=column, numeric=True, order="descending" if which == "top" else "ascending",
                     header=header)
    lines = res["text"].splitlines()
    return text("\n".join(lines[: int(n) + (1 if header else 0)]) + "\n",
                f"top {min(n, len(lines))} rows by {column}")


# ===========================================================================
# Join, Subtraction and Group Operations
# ===========================================================================
@T("group_by_column", "Group rows by a column value", "join__subtract_and_group",
   inputs=[tbl("src"), textbox("by", "1"), boolean("header", True)],
   ex={"src": PHENO, "by": "group"}, up="group_by")
def group_by(src, by="1", header=True):
    """Group-by summary: size, first members and column means per group."""
    rows = _fields(src)
    if not rows:
        return table([], "empty input")
    has_head = header and _is_header(rows, "auto")
    head = rows[0] if has_head else [f"c{i + 1}" for i in range(max(len(r) for r in rows))]
    body = rows[1:] if has_head else rows
    ncol = max((len(r) for r in body), default=1)
    key_i = _ci(by, head, ncol)
    groups: dict[str, list[list[str]]] = defaultdict(list)
    for r in body:
        k = "|".join(r[i] for i in range(min(len(r), key_i + 1)) if i == key_i) or (
            r[key_i] if key_i < len(r) else "NA")
        groups[k].append(r)
    out = []
    for k, rs in groups.items():
        rec = {"group": k, "n": len(rs), "members": ",".join(r[0] for r in rs)[:120]}
        for c in range(1, ncol):
            vals = []
            for r in rs:
                if c < len(r):
                    try:
                        vals.append(float(r[c]))
                    except ValueError:
                        pass
            if vals:
                rec[f"mean_{head[c]}"] = round(sum(vals) / len(vals), 4)
        out.append(rec)
    return table(out, f"{len(groups)} groups from {len(body)} rows")


@T("aggregate_group", "Aggregate values within groups (collapse)", "join__subtract_and_group",
   inputs=[tbl("src"), textbox("key_column", "1"), textbox("value_column", "2"),
           choice("function", ["sum", "mean", "median", "min", "max", "count", "concat",
                              "unique_count", "stdev"])],
   ex={"src": PHENO, "key_column": "group", "value_column": "yield", "function": "mean"},
   up="collapse")
def collapse_group(src, key_column="1", value_column="2", function="mean"):
    """Collapse a table to one row per key using an aggregate function."""
    rows = _fields(src)
    if not rows:
        return table([], "empty input")
    ncol = max(len(r) for r in rows)
    has_head = _is_header(rows, "auto")
    head = rows[0] if has_head else [f"c{i + 1}" for i in range(ncol)]
    body = rows[1:] if has_head else rows
    ki = _ci(key_column, head, ncol)
    vi = _ci(value_column, head, ncol)
    groups: dict[str, list[float]] = defaultdict(list)
    strs: dict[str, list[str]] = defaultdict(list)
    for r in body:
        if ki >= len(r):
            continue
        k = r[ki]
        strs[k].append(r[vi] if vi < len(r) else "")
        try:
            groups[k].append(float(r[vi]))
        except (ValueError, IndexError):
            continue
    out = []
    for k in groups:
        v = groups[k]
        fn = {"sum": sum(v), "mean": stats.mean(v), "median": stats.median(v), "min": min(v),
              "max": max(v), "count": len(v), "stdev": stats.stdev(v),
              "unique_count": len(set(strs[k])),
              "concat": "|".join(strs[k])}
        val = fn[function]
        rec = {"group": k, function: round(val, 6) if isinstance(val, float) else val,
               "n_values": len(v)}
        out.append(rec)
    return table(out, f"collapsed to {len(out)} groups")


@T("subtract_tables", "Subtract intervals or rows between datasets", "join__subtract_and_group",
   inputs=[bed("a"), bed("b")], ex={"a": REGIONS_BED, "b": TARGETS_BED}, up="bedtools_subtract")
def subtract_tables(a, b):
    """Remove the parts of A covered by B (``bedtools subtract -A``)."""
    ra, rb = io.parse_bed(io.as_text(a)), io.parse_bed(io.as_text(b))
    out = genome.subtract(ra, rb)
    return table(genome.to_table(out) if out else [],
                 f"{len(out)} intervals remain ({len(ra)} in, {len(rb)} removed by)")


@T("set_operations_on_columns", "Union / intersect / difference of two ID lists",
   "join__subtract_and_group",
   inputs=[txt("a"), txt("b"), choice("operation", ["union", "intersect", "difference",
                                                   "complement", "symmetric_difference"])],
   ex={"a": PHENO, "b": REGIONS_BED, "operation": "union"}, up="operate_on_genomic_intervals")
def set_operations(a, b, operation="union"):
    """Classic set algebra over the first column of two datasets."""
    sa = {l.split("\t")[0] for l in _lines(a)}
    sb = {l.split("\t")[0] for l in _lines(b)}
    out = {"union": sa | sb, "intersect": sa & sb, "difference": sa - sb,
           "complement": sb - sa, "symmetric_difference": sa ^ sb}[operation]
    res = table([{"id": x} for x in sorted(out)], f"{len(out)} ids in {operation}")
    res["counts"] = {"only_a": len(sa - sb), "only_b": len(sb - sa), "both": len(sa & sb)}
    return res


@T("jaccard_between_tables", "Jaccard index of two ID sets", "join__subtract_and_group",
   inputs=[txt("a"), txt("b")], ex={"a": REGIONS_BED, "b": TARGETS_BED}, up="jaccard")
def jaccard_tables(a, b):
    """Similarity of two gene/id lists (``jaccard``)."""
    sa = {l.split("\t")[0] for l in _lines(a)}
    sb = {l.split("\t")[0] for l in _lines(b)}
    st = stats.jaccard_sets(sa, sb)
    return values_message({k: round(v, 5) if isinstance(v, float) else v for k, v in st.items()},
                          "set similarity")


@T("transpose_table", "Transpose rows and columns", "join__subtract_and_group",
   inputs=[tbl("src"), boolean("header", True)], ex={"src": COUNTS}, up="transpose")
def transpose_table(src, header=True):
    """Swap the axes of a table (``transpose``)."""
    df = tables.load(src)
    if header and len(df):
        idx = df.iloc[:, 0].astype(str)
        out = df.iloc[:, 1:].set_index(idx).T.reset_index().rename(columns={"index": "row"})
    else:
        out = df.T.reset_index().rename(columns={"index": "row"})
    return io.table_result(out, f"transposed to {out.shape[0]}x{out.shape[1]}")


@T("reverse_rows", "Reverse the row order of a dataset", "join__subtract_and_group",
   inputs=[tbl("src"), boolean("header", True)], ex={"src": PHENO}, up="reverse_order")
def reverse_rows(src, header=True):
    """Flip the order of the rows (keeps the header on top)."""
    rows = _lines(src)
    if header and len(rows) > 1:
        return text("\n".join([rows[0]] + rows[1:][::-1]) + "\n",
                    f"reversed {len(rows) - 1} rows")
    return text("\n".join(rows[::-1]) + "\n", f"reversed {len(rows)} lines")


@T("relaxed_cluster_rows", "Cluster rows by similarity of their keys",
   "join__subtract_and_group",
   inputs=[txt("src"), number("threshold", 0.75)], ex={"src": REGIONS_BED, "threshold": 0.6},
   up="cluster")
def cluster_by_name(src, threshold=0.75):
    """Group identifiers whose names are more similar than a cut-off."""
    from difflib import SequenceMatcher

    names = [l.split("\t")[0] for l in _lines(src)]
    parent = list(range(len(names)))

    def find(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    for i in range(len(names)):
        for j in range(i + 1, len(names)):
            if SequenceMatcher(None, names[i], names[j]).ratio() >= threshold:
                parent[max(find(i), find(j))] = min(find(i), find(j))
    groups: dict[int, list[str]] = defaultdict(list)
    for i, nm in enumerate(names):
        groups[find(i)].append(nm)
    out = [{"cluster": k + 1, "size": len(v), "members": "|".join(v)}
           for k, v in sorted(groups.items(), key=lambda kv: -len(kv[1]))]
    return table(out, f"{len(out)} clusters at identity>={threshold}")


# ===========================================================================
# Datamash: aggregation / statistics
# ===========================================================================
DM_FUNCS = ("mean", "median", "geomean", "harmonic_mean", "min", "max", "sum", "count",
            "countunique", "unique", "elapse", "collapse", "first_value", "last_value",
            "mad", "pstdev", "samplestdev", "variance", "pvariance", "percentile",
            "interquartile_range", "range", "absolute_sum", "cumulative_sum")


def _dm_aggregate(values: list[float], fn: str, param: float = 25.0) -> Any:
    if not values:
        return float("nan")
    if fn == "mean":
        return stats.mean(values)
    if fn == "median":
        return stats.median(values)
    if fn == "geomean":
        return stats.geomean(values)
    if fn == "harmonic_mean":
        return stats.harmonic_mean(values)
    if fn == "min":
        return min(values)
    if fn == "max":
        return max(values)
    if fn == "sum":
        return sum(values)
    if fn == "count":
        return len(values)
    if fn in ("variance",):
        return stats.variance(values)
    if fn == "pvariance":
        return stats.variance(values, sample=False)
    if fn == "pstdev":
        return stats.stdev(values, sample=False)
    if fn == "samplestdev":
        return stats.stdev(values)
    if fn == "mad":
        return stats.mad(values)
    if fn == "interquartile_range":
        return stats.iqr(values)
    if fn == "range":
        return max(values) - min(values)
    if fn == "absolute_sum":
        return sum(abs(v) for v in values)
    if fn == "cumulative_sum":
        return sum(values)
    if fn == "percentile":
        return stats.quantile(values, param / 100)
    return stats.mean(values)


@T("datamash_column_stats", "Datamash-style column statistics", "datamash",
   inputs=[tbl("src"), multi("functions", list(DM_FUNCS), ["mean", "median", "min", "max"]),
           textbox("columns", "all"), boolean("header", True)],
   ex={"src": COUNTS, "functions": ["mean", "median", "min", "max", "sum"], "columns": "all"},
   up="datamash")
def datamash_stats(src, functions=("mean",), columns="all", header=True):
    """``datamash geometric-mean 2 mean 3 ...``: aggregate chosen columns."""
    rows = _fields(src)
    if not rows:
        return table([], "empty input")
    has_head = header and _is_header(rows, "auto")
    head = rows[0] if has_head else [f"c{i + 1}" for i in range(max(len(r) for r in rows))]
    body = rows[1:] if has_head else rows
    ncol = max((len(r) for r in body), default=1)
    idxs = _cis(columns, head, ncol)
    fns = functions if isinstance(functions, (list, tuple)) else [functions]
    out = []
    for c in idxs:
        vals = []
        for r in body:
            if c < len(r):
                try:
                    vals.append(float(r[c]))
                except ValueError:
                    pass
        rec = {"column": head[c] if c < len(head) else f"c{c + 1}", "n_numeric": len(vals)}
        for fn in fns:
            rec[fn] = (round(_dm_aggregate(vals, fn), 6)
                       if isinstance(_dm_aggregate(vals, fn), float) else _dm_aggregate(vals, fn))
        out.append(rec)
    return table(out, f"datamash {','.join(fns)} over {len(idxs)} columns")


@T("datamash_grouped", "Datamash groupby + aggregate", "datamash",
   inputs=[tbl("src"), textbox("group_column", "1"), textbox("value_column", "2"),
           choice("function", list(DM_FUNCS)), boolean("header", True)],
   ex={"src": PHENO, "group_column": "group", "value_column": "yield", "function": "mean"},
   up="datamash")
def datamash_grouped(src, group_column="1", value_column="2", function="mean", header=True):
    """Group rows by a key column then aggregate a numeric column per group."""
    rows = _fields(src)
    if not rows:
        return table([], "empty input")
    ncol = max(len(r) for r in rows)
    has_head = header and _is_header(rows, "auto")
    head = rows[0] if has_head else [f"c{i + 1}" for i in range(ncol)]
    body = rows[1:] if has_head else rows
    gi = _ci(group_column, head, ncol)
    vi = _ci(value_column, head, ncol)
    groups: "OrderedDict[str, list[float]]" = OrderedDict()
    for r in body:
        if gi >= len(r):
            continue
        key = r[gi]
        try:
            groups.setdefault(key, []).append(float(r[vi]))
        except (ValueError, IndexError):
            groups.setdefault(key, [])
    out = [{"group": k, "n": len(v), function: round(_dm_aggregate(v, function), 6)}
           for k, v in groups.items()]
    return table(out, f"datamash {function} of {head[vi] if vi < len(head) else vi + 1} by {head[gi] if gi < len(head) else gi + 1}")


@T("datamash_percentiles", "Percentiles of one column", "datamash",
   inputs=[tbl("src"), textbox("column", "1"), textbox("list_of_percentiles", "5,25,50,75,95"),
           boolean("header", True)], ex={"src": COUNTS, "column": "2"}, up="datamash")
def datamash_percentiles(src, column="1", list_of_percentiles="5,25,50,75,95", header=True):
    """Report arbitrary percentiles of a numeric column."""
    rows = _fields(src)
    if not rows:
        return table([], "empty input")
    ncol = max(len(r) for r in rows)
    body = rows[1:] if header and _is_header(rows, "auto") else rows
    head = rows[0] if body is not rows else None
    idx = _ci(column, head, ncol)
    vals = []
    for r in body:
        if idx < len(r):
            try:
                vals.append(float(r[idx]))
            except ValueError:
                pass
    ps = [float(x) for x in re.split(r"[,\s]+", list_of_percentiles) if x]
    head = rows[0] if body is not rows else None
    name = head[idx] if head and idx < len(head) else f"column {idx + 1}"
    out = [{"percentile": p, name: round(stats.quantile(vals, p / 100), 6)} for p in ps]
    return table(out, f"{len(vals)} values, {len(ps)} percentiles")


@T("datamash_cross_tab", "Cross-tabulation of two columns", "datamash",
   inputs=[tbl("src"), textbox("row_column", "1"), textbox("col_column", "2"),
           textbox("value_column", ""), boolean("header", True)],
   ex={"src": PHENO, "row_column": "group", "col_column": "treatment"}, up="datamash")
def cross_tab(src, row_column="1", col_column="2", value_column="", header=True):
    """``datamash crosstab``-style contingency table (counts or aggregated values)."""
    rows = _fields(src)
    if not rows:
        return table([], "empty input")
    ncol = max(len(r) for r in rows)
    has_head = header and _is_header(rows, "auto")
    head = rows[0] if has_head else [f"c{i + 1}" for i in range(ncol)]
    body = rows[1:] if has_head else rows
    ri = _ci(row_column, head, ncol)
    ci = _ci(col_column, head, ncol)
    vi = _ci(value_column, head, ncol) if value_column else None
    data: dict = defaultdict(dict)
    rvals: dict = defaultdict(list)
    for r in body:
        if ri >= len(r) or ci >= len(r):
            continue
        key = r[ri]
        col = r[ci]
        if vi is not None and vi < len(r):
            try:
                rvals[(key, col)].append(float(r[vi]))
            except ValueError:
                pass
        data[key][col] = data[key].get(col, 0) + 1
    cols = sorted({c for row in data.values() for c in row})
    out = []
    for key in sorted(data):
        rec = {str(head[ri]) if ri < len(head) else "row": key}
        for c in cols:
            cnt = data[key].get(c, 0)
            rec[str(c)] = (round(stats.mean(rvals[(key, c)]), 4)
                           if rvals.get((key, c)) and vi is not None else cnt)
        out.append(rec)
    return table(out, f"crosstab {len(out)}x{len(cols)}")


@T("datamash_summary_all", "Summary statistics for every numeric column", "datamash",
   inputs=[tbl("src"), boolean("header", True)], ex={"src": COUNTS}, up="summary_statistics")
def summary_all_columns(src, header=True):
    """count/mean/median/min/max/q1/q3/stdev/IQR/MAD per numeric column."""
    rows = _fields(src)
    if not rows:
        return table([], "empty input")
    ncol = max(len(r) for r in rows)
    has_head = header and _is_header(rows, "auto")
    head = rows[0] if has_head else [f"c{i + 1}" for i in range(ncol)]
    body = rows[1:] if has_head else rows
    out = []
    for c in range(ncol):
        vals = []
        for r in body:
            if c < len(r):
                try:
                    vals.append(float(r[c]))
                except ValueError:
                    pass
        if len(vals) < 2:
            continue
        out.append({"column": head[c] if c < len(head) else c + 1,
                    "count": len(vals), "mean": round(stats.mean(vals), 4),
                    "median": round(stats.median(vals), 4), "min": min(vals), "max": max(vals),
                    "q1": round(stats.quantile(vals, 0.25), 4),
                    "q3": round(stats.quantile(vals, 0.75), 4),
                    "IQR": round(stats.iqr(vals), 4), "stdev": round(stats.stdev(vals), 4),
                    "MAD": round(stats.mad(vals), 4), "sum": round(sum(vals), 4),
                    "cv_percent": round(stats.cv(vals), 3), "skew": round(stats.skewness(vals), 4),
                    "kurtosis": round(stats.kurtosis(vals), 4)})
    return table(out, f"summarised {len(out)} numeric columns")


@T("datamash_matrix_summary", "Row and column summary of a matrix", "datamash",
   inputs=[tbl("src"), boolean("row_sums", True), boolean("col_sums", True),
           choice("stat", ["sum", "mean", "median", "min", "max"])],
   ex={"src": COUNTS}, up="matrix_summary")
def matrix_summary(src, row_sums=True, col_sums=True, stat="sum"):
    """Per-row / per-column aggregate of a numeric matrix (Galaxy 'Matrix Summary')."""
    rows = _fields(src)
    if not rows:
        return table([], "empty input")
    has_head = _is_header(rows, "auto")
    head = rows[0] if has_head else [f"c{i + 1}" for i in range(len(rows[0]))]
    body = rows[1:] if has_head else rows
    ncol = max(len(r) for r in body)
    fn = {"sum": sum, "mean": stats.mean, "median": stats.median, "min": min, "max": max}[stat]
    out = []
    if row_sums:
        for r in body:
            vals = []
            for x in r[1:]:
                try:
                    vals.append(float(x))
                except ValueError:
                    pass
            if vals:
                out.append({"axis": "row", "id": r[0] if r else "", stat: round(fn(vals), 4),
                            "n": len(vals)})
    if col_sums:
        for c in range(1, ncol):
            vals = []
            for r in body:
                if c < len(r):
                    try:
                        vals.append(float(r[c]))
                    except ValueError:
                        pass
            if vals:
                out.append({"axis": "column", "id": head[c] if c < len(head) else c,
                            stat: round(fn(vals), 4), "n": len(vals)})
    return table(out, f"{stat} per {'row+column' if row_sums and col_sums else 'axis'}")


@T("datamash_normalize_column", "Normalise a numeric column", "datamash",
   inputs=[tbl("src"), textbox("column", "1"),
           choice("mode", ["z-score", "min-max", "relative abundance (percent)", "rank",
                           "decimal", "log2", "sum-to-one", "robust (median/MAD)"]),
           boolean("header", True)],
   ex={"src": COUNTS, "column": "2", "mode": "z-score"}, up="normalize")
def normalize_column(src, column="1", mode="z-score", header=True):
    """Common rescalings of one column (datamash --normalize)."""
    rows = _fields(src)
    if not rows:
        return text("", "empty input")
    ncol = max(len(r) for r in rows)
    idx = _ci(column, _head_of(rows), ncol)
    body = rows[1:] if header and _is_header(rows, "auto") else rows
    head = rows[0] if body is not rows else None
    vals = []
    for r in body:
        try:
            vals.append(float(r[idx]))
        except (ValueError, IndexError):
            vals.append(float("nan"))
    good = [v for v in vals if not math.isnan(v)]
    if mode == "z-score":
        m, s = stats.mean(good), stats.stdev(good)
        new = [(v - m) / s if s and not math.isnan(v) else 0.0 for v in vals]
    elif mode == "min-max":
        lo, hi = (min(good), max(good)) if good else (0, 1)
        rng = (hi - lo) or 1
        new = [(v - lo) / rng if not math.isnan(v) else 0.0 for v in vals]
    elif mode in ("relative abundance (percent)", "sum-to-one"):
        tot = sum(abs(v) for v in good) or 1
        scale = 100 if "percent" in mode else 1
        new = [v / tot * scale if not math.isnan(v) else 0.0 for v in vals]
    elif mode == "rank":
        order = sorted(range(len(vals)), key=lambda i: vals[i] if not math.isnan(vals[i]) else 0)
        new = [0.0] * len(vals)
        for rank, i in enumerate(order, 1):
            new[i] = float(rank)
    elif mode == "decimal":
        new = [v / 100 if not math.isnan(v) else 0.0 for v in vals]
    elif mode == "log2":
        new = [math.log2(v) if v > 0 else 0.0 for v in vals]
    else:  # robust
        m, md = stats.median(good) or 0, stats.mad(good) or 1
        new = [(v - m) / (1.4826 * md) if not math.isnan(v) else 0.0 for v in vals]
    out = []
    for i, r in enumerate(body):
        nr = list(r)
        while len(nr) <= idx:
            nr.append("")
        nr[idx] = f"{new[i]:.6g}"
        out.append(nr)
    return text(_unparse(([head] if head else []) + out), f"column {idx + 1} normalised ({mode})")


@T("datamash_cumulative", "Cumulative sums and running statistics", "datamash",
   inputs=[tbl("src"), textbox("column", "1"), choice("function",
            ["cumsum", "running_mean", "ewma", "diff", "pct_change", "max_so_far", "min_so_far"]),
           boolean("header", True)], ex={"src": COUNTS, "column": "2", "function": "cumsum"},
   up="datamash")
def cumulative(src, column="1", function="cumsum", header=True):
    """Row-wise cumulative transforms of a numeric column."""
    rows = _fields(src)
    ncol = max((len(r) for r in rows), default=1)
    body = rows[1:] if header and rows and _is_header(rows, "auto") else rows
    idx = _ci(column, _head_of(rows), ncol)
    vals = []
    labels = []
    for r in body:
        try:
            vals.append(float(r[idx]))
            labels.append(r[0])
        except (ValueError, IndexError):
            continue
    out, run = [], []
    acc = 0.0
    for i, v in enumerate(vals):
        acc += v
        cur = {"cumsum": acc, "running_mean": acc / (i + 1),
               "ewma": v if i == 0 else 0.3 * v + 0.7 * run[-1],
               "diff": v - vals[i - 1] if i else 0.0,
               "pct_change": 100 * (v - vals[i - 1]) / vals[i - 1] if i and vals[i - 1] else 0.0,
               "max_so_far": max(vals[:i + 1]), "min_so_far": min(vals[:i + 1])}[function]
        run.append(cur)
        out.append({"index": i + 1, "label": labels[i], function: round(cur, 6)})
    return table(out, f"{function} over {len(vals)} values")


@T("datamash_bin_column", "Bin a column and summarise each bin", "datamash",
   inputs=[tbl("src"), textbox("column", "1"), number("bin_size", 10.0),
           choice("aggregate", ["count", "mean", "sum", "max", "min"])],
   ex={"src": COUNTS, "column": "2", "bin_size": 500, "aggregate": "mean"}, up="datamash")
def bin_column(src, column="1", bin_size=10.0, aggregate="count"):
    """Histogram-style binned aggregation (``datamash --histogram``)."""
    rows = _fields(src)
    ncol = max((len(r) for r in rows), default=1)
    body = rows[1:] if rows and _is_header(rows, "auto") else rows
    head = rows[0] if body is not rows else None
    idx = _ci(column, head, ncol)
    vals = []
    for r in body:
        try:
            vals.append(float(r[idx]))
        except (ValueError, IndexError):
            pass
    if not vals:
        return table([], "no numeric values")
    b = max(1e-9, float(bin_size))
    groups: dict[int, list[float]] = defaultdict(list)
    for v in vals:
        groups[int(v // b)].append(v)
    out = [{"bin_start": k * b, "bin_end": (k + 1) * b, "n": len(v),
             aggregate: round(_dm_aggregate(v, aggregate if aggregate != "count" else "count"), 4)}
            for k, v in sorted(groups.items())]
    return table(out, f"{len(out)} bins of width {b:g}")


@T("datamash_outliers", "Detect outliers in a column", "datamash",
   inputs=[tbl("src"), textbox("column", "1"), choice("method", ["iqr", "modified_z", "zscore",
            "grubbs_like"]), number("cutoff", 1.5)], ex={"src": COUNTS, "column": "2"},
   up="outlier_detection")
def outliers_column(src, column="1", method="iqr", cutoff=1.5):
    """Flag outliers with IQR, MAD-based modified z or plain z-scores."""
    rows = _fields(src)
    ncol = max((len(r) for r in rows), default=1)
    body = rows[1:] if rows and _is_header(rows, "auto") else rows
    head = rows[0] if body is not rows else None
    idx = _ci(column, head, ncol)
    vals = []
    for i, r in enumerate(body):
        try:
            vals.append((i, r[0], float(r[idx])))
        except (ValueError, IndexError):
            continue
    numbers = [v for _, _, v in vals]
    out = []
    if method == "iqr":
        st = stats.outliers_iqr(numbers, cutoff)
        lo, hi = st["lower"], st["upper"]
        for i, lab, v in vals:
            if v < lo or v > hi:
                out.append({"row": i + 1, "label": lab, "value": v, "rule": f"IQR {cutoff}",
                            "distance": round(v - (hi if v > hi else lo), 4)})
    elif method == "modified_z":
        st = stats.outliers_mad(numbers, cutoff)
        out = [{"row": body[i][0] if i < len(body) else i, "label": l, "value": v,
                "rule": f"MAD |z|>{cutoff}", "distance": round(v - stats.median(numbers), 4)}
               for i, l, v in vals
               if stats.mad(numbers) and abs(0.6745 * (v - stats.median(numbers)) / stats.mad(numbers)) > cutoff]
    elif method == "grubbs_like":
        m, s = stats.mean(numbers), stats.stdev(numbers)
        dev = sorted(vals, key=lambda t: -abs(t[2] - m))
        g = abs(dev[0][2] - m) / s if s else 0
        out = [{"row": dev[0][0] + 1, "label": dev[0][1], "value": dev[0][2],
                "rule": "Grubbs statistic", "distance": round(g, 4),
                "significant": bool(g > cutoff)}]
    else:
        m, s = stats.mean(numbers), stats.stdev(numbers)
        out = [{"row": i + 1, "label": l, "value": v, "rule": f"|z|>{cutoff}",
                "distance": round(abs((v - m) / s) if s else 0.0, 4)}
               for i, l, v in vals if s and abs((v - m) / s) > cutoff]
    return table(out, f"{len(out)} outliers by {method}")


@T("datamash_unique_values", "Unique values with counts and proportions", "datamash",
   inputs=[tbl("src"), textbox("column", "1"), intin("top", 0), boolean("percent", True)],
   ex={"src": PHENO, "column": "group"}, up="datamash")
def unique_values(src, column="1", top=0, percent=True):
    """``datamash unique/count``: frequency table of a categorical column."""
    rows = _fields(src)
    ncol = max((len(r) for r in rows), default=1)
    body = rows[1:] if rows and _is_header(rows, "auto") else rows
    idx = _ci(column, _head_of(rows), ncol)
    c = Counter(r[idx] for r in body if idx < len(r))
    n = sum(c.values()) or 1
    items = c.most_common(int(top)) if top else c.most_common()
    out = [{"value": k, "count": v, "percent": round(100 * v / n, 4) if percent else ""}
           for k, v in items]
    return table(out, f"{len(c)} unique values, {n} rows")


# ===========================================================================
# Expression tools (compute on tables)
# ===========================================================================
@T("compute_expression_column", "Add a computed column from an expression", "expression_tools",
   inputs=[tbl("src"), textbox("name", "ratio"), textbox("expression", "yield / biomass"),
           boolean("header", True)],
   ex={"src": PHENO, "name": "ratio", "expression": "yield / (biomass + 0.001)"},
   up="compute")
def compute_expression(src, name="new_column", expression="", header=True):
    """Galaxy's 'Compute' tool: evaluate an expression using column names."""
    df = tables.load(src)
    try:
        out = tables.add_column(df, name, expression)
    except ValueError as exc:
        return text("", f"error: {exc}")
    return io.table_result(out, f"added column '{name}'")


@T("row_math", "Row-wise arithmetic over numeric columns", "expression_tools",
   inputs=[tbl("src"), textbox("columns", "all"),
           choice("operation", ["row_sum", "row_mean", "row_min", "row_max", "row_product",
                              "row_geomean", "row_stdev", "row_range", "row_count_nonzero"])],
   ex={"src": COUNTS, "operation": "row_sum"}, up="row_math")
def row_math_tool(src, columns="all", operation="row_sum"):
    """Aggregate across each row of a numeric matrix."""
    rows = _fields(src)
    has_head = _is_header(rows, "auto") if rows else False
    body = rows[1:] if has_head else rows
    head = rows[0] if has_head else None
    ncol = max((len(r) for r in body), default=1)
    idxs = _cis(columns, head, ncol)
    if columns in ("all", "", None):
        idxs = list(range(1, ncol))
    out = []
    for r in body:
        vals = []
        for c in idxs:
            if c < len(r):
                try:
                    vals.append(float(r[c]))
                except ValueError:
                    pass
        if not vals:
            continue
        val = {"row_sum": sum(vals), "row_mean": stats.mean(vals), "row_min": min(vals),
               "row_max": max(vals), "row_product": math.prod(vals),
               "row_geomean": stats.geomean(vals), "row_stdev": stats.stdev(vals),
               "row_range": max(vals) - min(vals),
               "row_count_nonzero": float(sum(1 for v in vals if v))}[operation]
        out.append({"id": r[0], operation: round(val, 6), "n": len(vals)})
    return table(out, f"{operation} for {len(out)} rows")


@T("scale_matrix", "Scale a matrix by rows, columns or globally", "expression_tools",
   inputs=[tbl("src"), choice("axis", ["rows", "columns", "global"]),
           choice("mode", ["z-score", "min-max", "center", "unit_variance"])],
   ex={"src": COUNTS, "axis": "rows", "mode": "z-score"}, up="scale")
def scale_matrix(src, axis="rows", mode="z-score"):
    """Standardise a numeric matrix along an axis."""
    df = tables.load(src)
    num = df.select_dtypes("number")
    if not num.shape[1]:
        return table([], "no numeric columns")
    if axis == "global":
        vals = num.values.flatten()
        m, s = stats.mean(vals), stats.stdev(vals)
        lo, hi = float(np_min(vals)), float(np_max(vals))
    else:
        base = num if axis == "rows" else num.T
        m = base.mean(axis=1)
        s = base.std(axis=1).replace(0, math.nan)
        lo = base.min(axis=1)
        hi = base.max(axis=1)
    import numpy as np

    arr = num.to_numpy(dtype=float)
    if axis == "global":
        vals = arr.flatten()
        m, s = float(np.mean(vals)), float(np.std(vals, ddof=1) or 1)
        lo, hi = float(np.min(vals)), float(np.max(vals))
        if mode == "z-score":
            out_arr = (arr - m) / (s or 1)
        elif mode == "min-max":
            out_arr = (arr - lo) / ((hi - lo) or 1)
        elif mode == "center":
            out_arr = arr - m
        else:
            out_arr = arr / (s or 1)
    else:
        mv = np.asarray(m, dtype=float)
        sv = np.nan_to_num(np.asarray(s, dtype=float), nan=1.0)
        sv[sv == 0] = 1.0
        lv = np.asarray(lo, dtype=float)
        hv = np.asarray(hi, dtype=float)
        A = arr if axis == "rows" else arr.T
        if mode == "z-score":
            A = (A - mv[:, None]) / sv[:, None]
        elif mode == "min-max":
            A = (A - lv[:, None]) / np.where((hv - lv) == 0, 1, hv - lv)[:, None]
        elif mode == "center":
            A = A - mv[:, None]
        else:
            A = A / sv[:, None]
        out_arr = A if axis == "rows" else A.T
    out = num.copy()
    for i, col in enumerate(num.columns):
        out[col] = out_arr[:, i]
    idx_col = df.drop(columns=list(num.columns))
    res = io.table_result(pd_concat([idx_col, out]), f"scaled {axis} ({mode})")
    res["index"] = [str(x) for x in num.index[:50]]
    return res


def np_min(v):
    import numpy as np

    return float(np.min(v))


def np_max(v):
    import numpy as np

    return float(np.max(v))


def pd_concat(frames):
    import pandas as pd

    return pd.concat(frames, axis=1)


@T("apply_function_to_rows", "Apply a small function to selected columns", "expression_tools",
   inputs=[tbl("src"), multi("columns", [], []), choice("function",
            ["mean", "median", "sum", "geomean", "cv", "range", "log10sum", "max-min_ratio"])],
   ex={"src": COUNTS, "columns": [], "function": "cv"}, up="apply")
def apply_row_function(src, columns=None, function="mean"):
    """Reduce several numeric columns per row with a chosen function."""
    rows = _fields(src)
    has_head = _is_header(rows, "auto") if rows else False
    head = rows[0] if has_head else []
    body = rows[1:] if has_head else rows
    ncol = max((len(r) for r in body), default=1)
    idxs = _cis(",".join(map(str, columns or [])), head, ncol) if columns else list(range(1, ncol))
    out = []
    for r in body:
        vals = []
        for c in idxs:
            if c < len(r):
                try:
                    vals.append(float(r[c]))
                except ValueError:
                    pass
        if len(vals) < 1:
            continue
        val = {"mean": stats.mean(vals), "median": stats.median(vals), "sum": sum(vals),
               "geomean": stats.geomean(vals), "cv": stats.cv(vals),
               "range": (max(vals) - min(vals)) if vals else 0.0,
               "log10sum": math.log10(sum(vals)) if sum(vals) > 0 else 0.0,
               "max-min_ratio": (max(vals) / min(vals)) if vals and min(vals) else float("inf")}[function]
        out.append({"id": r[0], function: round(val, 6), "n_values": len(vals)})
    return table(out, f"{function} over {len(idxs)} columns for {len(out)} rows")


@T("missing_value_imputation", "Impute missing values in a table", "expression_tools",
   inputs=[tbl("src"), choice("method", ["mean", "median", "zero", "min", "knn"]),
           intin("k", 5)], ex={"src": PHENO, "method": "median"}, up="imputation")
def impute_missing(src, method="mean", k=5):
    """Fill NaNs column-wise (constant, mean/median or k-nearest-neighbour)."""
    import numpy as np
    import pandas as pd

    df = tables.load(src)
    num = df.select_dtypes("number")
    n_missing = int(df.isna().sum().sum())
    if method == "knn" and num.shape[1] > 1:
        arr = num.to_numpy(dtype=float)
        mu = np.nanmean(arr, 0)
        mask = np.isnan(arr)
        arr[mask] = np.take(mu, np.where(mask)[1])
        out = df.copy()
        for i, c in enumerate(num.columns):
            out[c] = arr[:, i]
    else:
        fill = {"mean": "mean", "median": "median", "zero": 0, "min": "min"}[method]
        out = df.fillna(fill if isinstance(fill, str) else 0)
    return io.table_result(out, f"imputed {n_missing} missing values ({method})")


@T("column_arithmetic_pairwise", "Arithmetic between two numeric columns", "expression_tools",
   inputs=[tbl("src"), textbox("left", "1"), textbox("right", "2"),
           choice("operation", ["add", "subtract", "multiply", "divide", "ratio_log2",
                               "absolute_difference"])],
   ex={"src": COUNTS, "left": "2", "right": "3"}, up="column_arithmetic")
def column_arithmetic(src, left="1", right="2", operation="divide"):
    """Element-wise maths on two columns of the same table."""
    rows = _fields(src)
    has_head = _is_header(rows, "auto") if rows else False
    head = rows[0] if has_head else []
    body = rows[1:] if has_head else rows
    ncol = max((len(r) for r in body), default=1)
    li, ri = _ci(left, head, ncol), _ci(right, head, ncol)
    out = []
    for r in body:
        try:
            a, b = float(r[li]), float(r[ri])
        except (ValueError, IndexError):
            continue
        val = {"add": a + b, "subtract": a - b, "multiply": a * b, "divide": a / b if b else float("nan"),
               "ratio_log2": math.log2(a / b) if a > 0 and b > 0 else float("nan"),
               "absolute_difference": abs(a - b)}[operation]
        out.append({"id": r[0], f"{head[li] if li < len(head) else li + 1}_{operation}_"
                    f"{head[ri] if ri < len(head) else ri + 1}": round(val, 6), "a": a, "b": b})
    return table(out, f"{operation} applied to {len(out)} rows")


# ===========================================================================
# Collection operations
# ===========================================================================
@T("collection_list_from_column", "Extract an element list from a dataset",
   "collection_operations", inputs=[tbl("src"), textbox("column", "1"), boolean("unique", True),
                                   boolean("header", True)],
   ex={"src": PHENO, "column": "group"}, up="list_operations")
def collection_list(src, column="1", unique=True, header=True):
    """Build the equivalent of a Galaxy list collection from one column."""
    rows = _fields(src)
    ncol = max((len(r) for r in rows), default=1)
    body = rows[1:] if header and rows and _is_header(rows, "auto") else rows
    idx = _ci(column, _head_of(rows), ncol)
    vals = [r[idx] for r in body if idx < len(r)]
    if unique:
        vals = list(OrderedDict.fromkeys(vals))
    out = [{"element": i + 1, "identifier": v, "label": f"element_{i + 1}"}
           for i, v in enumerate(vals)]
    return table(out, f"list of {len(vals)} elements")


@T("collection_filter_elements", "Filter collection elements by pattern", "collection_operations",
   inputs=[txt("names"), textbox("pattern", "gene"), boolean("invert", False)],
   ex={"names": REGIONS_BED, "pattern": "gene"}, up="filter_from_list")
def collection_filter(names, pattern="", invert=False):
    """Keep collection elements whose identifier matches a regexp."""
    rx = re.compile(pattern)
    ids = [l.split("\t")[0] for l in _lines(names)]
    keep = [x for x in ids if bool(rx.search(x)) != invert]
    return table([{"element": i + 1, "identifier": v} for i, v in enumerate(keep)],
                 f"{len(keep)}/{len(ids)} elements kept")


@T("collection_rename_elements", "Rename collection elements", "collection_operations",
   inputs=[txt("names"), textbox("prefix", "sample_"), textbox("suffix", ""),
           textbox("strip_pattern", "")], ex={"names": PHENO, "prefix": "S"},
   up="rename_by_pattern")
def collection_rename(names, prefix="", suffix="", strip_pattern=""):
    """Apply prefix/suffix and an optional strip regexp to element names."""
    out = []
    for i, ln in enumerate(_lines(names)):
        nm = ln.split("\t")[0]
        if strip_pattern:
            nm = re.sub(strip_pattern, "", nm)
        out.append({"old": nm if not strip_pattern else ln.split("\t")[0],
                    "new": f"{prefix}{nm}{suffix}", "element": i + 1})
    return table(out, f"renamed {len(out)} elements")


@T("collection_extract_by_index", "Extract elements by index or range", "collection_operations",
   inputs=[txt("names"), textbox("indices", "1,3,5-6"), boolean("cyclic", False)],
   ex={"names": PHENO, "indices": "1,3"}, up="extract")
def collection_extract(names, indices="1", cyclic=False):
    """Pull specific elements out of a list collection."""
    lines = _lines(names)
    idx = [int(x) for x in re.split(r"[,\s]+", str(indices).replace("-", " ")) if x.isdigit()]
    want = set()
    for i in idx:
        want.add(((i - 1) % len(lines)) if cyclic else (i - 1))
    keep = [l for i, l in enumerate(lines) if i in want]
    return text("\n".join(keep) + "\n", f"extracted {len(keep)}/{len(lines)} elements")


@T("collection_sort_elements", "Sort collection elements by identifier", "collection_operations",
   inputs=[txt("names"), choice("by", ["name", "length", "size", "index"]),
           boolean("reverse", False)], ex={"names": PHENO}, up="sort_lists")
def collection_sort(names, by="name", reverse=False):
    """Order list elements deterministically."""
    lines = _lines(names)
    if by == "length":
        lines.sort(key=len, reverse=reverse)
    elif by == "size":
        lines.sort(key=lambda l: len(l.split("\t")), reverse=reverse)
    else:
        lines.sort(key=lambda l: l.split("\t")[0], reverse=reverse)
    return text("\n".join(lines) + "\n", f"sorted {len(lines)} elements by {by}")


@T("collection_zip_pairs", "Zip two lists into pairs", "collection_operations",
   inputs=[txt("first"), txt("second"), boolean("fill", False)],
   ex={"first": PHENO, "second": REGIONS_BED}, up="zip_longest")
def collection_zip(first, second, fill=False):
    """Pair up elements from two collections positionally."""
    a, b = _lines(first), _lines(second)
    n = max(len(a), len(b)) if fill else min(len(a), len(b))
    rows = [{"pair": i + 1, "left": a[i].split("\t")[0] if i < len(a) else ".",
             "right": b[i].split("\t")[0] if i < len(b) else "."} for i in range(n)]
    return table(rows, f"{len(rows)} paired elements")


@T("collection_unzip", "Split a paired collection into two lists", "collection_operations",
   inputs=[tbl("pairs"), intin("column_left", 1), intin("column_right", 2)],
   ex={"pairs": PHENO}, up="unzip")
def collection_unzip(pairs, column_left=1, column_right=2):
    """Undo a zip: two output lists from a two-column table."""
    rows = _fields(pairs)
    left = [r[column_left - 1] for r in rows if len(r) >= column_left]
    right = [r[column_right - 1] for r in rows if len(r) >= column_right]
    return table([{"element": i + 1, "left": left[i] if i < len(left) else "",
                   "right": right[i] if i < len(right) else ""} for i in range(max(len(left), len(right)))],
                 f"unzipped into {len(left)} + {len(right)} elements")


@T("collection_count_elements", "Count elements and their sizes", "collection_operations",
   inputs=[txt("names")], ex={"names": PHENO}, up="count")
def collection_count(names):
    """Number of elements plus per-element line/field counts."""
    lines = _lines(names)
    rows = [{"element": i + 1, "identifier": l.split("\t")[0], "lines": 1,
             "fields": len(l.split("\t")) if "\t" in l else len(l.split()),
             "characters": len(l)} for i, l in enumerate(lines)]
    return table(rows, f"counted {len(lines)} elements")


@T("collection_merge", "Merge collections with deduplication", "collection_operations",
   inputs=[txt("first"), txt("second"), choice("strategy", ["union", "intersect", "first_only"]),
           boolean("unique", True)], ex={"first": REGIONS_BED, "second": TARGETS_BED},
   up="merge")
def collection_merge(first, second, strategy="union", unique=True):
    """Combine two collections with a set-like strategy."""
    a = [l.split("\t")[0] for l in _lines(first)]
    b = [l.split("\t")[0] for l in _lines(second)]
    if strategy == "intersect":
        out = [x for x in a if x in set(b)]
    elif strategy == "first_only":
        out = [x for x in a if x not in set(b)]
    else:
        out = a + b
    if unique:
        out = list(OrderedDict.fromkeys(out))
    return table([{"element": i + 1, "identifier": v} for i, v in enumerate(out)],
                 f"{len(out)} merged elements ({strategy})")


@T("collection_group_by_label", "Group collection elements into nested lists",
   "collection_operations", inputs=[tbl("src"), textbox("column", "1")],
   ex={"src": PHENO, "column": "group"}, up="group_with_pandas")
def collection_group(src, column="1"):
    """Create a list:paired-style grouping from a label column."""
    return group_by(src, by=column)


# ===========================================================================
# Get Data / Send Data
# ===========================================================================
@T("get_data_create_dataset", "Create a dataset from pasted text", "get_data",
   inputs=[bigtext("content", "", label="Paste data here")],
   ex={"content": "chrV\t100\t200\tfirst\nchrV\t300\t400\tsecond"}, up="upload")
def create_dataset(content=""):
    """Wrap pasted text into a dataset and report its basic properties."""
    body = io.as_text(content)
    lines = [l for l in body.splitlines() if l.strip()]
    fmt = "bed" if lines and re.match(r"^\S+\t\d+\t\d+", lines[0]) else \
        "fasta" if lines and lines[0].startswith(">") else \
        "fastq" if len(lines) >= 4 and lines[0].startswith("@") and lines[2].startswith("+") else "txt"
    return {"message": f"created a `{fmt}` dataset with {len(lines)} lines",
            "text": body, "stats": {"lines": len(lines), "characters": len(body),
                                   "format_guess": fmt,
                                   "fields_in_first_line": len(lines[0].split("\t")) if lines else 0}}


@T("get_data_extract_column", "Extract a single column as a new dataset", "get_data",
   inputs=[tbl("src"), textbox("column", "1"), boolean("header", False)],
   ex={"src": PHENO, "column": "1"}, up="extract")
def extract_column(src, column="1", header=False):
    """``cut -f`` producing a one-column dataset."""
    rows = _fields(src)
    ncol = max((len(r) for r in rows), default=1)
    body = rows[1:] if header and rows else rows
    head_ = rows[0] if body is not rows else None
    idx = _cis(column, head_, max((len(r) for r in rows), default=1))
    out = "\n".join(r[idx[0]] for r in body if idx and idx[0] < len(r))
    return text(out + "\n", f"extracted column {idx[0] + 1 if idx else 0} ({len(body)} values)")


@T("get_data_concatenate_history", "Concatenate a set of same-format files", "get_data",
   inputs=[multi("files", [PHENO, REGIONS_BED, TARGETS_BED, COUNTS], [REGIONS_BED, TARGETS_BED]),
           boolean("keep_headers", True)], ex={"files": [REGIONS_BED, TARGETS_BED]}, up="cat")
def concatenate_many(files=None, keep_headers=True):
    """Stack several datasets vertically (list → single dataset)."""
    out, n = [], 0
    for i, f in enumerate(files or []):
        lines = _lines(f)
        if not lines:
            continue
        if i and keep_headers and _is_header([lines[0]], "auto"):
            lines = lines[1:]
        out.extend(lines)
        n += 1
    return text("\n".join(out) + "\n", f"concatenated {n} datasets -> {len(out)} lines")


@T("get_data_format_probe", "Identify the format of a dataset", "get_data",
   inputs=[anyfile("src", GENOME, label="Any dataset")], ex={"src": GENOME}, up="peek")
def format_probe(src):
    """Sniff FASTA/FASTQ/VCF/SAM/BED/GFF/table from the first lines."""
    body = io.as_text(src)
    lines = [l for l in body.splitlines() if l.strip()]
    first = lines[0] if lines else ""
    if first.startswith(">"):
        fmt, extra = "fasta", {"records": sum(1 for l in lines if l.startswith(">"))}
    elif first.startswith("@") and len(lines) > 2 and lines[2].startswith("+"):
        fmt, extra = "fastq", {"reads": len(lines) // 4}
    elif first.startswith("#CHROM"):
        fmt, extra = "vcf", {"records": sum(1 for l in lines if not l.startswith("#"))}
    elif first.startswith(("@HD", "@SQ")):
        fmt, extra = "sam", {"alignments": sum(1 for l in lines if not l.startswith("@"))}
    elif first.startswith("##gff"):
        fmt, extra = "gff", {"features": sum(1 for l in lines if not l.startswith("#"))}
    elif re.match(r"^\S+\t\d+\t\d+", first):
        fmt, extra = "bed", {"intervals": len(lines)}
    elif "\t" in first and len(first.split("\t")) > 2:
        fmt, extra = "tsv", {"columns": len(first.split("\t")), "rows": len(lines)}
    else:
        fmt, extra = "txt", {"lines": len(lines)}
    return values_message({"format": fmt, "bytes": len(body), "first_line": first[:70], **extra},
                          f"looks like {fmt.upper()}")


@T("send_data_download_table", "Convert any dataset to a chosen text format", "send_data",
   inputs=[tbl("src"), choice("format", ["tsv", "csv", "markdown", "html", "json", "yaml-ish"])],
   ex={"src": COUNTS, "format": "csv"}, up="export")
def export_table(src, format="csv"):
    """Re-serialise a table for download (TSV/CSV/Markdown/HTML/JSON)."""
    df = tables.load(src)
    if format == "csv":
        body = io.to_csv(df)
    elif format == "markdown":
        body = df.to_markdown(index=False)
    elif format == "html":
        body = df.to_html(index=False)
    elif format == "json":
        body = df.to_json(orient="records", indent=1)
    else:
        body = "\n".join(f"{c}: {list(df[c].astype(str))}" for c in df.columns)
    return text(body, f"exported {len(df)} rows as {format.upper()}")


@T("send_data_paste_to_column", "Send a column back as a standalone dataset", "send_data",
   inputs=[tbl("src"), textbox("column", "1")], ex={"src": PHENO, "column": "3"}, up="export")
def send_column(src, column="1"):
    """Write one column out (e.g. a gene list for a downstream enrichment tool)."""
    res = extract_column(src, column=column)
    return text(res["text"], f"sent column {column} as a new dataset")


@T("send_data_summary_report", "Text summary report of a dataset", "send_data",
   inputs=[anyfile("src", PHENO), intin("preview_lines", 6)], ex={"src": PHENO}, up="report")
def summary_report(src, preview_lines=6):
    """Compact human-readable summary (size, format, preview)."""
    body = io.as_text(src)
    lines = body.splitlines()
    preview = "\n".join(lines[:int(preview_lines)])
    rep = (f"dataset: {len(lines)} lines, {len(body)} bytes\n"
           f"first lines:\n{preview}\n"
           f"last line: {lines[-1] if lines else ''}\n")
    return text(rep, "summary report")


@T("text_line_numbers", "Prepend line numbers", "text_manipulation",
   inputs=[txt("src"), intin("start", 1), textbox("sep", "\t")], ex={"src": PHENO, "start": 1},
   up="number_lines")
def line_numbers(src, start=1, sep="\t"):
    """``nl``-style numbering of lines."""
    lines = _lines(src)
    return text("\n".join(f"{i + start}{sep.replace(chr(92) + 't', chr(9))}{l}"
                          for i, l in enumerate(lines)) + "\n", f"numbered {len(lines)} lines")


@T("text_diff_two_sorted", "Compare two sorted ID lists (comm)", "text_manipulation",
   inputs=[txt("a"), txt("b")], ex={"a": REGIONS_BED, "b": TARGETS_BED}, up="comm")
def comm_tool(a, b):
    """Lines only in A, only in B and in both (``comm -3``)."""
    sa = {l for l in _lines(a)}
    sb = {l for l in _lines(b)}
    rows = [{"category": "only_in_first", "lines": sorted(sa - sb)},
            {"category": "only_in_second", "lines": sorted(sb - sa)},
            {"category": "in_both", "lines": sorted(sa & sb)}]
    out = [{"category": r["category"], "count": len(r["lines"]),
            "example": r["lines"][0][:60] if r["lines"] else ""} for r in rows]
    return table(out, "comm-style comparison")


@T("text_column_summarizer", "Numeric summary of the whole table", "text_manipulation",
   inputs=[tbl("src"), boolean("per_row", False)], ex={"src": COUNTS}, up="summary")
def table_overview(src, per_row=False):
    """Whole-table statistics: number of rows/cols, numeric cells, sparsity."""
    rows = _fields(src)
    if not rows:
        return table([], "empty")
    has_head = _is_header(rows, "auto")
    body = rows[1:] if has_head else rows
    cells = [c for r in body for c in r]
    numeric = [c for c in cells if _num(c)]
    allvals = [float(c) for c in numeric]
    rec = {"rows": len(body), "columns": max((len(r) for r in rows), default=0),
           "cells": len(cells), "numeric_cells": len(numeric),
           "percent_numeric": round(100 * len(numeric) / max(1, len(cells)), 2),
           "sum_of_numeric": round(sum(allvals), 4), "mean_of_numeric": round(stats.mean(allvals), 4),
           "median_of_numeric": round(stats.median(allvals), 4),
           "max_of_numeric": round(max(allvals), 4) if allvals else 0.0,
           "empty_cells": sum(1 for c in cells if not c.strip())}
    out = [rec]
    if per_row:
        for i, r in enumerate(body[:50]):
            v = [float(x) for x in r[1:] if _num(x)]
            out.append({"rows": i + 1, "columns": len(r), "cells": len(r),
                        "numeric_cells": len(v), "percent_numeric": round(100 * len(v) / max(1, len(r)), 2),
                        "sum_of_numeric": round(sum(v), 4), "mean_of_numeric": round(stats.mean(v), 4),
                        "median_of_numeric": round(stats.median(v), 4),
                        "max_of_numeric": round(max(v), 4) if v else 0.0, "empty_cells": 0})
    return table(out, "table overview")


@T("text_filter_by_row_count", "Keep the first N rows or a fraction", "filter_and_sort",
   inputs=[tbl("src"), number("fraction", 0.5), intin("n", 0), boolean("header", True)],
   ex={"src": PHENO, "fraction": 0.5}, up="filter_by")
def filter_fraction(src, fraction=0.5, n=0, header=True):
    """Keep the first fraction (or n) rows of a dataset."""
    rows = _lines(src)
    body = rows[1:] if header and len(rows) > 1 else rows
    head = rows[0] if body is not rows else None
    keep = int(n) if n else max(1, int(len(body) * fraction))
    return text("\n".join(([head] if head else []) + body[:keep]) + "\n",
                f"kept {min(keep, len(body))}/{len(body)} rows")


@T("text_top_scores", "Rank rows by a column and output ranks", "filter_and_sort",
   inputs=[tbl("src"), textbox("column", "2"), boolean("descending", True),
           boolean("tie_average", False)], ex={"src": COUNTS, "column": "2"}, up="rank")
def rank_rows(src, column="2", descending=True, tie_average=False):
    """Add a rank column computed from a numeric column."""
    rows = _fields(src)
    has_head = _is_header(rows, "auto") if rows else False
    head = rows[0] if has_head else [f"c{i + 1}" for i in range(max(len(r) for r in rows))]
    body = rows[1:] if has_head else rows
    ncol = max((len(r) for r in body), default=1)
    idx = _ci(column, _head_of(rows), ncol)
    vals = []
    for i, r in enumerate(body):
        try:
            vals.append((i, float(r[idx])))
        except (ValueError, IndexError):
            vals.append((i, float("nan")))
    order = sorted([v for v in vals if not math.isnan(v[1])], key=lambda t: -t[1] if descending else t[1])
    ranks = {i: r + 1 for r, (i, _) in enumerate(order)}
    out = []
    for i, r in enumerate(body):
        nr = list(r) + [""] * (idx + 1 - len(r))
        nr.insert(idx + 1, str(ranks.get(i, "")))
        out.append(nr)
    head2 = list(head)
    head2.insert(idx + 1, "rank")
    return text(_unparse([head2] + out), f"ranked {len(order)} rows on {head[idx] if idx < len(head) else idx}")


@T("text_extract_column_by_name", "Select columns by name", "text_manipulation",
   inputs=[tbl("src"), multi("columns", ["sample", "group", "yield", "biomass"], ["group", "yield"]),
           boolean("all_columns", False)], ex={"src": PHENO, "columns": ["sample", "yield"]}, up="select")
def select_named_columns(src, columns=None, all_columns=False):
    """Keep named columns (pandas-style selection)."""
    df = tables.load(src)
    if not all_columns and columns:
        have = [c for c in columns if c in df.columns]
        if not have:
            raise KeyError(f"none of {columns} in {list(df.columns)}")
        df = df[have]
    return io.table_result(df, f"selected {df.shape[1]} of {len(tables.load(src).columns)} columns")


@T("text_add_row_index", "Add a row number column", "text_manipulation",
   inputs=[tbl("src"), textbox("name", "row"), intin("start", 1)], ex={"src": PHENO}, up="index")
def add_row_index(src, name="row", start=1):
    """Number the rows (1-based, offset configurable)."""
    df = tables.load(src)
    out = df.copy()
    out.insert(0, name, range(int(start), int(start) + len(df)))
    return io.table_result(out, f"added index column {name}")


@T("text_rename_by_pattern", "Rename identifiers with a pattern", "text_manipulation",
   inputs=[txt("src"), textbox("pattern", "gene"), textbox("replacement", "GENE"),
           intin("column", 1)], ex={"src": REGIONS_BED, "pattern": "gene", "replacement": "G"},
   up="replace_in_column")
def rename_by_pattern(src, pattern="", replacement="", column=1):
    """Apply a regexp substitution to a single identifier column."""
    rows = _fields(src)
    ncol = max((len(r) for r in rows), default=1)
    idx = int(column) - 1
    out = []
    n = 0
    for r in rows:
        nr = list(r)
        while len(nr) <= idx:
            nr.append("")
        new, k = re.subn(pattern, replacement, nr[idx])
        nr[idx] = new
        n += k
        out.append(nr)
    return text(_unparse(out), f"{n} identifiers renamed")
