"""File format readers/writers and path resolution.

Every tool receives *text or a path*: :func:`as_text` transparently reads a file
(a real path, a path relative to the repository root, or a Streamlit upload) so
tools can be unit-tested with inline strings and used in the UI with uploads.

Supported lightweight formats: FASTA, FASTQ (Sanger/Illumina), BED (3-12),
GFF/GTF, VCF 4.x, SAM, narrowPeak/bedGraph, pairs, matrix tables, GMT, Newick,
PFM/PWM and delimited tables.  ``.gz`` inputs are transparently uncompressed.
"""

from __future__ import annotations

import bz2
import gzip
import io as _io
import os
import re
import tempfile
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Iterable

REPO_ROOT = Path(__file__).resolve().parents[2]
EXAMPLE_DIR = REPO_ROOT / "data" / "examples"


def example_path(name: str) -> str:
    p = EXAMPLE_DIR / name
    return str(p)


def resolve(src: Any) -> str:
    """Return a readable local path for a path/uploaded-file object."""
    if src is None:
        raise ValueError("no file given")
    if hasattr(src, "getvalue"):  # Streamlit UploadedFile / StringIO
        data = src.getvalue()
        suffix = os.path.splitext(getattr(src, "name", "file.txt"))[1]
        tmp = Path(tempfile.mkdtemp(prefix="chroma_")) / f"upload{suffix}"
        tmp.write_bytes(data if isinstance(data, bytes) else data.encode())
        return str(tmp)
    if isinstance(src, (bytes, bytearray)):
        tmp = Path(tempfile.mkdtemp(prefix="chroma_")) / "upload.txt"
        tmp.write_bytes(bytes(src))
        return str(tmp)
    s = str(src).strip()
    if s.startswith("@"):  # @examples/genome.fa
        return example_path(s[1:].split("/", 1)[-1])
    p = Path(s).expanduser()
    if not p.exists():
        cand = REPO_ROOT / s
        if cand.exists():
            return str(cand)
        cand2 = EXAMPLE_DIR / p.name
        if cand2.exists():
            return str(cand2)
        raise FileNotFoundError(f"no such file: {s}")
    return str(p)


def open_text(src: Any):
    p = resolve(src) if not isinstance(src, str) or looks_like_path(src) else None
    if p is None:
        return _io.StringIO(src)
    if p.endswith(".gz"):
        return _io.TextIOWrapper(gzip.open(p, "rb"), encoding="utf-8", errors="replace")
    if p.endswith(".bz2"):
        return _io.TextIOWrapper(bz2.open(p, "rb"), encoding="utf-8", errors="replace")
    return open(p, encoding="utf-8", errors="replace")


PATH_RE = re.compile(r"\.(fa|fasta|faa|fna|seq|gb|gbk|fastq|fq|fqgz|bed|bedgraph|broadpeak|narrowpeak|gff|gff3|gtf|vcf|sam|bam|csv|tsv|txt|table|matrix|counts|gmt|newick|nwk|nw|pfm|pwm|icm|meme|pairs|interact|mcool|h5|json|md|yaml|yml|pep|proteins|id|lst|list|loci|pe|wig|bed6)\b", re.I)


def looks_like_path(s: Any) -> bool:
    """True for single-line strings that point at a readable file."""
    if not isinstance(s, str):
        return False
    t = s.strip()
    if not t or "\n" in t:
        return False
    if t.startswith("@"):
        return True
    p = Path(t).expanduser()
    if p.exists() and p.is_file():
        return True
    if (REPO_ROOT / t).exists():
        return True
    return bool(PATH_RE.search(t)) and (EXAMPLE_DIR / p.name).exists()


def as_text(src: Any) -> str:
    """Text payload of a string, path, upload or bytes object."""
    if src is None:
        return ""
    if isinstance(src, (bytes, bytearray)):
        return bytes(src).decode("utf-8", "replace")
    if hasattr(src, "read"):
        data = src.read()
        return data.decode("utf-8", "replace") if isinstance(data, bytes) else str(data)
    if hasattr(src, "getvalue"):
        data = src.getvalue()
        return data.decode("utf-8", "replace") if isinstance(data, bytes) else str(data)
    s = str(src)
    if looks_like_path(s):
        try:
            with open_text(s) as fh:
                return fh.read()
        except OSError:
            return s
    return s


def lines(src: Any, strip: bool = True) -> list[str]:
    out = as_text(src).splitlines()
    return [ln.strip() for ln in out] if strip else out


def read_table(src: Any, sep: str | None = None, header: Any = "infer",
               names: list[str] | None = None, comment: str = "#") -> Any:
    """Delimited table -> :class:`pandas.DataFrame`, numeric columns coerced."""
    import pandas as pd

    txt = as_text(src)
    if not txt.strip():
        return pd.DataFrame()
    body = [ln for ln in txt.splitlines() if ln.strip() and not ln.startswith("#")]
    if sep is None:
        if "\t" in body[0]:
            sep = "\t"
        elif "," in body[0]:
            sep = ","
        else:
            sep = None  # whitespace
    first = body[0].lstrip("#").strip()
    ncols = len(first.split(sep)) if sep else len(first.split())
    if header is None:
        header = None
    else:
        try:
            f = float(first.split(sep or None)[0]) if sep else float(first.split()[0])
            del f
            header = None
        except (ValueError, IndexError):
            header = 0
    if header is None and names is None:
        names = [f"col{i + 1}" for i in range(ncols)]
    df = pd.read_csv(
        _io.StringIO("\n".join(ln.lstrip("#") for ln in body)),
        sep=sep or r"\s+",
        header=header,
        names=None if header == 0 else names,
        engine="python",
    )
    for c in df.columns:
        conv = pd.to_numeric(df[c], errors="coerce")
        if conv.notna().sum() >= max(1, int(0.6 * len(df))):
            df[c] = conv
    return df


def write_table(df: Any, sep: str = "\t", index: bool = False) -> str:
    import pandas as pd

    if isinstance(df, pd.DataFrame):
        return df.to_csv(sep=sep, index=index).rstrip("\n")
    if isinstance(df, dict):
        return pd.DataFrame(df).to_csv(sep=sep, index=False).rstrip("\n")
    return str(df)


def to_csv(df: Any, index: bool = False) -> str:
    return write_table(df, sep=",", index=index)


def to_tsv(df: Any, index: bool = False) -> str:
    return write_table(df, sep="\t", index=index)


# ---------------------------------------------------------------------------
# FASTA
# ---------------------------------------------------------------------------
@dataclass
class Seq:
    """A named sequence (a FASTA record)."""

    id: str
    seq: str
    desc: str = ""

    @property
    def header(self) -> str:
        return f"{self.id} {self.desc}".strip()

    @property
    def description(self) -> str:
        return self.desc

    def __len__(self) -> int:
        return len(self.seq)

    def __repr__(self) -> str:  # pragma: no cover
        return f"Seq({self.id!r}, len={len(self.seq)})"


def parse_fasta(text: str) -> list[Seq]:
    recs: list[Seq] = []
    header, buf = None, []
    for ln in as_text(text).splitlines():
        if ln.startswith(">"):
            if header is not None:
                recs.append(_mkseq(header, "".join(buf)))
            header, buf = ln[1:].strip(), []
        elif header is not None:
            buf.append(ln.strip())
    if header is not None:
        recs.append(_mkseq(header, "".join(buf)))
    return recs


def _mkseq(header: str, seq: str) -> Seq:
    parts = header.split(None, 1)
    return Seq(parts[0], seq.strip(), parts[1] if len(parts) > 1 else "")


def write_fasta(recs: Iterable[Seq | tuple[str, str]], width: int = 60,
                format: str = "fasta") -> str:
    out = []
    for r in recs:
        if isinstance(r, Seq):
            head, s = r.header, r.seq
        else:
            head, s = str(r[0]), str(r[1])
        if not s:
            continue
        if width > 0:
            body = "\n".join(s[i: i + width] for i in range(0, len(s), width))
        else:
            body = s
        if format == "fasta":
            out.append(f">{head}\n{body}")
        elif format == "faqual":
            out.append(f";{head}\n{body}\n")
        elif format == "solexa" or format == "fastq":
            out.append(f"@{head}\n{s}\n+\n{'I' * len(s)}")
        else:
            out.append(f">{head}\n{body}")
    return "\n".join(out) + "\n" if out else ""


def seqs_or_text(src: Any) -> list[Seq]:
    """FASTA records from a text/path blob; falls back to a single record."""
    txt = as_text(src)
    recs = parse_fasta(txt)
    if recs:
        return recs
    clean = "".join(txt.split())
    return [Seq("sequence_1", clean, "pasted sequence")] if clean else []


# ---------------------------------------------------------------------------
# FASTQ
# ---------------------------------------------------------------------------
@dataclass
class Read:
    id: str
    seq: str
    qual: str
    desc: str = ""

    @property
    def header(self) -> str:
        return f"{self.id} {self.desc}".strip()

    @property
    def phred(self) -> list[int]:
        return [ord(c) - 33 for c in self.qual]

    def __len__(self) -> int:
        return len(self.seq)


def parse_fastq(text: str, offset: int | None = None) -> list[Read]:
    ln = [l for l in as_text(text).splitlines() if l.strip()]
    reads: list[Read] = []
    i = 0
    while i + 3 < len(ln) + 1 and i + 3 < len(ln):
        head = ln[i]
        if not head.startswith("@"):
            i += 1
            continue
        seq, plus, qual = ln[i + 1], ln[i + 2], ln[i + 3]
        off = offset if offset is not None else _guess_offset(qual)
        if qual and off != 33:
            qual = "".join(chr(min(93, max(33, ord(c) - off + 33))) for c in qual)
        parts = head[1:].split(None, 1)
        reads.append(Read(parts[0], seq.strip(), qual.strip(),
                          parts[1] if len(parts) > 1 else ""))
        i += 4
    return reads


def _guess_offset(qual: str) -> int:
    """Illumina-1.3+ (33) vs Illumina-1.0/Solexa (64) from character range."""
    if not qual:
        return 33
    mx = max(ord(c) for c in qual)
    mn = min(ord(c) for c in qual)
    return 33 if mn >= 33 and mx <= 74 else (64 if mx > 74 else 33)


def write_fastq(reads: Iterable[Read], width: int = 0) -> str:
    out = []
    for r in reads:
        out.append(f"@{r.header}\n{r.seq}\n+\n{r.qual}")
    return "\n".join(out) + "\n" if out else ""


def reads_or_text(src: Any) -> list[Read]:
    txt = as_text(src)
    rd = parse_fastq(txt)
    if rd:
        return rd
    recs = parse_fasta(txt)
    if recs:
        return [Read(r.id, r.seq, "I" * len(r.seq), r.desc) for r in recs]
    s = "".join(txt.split())
    if s and all(c in "ACGTNacgtn" for c in s[:20]):
        return [Read("read_1", s.upper(), "I" * len(s))]
    return rd


# ---------------------------------------------------------------------------
# BED / intervals
# ---------------------------------------------------------------------------
@dataclass
class Interval:
    chrom: str
    start: int
    end: int
    name: str = "."
    score: float = 0.0
    strand: str = "."
    extra: list[str] = field(default_factory=list)

    @property
    def length(self) -> int:
        return self.end - self.start

    def as_fields(self) -> list[str]:
        return [self.chrom, str(self.start), str(self.end), self.name,
                _fmt_num(self.score), self.strand, *self.extra]

    def key(self) -> tuple:
        return (self.chrom, self.start, self.end, self.name)


def _fmt_num(x: Any) -> str:
    """Format a numeric field, writing '.' for missing / NaN / infinite values."""
    try:
        f = float(x)
    except (TypeError, ValueError):
        return str(x)
    if f != f or f == float("inf") or f == float("-inf"):
        return "."
    return str(int(f)) if f == int(f) else f"{f:g}"


def parse_bed(text: str) -> list[Interval]:
    out: list[Interval] = []
    for ln in as_text(text).splitlines():
        ln = ln.strip()
        if not ln or ln.startswith(("#", "track", "browser")):
            continue
        p = ln.split("\t")
        if len(p) < 3:
            p = ln.split()
        if len(p) < 3:
            continue
        try:
            s, e = int(float(p[1])), int(float(p[2]))
        except ValueError:
            continue
        out.append(Interval(p[0], min(s, e), max(s, e),
                            p[3] if len(p) > 3 else ".",
                            float(p[4]) if len(p) > 4 and _isnum(p[4]) else 0.0,
                            p[5] if len(p) > 5 else ".",
                            p[6:]))
    return out


def write_bed(intervals: Iterable[Interval], bed12: bool = True) -> str:
    out = []
    for iv in intervals:
        f = iv.as_fields()
        if not bed12:
            f = f[:min(len(f), 6)]
        # keep a numeric score field, drop empties
        while len(f) > 6 and f[-1] in ("", "."):
            f.pop()
        out.append("\t".join(str(x) for x in f))
    return "\n".join(out) + "\n" if out else ""


def intervals_or_bed(src: Any) -> list[Interval]:
    txt = as_text(src)
    ivs = parse_bed(txt)
    if ivs:
        return ivs
    # accept "chr:start-end[+/-]" style region strings
    for ln in txt.splitlines():
        m = re.match(r"^(\S+?):(\d+)-(\d+)([+-])?$", ln.strip())
        if m:
            ivs.append(Interval(m.group(1), int(m.group(2)), int(m.group(3)),
                                ".", 0, m.group(4) or "."))
    return ivs


def parse_bedgraph(text: str) -> list[tuple[str, int, int, float]]:
    out = []
    for ln in as_text(text).splitlines():
        if not ln.strip() or ln.startswith(("#", "track", "variableStep", "fixedStep")):
            continue
        p = ln.split()
        if len(p) < 4:
            continue
        try:
            out.append((p[0], int(float(p[1])), int(float(p[2])), float(p[3])))
        except ValueError:  # not a bedGraph/WIG line (e.g. a genome size table)
            continue
    return out


def write_bedgraph(vals: Iterable[tuple[str, int, int, float]]) -> str:
    return "\n".join(f"{c}\t{s}\t{e}\t{_fmt_num(v)}" for c, s, e, v in vals) + "\n"


def parse_wig(text: str) -> list[tuple[str, int, int, float]]:
    out, cur = [], None
    for ln in as_text(text).splitlines():
        ln = ln.strip()
        m = re.match(r"^(variableStep|fixedStep)\s+(\S+)\s+start=(\d+)(?:\s+step=(\d+))?", ln)
        if m:
            cur = [m.group(2), int(m.group(3)), int(m.group(4) or 1)]
            continue
        if not ln or ln.startswith(("#", "track", "browser")):
            continue
        p = ln.split()
        if cur:
            pos = cur[1]
            if len(p) >= 2 and _isnum(p[0]):
                pos = int(float(p[0]))
                out.append((cur[0], pos, pos + 1, float(p[1])))
            elif len(p) >= 1 and _isnum(p[0]):
                out.append((cur[0], pos, pos + 1, float(p[0])))
            cur[1] = pos + cur[2]
        elif len(p) >= 4:
            out.append((p[0], int(float(p[1])), int(float(p[2])), float(p[3])))
    return out


# ---------------------------------------------------------------------------
# GFF / GTF
# ---------------------------------------------------------------------------
@dataclass
class GffRow:
    seqid: str
    source: str
    type: str
    start: int
    end: int
    score: float
    strand: str
    phase: str
    attrs: dict[str, str] = field(default_factory=dict)

    @property
    def length(self) -> int:
        return self.end - self.start + 1

    def get(self, key: str, default: str = "") -> str:
        return self.attrs.get(key, default)

    def as_line(self) -> str:
        return "\t".join([
            self.seqid, self.source, self.type, str(self.start), str(self.end),
            _fmt_num(self.score), self.strand, str(self.phase),
            ";".join(f"{k}={v}" for k, v in self.attrs.items())])


def parse_gff(text: str) -> list[GffRow]:
    out: list[GffRow] = []
    for ln in as_text(text).splitlines():
        if not ln.strip() or ln.startswith("#"):
            continue
        p = ln.split("\t")
        if len(p) < 8:
            p = ln.split()
        if len(p) < 8:
            continue
        attrs: dict[str, str] = {}
        for kv in p[8].split(";") if len(p) > 8 else []:
            kv = kv.strip()
            if not kv:
                continue
            if "=" in kv:
                k, v = kv.split("=", 1)
                attrs[k.strip()] = v.strip()
            else:  # GTF style: key "value";
                m = re.match(r'(\S+)\s+"?([^"]*)"?', kv)
                if m:
                    attrs[m.group(1)] = m.group(2)
        try:
            s, e = int(float(p[3])), int(float(p[4]))
        except ValueError:
            continue
        out.append(GffRow(p[0], p[1], p[2], min(s, e), max(s, e),
                          float(p[5]) if _isnum(p[5]) else float("nan"),
                          p[6], p[7], attrs))
    return out


def write_gff(rows: Iterable[GffRow]) -> str:
    return "##gff-version 3\n" + "\n".join(r.as_line() for r in rows) + "\n"


def gff_to_bed12(rows: Iterable[GffRow]) -> list[Interval]:
    """Collapse transcript rows into BED12 with exon blocks."""
    by_tx: dict[str, list[GffRow]] = {}
    for r in rows:
        if r.type.lower() in ("mrna", "transcript"):
            key = r.get("ID") or r.get("transcript_id") or r.get("Name") or "."
            by_tx.setdefault(key, []).append(r)
    out: list[Interval] = []
    exons: dict[str, list[GffRow]] = {}
    for r in rows:
        if r.type.lower() in ("exon", "cds"):
            parent = r.get("Parent") or r.get("transcript_id") or r.get("gene_id") or "."
            exons.setdefault(parent, []).append(r)
    for tx, rs in by_tx.items() or ():
        for r in rs:
            ex = sorted(exons.get(tx, []), key=lambda x: x.start)
            if not ex:
                ex = [r]
            start, end = r.start - 1, r.end
            blk = [f"{e.end - e.start + 1}" for e in ex]
            offs = [f"{e.start - start}" for e in ex]
            out.append(Interval(r.seqid, start, end, tx, 0, r.strand, [
                str(r.start - 1), str(r.end), "0", str(len(ex)),
                ",".join(blk) + ",", ",".join(offs) + ","]))
    return out


# ---------------------------------------------------------------------------
# VCF
# ---------------------------------------------------------------------------
@dataclass
class VcfRecord:
    chrom: str
    pos: int
    id: str
    ref: str
    alt: str
    qual: float
    filter: str
    info: dict[str, str] = field(default_factory=dict)
    fmt_keys: list[str] = field(default_factory=list)
    samples: list[dict[str, str]] = field(default_factory=list)

    @property
    def is_snv(self) -> bool:
        return len(self.ref) == 1 and all(len(a) == 1 for a in self.alt.split(","))

    @property
    def alts(self) -> list[str]:
        return [a for a in self.alt.split(",") if a and a != "."]

    @property
    def n_alt(self) -> int:
        return len(self.alts)

    @property
    def is_indel(self) -> bool:
        for a in self.alts:
            if len(a) != len(self.ref):
                return True
        return False

    @property
    def is_transition(self) -> bool:
        if not self.is_snv:
            return False
        ts = {("A", "G"), ("G", "A"), ("C", "T"), ("T", "C")}
        return any((self.ref, a) in ts for a in self.alts)

    def gt(self, i: int = 0) -> str:
        return self.samples[i].get("GT", "./.") if i < len(self.samples) else "./."

    def field(self, i: int, key: str) -> str:
        return self.samples[i].get(key, ".") if i < len(self.samples) else "."

    def info_num(self, key: str, default: float = float("nan")) -> float:
        v = self.info.get(key)
        try:
            return float(v.split(",")[0]) if v else default
        except (ValueError, AttributeError):
            return default

    def as_line(self) -> str:
        return "\t".join([
            self.chrom, str(self.pos), self.id, self.ref, self.alt,
            _fmt_num(self.qual) if self.qual == self.qual else ".",
            self.filter,
            ";".join(f"{k}={v}" if v != "" else k for k, v in self.info.items()) or ".",
        ] + (["."] + self.fmt_keys + [":"])[0:0] + (
            [":".join([self.fmt_keys[0]] if self.fmt_keys else [])] if False else []) + [])

    def format_line(self) -> str:
        head = [self.chrom, str(self.pos), self.id, self.ref, self.alt,
                _fmt_num(self.qual) if self.qual == self.qual else ".", self.filter,
                ";".join(f"{k}={v}" if v != "" else k for k, v in self.info.items()) or "."]
        if not self.fmt_keys:
            return "\t".join(head)
        cols = [head[i] if i < len(head) else "." for i in range(len(head))]
        return "\t".join(cols + [":".join(self.fmt_keys)] + [
            ":".join(s.get(k, ".") for k in self.fmt_keys) for s in self.samples])


@dataclass
class Vcf:
    header: list[str]
    samples: list[str]
    records: list[VcfRecord]

    def as_text(self) -> str:
        hdr = [h for h in self.header]
        if self.samples:
            hdr = [h for h in hdr if not h.startswith("#CHROM")]
            hdr.append("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"
                       + "\t".join(self.samples))
        lines_ = hdr + [r.format_line() for r in self.records]
        return "\n".join(lines_) + "\n"


def parse_vcf(text: str) -> Vcf:
    header: list[str] = []
    samples: list[str] = []
    recs: list[VcfRecord] = []
    for ln in as_text(text).splitlines():
        if ln.startswith("##"):
            header.append(ln.rstrip())
            continue
        if ln.startswith("#CHROM"):
            header.append(ln.rstrip())
            samples = ln.rstrip("\n").split("\t")[9:]
            continue
        if not ln.strip() or ln.startswith("#"):
            continue
        p = ln.rstrip("\n").split("\t")
        if len(p) < 8:
            p = ln.split()[:8] + ln.split()[8:]
        while len(p) < 8:
            p.append(".")
        info: dict[str, str] = {}
        if p[7] != ".":
            for kv in p[7].split(";"):
                if not kv:
                    continue
                if "=" in kv:
                    k, v = kv.split("=", 1)
                    info[k] = v
                else:
                    info[kv] = ""
        keys = p[8].split(":") if len(p) > 8 and p[8] not in (".", "") else []
        samps = []
        for cell in p[9:]:
            vals = cell.split(":")
            samps.append({k: (vals[i] if i < len(vals) else ".")
                          for i, k in enumerate(keys)})
        recs.append(VcfRecord(p[0], int(float(p[1])), p[2], p[3], p[4],
                              float(p[5]) if _isnum(p[5]) else float("nan"),
                              p[6], info, keys, samps))
    return Vcf(header, samples, recs)


def write_vcf(vcf: Vcf) -> str:
    return vcf.as_text()


def vcf_or_text(src: Any) -> Vcf:
    if isinstance(src, Vcf):
        return src
    return parse_vcf(src)


# ---------------------------------------------------------------------------
# SAM / light BAM surrogate
# ---------------------------------------------------------------------------
@dataclass
class Aln:
    qname: str
    flag: int
    rname: str
    pos: int
    mapq: int
    cigar: str
    rnext: str = "="
    pnext: int = 0
    tlen: int = 0
    seq: str = ""
    qual: str = ""
    tags: dict[str, str] = field(default_factory=dict)

    @property
    def reverse(self) -> bool:
        return bool(self.flag & 0x10)

    @property
    def secondary(self) -> bool:
        return bool(self.flag & 0x100)

    @property
    def supplementary(self) -> bool:
        return bool(self.flag & 0x800)

    @property
    def duplicate(self) -> bool:
        return bool(self.flag & 0x400)

    @property
    def unmapped(self) -> bool:
        return bool(self.flag & 0x4)

    @property
    def mapped(self) -> bool:
        return not self.unmapped

    @property
    def read1(self) -> bool:
        return bool(self.flag & 0x40)

    @property
    def read2(self) -> bool:
        return bool(self.flag & 0x80)

    @property
    def proper_pair(self) -> bool:
        return bool(self.flag & 0x2)

    @property
    def ref_end(self) -> int:
        return self.pos + ref_length(self.cigar)

    @property
    def query_length(self) -> int:
        return query_length(self.cigar)

    def as_line(self) -> str:
        t = [f"{k}:{v}" for k, v in self.tags.items()]
        return "\t".join([self.qname, str(self.flag), self.rname, str(self.pos),
                          str(self.mapq), self.cigar, self.rnext, str(self.pnext),
                          str(self.tlen), self.seq, self.qual or "*"] + t)

    def __getitem__(self, key):
        """Dict-style access so core helpers accept alignments either way."""
        if key in self.__dict__ or hasattr(type(self), key):
            return getattr(self, key)
        if key in self.tags:
            return self.tags[key]
        if key in ("read", "id"):
            return self.qname
        raise KeyError(key)

    def get(self, key, default=None):
        try:
            return self[key]
        except (KeyError, AttributeError):
            return default



CIGAR_RE = re.compile(r"(\d+)([MIDNSHP=X])")


def cigar_ops(cigar: str) -> list[tuple[int, str]]:
    return [(int(n), op) for n, op in CIGAR_RE.findall(cigar or "")]


def ref_length(cigar: str) -> int:
    return sum(n for n, op in cigar_ops(cigar) if op in "MDN=X")


def query_length(cigar: str) -> int:
    return sum(n for n, op in cigar_ops(cigar) if op in "MIS=X")


def parse_sam(text: str) -> tuple[list[str], list[Aln]]:
    header, alns = [], []
    for ln in as_text(text).splitlines():
        if not ln.strip():
            continue
        if ln.startswith("@"):
            header.append(ln.rstrip())
            continue
        p = ln.rstrip("\n").split("\t")
        if len(p) < 11:
            continue
        tags = {}
        for t in p[11:]:
            sp = t.split(":")
            if len(sp) >= 3:
                tags[sp[0]] = ":".join(sp[2:])
        try:
            alns.append(Aln(p[0], int(p[1]), p[2], int(p[3]), int(p[4]), p[5],
                            p[6], int(p[7]), int(p[8]), p[9],
                            "" if p[10] == "*" else p[10], tags))
        except ValueError:
            continue
    return header, alns


def write_sam(header: list[str], alns: Iterable[Aln]) -> str:
    hdr = list(header) or []
    if not any(h.startswith("@SQ") for h in hdr):
        pass
    if not any(h.startswith("@HD") for h in hdr):
        hdr.insert(0, "@HD\tVN:1.6\tSO:coordinate")
    return "\n".join(hdr + [a.as_line() for a in alns]) + "\n"


def sam_ref_names(header: list[str]) -> list[tuple[str, int]]:
    out = []
    for h in header:
        if h.startswith("@SQ"):
            d = dict(kv.split(":", 1) for kv in h.split("\t")[1:] if ":" in kv)
            out.append((d.get("SN", ""), int(d.get("LN", 0) or 0)))
    return out


def sam_header(refs: Iterable[tuple[str, int]], rg: list[dict] | None = None,
               pg: list[dict] | None = None, sort: str = "coordinate") -> list[str]:
    hdr = [f"@HD\tVN:1.6\tSO:{sort}"]
    for name, ln in refs:
        hdr.append(f"@SQ\tSN:{name}\tLN:{ln}")
    for r in rg or []:
        hdr.append("@RG\t" + "\t".join(f"{k}:{v}" for k, v in r.items()))
    for p in pg or []:
        hdr.append("@PG\t" + "\t".join(f"{k}:{v}" for k, v in p.items()))
    return hdr


# ---------------------------------------------------------------------------
# misc formats
# ---------------------------------------------------------------------------
def parse_matrix(text: str, sep: str | None = None) -> tuple[list[str], list[str], list[list[float]]]:
    """Row/column labelled matrix (e.g. expression counts) -> (rows, cols, data)."""
    body = [ln for ln in as_text(text).splitlines() if ln.strip()]
    if not body:
        return [], [], []
    if sep is None:
        sep = "\t" if "\t" in body[0] else ("," if "," in body[0] else None)
    split = (lambda s: s.split(sep)) if sep else (lambda s: s.split())
    cols = split(body[0].lstrip("#").strip())[1:]
    rnames, data = [], []
    for ln in body[1:]:
        p = split(ln)
        if not p or p[0].startswith("#"):
            continue
        rnames.append(p[0])
        row = []
        for x in p[1:]:
            try:
                row.append(float(x))
            except ValueError:
                row.append(float("nan"))
        while len(row) < len(cols):
            row.append(float("nan"))
        data.append(row[:len(cols)])
    return rnames, cols, data


def write_matrix(rnames: list[str], cols: list[str], data: list[list[float]],
                 row_label: str = "row", sep: str = "\t") -> str:
    lines_ = [sep.join([row_label] + list(cols))]
    for r, row in zip(rnames, data):
        lines_.append(sep.join([r] + [_fmt_num(v) for v in row]))
    return "\n".join(lines_) + "\n"


def parse_gmt(text: str) -> dict[str, list[str]]:
    sets: dict[str, list[str]] = {}
    for ln in as_text(text).splitlines():
        if not ln.strip():
            continue
        p = ln.rstrip("\n").split("\t")
        if len(p) < 2:
            p = ln.split()
        if not p:
            continue
        sets[p[0]] = [x for x in p[2:] if x] or [x for x in p[1:] if x and not x.startswith(("MIR", "WP"))]
    return sets


def write_gmt(sets: dict[str, list[str]], descriptions: dict[str, str] | None = None) -> str:
    out = []
    for name, members in sets.items():
        desc = (descriptions or {}).get(name, "")
        out.append("\t".join([name, desc] + list(members)))
    return "\n".join(out) + "\n"


def parse_pfm(text: str) -> dict[str, list[float]]:
    """JASPAR/MEME-ish PFM: four rows of counts (or one row per position)."""
    rows: list[list[float]] = []
    order = "ACGT"
    for ln in as_text(text).splitlines():
        ln = ln.strip()
        if not ln or ln[0] in ">#":
            continue
        p = ln.split()
        if len(p) >= 2 and p[0].upper() in order:
            vals = [float(x) for x in p[1:]]
            rows.append(vals)
            order = order  # rows are A,C,G,T
            continue
        if len(p) in (4, 5) and all(_isnum(x) for x in p[:4]):
            rows.append([float(x) for x in p[:4]])
    if not rows:
        return {}
    if len(rows) == 4:  # transpose so rows = positions
        n = len(rows[0])
        return {order[i]: [rows[i][j] for j in range(n)] for i in range(4)}
    return {order[i]: [r[i] for r in rows] for i in range(4)}


def parse_pairs(text: str) -> list[dict]:
    """Hi-C ``pairs`` format."""
    out = []
    names: list[str] = []
    cols: list[str] = []
    for ln in as_text(text).splitlines():
        if ln.startswith("#"):
            if "columns=" in ln:
                names = ln.split("columns=", 1)[1].split()
            continue
        if not ln.strip():
            continue
        p = ln.split()
        cols = cols or names
        ncol = len(cols) or 7
        if len(p) > ncol and len(p) == ncol + 1 and p[1].isdigit():
            p = [p[0]] + p[2:]  # pairs 2.1 with an extra read-count column
        rec = dict(zip(cols or ["id", "chrom1", "pos1", "chrom2", "pos2", "strand1",
                                "strand2"], p[:ncol]))
        out.append(rec)
    return out


# ---------------------------------------------------------------------------
# trees
def parse_newick(text: str) -> dict:
    """Newick -> nested dict ``{name, length, children}`` (pure python, no Biopython)."""
    s = as_text(text)
    s = ";".join(ln.strip() for ln in s.splitlines() if ln.strip() and not ln.startswith("#")).rstrip(";")
    tree = _newick_parse(s)
    if not tree["name"] and len(tree["children"]) == 1:
        return tree["children"][0]
    return tree


def _newick_parse(s: str) -> dict:
    """Recursive-descent Newick reader (labels, branch lengths, bootstrap values)."""
    s = re.sub(r"\[[^\]]*\]", "", s)          # drop [and,,] annotations
    node, _ = _newick_subtree(s, 0)
    return node


def _newick_subtree(s: str, i: int) -> tuple[dict, int]:
    node: dict[str, Any] = {"name": "", "length": None, "children": []}
    if i < len(s) and s[i] == "(":
        i += 1
        while i < len(s):
            child, i = _newick_subtree(s, i)
            node["children"].append(child)
            while i < len(s) and s[i].isspace():
                i += 1
            if i < len(s) and s[i] == ",":
                i += 1
                continue
            if i < len(s) and s[i] == ")":
                i += 1
            break
    name = ""
    while i < len(s) and s[i] not in ",();:":
        name += s[i]
        i += 1
    node["name"] = name.strip().strip("'")
    if i < len(s) and s[i] == ":":
        i += 1
        num = ""
        while i < len(s) and s[i] not in ",();":
            num += s[i]
            i += 1
        try:
            node["length"] = float(num.strip())
        except ValueError:
            pass
    return node, i


def write_newick(node: dict) -> str:
    def rec(n: dict) -> str:
        s = ""
        if n.get("children"):
            s = "(" + ",".join(rec(c) for c in n["children"]) + ")"
        s += (n.get("name") or "")
        if n.get("length") is not None:
            s += f":{n['length']:g}"
        return s

    return rec(node) + ";"


def newick_leaves(node: dict) -> list[str]:
    if not node.get("children"):
        return [node.get("name") or ""]
    out_: list[str] = []
    for c in node["children"]:
        out_.extend(newick_leaves(c))
    return out_


def out(message: str = "", **extra: Any) -> dict:
    """Standard tool result payload: a message plus optional table/image/files."""
    payload: dict[str, Any] = {}
    if message:
        payload["message"] = message
    payload.update({k: v for k, v in extra.items() if v is not None})
    return payload


def table_result(df: Any, message: str = "", name: str = "table.tsv") -> dict:
    import pandas as pd

    if isinstance(df, pd.Series):
        df = df.to_frame()
    if not isinstance(df, pd.DataFrame):
        df = pd.DataFrame(df)
    res: dict[str, Any] = {"table": df, "filename": name}
    if message:
        res["message"] = message
    res["rows"] = int(len(df))
    res["text"] = to_tsv(df)
    res.setdefault("stats", {})
    return res


def text_result(text: str, message: str = "") -> dict:
    res = {"text": text}
    if message:
        res["message"] = message
    return res


def image_result(fig: Any, message: str = "", name: str = "figure.png") -> dict:
    data = figure_bytes(fig)
    close(fig)
    res: dict[str, Any] = {"image": data, "filename": name}
    if message:
        res["message"] = message
    return res


def figure_bytes(fig: Any) -> bytes:
    buf = _io.BytesIO()
    fig.savefig(buf, format="png", dpi=110, bbox_inches="tight")
    return buf.getvalue()


def close(fig: Any) -> None:
    try:
        import matplotlib.pyplot as plt

        plt.close(fig)
    except Exception:  # pragma: no cover
        pass


def save_download(text: str, suffix: str = "txt") -> str:
    p = Path(tempfile.mkdtemp(prefix="chroma_out_")) / f"result.{suffix}"
    p.write_text(text, encoding="utf-8")
    return str(p)


def read_any(src: Any) -> Any:
    """Best-effort parse of a dataset into a python object (records or table)."""
    txt = as_text(src)
    first = next((l for l in txt.splitlines() if l.strip()), "")
    if first.startswith(">"):
        return parse_fasta(txt)
    if first.startswith("@") and not first.startswith("@HD") and not first.startswith("@r"):
        return parse_fastq(txt)
    if first.startswith("#CHROM"):
        return parse_vcf(txt)
    if first.startswith(("@SQ", "@HD")):
        return parse_sam(txt)[1]
    if first.startswith("##gff"):
        return parse_gff(txt)
    if "\t" in first and len(first.split("\t")) >= 3 and _isnum(first.split("\t")[1]):
        return parse_bed(txt)
    return read_table(txt)


def _isnum(x: Any) -> bool:
    try:
        float(x)
        return True
    except (TypeError, ValueError):
        return False


def write_rows(rows: Iterable[Iterable[Any]], sep: str = "\t") -> str:
    return "\n".join(sep.join(str(x) for x in r) for r in rows) + "\n"
