"""FASTA / FASTQ sequence tools (Galaxy 'FASTA/FASTQ', seqtk, format conversion)."""

from __future__ import annotations

import math
import random
import re
from collections import Counter, OrderedDict, defaultdict

from chroma_titan.core import io, seq, stats, tables
from chroma_titan.tools._common import *  # noqa: F401,F403
from chroma_titan.tools._common import GENES_FA, READS_SINGLE, ALN_SAM

FASTA_SECTIONS = "fasta_fastq"


def recs(src):
    return io.parse_fasta(io.as_text(src))


def _write(rs, width=60, fmt="fasta"):
    return {"text": io.write_fasta(rs, width=width, format=fmt),
            "filename": f"output.{fmt}", "message": f"{len(rs)} sequences"}


def _ids(rs):
    return [r.id for r in rs]


# ===========================================================================
# basic FASTA manipulation
# ===========================================================================
@T("fasta_filter_by_length", "Filter FASTA sequences by length", FASTA_SECTIONS, "fasta",
   [fa("src"), intin("min_len", 100), intin("max_len", 0), boolean("complement", False)],
   ex={"src": GENES_FA, "min_len": 150, "max_len": 0}, up="iuc/seqtk",
   summary="Keep sequences whose length falls inside a window.")
def fasta_filter_length(src, min_len=100, max_len=0, complement=False):
    """Sequence-length filter (``Fasta Filter``)."""
    out = []
    for r in recs(src):
        n = len(r.seq)
        if n >= int(min_len) and (not max_len or n <= int(max_len)):
            out.append(io.Seq(r.id, seq.reverse_complement(r.seq) if complement else r.seq, r.desc))
    res = _write(out)
    res["message"] = f"{len(out)}/{len(recs(src))} sequences kept"
    return res


@T("fasta_head", "First N sequences", FASTA_SECTIONS, "fasta", [fa("src"), intin("n", 5)],
   ex={"src": GENES_FA, "n": 3}, up="fasta_first_sequences", summary="Extract the first n records.")
def fasta_head(src, n=5):
    """``head`` for FASTA records."""
    return _write(recs(src)[: int(n)])


@T("fasta_tail", "Last N sequences", FASTA_SECTIONS, "fasta", [fa("src"), intin("n", 5)],
   ex={"src": GENES_FA, "n": 2}, up="fasta_last_sequences", summary="Extract the last n records.")
def fasta_tail(src, n=5):
    """``tail`` for FASTA records."""
    return _write(recs(src)[-int(n):])


@T("fasta_get_ids", "Extract sequence identifiers", FASTA_SECTIONS, "text",
   [fa("src"), choice("field", ["id", "description", "id_and_description", "index"]),
    boolean("unique", False)], ex={"src": GENES_FA}, up="fasta_to_tabulate",
   summary="List the names of all records in a FASTA file.")
def fasta_identifiers(src, field="id", unique=False):
    """Dump the defline fields of every record."""
    rs = recs(src)
    vals = [r.id for r in rs] if field == "id" else \
        [r.desc for r in rs] if field == "description" else \
        [f"{i}" for i in range(len(rs))] if field == "index" else \
        [f"{r.id} {r.desc}".strip() for r in rs]
    if unique:
        vals = list(OrderedDict.fromkeys(vals))
    return text("\n".join(vals) + "\n", f"{len(vals)} identifiers")


@T("fasta_extract_by_ids", "Fetch sequences by identifiers", FASTA_SECTIONS, "fasta",
   [fa("src"), bigtext("identifiers", "", help="one id per line"), boolean("invert", False),
    boolean("allow_prefix", True)], ex={"src": GENES_FA, "identifiers": "gene_1\ngene_3"},
   up="iuc/seqkit", summary="Extract (or drop) records whose id matches a list.")
def fasta_extract_by_ids(src, identifiers="", invert=False, allow_prefix=True):
    """``seqtk subseq``-like selection by name."""
    wanted = [l.strip() for l in _lines_like(identifiers)]
    wset = set(wanted)
    rs = recs(src)
    def hit(r):
        return r.id in wset or any(r.id.startswith(w) for w in wanted if allow_prefix)
    out = [r for r in rs if hit(r) != invert]
    res = _write(out)
    res["message"] = f"{len(out)}/{len(rs)} sequences fetched"
    return res


@T("fasta_sort", "Sort FASTA records", FASTA_SECTIONS, "fasta",
   [fa("src"), choice("by", ["id", "length", "gc", "alphabetical"]), boolean("reverse", False)],
   ex={"src": GENES_FA, "by": "length", "reverse": True}, up="fasta_order",
   summary="Order records by name, length or GC content.")
def fasta_sort_tool(src, by="id", reverse=False):
    """Sort the records of a FASTA file."""
    rs = recs(src)
    key = {"id": lambda r: r.id, "alphabetical": lambda r: r.id.lower(),
           "length": lambda r: len(r.seq), "gc": lambda r: seq.gc_content(r.seq)}[by]
    return _write(sorted(rs, key=key, reverse=reverse))


@T("fasta_gc_content", "GC content per sequence", FASTA_SECTIONS, "table",
   [fa("src"), intin("precision", 4)], ex={"src": GENES_FA}, up="seqtk",
   summary="Length, GC%, AT-skew and GC-skew of every record.")
def fasta_gc_content(src, precision=4):
    """Per-sequence GC content (EMBOSS ``newcpgreport``-like)."""
    out = [{"id": r.id, "length": len(r.seq), "gc_percent": round(100 * seq.gc_content(r.seq), precision),
            "gc_skew": round(seq.gc_skew(r.seq), precision), "at_skew": round(seq.at_skew(r.seq), precision),
            "a": r.seq.upper().count("A"), "c": r.seq.upper().count("C"),
            "g": r.seq.upper().count("G"), "t": r.seq.upper().count("T"),
            "n": r.seq.upper().count("N")} for r in recs(src)]
    return table(out, f"{len(out)} sequences summarised")


@T("fasta_count_nucleotides", "Count DNA nucleotides", FASTA_SECTIONS, "table",
   [fa("src"), boolean("per_sequence", True)], ex={"src": GENES_FA}, up="count_gatb",
   summary="A/C/G/T/N counts and percentages (Galaxy 'Count DNA nucleotides').")
def count_nucleotides(src, per_sequence=True):
    """Nucleotide composition, per record or for the whole file."""
    rs = recs(src)
    if per_sequence:
        out = []
        for r in rs:
            c = seq.counts(r.seq)
            tot = sum(c.values()) or 1
            out.append({"id": r.id, "A": c["A"], "C": c["C"], "G": c["G"], "T": c["T"], "N": c["N"],
                        "length": tot, "GC_percent": round(100 * (c["G"] + c["C"]) / tot, 4)})
        return table(out, "nucleotide counts per sequence")
    c = Counter()
    for r in rs:
        c.update(seq.counts(r.seq))
    tot = sum(c.values()) or 1
    out = [{"base": b, "count": c[b], "percent": round(100 * c[b] / tot, 4)} for b in "ACG TN".replace(" ", "")]
    return table(out, f"nucleotide counts over {len(rs)} sequences")


@T("fasta_nucleotide_frequencies", "Nucleotide frequencies", FASTA_SECTIONS, "table",
   [fa("src")], ex={"src": GENES_FA}, up="fasta-rc_cycle",
   summary="Frequency of every base in the concatenated alignment.")
def nuc_frequencies(src):
    """Base frequencies (Galaxy 'Compute Nucleotide Frequency')."""
    allseq = "".join(r.seq for r in recs(src))
    fr = seq.nucleotide_frequencies(allseq)
    return table([{"base": k, "frequency": round(v, 6), "percent": round(100 * v, 4)}
                  for k, v in sorted(fr.items())], f"{len(allseq)} bases")


@T("fasta_complement", "Reverse complement", FASTA_SECTIONS, "fasta",
   [fa("src"), choice("mode", ["revcomp", "complement", "reverse"])],
   ex={"src": GENES_FA, "mode": "revcomp"}, up="fastarevseq",
   summary="Return the reverse complement (or one half of it).")
def fasta_revcomp(src, mode="revcomp"):
    """Reverse-complement every record."""
    fn = {"revcomp": seq.reverse_complement, "complement": seq.complement, "reverse": lambda s: s[::-1]}[mode]
    return _write([io.Seq(r.id, fn(r.seq), r.desc) for r in recs(src)])


@T("fasta_uppercase", "Convert sequence case", FASTA_SECTIONS, "fasta",
   [fa("src"), choice("case", ["upper", "lower"])], ex={"src": GENES_FA, "case": "upper"},
   up="Change_case", summary="Upper/lower-case the residues (soft-masking convention).")
def fasta_case(src, case="upper"):
    """Case conversion of sequence characters."""
    fn = str.upper if case == "upper" else str.lower
    return _write([io.Seq(r.id, fn(r.seq), r.desc) for r in recs(src)])


@T("fasta_replace_ids", "Add or replace sequence descriptions", FASTA_SECTIONS, "fasta",
   [fa("src"), textbox("prefix", ""), textbox("suffix", ""), boolean("strip_description", True)],
   ex={"src": GENES_FA, "prefix": "chrV_"}, up="fasta_str_tie",
   summary="Rename records with a prefix/suffix and optional defline cleanup.")
def fasta_rename(src, prefix="", suffix="", strip_description=True):
    """Batch renaming of FASTA headers."""
    out = []
    for i, r in enumerate(recs(src), 1):
        nid = f"{prefix}{r.id}{suffix}"
        out.append(io.Seq(nid, r.seq, "" if strip_description else r.desc))
    return _write(out)


@T("fasta_wrap", "Change line wrapping", FASTA_SECTIONS, "fasta",
   [fa("src"), intin("width", 70), boolean("single_line", False)],
   ex={"src": GENES_FA, "width": 60}, up="fold", summary="Re-wrap sequence lines to a given width.")
def fasta_wrap(src, width=70, single_line=False):
    """Rewrap or unwrap FASTA lines."""
    return _write(recs(src), width=10 ** 6 if single_line else int(width))


@T("fasta_chop", "Fasta Chopper (split into fixed-size chunks)", FASTA_SECTIONS, "fasta",
   [fa("src"), intin("chunk_size", 200), boolean("discard_partial", False)],
   ex={"src": GENES_FA, "chunk_size": 150}, up="fasta_split_on_cases",
   summary="Split each sequence into equal-size pieces with numbered ids.")
def fasta_chop(src, chunk_size=200, discard_partial=False):
    """Galaxy 'Fasta Chopper'."""
    out = []
    for r in recs(src):
        n = int(chunk_size)
        parts = [r.seq[i:i + n] for i in range(0, len(r.seq), n)]
        if discard_partial and parts and len(parts[-1]) < n:
            parts = parts[:-1]
        out.extend(io.Seq(f"{r.id}_{k + 1}", p, r.desc) for k, p in enumerate(parts))
    return _write(out)


@T("fasta_split_on_cases", "Split FASTA on a column of a table", FASTA_SECTIONS, "table",
   [fa("src"), intin("n_parts", 2)], ex={"src": GENES_FA, "n_parts": 2}, up="fasta_split",
   summary="Distribute records into N roughly equal partitions.")
def fasta_split_parts(src, n_parts=2):
    """Round-robin split of records into n groups (reported as a table)."""
    rs = recs(src)
    n = max(1, int(n_parts))
    buckets = [[] for _ in range(n)]
    for i, r in enumerate(rs):
        buckets[i % n].append(r.id)
    return table([{"partition": i + 1, "n_sequences": len(b), "identifiers": "|".join(b)[:200]}
                  for i, b in enumerate(buckets)], f"split {len(rs)} sequences into {n} parts")


@T("fasta_split_by_size", "Split FASTA into files of a fixed size", FASTA_SECTIONS, "text",
   [fa("src"), intin("sequences_per_file", 3)], ex={"src": GENES_FA, "sequences_per_file": 3},
   up="fastasplit", summary="Chunk a multi-FASTA into groups of k sequences.")
def fasta_split_by_size(src, sequences_per_file=3):
    """Galaxy 'FASTASplit': report each chunk (as FASTA text per file)."""
    rs = recs(src)
    k = max(1, int(sequences_per_file))
    chunks = [rs[i:i + k] for i in range(0, len(rs), k)]
    body = "\n\n".join(f"# chunk {i + 1}\n{io.write_fasta(c)}" for i, c in enumerate(chunks))
    return text(body, f"{len(chunks)} chunks of up to {k} sequences")


@T("fasta_get_length", "Get lengths", FASTA_SECTIONS, "table",
   [fa("src"), boolean("cumulative", False)], ex={"src": GENES_FA}, up="fasta-get-metadata",
   summary="Length of every sequence, plus file totals.")
def fasta_lengths(src, cumulative=False):
    """Sequence lengths and cumulative sums."""
    rs = recs(src)
    out, tot = [], 0
    for r in rs:
        tot += len(r.seq)
        out.append({"id": r.id, "length": len(r.seq), "cumulative": tot,
                    "description": r.desc[:60]})
    return table(out, f"{len(rs)} sequences, {tot} bases total")


@T("fasta_statistics", "Sequence statistics of a FASTA file", FASTA_SECTIONS, "stats",
   [fa("src")], ex={"src": GENES_FA}, up="seqtk comp",
   summary="min/max/mean/median length, N50, GC and alphabet checks.")
def fasta_statistics(src):
    """Summary of a nucleotide or protein FASTA (like ``seqtk comp``/QC)."""
    rs = recs(src)
    lens = [len(r.seq) for r in rs]
    allseq = "".join(r.seq for r in rs)
    st = {"n_sequences": len(rs), "total_bases": sum(lens), "longest": max(lens) if lens else 0,
          "shortest": min(lens) if lens else 0, "mean_length": round(stats.mean(lens), 2) if lens else 0,
          "median_length": stats.median(lens) if lens else 0,
          "stdev_length": round(stats.stdev(lens), 2) if len(lens) > 1 else 0,
          "n50": seq.n50(lens), "l50": seq.l50(lens), "gc_percent": round(100 * seq.gc_content(allseq), 3),
          "alphabet": seq.alphabet_of(allseq), "n_bases": allseq.upper().count("N"),
          "ambiguous_fraction": round(sum(1 for c in allseq if c.upper() in "RYKMSWBDHVNX") / max(1, len(allseq)), 5)}
    return values_message(st, f"{len(rs)} sequences, {sum(lens)} residues")


@T("fasta_to_table", "FASTA to tabular", FASTA_SECTIONS, "table",
   [fa("src"), multi("fields", ["id", "description", "length", "sequence", "gc", "sha1", "first_base",
                               "last_base"], ["id", "length", "gc"]), intin("truncate_sequence", 0)],
   ex={"src": GENES_FA, "fields": ["id", "length", "gc", "sha1"]}, up="fasta_to_tabulate",
   summary="Convert records into a tab-delimited summary table.")
def fasta_to_table(src, fields=None, truncate_sequence=0):
    """Tabulate selected fields of a FASTA file."""
    fields = list(fields or ["id", "length"])
    out = []
    for r in recs(src):
        s = r.seq
        rec = {}
        for f in fields:
            if f == "id":
                rec["id"] = r.id
            elif f == "description":
                rec["description"] = r.desc
            elif f == "length":
                rec["length"] = len(s)
            elif f == "sequence":
                rec["sequence"] = s[: int(truncate_sequence)] if truncate_sequence else s
            elif f == "gc":
                rec["gc_percent"] = round(100 * seq.gc_content(s), 3)
            elif f == "sha1":
                rec["sha1"] = seq.sha1_of(s.upper())[:16]
            elif f == "first_base":
                rec["first_base"] = s[:1]
            elif f == "last_base":
                rec["last_base"] = s[-1:]
        out.append(rec)
    return table(out, f"{len(out)} rows")


@T("table_to_fasta", "Tabular to FASTA", "convert_formats", "fasta",
   [tbl("src"), textbox("id_column", "1"), textbox("sequence_column", "4"),
    choice("format", ["fasta", "fastq"]), textbox("quality", "I")],
   ex={"src": REGIONS_BED, "id_column": "4", "sequence_column": "3"}, up="table_to_fasta",
   summary="Build a FASTA (or FASTQ) file from two columns of a table.")
def table_to_fasta(src, id_column="1", sequence_column="4", format="fasta", quality="I"):
    """Two-column table → sequence file."""
    rows = [l.split("\t") if "\t" in l else l.split() for l in io.lines(src)]
    if not rows:
        return {"text": "", "message": "empty table"}
    ncol = max(len(r) for r in rows)
    hi = _ci(rows[0], ncol) if _looks_numeric_row(rows[0]) else None
    ii = int(id_column) - 1 if str(id_column).isdigit() else (rows[0].index(id_column) if hi is None else 0)
    si = int(sequence_column) - 1 if str(sequence_column).isdigit() else (
        rows[0].index(sequence_column) if hi is None else min(3, ncol - 1))
    body = rows[1:] if hi is not None else rows
    out = []
    for r in body:
        if max(ii, si) >= len(r):
            continue
        s = re.sub(r"[^A-Za-z]", "", r[si]).upper() or "NNNN"
        out.append(io.Seq(r[ii], s))
    if format == "fastq":
        return {"text": io.write_fastq([io.Read(x.id, x.seq, quality * len(x.seq)) for x in out]),
                "filename": "output.fastq", "message": f"{len(out)} reads"}
    return _write(out)


def _looks_numeric_row(row):
    return all(re.match(r"^-?\d+(\.\d+)?$", c.strip()) for c in row if c.strip())


@T("fasta_remove_gaps", "Remove gaps, dashes and unknown characters", FASTA_SECTIONS, "fasta",
   [fa("src"), multi("chars", ["-", ".", " ", "N", "*", "?", "X"], ["-", "."]),
    boolean("keep_case", True)], ex={"src": GENES_FA, "chars": ["-", ".", "N"]}, up="ungap",
   summary="Strip gap and filler characters from every sequence.")
def remove_gaps(src, chars=None, keep_case=True):
    """Delete the selected filler characters from all sequences."""
    drop = set("".join(chars or ["-", "."]))
    out = []
    for r in recs(src):
        s = "".join(c for c in r.seq if c not in drop and c.upper() not in drop)
        out.append(io.Seq(r.id, s if keep_case else s.upper(), r.desc))
    res = _write(out)
    res["message"] = f"removed {sum(len(r.seq) for r in recs(src)) - sum(len(x.seq) for x in out)} characters"
    return res


@T("fasta_clean_n", "Clean Ns from sequence ends", FASTA_SECTIONS, "fasta",
   [fa("src"), boolean("trim_internal", False), intin("max_internal_n", 5)],
   ex={"src": GENES_FA, "trim_internal": False}, up="clean_n",
   summary="Trim leading/trailing runs of N (optionally mask internal runs).")
def clean_ns(src, trim_internal=False, max_internal_n=5):
    """Galaxy 'Clean Ns' tool."""
    out = []
    for r in recs(src):
        s = r.seq.upper()
        if trim_internal:
            s = re.sub("N{%d,}" % max(1, int(max_internal_n)), "", s)
        s = re.sub(r"^N+", "", s)
        s = re.sub(r"N+$", "", s)
        if s:
            out.append(io.Seq(r.id, s, r.desc))
    res = _write(out)
    res["message"] = f"{len(out)} cleaned sequences"
    return res


@T("fasta_trim", "Trim sequences (fixed or by quality)", FASTA_SECTIONS, "fasta",
   [fa("src"), intin("trim_5p", 0), intin("trim_3p", 0), intin("keep_from", 0), intin("keep_to", 0)],
   ex={"src": GENES_FA, "trim_5p": 3, "trim_3p": 5}, up="trim_sequences",
   summary="Remove bases from the ends or keep a fixed window.")
def fasta_trim(src, trim_5p=0, trim_3p=0, keep_from=0, keep_to=0):
    """Trim N bases off each end, or extract an interval."""
    out = []
    for r in recs(src):
        s = r.seq
        if keep_from or keep_to:
            s = s[int(keep_from):int(keep_to) or len(s)]
        else:
            s = s[int(trim_5p): len(s) - int(trim_3p) if trim_3p else len(s)]
        if s:
            out.append(io.Seq(r.id, s, r.desc))
    return _write(out)


@T("fasta_reverse_sequences", "Reverse sequence order and orientation", FASTA_SECTIONS, "fasta",
   [fa("src"), boolean("reverse_records", True), boolean("reverse_residues", False)],
   ex={"src": GENES_FA}, up="fasta_rc", summary="Reverse the record order and/or the residues.")
def fasta_reverse(src, reverse_records=True, reverse_residues=False):
    """Flip record order (like Galaxy 'Reverse Complement' variants)."""
    rs = recs(src)
    if reverse_residues:
        rs = [io.Seq(r.id, r.seq[::-1], r.desc) for r in rs]
    if reverse_records:
        rs = rs[::-1]
    return _write(rs)


@T("fasta_dedupe", "Remove duplicate sequences", FASTA_SECTIONS, "fasta",
   [fa("src"), choice("key", ["sequence", "id", "sequence_ignoring_case", "canonical"]),
    boolean("report", True)], ex={"src": GENES_FA}, up="rdist",
   summary="Collapse records that share the same sequence (or id).")
def fasta_dedupe(src, key="sequence", report=True):
    """Deduplicate a FASTA file."""
    def k(r):
        s = r.seq.upper()
        if key == "id":
            return r.id
        if key == "sequence_ignoring_case":
            return s
        if key == "canonical":
            return min(s, seq.reverse_complement(s))
        return r.seq
    seen, out = set(), []
    for r in recs(src):
        kk = k(r)
        if kk in seen:
            continue
        seen.add(kk)
        out.append(r)
    res = _write(out)
    res["message"] = f"{len(out)} unique sequences (removed {len(recs(src)) - len(out)})"
    return res


@T("fasta_dna_to_protein", "DNA to protein", FASTA_SECTIONS, "fasta",
   [seqin("src", GENES_FA, label="Nucleotide FASTA"), intin("frame", 0),
    choice("table", ["Standard", "Vertebrate Mitochondrial", "Yeast Mitochondrial", "Bacterial"]),
    boolean("to_stop", True)], ex={"src": GENES_FA, "frame": 0, "to_stop": True}, up="mgrna2nuc",
   summary="Translate a nucleotide FASTA in one reading frame.")
def dna_to_protein(src, frame=0, table="Standard", to_stop=True):
    """Translation of every record (EMBOSS ``dna_to_pro`` style)."""
    out = []
    for r in recs(src):
        p = seq.translate(r.seq, frame=int(frame), to_stop=to_stop)
        if p:
            out.append(io.Seq(r.id, p, f"translated frame {int(frame) + 1}"))
    return _write(out)


@T("fasta_translate", "Translate (all six frames)", FASTA_SECTIONS, "table",
   [seqin("src", GENES_FA, label="Nucleotide FASTA"), choice("mode", ["best_orf", "all_frames"]),
    boolean("to_stop", True)], ex={"src": GENES_FA, "mode": "best_orf"}, up="translate_tool",
   summary="Report the translation(s) and longest ORF per sequence.")
def fasta_translate(src, mode="best_orf", to_stop=True):
    """Six-frame translation table with protein lengths."""
    out = []
    for r in recs(src):
        frames = seq.translate_all_frames(r.seq, to_stop=to_stop)
        if not frames:
            continue
        frames = [{**f, "length_nt": len(f["protein"]) * 3,
                   "strand": "+" if f["frame"] > 0 else "-"} for f in frames]
        best = max(frames, key=lambda d: len(d["protein"]))
        if mode == "best_orf":
            out.append({"id": r.id, "frame": best["frame"], "strand": best["strand"],
                        "length_nt": best["length_nt"], "protein_length": len(best["protein"]),
                        "protein": best["protein"][:120]})
        else:
            for f in frames:
                out.append({"id": r.id, "frame": f["frame"], "strand": f["strand"],
                            "length_nt": f["length_nt"], "protein_length": len(f["protein"]),
                            "protein": f["protein"][:120]})
    return table(out, f"translations for {len(recs(src))} sequences")


@T("fasta_get_orfs", "Get ORFs", FASTA_SECTIONS, "fasta",
   [seqin("src", GENES_FA, label="Nucleotide FASTA"), intin("min_length", 90),
    textbox("start_codons", "ATG"), boolean("both_strands", True), boolean("table", False)],
   ex={"src": GENES_FA, "min_length": 90}, up="getorf",
   summary="Find open reading frames (EMBOSS ``getorf``).")
def get_orfs(src, min_length=90, start_codons="ATG", both_strands=True, table=False):
    """Report ORFs as FASTA records or as a table."""
    out = []
    for r in recs(src):
        orfs = seq.find_orfs(r.seq, min_len=int(min_length), required_start=start_codons or "ATG",
                             partial=False)
        for o in orfs:
            out.append({"id": f"{r.id}_{o['strand']}_{o['start'] + 1}", "parent": r.id,
                        "strand": o["strand"], "frame": o["frame"], "start": o["start"],
                        "end": o["end"], "length": o["length"], "aa_length": o["aa_length"],
                        "start_codon": o.get("start_codon", ""), "stop_codon": o.get("stop_codon", ""),
                        "protein": o["protein"]})
    if table:
        return table([{k: v for k, v in o.items() if k != "protein"} for o in out],
                     f"{len(out)} ORFs")
    return _write([io.Seq(o["id"], o["protein"],
                          f"{o['strand']} {o['start'] + 1}-{o['end']} nt={o['length']}") for o in out])


@T("fasta_longest_orf", "Longest ORF as protein", FASTA_SECTIONS, "text",
   [seqin("src", GENES_FA, label="Nucleotide FASTA")], ex={"src": GENES_FA}, up="longest_orf",
   summary="Print the longest translation found in the file.")
def longest_orf(src):
    """Single longest ORF protein of a nucleotide dataset."""
    p = seq.longest_orf_protein(src)
    return text((p or "") + "\n", f"longest ORF protein: {len(p or '')} aa")


@T("fasta_get_hsp", "Consensus sequence of a multiple alignment (FASTA)", FASTA_SECTIONS, "text",
   [fa("src"), choice("mode", ["majority", "first", "ambiguous"])], ex={"src": MSA, "mode": "majority"},
   up="new_consensus", summary="Collapse aligned sequences into a consensus string.")
def consensus_from_alignment(src, mode="majority"):
    """Galaxy 'New consensus sequence' (majority / first / IUPAC)."""
    rs = recs(src)
    seqs = [r.seq for r in rs]
    if not seqs:
        return text("", "no sequences")
    cons = seq.consensus(seqs, mode="first" if mode == "first" else "majority") \
        if mode != "ambiguous" else seq.consensus_iupac(seqs)
    return text(cons + "\n", f"consensus of {len(seqs)} sequences ({mode})")


@T("fasta_nucleotide_percentage", "Nucleotide percentage plot data", FASTA_SECTIONS, "table",
   [fa("src"), intin("window", 200), intin("step", 100)], ex={"src": GENES_FA, "window": 300},
   up="gc_percent", summary="Sliding-window GC and nucleotide composition.")
def sliding_composition(src, window=200, step=100):
    """GC% along each sequence (windowed)."""
    out = []
    for r in recs(src):
        for w in seq.sliding_windows(r.seq, size=int(window), step=int(step) or None,
                                     fn=lambda s: round(100 * seq.gc_content(s), 3)):
            out.append({"id": r.id, "start": w["start"], "end": w["end"], "gc_percent": w["value"]})
    res = table(out, f"{len(out)} windows of {window} bp")
    if out:
        from chroma_titan.core import plot

        res["figure"] = plot.line({"GC %": [w["gc_percent"] for w in out[:500]]},
                                  x=[w["start"] for w in out[:500]],
                                  title="GC% along sequences", xlabel="position", ylabel="GC %")
    return res


@T("fasta_cpg_report", "CpG report (newcpgreport)", FASTA_SECTIONS, "table",
   [fa("src"), intin("window", 200), intin("step", 100), number("threshold", 0.6)],
   ex={"src": GENES_FA, "window": 200}, up="newcpgreport",
   summary="Windows with high CpG content and observed/expected CpG ratio.")
def cpg_report(src, window=200, step=100, threshold=0.6):
    """Islands of CpG richness per sequence (EMBOSS ``newcpgreport``)."""
    out = []
    for r in recs(src):
        for w in seq.sliding_windows(r.seq.upper(), size=int(window), step=int(step) or None,
                                     fn=lambda s: seq.cpg_obs_exp(s)):
            oe = w["value"].get("ratio", 0.0)
            gc = w["value"].get("gc", 0.0)
            if oe >= threshold * 100 and gc >= 40:
                out.append({"id": r.id, "start": w["start"], "end": w["end"],
                            "gc_percent": round(gc, 3), "cpg_obs_exp": round(oe, 4),
                            "cpg_count": w["value"].get("cg", 0)})
    return table(out, f"{len(out)} CpG-rich windows (obs/exp >= {threshold})")


@T("fasta_fasta_stats_of_records", "Enolpha/Phred quality conversion of FASTA", FASTA_SECTIONS, "text",
   [fa("src"), choice("direction", ["phred_to_fasta", "fasta_to_phred"]), textbox("qual", "!!!!!")],
   ex={"src": GENES_FA, "direction": "fasta_to_phred"}, up="fasta_to_fastq",
   summary="Attach a uniform quality string to FASTA (or drop qualities).")
def enolpha(src, direction="fasta_to_phred", qual="!!!!!"):
    """Convert between FASTA and FASTQ using a constant quality string."""
    rs = recs(src)
    if direction == "fasta_to_phred":
        q = qual if len(qual) else "!"
        return {"text": io.write_fastq([io.Read(r.id, r.seq, q * len(r.seq)) for r in rs]),
                "filename": "output.fastq", "message": f"{len(rs)} reads with quality {q!r}"}
    reads = io.parse_fastq(io.as_text(src))
    return {"text": io.write_fasta([io.Seq(x.id, x.seq) for x in reads]),
            "filename": "output.fasta", "message": f"{len(reads)} sequences without qualities"}


@T("fasta_barcode_split", "Barcode Split (demultiplex by inline barcodes)", FASTA_SECTIONS, "table",
   [fq("src", READS_SINGLE), bigtext("barcodes", "ACGT\nTGCA"), choice("where", ["start", "end"]),
    intin("mismatches", 0)], ex={"src": READS_SINGLE, "barcodes": "ACGT\nTGCA", "where": "start"},
   up="bctools_barcode", summary="Assign reads to samples by a 5'/3' barcode prefix.")
def barcode_split(src, barcodes="", where="start", mismatches=0):
    """Count and report reads per barcode bin."""
    reads = io.parse_fastq(io.as_text(src))
    bars = [b.strip().upper() for b in _lines_like(barcodes) if b.strip()]
    bins: dict[str, list] = defaultdict(list)
    unmatched = 0
    for rd in reads:
        s = rd.seq.upper()
        piece = s[: max(len(b) for b in bars)] if bars else ""
        found = None
        for b in bars:
            cand = s[:len(b)] if where == "start" else s[-len(b):]
            if sum(1 for x, y in zip(cand, b) if x != y) <= int(mismatches):
                found = b
                break
        if found:
            bins[found].append(rd)
        else:
            unmatched += 1
    out = [{"barcode": k, "reads": len(v), "mean_length": round(stats.mean([len(x.seq) for x in v]), 1),
            "first_read": v[0].id} for k, v in bins.items()]
    out.append({"barcode": "unassigned", "reads": unmatched, "mean_length": "", "first_read": ""})
    return table(out, f"demultiplexed {len(reads)} reads into {len(bins)} bins")


@T("fasta_random_subset", "Randomly choose N sequences", FASTA_SECTIONS, "fasta",
   [fa("src"), intin("n", 3), intin("seed", 1), number("fraction", 0.0), boolean("keep_order", False)],
   ex={"src": GENES_FA, "n": 3, "seed": 2}, up="random-choose", summary="Draw a random sample of records.")
def random_sequences(src, n=3, seed=1, fraction=0.0, keep_order=False):
    """Reservoir-style random subsetting of a FASTA file."""
    rs = recs(src)
    k = int(n) if n else int(len(rs) * fraction) if fraction else len(rs)
    k = max(1, min(k, len(rs)))
    picked = random.Random(seed).sample(rs, k)
    if keep_order:
        order = {r.id: i for i, r in enumerate(rs)}
        picked.sort(key=lambda r: order[r.id])
    res = _write(picked)
    res["message"] = f"{k} random sequences of {len(rs)}"
    return res


@T("fasta_extract_region", "Extract sequence by region", "fetch_sequences_alignments", "fasta",
   [fa("src"), textbox("regions", "chrV:100-200"), choice("strand", ["both", "plus", "minus"]),
    boolean("flank", False), intin("extend", 0)],
   ex={"src": GENOME, "regions": "chrV:100-200\nchrV:1000-1100"}, up="bedtools_getfasta",
   summary="Pull chrom:start-end slices out of a genome FASTA (``getfasta``).")
def extract_region_tool(src, regions="", strand="both", flank=False, extend=0):
    """Galaxy 'Fetch sequences align only to genomic regions'."""
    rs = recs(src)
    index = {r.id.upper(): r.seq.upper() for r in rs}
    out = []
    for line in _lines_like(regions):
        parts = line.split("\t")
        spec = parts[0].strip()
        m = re.match(r"^(\S+?)(?::(\d+)..?(\d+))?(?:\(([-+])\))?$", spec)
        if not m:
            continue
        chrom, a, b, st = m.group(1), m.group(2), m.group(3), m.group(4)
        if chrom.upper() not in index:
            continue
        s = index[chrom.upper()]
        lo = int(a) - 1 if a else 0
        hi = int(b) if b else len(s)
        if flank:
            lo = max(0, lo - int(extend))
            hi = min(len(s), hi + int(extend))
        piece = s[lo:hi]
        use_strand = st or ("+" if strand != "minus" else "-")
        if use_strand == "-":
            piece = seq.reverse_complement(piece)
        name = f"{chrom}:{lo + 1}-{hi}({use_strand})"
        out.append(io.Seq(name, piece))
        if strand == "both" and not st and len(piece) > 10:
            out.append(io.Seq(f"{chrom}:{lo + 1}-{hi}(-)", seq.reverse_complement(piece)))
    res = _write(out)
    res["message"] = f"{len(out)} regions extracted"
    return res


def _lines_like(x) -> list[str]:
    return [l for l in io.as_text(x).splitlines() if l.strip()]


@T("fasta_search_loop", "Find sequences matching a pattern", FASTA_SECTIONS, "table",
   [fa("src"), textbox("pattern", "GATC"), choice("search_in", ["sequence", "description", "id"]),
    boolean("reverse_complement", True), intin("max_hits", 50)],
   ex={"src": GENES_FA, "pattern": "ATG"}, up="fasta_search_loop",
   summary="Search a motif or regexp through the file and report hit positions.")
def search_sequences(src, pattern="", search_in="sequence", reverse_complement=True, max_hits=50):
    """Regexp search of a FASTA file with hit coordinates."""
    rx = re.compile(seq.iupac_to_regex(pattern) if search_in == "sequence" else pattern, re.I)
    out = []
    for r in recs(src):
        if search_in == "sequence":
            hits = [m.start() for m in rx.finditer(r.seq)]
            rc = seq.reverse_complement(r.seq) if reverse_complement else ""
            rhits = [len(r.seq) - m.end() for m in rx.finditer(rc)] if rc else []
            for h in hits[: max_hits]:
                out.append({"id": r.id, "strand": "+", "position": h, "matched": r.seq[h:h + len(pattern)]})
            for h in rhits[: max_hits]:
                out.append({"id": r.id, "strand": "-", "position": max(0, h), "matched": pattern})
        else:
            txt_ = r.desc if search_in == "description" else r.id
            if rx.search(txt_):
                out.append({"id": r.id, "strand": ".", "position": 0, "matched": txt_[:60]})
    return table(out[:500], f"{len(out)} hits in {len(recs(src))} sequences")


@T("fasta_kmer_count", "Nucleotide frequencies (k-mers)", FASTA_SECTIONS, "table",
   [fa("src"), intin("k", 4), boolean("canonical", True), intin("top", 20),
    boolean("normalize", True)], ex={"src": GENES_FA, "k": 4, "top": 15}, up="seqtk_seq",
   summary="Count k-mers over the file (oligonucleotide frequencies).")
def kmer_frequency(src, k=4, canonical=True, top=20, normalize=True):
    """Galaxy 'NucFreq'/'Kmer frequency' equivalent."""
    allseq = "".join(r.seq.upper() for r in recs(src))
    c = seq.kmer_counts(allseq, k=int(k), canonical=canonical, min_count=1)
    tot = sum(c.values()) or 1
    items = c.most_common(int(top)) if top else c.most_common()
    out = [{"kmer": km, "count": n, "frequency": round(n / tot, 6),
            "expected_uniform": round(1 / 4 ** int(k), 6)} for km, n in items]
    if normalize:
        mx = max((o["count"] for o in out), default=1) or 1
        for o in out:
            o["relative"] = round(o["count"] / mx, 4)
    return table(out, f"top {len(out)} {k}-mers of {len(allseq)} bases")


@T("fasta_profile_vector", "Profile vector (pseudo-monomer composition)", FASTA_SECTIONS, "table",
   [fa("src"), intin("k", 4), boolean("revcomp_aware", True)], ex={"src": GENES_FA, "k": 3},
   up="profile_vector", summary="Composition vector for sequence similarity work.")
def profile_vector(src, k=4, revcomp_aware=True):
    """Fixed-length k-mer profile for each record (useful for binning)."""
    out = []
    for r in recs(src):
        v = seq.profile_vector(r.seq.upper(), k=int(k))
        out.append({"id": r.id, "length": len(r.seq), "n_values": len(v),
                    "sum": round(sum(v), 4), "max": round(max(v), 4) if v else 0.0,
                    "profile": ",".join(f"{x:.4f}" for x in v[:64])})
    return table(out, f"profiles (k={k}) for {len(out)} sequences")


@T("fasta_mash_compare", "MinHash sketch and distance matrix", FASTA_SECTIONS, "table",
   [fa("src"), intin("k", 21), intin("num", 100), boolean("distance_matrix", True)],
   ex={"src": GENES_FA, "k": 15, "num": 50}, up="tools-mash",
   summary="Sketch each sequence and estimate pairwise Mash distances.")
def mash_compare(src, k=21, num=100, distance_matrix=True):
    """MinHash/Mash-style comparison without external tools."""
    rs = recs(src)
    sketches = {r.id: seq.minhash(r.seq.upper(), k=int(k), num=int(num)) for r in rs}
    out = []
    for a in rs:
        for b in rs:
            if a.id >= b.id:
                continue
            j = seq.minhash_jaccard(sketches[a.id], sketches[b.id])
            d = seq.mash_distance(sketches[a.id], sketches[b.id], k=int(k))
            out.append({"sequence_1": a.id, "sequence_2": b.id, "jaccard": round(j, 4),
                        "mash_distance": round(d, 4), "shared_hashes": int(round(j * int(num)))})
    return table(out, f"{len(rs)} sketches, {len(out)} pairs (k={k})")


@T("fasta_md5", "Sequence MD5 checksums", FASTA_SECTIONS, "table",
   [fa("src"), choice("algorithm", ["md5", "sha1"])], ex={"src": GENES_FA}, up="seqtk_md5",
   summary="Hash of every sequence (identical data → identical hash).")
def sequence_md5(src, algorithm="md5"):
    """Checksum each record for de-duplication audits."""
    out = [{"id": r.id, "length": len(r.seq), "hash": seq.sha1_of(r.seq.upper())
            if algorithm == "sha1" else seq.sequence_md5_id(r.seq.upper())} for r in recs(src)]
    return table(out, f"{len(out)} checksums")


@T("fasta_encode_aa_props", "Encode amino-acid properties", FASTA_SECTIONS, "table",
   [fa("src", PROTEINS), choice("property", ["hydrophobicity", "charge", "volume", "flexibility",
                                            "polarity"]), intin("window", 5)],
   ex={"src": PROTEINS, "property": "hydrophobicity", "window": 5}, up="aa_prop",
   summary="Numeric per-residue property track (Galaxy 'Encode Aminoacid Properties').")
def encode_aa_properties(src, property="hydrophobicity", window=5):
    """Kyte-Doolittle-like sliding property profile for each protein."""
    from chroma_titan.core import protein

    scales = {"hydrophobicity": {"A": 1.8, "R": -4.5, "N": -3.5, "D": -3.5, "C": 2.5, "Q": -3.5,
                                 "E": -3.5, "G": -0.4, "H": -3.2, "I": 4.5, "L": 3.8, "K": -3.9,
                                 "M": 1.9, "F": 2.8, "P": -1.6, "S": -0.8, "T": -0.7, "W": -0.9,
                                 "Y": -1.3, "V": 4.2},
              "charge": {"D": -1, "E": -1, "K": 1, "R": 1, "H": 0.1},
              "volume": {"G": 60, "A": 89, "S": 89, "C": 109, "D": 111, "P": 113, "N": 114,
                         "T": 116, "E": 138, "V": 140, "M": 163, "K": 169, "I": 167, "L": 167,
                         "R": 174, "F": 190, "Y": 194, "W": 228, "H": 153, "Q": 144},
              "flexibility": {"G": 1, "P": 0.5, "S": 0.75, "A": 0.35, "V": 0.3},
              "polarity": {"N": 1, "Q": 1, "S": 1, "T": 1, "Y": 1, "C": 1, "K": 1, "R": 1, "D": 1,
                           "E": 1, "H": 1}}
    table_ = scales.get(property, scales["hydrophobicity"])
    out = []
    for r in recs(src):
        vals = [table_.get(c.upper(), 0.0) for c in r.seq]
        w = max(1, int(window))
        prof = [round(stats.mean(vals[i:i + w]), 3) for i in range(0, max(1, len(vals) - w + 1))]
        out.append({"id": r.id, "length": len(r.seq), "mean": round(stats.mean(vals), 3),
                    "max": round(max(vals), 3) if vals else 0.0,
                    "profile": ",".join(map(str, prof[:80]))})
    return table(out, f"{property} profiles ({len(out)} sequences)")


@T("fasta_string_tie", "Concatenate sequences (String Tying)", FASTA_SECTIONS, "fasta",
   [fa("src"), boolean("revcomp_second", False), intin("gap_nt", 0), textbox("gap_char", "N")],
   ex={"src": GENES_FA, "gap_nt": 10}, up="string_tie",
   summary="Tie records end to end into one sequence with a spacer.")
def string_tie(src, revcomp_second=False, gap_nt=0, gap_char="N"):
    """Galaxy 'Fasta String Tying'."""
    rs = recs(src)
    parts = []
    for i, r in enumerate(rs):
        if i and gap_nt:
            parts.append(gap_char * int(gap_nt))
        parts.append(seq.reverse_complement(r.seq) if (revcomp_second and i % 2) else r.seq)
    return _write([io.Seq("tied", "".join(parts), f"{len(rs)} sequences joined")])


@T("fasta_filter_on_description", "Filter FASTA by description text", FASTA_SECTIONS, "fasta",
   [fa("src"), textbox("pattern", "gene"), boolean("invert", False)], ex={"src": GENES_FA},
   up="fasta-rc", summary="Keep records whose defline matches a regexp.")
def filter_on_description(src, pattern="", invert=False):
    """Defline-based filtering."""
    rx = re.compile(pattern)
    out = [r for r in recs(src) if bool(rx.search(r.desc or r.id)) != invert]
    res = _write(out)
    res["message"] = f"{len(out)} sequences match /{pattern}/"
    return res


@T("fasta_to_fastq_convert", "Convert between FASTA and FASTQ", "convert_formats", "text",
   [fq("src", READS_SINGLE), choice("quality_mode", ["constant", "index_derived", "zero"]),
    textbox("qual", "I")], ex={"src": READS_SINGLE, "quality_mode": "constant"}, up="fastq_to_fasta",
   summary="Drop or synthesise the quality lines when converting.")
def fasta_fastq_convert(src, quality_mode="constant", qual="I"):
    """FASTQ → FASTA and back (with fake qualities)."""
    body = io.as_text(src).lstrip()
    if body.startswith("@"):
        reads = io.parse_fastq(body)
        return {"text": io.write_fasta([io.Seq(r.id, r.seq) for r in reads]),
                "filename": "output.fasta", "message": f"{len(reads)} reads → FASTA"}
    rs = io.parse_fasta(body)
    q = []
    for r in rs:
        if quality_mode == "constant":
            q.append((qual or "I") * len(r.seq))
        elif quality_mode == "zero":
            q.append("!" * len(r.seq))
        else:
            q.append("".join(chr(33 + min(40, (ord(c) % 40))) for c in r.seq))
    return {"text": io.write_fastq([io.Read(r.id, r.seq, qq) for r, qq in zip(rs, q)]),
            "filename": "output.fastq", "message": f"{len(rs)} sequences → FASTQ"}


@T("convert_bed_to_fasta", "BED intervals to FASTA (genome required)", "convert_formats", "fasta",
   [bed("src"), fa("genome", GENOME), boolean("name_from_score", False)],
   ex={"src": REGIONS_BED, "genome": GENOME}, up="bedtools_getfasta",
   summary="Extract the sequence of every interval from a genome FASTA.")
def bed_to_fasta(src, genome, name_from_score=False):
    """BED → FASTA."""
    idx = {r.id.upper(): r.seq.upper() for r in io.parse_fasta(io.as_text(genome))}
    out = []
    for iv in io.parse_bed(io.as_text(src)):
        s = idx.get(iv.chrom.upper(), "")
        piece = s[iv.start:iv.end]
        if iv.strand == "-":
            piece = seq.reverse_complement(piece)
        out.append(io.Seq(iv.name if not name_from_score else f"{iv.chrom}:{iv.start}-{iv.end}",
                          piece, f"{iv.chrom}:{iv.start + 1}-{iv.end}{iv.strand}"))
    res = _write(out)
    res["message"] = f"{len(out)} intervals extracted"
    return res


@T("convert_gff_to_bed", "Convert GFF3 to BED12", "convert_formats", "bed",
   [gff("src"), textbox("feature", "mRNA"), boolean("merge_exons", True)],
   ex={"src": ANNOT_GFF, "feature": "mRNA"}, up="gff2bed",
   summary="Turn annotated gene/exon features into BED blocks.")
def gff_to_bed_tool(src, feature="mRNA", merge_exons=True):
    """GFF3 → BED12 via exon blocks."""
    rows = io.parse_gff(io.as_text(src))
    if feature:
        rows = [g for g in rows if g.type == feature or (feature in ("mRNA", "gene") and g.type in ("mRNA", "gene", "exon", "CDS"))]
    ivs = io.gff_to_bed12(rows) if merge_exons else \
        [io.Interval(g.seqid, g.start - 1, g.end, g.attrs.get("ID", g.attrs.get("gene_id", "")),
                     0, g.strand) for g in rows]
    return {"text": io.write_bed(ivs), "filename": "output.bed",
            "message": f"{len(ivs)} BED features"}


@T("convert_bed_to_gff", "Convert BED to GFF3", "convert_formats", "gff",
   [bed("src"), textbox("source", "chroma_titan"), textbox("feature", "region")],
   ex={"src": REGIONS_BED, "feature": "exon"}, up="bed2gff", summary="BED → GFF3 with name attributes.")
def bed_to_gff_tool(src, source="chroma_titan", feature="region"):
    """BED → GFF3."""
    from chroma_titan.core import genome as _g

    rows = _g.bed_to_gff(io.parse_bed(io.as_text(src)), source=source)
    if feature:
        rows = [io.GffRow(g.seqid, g.source, feature, g.start, g.end, g.score, g.strand, g.phase,
                          g.attrs) for g in rows]
    return {"text": io.write_gff(rows), "filename": "output.gff3",
            "message": f"{len(rows)} GFF features"}


@T("convert_table_to_bed", "Tabular to BED", "convert_formats", "bed",
   [tbl("src"), textbox("chrom_col", "1"), textbox("start_col", "2"), textbox("end_col", "3"),
    textbox("name_col", "4"), textbox("strand_col", "6"), boolean("start_is_1_based", False)],
   ex={"src": REGIONS_BED}, up="table_to_bed", summary="Map arbitrary table columns onto BED fields.")
def table_to_bed(src, chrom_col="1", start_col="2", end_col="3", name_col="4", strand_col="6",
                 start_is_1_based=False):
    """Generic table → BED."""
    rows = [l.split("\t") for l in io.lines(src)]
    if not rows:
        return {"text": "", "message": "empty table"}
    ncol = max(len(r) for r in rows)
    has_head = not all(len(r) > 1 and r[1].strip().isdigit() for r in rows[:2])
    body = rows[1:] if has_head else rows

    def col(spec):
        s = str(spec).strip()
        if has_head and s in rows[0]:
            return rows[0].index(s)
        return int(s) - 1 if s.isdigit() else 0
    ci, si, ei, ni, sti = map(col, (chrom_col, start_col, end_col, name_col, strand_col))
    out = []
    for r in body:
        try:
            st = int(float(r[si])) - (1 if start_is_1_based else 0)
            en = int(float(r[ei]))
        except (ValueError, IndexError):
            continue
        out.append(io.Interval(r[ci] if ci < len(r) else ".", st, en,
                               r[ni] if ni < len(r) else ".", 0,
                               (r[sti] if sti < len(r) else ".")[:1] or "."))
    return {"text": io.write_bed(out, bed12=False), "filename": "output.bed",
            "message": f"{len(out)} BED intervals"}


@T("convert_newick_to_table", "Newick tree to branch table", "convert_formats", "table",
   [anyfile("src", TREE, label="Newick file")], ex={"src": TREE}, up="newick2table",
   summary="Flatten a tree into a node/parent/branch-length table.")
def newick_to_table(src):
    """Newick → tabular edge list."""
    tree = io.parse_newick(io.as_text(src))
    edges = []

    def walk(node, parent="ROOT", depth=0):
        for ch in node.get("children", []):
            edges.append({"parent": parent, "child": ch.get("name") or "", "length": ch.get("length"),
                          "depth": depth + 1, "n_leaves": count_leaves(ch)})
            walk(ch, ch.get("name") or parent, depth + 1)

    def count_leaves(node):
        if not node.get("children"):
            return 1
        return sum(count_leaves(c) for c in node["children"])

    walk(tree)
    return table(edges, f"{len(edges)} branches")


@T("convert_table_to_newick", "Table of branches to Newick", "convert_formats", "text",
   [tbl("src"), textbox("parent_col", "parent"), textbox("child_col", "child"),
    textbox("length_col", "length")], ex={"src": REGIONS_BED}, up="table2newick",
   summary="Rebuild a Newick string from a parent/child edge list.")
def table_to_newick(src, parent_col="parent", child_col="child", length_col="length"):
    """Edge table → Newick."""
    rows = [l.split("\t") for l in io.lines(src)]
    if not rows:
        return text("", "empty table")
    head = rows[0]
    if parent_col not in head:
        return text("", "table must have parent/child/length columns")
    pi, ci = head.index(parent_col), head.index(child_col)
    li = head.index(length_col) if length_col in head else None
    kids: dict[str, list[tuple[str, float]]] = defaultdict(list)
    for r in rows[1:]:
        if len(r) <= max(pi, ci):
            continue
        try:
            ln = float(r[li]) if li is not None and li < len(r) and r[li] else 0.0
        except ValueError:
            ln = 0.0
        kids[r[pi]].append((r[ci], ln))

    def build(name):
        ch = kids.get(name, [])
        if not ch:
            return name
        inner = ",".join(f"{build(n)}:{l:g}" for n, l in ch)
        return f"({inner}){name if name != 'ROOT' else ''}"

    root = "ROOT" if "ROOT" in kids else next(iter(kids))
    return text(build(root) + ";\n", f"tree with {len(rows) - 1} edges")


@T("convert_matrix_to_dist", "Distance matrix formats (square ↔ phylip/lower)",
   "convert_formats", "text", [anyfile("src", MSA), choice("to", ["phylip_lower", "square", "vector"])],
   ex={"src": MSA, "to": "square"}, up="convert_matrix",
   summary="Rewrite a pairwise distance matrix between the usual layouts.")
def convert_matrix(src, to="square"):
    """Matrix reformatting for PHYLIP-style distance input."""
    from chroma_titan.core import align

    rs = io.parse_fasta(io.as_text(src))
    D = align.distance_matrix([r.seq for r in rs], model="identity")
    names = [r.id for r in rs]
    if to == "square":
        w = max(len(n) for n in names) + 1
        body = ["\n".join([f"{len(names):>6}"] +
                          [n.ljust(w) + " ".join(f"{v:.4f}" for v in row) for n, row in zip(names, D)])]
        return text(body[0] + "\n", "square distance matrix")
    if to == "vector":
        rows = [{"sequence_1": names[i], "sequence_2": names[j], "distance": round(D[i][j], 5)}
                for i in range(len(names)) for j in range(i + 1, len(names))]
        return table(rows, f"{len(rows)} pairwise distances")
    out = ["\n".join(f"{names[i]:<16}" + " ".join(f"{D[i][j]:.4f}" for j in range(i))
                     for i in range(1, len(names)))]
    return text(f"{len(names)}\n" + out[0] + "\n", "lower-triangle matrix")


# ===========================================================================
# seqtk toolkit (Galaxy 'seqtk' section)
# ===========================================================================
@T("seqtk_seq_subsample", "seqtk sample", "seqtk", "text",
   [fq("src", READS_SINGLE), number("fraction", 0.25), intin("n", 0), intin("seed", 100)],
   ex={"src": READS_SINGLE, "fraction": 0.25, "seed": 7}, up="iuc/seqtk",
   summary="Randomly subsample reads/sequences with a fixed seed.")
def seqtk_sample(src, fraction=0.25, n=0, seed=100):
    """``seqtk sample -s<seed> file <frac|n>``."""
    body = io.as_text(src).lstrip()
    is_fastq = body.startswith("@")
    items = io.parse_fastq(body) if is_fastq else io.parse_fasta(body)
    k = int(n) if n else max(1, int(len(items) * fraction))
    picked = random.Random(seed).sample(items, min(k, len(items)))
    txt_ = io.write_fastq(picked) if is_fastq else io.write_fasta(picked)
    return {"text": txt_, "filename": "output.fastq" if is_fastq else "output.fasta",
            "message": f"sampled {len(picked)}/{len(items)}"}


@T("seqtk_comp", "seqtk comp", "seqtk", "table", [fa("src")], ex={"src": GENES_FA}, up="iuc/seqtk",
   summary="Per-record length, GC and per-base counts (``seqtk comp``).")
def seqtk_comp(src):
    """``seqtk comp``: one row per sequence with base counts."""
    out = []
    for r in recs(src):
        c = seq.counts(r.seq)
        out.append({"id": r.id, "length": len(r.seq), "GC": round(100 * seq.gc_content(r.seq), 3),
                    "A": c["A"], "C": c["C"], "G": c["G"], "T": c["T"], "other": c["N"]})
    return table(out, f"{len(out)} records")


@T("seqtk_trim_qual", "seqtk trimq", "seqtk", "text", [fq("src"), intin("min_quality", 20),
                                                     boolean("trim_two_sided", True)],
   ex={"src": READS_SINGLE, "min_quality": 20}, up="iuc/seqtk",
   summary="Trim low quality tails from FASTQ reads.")
def seqtk_trimq(src, min_quality=20, trim_two_sided=True):
    """``seqtk trimq -c <q>``."""
    reads = io.parse_fastq(io.as_text(src))
    out = []
    for rd in reads:
        s, q = rd.seq, rd.qual
        while q and ord(q[0]) - 33 < min_quality:
            s, q = s[1:], q[1:]
        while q and ord(q[-1]) - 33 < min_quality:
            s, q = s[:-1], q[:-1]
        if len(s) >= 20:
            out.append(io.Read(rd.id, s, q, rd.desc))
    return {"text": io.write_fastq(out), "filename": "trimmed.fastq",
            "message": f"{len(out)}/{len(reads)} reads kept"}


@T("seqtk_listheads", "seqtk listh (headers)", "seqtk", "text", [fa("src")], ex={"src": GENES_FA},
   up="iuc/seqtk", summary="List all sequence headers.")
def seqtk_listheads(src):
    """``seqtk listh``: deflines only."""
    lines = [f">{r.id} {r.desc}".strip() for r in recs(src)]
    return text("\n".join(lines) + "\n", f"{len(lines)} headers")


@T("seqtk_fixhead", "seqtk fixhead (repair deflines)", "seqtk", "fasta",
   [fa("src"), textbox("prefix", "seq"), intin("start", 1)], ex={"src": GENES_FA}, up="iuc/seqtk",
   summary="Replace duplicate or invalid ids with sequential names.")
def seqtk_fixhead(src, prefix="seq", start=1):
    """Make ids unique by renaming sequentially."""
    out = [io.Seq(f"{prefix}{i}", r.seq, r.desc) for i, r in enumerate(recs(src), int(start))]
    return _write(out)


@T("seqtk_mut", "seqtk mutfa (apply mutations)", "seqtk", "fasta",
   [fa("src"), bigtext("mutations", "", help="id\\tposition\\told\\tnew")],
   ex={"src": GENES_FA, "mutations": ""}, up="iuc/seqtk",
   summary="Introduce substitutions/indels described in a mutation list.")
def seqtk_mut(src, mutations=""):
    """``seqtk mutfa``-style in-silico mutagenesis."""
    rs = {r.id: list(r.seq) for r in recs(src)}
    n = 0
    for line in _lines_like(mutations):
        parts = line.split("\t")
        if len(parts) < 2 or parts[0] not in rs:
            continue
        pos = int(float(parts[1]))
        new = parts[3] if len(parts) > 3 else parts[2]
        rs[parts[0]][pos - 1:pos] = list(new)
        n += 1
    out = [io.Seq(k, "".join(v)) for k, v in rs.items()]
    res = _write(out)
    res["message"] = f"{n} mutations applied"
    return res


@T("seqtk_dna2ca", "seqtk dna2ca (codon conversion)", "seqtk", "table",
   [fa("src"), choice("mode", ["count_codons", "gc3", "cai"])], ex={"src": GENES_FA, "mode": "gc3"},
   up="iuc/seqtk", summary="Codon-level statistics per sequence.")
def seqtk_codon_stats(src, mode="gc3"):
    """Codon usage, GC3 or CAI per record."""
    out = []
    for r in recs(src):
        if mode == "gc3":
            out.append({"id": r.id, "gc3": round(seq.gc3(r.seq), 4)})
        elif mode == "cai":
            out.append({"id": r.id, "cai": round(seq.cai(r.seq), 4), "enc": round(
                seq.effective_number_of_codons(r.seq), 2)})
        else:
            cu = seq.codon_usage(r.seq)
            out.append({"id": r.id, "n_codons": sum(cu.values()), "unique_codons": len(cu),
                        "most_common": (cu.most_common(1)[0][0] if cu else ""),
                        "table": "|".join(f"{k}:{v}" for k, v in cu.most_common(8))})
    return table(out, f"{mode} for {len(out)} sequences")


@T("seqtk_locut", "Low-complexity filter (locut/DUST)", "seqtk", "table",
   [fa("src"), number("dust_cutoff", 2.0), intin("window", 64)], ex={"src": GENES_FA, "dust_cutoff": 3.0},
   up="iuc/dust", summary="Flag low complexity sequences/windows like dustmasker.")
def locut(src, dust_cutoff=2.0, window=64):
    """Report DUST scores and low-complexity windows."""
    out = []
    for r in recs(src):
        d = seq.dust_score(r.seq, window=int(window))
        lc = seq.low_complexity(r.seq, window=int(window), level=2.4)
        out.append({"id": r.id, "length": len(r.seq), "dust": round(d, 3),
                    "low_complexity": bool(d < dust_cutoff), "n_lc_windows": len(lc),
                    "lc_fraction": round(sum(x["end"] - x["start"] for x in lc) / max(1, len(r.seq)), 4)})
    return table(out, f"{sum(1 for o in out if o['low_complexity'])} low-complexity sequences")


@T("seqtk_maskedcopy", "Mask low complexity regions", "seqtk", "fasta",
   [fa("src"), number("level", 2.4), intin("window", 64), textbox("mask_char", "N")],
   ex={"src": GENES_FA, "level": 2.4}, up="iuc/dustmasker",
   summary="Replace dust-masked segments with N (``maskedcopy -m``).")
def masked_copy(src, level=2.4, window=64, mask_char="N"):
    """Soft/hard masking of low complexity sequence."""
    out = [io.Seq(r.id, seq.mask_low_complexity(r.seq, window=int(window), level=level,
                                                 char=mask_char or "N"), r.desc) for r in recs(src)]
    res = _write(out)
    res["message"] = f"masked {sum(1 for a, b in zip(recs(src), out) if a.seq != b.seq)} sequences"
    return res


@T("seqtk_gapidx", "Report runs of N (gap index)", "seqtk", "table",
   [fa("src"), intin("min_len", 1)], ex={"src": GENOME, "min_len": 10}, up="iuc/seqtk",
   summary="List the coordinates of every gap in the sequences.")
def gap_index(src, min_len=1):
    """``seqtk gapidx``: N runs with 1-based coordinates."""
    out = []
    for r in recs(src):
        for blk in seq.find_n_blocks(r.seq, min_len=int(min_len)):
            out.append({"id": r.id, "start": blk["start"] + 1, "end": blk["end"],
                        "length": blk["end"] - blk["start"]})
    return table(out, f"{len(out)} gap blocks")


# ===========================================================================
# FASTQ quality-control-adjacent conversions used by the FASTA/FASTQ section
# ===========================================================================
@T("fastq_to_fasta", "FASTQ to FASTA", FASTA_SECTIONS, "fasta", [fq("src"), boolean("keep_first_word", True)],
   ex={"src": READS_SINGLE}, up="iuc/fastx_toolkit", summary="Drop quality lines and keep the sequence.")
def fastq_to_fasta(src, keep_first_word=True):
    """Lossless sequence extraction from FASTQ."""
    reads = io.parse_fastq(io.as_text(src))
    out = [io.Seq(r.id.split()[0] if keep_first_word else r.id, r.seq,
                  " ".join(r.id.split()[1:])) for r in reads]
    return _write(out)


@T("fastq_filter", "Filter FASTQ reads", FASTA_SECTIONS, "text",
   [fq("src"), intin("min_length", 50), intin("max_length", 0), number("min_mean_qual", 0.0),
    number("max_n_fraction", 1.0), textbox("contains", ""), boolean("remove_gaps", False)],
   ex={"src": READS_SINGLE, "min_length": 60, "min_mean_qual": 20}, up="iuc/fastx_toolkit",
   summary="Length / quality / N-content / substring filtering of reads.")
def fastq_filter_tool(src, min_length=0, max_length=0, min_mean_qual=0.0, max_n_fraction=1.0,
                      contains="", remove_gaps=False):
    """Compound read filter."""
    reads = io.parse_fastq(io.as_text(src))
    out = []
    for r in reads:
        s = r.seq.replace("-", "") if remove_gaps else r.seq
        if len(s) < min_length or (max_length and len(s) > max_length):
            continue
        if min_mean_qual and stats.mean([ord(c) - 33 for c in r.qual]) < min_mean_qual:
            continue
        nfrac = s.upper().count("N") / max(1, len(s))
        if nfrac > max_n_fraction:
            continue
        if contains and contains.upper() not in s.upper():
            continue
        out.append(io.Read(r.id, s, r.qual, r.desc))
    return {"text": io.write_fastq(out), "filename": "filtered.fastq",
            "message": f"{len(out)}/{len(reads)} reads pass the filter"}


@T("fastq_trim_length", "Truncate reads", FASTA_SECTIONS, "text",
   [fq("src"), intin("length", 0), intin("left_trim", 0), intin("right_trim", 0),
    boolean("pad_short_reads", False)], ex={"src": READS_SINGLE, "length": 60}, up="fastx_trimmer",
   summary="Fixed truncation and end trimming of FASTQ reads.")
def fastq_trim_length(src, length=0, left_trim=0, right_trim=0, pad_short_reads=False):
    """``fastx_trimmer`` behaviour."""
    reads = io.parse_fastq(io.as_text(src))
    out = []
    for r in reads:
        s, q = r.seq[int(left_trim):len(r.seq) - int(right_trim) if right_trim else len(r.seq)], \
            r.qual[int(left_trim):len(r.qual) - int(right_trim) if right_trim else len(r.qual)]
        if length:
            if pad_short_reads and len(s) < length:
                s = s + "N" * (length - len(s))
                q = q + "!" * (length - len(q))
            s, q = s[: int(length)], q[: int(length)]
        out.append(io.Read(r.id, s, q, r.desc))
    return {"text": io.write_fastq(out), "filename": "trimmed.fastq",
            "message": f"trimmed {len(out)} reads"}


@T("fastq_reverse_complement", "Reverse complement reads", FASTA_SECTIONS, "text",
   [fq("src")], ex={"src": READS_SINGLE}, up="fastx_revcomp",
   summary="Reverse-complement both sequence and quality strings.")
def fastq_revcomp(src):
    """Read-level reverse complement."""
    reads = io.parse_fastq(io.as_text(src))
    out = [io.Read(r.id, *seq.reverse_complement_fastq(r.seq, r.qual), r.desc) for r in reads]
    return {"text": io.write_fastq(out), "filename": "revcomp.fastq", "message": f"{len(out)} reads"}


@T("fastq_sort", "Sort FASTQ by read name or quality", FASTA_SECTIONS, "text",
   [fq("src"), choice("by", ["name", "quality", "length", "sequence"]), boolean("reverse", False)],
   ex={"src": READS_SINGLE, "by": "quality"}, up="fastq_sort",
   summary="Deterministic ordering of reads.")
def fastq_sort(src, by="name", reverse=False):
    """Order reads (paired-end workflows need name sorting)."""
    reads = io.parse_fastq(io.as_text(src))
    key = {"name": lambda r: r.id, "quality": lambda r: stats.mean([ord(c) - 33 for c in r.qual]),
           "length": lambda r: len(r.seq), "sequence": lambda r: r.seq}[by]
    return {"text": io.write_fastq(sorted(reads, key=key, reverse=reverse)),
            "filename": "sorted.fastq", "message": f"sorted {len(reads)} reads by {by}"}


@T("fastq_split_pairs", "Split paired-end FASTQ into two files", FASTA_SECTIONS, "text",
   [fq("src"), textbox("suffix_1", "_1"), textbox("suffix_2", "_2"), boolean("by_name", True)],
   ex={"src": READS_1, "suffix_1": "_1"}, up="fastx",
   summary="Interleaved → separate R1/R2 files (reported as a table).")
def split_pairs(src, suffix_1="_1", suffix_2="_2", by_name=True):
    """Split interleaved reads back into mates."""
    reads = io.parse_fastq(io.as_text(src))
    r1, r2, un = [], [], []
    for r in reads:
        if by_name and r.id.endswith(suffix_2):
            r2.append(r)
        elif by_name and r.id.endswith(suffix_1):
            r1.append(r)
        else:
            (r1 if len(r1) <= len(r2) else r2).append(r)
    return table([{"file": "mates_1", "reads": len(r1), "example": r1[0].id if r1 else ""},
                  {"file": "mates_2", "reads": len(r2), "example": r2[0].id if r2 else ""},
                  {"file": "unpaired", "reads": len(un), "example": ""}],
                 f"{len(r1)} R1 + {len(r2)} R2 reads")


@T("fastq_join_paired", "Join paired reads (overlap merging)", FASTA_SECTIONS, "text",
   [fq("src", READS_1), fq("src2", READS_2), intin("min_overlap", 11), intin("max_mismatches", 2)],
   ex={"src": READS_1, "src2": READS_2, "min_overlap": 11}, up="fastq_join",
   summary="Merge R1/R2 using the exact-overlap algorithm (``fastq_join``/PEAR-like).")
def join_paired_fastq(src, src2, min_overlap=11, max_mismatches=2):
    """Paired-end read merging with quality consensus."""
    from chroma_titan.core import align

    r1 = io.parse_fastq(io.as_text(src))
    r2 = io.parse_fastq(io.as_text(src2))
    res = align.merge_pairs([{"id": x.id, "seq": x.seq, "qual": x.qual} for x in r1],
                            [{"id": x.id, "seq": x.seq, "qual": x.qual} for x in r2],
                            min_overlap=int(min_overlap), max_mismatches=int(max_mismatches))
    merged = res.get("merged", []) if isinstance(res, dict) else []
    lines = []
    for m in merged:
        lines.append(f"@{m['id']}\n{m['seq']}\n+\n{m.get('qual', 'I' * len(m['seq']))}")
    stats_ = {k: v for k, v in res.items() if isinstance(v, (int, float))}
    return {"text": "\n".join(lines) + ("\n" if lines else ""), "filename": "merged.fastq",
            "message": f"merged {len(merged)}/{len(r1)} read pairs", "stats": stats_}


@T("fastq_convert_illumina18", "Convert FASTQ quality encoding", "convert_formats", "text",
   [fq("src"), choice("from_offset", [33, 64]), choice("to_offset", [33, 64]),
    intin("max_quality", 41)], ex={"src": READS_SINGLE, "from_offset": 33, "to_offset": 64},
   up=" fastaq", summary="Shift ASCII quality offsets (Illumina 1.3/1.5/1.8, Solexa).")
def convert_quality_encoding(src, from_offset=33, to_offset=64, max_quality=41):
    """Re-encode Phred qualities between offsets."""
    reads = io.parse_fastq(io.as_text(src), offset=int(from_offset))
    shift = int(to_offset) - int(from_offset)
    out = []
    for r in reads:
        q = "".join(chr(max(int(to_offset), min(int(max_quality) + int(to_offset), ord(c) + shift)))
                    for c in r.qual)
        out.append(io.Read(r.id, r.seq, q, r.desc))
    return {"text": io.write_fastq(out), "filename": f"offset{to_offset}.fastq",
            "message": f"re-encoded {len(out)} reads from Q{from_offset} to Q{to_offset}"}


@T("fastq_to_sam_like", "Convert FASTQ to an unmapped SAM", "convert_formats", "text",
   [fq("src"), fa("genome", GENOME), intin("expected_length", 0)], ex={"src": READS_SINGLE, "genome": GENOME},
   up="fastq_to_sam", summary="Emit a SAM header for the reference plus unmapped records.")
def fastq_to_unmapped_sam(src, genome, expected_length=0):
    """FASTQ → unmapped SAM (useful before an aligner)."""
    reads = io.parse_fastq(io.as_text(src))
    refs = [(r.id, len(r.seq)) for r in io.parse_fasta(io.as_text(genome))]
    lines = ["@HD\tVN:1.6\tSO:unsorted"] + [f"@RG\tID:1\tSM:sample\tPL:ILLUMINA"] + \
        [f"@SQ\tSN:{n}\tLN:{l}" for n, l in refs]
    for r in reads:
        lines.append("\t".join([r.id.split()[0], "4", "*", "0", "0", "*", "=", "0", str(len(r.seq)),
                               r.seq, r.qual]))
    return {"text": "\n".join(lines) + "\n", "filename": "unmapped.sam",
            "message": f"{len(reads)} unmapped alignments, {len(refs)} references"}
