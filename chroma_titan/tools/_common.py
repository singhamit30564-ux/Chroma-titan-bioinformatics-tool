"""Shared declarations for tool modules (example files, input widgets, results)."""

from __future__ import annotations

from typing import Any, Callable

from chroma_titan.core import genome, io, stats, tables
from chroma_titan.registry import Input, Tool, tool

# ---------------------------------------------------------------------------
# bundled example datasets (see scripts/make_examples.py)
# ---------------------------------------------------------------------------
EX = "examples/"
GENOME = EX + "genome.fa"
GENOME2 = EX + "genome2.fa"
GENES_FA = EX + "genes.fa"
REGIONS_BED = EX + "regions.bed"
TARGETS_BED = EX + "targets.bed"
ANNOT_GFF = EX + "annotation.gff"
VARIANTS_VCF = EX + "variants.vcf"
READS_1 = EX + "reads_1.fastq"
READS_2 = EX + "reads_2.fastq"
READS_SINGLE = EX + "reads_single.fastq"
ALN_SAM = EX + "alignments.sam"
COUNTS = EX + "counts.tsv"
GMT = EX + "genesets.gmt"
PFM = EX + "motif.pfm"
JASPAR = EX + "motifs_jaspar.txt"
PROTEINS = EX + "proteins.faa"
REF_PROTEINS = EX + "reference_proteins.faa"
MSA = EX + "msa.fasta"
TREE = EX + "tree.nwk"
PAIRS = EX + "contacts.pairs"
BEDGRAPH = EX + "coverage.bedgraph"
IMAGE = EX + "cells.pgm"
PEPTIDES = EX + "peptides.tsv"
PHENO = EX + "phenotypes.tsv"
OBO = EX + "phenotype.obo"
STR = EX + "str_profile.tsv"
CXCG = EX + "methylation.cxcg"
GENOME_TXT = EX + "genome.txt"
TAXMAP = EX + "taxmap.tsv"

SHORT_DNA = "ATGAAAGCTTTGCGATCGATCGATCGGCTAAGCATCGATCGATCGATTAA"
SHORT_PROT = "MKALLRVIWRQSMDDAEVPLQR"
SHORT_REGION = "chrV:401-600"

# ---------------------------------------------------------------------------
# input widget factories
# ---------------------------------------------------------------------------
def fa(name: str = "src", default: str = GENOME, label: str = "FASTA file or pasted FASTA",
       help: str = "") -> Input:
    return Input(name, kind="file", fmt="fasta", default=default, label=label, help=help)


def fq(name: str = "src", default: str = READS_SINGLE, label: str = "FASTQ file",
       help: str = "") -> Input:
    return Input(name, kind="file", fmt="fastq", default=default, label=label, help=help)


def bed(name: str = "src", default: str = REGIONS_BED, label: str = "BED file",
        help: str = "") -> Input:
    return Input(name, kind="file", fmt="bed", default=default, label=label, help=help)


def gff(name: str = "src", default: str = ANNOT_GFF, label: str = "GFF3/GTF file",
         help: str = "") -> Input:
    return Input(name, kind="file", fmt="gff", default=default, label=label, help=help)


def vcf(name: str = "src", default: str = VARIANTS_VCF, label: str = "VCF file",
        help: str = "") -> Input:
    return Input(name, kind="file", fmt="vcf", default=default, label=label, help=help)


def sam(name: str = "src", default: str = ALN_SAM, label: str = "SAM file",
        help: str = "") -> Input:
    return Input(name, kind="file", fmt="sam", default=default, label=label, help=help)


def tbl(name: str = "src", default: str = COUNTS, label: str = "Table (TSV/CSV)",
        help: str = "") -> Input:
    return Input(name, kind="file", fmt="tsv", default=default, label=label, help=help)


def txt(name: str = "src", default: str = "", label: str = "Text file", help: str = "") -> Input:
    return Input(name, kind="file", fmt="txt", default=default, label=label, help=help)


def img(name: str = "src", default: str = IMAGE, label: str = "Image (PGM/PPM)",
         help: str = "") -> Input:
    return Input(name, kind="file", fmt="pgm", default=default, label=label, help=help)


def seqin(name: str = "seq", default: str = SHORT_DNA, label: str = "Sequence",
          help: str = "") -> Input:
    return Input(name, kind="code", default=default, label=label, help=help)


def anyfile(name: str = "src", default: str = GENOME, fmt: str = "", label: str = "File",
            help: str = "") -> Input:
    return Input(name, kind="file", fmt=fmt, default=default, label=label, help=help)


def multi(name: str, options: list | tuple, default: list | None = None,
          label: str = "", help: str = "") -> Input:
    return Input(name, kind="multi", options=tuple(options), default=list(default or options),
                 label=label or name.replace("_", " ").title(), help=help)


def choice(name: str, options: list | tuple, default: Any = None, label: str = "",
           help: str = "") -> Input:
    return Input(name, kind="choice", options=tuple(options),
                 default=default if default is not None else options[0],
                 label=label or name.replace("_", " ").title(), help=help)


def number(name: str, default: float = 1.0, label: str = "", help: str = "",
           min: float | None = None, max: float | None = None, step: float | None = None) -> Input:
    return Input(name, kind="number", default=default, label=label or name.replace("_", " ").title(),
                 help=help, min=min, max=max, step=step)


def intin(name: str, default: int = 1, label: str = "", help: str = "",
          min: float | None = None, max: float | None = None) -> Input:
    return Input(name, kind="int", default=default, label=label or name.replace("_", " ").title(),
                 help=help, min=min, max=max, step=1)


def boolean(name: str, default: bool = False, label: str = "", help: str = "") -> Input:
    return Input(name, kind="bool", default=default, label=label or name.replace("_", " ").title(),
                 help=help)


def textbox(name: str, default: str = "", label: str = "", help: str = "") -> Input:
    return Input(name, kind="text", default=default, label=label or name.replace("_", " ").title(),
                 help=help)


def bigtext(name: str, default: str = "", label: str = "", help: str = "") -> Input:
    return Input(name, kind="code", default=default, label=label or name.replace("_", " ").title(),
                 help=help)


# ---------------------------------------------------------------------------
# result helpers
# ---------------------------------------------------------------------------
def table(rows: Any, message: str = "", name: str = "result.tsv",
          columns: list[str] | None = None) -> dict:
    df = tables.to_df(rows, columns) if not hasattr(rows, "to_csv") else rows
    if hasattr(rows, "to_csv") and columns and hasattr(df, "columns"):
        pass
    return io.table_result(df, message, name)


def figure(fig: Any, message: str = "", name: str = "figure.png") -> dict:
    return io.image_result(fig, message, name)


def text(body: str, message: str = "", name: str = "result.txt") -> dict:
    return {"text": body, "message": message, "filename": name} if message else \
        {"text": body, "filename": name}


def values_message(stats_dict: dict, message: str = "") -> dict:
    out: dict[str, Any] = {"stats": stats_dict}
    if message:
        out["message"] = message
    return out


def T(tool_id: str, name: str, section: str, output: str = "table",
      inputs: list[Input] | None = None, ex: dict | None = None, up: str = "",
      tags: tuple[str, ...] = (), summary: str = "", version: str = "1.0.0") -> Callable:
    """Short alias for :func:`chroma_titan.registry.tool`."""
    return tool(id=tool_id, name=name, section=section, output=output,
                inputs=tuple(inputs or ()), example=ex or {}, upstream=up, tags=tags,
                summary=summary, version=version)


# Common column orders
BED_COLS = ["chrom", "start", "end", "name", "score", "strand"]
GFF_COLS = ["seqid", "source", "type", "start", "end", "score", "strand", "phase"]
COLS = {"bed": BED_COLS}

__all__ = [
    "T", "Input", "Tool", "io", "genome", "stats", "tables",
    "fa", "fq", "bed", "gff", "vcf", "sam", "tbl", "txt", "img", "seqin", "anyfile",
    "multi", "choice", "number", "intin", "boolean", "textbox", "bigtext",
    "table", "figure", "text", "values_message",
    "GENOME", "GENOME2", "GENES_FA", "REGIONS_BED", "TARGETS_BED", "ANNOT_GFF", "VARIANTS_VCF",
    "READS_1", "READS_2", "READS_SINGLE", "ALN_SAM", "COUNTS", "GMT", "PFM", "JASPAR",
    "PROTEINS", "REF_PROTEINS", "MSA", "TREE", "PAIRS", "BEDGRAPH", "IMAGE", "PEPTIDES",
    "PHENO", "OBO", "STR", "CXCG", "GENOME_TXT", "TAXMAP", "SHORT_DNA", "SHORT_PROT",
    "SHORT_REGION", "BED_COLS", "GFF_COLS",
]
