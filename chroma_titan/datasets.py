"""The bundled example datasets (``data/examples``) as a browsable catalogue.

Every tool in the registry points at one of these files by default, which is what
makes the web app usable without uploading anything -- exactly the role the
"Send data to Galaxy" / bundled-history datasets play on usegalaxy.org.

>>> from chroma_titan import datasets
>>> datasets.count()
30
>>> datasets.load("counts.tsv").shape
(12, 7)
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

from chroma_titan.core import io

EXAMPLE_DIR: Path = io.EXAMPLE_DIR

#: name -> (format, panel section it feeds, one line description)
CATALOGUE: dict[str, tuple[str, str, str]] = {
    "genome.fa": ("FASTA", "fasta_fastq", "Small diploid genome (5 chromosomes) used as the reference"),
    "genome2.fa": ("FASTA", "fasta_fastq", "Second assembly of the same genome, for assembly-vs-reference comparisons"),
    "genome.txt": ("text", "assembly", "Plain-text contig listing used by assembly and scaffolding tools"),
    "genes.fa": ("FASTA", "fasta_fastq", "CDS sequences with descriptions, for ORF and codon tools"),
    "proteins.faa": ("FASTA", "proteomics", "Six protein chains with signal peptides, TM helices and motifs"),
    "reference_proteins.faa": ("FASTA", "annotation", "Reference proteome used by search/identity tools"),
    "reads_single.fastq": ("FASTQ", "fastq_quality_control", "Single-end reads with quality strings"),
    "reads_1.fastq": ("FASTQ", "fastq_quality_control", "Paired-end mates (forward)"),
    "reads_2.fastq": ("FASTQ", "fastq_quality_control", "Paired-end mates (reverse)"),
    "alignments.sam": ("SAM", "sam_bam", "Text SAM alignment of the reads onto genome.fa"),
    "regions.bed": ("BED", "bed", "Five genomic intervals with scores and strands"),
    "targets.bed": ("BED", "bed", "Candidate targets, used for intersect/subtract examples"),
    "two_bed.bed": ("BED", "bed", "A second interval set for set-operation demos"),
    "annotation.gff": ("GFF3", "annotation", "Gene/mRNA/exon annotation of genome.fa"),
    "variants.vcf": ("VCF", "vcf_bcf", "Short variants with FORMAT fields and two samples"),
    "counts.tsv": ("TSV", "rna_seq", "RNA-seq counts: 12 genes x 6 samples"),
    "phenotypes.tsv": ("TSV", "phenotype_association", "Sample sheet: group, treatment, yield, resistance, infection"),
    "genesets.gmt": ("GMT", "annotation", "Three gene sets for enrichment tools"),
    "phenotype.obo": ("OBO", "annotation", "Small ontology of phenotype terms"),
    "peptides.tsv": ("TSV", "proteomics", "Mass-spectrometry peptide table with intensities and runtimes"),
    "str_profile.tsv": ("TSV", "str_fm__microsatellite_analysis", "Microsatellite allele calls per locus and sample"),
    "taxmap.tsv": ("TSV", "metagenomic_analysis", "Accession to rank-7 lineage map ( SILVA style )"),
    "methylation.cxcg": ("cxCG", "epigenetics", "CpG/CHG/CHH methylation counts per cytosine"),
    "contacts.pairs": ("Hi-C pairs", "chromosome_conformation", "Restriction-site Hi-C contacts in pairs format"),
    "coverage.bedgraph": ("bedGraph", "deeptools", "Signal track (coverage) over the genome"),
    "motif.pfm": ("PFM", "motif", "Position frequency matrix in MEME/JASPAR style"),
    "motifs_jaspar.txt": ("JASPAR", "motif", "Two JASPAR-format matrices for comparison tools"),
    "msa.fasta": ("FASTA", "multiple_alignments", "Five aligned sequences for tree and MSA tools"),
    "tree.nwk": ("Newick", "phylogenetics", "Four-taxon rooted tree with branch lengths"),
    "cells.pgm": ("PGM", "imaging", "48x48 fluorescence-style microscopy image with four cells"),
}


def names() -> list[str]:
    """Every dataset known to the catalogue, in catalogue order."""
    return list(CATALOGUE)


def on_disk() -> list[str]:
    """Files actually present in ``data/examples``."""
    return sorted(p.name for p in EXAMPLE_DIR.glob("*") if p.is_file())


def count() -> int:
    return len(CATALOGUE)


def path(name: str) -> str:
    """Absolute path of a dataset (accepts bare names and ``examples/x`` forms)."""
    return io.resolve(_clean(name))


def size(name: str) -> int:
    return Path(path(name)).stat().st_size


def format_of(name: str) -> str:
    fmt, _sec, _desc = _info(name)
    return fmt


def section_for(name: str) -> str:
    _fmt, sec, _desc = _info(name)
    return sec


def description(name: str) -> str:
    _fmt, _sec, desc = _info(name)
    return desc


def entries() -> list[dict[str, Any]]:
    """Catalogue rows with paths and sizes, for the web app's dataset browser."""
    out = []
    for nm, (fmt, sec, desc) in CATALOGUE.items():
        p = EXAMPLE_DIR / nm
        out.append({"name": nm, "format": fmt, "section": sec, "description": desc,
                    "path": f"examples/{nm}", "bytes": p.stat().st_size if p.exists() else 0,
                    "present": p.exists()})
    return out


def table_rows() -> list[dict[str, Any]]:
    """``entries()`` with the number of lines/records, ready for a dataframe."""
    rows = []
    for e in entries():
        rec = dict(e)
        try:
            rec["lines"] = len(io.as_text(_clean(e["name"])).splitlines())
        except OSError:  # pragma: no cover - missing fixture
            rec["lines"] = 0
        rec["bytes"] = f"{e['bytes']:,}" if e["bytes"] else "-"
        rec.pop("present", None)
        rows.append(rec)
    return rows


def preview(name: str, rows: int = 8) -> str:
    """First ``rows`` lines, wrapped so they survive a browser window."""
    text = io.as_text(_clean(name))
    keep = [ln if len(ln) <= 160 else ln[:157] + "..." for ln in text.splitlines()[: int(rows)]]
    if len(text.splitlines()) > int(rows):
        keep.append(f"... {len(text.splitlines()) - len(keep)} more lines")
    return "\n".join(keep)


def load(name: str):
    """Parse a dataset into the natural Python object for its format."""
    src = _clean(name)
    low = name.lower()
    if low.endswith((".tsv", ".csv")):
        return io.read_table(src)
    if low.endswith((".fa", ".fasta", ".faa")):
        recs = io.parse_fasta(io.as_text(src))
        import pandas as pd

        return pd.DataFrame([{"id": r.id, "description": r.desc, "length": len(r.seq),
                             "sequence": r.seq} for r in recs])
    if low.endswith(".fastq"):
        recs = io.parse_fastq(io.as_text(src))
        import pandas as pd

        return pd.DataFrame([{"id": r.id, "length": len(r.seq), "mean_quality":
                             round(sum(ord(c) - 33 for c in r.qual) / max(1, len(r.qual)), 2)}
                            for r in recs])
    if low.endswith(".bed"):
        import pandas as pd

        ivs = io.parse_bed(io.as_text(src))
        return pd.DataFrame([{"chrom": i.chrom, "start": i.start, "end": i.end, "name": i.name,
                             "score": i.score, "strand": i.strand} for i in ivs])
    if low.endswith(".vcf"):
        import pandas as pd

        vcf = io.parse_vcf(io.as_text(src))
        return pd.DataFrame([{"chrom": v.chrom, "pos": v.pos, "id": v.id, "ref": v.ref,
                             "alt": v.alt, "qual": v.qual, "filter": v.filter} for v in vcf.records])
    if low.endswith(".gff"):
        import pandas as pd

        rows = io.parse_gff(io.as_text(src))
        return pd.DataFrame([{"seqid": g.seqid, "source": g.source, "type": g.type, "start": g.start,
                             "end": g.end, "score": g.score, "strand": g.strand} for g in rows])
    if low.endswith(".sam"):
        import pandas as pd

        return pd.DataFrame(io.parse_sam(io.as_text(src))[:200])
    if low.endswith((".pgm", ".ppm")):
        from chroma_titan.core import image

        return image.read_image(src)
    if low.endswith((".nwk", ".newick", ".tree")):
        return io.as_text(src).strip()
    return io.as_text(src)


def tools_using(name: str, limit: int = 40) -> list[str]:
    """Tool ids whose bundled example points at this dataset."""
    from chroma_titan import registry

    base = _clean(name).rsplit("/", 1)[-1]
    out = []
    for t in registry.all_tools():
        for _k, v in t.example.items():
            if isinstance(v, str) and base in v:
                out.append(t.id)
                break
    return out[:limit]


def _clean(name: str) -> str:
    nm = str(name).strip()
    if nm.startswith("examples/"):
        return nm
    return f"examples/{nm}"


def _info(name: str) -> tuple[str, str, str]:
    nm = _clean(name).rsplit("/", 1)[-1]
    if nm in CATALOGUE:
        return CATALOGUE[nm]
    suffix = nm.rsplit(".", 1)[-1].upper() if "." in nm else "text"
    return suffix, "get_data", f"Unlisted example file ({nm})"


__all__ = ["CATALOGUE", "EXAMPLE_DIR", "names", "on_disk", "count", "path", "size", "format_of",
           "section_for", "description", "entries", "table_rows", "preview", "load", "tools_using"]
