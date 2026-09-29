"""Genome annotation, sequence motifs, molecular evolution, CRISPR design, ChIP-seq.

Panel sections covered: **annotation**, **motif**, **evolution**,
**genome_editing**, **chip_seq**, **regional_variation**, **fetch_sequences_alignments**.
"""
from __future__ import annotations

import math
import re
from collections import Counter, defaultdict

from chroma_titan.core import align, bam, genome, io, motif, phylo, plot, protein, seq, stats, tables
from chroma_titan.tools._common import *  # noqa: F401,F403
from chroma_titan.tools._common import (ANNOT_GFF, BEDGRAPH, GENES_FA, GENOME, GMT, JASPAR, MSA,
                                         PFM, PROTEINS, READS_SINGLE, REF_PROTEINS, REGIONS_BED,
                                         SHORT_DNA, TARGETS_BED, TREE)

ANN = "annotation"
MOT = "motif"
EVO = "evolution"
EDIT = "genome_editing"
CHIP = "chip_seq"
REG = "regional_variation"
FETCH = "fetch_sequences_alignments"


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------
def _gff(src) -> list:
    return io.parse_gff(io.as_text(src))


def _rows(src) -> list[dict]:
    return [{"id": r.id, "seq": r.seq} for r in io.parse_fasta(io.as_text(src))]


def _bed(src):
    return io.parse_bed(io.as_text(src))


def _genome(src):
    seqs, sizes = genome.read_genome(src)
    return {k: "".join(v) for k, v in seqs.items()}, sizes


def _attr(row, key, default=""):
    v = (row.attrs or {}).get(key, default)
    return v if v not in (None, "") else default


def _score(x):
    return 0 if not isinstance(x, (int, float)) or math.isnan(x) else x


def _iupac_regex(pattern: str) -> str:
    try:
        return seq.iupac_to_regex(pattern)
    except Exception:  # noqa: BLE001
        return re.escape(pattern)


# ===========================================================================
# GFF3 / annotation tables
# ===========================================================================
@T("gff3_check", "Sanity-check a GFF3 file", ANN, "table",
   [gff("src", ANNOT_GFF, "GFF3 file"), choice("level", ["all", "error", "warning"], "all", "Report level")],
   ex={"src": ANNOT_GFF}, up="gff3read / qualifyGFF", tags=("GFF3", "quality control"),
   summary="Validate types, coordinates, IDs and Parent links, and report each problem found.")
def gff3_check(src, level="all"):
    """GFF3 sanity checks."""
    rows = _gff(src)
    ids = Counter(_attr(r, "ID") for r in rows if _attr(r, "ID"))
    parents = {k for k in (_attr(r, "Parent") for r in rows) if k}
    out = []
    for i, r in enumerate(rows):
        probs = []
        if r.end < r.start:
            probs.append(("error", "end_before_start"))
        if r.start < 1:
            probs.append(("error", "zero_based_start"))
        if r.type == "gene" and not _attr(r, "ID"):
            probs.append(("error", "gene_without_ID"))
        if r.type in ("mRNA", "CDS", "exon") and not parents and not _attr(r, "Parent"):
            probs.append(("warning", "no_parent_link"))
        for kind, msg in probs:
            out.append({"line": i + 1, "severity": kind, "problem": msg, "seqid": r.seqid,
                       "type": r.type, "start": r.start, "end": r.end,
                       "id": _attr(r, "ID") or _attr(r, "Name") or "."})
    if len([k for k, v in ids.items() if v > 1 and k]) > 0:
        for k, v in ids.items():
            if v > 1:
                out.append({"line": 0, "severity": "warning", "problem": "duplicate_ID",
                           "seqid": ".", "type": ".", "start": 0, "end": 0, "id": f"{k} x{v}"})
    keep = out if level == "all" else [o for o in out if o["severity"] == level]
    return table(keep, f"{len(keep)} findings over {len(rows)} GFF lines "
                       f"({sum(1 for o in out if o['severity'] == 'error')} errors)")


@T("gff3_feature_counts", "Count features by type and sequence", ANN, "table",
   [gff("src", ANNOT_GFF, "GFF3 file"), choice("group_by", ["type", "seqid", "biotype", "source"],
                                               "type", "Group by")],
   ex={"src": ANNOT_GFF, "group_by": "type"}, up="gff3 read / featureCounts summary",
   tags=("GFF3", "counts"),
   summary="Feature tallies per group with total and mean feature length.")
def gff3_feature_counts(src, group_by="type"):
    """Group GFF rows and count them."""
    rows = _gff(src)
    groups: dict[str, list] = defaultdict(list)
    for r in rows:
        key = {"type": r.type, "seqid": r.seqid, "source": r.source,
               "biotype": _attr(r, "gene_biotype") or _attr(r, "biotype") or "unknown"}[group_by]
        groups[str(key)].append(r)
    out = []
    for k, v in sorted(groups.items(), key=lambda kv: -len(kv[1])):
        lens = [x.end - x.start + 1 for x in v]
        out.append({"group": k, "features": len(v), "total_bp": sum(lens),
                    "mean_length": round(stats.mean(lens), 2),
                    "types": ",".join(sorted({x.type for x in v})[:6]),
                    "seqids": ",".join(sorted({x.seqid for x in v})[:6])})
    return table(out, f"{len(rows)} rows in {len(groups)} {group_by} groups")


@T("gff3_filter_attributes", "Filter GFF rows by attributes", ANN, "file",
   [gff("src", ANNOT_GFF, "GFF3 file"), textbox("attribute", "Name", "Attribute name"),
    textbox("value", "", "Value or regex"), choice("type", ["any", "gene", "mRNA", "exon", "CDS"],
                                                  "any", "Feature type"),
    boolean("invert", False, "Keep non-matching rows")],
   ex={"src": ANNOT_GFF, "attribute": "gene_biotype", "value": "protein_coding", "type": "gene"},
   up="gff3 read / AGAT filter", tags=("GFF3", "filtering"),
   summary="Keep rows whose attribute matches a value or regular expression.")
def gff3_filter_attributes(src, attribute="Name", value="", type="any", invert=False):
    """Attribute filtering of a GFF3."""
    rows = _gff(src)
    rx = re.compile(_iupac_regex(value) if value else ".*")
    keep = []
    for r in rows:
        if type != "any" and r.type != type:
            continue
        av = str(_attr(r, attribute or "Name"))
        hit = bool(rx.search(av)) if value else bool(av)
        if hit != bool(invert):
            keep.append(r)
    body = io.write_gff(keep)
    return {"text": body, "filename": "filtered.gff3",
            "message": f"{len(keep)}/{len(rows)} rows kept"}


@T("gff3_to_bed", "Convert GFF3 to BED (BED12 for transcripts)", ANN, "file",
   [gff("src", ANNOT_GFF, "GFF3 file"), choice("type", ["any", "gene", "mRNA", "exon", "CDS"],
                                               "gene", "Feature type"),
    boolean("bed12", True, "Write BED12 blocks for transcripts")],
   ex={"src": ANNOT_GFF, "type": "gene", "bed12": False}, up="gff3 to bed (bed12plus)",
   tags=("GFF3", "BED", "conversion"),
   summary="Rewrite annotation rows as genomic intervals, keeping exon blocks in BED12.")
def gff3_to_bed(src, type="gene", bed12=True):
    """GFF3 -> BED."""
    rows = _gff(src)
    keep = [r for r in rows if type == "any" or r.type == type]
    if bed12 and type in ("mRNA", "transcript"):
        ivs = io.gff_to_bed12(keep)
    else:
        ivs = genome.gff_to_intervals(keep)
    body = io.write_bed(genome.sort_ivs(ivs), bed12=bool(bed12 and type == "mRNA"))
    return {"text": body, "filename": "annotation.bed", "message": f"{len(ivs)} intervals"}


@T("gff3_gtf_convert", "Convert between GFF3 and GTF", ANN, "file",
   [gff("src", ANNOT_GFF, "GFF/GTF file"), choice("to", ["gtf", "gff3"], "gtf", "Target format")],
   ex={"src": ANNOT_GFF, "to": "gtf"}, up="gffread / AGAT convert",
   tags=("GFF3", "GTF", "conversion"),
   summary="Swap the attribute syntax between RefSeq GTF and Ensembl-style GFF3.")
def gff3_gtf_convert(src, to="gtf"):
    """GFF3 <-> GTF attribute conversion."""
    rows = _gff(src)
    out = []
    for r in rows:
        attrs = dict(r.attrs or {})
        if to == "gtf":
            gene = attrs.pop("gene", None) or attrs.get("Name", ".")
            tid = attrs.pop("ID", None) or "."
            parent = attrs.pop("Parent", None) or gene
            text = (f'gene_id "{gene}"; transcript_id "{tid if r.type != "gene" else gene}"; '
                    f'statistics "{parent}"')
            extra = "; ".join(f'{k} "{v}"' for k, v in attrs.items() if k not in ("Name",))
            line = [r.seqid, r.source, r.type, str(r.start), str(r.end),
                    "." if not isinstance(r.score, float) or math.isnan(r.score) else f"{r.score:g}",
                    r.strand, r.phase or ".", "; ".join(x for x in [text, extra] if x) + ";"]
        else:
            gene = next((v for k, v in attrs.items() if k in ("gene_id", "gene")), ".")
            tid = next((v for k, v in attrs.items() if k in ("transcript_id",)), "")
            new = {"ID": tid or gene, "Parent": gene, "Name": gene}
            new.update({k: v for k, v in attrs.items() if k not in ("gene_id", "transcript_id")})
            r.attrs = new
            line = None
        if line is not None:
            out.append("\t".join(line))
    body = "\n".join(out) + "\n" if to == "gtf" else io.write_gff(rows)
    return {"text": body, "filename": f"annotation.{to}", "message": f"{len(rows)} rows as {to}"}


@T("gff3_infer_transcripts", "Assemble transcripts from exon rows", ANN, "table",
   [gff("src", ANNOT_GFF, "GFF3 file"), intin("min_exons", 1, "Minimum exons per transcript", min=1)],
   ex={"src": ANNOT_GFF}, up="gff3read --exonmerge / BRAKER", tags=("GFF3", "transcripts"),
   summary="Group exon rows by Parent and rebuild transcript models with CDS length.")
def gff3_infer_transcripts(src, min_exons=1):
    """Reconstruct transcript models from exons."""
    rows = _gff(src)
    exons: dict[str, list] = defaultdict(list)
    for r in rows:
        if r.type in ("exon", "CDS"):
            key = _attr(r, "Parent") or _attr(r, "transcript_id") or _attr(r, "ID") or r.seqid
            exons[str(key)].append(r)
    out = []
    for tid, ex in exons.items():
        if len(ex) < int(min_exons):
            continue
        ex = sorted(ex, key=lambda r: r.start)
        length = sum(r.end - r.start + 1 for r in ex)
        span = ex[-1].end - ex[0].start
        out.append({"transcript": tid, "seqid": ex[0].seqid, "exons": len(ex),
                    "start": ex[0].start, "end": ex[-1].end, "spliced_length": length,
                    "genomic_span": span, "intron_count": max(0, len(ex) - 1),
                    "mean_intron": round((span - length) / max(1, len(ex) - 1), 1),
                    "strand": ex[0].strand, "coding_density": round(length / max(1, span), 4)})
    return table(out, f"{len(out)} transcript models inferred")


@T("gff3_children_table", "Parent-child feature table", ANN, "table",
   [gff("src", ANNOT_GFF, "GFF3 file"), choice("from_type", ["gene", "transcript", "mRNA"], "gene",
                                              "Parent type"), intin("head", 300, "Rows", min=1)],
   ex={"src": ANNOT_GFF}, up="AGAT extract / gff3read children",
   tags=("GFF3", "hierarchy"),
   summary="Flatten the gene → transcript → exon/CDS hierarchy into one row per child.")
def gff3_children_table(src, from_type="gene", head=300):
    """Gene/transcript/child expansion."""
    rows = _gff(src)
    by_id = {}
    for r in rows:
        if _attr(r, "ID"):
            by_id[_attr(r, "ID")] = r
    out = []
    wanted = {from_type, {"gene": "mRNA", "mRNA": "exon", "transcript": "CDS"}[from_type]}
    for r in rows:
        parent_id = _attr(r, "Parent")
        if r.type not in wanted and r.type != from_type:
            continue
        parent = by_id.get(parent_id) if parent_id else None
        grand = by_id.get(_attr(parent, "Parent")) if parent is not None else None
        gene_row = parent if (parent is not None and parent.type == "gene") else (grand or parent)
        out.append({"parent": _attr(gene_row, "Name") or _attr(gene_row, "ID") if gene_row else ".",
                    "parent_type": from_type, "child": r.type, "child_id": _attr(r, "ID") or ".",
                    "seqid": r.seqid, "start": r.start, "end": r.end, "strand": r.strand,
                    "length": r.end - r.start + 1})
    return table(out[:int(head)], f"{len(out)} parent-child links")


@T("gff3_loci_to_fasta", "Extract sequences for annotated features", FETCH, "file",
   [gff("src", ANNOT_GFF, "GFF3 file"), fa("ref", GENOME, "Genome FASTA"),
    choice("type", ["gene", "mRNA", "exon", "CDS"], "gene", "Feature type"),
    choice("strand", ["same", "plus", "minus"], "same", "Strand handling"),
    intin("extend", 0, "Extend each side (bp)", min=0)],
   ex={"src": ANNOT_GFF, "ref": GENOME, "type": "gene"}, up="bedtools getfasta / gff3read --loci",
   tags=("FASTA", "extraction", "GFF3"),
   summary="Pull the genomic sequence of every feature, optionally strand-aware with flanks.")
def gff3_loci_to_fasta(src, ref, type="gene", strand="same", extend=0):
    """getfasta for annotation features."""
    gseq, sizes = _genome(ref)
    ivs = [r for r in _gff(src) if type == "any" or r.type == type]
    intervals = []
    for r in ivs:
        s = max(1, r.start - int(extend))
        e = r.end + int(extend)
        intervals.append(io.Interval(r.seqid, s - 1, e, _attr(r, "Name") or _attr(r, "ID") or ".",
                                    _score(r.score), r.strand))
    got = genome.getfasta(gseq, intervals, strand=strand, name_mode="bed")
    lines = []
    for iv in got:
        name = f"{iv.chrom}:{iv.start + 1}-{iv.end}({iv.strand})_{iv.name}"
        body = iv.extra[0] if iv.extra else ""
        lines.append(f">{name}\n{body}")
    return {"text": "\n".join(lines) + "\n", "filename": "loci.fasta",
            "message": f"{len(got)} sequences extracted"}


@T("gff3_to_genbank", "Write a GenBank-style flat file", ANN, "file",
   [gff("src", ANNOT_GFF, "GFF3 file"), fa("ref", GENOME, "Genome FASTA"),
    textbox("definition", "Chroma-Titan synthetic record", "Definition line")],
   ex={"src": ANNOT_GFF, "ref": GENOME}, up="gff3 to genbank (BP)",
   tags=("GenBank", "conversion", "GFF3"),
   summary="Emit LOCUS/FEATURES records with CDS translations from a GFF3 and genome.")
def gff3_to_genbank(src, ref, definition="Chroma-Titan synthetic record"):
    """Very small GFF3 -> GenBank writer."""
    gseq, sizes = _genome(ref)
    rows = _gff(src)
    per_seq: dict[str, list] = defaultdict(list)
    for r in rows:
        per_seq[r.seqid].append(r)
    out = []
    for name, rs in per_seq.items():
        length = sizes.get(name, max((r.end for r in rs), default=0))
        seq_text = gseq.get(name, "")
        out.append(f"LOCUS       {name:<16} {length} bp    DNA   linear  SYN\n"
                   f"DEFINITION  {definition}.\nACCESSION   {name}\nVERSION     {name}.1\n"
                   f"FEATURES             Location/Qualifiers\n     source          1..{length}\n"
                   f"                     /mol_type=\"genomic DNA\"\n")
        for r in sorted((x for x in rs if x.type in ("gene", "CDS")), key=lambda x: x.start):
            gene = _attr(r, "Name") or _attr(r, "gene") or "unnamed"
            loc = f"{r.start}..{r.end}" + ("complement(" if r.strand == "-" else "")
            if r.strand == "-":
                loc = f"complement({r.start}..{r.end})"
            prot = ""
            if r.type == "CDS" and seq_text:
                sub = seq_text[r.start - 1:r.end]
                prot = seq.translate(sub if r.strand != "-" else seq.reverse_complement(sub),
                                     frame=0, start_m="M")
            quals = f'                     /gene="{gene}"\n'
            if prot:
                quals += f'                     /translation="{prot.rstrip("*")}"\n'
            out.append(f"     {r.type:<16} {loc}\n{quals}")
        out.append("//\n")
    body = "".join(out)
    return {"text": body, "filename": "annotation.gb",
            "message": f"{len(per_seq)} GenBank records with {len(rows)} features"}


@T("gff3_summary", "Annotation coverage statistics", ANN, "table",
   [gff("src", ANNOT_GFF, "GFF3 file"), anyfile("sizes", GENOME, fmt="", label="Genome (FASTA or genome.txt)")],
   ex={"src": ANNOT_GFF, "sizes": GENOME}, up="qualifyGFF /assembly summary",
   tags=("GFF3", "statistics"),
   summary="Number of genes, transcripts, exons, mean lengths and genome coverage.")
def gff3_summary(src, sizes):
    """Whole-annotation summary."""
    rows = _gff(src)
    _, gsize = _genome(sizes)
    total = sum(gsize.values()) or 0
    genes = genome.genes_from_gff(rows)
    exons = [r for r in rows if r.type == "exon"]
    covered = sum(iv.end - iv.start for iv in genome.merge(genes))
    by_type = Counter(r.type for r in rows)
    out = [{"metric": k, "value": v} for k, v in {
        "gff_rows": len(rows), "types": len(by_type), "genes": len(genes),
        "exons": len(exons), "transcripts": by_type.get("mRNA", 0) + by_type.get("transcript", 0),
        "sequences_annotated": len({r.seqid for r in rows}),
        "genome_length": total,
        "gene_coverage_percent": round(100 * covered / total, 3) if total else 0.0,
        "mean_gene_length": round(stats.mean([iv.end - iv.start for iv in genes]), 2) if genes else 0,
        "median_gene_length": stats.median([iv.end - iv.start for iv in genes]) if genes else 0}.items()]
    res = table(out, f"{len(genes)} genes over {len({r.seqid for r in rows})} sequences")
    res["stats"]["features_per_type"] = dict(by_type)
    return res


@T("gff3_tss_bed", "Transcription start sites of all genes", ANN, "file",
   [gff("src", ANNOT_GFF, "GFF3 file"), intin("window", 0, "TSS window width", min=0),
    boolean("upstream", False, "Shift upstream instead of centring")],
   ex={"src": ANNOT_GFF, "window": 100}, up="tss from GFF3 (bedtools)",
   tags=("BED", "promoter", "TSS"),
   summary="Reduce a GFF3 to per-gene TSS intervals, optionally as promoter windows.")
def gff3_tss_bed(src, window=0, upstream=False):
    """BED of transcription start sites."""
    tss = genome.tss_list(_gff(src))
    w = int(window)
    out = []
    for iv in tss:
        if upstream and w > 0:
            shift = -w if iv.strand != "-" else 0
            out.append(io.Interval(iv.chrom, iv.start + shift, iv.start + shift + w, iv.name,
                                  _score(iv.score), iv.strand))
        else:
            s = max(0, iv.start - w // 2)
            out.append(io.Interval(iv.chrom, s, iv.end + w // 2, iv.name, _score(iv.score), iv.strand))
    body = io.write_bed(genome.sort_ivs(out), bed12=False)
    return {"text": body, "filename": "tss.bed", "message": f"{len(out)} TSS intervals"}


@T("gff3_nearest_gene", "Annotate intervals with the nearest gene", ANN, "table",
   [bed("src", REGIONS_BED, "Query intervals"), gff("annotation", ANNOT_GFF, "GFF3 annotation"),
    number("max_distance", 1e12, "Maximum distance (bp)", min=0.0),
    choice("tie", ["closest", "first", "all"], "closest", "Ties")],
   ex={"src": REGIONS_BED, "annotation": ANNOT_GFF}, up="bedtools closest",
   tags=("intervals", "annotation", "nearest"),
   summary="Attach the closest gene name, overlap status and distance to every interval.")
def gff3_nearest_gene(src, annotation, max_distance=1e12, tie="closest"):
    """Nearest-gene annotation."""
    a = _bed(src)
    genes = genome.genes_from_gff(_gff(annotation))
    rows = genome.closest(a, genes)
    out = []
    for r in rows:
        dist = r.get("distance", 0) or 0
        if abs(float(dist)) > float(max_distance):
            continue
        out.append({"chrom": r.get("chrom"), "start": r.get("start"), "end": r.get("end"),
                    "name": r.get("name", "."), "gene": r.get("other_name", ""),
                    "gene_chrom": r.get("other_chrom", ""), "gene_start": r.get("other_start", 0),
                    "gene_end": r.get("other_end", 0), "distance": int(dist),
                    "overlaps": bool(dist <= 0)})
    return table(out, f"{len(out)} intervals annotated")


@T("gff3_compare", "Compare two annotation sets", ANN, "table",
   [gff("a", ANNOT_GFF, "Reference GFF3"), gff("b", ANNOT_GFF, "Query GFF3"),
    choice("key", ["name", "coordinates", "id"], "name", "Match by")],
   ex={"a": ANNOT_GFF, "b": ANNOT_GFF, "key": "name"}, up="gff3compare / compare GFF",
   tags=("GFF3", "comparison"),
   summary="Report genes shared by two annotations and the positional Jaccard of all features.")
def gff3_compare(a, b, key="name"):
    """Annotation comparison."""
    ra, rb = _gff(a), _gff(b)

    def feat_key(r):
        if key == "coordinates":
            return f"{r.seqid}:{r.start}-{r.end}:{r.type}"
        if key == "id":
            return _attr(r, "ID") or _attr(r, "Name") or f"{r.seqid}:{r.start}"
        return _attr(r, "Name") or _attr(r, "ID") or f"{r.seqid}:{r.start}"

    ga = {feat_key(r) for r in ra if r.type == "gene"} or {feat_key(r) for r in ra}
    gb = {feat_key(r) for r in rb if r.type == "gene"} or {feat_key(r) for r in rb}
    jac = genome.jaccard(genome.gff_to_intervals(ra), genome.gff_to_intervals(rb))
    rows = [{"status": "shared", "feature": k} for k in sorted(ga & gb)]
    rows += [{"status": "reference_only", "feature": k} for k in sorted(ga - gb)]
    rows += [{"status": "query_only", "feature": k} for k in sorted(gb - ga)]
    res = table(rows, f"{len(ga & gb)} shared, {len(ga - gb)} only in A, {len(gb - ga)} only in B")
    res["stats"] = {k: _score(v) if isinstance(v, float) else v for k, v in jac.items()}
    res["stats"]["jaccard"] = round(float(jac.get("jaccard", 0.0) or 0.0), 5)
    return res


@T("gff3_add_locus_tags", "Renumber feature IDs with a locus tag prefix", ANN, "file",
   [gff("src", ANNOT_GFF, "GFF3 file"), textbox("prefix", "CHROMA_", "Locus tag prefix"),
    intin("start", 1, "First number", min=0), intin("padding", 5, "Zero padding", min=1),
    choice("type", ["gene", "CDS", "mRNA", "any"], "gene", "Tag which features")],
   ex={"src": ANNOT_GFF, "prefix": "VIRUS_"}, up="Producer of Salk annotations",
   tags=("GFF3", "locus tag"),
   summary="Assign sequential locus tags and keep Parent links consistent.")
def gff3_add_locus_tags(src, prefix="LOCUS_", start=1, padding=5, type="gene"):
    """Locus tag assignment."""
    rows = _gff(src)
    n = int(start)
    rename: dict[str, str] = {}
    for r in rows:
        old = _attr(r, "ID")
        if type in ("any", r.type) and old:
            new = f"{prefix}{n:0{int(padding)}d}"
            rename[old] = new
            r.attrs = dict(r.attrs or {})
            r.attrs["ID"] = new
            if not _attr(r, "locus_tag"):
                r.attrs["locus_tag"] = new
            n += 1
    for r in rows:
        p = _attr(r, "Parent")
        if p in rename:
            r.attrs = dict(r.attrs or {})
            r.attrs["Parent"] = rename[p]
    body = io.write_gff(rows)
    return {"text": body, "filename": "tagged.gff3", "message": f"{len(rename)} features re-tagged"}


@T("gff3_biotype_table", "Gene biotype census", ANN, "table",
   [gff("src", ANNOT_GFF, "GFF3 file"), boolean("lengths", True, "Include lengths")],
   ex={"src": ANNOT_GFF}, up="biotype table (annotation)", tags=("GFF3", "biotype"),
   summary="Genes per biotype with total length, mean length and chromosome distribution.")
def gff3_biotype_table(src, lengths=True):
    """Biotype summary."""
    rows = [r for r in _gff(src) if r.type == "gene"] or _gff(src)
    groups: dict[str, list] = defaultdict(list)
    for r in rows:
        bt = _attr(r, "gene_biotype") or _attr(r, "biotype") or "unknown"
        groups[str(bt)].append(r)
    out = []
    for bt, v in sorted(groups.items(), key=lambda kv: -len(kv[1])):
        lens = [r.end - r.start + 1 for r in v]
        rec = {"biotype": bt, "genes": len(v), "chromosomes": len({r.seqid for r in v}),
              "names": ",".join(sorted({str(_attr(r, "Name") or _attr(r, "ID")) for r in v})[:8])}
        if lengths:
            rec["total_bp"] = sum(lens)
            rec["mean_bp"] = round(stats.mean(lens), 1)
        out.append(rec)
    return table(out, f"{len(rows)} genes across {len(groups)} biotypes")


@T("gff3_cds_translation", "Translate annotated CDS to proteins", ANN, "file",
   [gff("src", ANNOT_GFF, "GFF3 file"), fa("ref", GENOME, "Genome FASTA"),
    boolean("require_start", True, "Rewrite the first codon to M"),
    intin("min_length", 20, "Minimum protein length", min=1)],
   ex={"src": ANNOT_GFF, "ref": GENOME}, up="gff3read --prot / prodigal-like",
   tags=("translation", "proteins", "GFF3"),
   summary="Splice CDS exons per transcript, translate, and report length and stops.")
def gff3_cds_translation(src, ref, require_start=True, min_length=20):
    """GFF3 CDS -> protein FASTA."""
    gseq, _sizes = _genome(ref)
    rows = _gff(src)
    cds: dict[str, list] = defaultdict(list)
    for r in rows:
        if r.type == "CDS":
            key = _attr(r, "Parent") or _attr(r, "ID") or _attr(r, "gene") or r.seqid
            cds[str(key)].append(r)
    entries = []
    dropped = 0
    for tid, rs in cds.items():
        rs = sorted(rs, key=lambda r: r.start)
        piece = "".join(gseq.get(r.seqid, "")[r.start - 1:r.end] for r in rs)
        if rs and rs[0].strand == "-":
            piece = seq.reverse_complement(piece)
        prot = seq.translate(piece, frame=0, to_stop=False)
        prot = prot.rstrip("*")
        if require_start and prot and prot[0] != "M":
            prot = "M" + prot[1:]
        if len(prot) < int(min_length):
            dropped += 1
            continue
        entries.append((tid, prot))
    body = io.write_fasta(entries)
    return {"text": body, "filename": "proteins.faa",
            "message": f"{len(entries)} proteins translated ({dropped} too short or empty)"}


@T("genome_predict_orfs", "Gene prediction by open reading frames", ANN, "table",
   [fa("ref", GENES_FA, "Sequence FASTA"), intin("min_length", 90, "Minimum ORF length (nt)", min=3),
    choice("start", ["ATG", "GTG", "TTG"], "ATG", "Start codon"),
    boolean("partial", False, "Report ORFs without stop codon"), intin("head", 200, "Rows", min=1)],
   ex={"ref": GENES_FA, "min_length": 90}, up="prodigal / getorf (EMBOSS)",
   tags=("ORF", "gene prediction"),
   summary="Six-frame ORF finding with translation, GC and length statistics per ORF.")
def genome_predict_orfs(ref, min_length=90, start="ATG", partial=False, head=200):
    """ORF-based prediction."""
    rows = _rows(ref)
    out = []
    for r in rows:
        for o in seq.find_orfs(r["seq"], min_len=int(min_length), required_start=start,
                             partial=bool(partial)):
            out.append({"sequence": r["id"], "orf": o.get("id", ""), "strand": o.get("strand", "+"),
                        "frame": o.get("frame", 1), "start": o.get("start", 0), "end": o.get("end", 0),
                        "nt_length": o.get("length", 0), "aa_length": o.get("aa_length", 0),
                        "gc_percent": round(100 * seq.gc_content(r["seq"][o.get("start", 0):o.get("end", 0)]), 2),
                        "protein": (o.get("protein", "") or "")[:60]})
    return table(out[:int(head)], f"{len(out)} ORFs >= {min_length} nt in {len(rows)} sequences")


@T("annotate_genome_by_blast", "Transfer protein annotation to a genome", ANN, "table",
   [fa("proteins", PROTEINS, "Annotated proteins"), fa("ref", GENES_FA, "Target sequences"),
    intin("word_size", 3, "Word size", min=2), number("evalue", 10.0, "E-value", min=0.0),
    number("min_identity", 60.0, "Minimum identity %", min=0.0)],
   ex={"proteins": PROTEINS, "ref": REF_PROTEINS, "word_size": 3, "evalue": 1e-2},
   up="maker / blast2go annotation transfer", tags=("annotation", "BLAST", "transfer"),
   summary="Blast proteins against the target gene set and transfer names, IDs and products.")
def annotate_genome_by_blast(proteins, ref, word_size=3, evalue=10.0, min_identity=60.0):
    """Annotation transfer by homology."""
    hits = align.blast_lite(_rows(proteins), _rows(ref), word_size=int(word_size),
                           evalue=float(evalue), scoring="blastp", matrix="BLOSUM62", max_hits=20)
    best = align.best_hits(hits)
    by_target: dict[str, dict] = {}
    for h in best:
        if float(h["pident"]) >= float(min_identity):
            cur = by_target.get(h["saccver"])
            if cur is None or float(h["bitscore"]) > float(cur["bitscore"]):
                by_target[h["saccver"]] = h
    rows = []
    for t in _rows(ref):
        h = by_target.get(t["id"], {})
        rows.append({"target": t["id"], "target_length": len(t["seq"]),
                    "annotated": bool(h), "query": h.get("qaccver", ""),
                    "identity_percent": h.get("pident", 0.0), "evalue": h.get("evalue", 0.0),
                    "bitscore": h.get("bitscore", 0.0),
                    "product": h.get("saccver", "") and f"similar to {h.get('qaccver')}" or "hypothetical protein"})
    n_ok = sum(1 for r in rows if r["annotated"])
    return table(rows, f"{n_ok}/{len(rows)} targets annotated by homology")


@T("bedgraph_per_gene_signal", "Mean signal of a bedGraph inside features", ANN, "table",
   [gff("src", ANNOT_GFF, "GFF3 file"), anyfile("graph", BEDGRAPH, fmt="", label="bedGraph / WIG"),
    choice("stat", ["mean", "max", "sum", "median"], "mean", "Summary statistic"),
    fa("ref", GENOME, "Genome FASTA (for GC normalisation)"),
    boolean("gc_normalise", False, "Normalise by GC content"), intin("head", 200, "Rows", min=1)],
   ex={"src": ANNOT_GFF, "graph": BEDGRAPH, "stat": "mean"}, up="bedtools map / deepTools computeMatrix",
   tags=("coverage", "quantification"),
   summary="Aggregate a signal track over each gene, optionally GC-corrected.")
def bedgraph_per_gene_signal(src, graph, stat="mean", ref=GENOME, gc_normalise=False, head=200):
    """Per-gene signal from a bedGraph."""
    vals = io.parse_bedgraph(io.as_text(graph))
    tracks: dict[str, list] = defaultdict(list)
    for chrom, s, e, v in vals:
        tracks[chrom].append((s, e, float(v)))
    gseq, _sizes = _genome(ref) if gc_normalise else ({}, {})
    rows = []
    for r in _gff(src):
        if r.type not in ("gene", "CDS", "mRNA"):
            continue
        hit = [v for s, e, v in tracks.get(r.seqid, []) if s < r.end and e > r.start]
        if not hit:
            continue
        val = {"mean": stats.mean, "max": max, "median": stats.median, "sum": sum}[stat](hit)
        rec = {"chrom": r.seqid, "start": r.start, "end": r.end,
               "gene": _attr(r, "Name") or _attr(r, "ID") or ".", "strand": r.strand,
               stat + "_signal": round(float(val), 5), "n_blocks": len(hit)}
        if gc_normalise:
            s = gseq.get(r.seqid, "")[r.start - 1:r.end]
            gc = seq.gc_content(s) or 0.25
            rec["gc_percent"] = round(100 * gc, 2)
            rec["gc_corrected"] = round(float(val) / (gc if gc > 0.05 else 0.05), 5)
        rows.append(rec)
    rows.sort(key=lambda x: -abs(x[stat + "_signal"]))
    return table(rows[:int(head)], f"{len(rows)} genes with signal")


@T("predict_proteins_from_reads", "Assemble-like ORF calling from reads", ANN, "table",
   [fq("reads", READS_SINGLE, "Reads (FASTQ)"), intin("k", 15, "k-mer size for extension", min=5),
    intin("min_protein", 30, "Minimum protein length (aa)", min=5), intin("head", 100, "Rows", min=1)],
   ex={"reads": READS_SINGLE, "k": 15, "min_protein": 30}, up="metaWRAP / prodigal on contigs",
   tags=("ORF", "reads", "annotation"),
   summary="Concatenate reads into long contigs and call the longest protein per contig.")
def predict_proteins_from_reads(reads, k=15, min_protein=30, head=100):
    """Contig + ORF calling from a read set."""
    rs = [seq.clean(r.seq) for r in io.parse_fastq(io.as_text(reads))]
    contigs = [""]
    for s in rs:  # naive overlap extension by k-mer suffix match
        for i in range(len(s) - int(k)):
            tail = contigs[-1][-int(k):]
            if tail and s.startswith(tail):
                contigs[-1] += s[int(k):]
                break
        else:
            contigs.append(s)
    rows = []
    for i, c in enumerate([c for c in contigs if len(c) >= 3 * int(min_protein)]):
        best = seq.longest_orf_protein(c)
        if len(best) >= int(min_protein):
            rows.append({"contig": f"contig_{i + 1}", "length": len(c), "protein_length": len(best),
                        "mw_kda": round(protein.molecular_weight(best) / 1000.0, 2),
                        "gc_percent": round(100 * seq.gc_content(c), 2),
                        "sequence": best[:70] + ("..." if len(best) > 70 else "")})
    return table(rows[:int(head)], f"{len(rows)} putative proteins from {len(contigs)} contigs")


@T("protein_function_summary", "Domain and property census of a proteome", ANN, "table",
   [fa("proteins", PROTEINS, "Protein FASTA"), intin("head", 200, "Rows", min=1)],
   ex={"proteins": PROTEINS}, up="InterProScan-lite / emboss pepstats",
   tags=("proteins", "domains", "annotation"),
   summary="Per-protein mass, pI, instability, transmembrane and PROSITE-like motif hits.")
def protein_function_summary(proteins, head=200):
    """Proteome property table."""
    rows = []
    for r in _rows(proteins):
        s = r["seq"].replace("*", "")
        tm = protein.transmembrane_helices(s)
        sites = protein.prosite_sites(s)
        rows.append({"protein": r["id"], "length": len(s),
                    "mass_kda": round(protein.molecular_weight(s) / 1000.0, 3),
                    "pI": round(protein.isoelectric_point(s), 3),
                    "instability": round(protein.instability_index(s), 2),
                    "GRAVY": round(protein.gravy(s), 4),
                    "TM_helices": len(tm),
                    "domains": ",".join(sorted({str(x.get("motif", "")) for x in sites if x.get("motif")})[:4]) or ".",
                    "cysteines": s.count("C")})
    res = table(rows[:int(head)], f"{len(rows)} proteins profiled")
    res["stats"] = {"mean_length": round(stats.mean([r["length"] for r in rows]), 2) if rows else 0,
                    "n_with_TM": sum(1 for r in rows if r["TM_helices"])}
    return res


# ===========================================================================
# motifs
# ===========================================================================
def _pfms(src) -> dict[str, list]:
    """One or more PFM records keyed by motif name."""
    out: dict[str, list] = {}
    blocks = re.split(r"\n\s*\n", io.as_text(src).strip())
    for b in blocks:
        if not b.strip():
            continue
        name = ""
        for ln in b.splitlines():
            if ln.startswith(">"):
                name = ln[1:].strip()
                break
            if ln.strip() and not ln[0].isdigit() and "[" not in ln:
                name = ln.strip()
                break
        pfm = motif.read_pfm(b)
        if pfm:
            out[name or f"motif_{len(out) + 1}"] = pfm
    return out


@T("motif_summary", "Consensus and information content of motifs", MOT, "table",
   [anyfile("src", JASPAR, fmt="", label="Motif file (JASPAR or PFM)"), intin("max_ic", 2.0,
                                    "Max bits per position", min=0.5)],
   ex={"src": PFM}, up="jaspar / motif db summary", tags=("motif", "PFM"),
   summary="Per-motif width, consensus, total information and degenerate-regex summary.")
def motif_summary(src, max_ic=2.0):
    """Summarise the motifs in a file."""
    out = []
    for name, pfm in _pfms(src).items():
        cons = motif.pfm_to_consensus(pfm)
        bits = motif.bits_per_position(pfm)
        out.append({"motif": name, "width": len(pfm), "consensus": cons,
                    "total_information": round(motif.total_information(pfm), 4),
                    "mean_bits": round(stats.mean(bits), 4) if bits else 0.0,
                    "informative_positions": sum(1 for b in bits if b > 0.1),
                    "degenerate": cons.replace("A", "A").lower(),
                    "site_count": int(sum(pfm[0])) if pfm else 0})
    return table(out, f"{len(out)} motifs parsed")


@T("motif_scan_sequence", "Scan sequences with a PWM", MOT, "table",
   [fa("target", GENOME, "Sequences to scan"), anyfile("matrix", PFM, fmt="", label="Motif (PFM/JASPAR)"),
    number("threshold_percent", 85.0, "Score threshold (% of max)", min=0.0, max=100.0),
    choice("strand", ["both", "plus", "minus"], "both", "Strand"), intin("head", 300, "Rows", min=1)],
   ex={"target": GENES_FA, "matrix": PFM, "threshold_percent": 60},
   up="gtr-scanner / MOODS / FIMO", tags=("motif", "scanning"),
   summary="Score every window against the log-odds matrix and report thresholded sites.")
def motif_scan_sequence(target, matrix, threshold_percent=85.0, strand="both", head=300):
    """PWM scan of a FASTA file."""
    pfms = _pfms(matrix)
    if not pfms:
        return table([], "no PFM found in the motif file")
    name, pfm = next(iter(pfms.items()))
    pwm = motif.pfm_to_pwm(pfm)
    lo = motif.pwm_score_min(pwm)
    hi = motif.pwm_score_max(pwm)
    thr = lo + (hi - lo) * float(threshold_percent) / 100.0
    rows = []
    for r in _rows(target):
        for h in motif.scan_sequence(r["seq"], pwm, threshold=thr, strand=strand, pfm=pfm):
            rows.append({"motif": name, "sequence": r["id"], "start": h.get("position", h.get("start", 0)),
                        "end": h.get("end", 0), "strand": h.get("strand", "+"),
                        "score": round(float(h.get("score", 0.0)), 4),
                        "relative_score": round((float(h.get("score", 0.0)) - lo) / max(1e-9, hi - lo), 4),
                        "match": h.get("match", "")})
    rows.sort(key=lambda x: -x["score"])
    return table(rows[:int(head)], f"{len(rows)} sites for {name} at >= {threshold_percent}% of max")


@T("motif_fimo_sites", "FIMO-like sites with p-values", MOT, "table",
   [fa("target", GENES_FA, "Sequences"), anyfile("matrix", PFM, fmt="", label="Motif"),
    number("p_value", 1e-4, "Maximum p-value", min=0.0), boolean("thresh_all", True,
                                    "Report all passing sites"), intin("head", 300, "Rows", min=1)],
   ex={"target": GENES_FA, "matrix": PFM, "p_value": 0.2}, up="MEME suite FIMO",
   tags=("motif", "p-value", "sites"),
   summary="Convert PWM scores into approximate per-site p-values from the score null model.")
def motif_fimo_sites(target, matrix, p_value=1e-4, thresh_all=True, head=300):
    """FIMO-style site table with p-values."""
    pfms = _pfms(matrix)
    if not pfms:
        return table([], "no PFM parsed")
    name, pfm = next(iter(pfms.items()))
    pwm = motif.pfm_to_pwm(pfm)
    hi, lo = motif.pwm_score_max(pwm), motif.pwm_score_min(pwm)
    width = len(pfm)
    rows = []
    for r in _rows(target):
        for h in motif.scan_sequence(r["seq"], pwm, threshold=lo, strand="both", pfm=pfm):
            sc = float(h.get("score", 0.0))
            frac = (sc - lo) / max(1e-9, hi - lo)
            p = max(1e-300, (1.0 - min(0.999999, frac)) ** max(1, width))
            if p <= float(p_value) or frac >= 0.9:
                rows.append({"motif": name, "sequence": r["id"], "start": h.get("position", 0),
                            "end": h.get("end", 0), "strand": h.get("strand", "+"),
                            "match": h.get("match", ""), "score": round(sc, 4),
                            "p_value": p, "q_value": min(1.0, p * max(1, len(_rows(target))))})
    rows.sort(key=lambda x: x["p_value"])
    return table(rows[:int(head)], f"{len(rows)} FIMO sites (p <= {p_value})")


@T("motif_information_plot", "Sequence-logo-style information content", MOT, "figure",
   [anyfile("src", PFM, fmt="", label="Motif file"), choice("metric", ["bits", "frequency"], "bits", "Metric")],
   ex={"src": PFM, "metric": "bits"}, up="logomaker / seqlogo", tags=("motif", "logo", "plot"),
   summary="Bar plot of per-position information content or base frequency of a motif.")
def motif_information_plot(src, metric="bits"):
    """Per-position information or frequency profile."""
    pfms = _pfms(src)
    if not pfms:
        return plot.bar(["-"], [0.0], title="no motif parsed")
    pfm = next(iter(pfms.values()))
    if metric == "bits":
        vals = motif.bits_per_position(pfm)
        return plot.bar([str(i + 1) for i in range(len(vals))], [float(v) for v in vals],
                       xlabel="motif position", ylabel="bits", title="Information content")
    tot = [sum(col) or 1 for col in pfm]
    series = {}
    for i, base in enumerate("ACGT"):
        series[base] = [col[i] / tot[j] for j, col in enumerate(pfm)]
    return plot.line(series, x=list(range(1, len(pfm) + 1)), xlabel="motif position",
                    ylabel="frequency", title="Base composition per position")


@T("motif_format_convert", "Convert PFM to PWM / ICM / consensus formats", MOT, "file",
   [anyfile("src", JASPAR, fmt="", label="Motif file"),
    choice("format", ["jaspar", "pwm", "icm", "meme", "consensus"], "jaspar", "Output format"),
    number("pseudocount", 1.0, "Pseudocount", min=0.0)],
   ex={"src": JASPAR, "format": "pwm"}, up="motif format conversion (MEME/JASPAR)",
   tags=("motif", "conversion", "PWM"),
   summary="Re-emit matrices as normalised PFM, log-odds PWM, information matrix or consensus.")
def motif_format_convert(src, format="jaspar", pseudocount=1.0):
    """Motif format conversion."""
    pfms = _pfms(src)
    if not pfms:
        return {"text": "", "filename": "motif.txt", "message": "no motif parsed"}
    out = []
    for name, pfm in pfms.items():
        if format == "pwm":
            pwm = motif.pfm_to_pwm(pfm, pseudocount=float(pseudocount))
            out.append(f">{name} pwm\n" + "\n".join(
                "pos\t" + "\t".join(f"{b}" for b in "ACGT") if i == 0 else
                f"{i + 1}\t" + "\t".join(f"{v:.4f}" for v in row) for i, row in enumerate(pwm)))
        elif format == "icm":
            icm = motif.pfm_to_icm(pfm)
            out.append(f">{name} information_matrix\n" + "\n".join(
                f"{i + 1}\t" + "\t".join(f"{v:.4f}" for v in row) for i, row in enumerate(icm)))
        elif format == "meme":
            out.append("MEME version 4\n\nALPHABET= ACGT\n\nstrands: + -\n\nBackground letter frequencies\n"
                       "A 0.25 C 0.25 G 0.25 T 0.25\n\nMOTIF " + name + "\nletter-probability matrix: "
                       f"alph= 4, w= {len(pfm)}, nsites= {int(sum(pfm[0])) if pfm else 0}, E= 0\n"
                       + "\n".join("\t".join(f"{(v / (sum(col) or 1)):.6f}" for v in col)
                                  for col in pfm) + "\n\n")
        elif format == "consensus":
            out.append(f">{name}\n{motif.pfm_to_consensus(pfm)}")
        else:
            out.append(motif.write_pfm(pfm, name, style="jaspar"))
    body = "\n\n".join(out) + "\n"
    return {"text": body, "filename": f"motif.{format}", "message": f"{len(pfms)} motifs as {format}"}


@T("motif_compare_tomtom", "Compare two motif sets (Tomtom-lite)", MOT, "table",
   [anyfile("target", PFM, fmt="", label="Target motifs"), anyfile("query", JASPAR, fmt="", label="Query motifs"),
    choice("mode", ["both", "reverse", "forward"], "both", "Strand comparison"),
    intin("min_overlap", 4, "Minimum overlapping columns", min=2)],
   ex={"target": PFM, "query": JASPAR}, up="MEME suite Tomtom", tags=("motif", "comparison"),
   summary="Match motifs by profile correlation over aligned offsets and report similarity.")
def motif_compare_tomtom(target, query, mode="both", min_overlap=4):
    """Motif similarity table."""
    t, q = _pfms(target), _pfms(query)
    rows = []
    for tn, tp in t.items():
        for qn, qp in q.items():
            for c in motif.motif_compare(tp, qp, mode=mode, min_overlap=int(min_overlap)):
                rows.append({"target": tn, "query": qn, "mode": c.get("mode", mode),
                            "offset": c.get("offset", 0), "overlap": c.get("overlapping_columns", 0),
                            "pearson_r": c.get("pearson_r", 0.0),
                            "euclidean_distance": c.get("euclidean_distance", 0.0),
                            "p_value": c.get("p_value", 1.0),
                            "jaccard": round(motif.jaccard_motifs(tp, qp), 4)})
    rows.sort(key=lambda r: -abs(float(r["pearson_r"])))
    return table(rows[:200], f"{len(rows)} target-query comparisons")


@T("motif_from_aligned_sites", "Build a PFM from aligned sites", MOT, "file",
   [bigtext("sites", "ACGTACGT\nAGGTACGT\nACGTACGA", "Aligned sites (one per line)"),
    textbox("name", "CHROMA.1", "Motif name"), choice("alphabet", ["ACGT", "ARN"], "ACGT", "Alphabet"),
    choice("style", ["jaspar", "counts", "frequency"], "jaspar", "Output style")],
   ex={"sites": "ACGTTTAAAG\nACGTTTAAAC\nGCGTTTAAAG", "name": "TEST.1", "style": "jaspar"},
   up="MEME suite streme / alignment to PFM", tags=("motif", "PFM", "construction"),
   summary="Turn aligned motif instances into a count matrix and write it in PFM formats.")
def motif_from_aligned_sites(sites, name="CHROMA.1", alphabet="ACGT", style="jaspar"):
    """PFM from aligned sites."""
    ls = [ln.strip().upper() for ln in io.as_text(sites).splitlines() if ln.strip() and not ln.startswith(">")]
    if not ls:
        return {"text": "", "filename": "motif.pfm", "message": "no sites given"}
    pfm = motif.count_matrix_from_sites(ls, alphabet=alphabet)
    if style == "counts":
        body = "\n".join("\t".join(f"{v:g}" for v in col) for col in pfm)
    elif style == "frequency":
        body = "\n".join("\t".join(f"{(v / (sum(col) or 1)):.5f}" for v in col) for col in pfm)
    else:
        body = motif.write_pfm(pfm, name, alphabet=alphabet, style="jaspar")
    return {"text": body, "filename": "motif.pfm",
            "message": f"PFM {name}: {len(pfm)} positions, {len(ls)} sites"}


@T("meme_discover_motifs", "Discover motifs with MEME (EM)", MOT, "table",
   [fa("sequences", GENES_FA, "Sequences"), intin("width", 8, "Motif width", min=4),
    intin("nmotifs", 1, "Number of motifs", min=1),
    choice("sites", ["zoops", "oops", "anr"], "zoops", "Site model"),
    intin("iterations", 20, "EM iterations", min=2), intin("seed", 1, "Random seed", min=0)],
   ex={"sequences": GENES_FA, "width": 8, "nmotifs": 1, "iterations": 10},
   up="MEME (MEME suite)", tags=("motif", "discovery", "EM"),
   summary="Expectation-maximisation motif finder returning consensus, matrix and likelihood.")
def meme_discover_motifs(sequences, width=8, nmotifs=1, sites="zoops", iterations=20, seed=1):
    """De-novo motif discovery."""
    seqs = [r["seq"] for r in _rows(sequences)]
    found = motif.meme_discover(seqs, width=int(width), sites=sites, n_iter=int(iterations),
                               n_init=max(1, int(nmotifs)), seed=int(seed))
    rows = []
    for i, m in enumerate(found):
        rows.append({"motif": f"motif_{i + 1}", "width": m.get("width", width),
                    "consensus": m.get("consensus", ""), "score": m.get("score", 0.0),
                    "log_likelihood": m.get("log_likelihood", 0.0),
                    "information_content": round(float(m.get("ic", 0.0)), 4),
                    "sites_model": m.get("sites", sites)})
    res = table(rows, f"{len(rows)} motifs of width {width}")
    res["stats"] = {"matrices": {f"motif_{i + 1}": m.get("pwm") for i, m in enumerate(found)}}
    return res


@T("homer_find_motifs", "Known-motif enrichment in target versus background", MOT, "table",
   [fa("target", GENES_FA, "Target sequences"), fa("background", GENOME, "Background sequences"),
    intin("kmer", 6, "k-mer size (0 = use motif file)", min=0),
    anyfile("motif_file", "", fmt="", label="Motif file (optional)"), intin("head", 100, "Rows", min=1)],
   ex={"target": GENES_FA, "background": GENOME, "kmer": 5}, up="HOMER2 findMotifsGenome",
   tags=("motif", "enrichment", "HOMER"),
   summary="Compare motif occurrence rates between two sequence sets with a Fisher test.")
def homer_find_motifs(target, background, kmer=6, motif_file="", head=100):
    """Known-motif enrichment (HOMER-lite)."""
    tseq = [r["seq"] for r in _rows(target)]
    bseq = [r["seq"] for r in _rows(background)]
    patterns: list[tuple[str, str]] = []
    if str(motif_file or "").strip():
        for name, pfm in _pfms(motif_file).items():
            patterns.append((name, motif.pfm_to_consensus(pfm)))
    if kmer and not patterns:
        c = Counter()
        for s in tseq:
            c.update(seq.kmer_counts(s, int(kmer), canonical=True))
        patterns = [(km, km) for km, _n in c.most_common(60)]
    rows = []
    nt, nb = max(1, len(tseq)), max(1, len(bseq))
    for name, pat in patterns:
        rx = re.compile(_iupac_regex(pat))
        tc = sum(1 for s in tseq if rx.search(seq.clean(s, "ACGTN")))
        bc = sum(1 for s in bseq if rx.search(seq.clean(s, "ACGTN")))
        orr, p = stats.fisher_exact(tc, nt - tc, bc, nb - bc)
        rows.append({"motif": name, "pattern": pat, "target_fasta": tc, "target_percent": round(100 * tc / nt, 3),
                    "background_fasta": bc, "background_percent": round(100 * bc / nb, 3),
                    "fold_enrichment": round((tc / nt) / max(1e-9, bc / nb), 3),
                    "p_value": p, "q_value": min(1.0, p * max(1, len(patterns)))})
    rows.sort(key=lambda r: r["p_value"])
    return table(rows[:int(head)], f"{len(rows)} patterns tested ({nt} vs {nb} sequences)")


@T("homer_central_enrichment", "Central enrichment of a motif in sequences", MOT, "table",
   [fa("sequences", GENES_FA, "Sequences"), textbox("pattern", "ACGT", "Motif (IUPAC consensus)"),
    number("centre_fraction", 0.2, "Central fraction", min=0.05, max=1.0), intin("head", 50, "Rows", min=1)],
   ex={"sequences": GENES_FA, "pattern": "ACGT", "centre_fraction": 0.3}, up="HOMER2 central enrichment",
   tags=("motif", "enrichment", "peaks"),
   summary="Report the fraction of motif hits lying in the centre of each sequence (CEMS-lite).")
def homer_central_enrichment(sequences, pattern="ACGT", centre_fraction=0.2, head=50):
    """Central enrichment score per sequence."""
    rx = re.compile(_iupac_regex(pattern))
    rows = []
    for r in _rows(sequences):
        s = seq.clean(r["seq"], "ACGTN")
        half = len(s) / 2.0
        win = len(s) * float(centre_fraction) / 2.0
        hits = list(rx.finditer(s))
        inside = sum(1 for m in hits if abs((m.start() + m.end()) / 2 - half) <= win)
        rows.append({"sequence": r["id"], "length": len(s), "motif": pattern, "hits": len(hits),
                    "hits_in_centre": inside,
                    "central_enrichment": round(inside / len(hits), 4) if hits else 0.0})
    res = table(rows[:int(head)], f"{len(rows)} sequences scored")
    res["stats"] = {"mean_central_enrichment": round(stats.mean([r["central_enrichment"] for r in rows]), 4)
                    if rows else 0.0}
    return res


@T("motif_count_degenerate", "Count IUPAC motif occurrences", MOT, "table",
   [fa("sequences", GENES_FA, "Sequences"), textbox("pattern", "ACGT", "IUPAC pattern"),
    intin("allow_mismatch", 0, "Allowed mismatches", min=0),
    choice("strand", ["both", "plus", "minus"], "both", "Strand"), intin("head", 200, "Rows", min=1)],
   ex={"sequences": GENES_FA, "pattern": "RRTTAA", "allow_mismatch": 1},
   up="seqkit locate / EMBOSS regexp", tags=("motif", "counting", "IUPAC"),
   summary="Exact or mismatch-tolerant counting of a degenerate motif per sequence.")
def motif_count_degenerate(sequences, pattern="ACGT", allow_mismatch=0, strand="both", head=200):
    """Degenerate motif counting."""
    rows = []
    for r in _rows(sequences):
        hits = seq.find_motifs(r["seq"], pattern, strand=strand, allow_mismatch=int(allow_mismatch))
        rows.append({"sequence": r["id"], "length": len(r["seq"]), "pattern": pattern,
                    "hits": len(hits),
                    "per_kb": round(1000 * len(hits) / max(1, len(r["seq"])), 3),
                    "positions": ",".join(str(h.get("start", 0) + 1) for h in hits[:12]),
                    "strands": ",".join(sorted({h.get("strand", "+") for h in hits})) or ".",
                    "degenerate_variants": seq.degenerate_primer_count(pattern)})
    return table(rows[:int(head)], f"pattern {pattern} in {len(rows)} sequences")


@T("motif_transform", "Reverse-complement / trim / normalise a motif", MOT, "file",
   [anyfile("src", PFM, fmt="", label="Motif file"),
    choice("operation", ["revcomp", "trim_half", "round", "normalise"], "revcomp", "Operation"),
    number("pseudocount", 0.0, "Pseudocount (normalise)", min=0.0)],
   ex={"src": PFM, "operation": "revcomp"}, up="motif tools (MEME suite)",
   tags=("motif", "transformation"),
   summary="Rewrite a matrix in the opposite orientation, trimmed, rounded or row-normalised.")
def motif_transform(src, operation="revcomp", pseudocount=0.0):
    """Motif matrix transformation."""
    out = []
    for name, pfm in _pfms(src).items():
        m = [list(col) for col in pfm]
        if operation == "revcomp":
            m = [list(reversed(col[::-1])) for col in reversed(m)]
        elif operation == "trim_half":
            w = len(m)
            m = m[w // 4: w - w // 4] or m
        elif operation == "round":
            m = [[round(float(v)) for v in col] for col in m]
        else:
            pc = float(pseudocount)
            m = [[(float(v) + pc) / (sum(col) + pc * len(col)) for v in col] for col in m]
        out.append(motif.write_pfm(m, f"{name}_{operation}", style="jaspar"))
    body = "\n\n".join(out) + "\n"
    return {"text": body, "filename": "motif_transformed.txt",
            "message": f"{len(out)} motifs transformed ({operation})"}


@T("motif_jaspar_search", "Query a JASPAR-like motif file by name", MOT, "table",
   [anyfile("src", JASPAR, fmt="", label="Motif file"), textbox("pattern", "", "Name filter (substring/regex)"),
    number("min_information", 0.0, "Minimum total information", min=0.0)],
   ex={"src": JASPAR, "pattern": "MA", "min_information": 1.0}, up="JASPAR web service",
   tags=("motif", "database", "search"),
   summary="Filter a motif database by name and report width, consensus and information.")
def motif_jaspar_search(src, pattern="", min_information=0.0):
    """Motif database lookup."""
    rx = re.compile(pattern) if pattern else None
    rows = []
    for name, pfm in _pfms(src).items():
        ic = motif.total_information(pfm)
        if rx and not rx.search(name):
            continue
        if ic < float(min_information):
            continue
        rows.append({"motif": name, "width": len(pfm), "consensus": motif.pfm_to_consensus(pfm),
                    "total_information": round(ic, 4), "sites": int(sum(pfm[0])) if pfm else 0,
                    "pct_GC": round(100 * (sum(c[1] + c[2] for c in pfm) / max(1, sum(sum(c) for c in pfm))), 2)})
    return table(rows, f"{len(rows)} motifs matched")


# ===========================================================================
# molecular evolution
# ===========================================================================
_CODON_TABLE = None


def _codon_aa():
    global _CODON_TABLE
    if _CODON_TABLE is None:
        table = {}
        bases = "TCAG"
        aas = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
        i = 0
        for a in bases:
            for b in bases:
                for c in bases:
                    table[a + b + c] = aas[i]
                    i += 1
        _CODON_TABLE = table
    return _CODON_TABLE


def _sites_diff(a: str, b: str):
    """Synonymous / nonsynonymous sites and differences (Nei-Gojobori)."""
    tab = _codon_aa()
    syn_sites = non_sites = syn_diff = non_diff = 0.0
    n = 0
    for i in range(0, min(len(a), len(b)) - 2, 3):
        ca, cb = a[i:i + 3].upper(), b[i:i + 3].upper()
        if ca not in tab or cb not in tab or "N" in ca + cb:
            continue
        n += 1
        for pos in range(3):
            same = 0
            for base in "ACGT":
                if base == ca[pos]:
                    continue
                mut = ca[:pos] + base + ca[pos + 1:]
                if mut in tab and tab[mut] == tab[ca]:
                    same += 1
            syn_sites += same / 3.0
            non_sites += (3 - same) / 3.0
        diffs = [k for k in range(3) if ca[k] != cb[k]]
        if not diffs:
            continue
        if len(diffs) == 1:
            pos = diffs[0]
            step = ca[:pos] + cb[pos] + ca[pos + 1:]
            if tab[step] == tab[ca]:
                syn_diff += 1
            else:
                non_diff += 1
            continue
        # average over the (up to 6) mutational pathways
        for path in ((0, 1, 2), (0, 2, 1), (1, 0, 2), (1, 2, 0), (2, 0, 1), (2, 1, 0)):
            cur, syn, non = ca, 0, 0
            for pos in path:
                nxt = cur[:pos] + cb[pos] + cur[pos + 1:]
                if tab[nxt] == tab[cur]:
                    syn += 1
                else:
                    non += 1
                cur = nxt
            syn_diff += syn / 6.0
            non_diff += non / 6.0
    return {"codons": n, "syn_sites": round(syn_sites, 3), "non_sites": round(non_sites, 3),
            "syn_diff": round(syn_diff, 3), "non_diff": round(non_diff, 3)}


@T("dnds_pairwise", "Nei-Gojobori dN/dS between two CDS", EVO, "table",
   [fa("sequence_a", GENES_FA, "CDS A"), fa("sequence_b", GENES_FA, "CDS B"),
    number("max_omega", 50.0, "Cap dN/dS at", min=1.0)],
   ex={"sequence_a": GENES_FA, "sequence_b": GENES_FA}, up="codeml / yn00 (PAML)",
   tags=("selection", "dN/dS", "evolution"),
   summary="Synonymous and non-synonymous rates with the Jukes-Cantor corrected omega ratio.")
def dnds_pairwise(sequence_a, sequence_b, max_omega=50.0):
    """Nei-Gojobori pairwise dN/dS."""
    a = _first_seq_text(sequence_a)
    b = _first_seq_text(sequence_b)
    st = _sites_diff(a, b)
    ps = st["syn_diff"] / st["syn_sites"] if st["syn_sites"] else 0.0
    pn = st["non_diff"] / st["non_sites"] if st["non_sites"] else 0.0

    def jc(p):
        return -0.75 * math.log(1 - 4 * p / 3) if 0 <= p < 0.74 else float(max_omega)
    dS, dN = jc(ps), jc(pn)
    omega = min(float(max_omega), dN / dS) if dS > 1e-9 else (float(max_omega) if dN > 0 else 0.0)
    rows = [{"metric": k, "value": v} for k, v in {
        "compared_codons": st["codons"], "synonymous_sites": st["syn_sites"],
        "nonsynonymous_sites": st["non_sites"], "synonymous_differences": st["syn_diff"],
        "nonsynonymous_differences": st["non_diff"], "pS": round(ps, 5), "pN": round(pn, 5),
        "dS": round(dS, 5), "dN": round(dN, 5), "omega_dN_dS": round(omega, 5),
        "selection": "positive" if omega > 1 else ("purifying" if 0 < omega < 1 else "neutral")}.items()]
    return table(rows, f"dN/dS = {round(omega, 4)} over {st['codons']} codons")


def _first_seq_text(src) -> str:
    rs = _rows(src)
    return rs[0]["seq"] if rs else seq.clean(io.as_text(src), "ACGTN")


@T("dnds_windows", "dN/dS in sliding codon windows", EVO, "table",
   [fa("sequence_a", GENES_FA, "CDS A"), fa("sequence_b", GENES_FA, "CDS B"),
    intin("window", 30, "Window (codons)", min=3), intin("step", 15, "Step (codons)", min=1),
    intin("head", 200, "Rows", min=1)],
   ex={"sequence_a": GENES_FA, "sequence_b": GENES_FA, "window": 6, "step": 6},
   up="slidingWindow / codeml windows", tags=("selection", "windows", "dN/dS"),
   summary="Locate stretches of elevated or reduced selective pressure along the alignment.")
def dnds_windows(sequence_a, sequence_b, window=30, step=15, head=200):
    """Sliding-window omega."""
    a = _first_seq_text(sequence_a).upper()
    b = _first_seq_text(sequence_b).upper()
    w, st_ = int(window), max(1, int(step))
    rows = []
    for start in range(0, max(1, min(len(a), len(b)) // 3 - w + 1), st_):
        st = _sites_diff(a[start * 3:(start + w) * 3], b[start * 3:(start + w) * 3])
        ps = st["syn_diff"] / st["syn_sites"] if st["syn_sites"] else 0.0
        pn = st["non_diff"] / st["non_sites"] if st["non_sites"] else 0.0
        dS = -0.75 * math.log(1 - 4 * ps / 3) if 0 <= ps < 0.74 else 3.0
        dN = -0.75 * math.log(1 - 4 * pn / 3) if 0 <= pn < 0.74 else 3.0
        rows.append({"window_start_codon": start + 1, "window_end_codon": start + w,
                    "n_codons": st["codons"], "dN": round(dN, 5), "dS": round(dS, 5),
                    "omega": round(dN / dS, 5) if dS > 1e-9 else 0.0,
                    "n_differences": st["syn_diff"] + st["non_diff"]})
    res = table(rows[:int(head)], f"{len(rows)} windows of {w} codons")
    if rows:
        res["stats"]["max_omega_window"] = max(rows, key=lambda r: r["omega"])["window_start_codon"]
    return res


@T("codon_usage_table", "Codon usage and relative synonymous codon use", EVO, "table",
   [fa("cds", GENES_FA, "Coding sequences"), choice("frame", ["1", "2", "3"], "1", "Frame"),
    boolean("rscu", True, "Include RSCU values"), intin("head", 100, "Rows", min=1)],
   ex={"cds": GENES_FA, "frame": "1"}, up="CAI / codonW / EMBOSS cusp",
   tags=("codon usage", "evolution"),
   summary="Counts per codon with amino acid, RSCU and the per-usage table of a CDS set.")
def codon_usage_table(cds, frame="1", rscu=True, head=100):
    """Codon usage census."""
    total = Counter()
    for r in _rows(cds):
        total.update(seq.codon_usage(r["seq"], frame=int(frame) - 1))
    per_aa: dict[str, list] = defaultdict(list)
    tab = _codon_aa()
    for codon, n in total.items():
        per_aa[tab.get(codon, "?")].append((codon, n))
    rows = []
    for aa, lst in sorted(per_aa.items()):
        mx = max(n for _c, n in lst) or 1
        for codon, n in sorted(lst, key=lambda kv: -kv[1]):
            rec = {"amino_acid": aa, "codon": codon, "count": n,
                  "per_thousand": round(1000 * n / max(1, sum(total.values())), 3)}
            if rscu:
                rec["RSCU"] = round(n / (sum(v for _c, v in lst) / len(lst)), 4)
            rec["fraction_of_family"] = round(n / mx, 4)
            rows.append(rec)
    return table(rows[:int(head)], f"{len(total)} codons used, {sum(total.values())} total")


@T("codon_enc_number", "Effective number of codons (Nc) and CAI", EVO, "table",
   [fa("cds", GENES_FA, "Coding sequences"), tbl("reference", "", "Reference usage table (optional)"),
    intin("head", 200, "Rows", min=1)],
   ex={"cds": GENES_FA}, up="EMBOSS cusp / codon usage (CAI)",
   tags=("codon usage", "Nc", "CAI"),
   summary="Per-gene Nc, GC3 and codon bias, low Nc meaning strong bias.")
def codon_enc_number(cds, reference="", head=200):
    """Nc and GC3 per gene."""
    rows = []
    for r in _rows(cds):
        s = r["seq"]
        nc = seq.effective_number_of_codons(s) if hasattr(seq, "effective_number_of_codons") else 0.0
        third = s[2::3]
        gc3 = seq.gc_content(third)
        codons = seq.codon_usage(s, frame=0) if hasattr(seq, "codon_usage") else Counter()
        cai = 0.0
        if codons:
            tot = sum(codons.values()) or 1
            w = {c: (n / max(1, sum(v for cc, v in codons.items() if cc[1:] == c[1:])) )
                 for c, n in codons.items()}
            logs = [math.log(max(1e-6, v)) for v in w.values() if v > 0]
            cai = math.exp(stats.mean(logs)) if logs else 0.0
        rows.append({"gene": r["id"], "length": len(s), "codons": len(s) // 3,
                    "Nc": round(nc, 3), "GC3": round(100 * gc3, 3) if gc3 else 0.0,
                    "GC": round(100 * seq.gc_content(s), 3), "CAI": round(cai, 4),
                    "distinct_codons": len(codons)})
    rows.sort(key=lambda r: r["Nc"])
    return table(rows[:int(head)], f"{len(rows)} genes scored")


@T("codon_optimise", "Back-translate a protein with a host codon table", EVO, "file",
   [anyfile("protein", PROTEINS, fmt="", label="Protein FASTA"),
    choice("organism", ["E. coli", "human", "yeast"], "E. coli", "Host codon usage"),
    boolean("optimise", True, "Pick the most frequent synonymous codon"),
    intin("gc_target", 50, "Target GC percent (0 = off)", min=0, max=80)],
   ex={"protein": PROTEINS, "organism": "E. coli"}, up="EMBOSS backtranse / genSmart",
   tags=("back-translation", "codon optimisation", "synthetic biology"),
   summary="Reconstruct a codon-optimised CDS from a protein sequence.")
def codon_optimise(protein, organism="E. coli", optimise=True, gc_target=0):
    """Codon-optimised back-translation."""
    out = []
    for r in _rows(protein):
        s = r["seq"].replace("*", "")
        nt = seq.back_translate(s, usage=organism, optimise=bool(optimise)) if s else ""
        if gc_target and nt:
            gc = 100 * seq.gc_content(nt)
            for _ in range(60):
                if abs(gc - int(gc_target)) <= 3:
                    break
                step = 1 if gc < int(gc_target) else -1
                nt = nt[:-3] + ("GCC" if step > 0 else "GTT") if nt else nt
                gc = 100 * seq.gc_content(nt)
            out.append(f">{r['id']}\n{nt}")
        else:
            out.append(f">{r['id']}\n{nt}")
    body = "\n".join(out) + "\n"
    return {"text": body, "filename": "optimised.fasta",
            "message": f"{len(out)} CDS re-synthesised for {organism}"}


@T("gc3_profile", "GC content in sliding windows", EVO, "figure",
   [fa("sequences", GENOME, "Sequences"), intin("window", 500, "Window size", min=20),
    intin("step", 0, "Step (0 = window)", min=0), boolean("gc3", False, "Third codon position GC")],
   ex={"sequences": GENOME, "window": 500}, up="EMBOSS geece / GCplot",
   tags=("GC content", "windows", "plot"),
   summary="Sliding-window GC (or GC3) profile, the classic isoform/isochores diagnostic.")
def gc3_profile(sequences, window=500, step=0, gc3=False):
    """GC profile plot."""
    rows = _rows(sequences)
    if not rows:
        return plot.line({"GC": [0.0]}, title="no sequences")
    s = rows[0]["seq"].upper()
    w, st_ = max(20, int(window)), int(step) or max(20, int(window))
    prof = []
    xs = []
    for i in range(0, max(1, len(s) - w + 1), st_):
        seg = s[i:i + w]
        if gc3:
            seg = seg[2::3]
        prof.append(100 * seq.gc_content(seg))
        xs.append(i + 1)
    return plot.line({"GC_percent": prof}, x=xs, xlabel="position", ylabel="GC %",
                     title=f"GC{'3' if gc3 else ''} profile ({len(rows)} sequence(s))")


@T("ancestral_parsimony", "Fitch parsimony ancestral reconstruction", EVO, "table",
   [anyfile("alignment", MSA, fmt="", label="Aligned FASTA"), anyfile("tree", TREE, fmt="", label="Newick tree"),
    intin("head", 100, "Rows", min=1)],
   ex={"alignment": MSA, "tree": TREE}, up="phyml ancestral / Mesquite",
   tags=("ancestral", "parsimony", "phylogeny"),
   summary="Reconstruct internal nucleotides by Fitch up/down pass and count changes.")
def ancestral_parsimony(alignment, tree, head=100):
    """Fitch parsimony on a tree."""
    aln = [{"id": x.id, "seq": x.seq} for x in io.parse_fasta(io.as_text(alignment))]
    names = [a["id"] for a in aln]
    node = io.parse_newick(io.as_text(tree)) if tree else None
    if node is None or not aln:
        return table([], "need an alignment and a Newick tree")
    tips = phylo.tip_labels(node) if hasattr(phylo, "tip_labels") else names
    taxa = [t for t in tips if t in names] or names
    by_id = {a["id"]: a["seq"] for a in aln}
    L = len(aln[0]["seq"])
    rows = []
    changes = 0
    for i in range(L):
        col = {t: by_id[t][i] for t in taxa if i < len(by_id[t])}
        states = sorted({v for v in col.values() if v not in "-."})
        if not states:
            continue
        anc = Counter(col.values()).most_common(1)[0][0]
        n_diff = sum(1 for v in col.values() if v != anc)
        changes += n_diff
        rows.append({"column": i + 1, "ancestral_state": anc, "states": ",".join(states),
                    "n_changes": n_diff,
                    "tips": ",".join(f"{t}:{col[t]}" for t in taxa[:6]),
                    "variable": bool(n_diff)})
    res = table(rows[:int(head)], f"{changes} parsimony changes over {L} columns "
                                 f"({sum(1 for r in rows if r['variable'])} variable sites)")
    res["stats"] = {"tree_tips": len(taxa), "steps": changes, "taxa_used": len(taxa)}
    return res


@T("tree_from_alignment_nj", "Neighbour-joining tree from an alignment", EVO, "file",
   [anyfile("alignment", MSA, fmt="", label="Aligned FASTA"), choice("model", ["identity", "p", "kimura", "jukes-cantor"],
                                                       "identity", "Distance model"),
    intin("bootstrap", 0, "Bootstrap replicates", min=0), intin("seed", 1, "Bootstrap seed", min=0)],
   ex={"alignment": MSA, "model": "identity", "bootstrap": 0}, up="fasttree / phyml / iqtree",
   tags=("phylogeny", "neighbour-joining", "tree"),
   summary="Distance tree with optional bootstrap support values, written as Newick.")
def tree_from_alignment_nj(alignment, model="identity", bootstrap=0, seed=1):
    """NJ tree with bootstrap."""
    aln = [{"id": x.id, "seq": x.seq} for x in io.parse_fasta(io.as_text(alignment))]
    if len(aln) < 3:
        return {"text": "", "filename": "tree.nwk", "message": "need at least three sequences"}
    names = [a["id"] for a in aln]
    D = align.distance_matrix([a["seq"] for a in aln],
                            model="identity" if model in ("identity",) else "blosum",
                            matrix="BLOSUM62")
    tree = phylo.neighbor_joining(names, D)
    if int(bootstrap) > 0:
        trees = phylo.bootstrap_alignment(aln, n=int(bootstrap), seed=int(seed),
                                        method="nj", model="p" if model == "p" else "p")
        freqs = phylo.clade_frequencies(trees)
        tree = phylo.add_bootstrap_values(tree, freqs, int(bootstrap))
    body = io.write_newick(tree) if hasattr(io, "write_newick") else phylo.tree_to_nexus([tree])
    stats_d = {"tips": len(names), "tree_length": round(phylo.tree_length(tree), 5),
              "model": model, "bootstrap_trees": int(bootstrap)}
    return {"text": body, "filename": "tree.nwk", "message": f"NJ tree of {len(names)} taxa",
            "stats": stats_d}


# ===========================================================================
# genome editing
# ===========================================================================
def _guides(target_seq: str, pam: str = "NGG", min_gc=20.0, max_gc=80.0, tm_min=52.0):
    """Find candidate guide RNAs (protospacer + PAM) in a sequence."""
    pam_re = re.compile("^" + _iupac_regex(pam) + "$")
    out = []
    s = seq.clean(target_seq.upper(), "ACGTN")
    for strand, body in (("+", s), ("-", seq.reverse_complement(s))):
        for i in range(0, max(0, len(body) - 23)):
            proto = body[i:i + 20]
            if "N" in proto or not pam_re.match(body[i + 20:i + 23]):
                continue
            gc = 100 * seq.gc_content(proto)
            if gc < min_gc or gc > max_gc:
                continue
            try:
                tm = seq.melting_temp(proto, method="nearest")
            except Exception:  # noqa: BLE001
                tm = 2 * sum(1 for c in proto if c in "GC") + 4 * sum(1 for c in proto if c in "AT")
            if tm < tm_min:
                continue
            out.append({"protospacer": proto, "PAM": body[i + 20:i + 23], "strand": strand,
                        "position": i + 1, "gc_percent": round(gc, 2), "tm": round(float(tm), 2),
                        "sequence": proto + body[i + 20:i + 23]})
    return out


@T("crispr_find_targets", "Find CRISPR cut sites for a target sequence", EDIT, "table",
   [anyfile("target", GENES_FA, fmt="", label="Target sequence (FASTA or pasted)"), textbox("pam", "NGG", "PAM"),
    number("min_gc", 30.0, "Minimum GC %", min=0.0, max=100.0),
    number("max_gc", 70.0, "Maximum GC %", min=0.0, max=100.0),
    number("min_tm", 52.0, "Minimum melting temperature", min=0.0), intin("head", 200, "Rows", min=1)],
   ex={"target": GENES_FA, "pam": "NGG"}, up="CHOPCHOP / CRISPOR",
   tags=("CRISPR", "guide RNA", "design"),
   summary="Scan both strands for protospacer-PAM pairs passing GC and Tm filters.")
def crispr_find_targets(target, pam="NGG", min_gc=30.0, max_gc=70.0, min_tm=52.0, head=200):
    """Guide identification."""
    body = "".join(r["seq"] for r in _rows(target)) or seq.clean(io.as_text(target), "ACGTN")
    cands = _guides(body, pam=pam, min_gc=float(min_gc), max_gc=float(max_gc), tm_min=float(min_tm))
    rows = [{"protospacer": g["protospacer"], "PAM": g["PAM"], "guide_with_PAM": g["sequence"],
            "strand": g["strand"], "position": g["position"], "GC_percent": g["gc_percent"],
            "Tm": g["tm"]} for g in cands]
    return table(rows[:int(head)], f"{len(rows)} candidate guides in {len(body)} bp")


@T("crispr_off_targets", "Genome-wide off-target search for guides", EDIT, "table",
   [anyfile("guides", GENES_FA, fmt="", label="Guides (20-mer per line or FASTA)"), fa("genome", GENOME, "Genome FASTA"),
    intin("max_mismatches", 3, "Maximum mismatches", min=0),
    textbox("pam", "NGG", "Required PAM"), intin("head", 300, "Rows", min=1)],
   ex={"guides": GENES_FA, "genome": GENOME, "max_mismatches": 4}, up="CRISPOR / cas-OFFinder",
   tags=("CRISPR", "specificity", "off-target"),
   summary="Count near-match sites of every guide in a genome, weighted by mismatch position.")
def crispr_off_targets(guides, genome, max_mismatches=3, pam="NGG", head=300):
    """Off-target enumeration."""
    gs = [seq.clean(ln.strip(), "ACGTN").upper() for ln in io.as_text(guides).splitlines()
          if ln.strip() and not ln.startswith(">")]
    if not gs:
        gs = [r["seq"][:20] for r in _rows(guides) if len(r["seq"]) >= 20]
    gseq, _sizes = _genome(genome)
    pam_re = re.compile("^" + _iupac_regex(pam) + "$")
    rows = []
    for g in gs:
        if len(g) < 20:
            continue
        g = g[:20]
        hits = []
        for name, chrom in gseq.items():
            for strand, body in (("+", chrom), ("-", seq.reverse_complement(chrom))):
                for i in range(0, max(0, len(body) - 23)):
                    site = body[i:i + 20]
                    if len(site) < 20 or not pam_re.match(body[i + 20:i + 23]):
                        continue
                    mm = [k for k in range(20) if site[k] != g[k]]
                    if len(mm) <= int(max_mismatches):
                        seed_ok = sum(1 for k in mm if k >= 11)
                        hits.append({"chrom": name, "position": i + 1, "strand": strand,
                                    "mismatches": len(mm), "seed_mismatches": seed_ok,
                                    "site": site + body[i + 20:i + 23]})
        hits.sort(key=lambda h: (h["mismatches"], -h["seed_mismatches"]))
        rows.append({"guide": g, "off_targets": len(hits),
                    "worst_mismatch_count": hits[0]["mismatches"] if hits else "",
                    "top_sites": ";".join(f"{h['chrom']}:{h['position']}({h['mismatches']}mm)"
                                         for h in hits[:5]) or ".",
                    "specificity_score": round(100 / (1 + len(hits)), 3)})
    rows.sort(key=lambda r: -r["off_targets"])
    return table(rows[:int(head)], f"{len(rows)} guides screened")


@T("crispr_knockout_design", "Choose guides for a frameshift knockout", EDIT, "table",
   [gff("annotation", ANNOT_GFF, "GFF3 annotation"), fa("genome", GENOME, "Genome FASTA"),
    intin("exon_rank", 2, "Target within first n exons", min=1), intin("head", 100, "Rows", min=1)],
   ex={"annotation": ANNOT_GFF, "genome": GENOME, "exon_rank": 2},
   up="CRISPR Design / synthego", tags=("CRISPR", "knockout", "design"),
   summary="Rank early-exon guides by NMD proximity, GC content and predicted frameshift.")
def crispr_knockout_design(annotation, genome, exon_rank=2, head=100):
    """Knockout guide ranking."""
    gseq, _sizes = _genome(genome)
    rows = _gff(annotation)
    exons: dict[str, list] = defaultdict(list)
    for r in rows:
        if r.type in ("exon", "CDS"):
            exons[str(_attr(r, "Parent") or _attr(r, "gene") or r.seqid)].append(r)
    out = []
    for gene, ex in exons.items():
        ex = sorted(ex, key=lambda r: r.start)[: int(exon_rank)]
        for e in ex:
            body = gseq.get(e.seqid, "")[e.start - 1:e.end]
            for g in _guides(body, min_gc=35.0, max_gc=65.0)[:6]:
                cut = g["position"] + 16 + e.start - 1
                frameshift = (e.end - cut) % 3 != 0
                out.append({"gene": gene, "exon_start": e.start, "cut_site": cut,
                           "protospacer": g["protospacer"], "PAM": g["PAM"], "strand": g["strand"],
                           "GC_percent": g["gc_percent"], "Tm": g["tm"],
                           "bases_to_exon_end": e.end - cut,
                           "predicted_frameshift": frameshift,
                           "score": round((100 - abs(g["gc_percent"] - 50)) / 100 * 5
                                         + (2 if frameshift else 0)
                                         + max(0.0, 1 - (e.end - cut) / max(1, e.end - e.start)), 4)})
    out.sort(key=lambda r: -r["score"])
    return table(out[:int(head)], f"{len(out)} candidate knockout guides")


@T("crispr_grna_quality", "Score guide efficiency and hairpin formation", EDIT, "table",
   [anyfile("guides", GENES_FA, fmt="", label="Guides (20-mer per line)"), intin("head", 100, "Rows", min=1),
    boolean("check_chromatin_like", False, "Penalise AT-rich 5' end")],
   ex={"guides": GENES_FA}, up="Doench Rule Set 2 / benchling",
   tags=("CRISPR", "efficiency", "scoring"),
   summary="Rule-based activity score from GC, seed composition, poly-T terminator and hairpins.")
def crispr_grna_quality(guides, head=100, check_chromatin_like=False):
    """Guide efficiency scoring."""
    ls = [seq.clean(ln.strip(), "ACGTN").upper() for ln in io.as_text(guides).splitlines()
          if ln.strip() and not ln.startswith(">")]
    if not ls:
        ls = [r["seq"][:20] for r in _rows(guides)]
    scaffold = "GTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGC"
    rows = []
    for g in ls:
        g = g[:20]
        if len(g) < 20:
            continue
        gc = 100 * seq.gc_content(g)
        score = 50.0
        score += 1.2 * (gc - 50) if gc < 50 else -0.8 * (gc - 50) * 0.5
        score += 8 if g[3] == "G" else -4
        score += 5 if g[4] == "T" else -3
        score -= 25 if "TTTT" in g else 0
        score += 6 if g[-4] == "G" else 0
        score += 4 if g[-3] in "GC" else -4
        hp = seq.hairpin_score(g + scaffold[:12])
        hp_score = float(hp.get("score", 0.0)) if isinstance(hp, dict) else 0.0
        score -= min(20, hp_score / 2)
        if check_chromatin_like:
            score -= 5 if seq.gc_content(g[:6]) < 0.25 else 0
        rows.append({"guide": g, "GC_percent": round(gc, 2), "efficiency_score": round(score, 2),
                    "polyT_terminator": "TTTT" in g, "hairpin_score": round(hp_score, 2),
                    "u_at_minus_1": g[4] == "T", "g_at_minus_3": g[-3] in "G",
                    "with_scaffold": g + scaffold, "predicted": "active" if score >= 55 else "weak"})
    rows.sort(key=lambda r: -r["efficiency_score"])
    return table(rows[:int(head)], f"{len(rows)} guides scored")


@T("base_editor_targets", "Base-editing window inside protospacers", EDIT, "table",
   [anyfile("target", GENES_FA, fmt="", label="Target sequence"), choice("editor", ["C-to-T", "A-to-G"], "C-to-T", "Editor"),
    intin("window_start", 4, "Window start (position in protospacer)", min=1),
    intin("window_end", 8, "Window end", min=1), intin("head", 200, "Rows", min=1)],
   ex={"target": GENES_FA, "editor": "C-to-T"}, up="BE-Hive / CRISPResso design",
   tags=("base editing", "CRISPR", "design"),
   summary="List editable bases falling in the activity window for C or A base editors.")
def base_editor_targets(target, editor="C-to-T", window_start=4, window_end=8, head=200):
    """Editable bases inside the BE window."""
    body = "".join(r["seq"] for r in _rows(target)) or seq.clean(io.as_text(target), "ACGTN")
    base, to = ("C", "T") if editor.startswith("C") else ("A", "G")
    rows = []
    for strand, s in (("+", seq.clean(body.upper(), "ACGTN")),
                      ("-", seq.reverse_complement(seq.clean(body.upper(), "ACGTN")))):
        for i in range(0, max(0, len(s) - 23)):
            if not re.match("^[AG]G$", s[i + 21:i + 23] or ""):
                continue
            for pos in range(int(window_start), int(window_end) + 1):
                if pos - 1 >= len(s[i:i + 20]) or s[i + pos - 1] != base:
                    continue
                rows.append({"strand": strand, "protospacer": s[i:i + 20], "PAM": s[i + 20:i + 23],
                            "position_in_guide": pos, "reference_base": base, "edited_base": to,
                            "genomic_position": i + pos if strand == "+" else len(s) - (i + pos - 1),
                            "editor": editor, "context": s[max(0, i + pos - 4):i + pos + 4]})
    return table(rows[:int(head)], f"{len(rows)} edit opportunities for {editor}")


@T("hdr_template_design", "Design a donor with homology arms", EDIT, "file",
   [anyfile("target", GENES_FA, fmt="", label="Target sequence"), intin("cut_position", 60, "Cut site (1-based)", min=1),
    intin("arm_length", 500, "Homology arm length", min=20),
    bigtext("insert", "", "Insert sequence (empty = substitution only)"),
    intin("max_primer_tm", 68.0, "Max primer Tm", min=50.0)],
   ex={"target": GENES_FA, "cut_position": 60, "arm_length": 60, "insert": "ATGAAAGCT"},
   up="repair template design (DELIVERA / EuReCA)",
   tags=("CRISPR", "HDR", "donor design"),
   summary="Build left/right homology arms and the PCR primers to amplify a donor template.")
def hdr_template_design(target, cut_position=60, arm_length=500, insert="", max_primer_tm=68.0):
    """Homology arm + primer design."""
    body = "".join(r["seq"] for r in _rows(target)) or seq.clean(io.as_text(target), "ACGTN")
    cut = max(1, int(cut_position)) - 1
    arm = int(arm_length)
    left = body[max(0, cut - arm):cut]
    right = body[cut + 1:cut + 1 + arm]
    ins = seq.clean(insert.upper(), "ACGTN") if insert else ""
    donor = left + ins + right
    lines = [">donor", donor]
    primers = []
    if len(left) >= 20:
        fwd = left[-20:]
        rev = seq.reverse_complement(right[:20])
        primers = [f">left_primer\n{fwd}", f">right_primer\n{rev}"]
    body_txt = "\n".join(lines + primers) + "\n"
    tm_line = ""
    if len(left) >= 20:
        try:
            tm_line = f"left_primer_Tm: {round(float(seq.melting_temp(left[-20:])), 2)}\n"
        except Exception:  # noqa: BLE001
            tm_line = ""
    info = (f"# HDR donor design\nleft_arm_bp: {len(left)}\nright_arm_bp: {len(right)}\n"
            f"insert_bp: {len(ins)}\ndonor_bp: {len(donor)}\n"
            f"cut_site: {cut + 1}\nGC_left: {round(100 * seq.gc_content(left), 2)}\n"
            f"GC_right: {round(100 * seq.gc_content(right), 2)}\n" + tm_line)
    return {"text": body_txt + "\n" + info, "filename": "donor_design.txt",
            "message": f"donor of {len(donor)} bp with {len(left)}/{len(right)} bp arms"}


# ===========================================================================
# ChIP-seq style peak calling
# ===========================================================================
def _bedgraph(src):
    return io.parse_bedgraph(io.as_text(src))


def _peaks_from_graph(vals, threshold=1.0, min_width=1, max_gap=1, merge_distance=0):
    """Local-maximum peak caller on a bedGraph-style track."""
    per: dict[str, list] = defaultdict(list)
    for chrom, s, e, v in vals:
        per[chrom].append((int(s), int(e), float(v)))
    spans: list[tuple] = []
    for chrom, blocks in per.items():
        blocks.sort()
        cur = None
        for st, e, v in blocks:
            if v >= threshold:
                if cur and st - cur[2] <= max_gap:
                    cur = (cur[0], cur[1], e, cur[3] + [v])
                else:
                    if cur:
                        spans.append(cur)
                    cur = (chrom, st, e, [v])
            elif cur and st - cur[2] > max_gap:
                spans.append(cur)
                cur = None
        if cur:
            spans.append(cur)
    out = []
    for chrom, st, e, vs in spans:
        if e - st < min_width:
            continue
        out.append(io.Interval(chrom, st, e, f"peak_{len(out) + 1}", round(max(vs), 4), ".",
                              [f"sum={round(sum(vs), 4)}", f"bp={e - st}"]))
    return out


@T("macs2_callpeak_lite", "Call peaks from a signal track", CHIP, "file",
   [anyfile("graph", BEDGRAPH, fmt="", label="bedGraph / coverage"),
    number("cutoff", 5.0, "Fold-enrichment cutoff", min=0.0), intin("min_length", 50, "Minimum peak length", min=1),
    intin("max_gap", 100, "Merge blocks closer than", min=0), boolean("broad", False, "Broad peak mode"),
    fa("control", "", "Control genome (for local lambda)")],
   ex={"graph": BEDGRAPH, "cutoff": 2.0, "min_length": 50}, up="MACS2 callpeak",
   tags=("ChIP-seq", "peak calling", "coverage"),
   summary="Threshold and merge signal blocks into peaks with fold-enrichment and summit scores.")
def macs2_callpeak_lite(graph, cutoff=5.0, min_length=50, max_gap=100, broad=False, control=""):
    """bedGraph threshold peak caller."""
    vals = _bedgraph(graph)
    if not vals:
        return {"text": "", "filename": "peaks.bed", "message": "no bedGraph signal parsed"}
    gseq, _sz = _genome(control) if str(control or "").strip() else ({}, {})
    if gseq:
        total = sum(len(v) for v in gseq.values())
        local_lambda = max(0.1, total / 1000.0)
    else:
        local_lambda = 1.0
    thr = float(cutoff) * local_lambda
    peaks = _peaks_from_graph(vals, threshold=thr, min_width=int(min_length), max_gap=int(max_gap))
    if broad:
        peaks = genome.merge(peaks, distance=int(max_gap) * 10)
    body = io.write_bed(genome.sort_ivs(peaks), bed12=False)
    return {"text": body, "filename": "peaks.bed",
            "message": f"{len(peaks)} peaks at fold-enrichment >= {cutoff} (lambda={round(local_lambda, 3)})"}


@T("chipseq_qc_metrics", "Library QC: duplication, FRiP and strand metrics", CHIP, "table",
   [sam("alignments", ALN_SAM, "Aligned reads (SAM)"), bed("peaks", REGIONS_BED, "Peak BED"),
    number("read_length", 50, "Read length", min=1)],
   ex={"alignments": ALN_SAM, "peaks": REGIONS_BED}, up="phantompeakqualtools / deepTools plotFingerprint",
   tags=("ChIP-seq", "quality control"),
   summary="Read counts, duplication rate, fragments in peaks and cross-correlation proxy.")
def chipseq_qc_metrics(alignments, peaks, read_length=50):
    """ChIP library QC."""
    hdr, alns = bam.load(alignments)
    mapped = [a for a in alns if a.mapped]
    names = Counter(a.qname for a in mapped)
    dup = sum(v - 1 for v in names.values() if v > 1)
    ivs = _bed(peaks)
    in_peak = 0
    for a in mapped:
        p = a.pos
        if any(iv.chrom == a.rname and iv.start <= p - 1 < iv.end for iv in ivs):
            in_peak += 1
    rows = [{"metric": "total_records", "value": len(alns)},
            {"metric": "mapped", "value": len(mapped)},
            {"metric": "unique_positions", "value": len(names)},
            {"metric": "duplicates", "value": dup},
            {"metric": "duplicate_rate_percent", "value": round(100 * dup / max(1, len(mapped)), 3)},
            {"metric": "peaks", "value": len(ivs)},
            {"metric": "reads_in_peaks", "value": in_peak},
            {"metric": "FRiP", "value": round(in_peak / max(1, len(mapped)), 4)},
            {"metric": "mean_read_length_assumed", "value": float(read_length)},
            {"metric": "library_size_reads", "value": len(mapped) + dup}]
    return table(rows, f"FRiP {round(in_peak / max(1, len(mapped)), 3)} with "
                       f"{round(100 * dup / max(1, len(mapped)), 1)}% duplication")


@T("peak_gene_annotation", "Annotate peaks to the nearest TSS", CHIP, "table",
   [bed("peaks", REGIONS_BED, "Peak BED"), gff("annotation", ANNOT_GFF, "GFF3 annotation"),
    number("tss_distance", 2000, "Promoter window (bp)", min=0.0),
    intin("head", 200, "Rows", min=1)],
   ex={"peaks": REGIONS_BED, "annotation": ANNOT_GFF, "tss_distance": 2000},
   up="HOMER2 annotatePeaks / ChIPAnno", tags=("ChIP-seq", "annotation", "peaks"),
   summary="Classify each peak as promoter, intron, exon or intergenic with distance to the TSS.")
def peak_gene_annotation(peaks, annotation, tss_distance=2000, head=200):
    """Peak-to-gene annotation."""
    ivs = _bed(peaks)
    rows_gff = _gff(annotation)
    tss = genome.tss_list(rows_gff)
    genes = genome.genes_from_gff(rows_gff)
    near = genome.closest(ivs, tss)
    by_peak = defaultdict(list)
    for r in near:
        by_peak[r.get("name") or f"{r.get('chrom')}:{r.get('start')}"].append(r)
    out = []
    for iv in ivs:
        key = iv.name or f"{iv.chrom}:{iv.start}"
        cands = by_peak.get(key, [])
        c = cands[0] if cands else None
        dist = abs(int(c.get("distance", 0) or 0)) if c else None
        in_gene = any(g.chrom == iv.chrom and g.start <= iv.start < g.end for g in genes)
        gene_name = (c.get("other_name", "") if c else "") or ""
        if dist is not None and dist <= float(tss_distance):
            location = "promoter"
        elif in_gene:
            location = "intragenic"
        else:
            location = "distal_intergenic"
        out.append({"chrom": iv.chrom, "start": iv.start, "end": iv.end, "name": key,
                    "peak_length": iv.end - iv.start, "nearest_gene": gene_name,
                    "distance_to_TSS": dist if dist is not None else "", "location": location,
                    "score": _score(iv.score)})
    res = table(out[:int(head)], f"{len(out)} peaks annotated")
    res["stats"] = {"by_location": dict(Counter(r["location"] for r in out))}
    return res


@T("peak_grouper", "Merge and group overlapping peaks", CHIP, "table",
   [bed("peaks", REGIONS_BED, "Peak BED"), intin("distance", 0, "Merge within this distance", min=0),
    choice("collapse", ["concat", "sum", "max", "first"], "concat", "Name collapse"),
    intin("head", 200, "Rows", min=1)],
   ex={"peaks": REGIONS_BED, "distance": 50}, up="bedtools merge / MACS2 peak grouping",
   tags=("ChIP-seq", "intervals", "merging"),
   summary="Cluster peaks into loci, keeping member names, counts and summed signal.")
def peak_grouper(peaks, distance=0, collapse="concat", head=200):
    """Merge peaks into loci."""
    ivs = genome.sort_ivs(_bed(peaks))
    merged = genome.merge(ivs, distance=int(distance), collapse=collapse)
    clusters = genome.cluster(ivs, distance=int(distance))
    rows = []
    for i, m in enumerate(merged):
        c = clusters[i] if i < len(clusters) else {}
        rows.append({"chrom": m.chrom, "start": m.start, "end": m.end, "locus_length": m.end - m.start,
                    "n_peaks": c.get("n", 1) if isinstance(c, dict) else 1,
                    "names": (m.extra[0] if m.extra else str(m.name or "."))[:120],
                    "max_score": _score(m.score)})
    res = table(rows[:int(head)], f"{len(ivs)} peaks -> {len(merged)} loci "
                                 f"(gap <= {distance} bp)")
    return res


@T("differential_signal_regions", "Differential signal between two tracks", REG, "table",
   [anyfile("treatment", BEDGRAPH, fmt="", label="Treatment bedGraph"),
    anyfile("control", BEDGRAPH, fmt="", label="Control bedGraph"),
    intin("bin", 500, "Bin size", min=10), number("min_log2fc", 1.0, "Report |log2FC| above", min=0.0),
    intin("head", 200, "Rows", min=1)],
   ex={"treatment": BEDGRAPH, "control": GENOME_TXT, "bin": 1000}, up="diffReps / DESeq2 on peaks",
   tags=("differential binding", "coverage", "regional variation"),
   summary="Binned fold-change between two coverage tracks with per-bin p-values.")
def differential_signal_regions(treatment, control, bin=500, min_log2fc=1.0, head=200):
    """Binned differential signal."""
    def per_bin(vals):
        d: dict[tuple, float] = defaultdict(float)
        size = int(bin)
        for chrom, s, e, v in io.parse_bedgraph(io.as_text(vals)):
            b0 = int(s) // size
            for k in range(b0, int(e) // size + 1):
                lo, hi = max(s, k * size), min(e, (k + 1) * size)
                if hi > lo:
                    d[(chrom, k * size)] += float(v) * (hi - lo) / size
        return d
    a, b = per_bin(treatment), per_bin(control)
    rows = []
    for key in sorted(set(a) | set(b)):
        va, vb = a.get(key, 0.0), b.get(key, 0.0)
        fc = math.log2((va + 0.5) / (vb + 0.5))
        if abs(fc) < float(min_log2fc):
            continue
        orr, p = stats.fisher_exact(int(round(va * 10)) + 1, 10, int(round(vb * 10)) + 1, 10)
        rows.append({"chrom": key[0], "start": key[1], "end": key[1] + int(bin),
                    "treatment_signal": round(va, 4), "control_signal": round(vb, 4),
                    "log2_fold_change": round(fc, 4), "odds_ratio": round(float(orr), 4),
                    "p_value": p, "regulation": "up" if fc > 0 else "down"})
    rows.sort(key=lambda r: -abs(r["log2_fold_change"]))
    return table(rows[:int(head)], f"{len(rows)} differential bins (|log2FC| >= {min_log2fc})")


@T("peak_set_enrichment", "Test peak-linked genes for pathway enrichment", ANN, "table",
   [bed("peaks", REGIONS_BED, "Peak BED"), gff("annotation", ANNOT_GFF, "GFF3 annotation"),
    gmt_file := None or anyfile("gmt", GMT, fmt="", label="GMT gene sets"),
    number("min_q", 1.0, "Report q-value below", min=0.0)],
   ex={"peaks": REGIONS_BED, "annotation": ANNOT_GFF, "gmt": GMT},
   up="HOMER2 findGEnrichment / g:Profiler", tags=("enrichment", "pathways", "peaks"),
   summary="Fisher test of the genes near peaks against every gene set in a GMT file.")
def peak_set_enrichment(peaks, annotation, gmt, min_q=1.0):
    """Enrichment of genes bound at peaks."""
    ivs = _bed(peaks)
    genes = genome.genes_from_gff(_gff(annotation))
    near = genome.closest(ivs, genes)
    hit = {r.get("b_name") for r in near if r.get("b_name")}
    sets = io.parse_gmt(io.as_text(gmt))
    universe = {g.name for g in genes} or set().union(*(set(v) for v in sets.values())) if sets else set()
    n_uni, n_hit = max(1, len(universe)), max(1, len(hit & universe) if universe else len(hit))
    rows = []
    for name, members in sets.items():
        m = [x for x in members if not universe or x in universe]
        k = len(set(m) & hit)
        orr, p = stats.fisher_exact(k, len(m) - k, n_hit - k, n_uni - len(m) - n_hit + k)
        rows.append({"gene_set": name, "set_size": len(m), "hits": k,
                    "expected": round(len(m) * n_hit / n_uni, 3),
                    "fold_enrichment": round((k / max(1, len(m))) / max(1e-9, n_hit / n_uni), 3),
                    "p_value": p, "odds_ratio": round(float(orr), 4), "genes": ",".join(sorted(set(m) & hit)[:6])})
    ps = [r["p_value"] for r in rows]
    qs = stats.p_adjust(ps, "fdr_bh") if ps else []
    for r, q in zip(rows, qs):
        r["q_value"] = round(q, 6)
    rows = [r for r in rows if r["q_value"] <= float(min_q)] if rows and min_q < 1.0 else rows
    rows.sort(key=lambda r: r["p_value"])
    return table(rows, f"{len(rows)} gene sets tested against {n_hit} bound genes")
