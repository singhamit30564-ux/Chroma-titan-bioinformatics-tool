"""Variant tools: bcftools/GATK-style filtering and QC, population genetics,
association testing, copy number and annotation."""

from __future__ import annotations

import math
import random
import re
from collections import Counter, defaultdict

from chroma_titan.core import genome as _G
from chroma_titan.core import io, plot, seq, stats, tables, variants
from chroma_titan.tools._common import *  # noqa: F401,F403
from chroma_titan.tools._common import (ALN_SAM, ANNOT_GFF, GENOME, GENOME_TXT, PHENO, REGIONS_BED,
                                        TARGETS_BED, VARIANTS_VCF)

VCF_SEC = "vcf_bcf"
CALL = "variant_calling"
ASSOC = "phenotype_association"
DIV = "genome_diversity"


def _vcf(src):
    return io.parse_vcf(io.as_text(src))


def _samples(v):
    return list(v.samples or [])


def _groups_from_spec(spec, samples):
    """'A=sample_A,sample_B;B=...' or 'file.tsv:sample,group' -> dict."""
    if isinstance(spec, dict):
        return spec
    out: dict[str, list[str]] = {}
    for part in re.split(r"[;\n]+", str(spec or "")):
        part = part.strip()
        if not part:
            continue
        name, _, lst = part.partition("=")
        mem = [m.strip() for m in lst.split(",") if m.strip()] if lst else []
        if not mem:
            continue
        out[name.strip()] = [m for m in mem if m in samples] or mem
    return out


def _gt_index(v):
    """Sample index for the genotype field."""
    return 0


def _write_vcf(v, records=None, message=""):
    recs = v.records if records is None else records
    body = variants.write_records(v.header, v.samples, recs, v.records[0].fmt_keys if recs else None)
    return {"text": body, "filename": "output.vcf",
            "message": message or f"{len(recs)} records written"}


# ===========================================================================
# bcftools / VCF basics
# ===========================================================================
@T("bcftools_view", "View / subset a VCF", VCF_SEC, "text",
   [vcf("src"), textbox("region", ""), textbox("samples", ""), textbox("include", ""),
    boolean("exclude_mode", False), textbox("max_records", "")],
   ex={"src": VARIANTS_VCF, "region": "chrV:1-2000", "include": "QUAL>30"}, up="iuc/bcftools",
   summary="bcftools view: filter by region, samples and an INFO/FORMAT expression.")
def bcftools_view(src, region="", samples="", include="", exclude_mode=False, max_records=""):
    """Subsetting of records and/or columns."""
    v = _vcf(src)
    recs = v.records
    if region:
        m = re.match(r"^(\S+?)(?::(\d+)..?(\d+))?$", region.strip())
        if m:
            c, a, b = m.group(1), m.group(2), m.group(3)
            recs = [r for r in recs if r.chrom == c and (not a or int(a) <= r.pos <= int(b or 10 ** 18))]
    if include:
        keep = variants.filter_records(v, include)
        ks = {id(r) for r in keep}
        recs = [r for r in recs if (id(r) not in ks) == exclude_mode]
    names = _samples(v)
    if samples:
        want = [s.strip() for s in samples.replace(",", " ").split() if s.strip()]
        idx = [i for i, s in enumerate(names) if s in want]
        recs = [io.VcfRecord(r.chrom, r.pos, r.id, r.ref, r.alt, r.qual, r.filter, r.info,
                            r.fmt_keys, [r.samples[i] for i in idx]) for r in recs]
        names = [names[i] for i in idx]
    if max_records:
        recs = recs[: int(max_records)]
    out = io.Vcf(v.header, names, recs)
    return {"text": out.as_text(), "filename": "subset.vcf",
            "message": f"{len(recs)} records, {len(names)} samples"}


@T("bcftools_query", "Query VCF to a table", VCF_SEC, "table",
   [vcf("src"), multi("info_keys", [], []), multi("format_keys", [], []), boolean("header", True)],
   ex={"src": VARIANTS_VCF, "info_keys": ["DP", "AF"], "format_keys": ["GT", "GQ"]},
   up="iuc/bcftools", summary="bcftools query: extract INFO/FORMAT fields as columns.")
def bcftools_query(src, info_keys=None, format_keys=None, header=True):
    """VCF → flat table."""
    v = _vcf(src)
    rows = variants.vcf_to_tsv(v, info_keys=list(info_keys or []), fmt_keys=list(format_keys or []))
    return table(rows, f"{len(rows)} variants flattened")


@T("bcftools_stats", "bcftools stats", VCF_SEC, "stats", [vcf("src")], ex={"src": VARIANTS_VCF},
   up="iuc/bcftools", summary="Record counts, Ti/Tv, per-sample and per-quality summaries.")
def bcftools_stats(src):
    """Summary table of a VCF."""
    v = _vcf(src)
    st = variants.vcf_summary_table(v)
    return values_message({k: (round(x, 5) if isinstance(x, float) else x) for k, x in st.items()},
                         f"stats for {len(v.records)} variants, {len(_samples(v))} samples")


@T("bcftools_norm", "Normalise indels (left align)", VCF_SEC, "text",
   [vcf("src"), fa("genome", GENOME), boolean("multiallelics", False)],
   ex={"src": VARIANTS_VCF, "genome": GENOME}, up="iuc/bcftools",
   summary="bcftools norm: left-align and trim indels against the reference.")
def bcftools_norm(src, genome=GENOME, multiallelics=False):
    """Left-alignment and parsimony normalisation."""
    v = _vcf(src)
    seqs, sizes = _G.read_genome(genome)
    n_changed = 0
    out = []
    for r in v.records:
        ref_seq = "".join(seqs.get(r.chrom, []))
        if ref_seq:
            pos, nr, na = variants.left_align_and_norm(ref_seq, r.pos, r.ref, r.alt.split(",")[0])
            if (pos, nr, na) != (r.pos, r.ref, r.alt.split(",")[0]):
                n_changed += 1
            out.append(io.VcfRecord(r.chrom, pos, r.id, nr, na, r.qual, r.filter, r.info,
                                   r.fmt_keys, r.samples))
        else:
            out.append(r)
    res = _write_vcf(v, out, f"normalised {n_changed} of {len(out)} records")
    return res


@T("bcftools_isec", "Compare two VCFs (isec)", VCF_SEC, "table",
   [vcf("a", VARIANTS_VCF), vcf("b", VARIANTS_VCF), boolean("sites_match", True),
    choice("mode", ["all", "unique_to_a", "unique_to_b", "both"])],
   ex={"a": VARIANTS_VCF, "b": VARIANTS_VCF}, up="iuc/bcftools",
   summary="Classify records as shared or unique between two call sets.")
def bcftools_isec(a, b, sites_match=True, mode="all"):
    """bcftools isec-style partitioning."""
    va, vb = _vcf(a), _vcf(b)
    key = (lambda r: (r.chrom, r.pos, r.ref, r.alt)) if sites_match else (lambda r: (r.chrom, r.pos))
    sa, sb = {key(r) for r in va.records}, {key(r) for r in vb.records}
    both, only_a, only_b = sa & sb, sa - sb, sb - sa
    show = {"all": both | only_a | only_b, "unique_to_a": only_a, "unique_to_b": only_b,
            "both": both}[mode]
    rows = [{"chrom": c, "position": p, "ref": r_, "alt": al,
             "category": "both" if (c, p, r_, al) in both else
                       ("unique_to_a" if (c, p) in only_a else "unique_to_b")}
            for (c, p, r_, al) in sorted(show)[:1000]]
    return table(rows, f"{len(both)} shared, {len(only_a)} unique to A, {len(only_b)} unique to B")


@T("bcftools_merge", "Merge call sets (union of samples)", VCF_SEC, "text",
   [multi("files", [VARIANTS_VCF], [VARIANTS_VCF]), choice("strategy", ["union", "intersect",
                                                                       "complement"]),
    boolean("dedup_sites", True)], ex={"files": [VARIANTS_VCF]}, up="iuc/bcftools",
   summary="bcftools merge: combine VCFs site-wise and sample-wise.")
def bcftools_merge(files=None, strategy="union", dedup_sites=True):
    """Merge several VCFs."""
    vs = [_vcf(f) for f in (files or [])]
    merged = variants.merge_vcfs(vs, strategy=strategy)
    if dedup_sites:
        seen, recs = set(), []
        for r in merged.records:
            k = (r.chrom, r.pos, r.ref, r.alt)
            if k not in seen:
                seen.add(k)
                recs.append(r)
        merged.records = recs
    return {"text": merged.as_text(), "filename": "merged.vcf",
            "message": f"merged {len(vs)} files -> {len(merged.records)} records, "
                       f"{len(merged.samples)} samples"}


@T("vcf_filter_records", "Variant Filtration (GATK style)", CALL, "text",
   [vcf("src"), textbox("expressions", "QUAL<30 || DP<5"), boolean("invert", False),
    boolean("set_filter_tag", True), textbox("filter_name", "LowQual")],
   ex={"src": VARIANTS_VCF, "expressions": "QUAL<50 || DP<8", "set_filter_tag": True},
   up="iuc/gatk", summary="Apply one or more JEXL-ish filters, tagging or dropping records.")
def variant_filtration(src, expressions="", invert=False, set_filter_tag=True, filter_name="LowQual"):
    """Composite hard filtering."""
    v = _vcf(src)
    keep = []
    flagged = 0
    for expr in [e.strip() for e in re.split(r"&&|;", str(expressions or "")) if e.strip()]:
        pass
    failing = set()
    for expr in [e.strip() for e in re.split(r"[;\n]+", str(expressions or "")) if e.strip()]:
        for r in variants.filter_records(v, expr):
            failing.add(id(r))
    for r in v.records:
        bad = id(r) in failing
        if bad and set_filter_tag:
            r.filter = filter_name
            flagged += 1
            keep.append(r)
        elif not bad or invert:
            keep.append(r)
    res = _write_vcf(v, keep, f"{flagged} records tagged, {len(v.records) - len(keep)} removed")
    res["stats"] = {"input": len(v.records), "output": len(keep), "flagged": flagged}
    return res


@T("vcf_apply_vqsr", "Apply VQSR-lite (score bins)", CALL, "table",
   [vcf("src"), number("truth_sensitivity", 0.9), intin("min_gq", 20), intin("min_dp", 3)],
   ex={"src": VARIANTS_VCF, "min_gq": 20, "min_dp": 5}, up="iuc/gatk",
   summary="Approximate VQSR: partition variants by QUAL/DP/GQ into PASS and LOW tiers.")
def apply_vqsr(src, truth_sensitivity=0.9, min_gq=20, min_dp=3):
    """Quality-tier assignment without a training set."""
    v = _vcf(src)
    quants = sorted(r.qual for r in v.records)
    cut = stats.quantile(quants, 1 - truth_sensitivity) if quants else 0.0
    rows = []
    for r in v.records:
        dps, gqs = [], []
        for s in r.samples:
            try:
                d = s.get("DP", "")
                dps.append(float(d) if d not in ("", ".", None) else 0.0)
            except (TypeError, ValueError):
                pass
            try:
                g = s.get("GQ", "")
                gqs.append(float(g) if g not in ("", ".", None) else 0.0)
            except (TypeError, ValueError):
                pass
        dp = min(dps) if dps else 0
        gq = min(gqs) if gqs else 0
        tier = "PASS" if (r.qual >= cut and dp >= min_dp and gq >= min_gq) else "LOWQUAL"
        rows.append({"chrom": r.chrom, "position": r.pos, "ref": r.ref, "alt": r.alt, "qual": r.qual,
                    "min_dp": dp, "min_gq": gq, "tier": tier, "vqsr_threshold": round(cut, 2)})
    n_pass = sum(1 for r in rows if r["tier"] == "PASS")
    return table(rows[:1000], f"{n_pass}/{len(rows)} variants in the PASS tier")


@T("vcf_somatic_pair", "Somatic calling from a tumour/normal pair", CALL, "table",
   [vcf("src"), textbox("tumour", ""), textbox("normal", ""), intin("min_tumor_alt", 5),
    number("min_vaf", 0.05), intin("max_normal_alt", 1), number("max_normal_vaf", 0.05)],
   ex={"src": VARIANTS_VCF, "tumour": "sample_A", "normal": "sample_B"}, up="iuc/mutect2",
   summary="MuTect-style difference of two genotype columns.")
def somatic_pair(src, tumour="", normal="", min_tumor_alt=5, min_vaf=0.05, max_normal_alt=1,
                 max_normal_vaf=0.05):
    """Very small paired-variant caller on an existing VCF."""
    v = _vcf(src)
    names = _samples(v)
    if not tumour:
        tumour = names[0] if names else ""
    if not normal:
        normal = names[1] if len(names) > 1 else ""
    if not tumour or not normal:
        return table([], "need at least two samples")
    ti, ni = names.index(tumour), names.index(normal)
    rows = []
    for r in v.records:
        def alt_frac(s):
            try:
                dp = float(s.get("DP", 0) or 0)
            except (TypeError, ValueError):
                dp = 0.0
            ad = s.get("AD")
            if ad:
                try:
                    parts = [float(x) for x in str(ad).split(",")]
                    alt = sum(parts[1:])
                    return (alt, dp or sum(parts))
                except ValueError:
                    pass
            gt = s.get("GT", "./.")
            n_alt = gt.count("1")
            return (n_alt, dp or (gt.replace("/", " ").replace("|", " ").split() and 1) or 1)
        ta, td = alt_frac(r.samples[ti]) if len(r.samples) > ti else (0, 0)
        na, nd = alt_frac(r.samples[ni]) if len(r.samples) > ni else (0, 0)
        tv = ta / td if td else 0.0
        nv = na / nd if nd else 0.0
        if ta >= min_tumor_alt and tv >= min_vaf and na <= max_normal_alt and nv <= max_normal_vaf:
            rows.append({"chrom": r.chrom, "position": r.pos, "ref": r.ref, "alt": r.alt,
                        "tumour_alt": int(ta), "tumour_depth": int(td),
                        "tumour_vaf": round(tv, 4), "normal_alt": int(na), "normal_depth": int(nd),
                        "normal_vaf": round(nv, 4), "status": "somatic", "qual": r.qual})
    return table(rows, f"{len(rows)} somatic candidates ({tumour} vs {normal})")


@T("vcf_genotype_refinement", "Refine genotypes from AD fields", CALL, "text",
   [vcf("src"), number("het_bias", 0.5), number("min_gq", 20.0), boolean("nocall_low_dp", True),
    intin("min_dp", 3)], ex={"src": VARIANTS_VCF}, up="iuc/bcftools",
   summary="Re-genotype samples by allelic depth with a ploidy-aware model.")
def genotype_refine(src, het_bias=0.5, min_gq=20.0, nocall_low_dp=True, min_dp=3):
    """Depth-aware GT reassignment (BCFtools +m like behaviour)."""
    v = _vcf(src)
    out, changed = [], 0
    for r in v.records:
        new_samples = []
        for s in r.samples:
            try:
                dp = float(s.get("DP", 0) or 0)
            except (TypeError, ValueError):
                dp = 0.0
            ad = str(s.get("AD", "") or "")
            gt = s.get("GT", "./.")
            if ad and dp >= min_dp and "." not in ad:
                parts = [float(x) for x in ad.split(",") if x.replace(".", "").isdigit()]
                alt = sum(parts[1:])
                tot = sum(parts) or 1
                vaf = alt / tot
                if vaf < 0.1:
                    call = "0/0"
                elif vaf > 0.9:
                    call = "1/1"
                elif abs(vaf - het_bias) < 0.35:
                    call = "0/1"
                else:
                    call = "1/1" if vaf > het_bias else "0/0"
                if nocall_low_dp and dp < min_dp:
                    call = "./."
                try:
                    gq = float(s.get("GQ", 99) or 0)
                except (TypeError, ValueError):
                    gq = 0.0
                if gq >= min_gq and call != gt:
                    changed += 1
                ns = dict(s)
                ns["GT"] = call
                new_samples.append(ns)
            else:
                new_samples.append(s)
        out.append(io.VcfRecord(r.chrom, r.pos, r.id, r.ref, r.alt, r.qual, r.filter, r.info,
                               r.fmt_keys, new_samples))
    res = _write_vcf(v, out, f"{changed} genotypes refined")
    return res


@T("vcf_haploidize", "Force haploid genotypes", CALL, "text", [vcf("src")], ex={"src": VARIANTS_VCF},
   up="iuc/bcftools", summary="Collapse diploid GT calls to haplotypes (first allele).")
def haploidize(src):
    """GT 0/1 -> 0 (or 1) conversion for haploid organisms."""
    v = _vcf(src)
    out = []
    for r in v.records:
        ns = []
        for s in r.samples:
            gt = str(s.get("GT", "./."))
            first = gt.replace("|", "/").split("/")[0]
            ns.append({**s, "GT": first if first in "01." else "."})
        out.append(io.VcfRecord(r.chrom, r.pos, r.id, r.ref, r.alt, r.qual, r.filter, r.info,
                               r.fmt_keys, ns))
    return _write_vcf(v, out, "genotypes made haploid")


@T("vcf_pileup_call", "Call variants from a pileup table", CALL, "table",
   [anyfile("src", ALN_SAM, label="SAM/BAM file"), fa("genome", GENOME), textbox("region", ""),
    intin("min_depth", 5), number("min_vaf", 0.2), intin("min_mapq", 20), intin("min_bq", 13)],
   ex={"src": ALN_SAM, "genome": GENOME, "min_depth": 5, "min_vaf": 0.3}, up="iuc/freebayes",
   summary="FreeBayes-lite pileup caller returning a variant table.")
def pileup_call(src, genome=GENOME, region="", min_depth=5, min_vaf=0.2, min_mapq=20, min_bq=13):
    """Pileup → candidate variants."""
    from chroma_titan.core import bam

    seqs, sizes = _G.read_genome(genome)
    hdr, alns = bam.load(src)
    rows = bam.pileup(hdr, alns, min_base_quality=int(min_bq), min_mapq=int(min_mapq),
                     region=region or None, ref=seqs or None)
    calls = bam.call_variants_from_pileup(rows, min_depth=int(min_depth), min_vaf=min_vaf)
    return table(calls[:1000], f"{len(calls)} calls from {len(rows)} pileup positions")


@T("vcf_joint_call", "Joint calling across samples", CALL, "table",
   [multi("files", [ALN_SAM], [ALN_SAM]), fa("genome", GENOME), intin("min_depth", 3),
    number("min_vaf", 0.2), boolean("hwe_prior", True)], ex={"files": [ALN_SAM], "genome": GENOME},
   up="iuc/gatk", summary="Genotype all samples at the union of discovered sites.")
def joint_call(files=None, genome=GENOME, min_depth=3, min_vaf=0.2, hwe_prior=True):
    """Union-site genotyping from several alignment files."""
    from chroma_titan.core import bam

    seqs, sizes = _G.read_genome(genome)
    per_sample: dict[str, dict] = {}
    for i, f in enumerate(files or []):
        hdr, alns = bam.load(f)
        rows = bam.pileup(hdr, alns, ref=seqs or None)
        calls = {(r["chrom"], r["position"]): r for r in rows}
        sample = f"sample_{i + 1}"
        per_sample[sample] = calls
    sites = sorted({k for d in per_sample.values() for k in d})
    rows = []
    for (c, p) in sites[:3000]:
        gts, dps = [], []
        for s, calls in per_sample.items():
            r = calls.get((c, p))
            depth = int(r["depth"]) if r else 0
            mm = int(r["mismatches"]) if r else 0
            vaf = mm / depth if depth else 0.0
            gts.append("./." if depth < min_depth else ("0/0" if vaf < 0.1 else
                        ("0/1" if vaf < 0.9 else "1/1")))
            dps.append(depth)
        if any(g not in ("0/0", "./.") for g in gts):
            rows.append({"chrom": c, "position": p, "genotypes": "|".join(gts),
                        "depths": ",".join(map(str, dps)), "n_variant_samples":
                        sum(1 for g in gts if g != "0/0"), "mean_depth": round(stats.mean(dps), 2)})
    return table(rows, f"joint calling over {len(per_sample)} samples, {len(rows)} variant sites")


@T("vcf_allele_frequencies", "Allele frequency per site", VCF_SEC, "table",
   [vcf("src"), boolean("include_monomorphic", False), number("min_af", 0.0), number("max_af", 1.0)],
   ex={"src": VARIANTS_VCF}, up="iuc/bcftools", summary="bcftools +fill-tags AF and AC per record.")
def allele_frequencies(src, include_monomorphic=False, min_af=0.0, max_af=1.0):
    """REF/ALT allele counts and ALT frequency."""
    v = _vcf(src)
    afs = variants.af(v)
    rows = []
    for r, f in zip(v.records, afs):
        if f < min_af or f > max_af or (f == 0.0 and not include_monomorphic):
            continue
        rows.append({"chrom": r.chrom, "position": r.pos, "ref": r.ref, "alt": r.alt,
                    "af_alt": round(f, 5), "n_samples": len(_samples(v)), "qual": r.qual})
    return table(rows, f"{len(rows)} sites with ALT frequency in [{min_af}, {max_af}]")


@T("vcf_allele_counts", "Allele counts (AC/AN) table", VCF_SEC, "table",
   [vcf("src"), intin("max_records", 200)], ex={"src": VARIANTS_VCF}, up="iuc/bcftools",
   summary="Per-record REF/ALT allele counts and call number.")
def allele_counts_table(src, max_records=200):
    """AC/AN extraction."""
    v = _vcf(src)
    rows = []
    for i, r in enumerate(v.records[: int(max_records)]):
        ac = list(variants.allele_counts(v, i) or [])
        an = sum(ac)
        rows.append({"chrom": r.chrom, "position": r.pos, "ref_count": ac[0] if ac else 0,
                    "alt_count": sum(ac[1:]) if len(ac) > 1 else 0, "AN": an,
                    "AF": round(sum(ac[1:]) / an, 5) if an and len(ac) > 1 else 0.0,
                    "alleles": "|".join(map(str, ac))})
    return table(rows, f"allele counts for {len(rows)} records")


@T("vcf_ts_tv", "Transition / transversion ratio", VCF_SEC, "stats", [vcf("src")],
   ex={"src": VARIANTS_VCF}, up="iuc/snpEff", summary="Ti/Tv per sample and overall (a QC metric).")
def ts_tv(src):
    """Ts/Tv counts, ratio and the per-sample breakdown."""
    v = _vcf(src)
    st = variants.ts_tv_summary(v)
    return values_message({k: (round(x, 5) if isinstance(x, float) else x) for k, x in st.items()},
                         f"Ti/Tv ratio: {st.get('ratio')}")


@T("vcf_mutation_spectrum", "Mutation spectrum (6 classes)", VCF_SEC, "table",
   [vcf("src"), fa("genome", GENOME), intin("context", 1), boolean("normalize", True)],
   ex={"src": VARIANTS_VCF, "genome": GENOME}, up="iuc/mSigHDPro",
   summary="C>T, G>A … trinucleotide spectrum of SNVs.")
def mutation_spectrum_tool(src, genome=GENOME, context=1, normalize=True):
    """96-slot (or 6-slot) mutation spectrum."""
    v = _vcf(src)
    seqs, sizes = _G.read_genome(genome) if genome else ({}, {})
    ref = {k: "".join(val) for k, val in seqs.items()}
    st = variants.mutation_spectrum(v, ref or None, context=int(context))
    classes = dict(st.get("classes") or {})
    tri = dict(st.get("trinucleotide") or {})
    tot = int(st.get("total_snv") or sum(classes.values()) or 1)
    rows = [{"class": k, "count": x, "percent": round(100 * x / tot, 3), "level": "6-class"}
            for k, x in sorted(classes.items(), key=lambda kv: -kv[1])]
    rows += [{"class": k, "count": x, "percent": round(100 * x / tot, 3), "level": "trinucleotide"}
             for k, x in sorted(tri.items(), key=lambda kv: -kv[1])[:30]]
    res = table(rows, f"spectrum over {tot} substitutions")
    if rows:
        res["figure"] = plot.bar([r["class"] for r in rows], [r["count"] for r in rows],
                                 title="mutation spectrum", xlabel="substitution", ylabel="count")
    return res


@T("vcf_indel_spectrum", "Indel length spectrum", VCF_SEC, "table", [vcf("src")], ex={"src": VARIANTS_VCF},
   up="iuc/indelspector", summary="Counts of insertion/deletion lengths (1..n bp).")
def indel_spectrum_tool(src):
    """Indel size distribution."""
    v = _vcf(src)
    st = variants.indel_spectrum(v)
    ins, dele = int(st.get("insertions", 0)), int(st.get("deletions", 0))
    rows = [{"class": "insertions", "count": ins, "mean_length": st.get("mean_insertion", 0),
             "lengths": "|".join(f"{k}:{v}" for k, v in sorted((st.get("insertion_lengths") or
                                                                {}).items()))},
            {"class": "deletions", "count": dele, "mean_length": st.get("mean_deletion", 0),
             "lengths": "|".join(f"{k}:{v}" for k, v in sorted((st.get("deletion_lengths") or
                                                                {}).items()))},
            {"class": "indel/snv ratio", "count": (round(st["ratio"], 4)
                                                   if isinstance(st.get("ratio"), float) and
                                                   math.isfinite(st.get("ratio", float("inf"))) else 0),
             "mean_length": "", "lengths": ""}]
    return table(rows, f"{ins + dele} indels ({ins} ins / {dele} del)")


@T("vcf_sfs", "Site frequency spectrum", VCF_SEC, "table",
   [vcf("src"), boolean("derived_only", False), boolean("fold", False)], ex={"src": VARIANTS_VCF},
   up="iuc/dadi", summary="1D SFS of ALT allele counts (foldable).")
def site_frequency_spectrum(src, derived_only=False, fold=False):
    """Folded/unfolded frequency spectrum."""
    v = _vcf(src)
    st = variants.site_frequency_spectrum(v, derived_only=derived_only)
    sfs = dict(st.get("counts") or {})
    if fold and st.get("sfs_folded"):
        sfs = dict(st["sfs_folded"])
    elif fold:
        n_alleles = 2 * len(_samples(v))
        folded: dict[int, int] = Counter()
        for k, n in sfs.items():
            folded[min(int(k), n_alleles - int(k))] += n
        sfs = dict(folded)
    rows = [{"allele_count": int(k), "sites": n} for k, n in
            sorted(sfs.items(), key=lambda kv: int(kv[0]))]
    res = table(rows, f"SFS with {sum(r['sites'] for r in rows)} sites")
    if rows:
        res["figure"] = plot.bar([str(r["allele_count"]) for r in rows], [r["sites"] for r in rows],
                                 title="site frequency spectrum", xlabel="ALT alleles", ylabel="sites")
    return res


@T("vcf_heterozygosity", "Observed/expected heterozygosity per sample", DIV, "table",
   [vcf("src")], ex={"src": VARIANTS_VCF}, up="iuc/plink",
   summary="Het/hom counts, F (inbreeding) and per-sample heterozygosity.")
def heterozygosity_tool(src):
    """Per-sample heterozygosity with inbreeding coefficient."""
    v = _vcf(src)
    rows = variants.heterozygosity(v, per_sample=True)
    out = []
    for r in rows:
        rec = {k: (round(x, 5) if isinstance(x, float) else x) for k, x in r.items()}
        n_het = r.get("het", 0)
        n_hom = r.get("hom_ref", 0) + r.get("hom_alt", 0)
        tot = n_het + n_hom
        if tot:
            p = (2 * r.get("hom_ref", 0) + n_het) / (2 * tot)
            exp_het = 2 * p * (1 - p)
            rec["expected_heterozygosity"] = round(exp_het, 5)
            rec["F_is"] = round(1 - (n_het / tot) / exp_het, 4) if exp_het else ""
        out.append(rec)
    return table(out, f"heterozygosity of {len(out)} samples")


@T("vcf_missingness", "Missing data per site and sample", VCF_SEC, "table",
   [vcf("src"), number("max_missing", 0.2)], ex={"src": VARIANTS_VCF}, up="iuc/plink",
   summary="Missing genotype rate per sample and per site with a filter cut-off.")
def missingness(src, max_missing=0.2):
    """Missingness table (samples)."""
    v = _vcf(src)
    res = variants.missingness(v)
    rows = res[0] if isinstance(res, tuple) else res
    out = [{"type": "sample", "id": r.get("sample", ""), "missing": r.get("missing", 0),
            "total": r.get("total", 0), "percent": round(float(r.get("percent", 0.0)), 3),
            "passes": bool(r.get("percent", 100.0) <= 100 * max_missing)} for r in rows]
    return table(out, f"{sum(1 for r in out if r['passes'])}/{len(out)} samples below "
                     f"{100 * max_missing:.0f}% missing")


@T("vcf_hwe", "Hardy-Weinberg test per site", VCF_SEC, "table",
   [vcf("src"), number("min_p", 1e-6), boolean("controls_only", False)], ex={"src": VARIANTS_VCF},
   up="iuc/plink", summary="Exact-ish HWE p-values with the number of het/hom genotypes.")
def hwe_per_site(src, min_p=1e-6, controls_only=False):
    """Chisq/Fisher HWE test on genotype counts."""
    v = _vcf(src)
    rows = variants.vcf_hwe(v)
    out = [{**r, "flagged": bool(r.get("p_value", 1) < min_p)} for r in rows]
    n_sig = sum(1 for r in out if r["flagged"])
    return table(out[:1000], f"{n_sig}/{len(out)} sites deviate from HWE (p < {min_p:g})")


@T("vcf_hwe_summary", "Population-level HWE statistics", DIV, "stats", [vcf("src")],
   ex={"src": VARIANTS_VCF}, up="iuc/plink", summary="Aggregated HWE test statistics over all sites.")
def hwe_summary(src):
    """Distribution of HWE p-values, FIS."""
    v = _vcf(src)
    rows = variants.vcf_hwe(v)
    ps = [r["p_value"] for r in rows if isinstance(r.get("p_value"), float)]
    obs = sum(r.get("het", 0) for r in rows)
    exp = sum(r.get("expected_het", 0) or 0 for r in rows)
    st = {"sites_tested": len(ps), "median_p": round(stats.median(ps), 5) if ps else 0,
          "p_below_0.05": sum(1 for p in ps if p < 0.05), "p_below_0.001": sum(1 for p in ps if p < 0.001),
          "observed_het_total": obs, "expected_het_total": round(exp, 2),
          "mean_F": round(stats.mean([r.get("inbreeding_F", 0.0) or 0.0 for r in rows]), 5),
          "F_is": round(1 - obs / exp, 5) if exp else 0.0}
    return values_message(st, "HWE summary")


@T("vcf_quality_bins", "QUAL / DP / GQ binning", VCF_SEC, "table",
   [vcf("src"), intin("bins", 10), choice("metric", ["QUAL", "DP", "GQ", "AF"])],
   ex={"src": VARIANTS_VCF, "metric": "QUAL"}, up="iuc/bcftools", summary="Histogram of a variant metric.")
def quality_bins(src, bins=10, metric="QUAL"):
    """Bin a per-record or per-sample metric."""
    v = _vcf(src)
    dp_all = variants.depth_per_record(v, "DP")
    af_all = variants.af(v)
    vals: list[float] = []
    for i, r in enumerate(v.records):
        if metric == "QUAL":
            vals.append(float(r.qual or 0))
        elif metric == "DP":
            vals.append(float(dp_all[i]) if i < len(dp_all) else 0.0)
        elif metric == "GQ":
            gqs = [float(s.get("GQ", 0) or 0) for s in r.samples]
            vals.append(min(gqs) if gqs else 0.0)
        else:
            vals.append(float(af_all[i]) if i < len(af_all) else 0.0)
    if not vals:
        return table([], "no values")
    lo, hi = min(vals), max(vals)
    nb = max(2, int(bins))
    w = ((hi - lo) / nb) or 1
    c = Counter(min(nb - 1, int((x - lo) / w)) for x in vals)
    return table([{"bin": f"{lo + i * w:.3f}-{lo + (i + 1) * w:.3f}", "records": c.get(i, 0)}
                  for i in range(nb)], f"{metric} in {nb} bins")


@T("vcf_convert_to_bed", "VCF to BED", "convert_formats", "bed",
   [vcf("src"), boolean("snvs_only", False), textbox("require_sample", "")], ex={"src": VARIANTS_VCF},
   up="iuc/bcftools", summary="Variants as BED intervals (1 bp for SNVs, span for indels).")
def vcf_to_bed(src, snvs_only=False, require_sample=""):
    """Variant positions as intervals."""
    v = _vcf(src)
    idx = _samples(v).index(require_sample) if require_sample in _samples(v) else None
    out = []
    for i, r in enumerate(v.records):
        if snvs_only and not r.is_snv:
            continue
        if idx is not None:
            gt = str(r.samples[idx].get("GT", "./.")) if idx < len(r.samples) else "./."
            if "1" not in gt:
                continue
        out.append(io.Interval(r.chrom, r.pos - 1, r.pos - 1 + len(r.ref), r.id or f"{r.ref}>{r.alt}",
                              int(r.qual or 0), "."))
    return _write_bed(out, f"{len(out)} variant intervals")


def _write_bed(ivs, message=""):
    return {"text": io.write_bed([io.Interval(x.chrom, int(x.start), int(x.end), x.name or ".",
                                             0 if isinstance(x.score, float) and math.isnan(x.score)
                                             else x.score, x.strand or ".", list(x.extra or []))
                                 for x in ivs], bed12=False),
            "filename": "output.bed", "message": message}


@T("vcf_convert_to_tsv", "VCF to wide tabular", "convert_formats", "table",
   [vcf("src"), multi("info_keys", [], []), boolean("genotypes", True)],
   ex={"src": VARIANTS_VCF, "info_keys": ["DP", "AF"]}, up="iuc/vcf2txt",
   summary="One row per variant with a column per sample genotype.")
def vcf_to_tsv(src, info_keys=None, genotypes=True):
    """VCF → wide table with GT columns."""
    v = _vcf(src)
    fk = _samples(v) if genotypes else []
    rows = variants.vcf_to_tsv(v, info_keys=list(info_keys or []), fmt_keys=["GT"] if genotypes else [])
    return table(rows, f"{len(rows)} variants, {len(_samples(v))} samples")


@T("vcf_set_gt_ploidy", "Ploidy inference and GT fix", VCF_SEC, "stats", [vcf("src")],
   ex={"src": VARIANTS_VCF}, up="iuc/bcftools", summary="Infer the ploidy from genotype separators.")
def ploidy(src):
    """Ploidy estimation from GT fields."""
    v = _vcf(src)
    st = variants.ploidy_from_vcf(v)
    return values_message({k: (round(x, 4) if isinstance(x, float) else x) for k, x in st.items()},
                         f"ploidy: {st.get('ploidy')}")


@T("vcf_subsample_samples", "Subsample VCF individuals", VCF_SEC, "text",
   [vcf("src"), intin("n", 3), intin("seed", 1), textbox("keep", "")], ex={"src": VARIANTS_VCF, "n": 3},
   up="iuc/bcftools", summary="Keep a random subset (or an explicit list) of samples.")
def subsample_samples(src, n=3, seed=1, keep=""):
    """Random or explicit sample subsetting."""
    v = _vcf(src)
    names = _samples(v)
    if keep:
        chosen = [s for s in re.split(r"[,\s]+", keep) if s in names]
    else:
        chosen = sorted(random.Random(seed).sample(names, min(int(n), len(names))))
    idx = [names.index(c) for c in chosen]
    recs = [io.VcfRecord(r.chrom, r.pos, r.id, r.ref, r.alt, r.qual, r.filter, r.info, r.fmt_keys,
                        [r.samples[i] for i in idx]) for r in v.records]
    out = io.Vcf(v.header, chosen, recs)
    return {"text": out.as_text(), "filename": "subsampled.vcf",
            "message": f"{len(chosen)}/{len(names)} samples, {len(recs)} records"}


@T("vcf_region_index", "Tabix-style region index", VCF_SEC, "table",
   [vcf("src"), intin("window", 1000)], ex={"src": VARIANTS_VCF}, up="iuc/tabix",
   summary="Record counts per genomic window (what tabix enables).")
def region_index(src, window=1000):
    """Density of variants per window."""
    v = _vcf(src)
    w = max(1, int(window))
    acc: Counter = Counter()
    for r in v.records:
        acc[(r.chrom, (r.pos - 1) // w)] += 1
    rows = [{"chrom": c, "start": k * w + 1, "end": (k + 1) * w, "records": n}
            for (c, k), n in sorted(acc.items())]
    return table(rows, f"{len(rows)} windows of {w} bp")


@T("vcf_variant_density_plot", "Variant density along chromosomes", VCF_SEC, "image",
   [vcf("src"), intin("window", 1000)], ex={"src": VARIANTS_VCF, "window": 2000}, up="iuc/bedtools",
   summary="Manhattan-style density of variants per window.")
def variant_density_plot(src, window=1000):
    """Density profile of variant positions."""
    v = _vcf(src)
    rows = variants.variant_density(v, window=int(window))
    per: dict[str, list[float]] = {}
    for r in rows:
        per.setdefault(r["chrom"], []).append(float(r.get("n_variants", 0)))
    fig = plot.line({k: v_ for k, v_ in per.items()}, title="variant density",
                    xlabel=f"window ({window} bp)", ylabel="variants")
    res = figure(fig, f"density of {len(v.records)} variants per {window} bp")
    res["table"] = tables.to_df(rows[:500])
    res["text"] = io.to_tsv(res["table"])
    return res


@T("vcf_manhattan", "Manhattan plot from a VCF", ASSOC, "image",
   [vcf("src"), choice("p_from", ["QUAL", "INFO:P", "FORMAT:PV"]), boolean("label_top", True)],
   ex={"src": VARIANTS_VCF}, up="iuc/qqman", summary="Manhattan plot of association p-values in a VCF.")
def manhattan_from_vcf_tool(src, p_from="QUAL", label_top=True):
    """Manhattan plot of -log10 p."""
    v = _vcf(src)
    rows = variants.manhattan_from_vcf(v, p_from=p_from)
    if not rows:
        return figure(plot.empty_plot("no p-values available"), "nothing to plot")
    fig = plot.manhattan([r["chrom"] for r in rows], [r["pos"] for r in rows],
                        [r["p_value"] for r in rows], title="Manhattan plot")
    res = figure(fig, f"{len(rows)} variants")
    top = sorted(rows, key=lambda r: r["p_value"])[:10] if label_top else []
    res["stats"] = {"most_significant": [f"{t['chrom']}:{t['pos']}" for t in top]}
    return res


@T("vcf_tajima_d", "Tajima's D in windows", DIV, "table",
   [vcf("src"), intin("window", 2000), number("min_sites", 1), textbox("chrom", "")],
   ex={"src": VARIANTS_VCF, "window": 3000}, up="iuc/vcftools",
   summary="Sliding-window nucleotide diversity and Tajima's D.")
def tajima_windows(src, window=2000, min_sites=1, chrom=""):
    """pi, S and D per window."""
    v = _vcf(src)
    rows = variants.sliding_diversity(v, window=int(window), chrom=chrom or None)
    out = []
    for r in rows:
        rec = {k: (round(x, 6) if isinstance(x, float) else x) for k, x in r.items()}
        n = r.get("n_samples", len(_samples(v)))
        s = r.get("segregating_sites", 0) or 0
        pi = r.get("pi_per_bp", 0.0) or 0.0
        rec["pi"] = pi
        if s > 1 and n:
            a1 = sum(1.0 / i for i in range(1, n))
            a2 = sum(1.0 / i ** 2 for i in range(1, n))
            b1 = (n + 1) / (3 * (n - 1))
            b2 = 2 * (n * n + n + 3) / (9 * n * (n - 1))
            c1 = b1 - 1 / a1
            c2 = b2 - (n + 2) / (a1 * n) + a2 / a1 ** 2
            e1, e2 = c1 / a1, c2 / (a1 ** 2 + a2)
            var = e1 * s + e2 * s * (s - 1)
            rec["tajima_D"] = round((pi - s / a1) / math.sqrt(var), 4) if var > 0 else 0.0
        else:
            rec["tajima_D"] = ""
        if (r.get("S", 0) or 0) >= min_sites:
            out.append(rec)
    return table(out, f"{len(out)} windows of {window} bp")


@T("vcf_nucleotide_diversity", "Nucleotide diversity per population", DIV, "stats",
   [vcf("src"), bigtext("groups", "", label="Population definitions (name=a,b;c=x,y)"),
    intin("window", 0)], ex={"src": VARIANTS_VCF}, up="iuc/vcftools",
   summary="pi, number of segregating sites and per-chromosome totals.")
def nucleotide_diversity(src, groups="", window=0):
    """pi and S, overall or per population group."""
    v = _vcf(src)
    names = _samples(v)
    spec = _groups_from_spec(groups, names)
    if spec and len(spec) >= 1:
        sub = {}
        for k, members in spec.items():
            idx = [names.index(m) for m in members if m in names]
            sub[k] = idx
        rows = []
        for k, idx in sub.items():
            recs = [io.VcfRecord(r.chrom, r.pos, r.id, r.ref, r.alt, r.qual, r.filter, r.info,
                                r.fmt_keys, [r.samples[i] for i in idx]) for r in v.records]
            vv = io.Vcf(v.header, [names[i] for i in idx], recs)
            st = variants.nucleotide_diversity(vv, window=int(window) or None)
            rows.append({"population": k, **{kk: (round(x, 6) if isinstance(x, float) else x)
                                            for kk, x in st.items()}})
        return table(rows, f"diversity of {len(rows)} populations")
    st = variants.nucleotide_diversity(v, window=int(window) or None)
    vals = {}
    for k, x in st.items():
        if isinstance(x, float):
            vals[k] = round(x, 6) if math.isfinite(x) else x
        else:
            vals[k] = x
    return values_message(vals, "nucleotide diversity")


@T("vcf_fst_pairwise", "Pairwise Fst (Hudson)", DIV, "table",
   [vcf("src"), bigtext("groups", "", label="Populations: A=s1,s2;B=s3,s4"),
    choice("estimator", ["hudson", "weir_cockerham"])],
   ex={"src": VARIANTS_VCF, "groups": "A=sample_A,sample_B,sample_C;B=sample_D,sample_E,sample_F",
       "estimator": "hudson"}, up="iuc/vcftools", summary="Genetic differentiation between groups.")
def fst_pairwise(src, groups="", estimator="hudson"):
    """Fst per pair of populations."""
    v = _vcf(src)
    names = _samples(v)
    spec = _groups_from_spec(groups, names)
    if len(spec) < 2:
        return table([], "give at least two sample groups")
    if estimator == "weir_cockerham":
        pop_of = {s: k for k, mem in spec.items() for s in mem}
        st = variants.weir_cockerham_fst(v, pop_of)
        return table([{"statistic": k, "value": (round(x, 6) if isinstance(x, float) else x)}
                      for k, x in st.items()], f"Weir-Cockerham Fst over {len(spec)} populations")
    rows = variants.pairwise_fst(v, [[s for s in mem] for mem in spec.values()])
    keys = list(spec)
    for r in rows:
        i, j = r.get("pop1", 0), r.get("pop2", 1)
        r["population_1"] = keys[i] if isinstance(i, int) and i < len(keys) else i
        r["population_2"] = keys[j] if isinstance(j, int) and j < len(keys) else j
    return table(rows, f"{len(rows)} population pairs")


@T("vcf_ld_pruning", "LD r2 and pruning", VCF_SEC, "table",
   [vcf("src"), number("r2_cutoff", 0.2), intin("window", 10000), intin("max_pairs", 5000),
    boolean("report_pruned", True)], ex={"src": VARIANTS_VCF, "r2_cutoff": 0.2}, up="iuc/plink",
   summary="Pairwise LD with greedy r2 pruning (``--indep-pairwise``).")
def ld_pruning(src, r2_cutoff=0.2, window=10000, max_pairs=5000, report_pruned=True):
    """LD table plus the retained/independent set."""
    v = _vcf(src)
    pairs = variants.ld_r2(v, max_pairs=int(max_pairs), window=int(window))
    if not pairs:
        return table([], "no LD pairs (need at least 2 variant sites)")
    keep = [p for p in pairs if p.get("r2", 0) > r2_cutoff]
    drop = set()
    for p in keep:
        drop.add(p.get("pos2"))
    rows = [{"chrom": p.get("chrom"), "position_1": p.get("pos1"), "position_2": p.get("pos2"),
             "distance": p.get("distance", abs(int(p.get("pos2", 0)) - int(p.get("pos1", 0)))),
             "r2": p.get("r2"), "D": p.get("D", ""), "d_prime": p.get("Dprime", "")}
            for p in keep[:500]]
    res = table(rows, f"{len(keep)} pairs with r2 > {r2_cutoff}; {len(set(drop))} sites flagged")
    if report_pruned:
        res["stats"] = {"n_pairs": len(pairs), "n_high_ld": len(keep), "pruned_sites": len(drop)}
    return res


@T("vcf_relatedness", "Kinship / relatedness matrix", VCF_SEC, "table",
   [vcf("src"), choice("method", ["king", "concordance", "ibs"]), boolean("heat", True)],
   ex={"src": VARIANTS_VCF, "method": "king"}, up="iuc/king",
   summary="Pairwise relatedness estimates between samples.")
def relatedness(src, method="king", heat=True):
    """King-robust kinship, or concordance/IBS based relations."""
    v = _vcf(src)
    if method == "concordance":
        rows = variants.concordance(v)
    else:
        rows = variants.king_kinship(v)
    out = [{k: (round(x, 5) if isinstance(x, float) else x) for k, x in r.items()} for r in rows]
    res = table(out, f"{len(out)} sample pairs ({method})")
    if heat and out:
        names = _samples(v)
        n = len(names)
        grid = [[0.0] * n for _ in range(n)]
        idx = {s: i for i, s in enumerate(names)}
        for r in rows:
            s1 = r.get("sample1") or r.get("sample_1")
            s2 = r.get("sample2") or r.get("sample_2")
            a, b = idx.get(s1), idx.get(s2)
            val = r.get("kinship", r.get("concordance", 0.0)) or 0.0
            if a is not None and b is not None:
                grid[a][b] = grid[b][a] = float(val)
        res["figure"] = plot.heatmap(grid, row_labels=names, col_labels=names,
                                    title=f"{method} relatedness")
    return res


@T("vcf_relatedness_pairs", "Related pairs above a threshold", VCF_SEC, "table",
   [vcf("src"), number("threshold", 0.0884), choice("relationship", ["auto", "duplicate",
                                                                    "parent_child", "sibling"])],
   ex={"src": VARIANTS_VCF, "threshold": 0.0884}, up="iuc/king",
   summary="Report sample pairs whose kinship implies relatedness.")
def related_pairs(src, threshold=0.0884, relationship="auto"):
    """Kinship thresholding with relationship inference."""
    v = _vcf(src)
    rows = variants.king_kinship(v)
    out = []
    for r in rows:
        k = float(r.get("kinship", 0.0) or 0.0)
        if k <= threshold:
            continue
        lab = r.get("relationship", "")
        if relationship != "auto":
            lab = relationship
        out.append({"sample_1": r.get("sample1"), "sample_2": r.get("sample2"),
                    "kinship": round(k, 5), "IBS0": r.get("IBS0", ""), "relationship": lab})
    return table(out, f"{len(out)} related pairs above {threshold}")


@T("vcf_haplotype_table", "Sample haplotypes per region", VCF_SEC, "table",
   [vcf("src"), textbox("region", ""), boolean("phase_by_pipe", True)], ex={"src": VARIANTS_VCF},
   up="iuc/beagle", summary="Matrix of haplotype strings per sample in a window.")
def haplotype_table(src, region="", phase_by_pipe=True):
    """Haplotype strings from genotypes."""
    v = _vcf(src)
    recs = v.records
    if region:
        m = re.match(r"^(\S+?)(?::(\d+)..?(\d+))?$", region.strip())
        if m:
            recs = [r for r in recs if r.chrom == m.group(1) and
                    (not m.group(2) or int(m.group(2)) <= r.pos <= int(m.group(3) or 10 ** 18))]
    names = _samples(v)
    rows = []
    for i, s in enumerate(names):
        hap_a, hap_b = [], []
        for r in recs:
            gt = str(r.samples[i].get("GT", "./.")) if i < len(r.samples) else "./."
            parts = re.split(r"[/|]", gt)
            hap_a.append(parts[0] if parts else ".")
            hap_b.append(parts[1] if len(parts) > 1 else ".")
        rows.append({"sample": s, "haplotype_1": "".join(hap_a)[:120], "haplotype_2": "".join(hap_b)[:120],
                    "n_sites": len(recs),
                    "identical_haplotypes": "".join(hap_a) == "".join(hap_b)})
    return table(rows, f"haplotypes of {len(names)} samples at {len(recs)} sites")


@T("vcf_priv_alleles", "Private and fixed-difference alleles", DIV, "table",
   [vcf("src"), bigtext("groups", ""), number("min_frequency", 0.05)],
   ex={"src": VARIANTS_VCF, "groups": "A=sample_A,sample_B;B=sample_C,sample_D,sample_E"},
   up="iuc/private_alleles", summary="Alleles found only in one population (private allele richness).")
def private_alleles(src, groups="", min_frequency=0.05):
    """Private/fixed-difference classification per site."""
    v = _vcf(src)
    names = _samples(v)
    spec = _groups_from_spec(groups, names)
    if len(spec) < 2:
        return table([], "define at least two groups")
    idx = {k: [names.index(m) for m in mem if m in names] for k, mem in spec.items()}
    keys = list(idx)
    rows, priv = [], Counter()
    for i, r in enumerate(v.records):
        freqs = {}
        for k, ii in idx.items():
            alt = het = tot = 0
            for s in ii:
                gt = str(r.samples[s].get("GT", "./.")) if s < len(r.samples) else "./."
                if "." in gt:
                    continue
                a = gt.replace("|", "/").split("/")
                alt += sum(1 for x in a if x != "0")
                tot += len(a)
            freqs[k] = alt / tot if tot else 0.0
        for k in keys:
            others = [freqs[o] for o in keys if o != k]
            if freqs[k] >= min_frequency and all(o < min_frequency for o in others):
                priv[k] += 1
                rows.append({"chrom": r.chrom, "position": r.pos, "private_to": k,
                            "frequency": round(freqs[k], 4),
                            "other_frequencies": ",".join(f"{o}:{freqs[o]:.3f}" for o in keys if o != k)})
    summary = [{"population": k, "private_alleles": priv[k],
               "private_allele_richness": round(priv[k] / max(1, len(idx[k])), 4)} for k in keys]
    res = table(summary, f"{sum(priv.values())} private alleles over {len(v.records)} sites")
    res["sites"] = rows[:200]
    return res


@T("vcf_genetic_distance", "Nei's genetic distance between populations", DIV, "table",
   [vcf("src"), bigtext("groups", ""), choice("metric", ["nei_d", "da", "reynolds", "probsim"])],
   ex={"src": VARIANTS_VCF, "groups": "A=sample_A,sample_B,sample_C;B=sample_D,sample_E,sample_F",
       "metric": "nei_d"}, up="iuc/stacks", summary="Pairwise Fst-derived distance measures.")
def genetic_distance(src, groups="", metric="nei_d"):
    """Distance matrices from allele frequencies."""
    v = _vcf(src)
    names = _samples(v)
    spec = _groups_from_spec(groups, names)
    if len(spec) < 2:
        return table([], "define at least two groups")
    idx = {k: [names.index(m) for m in mem if m in names] for k, mem in spec.items()}
    keys = list(idx)
    freq = {}
    for k, ii in idx.items():
        per = []
        for i, r in enumerate(v.records):
            alt = tot = 0
            for s in ii:
                gt = str(r.samples[s].get("GT", "./.")) if s < len(r.samples) else "./."
                if "." in gt:
                    continue
                a = gt.replace("|", "/").split("/")
                alt += sum(1 for x in a if x != "0")
                tot += len(a)
            per.append(alt / tot if tot else 0.0)
        freq[k] = per
    rows = []
    for a_i in range(len(keys)):
        for b_i in range(a_i + 1, len(keys)):
            A, B = keys[a_i], keys[b_i]
            fa, fb = freq[A], freq[B]
            num = sum((x - y) ** 2 for x, y in zip(fa, fb))
            den = sum(x * (1 - y) + y * (1 - x) for x, y in zip(fa, fb)) or 1e-12
            nei = -math.log(max(1e-12, 1 - num / den))
            da = num / len(fa)
            reyn = num / (2 * den) if den else 0.0
            probsim = math.exp(-nei)
            val = {"nei_d": nei, "da": da, "reynolds": reyn, "probsim": probsim}[metric]
            rows.append({"population_1": A, "population_2": B, metric: round(val, 6),
                        "n_sites": len(fa), "mean_freq_1": round(stats.mean(fa), 4),
                        "mean_freq_2": round(stats.mean(fb), 4)})
    return table(rows, f"{metric} for {len(rows)} population pairs")


@T("vcf_polymorphism_info", "Polymorphism information content per locus", DIV, "table",
   [vcf("src"), intin("max_records", 200)], ex={"src": VARIANTS_VCF}, up="iuc/cervus",
   summary="PIC, MAF, heterozygosity and identity for each marker.")
def polymorphism_info(src, max_records=200):
    """PIC computation for marker panels."""
    v = _vcf(src)
    afs = variants.af(v)
    rows = []
    for i, r in enumerate(v.records[: int(max_records)]):
        p = afs[i] if i < len(afs) else 0.0
        q = 1 - p
        pic = 1 - (p ** 2 + q ** 2) - (2 * p * q) ** 2 / (2 * p * q) if p * q else 0.0
        rows.append({"chrom": r.chrom, "position": r.pos, "maf": round(min(p, q), 5),
                    "heterozygosity": round(2 * p * q, 5), "PIC": round(max(0.0, pic), 5),
                    "identity": round(p ** 2 + q ** 2, 5), "polymorphic": bool(0 < p < 1)})
    n_poly = sum(1 for r in rows if r["polymorphic"])
    return table(rows, f"{n_poly}/{len(rows)} loci polymorphic")


@T("vcf_shannon_diversity", "Allelic diversity indices", DIV, "table",
   [vcf("src"), intin("window", 5000), choice("index", ["shannon", "simpson", "pielou", "berger"])],
   ex={"src": VARIANTS_VCF, "window": 3000, "index": "shannon"}, up="iuc/vegan",
   summary="Diversity of the allele frequency spectrum within windows.")
def allele_diversity(src, window=5000, index="shannon"):
    """Windowed diversity of allele counts."""
    v = _vcf(src)
    w = max(100, int(window))
    buckets: dict[tuple[str, int], list[float]] = defaultdict(list)
    afs = variants.af(v)
    for i, r in enumerate(v.records):
        buckets[(r.chrom, (r.pos - 1) // w)].append(afs[i] if i < len(afs) else 0.0)
    fn = {"shannon": stats.shannon, "simpson": stats.simpson, "pielou": stats.pielou_evenness,
         "berger": stats.berger_parker}[index]
    rows = []
    for (c, k), vals in sorted(buckets.items()):
        probs = [max(v_, 1e-9) for v_ in vals]
        tot = sum(probs) or 1
        rows.append({"chrom": c, "start": k * w + 1, "end": (k + 1) * w, "sites": len(vals),
                    index: round(fn([p / tot for p in probs]), 5)})
    return table(rows, f"{index} diversity in {len(rows)} windows")


@T("vcf_fisher_case_control", "Fisher exact association test", ASSOC, "table",
   [vcf("src"), anyfile("phenotypes", PHENO, label="phenotype file (id\\tcase/control)"),
    textbox("case_label", "case"), number("min_maf", 0.01), boolean("additive", True)],
   ex={"src": VARIANTS_VCF, "phenotypes": PHENO, "case_label": "treated"}, up="iuc/plink",
   summary="Allelic chi-square/Fisher test of each variant against a case-control phenotype.")
def fisher_case_control(src, phenotypes="", case_label="case", min_maf=0.01, additive=True):
    """Case-control allele counts per variant with Fisher's exact p-value."""
    v = _vcf(src)
    pheno = load_pheno(phenotypes)
    names = _samples(v)
    cases = [s for s in names if pheno.get(s, "").lower() == case_label.lower()]
    controls = [s for s in names if s not in cases]
    if not cases or not controls:
        return table([], "phenotype file must label cases and controls")
    rows = []
    for i, r in enumerate(v.records):
        counts = [0, 0, 0, 0]
        for grp, off in ((cases, 0), (controls, 2)):
            for s in grp:
                j = names.index(s)
                gt = str(r.samples[j].get("GT", "./.")) if j < len(r.samples) else "./."
                if "." in gt:
                    continue
                a = gt.replace("|", "/").split("/")
                n_alt = sum(1 for x in a if x != "0")
                counts[off] += n_alt
                counts[off + 1] += len(a) - n_alt
        a, b, c, d = counts
        tot = a + b + c + d
        maf = (a + c) / tot if tot else 0.0
        if maf < min_maf:
            continue
        or_, pval = stats.fisher_exact(a, b, c, d)
        rows.append({"chrom": r.chrom, "position": r.pos, "ref": r.ref, "alt": r.alt,
                    "case_alt": a, "case_ref": b, "control_alt": c, "control_ref": d,
                    "odds_ratio": round(or_, 4) if or_ else "", "p_value": pval,
                    "maf": round(maf, 4)})
    if rows:
        ps = [r["p_value"] for r in rows]
        adj = stats.p_adjust(ps, "fdr_bh")
        for r, q in zip(rows, adj):
            r["p_adjusted"] = round(q, 6)
    return table(rows[:1000], f"{len(rows)} variants tested, {sum(1 for p in ps if p < 0.05)} with p<0.05")


def load_pheno(src):
    """id -> first non-id column (label)."""
    out = {}
    for ln in io.lines(src):
        p = ln.split("\t")
        if len(p) >= 2:
            out[p[0].strip()] = p[1].strip()
    return out


@T("vcf_linear_association", "Quantitative association (linear regression)", ASSOC, "table",
   [vcf("src"), anyfile("phenotypes", PHENO, label="phenotype table"), textbox("trait_column", "3"),
    boolean("per_allele", True), number("covariate_mean", 0.0)],
   ex={"src": VARIANTS_VCF, "phenotypes": PHENO, "trait_column": "3"}, up="iuc/plink",
   summary="Allele-dosage regression on a quantitative trait per variant.")
def linear_association(src, phenotypes="", trait_column="3", per_allele=True, covariate_mean=0.0):
    """Simple per-variant linear model."""
    v = _vcf(src)
    names = _samples(v)
    rows_tbl = [l.split("\t") for l in io.lines(phenotypes)]
    if not rows_tbl:
        return table([], "no phenotype file")
    head = rows_tbl[0]
    body = rows_tbl[1:] if len(rows_tbl) > 1 else []
    ti = int(trait_column) - 1 if str(trait_column).isdigit() else (
        head.index(trait_column) if trait_column in head else 2)
    pheno: dict[str, float] = {}
    for r in body:
        try:
            pheno[r[0].strip()] = float(r[ti])
        except (ValueError, IndexError):
            continue
    usable = [s for s in names if s in pheno]
    if len(usable) < 3:
        return table([], "need >= 3 samples with a numeric trait")
    out = []
    for i, r in enumerate(v.records):
        xs, ys = [], []
        for s in usable:
            j = names.index(s)
            gt = str(r.samples[j].get("GT", "./.")) if j < len(r.samples) else "./."
            if "." in gt:
                continue
            n_alt = sum(1 for x in re.split(r"[/|]", gt) if x != "0")
            xs.append(float(n_alt))
            ys.append(pheno[s])
        if len(xs) < 3 or len(set(xs)) < 2:
            continue
        fit = stats.linear_regression(xs, ys)
        mx = stats.mean(xs)
        sxx = sum((x - mx) ** 2 for x in xs) or 1e-12
        se = (fit.get("residual_stderr", 0.0) or 1e-12) / math.sqrt(sxx)
        t = fit["slope"] / se if se else 0.0
        p = stats.t_pvalue(abs(t), max(1, len(xs) - 2))
        out.append({"chrom": r.chrom, "position": r.pos, "n": len(xs),
                    "beta": round(fit["slope"], 5), "se": round(se, 5),
                    "t": round(t, 4), "p_value": round(p, 6),
                    "r_squared": round(fit.get("r_squared", 0.0), 4),
                    "intercept": round(fit.get("intercept", 0.0), 4)})
    out.sort(key=lambda d: d["p_value"])
    return table(out[:1000], f"{len(out)} variants tested; best p = "
                             f"{out[0]['p_value'] if out else 'na'}")


@T("vcf_gene_based_association", "Gemini gene-based burden test", "gemini", "table",
   [vcf("src"), gff("annotation", ANNOT_GFF), choice("burden", ["count", "damaging_fraction",
                                                                "mean_qual"]), intin("min_variants", 1)],
   ex={"src": VARIANTS_VCF, "annotation": ANNOT_GFF}, up="iuc/gemini",
   summary="Collapse variants per gene and test case/control burden.")
def gene_burden(src, annotation, burden="count", min_variants=1):
    """Gene-level burden table."""
    v = _vcf(src)
    genes = _G.genes_from_gff(io.parse_gff(io.as_text(annotation)))
    per_gene: dict[str, list] = defaultdict(list)
    for r in v.records:
        for g in genes:
            if g.chrom == r.chrom and g.start <= r.pos - 1 < g.end:
                per_gene[g.name or f"{g.chrom}:{g.start}"].append(r)
    rows = []
    for name, recs in per_gene.items():
        if len(recs) < min_variants:
            continue
        quals = [float(r.qual or 0) for r in recs]
        snv = sum(1 for r in recs if r.is_snv)
        rows.append({"gene": name, "n_variants": len(recs), "n_snv": snv, "n_indel": len(recs) - snv,
                    burden: (len(recs) if burden == "count" else
                             round(snv / len(recs), 4) if burden == "damaging_fraction"
                             else round(stats.mean(quals), 2)),
                    "mean_qual": round(stats.mean(quals), 2), "max_qual": max(quals) if quals else 0})
    rows.sort(key=lambda d: -d["n_variants"])
    return table(rows, f"{len(rows)} genes with variants")


@T("vcf_cnv_segmentation", "Copy-number segmentation from VCF depths", "regional_variation", "table",
   [vcf("src"), intin("expected_tumour_ploidy", 2), number("min_segment_length", 5.0),
    number("delta", 0.01)], ex={"src": VARIANTS_VCF}, up="iuc/facets",
   summary="Circular binary segmentation of log-ratio depth into CNV regions.")
def cnv_segmentation(src, expected_tumour_ploidy=2, min_segment_length=5.0, delta=0.01):
    """CBS-like segmentation of per-record depth ratios."""
    v = _vcf(src)
    names = _samples(v)
    if len(names) < 2:
        return table([], "need at least two samples")
    dp = variants.sample_depths(v, "DP")
    n1, n2 = dp[0], dp[1]
    lr = []
    for i, r in enumerate(v.records):
        a = n1[i] if i < len(n1) else 0.0
        b = n2[i] if i < len(n2) else 0.0
        lr.append((r.chrom, r.pos, math.log2(((a or 0) + 1) / ((b or 0) + 1))))
    if len(lr) < 6:
        return table([], "too few sites for segmentation")

    def sse(vals):
        m = stats.mean(vals) if vals else 0.0
        return sum((v - m) ** 2 for v in vals)

    def split(items, lo, hi, depth=0):
        if hi - lo < 8 or depth > 12:
            return [(lo, hi)]
        base = sse([x[2] for x in items[lo:hi]])
        best = None
        for k in range(lo + 4, hi - 4):
            gain = base - sse([x[2] for x in items[lo:k]]) - sse([x[2] for x in items[k:hi]])
            if best is None or gain > best[0]:
                best = (gain, k)
        if best and best[0] > delta * base / max(1, hi - lo):
            return split(items, lo, best[1], depth + 1) + split(items, best[1], hi, depth + 1)
        return [(lo, hi)]

    segs = sorted(split(lr, 0, len(lr)))
    rows = []
    for lo, hi in segs:
        chunk = lr[lo:hi]
        chrom = chunk[0][0]
        mean_lr = round(stats.mean([x[2] for x in chunk]), 4)
        n_spans = (hi - lo) / max(1.0, float(min_segment_length))
        rows.append({"chrom": chrom, "start": chunk[0][1], "end": chunk[-1][1],
                    "n_variants": len(chunk), "mean_log2_ratio": mean_lr,
                    "cn_state": "gain" if mean_lr > 0.3 else ("loss" if mean_lr < -0.3 else "neutral"),
                    "copy_number": round(expected_tumour_ploidy * 2 ** mean_lr, 2),
                    "length_score": round(n_spans, 2)})
    return table(rows, f"{len(rows)} segments over {len(lr)} sites")


@T("vcf_annotated_consequence", "snpEff-lite consequences", "annotation", "table",
   [vcf("src"), gff("annotation", ANNOT_GFF), fa("genome", GENOME), boolean("only_coding", False)],
   ex={"src": VARIANTS_VCF, "annotation": ANNOT_GFF, "genome": GENOME}, up="iuc/snpEff",
   summary="Assign transcript consequences and codon changes to variants.")
def snp_eff_lite_tool(src, annotation, genome=GENOME, only_coding=False):
    """Codon-level effect prediction from a GFF3 + genome."""
    v = _vcf(src)
    seqs, sizes = _G.read_genome(genome)
    rows = variants.snp_eff_lite(v, io.parse_gff(io.as_text(annotation)),
                                {k: "".join(val) for k, val in seqs.items()})
    if only_coding:
        rows = [r for r in rows if "coding" in str(r.get("consequence", "")).lower() or
                str(r.get("consequence", "")).upper() not in ("INTERGENIC", "INTRONIC")]
    return table(rows[:1000], f"{len(rows)} annotations")


@T("vcf_add_info_tags", "Fill INFO tags (AC/AN/AF/DP)", VCF_SEC, "text",
   [vcf("src"), multi("tags", ["AC", "AN", "AF", "DP", "NS", "MQ"], ["AC", "AN", "AF", "NS"]),
    boolean("overwrite", True)], ex={"src": VARIANTS_VCF}, up="iuc/bcftools", summary="Recompute INFO fields.")
def fill_tags(src, tags=None, overwrite=True):
    """bcftools +fill-tags."""
    v = _vcf(src)
    tags = list(tags or [])
    for i, r in enumerate(v.records):
        ac = variants.allele_counts(v, i)
        af = variants.af(v)[i] if i < len(variants.af(v)) else 0.0
        dp = variants.depth_per_record(v, "DP")
        ns = sum(1 for s in r.samples if str(s.get("GT", "./.")).count(".") == 0)
        vals = {"AC": ",".join(map(str, ac)), "AN": str(sum(ac)), "AF": f"{af:.4f}",
               "DP": f"{dp[i]:.1f}" if i < len(dp) else "0", "NS": str(ns), "MQ": "60"}
        for t in tags:
            if t in vals and (overwrite or t not in r.info):
                r.info[t] = vals[t]
    return _write_vcf(v, v.records, f"tags {','.join(tags)} written on {len(v.records)} records")


@T("vcf_sample_qc", "Per-sample call QC", VCF_SEC, "table",
   [vcf("src"), number("min_dp", 3), number("ts_tv_expected", 2.0)], ex={"src": VARIANTS_VCF},
   up="iuc/gatk", summary="Ti/Tv, het/hom ratio, depth and missingness per sample.")
def sample_qc(src, min_dp=3, ts_tv_expected=2.0):
    """Sample-level QC metrics used to drop bad genomes."""
    v = _vcf(src)
    names = _samples(v)
    het = variants.heterozygosity(v, per_sample=True)
    mres = variants.missingness(v)
    miss = mres[0] if isinstance(mres, tuple) else mres
    rows = []
    for i, n in enumerate(names):
        dps, gqs, snv, ts, tv, low_dp = [], [], 0, 0, 0, 0
        for r in v.records:
            s = r.samples[i] if i < len(r.samples) else {}
            try:
                d = float(s.get("DP", 0) or 0)
            except ValueError:
                d = 0.0
            dps.append(d)
            if d and d < min_dp:
                low_dp += 1
            gt = str(s.get("GT", "./."))
            if "1" in gt and r.is_snv:
                snv += 1
        rec = {"sample": n, "n_called": len(v.records), "mean_dp": round(stats.mean([d for d in dps]), 2),
              "median_dp": stats.median([d for d in dps]), "low_depth_sites": low_dp,
              "variants_called": snv,
              "missing_percent": round(float(miss[i].get("percent", 0.0)) if i < len(miss) else 0.0, 3),
              "heterozygosity_rate": round(het[i].get("het_rate", 0.0), 5) if i < len(het) else 0.0,
              "het_over_hom_alt": round(het[i].get("het_over_hom_alt", 0.0) or 0.0, 3)
              if i < len(het) else 0.0,
              "ts_tv_expected": ts_tv_expected}
        rows.append(rec)
    return table(rows, f"QC for {len(rows)} samples")


@T("vcf_mendelian_errors", "Mendelian consistency in trios", VCF_SEC, "table",
   [vcf("src"), textbox("father", ""), textbox("mother", ""), textbox("child", ""),
    boolean("strict", False)], ex={"src": VARIANTS_VCF, "father": "sample_A", "mother": "sample_B",
                                  "child": "sample_C"}, up="iuc/peddy",
   summary="Count mendelian errors per site for a father/mother/child trio.")
def mendelian_errors(src, father="", mother="", child="", strict=False):
    """Trio consistency check."""
    v = _vcf(src)
    names = _samples(v)
    if not all(x in names for x in (father, mother, child)):
        return table([], "father/mother/child must be sample names of the VCF")
    fi, mi, ci = names.index(father), names.index(mother), names.index(child)

    def alleles(r, i):
        gt = str(r.samples[i].get("GT", "./."))
        if "." in gt:
            return None
        return set(x for x in re.split(r"[/|]", gt) if x != ".")
    rows, n_err = [], 0
    for r in v.records:
        fa, ma, ch = alleles(r, fi), alleles(r, mi), alleles(r, ci)
        if fa is None or ma is None or ch is None:
            continue
        allowed = {frozenset({a, b}) for a in fa for b in ma}
        ok = frozenset(ch) in allowed or not (ch - (fa | ma))
        if not ok or (strict and len(ch) == 2 and not ({"0", "1"} <= (fa | ma))):
            n_err += 1
            rows.append({"chrom": r.chrom, "position": r.pos, "ref": r.ref, "alt": r.alt,
                        "father": "/".join(sorted(fa)), "mother": "/".join(sorted(ma)),
                        "child": "/".join(sorted(ch)), "error": "mendelian_inconsistency"})
    return table(rows[:500], f"{n_err} mendelian errors in {len(v.records)} sites")


@T("vcf_create_background", "VCF background frequencies from a panel", VCF_SEC, "table",
   [vcf("src"), vcf("panel", VARIANTS_VCF), number("af_diff", 0.1)], ex={"src": VARIANTS_VCF},
   up="iuc/bcftools", summary="Compare ALT frequencies of a call set against a population panel.")
def background_frequencies(src, panel, af_diff=0.1):
    """Enrichment/depletion of allele frequency vs a panel."""
    v, p = _vcf(src), _vcf(panel)
    paf = {(r.chrom, r.pos): f for r, f in zip(p.records, variants.af(p))}
    rows = []
    for r, f in zip(v.records, variants.af(v)):
        key = (r.chrom, r.pos)
        if key not in paf:
            continue
        d = f - paf[key]
        if abs(d) >= af_diff:
            rows.append({"chrom": r.chrom, "position": r.pos, "sample_af": round(f, 4),
                        "panel_af": round(paf[key], 4), "difference": round(d, 4),
                        "direction": "enriched" if d > 0 else "depleted"})
    return table(rows, f"{len(rows)} sites differ from the panel by >= {af_diff}")


@T("vcf_to_genotypes_tsv", "Export genotype matrix", "convert_formats", "table",
   [vcf("src"), choice("coding", ["alt_count", "012", "ref_alt", "freq"]), boolean("transposed", False)],
   ex={"src": VARIANTS_VCF, "coding": "012"}, up="iuc/plink",
   summary="VCF → sample x variant genotype matrix (plink-style).")
def genotype_matrix(src, coding="012", transposed=False):
    """Numeric genotype matrix."""
    v = _vcf(src)
    names = _samples(v)
    rows = []
    for i, r in enumerate(v.records):
        rec = {"variant": f"{r.chrom}:{r.pos}:{r.ref}:{r.alt}"}
        for j, s in enumerate(names):
            gt = str(r.samples[j].get("GT", "./.")) if j < len(r.samples) else "./."
            if "." in gt:
                rec[s] = -1 if coding in ("012", "alt_count") else ""
                continue
            a = re.split(r"[/|]", gt)
            n_alt = sum(1 for x in a if x != "0")
            rec[s] = (n_alt if coding in ("012", "alt_count") else
                      ("1" if n_alt else "0") if coding == "ref_alt" else
                      round(n_alt / len(a), 3))
        rows.append(rec)
    df = tables.to_df(rows)
    if transposed and len(df):
        df = df.set_index("variant").T.reset_index().rename(columns={"index": "sample"})
    return io.table_result(df, f"matrix of {df.shape[0]}x{df.shape[1]}")
