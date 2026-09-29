"""Variant/VCF engines: genotype parsing and population-genetic statistics."""

from __future__ import annotations

import math
import re
from collections import Counter, defaultdict
from typing import Sequence

from chroma_titan.core import stats
from chroma_titan.core.io import Vcf, VcfRecord, parse_vcf


def genotypes(vcf: Vcf, i: int) -> list[tuple[int, int]]:
    """Allele index pairs for sample ``i`` (missing -> (-1, -1))."""
    out = []
    for r in vcf.records:
        gt = r.gt(i)
        alleles = re.split(r"[/|]", gt)
        try:
            a, b = int(alleles[0]), int(alleles[1]) if len(alleles) > 1 else int(alleles[0])
        except (ValueError, IndexError):
            out.append((-1, -1))
            continue
        if a < 0 or b < 0:
            out.append((-1, -1))
        else:
            out.append((a, b))
    return out


def all_genotypes(vcf: Vcf) -> list[list[tuple[int, int]]]:
    return [genotypes(vcf, i) for i in range(len(vcf.samples))]


def allele_counts(vcf: Vcf, i: int) -> list[int]:
    return [a + b for a, b in genotypes(vcf, i) if a >= 0 and b >= 0]


def af(vcf: Vcf) -> list[float]:
    """ALT allele frequency per record (from GT)."""
    gts = all_genotypes(vcf)
    out = []
    for k in range(len(vcf.records)):
        tot = alt = 0
        for g in gts:
            a, b = g[k]
            if a < 0:
                continue
            tot += 2
            alt += (a == 1) + (b == 1)
        out.append(alt / tot if tot else float("nan"))
    return out


def depth_per_record(vcf: Vcf, tag: str = "DP") -> list[float]:
    out = []
    for r in vcf.records:
        if tag == "DP" and "DP" in r.info:
            out.append(r.info_num("DP"))
            continue
        vals = [float(r.field(i, tag)) for i in range(len(r.samples))
                if _isnum(r.field(i, tag))]
        out.append(sum(vals) if vals else (r.info_num("DP", 0.0)))
    return out


def sample_depths(vcf: Vcf, tag: str = "DP") -> list[list[float]]:
    out = []
    for i in range(len(vcf.samples)):
        row = []
        for r in vcf.records:
            v = r.field(i, tag)
            row.append(float(v) if _isnum(v) else 0.0)
        out.append(row)
    return out


# ---------------------------------------------------------------------------
# summary statistics
# ---------------------------------------------------------------------------
def ts_tv_summary(vcf: Vcf) -> dict:
    n_ts = n_tv = n_snv = n_indel = n_mnv = n_other = 0
    by_type: Counter = Counter()
    for r in vcf.records:
        if r.filter not in ("PASS", ".", ""):
            continue
        if r.is_indel:
            n_indel += 1
            by_type["INDEL"] += 1
            continue
        if not r.is_snv:
            n_mnv += 1
            by_type["MNV"] += 1
            continue
        n_snv += 1
        if r.is_transition:
            n_ts += 1
            by_type["Ti"] += 1
        else:
            n_tv += 1
            by_type["Tv"] += 1
        for a in r.alts:
            by_type[f"{r.ref}>{a}"] += 1
    return {"transitions": n_ts, "transversions": n_tv, "snvs": n_snv, "indels": n_indel,
            "mnvs": n_mnv, "ratio_ts_tv": round(n_ts / n_tv, 4) if n_tv else float("inf"),
            "ratio_indel_snv": round(n_indel / n_snv, 4) if n_snv else float("nan"),
            "total": len(vcf.records), "type_counts": dict(by_type)}


def mutation_spectrum(vcf: Vcf, genome: dict[str, str] | None = None,
                      context: int = 1) -> dict:
    """96-channel trinucleotide mutation spectrum + the 6 pyrimidine classes."""
    cnt: Counter = Counter()
    tri: Counter = Counter()
    for r in vcf.records:
        if not r.is_snv:
            continue
        for a in r.alts:
            pair = (r.ref, a)
            if pair in (("A", "G"), ("T", "C")):
                cnt["A>G"] += 1
            elif pair in (("C", "T"), ("G", "A")):
                cnt["C>T"] += 1
            elif pair in (("A", "C"), ("T", "G")):
                cnt["A>C"] += 1
            elif pair in (("C", "A"), ("G", "T")):
                cnt["C>A"] += 1
            elif pair in (("A", "T"), ("T", "A")):
                cnt["A>T"] += 1
            else:
                cnt["C>G"] += 1
        if genome:
            seqs = genome
            ch = seqs.get(r.chrom, "") if isinstance(seqs, dict) else ""
            if ch and 0 < r.pos <= len(ch) - 1:
                ctx = ch[r.pos - 2:r.pos + 1].upper()
                if len(ctx) == 3:
                    key = (ctx[0], r.ref, r.alts[0] if r.alts else "N", ctx[2])
                    tri[f"{key[0]}[{key[1]}>{key[2]}]{key[3]}"] += 1
    del genome
    return {"classes": dict(cnt), "trinucleotide": dict(tri),
            "total_snv": sum(cnt.values())}


def indel_spectrum(vcf: Vcf) -> dict:
    ins = delete = 0
    lengths_ins: list[int] = []
    lengths_del: list[int] = []
    homopolymers = 0
    for r in vcf.records:
        if r.is_snv or not r.is_indel:
            continue
        for a in r.alts:
            d = len(a) - len(r.ref)
            if d > 0:
                ins += 1
                lengths_ins.append(d)
            elif d < 0:
                delete += 1
                lengths_del.append(-d)
    return {"insertions": ins, "deletions": delete,
            "ratio": round(ins / delete, 4) if delete else float("inf"),
            "mean_insertion": round(stats.mean(lengths_ins), 3) if lengths_ins else 0,
            "mean_deletion": round(stats.mean(lengths_del), 3) if lengths_del else 0,
            "insertion_lengths": Counter(lengths_ins), "deletion_lengths": Counter(lengths_del),
            "homopolymer_nearby": homopolymers}


def site_frequency_spectrum(vcf: Vcf, derived_only: bool = False) -> dict:
    gts = all_genotypes(vcf)
    n = len(vcf.samples) or 1
    bins = Counter()
    for k, r in enumerate(vcf.records):
        tot = alt = 0
        for g in gts:
            a, b = g[k]
            if a < 0:
                continue
            tot += 2
            alt += (a > 0) + (b > 0)
        if not tot:
            continue
        d = alt if not derived_only else min(alt, tot - alt)
        bins[d] += 1
    maf = [min(v / (2 * n), 1 - v / (2 * n)) for v in bins]
    return {"counts": dict(sorted(bins.items())), "n_bins": len(bins),
            "singletons": bins.get(1, 0), "doubletons": bins.get(2, 0),
            "private_alleles": sum(v for k, v in bins.items() if k == 1),
            "min_af": round(min(maf), 6) if maf else 0.0, "sfs_folded": dict(sorted(bins.items()))}


def heterozygosity(vcf: Vcf, per_sample: bool = True) -> list[dict]:
    out = []
    for i, name in enumerate(vcf.samples or ["ALL"]):
        gts = genotypes(vcf, i) if vcf.samples else all_genotypes(vcf)[0] if vcf.samples else []
        het = hom_ref = hom_alt = miss = 0
        for a, b in gts:
            if a < 0:
                miss += 1
            elif a == b:
                if a == 0:
                    hom_ref += 1
                else:
                    hom_alt += 1
            else:
                het += 1
        n = max(1, len(gts))
        out.append({"sample": name, "heterozygous": het, "homozygous_ref": hom_ref,
                    "homozygous_alt": hom_alt, "missing": miss,
                    "het_rate": round(het / n, 6),
                    "het_over_hom_alt": round(het / hom_alt, 4) if hom_alt else float("inf"),
                    "n_sites": len(gts)})
    return out


def missingness(vcf: Vcf) -> list[dict]:
    rows = []
    for i, s in enumerate(vcf.samples):
        g = genotypes(vcf, i)
        miss = sum(1 for a, b in g if a < 0)
        rows.append({"sample": s, "missing": miss, "total": len(g),
                     "percent": round(100 * miss / max(1, len(g)), 3)})
    per_site = []
    for k, r in enumerate(vcf.records):
        cnt = sum(1 for i in range(len(vcf.samples)) if genotypes(vcf, i)[k][0] < 0)
        if cnt:
            per_site.append({"chrom": r.chrom, "pos": r.pos, "missing_samples": cnt})
    return rows, per_site


def hwe_test(freq: float, n_het: int, n_hom_ref: int, n_hom_alt: int) -> dict:
    """Weir-Smith exact test style chi-square HWE check (diploid)."""
    n = n_het + n_hom_ref + n_hom_alt
    if n == 0:
        return {"p_value": float("nan")}
    p = (2 * n_hom_ref + n_het) / (2 * n)
    e_ref, e_het, e_alt = p * p * n, 2 * p * (1 - p) * n, (1 - p) ** 2 * n
    chi = 0.0
    for obs, exp in ((n_hom_ref, e_ref), (n_het, e_het), (n_hom_alt, e_alt)):
        if exp:
            chi += (obs - exp) ** 2 / exp
    fis = 1 - (n_het / n) / (2 * p * (1 - p)) if p not in (0, 1) and n else float("nan")
    return {"expected_het": round(e_het, 3), "chi_square": round(chi, 4),
            "p_value": stats.chi2_pvalue(chi, 1), "inbreeding_F": round(fis, 4),
            "allele_freq": round(p, 4), "n": n}


def vcf_hwe(vcf: Vcf) -> list[dict]:
    out = []
    gts = all_genotypes(vcf)
    for k, r in enumerate(vcf.records):
        hr = sum(1 for g in gts if g[k] == (0, 0))
        het = sum(1 for g in gts if g[k][0] != g[k][1] and g[k][0] >= 0)
        ha = sum(1 for g in gts if g[k][0] == g[k][1] and g[k][0] > 0)
        res = hwe_test(0.5, het, hr, ha)
        out.append({"chrom": r.chrom, "pos": r.pos, "hom_ref": hr, "het": het,
                    "hom_alt": ha, **res})
    return out


def nucleotide_diversity(vcf: Vcf, window: int | None = None, genome_sizes: dict | None = None) -> dict:
    """pi, Watterson theta and Tajima's D over the whole file (or per window)."""
    gts = all_genotypes(vcf)
    n = len(gts)
    if not n:
        return {"pi": 0.0, "theta_w": 0.0, "tajima_D": 0.0, "S": 0}
    sepsites = 0
    pi_total = 0.0
    pos_list = []
    for k in range(len(vcf.records)):
        cnt = Counter()
        tot = 0
        for g in gts:
            a, b = g[k]
            if a < 0:
                continue
            cnt[a] += 1
            cnt[b] += 1
            tot += 2
        if tot < 2 or len(cnt) < 2:
            continue
        sepsites += 1
        pos_list.append(vcf.records[k].pos)
        pi_total += tot / (tot - 1) * (1 - sum((c / tot) ** 2 for c in cnt.values()))
    L = window or (max(genome_sizes.values()) if genome_sizes else 1) or 1
    m = max(2, 2 * n)  # Tajima D is defined over the number of haploid sequences
    a1 = sum(1.0 / i for i in range(1, m))
    a2 = sum(1.0 / i ** 2 for i in range(1, m))
    theta_w = sepsites / a1
    pi = pi_total
    b1, b2 = (m + 1) / (3 * (m - 1)), 2 * (m * m + m + 3) / (9 * m * (m - 1))
    c1, c2 = b1 - 1 / a1, b2 - (m + 2) / (a1 * m) + a2 / a1 ** 2
    e1, e2 = c1 / a1, c2 / (a1 ** 2 + a2)
    var = e1 * sepsites + e2 * sepsites * (sepsites - 1) if m > 1 else 0
    d = (pi - theta_w) / math.sqrt(var) if var > 0 else 0.0
    return {"pi_per_bp": round(pi / L, 8) if L else 0.0, "pi_total": round(pi, 4),
            "theta_w_per_bp": round(theta_w / L, 8) if L else 0.0,
            "segregating_sites": sepsites, "tajima_D": round(d, 4),
            "n_samples": n, "window": L, "mean_variant_spacing": round(
                stats.mean([b - a for a, b in zip(pos_list, pos_list[1:])]), 2)
            if len(pos_list) > 1 else 0.0}


def sliding_diversity(vcf: Vcf, window: int = 1000, step: int | None = None,
                      chrom: str | None = None) -> list[dict]:
    step = step or window
    rows = defaultdict(list)
    for r in vcf.records:
        rows[r.chrom].append(r)
    out = []
    for ch, recs in rows.items():
        if chrom and ch != chrom:
            continue
        recs.sort(key=lambda r: r.pos)
        maxpos = recs[-1].pos if recs else 0
        for start in range(1, maxpos, step):
            sub = Vcf(vcf.header, vcf.samples, [r for r in recs if start <= r.pos < start + window])
            if not sub.records:
                continue
            st = nucleotide_diversity(sub, window=window)
            out.append({"chrom": ch, "start": start, "end": start + window,
                        "n_variants": len(sub.records), **st})
    return out


def fst_hudson(p1: float, n1: int, p2: float, n2: int) -> dict:
    """Hudson's Fst estimator from allele frequencies and sample sizes."""
    if not (0 <= p1 <= 1 and 0 <= p2 <= 1) or n1 < 2 or n2 < 2:
        return {"numerator": 0.0, "denominator": 0.0, "fst": float("nan")}
    num = (p1 - p2) ** 2 - p1 * (1 - p1) / (n1 - 1) - p2 * (1 - p2) / (n2 - 1)
    den = p1 * (1 - p2) + p2 * (1 - p1)
    return {"numerator": num, "denominator": den,
            "fst": num / den if den else float("nan")}


def pairwise_fst(vcf: Vcf, groups: dict[str, list[str]] | list[list[str]]) -> list[dict]:
    samples = list(vcf.samples)
    if isinstance(groups, dict):
        named = list(groups.items())
        idx = {n: samples.index(n) for _, g in groups.items() for n in g if n in samples}
    else:
        named = [(f"group{i + 1}", g) for i, g in enumerate(groups)]
        idx = {n: samples.index(n) for g in groups for n in g if n in samples}
    gts = all_genotypes(vcf)
    out = []
    for i in range(len(named)):
        for j in range(i + 1, len(named)):
            a_names, b_names = named[i][1], named[j][1]
            a_idx = [idx[n] for n in a_names if n in idx]
            b_idx = [idx[n] for n in b_names if n in idx]
            num = den = 0.0
            n1, n2 = 2 * len(a_idx), 2 * len(b_idx)
            for k in range(len(vcf.records)):
                pa = _freq(gts, a_idx, k)
                pb = _freq(gts, b_idx, k)
                if pa is None or pb is None:
                    continue
                r = fst_hudson(pa, n1, pb, n2)
                if math.isnan(r["fst"]):
                    continue
                num += r["numerator"]
                den += r["denominator"]
            out.append({"population_1": named[i][0], "population_2": named[j][0],
                        "Fst": round(num / den, 6) if den else float("nan"),
                        "n_sites": len(vcf.records), "n1": n1, "n2": n2})
    return out


def _freq(gts, idx: list[int], k: int) -> float | None:
    tot = alt = 0
    for i in idx:
        a, b = gts[i][k]
        if a < 0:
            continue
        tot += 2
        alt += (a > 0) + (b > 0)
    return alt / tot if tot >= 2 else None


def weir_cockerham_fst(vcf: Vcf, pop_of_sample: dict[str, str]) -> dict:
    """Weir & Cockerham theta averaged over loci for named populations."""
    pops: dict[str, list[int]] = defaultdict(list)
    for i, sname in enumerate(vcf.samples):
        pops[pop_of_sample.get(sname, "all")].append(i)
    names = list(pops)
    if len(names) < 2:
        return {"theta": float("nan"), "populations": names, "n_loci": 0}
    gts = all_genotypes(vcf)
    a = b = c = 0.0
    n_loci = 0
    for k in range(len(vcf.records)):
        freqs, ns, hets = [], [], []
        for nm in names:
            idx = pops[nm]
            tot = alt = het = 0
            for i in idx:
                x, y = gts[i][k]
                if x < 0:
                    continue
                tot += 2
                alt += (x > 0) + (y > 0)
                het += (x != y)
            if tot:
                freqs.append(alt / tot)
                ns.append(tot)
                hets.append(het / len(idx))
        if len(freqs) < 2:
            continue
        n_loci += 1
        r_ = len(freqs)
        nbar = sum(ns) / r_
        nc = (r_ * nbar - sum(n * n for n in ns) / (r_ * nbar)) / max(1e-9, r_ - 1.0)
        pbar = sum(n * p for n, p in zip(ns, freqs)) / sum(ns)
        s2 = sum(n * (p - pbar) ** 2 for n, p in zip(ns, freqs)) / ((r_ - 1) * nbar)
        hbar = sum(hets) / r_
        a += nbar / nc * (s2 - 1 / (nbar - 1) * (pbar * (1 - pbar) - (r_ - 1) / r_ * s2
                                                   - 0.25 * hbar)) if nc else 0.0
        b += nbar / (nbar - 1) * (pbar * (1 - pbar) - (r_ - 1) / r_ * s2
                                 - (2 * nbar - 1) / (4 * nbar) * hbar) if nbar > 1 else 0.0
        c += 0.25 * hbar
    return {"theta": round(a / (a + b + c), 6) if (a + b + c) else float("nan"),
            "a": round(a, 6), "b": round(b, 6), "c": round(c, 6),
            "populations": names, "n_loci": n_loci}


def ld_r2(vcf: Vcf, sample_idx: int | None = None, max_pairs: int = 20000,
          window: int = 0) -> list[dict]:
    """Pairwise r^2 / D' from haplotype counts (phased or from GT dosage)."""
    haps = haplotypes(vcf)
    k = len(vcf.records)
    out = []
    for i in range(k):
        for j in range(i + 1, k):
            if window and abs(vcf.records[j].pos - vcf.records[i].pos) > window:
                continue
            if len(out) >= max_pairs:
                return out
            r = _r2(haps, i, j)
            out.append({"chrom": vcf.records[i].chrom, "pos1": vcf.records[i].pos,
                        "pos2": vcf.records[j].pos, "distance": vcf.records[j].pos - vcf.records[i].pos,
                        "r2": round(r["r2"], 5), "D": round(r["D"], 5),
                        "Dprime": round(r["Dprime"], 5)})
    return out


def haplotypes(vcf: Vcf) -> list[list[int]]:
    """Simple 0/1 haplotype matrix (phased if '|' present, else per-site dosage>0)."""
    n_hap = max(1, len(vcf.samples) * 2)
    mat = [[0] * len(vcf.records) for _ in range(n_hap)]
    for k, r in enumerate(vcf.records):
        h = 0
        for i in range(len(vcf.samples)):
            a, b = (genotypes(vcf, i)[k])
            if a >= 0:
                mat[h % n_hap][k] = 1 if a > 0 else 0
                h += 1
                mat[h % n_hap][k] = 1 if b > 0 else 0
                h += 1
    return mat


def _r2(haps: list[list[int]], i: int, j: int) -> dict:
    n = len(haps)
    if n < 2:
        return {"r2": 0.0, "D": 0.0, "Dprime": 0.0}
    p1 = sum(h[i] for h in haps) / n
    p2 = sum(h[j] for h in haps) / n
    p12 = sum(1 for h in haps if h[i] and h[j]) / n
    D = p12 - p1 * p2
    denom = math.sqrt(p1 * (1 - p1) * p2 * (1 - p2))
    r2 = (D * D / denom ** 2) if denom else 0.0
    dmax = min(p1 * (1 - p2), (1 - p1) * p2) if D > 0 else min(p1 * p2, (1 - p1) * (1 - p2))
    return {"r2": r2, "D": D, "Dprime": abs(D) / dmax if dmax else 0.0}


def king_kinship(vcf: Vcf) -> list[dict]:
    """KING-robust kinship / IBS0 for unrelated pairs (a PLINK --genome-lite)."""
    gts = [allele_counts(vcf, i) for i in range(len(vcf.samples))]
    out = []
    for i in range(len(gts)):
        for j in range(i + 1, len(gts)):
            a, b = gts[i], gts[j]
            n = min(len(a), len(b))
            het_a = sum(1 for x in a if x == 1)
            het_b = sum(1 for x in b if x == 1)
            N_Aa = sum(1 for x, y in zip(a[:n], b[:n]) if x == 1 and y in (0, 2))
            N_aA = sum(1 for x, y in zip(a[:n], b[:n]) if y == 1 and x in (0, 2))
            N_ib0 = sum(1 for x, y in zip(a[:n], b[:n]) if (x, y) in ((0, 2), (2, 0)))
            N_het_het = sum(1 for x, y in zip(a[:n], b[:n]) if x == 1 and y == 1)
            phi = (N_het_het * 2 - 2 * N_ib0) / (min(het_a, het_b) * 2) if min(het_a, het_b) else 0.0
            kinship = -0.5 if phi < -0.125 and N_ib0 > 0 else (phi / 2 if n else 0.0)
            out.append({"sample1": vcf.samples[i], "sample2": vcf.samples[j],
                        "kinship": round(max(-0.5, min(0.5, kinship)), 5),
                        "IBS0": N_ib0, "within_family_het_ratio": round(
                            N_het_het / het_a if het_a else 0.0, 4),
                        "n_sites": n,
                        "relationship": _infer_rel(max(-0.5, min(0.5, kinship)), N_ib0, n)})
    return out


def _infer_rel(k: float, ibs0: int, n: int) -> str:
    if n < 4:
        return "unknown"
    if k > 0.354:
        return "duplicate/MZ twin"
    if k > 0.177:
        return "1st degree"
    if k > 0.0884:
        return "2nd degree"
    if k > 0.0442:
        return "3rd degree"
    return "unrelated" if ibs0 > 0 else "unrelated/low coverage"


def concordance(vcf: Vcf) -> list[dict]:
    gts = [genotypes(vcf, i) for i in range(len(vcf.samples))]
    out = []
    for i in range(len(gts)):
        for j in range(i + 1, len(gts)):
            same = diff = miss = 0
            for a, b in zip(gts[i], gts[j]):
                if a[0] < 0 or b[0] < 0:
                    miss += 1
                elif a == b:
                    same += 1
                else:
                    diff += 1
            tot = same + diff
            out.append({"sample1": vcf.samples[i], "sample2": vcf.samples[j],
                        "concordance": round(same / tot, 5) if tot else float("nan"),
                        "discordant": diff, "concordant": same, "missing": miss})
    return out


def vaf_from_pileup(counts: list[tuple[int, int]], ref_base: str, alt_base: str) -> dict:
    ref_n = sum(c[0] for c in counts)
    alt_n = sum(c[1] for c in counts)
    tot = ref_n + alt_n
    return {"ref_count": ref_n, "alt_count": alt_n, "depth": tot,
            "vaf": round(alt_n / tot, 5) if tot else 0.0,
            "base": f"{ref_base}>{alt_base}"}


def genotype_likelihood(depth: int, alt: int, err: float = 0.01) -> dict:
    """Diploid genotype posteriors from read counts (uniform priors)."""
    def lik(g: int) -> float:
        p = 0.0
        for k in range(depth + 1):
            if k != alt:
                continue
            pa = 0.5 * g + err * (1 - 0.5 * g) if g < 2 else 1 - err
            pa = min(max(pa, 1e-9), 1 - 1e-9)
            p = math.comb(depth, k) * pa ** k * (1 - pa) ** (depth - k)
        return p

    l = [lik(g) for g in (0, 1, 2)]
    tot = sum(l) or 1.0
    post = [x / tot for x in l]
    gl = [10 * math.log10(max(1e-12, 1 - p)) for p in post]
    exp_alt = sum(g * p for g, p in zip((0, 1, 2), post)) / 2
    return {"PL": [round(x, 2) for x in gl], "GL": [round(math.log10(max(1e-12, x)), 4) for x in l],
            "posterior_0/0": round(post[0], 4), "posterior_0/1": round(post[1], 4),
            "posterior_1/1": round(post[2], 4), "GT": ["0/0", "0/1", "1/1"][int(round(
                sum(g * p for g, p in zip((0, 1, 2), post))))],
            "expected_alt_fraction": round(exp_alt, 4), "depth": depth, "alt": alt}


def left_align_and_norm(chrom_seq: str, pos: int, ref: str, alt: str) -> tuple[int, str, str]:
    """VCF normalisation: trim common suffix then shift left (bcftools norm)."""
    while len(ref) > 1 and len(alt) > 1 and ref[-1] == alt[-1]:
        ref, alt = ref[:-1], alt[:-1]
    while len(ref) != len(alt) and pos > 1:
        prev = chrom_seq[pos - 2] if pos - 2 < len(chrom_seq) else ""
        if not prev or (len(ref) > 1 and ref[0] == prev and alt[0] == ref[0]) or \
                (len(alt) > 1 and alt[0] == prev and (len(ref) == 1 or ref[0] == prev)):
            pass
        if len(alt) > len(ref) and alt[0] == prev and len(ref) >= 1:
            ref, alt = prev + ref, alt[1:]
            pos -= 1
        elif len(ref) > len(alt) and ref[0] == prev:
            ref, alt = ref[1:], prev + alt
            pos -= 1
        else:
            break
    return pos, ref, alt


def snp_eff_lite(vcf: Vcf, gff_rows, genome: dict[str, str] | None = None) -> list[dict]:
    """Annotated consequence by position relative to gene/exon/CDS features."""
    from chroma_titan.core.io import GffRow

    genes = [r for r in gff_rows if r.type.lower() in ("gene", "mrna", "transcript")]
    cds = [r for r in gff_rows if r.type.lower() in ("cds",)]
    rows = []
    for rec in vcf.records:
        ann = {"chrom": rec.chrom, "pos": rec.pos, "ref": rec.ref, "alt": rec.alt,
               "qual": rec.qual, "filter": rec.filter, "gene": ".", "feature": "intergenic",
               "consequence": "intergenic_variant", "codon_change": ".",
               "aa_change": ".", "distance_to_gene": ""}
        best = None
        for g in genes:
            if g.seqid != rec.chrom:
                continue
            if g.start - 2000 <= rec.pos <= g.end + 2000:
                inside = g.start <= rec.pos <= g.end
                dist = 0 if inside else min(abs(rec.pos - g.start), abs(rec.pos - g.end))
                if best is None or dist < best[0]:
                    best = (dist, g)
        if best:
            dist, g = best
            ann["gene"] = g.get("Name") or g.get("gene_id") or g.get("ID") or "."
            ann["distance_to_gene"] = dist
            coding = [c for c in cds if c.seqid == g.seqid and c.get("Parent") == g.get("ID")] or \
                     [c for c in cds if c.seqid == g.seqid and g.start - 1 <= c.start - 1 <= g.end]
            if coding and g.start <= rec.pos <= g.end:
                ann["feature"] = "exon"
                in_cds = next((c for c in coding if c.start <= rec.pos <= c.end), None)
                if rec.is_snv and in_cds and genome:
                    ann["consequence"] = _codon_effect(genome.get(g.seqid, ""), in_cds, g,
                                                       coding, rec)
                    if ann["consequence"] != "non_coding_exon_variant":
                        ann["feature"] = "CDS"
                else:
                    ann["consequence"] = "non_coding_exon_variant"
            elif dist == 0:
                ann["feature"] = "intron"
                ann["consequence"] = "intron_variant"
            elif dist <= 500:
                ann["feature"] = "upstream" if g.strand == "+" else "downstream"
                ann["consequence"] = f"{ann['feature']}_variant"
            else:
                ann["consequence"] = "upstream_gene_variant"
        rows.append(ann)
    return rows


def _codon_effect(genome_seq: str, cds: GffRow, gene: GffRow, coding: list, rec) -> str:
    from chroma_titan.core.seq import CODON_TABLE

    if not genome_seq:
        return "missense_variant"
    offset = rec.pos - cds.start
    codon_i = offset // 3
    if (rec.pos - cds.start) % 3:
        pos_in_codon = (rec.pos - cds.start) % 3
    else:
        pos_in_codon = 0
    start = cds.start - 1 + codon_i * 3
    codon = genome_seq[start:start + 3]
    if len(codon) < 3:
        return "frameshift_variant"
    new = list(codon)
    new[pos_in_codon] = rec.alt[0]
    new_codon = "".join(new)
    ref_aa, alt_aa = CODON_TABLE.get(codon, "X"), CODON_TABLE.get(new_codon, "X")
    cds_start = cds.start - 1
    del cds_start, gene, coding
    if ref_aa == alt_aa:
        return "synonymous_variant" if rec.is_snv else "inframe_deletion"
    if alt_aa == "*" or ref_aa == "*":
        return "stop_gained" if alt_aa == "*" else "stop_lost"
    if codon_i == 0 and (alt_aa == "M" or ref_aa == "M"):
        return "start_lost"
    return "missense_variant"


def variant_density(vcf: Vcf, window: int = 1000) -> list[dict]:
    by = defaultdict(Counter)
    for r in vcf.records:
        by[r.chrom][r.pos // window] += 1
    out = []
    for ch, cnt in by.items():
        for b in range(min(cnt), max(cnt) + 1):
            out.append({"chrom": ch, "start": b * window, "end": (b + 1) * window,
                        "n_variants": cnt.get(b, 0)})
    return out


def manhattan_from_vcf(vcf: Vcf, p_from: str = "QUAL") -> list[dict]:
    rows = []
    for r in vcf.records:
        p = 10 ** (-r.qual / 10) if p_from == "QUAL" and r.qual == r.qual else r.info_num("PV", float("nan"))
        rows.append({"chrom": r.chrom, "pos": r.pos, "p_value": min(1.0, max(1e-300, p)),
                     "variant": f"{r.ref}>{r.alt}", "qual": r.qual})
    return rows


def filter_records(vcf: Vcf, expression: str) -> list[VcfRecord]:
    """Filter with a small expression language over INFO/format fields."""
    ok = []
    for r in vcf.records:
        ctx = {"qual": r.qual if r.qual == r.qual else 0.0, "QUAL": ctx_qual(r),
               "DP": r.info_num("DP", 0), "MQ": r.info_num("MQ", 0), "AF": r.info_num("AF", 0),
               "AC": r.info_num("AC", 0), "filter": r.filter, "ref": r.ref, "alt": r.alt,
               "chrom": r.chrom, "pos": r.pos, "is_snv": r.is_snv, "is_indel": r.is_indel,
               "is_ts": r.is_transition, "n_alt": r.n_alt, "samples": len(r.samples)}
        for i, s in enumerate(r.samples):
            for k, v in s.items():
                if _isnum(v):
                    ctx[f"{k}_{i}"] = float(v)
        try:
            if eval(expression, {"__builtins__": {}}, ctx):  # noqa: S307
                ok.append(r)
        except Exception:  # noqa: BLE001
            continue
    return ok


def ctx_qual(r: VcfRecord) -> float:
    return r.qual if r.qual == r.qual else 0.0


def _isnum(x) -> bool:
    try:
        float(x)
        return True
    except (TypeError, ValueError):
        return False


def write_records(header: list[str], samples: list[str], records: list[VcfRecord],
                  fmt_keys: list[str] | None = None) -> str:
    out = list(header) or ["##fileformat=VCFv4.2"]
    keys = fmt_keys or (records[0].fmt_keys if records and records[0].fmt_keys else ["GT"])
    if samples:
        out.append("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(samples))
    for r in records:
        base = [r.chrom, str(r.pos), r.id, r.ref, r.alt,
                (f"{r.qual:g}" if r.qual == r.qual else "."), r.filter,
                ";".join(f"{k}={v}" if v != "" else k for k, v in r.info.items()) or "."]
        if samples:
            base.append(":".join(keys))
            for s in r.samples:
                base.append(":".join(s.get(k, ".") for k in keys))
        out.append("\t".join(base))
    return "\n".join(out) + "\n"


def vcf_to_tsv(vcf: Vcf, info_keys: list[str] | None = None,
              fmt_keys: list[str] | None = None) -> list[dict]:
    rows = []
    for r in vcf.records:
        rec = {"CHROM": r.chrom, "POS": r.pos, "ID": r.id, "REF": r.ref, "ALT": r.alt,
               "QUAL": r.qual, "FILTER": r.filter}
        for k in info_keys or sorted(r.info):
            v = r.info.get(k, ".")
            rec[f"INFO_{k}"] = float(v) if v and _isnum(v) else v
        for i, s in enumerate(r.samples):
            for k in (fmt_keys or r.fmt_keys or []):
                v = s.get(k, ".")
                rec[f"{k}:{i}"] = float(v) if v and _isnum(v) else v
        rows.append(rec)
    return rows


def vcf_summary_table(vcf: Vcf) -> dict:
    per_chrom = Counter(r.chrom for r in vcf.records)
    return {"n_variants": len(vcf.records), "n_samples": len(vcf.samples),
            "chromosomes": dict(per_chrom),
            "ti_tv": ts_tv_summary(vcf), "indel": indel_spectrum(vcf),
            "novel_vs_known": sum(1 for r in vcf.records if r.id in (".", "")),
            "quality_mean": round(stats.mean([r.qual for r in vcf.records if r.qual == r.qual]), 3),
            "filter_counts": dict(Counter(r.filter for r in vcf.records))}


def relatedness_matrix(vcf: Vcf) -> list[list[float]]:
    """GRM-style IBS relatedness (mean of genotype concordance - 0.5)."""
    gts = [allele_counts(vcf, i) for i in range(len(vcf.samples))]
    n = len(gts)
    out = [[0.0] * n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            a, b = gts[i], gts[j]
            m = min(len(a), len(b))
            if not m:
                continue
            out[i][j] = round(sum(1 - abs(a[k] - b[k]) / 2 for k in range(m)) / m - 0.5, 4)
    return out


def ploidy_from_vcf(vcf: Vcf) -> dict:
    ratios = []
    for r in vcf.records:
        for i in range(len(r.samples)):
            v = r.field(i, "AD")
            if v and "/" in v:
                try:
                    a, b = (float(x) for x in v.split("/")[:2])
                    if a + b:
                        ratios.append(b / (a + b))
                except ValueError:
                    continue
    return {"n_observed": len(ratios), "mean_alt_fraction": round(stats.mean(ratios), 4),
            "heterozygosity_peaks": sorted({round(x, 2) for x in ratios if x})[:10],
            "inferred_ploidy": 2 if 0.35 < stats.mean(ratios) < 0.65 else
            (1 if stats.mean(ratios) > 0.85 else "variable")}


def merge_vcfs(vcfs: list[Vcf], strategy: str = "union") -> Vcf:
    seen = {}
    for v in vcfs:
        for r in v.records:
            key = (r.chrom, r.pos, r.ref, r.alt)
            if key not in seen or (strategy == "first" and seen[key].qual < r.qual):
                seen[key] = r
            elif strategy == "best" and r.qual > seen[key].qual:
                seen[key] = r
    recs = sorted(seen.values(), key=lambda r: (r.chrom, r.pos))
    samples = list(vcfs[0].samples) if vcfs else []
    return Vcf(["##fileformat=VCFv4.2", "##merged_by=chroma-titan"], samples, recs)


def subsample(vcf: Vcf, n: int = 10, seed: int = 1, frac: float = 0.0) -> Vcf:
    import random

    rng = random.Random(seed)
    recs = vcf.records
    if frac:
        n = max(1, int(len(recs) * frac))
    if n >= len(recs):
        return Vcf(vcf.header, vcf.samples, list(recs))
    return Vcf(vcf.header, vcf.samples, sorted(rng.sample(recs, n), key=lambda r: (r.chrom, r.pos)))
