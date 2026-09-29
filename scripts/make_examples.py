"""Generate the bundled example datasets under ``data/examples/``.

Everything is deterministic (fixed seed) and self-consistent: the genome, its
annotation, the variants called against it, the simulated reads and their
alignments all describe the same tiny "organism", so cross-format tools
(pileup calling, coverage over BED, variant annotation...) produce meaningful
numbers in the demos and in the test-suite.

Run with::

    python scripts/make_examples.py
"""

from __future__ import annotations

import math
import random
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "data" / "examples"

ADAPTER = "AGATCGGAAGAGC"
CHROM1, CHROM2 = "chrV", "chrM"
SIZE1, SIZE2 = 6000, 3200


def rng(seed: int = 20260929) -> random.Random:
    return random.Random(seed)


def make_genome(r: random.Random) -> dict[str, str]:
    """Random genome, then overwrite gene bodies with coding-like sequence."""
    g = {CHROM1: "".join(r.choice("ACGT") for _ in range(SIZE1)),
          CHROM2: "".join(r.choice("ACGT") for _ in range(SIZE2))}
    return g


GENES = [
    # chrom, start(1-based), n_codons, strand, name, biotype
    (CHROM1, 401, 180, "+", "geneA", "protein_coding"),
    (CHROM1, 1501, 240, "-", "geneB", "protein_coding"),
    (CHROM1, 2601, 120, "+", "geneC", "protein_coding"),
    (CHROM1, 3901, 300, "+", "geneD", "lncRNA"),
    (CHROM1, 4601, 150, "-", "geneE", "protein_coding"),
    (CHROM2, 201, 210, "+", "viralRep", "protein_coding"),
    (CHROM2, 1201, 90, "-", "viralCap", "protein_coding"),
    (CHROM2, 2001, 140, "+", "viralEnv", "protein_coding"),
]
CODONS = ["GCT", "GCC", "GCA", "GCG", "TGT", "TGC", "GAT", "GAC", "AAA", "AAG", "TTC",
          "TTT", "GGT", "GGC", "CAT", "CAC", "ATC", "ATT", "ATG", "TTA", "CTG", "AAT",
          "CAA", "CAG", "TAT", "TAC", "TCC", "AGC", "ACT", "ACC", "GTT", "GTC", "TGG",
          "TAA", "TAG", "TGA"]
SENSE = [c for c in CODONS if c not in ("TAA", "TAG", "TGA")]
STOP = ("TAA", "TAG", "TGA")


def coding_seq(r: random.Random, n_codons: int, gc: float = 0.5) -> str:
    out = ["ATG"]
    for _ in range(n_codons - 2):
        if r.random() < gc:
            out.append(r.choice([c for c in SENSE if (c.count("G") + c.count("C")) >= 2]))
        else:
            out.append(r.choice(SENSE))
    out.append(r.choice(STOP))
    return "".join(out)


def revcomp(s: str) -> str:
    return s.translate(str.maketrans("ACGTN", "TGCAN"))[::-1]


def place_genes(g: dict[str, str], r: random.Random) -> tuple[dict, list]:
    from chroma_titan.core.seq import gc_content

    annot = []
    for i, (chrom, start, n_codons, strand, name, biotype) in enumerate(GENES, 1):
        cds_len = n_codons * 3
        if biotype == "lncRNA":
            body = "".join(r.choice("ACGT") for _ in range(cds_len))
        elif gc_content(g[chrom][start - 1:start - 1 + cds_len]) < 0:
            body = coding_seq(r, n_codons)
        else:
            body = coding_seq(r, n_codons, gc=0.35 + 0.3 * (i % 3) / 2)
        if strand == "-":
            body = revcomp(body)
        seq = list(g[chrom])
        seq[start - 1:start - 1 + cds_len] = list(body)
        g[chrom] = "".join(seq)
        end = start + cds_len - 1
        # two exons for coding genes (5'UTR-ish split), single exon otherwise
        if biotype == "protein_coding" and cds_len > 300:
            cut = start + cds_len // 3
            exons = [(start, cut - 1), (cut, end)]
        else:
            exons = [(start, end)]
        annot.append({"i": i, "chrom": chrom, "start": start, "end": end, "strand": strand,
                      "name": name, "biotype": biotype, "exons": exons, "cds": cds_len,
                      "gene_id": f"GENE{i:03d}", "transcript": f"TX{i:03d}.1"})
    return g, annot


def write_fasta(path: Path, g: dict[str, str], width: int = 60) -> None:
    lines = []
    for name, seq in g.items():
        lines.append(f">{name} synthetic example genome chromosome")
        for i in range(0, len(seq), width):
            lines.append(seq[i:i + width])
    path.write_text("\n".join(lines) + "\n")


def write_gff(path: Path, annot: list) -> None:
    rows = ["##gff-version 3",
            "##sequence-region chrV 1 6000",
            "##sequence-region chrM 1 3200"]
    for a in annot:
        rows.append(f"{a['chrom']}\tchroma\tgene\t{a['start']}\t{a['end']}\t.\t{a['strand']}\t.\t"
                    f"ID={a['gene_id']};Name={a['name']};gene_biotype={a['biotype']}")
        rows.append(f"{a['chrom']}\tchroma\ttranscript\t{a['start']}\t{a['end']}\t.\t{a['strand']}\t.\t"
                    f"ID={a['transcript']};Parent={a['gene_id']};Name={a['name']}")
        for j, (s, e) in enumerate(a["exons"], 1):
            rows.append(f"{a['chrom']}\tchroma\texon\t{s}\t{e}\t.\t{a['strand']}\t.\t"
                        f"ID=exon_{a['transcript']}_{j};Parent={a['transcript']}")
        if a["biotype"] == "protein_coding":
            rows.append(f"{a['chrom']}\tchroma\tCDS\t{a['start']}\t{a['end']}\t.\t{a['strand']}\t0\t"
                        f"ID=cds_{a['transcript']};Parent={a['transcript']};gene={a['name']}")
        rows.append(f"{a['chrom']}\tchroma\tfive_prime_UTR\t{a['start'] - 60}\t{a['start'] - 1}\t."
                    f"\t{a['strand']}\t.\tID=utr5_{a['transcript']};Parent={a['transcript']}")
    path.write_text("\n".join(rows) + "\n")


def make_variants(g: dict[str, str], annot: list, r: random.Random) -> list:
    """Spiked variants: mostly inside genes so annotation tools find effects."""
    variants = []
    n = 0
    for a in annot:
        span = a["end"] - a["start"]
        for _ in range(6):
            pos = a["start"] + r.randrange(max(1, span - 1))
            ref = g[a["chrom"]][pos - 1]
            n += 1
            kind = r.random()
            if kind < 0.62:
                alt = r.choice([b for b in "ACGT" if b != ref])
            elif kind < 0.85:
                alt = ref + "".join(r.choice("ACGT") for _ in range(r.choice([1, 2, 4])))
            else:
                alt = ref[: max(1, len(ref))]
                del_len = r.choice([1, 2, 3])
                alt = ref + "" if del_len == 0 else ref + g[a["chrom"]][pos:pos + del_len]
                alt = ""
                variants.append((a["chrom"], pos, f"chrV_{n}", ref, "", "DEL", del_len))
                continue
            variants.append((a["chrom"], pos, f"var_{n}", ref, alt, "SNP" if len(ref) == len(alt) == 1 else "INDEL", 0))
    for _ in range(26):
        chrom = r.choice([CHROM1, CHROM2])
        pos = r.randrange(1, len(g[chrom]))
        ref = g[chrom][pos - 1]
        alt = r.choice([b for b in "ACGT" if b != ref])
        n += 1
        variants.append((chrom, pos, f"var_{n}", ref, alt, "SNP", 0))
    variants.sort(key=lambda v: (v[0], v[1]))
    return variants


def write_reads(g: dict[str, str], variants: list, r: random.Random):
    """Simulated 2x150bp paired-end reads with errors, adapters and duplicates."""
    mut = {}
    for chrom, pos, _id, ref, alt, kind, dlen in variants:
        if kind == "SNP":
            mut[(chrom, pos)] = alt
    seqs = {c: list(g[c]) for c in g}
    for (c, p), alt in mut.items():
        seqs[c][p - 1] = alt
    seqs = {c: "".join(v) for c, v in seqs.items()}
    frag_len, read_len = 320, 150
    n_pairs = 260
    recs = []
    for i in range(n_pairs):
        chrom = r.choice([CHROM1, CHROM1, CHROM2])
        s = r.randrange(0, max(1, len(seqs[chrom]) - frag_len))
        frag = seqs[chrom][s:s + frag_len]
        fwd = frag[:read_len]
        rev = revcomp(frag[-read_len:])
        def err(seq, rate=0.006):
            out = list(seq)
            for j in range(len(out)):
                if r.random() < rate:
                    out[j] = r.choice([b for b in "ACGT" if b != out[j]])
            return "".join(out)
        def qual(seq, base=36, drop=0.0):
            q = []
            for j in range(len(seq)):
                v = base - drop * j / max(1, len(seq)) - (r.random() * 5 if r.random() < 0.25 else 0)
                q.append(chr(33 + max(2, min(41, int(round(v))))))
            return "".join(q)
        recs.append((f"read_{i + 1}/1", err(fwd), qual(fwd, 38, 8)))
        recs.append((f"read_{i + 1}/2", err(rev), qual(rev, 36, 12)))
        if i < 14:  # adapter read-through
            nm = recs[-1][1]
            recs[-1] = (recs[-1][0], nm[:110] + ADAPTER, recs[-1][2])
        if i in (7, 33, 61):  # poly-G / low quality tail artefacts
            recs[-1] = (recs[-1][0], recs[-1][1][:70] + "G" * 30, recs[-1][2][:70] + "5" * 30)
    dups = recs[:18]
    recs.extend([(f"{name}", sq, q) for name, sq, q in dups])
    n = int(0.05 * len(recs))
    for _ in range(n):
        recs.append((f"junk_{r.randrange(9999)}/{r.choice([1, 2])}",
                     "".join(r.choice("ACGT") for _ in range(read_len)),
                     "".join(r.choice("BCEHI") for _ in range(read_len))))
    for _ in range(6):
        recs.append((f"shorty_{r.randrange(99)}/{r.choice([1, 2])}",
                     "".join(r.choice("ACGT") for _ in range(r.randrange(8, 25))), "I" * 20))
    r1 = [x for x in recs if x[0].endswith("/1")]
    r2 = [x for x in recs if x[0].endswith("/2")]
    return seqs, r1, r2, mut


def write_fastq(path: Path, recs: list) -> None:
    lines = []
    for name, seq, qual in recs:
        lines.append(f"@{name} synthetic_simulated")
        lines.append(seq)
        lines.append("+")
        lines.append(qual)
    path.write_text("\n".join(lines) + "\n")


def write_sam(refs: dict[str, str], r1: list, r2: list, r: random.Random, mut) -> str:
    """Map the simulated reads back with a trivial exact seed search (self-consistent)."""
    lines = ["@HD\tVN:1.6\tSO:coordinate", "@PG\tID:chroma-sim\tPN:chroma-sim\tVN:1.0"]
    for name, seq in refs.items():
        lines.append(f"@SQ\tSN:{name}\tLN:{len(seq)}")
    lines.append("@RG\tID:flowcell1\tSM:sample_A\tPL:ILLUMINA\tLB:lib1")
    by_name = {}
    for name, seq, qual in list(r1) + list(r2):
        by_name.setdefault(name.rsplit("/", 1)[0], {})[name[-1]] = (seq, qual)
    alns = []
    idx = {}
    K = 24
    for chrom, seq in refs.items():
        d = {}
        for i in range(len(seq) - K + 1):
            d.setdefault(seq[i:i + K], []).append(i)
        idx[chrom] = d
    readno = 0
    for pair, mates in by_name.items():
        for which in ("1", "2"):
            if which not in mates:
                continue
            readno += 1
            seq, qual = mates[which]
            hits = []
            for off in (0, len(seq) - K):
                km = seq[off:off + K]
                for chrom, d in idx.items():
                    for pos in d.get(km, ()):
                        cand = pos - off
                        if 0 <= cand <= len(refs[chrom]) - len(seq):
                            hits.append((chrom, cand))
            if not hits:
                alns.append((pair + "/" + which, 4 | (64 if which == "1" else 128) | 1, "*", 0, 0,
                             f"{len(seq)}S", "*", 0, 0, seq, qual, None))
                continue
            chrom, pos = hits[0]
            ref_sub = refs[chrom][pos:pos + len(seq)]
            mm = [i for i in range(min(len(seq), len(ref_sub))) if seq[i] != ref_sub[i]]
            nm = len(mm)
            flag = (99 if which == "1" else 163) | (1024 if readno <= 36 else 0)
            clip = ""
            cigar = f"{len(seq)}M"
            if len(seq) < 100:
                clip = f"{len(seq) - 20}M20S" if which == "1" else f"20S{len(seq) - 20}M"
                cigar = clip
            alns.append((pair + "/" + which, flag, chrom, pos + 1, 60 if len(hits) == 1 else 3,
                         cigar, "=", pos + 1 + 320, 320, seq, qual, nm))
    alns.sort(key=lambda a: (a[2], a[3]))
    for a in alns:
        (qname, flag, rname, pos, mapq, cigar, rnext, pnext, tlen, seq, qual, nm) = a
        extra = f"\tNM:i:{nm}" if nm is not None else ""
        extra += "\tRG:Z:flowcell1"
        lines.append("\t".join([qname, str(flag), rname, str(pos), str(mapq), cigar, rnext,
                                str(pnext), str(tlen), seq, qual]) + extra)
    return "\n".join(lines) + "\n"


def write_vcf(path: Path, variants: list, samples: list[str], r: random.Random) -> None:
    lines = ["##fileformat=VCFv4.2", "##source=chroma-titan-simulator",
             "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total depth\">",
             "##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele freq\">",
             '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
             '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth">',
             '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths">',
             '##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality">',
             "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(samples)]
    for chrom, pos, vid, ref, alt, kind, dlen in variants:
        if kind == "DEL":
            ref = ref + (alt or "A") * 0
            ref = "A"
            alt = ""
        if kind == "DEL":
            continue
        qual = round(r.uniform(20, 250), 2)
        filt = "PASS" if qual > 60 else "LowQual"
        dp = r.randrange(8, 60)
        af = round(r.uniform(0.05, 0.6), 3)
        cells = []
        for _ in samples:
            g = r.choices(["0/0", "0/1", "1/1", "./."], weights=[5, 4, 2, 1])[0]
            d = r.randrange(4, 40)
            if g == "0/0":
                ad = f"{d},0"
            elif g == "1/1":
                ad = f"0,{d}"
            elif g == "./.":
                ad = ".,."
            else:
                a = r.randrange(1, d)
                ad = f"{d - a},{a}"
            cells.append(f"{g}:{d}:{ad}:{r.randrange(10, 99)}")
        lines.append("\t".join([chrom, str(pos), vid, ref, alt or ".", str(qual), filt,
                               f"DP={dp};AF={af}", "GT:DP:AD:GQ"] + cells))
    path.write_text("\n".join(lines) + "\n")


def write_bed(path: Path, annot: list, r: random.Random) -> None:
    rows = []
    for a in annot:
        rows.append(f"{a['chrom']}\t{a['start'] - 1}\t{a['end']}\t{a['name']}\t"
                    f"{r.randrange(100, 900)}\t{a['strand']}")
    for i in range(12):
        chrom = r.choice([CHROM1, CHROM2])
        s = r.randrange(0, 3000)
        rows.append(f"{chrom}\t{s}\t{s + r.randrange(80, 600)}\tpeak_{i + 1}\t{r.randrange(50, 999)}\t.")
    path.write_text("\n".join(sorted(rows)) + "\n")


def write_targets(path: Path, r: random.Random) -> None:
    rows = []
    for i in range(14):
        chrom = r.choice([CHROM1, CHROM2])
        s = r.randrange(0, 2800)
        rows.append(f"{chrom}\t{s}\t{s + r.randrange(150, 900)}\tsite_{i + 1}")
    path.write_text("\n".join(sorted(set(rows))) + "\n")


def write_counts(path: Path, annot: list, samples: list[str], r: random.Random) -> None:
    groups = {"control": samples[:3], "treated": samples[3:]}
    rows = []
    for a in annot:
        base = r.randrange(60, 3000)
        vals = []
        fold = 1.0 if a["biotype"] != "protein_coding" else r.choice([0.4, 0.55, 2.4, 3.1, 1.0])
        for s in samples:
            v = base * (fold if s in groups["treated"] else 1.0)
            vals.append(max(0, int(round(v * r.uniform(0.75, 1.25)))))
        rows.append([a["name"]] + vals)
    extra = [["geneF", *[r.randrange(0, 30) for _ in samples]],
             ["geneG", *[r.randrange(200, 900) for _ in samples]],
             ["spike_in", *[500 for _ in samples]],
             ["geneH_zero", *[0 for _ in samples]]]
    rows.extend(extra)
    path.write_text("\t".join(["gene"] + samples) + "\n" +
                    "\n".join("\t".join([str(x) for x in r_]) for r_ in rows) + "\n")


def write_gmt(path: Path, annot: list) -> None:
    names = [a["name"] for a in annot]
    sets = {"REPLICATION": [n for n in names if "Rep" in n or n in ("geneA", "geneC")] + names[:2],
            "CAPSID_ASSEMBLY": [n for n in names if "Cap" in n or "Env" in n],
            "LNC_REGIONS": [n for n in names if "D" in n or "E" in n],
            "HIGH_EXPRESSED": names[:4]}
    path.write_text("".join(f"{k}\tpathway_{i}\t" + "\t".join(v) + "\n"
                           for i, (k, v) in enumerate(sets.items(), 1)))


def write_motif(path: Path) -> None:
    path.write_text(""">MA0001.1 ExampleFactor
A [ 8 0 0 0 12 0 2 ]
C [ 1 0 11 0 0 12 0 ]
G [ 1 12 0 11 0 0 0 ]
T [ 2 0 1 1 0 0 10 ]
""")


def write_motifs_jaspar(path: Path) -> None:
    path.write_text("""MYB
A  [   8    0    0    0   12    0    2 ]
C  [   1    0   11    0    0   12    0 ]
G  [   1   12    0   11    0    0    0 ]
T  [   2    0    1    1    0    0   10 ]

MYC
A  [   0    1    0    2   10    1 ]
C  [   0    0    0    0    0    2 ]
G  [  11    0   11    8    0    0 ]
T  [   1   11    1    2    2    1 ]
""")


def write_proteins(path: Path, g: dict[str, str], annot: list) -> None:
    from chroma_titan.core.seq import translate

    lines = []
    for a in annot:
        if a["biotype"] != "protein_coding":
            continue
        seq = g[a["chrom"]][a["start"] - 1:a["end"]]
        if a["strand"] == "-":
            seq = revcomp(seq)
        prot = translate(seq)
        lines.append(f">{a['name']} protein_coding length={len(prot)}")
        for i in range(0, len(prot), 60):
            lines.append(prot[i:i + 60])
    lines.append(">kinase_like_control SH3-like domain sample")
    lines.append("MAAAKLIFVKEGKTVYCRQVLEAGDELPLSLAAGVIVSEKRPTKQFRG")
    path.write_text("\n".join(lines) + "\n")


def write_msa(path: Path, g: dict[str, str], annot: list, r: random.Random) -> None:
    """A gappy alignment derived from the coding sequences of the paralogue genes."""
    from chroma_titan.core.seq import translate

    seqs = []
    for a in annot[:5]:
        if a["biotype"] != "protein_coding":
            continue
        s = g[a["chrom"]][a["start"] - 1:a["end"]]
        if a["strand"] == "-":
            s = revcomp(s)
        p = translate(s)
        seqs.append(p[:60])
    L = max(len(s) for s in seqs)
    out = []
    for i, s in enumerate(seqs):
        row = list(s + "-" * (L - len(s)))
        for _ in range(int(0.06 * L)):
            row[r.randrange(L)] = r.choice("ACDEFGHIKLMNPQRSTVWY")
        out.append(f">seq_{i + 1}\n{''.join(row)}")
    path.write_text("\n".join(out) + "\n")


def write_tree(path: Path) -> None:
    path.write_text("((seq_1:0.12,seq_2:0.08)AB:0.21,(seq_3:0.17,seq_4:0.15)CD:0.09)root;\n")


def write_pairs(path: Path, r: random.Random) -> None:
    lines = ["#columns=	fragment	name	chrom1	pos1	chrom2	pos2	strand1	strand2"]
    for i in range(320):
        p1 = r.randrange(0, SIZE1 - 1) * 4
        same = r.random() > 0.25
        p2 = max(0, p1 + int(r.gauss(0, 12000))) if same else r.randrange(0, SIZE1 - 1) * 4
        lines.append(f"pair{i}\tMboI\t{CHROM1}\t{p1}\t{CHROM1 if same else CHROM2}\t"
                     f"{min(p2, SIZE1 - 2 if same else SIZE2 - 2)}\t+\t-")
    path.write_text("\n".join(lines) + "\n")


def write_bedgraph(path: Path, annot: list, r: random.Random) -> None:
    rows = []
    for a in annot:
        mid = (a["start"] + a["end"]) // 2
        for w in range(-600, 601, 100):
            s = mid + w
            d = 40 * math.exp(-(w / 260.0) ** 2) + r.uniform(0, 3)
            rows.append((a["chrom"], max(0, s), max(1, s + 100), round(d, 3)))
    for chrom, size in ((CHROM1, SIZE1), (CHROM2, SIZE2)):
        for s in range(0, size, 500):
            rows.append((chrom, s, s + 500, round(r.uniform(0, 6), 3)))
    rows.sort(key=lambda x: (x[0], x[1]))
    path.write_text("\n".join(f"{c}\t{s}\t{e}\t{v}" for c, s, e, v in rows) + "\n")


def write_pgm(path: Path, r: random.Random) -> None:
    """A 48x48 synthetic fluorescence micrograph: 4 nuclei + background."""
    h = w = 48
    img = [[max(0.0, r.gauss(30, 8)) for _ in range(w)] for _ in range(h)]
    blobs = [(10, 12, 6, 200), (14, 30, 5, 160), (30, 20, 7, 230), (36, 38, 4, 120)]
    for cy, cx, rad, peak in blobs:
        for y in range(h):
            for x in range(w):
                d2 = (y - cy) ** 2 + (x - cx) ** 2
                if d2 <= rad * rad:
                    img[y][x] = peak * math.exp(-d2 / (rad * rad * 0.7)) + r.gauss(0, 6)
    flat = [min(255, max(0, int(round(v)))) for row in img for v in row]
    text = f"P2\n{w} {h}\n255\n"
    for y in range(h):
        text += " ".join(str(v) for v in flat[y * w:(y + 1) * w]) + "\n"
    path.write_text(text)


def write_peptides(path: Path, r: random.Random) -> None:
    from chroma_titan.core.protein import peptide_mz

    rows = ["peptide\tmodified_sequence\tcharge\traw_intensity\truntime\tprotein"]
    proteins = ["geneA", "geneB", "geneC", "geneE", "viralRep"]
    for i in range(40):
        pep = "".join(r.choice("ACDEFGHIKLMNPQRSTVWY") for _ in range(r.randrange(7, 15)))
        mz = round(peptide_mz(pep, 2), 4)
        rows.append(f"{pep}\t{pep}\t2\t{r.randrange(5000, 990000)}\t{round(r.uniform(8, 62), 3)}\t{r.choice(proteins)}")
    path.write_text("\n".join(rows) + "\n")


def write_pheno(path: Path, r: random.Random, samples: list[str]) -> None:
    rows = ["sample\tgroup\tyield\tresistance_score\tinfected\tbiomass\ttreatment"]
    for s in samples:
        grp = "control" if samples.index(s) < 3 else "treated"
        base = 10 if grp == "treated" else 6
        rows.append(f"{s}\t{grp}\t{round(r.gauss(base, 1.1), 3)}\t{r.randrange(0, 100)}\t"
                    f"{r.choice([0, 1])}\t{round(abs(r.gauss(2.2, 0.4)), 3)}\t"
                    f"{'A' if grp == 'control' else 'B'}")
    path.write_text("\n".join(rows) + "\n")


def write_ontology(path: Path) -> None:
    path.write_text("""##obo-format=owl/functional
ontology: chroma-phenotype-example

[Term]
id: GO:0006260
name: DNA replication
is_a: GO:0006259 ! DNA metabolic process

[Term]
id: GO:0006259
name: DNA metabolic process
is_a: GO:0044238 ! cellular metabolic process

[Term]
id: GO:0019039
name: viral capsid assembly
is_a: GO:0019031 ! virion maturation

[Term]
id: GO:0019031
name: virion maturation

[Term]
id: GO:0003677
name: DNA binding
is_a: GO:0003674 ! molecular_function
""")


def write_str(path: Path, g: dict[str, str], r: random.Random) -> None:
    """STR profile table (locus x sample allele sizes)."""
    loci = ["VSTR1", "VSTR2", "MSTR1", "MSTR2", "VSTR3"]
    samples = ["sample_A", "sample_B", "sample_C", "sample_D"]
    rows = ["Locus\tMotif\tChr\tStart\t" + "\t".join(samples)]
    for i, loc in enumerate(loci):
        chrom = CHROM1 if loc.startswith("V") else CHROM2
        motif = r.choice(["CA", "GATA", "CAG", "TCTT"])
        start = r.randrange(100, 2000)
        alleles = []
        for _ in samples:
            n = r.randrange(6, 22)
            alleles.append(f"{start + n * len(motif)},{start + (n + r.choice([0, 1, 2])) * len(motif)}")
        rows.append(f"{loc}\t{motif}\t{chrom}\t{start}\t" + "\t".join(alleles))
    path.write_text("\n".join(rows) + "\n")
    del g


def write_methylation(path: Path, g: dict[str, str], r: random.Random) -> None:
    """Bismark-style CX_report lines."""
    lines = ["# CHROM\tPOS\tSTRAND\tCONTEXT\tCYTOSINE_FRACTION\tCHH_COUNT\t"
             "H_C_COUNT\tM_C_COUNT\tX_COUNT"]
    seq = g[CHROM1]
    for i in range(1, 3000):
        if seq[i - 1] == "C" and i < len(seq) and seq[i] == "G":
            ctx = "CG"
        elif i + 2 < len(seq) and seq[i + 1] != "G":
            ctx = "CHG" if seq[i + 1] != "A" else "CHH"
        else:
            continue
        x = r.randrange(3, 30)
        m = int(x * r.random())
        lines.append(f"{CHROM1}\t{i}\t.\t{ctx}\t{round(m / x, 3)}\t{x - m}\t0\t{m}\t{x}")
    path.write_text("\n".join(lines[:400]) + "\n")


def write_genome_file(path: Path) -> None:
    path.write_text(f"{CHROM1}\t{SIZE1}\tchromosome V (synthetic)\n{CHROM2}\t{SIZE2}\tchromosome M (synthetic)\n")


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    import sys

    sys.path.insert(0, str(ROOT))
    r = rng()
    g = make_genome(r)
    g, annot = place_genes(g, r)
    write_fasta(OUT / "genome.fa", g)
    (OUT / "genome2.fa").write_text(f">{CHROM2} second chromosome only\n"
                                     + "\n".join(g[CHROM2][i:i + 60] for i in range(0, len(g[CHROM2]), 60)) + "\n")
    (OUT / "genes.fa").write_text("".join(
        f">{a['name']} gene on {a['chrom']}:{a['start']}-{a['end']}{a['strand']}\n"
        + "\n".join((g[a['chrom']][a['start'] - 1:a['end']] if a['strand'] == "+"
                     else revcomp(g[a['chrom']][a['start'] - 1:a['end']]))[i:i + 60]
                    for i in range(0, a['end'] - a['start'] + 1, 60)) + "\n"
        for a in annot))
    write_gff(OUT / "annotation.gff", annot)
    write_bed(OUT / "regions.bed", annot, r)
    write_targets(OUT / "targets.bed", r)
    variants = make_variants(g, annot, r)
    refs, r1, r2, mut = write_reads(g, variants, r)
    write_fastq(OUT / "reads_1.fastq", r1)
    write_fastq(OUT / "reads_2.fastq", r2)
    write_fastq(OUT / "reads_single.fastq", (r1 + r2)[:120])
    (OUT / "alignments.sam").write_text(write_sam(refs, r1, r2, r, mut))
    samples = ["sample_A", "sample_B", "sample_C", "sample_D", "sample_E", "sample_F"]
    write_vcf(OUT / "variants.vcf", variants, samples, r)
    write_counts(OUT / "counts.tsv", annot, samples, r)
    write_gmt(OUT / "genesets.gmt", annot)
    write_motif(OUT / "motif.pfm")
    write_motifs_jaspar(OUT / "motifs_jaspar.txt")
    write_proteins(OUT / "proteins.faa", g, annot)
    write_msa(OUT / "msa.fasta", g, annot, r)
    write_tree(OUT / "tree.nwk")
    write_pairs(OUT / "contacts.pairs", r)
    write_bedgraph(OUT / "coverage.bedgraph", annot, r)
    write_pgm(OUT / "cells.pgm", r)
    write_peptides(OUT / "peptides.tsv", r)
    write_pheno(OUT / "phenotypes.tsv", r, samples)
    write_ontology(OUT / "phenotype.obo")
    write_str(OUT / "str_profile.tsv", g, r)
    write_methylation(OUT / "methylation.cxcg", g, r)
    write_genome_file(OUT / "genome.txt")
    (OUT / "two_bed.bed").write_text((OUT / "regions.bed").read_text().split("\n")[0] + "\n"
                                     + (OUT / "targets.bed").read_text().split("\n")[0] + "\n")
    (OUT / "taxmap.tsv").write_text("accession\tlineage\n" + "\n".join(
        f"ACC{i:03d}\tk__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Enterobacterales;"
        f"f__Enterobacteriaceae;g__Escherichia;s__Escherichia coli str{i}"
        for i in range(1, 26)))
    (OUT / "reference_proteins.faa").write_text(
        (OUT / "proteins.faa").read_text().replace("geneA", "ACC001")
        .replace("geneB", "ACC002").replace("geneC", "ACC003"))
    counts = {p.name: len(p.read_text().splitlines()) for p in sorted(OUT.iterdir())}
    print(f"wrote {len(counts)} example files to {OUT}")
    for k, v in counts.items():
        print(f"  {k:26s} {v:7d} lines")


if __name__ == "__main__":
    main()
