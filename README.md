# Chroma Titan — a Galaxy-style bioinformatics tool panel

**710 tools, 77 panel sections, zero infrastructure.**

Chroma Titan reproduces the tool panel of [usegalaxy.org](https://usegalaxy.org/) as a
self-contained Python project: a Streamlit app that looks and behaves like Galaxy's
"tool panel → form → run → view output" workflow, the same tools callable from a CLI,
and a pytest suite that executes every single one of them.

Unlike a real Galaxy deployment, a tool here is a **plain Python function** rather than a
`<tool>` XML wrapper plus a job script. There is no job queue, no conda resolvers, no
`samtools` binary and no scipy — the engines in `chroma_titan/core/` are implemented with
numpy/pandas/matplotlib so the whole thing runs in a browser tab, a Colab notebook or a
`python -m chroma_titan` call. It is meant for teaching and for exploring what a Galaxy
installation actually offers, tool by tool.

```bash
pip install -r requirements.txt
streamlit run app.py                 # the panel, in your browser
python -m chroma_titan check         # run the bundled example of all 710 tools
python -m chroma_titan sections      # tool counts per panel section
```

## The panel

`chroma_titan/panel.py` encodes the taxonomy that usegalaxy.org is deployed with (the
groups and section ids come from `panel_views/all_tools.yml` in
[`galaxyproject/usegalaxy-playbook`](https://github.com/galaxyproject/usegalaxy-playbook)),
including the specialised panel views — *Microbiology (MicroGalaxy)*, *Single Cell*,
*Imaging*, *Vertebrate Genomes Project*, *Biodiversity Genomics*. Every section of the
"All Tools" view is populated:

| Panel group | Sections | Tools |
|---|---|---|
| Get Data and Collections | 4 | 38 |
| General Text Tools | 4 | 65 |
| Genomic File Manipulation | 8 | 200 |
| Common Genomics Tools | 2 | 14 |
| Genomics Analysis | 18 | 148 |
| Statistics and Visualization | 4 | 43 |
| Genomics Toolkits | 22 | 109 |
| Domain Tools | 15 | 93 |
| **total** | **77** | **710** |

So you will find `bedtools_*`, `bcftools_*`, `samtools_*`, `seqtk_*`, `picard_*`,
`deeptools_*`, `emboss_*`, `mothur_*`, `qiime2_*`, `scanpy_*`, `seurat_*`, `monocle3_*`,
`presto_*`, `hicexplorer_*`, `rseqc_*`, `gemini_*`, `bbtools_*`, … next to the plain
"Filter and Sort", "Join, Subtraction and Group", "Phylogenetics", "Protein Structure",
"Proteomics", "Metabolomics", "ChemicalToolBox", "Pharmacology", "Multi-Omics",
"Imaging", "Climate analysis", "GIS data handling" and satellite-index sections.

## Using the app

* **Sidebar** – pick a *panel view*, then narrow with `Group → Section → Tool`, or type a
  search term (`coverage`, `motif`, `NDVI`, `UMAP`, …). Recently used tools stay pinned at
  the bottom of the sidebar.
* **Tool form** – generated from the tool's declared inputs: FASTA/FASTQ/BED/GFF/VCF/SAM,
  tables, images and Hi-C pairs can come from the bundled example data, from an upload, or
  be pasted in directly. Numbers, sliders, tick-boxes, pick-lists and column multi-selects
  map to the matching Galaxy parameter types.
* **Run** – the output panel shows the job message, the result table (with a CSV download),
  figures, text/sequence output, and a *Statistics and extra datasets* section for whatever
  side-data the tool returned (matrix, newick, per-bin lists, …). Failures render as an
  expanded traceback, like a Galaxy job error box.
* **All tools** – a filterable table of the whole registry, with a CSV download and a
  jump-to-form selector. **Datasets** – the 30 bundled files with previews and downloads.

## Using the CLI

```bash
python -m chroma_titan search "intersect"
python -m chroma_titan info  bedtools_intersect
python -m chroma_titan run   fasta_gc_content --src=examples/genes.fa
python -m chroma_titan run   bedtools_merge --src examples/regions.bed --d 100 --json
python -m chroma_titan docs                    # regenerate docs/TOOLS.md
python -m chroma_titan check -v                # self-test: 710 ok, 0 failed
python -m chroma_titan examples                # list the bundled data
```

## Writing a tool

One decorated function per tool, grouped by the panel section it belongs to:

```python
from chroma_titan.tools._common import T, fa, io, table, GENES_FA

@T("my_gc_summary", "GC summary per record", "fasta_fastq", "table",
   [fa("src", GENES_FA, "FASTA file"), number("min_gc", 40.0, "Report above (%)")],
   ex={"src": GENES_FA, "min_gc": 40.0}, up="iuc/fasta-gc",
   tags=("gc", "qc"), summary="Report GC content per sequence and drop the low ones.")
def my_gc_summary(src, min_gc=40.0):
    """GC content per record, filtered."""
    rows = [{"id": r.id, "length": len(r.seq), "gc": round(seq_gc(r.seq), 3)}
            for r in io.parse_fasta(io.as_text(src))]
    keep = [r for r in rows if r["gc"] >= float(min_gc)]
    return table(keep, f"{len(keep)}/{len(rows)} records at or above {min_gc}% GC")
```

Rules the registry enforces, and that `tests/test_tools.py` checks for every tool:

* the tool id is unique, snake_case, and is how `registry.run(id, **values)` dispatches;
* `section` must be a real `usegalaxy.org` section id from `chroma_titan.panel`;
* every parameter needs an `Input` (or a default, from which the widget spec is inferred);
* `ex` is the runnable example — `python -m chroma_titan check` executes it for all tools,
  so a tool without a working example cannot be merged;
* results are built with `table() / figure() / text() / values_message()`; anything put in
  `result["stats"]` shows up in the app's "Statistics and extra datasets" expander and in
  the CLI's `--json` output.

## Layout

```
app.py                      Streamlit front end: sidebar panel + tool forms + results
chroma_titan/
  registry.py               @tool decorator, Input dataclass, search, coercion, dispatch
  panel.py                  usegalaxy.org groups, sections, panel views
  datasets.py               catalogue + loaders/previews for data/examples
  cli.py                    list / search / sections / info / run / examples / check / docs
  core/                     engines: io, seq, genome, tables, variants, align, motif,
                            phylo, protein, stats, ml, bam, image, plot
  tools/                    the tools, one module per area of the panel
    text.py fasta.py fastq_qc.py bed.py sam.py vcf.py alignmap.py annotation.py
    rnaseq.py statsml.py omics_misc.py prochem.py phyloimag.py
data/examples/              30 small example files (genome, reads, BED/GFF/VCF/SAM,
                            counts, phenotypes, Hi-C pairs, PGM image, trees, motifs …)
scripts/make_examples.py    regenerates every bundled example deterministically
tests/                      test_tools.py (registry, panel, CLI, all examples),
                            test_core.py (engines), test_app.py (Streamlit AppTest smoke test)
docs/TOOLS.md               generated reference: every section, tool and parameter
```

| module | tools | sections covered |
|---|---|---|
| `text.py` | 88 | 8 |
| `fasta.py` | 71 | 4 |
| `omics_misc.py` | 64 | 15 |
| `annotation.py` | 54 | 7 |
| `alignmap.py` | 53 | 7 |
| `vcf.py` | 52 | 8 |
| `sam.py` | 51 | 4 |
| `bed.py` | 50 | 6 |
| `prochem.py` | 49 | 10 |
| `phyloimag.py` | 48 | 7 |
| `rnaseq.py` | 44 | 8 |
| `fastq_qc.py` | 43 | 3 |
| `statsml.py` | 43 | 4 |

## Tests

```bash
pytest                    # 84 fast tests: registry invariants, panel coverage, engines, app
pytest -m slow            # + executes the example of every registered tool
python -m chroma_titan check
```

The fast suite verifies structural things that keep the UI honest — unique ids, real panel
sections, no section left empty, every tool has a runnable example, widget defaults are
renderable (a default outside `min/max`, a choice missing from its options, a missing
bundled file all fail the test). `tests/test_app.py` boots the Streamlit app with
`streamlit.testing.v1.AppTest`, jumps to eight representative tools, clicks *Run tool* and
asserts the results render.

## Honest limitations

* Tools are **teaching implementations**. A `bcftools_view` or `deseq_lite` here implements
  the observable behaviour (filters, column semantics, ranking, the numbers that matter) in
  a few dozen lines of pandas/numpy, not the upstream algorithm; p-values come from
  `chroma_titan/core/stats.py` (no scipy), so treat results as illustrative.
* SAM/BAM input is read as **text SAM** (`chroma_titan.core.io.parse_sam`); binary BAM, CRAM
  and indexed queries are out of scope. Big-wig/BigBed are handled as bedGraph.
* Coordinates are 0-based BED / 1-based GFF·VCF, matching Galaxy's own conventions.
* No history, no jobs, no user accounts, no data transfer between tools — that is exactly
  what a real Galaxy server adds on top of this kind of registry.

## Licence

Apache License 2.0 — see [LICENSE](LICENSE). Galaxy and usegalaxy.org are trademarks of the
Galaxy Project / Biotech Partners; this repository is an independent educational
re-implementation of the tool panel and is not affiliated with them.
