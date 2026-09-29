"""The tool-panel taxonomy of https://usegalaxy.org/.

The section ids, the section groups ("labels") and the panel *views* below were
read from the Galaxy configuration that the public usegalaxy.org server is
deployed with (``panel_views/all_tools.yml`` in the
``galaxyproject/usegalaxy-playbook`` repository).  ``chroma_titan`` groups its
tools with exactly the same structure, so the sidebar of the web app behaves
like the Galaxy tool panel.
"""

from __future__ import annotations

#: pretty label for every Galaxy panel section id
SECTION_LABELS: dict[str, str] = {
    "get_data": "Get Data",
    "send_data": "Send Data",
    "collection_operations": "Collection Operations",
    "expression_tools": "Expression Tools",
    "text_manipulation": "Text Manipulation",
    "filter_and_sort": "Filter and Sort",
    "join__subtract_and_group": "Join, Subtraction and Group Operations",
    "datamash": "Datamash: Aggregation and Statistics",
    "fasta_fastq": "FASTA/FASTQ",
    "fastq_quality_control": "FASTQ Quality Control",
    "sam_bam": "SAM/BAM",
    "bed": "BED",
    "vcf_bcf": "VCF/BCF",
    "nanopore": "Nanopore",
    "convert_formats": "Convert Formats",
    "lift_over": "Lift-Over",
    "operate_on_genomic_intervals": "Operate on Genomic Intervals",
    "fetch_sequences_alignments": "Fetch Sequences and Alignments",
    "assembly": "Assembly",
    "annotation": "Annotation",
    "mapping": "Mapping",
    "variant_calling": "Variant Calling",
    "chip_seq": "Chromatin Accessibility (ChIP-seq/ATAC)",
    "rna_seq": "RNA Sequencing (RNA-Seq)",
    "multiple_alignments": "Multiple Alignments",
    "phenotype_association": "Phenotype Association",
    "phylogenetics": "Phylogenetics",
    "evolution": "Evolution",
    "regional_variation": "Regional Variation",
    "chromosome_conformation": "Chromosome Conformation (Hi-C)",
    "transposon_insertion_sequencing": "Transposon Insertion Sequencing",
    "protein_modeling": "Protein Structure",
    "genome_editing": "Genome Editing (CRISPR)",
    "biodiversity_data_exploration": "Biodiversity Data Exploration",
    "sequence_contamination_filtering": "Sequence Contamination Filtering",
    "genome_diversity": "Genome Diversity",
    "statistics": "Statistics",
    "machine_learning": "Machine Learning",
    "graph_display_data": "Graph/Display Data",
    "interactive_tools": "Interactive Tools (Visualization)",
    "str_fm__microsatellite_analysis": "STR-FM: Microsatellite Analysis",
    "mothur": "Mothur (16S rRNA)",
    "qiime2": "QIIME 2",
    "picard": "Picard",
    "deeptools": "deepTools",
    "emboss": "EMBOSS",
    "ncbi_blast_": "NCBI BLAST+",
    "rseqc": "RSeQC",
    "mimodd": "MiModD",
    "hca-scanpy": "Scanpy (Single Cell)",
    "hicexplorer": "HiCExplorer",
    "du_novo": "Du Novo (Duplex Consensus)",
    "seqtk": "SeqTK",
    "monocle3": "Monocle3",
    "gemini": "Gemini",
    "bbtools": "BBTools",
    "iwtomics": "IWTomics",
    "presto": "PreSTO",
    "planttribes": "PlantTribes",
    "motif": "Motif Tools",
    "sccaf": "sccaf",
    "seurat": "Seurat",
    "virology": "Virology",
    "metagenomic_analysis": "Metagenomic Analysis",
    "single_cell": "Single Cell Analysis",
    "import/manipulate_sc_data": "Import and Manipulate Single Cell Data",
    "imaging": "Imaging",
    "proteomics": "Proteomics",
    "metabolomics": "Metabolomics",
    "chemicaltoolbox": "ChemicalToolBox",
    "pharmacology": "Pharmacology",
    "multiomics": "Multi-Omics",
    "climate_analysis": "Climate Analysis",
    "gis_data_handling": "GIS Data Handling",
    "epigenetics": "Epigenetics",
    "spatial": "Spatial",
    "compute_indicators_for_satellite_remote_sensing": "Indicators for Satellite Remote Sensing",
}

#: the "All Tools" panel view, verbatim structure from usegalaxy.org
ALL_TOOLS_GROUPS: list[tuple[str, list[str]]] = [
    ("Get Data and Collections", [
        "get_data", "send_data", "collection_operations", "expression_tools"]),
    ("General Text Tools", [
        "text_manipulation", "filter_and_sort", "join__subtract_and_group", "datamash"]),
    ("Genomic File Manipulation", [
        "fasta_fastq", "fastq_quality_control", "sam_bam", "bed", "vcf_bcf",
        "nanopore", "convert_formats", "lift_over"]),
    ("Common Genomics Tools", [
        "operate_on_genomic_intervals", "fetch_sequences_alignments"]),
    ("Genomics Analysis", [
        "assembly", "annotation", "mapping", "variant_calling", "chip_seq", "rna_seq",
        "multiple_alignments", "phenotype_association", "phylogenetics", "evolution",
        "regional_variation", "chromosome_conformation",
        "transposon_insertion_sequencing", "protein_modeling", "genome_editing",
        "biodiversity_data_exploration", "sequence_contamination_filtering",
        "genome_diversity"]),
    ("Statistics and Visualization", [
        "statistics", "machine_learning", "graph_display_data", "interactive_tools"]),
    ("Genomics Toolkits", [
        "str_fm__microsatellite_analysis", "mothur", "qiime2", "picard", "deeptools",
        "emboss", "ncbi_blast_", "rseqc", "hca-scanpy", "hicexplorer", "seqtk",
        "bbtools", "motif", "gemini", "du_novo", "iwtomics", "presto", "planttribes",
        "sccaf", "seurat", "monocle3", "mimodd"]),
    ("Domain Tools", [
        "virology", "metagenomic_analysis", "single_cell", "import/manipulate_sc_data",
        "imaging", "proteomics", "metabolomics", "chemicaltoolbox", "pharmacology",
        "multiomics", "climate_analysis", "gis_data_handling", "epigenetics", "spatial",
        "compute_indicators_for_satellite_remote_sensing"]),
]

#: extra panel views offered by usegalaxy.org (Microbiology, Single Cell, ...).
PANEL_VIEWS: dict[str, dict] = {
    "all_tools": {
        "name": "All Tools",
        "description": "Every tool in the panel, grouped the way Galaxy groups them.",
        "groups": None,  # None == all groups
    },
    "microgalaxy": {
        "name": "Microbiology (MicroGalaxy)",
        "description": "Reduced panel for microbial genomics and metagenomics.",
        "groups": {
            "General Text Tools": ["text_manipulation", "filter_and_sort",
                                   "join__subtract_and_group"],
            "Genomic File Manipulation": ["convert_formats", "fasta_fastq",
                                         "fastq_quality_control", "sam_bam", "bed"],
            "Common Genomics Tools": ["operate_on_genomic_intervals"],
            "Genomics Analysis": ["annotation", "multiple_alignments", "assembly",
                                  "mapping"],
            "Genomics Toolkits": ["picard", "deeptools", "ncbi_blast_", "mothur",
                                  "qiime2"],
            "Metagenomics": ["metagenomic_analysis", "epigenetics"],
            "Metabolomics": ["metabolomics"],
            "Statistics and Visualization": ["graph_display_data", "interactive_tools",
                                             "statistics"],
        },
    },
    "singlecell": {
        "name": "Single Cell",
        "description": "Scanpy/Seurat/Monocle3 flavoured single-cell panel.",
        "groups": {
            "General Tools": ["text_manipulation", "filter_and_sort", "fasta_fastq",
                              "fastq_quality_control", "bed", "graph_display_data",
                              "get_data"],
            "Single Cell Tools": ["single_cell", "hca-scanpy", "seurat", "monocle3",
                                  "import/manipulate_sc_data", "spatial"],
            "Statistics and Visualization": ["statistics", "machine_learning",
                                            "interactive_tools"],
        },
    },
    "imaging": {
        "name": "Imaging",
        "description": "Bio-image analysis panel (2D/3D image processing).",
        "groups": {
            "General Text Tools": ["text_manipulation", "filter_and_sort", "datamash"],
            "Imaging": ["imaging"],
            "Statistics and Visualization": ["statistics", "machine_learning",
                                             "graph_display_data", "interactive_tools"],
        },
    },
    "biodiversity": {
        "name": "Biodiversity Genomics",
        "description": "GBF/biodiversity genome workflow panel.",
        "groups": {
            "General Text Tools": ["text_manipulation", "filter_and_sort",
                                   "join__subtract_and_group"],
            "Genomic File Manipulation": ["convert_formats", "fasta_fastq",
                                         "fastq_quality_control", "sam_bam", "bed"],
            "Common Genomics Tools": ["operate_on_genomic_intervals",
                                      "fetch_sequences_alignments"],
            "Genomics Analysis": ["assembly", "annotation", "multiple_alignments",
                                  "mapping", "phenotype_association", "evolution",
                                  "biodiversity_data_exploration",
                                  "sequence_contamination_filtering",
                                  "genome_diversity", "phylogenetics"],
            "Genomics Toolkits": ["picard", "deeptools", "ncbi_blast_"],
            "Domain Tools": ["metagenomic_analysis", "genome_editing"],
        },
    },
}

# ---------------------------------------------------------------------------
_SECTION_TO_GROUP: dict[str, str] = {}
for _label, _secs in ALL_TOOLS_GROUPS:
    for _s in _secs:
        _SECTION_TO_GROUP.setdefault(_s, _label)


def group_of(section: str) -> str:
    """Panel label a section lives under on usegalaxy.org."""
    return _SECTION_TO_GROUP.get(section, "Other Tools")


def section_label(section: str) -> str:
    return SECTION_LABELS.get(section, section.replace("_", " ").title())


def ordered_sections(group: str | None = None) -> list[str]:
    """All section ids, in panel order (optionally restricted to one group)."""
    out: list[str] = []
    for label, secs in ALL_TOOLS_GROUPS:
        if group and label != group:
            continue
        for s in secs:
            if s not in out:
                out.append(s)
    if group is None:
        from chroma_titan.registry import REGISTRY

        extra = sorted({t.section for t in REGISTRY.values()} - set(out))
        out.extend(extra)
    return out


def panel_layout(view: str = "all_tools") -> list[tuple[str, list[tuple[str, list[str]]]]]:
    """``[(group label, [(section label, [tool ids])])]`` for a panel view."""
    from chroma_titan.registry import REGISTRY

    spec = PANEL_VIEWS.get(view) or PANEL_VIEWS["all_tools"]
    groups: list[tuple[str, list[tuple[str, list[str]]]]] = []
    source = spec["groups"]
    if source is None:
        allowed = {s for _, secs in ALL_TOOLS_GROUPS for s in secs}
        pairs = [(label, secs) for label, secs in ALL_TOOLS_GROUPS]
        pairs.append(("Other Tools", ordered_sections() and [
            s for s in ordered_sections() if s not in allowed]))
    else:
        pairs = list(source.items())
        allowed = None
    for label, secs in pairs:
        rows: list[tuple[str, list[str]]] = []
        for sec in secs:
            ids = [t.id for t in REGISTRY.values() if t.section == sec]
            if ids:
                rows.append((section_label(sec), sorted(ids)))
        if rows:
            groups.append((label, rows))
    # dedupe tools that appear in several sections of a view, keeping first slot
    seen: set[str] = set()
    deduped: list[tuple[str, list[tuple[str, list[str]]]]] = []
    for label, rows in groups:
        newrows = []
        for seclabel, ids in rows:
            keep = [i for i in ids if i not in seen]
            seen.update(keep)
            if keep:
                newrows.append((seclabel, keep))
        if newrows:
            deduped.append((label, newrows))
    return deduped


def view_ids() -> list[str]:
    return list(PANEL_VIEWS)
