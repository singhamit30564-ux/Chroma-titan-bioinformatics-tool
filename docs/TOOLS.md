# Chroma Titan tool reference

710 tools, grouped like the tool panel of [usegalaxy.org](https://usegalaxy.org/).
This file is generated: `python -m chroma_titan docs`.

| Section | Group | Tools |
|---|---|---|
| [`get_data`](#get_data) | Get Data and Collections | 4 |
| [`send_data`](#send_data) | Get Data and Collections | 3 |
| [`collection_operations`](#collection_operations) | Get Data and Collections | 10 |
| [`expression_tools`](#expression_tools) | Get Data and Collections | 21 |
| [`text_manipulation`](#text_manipulation) | General Text Tools | 30 |
| [`filter_and_sort`](#filter_and_sort) | General Text Tools | 16 |
| [`join__subtract_and_group`](#join__subtract_and_group) | General Text Tools | 8 |
| [`datamash`](#datamash) | General Text Tools | 11 |
| [`fasta_fastq`](#fasta_fastq) | Genomic File Manipulation | 49 |
| [`fastq_quality_control`](#fastq_quality_control) | Genomic File Manipulation | 30 |
| [`sam_bam`](#sam_bam) | Genomic File Manipulation | 33 |
| [`bed`](#bed) | Genomic File Manipulation | 32 |
| [`vcf_bcf`](#vcf_bcf) | Genomic File Manipulation | 27 |
| [`nanopore`](#nanopore) | Genomic File Manipulation | 11 |
| [`convert_formats`](#convert_formats) | Genomic File Manipulation | 14 |
| [`lift_over`](#lift_over) | Genomic File Manipulation | 4 |
| [`operate_on_genomic_intervals`](#operate_on_genomic_intervals) | Common Genomics Tools | 9 |
| [`fetch_sequences_alignments`](#fetch_sequences_alignments) | Common Genomics Tools | 5 |
| [`assembly`](#assembly) | Genomics Analysis | 9 |
| [`annotation`](#annotation) | Genomics Analysis | 24 |
| [`mapping`](#mapping) | Genomics Analysis | 12 |
| [`variant_calling`](#variant_calling) | Genomics Analysis | 7 |
| [`chip_seq`](#chip_seq) | Genomics Analysis | 5 |
| [`rna_seq`](#rna_seq) | Genomics Analysis | 8 |
| [`multiple_alignments`](#multiple_alignments) | Genomics Analysis | 15 |
| [`phenotype_association`](#phenotype_association) | Genomics Analysis | 3 |
| [`phylogenetics`](#phylogenetics) | Genomics Analysis | 12 |
| [`evolution`](#evolution) | Genomics Analysis | 8 |
| [`regional_variation`](#regional_variation) | Genomics Analysis | 2 |
| [`chromosome_conformation`](#chromosome_conformation) | Genomics Analysis | 6 |
| [`transposon_insertion_sequencing`](#transposon_insertion_sequencing) | Genomics Analysis | 3 |
| [`protein_modeling`](#protein_modeling) | Genomics Analysis | 11 |
| [`genome_editing`](#genome_editing) | Genomics Analysis | 6 |
| [`biodiversity_data_exploration`](#biodiversity_data_exploration) | Genomics Analysis | 4 |
| [`sequence_contamination_filtering`](#sequence_contamination_filtering) | Genomics Analysis | 4 |
| [`genome_diversity`](#genome_diversity) | Genomics Analysis | 9 |
| [`statistics`](#statistics) | Statistics and Visualization | 21 |
| [`machine_learning`](#machine_learning) | Statistics and Visualization | 12 |
| [`graph_display_data`](#graph_display_data) | Statistics and Visualization | 6 |
| [`interactive_tools`](#interactive_tools) | Statistics and Visualization | 4 |
| [`str_fm__microsatellite_analysis`](#str_fm__microsatellite_analysis) | Genomics Toolkits | 5 |
| [`mothur`](#mothur) | Genomics Toolkits | 7 |
| [`qiime2`](#qiime2) | Genomics Toolkits | 6 |
| [`picard`](#picard) | Genomics Toolkits | 5 |
| [`deeptools`](#deeptools) | Genomics Toolkits | 8 |
| [`emboss`](#emboss) | Genomics Toolkits | 7 |
| [`ncbi_blast_`](#ncbi_blast_) | Genomics Toolkits | 8 |
| [`rseqc`](#rseqc) | Genomics Toolkits | 10 |
| [`hca-scanpy`](#hca-scanpy) | Genomics Toolkits | 4 |
| [`hicexplorer`](#hicexplorer) | Genomics Toolkits | 3 |
| [`seqtk`](#seqtk) | Genomics Toolkits | 10 |
| [`bbtools`](#bbtools) | Genomics Toolkits | 8 |
| [`motif`](#motif) | Genomics Toolkits | 13 |
| [`gemini`](#gemini) | Genomics Toolkits | 3 |
| [`du_novo`](#du_novo) | Genomics Toolkits | 1 |
| [`iwtomics`](#iwtomics) | Genomics Toolkits | 1 |
| [`presto`](#presto) | Genomics Toolkits | 2 |
| [`planttribes`](#planttribes) | Genomics Toolkits | 2 |
| [`sccaf`](#sccaf) | Genomics Toolkits | 1 |
| [`seurat`](#seurat) | Genomics Toolkits | 2 |
| [`monocle3`](#monocle3) | Genomics Toolkits | 2 |
| [`mimodd`](#mimodd) | Genomics Toolkits | 1 |
| [`virology`](#virology) | Domain Tools | 4 |
| [`metagenomic_analysis`](#metagenomic_analysis) | Domain Tools | 13 |
| [`single_cell`](#single_cell) | Domain Tools | 8 |
| [`import/manipulate_sc_data`](#import/manipulate_sc_data) | Domain Tools | 5 |
| [`imaging`](#imaging) | Domain Tools | 8 |
| [`proteomics`](#proteomics) | Domain Tools | 12 |
| [`metabolomics`](#metabolomics) | Domain Tools | 7 |
| [`chemicaltoolbox`](#chemicaltoolbox) | Domain Tools | 7 |
| [`pharmacology`](#pharmacology) | Domain Tools | 5 |
| [`multiomics`](#multiomics) | Domain Tools | 5 |
| [`climate_analysis`](#climate_analysis) | Domain Tools | 3 |
| [`gis_data_handling`](#gis_data_handling) | Domain Tools | 4 |
| [`epigenetics`](#epigenetics) | Domain Tools | 6 |
| [`spatial`](#spatial) | Domain Tools | 2 |
| [`compute_indicators_for_satellite_remote_sensing`](#compute_indicators_for_satellite_remote_sensing) | Domain Tools | 4 |
| **total** |  | **710** |

## get_data

*Get Data* — in Galaxy group *Get Data and Collections*.

| Tool | What it does | Upstream Galaxy tool |
|---|---|---|
| [`get_data_concatenate_history`](#get_data_concatenate_history) | Stack several datasets vertically (list → single dataset). | `cat` |
| [`get_data_create_dataset`](#get_data_create_dataset) | Wrap pasted text into a dataset and report its basic properties. | `upload` |
| [`get_data_extract_column`](#get_data_extract_column) | ``cut -f`` producing a one-column dataset. | `extract` |
| [`get_data_format_probe`](#get_data_format_probe) | Sniff FASTA/FASTQ/VCF/SAM/BED/GFF/table from the first lines. | `peek` |

### get_data_concatenate_history

**Concatenate a set of same-format files** — Stack several datasets vertically (list → single dataset).

| parameter | kind | default | options |
|---|---|---|---|
| `files` | multi | `['examples/regions.bed', 'examples/targets.bed']` | examples/phenotypes.tsv, examples/regions.bed, examples/targets.bed, examples/counts.tsv |
| `keep_headers` | bool | `True` | — |

Galaxy Tool Shed: `cat`

### get_data_create_dataset

**Create a dataset from pasted text** — Wrap pasted text into a dataset and report its basic properties.

| parameter | kind | default | options |
|---|---|---|---|
| `content` | code | `` | — |

Galaxy Tool Shed: `upload`

### get_data_extract_column

**Extract a single column as a new dataset** — ``cut -f`` producing a one-column dataset.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/counts.tsv` | — |
| `column` | text | `1` | — |
| `header` | bool | `False` | — |

Galaxy Tool Shed: `extract`

### get_data_format_probe

**Identify the format of a dataset** — Sniff FASTA/FASTQ/VCF/SAM/BED/GFF/table from the first lines.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/genome.fa` | — |

Galaxy Tool Shed: `peek`

## send_data

*Send Data* — in Galaxy group *Get Data and Collections*.

| Tool | What it does | Upstream Galaxy tool |
|---|---|---|
| [`send_data_download_table`](#send_data_download_table) | Re-serialise a table for download (TSV/CSV/Markdown/HTML/JSON). | `export` |
| [`send_data_paste_to_column`](#send_data_paste_to_column) | Write one column out (e.g. a gene list for a downstream enrichment tool). | `export` |
| [`send_data_summary_report`](#send_data_summary_report) | Compact human-readable summary (size, format, preview). | `report` |

### send_data_download_table

**Convert any dataset to a chosen text format** — Re-serialise a table for download (TSV/CSV/Markdown/HTML/JSON).

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/counts.tsv` | — |
| `format` | choice | `tsv` | tsv, csv, markdown, html, json, yaml-ish |

Galaxy Tool Shed: `export`

### send_data_paste_to_column

**Send a column back as a standalone dataset** — Write one column out (e.g. a gene list for a downstream enrichment tool).

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/counts.tsv` | — |
| `column` | text | `1` | — |

Galaxy Tool Shed: `export`

### send_data_summary_report

**Text summary report of a dataset** — Compact human-readable summary (size, format, preview).

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/phenotypes.tsv` | — |
| `preview_lines` | int | `6` | — |

Galaxy Tool Shed: `report`

## collection_operations

*Collection Operations* — in Galaxy group *Get Data and Collections*.

| Tool | What it does | Upstream Galaxy tool |
|---|---|---|
| [`collection_count_elements`](#collection_count_elements) | Number of elements plus per-element line/field counts. | `count` |
| [`collection_extract_by_index`](#collection_extract_by_index) | Pull specific elements out of a list collection. | `extract` |
| [`collection_filter_elements`](#collection_filter_elements) | Keep collection elements whose identifier matches a regexp. | `filter_from_list` |
| [`collection_group_by_label`](#collection_group_by_label) | Create a list:paired-style grouping from a label column. | `group_with_pandas` |
| [`collection_list_from_column`](#collection_list_from_column) | Build the equivalent of a Galaxy list collection from one column. | `list_operations` |
| [`collection_merge`](#collection_merge) | Combine two collections with a set-like strategy. | `merge` |
| [`collection_rename_elements`](#collection_rename_elements) | Apply prefix/suffix and an optional strip regexp to element names. | `rename_by_pattern` |
| [`collection_sort_elements`](#collection_sort_elements) | Order list elements deterministically. | `sort_lists` |
| [`collection_unzip`](#collection_unzip) | Undo a zip: two output lists from a two-column table. | `unzip` |
| [`collection_zip_pairs`](#collection_zip_pairs) | Pair up elements from two collections positionally. | `zip_longest` |

### collection_count_elements

**Count elements and their sizes** — Number of elements plus per-element line/field counts.

| parameter | kind | default | options |
|---|---|---|---|
| `names` | file | `` | — |

Galaxy Tool Shed: `count`

### collection_extract_by_index

**Extract elements by index or range** — Pull specific elements out of a list collection.

| parameter | kind | default | options |
|---|---|---|---|
| `names` | file | `` | — |
| `indices` | text | `1,3,5-6` | — |
| `cyclic` | bool | `False` | — |

Galaxy Tool Shed: `extract`

### collection_filter_elements

**Filter collection elements by pattern** — Keep collection elements whose identifier matches a regexp.

| parameter | kind | default | options |
|---|---|---|---|
| `names` | file | `` | — |
| `pattern` | text | `gene` | — |
| `invert` | bool | `False` | — |

Galaxy Tool Shed: `filter_from_list`

### collection_group_by_label

**Group collection elements into nested lists** — Create a list:paired-style grouping from a label column.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/counts.tsv` | — |
| `column` | text | `1` | — |

Galaxy Tool Shed: `group_with_pandas`

### collection_list_from_column

**Extract an element list from a dataset** — Build the equivalent of a Galaxy list collection from one column.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/counts.tsv` | — |
| `column` | text | `1` | — |
| `unique` | bool | `True` | — |
| `header` | bool | `True` | — |

Galaxy Tool Shed: `list_operations`

### collection_merge

**Merge collections with deduplication** — Combine two collections with a set-like strategy.

| parameter | kind | default | options |
|---|---|---|---|
| `first` | file | `` | — |
| `second` | file | `` | — |
| `strategy` | choice | `union` | union, intersect, first_only |
| `unique` | bool | `True` | — |

Galaxy Tool Shed: `merge`

### collection_rename_elements

**Rename collection elements** — Apply prefix/suffix and an optional strip regexp to element names.

| parameter | kind | default | options |
|---|---|---|---|
| `names` | file | `` | — |
| `prefix` | text | `sample_` | — |
| `suffix` | text | `` | — |
| `strip_pattern` | text | `` | — |

Galaxy Tool Shed: `rename_by_pattern`

### collection_sort_elements

**Sort collection elements by identifier** — Order list elements deterministically.

| parameter | kind | default | options |
|---|---|---|---|
| `names` | file | `` | — |
| `by` | choice | `name` | name, length, size, index |
| `reverse` | bool | `False` | — |

Galaxy Tool Shed: `sort_lists`

### collection_unzip

**Split a paired collection into two lists** — Undo a zip: two output lists from a two-column table.

| parameter | kind | default | options |
|---|---|---|---|
| `pairs` | file | `examples/counts.tsv` | — |
| `column_left` | int | `1` | — |
| `column_right` | int | `2` | — |

Galaxy Tool Shed: `unzip`

### collection_zip_pairs

**Zip two lists into pairs** — Pair up elements from two collections positionally.

| parameter | kind | default | options |
|---|---|---|---|
| `first` | file | `` | — |
| `second` | file | `` | — |
| `fill` | bool | `False` | — |

Galaxy Tool Shed: `zip_longest`

## expression_tools

*Expression Tools* — in Galaxy group *Get Data and Collections*.

| Tool | What it does | Upstream Galaxy tool |
|---|---|---|
| [`apply_function_to_rows`](#apply_function_to_rows) | Reduce several numeric columns per row with a chosen function. | `apply` |
| [`column_arithmetic_pairwise`](#column_arithmetic_pairwise) | Element-wise maths on two columns of the same table. | `column_arithmetic` |
| [`compute_expression_column`](#compute_expression_column) | Galaxy's 'Compute' tool: evaluate an expression using column names. | `compute` |
| [`counts_add_calc_column`](#counts_add_calc_column) | Evaluate an arithmetic expression over the numeric columns of a table. | `Galaxy Add/redo calculation` |
| [`counts_clustering_heatmap`](#counts_clustering_heatmap) | Order rows by hierarchical clustering of the most variable genes and draw the matrix. | `pheatmap / ComplexHeatmap` |
| [`counts_cpm`](#counts_cpm) | Library-size normalisation to counts per million, optionally log transformed. | `edgeR cpm / DESeq2 vst` |
| [`counts_dedupe_merge`](#counts_dedupe_merge) | Collapse sample columns sharing a prefix (replicates) into one value per feature. | `collate / sum columns` |
| [`counts_feature_stats`](#counts_feature_stats) | Mean, median, spread, CV and detection rate for every gene. | `row stats (base R)` |
| [`counts_filter_low_expression`](#counts_filter_low_expression) | Drop features that are not expressed well enough for statistical testing. | `DESeq2 row filtering` |
| [`counts_fold_change`](#counts_fold_change) | Per-feature log2 fold change between two samples, columns or sample groups. | `DESeq2 results / edgeR` |
| [`counts_log_transform`](#counts_log_transform) | Variance-stabilising log transform with a configurable pseudocount. | `log2 transform (limma)` |
| [`counts_pca_samples`](#counts_pca_samples) | Project libraries onto principal components to reveal replicate structure. | `DESeq2 plotPCA / prcomp` |
| [`counts_quantile_normalize`](#counts_quantile_normalize) | Force every sample onto the same distribution before comparing groups. | `limma normalizeQuantiles` |
| [`counts_sample_correlation`](#counts_sample_correlation) | Heatmap of pairwise correlations between libraries, the batch-effect smoke test. | `DESeq2 sampleDistance / corrplot` |
| [`counts_summary`](#counts_summary) | Library sizes, genes detected, zero counts and the most abundant features. | `HTSeq-count / featureCounts output QC` |
| [`counts_tidy_long`](#counts_tidy_long) | Unroll the matrix into feature / sample / value rows for plotting tools. | `tidyr pivot_longer` |
| [`counts_tpm`](#counts_tpm) | Length-corrected transcripts per million with a FASTA of transcript sequences. | `salmon / featureCounts + TPM` |
| [`counts_zscore`](#counts_zscore) | Centre and scale each row (or column) so patterns are visible in a heatmap. | `pheatmap / row z-score` |
| [`missing_value_imputation`](#missing_value_imputation) | Fill NaNs column-wise (constant, mean/median or k-nearest-neighbour). | `imputation` |
| [`row_math`](#row_math) | Aggregate across each row of a numeric matrix. | `row_math` |
| [`scale_matrix`](#scale_matrix) | Standardise a numeric matrix along an axis. | `scale` |

### apply_function_to_rows

**Apply a small function to selected columns** — Reduce several numeric columns per row with a chosen function.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/counts.tsv` | — |
| `columns` | multi | `[]` | — |
| `function` | choice | `mean` | mean, median, sum, geomean, cv, range, log10sum, max-min_ratio |

Galaxy Tool Shed: `apply`

### column_arithmetic_pairwise

**Arithmetic between two numeric columns** — Element-wise maths on two columns of the same table.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/counts.tsv` | — |
| `left` | text | `1` | — |
| `right` | text | `2` | — |
| `operation` | choice | `add` | add, subtract, multiply, divide, ratio_log2, absolute_difference |

Galaxy Tool Shed: `column_arithmetic`

### compute_expression_column

**Add a computed column from an expression** — Galaxy's 'Compute' tool: evaluate an expression using column names.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/counts.tsv` | — |
| `name` | text | `ratio` | — |
| `expression` | text | `yield / biomass` | — |
| `header` | bool | `True` | — |

Galaxy Tool Shed: `compute`

### counts_add_calc_column

**Add a computed column to a table** — Evaluate an arithmetic expression over the numeric columns of a table.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/counts.tsv` | — |
| `name` | text | `ratio` | — |
| `expression` | text | `(sample_A + sample_B) / sample_D` | — |
| `digits` | int | `4` | — |

Galaxy Tool Shed: `Galaxy Add/redo calculation`

### counts_clustering_heatmap

**Clustered heatmap of features** — Order rows by hierarchical clustering of the most variable genes and draw the matrix.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/counts.tsv` | — |
| `top_features` | int | `20` | — |
| `method` | choice | `average` | average, complete, single, ward |
| `zscore_rows` | bool | `True` | — |

Galaxy Tool Shed: `pheatmap / ComplexHeatmap`

### counts_cpm

**Normalise counts to CPM** — Library-size normalisation to counts per million, optionally log transformed.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/counts.tsv` | — |
| `log` | bool | `False` | — |
| `digits` | int | `4` | — |

Galaxy Tool Shed: `edgeR cpm / DESeq2 vst`

### counts_dedupe_merge

**Merge replicate columns** — Collapse sample columns sharing a prefix (replicates) into one value per feature.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/counts.tsv` | — |
| `suffix` | text | `_rep` | — |
| `agg` | choice | `mean` | mean, sum, median, max |
| `digits` | int | `3` | — |

Galaxy Tool Shed: `collate / sum columns`

### counts_feature_stats

**Descriptive statistics per feature** — Mean, median, spread, CV and detection rate for every gene.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/counts.tsv` | — |
| `head` | int | `200` | — |
| `cv_sort` | bool | `True` | — |

Galaxy Tool Shed: `row stats (base R)`

### counts_filter_low_expression

**Filter low-count features** — Drop features that are not expressed well enough for statistical testing.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/counts.tsv` | — |
| `min_counts` | int | `10` | — |
| `min_samples` | int | `1` | — |
| `min_cpm` | number | `0.0` | — |

Galaxy Tool Shed: `DESeq2 row filtering`

### counts_fold_change

**Fold change between two columns** — Per-feature log2 fold change between two samples, columns or sample groups.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/counts.tsv` | — |
| `column_a` | text | `` | — |
| `column_b` | text | `` | — |
| `pseudocount` | number | `1.0` | — |
| `min_abs_log2fc` | number | `0.0` | — |

Galaxy Tool Shed: `DESeq2 results / edgeR`

### counts_log_transform

**Log transform a matrix** — Variance-stabilising log transform with a configurable pseudocount.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/counts.tsv` | — |
| `base` | choice | `2` | 2, 10, e |
| `pseudocount` | int | `1` | — |

Galaxy Tool Shed: `log2 transform (limma)`

### counts_pca_samples

**PCA of samples** — Project libraries onto principal components to reveal replicate structure.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/counts.tsv` | — |
| `components` | int | `2` | — |
| `scale` | choice | `log` | none, zscore, log |

Galaxy Tool Shed: `DESeq2 plotPCA / prcomp`

### counts_quantile_normalize

**Quantile normalisation of a matrix** — Force every sample onto the same distribution before comparing groups.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/counts.tsv` | — |
| `digits` | int | `3` | — |

Galaxy Tool Shed: `limma normalizeQuantiles`

### counts_sample_correlation

**Sample-to-sample correlation matrix** — Heatmap of pairwise correlations between libraries, the batch-effect smoke test.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/counts.tsv` | — |
| `method` | choice | `pearson` | pearson, spearman |
| `log` | bool | `True` | — |

Galaxy Tool Shed: `DESeq2 sampleDistance / corrplot`

### counts_summary

**Summary of a counts matrix** — Library sizes, genes detected, zero counts and the most abundant features.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/counts.tsv` | — |
| `top` | int | `12` | — |

Galaxy Tool Shed: `HTSeq-count / featureCounts output QC`

### counts_tidy_long

**Long format of an expression matrix** — Unroll the matrix into feature / sample / value rows for plotting tools.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/counts.tsv` | — |
| `id_name` | text | `feature` | — |
| `head` | int | `500` | — |

Galaxy Tool Shed: `tidyr pivot_longer`

### counts_tpm

**Normalise to TPM using transcript lengths** — Length-corrected transcripts per million with a FASTA of transcript sequences.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/counts.tsv` | — |
| `lengths` | file | `examples/genes.fa` | — |
| `digits` | int | `2` | — |

Galaxy Tool Shed: `salmon / featureCounts + TPM`

### counts_zscore

**Row z-scores of a matrix** — Centre and scale each row (or column) so patterns are visible in a heatmap.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/counts.tsv` | — |
| `axis` | int | `1` | — |
| `digits` | int | `3` | — |

Galaxy Tool Shed: `pheatmap / row z-score`

### missing_value_imputation

**Impute missing values in a table** — Fill NaNs column-wise (constant, mean/median or k-nearest-neighbour).

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/counts.tsv` | — |
| `method` | choice | `mean` | mean, median, zero, min, knn |
| `k` | int | `5` | — |

Galaxy Tool Shed: `imputation`

### row_math

**Row-wise arithmetic over numeric columns** — Aggregate across each row of a numeric matrix.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/counts.tsv` | — |
| `columns` | text | `all` | — |
| `operation` | choice | `row_sum` | row_sum, row_mean, row_min, row_max, row_product, row_geomean, row_stdev, row_range, row_count_nonzero |

Galaxy Tool Shed: `row_math`

### scale_matrix

**Scale a matrix by rows, columns or globally** — Standardise a numeric matrix along an axis.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/counts.tsv` | — |
| `axis` | choice | `rows` | rows, columns, global |
| `mode` | choice | `z-score` | z-score, min-max, center, unit_variance |

Galaxy Tool Shed: `scale`

## text_manipulation

*Text Manipulation* — in Galaxy group *General Text Tools*.

| Tool | What it does | Upstream Galaxy tool |
|---|---|---|
| [`text_add_row_index`](#text_add_row_index) | Number the rows (1-based, offset configurable). | `index` |
| [`text_change_case`](#text_change_case) | Convert text to upper/lower/title case (Galaxy 'Change Case'). | `change_case` |
| [`text_column_formatter`](#text_column_formatter) | Pretty-print a delimited table with aligned columns. | `table_formatter` |
| [`text_column_summarizer`](#text_column_summarizer) | Whole-table statistics: number of rows/cols, numeric cells, sparsity. | `summary` |
| [`text_compare_two_files`](#text_compare_two_files) | ``cmp``-style diff: identical lines, only-in-A and only-in-B counts. | `compare` |
| [`text_concatenate`](#text_concatenate) | Stack two text datasets vertically. | `concat` |
| [`text_cut_columns`](#text_cut_columns) | Keep the selected columns of a delimited table (POSIX ``cut``). | `text_processing` |
| [`text_diff_two_sorted`](#text_diff_two_sorted) | Lines only in A, only in B and in both (``comm -3``). | `comm` |
| [`text_extract_column_by_name`](#text_extract_column_by_name) | Keep named columns (pandas-style selection). | `select` |
| [`text_extract_with_regex`](#text_extract_with_regex) | ``grep -oP`` style extraction of capture groups. | `extract` |
| [`text_find_most_dispensable_column`](#text_find_most_dispensable_column) | Rank columns by how redundant they are (constant or duplicate-heavy). | `dispensable` |
| [`text_grep`](#text_grep) | ``grep``: keep (or drop with -v) lines matching a regular expression. | `text_processing` |
| [`text_head`](#text_head) | ``head -n``: the first n lines of a text file. | `head` |
| [`text_histogram_of_column`](#text_histogram_of_column) | Bin the values of one column (``datamash histogram``). | `histogram` |
| [`text_insert_lines`](#text_insert_lines) | Insert a fixed line every ``interval`` lines. | `insert_lines` |
| [`text_join_files`](#text_join_files) | Equi-join on one key column (``join``/SQL inner-outer). | `join` |
| [`text_join_with_loose_cutoff`](#text_join_with_loose_cutoff) | Join two tables on their first column allowing up to ``cutoff`` mismatches. | `loose_join` |
| [`text_line_numbers`](#text_line_numbers) | ``nl``-style numbering of lines. | `number_lines` |
| [`text_paste_columns`](#text_paste_columns) | Concatenate two files line by line (``paste``). | `text_processing` |
| [`text_random_lines`](#text_random_lines) | ``shuf``: shuffle lines, or draw a random subset. | `shuffle` |
| [`text_rename_by_pattern`](#text_rename_by_pattern) | Apply a regexp substitution to a single identifier column. | `replace_in_column` |
| [`text_replace`](#text_replace) | ``sed``-style substitution over every line of a text dataset. | `text_processing` |
| [`text_sort`](#text_sort) | ``sort -k``: lexicographic or numeric sort of lines / a column. | `sort` |
| [`text_split_on_column`](#text_split_on_column) | Group lines by the value of one column (like ``split`` per category). | `split` |
| [`text_subtract`](#text_subtract) | Remove every line of ``main`` that also occurs in ``sub``. | `subtract` |
| [`text_tail`](#text_tail) | ``tail -n``: the last n lines, optionally keeping the header. | `tail` |
| [`text_tr`](#text_tr) | ``tr``: translate set1→set2, delete characters, or squeeze repeats. | `text_processing` |
| [`text_uniq`](#text_uniq) | ``uniq -c``: collapse repeated lines and count them. | `uniq` |
| [`text_word_counter`](#text_word_counter) | ``wc``: lines, words and characters, plus the top words. | `word_counter` |
| [`text_wrap`](#text_wrap) | ``fold``-style line wrapping. | `textwrapper` |

### text_add_row_index

**Add a row number column** — Number the rows (1-based, offset configurable).

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/counts.tsv` | — |
| `name` | text | `row` | — |
| `start` | int | `1` | — |

Galaxy Tool Shed: `index`

### text_change_case

**Change case of a text dataset** — Convert text to upper/lower/title case (Galaxy 'Change Case').

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `` | — |
| `mode` | choice | `upper` | upper, lower, title, sentence, swap, inverse |

Galaxy Tool Shed: `change_case`

### text_column_formatter

**Reformat a table to aligned fixed width** — Pretty-print a delimited table with aligned columns.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/counts.tsv` | — |
| `pad` | int | `12` | — |
| `header` | bool | `True` | — |

Galaxy Tool Shed: `table_formatter`

### text_column_summarizer

**Numeric summary of the whole table** — Whole-table statistics: number of rows/cols, numeric cells, sparsity.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/counts.tsv` | — |
| `per_row` | bool | `False` | — |

Galaxy Tool Shed: `summary`

### text_compare_two_files

**Compare two files line by line** — ``cmp``-style diff: identical lines, only-in-A and only-in-B counts.

| parameter | kind | default | options |
|---|---|---|---|
| `first` | file | `` | — |
| `second` | file | `` | — |

Galaxy Tool Shed: `compare`

### text_concatenate

**Concatenate datasets end to end** — Stack two text datasets vertically.

| parameter | kind | default | options |
|---|---|---|---|
| `first` | file | `` | — |
| `second` | file | `` | — |
| `blank_line` | bool | `False` | — |

Galaxy Tool Shed: `concat`

### text_cut_columns

**Cut columns from a table** — Keep the selected columns of a delimited table (POSIX ``cut``).

| parameter | kind | default | options |
|---|---|---|---|
| `columns` | text | `1,3` | — |
| `src` | file | `examples/counts.tsv` | — |

Galaxy Tool Shed: `text_processing`

### text_diff_two_sorted

**Compare two sorted ID lists (comm)** — Lines only in A, only in B and in both (``comm -3``).

| parameter | kind | default | options |
|---|---|---|---|
| `a` | file | `` | — |
| `b` | file | `` | — |

Galaxy Tool Shed: `comm`

### text_extract_column_by_name

**Select columns by name** — Keep named columns (pandas-style selection).

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/counts.tsv` | — |
| `columns` | multi | `['group', 'yield']` | sample, group, yield, biomass |
| `all_columns` | bool | `False` | — |

Galaxy Tool Shed: `select`

### text_extract_with_regex

**Extract fields with a regular expression** — ``grep -oP`` style extraction of capture groups.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `` | — |
| `pattern` | text | `(\w+)\t(\d+)` | — |
| `groups` | bool | `True` | — |

Galaxy Tool Shed: `extract`

### text_find_most_dispensable_column

**Find the most dispensable columns** — Rank columns by how redundant they are (constant or duplicate-heavy).

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/counts.tsv` | — |
| `header` | bool | `True` | — |

Galaxy Tool Shed: `dispensable`

### text_grep

**Filter lines matching a pattern** — ``grep``: keep (or drop with -v) lines matching a regular expression.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `` | — |
| `pattern` | text | `gene` | — |
| `invert` | bool | `False` | — |
| `ignore_case` | bool | `True` | — |
| `whole_line` | bool | `False` | — |
| `max_results` | int | `0` | — |

Galaxy Tool Shed: `text_processing`

### text_head

**Select first N lines** — ``head -n``: the first n lines of a text file.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `` | — |
| `n` | int | `10` | — |

Galaxy Tool Shed: `head`

### text_histogram_of_column

**Histogram of a numeric column** — Bin the values of one column (``datamash histogram``).

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/counts.tsv` | — |
| `column` | int | `2` | — |
| `bins` | int | `20` | — |
| `header` | bool | `True` | — |

Galaxy Tool Shed: `histogram`

### text_insert_lines

**Insert lines before/after every Nth line** — Insert a fixed line every ``interval`` lines.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `` | — |
| `interval` | int | `2` | — |
| `text_line` | text | `## separator` | — |
| `position` | choice | `before` | before, after |

Galaxy Tool Shed: `insert_lines`

### text_join_files

**Join two delimited files on a column** — Equi-join on one key column (``join``/SQL inner-outer).

| parameter | kind | default | options |
|---|---|---|---|
| `left` | file | `examples/counts.tsv` | — |
| `right` | file | `examples/counts.tsv` | — |
| `left_key` | int | `0` | — |
| `right_key` | int | `0` | — |
| `mode` | choice | `inner` | inner, outer, left outer |
| `header` | bool | `True` | — |

Galaxy Tool Shed: `join`

### text_join_with_loose_cutoff

**Relaxed join on a prefix** — Join two tables on their first column allowing up to ``cutoff`` mismatches.

| parameter | kind | default | options |
|---|---|---|---|
| `left` | file | `examples/counts.tsv` | — |
| `right` | file | `examples/counts.tsv` | — |
| `prefix_length` | int | `8` | — |
| `cutoff` | int | `2` | — |
| `ignore_case` | bool | `True` | — |

Galaxy Tool Shed: `loose_join`

### text_line_numbers

**Prepend line numbers** — ``nl``-style numbering of lines.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `` | — |
| `start` | int | `1` | — |
| `sep` | text | `	` | — |

Galaxy Tool Shed: `number_lines`

### text_paste_columns

**Paste tool (join files side by side)** — Concatenate two files line by line (``paste``).

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/counts.tsv` | — |
| `other` | file | `examples/counts.tsv` | — |
| `delimiter` | choice | `tab` | tab, comma, space, newline |
| `whole_line` | bool | `True` | — |

Galaxy Tool Shed: `text_processing`

### text_random_lines

**Shuffle or randomly sample lines** — ``shuf``: shuffle lines, or draw a random subset.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `` | — |
| `n` | int | `0` | — |
| `fraction` | number | `1.0` | — |
| `seed` | int | `1` | — |
| `header_kept` | bool | `False` | — |

Galaxy Tool Shed: `shuffle`

### text_rename_by_pattern

**Rename identifiers with a pattern** — Apply a regexp substitution to a single identifier column.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `` | — |
| `pattern` | text | `gene` | — |
| `replacement` | text | `GENE` | — |
| `column` | int | `1` | — |

Galaxy Tool Shed: `replace_in_column`

### text_replace

**Find and replace (regular expressions)** — ``sed``-style substitution over every line of a text dataset.

| parameter | kind | default | options |
|---|---|---|---|
| `pattern` | text | `` | — |
| `replacement` | text | `` | — |
| `src` | file | `` | — |
| `regex` | bool | `True` | — |
| `case_insensitive` | bool | `False` | — |

Galaxy Tool Shed: `text_processing`

### text_sort

**Sort a text file line by line** — ``sort -k``: lexicographic or numeric sort of lines / a column.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `` | — |
| `column` | int | `-1` | — |
| `numeric` | bool | `False` | — |
| `reverse` | bool | `False` | — |
| `unique` | bool | `False` | — |
| `header` | bool | `False` | — |

Galaxy Tool Shed: `sort`

### text_split_on_column

**Split a table into chunks by column value** — Group lines by the value of one column (like ``split`` per category).

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/counts.tsv` | — |
| `column` | int | `1` | — |
| `max_groups` | int | `20` | — |

Galaxy Tool Shed: `split`

### text_subtract

**Subtract one text file from another** — Remove every line of ``main`` that also occurs in ``sub``.

| parameter | kind | default | options |
|---|---|---|---|
| `main` | file | `` | — |
| `sub` | file | `` | — |
| `whole_line` | bool | `True` | — |

Galaxy Tool Shed: `subtract`

### text_tail

**Select last N lines** — ``tail -n``: the last n lines, optionally keeping the header.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `` | — |
| `n` | int | `10` | — |
| `skip_header` | bool | `False` | — |

Galaxy Tool Shed: `tail`

### text_tr

**Translate or delete characters** — ``tr``: translate set1→set2, delete characters, or squeeze repeats.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `` | — |
| `set1` | text | `ACGT` | — |
| `set2` | text | `TGCA` | — |
| `mode` | choice | `translate` | translate, delete, squeeze, complement |

Galaxy Tool Shed: `text_processing`

### text_uniq

**Report unique/adjacent duplicate lines** — ``uniq -c``: collapse repeated lines and count them.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `` | — |
| `count` | bool | `True` | — |
| `ignore_case` | bool | `False` | — |
| `adjacent_only` | bool | `False` | — |

Galaxy Tool Shed: `uniq`

### text_word_counter

**Word/character/line counter** — ``wc``: lines, words and characters, plus the top words.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `` | — |

Galaxy Tool Shed: `word_counter`

### text_wrap

**Wrap text to fixed width** — ``fold``-style line wrapping.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `` | — |
| `width` | int | `60` | — |
| `indent` | int | `0` | — |
| `newlines` | bool | `False` | — |

Galaxy Tool Shed: `textwrapper`

## filter_and_sort

*Filter and Sort* — in Galaxy group *General Text Tools*.

| Tool | What it does | Upstream Galaxy tool |
|---|---|---|
| [`filter_by_expression`](#filter_by_expression) | Row filter with a pandas ``query`` expression (safe, no eval of code). | `table_filter` |
| [`filter_by_ids`](#filter_by_ids) | Keep rows whose key appears in a list (or does not, with invert). | `filter_from_file` |
| [`filter_first_n`](#filter_first_n) | Keep the first or last n lines. | `head` |
| [`filter_last_column`](#filter_last_column) | Compare the last column against a cut-off. | `filter_last_column` |
| [`filter_missing_values`](#filter_missing_values) | Drop rows/columns whose missing-value fraction exceeds the cut-off. | `filter` |
| [`filter_non_numeric_rows`](#filter_non_numeric_rows) | Drop rows where fewer than ``min_numeric_fraction``% of the fields are numeric. | `filter_rows` |
| [`filter_top_n_by_column`](#filter_top_n_by_column) | Rank rows by a column and keep the best/worst N. | `top_scores` |
| [`filter_with_regexp`](#filter_with_regexp) | ``awk``-like regexp filter on one field. | `filter_with_regexp` |
| [`first_and_last_column_stats`](#first_and_last_column_stats) | Summary stats for the outer columns of a numeric table. | `first_last` |
| [`remove_duplicate_rows`](#remove_duplicate_rows) | Collapse rows that share the same key columns. | `unique` |
| [`sample_rows`](#sample_rows) | Random / systematic / first-n sampling of table rows. | `random_lines` |
| [`shuffle_table_rows`](#shuffle_table_rows) | Deterministic row shuffling with a seed. | `random_shuffle` |
| [`sort_genomic_ranges`](#sort_genomic_ranges) | BED-aware sort: chromosome then start (or interval size). | `sort_bed` |
| [`sort_table_column`](#sort_table_column) | Sort rows by a named column or an index. | `sort` |
| [`text_filter_by_row_count`](#text_filter_by_row_count) | Keep the first fraction (or n) rows of a dataset. | `filter_by` |
| [`text_top_scores`](#text_top_scores) | Add a rank column computed from a numeric column. | `rank` |

### filter_by_expression

**Filter table rows with an expression** — Row filter with a pandas ``query`` expression (safe, no eval of code).

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/counts.tsv` | — |
| `expression` | text | `yield > 7` | — |
| `header` | bool | `True` | — |

Galaxy Tool Shed: `table_filter`

### filter_by_ids

**Filter lines matching identifiers from another dataset** — Keep rows whose key appears in a list (or does not, with invert).

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/counts.tsv` | — |
| `ids` | file | `` | — |
| `column` | int | `1` | — |
| `invert` | bool | `False` | — |

Galaxy Tool Shed: `filter_from_file`

### filter_first_n

**Filter first/last lines of a dataset** — Keep the first or last n lines.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/counts.tsv` | — |
| `which` | choice | `first` | first, last |
| `n` | int | `5` | — |

Galaxy Tool Shed: `head`

### filter_last_column

**Filter rows by the value in the last column** — Compare the last column against a cut-off.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/counts.tsv` | — |
| `threshold` | number | `5.0` | — |
| `op` | choice | `>` | >, >=, <, <=, ==, != |

Galaxy Tool Shed: `filter_last_column`

### filter_missing_values

**Remove rows or columns with missing data** — Drop rows/columns whose missing-value fraction exceeds the cut-off.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/counts.tsv` | — |
| `axis` | choice | `rows` | rows, columns |
| `max_missing_fraction` | number | `0.0` | — |

Galaxy Tool Shed: `filter`

### filter_non_numeric_rows

**Remove rows with non-numeric values** — Drop rows where fewer than ``min_numeric_fraction``% of the fields are numeric.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/counts.tsv` | — |
| `min_numeric_fraction` | int | `100` | — |

Galaxy Tool Shed: `filter_rows`

### filter_top_n_by_column

**Top/bottom N rows by a column** — Rank rows by a column and keep the best/worst N.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/counts.tsv` | — |
| `column` | text | `1` | — |
| `n` | int | `10` | — |
| `which` | choice | `top` | top, bottom |
| `header` | bool | `True` | — |

Galaxy Tool Shed: `top_scores`

### filter_with_regexp

**Filter rows whose key matches a pattern** — ``awk``-like regexp filter on one field.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/counts.tsv` | — |
| `pattern` | text | `gene` | — |
| `field` | choice | `first` | first, last, any |
| `invert` | bool | `False` | — |

Galaxy Tool Shed: `filter_with_regexp`

### first_and_last_column_stats

**Statistics of the first and last columns** — Summary stats for the outer columns of a numeric table.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/counts.tsv` | — |

Galaxy Tool Shed: `first_last`

### remove_duplicate_rows

**Unique rows / detect duplicates** — Collapse rows that share the same key columns.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/counts.tsv` | — |
| `key_columns` | int | `1` | — |
| `header` | bool | `True` | — |
| `keep` | choice | `first` | first, last, none |

Galaxy Tool Shed: `unique`

### sample_rows

**Take a random subset of rows** — Random / systematic / first-n sampling of table rows.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/counts.tsv` | — |
| `n` | int | `5` | — |
| `seed` | int | `1` | — |
| `header` | bool | `True` | — |
| `strategy` | choice | `random` | random, systematic, first, last |

Galaxy Tool Shed: `random_lines`

### shuffle_table_rows

**Randomly shuffle the order of rows** — Deterministic row shuffling with a seed.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/counts.tsv` | — |
| `seed` | int | `42` | — |
| `header` | bool | `True` | — |

Galaxy Tool Shed: `random_shuffle`

### sort_genomic_ranges

**Sort ranges by chromosomal location** — BED-aware sort: chromosome then start (or interval size).

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/regions.bed` | — |
| `order` | choice | `asc` | asc, desc |
| `by_size` | bool | `False` | — |

Galaxy Tool Shed: `sort_bed`

### sort_table_column

**Sort by column** — Sort rows by a named column or an index.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/counts.tsv` | — |
| `column` | text | `1` | — |
| `order` | choice | `ascending` | ascending, descending |
| `numeric` | bool | `True` | — |
| `header` | bool | `True` | — |

Galaxy Tool Shed: `sort`

### text_filter_by_row_count

**Keep the first N rows or a fraction** — Keep the first fraction (or n) rows of a dataset.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/counts.tsv` | — |
| `fraction` | number | `0.5` | — |
| `n` | int | `0` | — |
| `header` | bool | `True` | — |

Galaxy Tool Shed: `filter_by`

### text_top_scores

**Rank rows by a column and output ranks** — Add a rank column computed from a numeric column.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/counts.tsv` | — |
| `column` | text | `2` | — |
| `descending` | bool | `True` | — |
| `tie_average` | bool | `False` | — |

Galaxy Tool Shed: `rank`

## join__subtract_and_group

*Join, Subtraction and Group Operations* — in Galaxy group *General Text Tools*.

| Tool | What it does | Upstream Galaxy tool |
|---|---|---|
| [`aggregate_group`](#aggregate_group) | Collapse a table to one row per key using an aggregate function. | `collapse` |
| [`group_by_column`](#group_by_column) | Group-by summary: size, first members and column means per group. | `group_by` |
| [`jaccard_between_tables`](#jaccard_between_tables) | Similarity of two gene/id lists (``jaccard``). | `jaccard` |
| [`relaxed_cluster_rows`](#relaxed_cluster_rows) | Group identifiers whose names are more similar than a cut-off. | `cluster` |
| [`reverse_rows`](#reverse_rows) | Flip the order of the rows (keeps the header on top). | `reverse_order` |
| [`set_operations_on_columns`](#set_operations_on_columns) | Classic set algebra over the first column of two datasets. | `operate_on_genomic_intervals` |
| [`subtract_tables`](#subtract_tables) | Remove the parts of A covered by B (``bedtools subtract -A``). | `bedtools_subtract` |
| [`transpose_table`](#transpose_table) | Swap the axes of a table (``transpose``). | `transpose` |

### aggregate_group

**Aggregate values within groups (collapse)** — Collapse a table to one row per key using an aggregate function.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/counts.tsv` | — |
| `key_column` | text | `1` | — |
| `value_column` | text | `2` | — |
| `function` | choice | `sum` | sum, mean, median, min, max, count, concat, unique_count, stdev |

Galaxy Tool Shed: `collapse`

### group_by_column

**Group rows by a column value** — Group-by summary: size, first members and column means per group.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/counts.tsv` | — |
| `by` | text | `1` | — |
| `header` | bool | `True` | — |

Galaxy Tool Shed: `group_by`

### jaccard_between_tables

**Jaccard index of two ID sets** — Similarity of two gene/id lists (``jaccard``).

| parameter | kind | default | options |
|---|---|---|---|
| `a` | file | `` | — |
| `b` | file | `` | — |

Galaxy Tool Shed: `jaccard`

### relaxed_cluster_rows

**Cluster rows by similarity of their keys** — Group identifiers whose names are more similar than a cut-off.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `` | — |
| `threshold` | number | `0.75` | — |

Galaxy Tool Shed: `cluster`

### reverse_rows

**Reverse the row order of a dataset** — Flip the order of the rows (keeps the header on top).

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/counts.tsv` | — |
| `header` | bool | `True` | — |

Galaxy Tool Shed: `reverse_order`

### set_operations_on_columns

**Union / intersect / difference of two ID lists** — Classic set algebra over the first column of two datasets.

| parameter | kind | default | options |
|---|---|---|---|
| `a` | file | `` | — |
| `b` | file | `` | — |
| `operation` | choice | `union` | union, intersect, difference, complement, symmetric_difference |

Galaxy Tool Shed: `operate_on_genomic_intervals`

### subtract_tables

**Subtract intervals or rows between datasets** — Remove the parts of A covered by B (``bedtools subtract -A``).

| parameter | kind | default | options |
|---|---|---|---|
| `a` | file | `examples/regions.bed` | — |
| `b` | file | `examples/regions.bed` | — |

Galaxy Tool Shed: `bedtools_subtract`

### transpose_table

**Transpose rows and columns** — Swap the axes of a table (``transpose``).

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/counts.tsv` | — |
| `header` | bool | `True` | — |

Galaxy Tool Shed: `transpose`

## datamash

*Datamash: Aggregation and Statistics* — in Galaxy group *General Text Tools*.

| Tool | What it does | Upstream Galaxy tool |
|---|---|---|
| [`datamash_bin_column`](#datamash_bin_column) | Histogram-style binned aggregation (``datamash --histogram``). | `datamash` |
| [`datamash_column_stats`](#datamash_column_stats) | ``datamash geometric-mean 2 mean 3 ...``: aggregate chosen columns. | `datamash` |
| [`datamash_cross_tab`](#datamash_cross_tab) | ``datamash crosstab``-style contingency table (counts or aggregated values). | `datamash` |
| [`datamash_cumulative`](#datamash_cumulative) | Row-wise cumulative transforms of a numeric column. | `datamash` |
| [`datamash_grouped`](#datamash_grouped) | Group rows by a key column then aggregate a numeric column per group. | `datamash` |
| [`datamash_matrix_summary`](#datamash_matrix_summary) | Per-row / per-column aggregate of a numeric matrix (Galaxy 'Matrix Summary'). | `matrix_summary` |
| [`datamash_normalize_column`](#datamash_normalize_column) | Common rescalings of one column (datamash --normalize). | `normalize` |
| [`datamash_outliers`](#datamash_outliers) | Flag outliers with IQR, MAD-based modified z or plain z-scores. | `outlier_detection` |
| [`datamash_percentiles`](#datamash_percentiles) | Report arbitrary percentiles of a numeric column. | `datamash` |
| [`datamash_summary_all`](#datamash_summary_all) | count/mean/median/min/max/q1/q3/stdev/IQR/MAD per numeric column. | `summary_statistics` |
| [`datamash_unique_values`](#datamash_unique_values) | ``datamash unique/count``: frequency table of a categorical column. | `datamash` |

### datamash_bin_column

**Bin a column and summarise each bin** — Histogram-style binned aggregation (``datamash --histogram``).

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/counts.tsv` | — |
| `column` | text | `1` | — |
| `bin_size` | number | `10.0` | — |
| `aggregate` | choice | `count` | count, mean, sum, max, min |

Galaxy Tool Shed: `datamash`

### datamash_column_stats

**Datamash-style column statistics** — ``datamash geometric-mean 2 mean 3 ...``: aggregate chosen columns.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/counts.tsv` | — |
| `functions` | multi | `['mean', 'median', 'min', 'max']` | mean, median, geomean, harmonic_mean, min, max, sum, count, countunique, unique, elapse, collapse, first_value, last_value, mad, pstdev, samplestdev, variance, pvariance, percentile, interquartile_range, range, absolute_sum, cumulative_sum |
| `columns` | text | `all` | — |
| `header` | bool | `True` | — |

Galaxy Tool Shed: `datamash`

### datamash_cross_tab

**Cross-tabulation of two columns** — ``datamash crosstab``-style contingency table (counts or aggregated values).

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/counts.tsv` | — |
| `row_column` | text | `1` | — |
| `col_column` | text | `2` | — |
| `value_column` | text | `` | — |
| `header` | bool | `True` | — |

Galaxy Tool Shed: `datamash`

### datamash_cumulative

**Cumulative sums and running statistics** — Row-wise cumulative transforms of a numeric column.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/counts.tsv` | — |
| `column` | text | `1` | — |
| `function` | choice | `cumsum` | cumsum, running_mean, ewma, diff, pct_change, max_so_far, min_so_far |
| `header` | bool | `True` | — |

Galaxy Tool Shed: `datamash`

### datamash_grouped

**Datamash groupby + aggregate** — Group rows by a key column then aggregate a numeric column per group.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/counts.tsv` | — |
| `group_column` | text | `1` | — |
| `value_column` | text | `2` | — |
| `function` | choice | `mean` | mean, median, geomean, harmonic_mean, min, max, sum, count, countunique, unique, elapse, collapse, first_value, last_value, mad, pstdev, samplestdev, variance, pvariance, percentile, interquartile_range, range, absolute_sum, cumulative_sum |
| `header` | bool | `True` | — |

Galaxy Tool Shed: `datamash`

### datamash_matrix_summary

**Row and column summary of a matrix** — Per-row / per-column aggregate of a numeric matrix (Galaxy 'Matrix Summary').

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/counts.tsv` | — |
| `row_sums` | bool | `True` | — |
| `col_sums` | bool | `True` | — |
| `stat` | choice | `sum` | sum, mean, median, min, max |

Galaxy Tool Shed: `matrix_summary`

### datamash_normalize_column

**Normalise a numeric column** — Common rescalings of one column (datamash --normalize).

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/counts.tsv` | — |
| `column` | text | `1` | — |
| `mode` | choice | `z-score` | z-score, min-max, relative abundance (percent), rank, decimal, log2, sum-to-one, robust (median/MAD) |
| `header` | bool | `True` | — |

Galaxy Tool Shed: `normalize`

### datamash_outliers

**Detect outliers in a column** — Flag outliers with IQR, MAD-based modified z or plain z-scores.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/counts.tsv` | — |
| `column` | text | `1` | — |
| `method` | choice | `iqr` | iqr, modified_z, zscore, grubbs_like |
| `cutoff` | number | `1.5` | — |

Galaxy Tool Shed: `outlier_detection`

### datamash_percentiles

**Percentiles of one column** — Report arbitrary percentiles of a numeric column.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/counts.tsv` | — |
| `column` | text | `1` | — |
| `list_of_percentiles` | text | `5,25,50,75,95` | — |
| `header` | bool | `True` | — |

Galaxy Tool Shed: `datamash`

### datamash_summary_all

**Summary statistics for every numeric column** — count/mean/median/min/max/q1/q3/stdev/IQR/MAD per numeric column.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/counts.tsv` | — |
| `header` | bool | `True` | — |

Galaxy Tool Shed: `summary_statistics`

### datamash_unique_values

**Unique values with counts and proportions** — ``datamash unique/count``: frequency table of a categorical column.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/counts.tsv` | — |
| `column` | text | `1` | — |
| `top` | int | `0` | — |
| `percent` | bool | `True` | — |

Galaxy Tool Shed: `datamash`

## fasta_fastq

*FASTA/FASTQ* — in Galaxy group *Genomic File Manipulation*.

| Tool | What it does | Upstream Galaxy tool |
|---|---|---|
| [`fasta_barcode_split`](#fasta_barcode_split) | Assign reads to samples by a 5'/3' barcode prefix. | `bctools_barcode` |
| [`fasta_chop`](#fasta_chop) | Split each sequence into equal-size pieces with numbered ids. | `fasta_split_on_cases` |
| [`fasta_clean_n`](#fasta_clean_n) | Trim leading/trailing runs of N (optionally mask internal runs). | `clean_n` |
| [`fasta_complement`](#fasta_complement) | Return the reverse complement (or one half of it). | `fastarevseq` |
| [`fasta_count_nucleotides`](#fasta_count_nucleotides) | A/C/G/T/N counts and percentages (Galaxy 'Count DNA nucleotides'). | `count_gatb` |
| [`fasta_cpg_report`](#fasta_cpg_report) | Windows with high CpG content and observed/expected CpG ratio. | `newcpgreport` |
| [`fasta_dedupe`](#fasta_dedupe) | Collapse records that share the same sequence (or id). | `rdist` |
| [`fasta_dna_to_protein`](#fasta_dna_to_protein) | Translate a nucleotide FASTA in one reading frame. | `mgrna2nuc` |
| [`fasta_encode_aa_props`](#fasta_encode_aa_props) | Numeric per-residue property track (Galaxy 'Encode Aminoacid Properties'). | `aa_prop` |
| [`fasta_extract_by_ids`](#fasta_extract_by_ids) | Extract (or drop) records whose id matches a list. | `iuc/seqkit` |
| [`fasta_fasta_stats_of_records`](#fasta_fasta_stats_of_records) | Attach a uniform quality string to FASTA (or drop qualities). | `fasta_to_fastq` |
| [`fasta_filter_by_length`](#fasta_filter_by_length) | Keep sequences whose length falls inside a window. | `iuc/seqtk` |
| [`fasta_filter_on_description`](#fasta_filter_on_description) | Keep records whose defline matches a regexp. | `fasta-rc` |
| [`fasta_gc_content`](#fasta_gc_content) | Length, GC%, AT-skew and GC-skew of every record. | `seqtk` |
| [`fasta_get_hsp`](#fasta_get_hsp) | Collapse aligned sequences into a consensus string. | `new_consensus` |
| [`fasta_get_ids`](#fasta_get_ids) | List the names of all records in a FASTA file. | `fasta_to_tabulate` |
| [`fasta_get_length`](#fasta_get_length) | Length of every sequence, plus file totals. | `fasta-get-metadata` |
| [`fasta_get_orfs`](#fasta_get_orfs) | Find open reading frames (EMBOSS ``getorf``). | `getorf` |
| [`fasta_head`](#fasta_head) | Extract the first n records. | `fasta_first_sequences` |
| [`fasta_kmer_count`](#fasta_kmer_count) | Count k-mers over the file (oligonucleotide frequencies). | `seqtk_seq` |
| [`fasta_longest_orf`](#fasta_longest_orf) | Print the longest translation found in the file. | `longest_orf` |
| [`fasta_mash_compare`](#fasta_mash_compare) | Sketch each sequence and estimate pairwise Mash distances. | `tools-mash` |
| [`fasta_md5`](#fasta_md5) | Hash of every sequence (identical data → identical hash). | `seqtk_md5` |
| [`fasta_nucleotide_frequencies`](#fasta_nucleotide_frequencies) | Frequency of every base in the concatenated alignment. | `fasta-rc_cycle` |
| [`fasta_nucleotide_percentage`](#fasta_nucleotide_percentage) | Sliding-window GC and nucleotide composition. | `gc_percent` |
| [`fasta_profile_vector`](#fasta_profile_vector) | Composition vector for sequence similarity work. | `profile_vector` |
| [`fasta_random_subset`](#fasta_random_subset) | Draw a random sample of records. | `random-choose` |
| [`fasta_remove_gaps`](#fasta_remove_gaps) | Strip gap and filler characters from every sequence. | `ungap` |
| [`fasta_replace_ids`](#fasta_replace_ids) | Rename records with a prefix/suffix and optional defline cleanup. | `fasta_str_tie` |
| [`fasta_reverse_sequences`](#fasta_reverse_sequences) | Reverse the record order and/or the residues. | `fasta_rc` |
| [`fasta_search_loop`](#fasta_search_loop) | Search a motif or regexp through the file and report hit positions. | `fasta_search_loop` |
| [`fasta_sort`](#fasta_sort) | Order records by name, length or GC content. | `fasta_order` |
| [`fasta_split_by_size`](#fasta_split_by_size) | Chunk a multi-FASTA into groups of k sequences. | `fastasplit` |
| [`fasta_split_on_cases`](#fasta_split_on_cases) | Distribute records into N roughly equal partitions. | `fasta_split` |
| [`fasta_statistics`](#fasta_statistics) | min/max/mean/median length, N50, GC and alphabet checks. | `seqtk comp` |
| [`fasta_string_tie`](#fasta_string_tie) | Tie records end to end into one sequence with a spacer. | `string_tie` |
| [`fasta_tail`](#fasta_tail) | Extract the last n records. | `fasta_last_sequences` |
| [`fasta_to_table`](#fasta_to_table) | Convert records into a tab-delimited summary table. | `fasta_to_tabulate` |
| [`fasta_translate`](#fasta_translate) | Report the translation(s) and longest ORF per sequence. | `translate_tool` |
| [`fasta_trim`](#fasta_trim) | Remove bases from the ends or keep a fixed window. | `trim_sequences` |
| [`fasta_uppercase`](#fasta_uppercase) | Upper/lower-case the residues (soft-masking convention). | `Change_case` |
| [`fasta_wrap`](#fasta_wrap) | Re-wrap sequence lines to a given width. | `fold` |
| [`fastq_filter`](#fastq_filter) | Length / quality / N-content / substring filtering of reads. | `iuc/fastx_toolkit` |
| [`fastq_join_paired`](#fastq_join_paired) | Merge R1/R2 using the exact-overlap algorithm (``fastq_join``/PEAR-like). | `fastq_join` |
| [`fastq_reverse_complement`](#fastq_reverse_complement) | Reverse-complement both sequence and quality strings. | `fastx_revcomp` |
| [`fastq_sort`](#fastq_sort) | Deterministic ordering of reads. | `fastq_sort` |
| [`fastq_split_pairs`](#fastq_split_pairs) | Interleaved → separate R1/R2 files (reported as a table). | `fastx` |
| [`fastq_to_fasta`](#fastq_to_fasta) | Drop quality lines and keep the sequence. | `iuc/fastx_toolkit` |
| [`fastq_trim_length`](#fastq_trim_length) | Fixed truncation and end trimming of FASTQ reads. | `fastx_trimmer` |

### fasta_barcode_split

**Barcode Split (demultiplex by inline barcodes)** — Assign reads to samples by a 5'/3' barcode prefix.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/reads_single.fastq` | — |
| `barcodes` | code | `ACGT
TGCA` | — |
| `where` | choice | `start` | start, end |
| `mismatches` | int | `0` | — |

Galaxy Tool Shed: `bctools_barcode`

### fasta_chop

**Fasta Chopper (split into fixed-size chunks)** — Split each sequence into equal-size pieces with numbered ids.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/genome.fa` | — |
| `chunk_size` | int | `200` | — |
| `discard_partial` | bool | `False` | — |

Galaxy Tool Shed: `fasta_split_on_cases`

### fasta_clean_n

**Clean Ns from sequence ends** — Trim leading/trailing runs of N (optionally mask internal runs).

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/genome.fa` | — |
| `trim_internal` | bool | `False` | — |
| `max_internal_n` | int | `5` | — |

Galaxy Tool Shed: `clean_n`

### fasta_complement

**Reverse complement** — Return the reverse complement (or one half of it).

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/genome.fa` | — |
| `mode` | choice | `revcomp` | revcomp, complement, reverse |

Galaxy Tool Shed: `fastarevseq`

### fasta_count_nucleotides

**Count DNA nucleotides** — A/C/G/T/N counts and percentages (Galaxy 'Count DNA nucleotides').

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/genome.fa` | — |
| `per_sequence` | bool | `True` | — |

Galaxy Tool Shed: `count_gatb`

### fasta_cpg_report

**CpG report (newcpgreport)** — Windows with high CpG content and observed/expected CpG ratio.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/genome.fa` | — |
| `window` | int | `200` | — |
| `step` | int | `100` | — |
| `threshold` | number | `0.6` | — |

Galaxy Tool Shed: `newcpgreport`

### fasta_dedupe

**Remove duplicate sequences** — Collapse records that share the same sequence (or id).

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/genome.fa` | — |
| `key` | choice | `sequence` | sequence, id, sequence_ignoring_case, canonical |
| `report` | bool | `True` | — |

Galaxy Tool Shed: `rdist`

### fasta_dna_to_protein

**DNA to protein** — Translate a nucleotide FASTA in one reading frame.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | code | `examples/genes.fa` | — |
| `frame` | int | `0` | — |
| `table` | choice | `Standard` | Standard, Vertebrate Mitochondrial, Yeast Mitochondrial, Bacterial |
| `to_stop` | bool | `True` | — |

Galaxy Tool Shed: `mgrna2nuc`

### fasta_encode_aa_props

**Encode amino-acid properties** — Numeric per-residue property track (Galaxy 'Encode Aminoacid Properties').

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/proteins.faa` | — |
| `property` | choice | `hydrophobicity` | hydrophobicity, charge, volume, flexibility, polarity |
| `window` | int | `5` | — |

Galaxy Tool Shed: `aa_prop`

### fasta_extract_by_ids

**Fetch sequences by identifiers** — Extract (or drop) records whose id matches a list.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/genome.fa` | — |
| `identifiers` | code | `` | — |
| `invert` | bool | `False` | — |
| `allow_prefix` | bool | `True` | — |

Galaxy Tool Shed: `iuc/seqkit`

### fasta_fasta_stats_of_records

**Enolpha/Phred quality conversion of FASTA** — Attach a uniform quality string to FASTA (or drop qualities).

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/genome.fa` | — |
| `direction` | choice | `phred_to_fasta` | phred_to_fasta, fasta_to_phred |
| `qual` | text | `!!!!!` | — |

Galaxy Tool Shed: `fasta_to_fastq`

### fasta_filter_by_length

**Filter FASTA sequences by length** — Keep sequences whose length falls inside a window.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/genome.fa` | — |
| `min_len` | int | `100` | — |
| `max_len` | int | `0` | — |
| `complement` | bool | `False` | — |

Galaxy Tool Shed: `iuc/seqtk`

### fasta_filter_on_description

**Filter FASTA by description text** — Keep records whose defline matches a regexp.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/genome.fa` | — |
| `pattern` | text | `gene` | — |
| `invert` | bool | `False` | — |

Galaxy Tool Shed: `fasta-rc`

### fasta_gc_content

**GC content per sequence** — Length, GC%, AT-skew and GC-skew of every record.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/genome.fa` | — |
| `precision` | int | `4` | — |

Galaxy Tool Shed: `seqtk`

### fasta_get_hsp

**Consensus sequence of a multiple alignment (FASTA)** — Collapse aligned sequences into a consensus string.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/genome.fa` | — |
| `mode` | choice | `majority` | majority, first, ambiguous |

Galaxy Tool Shed: `new_consensus`

### fasta_get_ids

**Extract sequence identifiers** — List the names of all records in a FASTA file.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/genome.fa` | — |
| `field` | choice | `id` | id, description, id_and_description, index |
| `unique` | bool | `False` | — |

Galaxy Tool Shed: `fasta_to_tabulate`

### fasta_get_length

**Get lengths** — Length of every sequence, plus file totals.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/genome.fa` | — |
| `cumulative` | bool | `False` | — |

Galaxy Tool Shed: `fasta-get-metadata`

### fasta_get_orfs

**Get ORFs** — Find open reading frames (EMBOSS ``getorf``).

| parameter | kind | default | options |
|---|---|---|---|
| `src` | code | `examples/genes.fa` | — |
| `min_length` | int | `90` | — |
| `start_codons` | text | `ATG` | — |
| `both_strands` | bool | `True` | — |
| `table` | bool | `False` | — |

Galaxy Tool Shed: `getorf`

### fasta_head

**First N sequences** — Extract the first n records.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/genome.fa` | — |
| `n` | int | `5` | — |

Galaxy Tool Shed: `fasta_first_sequences`

### fasta_kmer_count

**Nucleotide frequencies (k-mers)** — Count k-mers over the file (oligonucleotide frequencies).

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/genome.fa` | — |
| `k` | int | `4` | — |
| `canonical` | bool | `True` | — |
| `top` | int | `20` | — |
| `normalize` | bool | `True` | — |

Galaxy Tool Shed: `seqtk_seq`

### fasta_longest_orf

**Longest ORF as protein** — Print the longest translation found in the file.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | code | `examples/genes.fa` | — |

Galaxy Tool Shed: `longest_orf`

### fasta_mash_compare

**MinHash sketch and distance matrix** — Sketch each sequence and estimate pairwise Mash distances.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/genome.fa` | — |
| `k` | int | `21` | — |
| `num` | int | `100` | — |
| `distance_matrix` | bool | `True` | — |

Galaxy Tool Shed: `tools-mash`

### fasta_md5

**Sequence MD5 checksums** — Hash of every sequence (identical data → identical hash).

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/genome.fa` | — |
| `algorithm` | choice | `md5` | md5, sha1 |

Galaxy Tool Shed: `seqtk_md5`

### fasta_nucleotide_frequencies

**Nucleotide frequencies** — Frequency of every base in the concatenated alignment.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/genome.fa` | — |

Galaxy Tool Shed: `fasta-rc_cycle`

### fasta_nucleotide_percentage

**Nucleotide percentage plot data** — Sliding-window GC and nucleotide composition.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/genome.fa` | — |
| `window` | int | `200` | — |
| `step` | int | `100` | — |

Galaxy Tool Shed: `gc_percent`

### fasta_profile_vector

**Profile vector (pseudo-monomer composition)** — Composition vector for sequence similarity work.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/genome.fa` | — |
| `k` | int | `4` | — |
| `revcomp_aware` | bool | `True` | — |

Galaxy Tool Shed: `profile_vector`

### fasta_random_subset

**Randomly choose N sequences** — Draw a random sample of records.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/genome.fa` | — |
| `n` | int | `3` | — |
| `seed` | int | `1` | — |
| `fraction` | number | `0.0` | — |
| `keep_order` | bool | `False` | — |

Galaxy Tool Shed: `random-choose`

### fasta_remove_gaps

**Remove gaps, dashes and unknown characters** — Strip gap and filler characters from every sequence.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/genome.fa` | — |
| `chars` | multi | `['-', '.']` | -, .,  , N, *, ?, X |
| `keep_case` | bool | `True` | — |

Galaxy Tool Shed: `ungap`

### fasta_replace_ids

**Add or replace sequence descriptions** — Rename records with a prefix/suffix and optional defline cleanup.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/genome.fa` | — |
| `prefix` | text | `` | — |
| `suffix` | text | `` | — |
| `strip_description` | bool | `True` | — |

Galaxy Tool Shed: `fasta_str_tie`

### fasta_reverse_sequences

**Reverse sequence order and orientation** — Reverse the record order and/or the residues.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/genome.fa` | — |
| `reverse_records` | bool | `True` | — |
| `reverse_residues` | bool | `False` | — |

Galaxy Tool Shed: `fasta_rc`

### fasta_search_loop

**Find sequences matching a pattern** — Search a motif or regexp through the file and report hit positions.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/genome.fa` | — |
| `pattern` | text | `GATC` | — |
| `search_in` | choice | `sequence` | sequence, description, id |
| `reverse_complement` | bool | `True` | — |
| `max_hits` | int | `50` | — |

Galaxy Tool Shed: `fasta_search_loop`

### fasta_sort

**Sort FASTA records** — Order records by name, length or GC content.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/genome.fa` | — |
| `by` | choice | `id` | id, length, gc, alphabetical |
| `reverse` | bool | `False` | — |

Galaxy Tool Shed: `fasta_order`

### fasta_split_by_size

**Split FASTA into files of a fixed size** — Chunk a multi-FASTA into groups of k sequences.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/genome.fa` | — |
| `sequences_per_file` | int | `3` | — |

Galaxy Tool Shed: `fastasplit`

### fasta_split_on_cases

**Split FASTA on a column of a table** — Distribute records into N roughly equal partitions.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/genome.fa` | — |
| `n_parts` | int | `2` | — |

Galaxy Tool Shed: `fasta_split`

### fasta_statistics

**Sequence statistics of a FASTA file** — min/max/mean/median length, N50, GC and alphabet checks.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/genome.fa` | — |

Galaxy Tool Shed: `seqtk comp`

### fasta_string_tie

**Concatenate sequences (String Tying)** — Tie records end to end into one sequence with a spacer.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/genome.fa` | — |
| `revcomp_second` | bool | `False` | — |
| `gap_nt` | int | `0` | — |
| `gap_char` | text | `N` | — |

Galaxy Tool Shed: `string_tie`

### fasta_tail

**Last N sequences** — Extract the last n records.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/genome.fa` | — |
| `n` | int | `5` | — |

Galaxy Tool Shed: `fasta_last_sequences`

### fasta_to_table

**FASTA to tabular** — Convert records into a tab-delimited summary table.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/genome.fa` | — |
| `fields` | multi | `['id', 'length', 'gc']` | id, description, length, sequence, gc, sha1, first_base, last_base |
| `truncate_sequence` | int | `0` | — |

Galaxy Tool Shed: `fasta_to_tabulate`

### fasta_translate

**Translate (all six frames)** — Report the translation(s) and longest ORF per sequence.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | code | `examples/genes.fa` | — |
| `mode` | choice | `best_orf` | best_orf, all_frames |
| `to_stop` | bool | `True` | — |

Galaxy Tool Shed: `translate_tool`

### fasta_trim

**Trim sequences (fixed or by quality)** — Remove bases from the ends or keep a fixed window.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/genome.fa` | — |
| `trim_5p` | int | `0` | — |
| `trim_3p` | int | `0` | — |
| `keep_from` | int | `0` | — |
| `keep_to` | int | `0` | — |

Galaxy Tool Shed: `trim_sequences`

### fasta_uppercase

**Convert sequence case** — Upper/lower-case the residues (soft-masking convention).

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/genome.fa` | — |
| `case` | choice | `upper` | upper, lower |

Galaxy Tool Shed: `Change_case`

### fasta_wrap

**Change line wrapping** — Re-wrap sequence lines to a given width.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/genome.fa` | — |
| `width` | int | `70` | — |
| `single_line` | bool | `False` | — |

Galaxy Tool Shed: `fold`

### fastq_filter

**Filter FASTQ reads** — Length / quality / N-content / substring filtering of reads.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/reads_single.fastq` | — |
| `min_length` | int | `50` | — |
| `max_length` | int | `0` | — |
| `min_mean_qual` | number | `0.0` | — |
| `max_n_fraction` | number | `1.0` | — |
| `contains` | text | `` | — |
| `remove_gaps` | bool | `False` | — |

Galaxy Tool Shed: `iuc/fastx_toolkit`

### fastq_join_paired

**Join paired reads (overlap merging)** — Merge R1/R2 using the exact-overlap algorithm (``fastq_join``/PEAR-like).

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/reads_1.fastq` | — |
| `src2` | file | `examples/reads_2.fastq` | — |
| `min_overlap` | int | `11` | — |
| `max_mismatches` | int | `2` | — |

Galaxy Tool Shed: `fastq_join`

### fastq_reverse_complement

**Reverse complement reads** — Reverse-complement both sequence and quality strings.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/reads_single.fastq` | — |

Galaxy Tool Shed: `fastx_revcomp`

### fastq_sort

**Sort FASTQ by read name or quality** — Deterministic ordering of reads.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/reads_single.fastq` | — |
| `by` | choice | `name` | name, quality, length, sequence |
| `reverse` | bool | `False` | — |

Galaxy Tool Shed: `fastq_sort`

### fastq_split_pairs

**Split paired-end FASTQ into two files** — Interleaved → separate R1/R2 files (reported as a table).

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/reads_single.fastq` | — |
| `suffix_1` | text | `_1` | — |
| `suffix_2` | text | `_2` | — |
| `by_name` | bool | `True` | — |

Galaxy Tool Shed: `fastx`

### fastq_to_fasta

**FASTQ to FASTA** — Drop quality lines and keep the sequence.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/reads_single.fastq` | — |
| `keep_first_word` | bool | `True` | — |

Galaxy Tool Shed: `iuc/fastx_toolkit`

### fastq_trim_length

**Truncate reads** — Fixed truncation and end trimming of FASTQ reads.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/reads_single.fastq` | — |
| `length` | int | `0` | — |
| `left_trim` | int | `0` | — |
| `right_trim` | int | `0` | — |
| `pad_short_reads` | bool | `False` | — |

Galaxy Tool Shed: `fastx_trimmer`

## fastq_quality_control

*FASTQ Quality Control* — in Galaxy group *Genomic File Manipulation*.

| Tool | What it does | Upstream Galaxy tool |
|---|---|---|
| [`fastq_bbsplit_screen`](#fastq_bbsplit_screen) | Assign each read to the best reference by k-mer seeding. | `iuc/bbsplit` |
| [`fastq_collapser`](#fastq_collapser) | Unique sequences with abundance and mean quality (``fastx_collapser``). | `iuc/fastx_toolkit` |
| [`fastq_cutadapt`](#fastq_cutadapt) | Find the best 3' adapter match and clip everything from there. | `iuc/cutadapt` |
| [`fastq_duplication_levels`](#fastq_duplication_levels) | Percentage of reads seen 1x, 2x, … >N times (PCR duplicates). | `iuc/fastqc` |
| [`fastq_encoding`](#fastq_encoding) | Guess Phred+33/+64/Solexa+66 from ASCII ranges. | `fastq-encoding-detect` |
| [`fastq_error_rate`](#fastq_error_rate) | Convert Phred scores into an expected substitution error rate. | `iuc/fastp` |
| [`fastq_fastp`](#fastq_fastp) | Adapter trim + quality filter + duplication stats in one pass. | `iuc/fastp` |
| [`fastq_fastx_toolkit_qc`](#fastq_fastx_toolkit_qc) | ``fastx_quality_stats`` per-cycle statistics for one metric. | `iuc/fastx_toolkit` |
| [`fastq_gc_distribution`](#fastq_gc_distribution) | GC% histogram of reads; a bimodal plot means contamination. | `fastqgc` |
| [`fastq_kmer_content`](#fastq_kmer_content) | K-mer frequencies at the start of reads (adapter signature). | `iuc/fastqc` |
| [`fastq_mott_trimmer`](#fastq_mott_trimmer) | Sliding-window quality trimming as used by prinseq/dynamictrim. | `dynamictrim` |
| [`fastq_n_stats`](#fastq_n_stats) | Number and percentage of Ns at each cycle (or a single summary). | `fastx_n_statistics` |
| [`fastq_normalise`](#fastq_normalise) | Discard reads whose median k-mer coverage exceeds a target (BBNorm-like). | `iuc/bbnorm` |
| [`fastq_overrepresented`](#fastq_overrepresented) | Sequences appearing more often than a threshold. | `iuc/fastqc` |
| [`fastq_per_base_quality`](#fastq_per_base_quality) | Mean/median/quartiles of the quality score at every read position. | `iuc/fastqc` |
| [`fastq_per_sequence_quality`](#fastq_per_sequence_quality) | Histogram of the mean quality of each read. | `iuc/fastqc` |
| [`fastq_polyx_trim`](#fastq_polyx_trim) | Remove homopolymer A/T tails (RNA-seq artefacts). | `iuc/trim_galore` |
| [`fastq_prinseq`](#fastq_prinseq) | Prinseq++-style combined quality/complexity filters. | `du_novo/prinseq` |
| [`fastq_qc_duplicate_report`](#fastq_qc_duplicate_report) | Library complexity: unique vs duplicated read pairs. | `iuc/fastp` |
| [`fastq_quality_hist`](#fastq_quality_hist) | Counts of each Phred score in the whole file. | `fastx_quality_statistics` |
| [`fastq_randomize`](#fastq_randomize) | Shuffle read order to remove positional bias. | `iuc/bbtools` |
| [`fastq_read_length_distribution`](#fastq_read_length_distribution) | Histogram of read lengths (with optional binning). | `fastq_len` |
| [`fastq_rrbs_trim`](#fastq_rrbs_trim) | Cut at the MspI/CTAG site (``trim_galore --rrbs``). | `iuc/trim_galore` |
| [`fastq_sliding_trim`](#fastq_sliding_trim) | Trim bases below a Phred cut-off from both ends, drop short reads. | `iuc/trim_galore` |
| [`fastq_sortbylength`](#fastq_sortbylength) | Deterministic read ordering (SortByLength/SortBySize). | `iuc/fastx_toolkit` |
| [`fastq_subsample`](#fastq_subsample) | Draw a random subset of reads (mates kept in sync). | `iuc/bbtools` |
| [`fastq_tile_stats`](#fastq_tile_stats) | Aggregate quality by the flowcell/lane field of the read name when present. | `fastq_stats` |
| [`fastq_trim_by_length`](#fastq_trim_by_length) | Keep or drop reads outside a length window. | `fastqlenfilter` |
| [`fastq_trim_ns`](#fastq_trim_ns) | Remove N bases from 5' and 3' ends. | `du_novo/prinseq` |
| [`fastqc_summary`](#fastqc_summary) | Per-base quality, GC, N content, duplication and adapter flags. | `iuc/fastqc` |

### fastq_bbsplit_screen

**Read screening against reference sequences** — Assign each read to the best reference by k-mer seeding.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/reads_single.fastq` | — |
| `references` | file | `examples/genome.fa` | — |
| `k` | int | `21` | — |
| `trim` | bool | `False` | — |
| `mismatches` | int | `3` | — |

Galaxy Tool Shed: `iuc/bbsplit`

### fastq_collapser

**Collapse duplicate sequences** — Unique sequences with abundance and mean quality (``fastx_collapser``).

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/reads_single.fastq` | — |
| `min_count` | int | `1` | — |

Galaxy Tool Shed: `iuc/fastx_toolkit`

### fastq_cutadapt

**Cutadapt-style adapter trimming** — Find the best 3' adapter match and clip everything from there.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/reads_single.fastq` | — |
| `adapter` | text | `AGATCGGAAGAGC` | — |
| `min_overlap` | int | `3` | — |
| `errors` | number | `0.1` | — |
| `trim_nextera` | bool | `False` | — |
| `minimum_length` | int | `15` | — |

Galaxy Tool Shed: `iuc/cutadapt`

### fastq_duplication_levels

**Sequence duplication levels** — Percentage of reads seen 1x, 2x, … >N times (PCR duplicates).

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/reads_single.fastq` | — |
| `max_level` | int | `10` | — |

Galaxy Tool Shed: `iuc/fastqc`

### fastq_encoding

**Detect the quality encoding** — Guess Phred+33/+64/Solexa+66 from ASCII ranges.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/reads_single.fastq` | — |

Galaxy Tool Shed: `fastq-encoding-detect`

### fastq_error_rate

**Estimated error rate from qualities** — Convert Phred scores into an expected substitution error rate.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/reads_single.fastq` | — |

Galaxy Tool Shed: `iuc/fastp`

### fastq_fastp

**fastp-like all-in-one QC** — Adapter trim + quality filter + duplication stats in one pass.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/reads_single.fastq` | — |
| `src2` | file | `` | — |
| `qualifier` | int | `15` | — |
| `unqualified_percent_limit` | int | `40` | — |
| `length_required` | int | `15` | — |
| `adapter_sequence` | text | `` | — |
| `dedup` | bool | `False` | — |
| `cut_window_size` | int | `4` | — |
| `cut_mean_quality` | int | `20` | — |

Galaxy Tool Shed: `iuc/fastp`

### fastq_fastx_toolkit_qc

**FASTX-Toolkit per-cycle statistics** — ``fastx_quality_stats`` per-cycle statistics for one metric.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/reads_single.fastq` | — |
| `stat` | choice | `quality` | quality, gc, length, maxquality, minquality, ap, cp |

Galaxy Tool Shed: `iuc/fastx_toolkit`

### fastq_gc_distribution

**Per read GC content distribution** — GC% histogram of reads; a bimodal plot means contamination.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/reads_single.fastq` | — |
| `bins` | int | `20` | — |

Galaxy Tool Shed: `fastqgc`

### fastq_kmer_content

**Overrepresented k-mers** — K-mer frequencies at the start of reads (adapter signature).

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/reads_single.fastq` | — |
| `k` | int | `5` | — |
| `top` | int | `10` | — |

Galaxy Tool Shed: `iuc/fastqc`

### fastq_mott_trimmer

**Trim sequences by quality (Mott's algorithm)** — Sliding-window quality trimming as used by prinseq/dynamictrim.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/reads_single.fastq` | — |
| `threshold` | number | `5.0` | — |
| `window` | number | `0.5` | — |
| `min_length` | int | `20` | — |
| `which_end` | choice | `five_prime` | five_prime, three_prime, both |

Galaxy Tool Shed: `dynamictrim`

### fastq_n_stats

**N content per read position** — Number and percentage of Ns at each cycle (or a single summary).

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/reads_single.fastq` | — |
| `summary` | bool | `False` | — |

Galaxy Tool Shed: `fastx_n_statistics`

### fastq_normalise

**Digital normalisation of k-mer coverage** — Discard reads whose median k-mer coverage exceeds a target (BBNorm-like).

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/reads_single.fastq` | — |
| `k` | int | `21` | — |
| `target_depth` | number | `5.0` | — |

Galaxy Tool Shed: `iuc/bbnorm`

### fastq_overrepresented

**Overrepresented sequences** — Sequences appearing more often than a threshold.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/reads_single.fastq` | — |
| `threshold_percent` | number | `0.1` | — |
| `top` | int | `20` | — |

Galaxy Tool Shed: `iuc/fastqc`

### fastq_per_base_quality

**Per base sequence quality** — Mean/median/quartiles of the quality score at every read position.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/reads_single.fastq` | — |

Galaxy Tool Shed: `iuc/fastqc`

### fastq_per_sequence_quality

**Per sequence quality scores** — Histogram of the mean quality of each read.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/reads_single.fastq` | — |
| `bins` | int | `15` | — |

Galaxy Tool Shed: `iuc/fastqc`

### fastq_polyx_trim

**Trim poly-A / poly-T tails** — Remove homopolymer A/T tails (RNA-seq artefacts).

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/reads_single.fastq` | — |
| `tail` | choice | `polyA` | polyA, polyT, both |
| `min_run` | int | `5` | — |

Galaxy Tool Shed: `iuc/trim_galore`

### fastq_prinseq

**Prinseq-like multi-criteria filter** — Prinseq++-style combined quality/complexity filters.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/reads_single.fastq` | — |
| `min_len` | int | `50` | — |
| `max_ns_frac` | number | `1.0` | — |
| `min_qual_mean` | int | `0` | — |
| `trim_qual_left` | number | `0.0` | — |
| `trim_qual_right` | number | `0.0` | — |
| `derep` | int | `0` | — |
| `dust_threshold` | number | `0.0` | — |
| `max_homopolymer` | int | `0` | — |

Galaxy Tool Shed: `du_novo/prinseq`

### fastq_qc_duplicate_report

**Duplicate rate from a read set** — Library complexity: unique vs duplicated read pairs.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/reads_1.fastq` | — |
| `src2` | file | `` | — |
| `kmer_for_key` | int | `0` | — |

Galaxy Tool Shed: `iuc/fastp`

### fastq_quality_hist

**Quality score distribution** — Counts of each Phred score in the whole file.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/reads_single.fastq` | — |

Galaxy Tool Shed: `fastx_quality_statistics`

### fastq_randomize

**Randomly shuffle reads** — Shuffle read order to remove positional bias.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/reads_single.fastq` | — |
| `seed` | int | `42` | — |

Galaxy Tool Shed: `iuc/bbtools`

### fastq_read_length_distribution

**Read length distribution** — Histogram of read lengths (with optional binning).

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/reads_single.fastq` | — |
| `bins` | int | `0` | — |

Galaxy Tool Shed: `fastq_len`

### fastq_rrbs_trim

**Trim for RRBS (MspI site)** — Cut at the MspI/CTAG site (``trim_galore --rrbs``).

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/reads_single.fastq` | — |
| `site` | text | `CTAG` | — |
| `trim_site` | int | `0` | — |

Galaxy Tool Shed: `iuc/trim_galore`

### fastq_sliding_trim

**Trim Galore! style quality/adapter trim** — Trim bases below a Phred cut-off from both ends, drop short reads.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/reads_single.fastq` | — |
| `quality_cutoff` | int | `20` | — |
| `length_min` | int | `25` | — |
| `trim_ns` | bool | `True` | — |
| `retain_untrimmed` | bool | `False` | — |

Galaxy Tool Shed: `iuc/trim_galore`

### fastq_sortbylength

**Sort reads by length or name** — Deterministic read ordering (SortByLength/SortBySize).

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/reads_single.fastq` | — |
| `by` | choice | `length` | length, name, quality, sequence |
| `largest_first` | bool | `True` | — |

Galaxy Tool Shed: `iuc/fastx_toolkit`

### fastq_subsample

**Subsample reads (bbduk style)** — Draw a random subset of reads (mates kept in sync).

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/reads_single.fastq` | — |
| `n` | int | `0` | — |
| `fraction` | number | `0.1` | — |
| `seed` | int | `1` | — |
| `src2` | file | `` | — |

Galaxy Tool Shed: `iuc/bbtools`

### fastq_tile_stats

**Per-lane / per-tile quality summary** — Aggregate quality by the flowcell/lane field of the read name when present.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/reads_single.fastq` | — |

Galaxy Tool Shed: `fastq_stats`

### fastq_trim_by_length

**Filter reads by length** — Keep or drop reads outside a length window.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/reads_single.fastq` | — |
| `min_len` | int | `0` | — |
| `max_len` | int | `0` | — |
| `mode` | choice | `keep` | keep, discard |

Galaxy Tool Shed: `fastqlenfilter`

### fastq_trim_ns

**Trim Ns from read ends** — Remove N bases from 5' and 3' ends.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/reads_single.fastq` | — |
| `min_length` | int | `0` | — |

Galaxy Tool Shed: `du_novo/prinseq`

### fastqc_summary

**FastQC-style summary report** — Per-base quality, GC, N content, duplication and adapter flags.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/reads_single.fastq` | — |
| `top_sequences` | int | `5` | — |
| `contaminant_k` | int | `5` | — |

Galaxy Tool Shed: `iuc/fastqc`

## sam_bam

*SAM/BAM* — in Galaxy group *Genomic File Manipulation*.

| Tool | What it does | Upstream Galaxy tool |
|---|---|---|
| [`samtools_addreplacerg`](#samtools_addreplacerg) | Stamp every alignment with a new @RG line and RG tag. | `iuc/samtools` |
| [`samtools_bam_sanity`](#samtools_bam_sanity) | Traffic-light table of mapping metrics against thresholds. | `iuc/multiqc` |
| [`samtools_bedgraph`](#samtools_bedgraph) | Depth per base interval as bedGraph (``bedtools genomecov -bg``). | `iuc/samtools` |
| [`samtools_call_variants`](#samtools_call_variants) | Simple pileup genotype caller producing a variant table. | `iuc/samtools` |
| [`samtools_calmd`](#samtools_calmd) | Derive MD/NM-like mismatch strings against the reference. | `iuc/samtools` |
| [`samtools_cigar_summary`](#samtools_cigar_summary) | Totals of M/I/D/S/H/N/P operations in the alignment. | `iuc/samtools` |
| [`samtools_consensus`](#samtools_consensus) | Build a consensus FASTA of the covered regions (``samtools consensus``). | `iuc/samtools` |
| [`samtools_coverage_bed`](#samtools_coverage_bed) | bedtools genomecov -bga style covered/uncovered blocks. | `iuc/bedtools` |
| [`samtools_coverage_histogram`](#samtools_coverage_histogram) | Mean depth in fixed-size bins (``plotCoverage`` data). | `iuc/samtools` |
| [`samtools_depth`](#samtools_depth) | Per-base depth (``samtools depth``) with quality filters. | `iuc/samtools` |
| [`samtools_downsample`](#samtools_downsample) | Randomly keep a fraction of reads. | `iuc/samtools` |
| [`samtools_fastq`](#samtools_fastq) | Convert alignments back into a FASTQ file. | `iuc/samtools` |
| [`samtools_fixmate`](#samtools_fixmate) | Populate RNEXT/PNEXT/TLEN for read pairs. | `iuc/samtools` |
| [`samtools_flagstat`](#samtools_flagstat) | samtools flagstat: counts by flag category. | `iuc/samtools` |
| [`samtools_idxstats`](#samtools_idxstats) | Per-reference mapped/unmapped counts and lengths. | `iuc/samtools` |
| [`samtools_insert_size_hist`](#samtools_insert_size_hist) | Histogram of fragment lengths with mean/SD/3xSD cut-off. | `iuc/samtools` |
| [`samtools_mapq_hist`](#samtools_mapq_hist) | Counts per MAPQ bin. | `iuc/samtools` |
| [`samtools_markdup`](#samtools_markdup) | Identify optical-agnostic duplicates by outer coordinates. | `iuc/samtools` |
| [`samtools_merge`](#samtools_merge) | ``samtools merge``: concatenate headers and alignments. | `iuc/samtools` |
| [`samtools_mismatch_profile`](#samtools_mismatch_profile) | Error rate at each read position (mismatch profile). | `iuc/samtools` |
| [`samtools_nm_distribution`](#samtools_nm_distribution) | Number of alignments per NM value. | `iuc/samtools` |
| [`samtools_pileup`](#samtools_pileup) | Per-base pileup strings (``samtools mpileup``). | — |
| [`samtools_read_bam_nx`](#samtools_read_bam_nx) | One row per alignment: length, NM, clipping, insert size. | `iuc/sambamba` |
| [`samtools_reheader`](#samtools_reheader) | Swap the header of a SAM file (``samtools reheader``). | `iuc/samtools` |
| [`samtools_rmdup_position`](#samtools_rmdup_position) | Drop PCR duplicates sharing the same alignment start. | `iuc/samtools` |
| [`samtools_softclip_stats`](#samtools_softclip_stats) | Clipped bases per read end (breakpoint signal). | `iuc/samtools` |
| [`samtools_sort`](#samtools_sort) | Reorder the alignment records. | `iuc/samtools` |
| [`samtools_stats`](#samtools_stats) | SNPs/indels, insert size, read lengths, error rate. | `iuc/samtools` |
| [`samtools_tag_summary`](#samtools_tag_summary) | Value counts for a SAM tag (RG, NM, AS, …). | `iuc/samtools` |
| [`samtools_tlen_scatter`](#samtools_tlen_scatter) | Per-pair TLEN track (for library QC). | `iuc/samtools` |
| [`samtools_vcf_from_pileup`](#samtools_vcf_from_pileup) | Write a VCF from pileup-based calls. | `iuc/samtools` |
| [`samtools_view`](#samtools_view) | samtools view with region, flag, MAPQ and pairing filters. | `iuc/samtools` |
| [`samtools_view_count`](#samtools_view_count) | ``samtools view -c`` grouped by a chosen key. | `iuc/samtools` |

### samtools_addreplacerg

**Add or replace read groups** — Stamp every alignment with a new @RG line and RG tag.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/alignments.sam` | — |
| `id` | text | `run1` | — |
| `sample` | text | `SAMPLE_X` | — |
| `library` | text | `lib1` | — |
| `platform` | text | `ILLUMINA` | — |

Galaxy Tool Shed: `iuc/samtools`

### samtools_bam_sanity

**Alignment sanity checks (NGS QC)** — Traffic-light table of mapping metrics against thresholds.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/alignments.sam` | — |
| `max_unmapped_fraction` | number | `0.5` | — |
| `min_proper_pair_fraction` | number | `0.5` | — |
| `max_duplicate_fraction` | number | `0.6` | — |

Galaxy Tool Shed: `iuc/multiqc`

### samtools_bedgraph

**Genome coverage to bedGraph** — Depth per base interval as bedGraph (``bedtools genomecov -bg``).

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/alignments.sam` | — |
| `scale` | number | `1.0` | — |
| `fragment` | bool | `True` | — |
| `bin_size` | int | `0` | — |
| `min_depth` | number | `0.0` | — |

Galaxy Tool Shed: `iuc/samtools`

### samtools_call_variants

**Call variants from a pileup** — Simple pileup genotype caller producing a variant table.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/alignments.sam` | — |
| `genome` | file | `examples/genome.fa` | — |
| `min_depth` | int | `5` | — |
| `min_vaf` | number | `0.2` | — |
| `min_qual` | number | `20.0` | — |
| `het_fraction` | number | `0.75` | — |
| `region` | text | `` | — |

Galaxy Tool Shed: `iuc/samtools`

### samtools_calmd

**Compute MD tags** — Derive MD/NM-like mismatch strings against the reference.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/alignments.sam` | — |
| `genome` | file | `examples/genome.fa` | — |
| `write_tags` | bool | `False` | — |

Galaxy Tool Shed: `iuc/samtools`

### samtools_cigar_summary

**CIGAR summary** — Totals of M/I/D/S/H/N/P operations in the alignment.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/alignments.sam` | — |

Galaxy Tool Shed: `iuc/samtools`

### samtools_consensus

**Generate consensus from alignments** — Build a consensus FASTA of the covered regions (``samtools consensus``).

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/alignments.sam` | — |
| `genome` | file | `examples/genome.fa` | — |
| `caller` | text | `samtools` | — |
| `min_depth` | int | `3` | — |
| `report_variants` | bool | `True` | — |

Galaxy Tool Shed: `iuc/samtools`

### samtools_coverage_bed

**Intervals covered by at least N reads** — bedtools genomecov -bga style covered/uncovered blocks.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/alignments.sam` | — |
| `min_depth` | int | `1` | — |
| `max_depth` | number | `0.0` | — |
| `sort` | bool | `True` | — |

Galaxy Tool Shed: `iuc/bedtools`

### samtools_coverage_histogram

**Coverage histogram** — Mean depth in fixed-size bins (``plotCoverage`` data).

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/alignments.sam` | — |
| `binsize` | int | `500` | — |
| `per_base` | bool | `True` | — |

Galaxy Tool Shed: `iuc/samtools`

### samtools_depth

**Depth** — Per-base depth (``samtools depth``) with quality filters.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/alignments.sam` | — |
| `region` | text | `` | — |
| `min_base_quality` | int | `0` | — |
| `min_mapq` | int | `0` | — |
| `max_depth` | int | `0` | — |
| `all_sites` | bool | `False` | — |
| `head` | int | `200` | — |

Galaxy Tool Shed: `iuc/samtools`

### samtools_downsample

**Downsample alignments** — Randomly keep a fraction of reads.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/alignments.sam` | — |
| `fraction` | number | `0.5` | — |
| `n` | int | `0` | — |
| `seed` | int | `42` | — |

Galaxy Tool Shed: `iuc/samtools`

### samtools_fastq

**Extract reads as FASTQ** — Convert alignments back into a FASTQ file.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/alignments.sam` | — |
| `mapped_only` | bool | `True` | — |
| `which` | choice | `all` | all, read1, read2, unmapped |

Galaxy Tool Shed: `iuc/samtools`

### samtools_fixmate

**Fix mate information** — Populate RNEXT/PNEXT/TLEN for read pairs.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/alignments.sam` | — |
| `assume_s1` | bool | `False` | — |

Galaxy Tool Shed: `iuc/samtools`

### samtools_flagstat

**Flagstat** — samtools flagstat: counts by flag category.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/alignments.sam` | — |

Galaxy Tool Shed: `iuc/samtools`

### samtools_idxstats

**Idxstats** — Per-reference mapped/unmapped counts and lengths.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/alignments.sam` | — |

Galaxy Tool Shed: `iuc/samtools`

### samtools_insert_size_hist

**Insert size distribution** — Histogram of fragment lengths with mean/SD/3xSD cut-off.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/alignments.sam` | — |
| `bins` | int | `25` | — |

Galaxy Tool Shed: `iuc/samtools`

### samtools_mapq_hist

**Mapping quality histogram** — Counts per MAPQ bin.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/alignments.sam` | — |
| `bins` | int | `10` | — |

Galaxy Tool Shed: `iuc/samtools`

### samtools_markdup

**Mark duplicates** — Identify optical-agnostic duplicates by outer coordinates.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/alignments.sam` | — |
| `flow_mode` | bool | `False` | — |
| `output_text` | bool | `False` | — |
| `head` | int | `30` | — |

Galaxy Tool Shed: `iuc/samtools`

### samtools_merge

**Merge several SAM files** — ``samtools merge``: concatenate headers and alignments.

| parameter | kind | default | options |
|---|---|---|---|
| `files` | multi | `['examples/alignments.sam']` | examples/alignments.sam |
| `assume_sorted` | bool | `False` | — |
| `add_read_group` | bool | `False` | — |

Galaxy Tool Shed: `iuc/samtools`

### samtools_mismatch_profile

**Mismatch rate per cycle** — Error rate at each read position (mismatch profile).

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/alignments.sam` | — |
| `read_len` | int | `0` | — |
| `min_mapq` | int | `0` | — |

Galaxy Tool Shed: `iuc/samtools`

### samtools_nm_distribution

**Edit-distance (NM) distribution** — Number of alignments per NM value.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/alignments.sam` | — |

Galaxy Tool Shed: `iuc/samtools`

### samtools_pileup

**Mpileup text** — Per-base pileup strings (``samtools mpileup``).

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/alignments.sam` | — |
| `genome` | file | `examples/genome.fa` | — |
| `region` | text | `` | — |
| `min_base_quality` | int | `13` | — |
| `min_mapq` | int | `0` | — |
| `max_rows` | int | `60` | — |

### samtools_read_bam_nx

**Read-level stats table** — One row per alignment: length, NM, clipping, insert size.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/alignments.sam` | — |
| `max_rows` | int | `200` | — |
| `min_mapq` | int | `0` | — |

Galaxy Tool Shed: `iuc/sambamba`

### samtools_reheader

**Replace the SAM header** — Swap the header of a SAM file (``samtools reheader``).

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/alignments.sam` | — |
| `header` | code | `@HD	VN:1.6	SO:coordinate` | — |
| `keep_sq` | bool | `False` | — |

Galaxy Tool Shed: `iuc/samtools`

### samtools_rmdup_position

**Remove duplicate reads by position** — Drop PCR duplicates sharing the same alignment start.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/alignments.sam` | — |
| `keep_highest_mapq` | bool | `True` | — |
| `paired` | bool | `False` | — |

Galaxy Tool Shed: `iuc/samtools`

### samtools_softclip_stats

**Soft-clip statistics** — Clipped bases per read end (breakpoint signal).

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/alignments.sam` | — |
| `min_clip` | int | `5` | — |

Galaxy Tool Shed: `iuc/samtools`

### samtools_sort

**Sort SAM by coordinate or name** — Reorder the alignment records.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/alignments.sam` | — |
| `by` | choice | `coordinate` | coordinate, queryname, unsorted |

Galaxy Tool Shed: `iuc/samtools`

### samtools_stats

**Samtools stats** — SNPs/indels, insert size, read lengths, error rate.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/alignments.sam` | — |
| `genome` | file | `` | — |

Galaxy Tool Shed: `iuc/samtools`

### samtools_tag_summary

**Summarise an alignment tag** — Value counts for a SAM tag (RG, NM, AS, …).

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/alignments.sam` | — |
| `tag` | text | `NM` | — |

Galaxy Tool Shed: `iuc/samtools`

### samtools_tlen_scatter

**Fragment length vs position** — Per-pair TLEN track (for library QC).

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/alignments.sam` | — |
| `max_rows` | int | `400` | — |

Galaxy Tool Shed: `iuc/samtools`

### samtools_vcf_from_pileup

**Pileup to VCF** — Write a VCF from pileup-based calls.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/alignments.sam` | — |
| `genome` | file | `examples/genome.fa` | — |
| `sample` | text | `SAMPLE1` | — |
| `min_depth` | int | `5` | — |
| `min_vaf` | number | `0.25` | — |

Galaxy Tool Shed: `iuc/samtools`

### samtools_view

**View / filter alignments** — samtools view with region, flag, MAPQ and pairing filters.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/alignments.sam` | — |
| `region` | text | `` | — |
| `min_mapq` | int | `0` | — |
| `flag_filter` | int | `0` | — |
| `flag_require` | int | `0` | — |
| `read_type` | choice | `all` | all, mapped, unmapped, paired, proper |
| `drop_secondary` | bool | `True` | — |
| `read_group` | text | `` | — |
| `head` | int | `50` | — |

Galaxy Tool Shed: `iuc/samtools`

### samtools_view_count

**Count reads matching a filter** — ``samtools view -c`` grouped by a chosen key.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/alignments.sam` | — |
| `group_by` | choice | `reference` | reference, read1, flag, mapq, cigar, pair, name |

Galaxy Tool Shed: `iuc/samtools`

## bed

*BED* — in Galaxy group *Genomic File Manipulation*.

| Tool | What it does | Upstream Galaxy tool |
|---|---|---|
| [`bedtools_bamtobed`](#bedtools_bamtobed) | bedtools bamtobed: alignments as BED intervals. | `iuc/bedtools` |
| [`bedtools_bed12_blocks`](#bedtools_bed12_blocks) | Per-interval block summary (BED12 aware). | `iuc/bedtools` |
| [`bedtools_bed_graph_to_bag`](#bedtools_bed_graph_to_bag) | Re-bin a bedGraph into fixed-size windows with a summary statistic. | `iuc/bedtools` |
| [`bedtools_binnify`](#bedtools_binnify) | bedtools binnify: fixed-width bins covering every chromosome. | `iuc/bedtools` |
| [`bedtools_blacklist_filter`](#bedtools_blacklist_filter) | Drop intervals that overlap the blacklist by more than a fraction. | `iuc/bedtools` |
| [`bedtools_cluster`](#bedtools_cluster) | bedtools cluster: assign a cluster id to intervals close together. | `iuc/bedtools` |
| [`bedtools_complement`](#bedtools_complement) | bedtools complement: everything in the genome that is not covered. | `iuc/bedtools` |
| [`bedtools_coverage`](#bedtools_coverage) | bedtools coverage: overlap bases, counts and fraction per A interval. | `iuc/bedtools` |
| [`bedtools_distance`](#bedtools_distance) | bedtools distance: pairwise genomic distances between two files. | `iuc/bedtools` |
| [`bedtools_firstline`](#bedtools_firstline) | bedtools firstline: drop duplicated coordinate rows. | `iuc/bedtools` |
| [`bedtools_flank`](#bedtools_flank) | bedtools flank: upstream/downstream regions relative to the strand. | `iuc/bedtools` |
| [`bedtools_genomecov_binned`](#bedtools_genomecov_binned) | bedtools genomecov -bga style depth in fixed bins. | `iuc/bedtools` |
| [`bedtools_gff_to_bed`](#bedtools_gff_to_bed) | Convert GFF features to BED. | `iuc/gtf2bed` |
| [`bedtools_intersect_bp`](#bedtools_intersect_bp) | bedtools intersect -loj style base-pair totals for the two files. | `iuc/bedtools` |
| [`bedtools_intersect_fractions`](#bedtools_intersect_fractions) | Per-A overlap fraction and base counts (bedtools coverage -d style). | `iuc/bedtools` |
| [`bedtools_intersect_with_counts_by_name`](#bedtools_intersect_with_counts_by_name) | Count how many B intervals hit each named A feature. | `iuc/bedtools` |
| [`bedtools_links`](#bedtools_links) | bedtools links: pairs of intervals within a distance. | `iuc/bedtools` |
| [`bedtools_mask_from_genome`](#bedtools_mask_from_genome) | Return the genome with the given intervals lower-cased or N-masked. | `iuc/bedtools` |
| [`bedtools_merge`](#bedtools_merge) | bedtools merge with distance and column summaries. | `iuc/bedtools` |
| [`bedtools_merge_counts`](#bedtools_merge_counts) | Bin a BED file and count (optionally normalise) intervals per bin. | `iuc/featurecounts` |
| [`bedtools_negativeb`](#bedtools_negativeb) | bedtools negativeb: A intervals minus any that overlap B. | `iuc/bedtools` |
| [`bedtools_random`](#bedtools_random) | bedtools random: generate random intervals of a fixed width. | `iuc/bedtools` |
| [`bedtools_remove_overlaps`](#bedtools_remove_overlaps) | Drop intervals that overlap another interval by at least N bp. | `iuc/bedtools` |
| [`bedtools_restrict_to`](#bedtools_restrict_to) | Keep only intervals that overlap (or do not overlap) a restriction set. | `iuc/bedtools` |
| [`bedtools_scale`](#bedtools_scale) | bedtools scale: scale coordinates relative to genome length. | `iuc/bedtools` |
| [`bedtools_shift`](#bedtools_shift) | bedtools shift: move every interval. | `iuc/bedtools` |
| [`bedtools_shuffle`](#bedtools_shuffle) | bedtools shuffle: randomise interval positions inside the genome. | `iuc/bedtools` |
| [`bedtools_slop`](#bedtools_slop) | bedtools slop: expand intervals, optionally clipped to the genome. | `iuc/bedtools` |
| [`bedtools_summarize`](#bedtools_summarize) | Length statistics of a BED file (total, mean, N50, per-chrom). | `iuc/bedtools` |
| [`bedtools_tidy`](#bedtools_tidy) | bedtools tidy: fix coordinates and extra columns. | `iuc/bedtools` |
| [`bedtools_window_counts`](#bedtools_window_counts) | Number of intervals per fixed window (windowed counts). | `iuc/bedtools` |
| [`bedtools_zero_removal`](#bedtools_zero_removal) | Keep bedGraph segments above (or below) a value, as BED. | `iuc/bedtools` |

### bedtools_bamtobed

**BAM-to-BED** — bedtools bamtobed: alignments as BED intervals.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/alignments.sam` | — |
| `split` | bool | `True` | — |
| `min_mapq` | number | `20.0` | — |
| `flag_filter` | int | `3844` | — |

Galaxy Tool Shed: `iuc/bedtools`

### bedtools_bed12_blocks

**Explain a BED12 file (block stats)** — Per-interval block summary (BED12 aware).

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/regions.bed` | — |
| `as_gtf` | bool | `False` | — |

Galaxy Tool Shed: `iuc/bedtools`

### bedtools_bed_graph_to_bag

**BED graph to genome-wide bins** — Re-bin a bedGraph into fixed-size windows with a summary statistic.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/coverage.bedgraph` | — |
| `binsize` | int | `1000` | — |
| `stat` | choice | `mean` | mean, sum, max, min |

Galaxy Tool Shed: `iuc/bedtools`

### bedtools_binnify

**Binnify genome into equal bins** — bedtools binnify: fixed-width bins covering every chromosome.

| parameter | kind | default | options |
|---|---|---|---|
| `genome` | file | `examples/genome.txt` | — |
| `binsize` | int | `1000` | — |

Galaxy Tool Shed: `iuc/bedtools`

### bedtools_blacklist_filter

**Remove intervals overlapping a blacklist** — Drop intervals that overlap the blacklist by more than a fraction.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/regions.bed` | — |
| `blacklist` | file | `` | — |
| `max_frac` | number | `0.5` | — |

Galaxy Tool Shed: `iuc/bedtools`

### bedtools_cluster

**Cluster** — bedtools cluster: assign a cluster id to intervals close together.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/regions.bed` | — |
| `distance` | int | `1` | — |

Galaxy Tool Shed: `iuc/bedtools`

### bedtools_complement

**Complement** — bedtools complement: everything in the genome that is not covered.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/regions.bed` | — |
| `genome` | file | `examples/genome.txt` | — |

Galaxy Tool Shed: `iuc/bedtools`

### bedtools_coverage

**Coverage** — bedtools coverage: overlap bases, counts and fraction per A interval.

| parameter | kind | default | options |
|---|---|---|---|
| `a` | file | `examples/regions.bed` | — |
| `b` | file | `examples/regions.bed` | — |
| `mode` | choice | `bases` | bases, count, fraction |

Galaxy Tool Shed: `iuc/bedtools`

### bedtools_distance

**Distance** — bedtools distance: pairwise genomic distances between two files.

| parameter | kind | default | options |
|---|---|---|---|
| `a` | file | `examples/regions.bed` | — |
| `b` | file | `examples/regions.bed` | — |

Galaxy Tool Shed: `iuc/bedtools`

### bedtools_firstline

**First line per interval (deduplicate coordinates)** — bedtools firstline: drop duplicated coordinate rows.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/regions.bed` | — |

Galaxy Tool Shed: `iuc/bedtools`

### bedtools_flank

**Flank** — bedtools flank: upstream/downstream regions relative to the strand.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/regions.bed` | — |
| `left` | int | `100` | — |
| `right` | int | `100` | — |
| `both` | bool | `False` | — |

Galaxy Tool Shed: `iuc/bedtools`

### bedtools_genomecov_binned

**Genome coverage in bins** — bedtools genomecov -bga style depth in fixed bins.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/coverage.bedgraph` | — |
| `binsize` | int | `500` | — |
| `max_depth` | int | `0` | — |

Galaxy Tool Shed: `iuc/bedtools`

### bedtools_gff_to_bed

**GFF to BED (operate on features)** — Convert GFF features to BED.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/annotation.gff` | — |
| `feature` | text | `gene` | — |
| `merge_exons` | bool | `False` | — |

Galaxy Tool Shed: `iuc/gtf2bed`

### bedtools_intersect_bp

**Intersect bp** — bedtools intersect -loj style base-pair totals for the two files.

| parameter | kind | default | options |
|---|---|---|---|
| `a` | file | `examples/regions.bed` | — |
| `b` | file | `examples/regions.bed` | — |

Galaxy Tool Shed: `iuc/bedtools`

### bedtools_intersect_fractions

**Fraction of each interval overlapped** — Per-A overlap fraction and base counts (bedtools coverage -d style).

| parameter | kind | default | options |
|---|---|---|---|
| `a` | file | `examples/regions.bed` | — |
| `b` | file | `examples/regions.bed` | — |

Galaxy Tool Shed: `iuc/bedtools`

### bedtools_intersect_with_counts_by_name

**Count overlaps per feature name** — Count how many B intervals hit each named A feature.

| parameter | kind | default | options |
|---|---|---|---|
| `a` | file | `examples/regions.bed` | — |
| `b` | file | `examples/regions.bed` | — |
| `reciprocal` | bool | `False` | — |

Galaxy Tool Shed: `iuc/bedtools`

### bedtools_links

**Links** — bedtools links: pairs of intervals within a distance.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/regions.bed` | — |
| `distance` | int | `1` | — |

Galaxy Tool Shed: `iuc/bedtools`

### bedtools_mask_from_genome

**Mask genome regions (soft/hard)** — Return the genome with the given intervals lower-cased or N-masked.

| parameter | kind | default | options |
|---|---|---|---|
| `genome` | file | `examples/genome.fa` | — |
| `intervals` | file | `examples/regions.bed` | — |
| `mode` | choice | `soft` | soft, hard |

Galaxy Tool Shed: `iuc/bedtools`

### bedtools_merge

**Merge** — bedtools merge with distance and column summaries.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/regions.bed` | — |
| `distance` | int | `-1` | — |
| `strand` | bool | `False` | — |
| `collapse` | choice | `none` | none, concat, count, mean, min, max, sum, distinct |

Galaxy Tool Shed: `iuc/bedtools`

### bedtools_merge_counts

**Count intervals per merged bin** — Bin a BED file and count (optionally normalise) intervals per bin.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/regions.bed` | — |
| `size` | int | `1000` | — |
| `normalisation` | choice | `none` | none, per_kb, fraction_of_chrom |
| `genome` | file | `examples/genome.txt` | — |

Galaxy Tool Shed: `iuc/featurecounts`

### bedtools_negativeb

**Negative B (remove B from A, keep strand)** — bedtools negativeb: A intervals minus any that overlap B.

| parameter | kind | default | options |
|---|---|---|---|
| `a` | file | `examples/regions.bed` | — |
| `b` | file | `examples/regions.bed` | — |

Galaxy Tool Shed: `iuc/bedtools`

### bedtools_random

**Random** — bedtools random: generate random intervals of a fixed width.

| parameter | kind | default | options |
|---|---|---|---|
| `genome` | file | `examples/genome.txt` | — |
| `n` | int | `10` | — |
| `length` | int | `200` | — |
| `seed` | int | `1` | — |
| `chrom` | text | `` | — |

Galaxy Tool Shed: `iuc/bedtools`

### bedtools_remove_overlaps

**Remove overlapping intervals** — Drop intervals that overlap another interval by at least N bp.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/regions.bed` | — |
| `min_overlap` | int | `1` | — |

Galaxy Tool Shed: `iuc/bedtools`

### bedtools_restrict_to

**Restrict ranges to a reference list** — Keep only intervals that overlap (or do not overlap) a restriction set.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/regions.bed` | — |
| `restrict` | file | `examples/targets.bed` | — |
| `invert` | bool | `False` | — |

Galaxy Tool Shed: `iuc/bedtools`

### bedtools_scale

**Scale** — bedtools scale: scale coordinates relative to genome length.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/regions.bed` | — |
| `genome` | file | `examples/genome.txt` | — |
| `factor` | number | `1.5` | — |

Galaxy Tool Shed: `iuc/bedtools`

### bedtools_shift

**Shift** — bedtools shift: move every interval.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/regions.bed` | — |
| `by` | int | `100` | — |
| `strand_aware` | bool | `True` | — |

Galaxy Tool Shed: `iuc/bedtools`

### bedtools_shuffle

**Shuffle** — bedtools shuffle: randomise interval positions inside the genome.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/regions.bed` | — |
| `genome` | file | `examples/genome.txt` | — |
| `seed` | int | `42` | — |
| `exclude` | file | `` | — |
| `include` | bool | `False` | — |

Galaxy Tool Shed: `iuc/bedtools`

### bedtools_slop

**Slop** — bedtools slop: expand intervals, optionally clipped to the genome.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/regions.bed` | — |
| `bp` | int | `100` | — |
| `fraction` | number | `0.0` | — |
| `direction` | choice | `both` | both, left, right |
| `genome_bound` | bool | `False` | — |
| `genome` | file | `examples/genome.txt` | — |

Galaxy Tool Shed: `iuc/bedtools`

### bedtools_summarize

**Interval file summary** — Length statistics of a BED file (total, mean, N50, per-chrom).

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/regions.bed` | — |

Galaxy Tool Shed: `iuc/bedtools`

### bedtools_tidy

**Tidy** — bedtools tidy: fix coordinates and extra columns.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/regions.bed` | — |
| `max_extra` | int | `6` | — |

Galaxy Tool Shed: `iuc/bedtools`

### bedtools_window_counts

**Count intervals in sliding windows** — Number of intervals per fixed window (windowed counts).

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/regions.bed` | — |
| `size` | int | `1000` | — |
| `chrom` | text | `` | — |

Galaxy Tool Shed: `iuc/bedtools`

### bedtools_zero_removal

**Remove zero-coverage regions** — Keep bedGraph segments above (or below) a value, as BED.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/coverage.bedgraph` | — |
| `threshold` | number | `0.0` | — |
| `below` | bool | `True` | — |

Galaxy Tool Shed: `iuc/bedtools`

## vcf_bcf

*VCF/BCF* — in Galaxy group *Genomic File Manipulation*.

| Tool | What it does | Upstream Galaxy tool |
|---|---|---|
| [`bcftools_isec`](#bcftools_isec) | Classify records as shared or unique between two call sets. | `iuc/bcftools` |
| [`bcftools_merge`](#bcftools_merge) | bcftools merge: combine VCFs site-wise and sample-wise. | `iuc/bcftools` |
| [`bcftools_norm`](#bcftools_norm) | bcftools norm: left-align and trim indels against the reference. | `iuc/bcftools` |
| [`bcftools_query`](#bcftools_query) | bcftools query: extract INFO/FORMAT fields as columns. | `iuc/bcftools` |
| [`bcftools_stats`](#bcftools_stats) | Record counts, Ti/Tv, per-sample and per-quality summaries. | `iuc/bcftools` |
| [`bcftools_view`](#bcftools_view) | bcftools view: filter by region, samples and an INFO/FORMAT expression. | `iuc/bcftools` |
| [`vcf_add_info_tags`](#vcf_add_info_tags) | Recompute INFO fields. | `iuc/bcftools` |
| [`vcf_allele_counts`](#vcf_allele_counts) | Per-record REF/ALT allele counts and call number. | `iuc/bcftools` |
| [`vcf_allele_frequencies`](#vcf_allele_frequencies) | bcftools +fill-tags AF and AC per record. | `iuc/bcftools` |
| [`vcf_create_background`](#vcf_create_background) | Compare ALT frequencies of a call set against a population panel. | `iuc/bcftools` |
| [`vcf_haplotype_table`](#vcf_haplotype_table) | Matrix of haplotype strings per sample in a window. | `iuc/beagle` |
| [`vcf_hwe`](#vcf_hwe) | Exact-ish HWE p-values with the number of het/hom genotypes. | `iuc/plink` |
| [`vcf_indel_spectrum`](#vcf_indel_spectrum) | Counts of insertion/deletion lengths (1..n bp). | `iuc/indelspector` |
| [`vcf_ld_pruning`](#vcf_ld_pruning) | Pairwise LD with greedy r2 pruning (``--indep-pairwise``). | `iuc/plink` |
| [`vcf_mendelian_errors`](#vcf_mendelian_errors) | Count mendelian errors per site for a father/mother/child trio. | `iuc/peddy` |
| [`vcf_missingness`](#vcf_missingness) | Missing genotype rate per sample and per site with a filter cut-off. | `iuc/plink` |
| [`vcf_mutation_spectrum`](#vcf_mutation_spectrum) | C>T, G>A … trinucleotide spectrum of SNVs. | `iuc/mSigHDPro` |
| [`vcf_quality_bins`](#vcf_quality_bins) | Histogram of a variant metric. | `iuc/bcftools` |
| [`vcf_region_index`](#vcf_region_index) | Record counts per genomic window (what tabix enables). | `iuc/tabix` |
| [`vcf_relatedness`](#vcf_relatedness) | Pairwise relatedness estimates between samples. | `iuc/king` |
| [`vcf_relatedness_pairs`](#vcf_relatedness_pairs) | Report sample pairs whose kinship implies relatedness. | `iuc/king` |
| [`vcf_sample_qc`](#vcf_sample_qc) | Ti/Tv, het/hom ratio, depth and missingness per sample. | `iuc/gatk` |
| [`vcf_set_gt_ploidy`](#vcf_set_gt_ploidy) | Infer the ploidy from genotype separators. | `iuc/bcftools` |
| [`vcf_sfs`](#vcf_sfs) | 1D SFS of ALT allele counts (foldable). | `iuc/dadi` |
| [`vcf_subsample_samples`](#vcf_subsample_samples) | Keep a random subset (or an explicit list) of samples. | `iuc/bcftools` |
| [`vcf_ts_tv`](#vcf_ts_tv) | Ti/Tv per sample and overall (a QC metric). | `iuc/snpEff` |
| [`vcf_variant_density_plot`](#vcf_variant_density_plot) | Manhattan-style density of variants per window. | `iuc/bedtools` |

### bcftools_isec

**Compare two VCFs (isec)** — Classify records as shared or unique between two call sets.

| parameter | kind | default | options |
|---|---|---|---|
| `a` | file | `examples/variants.vcf` | — |
| `b` | file | `examples/variants.vcf` | — |
| `sites_match` | bool | `True` | — |
| `mode` | choice | `all` | all, unique_to_a, unique_to_b, both |

Galaxy Tool Shed: `iuc/bcftools`

### bcftools_merge

**Merge call sets (union of samples)** — bcftools merge: combine VCFs site-wise and sample-wise.

| parameter | kind | default | options |
|---|---|---|---|
| `files` | multi | `['examples/variants.vcf']` | examples/variants.vcf |
| `strategy` | choice | `union` | union, intersect, complement |
| `dedup_sites` | bool | `True` | — |

Galaxy Tool Shed: `iuc/bcftools`

### bcftools_norm

**Normalise indels (left align)** — bcftools norm: left-align and trim indels against the reference.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/variants.vcf` | — |
| `genome` | file | `examples/genome.fa` | — |
| `multiallelics` | bool | `False` | — |

Galaxy Tool Shed: `iuc/bcftools`

### bcftools_query

**Query VCF to a table** — bcftools query: extract INFO/FORMAT fields as columns.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/variants.vcf` | — |
| `info_keys` | multi | `[]` | — |
| `format_keys` | multi | `[]` | — |
| `header` | bool | `True` | — |

Galaxy Tool Shed: `iuc/bcftools`

### bcftools_stats

**bcftools stats** — Record counts, Ti/Tv, per-sample and per-quality summaries.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/variants.vcf` | — |

Galaxy Tool Shed: `iuc/bcftools`

### bcftools_view

**View / subset a VCF** — bcftools view: filter by region, samples and an INFO/FORMAT expression.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/variants.vcf` | — |
| `region` | text | `` | — |
| `samples` | text | `` | — |
| `include` | text | `` | — |
| `exclude_mode` | bool | `False` | — |
| `max_records` | text | `` | — |

Galaxy Tool Shed: `iuc/bcftools`

### vcf_add_info_tags

**Fill INFO tags (AC/AN/AF/DP)** — Recompute INFO fields.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/variants.vcf` | — |
| `tags` | multi | `['AC', 'AN', 'AF', 'NS']` | AC, AN, AF, DP, NS, MQ |
| `overwrite` | bool | `True` | — |

Galaxy Tool Shed: `iuc/bcftools`

### vcf_allele_counts

**Allele counts (AC/AN) table** — Per-record REF/ALT allele counts and call number.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/variants.vcf` | — |
| `max_records` | int | `200` | — |

Galaxy Tool Shed: `iuc/bcftools`

### vcf_allele_frequencies

**Allele frequency per site** — bcftools +fill-tags AF and AC per record.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/variants.vcf` | — |
| `include_monomorphic` | bool | `False` | — |
| `min_af` | number | `0.0` | — |
| `max_af` | number | `1.0` | — |

Galaxy Tool Shed: `iuc/bcftools`

### vcf_create_background

**VCF background frequencies from a panel** — Compare ALT frequencies of a call set against a population panel.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/variants.vcf` | — |
| `panel` | file | `examples/variants.vcf` | — |
| `af_diff` | number | `0.1` | — |

Galaxy Tool Shed: `iuc/bcftools`

### vcf_haplotype_table

**Sample haplotypes per region** — Matrix of haplotype strings per sample in a window.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/variants.vcf` | — |
| `region` | text | `` | — |
| `phase_by_pipe` | bool | `True` | — |

Galaxy Tool Shed: `iuc/beagle`

### vcf_hwe

**Hardy-Weinberg test per site** — Exact-ish HWE p-values with the number of het/hom genotypes.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/variants.vcf` | — |
| `min_p` | number | `1e-06` | — |
| `controls_only` | bool | `False` | — |

Galaxy Tool Shed: `iuc/plink`

### vcf_indel_spectrum

**Indel length spectrum** — Counts of insertion/deletion lengths (1..n bp).

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/variants.vcf` | — |

Galaxy Tool Shed: `iuc/indelspector`

### vcf_ld_pruning

**LD r2 and pruning** — Pairwise LD with greedy r2 pruning (``--indep-pairwise``).

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/variants.vcf` | — |
| `r2_cutoff` | number | `0.2` | — |
| `window` | int | `10000` | — |
| `max_pairs` | int | `5000` | — |
| `report_pruned` | bool | `True` | — |

Galaxy Tool Shed: `iuc/plink`

### vcf_mendelian_errors

**Mendelian consistency in trios** — Count mendelian errors per site for a father/mother/child trio.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/variants.vcf` | — |
| `father` | text | `` | — |
| `mother` | text | `` | — |
| `child` | text | `` | — |
| `strict` | bool | `False` | — |

Galaxy Tool Shed: `iuc/peddy`

### vcf_missingness

**Missing data per site and sample** — Missing genotype rate per sample and per site with a filter cut-off.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/variants.vcf` | — |
| `max_missing` | number | `0.2` | — |

Galaxy Tool Shed: `iuc/plink`

### vcf_mutation_spectrum

**Mutation spectrum (6 classes)** — C>T, G>A … trinucleotide spectrum of SNVs.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/variants.vcf` | — |
| `genome` | file | `examples/genome.fa` | — |
| `context` | int | `1` | — |
| `normalize` | bool | `True` | — |

Galaxy Tool Shed: `iuc/mSigHDPro`

### vcf_quality_bins

**QUAL / DP / GQ binning** — Histogram of a variant metric.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/variants.vcf` | — |
| `bins` | int | `10` | — |
| `metric` | choice | `QUAL` | QUAL, DP, GQ, AF |

Galaxy Tool Shed: `iuc/bcftools`

### vcf_region_index

**Tabix-style region index** — Record counts per genomic window (what tabix enables).

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/variants.vcf` | — |
| `window` | int | `1000` | — |

Galaxy Tool Shed: `iuc/tabix`

### vcf_relatedness

**Kinship / relatedness matrix** — Pairwise relatedness estimates between samples.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/variants.vcf` | — |
| `method` | choice | `king` | king, concordance, ibs |
| `heat` | bool | `True` | — |

Galaxy Tool Shed: `iuc/king`

### vcf_relatedness_pairs

**Related pairs above a threshold** — Report sample pairs whose kinship implies relatedness.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/variants.vcf` | — |
| `threshold` | number | `0.0884` | — |
| `relationship` | choice | `auto` | auto, duplicate, parent_child, sibling |

Galaxy Tool Shed: `iuc/king`

### vcf_sample_qc

**Per-sample call QC** — Ti/Tv, het/hom ratio, depth and missingness per sample.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/variants.vcf` | — |
| `min_dp` | number | `3` | — |
| `ts_tv_expected` | number | `2.0` | — |

Galaxy Tool Shed: `iuc/gatk`

### vcf_set_gt_ploidy

**Ploidy inference and GT fix** — Infer the ploidy from genotype separators.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/variants.vcf` | — |

Galaxy Tool Shed: `iuc/bcftools`

### vcf_sfs

**Site frequency spectrum** — 1D SFS of ALT allele counts (foldable).

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/variants.vcf` | — |
| `derived_only` | bool | `False` | — |
| `fold` | bool | `False` | — |

Galaxy Tool Shed: `iuc/dadi`

### vcf_subsample_samples

**Subsample VCF individuals** — Keep a random subset (or an explicit list) of samples.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/variants.vcf` | — |
| `n` | int | `3` | — |
| `seed` | int | `1` | — |
| `keep` | text | `` | — |

Galaxy Tool Shed: `iuc/bcftools`

### vcf_ts_tv

**Transition / transversion ratio** — Ti/Tv per sample and overall (a QC metric).

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/variants.vcf` | — |

Galaxy Tool Shed: `iuc/snpEff`

### vcf_variant_density_plot

**Variant density along chromosomes** — Manhattan-style density of variants per window.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/variants.vcf` | — |
| `window` | int | `1000` | — |

Galaxy Tool Shed: `iuc/bedtools`

## nanopore

*Nanopore* — in Galaxy group *Genomic File Manipulation*.

| Tool | What it does | Upstream Galaxy tool |
|---|---|---|
| [`longread_align`](#longread_align) | Minimap2-style long-read mapping with tolerant seed extension. | `minimap2` |
| [`nanopore_fast5_proxy`](#nanopore_fast5_proxy) | Quality moving average along each read (NanoPack-style diagnostics). | `iuc/nanopack` |
| [`nanopore_filter`](#nanopore_filter) | ``filtlong``-style selection of reads. | `iuc/filtlong` |
| [`nanopore_kmer_quals`](#nanopore_kmer_quals) | Mean quality per k-mer context (homopolymer/systematic errors). | `iuc/nanopore` |
| [`nanopore_porechop`](#nanopore_porechop) | Detect adapter matches and report the trimmed length of each read. | `iuc/porechop` |
| [`nanopore_quality_by_length`](#nanopore_quality_by_length) | Binned mean quality as a function of read length. | `iuc/nanoq` |
| [`nanopore_quality_profile`](#nanopore_quality_profile) | Mean base quality at every cycle of the run. | `porechop / NanoPlot` |
| [`nanopore_read_table`](#nanopore_read_table) | One row per read: length, quality, GC, Ns, homopolymers. | `iuc/poretools` |
| [`nanopore_run_stats`](#nanopore_run_stats) | Total bases, N50, mean quality and Q20 fraction of a long-read run. | `nanoplot / poretools` |
| [`nanopore_stats`](#nanopore_stats) | Read count, total bases, N50 and quality of long reads. | `iuc/nanoq` |
| [`nanopore_trim_adapters`](#nanopore_trim_adapters) | Remove adapter occurrences and quality-trim read ends, keeping long enough reads. | `porechop / fastp` |

### longread_align

**Map long reads onto a reference** — Minimap2-style long-read mapping with tolerant seed extension.

| parameter | kind | default | options |
|---|---|---|---|
| `reads` | file | `examples/reads_single.fastq` | — |
| `ref` | file | `examples/genome.fa` | — |
| `k` | int | `15` | — |
| `max_mismatches` | int | `12` | — |
| `head` | int | `100` | — |

Galaxy Tool Shed: `minimap2`

### nanopore_fast5_proxy

**Signal-level proxy statistics from FASTQ** — Quality moving average along each read (NanoPack-style diagnostics).

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/reads_single.fastq` | — |
| `window` | int | `25` | — |

Galaxy Tool Shed: `iuc/nanopack`

### nanopore_filter

**Filter long reads by length and quality** — ``filtlong``-style selection of reads.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/reads_single.fastq` | — |
| `min_length` | int | `100` | — |
| `min_mean_quality` | number | `7.0` | — |
| `max_length` | int | `0` | — |

Galaxy Tool Shed: `iuc/filtlong`

### nanopore_kmer_quals

**Nanopore k-mer quality bias** — Mean quality per k-mer context (homopolymer/systematic errors).

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/reads_single.fastq` | — |
| `k` | int | `5` | — |
| `top` | int | `15` | — |

Galaxy Tool Shed: `iuc/nanopore`

### nanopore_porechop

**Porechop-style adapter splitting** — Detect adapter matches and report the trimmed length of each read.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/reads_single.fastq` | — |
| `adapter_5p` | text | `AATGTACTTCGTTCAGTTACGTATTGCT` | — |
| `adapter_3p` | text | `GCAATACGTAACTGAACGAAGT` | — |
| `min_subseq` | int | `50` | — |

Galaxy Tool Shed: `iuc/porechop`

### nanopore_quality_by_length

**Quality vs read length** — Binned mean quality as a function of read length.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/reads_single.fastq` | — |
| `bins` | int | `10` | — |

Galaxy Tool Shed: `iuc/nanoq`

### nanopore_quality_profile

**Per-cycle quality profile** — Mean base quality at every cycle of the run.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/reads_single.fastq` | — |
| `min_length` | int | `0` | — |

Galaxy Tool Shed: `porechop / NanoPlot`

### nanopore_read_table

**Read-level summary table** — One row per read: length, quality, GC, Ns, homopolymers.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/reads_single.fastq` | — |
| `max_rows` | int | `500` | — |

Galaxy Tool Shed: `iuc/poretools`

### nanopore_run_stats

**Long-read length and quality summary** — Total bases, N50, mean quality and Q20 fraction of a long-read run.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/reads_single.fastq` | — |
| `min_length` | int | `0` | — |
| `bins` | int | `8` | — |

Galaxy Tool Shed: `nanoplot / poretools`

### nanopore_stats

**Nanopore read statistics** — Read count, total bases, N50 and quality of long reads.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/reads_single.fastq` | — |
| `min_length` | int | `200` | — |

Galaxy Tool Shed: `iuc/nanoq`

### nanopore_trim_adapters

**Trim adapters and low-quality tails** — Remove adapter occurrences and quality-trim read ends, keeping long enough reads.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/reads_single.fastq` | — |
| `adapters` | code | `AGATCGGAAGAGC` | — |
| `quality_cut` | int | `7` | — |
| `min_length` | int | `50` | — |

Galaxy Tool Shed: `porechop / fastp`

## convert_formats

*Convert Formats* — in Galaxy group *Genomic File Manipulation*.

| Tool | What it does | Upstream Galaxy tool |
|---|---|---|
| [`convert_bed_to_fasta`](#convert_bed_to_fasta) | Extract the sequence of every interval from a genome FASTA. | `bedtools_getfasta` |
| [`convert_bed_to_gff`](#convert_bed_to_gff) | BED → GFF3 with name attributes. | `bed2gff` |
| [`convert_gff_to_bed`](#convert_gff_to_bed) | Turn annotated gene/exon features into BED blocks. | `gff2bed` |
| [`convert_matrix_to_dist`](#convert_matrix_to_dist) | Rewrite a pairwise distance matrix between the usual layouts. | `convert_matrix` |
| [`convert_newick_to_table`](#convert_newick_to_table) | Flatten a tree into a node/parent/branch-length table. | `newick2table` |
| [`convert_table_to_bed`](#convert_table_to_bed) | Map arbitrary table columns onto BED fields. | `table_to_bed` |
| [`convert_table_to_newick`](#convert_table_to_newick) | Rebuild a Newick string from a parent/child edge list. | `table2newick` |
| [`fasta_to_fastq_convert`](#fasta_to_fastq_convert) | Drop or synthesise the quality lines when converting. | `fastq_to_fasta` |
| [`fastq_convert_illumina18`](#fastq_convert_illumina18) | Shift ASCII quality offsets (Illumina 1.3/1.5/1.8, Solexa). | ` fastaq` |
| [`fastq_to_sam_like`](#fastq_to_sam_like) | Emit a SAM header for the reference plus unmapped records. | `fastq_to_sam` |
| [`table_to_fasta`](#table_to_fasta) | Build a FASTA (or FASTQ) file from two columns of a table. | `table_to_fasta` |
| [`vcf_convert_to_bed`](#vcf_convert_to_bed) | Variants as BED intervals (1 bp for SNVs, span for indels). | `iuc/bcftools` |
| [`vcf_convert_to_tsv`](#vcf_convert_to_tsv) | One row per variant with a column per sample genotype. | `iuc/vcf2txt` |
| [`vcf_to_genotypes_tsv`](#vcf_to_genotypes_tsv) | VCF → sample x variant genotype matrix (plink-style). | `iuc/plink` |

### convert_bed_to_fasta

**BED intervals to FASTA (genome required)** — Extract the sequence of every interval from a genome FASTA.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/regions.bed` | — |
| `genome` | file | `examples/genome.fa` | — |
| `name_from_score` | bool | `False` | — |

Galaxy Tool Shed: `bedtools_getfasta`

### convert_bed_to_gff

**Convert BED to GFF3** — BED → GFF3 with name attributes.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/regions.bed` | — |
| `source` | text | `chroma_titan` | — |
| `feature` | text | `region` | — |

Galaxy Tool Shed: `bed2gff`

### convert_gff_to_bed

**Convert GFF3 to BED12** — Turn annotated gene/exon features into BED blocks.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/annotation.gff` | — |
| `feature` | text | `mRNA` | — |
| `merge_exons` | bool | `True` | — |

Galaxy Tool Shed: `gff2bed`

### convert_matrix_to_dist

**Distance matrix formats (square ↔ phylip/lower)** — Rewrite a pairwise distance matrix between the usual layouts.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/msa.fasta` | — |
| `to` | choice | `phylip_lower` | phylip_lower, square, vector |

Galaxy Tool Shed: `convert_matrix`

### convert_newick_to_table

**Newick tree to branch table** — Flatten a tree into a node/parent/branch-length table.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/tree.nwk` | — |

Galaxy Tool Shed: `newick2table`

### convert_table_to_bed

**Tabular to BED** — Map arbitrary table columns onto BED fields.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/counts.tsv` | — |
| `chrom_col` | text | `1` | — |
| `start_col` | text | `2` | — |
| `end_col` | text | `3` | — |
| `name_col` | text | `4` | — |
| `strand_col` | text | `6` | — |
| `start_is_1_based` | bool | `False` | — |

Galaxy Tool Shed: `table_to_bed`

### convert_table_to_newick

**Table of branches to Newick** — Rebuild a Newick string from a parent/child edge list.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/counts.tsv` | — |
| `parent_col` | text | `parent` | — |
| `child_col` | text | `child` | — |
| `length_col` | text | `length` | — |

Galaxy Tool Shed: `table2newick`

### fasta_to_fastq_convert

**Convert between FASTA and FASTQ** — Drop or synthesise the quality lines when converting.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/reads_single.fastq` | — |
| `quality_mode` | choice | `constant` | constant, index_derived, zero |
| `qual` | text | `I` | — |

Galaxy Tool Shed: `fastq_to_fasta`

### fastq_convert_illumina18

**Convert FASTQ quality encoding** — Shift ASCII quality offsets (Illumina 1.3/1.5/1.8, Solexa).

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/reads_single.fastq` | — |
| `from_offset` | choice | `33` | 33, 64 |
| `to_offset` | choice | `33` | 33, 64 |
| `max_quality` | int | `41` | — |

Galaxy Tool Shed: ` fastaq`

### fastq_to_sam_like

**Convert FASTQ to an unmapped SAM** — Emit a SAM header for the reference plus unmapped records.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/reads_single.fastq` | — |
| `genome` | file | `examples/genome.fa` | — |
| `expected_length` | int | `0` | — |

Galaxy Tool Shed: `fastq_to_sam`

### table_to_fasta

**Tabular to FASTA** — Build a FASTA (or FASTQ) file from two columns of a table.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/counts.tsv` | — |
| `id_column` | text | `1` | — |
| `sequence_column` | text | `4` | — |
| `format` | choice | `fasta` | fasta, fastq |
| `quality` | text | `I` | — |

Galaxy Tool Shed: `table_to_fasta`

### vcf_convert_to_bed

**VCF to BED** — Variants as BED intervals (1 bp for SNVs, span for indels).

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/variants.vcf` | — |
| `snvs_only` | bool | `False` | — |
| `require_sample` | text | `` | — |

Galaxy Tool Shed: `iuc/bcftools`

### vcf_convert_to_tsv

**VCF to wide tabular** — One row per variant with a column per sample genotype.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/variants.vcf` | — |
| `info_keys` | multi | `[]` | — |
| `genotypes` | bool | `True` | — |

Galaxy Tool Shed: `iuc/vcf2txt`

### vcf_to_genotypes_tsv

**Export genotype matrix** — VCF → sample x variant genotype matrix (plink-style).

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/variants.vcf` | — |
| `coding` | choice | `alt_count` | alt_count, 012, ref_alt, freq |
| `transposed` | bool | `False` | — |

Galaxy Tool Shed: `iuc/plink`

## lift_over

*Lift-Over* — in Galaxy group *Genomic File Manipulation*.

| Tool | What it does | Upstream Galaxy tool |
|---|---|---|
| [`lift_over_chain`](#lift_over_chain) | Re-map interval coordinates between assemblies using aligned chain blocks. | `UCSC liftOver` |
| [`liftover_chain_bed`](#liftover_chain_bed) | Move intervals between assemblies with a chain or offset. | `iuc/ucsc_liftover` |
| [`liftover_chain_summary`](#liftover_chain_summary) | Report chain blocks per chromosome pair. | `iuc/ucsc_liftover` |
| [`liftover_pdb_chain`](#liftover_pdb_chain) | Shift variant coordinates between assemblies and report per-record status. | `iuc/ucsc_liftover` |

### lift_over_chain

**Lift BED intervals through a chain file** — Re-map interval coordinates between assemblies using aligned chain blocks.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/regions.bed` | — |
| `chain` | file | `` | — |
| `mode` | choice | `smart` | smart, start, end |
| `fail` | choice | `keep` | drop, keep |

Galaxy Tool Shed: `UCSC liftOver`

### liftover_chain_bed

**Lift Over BED with a chain file** — Move intervals between assemblies with a chain or offset.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/regions.bed` | — |
| `chain` | file | `` | — |
| `mode` | choice | `dotter` | dotter, ucsc_chain, manual_offset |
| `offset` | int | `0` | — |
| `fail_to_original` | bool | `True` | — |

Galaxy Tool Shed: `iuc/ucsc_liftover`

### liftover_chain_summary

**Chain QC: coverage of the query assembly** — Report chain blocks per chromosome pair.

| parameter | kind | default | options |
|---|---|---|---|
| `chain` | file | `` | — |
| `top` | int | `20` | — |

Galaxy Tool Shed: `iuc/ucsc_liftover`

### liftover_pdb_chain

**LiftOver positions (VCF style)** — Shift variant coordinates between assemblies and report per-record status.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/variants.vcf` | — |
| `chain` | file | `` | — |
| `offset` | int | `0` | — |
| `min_mapping_rate` | number | `0.5` | — |

Galaxy Tool Shed: `iuc/ucsc_liftover`

## operate_on_genomic_intervals

*Operate on Genomic Intervals* — in Galaxy group *Common Genomics Tools*.

| Tool | What it does | Upstream Galaxy tool |
|---|---|---|
| [`bedtools_closest`](#bedtools_closest) | bedtools closest: the nearest B interval for each A interval. | `iuc/bedtools` |
| [`bedtools_intersect`](#bedtools_intersect) | bedtools intersect: overlaps between two interval files. | `iuc/bedtools` |
| [`bedtools_intersect_intervals_only`](#bedtools_intersect_intervals_only) | bedtools intersect -u: A intervals that overlap any B interval. | `iuc/bedtools` |
| [`bedtools_intersect_v`](#bedtools_intersect_v) | bedtools intersect -v: A intervals with no overlap in B. | `iuc/bedtools` |
| [`bedtools_jaccard`](#bedtools_jaccard) | bedtools jaccard plus the Fisher/ t-tests on interval overlap. | `iuc/bedtools` |
| [`bedtools_multiinter`](#bedtools_multiinter) | bedtools multiinter: which files cover each genomic bin. | `iuc/bedtools` |
| [`bedtools_subtract`](#bedtools_subtract) | bedtools subtract: remove the parts of A covered by B. | `iuc/bedtools` |
| [`bedtools_union_bed_graphs`](#bedtools_union_bed_graphs) | bedtools unionbedg-style coverage matrix over all intervals. | `iuc/bedtools` |
| [`bedtools_window`](#bedtools_window) | bedtools window: B intervals within W bp of A intervals. | `iuc/bedtools` |

### bedtools_closest

**Closest** — bedtools closest: the nearest B interval for each A interval.

| parameter | kind | default | options |
|---|---|---|---|
| `a` | file | `examples/regions.bed` | — |
| `b` | file | `examples/regions.bed` | — |
| `distance_only` | bool | `False` | — |

Galaxy Tool Shed: `iuc/bedtools`

### bedtools_intersect

**Intersect** — bedtools intersect: overlaps between two interval files.

| parameter | kind | default | options |
|---|---|---|---|
| `a` | file | `examples/regions.bed` | — |
| `b` | file | `examples/regions.bed` | — |
| `wa` | choice | `none` | none, wa, wb, both |
| `strand` | bool | `False` | — |
| `report` | choice | `rows` | rows, unique_a, count_a, invert |

Galaxy Tool Shed: `iuc/bedtools`

### bedtools_intersect_intervals_only

**Intersect intervals (geometry only)** — bedtools intersect -u: A intervals that overlap any B interval.

| parameter | kind | default | options |
|---|---|---|---|
| `a` | file | `examples/regions.bed` | — |
| `b` | file | `examples/regions.bed` | — |
| `strand` | bool | `False` | — |

Galaxy Tool Shed: `iuc/bedtools`

### bedtools_intersect_v

**Intersect -v (only non-overlapping)** — bedtools intersect -v: A intervals with no overlap in B.

| parameter | kind | default | options |
|---|---|---|---|
| `a` | file | `examples/regions.bed` | — |
| `b` | file | `examples/regions.bed` | — |
| `strand` | bool | `False` | — |

Galaxy Tool Shed: `iuc/bedtools`

### bedtools_jaccard

**Jaccard and Fisher tests** — bedtools jaccard plus the Fisher/ t-tests on interval overlap.

| parameter | kind | default | options |
|---|---|---|---|
| `a` | file | `examples/regions.bed` | — |
| `b` | file | `examples/regions.bed` | — |

Galaxy Tool Shed: `iuc/bedtools`

### bedtools_multiinter

**MultiIntersect** — bedtools multiinter: which files cover each genomic bin.

| parameter | kind | default | options |
|---|---|---|---|
| `files` | multi | `['examples/regions.bed', 'examples/targets.bed']` | examples/regions.bed, examples/targets.bed |
| `header` | bool | `True` | — |

Galaxy Tool Shed: `iuc/bedtools`

### bedtools_subtract

**Subtract** — bedtools subtract: remove the parts of A covered by B.

| parameter | kind | default | options |
|---|---|---|---|
| `a` | file | `examples/regions.bed` | — |
| `b` | file | `examples/regions.bed` | — |
| `whole_A` | bool | `False` | — |

Galaxy Tool Shed: `iuc/bedtools`

### bedtools_union_bed_graphs

**Union of multiple interval files** — bedtools unionbedg-style coverage matrix over all intervals.

| parameter | kind | default | options |
|---|---|---|---|
| `files` | multi | `['examples/regions.bed', 'examples/targets.bed']` | examples/regions.bed, examples/targets.bed |

Galaxy Tool Shed: `iuc/bedtools`

### bedtools_window

**Window** — bedtools window: B intervals within W bp of A intervals.

| parameter | kind | default | options |
|---|---|---|---|
| `a` | file | `examples/regions.bed` | — |
| `b` | file | `examples/regions.bed` | — |
| `window` | int | `500` | — |

Galaxy Tool Shed: `iuc/bedtools`

## fetch_sequences_alignments

*Fetch Sequences and Alignments* — in Galaxy group *Common Genomics Tools*.

| Tool | What it does | Upstream Galaxy tool |
|---|---|---|
| [`bedtools_getfasta`](#bedtools_getfasta) | bedtools getfasta: extract the sequence of each interval. | `iuc/bedtools` |
| [`bedtools_getfasta_split`](#bedtools_getfasta_split) | bedtools getfasta -split style extraction with padding. | `iuc/bedtools` |
| [`bedtools_nuc_content`](#bedtools_nuc_content) | Per-interval nucleotide composition (GC%, CpG obs/exp). | `iuc/bedtools` |
| [`fasta_extract_region`](#fasta_extract_region) | Pull chrom:start-end slices out of a genome FASTA (``getfasta``). | `bedtools_getfasta` |
| [`gff3_loci_to_fasta`](#gff3_loci_to_fasta) | Pull the genomic sequence of every feature, optionally strand-aware with flanks. | `bedtools getfasta / gff3read --loci` |

### bedtools_getfasta

**Get Fasta (BED intervals to sequence)** — bedtools getfasta: extract the sequence of each interval.

| parameter | kind | default | options |
|---|---|---|---|
| `genome` | file | `examples/genome.fa` | — |
| `intervals` | file | `examples/regions.bed` | — |
| `strand` | choice | `same` | same, plus, minus |
| `name_mode` | choice | `bed` | bed, sequence, chrom_only |
| `split` | bool | `False` | — |

Galaxy Tool Shed: `iuc/bedtools`

### bedtools_getfasta_split

**Split BED12 blocks into FASTA** — bedtools getfasta -split style extraction with padding.

| parameter | kind | default | options |
|---|---|---|---|
| `genome` | file | `examples/genome.fa` | — |
| `intervals` | file | `examples/regions.bed` | — |
| `padding` | int | `0` | — |

Galaxy Tool Shed: `iuc/bedtools`

### bedtools_nuc_content

**Nuc content of intervals** — Per-interval nucleotide composition (GC%, CpG obs/exp).

| parameter | kind | default | options |
|---|---|---|---|
| `genome` | file | `examples/genome.fa` | — |
| `intervals` | file | `examples/regions.bed` | — |

Galaxy Tool Shed: `iuc/bedtools`

### fasta_extract_region

**Extract sequence by region** — Pull chrom:start-end slices out of a genome FASTA (``getfasta``).

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/genome.fa` | — |
| `regions` | text | `chrV:100-200` | — |
| `strand` | choice | `both` | both, plus, minus |
| `flank` | bool | `False` | — |
| `extend` | int | `0` | — |

Galaxy Tool Shed: `bedtools_getfasta`

### gff3_loci_to_fasta

**Extract sequences for annotated features** — Pull the genomic sequence of every feature, optionally strand-aware with flanks.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/annotation.gff` | — |
| `ref` | file | `examples/genome.fa` | — |
| `type` | choice | `gene` | gene, mRNA, exon, CDS |
| `strand` | choice | `same` | same, plus, minus |
| `extend` | int | `0` | — |

Galaxy Tool Shed: `bedtools getfasta / gff3read --loci`

## assembly

*Assembly* — in Galaxy group *Genomics Analysis*.

| Tool | What it does | Upstream Galaxy tool |
|---|---|---|
| [`assembly_completeness`](#assembly_completeness) | BLAST each expected gene against the assembly and score the recovered fraction. | `BUSCO-lite` |
| [`assembly_contigs_fasta`](#assembly_contigs_fasta) | Export assembled contigs in FASTA with length and coverage in the headers. | `spades contigs.fasta` |
| [`assembly_metrics`](#assembly_metrics) | Assembly summary: sizes, N50/L50, GC and a contig length histogram. | `quast / asmstats` |
| [`contig_nuc_content`](#contig_nuc_content) | Per-contig GC, N content, ambiguity and longest homopolymer. | `seqkit stats / emboss geece` |
| [`contig_read_coverage`](#contig_read_coverage) | Per-contig mapped-read counts, mean depth and zero-coverage fraction. | `bwa + samtools depth` |
| [`debruijn_unitigs`](#debruijn_unitigs) | Build a de Bruijn graph of read k-mers and emit maximal non-branching unitigs. | `spades / velvet` |
| [`greedy_overlap_assembly`](#greedy_overlap_assembly) | Repeatedly join the read pair with the largest exact overlap to build contigs. | `cap3 / phrap` |
| [`kmer_spectrum`](#kmer_spectrum) | Histogram of k-mer multiplicities - the error peak and genomic peak used to pick k. | `kmergenie / jellyfish hist` |
| [`scaffold_with_pairs`](#scaffold_with_pairs) | Count read pairs bridging contig pairs to infer scaffold adjacency. | `scaffolds / OPERA` |

### assembly_completeness

**Search conserved genes in an assembly** — BLAST each expected gene against the assembly and score the recovered fraction.

| parameter | kind | default | options |
|---|---|---|---|
| `contigs` | file | `examples/genome.fa` | — |
| `genes` | file | `examples/genes.fa` | — |
| `word_size` | int | `11` | — |
| `min_identity` | number | `80.0` | — |

Galaxy Tool Shed: `BUSCO-lite`

### assembly_contigs_fasta

**Write the unitig assembly as FASTA** — Export assembled contigs in FASTA with length and coverage in the headers.

| parameter | kind | default | options |
|---|---|---|---|
| `reads` | file | `examples/reads_single.fastq` | — |
| `k` | int | `15` | — |
| `min_length` | int | `40` | — |
| `orientation` | bool | `False` | — |

Galaxy Tool Shed: `spades contigs.fasta`

### assembly_metrics

**N50, L50 and contig size distribution** — Assembly summary: sizes, N50/L50, GC and a contig length histogram.

| parameter | kind | default | options |
|---|---|---|---|
| `contigs` | file | `examples/genome.fa` | — |
| `bins` | int | `8` | — |
| `min_length` | int | `0` | — |

Galaxy Tool Shed: `quast / asmstats`

### contig_nuc_content

**Nucleotide composition of each contig** — Per-contig GC, N content, ambiguity and longest homopolymer.

| parameter | kind | default | options |
|---|---|---|---|
| `contigs` | file | `examples/genome.fa` | — |
| `window` | int | `0` | — |

Galaxy Tool Shed: `seqkit stats / emboss geece`

### contig_read_coverage

**Map reads back onto contigs for coverage** — Per-contig mapped-read counts, mean depth and zero-coverage fraction.

| parameter | kind | default | options |
|---|---|---|---|
| `reads` | file | `examples/reads_single.fastq` | — |
| `contigs` | file | `examples/genome.fa` | — |
| `k` | int | `21` | — |
| `bins` | int | `10` | — |
| `head` | int | `100` | — |

Galaxy Tool Shed: `bwa + samtools depth`

### debruijn_unitigs

**De Bruijn graph unitigs from reads** — Build a de Bruijn graph of read k-mers and emit maximal non-branching unitigs.

| parameter | kind | default | options |
|---|---|---|---|
| `reads` | file | `examples/reads_single.fastq` | — |
| `k` | int | `15` | — |
| `min_coverage` | int | `1` | — |
| `min_length` | int | `40` | — |
| `head` | int | `200` | — |

Galaxy Tool Shed: `spades / velvet`

### greedy_overlap_assembly

**Greedy overlap-layout-consensus assembly** — Repeatedly join the read pair with the largest exact overlap to build contigs.

| parameter | kind | default | options |
|---|---|---|---|
| `reads` | file | `examples/reads_single.fastq` | — |
| `min_overlap` | int | `20` | — |
| `head` | int | `100` | — |

Galaxy Tool Shed: `cap3 / phrap`

### kmer_spectrum

**Read k-mer frequency spectrum** — Histogram of k-mer multiplicities - the error peak and genomic peak used to pick k.

| parameter | kind | default | options |
|---|---|---|---|
| `reads` | file | `examples/reads_single.fastq` | — |
| `k` | int | `17` | — |
| `max_count` | int | `20` | — |
| `canonical` | bool | `True` | — |

Galaxy Tool Shed: `kmergenie / jellyfish hist`

### scaffold_with_pairs

**Link contigs with paired reads** — Count read pairs bridging contig pairs to infer scaffold adjacency.

| parameter | kind | default | options |
|---|---|---|---|
| `reads1` | file | `examples/reads_1.fastq` | — |
| `reads2` | file | `examples/reads_2.fastq` | — |
| `contigs` | file | `examples/genome.fa` | — |
| `k` | int | `21` | — |
| `min_links` | int | `1` | — |

Galaxy Tool Shed: `scaffolds / OPERA`

## annotation

*Annotation* — in Galaxy group *Genomics Analysis*.

| Tool | What it does | Upstream Galaxy tool |
|---|---|---|
| [`annotate_genome_by_blast`](#annotate_genome_by_blast) | Blast proteins against the target gene set and transfer names, IDs and products. | `maker / blast2go annotation transfer` |
| [`bedgraph_per_gene_signal`](#bedgraph_per_gene_signal) | Aggregate a signal track over each gene, optionally GC-corrected. | `bedtools map / deepTools computeMatrix` |
| [`bedtools_genes_from_gff`](#bedtools_genes_from_gff) | Extract gene bodies (or TSSs) from a GTF/GFF3. | `iuc/gtf_sort` |
| [`bedtools_tss_list`](#bedtools_tss_list) | Promoter windows around every transcript start site. | `iuc/homer` |
| [`genome_predict_orfs`](#genome_predict_orfs) | Six-frame ORF finding with translation, GC and length statistics per ORF. | `prodigal / getorf (EMBOSS)` |
| [`gff3_add_locus_tags`](#gff3_add_locus_tags) | Assign sequential locus tags and keep Parent links consistent. | `Producer of Salk annotations` |
| [`gff3_biotype_table`](#gff3_biotype_table) | Genes per biotype with total length, mean length and chromosome distribution. | `biotype table (annotation)` |
| [`gff3_cds_translation`](#gff3_cds_translation) | Splice CDS exons per transcript, translate, and report length and stops. | `gff3read --prot / prodigal-like` |
| [`gff3_check`](#gff3_check) | Validate types, coordinates, IDs and Parent links, and report each problem found. | `gff3read / qualifyGFF` |
| [`gff3_children_table`](#gff3_children_table) | Flatten the gene → transcript → exon/CDS hierarchy into one row per child. | `AGAT extract / gff3read children` |
| [`gff3_compare`](#gff3_compare) | Report genes shared by two annotations and the positional Jaccard of all features. | `gff3compare / compare GFF` |
| [`gff3_feature_counts`](#gff3_feature_counts) | Feature tallies per group with total and mean feature length. | `gff3 read / featureCounts summary` |
| [`gff3_filter_attributes`](#gff3_filter_attributes) | Keep rows whose attribute matches a value or regular expression. | `gff3 read / AGAT filter` |
| [`gff3_gtf_convert`](#gff3_gtf_convert) | Swap the attribute syntax between RefSeq GTF and Ensembl-style GFF3. | `gffread / AGAT convert` |
| [`gff3_infer_transcripts`](#gff3_infer_transcripts) | Group exon rows by Parent and rebuild transcript models with CDS length. | `gff3read --exonmerge / BRAKER` |
| [`gff3_nearest_gene`](#gff3_nearest_gene) | Attach the closest gene name, overlap status and distance to every interval. | `bedtools closest` |
| [`gff3_summary`](#gff3_summary) | Number of genes, transcripts, exons, mean lengths and genome coverage. | `qualifyGFF /assembly summary` |
| [`gff3_to_bed`](#gff3_to_bed) | Rewrite annotation rows as genomic intervals, keeping exon blocks in BED12. | `gff3 to bed (bed12plus)` |
| [`gff3_to_genbank`](#gff3_to_genbank) | Emit LOCUS/FEATURES records with CDS translations from a GFF3 and genome. | `gff3 to genbank (BP)` |
| [`gff3_tss_bed`](#gff3_tss_bed) | Reduce a GFF3 to per-gene TSS intervals, optionally as promoter windows. | `tss from GFF3 (bedtools)` |
| [`peak_set_enrichment`](#peak_set_enrichment) | Fisher test of the genes near peaks against every gene set in a GMT file. | `HOMER2 findGEnrichment / g:Profiler` |
| [`predict_proteins_from_reads`](#predict_proteins_from_reads) | Concatenate reads into long contigs and call the longest protein per contig. | `metaWRAP / prodigal on contigs` |
| [`protein_function_summary`](#protein_function_summary) | Per-protein mass, pI, instability, transmembrane and PROSITE-like motif hits. | `InterProScan-lite / emboss pepstats` |
| [`vcf_annotated_consequence`](#vcf_annotated_consequence) | Assign transcript consequences and codon changes to variants. | `iuc/snpEff` |

### annotate_genome_by_blast

**Transfer protein annotation to a genome** — Blast proteins against the target gene set and transfer names, IDs and products.

| parameter | kind | default | options |
|---|---|---|---|
| `proteins` | file | `examples/proteins.faa` | — |
| `ref` | file | `examples/genes.fa` | — |
| `word_size` | int | `3` | — |
| `evalue` | number | `10.0` | — |
| `min_identity` | number | `60.0` | — |

Galaxy Tool Shed: `maker / blast2go annotation transfer`

### bedgraph_per_gene_signal

**Mean signal of a bedGraph inside features** — Aggregate a signal track over each gene, optionally GC-corrected.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/annotation.gff` | — |
| `graph` | file | `examples/coverage.bedgraph` | — |
| `stat` | choice | `mean` | mean, max, sum, median |
| `ref` | file | `examples/genome.fa` | — |
| `gc_normalise` | bool | `False` | — |
| `head` | int | `200` | — |

Galaxy Tool Shed: `bedtools map / deepTools computeMatrix`

### bedtools_genes_from_gff

**Get gene intervals from an annotation** — Extract gene bodies (or TSSs) from a GTF/GFF3.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/annotation.gff` | — |
| `biotype` | text | `` | — |
| `tss_only` | bool | `False` | — |
| `extend` | int | `0` | — |

Galaxy Tool Shed: `iuc/gtf_sort`

### bedtools_tss_list

**TSS list from annotation** — Promoter windows around every transcript start site.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/annotation.gff` | — |
| `upstream` | int | `1000` | — |
| `downstream` | int | `1000` | — |

Galaxy Tool Shed: `iuc/homer`

### genome_predict_orfs

**Gene prediction by open reading frames** — Six-frame ORF finding with translation, GC and length statistics per ORF.

| parameter | kind | default | options |
|---|---|---|---|
| `ref` | file | `examples/genes.fa` | — |
| `min_length` | int | `90` | — |
| `start` | choice | `ATG` | ATG, GTG, TTG |
| `partial` | bool | `False` | — |
| `head` | int | `200` | — |

Galaxy Tool Shed: `prodigal / getorf (EMBOSS)`

### gff3_add_locus_tags

**Renumber feature IDs with a locus tag prefix** — Assign sequential locus tags and keep Parent links consistent.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/annotation.gff` | — |
| `prefix` | text | `CHROMA_` | — |
| `start` | int | `1` | — |
| `padding` | int | `5` | — |
| `type` | choice | `gene` | gene, CDS, mRNA, any |

Galaxy Tool Shed: `Producer of Salk annotations`

### gff3_biotype_table

**Gene biotype census** — Genes per biotype with total length, mean length and chromosome distribution.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/annotation.gff` | — |
| `lengths` | bool | `True` | — |

Galaxy Tool Shed: `biotype table (annotation)`

### gff3_cds_translation

**Translate annotated CDS to proteins** — Splice CDS exons per transcript, translate, and report length and stops.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/annotation.gff` | — |
| `ref` | file | `examples/genome.fa` | — |
| `require_start` | bool | `True` | — |
| `min_length` | int | `20` | — |

Galaxy Tool Shed: `gff3read --prot / prodigal-like`

### gff3_check

**Sanity-check a GFF3 file** — Validate types, coordinates, IDs and Parent links, and report each problem found.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/annotation.gff` | — |
| `level` | choice | `all` | all, error, warning |

Galaxy Tool Shed: `gff3read / qualifyGFF`

### gff3_children_table

**Parent-child feature table** — Flatten the gene → transcript → exon/CDS hierarchy into one row per child.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/annotation.gff` | — |
| `from_type` | choice | `gene` | gene, transcript, mRNA |
| `head` | int | `300` | — |

Galaxy Tool Shed: `AGAT extract / gff3read children`

### gff3_compare

**Compare two annotation sets** — Report genes shared by two annotations and the positional Jaccard of all features.

| parameter | kind | default | options |
|---|---|---|---|
| `a` | file | `examples/annotation.gff` | — |
| `b` | file | `examples/annotation.gff` | — |
| `key` | choice | `name` | name, coordinates, id |

Galaxy Tool Shed: `gff3compare / compare GFF`

### gff3_feature_counts

**Count features by type and sequence** — Feature tallies per group with total and mean feature length.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/annotation.gff` | — |
| `group_by` | choice | `type` | type, seqid, biotype, source |

Galaxy Tool Shed: `gff3 read / featureCounts summary`

### gff3_filter_attributes

**Filter GFF rows by attributes** — Keep rows whose attribute matches a value or regular expression.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/annotation.gff` | — |
| `attribute` | text | `Name` | — |
| `value` | text | `` | — |
| `type` | choice | `any` | any, gene, mRNA, exon, CDS |
| `invert` | bool | `False` | — |

Galaxy Tool Shed: `gff3 read / AGAT filter`

### gff3_gtf_convert

**Convert between GFF3 and GTF** — Swap the attribute syntax between RefSeq GTF and Ensembl-style GFF3.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/annotation.gff` | — |
| `to` | choice | `gtf` | gtf, gff3 |

Galaxy Tool Shed: `gffread / AGAT convert`

### gff3_infer_transcripts

**Assemble transcripts from exon rows** — Group exon rows by Parent and rebuild transcript models with CDS length.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/annotation.gff` | — |
| `min_exons` | int | `1` | — |

Galaxy Tool Shed: `gff3read --exonmerge / BRAKER`

### gff3_nearest_gene

**Annotate intervals with the nearest gene** — Attach the closest gene name, overlap status and distance to every interval.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/regions.bed` | — |
| `annotation` | file | `examples/annotation.gff` | — |
| `max_distance` | number | `1000000000000.0` | — |
| `tie` | choice | `closest` | closest, first, all |

Galaxy Tool Shed: `bedtools closest`

### gff3_summary

**Annotation coverage statistics** — Number of genes, transcripts, exons, mean lengths and genome coverage.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/annotation.gff` | — |
| `sizes` | file | `examples/genome.fa` | — |

Galaxy Tool Shed: `qualifyGFF /assembly summary`

### gff3_to_bed

**Convert GFF3 to BED (BED12 for transcripts)** — Rewrite annotation rows as genomic intervals, keeping exon blocks in BED12.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/annotation.gff` | — |
| `type` | choice | `gene` | any, gene, mRNA, exon, CDS |
| `bed12` | bool | `True` | — |

Galaxy Tool Shed: `gff3 to bed (bed12plus)`

### gff3_to_genbank

**Write a GenBank-style flat file** — Emit LOCUS/FEATURES records with CDS translations from a GFF3 and genome.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/annotation.gff` | — |
| `ref` | file | `examples/genome.fa` | — |
| `definition` | text | `Chroma-Titan synthetic record` | — |

Galaxy Tool Shed: `gff3 to genbank (BP)`

### gff3_tss_bed

**Transcription start sites of all genes** — Reduce a GFF3 to per-gene TSS intervals, optionally as promoter windows.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/annotation.gff` | — |
| `window` | int | `0` | — |
| `upstream` | bool | `False` | — |

Galaxy Tool Shed: `tss from GFF3 (bedtools)`

### peak_set_enrichment

**Test peak-linked genes for pathway enrichment** — Fisher test of the genes near peaks against every gene set in a GMT file.

| parameter | kind | default | options |
|---|---|---|---|
| `peaks` | file | `examples/regions.bed` | — |
| `annotation` | file | `examples/annotation.gff` | — |
| `gmt` | file | `examples/genesets.gmt` | — |
| `min_q` | number | `1.0` | — |

Galaxy Tool Shed: `HOMER2 findGEnrichment / g:Profiler`

### predict_proteins_from_reads

**Assemble-like ORF calling from reads** — Concatenate reads into long contigs and call the longest protein per contig.

| parameter | kind | default | options |
|---|---|---|---|
| `reads` | file | `examples/reads_single.fastq` | — |
| `k` | int | `15` | — |
| `min_protein` | int | `30` | — |
| `head` | int | `100` | — |

Galaxy Tool Shed: `metaWRAP / prodigal on contigs`

### protein_function_summary

**Domain and property census of a proteome** — Per-protein mass, pI, instability, transmembrane and PROSITE-like motif hits.

| parameter | kind | default | options |
|---|---|---|---|
| `proteins` | file | `examples/proteins.faa` | — |
| `head` | int | `200` | — |

Galaxy Tool Shed: `InterProScan-lite / emboss pepstats`

### vcf_annotated_consequence

**snpEff-lite consequences** — Assign transcript consequences and codon changes to variants.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/variants.vcf` | — |
| `annotation` | file | `examples/annotation.gff` | — |
| `genome` | file | `examples/genome.fa` | — |
| `only_coding` | bool | `False` | — |

Galaxy Tool Shed: `iuc/snpEff`

## mapping

*Mapping* — in Galaxy group *Genomics Analysis*.

| Tool | What it does | Upstream Galaxy tool |
|---|---|---|
| [`bowtie2_end_to_end`](#bowtie2_end_to_end) | Full-length read alignment table with edit distance, MAPQ and per-read score. | `bowtie2 --end-to-end` |
| [`bwa_index_summary`](#bwa_index_summary) | Distinct k-mers, repeated k-mers and the most frequent seeds of an index. | `bwa index` |
| [`bwa_mem_lite`](#bwa_mem_lite) | Seed-and-extend mapping of FASTQ reads onto a FASTA reference, SAM output. | `bwa mem` |
| [`hisat2_spliced`](#hisat2_spliced) | Local mapping with clipped ends; clipped bases are reported as putative junctions. | `hisat2` |
| [`identity_histogram`](#identity_histogram) | Percent-identity histogram of every alignment in the file. | `samtools stats / plot` |
| [`insert_size_stats`](#insert_size_stats) | Median, mean, standard deviation and outlier counts of the fragment insert sizes. | `picard CollectsInsertSizeMetrics` |
| [`kmer_pseudoalign`](#kmer_pseudoalign) | Assign every read to the transcript sharing the most k-mers and tabulate counts. | `kallisto / salmon` |
| [`mapq_filter`](#mapq_filter) | Write a SAM containing only alignments above a MAPQ threshold. | `samtools view -q` |
| [`merge_pairs_lite`](#merge_pairs_lite) | Detect mate overlaps, merge into one read and report the merging rate. | `fastq-join / vsearch mergepairs` |
| [`samtools_stats_lite`](#samtools_stats_lite) | Read totals, mapping rate, mismatch and MAPQ means for an alignment file. | `samtools stats` |
| [`softclip_report`](#softclip_report) | List alignments whose clipped tails suggest chimeric reads or junctions. | `samtools view / clipping analysis` |
| [`unmapped_to_fastq`](#unmapped_to_fastq) | Recover unmapped sequences from an alignment as FASTQ for re-analysis. | `samtools fastq -f 4` |

### bowtie2_end_to_end

**End-to-end mapping with Bowtie2 scoring** — Full-length read alignment table with edit distance, MAPQ and per-read score.

| parameter | kind | default | options |
|---|---|---|---|
| `reads` | file | `examples/reads_single.fastq` | — |
| `ref` | file | `examples/genome.fa` | — |
| `k` | int | `31` | — |
| `max_mismatches` | int | `4` | — |
| `head` | int | `200` | — |

Galaxy Tool Shed: `bowtie2 --end-to-end`

### bwa_index_summary

**Summarise a k-mer index of the reference** — Distinct k-mers, repeated k-mers and the most frequent seeds of an index.

| parameter | kind | default | options |
|---|---|---|---|
| `ref` | file | `examples/genome.fa` | — |
| `k` | int | `11` | — |
| `top` | int | `20` | — |

Galaxy Tool Shed: `bwa index`

### bwa_mem_lite

**Map reads to a reference (BWA-MEM style)** — Seed-and-extend mapping of FASTQ reads onto a FASTA reference, SAM output.

| parameter | kind | default | options |
|---|---|---|---|
| `reads` | file | `examples/reads_single.fastq` | — |
| `ref` | file | `examples/genome.fa` | — |
| `k` | int | `21` | — |
| `max_mismatches` | int | `3` | — |
| `seed_stride` | int | `1` | — |
| `min_mapq` | int | `0` | — |
| `local` | bool | `False` | — |

Galaxy Tool Shed: `bwa mem`

### hisat2_spliced

**Spliced / soft-clipped mapping (HISAT2 style)** — Local mapping with clipped ends; clipped bases are reported as putative junctions.

| parameter | kind | default | options |
|---|---|---|---|
| `reads` | file | `examples/reads_1.fastq` | — |
| `ref` | file | `examples/genes.fa` | — |
| `k` | int | `17` | — |
| `max_mismatches` | int | `5` | — |
| `head` | int | `200` | — |

Galaxy Tool Shed: `hisat2`

### identity_histogram

**Distribution of alignment identity** — Percent-identity histogram of every alignment in the file.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/alignments.sam` | — |
| `bin_size` | int | `5` | — |
| `source` | choice | `sam` | sam, reads |

Galaxy Tool Shed: `samtools stats / plot`

### insert_size_stats

**Insert size (TLEN) distribution** — Median, mean, standard deviation and outlier counts of the fragment insert sizes.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/alignments.sam` | — |
| `bins` | int | `10` | — |
| `sd_cut` | int | `3.0` | — |

Galaxy Tool Shed: `picard CollectsInsertSizeMetrics`

### kmer_pseudoalign

**k-mer pseudoalignment to transcripts** — Assign every read to the transcript sharing the most k-mers and tabulate counts.

| parameter | kind | default | options |
|---|---|---|---|
| `reads` | file | `examples/reads_1.fastq` | — |
| `transcripts` | file | `examples/genes.fa` | — |
| `k` | int | `21` | — |
| `min_frac` | number | `0.4` | — |
| `head` | int | `100` | — |

Galaxy Tool Shed: `kallisto / salmon`

### mapq_filter

**Filter alignments by MAPQ** — Write a SAM containing only alignments above a MAPQ threshold.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/alignments.sam` | — |
| `min_mapq` | int | `30` | — |
| `drop_unmapped` | bool | `False` | — |

Galaxy Tool Shed: `samtools view -q`

### merge_pairs_lite

**Merge overlapping read pairs** — Detect mate overlaps, merge into one read and report the merging rate.

| parameter | kind | default | options |
|---|---|---|---|
| `reads1` | file | `examples/reads_1.fastq` | — |
| `reads2` | file | `examples/reads_2.fastq` | — |
| `min_overlap` | int | `11` | — |
| `max_mismatches` | int | `2` | — |
| `head` | int | `100` | — |

Galaxy Tool Shed: `fastq-join / vsearch mergepairs`

### samtools_stats_lite

**Mapping statistics of a SAM file** — Read totals, mapping rate, mismatch and MAPQ means for an alignment file.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/alignments.sam` | — |

Galaxy Tool Shed: `samtools stats`

### softclip_report

**Soft-clipped bases per alignment** — List alignments whose clipped tails suggest chimeric reads or junctions.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/alignments.sam` | — |
| `min_clip` | int | `10` | — |
| `head` | int | `100` | — |

Galaxy Tool Shed: `samtools view / clipping analysis`

### unmapped_to_fastq

**Extract unmapped reads** — Recover unmapped sequences from an alignment as FASTQ for re-analysis.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/alignments.sam` | — |
| `singletons` | bool | `True` | — |

Galaxy Tool Shed: `samtools fastq -f 4`

## variant_calling

*Variant Calling* — in Galaxy group *Genomics Analysis*.

| Tool | What it does | Upstream Galaxy tool |
|---|---|---|
| [`vcf_apply_vqsr`](#vcf_apply_vqsr) | Approximate VQSR: partition variants by QUAL/DP/GQ into PASS and LOW tiers. | `iuc/gatk` |
| [`vcf_filter_records`](#vcf_filter_records) | Apply one or more JEXL-ish filters, tagging or dropping records. | `iuc/gatk` |
| [`vcf_genotype_refinement`](#vcf_genotype_refinement) | Re-genotype samples by allelic depth with a ploidy-aware model. | `iuc/bcftools` |
| [`vcf_haploidize`](#vcf_haploidize) | Collapse diploid GT calls to haplotypes (first allele). | `iuc/bcftools` |
| [`vcf_joint_call`](#vcf_joint_call) | Genotype all samples at the union of discovered sites. | `iuc/gatk` |
| [`vcf_pileup_call`](#vcf_pileup_call) | FreeBayes-lite pileup caller returning a variant table. | `iuc/freebayes` |
| [`vcf_somatic_pair`](#vcf_somatic_pair) | MuTect-style difference of two genotype columns. | `iuc/mutect2` |

### vcf_apply_vqsr

**Apply VQSR-lite (score bins)** — Approximate VQSR: partition variants by QUAL/DP/GQ into PASS and LOW tiers.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/variants.vcf` | — |
| `truth_sensitivity` | number | `0.9` | — |
| `min_gq` | int | `20` | — |
| `min_dp` | int | `3` | — |

Galaxy Tool Shed: `iuc/gatk`

### vcf_filter_records

**Variant Filtration (GATK style)** — Apply one or more JEXL-ish filters, tagging or dropping records.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/variants.vcf` | — |
| `expressions` | text | `QUAL<30 || DP<5` | — |
| `invert` | bool | `False` | — |
| `set_filter_tag` | bool | `True` | — |
| `filter_name` | text | `LowQual` | — |

Galaxy Tool Shed: `iuc/gatk`

### vcf_genotype_refinement

**Refine genotypes from AD fields** — Re-genotype samples by allelic depth with a ploidy-aware model.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/variants.vcf` | — |
| `het_bias` | number | `0.5` | — |
| `min_gq` | number | `20.0` | — |
| `nocall_low_dp` | bool | `True` | — |
| `min_dp` | int | `3` | — |

Galaxy Tool Shed: `iuc/bcftools`

### vcf_haploidize

**Force haploid genotypes** — Collapse diploid GT calls to haplotypes (first allele).

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/variants.vcf` | — |

Galaxy Tool Shed: `iuc/bcftools`

### vcf_joint_call

**Joint calling across samples** — Genotype all samples at the union of discovered sites.

| parameter | kind | default | options |
|---|---|---|---|
| `files` | multi | `['examples/alignments.sam']` | examples/alignments.sam |
| `genome` | file | `examples/genome.fa` | — |
| `min_depth` | int | `3` | — |
| `min_vaf` | number | `0.2` | — |
| `hwe_prior` | bool | `True` | — |

Galaxy Tool Shed: `iuc/gatk`

### vcf_pileup_call

**Call variants from a pileup table** — FreeBayes-lite pileup caller returning a variant table.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/alignments.sam` | — |
| `genome` | file | `examples/genome.fa` | — |
| `region` | text | `` | — |
| `min_depth` | int | `5` | — |
| `min_vaf` | number | `0.2` | — |
| `min_mapq` | int | `20` | — |
| `min_bq` | int | `13` | — |

Galaxy Tool Shed: `iuc/freebayes`

### vcf_somatic_pair

**Somatic calling from a tumour/normal pair** — MuTect-style difference of two genotype columns.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/variants.vcf` | — |
| `tumour` | text | `` | — |
| `normal` | text | `` | — |
| `min_tumor_alt` | int | `5` | — |
| `min_vaf` | number | `0.05` | — |
| `max_normal_alt` | int | `1` | — |
| `max_normal_vaf` | number | `0.05` | — |

Galaxy Tool Shed: `iuc/mutect2`

## chip_seq

*Chromatin Accessibility (ChIP-seq/ATAC)* — in Galaxy group *Genomics Analysis*.

| Tool | What it does | Upstream Galaxy tool |
|---|---|---|
| [`bedtools_metagene_profile`](#bedtools_metagene_profile) | Density of A intervals around TSSs in fixed bins. | `iuc/deepTools` |
| [`chipseq_qc_metrics`](#chipseq_qc_metrics) | Read counts, duplication rate, fragments in peaks and cross-correlation proxy. | `phantompeakqualtools / deepTools plotFingerprint` |
| [`macs2_callpeak_lite`](#macs2_callpeak_lite) | Threshold and merge signal blocks into peaks with fold-enrichment and summit scores. | `MACS2 callpeak` |
| [`peak_gene_annotation`](#peak_gene_annotation) | Classify each peak as promoter, intron, exon or intergenic with distance to the TSS. | `HOMER2 annotatePeaks / ChIPAnno` |
| [`peak_grouper`](#peak_grouper) | Cluster peaks into loci, keeping member names, counts and summed signal. | `bedtools merge / MACS2 peak grouping` |

### bedtools_metagene_profile

**Metagene TSS profile of intervals** — Density of A intervals around TSSs in fixed bins.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/regions.bed` | — |
| `annotation` | file | `examples/annotation.gff` | — |
| `bin_size` | int | `100` | — |
| `flank` | int | `2000` | — |

Galaxy Tool Shed: `iuc/deepTools`

### chipseq_qc_metrics

**Library QC: duplication, FRiP and strand metrics** — Read counts, duplication rate, fragments in peaks and cross-correlation proxy.

| parameter | kind | default | options |
|---|---|---|---|
| `alignments` | file | `examples/alignments.sam` | — |
| `peaks` | file | `examples/regions.bed` | — |
| `read_length` | number | `50` | — |

Galaxy Tool Shed: `phantompeakqualtools / deepTools plotFingerprint`

### macs2_callpeak_lite

**Call peaks from a signal track** — Threshold and merge signal blocks into peaks with fold-enrichment and summit scores.

| parameter | kind | default | options |
|---|---|---|---|
| `graph` | file | `examples/coverage.bedgraph` | — |
| `cutoff` | number | `5.0` | — |
| `min_length` | int | `50` | — |
| `max_gap` | int | `100` | — |
| `broad` | bool | `False` | — |
| `control` | file | `` | — |

Galaxy Tool Shed: `MACS2 callpeak`

### peak_gene_annotation

**Annotate peaks to the nearest TSS** — Classify each peak as promoter, intron, exon or intergenic with distance to the TSS.

| parameter | kind | default | options |
|---|---|---|---|
| `peaks` | file | `examples/regions.bed` | — |
| `annotation` | file | `examples/annotation.gff` | — |
| `tss_distance` | number | `2000` | — |
| `head` | int | `200` | — |

Galaxy Tool Shed: `HOMER2 annotatePeaks / ChIPAnno`

### peak_grouper

**Merge and group overlapping peaks** — Cluster peaks into loci, keeping member names, counts and summed signal.

| parameter | kind | default | options |
|---|---|---|---|
| `peaks` | file | `examples/regions.bed` | — |
| `distance` | int | `0` | — |
| `collapse` | choice | `concat` | concat, sum, max, first |
| `head` | int | `200` | — |

Galaxy Tool Shed: `bedtools merge / MACS2 peak grouping`

## rna_seq

*RNA Sequencing (RNA-Seq)* — in Galaxy group *Genomics Analysis*.

| Tool | What it does | Upstream Galaxy tool |
|---|---|---|
| [`deseq_lite`](#deseq_lite) | Per-gene log2 fold change, test statistic, raw and adjusted p-values. | `DESeq2 / edgeR` |
| [`differential_fisher_on_off`](#differential_fisher_on_off) | Test on/off patterns between groups when only detection matters. | `Fisher exact on/off (edgeR-ancestral)` |
| [`gene_set_counts`](#gene_set_counts) | Score each pathway by summarising the expression of its member genes. | `GSVA / roasters` |
| [`gsea_lite`](#gsea_lite) | Running-sum enrichment score with a permutation null for every gene set. | `GSEA / fgsea` |
| [`limma_moderated_t`](#limma_moderated_t) | Shrink the gene-wise variance toward a pooled estimate and recompute t statistics. | `limma eBayes` |
| [`ma_plot`](#ma_plot) | Detect intensity-dependent bias by plotting fold change against abundance. | `edgeR plotMD / limma` |
| [`permutation_group_test`](#permutation_group_test) | Distribution-based p-values for group contrasts when sample numbers are tiny. | `permutation test / signal` |
| [`volcano_plot`](#volcano_plot) | Plot log2 fold change against -log10 p-value with significance thresholds. | `EnhancedVolcano / DESeq2 results plot` |

### deseq_lite

**Differential expression (negative-binomial-free DESeq2 style)** — Per-gene log2 fold change, test statistic, raw and adjusted p-values.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/counts.tsv` | — |
| `groups` | text | `` | — |
| `phenotype_file` | file | `` | — |
| `test` | choice | `wald` | wald, welch, exact |
| `adjust` | choice | `fdr_bh` | fdr_bh, bonferroni, none |
| `alpha` | number | `0.05` | — |

Galaxy Tool Shed: `DESeq2 / edgeR`

### differential_fisher_on_off

**Fisher test of presence/absence per feature** — Test on/off patterns between groups when only detection matters.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/counts.tsv` | — |
| `groups` | text | `` | — |
| `cutoff` | int | `1` | — |

Galaxy Tool Shed: `Fisher exact on/off (edgeR-ancestral)`

### gene_set_counts

**Aggregate expression over GMT gene sets** — Score each pathway by summarising the expression of its member genes.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/counts.tsv` | — |
| `sets` | file | `examples/genesets.gmt` | — |
| `agg` | choice | `mean` | mean, sum, median, max |
| `zscore` | bool | `True` | — |

Galaxy Tool Shed: `GSVA / roasters`

### gsea_lite

**Rank-based GSEA of a gene set** — Running-sum enrichment score with a permutation null for every gene set.

| parameter | kind | default | options |
|---|---|---|---|
| `results` | file | `examples/counts.tsv` | — |
| `rank_column` | text | `` | — |
| `sets` | file | `examples/genesets.gmt` | — |
| `permutations` | int | `200` | — |
| `min_size` | number | `2` | — |

Galaxy Tool Shed: `GSEA / fgsea`

### limma_moderated_t

**limma-style moderated t-test** — Shrink the gene-wise variance toward a pooled estimate and recompute t statistics.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/counts.tsv` | — |
| `groups` | text | `` | — |
| `prior_df` | number | `3.0` | — |
| `alpha` | number | `0.05` | — |

Galaxy Tool Shed: `limma eBayes`

### ma_plot

**MA plot (fold change versus mean expression)** — Detect intensity-dependent bias by plotting fold change against abundance.

| parameter | kind | default | options |
|---|---|---|---|
| `results` | file | `examples/counts.tsv` | — |
| `mean_column` | text | `` | — |
| `log2fc_column` | text | `` | — |

Galaxy Tool Shed: `edgeR plotMD / limma`

### permutation_group_test

**Permutation test between two sample groups** — Distribution-based p-values for group contrasts when sample numbers are tiny.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/counts.tsv` | — |
| `groups` | text | `` | — |
| `permutations` | int | `500` | — |
| `stat` | choice | `mean_difference` | mean_difference, t, sum |
| `seed` | int | `1` | — |

Galaxy Tool Shed: `permutation test / signal`

### volcano_plot

**Volcano plot of differential results** — Plot log2 fold change against -log10 p-value with significance thresholds.

| parameter | kind | default | options |
|---|---|---|---|
| `results` | file | `examples/counts.tsv` | — |
| `log2fc_column` | text | `` | — |
| `p_column` | text | `` | — |
| `lfc_cut` | number | `1.0` | — |
| `p_cut` | number | `0.05` | — |
| `label_column` | text | `` | — |

Galaxy Tool Shed: `EnhancedVolcano / DESeq2 results plot`

## multiple_alignments

*Multiple Alignments* — in Galaxy group *Genomics Analysis*.

| Tool | What it does | Upstream Galaxy tool |
|---|---|---|
| [`clustal_omega_align`](#clustal_omega_align) | Alignment guided by a distance tree, output in Clustal-style columns. | `clustalo` |
| [`edit_distance_report`](#edit_distance_report) | Levenshtein and Hamming distances for every pair of sequences in a file. | `emboss distmat / edlib` |
| [`mafft_align`](#mafft_align) | Center-star progressive MSA with optional profile refinement. | `mafft` |
| [`msa_back_translate`](#msa_back_translate) | Re-insert codons into a protein alignment to obtain a codon-based alignment. | `pal2nal / rev-trans` |
| [`msa_column_frequencies`](#msa_column_frequencies) | Frequency table of every residue at every alignment column. | `seqkit stats on alignment` |
| [`msa_consensus`](#msa_consensus) | Collapse the alignment to one sequence using majority votes or IUPAC codes. | `emboss cons / samtools consensus` |
| [`msa_entropy_profile`](#msa_entropy_profile) | Shannon entropy per alignment column: conserved core versus variable loops. | `plot of per-column entropy` |
| [`msa_format_convert`](#msa_format_convert) | Re-serialise an alignment in the format expected by phylogenetics or HMM tools. | `alignment conversion` |
| [`msa_stats`](#msa_stats) | Conserved, variable, gappy and gap-only columns plus mean pairwise identity. | `aliview / msa summary` |
| [`msa_trim_columns`](#msa_trim_columns) | Remove poorly aligned columns before building a tree. | `trimal` |
| [`msa_variable_sites`](#msa_variable_sites) | Columns that vary between sequences, with entropy and taxon pattern. | `variable sites / seaview` |
| [`muscle_align`](#muscle_align) | Progressive alignment with refinement cycles, as muscle -refine does. | `muscle` |
| [`pairwise_align`](#pairwise_align) | Dynamic-programming alignment of two sequences with affine gaps and a matrix. | `needle / water (EMBOSS)` |
| [`pairwise_identity_matrix`](#pairwise_identity_matrix) | N x N matrix of pairwise identity or substitution-corrected distance. | `peptide-database / emboss proteity` |
| [`pam_kimura_distances`](#pam_kimura_distances) | Pairwise substitution distances corrected for multiple hits. | `phylip proml / distmat` |

### clustal_omega_align

**Multiple sequence alignment (Clustal Omega)** — Alignment guided by a distance tree, output in Clustal-style columns.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/msa.fasta` | — |
| `guide_tree` | choice | `upgma` | upgma, nj |

Galaxy Tool Shed: `clustalo`

### edit_distance_report

**Edit distance and Hamming distance between sequences** — Levenshtein and Hamming distances for every pair of sequences in a file.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/msa.fasta` | — |
| `head` | int | `200` | — |
| `upper_only` | bool | `True` | — |

Galaxy Tool Shed: `emboss distmat / edlib`

### mafft_align

**Multiple sequence alignment (MAFFT)** — Center-star progressive MSA with optional profile refinement.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/msa.fasta` | — |
| `algorithm` | choice | `auto` | auto, fft-ns-1, linsi, refine |
| `gap_open` | number | `-1.5` | — |
| `gap_extend` | number | `-0.5` | — |
| `maxiterate` | int | `0` | — |

Galaxy Tool Shed: `mafft`

### msa_back_translate

**Protein alignment to codon alignment** — Re-insert codons into a protein alignment to obtain a codon-based alignment.

| parameter | kind | default | options |
|---|---|---|---|
| `protein_alignment` | file | `examples/proteins.faa` | — |
| `nucleotides` | file | `examples/genes.fa` | — |

Galaxy Tool Shed: `pal2nal / rev-trans`

### msa_column_frequencies

**Per-column residue frequencies** — Frequency table of every residue at every alignment column.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/msa.fasta` | — |
| `start` | int | `1` | — |
| `end` | int | `0` | — |

Galaxy Tool Shed: `seqkit stats on alignment`

### msa_consensus

**Consensus sequence of an alignment** — Collapse the alignment to one sequence using majority votes or IUPAC codes.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/msa.fasta` | — |
| `mode` | choice | `majority` | majority, IUPAC, first |
| `threshold` | number | `0.5` | — |

Galaxy Tool Shed: `emboss cons / samtools consensus`

### msa_entropy_profile

**Conservation profile of an alignment** — Shannon entropy per alignment column: conserved core versus variable loops.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/msa.fasta` | — |
| `style` | choice | `line` | line, bar |
| `threshold` | number | `1.0` | — |

Galaxy Tool Shed: `plot of per-column entropy`

### msa_format_convert

**Convert alignment to Clustal / Stockholm / PHYLIP** — Re-serialise an alignment in the format expected by phylogenetics or HMM tools.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/msa.fasta` | — |
| `format` | choice | `clustal` | clustal, stockholm, phylip, fasta |

Galaxy Tool Shed: `alignment conversion`

### msa_stats

**Alignment column statistics** — Conserved, variable, gappy and gap-only columns plus mean pairwise identity.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/msa.fasta` | — |
| `gap_threshold` | int | `50` | — |

Galaxy Tool Shed: `aliview / msa summary`

### msa_trim_columns

**Trim gappy alignment columns (trimAl)** — Remove poorly aligned columns before building a tree.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/msa.fasta` | — |
| `mode` | choice | `gappyout` | gappyout, strict, automated1 |
| `threshold` | number | `0.5` | — |

Galaxy Tool Shed: `trimal`

### msa_variable_sites

**Variable columns of an alignment** — Columns that vary between sequences, with entropy and taxon pattern.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/msa.fasta` | — |
| `head` | int | `300` | — |
| `informative_only` | bool | `False` | — |

Galaxy Tool Shed: `variable sites / seaview`

### muscle_align

**Multiple sequence alignment (MUSCLE)** — Progressive alignment with refinement cycles, as muscle -refine does.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/msa.fasta` | — |
| `refinement` | int | `1` | — |

Galaxy Tool Shed: `muscle`

### pairwise_align

**Pairwise alignment (Needleman-Wunsch / Smith-Waterman)** — Dynamic-programming alignment of two sequences with affine gaps and a matrix.

| parameter | kind | default | options |
|---|---|---|---|
| `seq_a` | file | `examples/genes.fa` | — |
| `seq_b` | file | `examples/genes.fa` | — |
| `mode` | choice | `global` | global, local |
| `matrix` | choice | `` | , BLOSUM62, PAM250 |
| `gap_open` | number | `-1.5` | — |
| `gap_extend` | number | `-0.5` | — |
| `match` | int | `2` | — |
| `mismatch` | int | `-1` | — |

Galaxy Tool Shed: `needle / water (EMBOSS)`

### pairwise_identity_matrix

**Pairwise identity / distance matrix** — N x N matrix of pairwise identity or substitution-corrected distance.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/proteins.faa` | — |
| `model` | choice | `identity` | identity, distance, blosum |
| `heatmap` | bool | `True` | — |

Galaxy Tool Shed: `peptide-database / emboss proteity`

### pam_kimura_distances

**Substitution-corrected distances** — Pairwise substitution distances corrected for multiple hits.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/msa.fasta` | — |
| `model` | choice | `pam` | pam, kimura, uncorrected |
| `max_distance` | number | `3.0` | — |

Galaxy Tool Shed: `phylip proml / distmat`

## phenotype_association

*Phenotype Association* — in Galaxy group *Genomics Analysis*.

| Tool | What it does | Upstream Galaxy tool |
|---|---|---|
| [`vcf_fisher_case_control`](#vcf_fisher_case_control) | Allelic chi-square/Fisher test of each variant against a case-control phenotype. | `iuc/plink` |
| [`vcf_linear_association`](#vcf_linear_association) | Allele-dosage regression on a quantitative trait per variant. | `iuc/plink` |
| [`vcf_manhattan`](#vcf_manhattan) | Manhattan plot of association p-values in a VCF. | `iuc/qqman` |

### vcf_fisher_case_control

**Fisher exact association test** — Allelic chi-square/Fisher test of each variant against a case-control phenotype.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/variants.vcf` | — |
| `phenotypes` | file | `examples/phenotypes.tsv` | — |
| `case_label` | text | `case` | — |
| `min_maf` | number | `0.01` | — |
| `additive` | bool | `True` | — |

Galaxy Tool Shed: `iuc/plink`

### vcf_linear_association

**Quantitative association (linear regression)** — Allele-dosage regression on a quantitative trait per variant.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/variants.vcf` | — |
| `phenotypes` | file | `examples/phenotypes.tsv` | — |
| `trait_column` | text | `3` | — |
| `per_allele` | bool | `True` | — |
| `covariate_mean` | number | `0.0` | — |

Galaxy Tool Shed: `iuc/plink`

### vcf_manhattan

**Manhattan plot from a VCF** — Manhattan plot of association p-values in a VCF.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/variants.vcf` | — |
| `p_from` | choice | `QUAL` | QUAL, INFO:P, FORMAT:PV |
| `label_top` | bool | `True` | — |

Galaxy Tool Shed: `iuc/qqman`

## phylogenetics

*Phylogenetics* — in Galaxy group *Genomics Analysis*.

| Tool | What it does | Upstream Galaxy tool |
|---|---|---|
| [`phylo_bootstrap`](#phylo_bootstrap) | Resample alignment columns and annotate the consensus tree with support values. | `phangorn bootstrap / seqboot` |
| [`phylo_consensus_tree`](#phylo_consensus_tree) | Collapse a set of trees into one topology with clade frequencies. | `consense / APE consensus` |
| [`phylo_dating_check`](#phylo_dating_check) | Test whether root-to-tip divergence correlates with sampling time. | `TempEst / liTo root-to-tip regression` |
| [`phylo_distance_from_tree`](#phylo_distance_from_tree) | Pairwise tip distances implied by a tree, with an ultrametricity test. | `ape::cophenetic.phylo / diagonalize` |
| [`phylo_nj_from_alignment`](#phylo_nj_from_alignment) | Distance matrix plus neighbour-joining tree, written as Newick. | `FastTree / MEGA NJ / ape NJ` |
| [`phylo_pairwise_distances`](#phylo_pairwise_distances) | Uncorrected and model-corrected pairwise distances for every pair. | `ape::dist.dna / phylip protdist` |
| [`phylo_read_nexus`](#phylo_read_nexus) | Parse the trees stored in a NEXUS block and summarise them. | `phytools read.nexus / DendroPy` |
| [`phylo_robinson_foulds`](#phylo_robinson_foulds) | Count the splits that differ between two topologies. | `ape::RF / treedist` |
| [`phylo_simulate_evolution`](#phylo_simulate_evolution) | Generate an alignment by evolving a sequence down a topology. | `evolver / Seq-Gen` |
| [`phylo_tree_operations`](#phylo_tree_operations) | Apply a single topological operation and return the modified tree. | `ape root/collapse / ete3 prune` |
| [`phylo_tree_statistics`](#phylo_tree_statistics) | Tips, clades, tree length and a balance measure for each input tree. | `ape::NNI stats / treebalance` |
| [`phylo_upgma`](#phylo_upgma) | Ultrametric phenetic tree from an alignment distance matrix. | `ape::upgma / MEGA UPGMA` |

### phylo_bootstrap

**Non-parametric bootstrap support for a tree** — Resample alignment columns and annotate the consensus tree with support values.

| parameter | kind | default | options |
|---|---|---|---|
| `alignment_file` | file | `examples/msa.fasta` | — |
| `replicates` | int | `100` | — |
| `method` | choice | `nj` | nj, upgma |
| `seed` | int | `1` | — |
| `model` | choice | `p` | p, identity, k2p |

Galaxy Tool Shed: `phangorn bootstrap / seqboot`

### phylo_consensus_tree

**Majority-rule consensus of a tree set** — Collapse a set of trees into one topology with clade frequencies.

| parameter | kind | default | options |
|---|---|---|---|
| `trees` | file | `` | — |
| `cutoff` | number | `0.5` | — |
| `support_values` | bool | `True` | — |

Galaxy Tool Shed: `consense / APE consensus`

### phylo_dating_check

**Root-to-tip distances against sampling order** — Test whether root-to-tip divergence correlates with sampling time.

| parameter | kind | default | options |
|---|---|---|---|
| `newick` | file | `examples/tree.nwk` | — |
| `dates` | file | `examples/phenotypes.tsv` | — |
| `taxon_column` | text | `` | — |
| `date_column` | text | `` | — |
| `regression` | choice | `spearman` | pearson, spearman |

Galaxy Tool Shed: `TempEst / liTo root-to-tip regression`

### phylo_distance_from_tree

**Patristic and cophenetic distance matrices** — Pairwise tip distances implied by a tree, with an ultrametricity test.

| parameter | kind | default | options |
|---|---|---|---|
| `newick` | file | `examples/tree.nwk` | — |
| `measure` | choice | `patristic` | patristic, cophenetic |
| `heatmap` | bool | `True` | — |
| `ultrametric_cutoff` | number | `0.05` | — |

Galaxy Tool Shed: `ape::cophenetic.phylo / diagonalize`

### phylo_nj_from_alignment

**Neighbour-joining tree from an alignment** — Distance matrix plus neighbour-joining tree, written as Newick.

| parameter | kind | default | options |
|---|---|---|---|
| `alignment_file` | file | `examples/msa.fasta` | — |
| `model` | choice | `p` | p, identity, k2p, jc69, tn93, jaccard_kmer |
| `root_midpoint` | bool | `True` | — |
| `min_length` | int | `4` | — |

Galaxy Tool Shed: `FastTree / MEGA NJ / ape NJ`

### phylo_pairwise_distances

**Pairwise distances between sequences** — Uncorrected and model-corrected pairwise distances for every pair.

| parameter | kind | default | options |
|---|---|---|---|
| `alignment_file` | file | `examples/msa.fasta` | — |
| `model` | choice | `p` | p, identity, k2p, jc69, tn93, jaccard_kmer |
| `matrix` | bool | `True` | — |
| `head` | int | `60` | — |

Galaxy Tool Shed: `ape::dist.dna / phylip protdist`

### phylo_read_nexus

**Read a multi-tree NEXUS file** — Parse the trees stored in a NEXUS block and summarise them.

| parameter | kind | default | options |
|---|---|---|---|
| `nexus` | file | `#NEXUS
begin trees;
tree t1 = ((a:0.1,b:0.2)c:0.3,d:0.4);
tree t2 = ((a:0.1,b:0.3)c:0.2,d:0.4);
end;` | — |
| `consensus` | bool | `True` | — |

Galaxy Tool Shed: `phytools read.nexus / DendroPy`

### phylo_robinson_foulds

**Robinson-Foulds distance between two trees** — Count the splits that differ between two topologies.

| parameter | kind | default | options |
|---|---|---|---|
| `tree_a` | text | `((seq_1,seq_2),(seq_3,seq_4));` | — |
| `tree_b` | text | `((seq_1,seq_3),(seq_2,seq_4));` | — |
| `normalised` | bool | `True` | — |

Galaxy Tool Shed: `ape::RF / treedist`

### phylo_simulate_evolution

**Simulate sequences along a tree** — Generate an alignment by evolving a sequence down a topology.

| parameter | kind | default | options |
|---|---|---|---|
| `newick` | file | `examples/tree.nwk` | — |
| `length` | int | `300` | — |
| `rate` | number | `1.0` | — |
| `seed` | int | `4` | — |
| `random_tree_mode` | choice | `yule` | yule, pda |
| `ntaxa` | int | `4` | — |

Galaxy Tool Shed: `evolver / Seq-Gen`

### phylo_tree_operations

**Root, ladderize, prune or subset a tree** — Apply a single topological operation and return the modified tree.

| parameter | kind | default | options |
|---|---|---|---|
| `newick` | file | `examples/tree.nwk` | — |
| `operation` | choice | `midpoint_root` | midpoint_root, reroot_leaf, ladderize, prune, subtree, strip_lengths |
| `leaf` | text | `` | — |
| `keep` | text | `` | — |
| `report_stats` | bool | `True` | — |

Galaxy Tool Shed: `ape root/collapse / ete3 prune`

### phylo_tree_statistics

**Topology statistics of one or more trees** — Tips, clades, tree length and a balance measure for each input tree.

| parameter | kind | default | options |
|---|---|---|---|
| `trees` | file | `examples/tree.nwk` | — |
| `splits` | bool | `True` | — |
| `colless` | bool | `True` | — |

Galaxy Tool Shed: `ape::NNI stats / treebalance`

### phylo_upgma

**UPGMA / WPGMA distance clustering tree** — Ultrametric phenetic tree from an alignment distance matrix.

| parameter | kind | default | options |
|---|---|---|---|
| `alignment_file` | file | `examples/msa.fasta` | — |
| `weighted` | bool | `False` | — |

Galaxy Tool Shed: `ape::upgma / MEGA UPGMA`

## evolution

*Evolution* — in Galaxy group *Genomics Analysis*.

| Tool | What it does | Upstream Galaxy tool |
|---|---|---|
| [`ancestral_parsimony`](#ancestral_parsimony) | Reconstruct internal nucleotides by Fitch up/down pass and count changes. | `phyml ancestral / Mesquite` |
| [`codon_enc_number`](#codon_enc_number) | Per-gene Nc, GC3 and codon bias, low Nc meaning strong bias. | `EMBOSS cusp / codon usage (CAI)` |
| [`codon_optimise`](#codon_optimise) | Reconstruct a codon-optimised CDS from a protein sequence. | `EMBOSS backtranse / genSmart` |
| [`codon_usage_table`](#codon_usage_table) | Counts per codon with amino acid, RSCU and the per-usage table of a CDS set. | `CAI / codonW / EMBOSS cusp` |
| [`dnds_pairwise`](#dnds_pairwise) | Synonymous and non-synonymous rates with the Jukes-Cantor corrected omega ratio. | `codeml / yn00 (PAML)` |
| [`dnds_windows`](#dnds_windows) | Locate stretches of elevated or reduced selective pressure along the alignment. | `slidingWindow / codeml windows` |
| [`gc3_profile`](#gc3_profile) | Sliding-window GC (or GC3) profile, the classic isoform/isochores diagnostic. | `EMBOSS geece / GCplot` |
| [`tree_from_alignment_nj`](#tree_from_alignment_nj) | Distance tree with optional bootstrap support values, written as Newick. | `fasttree / phyml / iqtree` |

### ancestral_parsimony

**Fitch parsimony ancestral reconstruction** — Reconstruct internal nucleotides by Fitch up/down pass and count changes.

| parameter | kind | default | options |
|---|---|---|---|
| `alignment` | file | `examples/msa.fasta` | — |
| `tree` | file | `examples/tree.nwk` | — |
| `head` | int | `100` | — |

Galaxy Tool Shed: `phyml ancestral / Mesquite`

### codon_enc_number

**Effective number of codons (Nc) and CAI** — Per-gene Nc, GC3 and codon bias, low Nc meaning strong bias.

| parameter | kind | default | options |
|---|---|---|---|
| `cds` | file | `examples/genes.fa` | — |
| `reference` | file | `` | — |
| `head` | int | `200` | — |

Galaxy Tool Shed: `EMBOSS cusp / codon usage (CAI)`

### codon_optimise

**Back-translate a protein with a host codon table** — Reconstruct a codon-optimised CDS from a protein sequence.

| parameter | kind | default | options |
|---|---|---|---|
| `protein` | file | `examples/proteins.faa` | — |
| `organism` | choice | `E. coli` | E. coli, human, yeast |
| `optimise` | bool | `True` | — |
| `gc_target` | int | `50` | — |

Galaxy Tool Shed: `EMBOSS backtranse / genSmart`

### codon_usage_table

**Codon usage and relative synonymous codon use** — Counts per codon with amino acid, RSCU and the per-usage table of a CDS set.

| parameter | kind | default | options |
|---|---|---|---|
| `cds` | file | `examples/genes.fa` | — |
| `frame` | choice | `1` | 1, 2, 3 |
| `rscu` | bool | `True` | — |
| `head` | int | `100` | — |

Galaxy Tool Shed: `CAI / codonW / EMBOSS cusp`

### dnds_pairwise

**Nei-Gojobori dN/dS between two CDS** — Synonymous and non-synonymous rates with the Jukes-Cantor corrected omega ratio.

| parameter | kind | default | options |
|---|---|---|---|
| `sequence_a` | file | `examples/genes.fa` | — |
| `sequence_b` | file | `examples/genes.fa` | — |
| `max_omega` | number | `50.0` | — |

Galaxy Tool Shed: `codeml / yn00 (PAML)`

### dnds_windows

**dN/dS in sliding codon windows** — Locate stretches of elevated or reduced selective pressure along the alignment.

| parameter | kind | default | options |
|---|---|---|---|
| `sequence_a` | file | `examples/genes.fa` | — |
| `sequence_b` | file | `examples/genes.fa` | — |
| `window` | int | `30` | — |
| `step` | int | `15` | — |
| `head` | int | `200` | — |

Galaxy Tool Shed: `slidingWindow / codeml windows`

### gc3_profile

**GC content in sliding windows** — Sliding-window GC (or GC3) profile, the classic isoform/isochores diagnostic.

| parameter | kind | default | options |
|---|---|---|---|
| `sequences` | file | `examples/genome.fa` | — |
| `window` | int | `500` | — |
| `step` | int | `0` | — |
| `gc3` | bool | `False` | — |

Galaxy Tool Shed: `EMBOSS geece / GCplot`

### tree_from_alignment_nj

**Neighbour-joining tree from an alignment** — Distance tree with optional bootstrap support values, written as Newick.

| parameter | kind | default | options |
|---|---|---|---|
| `alignment` | file | `examples/msa.fasta` | — |
| `model` | choice | `identity` | identity, p, kimura, jukes-cantor |
| `bootstrap` | int | `0` | — |
| `seed` | int | `1` | — |

Galaxy Tool Shed: `fasttree / phyml / iqtree`

## regional_variation

*Regional Variation* — in Galaxy group *Genomics Analysis*.

| Tool | What it does | Upstream Galaxy tool |
|---|---|---|
| [`differential_signal_regions`](#differential_signal_regions) | Binned fold-change between two coverage tracks with per-bin p-values. | `diffReps / DESeq2 on peaks` |
| [`vcf_cnv_segmentation`](#vcf_cnv_segmentation) | Circular binary segmentation of log-ratio depth into CNV regions. | `iuc/facets` |

### differential_signal_regions

**Differential signal between two tracks** — Binned fold-change between two coverage tracks with per-bin p-values.

| parameter | kind | default | options |
|---|---|---|---|
| `treatment` | file | `examples/coverage.bedgraph` | — |
| `control` | file | `examples/coverage.bedgraph` | — |
| `bin` | int | `500` | — |
| `min_log2fc` | number | `1.0` | — |
| `head` | int | `200` | — |

Galaxy Tool Shed: `diffReps / DESeq2 on peaks`

### vcf_cnv_segmentation

**Copy-number segmentation from VCF depths** — Circular binary segmentation of log-ratio depth into CNV regions.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/variants.vcf` | — |
| `expected_tumour_ploidy` | int | `2` | — |
| `min_segment_length` | number | `5.0` | — |
| `delta` | number | `0.01` | — |

Galaxy Tool Shed: `iuc/facets`

## chromosome_conformation

*Chromosome Conformation (Hi-C)* — in Galaxy group *Genomics Analysis*.

| Tool | What it does | Upstream Galaxy tool |
|---|---|---|
| [`hic_call_tads`](#hic_call_tads) | Merge adjacent bins whose internal contacts dominate the flanks. | `insulation-domain-finding / HiCExplorer hicFindTADs` |
| [`hic_compartments_eigenvec`](#hic_compartments_eigenvec) | First eigenvector of the correlation matrix, signed by a genomic track. | `HiCExplorer hicCorrectMatrix / eigenvecsfinder` |
| [`hic_insulation_score`](#hic_insulation_score) | Insulation dip score per bin with the strongest boundaries reported. | `HiCExplorer insulation / cooltools insulation` |
| [`hic_matrix_normalisation`](#hic_matrix_normalisation) | Matrix balancing so every row carries the same total coverage. | `cooler balance / HiCExplorer matrix balance` |
| [`hic_observed_expected`](#hic_observed_expected) | Contact frequency as a function of genomic separation, in log space. | `HiCExplorer plot_detector / cooltools econtact` |
| [`hic_pileup`](#hic_pileup) | Average contacts around each feature, corrected for the decay curve. | `HiCExplorer hicPileup / fit-hic-composite-log` |

### hic_call_tads

**Call domains from a contact matrix** — Merge adjacent bins whose internal contacts dominate the flanks.

| parameter | kind | default | options |
|---|---|---|---|
| `pairs` | file | `examples/contacts.pairs` | — |
| `bin_size` | int | `500` | — |
| `min_ratio` | number | `1.2` | — |
| `min_bins` | int | `2` | — |

Galaxy Tool Shed: `insulation-domain-finding / HiCExplorer hicFindTADs`

### hic_compartments_eigenvec

**A/B compartment eigenvector** — First eigenvector of the correlation matrix, signed by a genomic track.

| parameter | kind | default | options |
|---|---|---|---|
| `pairs` | file | `examples/contacts.pairs` | — |
| `bin_size` | int | `600` | — |
| `eigenvector` | int | `1` | — |
| `correlate_with` | choice | `none` | none, gc_content, coverage |
| `signal_track` | file | `examples/coverage.bedgraph` | — |

Galaxy Tool Shed: `HiCExplorer hicCorrectMatrix / eigenvecsfinder`

### hic_insulation_score

**Insulation profile and boundary calls** — Insulation dip score per bin with the strongest boundaries reported.

| parameter | kind | default | options |
|---|---|---|---|
| `pairs` | file | `examples/contacts.pairs` | — |
| `bin_size` | int | `500` | — |
| `window` | int | `2` | — |
| `strength` | number | `0.2` | — |

Galaxy Tool Shed: `HiCExplorer insulation / cooltools insulation`

### hic_matrix_normalisation

**Balance and normalise a contact matrix** — Matrix balancing so every row carries the same total coverage.

| parameter | kind | default | options |
|---|---|---|---|
| `pairs` | file | `examples/contacts.pairs` | — |
| `bin_size` | int | `500` | — |
| `method` | choice | `VC_sqrt` | ICE, VC, VC_sqrt, KR, none |
| `iterations` | int | `10` | — |
| `matrix` | bool | `True` | — |

Galaxy Tool Shed: `cooler balance / HiCExplorer matrix balance`

### hic_observed_expected

**Distance-decay (observed / expected) curve** — Contact frequency as a function of genomic separation, in log space.

| parameter | kind | default | options |
|---|---|---|---|
| `pairs` | file | `examples/contacts.pairs` | — |
| `bins` | int | `60` | — |
| `bin_size` | int | `400` | — |
| `log_y` | bool | `True` | — |

Galaxy Tool Shed: `HiCExplorer plot_detector / cooltools econtact`

### hic_pileup

**Pile-up around a set of features** — Average contacts around each feature, corrected for the decay curve.

| parameter | kind | default | options |
|---|---|---|---|
| `pairs` | file | `examples/contacts.pairs` | — |
| `features` | file | `examples/regions.bed` | — |
| `flank` | int | `2000` | — |
| `resolution` | int | `20` | — |
| `normalisation` | choice | `distance_decay` | none, distance_decay |

Galaxy Tool Shed: `HiCExplorer hicPileup / fit-hic-composite-log`

## transposon_insertion_sequencing

*Transposon Insertion Sequencing* — in Galaxy group *Genomics Analysis*.

| Tool | What it does | Upstream Galaxy tool |
|---|---|---|
| [`tnseq_ess`](#tnseq_ess) | Combine fold change and insertion count into an essentiality classification. | `TRANSIT-PROC / BaDIeS` |
| [`tnseq_hgc_flux`](#tnseq_hgc_flux) | Normalise insertion counts per gene to obtain a simple flux estimate. | `TRANSIT-PROC / ESME flux` |
| [`tnseq_log2fc`](#tnseq_log2fc) | Normalised log2 fold changes with per-gene significance for screens. | `ESME / TRANSIT` |

### tnseq_ess

**Essentiality score per gene** — Combine fold change and insertion count into an essentiality classification.

| parameter | kind | default | options |
|---|---|---|---|
| `counts_a` | file | `examples/counts.tsv` | — |
| `counts_b` | file | `examples/counts.tsv` | — |
| `cut_low` | number | `-1.0` | — |
| `cut_high` | number | `1.0` | — |

Galaxy Tool Shed: `TRANSIT-PROC / BaDIeS`

### tnseq_hgc_flux

**Gene-level insertion counts (HGC flux)** — Normalise insertion counts per gene to obtain a simple flux estimate.

| parameter | kind | default | options |
|---|---|---|---|
| `alignment` | file | `examples/alignments.sam` | — |
| `annotation` | file | `examples/annotation.gff` | — |
| `min_insertions` | int | `0` | — |

Galaxy Tool Shed: `TRANSIT-PROC / ESME flux`

### tnseq_log2fc

**Log2 fold change between insertion libraries** — Normalised log2 fold changes with per-gene significance for screens.

| parameter | kind | default | options |
|---|---|---|---|
| `counts_a` | file | `examples/counts.tsv` | — |
| `counts_b` | file | `examples/counts.tsv` | — |
| `pseudocount` | number | `1.0` | — |
| `alpha` | number | `0.05` | — |

Galaxy Tool Shed: `ESME / TRANSIT`

## protein_modeling

*Protein Structure* — in Galaxy group *Genomics Analysis*.

| Tool | What it does | Upstream Galaxy tool |
|---|---|---|
| [`protein_align_identities`](#protein_align_identities) | Identity of every query/reference pair using scoring or k-mer overlap. | `mmseqs easy-search / blastp identity` |
| [`protein_charge_properties`](#protein_charge_properties) | Net charge and pI per chain, plus the pH window where it buffers best. | `IPC / EMBOSS charge` |
| [`protein_composition_table`](#protein_composition_table) | Composition of every residue class with global physicochemical indices. | `ProtComp / pepstats composition` |
| [`protein_disorder`](#protein_disorder) | Per-residue disorder probability with disordered region calls. | `IUPred / Disopred` |
| [`protein_domain_sites`](#protein_domain_sites) | Conserved pattern matches along each chain with their positions. | `InterProScan / PaSE` |
| [`protein_hydrophobicity_moment`](#protein_hydrophobicity_moment) | Amphipathicity of each window, the classic test for facial helices. | `WEB-based helical wheel / EMBOSS heel` |
| [`protein_mass_and_formula`](#protein_mass_and_formula) | Molecular weight, empirical formula, extinction and molar absorptivity. | `peptide mass / expasy compute pI` |
| [`protein_secondary_structure`](#protein_secondary_structure) | Helix, strand and coil assignment from propensities and hydropathy. | `PSIPRED / Jpred` |
| [`protein_signal_peptide`](#protein_signal_peptide) | N-region hydrophobicity and cleavage-site score for secretion prediction. | `SignalP / Phobius` |
| [`protein_structure_summary_figure`](#protein_structure_summary_figure) | Windowed feature track along a chain, for a quick structural overview. | `UCSF ChimeraX / PyMOL feature view` |
| [`protein_transmembrane`](#protein_transmembrane) | Predicted TM segments with length, centred hydropathy and orientation. | `TMHMM / Phobius` |

### protein_align_identities

**Pairwise protein identity matrix** — Identity of every query/reference pair using scoring or k-mer overlap.

| parameter | kind | default | options |
|---|---|---|---|
| `protein_file` | file | `examples/proteins.faa` | — |
| `reference_file` | file | `examples/reference_proteins.faa` | — |
| `alignment` | choice | `pairwise` | pairwise, hamming |
| `min_identity` | number | `0.0` | — |
| `head` | int | `100` | — |

Galaxy Tool Shed: `mmseqs easy-search / blastp identity`

### protein_charge_properties

**Charge, pI and buffer behaviour** — Net charge and pI per chain, plus the pH window where it buffers best.

| parameter | kind | default | options |
|---|---|---|---|
| `protein_file` | file | `examples/proteins.faa` | — |
| `pH` | number | `7.0` | — |
| `ionic_strength` | number | `0.15` | — |
| `profile` | bool | `True` | — |

Galaxy Tool Shed: `IPC / EMBOSS charge`

### protein_composition_table

**Amino-acid composition and physicochemistry** — Composition of every residue class with global physicochemical indices.

| parameter | kind | default | options |
|---|---|---|---|
| `protein_file` | file | `examples/proteins.faa` | — |
| `normalisation` | choice | `percent` | percent, count, per_1000 |
| `properties` | bool | `True` | — |

Galaxy Tool Shed: `ProtComp / pepstats composition`

### protein_disorder

**Disorder propensity profile** — Per-residue disorder probability with disordered region calls.

| parameter | kind | default | options |
|---|---|---|---|
| `protein_file` | file | `examples/proteins.faa` | — |
| `cutoff` | number | `0.5` | — |
| `smoothing` | int | `5` | — |
| `protein_id` | text | `` | — |

Galaxy Tool Shed: `IUPred / Disopred`

### protein_domain_sites

**Prosite-style motif and domain sites** — Conserved pattern matches along each chain with their positions.

| parameter | kind | default | options |
|---|---|---|---|
| `protein_file` | file | `examples/proteins.faa` | — |
| `head` | int | `400` | — |

Galaxy Tool Shed: `InterProScan / PaSE`

### protein_hydrophobicity_moment

**Hydrophobic moment of helical segments** — Amphipathicity of each window, the classic test for facial helices.

| parameter | kind | default | options |
|---|---|---|---|
| `protein_file` | file | `examples/proteins.faa` | — |
| `helix_period` | int | `3.6` | — |
| `window` | number | `18` | — |
| `head` | int | `200` | — |

Galaxy Tool Shed: `WEB-based helical wheel / EMBOSS heel`

### protein_mass_and_formula

**Mass, formula and extinction coefficients** — Molecular weight, empirical formula, extinction and molar absorptivity.

| parameter | kind | default | options |
|---|---|---|---|
| `protein_file` | file | `examples/proteins.faa` | — |
| `mode` | choice | `average` | average, mono |
| `pH` | number | `8.0` | — |
| `reduced` | bool | `True` | — |

Galaxy Tool Shed: `peptide mass / expasy compute pI`

### protein_secondary_structure

**Predicted secondary structure** — Helix, strand and coil assignment from propensities and hydropathy.

| parameter | kind | default | options |
|---|---|---|---|
| `protein_file` | file | `examples/proteins.faa` | — |
| `window` | int | `9` | — |
| `helix_cutoff` | number | `0.6` | — |
| `profile` | bool | `True` | — |

Galaxy Tool Shed: `PSIPRED / Jpred`

### protein_signal_peptide

**Signal peptide and targeting prediction** — N-region hydrophobicity and cleavage-site score for secretion prediction.

| parameter | kind | default | options |
|---|---|---|---|
| `protein_file` | file | `examples/proteins.faa` | — |
| `cutoff` | number | `0.5` | — |
| `n_window` | int | `15` | — |

Galaxy Tool Shed: `SignalP / Phobius`

### protein_structure_summary_figure

**Structure-feature summary figure** — Windowed feature track along a chain, for a quick structural overview.

| parameter | kind | default | options |
|---|---|---|---|
| `protein_file` | file | `examples/proteins.faa` | — |
| `window` | int | `15` | — |
| `tm_cutoff` | number | `1.4` | — |
| `track` | choice | `hydropathy` | hydropathy, charge, disorder |
| `protein_id` | text | `` | — |

Galaxy Tool Shed: `UCSF ChimeraX / PyMOL feature view`

### protein_transmembrane

**Transmembrane helix prediction** — Predicted TM segments with length, centred hydropathy and orientation.

| parameter | kind | default | options |
|---|---|---|---|
| `protein_file` | file | `examples/proteins.faa` | — |
| `window` | int | `19` | — |
| `cutoff` | number | `1.6` | — |
| `topology` | bool | `True` | — |

Galaxy Tool Shed: `TMHMM / Phobius`

## genome_editing

*Genome Editing (CRISPR)* — in Galaxy group *Genomics Analysis*.

| Tool | What it does | Upstream Galaxy tool |
|---|---|---|
| [`base_editor_targets`](#base_editor_targets) | List editable bases falling in the activity window for C or A base editors. | `BE-Hive / CRISPResso design` |
| [`crispr_find_targets`](#crispr_find_targets) | Scan both strands for protospacer-PAM pairs passing GC and Tm filters. | `CHOPCHOP / CRISPOR` |
| [`crispr_grna_quality`](#crispr_grna_quality) | Rule-based activity score from GC, seed composition, poly-T terminator and hairpins. | `Doench Rule Set 2 / benchling` |
| [`crispr_knockout_design`](#crispr_knockout_design) | Rank early-exon guides by NMD proximity, GC content and predicted frameshift. | `CRISPR Design / synthego` |
| [`crispr_off_targets`](#crispr_off_targets) | Count near-match sites of every guide in a genome, weighted by mismatch position. | `CRISPOR / cas-OFFinder` |
| [`hdr_template_design`](#hdr_template_design) | Build left/right homology arms and the PCR primers to amplify a donor template. | `repair template design (DELIVERA / EuReCA)` |

### base_editor_targets

**Base-editing window inside protospacers** — List editable bases falling in the activity window for C or A base editors.

| parameter | kind | default | options |
|---|---|---|---|
| `target` | file | `examples/genes.fa` | — |
| `editor` | choice | `C-to-T` | C-to-T, A-to-G |
| `window_start` | int | `4` | — |
| `window_end` | int | `8` | — |
| `head` | int | `200` | — |

Galaxy Tool Shed: `BE-Hive / CRISPResso design`

### crispr_find_targets

**Find CRISPR cut sites for a target sequence** — Scan both strands for protospacer-PAM pairs passing GC and Tm filters.

| parameter | kind | default | options |
|---|---|---|---|
| `target` | file | `examples/genes.fa` | — |
| `pam` | text | `NGG` | — |
| `min_gc` | number | `30.0` | — |
| `max_gc` | number | `70.0` | — |
| `min_tm` | number | `52.0` | — |
| `head` | int | `200` | — |

Galaxy Tool Shed: `CHOPCHOP / CRISPOR`

### crispr_grna_quality

**Score guide efficiency and hairpin formation** — Rule-based activity score from GC, seed composition, poly-T terminator and hairpins.

| parameter | kind | default | options |
|---|---|---|---|
| `guides` | file | `examples/genes.fa` | — |
| `head` | int | `100` | — |
| `check_chromatin_like` | bool | `False` | — |

Galaxy Tool Shed: `Doench Rule Set 2 / benchling`

### crispr_knockout_design

**Choose guides for a frameshift knockout** — Rank early-exon guides by NMD proximity, GC content and predicted frameshift.

| parameter | kind | default | options |
|---|---|---|---|
| `annotation` | file | `examples/annotation.gff` | — |
| `genome` | file | `examples/genome.fa` | — |
| `exon_rank` | int | `2` | — |
| `head` | int | `100` | — |

Galaxy Tool Shed: `CRISPR Design / synthego`

### crispr_off_targets

**Genome-wide off-target search for guides** — Count near-match sites of every guide in a genome, weighted by mismatch position.

| parameter | kind | default | options |
|---|---|---|---|
| `guides` | file | `examples/genes.fa` | — |
| `genome` | file | `examples/genome.fa` | — |
| `max_mismatches` | int | `3` | — |
| `pam` | text | `NGG` | — |
| `head` | int | `300` | — |

Galaxy Tool Shed: `CRISPOR / cas-OFFinder`

### hdr_template_design

**Design a donor with homology arms** — Build left/right homology arms and the PCR primers to amplify a donor template.

| parameter | kind | default | options |
|---|---|---|---|
| `target` | file | `examples/genes.fa` | — |
| `cut_position` | int | `60` | — |
| `arm_length` | int | `500` | — |
| `insert` | code | `` | — |
| `max_primer_tm` | int | `68.0` | — |

Galaxy Tool Shed: `repair template design (DELIVERA / EuReCA)`

## biodiversity_data_exploration

*Biodiversity Data Exploration* — in Galaxy group *Genomics Analysis*.

| Tool | What it does | Upstream Galaxy tool |
|---|---|---|
| [`biodiv_barcode_gap`](#biodiv_barcode_gap) | Compare within- and between-group barcode distances to find a barcoding gap. | `ABGD / spades` |
| [`biodiv_checklist`](#biodiv_checklist) | Unique names at a rank with how often each was observed. | `GBIF checklist / spocc` |
| [`biodiv_endemism_score`](#biodiv_endemism_score) | Score each taxon by how narrowly it is represented in the sample set. | `spocc / endemicity index` |
| [`biodiv_marker_match`](#biodiv_marker_match) | Rank reference matches per query and flag confident identifications. | `BLAST / IDTaxa` |

### biodiv_barcode_gap

**Distance-based species delimitation** — Compare within- and between-group barcode distances to find a barcoding gap.

| parameter | kind | default | options |
|---|---|---|---|
| `records` | file | `examples/genes.fa` | — |
| `threshold` | number | `0.03` | — |
| `head` | int | `60` | — |

Galaxy Tool Shed: `ABGD / spades`

### biodiv_checklist

**Species checklist from a taxonomy table** — Unique names at a rank with how often each was observed.

| parameter | kind | default | options |
|---|---|---|---|
| `taxonomy` | file | `examples/taxmap.tsv` | — |
| `rank` | choice | `species` | kingdom, phylum, class, order, family, genus, species |
| `counts` | bool | `True` | — |

Galaxy Tool Shed: `GBIF checklist / spocc`

### biodiv_endemism_score

**Range-size proxy from taxon record counts** — Score each taxon by how narrowly it is represented in the sample set.

| parameter | kind | default | options |
|---|---|---|---|
| `taxonomy` | file | `examples/taxmap.tsv` | — |
| `rank` | choice | `genus` | kingdom, phylum, class, order, family, genus, species |
| `head` | int | `60` | — |

Galaxy Tool Shed: `spocc / endemicity index`

### biodiv_marker_match

**Match barcodes to a reference library** — Rank reference matches per query and flag confident identifications.

| parameter | kind | default | options |
|---|---|---|---|
| `records` | file | `examples/genes.fa` | — |
| `reference` | file | `examples/genes.fa` | — |
| `min_identity` | number | `97.0` | — |
| `max_results` | int | `5` | — |

Galaxy Tool Shed: `BLAST / IDTaxa`

## sequence_contamination_filtering

*Sequence Contamination Filtering* — in Galaxy group *Genomics Analysis*.

| Tool | What it does | Upstream Galaxy tool |
|---|---|---|
| [`blast_lca_taxonomy`](#blast_lca_taxonomy) | Classify each query by the lowest common ancestor of its best BLAST hits. | `MEGAN / kraken LCA` |
| [`contamination_filter`](#contamination_filter) | Flag sequences with a high-identity hit in the contaminant database and filter them out. | `deconseq / blastn screen` |
| [`dereplicate_sequences`](#dereplicate_sequences) | Group sequences by checksum and keep one representative with an abundance count. | `vsearch --derep_fulllength` |
| [`mash_screen`](#mash_screen) | MinHash sketch of every sequence and its Mash distance to each target. | `mash dist` |

### blast_lca_taxonomy

**Assign lowest common ancestor taxonomy** — Classify each query by the lowest common ancestor of its best BLAST hits.

| parameter | kind | default | options |
|---|---|---|---|
| `query` | file | `examples/genes.fa` | — |
| `database` | file | `examples/genome.fa` | — |
| `taxmap` | code | `` | — |
| `word_size` | int | `11` | — |
| `top_hits` | int | `5` | — |

Galaxy Tool Shed: `MEGAN / kraken LCA`

### contamination_filter

**Remove contaminant / host sequences** — Flag sequences with a high-identity hit in the contaminant database and filter them out.

| parameter | kind | default | options |
|---|---|---|---|
| `reads` | file | `examples/genes.fa` | — |
| `host` | file | `examples/genome.fa` | — |
| `min_identity` | number | `90.0` | — |
| `word_size` | int | `11` | — |
| `keep` | choice | `clean` | clean, hits |

Galaxy Tool Shed: `deconseq / blastn screen`

### dereplicate_sequences

**Collapse identical sequences** — Group sequences by checksum and keep one representative with an abundance count.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/genes.fa` | — |
| `count_by` | choice | `abundance` | abundance, first |
| `min_size` | int | `1` | — |

Galaxy Tool Shed: `vsearch --derep_fulllength`

### mash_screen

**Mash / MinHash distance between sequence sets** — MinHash sketch of every sequence and its Mash distance to each target.

| parameter | kind | default | options |
|---|---|---|---|
| `query` | file | `examples/genome.fa` | — |
| `database` | file | `examples/genome.fa` | — |
| `k` | int | `21` | — |
| `num` | int | `500` | — |

Galaxy Tool Shed: `mash dist`

## genome_diversity

*Genome Diversity* — in Galaxy group *Genomics Analysis*.

| Tool | What it does | Upstream Galaxy tool |
|---|---|---|
| [`vcf_fst_pairwise`](#vcf_fst_pairwise) | Genetic differentiation between groups. | `iuc/vcftools` |
| [`vcf_genetic_distance`](#vcf_genetic_distance) | Pairwise Fst-derived distance measures. | `iuc/stacks` |
| [`vcf_heterozygosity`](#vcf_heterozygosity) | Het/hom counts, F (inbreeding) and per-sample heterozygosity. | `iuc/plink` |
| [`vcf_hwe_summary`](#vcf_hwe_summary) | Aggregated HWE test statistics over all sites. | `iuc/plink` |
| [`vcf_nucleotide_diversity`](#vcf_nucleotide_diversity) | pi, number of segregating sites and per-chromosome totals. | `iuc/vcftools` |
| [`vcf_polymorphism_info`](#vcf_polymorphism_info) | PIC, MAF, heterozygosity and identity for each marker. | `iuc/cervus` |
| [`vcf_priv_alleles`](#vcf_priv_alleles) | Alleles found only in one population (private allele richness). | `iuc/private_alleles` |
| [`vcf_shannon_diversity`](#vcf_shannon_diversity) | Diversity of the allele frequency spectrum within windows. | `iuc/vegan` |
| [`vcf_tajima_d`](#vcf_tajima_d) | Sliding-window nucleotide diversity and Tajima's D. | `iuc/vcftools` |

### vcf_fst_pairwise

**Pairwise Fst (Hudson)** — Genetic differentiation between groups.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/variants.vcf` | — |
| `groups` | code | `` | — |
| `estimator` | choice | `hudson` | hudson, weir_cockerham |

Galaxy Tool Shed: `iuc/vcftools`

### vcf_genetic_distance

**Nei's genetic distance between populations** — Pairwise Fst-derived distance measures.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/variants.vcf` | — |
| `groups` | code | `` | — |
| `metric` | choice | `nei_d` | nei_d, da, reynolds, probsim |

Galaxy Tool Shed: `iuc/stacks`

### vcf_heterozygosity

**Observed/expected heterozygosity per sample** — Het/hom counts, F (inbreeding) and per-sample heterozygosity.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/variants.vcf` | — |

Galaxy Tool Shed: `iuc/plink`

### vcf_hwe_summary

**Population-level HWE statistics** — Aggregated HWE test statistics over all sites.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/variants.vcf` | — |

Galaxy Tool Shed: `iuc/plink`

### vcf_nucleotide_diversity

**Nucleotide diversity per population** — pi, number of segregating sites and per-chromosome totals.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/variants.vcf` | — |
| `groups` | code | `` | — |
| `window` | int | `0` | — |

Galaxy Tool Shed: `iuc/vcftools`

### vcf_polymorphism_info

**Polymorphism information content per locus** — PIC, MAF, heterozygosity and identity for each marker.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/variants.vcf` | — |
| `max_records` | int | `200` | — |

Galaxy Tool Shed: `iuc/cervus`

### vcf_priv_alleles

**Private and fixed-difference alleles** — Alleles found only in one population (private allele richness).

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/variants.vcf` | — |
| `groups` | code | `` | — |
| `min_frequency` | number | `0.05` | — |

Galaxy Tool Shed: `iuc/private_alleles`

### vcf_shannon_diversity

**Allelic diversity indices** — Diversity of the allele frequency spectrum within windows.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/variants.vcf` | — |
| `window` | int | `5000` | — |
| `index` | choice | `shannon` | shannon, simpson, pielou, berger |

Galaxy Tool Shed: `iuc/vegan`

### vcf_tajima_d

**Tajima's D in windows** — Sliding-window nucleotide diversity and Tajima's D.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/variants.vcf` | — |
| `window` | int | `2000` | — |
| `min_sites` | number | `1` | — |
| `chrom` | text | `` | — |

Galaxy Tool Shed: `iuc/vcftools`

## statistics

*Statistics* — in Galaxy group *Statistics and Visualization*.

| Tool | What it does | Upstream Galaxy tool |
|---|---|---|
| [`stats_anova`](#stats_anova) | Multi-group comparison with group means, F or H statistic and p-value. | `aov / kruskal.test / car::leveneTest` |
| [`stats_auc_roc`](#stats_auc_roc) | Sensitivity/specificity trade-off curve with the area under it. | `pROC / roc()` |
| [`stats_bootstrap_ci`](#stats_bootstrap_ci) | Resampling confidence interval and standard error for a column statistic. | `boot / simpleboot` |
| [`stats_contingency`](#stats_contingency) | Test association in a 2x2 slice or the whole table, with expected counts. | `chisq.test / fisher.test` |
| [`stats_correlation_matrix`](#stats_correlation_matrix) | Pairwise correlation coefficients for every column pair as a heatmap. | `corrplot / Hmisc rcorr` |
| [`stats_correlation_pair`](#stats_correlation_pair) | Correlation coefficient, test statistic and p-value for one pair of columns. | `cor.test` |
| [`stats_describe_columns`](#stats_describe_columns) | n, mean, median, spread, quantiles and skewness for each numeric column. | `datamash describe / summary` |
| [`stats_distribution_plot`](#stats_distribution_plot) | Empirical distribution of a value, or of every column for comparison. | `hist / ggplot2 density` |
| [`stats_group_summary`](#stats_group_summary) | Group-wise aggregation with sample size, spread and standard error. | `datamash groupby / aggregate` |
| [`stats_mutual_information`](#stats_mutual_information) | Information each column shares with a target, for nonlinear feature ranking. | `sklearn mutual_info_regression` |
| [`stats_normality`](#stats_normality) | Skew, kurtosis, Jarque-Bera p-value and equality-of-variance diagnostics per column. | `shapiro / bartlett.test` |
| [`stats_outlier_detection`](#stats_outlier_detection) | Flag rows whose value deviates from the robust centre of the distribution. | `outlier detection (base R)` |
| [`stats_p_adjust_table`](#stats_p_adjust_table) | Adjust a column of p-values and flag the features that stay significant. | `p.adjust / BH` |
| [`stats_power_analysis`](#stats_power_analysis) | Required group size for a given effect, or the power of a planned experiment. | `pwr.t.test / power.t.test` |
| [`stats_qqplot`](#stats_qqplot) | Check normality of a metric visually - heavy tails bend away from the line. | `qqnorm / qqplot` |
| [`stats_rank_tests`](#stats_rank_tests) | Distribution-free group comparison for small or non-normal samples. | `wilcox.test / exactRankTests` |
| [`stats_regression`](#stats_regression) | Slope, intercept, R^2, residual standard error and the slope p-value. | `lm / statsmodels OLS` |
| [`stats_summary_of_bedgraph`](#stats_summary_of_bedgraph) | Signal statistics of a genomic track: mean, percentiles and zero fraction. | `deepTools plotCoverage / summary` |
| [`stats_time_series_trend`](#stats_time_series_trend) | EWMA smoothing, Durbin-Watson residual statistic and cumulative-sum change points. | `trend / cusum (statistics)` |
| [`stats_ttest`](#stats_ttest) | Student, Welch or paired t-test per feature with the mean difference and effect size. | `t.test / statsmodels` |
| [`stats_two_sample_ks`](#stats_two_sample_ks) | Maximum ECDF difference and its p-value between two samples. | `ks.test` |

### stats_anova

**One-way ANOVA and Kruskal-Wallis by group** — Multi-group comparison with group means, F or H statistic and p-value.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/phenotypes.tsv` | — |
| `value_column` | text | `` | — |
| `group_column` | text | `` | — |
| `test` | choice | `anova` | anova, kruskal, levene, bartlett |

Galaxy Tool Shed: `aov / kruskal.test / car::leveneTest`

### stats_auc_roc

**ROC curve and AUC for a classifier score** — Sensitivity/specificity trade-off curve with the area under it.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/counts.tsv` | — |
| `score_column` | text | `` | — |
| `truth_column` | text | `` | — |
| `positive` | text | `` | — |

Galaxy Tool Shed: `pROC / roc()`

### stats_bootstrap_ci

**Bootstrap confidence interval of a statistic** — Resampling confidence interval and standard error for a column statistic.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/counts.tsv` | — |
| `column` | text | `` | — |
| `stat` | choice | `mean` | mean, median, difference, ratio |
| `second_column` | text | `` | — |
| `resamples` | int | `1000` | — |
| `level` | number | `0.95` | — |
| `seed` | int | `1` | — |

Galaxy Tool Shed: `boot / simpleboot`

### stats_contingency

**Chi-square and Fisher tests on a contingency table** — Test association in a 2x2 slice or the whole table, with expected counts.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/counts.tsv` | — |
| `test` | choice | `chi_square` | chi_square, fisher, g_test |
| `row_a` | text | `` | — |
| `row_b` | text | `` | — |
| `col_a` | text | `` | — |
| `col_b` | text | `` | — |

Galaxy Tool Shed: `chisq.test / fisher.test`

### stats_correlation_matrix

**Correlation matrix of columns** — Pairwise correlation coefficients for every column pair as a heatmap.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/counts.tsv` | — |
| `method` | choice | `pearson` | pearson, spearman, kendall |
| `annot` | bool | `True` | — |
| `top` | int | `0` | — |

Galaxy Tool Shed: `corrplot / Hmisc rcorr`

### stats_correlation_pair

**Correlate two columns with a p-value** — Correlation coefficient, test statistic and p-value for one pair of columns.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/counts.tsv` | — |
| `column_x` | text | `` | — |
| `column_y` | text | `` | — |
| `method` | choice | `pearson` | pearson, spearman, kendall |

Galaxy Tool Shed: `cor.test`

### stats_describe_columns

**Descriptive statistics of every numeric column** — n, mean, median, spread, quantiles and skewness for each numeric column.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/counts.tsv` | — |
| `include_index` | bool | `False` | — |

Galaxy Tool Shed: `datamash describe / summary`

### stats_distribution_plot

**Histogram with density overlay** — Empirical distribution of a value, or of every column for comparison.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/counts.tsv` | — |
| `column` | text | `` | — |
| `bins` | int | `25` | — |
| `log_y` | bool | `False` | — |
| `multi` | bool | `False` | — |

Galaxy Tool Shed: `hist / ggplot2 density`

### stats_group_summary

**Summarise a value column by group** — Group-wise aggregation with sample size, spread and standard error.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/phenotypes.tsv` | — |
| `value_column` | text | `` | — |
| `group_column` | text | `` | — |
| `agg` | choice | `mean` | mean, median, sum, count, stdev, min, max |
| `sem` | bool | `True` | — |

Galaxy Tool Shed: `datamash groupby / aggregate`

### stats_mutual_information

**Mutual information feature scores** — Information each column shares with a target, for nonlinear feature ranking.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/counts.tsv` | — |
| `target_column` | text | `` | — |
| `bins` | int | `10` | — |
| `head` | int | `100` | — |

Galaxy Tool Shed: `sklearn mutual_info_regression`

### stats_normality

**Normality and variance homogeneity tests** — Skew, kurtosis, Jarque-Bera p-value and equality-of-variance diagnostics per column.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/counts.tsv` | — |
| `head` | int | `30` | — |
| `variance_tests` | bool | `True` | — |

Galaxy Tool Shed: `shapiro / bartlett.test`

### stats_outlier_detection

**IQR, MAD and z-score outliers** — Flag rows whose value deviates from the robust centre of the distribution.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/counts.tsv` | — |
| `column` | text | `` | — |
| `iqr_factor` | number | `1.5` | — |
| `mad_factor` | number | `3.5` | — |
| `z_threshold` | number | `3.0` | — |
| `method` | choice | `iqr` | iqr, mad, zscore |

Galaxy Tool Shed: `outlier detection (base R)`

### stats_p_adjust_table

**Multiple-testing correction of a p-value column** — Adjust a column of p-values and flag the features that stay significant.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/counts.tsv` | — |
| `p_column` | text | `` | — |
| `method` | choice | `fdr_bh` | fdr_bh, bonferroni, holm, holm_bonferroni |
| `alpha` | number | `0.05` | — |
| `head` | int | `300` | — |

Galaxy Tool Shed: `p.adjust / BH`

### stats_power_analysis

**Power / sample size calculation** — Required group size for a given effect, or the power of a planned experiment.

| parameter | kind | default | options |
|---|---|---|---|
| `effect_size` | number | `0.5` | — |
| `alpha` | number | `0.05` | — |
| `power` | number | `0.8` | — |
| `n` | int | `0` | — |

Galaxy Tool Shed: `pwr.t.test / power.t.test`

### stats_qqplot

**Quantile-quantile plot against the normal distribution** — Check normality of a metric visually - heavy tails bend away from the line.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/counts.tsv` | — |
| `column` | text | `` | — |
| `dist` | choice | `normal` | normal, uniform, exponential |

Galaxy Tool Shed: `qqnorm / qqplot`

### stats_rank_tests

**Mann-Whitney and Wilcoxon tests** — Distribution-free group comparison for small or non-normal samples.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/counts.tsv` | — |
| `group_a` | text | `` | — |
| `group_b` | text | `` | — |
| `test` | choice | `mann_whitney` | mann_whitney, wilcoxon |

Galaxy Tool Shed: `wilcox.test / exactRankTests`

### stats_regression

**Linear regression of y on x** — Slope, intercept, R^2, residual standard error and the slope p-value.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/counts.tsv` | — |
| `x_column` | text | `` | — |
| `y_column` | text | `` | — |
| `order` | int | `1` | — |
| `loess` | bool | `False` | — |

Galaxy Tool Shed: `lm / statsmodels OLS`

### stats_summary_of_bedgraph

**Distribution of a coverage track** — Signal statistics of a genomic track: mean, percentiles and zero fraction.

| parameter | kind | default | options |
|---|---|---|---|
| `graph` | file | `examples/coverage.bedgraph` | — |
| `bins` | int | `20` | — |
| `weight` | choice | `length` | none, length |

Galaxy Tool Shed: `deepTools plotCoverage / summary`

### stats_time_series_trend

**Trend, smoothing and autocorrelation** — EWMA smoothing, Durbin-Watson residual statistic and cumulative-sum change points.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/counts.tsv` | — |
| `value_column` | text | `` | — |
| `alpha` | number | `0.3` | — |
| `lag` | int | `1` | — |
| `cusum` | bool | `True` | — |

Galaxy Tool Shed: `trend / cusum (statistics)`

### stats_ttest

**Two-group t-test (and Welch)** — Student, Welch or paired t-test per feature with the mean difference and effect size.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/counts.tsv` | — |
| `group_a` | text | `` | — |
| `group_b` | text | `` | — |
| `variant` | choice | `student` | student, welch, paired |

Galaxy Tool Shed: `t.test / statsmodels`

### stats_two_sample_ks

**Kolmogorov-Smirnov test between columns** — Maximum ECDF difference and its p-value between two samples.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/counts.tsv` | — |
| `column_a` | text | `` | — |
| `column_b` | text | `` | — |

Galaxy Tool Shed: `ks.test`

## machine_learning

*Machine Learning* — in Galaxy group *Statistics and Visualization*.

| Tool | What it does | Upstream Galaxy tool |
|---|---|---|
| [`ml_anomaly_detection`](#ml_anomaly_detection) | Flag observations that are unusual in any feature using robust z-scores. | `anomaly detection (statistics)` |
| [`ml_classifier`](#ml_classifier) | Train a model on a feature matrix and report accuracy, confusion counts and importances. | `caret train / sklearn` |
| [`ml_cluster_validation`](#ml_cluster_validation) | Sweep the number of clusters and report silhouette and inertia for each k. | `NbClust / clusterCrit` |
| [`ml_dbscan`](#ml_dbscan) | Density-based clustering that leaves noise unassigned, with an eps diagnostic curve. | `cluster::dbscan / sklearn DBSCAN` |
| [`ml_embedding_plot`](#ml_embedding_plot) | Non-linear 2-D projection of a feature matrix, optionally coloured by group. | `Rtsne / uwot / cmdscale` |
| [`ml_feature_selection_cv`](#ml_feature_selection_cv) | How much each feature matters, measured by the error increase when shuffled. | `caret varImp / permutation importance` |
| [`ml_gmm_bic`](#ml_gmm_bic) | Fit mixtures of increasing complexity and pick the BIC optimum. | `mclust / GaussianMixture` |
| [`ml_hierarchical`](#ml_hierarchical) | Linkage dendrogram of the observations with the cut height reported. | `hclust / scipy linkage` |
| [`ml_kmeans`](#ml_kmeans) | Lloyd's algorithm with inertia, silhouette score and cluster sizes. | `kmeans / scikit-learn KMeans` |
| [`ml_pca`](#ml_pca) | Scores, loadings and explained variance of a principal component analysis. | `prcomp / sklearn PCA` |
| [`ml_regression`](#ml_regression) | Fit a continuous target and compare in-sample with cross-validated error. | `glmnet / xgboost / caret train` |
| [`ml_split_and_encode`](#ml_split_and_encode) | Split a table into training and test parts and encode its categorical column. | `caret createDataPartition / fastDummies` |

### ml_anomaly_detection

**z-score anomaly detection per feature** — Flag observations that are unusual in any feature using robust z-scores.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/counts.tsv` | — |
| `threshold` | number | `3.0` | — |
| `robust` | choice | `median` | median, mean |
| `head` | int | `200` | — |

Galaxy Tool Shed: `anomaly detection (statistics)`

### ml_classifier

**Train a classifier (tree, random forest, logistic, kNN, Naive Bayes)** — Train a model on a feature matrix and report accuracy, confusion counts and importances.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/counts.tsv` | — |
| `label_column` | text | `` | — |
| `model` | choice | `random_forest` | decision_tree, random_forest, logistic, naive_bayes, knn, svm |
| `test_percent` | int | `25` | — |
| `depth` | int | `3` | — |
| `seed` | int | `1` | — |

Galaxy Tool Shed: `caret train / sklearn`

### ml_cluster_validation

**Silhouette and gap-style cluster validation** — Sweep the number of clusters and report silhouette and inertia for each k.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/counts.tsv` | — |
| `k_min` | int | `2` | — |
| `k_max` | int | `6` | — |
| `metric` | choice | `both` | silhouette, inertia, both |
| `seed` | int | `1` | — |

Galaxy Tool Shed: `NbClust / clusterCrit`

### ml_dbscan

**DBSCAN density clustering with eps curve** — Density-based clustering that leaves noise unassigned, with an eps diagnostic curve.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/counts.tsv` | — |
| `eps` | number | `1.0` | — |
| `min_samples` | int | `3` | — |
| `standardize` | bool | `True` | — |

Galaxy Tool Shed: `cluster::dbscan / sklearn DBSCAN`

### ml_embedding_plot

**t-SNE / UMAP / MDS embedding** — Non-linear 2-D projection of a feature matrix, optionally coloured by group.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/counts.tsv` | — |
| `method` | choice | `umap` | tsne, umap, mds |
| `components` | int | `2` | — |
| `neighbours` | int | `5` | — |
| `label_column` | text | `` | — |

Galaxy Tool Shed: `Rtsne / uwot / cmdscale`

### ml_feature_selection_cv

**Permutation importance for a fitted model** — How much each feature matters, measured by the error increase when shuffled.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/counts.tsv` | — |
| `target_column` | text | `` | — |
| `model` | choice | `linear` | linear, ridge, tree |
| `repeats` | int | `10` | — |
| `seed` | int | `2` | — |

Galaxy Tool Shed: `caret varImp / permutation importance`

### ml_gmm_bic

**Gaussian mixture model selection by BIC** — Fit mixtures of increasing complexity and pick the BIC optimum.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/counts.tsv` | — |
| `max_components` | int | `4` | — |
| `standardize` | bool | `True` | — |

Galaxy Tool Shed: `mclust / GaussianMixture`

### ml_hierarchical

**Hierarchical clustering with cophenetic correlation** — Linkage dendrogram of the observations with the cut height reported.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/counts.tsv` | — |
| `method` | choice | `average` | average, complete, single, ward |
| `clusters` | int | `3` | — |
| `rows` | bool | `True` | — |

Galaxy Tool Shed: `hclust / scipy linkage`

### ml_kmeans

**k-means clustering of samples or features** — Lloyd's algorithm with inertia, silhouette score and cluster sizes.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/counts.tsv` | — |
| `k` | int | `3` | — |
| `iters` | int | `100` | — |
| `transpose` | choice | `rows` | rows, columns |
| `seed` | int | `1` | — |
| `standardize` | bool | `True` | — |

Galaxy Tool Shed: `kmeans / scikit-learn KMeans`

### ml_pca

**Principal component analysis** — Scores, loadings and explained variance of a principal component analysis.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/counts.tsv` | — |
| `components` | int | `2` | — |
| `standardize` | bool | `True` | — |
| `head` | int | `200` | — |

Galaxy Tool Shed: `prcomp / sklearn PCA`

### ml_regression

**Regression with cross-validation** — Fit a continuous target and compare in-sample with cross-validated error.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/counts.tsv` | — |
| `target_column` | text | `` | — |
| `model` | choice | `ridge` | linear, ridge, lasso, gradient_boosting, regression_tree |
| `folds` | int | `4` | — |
| `alpha` | number | `1.0` | — |
| `trees` | int | `40` | — |

Galaxy Tool Shed: `glmnet / xgboost / caret train`

### ml_split_and_encode

**Train/test split and one-hot preview** — Split a table into training and test parts and encode its categorical column.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/phenotypes.tsv` | — |
| `label_column` | text | `` | — |
| `test_fraction` | number | `0.25` | — |
| `seed` | int | `7` | — |
| `head` | int | `200` | — |

Galaxy Tool Shed: `caret createDataPartition / fastDummies`

## graph_display_data

*Graph/Display Data* — in Galaxy group *Statistics and Visualization*.

| Tool | What it does | Upstream Galaxy tool |
|---|---|---|
| [`network_betweenness`](#network_betweenness) | Approximate betweenness by sampling shortest paths between node pairs. | `igraph betweenness / networkx` |
| [`network_communities`](#network_communities) | Partition a graph into communities by iterative neighbour voting. | `igraph cluster_label_propagation` |
| [`network_flow_table`](#network_flow_table) | Aggregate transitions between two categorical columns into weighted flows. | `networkD3 sankey` |
| [`network_from_edge_list`](#network_from_edge_list) | Degree distribution, components, density and self-loop counts of a graph. | `igraph / networkx` |
| [`network_upset_plot`](#network_upset_plot) | Enumerate every combination of sets with its intersection size (UpSet table). | `UpSetR / ggupset` |
| [`network_venn`](#network_venn) | Shared and unique elements of two key columns with a Jaccard index. | `venneuler / ggVennDiagram` |

### network_betweenness

**Betweenness centrality of nodes** — Approximate betweenness by sampling shortest paths between node pairs.

| parameter | kind | default | options |
|---|---|---|---|
| `edges` | file | `examples/counts.tsv` | — |
| `samples` | int | `200` | — |
| `head` | int | `100` | — |

Galaxy Tool Shed: `igraph betweenness / networkx`

### network_communities

**Label-propagation communities** — Partition a graph into communities by iterative neighbour voting.

| parameter | kind | default | options |
|---|---|---|---|
| `edges` | file | `examples/counts.tsv` | — |
| `iterations` | int | `20` | — |
| `head` | int | `300` | — |

Galaxy Tool Shed: `igraph cluster_label_propagation`

### network_flow_table

**Sankey-style flow table between stages** — Aggregate transitions between two categorical columns into weighted flows.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/counts.tsv` | — |
| `from_column` | text | `` | — |
| `to_column` | text | `` | — |
| `head` | int | `200` | — |

Galaxy Tool Shed: `networkD3 sankey`

### network_from_edge_list

**Graph metrics from an edge list** — Degree distribution, components, density and self-loop counts of a graph.

| parameter | kind | default | options |
|---|---|---|---|
| `edges` | file | `examples/counts.tsv` | — |
| `directed` | bool | `False` | — |
| `head` | int | `200` | — |

Galaxy Tool Shed: `igraph / networkx`

### network_upset_plot

**Intersection sizes of gene or feature sets** — Enumerate every combination of sets with its intersection size (UpSet table).

| parameter | kind | default | options |
|---|---|---|---|
| `sets` | file | `examples/genesets.gmt` | — |
| `min_size` | int | `1` | — |
| `show_plot` | bool | `True` | — |

Galaxy Tool Shed: `UpSetR / ggupset`

### network_venn

**Two-set Venn comparison** — Shared and unique elements of two key columns with a Jaccard index.

| parameter | kind | default | options |
|---|---|---|---|
| `left` | file | `examples/counts.tsv` | — |
| `right` | file | `examples/counts.tsv` | — |
| `key_left` | text | `` | — |
| `key_right` | text | `` | — |

Galaxy Tool Shed: `venneuler / ggVennDiagram`

## interactive_tools

*Interactive Tools (Visualization)* — in Galaxy group *Statistics and Visualization*.

| Tool | What it does | Upstream Galaxy tool |
|---|---|---|
| [`interactive_heatmap_viewer`](#interactive_heatmap_viewer) | Display a matrix with value labels and row ordering by mean expression. | `ComplexHeatmap / plotly heatmap` |
| [`interactive_pivot_table`](#interactive_pivot_table) | Melt a wide table and pivot it back with an aggregation function. | `tidyr pivot_wider / reshape2` |
| [`interactive_scatter_explorer`](#interactive_scatter_explorer) | Scatter of two columns coloured by a third, plus the underlying point table. | `plotly / interactive scatter` |
| [`interactive_track_view`](#interactive_track_view) | Metaprofile of a signal track resampled across every annotated feature. | `IGV / deepTools plotProfile` |

### interactive_heatmap_viewer

**Annotated heatmap of a matrix** — Display a matrix with value labels and row ordering by mean expression.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/counts.tsv` | — |
| `zscore` | bool | `True` | — |
| `cmap` | choice | `coolwarm` | viridis, magma, coolwarm, Greys |
| `max_rows` | int | `40` | — |
| `max_cols` | int | `40` | — |

Galaxy Tool Shed: `ComplexHeatmap / plotly heatmap`

### interactive_pivot_table

**Pivot a table into a matrix** — Melt a wide table and pivot it back with an aggregation function.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/counts.tsv` | — |
| `index_column` | text | `` | — |
| `value_column` | text | `` | — |
| `agg` | choice | `mean` | mean, sum, count, max |

Galaxy Tool Shed: `tidyr pivot_wider / reshape2`

### interactive_scatter_explorer

**Scatter view with per-point table** — Scatter of two columns coloured by a third, plus the underlying point table.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/counts.tsv` | — |
| `x_column` | text | `` | — |
| `y_column` | text | `` | — |
| `colour_column` | text | `` | — |
| `label_column` | text | `` | — |
| `head` | int | `200` | — |

Galaxy Tool Shed: `plotly / interactive scatter`

### interactive_track_view

**Coverage track over annotated features** — Metaprofile of a signal track resampled across every annotated feature.

| parameter | kind | default | options |
|---|---|---|---|
| `graph` | file | `examples/coverage.bedgraph` | — |
| `annotation` | file | `examples/annotation.gff` | — |
| `bins` | int | `40` | — |
| `mode` | choice | `mean` | mean, max |

Galaxy Tool Shed: `IGV / deepTools plotProfile`

## str_fm__microsatellite_analysis

*STR-FM: Microsatellite Analysis* — in Galaxy group *Genomics Toolkits*.

| Tool | What it does | Upstream Galaxy tool |
|---|---|---|
| [`str_allele_table`](#str_allele_table) | Tidy an ampFLPSTR-style profile into one row per allele call. | `GeneMapper / STRAF` |
| [`str_expansion_test`](#str_expansion_test) | Flag samples whose allele sizes are outliers for the locus repeat structure. | `hipSTRyda / ExpansionHunter` |
| [`str_locus_statistics`](#str_locus_statistics) | Allele counts, observed and expected heterozygosity per microsatellite. | `STRAF / poppr` |
| [`str_profile_match`](#str_profile_match) | Count matching, off-by-one and discordant loci between two profiles. | `match.py / CODIS comparison` |
| [`str_relatedness`](#str_relatedness) | Estimate IBS sharing between samples to flag duplicates and relatives. | `King relatedness / LRmix` |

### str_allele_table

**Long-format STR allele table** — Tidy an ampFLPSTR-style profile into one row per allele call.

| parameter | kind | default | options |
|---|---|---|---|
| `profile` | file | `examples/str_profile.tsv` | — |
| `expand` | bool | `True` | — |

Galaxy Tool Shed: `GeneMapper / STRAF`

### str_expansion_test

**Test for repeat expansions at a locus** — Flag samples whose allele sizes are outliers for the locus repeat structure.

| parameter | kind | default | options |
|---|---|---|---|
| `profile` | file | `examples/str_profile.tsv` | — |
| `locus` | text | `` | — |
| `normal_mean` | number | `0.0` | — |
| `n_sds` | number | `3.0` | — |

Galaxy Tool Shed: `hipSTRyda / ExpansionHunter`

### str_locus_statistics

**Allelic diversity per STR locus** — Allele counts, observed and expected heterozygosity per microsatellite.

| parameter | kind | default | options |
|---|---|---|---|
| `profile` | file | `examples/str_profile.tsv` | — |
| `expected_heterozygosity` | bool | `True` | — |

Galaxy Tool Shed: `STRAF / poppr`

### str_profile_match

**Compare two STR profiles** — Count matching, off-by-one and discordant loci between two profiles.

| parameter | kind | default | options |
|---|---|---|---|
| `profile` | file | `examples/str_profile.tsv` | — |
| `sample_a` | text | `` | — |
| `sample_b` | text | `` | — |
| `tolerance` | int | `1` | — |

Galaxy Tool Shed: `match.py / CODIS comparison`

### str_relatedness

**Pairwise relatedness from shared alleles** — Estimate IBS sharing between samples to flag duplicates and relatives.

| parameter | kind | default | options |
|---|---|---|---|
| `profile` | file | `examples/str_profile.tsv` | — |
| `min_loci` | number | `1` | — |

Galaxy Tool Shed: `King relatedness / LRmix`

## mothur

*Mothur (16S rRNA)* — in Galaxy group *Genomics Toolkits*.

| Tool | What it does | Upstream Galaxy tool |
|---|---|---|
| [`mothur_classify_seqs`](#mothur_classify_seqs) | Best reference match per query with identity and the assigned lineage. | `mothur classify.seqs` |
| [`mothur_count_seqs`](#mothur_count_seqs) | Length distribution, GC and ambiguity statistics of a nucleotide library. | `mothur count.seqs` |
| [`mothur_dist_seqs`](#mothur_dist_seqs) | Every pair of sequences scored by an evolutionary or identity distance. | `mothur dist.seqs / phylip dist` |
| [`mothur_get_current`](#mothur_get_current) | Dereplicate a FASTA and report how often each unique sequence occurs. | `mothur get.current / count.seqs` |
| [`mothur_pre_cluster`](#mothur_pre_cluster) | Greedy agglomerative clustering of sequences into OTU-like groups. | `mothur pre.cluster / vsearch cluster` |
| [`mothur_screen_seqs`](#mothur_screen_seqs) | Drop sequences that are too short, too long, ambiguous or contain bad motifs. | `mothur screen.seqs` |
| [`mothur_trim_seqs`](#mothur_trim_seqs) | Cut every sequence to the same alignment window for comparable regions. | `mothur trim.seqs` |

### mothur_classify_seqs

**Classify sequences against a reference taxonomy** — Best reference match per query with identity and the assigned lineage.

| parameter | kind | default | options |
|---|---|---|---|
| `records` | file | `examples/genes.fa` | — |
| `reference` | file | `examples/genes.fa` | — |
| `taxonomy` | file | `examples/taxmap.tsv` | — |
| `cutoff` | number | `60.0` | — |
| `head` | int | `100` | — |

Galaxy Tool Shed: `mothur classify.seqs`

### mothur_count_seqs

**Sequence length and composition summary** — Length distribution, GC and ambiguity statistics of a nucleotide library.

| parameter | kind | default | options |
|---|---|---|---|
| `records` | file | `examples/genes.fa` | — |
| `bin` | int | `100` | — |

Galaxy Tool Shed: `mothur count.seqs`

### mothur_dist_seqs

**Pairwise sequence distances** — Every pair of sequences scored by an evolutionary or identity distance.

| parameter | kind | default | options |
|---|---|---|---|
| `records` | file | `examples/genes.fa` | — |
| `model` | choice | `identity` | identity, p, k2p, raw |
| `head` | int | `60` | — |

Galaxy Tool Shed: `mothur dist.seqs / phylip dist`

### mothur_get_current

**Unique sequences with abundances** — Dereplicate a FASTA and report how often each unique sequence occurs.

| parameter | kind | default | options |
|---|---|---|---|
| `records` | file | `examples/genes.fa` | — |
| `rename` | bool | `True` | — |
| `head` | int | `200` | — |

Galaxy Tool Shed: `mothur get.current / count.seqs`

### mothur_pre_cluster

**Cluster sequences at a distance cutoff** — Greedy agglomerative clustering of sequences into OTU-like groups.

| parameter | kind | default | options |
|---|---|---|---|
| `records` | file | `examples/genes.fa` | — |
| `distance` | number | `0.03` | — |
| `method` | choice | `average` | average, maximum, nearest |
| `head` | int | `200` | — |

Galaxy Tool Shed: `mothur pre.cluster / vsearch cluster`

### mothur_screen_seqs

**Screen sequences by length and ambiguity** — Drop sequences that are too short, too long, ambiguous or contain bad motifs.

| parameter | kind | default | options |
|---|---|---|---|
| `records` | file | `examples/genes.fa` | — |
| `min_length` | int | `10` | — |
| `max_length` | int | `100000` | — |
| `max_ambiguous` | int | `0` | — |
| `patterns` | text | `` | — |

Galaxy Tool Shed: `mothur screen.seqs`

### mothur_trim_seqs

**Trim sequences to a start and end position** — Cut every sequence to the same alignment window for comparable regions.

| parameter | kind | default | options |
|---|---|---|---|
| `records` | file | `examples/genes.fa` | — |
| `start` | int | `1` | — |
| `end` | int | `300` | — |
| `left_justify` | bool | `False` | — |

Galaxy Tool Shed: `mothur trim.seqs`

## qiime2

*QIIME 2* — in Galaxy group *Genomics Toolkits*.

| Tool | What it does | Upstream Galaxy tool |
|---|---|---|
| [`qiime_dereplicate`](#qiime_dereplicate) | Collapse identical sequences into features with a size column. | `qiime quality-control dereplicate` |
| [`qiime_fastq_stats`](#qiime_fastq_stats) | Read counts, length profile and mean quality of an amplicon library. | `qiime demux summarize / fastp` |
| [`qiime_filter_features`](#qiime_filter_features) | Drop rare or low-abundance features and report what survives. | `qiime feature-table filter-features` |
| [`qiime_merge_metadata`](#qiime_merge_metadata) | Attach one tabular annotation to another by a shared key column. | `qiime metadata join / phyloseq merge` |
| [`qiime_taxa_barplot`](#qiime_taxa_barplot) | Relative composition of every sample, grouped at a taxonomic rank. | `qiime taxa barplot` |
| [`qiime_trim_primers`](#qiime_trim_primers) | Locate and remove amplicon primers, reporting how many reads were affected. | `qiime cutadapt trim-pairs` |

### qiime_dereplicate

**Dereplicate and count feature sequences** — Collapse identical sequences into features with a size column.

| parameter | kind | default | options |
|---|---|---|---|
| `records` | file | `examples/genes.fa` | — |
| `relabelling` | bool | `True` | — |

Galaxy Tool Shed: `qiime quality-control dereplicate`

### qiime_fastq_stats

**Per-length and per-quality FASTQ summary** — Read counts, length profile and mean quality of an amplicon library.

| parameter | kind | default | options |
|---|---|---|---|
| `reads` | file | `examples/reads_single.fastq` | — |
| `quality_window` | int | `4` | — |

Galaxy Tool Shed: `qiime demux summarize / fastp`

### qiime_filter_features

**Filter a feature table by prevalence and depth** — Drop rare or low-abundance features and report what survives.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/counts.tsv` | — |
| `min_prevalence` | number | `0.5` | — |
| `min_total` | int | `10` | — |
| `min_per_sample` | int | `0` | — |

Galaxy Tool Shed: `qiime feature-table filter-features`

### qiime_merge_metadata

**Join a feature table with sample metadata** — Attach one tabular annotation to another by a shared key column.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/phenotypes.tsv` | — |
| `extra` | file | `examples/counts.tsv` | — |
| `left_key` | text | `` | — |
| `right_key` | text | `` | — |
| `how` | choice | `left` | left, inner |

Galaxy Tool Shed: `qiime metadata join / phyloseq merge`

### qiime_taxa_barplot

**Stacked bar plot of taxon composition** — Relative composition of every sample, grouped at a taxonomic rank.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/counts.tsv` | — |
| `taxonomy` | file | `examples/taxmap.tsv` | — |
| `rank` | choice | `phylum` | kingdom, phylum, class, order, family, genus, species |
| `top` | int | `5` | — |

Galaxy Tool Shed: `qiime taxa barplot`

### qiime_trim_primers

**Cut primer sequences from reads** — Locate and remove amplicon primers, reporting how many reads were affected.

| parameter | kind | default | options |
|---|---|---|---|
| `reads` | file | `examples/reads_1.fastq` | — |
| `forward_primer` | text | `` | — |
| `reverse_primer` | text | `` | — |
| `remove_primers` | bool | `True` | — |
| `failure` | choice | `discard` | keep, discard |

Galaxy Tool Shed: `qiime cutadapt trim-pairs`

## picard

*Picard* — in Galaxy group *Genomics Toolkits*.

| Tool | What it does | Upstream Galaxy tool |
|---|---|---|
| [`samtools_alignment_summary`](#samtools_alignment_summary) | Total/proper pairs, mismatch and indel rates per read. | `iuc/picard` |
| [`samtools_fingerprint`](#samtools_fingerprint) | CheckBias/fingerprint style binned read-count profile. | `iuc/picard` |
| [`samtools_gc_bias`](#samtools_gc_bias) | CollectGcBiasMetrics-style observed vs expected GC distribution. | `iuc/picard` |
| [`samtools_library_complexity`](#samtools_library_complexity) | Duplicate-based library size estimate (L = N ln(N/(N-M))). | `iuc/picard` |
| [`samtools_validate`](#samtools_validate) | Structural validation of SAM records (CIGAR vs length, flags). | `iuc/picard` |

### samtools_alignment_summary

**Picard AlignmentSummaryMetrics** — Total/proper pairs, mismatch and indel rates per read.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/alignments.sam` | — |
| `read_length` | int | `100` | — |

Galaxy Tool Shed: `iuc/picard`

### samtools_fingerprint

**Reads fingerprint (cross-sample check)** — CheckBias/fingerprint style binned read-count profile.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/alignments.sam` | — |
| `genome` | file | `examples/genome.fa` | — |
| `binsize` | int | `1000` | — |

Galaxy Tool Shed: `iuc/picard`

### samtools_gc_bias

**GC bias of the library** — CollectGcBiasMetrics-style observed vs expected GC distribution.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/alignments.sam` | — |
| `genome` | file | `examples/genome.fa` | — |
| `bins` | int | `10` | — |

Galaxy Tool Shed: `iuc/picard`

### samtools_library_complexity

**Picard EstimateLibraryComplexity** — Duplicate-based library size estimate (L = N ln(N/(N-M))).

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/alignments.sam` | — |
| `pairs` | int | `100000` | — |

Galaxy Tool Shed: `iuc/picard`

### samtools_validate

**Picard ValidateSamFile style checks** — Structural validation of SAM records (CIGAR vs length, flags).

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/alignments.sam` | — |
| `strict` | bool | `False` | — |
| `max_errors` | int | `20` | — |

Galaxy Tool Shed: `iuc/picard`

## deeptools

*deepTools* — in Galaxy group *Genomics Toolkits*.

| Tool | What it does | Upstream Galaxy tool |
|---|---|---|
| [`deeptools_bam_fragmentsize`](#deeptools_bam_fragmentsize) | Fragment-size statistics used by deepTools before extending reads. | `iuc/deepTools` |
| [`deeptools_bamcompare`](#deeptools_bamcompare) | Per-bin comparison of two alignment files. | `iuc/deepTools` |
| [`deeptools_bamcoverage`](#deeptools_bamcoverage) | deepTools bamCoverage: normalised bedGraph of read depth. | `iuc/deepTools` |
| [`deeptools_computeMatrix`](#deeptools_computeMatrix) | Matrix of coverage values per region, for heatmaps. | `iuc/deepTools` |
| [`deeptools_multibamsummary`](#deeptools_multibamsummary) | Count reads per genomic bin for several BAM files at once. | `iuc/deepTools` |
| [`deeptools_normalize_per_read`](#deeptools_normalize_per_read) | Rescale coverage values to counts per million mapped reads. | `iuc/deepTools` |
| [`deeptools_plot_heatmap`](#deeptools_plot_heatmap) | Render a coverage heatmap of sorted regions. | `iuc/deepTools` |
| [`deeptools_plotprofile`](#deeptools_plotprofile) | Average read density around a set of regions. | `iuc/deepTools` |

### deeptools_bam_fragmentsize

**bamPEFragmentSize summary** — Fragment-size statistics used by deepTools before extending reads.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/alignments.sam` | — |
| `max_fraction` | int | `2` | — |

Galaxy Tool Shed: `iuc/deepTools`

### deeptools_bamcompare

**bamCompare (log2 ratio of two samples)** — Per-bin comparison of two alignment files.

| parameter | kind | default | options |
|---|---|---|---|
| `a` | file | `examples/alignments.sam` | — |
| `b` | file | `examples/alignments.sam` | — |
| `operation` | choice | `log2` | log2, ratio, difference, sub |
| `bin_size` | int | `500` | — |
| `pseudocount` | number | `1.0` | — |

Galaxy Tool Shed: `iuc/deepTools`

### deeptools_bamcoverage

**bamCoverage (normalise depth)** — deepTools bamCoverage: normalised bedGraph of read depth.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/alignments.sam` | — |
| `normalisation` | choice | `none` | none, RPKM, FPKM, BAM, CPS |
| `scale_factor` | number | `1.0` | — |
| `bin_size` | int | `0` | — |
| `effective_genome_length` | int | `1` | — |
| `extend_reads` | number | `0.0` | — |
| `genome` | file | `` | — |

Galaxy Tool Shed: `iuc/deepTools`

### deeptools_computeMatrix

**computeMatrix (regions x bins matrix)** — Matrix of coverage values per region, for heatmaps.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/alignments.sam` | — |
| `regions` | file | `examples/regions.bed` | — |
| `bins` | int | `10` | — |
| `up` | int | `300` | — |
| `down` | int | `300` | — |
| `scale` | bool | `True` | — |

Galaxy Tool Shed: `iuc/deepTools`

### deeptools_multibamsummary

**multiBamSummary** — Count reads per genomic bin for several BAM files at once.

| parameter | kind | default | options |
|---|---|---|---|
| `files` | multi | `['examples/alignments.sam']` | examples/alignments.sam |
| `bin_size` | int | `1000` | — |
| `all_bins` | bool | `False` | — |

Galaxy Tool Shed: `iuc/deepTools`

### deeptools_normalize_per_read

**Normalise a bedGraph per million reads** — Rescale coverage values to counts per million mapped reads.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/coverage.bedgraph` | — |
| `total_reads` | int | `1000000` | — |
| `scale` | number | `1.0` | — |

Galaxy Tool Shed: `iuc/deepTools`

### deeptools_plot_heatmap

**plotHeatmap figure** — Render a coverage heatmap of sorted regions.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/alignments.sam` | — |
| `regions` | file | `examples/regions.bed` | — |
| `bins` | int | `24` | — |
| `up` | int | `400` | — |
| `down` | int | `400` | — |

Galaxy Tool Shed: `iuc/deepTools`

### deeptools_plotprofile

**plotProfile data (per-bin mean)** — Average read density around a set of regions.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/alignments.sam` | — |
| `regions` | file | `examples/regions.bed` | — |
| `flank` | int | `500` | — |
| `bins` | int | `20` | — |
| `per_interval` | bool | `False` | — |

Galaxy Tool Shed: `iuc/deepTools`

## emboss

*EMBOSS* — in Galaxy group *Genomics Toolkits*.

| Tool | What it does | Upstream Galaxy tool |
|---|---|---|
| [`emboss_einverted`](#emboss_einverted) | Direct and inverted repeats within a sequence, with their positions. | `EMBOSS einverted` |
| [`emboss_geece`](#emboss_geece) | Fraction of each sequence covered by ORFs on the forward and reverse strands. | `EMBOSS geece` |
| [`emboss_getorf`](#emboss_getorf) | All ORFs above a size threshold, written out as a nucleotide FASTA. | `EMBOSS getorf` |
| [`emboss_newcleft`](#emboss_newcleft) | Report enzymes whose cuts leave at least one very large fragment. | `EMBOSS newcleft` |
| [`emboss_pepstats`](#emboss_pepstats) | Molecular weight, pI, charge, instability and hydrophobicity per chain. | `EMBOSS pepstats / ProtParam` |
| [`emboss_transeq`](#emboss_transeq) | Six-frame translation with a length filter, as a FASTA result. | `EMBOSS transeq` |
| [`emboss_water`](#emboss_water) | Optimal alignment with identity, gaps and the raw score. | `EMBOSS water / needle` |

### emboss_einverted

**Find inverted (palindromic) repeats** — Direct and inverted repeats within a sequence, with their positions.

| parameter | kind | default | options |
|---|---|---|---|
| `records` | file | `examples/genes.fa` | — |
| `threshold` | int | `8` | — |
| `period` | int | `100000` | — |
| `head` | int | `200` | — |

Galaxy Tool Shed: `EMBOSS einverted`

### emboss_geece

**Coding potential of each DNA strand** — Fraction of each sequence covered by ORFs on the forward and reverse strands.

| parameter | kind | default | options |
|---|---|---|---|
| `records` | file | `examples/genes.fa` | — |
| `min_len` | int | `100` | — |
| `genetic_code` | choice | `Standard` | Standard, Vertebrate Mitochondrial |

Galaxy Tool Shed: `EMBOSS geece`

### emboss_getorf

**Extract open reading frames** — All ORFs above a size threshold, written out as a nucleotide FASTA.

| parameter | kind | default | options |
|---|---|---|---|
| `records` | file | `examples/genes.fa` | — |
| `min_len` | int | `90` | — |
| `frame` | choice | `all` | all, 1, 2, 3 |
| `required_start` | text | `ATG` | — |
| `partial` | bool | `True` | — |

Galaxy Tool Shed: `EMBOSS getorf`

### emboss_newcleft

**Find restriction enzymes that leave no cut** — Report enzymes whose cuts leave at least one very large fragment.

| parameter | kind | default | options |
|---|---|---|---|
| `records` | file | `examples/genes.fa` | — |
| `min_size` | int | `200` | — |
| `max_enzymes` | int | `25` | — |

Galaxy Tool Shed: `EMBOSS newcleft`

### emboss_pepstats

**Protein physicochemical statistics** — Molecular weight, pI, charge, instability and hydrophobicity per chain.

| parameter | kind | default | options |
|---|---|---|---|
| `proteins` | file | `examples/proteins.faa` | — |
| `pH` | number | `7.0` | — |
| `mode` | choice | `average` | average, mono |

Galaxy Tool Shed: `EMBOSS pepstats / ProtParam`

### emboss_transeq

**Translate nucleotides in six frames** — Six-frame translation with a length filter, as a FASTA result.

| parameter | kind | default | options |
|---|---|---|---|
| `records` | file | `examples/genes.fa` | — |
| `frame` | choice | `all` | 1, 2, 3, -1, -2, -3, all |
| `to_stop` | bool | `False` | — |
| `min_len` | int | `30` | — |

Galaxy Tool Shed: `EMBOSS transeq`

### emboss_water

**Local pairwise alignment of two sequences** — Optimal alignment with identity, gaps and the raw score.

| parameter | kind | default | options |
|---|---|---|---|
| `sequence_a` | text | `ATGAAAGCTTTGCGATCGATCGATCGGCTAAGCATCGATCGATCGATTAA` | — |
| `sequence_b` | text | `ATGAAAGCTTTGCGATCGATCGATCGGCTAAGCATCGATCGATCGATTAA` | — |
| `mode` | choice | `local` | local, global |
| `gap_open` | number | `-10.0` | — |
| `gap_extend` | number | `-0.5` | — |
| `matrix` | choice | `` | , BLOSUM62, PAM250 |

Galaxy Tool Shed: `EMBOSS water / needle`

## ncbi_blast_

*NCBI BLAST+* — in Galaxy group *Genomics Toolkits*.

| Tool | What it does | Upstream Galaxy tool |
|---|---|---|
| [`blast_best_hit_per_query`](#blast_best_hit_per_query) | Reduce all-vs-all hits to the single top-scoring alignment per query. | `blast + best-hit filter` |
| [`blast_identity_histogram`](#blast_identity_histogram) | Histogram of percent identity across the reported hits. | `BLAST histogram` |
| [`blast_make_database`](#blast_make_database) | Write a database descriptor with lengths, alphabets and k-mer index size. | `makeblastdb` |
| [`blast_tabular_normalise`](#blast_tabular_normalise) | Rename columns of a -outfmt 6 table and filter by identity and E-value. | `blast -outfmt 6` |
| [`blastn_search`](#blastn_search) | Word-hit BLASTN search with extension, bit scores and E-values. | `blastn` |
| [`blastp_search`](#blastp_search) | Protein-protein search scored with BLOSUM/PAM matrices. | `blastp` |
| [`blastx_search`](#blastx_search) | Translate the query in six frames and search a protein database. | `blastx` |
| [`tblastn_search`](#tblastn_search) | Search a protein query against the six translated frames of a nucleotide database. | `tblastn` |

### blast_best_hit_per_query

**Keep the best hit of every query** — Reduce all-vs-all hits to the single top-scoring alignment per query.

| parameter | kind | default | options |
|---|---|---|---|
| `query` | file | `examples/genes.fa` | — |
| `database` | file | `examples/genome.fa` | — |
| `by` | choice | `bitscore` | bitscore, evalue |
| `word_size` | int | `11` | — |
| `evalue` | number | `0.001` | — |

Galaxy Tool Shed: `blast + best-hit filter`

### blast_identity_histogram

**Identity distribution of BLAST hits** — Histogram of percent identity across the reported hits.

| parameter | kind | default | options |
|---|---|---|---|
| `query` | file | `examples/genes.fa` | — |
| `database` | file | `examples/genome.fa` | — |
| `bin_size` | int | `5` | — |
| `word_size` | int | `11` | — |

Galaxy Tool Shed: `BLAST histogram`

### blast_make_database

**Prepare a BLAST database** — Write a database descriptor with lengths, alphabets and k-mer index size.

| parameter | kind | default | options |
|---|---|---|---|
| `sequences` | file | `examples/genome.fa` | — |
| `word_size` | int | `11` | — |
| `protein` | bool | `False` | — |

Galaxy Tool Shed: `makeblastdb`

### blast_tabular_normalise

**Normalise a BLAST tabular output** — Rename columns of a -outfmt 6 table and filter by identity and E-value.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/counts.tsv` | — |
| `columns` | int | `12` | — |
| `min_pident` | number | `0.0` | — |
| `max_evalue` | number | `1e+100` | — |

Galaxy Tool Shed: `blast -outfmt 6`

### blastn_search

**BLASTN: nucleotide query versus nucleotide database** — Word-hit BLASTN search with extension, bit scores and E-values.

| parameter | kind | default | options |
|---|---|---|---|
| `query` | file | `examples/genes.fa` | — |
| `database` | file | `examples/genome.fa` | — |
| `word_size` | int | `11` | — |
| `evalue` | number | `1e-05` | — |
| `max_target_seqs` | int | `10` | — |
| `dust` | bool | `True` | — |
| `gap_open` | int | `-11` | — |

Galaxy Tool Shed: `blastn`

### blastp_search

**BLASTP: protein query versus protein database** — Protein-protein search scored with BLOSUM/PAM matrices.

| parameter | kind | default | options |
|---|---|---|---|
| `query` | file | `examples/proteins.faa` | — |
| `database` | file | `examples/reference_proteins.faa` | — |
| `matrix` | choice | `BLOSUM62` | BLOSUM62, BLOSUM45, PAM250 |
| `evalue` | number | `10.0` | — |
| `word_size` | int | `3` | — |
| `max_target_seqs` | int | `10` | — |

Galaxy Tool Shed: `blastp`

### blastx_search

**BLASTX: nucleotide query in six frames** — Translate the query in six frames and search a protein database.

| parameter | kind | default | options |
|---|---|---|---|
| `query` | file | `examples/genes.fa` | — |
| `database` | file | `examples/reference_proteins.faa` | — |
| `word_size` | int | `3` | — |
| `evalue` | number | `10.0` | — |
| `matrix` | choice | `BLOSUM62` | BLOSUM62, PAM250 |

Galaxy Tool Shed: `blastx`

### tblastn_search

**TBLASTN: protein query versus translated database** — Search a protein query against the six translated frames of a nucleotide database.

| parameter | kind | default | options |
|---|---|---|---|
| `query` | file | `examples/proteins.faa` | — |
| `database` | file | `examples/genome.fa` | — |
| `word_size` | int | `3` | — |
| `evalue` | number | `10.0` | — |

Galaxy Tool Shed: `tblastn`

## rseqc

*RSeQC* — in Galaxy group *Genomics Toolkits*.

| Tool | What it does | Upstream Galaxy tool |
|---|---|---|
| [`rna_seq_mapping_qc`](#rna_seq_mapping_qc) | Mapped, duplicate, low-MAPQ and clipping statistics for an RNA-seq library. | `RSeQC / hisat2 reporting` |
| [`rseqc_gene_body_coverage`](#rseqc_gene_body_coverage) | Read density along the normalised gene body: 3' bias shows up immediately. | `RSeQC gene_body_coverage` |
| [`rseqc_inner_distance`](#rseqc_inner_distance) | Fragment inner-distance histogram - the Fritsch-Hill criterion for strandedness. | `RSeQC inner_distance` |
| [`rseqc_junction_annotation`](#rseqc_junction_annotation) | Infer N-gaps in CIGARs (or introns between blocks) and check them against the annotation. | `RSeQC junction_annotation` |
| [`rseqc_read_distribution`](#rseqc_read_distribution) | Split alignments into CDS, exon-intron, upstream and downstream categories. | `RSeQC read_distribution` |
| [`samtools_feature_counts`](#samtools_feature_counts) | HTSeq-style counting of alignments per feature. | `iuc/featureCounts` |
| [`samtools_inner_distance`](#samtools_inner_distance) | Mean insert size inside each feature (library QC). | `iuc/inner_distance` |
| [`samtools_junctions`](#samtools_junctions) | Spliced junctions (N in CIGAR) with support counts. | `iuc/hisat2` |
| [`samtools_metagene`](#samtools_metagene) | Average read density across scaled gene bodies. | `iuc/RSeQC` |
| [`samtools_read_distribution`](#samtools_read_distribution) | RSeQC read_distribution: reads in CDS/5'UTR/3'UTR/introns/intergenic. | `iuc/RSeQC` |

### rna_seq_mapping_qc

**RNA-seq alignment QC metrics** — Mapped, duplicate, low-MAPQ and clipping statistics for an RNA-seq library.

| parameter | kind | default | options |
|---|---|---|---|
| `alignments` | file | `examples/alignments.sam` | — |
| `head` | int | `100` | — |

Galaxy Tool Shed: `RSeQC / hisat2 reporting`

### rseqc_gene_body_coverage

**Gene body coverage profile** — Read density along the normalised gene body: 3' bias shows up immediately.

| parameter | kind | default | options |
|---|---|---|---|
| `alignments` | file | `examples/alignments.sam` | — |
| `annotation` | file | `examples/annotation.gff` | — |
| `bins` | int | `20` | — |

Galaxy Tool Shed: `RSeQC gene_body_coverage`

### rseqc_inner_distance

**Inner distance / insert size distribution** — Fragment inner-distance histogram - the Fritsch-Hill criterion for strandedness.

| parameter | kind | default | options |
|---|---|---|---|
| `alignments` | file | `examples/alignments.sam` | — |
| `bins` | int | `25` | — |
| `max_size` | int | `500` | — |

Galaxy Tool Shed: `RSeQC inner_distance`

### rseqc_junction_annotation

**Splice junctions from alignments and annotation** — Infer N-gaps in CIGARs (or introns between blocks) and check them against the annotation.

| parameter | kind | default | options |
|---|---|---|---|
| `alignments` | file | `examples/alignments.sam` | — |
| `annotation` | file | `examples/annotation.gff` | — |
| `min_intron` | int | `20` | — |
| `head` | int | `200` | — |

Galaxy Tool Shed: `RSeQC junction_annotation`

### rseqc_read_distribution

**Where reads fall relative to features** — Split alignments into CDS, exon-intron, upstream and downstream categories.

| parameter | kind | default | options |
|---|---|---|---|
| `alignments` | file | `examples/alignments.sam` | — |
| `annotation` | file | `examples/annotation.gff` | — |
| `transposons` | file | `examples/regions.bed` | — |
| `head` | int | `100` | — |

Galaxy Tool Shed: `RSeQC read_distribution`

### samtools_feature_counts

**Assign reads to features** — HTSeq-style counting of alignments per feature.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/alignments.sam` | — |
| `annotation` | file | `examples/annotation.gff` | — |
| `mode` | choice | `intersection_nonempty` | intersection_nonempty, union, strict |
| `count_read_overlaps` | bool | `False` | — |
| `nonunique` | bool | `True` | — |

Galaxy Tool Shed: `iuc/featureCounts`

### samtools_inner_distance

**Inner distance (fragment size per bin)** — Mean insert size inside each feature (library QC).

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/alignments.sam` | — |
| `regions` | file | `examples/regions.bed` | — |
| `bin_size` | int | `500` | — |

Galaxy Tool Shed: `iuc/inner_distance`

### samtools_junctions

**Exon junctions from alignments** — Spliced junctions (N in CIGAR) with support counts.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/alignments.sam` | — |
| `min_count` | int | `1` | — |

Galaxy Tool Shed: `iuc/hisat2`

### samtools_metagene

**Metagene profile of alignments** — Average read density across scaled gene bodies.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/alignments.sam` | — |
| `annotation` | file | `examples/annotation.gff` | — |
| `profile_bins` | int | `40` | — |
| `up` | int | `500` | — |
| `down` | int | `500` | — |

Galaxy Tool Shed: `iuc/RSeQC`

### samtools_read_distribution

**Read distribution relative to features** — RSeQC read_distribution: reads in CDS/5'UTR/3'UTR/introns/intergenic.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/alignments.sam` | — |
| `annotation` | file | `examples/annotation.gff` | — |
| `bin_size` | int | `100` | — |

Galaxy Tool Shed: `iuc/RSeQC`

## hca-scanpy

*Scanpy (Single Cell)* — in Galaxy group *Genomics Toolkits*.

| Tool | What it does | Upstream Galaxy tool |
|---|---|---|
| [`scanpy_cell_cycle_score`](#scanpy_cell_cycle_score) | Score S and G2M signatures per cell and call a phase for each one. | `scanpy tl.score_genes / cell cycle scoring` |
| [`scanpy_cluster_knn_graph`](#scanpy_cluster_knn_graph) | Cluster cells from a shared-nearest-neighbour graph without external packages. | `scanpy tl.leiden / scran buildSNNGraph` |
| [`scanpy_highly_variable`](#scanpy_highly_variable) | Binned dispersion ranking used to pick informative genes before PCA. | `scanpy pp.highly_variable_genes` |
| [`scanpy_qc_metrics`](#scanpy_qc_metrics) | Genes, UMIs, mito and ribosomal fractions per cell with filter flags. | `scanpy qc.calculate_qc_metrics` |

### scanpy_cell_cycle_score

**Cell-cycle phase scores** — Score S and G2M signatures per cell and call a phase for each one.

| parameter | kind | default | options |
|---|---|---|---|
| `s_phase_genes` | file | `geneA
GeneB
GENEC` | — |
| `g2m_genes` | file | `geneD
geneE` | — |
| `counts` | file | `examples/counts.tsv` | — |
| `normalisation` | choice | `zscore` | zscore, none |
| `head` | int | `100` | — |

Galaxy Tool Shed: `scanpy tl.score_genes / cell cycle scoring`

### scanpy_cluster_knn_graph

**k-NN graph clustering (SNN + label propagation)** — Cluster cells from a shared-nearest-neighbour graph without external packages.

| parameter | kind | default | options |
|---|---|---|---|
| `counts` | file | `examples/counts.tsv` | — |
| `neighbours` | int | `5` | — |
| `reduction` | choice | `pca` | pca, raw |
| `resolution` | int | `2` | — |
| `standardize` | bool | `True` | — |
| `head` | int | `200` | — |

Galaxy Tool Shed: `scanpy tl.leiden / scran buildSNNGraph`

### scanpy_highly_variable

**Select highly variable genes** — Binned dispersion ranking used to pick informative genes before PCA.

| parameter | kind | default | options |
|---|---|---|---|
| `counts` | file | `examples/counts.tsv` | — |
| `flavour` | choice | `seurat` | seurat, cell_ranger, seurat_v3 |
| `n_bins` | int | `6` | — |
| `top_genes` | int | `10` | — |
| `show_plot` | bool | `True` | — |

Galaxy Tool Shed: `scanpy pp.highly_variable_genes`

### scanpy_qc_metrics

**Scanpy-style per-cell QC metrics** — Genes, UMIs, mito and ribosomal fractions per cell with filter flags.

| parameter | kind | default | options |
|---|---|---|---|
| `counts` | file | `examples/counts.tsv` | — |
| `mitochondrial_prefix` | text | `MT-` | — |
| `ribosomal_prefix` | text | `RP` | — |
| `mt_cut` | number | `20.0` | — |
| `min_genes` | int | `1` | — |
| `head` | int | `100` | — |

Galaxy Tool Shed: `scanpy qc.calculate_qc_metrics`

## hicexplorer

*HiCExplorer* — in Galaxy group *Genomics Toolkits*.

| Tool | What it does | Upstream Galaxy tool |
|---|---|---|
| [`hic_compartments`](#hic_compartments) | First eigenvector of the correlation matrix - the classic A/B compartment signal. | `eigenvector_c / HiCExplorer` |
| [`hic_contact_matrix`](#hic_contact_matrix) | Bin ligation events into a matrix and optionally apply vanilla-ICE normalisation. | `cooler dump / HiCExplorer` |
| [`hic_insulation_profile`](#hic_insulation_profile) | Boundary detection from the fraction of contacts crossing each bin. | `insulation4phipson / HiCExplorer insulation` |

### hic_compartments

**A/B compartment eigenvector** — First eigenvector of the correlation matrix - the classic A/B compartment signal.

| parameter | kind | default | options |
|---|---|---|---|
| `pairs` | file | `examples/contacts.pairs` | — |
| `bins` | int | `1000` | — |
| `chrom` | choice | `all` | all |
| `head` | int | `100` | — |

Galaxy Tool Shed: `eigenvector_c / HiCExplorer`

### hic_contact_matrix

**Contact matrix from a pairs file** — Bin ligation events into a matrix and optionally apply vanilla-ICE normalisation.

| parameter | kind | default | options |
|---|---|---|---|
| `pairs` | file | `examples/contacts.pairs` | — |
| `bins` | int | `1000` | — |
| `chrom` | text | `` | — |
| `normalise` | bool | `True` | — |
| `head` | int | `50` | — |

Galaxy Tool Shed: `cooler dump / HiCExplorer`

### hic_insulation_profile

**Insulation score along the genome** — Boundary detection from the fraction of contacts crossing each bin.

| parameter | kind | default | options |
|---|---|---|---|
| `pairs` | file | `examples/contacts.pairs` | — |
| `bins` | int | `1000` | — |
| `window` | int | `3` | — |

Galaxy Tool Shed: `insulation4phipson / HiCExplorer insulation`

## seqtk

*SeqTK* — in Galaxy group *Genomics Toolkits*.

| Tool | What it does | Upstream Galaxy tool |
|---|---|---|
| [`seqtk_comp`](#seqtk_comp) | Per-record length, GC and per-base counts (``seqtk comp``). | `iuc/seqtk` |
| [`seqtk_dna2ca`](#seqtk_dna2ca) | Codon-level statistics per sequence. | `iuc/seqtk` |
| [`seqtk_fixhead`](#seqtk_fixhead) | Replace duplicate or invalid ids with sequential names. | `iuc/seqtk` |
| [`seqtk_gapidx`](#seqtk_gapidx) | List the coordinates of every gap in the sequences. | `iuc/seqtk` |
| [`seqtk_listheads`](#seqtk_listheads) | List all sequence headers. | `iuc/seqtk` |
| [`seqtk_locut`](#seqtk_locut) | Flag low complexity sequences/windows like dustmasker. | `iuc/dust` |
| [`seqtk_maskedcopy`](#seqtk_maskedcopy) | Replace dust-masked segments with N (``maskedcopy -m``). | `iuc/dustmasker` |
| [`seqtk_mut`](#seqtk_mut) | Introduce substitutions/indels described in a mutation list. | `iuc/seqtk` |
| [`seqtk_seq_subsample`](#seqtk_seq_subsample) | Randomly subsample reads/sequences with a fixed seed. | `iuc/seqtk` |
| [`seqtk_trim_qual`](#seqtk_trim_qual) | Trim low quality tails from FASTQ reads. | `iuc/seqtk` |

### seqtk_comp

**seqtk comp** — Per-record length, GC and per-base counts (``seqtk comp``).

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/genome.fa` | — |

Galaxy Tool Shed: `iuc/seqtk`

### seqtk_dna2ca

**seqtk dna2ca (codon conversion)** — Codon-level statistics per sequence.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/genome.fa` | — |
| `mode` | choice | `count_codons` | count_codons, gc3, cai |

Galaxy Tool Shed: `iuc/seqtk`

### seqtk_fixhead

**seqtk fixhead (repair deflines)** — Replace duplicate or invalid ids with sequential names.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/genome.fa` | — |
| `prefix` | text | `seq` | — |
| `start` | int | `1` | — |

Galaxy Tool Shed: `iuc/seqtk`

### seqtk_gapidx

**Report runs of N (gap index)** — List the coordinates of every gap in the sequences.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/genome.fa` | — |
| `min_len` | int | `1` | — |

Galaxy Tool Shed: `iuc/seqtk`

### seqtk_listheads

**seqtk listh (headers)** — List all sequence headers.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/genome.fa` | — |

Galaxy Tool Shed: `iuc/seqtk`

### seqtk_locut

**Low-complexity filter (locut/DUST)** — Flag low complexity sequences/windows like dustmasker.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/genome.fa` | — |
| `dust_cutoff` | number | `2.0` | — |
| `window` | int | `64` | — |

Galaxy Tool Shed: `iuc/dust`

### seqtk_maskedcopy

**Mask low complexity regions** — Replace dust-masked segments with N (``maskedcopy -m``).

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/genome.fa` | — |
| `level` | number | `2.4` | — |
| `window` | int | `64` | — |
| `mask_char` | text | `N` | — |

Galaxy Tool Shed: `iuc/dustmasker`

### seqtk_mut

**seqtk mutfa (apply mutations)** — Introduce substitutions/indels described in a mutation list.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/genome.fa` | — |
| `mutations` | code | `` | — |

Galaxy Tool Shed: `iuc/seqtk`

### seqtk_seq_subsample

**seqtk sample** — Randomly subsample reads/sequences with a fixed seed.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/reads_single.fastq` | — |
| `fraction` | number | `0.25` | — |
| `n` | int | `0` | — |
| `seed` | int | `100` | — |

Galaxy Tool Shed: `iuc/seqtk`

### seqtk_trim_qual

**seqtk trimq** — Trim low quality tails from FASTQ reads.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/reads_single.fastq` | — |
| `min_quality` | int | `20` | — |
| `trim_two_sided` | bool | `True` | — |

Galaxy Tool Shed: `iuc/seqtk`

## bbtools

*BBTools* — in Galaxy group *Genomics Toolkits*.

| Tool | What it does | Upstream Galaxy tool |
|---|---|---|
| [`bbduk_trim`](#bbduk_trim) | Trim adapter matches and poor-quality tails, dropping short leftovers. | `BBMap bbduk.sh` |
| [`bbrep_repearness`](#bbrep_repearness) | Fraction of a read made of repeated k-mers, for filtering low-complexity data. | `BBMap bbrep` |
| [`fastq_bbduk_contaminant`](#fastq_bbduk_contaminant) | ``bbduk.sh ktrim=r k=…``: drop or trim reads containing contaminant k-mers. | `iuc/bbduk` |
| [`fastq_bbinsert_check`](#fastq_bbinsert_check) | Flag abnormal fragment lengths in a paired alignment file. | `iuc/bbtools` |
| [`fastq_bblearn_errors`](#fastq_bblearn_errors) | Empirical ambiguity rate per quality bin (``errorsfromfastq``). | `iuc/bbduk` |
| [`fastq_bbmerge_report`](#fastq_bbmerge_report) | ``bbmerge.sh``-style overlap merging statistics. | `iuc/bbmerge` |
| [`fastq_bbnorm_qc`](#fastq_bbnorm_qc) | Group reads by mean quality bin (BBNorm-style diagnostics). | `iuc/bbtools` |
| [`fastq_bbsim_shard`](#fastq_bbsim_shard) | ``shards.sh``: split a FASTQ into balanced random shards. | `iuc/bbtools` |

### bbduk_trim

**Quality and adapter trimming (kmer filter)** — Trim adapter matches and poor-quality tails, dropping short leftovers.

| parameter | kind | default | options |
|---|---|---|---|
| `reads` | file | `examples/reads_single.fastq` | — |
| `qual_cutoff` | number | `20.0` | — |
| `min_len` | int | `25` | — |
| `adapters` | text | `AGATCGGAAGAGC` | — |
| `trim_qlow` | bool | `True` | — |

Galaxy Tool Shed: `BBMap bbduk.sh`

### bbrep_repearness

**Repearness of each read** — Fraction of a read made of repeated k-mers, for filtering low-complexity data.

| parameter | kind | default | options |
|---|---|---|---|
| `reads` | file | `examples/reads_single.fastq` | — |
| `kmer` | int | `6` | — |
| `threshold` | number | `0.3` | — |
| `head` | int | `200` | — |

Galaxy Tool Shed: `BBMap bbrep`

### fastq_bbduk_contaminant

**Scrub contaminant k-mers** — ``bbduk.sh ktrim=r k=…``: drop or trim reads containing contaminant k-mers.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/reads_single.fastq` | — |
| `contaminants` | code | `GGGGGGGGGGGGGGGG
AAAAAAAAAAAAAAAAAA` | — |
| `k` | int | `16` | — |
| `mismatches` | int | `0` | — |
| `trim` | bool | `True` | — |

Galaxy Tool Shed: `iuc/bbduk`

### fastq_bbinsert_check

**Insert-size and pair sanity check** — Flag abnormal fragment lengths in a paired alignment file.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/alignments.sam` | — |
| `sd_cutoff` | number | `3.0` | — |

Galaxy Tool Shed: `iuc/bbtools`

### fastq_bblearn_errors

**Error model from qualities** — Empirical ambiguity rate per quality bin (``errorsfromfastq``).

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/reads_single.fastq` | — |
| `k` | int | `1` | — |

Galaxy Tool Shed: `iuc/bbduk`

### fastq_bbmerge_report

**Merge mate pairs and report rates** — ``bbmerge.sh``-style overlap merging statistics.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/reads_1.fastq` | — |
| `src2` | file | `examples/reads_2.fastq` | — |
| `min_overlap` | int | `11` | — |
| `max_insert` | int | `1000` | — |

Galaxy Tool Shed: `iuc/bbmerge`

### fastq_bbnorm_qc

**Quality-based read binning** — Group reads by mean quality bin (BBNorm-style diagnostics).

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/reads_single.fastq` | — |
| `bin_size` | int | `5` | — |

Galaxy Tool Shed: `iuc/bbtools`

### fastq_bbsim_shard

**Shard reads into N files** — ``shards.sh``: split a FASTQ into balanced random shards.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/reads_single.fastq` | — |
| `shards` | int | `4` | — |
| `seed` | int | `1` | — |

Galaxy Tool Shed: `iuc/bbtools`

## motif

*Motif Tools* — in Galaxy group *Genomics Toolkits*.

| Tool | What it does | Upstream Galaxy tool |
|---|---|---|
| [`homer_central_enrichment`](#homer_central_enrichment) | Report the fraction of motif hits lying in the centre of each sequence (CEMS-lite). | `HOMER2 central enrichment` |
| [`homer_find_motifs`](#homer_find_motifs) | Compare motif occurrence rates between two sequence sets with a Fisher test. | `HOMER2 findMotifsGenome` |
| [`meme_discover_motifs`](#meme_discover_motifs) | Expectation-maximisation motif finder returning consensus, matrix and likelihood. | `MEME (MEME suite)` |
| [`motif_compare_tomtom`](#motif_compare_tomtom) | Match motifs by profile correlation over aligned offsets and report similarity. | `MEME suite Tomtom` |
| [`motif_count_degenerate`](#motif_count_degenerate) | Exact or mismatch-tolerant counting of a degenerate motif per sequence. | `seqkit locate / EMBOSS regexp` |
| [`motif_fimo_sites`](#motif_fimo_sites) | Convert PWM scores into approximate per-site p-values from the score null model. | `MEME suite FIMO` |
| [`motif_format_convert`](#motif_format_convert) | Re-emit matrices as normalised PFM, log-odds PWM, information matrix or consensus. | `motif format conversion (MEME/JASPAR)` |
| [`motif_from_aligned_sites`](#motif_from_aligned_sites) | Turn aligned motif instances into a count matrix and write it in PFM formats. | `MEME suite streme / alignment to PFM` |
| [`motif_information_plot`](#motif_information_plot) | Bar plot of per-position information content or base frequency of a motif. | `logomaker / seqlogo` |
| [`motif_jaspar_search`](#motif_jaspar_search) | Filter a motif database by name and report width, consensus and information. | `JASPAR web service` |
| [`motif_scan_sequence`](#motif_scan_sequence) | Score every window against the log-odds matrix and report thresholded sites. | `gtr-scanner / MOODS / FIMO` |
| [`motif_summary`](#motif_summary) | Per-motif width, consensus, total information and degenerate-regex summary. | `jaspar / motif db summary` |
| [`motif_transform`](#motif_transform) | Rewrite a matrix in the opposite orientation, trimmed, rounded or row-normalised. | `motif tools (MEME suite)` |

### homer_central_enrichment

**Central enrichment of a motif in sequences** — Report the fraction of motif hits lying in the centre of each sequence (CEMS-lite).

| parameter | kind | default | options |
|---|---|---|---|
| `sequences` | file | `examples/genes.fa` | — |
| `pattern` | text | `ACGT` | — |
| `centre_fraction` | number | `0.2` | — |
| `head` | int | `50` | — |

Galaxy Tool Shed: `HOMER2 central enrichment`

### homer_find_motifs

**Known-motif enrichment in target versus background** — Compare motif occurrence rates between two sequence sets with a Fisher test.

| parameter | kind | default | options |
|---|---|---|---|
| `target` | file | `examples/genes.fa` | — |
| `background` | file | `examples/genome.fa` | — |
| `kmer` | int | `6` | — |
| `motif_file` | file | `` | — |
| `head` | int | `100` | — |

Galaxy Tool Shed: `HOMER2 findMotifsGenome`

### meme_discover_motifs

**Discover motifs with MEME (EM)** — Expectation-maximisation motif finder returning consensus, matrix and likelihood.

| parameter | kind | default | options |
|---|---|---|---|
| `sequences` | file | `examples/genes.fa` | — |
| `width` | int | `8` | — |
| `nmotifs` | int | `1` | — |
| `sites` | choice | `zoops` | zoops, oops, anr |
| `iterations` | int | `20` | — |
| `seed` | int | `1` | — |

Galaxy Tool Shed: `MEME (MEME suite)`

### motif_compare_tomtom

**Compare two motif sets (Tomtom-lite)** — Match motifs by profile correlation over aligned offsets and report similarity.

| parameter | kind | default | options |
|---|---|---|---|
| `target` | file | `examples/motif.pfm` | — |
| `query` | file | `examples/motifs_jaspar.txt` | — |
| `mode` | choice | `both` | both, reverse, forward |
| `min_overlap` | int | `4` | — |

Galaxy Tool Shed: `MEME suite Tomtom`

### motif_count_degenerate

**Count IUPAC motif occurrences** — Exact or mismatch-tolerant counting of a degenerate motif per sequence.

| parameter | kind | default | options |
|---|---|---|---|
| `sequences` | file | `examples/genes.fa` | — |
| `pattern` | text | `ACGT` | — |
| `allow_mismatch` | int | `0` | — |
| `strand` | choice | `both` | both, plus, minus |
| `head` | int | `200` | — |

Galaxy Tool Shed: `seqkit locate / EMBOSS regexp`

### motif_fimo_sites

**FIMO-like sites with p-values** — Convert PWM scores into approximate per-site p-values from the score null model.

| parameter | kind | default | options |
|---|---|---|---|
| `target` | file | `examples/genes.fa` | — |
| `matrix` | file | `examples/motif.pfm` | — |
| `p_value` | number | `0.0001` | — |
| `thresh_all` | bool | `True` | — |
| `head` | int | `300` | — |

Galaxy Tool Shed: `MEME suite FIMO`

### motif_format_convert

**Convert PFM to PWM / ICM / consensus formats** — Re-emit matrices as normalised PFM, log-odds PWM, information matrix or consensus.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/motifs_jaspar.txt` | — |
| `format` | choice | `jaspar` | jaspar, pwm, icm, meme, consensus |
| `pseudocount` | number | `1.0` | — |

Galaxy Tool Shed: `motif format conversion (MEME/JASPAR)`

### motif_from_aligned_sites

**Build a PFM from aligned sites** — Turn aligned motif instances into a count matrix and write it in PFM formats.

| parameter | kind | default | options |
|---|---|---|---|
| `sites` | code | `ACGTACGT
AGGTACGT
ACGTACGA` | — |
| `name` | text | `CHROMA.1` | — |
| `alphabet` | choice | `ACGT` | ACGT, ARN |
| `style` | choice | `jaspar` | jaspar, counts, frequency |

Galaxy Tool Shed: `MEME suite streme / alignment to PFM`

### motif_information_plot

**Sequence-logo-style information content** — Bar plot of per-position information content or base frequency of a motif.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/motif.pfm` | — |
| `metric` | choice | `bits` | bits, frequency |

Galaxy Tool Shed: `logomaker / seqlogo`

### motif_jaspar_search

**Query a JASPAR-like motif file by name** — Filter a motif database by name and report width, consensus and information.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/motifs_jaspar.txt` | — |
| `pattern` | text | `` | — |
| `min_information` | number | `0.0` | — |

Galaxy Tool Shed: `JASPAR web service`

### motif_scan_sequence

**Scan sequences with a PWM** — Score every window against the log-odds matrix and report thresholded sites.

| parameter | kind | default | options |
|---|---|---|---|
| `target` | file | `examples/genome.fa` | — |
| `matrix` | file | `examples/motif.pfm` | — |
| `threshold_percent` | number | `85.0` | — |
| `strand` | choice | `both` | both, plus, minus |
| `head` | int | `300` | — |

Galaxy Tool Shed: `gtr-scanner / MOODS / FIMO`

### motif_summary

**Consensus and information content of motifs** — Per-motif width, consensus, total information and degenerate-regex summary.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/motifs_jaspar.txt` | — |
| `max_ic` | int | `2.0` | — |

Galaxy Tool Shed: `jaspar / motif db summary`

### motif_transform

**Reverse-complement / trim / normalise a motif** — Rewrite a matrix in the opposite orientation, trimmed, rounded or row-normalised.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/motif.pfm` | — |
| `operation` | choice | `revcomp` | revcomp, trim_half, round, normalise |
| `pseudocount` | number | `0.0` | — |

Galaxy Tool Shed: `motif tools (MEME suite)`

## gemini

*Gemini* — in Galaxy group *Genomics Toolkits*.

| Tool | What it does | Upstream Galaxy tool |
|---|---|---|
| [`gemini_annotate_regions`](#gemini_annotate_regions) | Intersect a VCF with genes and summarise the coding consequences. | `gemini annotate / VEP` |
| [`gemini_chrom_state`](#gemini_chrom_state) | Bucket regions into repressed, poised and active by their signal level. | `chromHMM-lite / gemini chrom_state` |
| [`vcf_gene_based_association`](#vcf_gene_based_association) | Collapse variants per gene and test case/control burden. | `iuc/gemini` |

### gemini_annotate_regions

**Annotate variants falling in features** — Intersect a VCF with genes and summarise the coding consequences.

| parameter | kind | default | options |
|---|---|---|---|
| `variants` | file | `examples/variants.vcf` | — |
| `annotation` | file | `examples/annotation.gff` | — |
| `max_impact` | choice | `modifier` | high, medium, low, modifier |
| `gene_summary` | bool | `False` | — |

Galaxy Tool Shed: `gemini annotate / VEP`

### gemini_chrom_state

**Coverage-based chromatin state per interval** — Bucket regions into repressed, poised and active by their signal level.

| parameter | kind | default | options |
|---|---|---|---|
| `graph` | file | `examples/coverage.bedgraph` | — |
| `regions` | file | `examples/regions.bed` | — |
| `low` | number | `0.2` | — |
| `high` | number | `0.8` | — |

Galaxy Tool Shed: `chromHMM-lite / gemini chrom_state`

### vcf_gene_based_association

**Gemini gene-based burden test** — Collapse variants per gene and test case/control burden.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/variants.vcf` | — |
| `annotation` | file | `examples/annotation.gff` | — |
| `burden` | choice | `count` | count, damaging_fraction, mean_qual |
| `min_variants` | int | `1` | — |

Galaxy Tool Shed: `iuc/gemini`

## du_novo

*Du Novo (Duplex Consensus)* — in Galaxy group *Genomics Toolkits*.

| Tool | What it does | Upstream Galaxy tool |
|---|---|---|
| [`dunovo_duplex_consensus`](#dunovo_duplex_consensus) | Merge overlapping read pairs into one consensus with per-read quality stats. | `Duplex Consensus (du novo)` |

### dunovo_duplex_consensus

**Duplex consensus from paired reads** — Merge overlapping read pairs into one consensus with per-read quality stats.

| parameter | kind | default | options |
|---|---|---|---|
| `reads1` | file | `examples/reads_1.fastq` | — |
| `reads2` | file | `examples/reads_2.fastq` | — |
| `min_overlap` | int | `8` | — |
| `min_quality` | number | `20.0` | — |

Galaxy Tool Shed: `Duplex Consensus (du novo)`

## iwtomics

*IWTomics* — in Galaxy group *Genomics Toolkits*.

| Tool | What it does | Upstream Galaxy tool |
|---|---|---|
| [`iwtonomics_interval_test`](#iwtonomics_interval_test) | Smooth the profile along the feature order and find intervals differing between groups. | `IWTomics (iwtomics)` |

### iwtonomics_interval_test

**IWTomics-style interval test on profiles** — Smooth the profile along the feature order and find intervals differing between groups.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/counts.tsv` | — |
| `groups` | text | `` | — |
| `stat` | choice | `mean` | mean, max, sum |
| `neighbourhood` | int | `2` | — |
| `permutations` | int | `200` | — |

Galaxy Tool Shed: `IWTomics (iwtomics)`

## presto

*PreSTO* — in Galaxy group *Genomics Toolkits*.

| Tool | What it does | Upstream Galaxy tool |
|---|---|---|
| [`presto_aucell`](#presto_aucell) | Area under the recovery curve of each gene set within a cell's ranked genes. | `AUCell / presto:: aucell` |
| [`presto_marker_auc`](#presto_marker_auc) | Find genes that separate one cluster from the rest by area under the ROC curve. | `presto::marker_genes / AUCell markers` |

### presto_aucell

**AUCell-like gene-set activity scores** — Area under the recovery curve of each gene set within a cell's ranked genes.

| parameter | kind | default | options |
|---|---|---|---|
| `counts` | file | `examples/counts.tsv` | — |
| `gene_sets` | file | `examples/genesets.gmt` | — |
| `cut` | number | `0.05` | — |
| `rank_normalized` | bool | `True` | — |
| `head` | int | `100` | — |

Galaxy Tool Shed: `AUCell / presto:: aucell`

### presto_marker_auc

**Per-cluster marker ranking by AUC** — Find genes that separate one cluster from the rest by area under the ROC curve.

| parameter | kind | default | options |
|---|---|---|---|
| `counts` | file | `examples/counts.tsv` | — |
| `cluster_column` | text | `` | — |
| `metadata` | file | `examples/phenotypes.tsv` | — |
| `min_auc` | number | `0.7` | — |
| `head` | int | `100` | — |

Galaxy Tool Shed: `presto::marker_genes / AUCell markers`

## planttribes

*PlantTribes* — in Galaxy group *Genomics Toolkits*.

| Tool | What it does | Upstream Galaxy tool |
|---|---|---|
| [`planttribes_align_contigs`](#planttribes_align_contigs) | Simulated enzyme digestion of contigs into plausibly-sized loci. | `ipyrad / UNEAF1 demultiplex` |
| [`planttribes_enzyme_cuts`](#planttribes_enzyme_cuts) | Count how often each enzyme cuts every sequence in a RAD library. | `TASSEL / ipyrad --enzymes` |

### planttribes_align_contigs

**Digest contigs into RAD loci** — Simulated enzyme digestion of contigs into plausibly-sized loci.

| parameter | kind | default | options |
|---|---|---|---|
| `records` | file | `examples/genes.fa` | — |
| `enzyme` | choice | `EcoRI` | EcoRI, HindIII, BamHI, PstI |
| `min_len` | int | `30` | — |
| `sticky_end` | bool | `True` | — |

Galaxy Tool Shed: `ipyrad / UNEAF1 demultiplex`

### planttribes_enzyme_cuts

**Restriction enzyme cut frequencies** — Count how often each enzyme cuts every sequence in a RAD library.

| parameter | kind | default | options |
|---|---|---|---|
| `records` | file | `examples/genes.fa` | — |
| `enzymes` | multi | `['EcoRI', 'HindIII', 'SbfI']` | EcoRI, HindIII, BamHI, SbfI, EcoT22I, PstI, NsiI |
| `min_cut` | int | `0` | — |

Galaxy Tool Shed: `TASSEL / ipyrad --enzymes`

## sccaf

*sccaf* — in Galaxy group *Genomics Toolkits*.

| Tool | What it does | Upstream Galaxy tool |
|---|---|---|
| [`sccaf_split_contigs`](#sccaf_split_contigs) | Break long sequences into fixed-size fragments with sequential ids. | `sccaf split_contigs` |

### sccaf_split_contigs

**Split contigs into fragment files** — Break long sequences into fixed-size fragments with sequential ids.

| parameter | kind | default | options |
|---|---|---|---|
| `contigs` | file | `examples/genome.fa` | — |
| `max_bp` | int | `500` | — |
| `keep_order` | bool | `True` | — |

Galaxy Tool Shed: `sccaf split_contigs`

## seurat

*Seurat* — in Galaxy group *Genomics Toolkits*.

| Tool | What it does | Upstream Galaxy tool |
|---|---|---|
| [`seurat_project_onto_pca`](#seurat_project_onto_pca) | Score query cells on the reference PCA space without recomputing it. | `Seurat ProjectDim / Azimuth mapping` |
| [`seurat_summary_report`](#seurat_summary_report) | Cells, genes, layers and per-cluster medians of an object in one table. | `Seurat object summary` |

### seurat_project_onto_pca

**Project new data onto existing loadings** — Score query cells on the reference PCA space without recomputing it.

| parameter | kind | default | options |
|---|---|---|---|
| `reference` | file | `examples/counts.tsv` | — |
| `query` | file | `examples/counts.tsv` | — |
| `components` | int | `2` | — |
| `standardize` | bool | `True` | — |
| `head` | int | `100` | — |

Galaxy Tool Shed: `Seurat ProjectDim / Azimuth mapping`

### seurat_summary_report

**Seurat-style summary of an object** — Cells, genes, layers and per-cluster medians of an object in one table.

| parameter | kind | default | options |
|---|---|---|---|
| `counts` | file | `examples/counts.tsv` | — |
| `metadata` | file | `examples/phenotypes.tsv` | — |
| `cluster_column` | text | `` | — |
| `top_genes` | int | `3` | — |

Galaxy Tool Shed: `Seurat object summary`

## monocle3

*Monocle3* — in Galaxy group *Genomics Toolkits*.

| Tool | What it does | Upstream Galaxy tool |
|---|---|---|
| [`monocle3_correlated_genes`](#monocle3_correlated_genes) | Rank genes by how monotonically they change along a pseudotime vector. | `monocle3 graph_test / tradeSeq` |
| [`monocle3_pseudotime`](#monocle3_pseudotime) | Order single cells along a nearest-neighbour graph and report pseudotime. | `monocle3 learn_graph / order_cells` |

### monocle3_correlated_genes

**Genes correlated with pseudotime** — Rank genes by how monotonically they change along a pseudotime vector.

| parameter | kind | default | options |
|---|---|---|---|
| `counts` | file | `examples/counts.tsv` | — |
| `pseudotime_column` | text | `` | — |
| `method` | choice | `spearman` | spearman, pearson |
| `min_abs` | number | `0.0` | — |
| `head` | int | `100` | — |

Galaxy Tool Shed: `monocle3 graph_test / tradeSeq`

### monocle3_pseudotime

**Order cells along a principal graph** — Order single cells along a nearest-neighbour graph and report pseudotime.

| parameter | kind | default | options |
|---|---|---|---|
| `counts` | file | `examples/counts.tsv` | — |
| `neighbours` | int | `4` | — |
| `start` | choice | `first` | first, highest_loading, lowest_total |
| `components` | int | `2` | — |
| `backtrack` | bool | `True` | — |

Galaxy Tool Shed: `monocle3 learn_graph / order_cells`

## mimodd

*MiModD* — in Galaxy group *Genomics Toolkits*.

| Tool | What it does | Upstream Galaxy tool |
|---|---|---|
| [`mimodd_mirna_hairpins`](#mimodd_mirna_hairpins) | Stem-loop prediction with free energy and mature-arm extraction. | `MiModD / RNAfold` |

### mimodd_mirna_hairpins

**miRNA hairpin detection** — Stem-loop prediction with free energy and mature-arm extraction.

| parameter | kind | default | options |
|---|---|---|---|
| `records` | file | `examples/genes.fa` | — |
| `min_mfe` | number | `-3.0` | — |
| `min_stem` | int | `8` | — |
| `loop_min` | number | `3` | — |

Galaxy Tool Shed: `MiModD / RNAfold`

## virology

*Virology* — in Galaxy group *Domain Tools*.

| Tool | What it does | Upstream Galaxy tool |
|---|---|---|
| [`viro_consensus_from_variants`](#viro_consensus_from_variants) | Apply called variants to a reference to obtain a consensus sequence. | `bcftools consensus / iVar consensus` |
| [`viro_coverage_by_target`](#viro_coverage_by_target) | Mean and minimum depth over each target region of a viral panel. | `ivar trim / nCoV-2019 depth` |
| [`viro_drift_distances`](#viro_drift_distances) | Pairwise divergence between isolates plus a simple drift cluster label. | `Nextclade / phydrift` |
| [`viro_variant_load`](#viro_variant_load) | Per-site allele balance and the intra-host variant load of a sample. | `lofreq report / ivar variants` |

### viro_consensus_from_variants

**Build a consensus genome from a VCF** — Apply called variants to a reference to obtain a consensus sequence.

| parameter | kind | default | options |
|---|---|---|---|
| `reference` | file | `examples/genome.fa` | — |
| `variants` | file | `examples/variants.vcf` | — |
| `ambiguity` | choice | `majority` | majority, iupac, reference |
| `min_quality` | number | `0.0` | — |

Galaxy Tool Shed: `bcftools consensus / iVar consensus`

### viro_coverage_by_target

**Read coverage per genomic target** — Mean and minimum depth over each target region of a viral panel.

| parameter | kind | default | options |
|---|---|---|---|
| `alignment` | file | `examples/alignments.sam` | — |
| `targets` | file | `examples/regions.bed` | — |
| `min_depth` | number | `1.0` | — |

Galaxy Tool Shed: `ivar trim / nCoV-2019 depth`

### viro_drift_distances

**Antigenic drift matrix between isolates** — Pairwise divergence between isolates plus a simple drift cluster label.

| parameter | kind | default | options |
|---|---|---|---|
| `proteins` | file | `examples/proteins.faa` | — |
| `model` | choice | `identity` | identity, p, k2p |
| `cluster_cutoff` | number | `0.05` | — |

Galaxy Tool Shed: `Nextclade / phydrift`

### viro_variant_load

**Allele fractions and variant load** — Per-site allele balance and the intra-host variant load of a sample.

| parameter | kind | default | options |
|---|---|---|---|
| `variants` | file | `examples/variants.vcf` | — |
| `min_depth` | int | `0` | — |
| `min_fraction` | number | `0.0` | — |

Galaxy Tool Shed: `lofreq report / ivar variants`

## metagenomic_analysis

*Metagenomic Analysis* — in Galaxy group *Domain Tools*.

| Tool | What it does | Upstream Galaxy tool |
|---|---|---|
| [`meta_abundance_at_rank`](#meta_abundance_at_rank) | Sum sample abundances per taxon at a chosen rank. | `phylosef tax_glom / qiime agglomerate` |
| [`meta_alpha_diversity`](#meta_alpha_diversity) | Richness and evenness per sample, optionally rarefied to a common depth. | `vegan::diversity / phyloseq estimate richness` |
| [`meta_beta_brays`](#meta_beta_brays) | Sample-by-sample community dissimilarity matrix as a heatmap. | `vegan::vegdist / phyloseq distance` |
| [`meta_chimera_screen`](#meta_chimera_screen) | Detect sequences whose halves match two different parents better than themselves. | `vsearch --uchime-denovo / mothur chimera` |
| [`meta_count_table_summary`](#meta_count_table_summary) | Library sizes, feature counts and singleton fractions of an abundance table. | `qiime summarize table / phyloseq smry` |
| [`meta_denoise_asvs`](#meta_denoise_asvs) | Collapse sequences that differ by a few substitutions into one ASV. | `DADA2 / deblur / mothur pre.cluster` |
| [`meta_filter_by_lineage`](#meta_filter_by_lineage) | Subset a FASTA by matching lineage text, reporting what was kept. | `qiime taxa filter-table / phyloseq prune` |
| [`meta_kraken_report`](#meta_kraken_report) | Per-taxon percentages at every rank, like a kraken report. | `kraken-report / pavian` |
| [`meta_lca_of_hits`](#meta_lca_of_hits) | Assign each query to the lowest common ancestor of its best hits. | `MEGAN LCA / kraken2 classify` |
| [`meta_ordination_mds`](#meta_ordination_mds) | Two-dimensional ordination of community samples with a stress value. | `vegan::metaMDS / cmdscale` |
| [`meta_prevalence_filter`](#meta_prevalence_filter) | Compare feature prevalence between control and real samples with Fisher's test. | `decontam prevalence` |
| [`meta_rarefaction`](#meta_rarefaction) | Observed richness as a function of sequencing depth for every sample. | `vegan::specpool / iNEXT` |
| [`meta_taxa_summary_by_rank`](#meta_taxa_summary_by_rank) | How many sequences fall into each taxon at the chosen rank. | `kraken-report / phylosep tax_table` |

### meta_abundance_at_rank

**Abundance table aggregated to a taxonomic rank** — Sum sample abundances per taxon at a chosen rank.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/counts.tsv` | — |
| `taxonomy` | file | `examples/taxmap.tsv` | — |
| `rank` | choice | `genus` | kingdom, phylum, class, order, family, genus, species |
| `relative` | bool | `True` | — |

Galaxy Tool Shed: `phylosef tax_glom / qiime agglomerate`

### meta_alpha_diversity

**Within-sample diversity indices** — Richness and evenness per sample, optionally rarefied to a common depth.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/counts.tsv` | — |
| `index` | choice | `shannon` | shannon, simpson, inverse_simpson, chao1, gini, pielou |
| `head` | int | `50` | — |

Galaxy Tool Shed: `vegan::diversity / phyloseq estimate richness`

### meta_beta_brays

**Bray-Curtis dissimilarity between samples** — Sample-by-sample community dissimilarity matrix as a heatmap.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/counts.tsv` | — |
| `relativise` | bool | `True` | — |
| `metric` | choice | `bray_curtis` | bray_curtis, jaccard, euclidean |

Galaxy Tool Shed: `vegan::vegdist / phyloseq distance`

### meta_chimera_screen

**Flag chimeric sequences** — Detect sequences whose halves match two different parents better than themselves.

| parameter | kind | default | options |
|---|---|---|---|
| `records` | file | `examples/genes.fa` | — |
| `min_gain` | int | `2` | — |

Galaxy Tool Shed: `vsearch --uchime-denovo / mothur chimera`

### meta_count_table_summary

**OTU / feature table summary** — Library sizes, feature counts and singleton fractions of an abundance table.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/counts.tsv` | — |
| `per_sample` | bool | `True` | — |

Galaxy Tool Shed: `qiime summarize table / phyloseq smry`

### meta_denoise_asvs

**Denoise sequences by merging near-identical variants** — Collapse sequences that differ by a few substitutions into one ASV.

| parameter | kind | default | options |
|---|---|---|---|
| `records` | file | `examples/genes.fa` | — |
| `mismatches` | int | `1` | — |
| `min_size` | int | `1` | — |

Galaxy Tool Shed: `DADA2 / deblur / mothur pre.cluster`

### meta_filter_by_lineage

**Keep or drop sequences by taxonomic lineage** — Subset a FASTA by matching lineage text, reporting what was kept.

| parameter | kind | default | options |
|---|---|---|---|
| `records` | file | `examples/genes.fa` | — |
| `taxonomy` | file | `examples/taxmap.tsv` | — |
| `pattern` | text | `Proteobacteria` | — |
| `mode` | choice | `keep` | keep, remove |

Galaxy Tool Shed: `qiime taxa filter-table / phyloseq prune`

### meta_kraken_report

**Kraken-style hierarchical report** — Per-taxon percentages at every rank, like a kraken report.

| parameter | kind | default | options |
|---|---|---|---|
| `taxonomy` | file | `examples/taxmap.tsv` | — |
| `min_percent` | int | `0.0` | — |
| `cumulative` | bool | `True` | — |

Galaxy Tool Shed: `kraken-report / pavian`

### meta_lca_of_hits

**Lowest common ancestor from BLAST-like hits** — Assign each query to the lowest common ancestor of its best hits.

| parameter | kind | default | options |
|---|---|---|---|
| `query_records` | file | `examples/genes.fa` | — |
| `subject_records` | file | `examples/proteins.faa` | — |
| `taxonomy` | file | `examples/taxmap.tsv` | — |
| `max_hits` | int | `5` | — |
| `scoring` | choice | `blastn` | blastn, blastp |

Galaxy Tool Shed: `MEGAN LCA / kraken2 classify`

### meta_ordination_mds

**MDS / PCoA ordination of a distance matrix** — Two-dimensional ordination of community samples with a stress value.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/counts.tsv` | — |
| `metric` | choice | `bray_curtis` | bray_curtis, euclidean |
| `colour_column` | text | `` | — |

Galaxy Tool Shed: `vegan::metaMDS / cmdscale`

### meta_prevalence_filter

**Remove features seen only in negative controls** — Compare feature prevalence between control and real samples with Fisher's test.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/counts.tsv` | — |
| `blank_columns` | text | `` | — |
| `min_fold` | number | `2.0` | — |
| `action` | choice | `report` | report, keep |

Galaxy Tool Shed: `decontam prevalence`

### meta_rarefaction

**Rarefaction curves of observed features** — Observed richness as a function of sequencing depth for every sample.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/counts.tsv` | — |
| `steps` | int | `12` | — |

Galaxy Tool Shed: `vegan::specpool / iNEXT`

### meta_taxa_summary_by_rank

**Taxon counts at a lineage rank** — How many sequences fall into each taxon at the chosen rank.

| parameter | kind | default | options |
|---|---|---|---|
| `taxonomy` | file | `examples/taxmap.tsv` | — |
| `rank` | choice | `genus` | kingdom, phylum, class, order, family, genus, species, 0, 1, 2, 3, 4, 5, 6 |
| `head` | int | `50` | — |
| `bar` | bool | `True` | — |

Galaxy Tool Shed: `kraken-report / phylosep tax_table`

## single_cell

*Single Cell Analysis* — in Galaxy group *Domain Tools*.

| Tool | What it does | Upstream Galaxy tool |
|---|---|---|
| [`monocle_trajectory`](#monocle_trajectory) | Pseudotime from geodesic distance along the principal curve of the embedding. | `Monocle3 order_cells` |
| [`presto_wilcoxon_markers`](#presto_wilcoxon_markers) | Fast rank-sum ranking of features between two groups defined by a label column. | `presto wilcoxauc` |
| [`scanpy_score_genes`](#scanpy_score_genes) | Mean expression of a signature relative to random control genes per cell. | `scanpy score_genes` |
| [`seurat_find_clusters`](#seurat_find_clusters) | PCA, mutual k-NN graph and label propagation into clusters (Louvain-lite). | `Seurat FindNeighbors/FindClusters` |
| [`seurat_marker_genes`](#seurat_marker_genes) | Rank-sum test of each gene inside versus outside a cluster with fold changes. | `Seurat FindAllMarkers / presto` |
| [`seurat_normalise`](#seurat_normalise) | Normalise each cell to a total, log1p transform and return the expression matrix. | `Seurat NormalizeData / SCTransform` |
| [`seurat_qc_metrics`](#seurat_qc_metrics) | QC table per cell with total counts, detected genes and mito fraction plus a filter. | `Seurat dataQC / scDblFinder` |
| [`seurat_umap`](#seurat_umap) | Two-dimensional embedding of cells coloured by cluster index. | `Seurat RunUMAP / scanpy` |

### monocle_trajectory

**Order cells along a pseudotime trajectory** — Pseudotime from geodesic distance along the principal curve of the embedding.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/counts.tsv` | — |
| `components` | int | `2` | — |
| `root_cell` | text | `` | — |
| `head` | int | `200` | — |

Galaxy Tool Shed: `Monocle3 order_cells`

### presto_wilcoxon_markers

**Presto-style Wilcoxon marker ranking** — Fast rank-sum ranking of features between two groups defined by a label column.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/counts.tsv` | — |
| `group_column` | text | `` | — |
| `head` | int | `300` | — |

Galaxy Tool Shed: `presto wilcoxauc`

### scanpy_score_genes

**Score a gene set per cell** — Mean expression of a signature relative to random control genes per cell.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/counts.tsv` | — |
| `sets` | file | `examples/genesets.gmt` | — |
| `ctrl` | int | `50` | — |
| `scale` | bool | `True` | — |
| `head` | int | `200` | — |

Galaxy Tool Shed: `scanpy score_genes`

### seurat_find_clusters

**Cluster cells on a k-NN graph** — PCA, mutual k-NN graph and label propagation into clusters (Louvain-lite).

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/counts.tsv` | — |
| `pcs` | int | `2` | — |
| `resolution` | int | `3` | — |
| `neighbours` | int | `5` | — |
| `seed` | int | `1` | — |

Galaxy Tool Shed: `Seurat FindNeighbors/FindClusters`

### seurat_marker_genes

**Marker genes per cluster** — Rank-sum test of each gene inside versus outside a cluster with fold changes.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/counts.tsv` | — |
| `clusters` | file | `` | — |
| `min_auc` | number | `0.0` | — |
| `head` | int | `300` | — |
| `min_logfc` | number | `0.0` | — |

Galaxy Tool Shed: `Seurat FindAllMarkers / presto`

### seurat_normalise

**SCTransform-lite normalisation of a cell matrix** — Normalise each cell to a total, log1p transform and return the expression matrix.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/counts.tsv` | — |
| `scale_factor` | int | `10000` | — |
| `clip` | int | `0` | — |
| `head` | int | `200` | — |

Galaxy Tool Shed: `Seurat NormalizeData / SCTransform`

### seurat_qc_metrics

**Per-cell QC: counts, genes and mitochondrial share** — QC table per cell with total counts, detected genes and mito fraction plus a filter.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/counts.tsv` | — |
| `mito_prefix` | text | `MT-` | — |
| `min_genes` | number | `1.0` | — |
| `head` | int | `200` | — |

Galaxy Tool Shed: `Seurat dataQC / scDblFinder`

### seurat_umap

**UMAP / t-SNE embedding of cells** — Two-dimensional embedding of cells coloured by cluster index.

| parameter | kind | default | options |
|---|---|---|---|
| `src` | file | `examples/counts.tsv` | — |
| `method` | choice | `umap` | umap, tsne, pca |
| `neighbours` | int | `8` | — |
| `components` | int | `2` | — |

Galaxy Tool Shed: `Seurat RunUMAP / scanpy`

## import/manipulate_sc_data

*Import and Manipulate Single Cell Data* — in Galaxy group *Domain Tools*.

| Tool | What it does | Upstream Galaxy tool |
|---|---|---|
| [`sc_barcode_whitelist`](#sc_barcode_whitelist) | Group near-identical barcodes so low-count neighbours are merged in. | `umi_tools whitelist / kb-tools` |
| [`sc_collate_samples`](#sc_collate_samples) | Concatenate two count tables by cells or by features, keeping provenance. | `Seurat merge / SingleCellExperiment cbind` |
| [`sc_image_to_cells`](#sc_image_to_cells) | Segment bright objects in a micrograph and tabulate their intensities. | `CellProfiler / Scanpy image stats` |
| [`sc_matrix_market_preview`](#sc_matrix_market_preview) | Shape, density and per-row totals of a sparse count matrix. | `10x / Seurat ReadMM` |
| [`sc_matrix_to_long`](#sc_matrix_to_long) | Flatten a sparse count matrix into cell, feature and count columns. | `10x read10x / Seurat as.data.frame` |

### sc_barcode_whitelist

**Collapse barcodes that differ by one base** — Group near-identical barcodes so low-count neighbours are merged in.

| parameter | kind | default | options |
|---|---|---|---|
| `barcodes` | file | `ACGTACGT
ACGTACGA
TTTTGGGG
TTTTGGGA` | — |
| `min_count` | int | `1` | — |
| `mismatch` | int | `1` | — |

Galaxy Tool Shed: `umi_tools whitelist / kb-tools`

### sc_collate_samples

**Merge per-sample count tables into one object** — Concatenate two count tables by cells or by features, keeping provenance.

| parameter | kind | default | options |
|---|---|---|---|
| `table_a` | file | `examples/counts.tsv` | — |
| `table_b` | file | `examples/counts.tsv` | — |
| `sample_a` | text | `A` | — |
| `sample_b` | text | `B` | — |
| `mode` | choice | `cbind` | cbind, rbind |

Galaxy Tool Shed: `Seurat merge / SingleCellExperiment cbind`

### sc_image_to_cells

**Extract per-cell intensities from an image** — Segment bright objects in a micrograph and tabulate their intensities.

| parameter | kind | default | options |
|---|---|---|---|
| `image_file` | file | `examples/cells.pgm` | — |
| `threshold` | number | `0.0` | — |
| `min_area` | int | `4` | — |
| `intensity` | choice | `mean` | mean, max, sum |
| `head` | int | `100` | — |

Galaxy Tool Shed: `CellProfiler / Scanpy image stats`

### sc_matrix_market_preview

**Inspect a Matrix Market (mtx) file** — Shape, density and per-row totals of a sparse count matrix.

| parameter | kind | default | options |
|---|---|---|---|
| `matrix` | file | `%%MatrixMarket matrix coordinate real general
3 4 6
1 1 4.0
1 2 1.0
2 3 7.0
3 1 2.0
3 4 5.0
2 2 3.0` | — |
| `head` | int | `20` | — |

Galaxy Tool Shed: `10x / Seurat ReadMM`

### sc_matrix_to_long

**Sparse coordinates to a long table** — Flatten a sparse count matrix into cell, feature and count columns.

| parameter | kind | default | options |
|---|---|---|---|
| `matrix` | file | `%%MatrixMarket matrix coordinate real general
3 4 6
1 1 4.0
1 2 1.0
2 3 7.0
3 1 2.0
3 4 5.0
2 2 3.0` | — |
| `gene_file` | text | `` | — |
| `cell_file` | text | `` | — |
| `drop_zeros` | bool | `True` | — |
| `head` | int | `300` | — |

Galaxy Tool Shed: `10x read10x / Seurat as.data.frame`

## imaging

*Imaging* — in Galaxy group *Domain Tools*.

| Tool | What it does | Upstream Galaxy tool |
|---|---|---|
| [`image_filter_morphology`](#image_filter_morphology) | Apply one image filter and quantify how it changed the intensities. | `Fiji filters / skimage morphology` |
| [`image_measure_objects`](#image_measure_objects) | Area, perimeter, eccentricity and solidity for each detected object. | `CellProfiler MeasureObjectSizeShape` |
| [`image_projection_zstack`](#image_projection_zstack) | Collapse a stack into one plane and report how the statistics change. | `Fiji Z-project / colocalisation threshold` |
| [`image_quality_summary`](#image_quality_summary) | Mean, dynamic range, saturation and the intensity histogram of an image. | `Fiji / ImageJ summary stats` |
| [`image_radial_profile`](#image_radial_profile) | Mean intensity as a function of distance from the image centre. | `Fiji radial profile / CellProfiler` |
| [`image_texture_features`](#image_texture_features) | Contrast, homogeneity, energy and entropy at one or several scales. | `Ilastik texture / skimage features` |
| [`image_threshold_cells`](#image_threshold_cells) | Global threshold, morphological cleanup and an object count. | `Cellpose / ilastik all-pixels` |
| [`image_watershed_split`](#image_watershed_split) | Split clumped objects by flooding a smoothed intensity surface. | `CellProfiler SplitNuclei / watershed` |

### image_filter_morphology

**Filters and morphological operations** — Apply one image filter and quantify how it changed the intensities.

| parameter | kind | default | options |
|---|---|---|---|
| `image_file` | file | `examples/cells.pgm` | — |
| `operation` | choice | `gaussian` | gaussian, median, sobel, laplacian, dilate, erode, fill_holes, equalize, distance_transform, invert |
| `sigma` | number | `1.0` | — |
| `size` | int | `3` | — |
| `iterations` | int | `1` | — |
| `report_stats` | bool | `True` | — |

Galaxy Tool Shed: `Fiji filters / skimage morphology`

### image_measure_objects

**Morphometry of segmented objects** — Area, perimeter, eccentricity and solidity for each detected object.

| parameter | kind | default | options |
|---|---|---|---|
| `image_file` | file | `examples/cells.pgm` | — |
| `threshold` | number | `0.0` | — |
| `min_area` | int | `8` | — |
| `shape_metric` | choice | `solidity` | solidity, eccentricity, extent |
| `head` | int | `100` | — |

Galaxy Tool Shed: `CellProfiler MeasureObjectSizeShape`

### image_projection_zstack

**Maximum/mean intensity projection** — Collapse a stack into one plane and report how the statistics change.

| parameter | kind | default | options |
|---|---|---|---|
| `image_file` | file | `examples/cells.pgm` | — |
| `mode` | choice | `mean` | max, min, mean, sum |
| `axis` | int | `0` | — |
| `slices` | int | `4` | — |
| `figure` | bool | `True` | — |

Galaxy Tool Shed: `Fiji Z-project / colocalisation threshold`

### image_quality_summary

**Basic image statistics and histogram** — Mean, dynamic range, saturation and the intensity histogram of an image.

| parameter | kind | default | options |
|---|---|---|---|
| `image_file` | file | `examples/cells.pgm` | — |
| `bins` | int | `16` | — |
| `per_channel` | bool | `True` | — |

Galaxy Tool Shed: `Fiji / ImageJ summary stats`

### image_radial_profile

**Radial intensity profile** — Mean intensity as a function of distance from the image centre.

| parameter | kind | default | options |
|---|---|---|---|
| `image_file` | file | `examples/cells.pgm` | — |
| `bins` | int | `16` | — |
| `centre` | choice | `image` | image, brightest_object |
| `normalise` | bool | `True` | — |

Galaxy Tool Shed: `Fiji radial profile / CellProfiler`

### image_texture_features

**Texture descriptors (GLCM and coarseness)** — Contrast, homogeneity, energy and entropy at one or several scales.

| parameter | kind | default | options |
|---|---|---|---|
| `image_file` | file | `examples/cells.pgm` | — |
| `levels` | int | `8` | — |
| `distance` | int | `1` | — |
| `angle` | choice | `0.0` | 0.0, 0.7854, 1.5708 |
| `multi_scale` | bool | `True` | — |

Galaxy Tool Shed: `Ilastik texture / skimage features`

### image_threshold_cells

**Binarise and count objects** — Global threshold, morphological cleanup and an object count.

| parameter | kind | default | options |
|---|---|---|---|
| `image_file` | file | `examples/cells.pgm` | — |
| `method` | choice | `otsu` | otsu, manual, background |
| `value` | number | `0.5` | — |
| `invert` | bool | `False` | — |
| `min_area` | int | `5` | — |
| `show` | bool | `True` | — |

Galaxy Tool Shed: `Cellpose / ilastik all-pixels`

### image_watershed_split

**Watershed splitting of touching objects** — Split clumped objects by flooding a smoothed intensity surface.

| parameter | kind | default | options |
|---|---|---|---|
| `image_file` | file | `examples/cells.pgm` | — |
| `markers` | int | `4` | — |
| `smoothing` | number | `1.5` | — |
| `otsu_seeds` | bool | `True` | — |

Galaxy Tool Shed: `CellProfiler SplitNuclei / watershed`

## proteomics

*Proteomics* — in Galaxy group *Domain Tools*.

| Tool | What it does | Upstream Galaxy tool |
|---|---|---|
| [`proteo_differential_abundance`](#proteo_differential_abundance) | Per-protein fold change with a p-value and significance flags for two groups. | `DEP / OpenMS DiffEnrich` |
| [`proteo_digest_search`](#proteo_digest_search) | Digest proteins in silico and match each observed peptide by mass tolerance. | `MS-Fragger / PEAKS / OpenMS FidoAdapter` |
| [`proteo_disulfide_pairs`](#proteo_disulfide_pairs) | Cysteine inventory and plausible disulfide pairings by spacing and oxidation state. | `DisulfidePrediction / cysteine stats` |
| [`proteo_enzyme_specificity`](#proteo_enzyme_specificity) | Which protease best explains the observed peptide termini and specificity? | `ProteinLysis / digestion evaluation` |
| [`proteo_fdr_from_decoys`](#proteo_fdr_from_decoys) | Rank PSMs by score and report the decoy-based FDR at every threshold. | `Percolator / Prophet PSM FDR` |
| [`proteo_hydrophobicity_profile`](#proteo_hydrophobicity_profile) | Mean hydropathy along each chain with predicted transmembrane windows. | `Hydropathy (Kyte-Doolittle) / TMHMM plot` |
| [`proteo_normalise_median`](#proteo_normalise_median) | Scale every sample so protein intensities are comparable across runs. | `DEP normalisation / MaxQuant` |
| [`proteo_peptide_properties`](#proteo_peptide_properties) | Mass, m/z, pI, net charge, hydrophobicity and instability for each peptide. | `ProtParam / pyMSpec` |
| [`proteo_ptm_sites`](#proteo_ptm_sites) | Sequence-motif based prediction of glycosylation, phosphorylation and oxidation sites. | `Uniprot PTM scan / SeQuest PTM` |
| [`proteo_quant_summary`](#proteo_quant_summary) | Peptide counts, total and unique intensity per protein group. | `MaxQuant proteinGroups / FlashLFQ` |
| [`proteo_sequence_coverage`](#proteo_sequence_coverage) | Share of each protein explained by identified peptides, plus uncovered gaps. | `PeptideShaker / coverage plot` |
| [`proteo_theoretical_spectra`](#proteo_theoretical_spectra) | Digest and list every peptide with masses and m/z for the chosen charges. | `PeptideGenerator / msconvert theoretical` |

### proteo_differential_abundance

**Differential protein abundance between groups** — Per-protein fold change with a p-value and significance flags for two groups.

| parameter | kind | default | options |
|---|---|---|---|
| `abundance` | file | `examples/counts.tsv` | — |
| `groups` | text | `` | — |
| `test` | choice | `ttest` | ttest, welch, mannwhitney |
| `fc_cutoff` | number | `2.0` | — |
| `alpha` | number | `0.05` | — |

Galaxy Tool Shed: `DEP / OpenMS DiffEnrich`

### proteo_digest_search

**Match observed peptides to an in-silico digest** — Digest proteins in silico and match each observed peptide by mass tolerance.

| parameter | kind | default | options |
|---|---|---|---|
| `protein_file` | file | `examples/proteins.faa` | — |
| `spectrum_table` | file | `examples/peptides.tsv` | — |
| `enzyme` | choice | `trypsin` | trypsin, lys-c, arg-c, glu-c, chymotrypsin, none |
| `ppm_tolerance` | number | `50.0` | — |
| `missed_cleavages` | int | `1` | — |
| `min_mass` | number | `300.0` | — |

Galaxy Tool Shed: `MS-Fragger / PEAKS / OpenMS FidoAdapter`

### proteo_disulfide_pairs

**Cysteine and disulfide analysis** — Cysteine inventory and plausible disulfide pairings by spacing and oxidation state.

| parameter | kind | default | options |
|---|---|---|---|
| `protein_file` | file | `examples/proteins.faa` | — |
| `mode` | choice | `pairs` | all_cysteines, pairs |
| `min_spacing` | int | `12` | — |

Galaxy Tool Shed: `DisulfidePrediction / cysteine stats`

### proteo_enzyme_specificity

**Compare observed peptides across enzymes** — Which protease best explains the observed peptide termini and specificity?

| parameter | kind | default | options |
|---|---|---|---|
| `protein_file` | file | `examples/proteins.faa` | — |
| `spectrum_table` | file | `examples/peptides.tsv` | — |
| `missed_cleavages` | int | `2` | — |
| `head` | int | `20` | — |

Galaxy Tool Shed: `ProteinLysis / digestion evaluation`

### proteo_fdr_from_decoys

**Target-decoy false discovery rate** — Rank PSMs by score and report the decoy-based FDR at every threshold.

| parameter | kind | default | options |
|---|---|---|---|
| `spectrum_table` | file | `examples/peptides.tsv` | — |
| `score_column` | text | `` | — |
| `target_prefix` | text | `sp|` | — |
| `decoy_prefix` | text | `decoy_` | — |
| `fdr_cutoff` | number | `0.01` | — |

Galaxy Tool Shed: `Percolator / Prophet PSM FDR`

### proteo_hydrophobicity_profile

**Sliding hydrophobicity profile** — Mean hydropathy along each chain with predicted transmembrane windows.

| parameter | kind | default | options |
|---|---|---|---|
| `protein_file` | file | `examples/proteins.faa` | — |
| `window` | int | `19` | — |
| `scale` | choice | `kyte_doolittle` | kyte_doolittle, eisenger, wilson |
| `tm_cutoff` | number | `1.6` | — |
| `protein_id` | text | `` | — |

Galaxy Tool Shed: `Hydropathy (Kyte-Doolittle) / TMHMM plot`

### proteo_normalise_median

**Median-intensity normalisation of a protein table** — Scale every sample so protein intensities are comparable across runs.

| parameter | kind | default | options |
|---|---|---|---|
| `abundance` | file | `examples/counts.tsv` | — |
| `mode` | choice | `median` | median, total, rms, none |
| `log2` | number | `True` | — |

Galaxy Tool Shed: `DEP normalisation / MaxQuant`

### proteo_peptide_properties

**Physicochemical properties of peptides** — Mass, m/z, pI, net charge, hydrophobicity and instability for each peptide.

| parameter | kind | default | options |
|---|---|---|---|
| `spectrum_table` | file | `examples/peptides.tsv` | — |
| `sequence_column` | text | `` | — |
| `charge` | int | `2` | — |
| `with_mods` | bool | `False` | — |
| `head` | int | `200` | — |

Galaxy Tool Shed: `ProtParam / pyMSpec`

### proteo_ptm_sites

**Potential modification sites in proteins** — Sequence-motif based prediction of glycosylation, phosphorylation and oxidation sites.

| parameter | kind | default | options |
|---|---|---|---|
| `protein_file` | file | `examples/proteins.faa` | — |
| `modifications` | multi | `['N-glycosylation', 'phospho-ST', 'oxidation-M']` | N-glycosylation, phospho-ST, phospho-Y, acetylation-K, ubiquitination-K, oxidation-M, amidation-Cterm |
| `context` | bool | `True` | — |

Galaxy Tool Shed: `Uniprot PTM scan / SeQuest PTM`

### proteo_quant_summary

**Quantitation summary per protein** — Peptide counts, total and unique intensity per protein group.

| parameter | kind | default | options |
|---|---|---|---|
| `spectrum_table` | file | `examples/peptides.tsv` | — |
| `protein_column` | text | `` | — |
| `intensity_column` | text | `` | — |
| `share` | bool | `True` | — |
| `head` | int | `100` | — |

Galaxy Tool Shed: `MaxQuant proteinGroups / FlashLFQ`

### proteo_sequence_coverage

**Sequence coverage from identified peptides** — Share of each protein explained by identified peptides, plus uncovered gaps.

| parameter | kind | default | options |
|---|---|---|---|
| `protein_file` | file | `examples/proteins.faa` | — |
| `spectrum_table` | file | `examples/peptides.tsv` | — |
| `protein_column` | text | `` | — |
| `with_mods` | bool | `False` | — |

Galaxy Tool Shed: `PeptideShaker / coverage plot`

### proteo_theoretical_spectra

**Theoretical peptide m/z list** — Digest and list every peptide with masses and m/z for the chosen charges.

| parameter | kind | default | options |
|---|---|---|---|
| `protein_file` | file | `examples/proteins.faa` | — |
| `enzyme` | choice | `trypsin` | trypsin, lys-c, arg-c, glu-c, chymotrypsin, none |
| `min_len` | int | `7` | — |
| `max_len` | int | `30` | — |
| `charges` | multi | `['1', '2']` | 1, 2, 3 |
| `head` | int | `300` | — |

Galaxy Tool Shed: `PeptideGenerator / msconvert theoretical`

## metabolomics

*Metabolomics* — in Galaxy group *Domain Tools*.

| Tool | What it does | Upstream Galaxy tool |
|---|---|---|
| [`metab_adduct_annotate`](#metab_adduct_annotate) | Candidate neutral masses and formulas for each observed m/z within a ppm window. | `SIRIUS / MZmine isotope pattern` |
| [`metab_blank_subtraction`](#metab_blank_subtraction) | Remove background signal measured in blank injections from every sample. | `MetaboAnalyst remove blanks / bls` |
| [`metab_feature_summary`](#metab_feature_summary) | Intensity statistics, detection rate and dynamic range of every feature. | `XCMS / MZmine feature list` |
| [`metab_missing_values`](#metab_missing_values) | Per-feature missingness with an imputed value matrix for downstream stats. | `MetaboAnalyst imputation / minimol` |
| [`metab_normalise_methods`](#metab_normalise_methods) | Median, quantile, TIC or internal-standard scaling with before/after diagnostics. | `MetaboAnalyst normalisation / normq` |
| [`metab_qc_quality_control`](#metab_qc_quality_control) | Which features are stable in QC injections, with drift diagnostics. | `MetaboAnalyst quality control / osram` |
| [`metab_rt_alignment`](#metab_rt_alignment) | Estimate and remove a linear RT drift between two LC-MS runs. | `XCMS correctionRt / MZmine align` |

### metab_adduct_annotate

**Annotate m/z with adduct formulas** — Candidate neutral masses and formulas for each observed m/z within a ppm window.

| parameter | kind | default | options |
|---|---|---|---|
| `features` | file | `examples/counts.tsv` | — |
| `mz_column` | text | `` | — |
| `adducts` | multi | `['[M+H]+', '[M-H]-']` | [M+H]+, [M+Na]+, [M+K]+, [M-H]-, [M+NH4]+, [M+2H]2+, [M-H2O+H]+ |
| `ppm` | number | `10.0` | — |
| `head` | int | `200` | — |

Galaxy Tool Shed: `SIRIUS / MZmine isotope pattern`

### metab_blank_subtraction

**Subtract process blanks** — Remove background signal measured in blank injections from every sample.

| parameter | kind | default | options |
|---|---|---|---|
| `features` | file | `examples/counts.tsv` | — |
| `blanks` | file | `examples/counts.tsv` | — |
| `factor` | number | `1.0` | — |
| `floor_zero` | bool | `True` | — |
| `min_ratio` | number | `3.0` | — |

Galaxy Tool Shed: `MetaboAnalyst remove blanks / bls`

### metab_feature_summary

**Peak table feature summary** — Intensity statistics, detection rate and dynamic range of every feature.

| parameter | kind | default | options |
|---|---|---|---|
| `features` | file | `examples/counts.tsv` | — |
| `mz_column` | text | `` | — |
| `rt_column` | text | `` | — |
| `head` | int | `200` | — |

Galaxy Tool Shed: `XCMS / MZmine feature list`

### metab_missing_values

**Missing-value patterns and imputation** — Per-feature missingness with an imputed value matrix for downstream stats.

| parameter | kind | default | options |
|---|---|---|---|
| `features` | file | `examples/counts.tsv` | — |
| `group_column` | choice | `none` | none |
| `method` | choice | `half_min` | half_min, zero, knn, mean, none |
| `detection_cut` | number | `0.5` | — |
| `head` | int | `200` | — |

Galaxy Tool Shed: `MetaboAnalyst imputation / minimol`

### metab_normalise_methods

**Normalise a feature table** — Median, quantile, TIC or internal-standard scaling with before/after diagnostics.

| parameter | kind | default | options |
|---|---|---|---|
| `features` | file | `examples/counts.tsv` | — |
| `samples` | text | `` | — |
| `method` | choice | `quantile` | median, quantile, internal_standard, tic, log, none |
| `standard_feature` | text | `` | — |
| `head` | int | `200` | — |

Galaxy Tool Shed: `MetaboAnalyst normalisation / normq`

### metab_qc_quality_control

**QC-based performance metrics** — Which features are stable in QC injections, with drift diagnostics.

| parameter | kind | default | options |
|---|---|---|---|
| `features` | file | `examples/counts.tsv` | — |
| `qc_columns` | text | `` | — |
| `cv_threshold` | number | `30.0` | — |
| `drift_column` | text | `` | — |
| `head` | int | `200` | — |

Galaxy Tool Shed: `MetaboAnalyst quality control / osram`

### metab_rt_alignment

**Align retention times between two runs** — Estimate and remove a linear RT drift between two LC-MS runs.

| parameter | kind | default | options |
|---|---|---|---|
| `query` | file | `examples/peptides.tsv` | — |
| `reference` | file | `examples/peptides.tsv` | — |
| `query_rt` | text | `runtime` | — |
| `reference_rt` | text | `runtime` | — |
| `query_key` | text | `peptide` | — |
| `reference_key` | text | `peptide` | — |
| `max_shift` | number | `5.0` | — |

Galaxy Tool Shed: `XCMS correctionRt / MZmine align`

## chemicaltoolbox

*ChemicalToolBox* — in Galaxy group *Domain Tools*.

| Tool | What it does | Upstream Galaxy tool |
|---|---|---|
| [`chem_compound_hash`](#chem_compound_hash) | Deterministic connectivity hash so duplicated compounds can be collapsed. | `InChIKey generation` |
| [`chem_library_diversity`](#chem_library_diversity) | Collapse structures to their ring or formula scaffold and score diversity. | `Bemis-Murcko scaffolds / infomax` |
| [`chem_molecular_property_filter`](#chem_molecular_property_filter) | Apply drug-likeness windows to a SMILES list and emit the survivors. | `SwissTargetPrediction filtering / medchem rules` |
| [`chem_reactant_mass_balance`](#chem_reactant_mass_balance) | Compute mmol, mass and the limiting reagent for each reaction component. | `Synthia / Amide reaction calculator` |
| [`chem_similarity_screen`](#chem_similarity_screen) | Rank a compound library by fingerprint similarity to a query structure. | `RDKit similarity search / SwissScreen` |
| [`chem_smiles_descriptors`](#chem_smiles_descriptors) | Formula, mass, ring and H-bond descriptors computed directly from SMILES. | `RDKit Descriptors / cdk descriptors` |
| [`chem_substructure_search`](#chem_substructure_search) | Match a symbolic pattern against every structure in a library. | `RDKit SubstructSearch / DAYLIGHT SMARTS` |

### chem_compound_hash

**Canonical hash / InChIKey-style identifier** — Deterministic connectivity hash so duplicated compounds can be collapsed.

| parameter | kind | default | options |
|---|---|---|---|
| `library` | file | `CC(=O)Oc1ccccc1C(=O)O
CCOc1ccccc1` | — |
| `key_length` | int | `14` | — |
| `group_isomers` | bool | `True` | — |

Galaxy Tool Shed: `InChIKey generation`

### chem_library_diversity

**Scaffold diversity of a compound library** — Collapse structures to their ring or formula scaffold and score diversity.

| parameter | kind | default | options |
|---|---|---|---|
| `library` | file | `CC(=O)Oc1ccccc1C(=O)O
c1ccccc1O
c1ccncc1
C1CCNCC1
CC(=O)Nc1ccccc1` | — |
| `scaffold_definition` | choice | `rings` | rings, heavy_atom_graph, formula |
| `top` | int | `20` | — |

Galaxy Tool Shed: `Bemis-Murcko scaffolds / infomax`

### chem_molecular_property_filter

**Filter a compound list by property windows** — Apply drug-likeness windows to a SMILES list and emit the survivors.

| parameter | kind | default | options |
|---|---|---|---|
| `library` | file | `CC(=O)Oc1ccccc1C(=O)O
c1ccccc1O
CCN` | — |
| `min_mw` | number | `0.0` | — |
| `max_mw` | number | `1000.0` | — |
| `max_donors` | number | `10` | — |
| `max_acceptors` | number | `15` | — |
| `max_rings` | int | `10` | — |
| `min_logp` | number | `-5.0` | — |
| `max_logp` | number | `10.0` | — |
| `keep_only` | bool | `True` | — |

Galaxy Tool Shed: `SwissTargetPrediction filtering / medchem rules`

### chem_reactant_mass_balance

**Reaction mass balance from a stoichiometry table** — Compute mmol, mass and the limiting reagent for each reaction component.

| parameter | kind | default | options |
|---|---|---|---|
| `stoichiometry` | file | `examples/counts.tsv` | — |
| `compound_column` | text | `` | — |
| `mw_column` | text | `` | — |
| `equiv_column` | text | `` | — |
| `scale` | number | `1.0` | — |
| `limiting` | choice | `auto` | auto, first_row |

Galaxy Tool Shed: `Synthia / Amide reaction calculator`

### chem_similarity_screen

**Tanimoto similarity screen against a query** — Rank a compound library by fingerprint similarity to a query structure.

| parameter | kind | default | options |
|---|---|---|---|
| `query_smiles` | text | `CC(=O)Oc1ccccc1C(=O)O` | — |
| `library` | file | `CC(=O)Oc1ccccc1C(=O)O
c1ccccc1O
CN1C=NC2=C1C(=O)N(C)C(=O)N2C
CCO
CC(=O)O` | — |
| `head` | int | `50` | — |

Galaxy Tool Shed: `RDKit similarity search / SwissScreen`

### chem_smiles_descriptors

**Descriptors from SMILES strings** — Formula, mass, ring and H-bond descriptors computed directly from SMILES.

| parameter | kind | default | options |
|---|---|---|---|
| `smiles_list` | file | `CC(=O)Oc1ccccc1C(=O)O
CN1C=NC2=C1C(=O)N(C)C(=O)N2C
c1ccccc1` | — |
| `extra` | bool | `True` | — |
| `head` | int | `200` | — |

Galaxy Tool Shed: `RDKit Descriptors / cdk descriptors`

### chem_substructure_search

**Substructure and SMARTS-lite search** — Match a symbolic pattern against every structure in a library.

| parameter | kind | default | options |
|---|---|---|---|
| `library` | file | `CC(=O)Oc1ccccc1C(=O)O
c1ccccc1O
CCN
c1ccncc1` | — |
| `smarts` | text | `c1ccccc1` | — |
| `count_only` | bool | `False` | — |
| `head` | int | `200` | — |

Galaxy Tool Shed: `RDKit SubstructSearch / DAYLIGHT SMARTS`

## pharmacology

*Pharmacology* — in Galaxy group *Domain Tools*.

| Tool | What it does | Upstream Galaxy tool |
|---|---|---|
| [`pharm_bliss_synergy`](#pharm_bliss_synergy) | Compare the observed combination effect with the Bliss expectation per row. | `CombroScore / synergyfinder` |
| [`pharm_hill_curve`](#pharm_hill_curve) | Fit EC50/IC50, Hill slope and the curve parameters to a dose-response series. | `GraphPad Prism / drc::drm` |
| [`pharm_ic50_panel`](#pharm_ic50_panel) | Interpolate the dose that reaches a given inhibition for every compound. | `Prism curve fitting panel / activity-based profiling` |
| [`pharm_pk_compartment_model`](#pharm_pk_compartment_model) | Fit a simple IV bolus (or two-compartment) curve and report the rate constants. | `nlmixr / Phoenix WinNonlin model fit` |
| [`pharm_pk_noncompartmental`](#pharm_pk_noncompartmental) | Cmax, Tmax, AUC by the linear-log trapezoid rule, clearance and half-life. | `PKSolver / noncompartamental analysis` |

### pharm_bliss_synergy

**Bliss and Loewe combination scores** — Compare the observed combination effect with the Bliss expectation per row.

| parameter | kind | default | options |
|---|---|---|---|
| `combination` | file | `examples/counts.tsv` | — |
| `drug_a_column` | text | `` | — |
| `drug_b_column` | text | `` | — |
| `combination_column` | text | `` | — |
| `model` | choice | `both` | bliss, loewe, both |
| `cutoff` | number | `0.1` | — |

Galaxy Tool Shed: `CombroScore / synergyfinder`

### pharm_hill_curve

**Hill / dose-response fit** — Fit EC50/IC50, Hill slope and the curve parameters to a dose-response series.

| parameter | kind | default | options |
|---|---|---|---|
| `dose_response` | file | `examples/counts.tsv` | — |
| `dose_column` | text | `` | — |
| `response_column` | text | `` | — |
| `fit` | choice | `hill` | hill, logistic, linear |
| `top` | number | `0.0` | — |
| `bottom` | number | `0.0` | — |

Galaxy Tool Shed: `GraphPad Prism / drc::drm`

### pharm_ic50_panel

**Per-compound IC50 from a dose matrix** — Interpolate the dose that reaches a given inhibition for every compound.

| parameter | kind | default | options |
|---|---|---|---|
| `screen` | file | `examples/counts.tsv` | — |
| `doses` | text | `` | — |
| `inhibition_target` | number | `50.0` | — |
| `top` | number | `0.0` | — |
| `head` | int | `200` | — |

Galaxy Tool Shed: `Prism curve fitting panel / activity-based profiling`

### pharm_pk_compartment_model

**One- and two-compartment model fit** — Fit a simple IV bolus (or two-compartment) curve and report the rate constants.

| parameter | kind | default | options |
|---|---|---|---|
| `pk` | file | `examples/counts.tsv` | — |
| `time_column` | text | `` | — |
| `concentration_column` | text | `` | — |
| `model` | choice | `one` | one, two |
| `dose` | number | `100.0` | — |
| `iterations` | int | `400` | — |

Galaxy Tool Shed: `nlmixr / Phoenix WinNonlin model fit`

### pharm_pk_noncompartmental

**Non-compartmental PK analysis** — Cmax, Tmax, AUC by the linear-log trapezoid rule, clearance and half-life.

| parameter | kind | default | options |
|---|---|---|---|
| `pk` | file | `examples/counts.tsv` | — |
| `time_column` | text | `` | — |
| `concentration_column` | text | `` | — |
| `dose_route` | choice | `iv` | iv, extravascular |
| `dose` | number | `100.0` | — |
| `body_weight` | number | `70.0` | — |
| `lambda_z_points` | number | `3` | — |

Galaxy Tool Shed: `PKSolver / noncompartamental analysis`

## multiomics

*Multi-Omics* — in Galaxy group *Domain Tools*.

| Tool | What it does | Upstream Galaxy tool |
|---|---|---|
| [`multiomics_cross_correlation`](#multiomics_cross_correlation) | Correlation heatmap between the features of two omics layers. | `mixOmics correlogram / MOFA loadings` |
| [`multiomics_integrated_score`](#multiomics_integrated_score) | Combine two feature tables into one ranked multi-omics score per feature. | `mixOmics DIABLO score / iCluster` |
| [`multiomics_merge_tables`](#multiomics_merge_tables) | Merge expression with another data layer and correlate the shared columns. | `multi-omics integration (MOFA-lite)` |
| [`multiomics_rank_product`](#multiomics_rank_product) | Combine two ranked lists with the rank product to find consistently top features. | `RankProd / tie within replicates` |
| [`multiomics_set_overlap`](#multiomics_set_overlap) | Shared and unique feature names between two or three omics tables. | `mixOmics / clusterProfiler compareSets` |

### multiomics_cross_correlation

**Cross-omics correlation block** — Correlation heatmap between the features of two omics layers.

| parameter | kind | default | options |
|---|---|---|---|
| `table_a` | file | `examples/counts.tsv` | — |
| `table_b` | file | `examples/counts.tsv` | — |
| `method` | choice | `spearman` | pearson, spearman |
| `max_rows` | int | `12` | — |
| `max_cols` | int | `12` | — |

Galaxy Tool Shed: `mixOmics correlogram / MOFA loadings`

### multiomics_integrated_score

**Z-score integration across omics layers** — Combine two feature tables into one ranked multi-omics score per feature.

| parameter | kind | default | options |
|---|---|---|---|
| `table_a` | file | `examples/counts.tsv` | — |
| `table_b` | file | `examples/counts.tsv` | — |
| `weight_a` | choice | `1` | 1, 2, 3 |
| `weight_b` | choice | `1` | 1, 2, 3 |
| `rank_instead` | bool | `False` | — |
| `head` | int | `100` | — |

Galaxy Tool Shed: `mixOmics DIABLO score / iCluster`

### multiomics_merge_tables

**Join two omics tables by feature** — Merge expression with another data layer and correlate the shared columns.

| parameter | kind | default | options |
|---|---|---|---|
| `left` | file | `examples/counts.tsv` | — |
| `right` | file | `examples/phenotypes.tsv` | — |
| `key_left` | text | `` | — |
| `key_right` | text | `` | — |
| `how` | choice | `inner` | inner, left, right, outer |
| `min_abs_corr` | number | `0.0` | — |

Galaxy Tool Shed: `multi-omics integration (MOFA-lite)`

### multiomics_rank_product

**Rank product consensus across datasets** — Combine two ranked lists with the rank product to find consistently top features.

| parameter | kind | default | options |
|---|---|---|---|
| `left` | file | `examples/counts.tsv` | — |
| `right` | file | `examples/counts.tsv` | — |
| `score_left` | text | `` | — |
| `score_right` | text | `` | — |
| `key_left` | text | `` | — |
| `key_right` | text | `` | — |
| `head` | int | `200` | — |

Galaxy Tool Shed: `RankProd / tie within replicates`

### multiomics_set_overlap

**Overlap of feature sets across omics layers** — Shared and unique feature names between two or three omics tables.

| parameter | kind | default | options |
|---|---|---|---|
| `table_a` | file | `examples/counts.tsv` | — |
| `table_b` | file | `examples/counts.tsv` | — |
| `table_c` | file | `` | — |
| `key_column` | text | `` | — |
| `jaccard` | bool | `True` | — |

Galaxy Tool Shed: `mixOmics / clusterProfiler compareSets`

## climate_analysis

*Climate Analysis* — in Galaxy group *Domain Tools*.

| Tool | What it does | Upstream Galaxy tool |
|---|---|---|
| [`climate_growing_degree_days`](#climate_growing_degree_days) | Daily degree-day accumulation above a base temperature, with frost counts. | `climate GDD / degree-days (agronomy)` |
| [`climate_precipitation_indices`](#climate_precipitation_indices) | Wet-day counts, simple precipitation intensity and the longest dry spell. | `climdex / RClimdex` |
| [`climate_spatial_autocorrelation`](#climate_spatial_autocorrelation) | Spatial autocorrelation of a gridded field with a significance test. | `spdep::moran.test / PySAL` |

### climate_growing_degree_days

**Growing degree days and frost risk** — Daily degree-day accumulation above a base temperature, with frost counts.

| parameter | kind | default | options |
|---|---|---|---|
| `daily_weather` | file | `examples/counts.tsv` | — |
| `tmax_column` | text | `` | — |
| `tmin_column` | text | `` | — |
| `base` | number | `10.0` | — |
| `ceiling` | number | `30.0` | — |
| `frost` | bool | `True` | — |

Galaxy Tool Shed: `climate GDD / degree-days (agronomy)`

### climate_precipitation_indices

**Precipitation extremes and dry spells** — Wet-day counts, simple precipitation intensity and the longest dry spell.

| parameter | kind | default | options |
|---|---|---|---|
| `rainfall` | file | `examples/counts.tsv` | — |
| `column` | text | `` | — |
| `wet_day` | number | `1.0` | — |
| `heavy` | number | `20.0` | — |
| `max_dry` | int | `5` | — |

Galaxy Tool Shed: `climdex / RClimdex`

### climate_spatial_autocorrelation

**Moran's I of a spatial field** — Spatial autocorrelation of a gridded field with a significance test.

| parameter | kind | default | options |
|---|---|---|---|
| `raster` | file | `examples/cells.pgm` | — |
| `blocks` | int | `6` | — |
| `neighbours` | number | `1.0` | — |
| `standardise` | choice | `row` | row, none |

Galaxy Tool Shed: `spdep::moran.test / PySAL`

## gis_data_handling

*GIS Data Handling* — in Galaxy group *Domain Tools*.

| Tool | What it does | Upstream Galaxy tool |
|---|---|---|
| [`gis_coordinate_transform`](#gis_coordinate_transform) | Reproject lat/lon into Web-Mercator or a simple UTM-style grid. | `sf::st_transform / pyproj` |
| [`gis_nearest_neighbour`](#gis_nearest_neighbour) | Link each point to its k nearest neighbours with distances. | `sf nearest / geopandas sjoin_nearest` |
| [`gis_raster_zonal_stats`](#gis_raster_zonal_stats) | Summarise a raster inside each zone of a polygon table. | `raster zonal / exactextractr` |
| [`gis_station_distances`](#gis_station_distances) | Great-circle distances between every pair of coordinates in a table. | `sf::st_distance / PostGIS ST_Distance` |

### gis_coordinate_transform

**Convert and validate coordinates** — Reproject lat/lon into Web-Mercator or a simple UTM-style grid.

| parameter | kind | default | options |
|---|---|---|---|
| `coords` | file | `examples/phenotypes.tsv` | — |
| `lat_column` | text | `latitude` | — |
| `lon_column` | text | `longitude` | — |
| `target` | choice | `web_mercator` | decimal_degrees, _UTM_lite, web_mercator |
| `zone_width` | number | `6.0` | — |

Galaxy Tool Shed: `sf::st_transform / pyproj`

### gis_nearest_neighbour

**Nearest site for every record** — Link each point to its k nearest neighbours with distances.

| parameter | kind | default | options |
|---|---|---|---|
| `records` | file | `examples/phenotypes.tsv` | — |
| `name_column` | text | `` | — |
| `lat_column` | text | `latitude` | — |
| `lon_column` | text | `longitude` | — |
| `k` | int | `3` | — |
| `summary` | bool | `True` | — |

Galaxy Tool Shed: `sf nearest / geopandas sjoin_nearest`

### gis_raster_zonal_stats

**Zonal statistics of a raster-like image** — Summarise a raster inside each zone of a polygon table.

| parameter | kind | default | options |
|---|---|---|---|
| `raster` | file | `examples/cells.pgm` | — |
| `zones` | file | `examples/regions.bed` | — |
| `stat` | choice | `mean` | mean, sum, min, max, stdev, fraction_above |
| `cut` | number | `0.5` | — |
| `head` | int | `200` | — |

Galaxy Tool Shed: `raster zonal / exactextractr`

### gis_station_distances

**Distance matrix between sampling sites** — Great-circle distances between every pair of coordinates in a table.

| parameter | kind | default | options |
|---|---|---|---|
| `sites` | file | `examples/phenotypes.tsv` | — |
| `name_column` | text | `` | — |
| `lat_column` | text | `latitude` | — |
| `lon_column` | text | `longitude` | — |
| `head` | int | `200` | — |

Galaxy Tool Shed: `sf::st_distance / PostGIS ST_Distance`

## epigenetics

*Epigenetics* — in Galaxy group *Domain Tools*.

| Tool | What it does | Upstream Galaxy tool |
|---|---|---|
| [`epi_bisulfite_extractor`](#epi_bisulfite_extractor) | Clean per-cytosine table with context, fraction and call for each site. | `bismark --methylation_extractor` |
| [`epi_coverage_track_stats`](#epi_coverage_track_stats) | Reduce a signal track to one value per region for downstream modelling. | `computeMatrix / bigWigStatistics` |
| [`epi_cpg_islands`](#epi_cpg_islands) | Sliding-window detection of CpG-rich, GC-rich regions reported as BED. | `cpgislandfinder / newcpgreport (EMBOSS)` |
| [`epi_dmr_scan`](#epi_dmr_scan) | Call contiguous stretches of cytosines that cross a methylation threshold. | `DSS / methylKit dmr` |
| [`epi_methylation_context_summary`](#epi_methylation_context_summary) | Mean methylation and coverage of cytosines in each sequence context. | `bismark2bedGraph / MethylDackel` |
| [`epi_methylation_per_region`](#epi_methylation_per_region) | Aggregate cytosine-level methylation into per-region betas. | `methylprep / DSS merge` |

### epi_bisulfite_extractor

**Per-cytosine extraction table** — Clean per-cytosine table with context, fraction and call for each site.

| parameter | kind | default | options |
|---|---|---|---|
| `methyl` | file | `examples/methylation.cxcg` | — |
| `context` | choice | `all` | all, CpG, CHG, CHH |
| `min_fraction` | number | `0.0` | — |
| `max_fraction` | number | `1.0` | — |
| `head` | int | `300` | — |

Galaxy Tool Shed: `bismark --methylation_extractor`

### epi_coverage_track_stats

**Coverage track summary per region** — Reduce a signal track to one value per region for downstream modelling.

| parameter | kind | default | options |
|---|---|---|---|
| `graph` | file | `examples/coverage.bedgraph` | — |
| `regions` | file | `examples/regions.bed` | — |
| `mode` | choice | `mean` | mean, median, max, sum |

Galaxy Tool Shed: `computeMatrix / bigWigStatistics`

### epi_cpg_islands

**Find CpG islands in a genome** — Sliding-window detection of CpG-rich, GC-rich regions reported as BED.

| parameter | kind | default | options |
|---|---|---|---|
| `reference` | file | `examples/genome.fa` | — |
| `min_length` | int | `200` | — |
| `gc_fraction` | number | `0.5` | — |
| `obs_exp` | number | `0.6` | — |
| `window` | int | `100` | — |

Galaxy Tool Shed: `cpgislandfinder / newcpgreport (EMBOSS)`

### epi_dmr_scan

**Sliding-window DMR detection** — Call contiguous stretches of cytosines that cross a methylation threshold.

| parameter | kind | default | options |
|---|---|---|---|
| `methyl` | file | `examples/methylation.cxcg` | — |
| `window` | int | `50` | — |
| `step` | int | `25` | — |
| `min_diff` | number | `0.2` | — |
| `min_fraction` | number | `0.5` | — |
| `bed12` | bool | `False` | — |

Galaxy Tool Shed: `DSS / methylKit dmr`

### epi_methylation_context_summary

**Methylation by cytosine context** — Mean methylation and coverage of cytosines in each sequence context.

| parameter | kind | default | options |
|---|---|---|---|
| `methyl` | file | `examples/methylation.cxcg` | — |
| `context` | choice | `all` | all, CpG, CHG, CHH |
| `min_coverage` | number | `0.0` | — |

Galaxy Tool Shed: `bismark2bedGraph / MethylDackel`

### epi_methylation_per_region

**Methylation over annotated regions** — Aggregate cytosine-level methylation into per-region betas.

| parameter | kind | default | options |
|---|---|---|---|
| `methyl` | file | `examples/methylation.cxcg` | — |
| `regions` | file | `examples/regions.bed` | — |
| `min_fraction` | number | `0.0` | — |
| `head` | int | `200` | — |

Galaxy Tool Shed: `methylprep / DSS merge`

## spatial

*Spatial* — in Galaxy group *Domain Tools*.

| Tool | What it does | Upstream Galaxy tool |
|---|---|---|
| [`spatial_expression_overlay`](#spatial_expression_overlay) | Scatter the tissue coordinates with colour encoding marker expression. | `scanpy plot_spatial` |
| [`spatial_spot_clusters`](#spatial_spot_clusters) | Smooth spot expression over its spatial neighbours then cluster into domains. | `Scanpy spatial / Seurat SPOTlight` |

### spatial_expression_overlay

**Spot plot coloured by a marker** — Scatter the tissue coordinates with colour encoding marker expression.

| parameter | kind | default | options |
|---|---|---|---|
| `coordinates` | file | `examples/counts.tsv` | — |
| `marker` | text | `` | — |
| `x_column` | text | `` | — |
| `y_column` | text | `` | — |

Galaxy Tool Shed: `scanpy plot_spatial`

### spatial_spot_clusters

**Cluster spatial spots into domains** — Smooth spot expression over its spatial neighbours then cluster into domains.

| parameter | kind | default | options |
|---|---|---|---|
| `coordinates` | file | `examples/counts.tsv` | — |
| `k` | int | `3` | — |
| `neighbours` | int | `4` | — |
| `head` | int | `500` | — |

Galaxy Tool Shed: `Scanpy spatial / Seurat SPOTlight`

## compute_indicators_for_satellite_remote_sensing

*Indicators for Satellite Remote Sensing* — in Galaxy group *Domain Tools*.

| Tool | What it does | Upstream Galaxy tool |
|---|---|---|
| [`rs_image_change_detection`](#rs_image_change_detection) | Pixel-wise difference between two dates with area statistics. | `r.change.detect / CCDC` |
| [`rs_indices_time_series`](#rs_indices_time_series) | Per-scene index values with the seasonal trend and a break point. | `phenofit / TIMESAT` |
| [`rs_terrain_from_raster`](#rs_terrain_from_raster) | Slope, aspect and shading statistics for a DEM-like grid. | `r.slope.aspect / SAGA terrain analysis` |
| [`rs_vegetation_indices`](#rs_vegetation_indices) | Compute a normalised band-difference index and its per-class area. | `Sentinel-2 band math / QGIS raster calculator` |

### rs_image_change_detection

**Change detection between two images** — Pixel-wise difference between two dates with area statistics.

| parameter | kind | default | options |
|---|---|---|---|
| `before` | file | `examples/cells.pgm` | — |
| `after` | file | `examples/cells.pgm` | — |
| `threshold` | number | `10.0` | — |
| `mode` | choice | `difference` | difference, ratio, ndbi |
| `figure` | bool | `True` | — |

Galaxy Tool Shed: `r.change.detect / CCDC`

### rs_indices_time_series

**Index trend over a stack of scenes** — Per-scene index values with the seasonal trend and a break point.

| parameter | kind | default | options |
|---|---|---|---|
| `scenes` | file | `examples/counts.tsv` | — |
| `value_column` | text | `` | — |
| `timestep` | int | `16` | — |
| `trend` | choice | `both` | linear, sen_like, both |
| `break_fraction` | number | `0.2` | — |

Galaxy Tool Shed: `phenofit / TIMESAT`

### rs_terrain_from_raster

**Slope and aspect from an elevation grid** — Slope, aspect and shading statistics for a DEM-like grid.

| parameter | kind | default | options |
|---|---|---|---|
| `raster` | file | `examples/cells.pgm` | — |
| `cell_size` | number | `30.0` | — |
| `hillshade` | bool | `True` | — |
| `azimuth` | number | `315.0` | — |
| `altitude` | number | `45.0` | — |

Galaxy Tool Shed: `r.slope.aspect / SAGA terrain analysis`

### rs_vegetation_indices

**Vegetation indices from image bands** — Compute a normalised band-difference index and its per-class area.

| parameter | kind | default | options |
|---|---|---|---|
| `image_file` | file | `examples/cells.pgm` | — |
| `index` | choice | `NDVI` | NDVI, EVI, SAVI, NDWI, GNDVI |
| `soil` | number | `1.0` | — |
| `classify` | bool | `True` | — |
| `cut_low` | number | `0.2` | — |
| `cut_high` | number | `0.6` | — |

Galaxy Tool Shed: `Sentinel-2 band math / QGIS raster calculator`

