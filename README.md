# PC9-PC9OR Transcriptome Analysis

This repository contains an R pipeline for PC9 vs PC9-OR transcriptome analysis. It performs:

- raw-count import and sample-group detection;
- DESeq2 differential expression analysis;
- PCA, MA plot, volcano plot, and DEG heatmap generation;
- GO and KEGG pathway enrichment with `clusterProfiler`;
- CSV export of all genes, significant DEGs, up-regulated genes, down-regulated genes, normalized counts, and enrichment results.

## Quick start

Open R/RStudio or a terminal on Windows and run:

```r
Rscript scripts/rnaseq_pc9_pc9or_analysis.R \
  --input_dir "C:/Users/28907/Desktop/GEOPC9-PC9OR" \
  --organism human \
  --gene_id_type SYMBOL
```

If your dataset has a specific count matrix and metadata file, use:

```r
Rscript scripts/rnaseq_pc9_pc9or_analysis.R \
  --count_matrix "C:/Users/28907/Desktop/GEOPC9-PC9OR/counts.csv" \
  --metadata "C:/Users/28907/Desktop/GEOPC9-PC9OR/metadata.csv" \
  --output_dir "C:/Users/28907/Desktop/GEOPC9-PC9OR/results" \
  --organism human \
  --gene_id_type SYMBOL
```

## Input format

The count matrix should be a CSV, TSV, or TXT file with genes in rows and samples in columns. The first column should contain gene symbols, Ensembl IDs, or Entrez IDs.

Example:

| gene | PC9_1 | PC9_2 | PC9OR_1 | PC9OR_2 |
| --- | ---: | ---: | ---: | ---: |
| EGFR | 120 | 115 | 430 | 410 |
| ACTB | 900 | 870 | 860 | 890 |

The optional metadata file should contain `sample` and `group` columns:

| sample | group |
| --- | --- |
| PC9_1 | PC9 |
| PC9_2 | PC9 |
| PC9OR_1 | PC9OR |
| PC9OR_2 | PC9OR |

## Important options

- `--gene_id_type`: `SYMBOL`, `ENSEMBL`, or `ENTREZID`.
- `--organism`: `human` or `mouse`.
- `--padj_cutoff`: adjusted P-value cutoff; default is `0.05`.
- `--log2fc_cutoff`: absolute log2 fold-change cutoff; default is `1`.
- `--control_group`: default is `PC9`.
- `--treatment_group`: default is `PC9OR`.

## Output

By default, results are written to `C:/Users/28907/Desktop/GEOPC9-PC9OR/results`:

- `tables/DESeq2_all_genes_PC9OR_vs_PC9.csv`
- `tables/DEG_significant_PC9OR_vs_PC9.csv`
- `tables/DEG_up_in_PC9OR.csv`
- `tables/DEG_down_in_PC9OR.csv`
- `tables/DESeq2_normalized_counts.csv`
- `tables/*_GO_BP.csv`, `*_GO_CC.csv`, `*_GO_MF.csv`, `*_KEGG.csv`
- `plots/PCA_PC9_PC9OR.*`
- `plots/Volcano_PC9OR_vs_PC9.*`
- `plots/MA_plot_PC9OR_vs_PC9.pdf`
- `plots/Top_DEG_heatmap.pdf`
