#!/usr/bin/env Rscript

# PC9 vs PC9-OR transcriptome analysis pipeline
# Author: generated for the PC9-PC9OR dataset
# Purpose: differential expression analysis and pathway enrichment analysis
#
# Usage examples:
#   Rscript scripts/rnaseq_pc9_pc9or_analysis.R \
#     --input_dir "C:/Users/28907/Desktop/GEOPC9-PC9OR" \
#     --organism human \
#     --gene_id_type SYMBOL
#
#   Rscript scripts/rnaseq_pc9_pc9or_analysis.R \
#     --count_matrix "C:/Users/28907/Desktop/GEOPC9-PC9OR/counts.csv" \
#     --metadata "C:/Users/28907/Desktop/GEOPC9-PC9OR/metadata.csv" \
#     --output_dir "C:/Users/28907/Desktop/GEOPC9-PC9OR/results"
#
# Expected input:
#   1. A raw integer count matrix in CSV/TSV/TXT format. Rows are genes and columns
#      are samples. The first column should contain gene IDs or gene symbols.
#   2. Optional metadata CSV/TSV/TXT with at least two columns: sample and group.
#      group should contain PC9 and PC9OR/PC9-OR values. If metadata is not given,
#      sample groups are inferred from column names containing PC9OR/PC9-OR/OR and PC9.
#
# Output:
#   - DESeq2 differential expression tables
#   - PCA, volcano, MA, heatmap plots
#   - GO BP/CC/MF and KEGG enrichment tables and dotplots
#   - Up/down-regulated gene lists

options(stringsAsFactors = FALSE)

required_cran <- c("optparse", "ggplot2", "ggrepel", "pheatmap", "readr", "dplyr", "tibble", "ashr")
required_bioc <- c("DESeq2", "clusterProfiler", "org.Hs.eg.db", "org.Mm.eg.db", "AnnotationDbi", "enrichplot")

install_if_missing <- function(pkg, bioc = FALSE) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    message("Installing missing package: ", pkg)
    if (bioc) {
      if (!requireNamespace("BiocManager", quietly = TRUE)) {
        install.packages("BiocManager", repos = "https://cloud.r-project.org")
      }
      BiocManager::install(pkg, ask = FALSE, update = FALSE)
    } else {
      install.packages(pkg, repos = "https://cloud.r-project.org")
    }
  }
}

for (pkg in required_cran) install_if_missing(pkg, bioc = FALSE)
for (pkg in required_bioc) install_if_missing(pkg, bioc = TRUE)

suppressPackageStartupMessages({
  library(optparse)
  library(DESeq2)
  library(clusterProfiler)
  library(AnnotationDbi)
  library(enrichplot)
  library(ggplot2)
  library(ggrepel)
  library(pheatmap)
  library(readr)
  library(dplyr)
  library(tibble)
})

option_list <- list(
  make_option("--input_dir", type = "character", default = "C:/Users/28907/Desktop/GEOPC9-PC9OR",
              help = "Directory containing transcriptome count files [default: %default]"),
  make_option("--count_matrix", type = "character", default = NA,
              help = "Path to count matrix CSV/TSV/TXT. If omitted, detected automatically from input_dir."),
  make_option("--metadata", type = "character", default = NA,
              help = "Path to sample metadata with sample and group columns."),
  make_option("--output_dir", type = "character", default = NA,
              help = "Output directory [default: input_dir/results]"),
  make_option("--organism", type = "character", default = "human",
              help = "Organism: human or mouse [default: %default]"),
  make_option("--gene_id_type", type = "character", default = "SYMBOL",
              help = "Gene ID type in count matrix: SYMBOL, ENSEMBL, or ENTREZID [default: %default]"),
  make_option("--control_group", type = "character", default = "PC9",
              help = "Control group name [default: %default]"),
  make_option("--treatment_group", type = "character", default = "PC9OR",
              help = "Treatment/resistant group name [default: %default]"),
  make_option("--padj_cutoff", type = "double", default = 0.05,
              help = "Adjusted P-value cutoff [default: %default]"),
  make_option("--log2fc_cutoff", type = "double", default = 1,
              help = "Absolute log2 fold-change cutoff [default: %default]")
)
opt <- parse_args(OptionParser(option_list = option_list))

if (is.na(opt$output_dir)) {
  opt$output_dir <- file.path(opt$input_dir, "results")
}
dir.create(opt$output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(opt$output_dir, "plots"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(opt$output_dir, "tables"), recursive = TRUE, showWarnings = FALSE)

message("Input directory: ", opt$input_dir)
message("Output directory: ", opt$output_dir)

read_table_auto <- function(file) {
  ext <- tolower(tools::file_ext(file))
  if (ext == "csv") {
    readr::read_csv(file, show_col_types = FALSE)
  } else if (ext %in% c("tsv", "txt")) {
    readr::read_tsv(file, show_col_types = FALSE)
  } else {
    stop("Unsupported file extension: ", ext, ". Please use CSV, TSV, or TXT.")
  }
}

detect_count_matrix <- function(input_dir) {
  candidates <- list.files(input_dir, pattern = "\\.(csv|tsv|txt)$", full.names = TRUE, ignore.case = TRUE)
  if (length(candidates) == 0) stop("No CSV/TSV/TXT files found in ", input_dir)

  score_file <- function(file) {
    name <- tolower(basename(file))
    score <- 0
    if (grepl("count|counts|raw|matrix|expression|expr", name)) score <- score + 10
    if (grepl("metadata|meta|sample|clinical|phenotype", name)) score <- score - 20
    score
  }
  candidates[order(vapply(candidates, score_file, numeric(1)), decreasing = TRUE)][1]
}

clean_gene_ids <- function(ids) {
  ids <- as.character(ids)
  # Remove common Ensembl version suffixes, e.g. ENSG000001.13 -> ENSG000001
  sub("\\.[0-9]+$", "", ids)
}

if (is.na(opt$count_matrix)) {
  opt$count_matrix <- detect_count_matrix(opt$input_dir)
  message("Detected count matrix: ", opt$count_matrix)
}

count_df <- read_table_auto(opt$count_matrix)
if (ncol(count_df) < 3) stop("Count matrix must contain one gene column and at least two sample columns.")

gene_col <- names(count_df)[1]
count_mat <- as.data.frame(count_df)
rownames(count_mat) <- clean_gene_ids(count_mat[[gene_col]])
count_mat[[gene_col]] <- NULL
count_mat <- as.matrix(count_mat)
storage.mode(count_mat) <- "numeric"
count_mat[is.na(count_mat)] <- 0
count_mat <- round(count_mat)

# Remove duplicated gene IDs by keeping the row with the largest total count.
row_totals <- rowSums(count_mat)
count_mat <- count_mat[order(row_totals, decreasing = TRUE), , drop = FALSE]
count_mat <- count_mat[!duplicated(rownames(count_mat)), , drop = FALSE]

# Remove genes with almost no reads to improve statistical power.
keep <- rowSums(count_mat >= 10) >= 2
count_mat <- count_mat[keep, , drop = FALSE]
message("Genes retained after low-count filtering: ", nrow(count_mat))

read_metadata <- function(metadata_file, samples, control_group, treatment_group) {
  meta <- read_table_auto(metadata_file) |> as.data.frame()
  names(meta) <- tolower(names(meta))
  if (!all(c("sample", "group") %in% names(meta))) {
    stop("Metadata must contain columns named sample and group.")
  }
  meta <- meta[, c("sample", "group")]
  meta$sample <- as.character(meta$sample)
  meta$group <- as.character(meta$group)
  rownames(meta) <- meta$sample
  meta <- meta[samples, , drop = FALSE]
  if (any(is.na(meta$group))) stop("Metadata is missing group labels for: ", paste(samples[is.na(meta$group)], collapse = ", "))
  normalize_groups(meta, control_group, treatment_group)
}

normalize_groups <- function(meta, control_group, treatment_group) {
  raw <- toupper(gsub("[-_ ]", "", meta$group))
  control_key <- toupper(gsub("[-_ ]", "", control_group))
  treatment_key <- toupper(gsub("[-_ ]", "", treatment_group))
  raw[raw %in% c("PC9OR", "PC9ORG", "PC9RESISTANT", "OR", "RESISTANT")] <- treatment_key
  raw[raw %in% c("PC9", "PARENTAL", "SENSITIVE")] <- control_key
  if (!all(raw %in% c(control_key, treatment_key))) {
    stop("Only two groups are supported. Found unrecognized group labels: ", paste(unique(meta$group[!raw %in% c(control_key, treatment_key)]), collapse = ", "))
  }
  meta$group <- ifelse(raw == treatment_key, treatment_group, control_group)
  meta$group <- factor(meta$group, levels = c(control_group, treatment_group))
  meta
}

infer_metadata <- function(samples, control_group, treatment_group) {
  normalized <- toupper(gsub("[-_ ]", "", samples))
  group <- rep(NA_character_, length(samples))
  group[grepl("PC9OR|OR|RESIST", normalized)] <- treatment_group
  group[is.na(group) & grepl("PC9", normalized)] <- control_group
  if (any(is.na(group))) {
    stop("Could not infer group from sample names: ", paste(samples[is.na(group)], collapse = ", "),
         ". Please provide --metadata with sample and group columns.")
  }
  meta <- data.frame(sample = samples, group = group, row.names = samples)
  meta$group <- factor(meta$group, levels = c(control_group, treatment_group))
  meta
}

sample_names <- colnames(count_mat)
if (!is.na(opt$metadata)) {
  coldata <- read_metadata(opt$metadata, sample_names, opt$control_group, opt$treatment_group)
} else {
  coldata <- infer_metadata(sample_names, opt$control_group, opt$treatment_group)
}

write.csv(coldata, file.path(opt$output_dir, "tables", "sample_metadata_used.csv"), row.names = TRUE)
message("Group summary:")
print(table(coldata$group))

if (any(table(coldata$group) < 2)) {
  warning("One group has fewer than two samples. DESeq2 can run, but biological replicates are strongly recommended.")
}

# Differential expression with DESeq2. The contrast reports treatment vs control,
# so positive log2FoldChange means higher in PC9OR than PC9.
dds <- DESeqDataSetFromMatrix(countData = count_mat, colData = coldata, design = ~ group)
dds <- DESeq(dds)
res <- results(dds, contrast = c("group", opt$treatment_group, opt$control_group), alpha = opt$padj_cutoff)
res_lfc <- lfcShrink(dds, contrast = c("group", opt$treatment_group, opt$control_group), res = res, type = "ashr")

res_df <- as.data.frame(res_lfc) |>
  rownames_to_column("gene_id") |>
  arrange(padj)
res_df$direction <- "Not_significant"
res_df$direction[!is.na(res_df$padj) & res_df$padj < opt$padj_cutoff & res_df$log2FoldChange >= opt$log2fc_cutoff] <- "Up_in_PC9OR"
res_df$direction[!is.na(res_df$padj) & res_df$padj < opt$padj_cutoff & res_df$log2FoldChange <= -opt$log2fc_cutoff] <- "Down_in_PC9OR"

write.csv(res_df, file.path(opt$output_dir, "tables", "DESeq2_all_genes_PC9OR_vs_PC9.csv"), row.names = FALSE)
write.csv(filter(res_df, direction == "Up_in_PC9OR"), file.path(opt$output_dir, "tables", "DEG_up_in_PC9OR.csv"), row.names = FALSE)
write.csv(filter(res_df, direction == "Down_in_PC9OR"), file.path(opt$output_dir, "tables", "DEG_down_in_PC9OR.csv"), row.names = FALSE)
write.csv(filter(res_df, direction != "Not_significant"), file.path(opt$output_dir, "tables", "DEG_significant_PC9OR_vs_PC9.csv"), row.names = FALSE)

message("Significant DEG counts:")
print(table(res_df$direction))

# Quality-control and DEG plots.
vsd <- vst(dds, blind = FALSE)

pca_data <- plotPCA(vsd, intgroup = "group", returnData = TRUE)
percent_var <- round(100 * attr(pca_data, "percentVar"))
p_pca <- ggplot(pca_data, aes(PC1, PC2, color = group, label = name)) +
  geom_point(size = 4) +
  ggrepel::geom_text_repel(max.overlaps = Inf) +
  xlab(paste0("PC1: ", percent_var[1], "% variance")) +
  ylab(paste0("PC2: ", percent_var[2], "% variance")) +
  theme_bw(base_size = 12)
ggsave(file.path(opt$output_dir, "plots", "PCA_PC9_PC9OR.pdf"), p_pca, width = 7, height = 5)
ggsave(file.path(opt$output_dir, "plots", "PCA_PC9_PC9OR.png"), p_pca, width = 7, height = 5, dpi = 300)

pdf(file.path(opt$output_dir, "plots", "MA_plot_PC9OR_vs_PC9.pdf"), width = 7, height = 5)
plotMA(res_lfc, ylim = c(-6, 6), main = "PC9OR vs PC9")
dev.off()

volcano_df <- res_df |>
  mutate(log10_padj = -log10(padj),
         label = ifelse(direction != "Not_significant" & rank(padj, ties.method = "first") <= 20, gene_id, NA))
volcano_df$log10_padj[is.infinite(volcano_df$log10_padj)] <- max(volcano_df$log10_padj[is.finite(volcano_df$log10_padj)], na.rm = TRUE) + 1
p_volcano <- ggplot(volcano_df, aes(log2FoldChange, log10_padj, color = direction)) +
  geom_point(alpha = 0.7, size = 1.4, na.rm = TRUE) +
  geom_vline(xintercept = c(-opt$log2fc_cutoff, opt$log2fc_cutoff), linetype = "dashed") +
  geom_hline(yintercept = -log10(opt$padj_cutoff), linetype = "dashed") +
  ggrepel::geom_text_repel(aes(label = label), max.overlaps = Inf, na.rm = TRUE) +
  scale_color_manual(values = c("Down_in_PC9OR" = "#377EB8", "Not_significant" = "grey70", "Up_in_PC9OR" = "#E41A1C")) +
  labs(title = "PC9OR vs PC9", x = "log2 fold change", y = "-log10 adjusted P-value") +
  theme_bw(base_size = 12)
ggsave(file.path(opt$output_dir, "plots", "Volcano_PC9OR_vs_PC9.pdf"), p_volcano, width = 7, height = 6)
ggsave(file.path(opt$output_dir, "plots", "Volcano_PC9OR_vs_PC9.png"), p_volcano, width = 7, height = 6, dpi = 300)

sig_genes <- filter(res_df, direction != "Not_significant")$gene_id
if (length(sig_genes) >= 2) {
  top_genes <- head(filter(res_df, direction != "Not_significant")$gene_id, 50)
  heat_mat <- assay(vsd)[top_genes, , drop = FALSE]
  heat_mat <- heat_mat - rowMeans(heat_mat)
  pdf(file.path(opt$output_dir, "plots", "Top_DEG_heatmap.pdf"), width = 8, height = 10)
  pheatmap(heat_mat, annotation_col = data.frame(group = coldata$group, row.names = rownames(coldata)),
           show_rownames = TRUE, fontsize_row = 7, main = "Top DEGs: PC9OR vs PC9")
  dev.off()
}

# Gene ID conversion and enrichment analysis.
organism <- tolower(opt$organism)
if (organism == "human") {
  org_db <- org.Hs.eg.db
  kegg_org <- "hsa"
} else if (organism == "mouse") {
  org_db <- org.Mm.eg.db
  kegg_org <- "mmu"
} else {
  stop("Unsupported organism: ", opt$organism, ". Use human or mouse.")
}

gene_id_type <- toupper(opt$gene_id_type)
if (!gene_id_type %in% c("SYMBOL", "ENSEMBL", "ENTREZID")) {
  stop("--gene_id_type must be SYMBOL, ENSEMBL, or ENTREZID.")
}

map_to_entrez <- function(genes) {
  genes <- unique(clean_gene_ids(genes))
  if (length(genes) == 0) return(character(0))
  if (gene_id_type == "ENTREZID") return(unique(as.character(genes)))
  mapped <- AnnotationDbi::select(org_db, keys = genes, keytype = gene_id_type,
                                  columns = c("ENTREZID", "SYMBOL"))
  unique(na.omit(mapped$ENTREZID))
}

run_enrichment <- function(genes, prefix) {
  entrez <- map_to_entrez(genes)
  if (length(entrez) < 5) {
    warning(prefix, ": fewer than five genes mapped to ENTREZID; skipping enrichment.")
    return(invisible(NULL))
  }

  for (ont in c("BP", "CC", "MF")) {
    ego <- enrichGO(gene = entrez, OrgDb = org_db, keyType = "ENTREZID", ont = ont,
                    pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.2,
                    readable = TRUE)
    ego_df <- as.data.frame(ego)
    write.csv(ego_df, file.path(opt$output_dir, "tables", paste0(prefix, "_GO_", ont, ".csv")), row.names = FALSE)
    if (nrow(ego_df) > 0) {
      p <- dotplot(ego, showCategory = 20) + ggtitle(paste(prefix, "GO", ont))
      ggsave(file.path(opt$output_dir, "plots", paste0(prefix, "_GO_", ont, "_dotplot.pdf")), p, width = 9, height = 7)
    }
  }

  ekegg <- enrichKEGG(gene = entrez, organism = kegg_org, pAdjustMethod = "BH",
                      pvalueCutoff = 0.05, qvalueCutoff = 0.2)
  ekegg <- setReadable(ekegg, OrgDb = org_db, keyType = "ENTREZID")
  ekegg_df <- as.data.frame(ekegg)
  write.csv(ekegg_df, file.path(opt$output_dir, "tables", paste0(prefix, "_KEGG.csv")), row.names = FALSE)
  if (nrow(ekegg_df) > 0) {
    p <- dotplot(ekegg, showCategory = 20) + ggtitle(paste(prefix, "KEGG"))
    ggsave(file.path(opt$output_dir, "plots", paste0(prefix, "_KEGG_dotplot.pdf")), p, width = 9, height = 7)
  }
}

run_enrichment(filter(res_df, direction != "Not_significant")$gene_id, "DEG_all")
run_enrichment(filter(res_df, direction == "Up_in_PC9OR")$gene_id, "DEG_up_in_PC9OR")
run_enrichment(filter(res_df, direction == "Down_in_PC9OR")$gene_id, "DEG_down_in_PC9OR")

# Save normalized counts for downstream plotting or validation.
norm_counts <- counts(dds, normalized = TRUE) |>
  as.data.frame() |>
  rownames_to_column("gene_id")
write.csv(norm_counts, file.path(opt$output_dir, "tables", "DESeq2_normalized_counts.csv"), row.names = FALSE)

message("Analysis complete. Results saved to: ", opt$output_dir)
