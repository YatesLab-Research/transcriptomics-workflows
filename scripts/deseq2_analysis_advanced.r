# -------- Parse Command-Line Arguments --------
args <- commandArgs(trailingOnly = TRUE)

# Default values
input_counts <- "data/gene_counts.txt"
input_metadata <- "data/metadata.tsv"
output_dir <- "results/deseq2"
groupA <- "treatment"
groupB <- "control"

# Parse flags
for (i in seq(1, length(args), by = 2)) {
  if (args[i] == "--counts")    input_counts   <- args[i + 1]
  if (args[i] == "--metadata")  input_metadata <- args[i + 1]
  if (args[i] == "--outdir")    output_dir     <- args[i + 1]
  if (args[i] == "--groupA")    groupA         <- args[i + 1]
  if (args[i] == "--groupB")    groupB         <- args[i + 1]
}

plot_dir <- file.path(output_dir, "figures")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

# -------- Load Libraries --------
suppressPackageStartupMessages({
  library(DESeq2)
  library(pheatmap)
  library(ggplot2)
  library(EnhancedVolcano)
  library(RColorBrewer)
  library(biomaRt)
})

# -------- Load Data --------
counts <- read.delim(input_counts, comment.char = "#", row.names = 1)
counts <- counts[, 6:ncol(counts)]
colnames(counts) <- gsub(".sorted.bam", "", colnames(counts))

metadata <- read.table(input_metadata, header = TRUE, sep = "\t", row.names = 1)
metadata <- metadata[colnames(counts), ]
metadata$condition <- factor(metadata$condition)
metadata$batch <- factor(metadata$batch)

# -------- Validate Group Names --------
available_groups <- levels(metadata$condition)

if (!(groupA %in% available_groups) || !(groupB %in% available_groups)) {
  cat("❌ Invalid group name(s) provided.\n")
  cat("Available groups in metadata$condition:\n")
  cat(paste(" -", available_groups), sep = "\n")

  # Suggest valid pairwise combinations
  cat("\n💡 Valid contrast pairs (groupA vs groupB):\n")
  valid_pairs <- combn(available_groups, 2, simplify = FALSE)
  for (pair in valid_pairs) {
    cat(" -", pair[2], "vs", pair[1], "\n")
  }

  stop("Please provide valid groupA and groupB from the list above.")
}

message("✅ Group contrast set: ", groupA, " vs ", groupB)

# -------- DESeq2 with Batch Correction --------
dds <- DESeqDataSetFromMatrix(countData = counts, colData = metadata, design = ~ batch + condition)
dds <- DESeq(dds)

# -------- Contrast Analysis --------
contrast <- c("condition", groupA, groupB)
res <- results(dds, contrast = contrast, alpha = 0.05)
res <- lfcShrink(dds, contrast = contrast, type = "apeglm")

# Save results
write.csv(as.data.frame(res), file = file.path(output_dir, paste0("deseq2_results_", groupA, "_vs_", groupB, ".csv")))

# -------- DE Summary Statistics --------
res_df <- as.data.frame(res)
res_df$significance <- "not_significant"
res_df$significance[res_df$padj < 0.05 & res_df$log2FoldChange > 1] <- "upregulated"
res_df$significance[res_df$padj < 0.05 & res_df$log2FoldChange < -1] <- "downregulated"

summary_table <- table(res_df$significance)
write.table(summary_table, file = file.path(output_dir, "de_gene_summary.txt"), quote = FALSE, sep = "\t", col.names = NA)

# -------- Gene Annotation --------
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
annot <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol", "description"),
               filters = "ensembl_gene_id",
               values = rownames(res),
               mart = ensembl)

res_annot <- merge(res_df, annot, by.x = "row.names", by.y = "ensembl_gene_id", all.x = TRUE)
colnames(res_annot)[1] <- "Ensembl_ID"

write.csv(res_annot, file = file.path(output_dir, paste0("deseq2_results_annotated_", groupA, "_vs_", groupB, ".csv")), row.names = FALSE)

# -------- Visualization --------
cb_palette <- brewer.pal(8, "Set2")

save_plot <- function(name, expr, width = 7, height = 6, res = 300) {
  png(file.path(plot_dir, paste0(name, ".png")), width = width, height = height, units = "in", res = res)
  expr
  dev.off()

  pdf(file.path(plot_dir, paste0(name, ".pdf")), width = width, height = height)
  expr
  dev.off()
}

# MA Plot
save_plot("ma_plot", {
  plotMA(res, main = paste("MA Plot:", groupA, "vs", groupB), ylim = c(-5, 5), alpha = 0.05, cex = 0.8)
})

# Volcano Plot
save_plot("volcano_plot", {
  EnhancedVolcano(res,
    lab = rownames(res),
    x = 'log2FoldChange',
    y = 'pvalue',
    title = paste("Volcano Plot:", groupA, "vs", groupB),
    pCutoff = 0.05,
    FCcutoff = 1.5,
    pointSize = 2.5,
    labSize = 3.5,
    colAlpha = 0.75,
    col = c("grey40", cb_palette[2], cb_palette[1], cb_palette[3]),
    drawConnectors = TRUE
  )
})

# Heatmap
vsd <- vst(dds, blind = FALSE)
top_genes <- head(order(res$padj), 30)
save_plot("heatmap_top30", {
  pheatmap(assay(vsd)[top_genes, ],
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    show_rownames = TRUE,
    show_colnames = TRUE,
    annotation_col = metadata,
    fontsize_row = 8,
    fontsize_col = 10,
    color = colorRampPalette(rev(brewer.pal(n = 7, "RdYlBu")))(100),
    main = paste("Top 30 DE Genes:", groupA, "vs", groupB)
  )
})

# PCA
save_plot("pca_plot", {
  pca_data <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
  percentVar <- round(100 * attr(pca_data, "percentVar"))
  ggplot(pca_data, aes(PC1, PC2, color = condition)) +
    geom_point(size = 3) +
    scale_color_manual(values = cb_palette) +
    labs(
      title = "PCA of Samples",
      x = paste0("PC1 (", percentVar[1], "%)"),
      y = paste0("PC2 (", percentVar[2], "%)")
    ) +
    theme_minimal(base_size = 14)
})
