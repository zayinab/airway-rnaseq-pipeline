library(tximport)
library(DESeq2)
library(readr)
library(ggplot2)
library(pheatmap)
library(org.Hs.eg.db)

dir.create("results/figures", showWarnings = FALSE, recursive = TRUE)

meta <- read_csv("metadata/metadata.csv", show_col_types = FALSE)
meta$condition <- factor(meta$condition, levels = c("untreated", "treated"))
meta$cell_line <- factor(meta$cell_line)
files <- file.path("results/salmon", meta$sample, "quant.sf")
names(files) <- meta$sample
stopifnot(all(file.exists(files)))

tx2gene <- read_csv("data/tx2gene.csv", show_col_types = FALSE)
txi <- tximport(files, type = "salmon", tx2gene = tx2gene,
                ignoreTxVersion = TRUE)

dds <- DESeqDataSetFromTximport(txi, colData = as.data.frame(meta),
                                design = ~ cell_line + condition)
dds <- dds[rowSums(counts(dds)) >= 10, ]
dds <- DESeq(dds)
res <- results(dds, contrast = c("condition", "treated", "untreated"))

df <- as.data.frame(res)
df$gene_id <- rownames(df)
df$symbol <- mapIds(org.Hs.eg.db, keys = df$gene_id, keytype = "ENSEMBL",
                    column = "SYMBOL", multiVals = "first")
df <- df[order(df$padj), ]
write.csv(df, "results/deseq2_results.csv", row.names = FALSE)

cat("Genes tested:", nrow(df), "\n")
cat("Significant (padj < 0.05, |log2FC| > 1):",
    sum(df$padj < 0.05 & abs(df$log2FoldChange) > 1, na.rm = TRUE), "\n")

vsd <- vst(dds, blind = FALSE)
p <- plotPCA(vsd, intgroup = c("condition", "cell_line"))
ggsave("results/figures/pca.png", p, width = 6, height = 5, dpi = 150)

df$status <- ifelse(df$padj < 0.05 & df$log2FoldChange > 1, "Up",
             ifelse(df$padj < 0.05 & df$log2FoldChange < -1, "Down", "Not sig"))
v <- ggplot(na.omit(df[, c("log2FoldChange", "pvalue", "status")]),
            aes(log2FoldChange, -log10(pvalue), color = status)) +
  geom_point(alpha = 0.5, size = 1) +
  scale_color_manual(values = c(Up = "firebrick", Down = "steelblue",
                                "Not sig" = "grey70")) +
  theme_minimal() +
  labs(title = "Dexamethasone vs untreated", x = "log2 fold change")
ggsave("results/figures/volcano.png", v, width = 6, height = 5, dpi = 150)

top <- head(df$gene_id[!is.na(df$padj)], 30)
mat <- assay(vsd)[top, ]
rownames(mat) <- ifelse(is.na(df$symbol[match(top, df$gene_id)]), top,
                        df$symbol[match(top, df$gene_id)])
annot <- data.frame(condition = meta$condition, row.names = meta$sample)
pheatmap(mat, scale = "row", annotation_col = annot,
         filename = "results/figures/heatmap_top30.png",
         width = 7, height = 8)

known <- c("DUSP1", "KLF15", "PER1", "TSC22D3")
print(df[df$symbol %in% known, c("symbol", "log2FoldChange", "padj")])
