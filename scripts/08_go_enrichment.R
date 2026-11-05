library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)

dir.create("results/go", showWarnings = FALSE, recursive = TRUE)
dir.create("results/figures", showWarnings = FALSE, recursive = TRUE)

# 1. Load DESeq2 results
d <- read.csv("results/deseq2_results.csv")
d <- d[!is.na(d$padj), ]

# 2. Gene sets: all tested genes are the background (universe)
universe  <- d$gene_id
sig_up    <- d$gene_id[d$padj < 0.05 & d$log2FoldChange > 1]
sig_down  <- d$gene_id[d$padj < 0.05 & d$log2FoldChange < -1]
cat("Universe:", length(universe), " Up:", length(sig_up),
    " Down:", length(sig_down), "\n")

# 3. GO Biological Process enrichment for one gene set
run_go <- function(genes, label) {
  ego <- enrichGO(gene = genes, universe = universe,
                  OrgDb = org.Hs.eg.db, keyType = "ENSEMBL",
                  ont = "BP", pAdjustMethod = "BH",
                  pvalueCutoff = 0.05, qvalueCutoff = 0.05,
                  readable = TRUE)
  if (is.null(ego) || nrow(as.data.frame(ego)) == 0) {
    message("No enriched GO terms for: ", label)
    return(invisible(NULL))
  }
  write.csv(as.data.frame(ego),
            paste0("results/go/go_bp_", label, ".csv"), row.names = FALSE)
  p <- dotplot(ego, showCategory = 15) +
    ggtitle(paste("GO Biological Process:", label))
  ggsave(paste0("results/figures/go_dotplot_", label, ".png"),
         p, width = 8, height = 7, dpi = 150)
  cat("\nTop terms for", label, ":\n")
  print(head(as.data.frame(ego)[, c("Description", "Count", "p.adjust")], 8))
}

run_go(sig_up, "upregulated")
run_go(sig_down, "downregulated")
