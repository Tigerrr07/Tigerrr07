#!/usr/bin/env Rscript
# Load and summarize KO vs WT DESeq2 peak results.

# Use a toy/example project root to avoid exposing real filesystem paths.
proj <- "toy_data/FHL2_ATACseq"
out_dir <- file.path(proj, "downstream/deseq2_consensus")

merged_file <- file.path(out_dir, "KO_vs_WT_consensus_peaks_DESeq2.annotated.csv")
sig_file <- file.path(out_dir, "KO_vs_WT_consensus_peaks_DESeq2.sig.csv")

merged <- read.csv(merged_file, check.names = FALSE)
sig <- read.csv(sig_file, check.names = FALSE)
sig_up <- sig[sig$log2FoldChange > 0, , drop = FALSE]
sig_down <- sig[sig$log2FoldChange < 0, , drop = FALSE]

cat("Merged peaks:", nrow(merged), "\n")
cat("Significant peaks:", nrow(sig), "\n")
cat("KO-up peaks (log2FC > 0):", nrow(sig_up), "\n")
cat("KO-down peaks (log2FC < 0):", nrow(sig_down), "\n")
cat("Median |log2FC| in significant peaks:", round(median(abs(sig$log2FoldChange), na.rm = TRUE), 3), "\n")

cat("\nTop annotation categories in significant peaks:\n")
print(head(sort(table(sig$Annotation), decreasing = TRUE), 10))

cat("\nTop genes in significant peaks (nearest gene):\n")
print(head(sort(table(sig$`Gene Name`), decreasing = TRUE), 10))

write.csv(sig_up, file.path(out_dir, "KO_vs_WT_consensus_peaks_DESeq2.sig_up.csv"), row.names = FALSE, quote = FALSE)
write.csv(sig_down, file.path(out_dir, "KO_vs_WT_consensus_peaks_DESeq2.sig_down.csv"), row.names = FALSE, quote = FALSE)

genes_up <- sort(unique(sig_up$`Gene Name`[!is.na(sig_up$`Gene Name`) & sig_up$`Gene Name` != ""]))
genes_down <- sort(unique(sig_down$`Gene Name`[!is.na(sig_down$`Gene Name`) & sig_down$`Gene Name` != ""]))

writeLines(genes_up, file.path(out_dir, "gene_list_up.txt"))
writeLines(genes_down, file.path(out_dir, "gene_list_down.txt"))

cat("\nWrote files:\n")
cat("- KO_vs_WT_consensus_peaks_DESeq2.sig_up.csv\n")
cat("- KO_vs_WT_consensus_peaks_DESeq2.sig_down.csv\n")
cat("- gene_list_up.txt\n")
cat("- gene_list_down.txt\n")
