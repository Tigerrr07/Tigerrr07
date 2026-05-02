#!/usr/bin/env Rscript
# Consensus peaks: KO vs WT differential accessibility (DESeq2).
# Requires: BiocManager::install("DESeq2")

# module load gcc/12.3.0
# module load R/4.4.0
library(DESeq2)
library(ggplot2)
library(ggrepel)

# Use a toy/example project root to avoid exposing real filesystem paths.
proj <- "toy_data/FHL2_ATACseq"
fc_file <- file.path(
  proj,
  "results/bwa/merged_replicate/macs2/broad_peak/consensus",
  "consensus_peaks.mRp.clN.featureCounts.txt"
)
anno_file <- file.path(
  proj,
  "results/bwa/merged_replicate/macs2/broad_peak/consensus",
  "consensus_peaks.mRp.clN.annotatePeaks.txt"
)
out_dir <- file.path(proj, "downstream/deseq2_consensus")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

fc <- read.delim(fc_file, comment.char = "#", check.names = FALSE)

peak_id <- fc[["Geneid"]]
count_cols <- grep("\\.sorted\\.bam$", colnames(fc), value = TRUE)
if (length(count_cols) != 4L) {
  count_cols <- colnames(fc)[(ncol(fc) - 3):ncol(fc)]
}

counts <- as.matrix(fc[, count_cols, drop = FALSE])
storage.mode(counts) <- "integer"
colnames(counts) <- sub("\\.mLb\\.clN\\.sorted\\.bam$", "", colnames(counts))
rownames(counts) <- peak_id

samples <- colnames(counts)
condition <- ifelse(grepl("^WT_", samples), "WT", "KO")
coldata <- data.frame(
  condition = factor(condition, levels = c("WT", "KO")),
  row.names = samples
)
coldata <- coldata[colnames(counts), , drop = FALSE]

dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = coldata,
  design = ~ condition
)
dds <- dds[rowSums(counts(dds)) >= 10, ]
dds <- DESeq(dds)

res <- results(dds, contrast = c("condition", "KO", "WT"), alpha = 0.05)
res <- res[order(res$padj, na.last = TRUE), ]

res_df <- as.data.frame(res)
res_df$PeakID <- rownames(res_df)

anno <- read.delim(anno_file, check.names = FALSE)
names(anno)[1] <- "PeakID"

merged <- merge(res_df, anno, by = "PeakID", all.x = TRUE, sort = FALSE)

out_csv <- file.path(out_dir, "KO_vs_WT_consensus_peaks_DESeq2.annotated.csv")
write.csv(merged, out_csv, row.names = FALSE, quote = FALSE)

sig <- merged[!is.na(merged$padj) & merged$padj < 0.05 & abs(merged$log2FoldChange) >= 1, ]
write.csv(sig, file.path(out_dir, "KO_vs_WT_consensus_peaks_DESeq2.sig.csv"), row.names = FALSE, quote = FALSE)

vsd <- vst(dds, blind = FALSE)
pca_data <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
percent_var <- round(100 * attr(pca_data, "percentVar"))
pca_pdf <- file.path(out_dir, "KO_vs_WT_consensus_peaks_PCA.pdf")

pca_data$sample <- rownames(pca_data)

p <- ggplot(pca_data, aes(x = PC1, y = PC2, color = condition, label = sample)) +
  geom_point(size = 3) +
  geom_text_repel(
    size = 3,
    show.legend = FALSE,
    box.padding = 0.3,
    max.overlaps = Inf,
    segment.size = 0.2
  ) +
  scale_color_manual(values = c("WT" = "#d62728", "KO" = "#1f77b4")) +
  labs(
    x = paste0("PC1: ", percent_var[1], "% variance"),
    y = paste0("PC2: ", percent_var[2], "% variance")
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.title = element_blank(),
    panel.grid = element_blank()
  )

ggsave(filename = pca_pdf, plot = p, width = 3.5, height = 2.5)
saveRDS(dds, file.path(out_dir, "dds_consensus_peaks.rds"))
# qs::qsave(dds, file.path(out_dir, "dds_consensus_peaks.qs"))
message("Wrote: ", out_csv)
message("Significant peaks (padj<0.05, |log2FC|>=1): ", nrow(sig))
message("PCA plot: ", pca_pdf)
