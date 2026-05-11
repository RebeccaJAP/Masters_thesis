  # Load packages
library(ggplot2)
library(edgeR)
library(tidyr)
library(dplyr)
library(Seurat)
library(ggrepel)

  # Adjust the cutoffs if needed
padj_cutoff <- 0.05
sample_cutoff <- 5 # How many pools does a gene need to be DE to be of interest?


  # Get the filtered pseudobulk object
load("Pseudobulk_data/pseudobulk_obj_filtered.RData")


  # Add the pericyte category for each sample: is the proportion of pericytes high or low
  # cutoff limit between the categories defined in pericyte_high_low_categorizing.R
pseudobulk_obj@meta.data <- pseudobulk_obj@meta.data %>%
  mutate(per_category = case_when(Day_fixed == 20 & proportion_Per > 0.1 ~ "high",
                                  Day_fixed == 20 & proportion_Per <= 0.1 ~ "low"))


  # Get samples with at least one low-pericyte and one high-pericyte donor
data <- pseudobulk_obj@meta.data %>%
  filter(diff_method != "indv" & Day_fixed == 20) %>%
  group_by(SampleID) %>%
  mutate(n=n_distinct(per_category)) %>%
  ungroup() %>%
  filter(n > 1)

valid_samples <- unique(data$SampleID)


  # Get the relevant metadata and separate by samples
meta_by_sample <- lapply(valid_samples, function(s) {
  pseudobulk_obj@meta.data[pseudobulk_obj$SampleID == s,]
})

names(meta_by_sample) <- valid_samples


  # Get the relevant expression data and separate by samples
expr_by_sample <- lapply(valid_samples, function(s) {
  cols <- rownames(pseudobulk_obj@meta.data[pseudobulk_obj$SampleID == s & pseudobulk_obj$Day_fixed == 20,])
  pseudobulk_obj@assays$RNA$counts[, cols, drop = FALSE]
})

names(expr_by_sample) <- valid_samples


  # Calculate normalization factors
d0 <- lapply(expr_by_sample, function(mat) {
  calcNormFactors(DGEList(mat))
})


  # Create design matrix
mm_all <- c()

for (i in seq(1, length(valid_samples))) {
  per_category <- meta_by_sample[[i]]$per_category
  per_category <- factor(per_category, levels = c("low", "high"))
  mm <- model.matrix(~0 + per_category)
  mm_all[[i]] <- mm
}

names(mm_all) <- valid_samples


  # Run voom, contrast pericyte-high and pericyte-low, get DEGs
toptables <- c()

for (i in seq(1, length(valid_samples))) {
  y <- voom(d0[[i]], mm_all[[i]], plot = T)
  fit <- lmFit(y, mm_all[[i]])
  contr <- makeContrasts(per_categoryhigh - per_categorylow, levels = colnames(coef(fit)))
  tmp <- contrasts.fit(fit, contr)
  tmp <- eBayes(tmp)
  top.table <- topTable(tmp, sort.by = "P", n = Inf)
  toptables[[i]] <- top.table
}

names(toptables) <- valid_samples


  # Check how well the high and low pericyte samples separate
for (i in seq(1, length(valid_samples))) {
  
  mds <- plotMDS(d0[[i]], labels=meta_by_sample[[i]]$donor, plot = F)
  
  mds_df <- data.frame(
    Dim1 = mds$x,
    Dim2 = mds$y,
    Color = meta_by_sample[[i]]$per_category,
    Sample = meta_by_sample[[i]]$SampleID
  )
  
  plot <- ggplot(mds_df, aes(x = Dim1, y = Dim2, color = as.factor(Color))) +
    geom_point(size = 3, pch=19) +
    geom_text_repel(aes(label = meta_by_sample[[i]]$donor), size = 3, max.overlaps = Inf, box.padding = 0.5) +
    theme_minimal() +
    labs(color = "Category", title = paste0("MDS Plot for sample ", valid_samples[i])) +
    xlab(paste0(mds$axislabel[1], " ", mds$dim.plot[1], " (", round(100*mds$var.explained[1]), "%)")) +
    ylab(paste0(mds$axislabel[1], " ", mds$dim.plot[2], " (", round(100*mds$var.explained[2]), "%)"))
  
  print(plot)
  
}


  # Get the names of differentially expressed genes for each sample
signif_DE_genes <- lapply(toptables, function(sample) {
  sample$gene <- rownames(sample)
  sample$gene[which(sample$adj.P.Val < padj_cutoff)]
})


  # Limit the limma results to statistically significant ones
toptables_signif <- lapply(toptables, function(sample) {
  sample$gene <- rownames(sample)
  sample[which(sample$adj.P.Val < padj_cutoff),]
})


  # How many samples each gene is differentially expressed in?
gene_counts_limma <- sort(table(unlist(signif_DE_genes)), decreasing = T)

  # All significant DEGs in one list
all_de_genes <- names(gene_counts_limma)


  # Is gene differentailly expressed per pool? TRUE/FALSE value for each sample
binary_mat <- sapply(signif_DE_genes, function(g) all_de_genes %in% g)
rownames(binary_mat) <- all_de_genes


  # Make into dataframe and and a column for gene name
binary_df <- as.data.frame(binary_mat)
binary_df$gene <- rownames(binary_df)

  # Change into long format
binary_long <- pivot_longer(binary_df, cols = -gene, names_to = "dataset", values_to = "present")


  # Create a matrix of log fold change values
logfc_matrix <- matrix(NA, nrow = length(all_de_genes), ncol = length(toptables_signif),
                       dimnames = list(all_de_genes, names(toptables_signif)))

for (i in seq_along(toptables_signif)) {
  df <- toptables_signif[[i]]
  gene_logfc <- setNames(df$logFC, df$gene)
  logfc_matrix[names(gene_logfc), i] <- gene_logfc
}

  # Make into a dataframe and long dataframe
logfc_long <- as.data.frame(logfc_matrix)
logfc_long$gene <- rownames(logfc_matrix)
logfc_long <- pivot_longer(logfc_long, -gene, names_to = "dataset", values_to = "logFC")

heatmap_df <- left_join(binary_long, logfc_long, by = c("gene", "dataset"))

heatmap_df_cutoff <- heatmap_df %>%
  group_by(gene) %>%
  mutate(n_present = sum(present), dataset = sub("\\_d..", "", dataset)) %>%
  filter(n_present >= sample_cutoff) %>%
  arrange(desc(abs(logFC)))


  # Plot heatmap
  # Note: produces warning about missing points for the non-significant genes;
  # these are missing on purpose and the warning should be ignored

heatmap_d20_DGE <- ggplot((heatmap_df_cutoff %>% filter(n_present >= sample_cutoff)), aes(x = gene, y = dataset)) +
  geom_point(aes(color = logFC, size = abs(logFC))) +
  scale_color_gradientn(colors=c("#1600A3", "slateblue1", "white", "hotpink", "#D10069"), rescaler = ~ scales::rescale_mid(.x, mid = 0)) +
  theme_minimal() +
  labs(x = "Gene", y = "Pool", color = "log2FC", size = "Absolute value of log2FC") +
  ggtitle(paste0("DGE in pericyte-high vs. pericyte-low donors on day 20")) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size=8),
        title = element_text(size = 12))

  # Save the heatmap
pdf(file = "Gene_expression_analyses/DGE_and_GO_enrichment/Pericyte_high_low/Pericyte_high_low_DGE_heatmap.pdf",
    width = 8.27,
    height = 4.5)
plot(heatmap_d20_DGE)
dev.off()


  # Get genes for GO enrichment analysis:

  # Get genes which are significantly differentially expressed in at least 4 pools
  # and which are not significantly upregulated in any pool
downregulated_genes <- heatmap_df %>%
  group_by(gene) %>%
  mutate(n_present = sum(present), max_logFC = max(logFC, na.rm = T)) %>%
  filter(n_present >= 4 & max_logFC < 0) %>%
  distinct(gene)

  # Get genes which are significantly differentially expressed in at least 4 pools
  # and which are not significantly downregulated in any pool
upregulated_genes <- heatmap_df %>%
  group_by(gene) %>%
  mutate(n_present = sum(present), min_logFC = min(logFC, na.rm = T)) %>%
  filter(n_present >= 4 & min_logFC > 0) %>%
  distinct(gene)

  # Save the genes for further analysis
save(downregulated_genes,
     file = "Gene_expression_analyses/DGE_and_GO_enrichment/Pericyte_high_low/limma_4_down_pericyte_high_low.RData")
save(upregulated_genes,
     file = "Gene_expression_analyses/DGE_and_GO_enrichment/Pericyte_high_low/limma_4_up_pericyte_high_low.RData")
