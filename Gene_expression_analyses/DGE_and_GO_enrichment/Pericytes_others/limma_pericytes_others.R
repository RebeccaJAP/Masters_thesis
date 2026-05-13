  # Load packages
library(ComplexHeatmap)
library(rlist)
library(gridtext)
library(ggplot2)
library(edgeR)
library(tidyr)
library(dplyr)
library(Seurat)
library(ggrepel)


  # Set the timepoint of interest
day <- 20
#day <- 60

  # Adjust the cutoffs if needed
padj_cutoff <- 0.05

  # How many pools does a gene need to be DE to be of interest?
sample_cutoff <- ifelse(day == 20, 10, 4)

  # Get the filtered pseudobulk object where pericytes and non-pericytes are separated
load("Pseudobulk_data/pseudobulk_obj_pericytes_others_filtered.RData")

  # Get the pools with at least one pericyte sample and one non-pericyte sample
data <- pseudobulk_obj@meta.data %>%
  filter(diff_method != "indv" & Day_fixed == day) %>%
  group_by(SampleID) %>%
  mutate(n=n_distinct(cell_type_group)) %>%
  ungroup() %>%
  filter(n > 1)

valid_samples <- unique(data$SampleID)

  # Get the relevant metadata and separate by pools
meta_by_sample <- lapply(valid_samples, function(s) {
  pseudobulk_obj@meta.data[pseudobulk_obj$SampleID == s,]
})

names(meta_by_sample) <- valid_samples

  # Get the relevant expression data and separate by pools
expr_by_sample <- lapply(valid_samples, function(s) {
  cols <- rownames(pseudobulk_obj@meta.data[pseudobulk_obj$SampleID == s & pseudobulk_obj$Day_fixed == day,])
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
  cell_category_ <- meta_by_sample[[i]]$cell_type_group
  cell_category_ <- factor(cell_category_, levels = c("pericyte", "nonPericyte"))
  mm <- model.matrix(~0 + cell_category_)
  mm_all[[i]] <- mm
}

names(mm_all) <- valid_samples


  # Run voom, contrast pericytes and non-pericytes, get DEGs
toptables <- c()

for (i in seq(1, length(valid_samples))) {
  y <- voom(d0[[i]], mm_all[[i]], plot = T)
  fit <- lmFit(y, mm_all[[i]])
  contr <- makeContrasts(cell_category_pericyte - cell_category_nonPericyte, levels = colnames(coef(fit)))
  tmp <- contrasts.fit(fit, contr)
  tmp <- eBayes(tmp)
  top.table <- topTable(tmp, sort.by = "P", n = Inf)
  toptables[[i]] <- top.table
}

names(toptables) <- valid_samples


  # Check how well the pericyte and non-pericyte samples separate
for (i in seq(1, length(valid_samples))) {
  
  mds <- plotMDS(d0[[i]], labels=meta_by_sample[[i]]$donor, plot = F)
  
  mds_df <- data.frame(
    Dim1 = mds$x,
    Dim2 = mds$y,
    Color = meta_by_sample[[i]]$cell_type_group,
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


  # Get the names of differentially expressed genes for each pool
signif_DE_genes <- lapply(toptables, function(sample) {
  sample$gene <- rownames(sample)
  sample$gene[which(sample$adj.P.Val < padj_cutoff)]
})


  # Limit the limma results to statistically significant ones
toptables_signif <- lapply(toptables, function(sample) {
  sample$gene <- rownames(sample)
  sample[which(sample$adj.P.Val < padj_cutoff),]
})


  # How many pools each gene is differentially expressed in?
gene_counts_limma <- sort(table(unlist(signif_DE_genes)), decreasing = T)

  # All significant DEGs in one list
all_de_genes <- names(gene_counts_limma)


  # Is gene differentailly expressed per pool? TRUE/FALSE value for each pool
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
  mutate(n_present = sum(present),
         dataset = as.factor(sub("\\_d..", "", dataset)),
         abs_logFC = min(abs(logFC), na.rm = T)) %>%
  arrange(desc(abs_logFC))

  # Set a limit for how much a gene needs to be differentially expressed to be plotted
logFC_limit <- ifelse(day == 20, 0, 3)

  # Plot heatmap
  # Note: produces warning about missing points for the non-significant genes;
  # these are missing on purpose and the warning should be ignored

heatmap_DGE <- ggplot((heatmap_df_cutoff %>% filter(n_present >= sample_cutoff & abs_logFC > logFC_limit)), aes(x = gene, y = dataset)) +
  geom_point(aes(color = logFC, size = abs(logFC))) +
  scale_y_discrete(drop = FALSE) +
  scale_color_gradientn(colors=c("#1600A3", "slateblue1", "white", "hotpink", "#D10069"), rescaler = ~ scales::rescale_mid(.x, mid = 0)) +
  theme_minimal() +
  labs(x = "Gene", y = "Pool", color = "log2FC", size = "Absolute value of log2FC") +
  ggtitle(paste0("DGE in pericytes vs. non-pericyte cells on day ", day)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size=8),
        title = element_text(size = 12))

heatmap_DGE

  # Save the heatmap
pdf(file = paste0("Gene_expression_analyses/DGE_and_GO_enrichment/Pericytes_others/Pericytes_others_DGE_heatmap_d", day, ".pdf"),
    width = 8.27,
    height = 4.5)
plot(heatmap_DGE)
dev.off()


  # Get genes for GO enrichment analysis:

  # Set a limit for how many pools the gene needs to be significantly differentially expressed in
  # to be included in GO enrichment analysis
GO_limit <- ifelse(day == 20, 6, 3)

  # Get genes which are significantly differentially expressed in enough pools
  # and which are not significantly upregulated in any pool
downregulated_genes <- heatmap_df %>%
  group_by(gene) %>%
  mutate(n_present = sum(present), max_logFC = max(logFC, na.rm = T)) %>%
  filter(n_present >= GO_limit & max_logFC < 0) %>%
  distinct(gene)

  # Get genes which are significantly differentially expressed in enough pools
  # and which are not significantly downregulated in any pool
upregulated_genes <- heatmap_df %>%
  group_by(gene) %>%
  mutate(n_present = sum(present), min_logFC = min(logFC, na.rm = T)) %>%
  filter(n_present >= GO_limit & min_logFC > 0) %>%
  distinct(gene)

  # Save the genes for further analysis
save(downregulated_genes,
     file = paste0("Gene_expression_analyses/DGE_and_GO_enrichment/Pericytes_others/limma_down_pericyte_nonpericyte_d", day, ".RData"))
save(upregulated_genes,
     file = paste0("Gene_expression_analyses/DGE_and_GO_enrichment/Pericytes_others/limma_up_pericyte_nonpericyte_d", day, ".RData"))


                     
  # Create an upset plot:

  # Set color palette
palette_all <- list.reverse(c("slateblue1", "hotpink", "limegreen", "mediumorchid1","#B57457", "chocolate1","#1CC9BB",  "goldenrod1", "#E0405B", "olivedrab2"))

  # Remove the time point from pool name
colnames(binary_mat) <- sub("\\_d..", "", colnames(binary_mat))

  # Create and normalize a matrix with intersections of all differentially expressed genes between different pool combination
comb_mat <- make_comb_mat(binary_mat, mode = "intersect")
comb_mat <- normalize_comb_mat(comb_mat, full_comb_sets = T)

  # On day 20, smaller pool sets should not be included to keep the plot clear
if (day == 20) {
  min_pool_count <- 8
  comb_mat <- comb_mat[comb_degree(comb_mat) >= min_pool_count]
} else {
  min_pool_count <- 1
}

  # Get the combination sizes, i.e. how many differentially expressed genes per pool combination
cs = comb_size(comb_mat)

  # Get the set sizes, i.e. how many differentially expressed genes per pool
ss = set_size(comb_mat)

  # Create a data frame combining above information for ordering the combinations
order_df <- data.frame(name = comb_name(comb_mat), degree = comb_degree(comb_mat), size = cs)


  # Plot the upset plot
upset <- UpSet(comb_mat,
               pt_size = unit(2, "mm"),
               lwd = 1,
               column_title = paste0("Pericytes vs. other cell types, day ", day),
               comb_col = palette_all[comb_degree(comb_mat)],
               comb_order = order(-order_df$degree, order_df$size),
               bg_pt_col = "grey58",
               row_names_side = "left",
               row_names_gp = grid::gpar(fontsize = 9),
               bottom_annotation = HeatmapAnnotation("Number of pools" = factor(comb_degree(comb_mat),
                                                                                levels = as.character(min_pool_count:length(unique(data$SampleID)))),
                                                     col = list("Number of pools" = c("1" = palette_all[1],
                                                                                      "2" = palette_all[2],
                                                                                      "3" = palette_all[3],
                                                                                      "4" = palette_all[4],
                                                                                      "5"=palette_all[5],
                                                                                      "6"=palette_all[6],
                                                                                      "7"=palette_all[7],
                                                                                      "8" = palette_all[8],
                                                                                      "9" = palette_all[9],
                                                                                      "10" = palette_all[10])),
                                                     show_annotation_name = F),
               top_annotation = HeatmapAnnotation("Gene intersection size" = anno_barplot(cs,
                                                                                          gp = gpar(fill = palette_all[comb_degree(comb_mat)],
                                                                                                    col = palette_all[comb_degree(comb_mat)]), 
                                                                                          border = F,
                                                                                          add_numbers = T,
                                                                                          numbers_gp = gpar(fontsize = 8),
                                                                                          numbers_rot = 90,
                                                                                          height = unit(4, "cm")),
                                                  annotation_name_side = c("left"),
                                                  annotation_name_gp = gpar(fontsize = 9)),
               right_annotation = rowAnnotation("Differentially expressed \ngenes per pool" = anno_barplot(ss,
                                                                                                           gp = gpar(fill = "black",
                                                                                                                     col="black"),
                                                                                                           axis_param = list(
                                                                                                             labels = ss,
                                                                                                             labels_rot = 0),
                                                                                                           border = F,
                                                                                                           width = unit(4, "cm"),
                                                                                                           add_numbers = T,
                                                                                                           numbers_gp = gpar(fontsize = 8)),
                                                annotation_name_gp = gpar(fontsize = 9)))


upset

  # Save the upset plot
pdf(file = paste0("Gene_expression_analyses/DGE_and_GO_enrichment/Pericytes_others/Pericytes_others_DGE_upset_d", day, ".pdf"),
    width = 11,
    height = 5)
plot(upset)
dev.off()
