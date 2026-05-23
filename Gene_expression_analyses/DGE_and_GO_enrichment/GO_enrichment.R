  # Set time point
day <- 20
#day <- 60

  # Pericytes vs. others or pericyte-high vs. pericyte-low analysis?
#comparison <- "pericytes_others"
 comparison <- "high_low"

  # Load packages
library(clusterProfiler)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(dplyr)
library(Seurat)

  # Load the up- and down-regulated genes and the pseudobulk object
if (comparison == "pericytes_others") {
  load(paste0("Gene_expression_analyses/DGE_and_GO_enrichment/Pericytes_others/limma_up_pericyte_nonpericyte_d", day, ".RData"))
  load(paste0("Gene_expression_analyses/DGE_and_GO_enrichment/Pericytes_others/limma_down_pericyte_nonpericyte_d", day, ".RData"))
  load("Pseudobulk_data/pseudobulk_obj_pericytes_others_filtered.RData")
  title_end <- paste0("pericytes on day ", day)
  } else {
  load("Gene_expression_analyses/DGE_and_GO_enrichment/Pericyte_high_low/limma_up_pericyte_high_low.RData"
  load("Gene_expression_analyses/DGE_and_GO_enrichment/Pericyte_high_low/limma_down_pericyte_high_low.RData")
  load("Pseudobulk_data/pseudobulk_obj_filtered.RData")
  title_end <- "pericyte-high cell lines on day 20"
  }
  
  # Get all the genes in the pseudobulk object and set as the background
background <- unique(rownames(pseudobulk_obj@assays$RNA$counts))

  # Do the GO enrichment analysis for the upregulated genes
GO_results_up <- enrichGO(unlist(upregulated_genes),
                          OrgDb = "org.Hs.eg.db",
                          keyType = "SYMBOL",
                          ont = "BP",
                          universe = background,
                          pvalueCutoff = 0.05)

GO_df_up <- as.data.frame(GO_results_up)

  # Simplify the results, i.e. get less overlapping results
GO_up_simple <- simplify(
  GO_results_up,
  cutoff = 0.7,
  by = "p.adjust",
  select_fun = min,
  measure = "Wang",
  semData = NULL
)

GO_df_up_simple <- as.data.frame(GO_up_simple)


  # Do the GO enrichment analysis for the downregulated genes
GO_results_down <- enrichGO(unlist(downregulated_genes),
                          OrgDb = "org.Hs.eg.db",
                          keyType = "SYMBOL",
                          ont = "BP",
                          universe = background,
                          pvalueCutoff = 0.05)

GO_df_down <- as.data.frame(GO_results_down)

# Simplify the results, i.e. get less overlapping results
GO_down_simple <- simplify(
  GO_results_down,
  cutoff = 0.7,
  by = "p.adjust",
  select_fun = min,
  measure = "Wang",
  semData = NULL
)

GO_df_down_simple <- as.data.frame(GO_down_simple)

  # Plot the results
barplot_up_simple <- barplot(GO_up_simple %>%
                               arrange(desc(Count)),
                             showCategory = 25,
                             font.size = 10,
                             title = paste0("Biological processes related to upregulated genes in \n", title_end),
                             label_format = 100,
                             color = "p.adjust") +
  scale_fill_gradientn(colors=c("goldenrod1","hotpink"))


barplot_down_simple <- barplot(GO_down_simple %>%
                               arrange(desc(Count)),
                             showCategory = 25,
                             font.size = 10,
                             title = paste0("Biological processes related to downregulated genes in \n", title_end),
                             label_format = 100,
                             color = "p.adjust") +
  scale_fill_gradientn(colors=c("goldenrod1","hotpink"))

  # Set the filename depending on the analysis
if (comparison == "pericytes_others") {
  filename_up <- paste0("Gene_expression_analyses/DGE_and_GO_enrichment/Pericytes_others/GO_barplot_up_d", day, ".pdf")
  filename_down <- paste0("Gene_expression_analyses/DGE_and_GO_enrichment/Pericytes_others/GO_barplot_down_d", day, ".pdf")
} else {
  filename_up <- "Gene_expression_analyses/DGE_and_GO_enrichment/Pericyte_high_low/GO_barplot_up_d20.pdf"
  filename_down <- "Gene_expression_analyses/DGE_and_GO_enrichment/Pericyte_high_low/GO_barplot_down_d20.pdf"
}

  # Save the barplots
pdf(file = filename_up, paper = "a5r")
plot(barplot_up_simple)
dev.off()

  # No results for d20 pericytes vs. others -> do not try to save
if (comparison != "pericytes_others" | day != 20) {
  pdf(file = filename_down, paper = "a5r")
  plot(barplot_down_simple)
  dev.off()
}
