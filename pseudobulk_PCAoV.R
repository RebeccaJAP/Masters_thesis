  # Load packages
library(dplyr)
library(ggplot2)
library(oob)
library(Seurat)
library(ggbeeswarm)
library(cowplot)
library(ggpubr)


  # Get the necessary data and change names to avoid overwriting variables
  # Pericytes
load("Pseudobulk_data/pseudobulk_obj_pericytes_filtered.RData")
pseudobulk_obj_per <- pseudobulk_obj

  # Non-pericytes
load("Pseudobulk_data/pseudobulk_obj_others_filtered.RData")
pseudobulk_obj_others <- pseudobulk_obj

  # All cell types
load("Pseudobulk_data/pseudobulk_obj_filtered.RData")


	# Find the principal components
pca_all <- PCA(pseudobulk_obj@assays$RNA$scale.data)
pca_per <- PCA(pseudobulk_obj_per@assays$RNA$scale.data)
pca_others <- PCA(pseudobulk_obj_others@assays$RNA$scale.data)


	# Create the elbow plots

x_labels <- paste0("PC", seq(1,30))
  
scree_all <- barplotPercentVar(pca_all, returnGraph = TRUE, nPC = 30) + 
  ggtitle("Variance explained by each principal component",
          subtitle = "All cell types") +
  geom_col(fill="#EF9AA8") +
  scale_x_discrete(limits = x_labels) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5),
        plot.title = element_text(size = 12)) +
  coord_cartesian(ylim = c(0, 20))
      
scree_per <- barplotPercentVar(pca_per, returnGraph = TRUE, nPC = 30) + 
  ggtitle(NULL,
          subtitle = "Only pericytes") +
  geom_col(fill="#EF9AA8") +
  scale_x_discrete(limits = x_labels) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5),
        plot.title = element_text(size = 12)) +
  coord_cartesian(ylim = c(0, 20))
      
scree_others <- barplotPercentVar(pca_others, returnGraph = TRUE, nPC = 30) + 
  ggtitle(NULL,
          subtitle = "Non-pericyte cell types") +
  geom_col(fill="#EF9AA8") +
  scale_x_discrete(limits = x_labels) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5),
        plot.title = element_text(size = 12)) +
  coord_cartesian(ylim = c(0, 20))
      
      
  # Set the colors
colors_pcaov <- c("proportion_oRG" = "thistle2", "proportion_vRG" = "lightblue", "proportion_panRG_O" = "peachpuff",
                  "proportion_Per" = "lightgoldenrod", "proportion_PgS" = "#C0FAED", "Day_fixed" = "#E0405B",
                  "prop_patient" = "olivedrab2", "donors_per_pool" = "slateblue1", "diff_method" = "hotpink", "R" = "chocolate1",
                  "s" = "limegreen", "g2m" = "goldenrod1", "disease" = "tan", "gender" = "#1CC9BB", "pluri_novelty" = "#B57457",
                  "rnaseq.sendai_reads" = "mediumorchid1")
      

  # Define the function for PCAoV
  # Modified from function PCaov of the oob package
PCaov <- function(pca, colData, nComponent = 10){
  aovFormula <- paste0(" ~ ", paste0(rev(colnames(colData)), 
                                     collapse = " + "))
  modelData <- data.frame(colData, pca$x)
  nVar <- ncol(colData)
  anovPerPCdt <- do.call("rbind", lapply(colnames(pca$x[, seq_len(nComponent)]),
                                         function(PC) {
                                           model <- lm(formula(paste0(PC, aovFormula)), data = modelData)
                                           aov_results <- car::Anova(model, type = 2)
                                           pvals <- aov_results$`Pr(>F)`[seq_len(nVar)]
                                           sumsq <- aov_results$`Sum Sq`[seq_len(nVar)]
                                           sumSqPercent <- sumsq / sum(aov_results$`Sum Sq`, na.rm = TRUE) * 100
                                           data.frame(feature = rev(colnames(colData)), PC,
                                                      sumSqPercent, pval = pvals)
                                           }))
  anovPerPCdt$feature <- as.factor(anovPerPCdt$feature)
  anovPerPCdt$padj <- p.adjust(anovPerPCdt$pval)
  anovPerPCdt$PC <- paste0("PC", formatNumber2Character(substr(anovPerPCdt$PC,
                                                               3, nchar(anovPerPCdt$PC))))
  return(anovPerPCdt)
  }
      
      
  # Run and plot the PCAoV for the entire pseudobulk dataset
pca_aov_all <- PCaov(pca_all,
                     colData = pseudobulk_obj@meta.data[, c("pluri_novelty", "gender", "rnaseq.sendai_reads", "diff_method",
                                                            "R","s", "g2m", "Day_fixed", "disease", "donors_per_pool", "prop_patient",
                                                            "proportion_oRG", "proportion_vRG", "proportion_panRG_O", "proportion_Per",
                                                            "proportion_PgS")])
      
pca_aov_all$PC <- factor(pca_aov_all$PC, levels=unique(pca_aov_all$PC))
      
pca_aov_all$feature <- factor(pca_aov_all$feature, levels = c("proportion_oRG", "proportion_vRG", "proportion_panRG_O",
                                                              "proportion_Per", "proportion_PgS", "Day_fixed", "prop_patient",
                                                              "donors_per_pool", "diff_method", "R", "s", "g2m", "disease", "gender",
                                                              "pluri_novelty", "rnaseq.sendai_reads"))
      
pcaov_plot_all <- ggBorderedFactors(
  ggplot(pca_aov_all, aes(x = PC, y = sumSqPercent, fill = feature)) +
    geom_beeswarm(pch = 21, size = 2.5, cex = 1) +
    xlab("Principal component") +
    ylab("% sum of squares") +
    scale_fill_manual(values = colors_pcaov,
                      labels = c("oRG proportion", "vRG proportion", "panRG-O proportion", "Per proportion", "PgS proportion",
                                 "Time point", "Proportion of patient lines", "Pool size","Differentiation method", "R", "S score",
                                 "G2M score", "Disorder", "Gender", "Pluripotency novelty score", "Number of sendai \n positive reads")) +
    guides(fill = guide_legend(title = "Feature")) +
    ggtitle("Analysis of variance for principal components",
            subtitle = "All cell types") +
    theme(panel.grid.major.y = element_line(colour = "grey75"),
          panel.grid.minor.y = element_line(colour = "grey75"),
          panel.background = element_rect(fill = NA, colour = "black"),
          plot.title = element_text(size=12),
          )) +
  coord_cartesian(ylim = c(0, 70))

      
  # Run and plot the PCAoV for the pericytes
pca_aov_per <- PCaov(pca_per, colData = pseudobulk_obj_per@meta.data[, c("pluri_novelty", "gender", "rnaseq.sendai_reads",
                                                                         "diff_method", "R","s", "g2m", "Day_fixed", "disease",
                                                                         "donors_per_pool", "prop_patient")])
      
pca_aov_per$PC <- factor(pca_aov_per$PC, levels=unique(pca_aov_per$PC))
      
pca_aov_per$feature <- factor(pca_aov_per$feature, levels = c("Day_fixed", "prop_patient", "donors_per_pool", "diff_method", "R", "s",
                                                              "g2m", "disease", "gender", "pluri_novelty", "rnaseq.sendai_reads"))
      
pcaov_plot_per <- ggBorderedFactors(
  ggplot(pca_aov_per, aes(
    x = PC, y = sumSqPercent, fill = feature)) +
    geom_beeswarm(pch = 21, size = 2.5, cex = 1) +
    xlab("Principal component") +
    ylab("% sum of squares") +
    scale_fill_manual(values = colors_pcaov,
                      labels = c("Time point", "Proportion of patient lines", "Pool size", "Differentiation method", "R", "S score",
                                 "G2M score", "Disorder", "Gender", "Pluripotency novelty score", "Number of sendai \n positive reads")) +
    guides(fill = guide_legend(title = "Feature")) +
    ggtitle(NULL,
            subtitle = "Only pericytes") +
    theme(panel.grid.major.y = element_line(colour = "grey75"),
          panel.grid.minor.y = element_line(colour = "grey75"),
          panel.background = element_rect(fill = NA, colour = "black"),
          plot.title = element_text(size=12)
          )) +
  coord_cartesian(ylim = c(0, 70))
 

  # Run and plot the PCAoV for the pericytes
pca_aov_others <- PCaov(pca_others, colData = pseudobulk_obj_others@meta.data[, c("pluri_novelty", "gender", "rnaseq.sendai_reads",
                                                                                  "diff_method", "R","s", "g2m", "Day_fixed", "disease",
                                                                                  "donors_per_pool", "prop_patient")])

pca_aov_others$PC <- factor(pca_aov_others$PC, levels=unique(pca_aov_others$PC))
      
pca_aov_others$feature <- factor(pca_aov_others$feature, levels = c("Day_fixed", "prop_patient", "donors_per_pool", "diff_method", "R",
                                                                    "s", "g2m", "disease", "gender", "pluri_novelty",
                                                                    "rnaseq.sendai_reads"))
      
pcaov_plot_others <- ggBorderedFactors(
  ggplot(pca_aov_others, aes(x = PC, y = sumSqPercent, fill = feature)) +
    geom_beeswarm(pch = 21, size = 2.5, cex = 1) +
    xlab("Principal component") +
    ylab("% sum of squares") +
    scale_fill_manual(values = colors_pcaov,
                      labels = c("Time point", "Proportion of patient lines", "Pool size","Differentiation method", "R", "S score",
                                 "G2M score", "Disorder", "Gender", "Pluripotency novelty score", "Number of sendai \n positive reads")) +
    guides(fill = guide_legend(title = "Feature")) +
    ggtitle(NULL,
            subtitle = "Non-pericyte cell types") +
    theme(panel.grid.major.y = element_line(colour = "grey75"),
          panel.grid.minor.y = element_line(colour = "grey75"),
          panel.background = element_rect(fill = NA, colour = "black"),
          plot.title = element_text(size=12)
          )) +
  coord_cartesian(ylim = c(0, 70))
      
  
  # Create a shared legend
leg <- ggdraw(ggpubr::get_legend(pcaov_plot_all +
                           guides(fill = guide_legend(title = "Feature", nrow = 8)) +
                           theme(legend.key.size = unit(0.4, 'cm'),
                                 legend.title = element_text(size=10),
                                 legend.text = element_text(size=9),
                                 legend.key = element_blank(),
                                 plot.margin = margin(t = 0, r = 0, b = 0, l = 0.2, unit = "cm")),
                         "right"))

  # Plot all figures and the shared legend together
plotlist <- list(a = scree_all +
                   theme(axis.title.x=element_blank(),
                         plot.margin = margin(t = 0.1, r = 1.2, b = 0, l = 0.1, unit = "cm")),
                 b = pcaov_plot_all +
                   theme(axis.title.x=element_blank(),
                         plot.margin = margin(t = 0.1, r = 0.1, b = 0, l = 0.2, unit = "cm")) +
                   guides(fill = "none"),
                 c = scree_per +
                   theme(axis.title.x=element_blank(),
                         plot.margin = margin(t = 0, r = 1.2, b = 0, l = 0.1, unit = "cm")),
                 d = pcaov_plot_per +
                   theme(axis.title.x=element_blank(),
                         plot.margin = margin(t = 0, r = 0.1, b = 0, l = 0.2, unit = "cm")) +
                   guides(fill = "none"),
                 e = scree_others +
                   theme(plot.margin = margin(t = 0, r = 1.2, b = 0, l = 0.1, unit = "cm")),
                 f = pcaov_plot_others +
                   theme(plot.margin = margin(t = 0, r = 0.1, b = 0, l = 0.2, unit = "cm")) +
                   guides(fill = "none"),
                 g = leg)
      
pgs <- lapply(plotlist, function(plot) {ggplotGrob(plot)})

  # Get the height of each plot and find the maximum height
pg_heights <- lapply(pgs, function(pg) {pg$heights})
maxHeight = grid::unit.pmax(pg_heights$a, pg_heights$b, pg_heights$c, pg_heights$d, pg_heights$e, pg_heights$f)

  # Set the height of each plot the maximum found to ensure even heights
pgs_new <- lapply(pgs, function(pg) {
  pg$heights <- as.list(maxHeight)
  pg
  })
      
ggarrange(plotlist = list(pgs_new$a, pgs_new$b, pgs_new$c, pgs_new$d, pgs_new$e, pgs_new$f, NULL, pgs$g),
          ncol = 2,
          nrow = 4,
          labels = c("A", "B", "", "", "", "", "", ""),
          heights = c(1,1,1,0.55))
    
      
      
      
      
      
      
