load("pseudobulk_obj_whole_meta_per.RData")
pseudobulk_obj_per <- pseudobulk_obj

load("pseudobulk_obj_whole_meta_others.RData")
pseudobulk_obj_others <- pseudobulk_obj

load("pseudobulk_obj_whole_meta_filtered_oct.RData")


pca_all <- PCA(pseudobulk_obj@assays$RNA$scale.data)

pca_per <- PCA(pseudobulk_obj_per@assays$RNA$scale.data)

pca_others <- PCA(pseudobulk_obj_others@assays$RNA$scale.data)


# Scree plots:

x_labels <- paste0("PC", seq(1,30))
  
scree_all <- barplotPercentVar(pca_all, returnGraph = TRUE, nPC = 30) + 
        ggtitle("Variance explained by each principal component",
          #NULL,
          subtitle = "All cell types") +
      geom_col(fill="#EF9AA8") +
      scale_x_discrete(limits = x_labels) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5),
        plot.title = element_text(size = 12)) +
  coord_cartesian(ylim = c(0, 20))

scree_per <- barplotPercentVar(pca_per, returnGraph = TRUE, nPC = 30) + 
      ggtitle(#"Variance explained by each principal component",
        NULL,
        subtitle = "Only pericytes") +
      geom_col(fill="#EF9AA8") +
      scale_x_discrete(limits = x_labels) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5),
        plot.title = element_text(size = 12)) +
  coord_cartesian(ylim = c(0, 20))

scree_others <- barplotPercentVar(pca_others, returnGraph = TRUE, nPC = 30) + 
      ggtitle(#"Variance explained by each principal component",
          NULL,
          subtitle = "Non-pericyte cell types") +
      geom_col(fill="#EF9AA8") +
      scale_x_discrete(limits = x_labels) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5),
        plot.title = element_text(size = 12)) +
  coord_cartesian(ylim = c(0, 20))



# Colors
colors <- c("slateblue1", "hotpink", "limegreen", "mediumorchid1","#B57457", "chocolate1","#1CC9BB", "goldenrod1", "#E0405B", "olivedrab2", "lightblue", "pink", "lightgreen", "thistle2","tan", "peachpuff", "#C0FAED", "lightgoldenrod","#EF9AA8", "darkolivegreen1")

colors_16 <- c("slateblue1", "hotpink", "limegreen", "mediumorchid1","#B57457", "chocolate1","#1CC9BB", "goldenrod1", "#E0405B", "olivedrab2", "lightblue", "thistle2", "peachpuff", "tan", "#C0FAED", "lightgoldenrod","#EF9AA8")

colors_pcaov <- c("proportion_oRG" = "thistle2", "proportion_vRG" = "lightblue", "proportion_panRG_O" = "peachpuff", "proportion_Per" = "lightgoldenrod", "proportion_PgS" = "#C0FAED", "Day_fixed" = "#E0405B", "prop_patient" = "olivedrab2", "donors_per_pool" = "slateblue1", "diff_method" = "hotpink", "R" = "chocolate1", "s" = "limegreen", "g2m" = "goldenrod1", "disease" = "tan", "gender" = "#1CC9BB", "pluri_novelty" = "#B57457", "rnaseq.sendai_reads" = "mediumorchid1")


# PCAoV plots

## Whole data
pca_aov_all <- PCaov(pca_all, colData = pseudobulk_obj@meta.data[, c("pluri_novelty", "gender", "rnaseq.sendai_reads", "diff_method", "R","s", "g2m", "Day_fixed", "disease", "donors_per_pool", "prop_patient",
                                                             "proportion_oRG", "proportion_vRG",
                                                             "proportion_panRG_O", "proportion_Per",
                                                             "proportion_PgS")])

pca_aov_all$PC <- factor(pca_aov_all$PC, levels=unique(pca_aov_all$PC))

pca_aov_all$feature <- factor(pca_aov_all$feature, levels = c("proportion_oRG", "proportion_vRG", "proportion_panRG_O", "proportion_Per", "proportion_PgS", "Day_fixed", "prop_patient", "donors_per_pool", "diff_method", "R", "s", "g2m", "disease", "gender", "pluri_novelty", "rnaseq.sendai_reads"))



pcaov_plot_all <- ggBorderedFactors(
    ggplot(pca_aov_all, aes(
        x = PC, y = sumSqPercent, fill = feature)) +
    geom_beeswarm(pch = 21, size = 2.5, cex = 1) +
    xlab("Principal component") +
      ylab("% sum of squares") +
      scale_fill_manual(values = colors_pcaov,
                        labels = c("oRG proportion", "vRG proportion", "panRG-O proportion", "Per proportion", "PgS proportion", "Time point", "Proportion of patient lines", "Pool size","Differentiation method", "R", "S score", "G2M score", "Disorder", "Gender", "Pluripotency score", "Number of sendai \n positive reads")) +
      guides(fill = guide_legend(title = "Feature")) +
  ggtitle("Analysis of variance for principal components",
   # NULL,
    subtitle = "All cell types") +
    theme(panel.grid.major.y = element_line(colour = "grey75"),
          panel.grid.minor.y = element_line(colour = "grey75"),
          panel.background = element_rect(fill = NA, colour = "black"),
          plot.title = element_text(size=12),
       #   axis.text.x = element_text(angle = 45, vjust = 0.5)
          )) +
    coord_cartesian(ylim = c(0, 70))



## Pericytes
pca_aov_per <- PCaov(pca_per, colData = pseudobulk_obj_per@meta.data[, c("pluri_novelty", "gender", "rnaseq.sendai_reads", "diff_method", "R","s", "g2m", "Day_fixed", "disease", "donors_per_pool", "prop_patient")])

pca_aov_per$PC <- factor(pca_aov_per$PC, levels=unique(pca_aov_per$PC))

pca_aov_per$feature <- factor(pca_aov_per$feature, levels = c("Day_fixed", "prop_patient", "donors_per_pool", "diff_method", "R", "s", "g2m", "disease", "gender", "pluri_novelty", "rnaseq.sendai_reads"))

pcaov_plot_per <- ggBorderedFactors(
    ggplot(pca_aov_per, aes(
        x = PC, y = sumSqPercent, fill = feature)) +
    geom_beeswarm(pch = 21, size = 2.5, cex = 1) +
    xlab("Principal component") +
      ylab("% sum of squares") +
      scale_fill_manual(values = colors_pcaov,
                        labels = c("Time point", "Proportion of patient lines", "Pool size","Differentiation method", "R", "S score", "G2M score", "Disorder", "Gender", "Pluripotency score", "Number of sendai \n positive reads")) +
      guides(fill = guide_legend(title = "Feature")) +
  ggtitle(#"Analysis of variance for principal components",
    NULL,
    subtitle = "Only pericytes") +
    theme(panel.grid.major.y = element_line(colour = "grey75"),
          panel.grid.minor.y = element_line(colour = "grey75"),
          panel.background = element_rect(fill = NA, colour = "black"),
          plot.title = element_text(size=12),
       #   axis.text.x = element_text(angle = 45, vjust = 0.5)
          )) +
    coord_cartesian(ylim = c(0, 70))




## Non-pericytes
pca_aov_others <- PCaov(pca_others, colData = pseudobulk_obj_others@meta.data[, c("pluri_novelty", "gender", "rnaseq.sendai_reads", "diff_method", "R","s", "g2m", "Day_fixed", "disease", "donors_per_pool", "prop_patient")])

pca_aov_others$PC <- factor(pca_aov_others$PC, levels=unique(pca_aov_others$PC))

pca_aov_others$feature <- factor(pca_aov_others$feature, levels = c("Day_fixed", "prop_patient", "donors_per_pool", "diff_method", "R", "s", "g2m", "disease", "gender", "pluri_novelty", "rnaseq.sendai_reads"))

pcaov_plot_others <- ggBorderedFactors(
    ggplot(pca_aov_others, aes(
        x = PC, y = sumSqPercent, fill = feature)) +
    geom_beeswarm(pch = 21, size = 2.5, cex = 1) +
    xlab("Principal component") +
      ylab("% sum of squares") +
      scale_fill_manual(values = colors_pcaov,
                        labels = c("Time point", "Proportion of patient lines", "Pool size","Differentiation method", "R", "S score", "G2M score", "Disorder", "Gender", "Pluripotency score", "Number of sendai \n positive reads")) +
      guides(fill = guide_legend(title = "Feature")) +
  ggtitle(#"Analysis of variance for principal components",
    NULL,
    subtitle = "Non-pericyte cell types") +
    theme(panel.grid.major.y = element_line(colour = "grey75"),
          panel.grid.minor.y = element_line(colour = "grey75"),
          panel.background = element_rect(fill = NA, colour = "black"),
          plot.title = element_text(size=12),
         # axis.text.x = element_text(angle = 45, vjust = 0.5)
          )) +
    coord_cartesian(ylim = c(0, 70))



# Create a shared legend for all PCAoV plots
leg <- ggdraw(get_legend(pcaov_plot_all +
                           guides(fill = guide_legend(title = "Feature", nrow = 8)) +
                           theme(legend.key.size = unit(0.4, 'cm'),
                                 legend.title = element_text(size=10),
                                 legend.text = element_text(size=9),
                                 legend.key = element_blank(),
                              #   legend.background = element_rect(fill = "lavenderblush"),
                                plot.margin = margin(t = 0,
                             r = 0, 
                             b = 0, 
                             l = 0.2,
                             unit = "cm")),
                       #  "bottom"))
                       "right"))

# Arrange the plots and ensure same sizes

plotlist <- list(a = scree_all +
                            theme(axis.title.x=element_blank(), plot.margin = margin(t = 0.1, 
                             r = 1.2,  
                             b = 0,  
                             l = 0.1,  
                             unit = "cm"),
                           #  plot.background = element_rect(fill = 'lightblue')
                           ),
                          b = pcaov_plot_all +
                            theme(axis.title.x=element_blank(), plot.margin = margin(t = 0.1, 
                             r = 0.1,  
                             b = 0, 
                             l = 0.2, 
                             unit = "cm"),
                         #    plot.background = element_rect(fill = 'pink')
                         ) +
                     guides(fill = "none"),
                          c = scree_per +
                            theme(axis.title.x=element_blank(), plot.margin = margin(t = 0,  
                             r = 1.2, 
                             b = 0,  
                             l = 0.1, 
                             unit = "cm"),
                          #   plot.background = element_rect(fill = 'lavender')
                          ),
                          d = pcaov_plot_per +
                            theme(axis.title.x=element_blank(), plot.margin = margin(t = 0, 
                             r = 0.1,  
                             b = 0,  
                             l = 0.2, 
                             unit = "cm"),
                           #  plot.background = element_rect(fill = 'peachpuff')
                           ) +
                     guides(fill = "none"),
                          e = scree_others +
                            theme(plot.margin = margin(t = 0,
                             r = 1.2,  
                             b = 0,  
                             l = 0.1,
                             unit = "cm"),
                          #   plot.background = element_rect(fill = 'thistle1')
                             ),
                          f = pcaov_plot_others +
                            theme(plot.margin =margin(t = 0,
                             r = 0.1, 
                             b = 0, 
                             l = 0.2,
                             unit = "cm"),
                         #    plot.background = element_rect(fill = 'darkolivegreen1')
                         ) +
                     guides(fill = "none"),
                          g = leg)


pgs <- lapply(plotlist, function(plot) {ggplotGrob(plot)})
pg_heights <- lapply(pgs, function(pg) {pg$heights})

maxHeight = grid::unit.pmax(pg_heights$a, pg_heights$b, pg_heights$c, pg_heights$d, pg_heights$e, pg_heights$f)

pgs_new <- lapply(pgs, function(pg) {
  pg$heights <- as.list(maxHeight)
  pg
  })


pcaov_plots <- ggarrange(plotlist = list(pgs_new$a, pgs_new$b, pgs_new$c, pgs_new$d, pgs_new$e, pgs_new$f, NULL, pgs$g), ncol = 2, nrow = 4, labels = c("A", "B", "", "", "", "", "", ""), heights = c(1,1,1,0.55))

pdf(pcaov_plots, file = "	Pseudobulk_PCA_plots.pdf", paper = "a4")
