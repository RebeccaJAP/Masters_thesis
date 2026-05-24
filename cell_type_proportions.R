  # Load packages
library(ggplot2)
library(dplyr)

  # Load the full metadata
load("joined_meta.RData")

  # Compute the cell type proportions per pool, time point and donor
  # and prepare the labels
cell_type_proportion <- joined_meta %>%
  mutate(label = gsub("_with", " with", paste0(pool, "\n", "day ", Day_fixed), ignore.case = T)) %>%
  group_by(SampleID, pool, label, donor, cell_type_annot) %>%
  mutate(donor = case_when(grepl("HEL",donor) ~ donor, .default = gsub("^[^-]*-", "", donor))) %>%
  summarise(freq = n()) %>%
  mutate(proportion = freq/sum(freq), .keep="unused")

  # Assign colors to cell types
palette = c("Per" = "slateblue1",
            "vRG" = "peachpuff",
            "oRG" = "pink",
            "panRG-O" = "lightgoldenrod",
            "hRG-O" = "#EF9AA8",
            "LQ-RG-O" = "thistle2",
            
            "IP" = "#1CC9BB",
            "PgS" = "goldenrod1",
            "PgG2M" = "mediumorchid1",
            
            "ExDp1" = "darkolivegreen1",
            "ExPanNeu-O" = "forestgreen",
            "ExNeuNew-O" = "#C0FAED",
            "ExDp2" = "olivedrab2",
            "ExN" =  "lightgreen",
            "ExU-O" = "limegreen",
            
            "InCGE" = "#E0405B",
            "InMGE" = "hotpink",
            
            "End" = "lightblue",
            
            "AstroHindb-O" = "salmon2",
            "OPC" = "brown2",
            "Mic" = "chocolate1",
            
            
            "Unmapped" = "tan",
            "Others" = "#B57457"
            
)


  # Set the cell type order
cell_types <- c("Per", "vRG", "oRG", "panRG-O", "hRG-O", "LQ-RG-O", "IP", "PgS", "PgG2M", "ExDp1", "ExPanNeu-O", "ExNeuNew-O", "ExDp2", "ExN", "ExU-O", "InCGE", "InMGE", "End", "AstroHindb-O", "OPC", "Mic", "Unmapped", "Others")

  # Create the plot
cell_type_prop_plot <- ggplot(cell_type_proportion, aes(fill=factor(cell_type_annot, levels = cell_types), y=proportion, x=donor)) +
  geom_bar(position="stack", stat="identity") +
  scale_fill_manual(values = palette, breaks = cell_types) + 
  scale_x_discrete(guide=guide_axis(angle=90)) +
  facet_wrap(~label, scales = "free_x", labeller = label_wrap_gen(10)) +
  ggtitle("Proportions of cell types per pool, day, and donor") +
  theme(axis.text.x = element_text(size = 8) ,
        axis.text.y = element_text(size = 8),
        axis.title = element_text(size = 8),
        plot.title = element_text(size = 10),
        strip.text = element_text(size = 7.5, margin = margin(t = 0.5, r = 0.5, b = 0.5, l = 0.5)),
        legend.title = element_text(size = 8), 
        legend.text = element_text(size = 8),
        legend.key.size = unit(0.4, 'cm'),
        legend.margin = margin(t = 0, r = 0, b = 0, l = -8),
        legend.position.inside = c(.85, .03),
        panel.spacing.x = unit(0.1, "cm")) +
  guides(fill = guide_legend(position = "inside", ncol = 2, title.position = "left")) +
  labs(fill = "Cell type") +
  ylab("Proportion of donor cells") +
  xlab("Donor")

pdf(file=paste0("cell_type_proportions_plot.pdf"), paper = "a4")
plot(cell_type_prop_plot)
dev.off()
