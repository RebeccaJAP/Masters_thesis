  # Get the filtered pseudobulk object
load("Pseudobulk_data/pseudobulk_obj_filtered.RData")

  # Exclude individual lines and only look into day 20 data
data <- pseudobulk_obj@meta.data %>%
  filter(diff_method != "indv" & Day_fixed == 20)

  # Get the maximum pericyte proportion of the successful line on day 38
Sanger_ref <- pseudobulk_obj@meta.data %>%
  filter(SampleID == "OT_G31139991327") %>%
  select(SampleID, proportion_Per, Day_fixed)

Sanger_limit <- max(Sanger_ref$proportion_Per)

  # Get the combined mean proportion of cell types corresponding to pericytes in the fetal data of Braun et al.
Braun_pericyte_props <- read.table("Gene_expression_analyses/DGE_and_GO_enrichment/Pericyte_high_low/combined_proportions_per_day_fetal.csv", sep=",", header=T)
Braun_limit <- subset(Braun_pericyte_props, Days_fixed == 35)$X0

  # Set colors
colors <- c("hotpink", "#EF9AA8", "pink", "mediumorchid1", "#E0405B", "thistle2","chocolate1", "goldenrod1", "tan", "peachpuff")


  # Check how well each limit separates high- and low-pericyte donors
ggplot(data, aes(x = reorder(paste(pool, donor, sep="_"), -proportion_Per), y = proportion_Per)) +
  geom_bar(stat="identity", aes(fill=pool)) +
  scale_fill_manual(values = colors) +
  scale_x_discrete(labels = reorder(data$donor, -data$proportion_Per)) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, size = 5)) +
  xlab("Donor") +
  ylab("Pericyte proportion") +
  geom_abline(intercept = mean(data$proportion_Per), slope = 0, linewidth = 0.2) +
  geom_abline(intercept = Sanger_limit, slope = 0, linewidth = 0.2) +
  geom_abline(intercept = Braun_limit, slope = 0, linewidth = 0.2) +
  geom_abline(intercept = 0.1, slope = 0, linewidth = 0.2) +
  annotate("text", x = 80, y = Sanger_limit + 0.007, label = "Pool OT_G31139991327", size = 2) +
  annotate("text", x = 80, y = Braun_limit + 0.007, label = "Braun et al.", size = 2) +
  annotate("text", x = 80, y = mean(data$proportion_Per) + 0.007, label = "Mean proportion d20", size = 2) +
  annotate("text", x = 80, y = 0.107, label = "Selected limit", size = 2)


  # Check how the limits separate high- and low-pericyte donors withing each pool
ggplot(data, aes(x = paste(pool, donor, sep="_"), y = proportion_Per)) +
  geom_bar(stat="identity", aes(fill=pool)) +
  xlab("Donor") +
  scale_fill_manual(values = colors) +
  ylab("Pericyte proportion") +
  scale_x_discrete(labels = data$donor) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, size = 5)) +
  geom_abline(intercept = mean(data$proportion_Per), slope = 0, linewidth = 0.2) +
  geom_abline(intercept = Sanger_limit, slope = 0, linewidth = 0.2) +
  geom_abline(intercept = Braun_limit, slope = 0, linewidth = 0.2) +
  geom_abline(intercept = 0.1, slope = 0, linewidth = 0.2) +
  annotate("text", x = 80, y = Sanger_limit + 0.007, label = "Pool OT_G31139991327", size = 2) +
  annotate("text", x = 80, y = Braun_limit + 0.007, label = "Braun et al.", size = 2) +
  annotate("text", x = 80, y = mean(data$proportion_Per) + 0.007, label = "Mean proportion d20", size = 2) +
  annotate("text", x = 80, y = 0.107, label = "Selected limit", size = 2)
