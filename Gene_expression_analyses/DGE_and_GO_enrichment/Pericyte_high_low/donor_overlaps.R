  # Load packages
library(gt)
library(gtExtras)
library(purrr)
library(dplyr)
library(ggplot2)

  # Load the filtered pseudobulk object
load("Pseudobulk_data/pseudobulk_obj_filtered.RData")

  # Categorize the cell lines as either high- or low-pericyte
pseudobulk_obj@meta.data <- pseudobulk_obj@meta.data %>%
  mutate(per_category = case_when(Day_fixed == 20 & proportion_Per > 0.1 ~ "high",
                                  Day_fixed == 20 & proportion_Per <= 0.1 ~ "low"))


  # Get samples with at least one low-pericyte and one high-pericyte cell line
  # and filter out individual lines and later time points
data <- pseudobulk_obj@meta.data %>%
  filter(diff_method != "indv" & Day_fixed == 20) %>%
  group_by(SampleID) %>%
  mutate(n=n_distinct(per_category)) %>%
  ungroup() %>%
  filter(n > 1)

  # Get the SampleIDs of the pools of interest
valid_samples <- unique(data$SampleID)

  # For each donor, get the total number of cell lines, the number of pericyte-high cell lines,
  #and the number of pericyte-low cell lines
check_overlap <- data %>%
  filter(SampleID %in% valid_samples) %>%
  group_by(donor) %>%
  mutate(per_category = as.factor(per_category),
         samples_donor = n_distinct(SampleID),
         categories_donor = n_distinct(per_category)) %>%
  filter(samples_donor > 1) %>%
  group_by(donor, samples_donor, per_category, .drop=F) %>%
  summarise(categorywise_count_donor = n_distinct(SampleID)) %>%
  arrange(desc(samples_donor), donor) %>%
  dplyr::select(donor, samples_donor, per_category,
                categorywise_count_donor) %>%
  distinct()

  # Change the pericyte category to factors to ensure the order is kept
plot_data <- check_overlap %>% mutate(per_category = factor(per_category, levels = c("low", "high")))

  # Compute the category distribution per donor and the plot label placements
plot_data <- plot_data %>%
  arrange(donor, per_category) %>%
  group_by(donor) %>%
  mutate(
    proportion = categorywise_count_donor / samples_donor,
    cumulative = cumsum(proportion),
    midpoint = cumulative - proportion / 2) %>%
  ungroup()

  # Set colors
colors <- c("pink", "lightblue")

  # Create a function for plotting a barplot for each donor
plot_stacked_bar_donor <- function(donor_i) {
  plot_data %>%
    filter(donor == donor_i) %>%
    ggplot(aes(y = donor, x = proportion, fill = per_category)) +
    geom_bar(stat = "identity",
             width = 0.6,
             position = position_stack(reverse = TRUE)) +
    scale_fill_manual(values = c(colors[2], colors[1])) +
    geom_text(aes(x = midpoint,
                  label = categorywise_count_donor),
              color = "black",
              size = 40) +
    scale_x_continuous(labels = scales::percent) +
    labs(y = "Donor",
         x = "Proportion",
         fill = "Category") +
    guides(fill="none") +
    theme_void() 
}


  # Create the table column label
labels = c("Low-pericyte", "high-pericyte")

label_built <- if (length(labels) == 2) {
  lab_pal1 <- colors2[2]
  lab_pal2 <- colors2[1]
  lab1 <- labels[1]
  lab2 <- labels[2]
  glue::glue(
    "<span style='color:{lab_pal1}'><b>{lab1}</b></span>",
    " vs. ",
    "<span style='color:{lab_pal2}'><b>{lab2}</b></span>",
    " cell lines"
  ) %>%
    gt::html()
}

  # Create the final table
final_table <- check_overlap %>%
  select(-c(per_category, categorywise_count_donor)) %>%
  mutate(Distribution = donor) %>%
  distinct() %>%
  ungroup() %>%
  gt() %>%
  cols_label(donor = "Donor", samples_donor = "Number of donor cell lines", Distribution = gt::html(label_built)) %>%
  cols_align(align = c("center"),
             columns = c(-donor)) %>%
  tab_options(data_row.padding = px(0)) %>%
  tab_style(style=cell_text(weight = "bold"), locations = cells_column_labels(columns = everything())) %>%
  text_transform(
    locations = cells_body(columns = 'Distribution'),
    fn = function(column) {
      map(column, plot_stacked_bar_donor) |>
        ggplot_image(height = px(30), aspect_ratio = 3)
    })

final_table

final_table %>% gtsave("Donor_overlaps.html")


