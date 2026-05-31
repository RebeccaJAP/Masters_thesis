  # Download packages
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(stringr)
library(gginnards)

  # Get the data
load("Pseudobulk_data/pseudobulk_obj_filtered.RData")

  # Set the colors
day_cols <- c("slateblue1", "hotpink",  "olivedrab2", "goldenrod1", "chocolate1", "black")

donor_cols = c("HPSI0316i-aask_4" = "hotpink",    
         "HPSI0316i-ierp_4" = "lightblue",     
         "HPSI0414i-seru_1" = "chocolate1",
         "HPSI0316i-qeti_2" = "slateblue1",
         "HPSI0516i-oadp_4" = "#E0405B",
         "HPSI1016i-riwg_2" = "peachpuff",
         "HPSI0314i-hoik_1" = "pink",
         "HPSI0516i-oadp_5_G6C3" = "goldenrod1",
         "HPSI0414i-oikd_2" = "thistle2",
         "HEL_315.7" = "mediumorchid1",
         "HPSI0114i-kolf_2" = "darkolivegreen1",
         "HPSI1016i-livl_2" = "#C0FAED",
         "HPSI0614i-paab_3" = "tan",
         "HPSI0214i-kucg_2" =  "#EF9AA8",
         "HEL_316.5" = "#1CC9BB",
         "HPSI0414i-oaqd_2" = "lightgoldenrod",
         "HPSI1016i-eoxu_2" = "olivedrab2",
         "HEL_342.2" = "limegreen",            
         "HPSI0914i-suop_5" = "#1600A3",
         "HPSI1213i-tolg_4" = "blue3",     
         "HPSI0514i-letw_1" = "orchid3",
         "HPSI0913i-lise_1" = "deepskyblue",     
         "HPSI0714i-kute_5" = "brown",
         "HEL_314.1" = "maroon2",
         "HPSI1113i-bima_1" = "lightgreen",     
         "HPSI0114i-eipl_1" = "forestgreen",     
         "HPSI0214i-rayr_1" = "#B57457",     
         "HPSI0714i-pipw_5" = "firebrick2",     
         "HPSI0115i-paim_3" = "orange",
         "HEL_313.5" = "darkorchid",            
         "HPSI1113i-podx_1" = "darkolivegreen3",
         "HEL_318.3" = "darkgoldenrod3",            
         "HEL_317.1" = "violet",            
         "HEL_319.1" = "deeppink",            
         "HEL_344.2" = "darkslateblue",
         "HEL_312.3" = "turquoise2")

pseudobulk_obj@meta.data <- pseudobulk_obj@meta.data %>%
  group_by(SampleID) %>%
  mutate(max_R = max(R)) %>%
  ungroup()


  # Define a function for removing NAs and getting sample sizes per day for plotting
  # Input names of two numerical variables (assumed to be continuous)
sample_sizer <- function(independent_var, dependent_var) {
  
    # Remove observations where at least one variable of interest is NA
  data <- pseudobulk_obj@meta.data[!is.na(pseudobulk_obj@meta.data[[independent_var]]) &
                                     !is.na(pseudobulk_obj@meta.data[[dependent_var]]), ]

    # Remove duplicate values when pool-level variables and R are of interest
  if(dependent_var == "max_R") {
    data <- data %>% distinct_at(c("SampleID", "Day_fixed", "max_R", independent_var), .keep_all = T)
  }
  
    # Get the sample size per time point and create a label used for plots
  data <- data %>%
    group_by(Day_fixed) %>%
    mutate(N = n(),
           label = paste0("Day ", Day_fixed, " (N = ", N, ")")) %>%
    ungroup()
  
  return(data)
}

  
# Adjust the function above for one or more categorical variables
    # independent_var needs to be categorical
    # dependent_var can be either categorical or continuous
sample_sizer_cat <- function(independent_var, dependent_var) {
  
    # Remove observations where at least one variable of interest is NA
  data <- pseudobulk_obj@meta.data[!is.na(pseudobulk_obj@meta.data[[independent_var]]) &
                                     !is.na(pseudobulk_obj@meta.data[[dependent_var]]), ]
  
    # Get the sample size per group and create a label used for plots
  data <- data %>%
    group_by_at(independent_var) %>%
    mutate(N = n(),
           label = paste0("N = ", N)) %>%
    ungroup()
  
  if(dependent_var == "max_R") {
    data <- data %>% distinct_at(c("SampleID", "Day_fixed", "max_R", independent_var), .keep_all = T)
  }
  
  return(data)
}


  # Define a function for plotting when the independent variable is continuous
plot_cont <- function(independent_var, dependent_var) {
  
  data_plot <- sample_sizer(independent_var, dependent_var)
  
  p <- ggplot(data_plot, aes(x = .data[[independent_var]], y = .data[[dependent_var]])) +
    geom_point(alpha = 0.70, size = 2.5,
               aes(colour = as.factor(Day_fixed),
                   shape = factor(diff_method, levels = c("pool", "PM", "indv")))) +
    scale_color_manual(values = day_cols) +
    scale_shape_manual(values = c("pool" = 16,
                                  "PM" = 15,
                                  "indv" = 17),
                       labels = c("Cell village", "Post-mitotic pool", "Individual line")) +
    labs(colour = "Day", shape = "Culturing condition") +
    theme(axis.text.x = element_text(size = 8),
          plot.title = element_text(size = 10),
          legend.title = element_text(size=9),
          axis.title = element_text(size=8))
  
  return(p)
  }


  # Define a function for plotting when the independent variable is categorical
    # categories: Input a vector of strings containing the category names to adjust x labels
    # color_by: Should the dots be colored by donor ("don"), time point ("time") or both ("both")?
plot_cat <- function(independent_var, dependent_var, categories = NULL, jitter_width = 0.35, color_by = "time") {
  
  data_plot <- sample_sizer_cat(independent_var, dependent_var)
  
    # Control jitter
  jitter_constant <- position_jitter(width = jitter_width, seed=34)
  
  
  levels <- unique(data_plot[[independent_var]])
  
  if (is.null(categories)) {categories <- levels}
  
    # Create labels
  
  label_levels <- vector(length = length(categories))
  
  for (i in 1:length(categories)) {
    label_levels[i] <- paste(#str_wrap(str_to_sentence(categories[i]), width = 20),
      str_wrap(categories[i], width = 20),
                             data_plot$label[!duplicated(data_plot$label) & data_plot[[independent_var]] == levels[i]], sep = "\n")
  }
  
  p <- ggplot(data_plot, aes(x=factor(.data[[independent_var]], levels = levels),
                             y=.data[[dependent_var]],
                             fill=factor(.data[[independent_var]], levels = levels))) +
    geom_boxplot(fill = "white", outlier.shape = NA) +
    geom_jitter(position = jitter_constant, size = 2, aes(colour = as.factor(donor)), alpha = 0.75) +
    scale_color_manual(values = donor_cols) +
    scale_x_discrete(labels = label_levels) +
    guides(alpha = "none", size = "none", fill = "none") +
    labs(x=NULL, colour = "Donor") +
    theme(legend.title = element_text(size = 9),
          plot.title = element_text(size = 10),
          axis.title = element_text(size = 8),
          axis.text = element_text(size = 8),
          legend.position = "right",
          legend.key.size = unit(0.4, 'cm'),
          legend.text = element_text(size = 8),
          legend.key = element_blank(),
          legend.margin = margin(20, 1, 1, 1)
         )
  
  pp <- ggplot(data_plot, aes(x=factor(.data[[independent_var]], levels = levels),
                              y=.data[[dependent_var]],
                              fill=factor(.data[[independent_var]], levels = levels))) +
    geom_boxplot(fill = "white", outlier.shape = NA) +
    geom_jitter(position = jitter_constant, size = 2, aes(colour = as.factor(Day_fixed)), alpha = 0.75) +
    scale_color_manual(values = day_cols) +
    scale_x_discrete(labels = label_levels) +
    guides(alpha = "none", size = "none", fill = "none", colour = guide_legend(ncol = 1)) +
    labs(x = NULL, colour = "Day") +
    theme(legend.title = element_text(size = 9),
          plot.title = element_text(size = 10),
          axis.title = element_text(size = 8),
          axis.text = element_text(size = 8),
          legend.position = "right",
          legend.key.size = unit(0.4, 'cm'),
          legend.text = element_text(size = 8),
          legend.key = element_blank(),
          legend.margin = margin(1, 1, 1, 1))
  
  combined <- p | pp &
    plot_annotation(tag_levels = 'A') &
    theme(plot.tag = element_text(face = 'bold'),
          plot.tag.position = c(0.01, 1.01))
  
  if (color_by == "both") {return(combined)
  } else if (color_by == "don") {
    return(p)
  } else if (color_by == "time") {
    return(pp)
  }
  
}



  # Function for computing the correlation with continuous independent variables

cor_cont <- function(independent_var, dependent_var, cor_method = "kendall") {
  data_cor <- sample_sizer(independent_var, dependent_var)
  
    # Correlation with all datapoints included
  print(paste0("All days; N = ", nrow(data_cor)))
  print(cor.test(data_cor[[independent_var]], data_cor[[dependent_var]], method = cor_method))
  
    # Correlation per day
  for (day in c(20,40,60,70,80)) {
    data_day <- data_cor %>% filter(Day_fixed == day)
    print(paste0("Day: ", day, "; N = ", unique(data_day$N)))
    print(cor.test(data_day[[independent_var]], data_day[[dependent_var]], method = cor_method))
  }
}



cor_cat <- function(independent_var, dependent_var, eq_var = F) {
  data_cor <- sample_sizer_cat(independent_var, dependent_var)
  print(oneway.test(data_cor[[dependent_var]] ~ as.factor(data_cor[[independent_var]]), var.equal = eq_var))
}


########

  # Visualization

  # Create the plots for R number and each independent variable

g2m_R_plot <- plot_cont("z_g2m", "R") +
  xlab("Z-score of the G2M score") +
  ylab("R number") +
  ggtitle("G2M score")

s_R_plot <- plot_cont("z_s", "R") +
  xlab("Z-score of the S score") +
  ylab("R number") +
  ggtitle("S score")

novelty_R_plot <- plot_cont("pluri_novelty", "R") +
  xlab("Pluripotency novelty score") +
  ylab("R number") +
  ggtitle("Pluripotency novelty")

  
poolsize_R_plot <- plot_cont("donors_per_pool", "max_R") +
  xlab("Number of donors in pool") +
  ylab("R number") +
  ggtitle("Pool size") +
  xlim(0,18)


patient_prop_R_plot <- plot_cont("prop_patient", "max_R") +
  xlab("Proportion of patient lines") +
  ylab("R number") +
  ggtitle("Proportion of patient lines in pool") +
  xlim(0,1)

disease_R <- plot_cat("disease", "R", color_by = "time",
                      categories = c("no NDD", "Kabuki", "RND")) +
  labs(y = "R number") +
  guides(color = "none") +
  ggtitle("Donor disorder status")

diff_method_R <- plot_cat("diff_method", "R",
                          color_by = "time",
                          categories = c("cell village", "individual lines", "PM")) +
  labs(y = "R number") +
  guides(color="none") +
  ggtitle("Differentiation method")

diff_method_maxR <- plot_cat("diff_method", "max_R",
                          color_by = "time",
                          categories = c("cell village", "individual lines", "PM")) +
  labs(y = "R number") +
  guides(color="none") +
  ggtitle("Differentiation method")

sex_R <- plot_cat("gender", "R", color_by = "time") +
  labs(y = "R number") +
  guides(color="none") +
  ggtitle("Donor sex")

sex_R_donor <- plot_cat("gender", "R", color_by = "don") +
  labs(y = "R number") +
  ggtitle("Donor sex")

sendai_R <- plot_cat("sendai", "R", color_by = "time",
                     categories = c("Sendai negative", "Sendai positive")) +
  labs(y = "R number") +
  guides(color="none") +
  ggtitle("Sendai positivity")

sendai_R_donor <- delete_layers(sendai_R, "GeomPoint") + 
  geom_jitter(color = ifelse(sendai_R$data$donor == "HPSI0913i-lise_1", "slateblue1", "grey"),
              alpha = ifelse(sendai_R$data$donor == "HPSI0913i-lise_1",1, 0.5),
              position = position_jitter(width = 0.35, seed=34),
              size = 2) +
  ggtitle("Sendai positivity of donor \nHPSI0913i-lise_1")


R_cor_plots <- (free((g2m_R_plot / s_R_plot / poolsize_R_plot / patient_prop_R_plot), type = "label") |
    ((free(novelty_R_plot, type = "label") | sex_R) /
       (diff_method_R | disease_R) /
       (sendai_R | sendai_R_donor))) +
  plot_layout(guides = "collect", widths = c(1,2)) +
  plot_annotation(tag_levels = 'A') &
  theme(plot.tag = element_text(face = 'bold'),
        plot.tag.position = c(0.01, 0.98),
      #  legend.position = "bottom") &
      legend.position = "right") &
  scale_y_continuous(breaks = c(0, 4, 8, 12))

pdf(file=paste0("R_correlation_plots.pdf"), width = 10, height = 10)
plot(R_cor_plots)
dev.off()


    ## Create plots for pericyte proportion and each independent variable

  # Visualize the relationship between cell line imbalance and the pericyte proportion
R_peri_plot <- plot_cont("R", "proportion_Per") +
  ylab("Pericyte proportion") +
  xlab("R number") +
  ggtitle("R number")

g2m_peri_plot <- plot_cont("z_g2m", "proportion_Per") +
  xlab("Z-score of the G2M score") +
  ylab("Pericyte proportion") +
  ggtitle("G2M score")

s_peri_plot <- plot_cont("z_s", "proportion_Per") +
  xlab("Z-score of the S score") +
  ylab("Pericyte proportion") +
  ggtitle("S score")

novelty_peri_plot <- plot_cont("pluri_novelty", "proportion_Per") +
  xlab("Pluripotency novelty score") +
  ylab("Pericyte proportion") +
  ggtitle("Pluripotency novelty score")

poolsize_peri_plot <- plot_cont("donors_per_pool", "proportion_Per") +
  xlab("Number of donors in pool") +
  ylab("Pericyte proportion") +
  ggtitle("Pool size") +
  xlim(0,18)

patient_prop_peri_plot <- plot_cont("prop_patient", "proportion_Per") +
  xlab("Proportion of patient lines") +
  ylab("Pericyte proportion") +
  ggtitle("Proportion of patient lines") +
  xlim(0,1)

sex_peri <- plot_cat("gender", "proportion_Per") +
  ggtitle("Donor sex") +
  guides(color = "none") +
  labs(y = "Pericyte proportion")

sendai_peri <- plot_cat("sendai", "proportion_Per",
                       categories = c("Sendai negative", "Sendai positive")) +
  ggtitle("Sendai positivity") +
  guides(color = "none") +
  labs(y = "Pericyte proportion")

sendai_peri_donor <- delete_layers(sendai_peri, "GeomPoint") + 
  geom_jitter(color = ifelse(sendai_peri$data$donor == "HPSI0913i-lise_1", "slateblue1", "grey"),
              alpha = ifelse(sendai_peri$data$donor == "HPSI0913i-lise_1",1, 0.5),
              position = position_jitter(width = 0.35, seed=34),
              size = 2) +
  guides(color = "none") +
  ggtitle("Sendai positivity of donor \nHPSI0913i-lise_1")


disease_peri <- plot_cat("disease", "proportion_Per",
                          categories = c("no NDD", "Kabuki", "RND")) +
  ggtitle("Donor disorder status") +
  guides(color = "none") +
  labs(y = "Pericyte proportion")
  

diff_method_peri <- plot_cat("diff_method", "proportion_Per",
                              categories = c("cell village", "indvidual", "PM")) +
  ggtitle("Differentiation method") +
  guides(color = "none") +
  labs(y = "Pericyte proportion")


peri_cor_plots <- (free((R_peri_plot / g2m_peri_plot / s_peri_plot / poolsize_peri_plot / patient_prop_peri_plot), type = "label") |
    ((free(novelty_peri_plot, type = "label") | sex_peri) /
       (diff_method_peri | disease_peri) /
       (sendai_peri | sendai_peri_donor))) +
  plot_layout(guides = "collect", widths = c(1,2)) +
  plot_annotation(tag_levels = 'A') &
  theme(plot.tag = element_text(face = 'bold'),
        plot.tag.position = c(0.01, 0.98),
        legend.position = "right")

pdf(file=paste0("pericyte_correlation_plots.pdf"), width = 10, height = 10)
plot(peri_cor_plots)
dev.off()

########

  # Correlation computations for continuous variables (R)

cor_cont("z_g2m", "R")

cor_cont("z_s", "R")

cor_cont("donors_per_pool", "max_R")

cor_cont("prop_patient", "max_R")

cor_cont("pluri_novelty", "R")


  # Test if means differ between different levels of categorical independent variables (R)

cor_cat("disease", "R")

cor_cat("diff_method", "R")
cor_cat("diff_method", "max_R")

cor_cat("gender", "R")

cor_cat("sendai", "R")
cor_cat("sendai_approx", "R")


  # Correlation computations for continuous variables (proportion of pericytes)

cor_cont("R", "proportion_Per")
cor_cont("z_g2m", "proportion_Per")
cor_cont("z_s", "proportion_Per")

cor_cont("pluri_novelty", "proportion_Per")
cor_cont("donors_per_pool", "proportion_Per")
cor_cont("prop_patient", "proportion_Per")


  # Test if means differ between different levels of categorical independent variables (pericyte proportion)

cor_cat("gender", "proportion_Per")
cor_cat("sendai", "proportion_Per")
cor_cat("sendai_approx", "proportion_Per")
cor_cat("disease", "proportion_Per")
