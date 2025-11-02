library(tidyr)
library(dplyr)
library(ggplot2)

load("joined_meta.RData")

# Join day 0 data with the rest of the days
d0_meta <- joined_meta %>%
              distinct(pool, donor, proportion_d0, diff_method) %>%
              mutate(SampleID = paste0(pool, "_d0"),
                     Day_fixed = 0,
                     proportion_day_fixed = proportion_d0) %>%
              select(-proportion_d0)

joined_with_d0 <- bind_rows(joined_meta, d0_meta)

# Create all combinations of donor and Day_fixed per pool
prop_data_final <- joined_with_d0 %>%
                      distinct(SampleID, pool, donor, Day_fixed, proportion_day_fixed, diff_method) %>%
                      group_by(pool, diff_method) %>%
                      complete(donor, Day_fixed) %>%
                      mutate(proportion_day_fixed = case_when(is.na(proportion_day_fixed) ~ 0,
                                                              .default = proportion_day_fixed))

# Designate colors
cols = c("HPSI0316i-aask_4" = "hotpink",    
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


donor_imbalance_plot <- ggplot(prop_data_final,
                                    aes(x=Day_fixed, y=proportion_day_fixed, fill=donor)) + 
                        geom_area(color="white", linewidth=0.2) +
                        theme(axis.title = element_text(size = 9),
                              strip.text = element_text(size = 8,
                                  margin = margin(0.02,0.05,0.02,0.05, "cm")),
                              legend.key.size = unit(0.4, 'cm'),
                              legend.title = element_text(size = 9),
                              legend.text = element_text(size = 8),
                              legend.margin = margin(0.01, 0.01, 0.01, 0.01),
                              legend.position = "bottom"
                              ) +
                        guides(fill = guide_legend(title = "Donor")) +
                        scale_fill_manual(values = cols) + 
                        labs(x = "Day", y="Proportion of cells by donor") +
                        facet_wrap(~ pool)


# Save
pdf(donor_imbalance_plot, file = "donor_imbalance.pdf", paper = "a4")
