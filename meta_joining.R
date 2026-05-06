  # Load packages
library(Seurat)
library(dplyr)
library(ggplot2)
library(tidyr)


  # Load original meta data
load("meta.RData")

  # To ensure the final row names match the original ones
meta <- meta %>% mutate(orig_rownames = rownames(meta))


  # Load marker-based cell types and join with original
marker_annot <- read.csv("desired_columns.csv",
                         sep = ",",
                         header = T)

meta <- left_join(meta, marker_annot)


  # Choose the main annotation type used (both are kept regardless)
meta$cell_type_annot <- meta$twoPassAnnotation_clean
# meta$cell_type_annot <- meta$cell_type_pred_y


  # Add column for pool, i.e. leave day out and unify the Sanger pools' names
meta <- meta %>%
  mutate(pool = sub("\\_d..", "", SampleID)) %>%
  mutate(pool = case_when(pool == "OT_G311310130744" |
                            pool == "OT_G31139991327" ~ "OT_PM_SangerPools",
                          .default = pool))


  # Add the number of cells per donor, pool, and day
meta <- meta %>% group_by(SampleID, donor) %>% mutate(freq_donor = n()) %>% ungroup()
  
  # Get the first day, original donor count per pool & proportion of each donor's cells on the first day
meta <- meta %>%
  group_by(pool) %>%
  mutate(first_day = min(Day_fixed),
         donors_per_pool = n_distinct(donor),
         proportion_d0 = 1/donors_per_pool) %>%
  ungroup()

  # Compute the donor proportions and R for each time point
donor_props <- meta %>%
  distinct(SampleID, donor, Day_fixed, donors_per_pool, proportion_d0, freq_donor) %>%
  group_by(SampleID) %>%
  mutate(cells_per_SampleID = sum(freq_donor)) %>%
  ungroup() %>%
  mutate(proportion_day_fixed = freq_donor/cells_per_SampleID) %>%
  mutate(R = proportion_day_fixed / proportion_d0)

meta <- left_join(meta, donor_props)


    ## Adding the donor metadata
donor_meta <- read.table("hipsci.qc1_sample_info.20170927.tsv", sep="\t", header=T)

  # Make the donor meta easier to join with the original meta
colnames(donor_meta)[which(names(donor_meta) == "donor")] <- "donor_short"
colnames(donor_meta)[which(names(donor_meta) == "name")] <- "donor"

meta <- left_join(meta, donor_meta)


  # Remove completely empty columns
meta <- meta[sapply(meta, function(x) sum(is.na(x)) != length(x))]


## Add the sex of a donor from karyotype results when missing

# Check which donors are missing from donor_meta
unique(meta$donor[is.na(match(meta$donor, donor_meta$donor))])

  # Add sex information manually from karyotyping results
males <- c("HEL_312.3", "HEL_313.5", "HEL_314.1", "HEL_315.7", "HEL_318.3", "HEL_319.1",  "HEL_344.2",  "HPSI0516i-oadp_5_G6C3")
females <- c("HEL_316.5", "HEL_342.2", "HEL_317.1")

meta <- meta %>%
  mutate(gender = case_when(donor %in% males ~ "male",
                            donor %in% females ~ "female",
                            .default = gender))


    ## Proportion of patient lines per pool

# Add the proportion of donors with a known disease in the first known day to compare pools based on the original setting

patient_prop_orig <- meta %>%
  filter(disease != "Normal" & !is.na(disease) & Day_fixed == first_day) %>%
  group_by(pool) %>%
  mutate(n_patient=n_distinct(donor)) %>%
  mutate(prop_patient_orig = n_patient/donors_per_pool, .keep="used") %>%
  distinct() %>%
  select(pool, prop_patient_orig)

meta <- left_join(meta, patient_prop_orig)


  # Add the proportion of donors with a known disease for each day to see if the proportion of donors with a disease grows

patient_prop <- meta %>%
  filter(disease != "Normal" & !is.na(disease)) %>%
  group_by(SampleID) %>%
  mutate(n_patient = n_distinct(donor)) %>%
  mutate(prop_patient = n_patient/donors_per_pool, .keep="used") %>%
  distinct() %>%
  select(SampleID, prop_patient)

meta <- left_join(meta, patient_prop)


  # Discretize the number of sendai positive reads
meta <- meta %>%
  mutate(sendai = case_when(
    rnaseq.sendai_reads == 0 ~ 0,
    rnaseq.sendai_reads > 0 ~ 1),
    sendai_approx = case_when(
      rnaseq.sendai_reads < 10 ~ 0,
      rnaseq.sendai_reads >= 10 ~ 1)
  )


  # Define the top cell types (using the main annotation type)

  # Choose the number of top cell types to use based on an elbow plot
barplot(sort(table(meta$cell_type_annot)), las = 2)

  # The top 6 cell types are distinct from the others in their frequency
top_types <- meta %>%
  group_by(cell_type_annot) %>%
  summarize(n = n()) %>%
  arrange(desc(n)) %>%
  slice_head(n=6)

  # Calculate the cell type proportions for the top 6 types
cell_type_proportion <- meta %>%
  group_by(SampleID, donor, cell_type_annot) %>%
  dplyr::count(cell_type_annot, .drop=F) %>%
  ungroup(cell_type_annot) %>%
  mutate(cell_type_prop = n/sum(n), .keep = "unused") %>%
  pivot_wider(names_from = cell_type_annot,
              values_from = cell_type_prop,
              values_fill = 0,
              names_prefix = "proportion_") %>%
  select(SampleID, donor, contains(top_types$cell_type_annot))

colnames(cell_type_proportion) <- sub("-", "_", colnames(cell_type_proportion))

meta <- left_join(meta, cell_type_proportion)


  # Add the differentiation method (pooled, pooled PM, differentiated individually)
meta <- meta %>%
  mutate(diff_method = case_when(grepl("indv", pool, fixed = TRUE) ~ "indv",
                                 grepl("PM", pool, fixed = TRUE) ~ "PM",
                                 .default = "pool"))


    ## Donor meta + original meta confounders

  # Calculate confounder Z scores within pool and day
top_type_cols <- paste0("proportion_", sub("-", "_", top_types$cell_type_annot))

z_scores <- meta %>%
  group_by(SampleID) %>%
  distinct(SampleID, donor, across(all_of(top_type_cols)), R, pluri_novelty) %>%
  mutate(across(c(all_of(top_type_cols), R, pluri_novelty), .fns = c(sd = sd, mean = mean))) %>%
  ungroup() %>%
  mutate(across(ends_with("_mean"), 
                .fns = ~ (get(sub("_mean$", "", cur_column())) - .) /
                  get(sub("_mean$", "_sd", cur_column())),
                .names = "z_{col}")) %>%
  rename_with(~ sub("_mean", "", .), starts_with("z_"))

z_scores <- z_scores %>% select(-contains(c("_mean", "_sd")))

meta <- left_join(meta, z_scores)

  # Change the object name
joined_meta <- meta

  # Add differently formatted columns for SampleID & donor to enable joining with pseudobulk data
  # Add a column for pericyte/non-pericyte category
joined_meta <- joined_meta %>%
  mutate(new_SampleID =  gsub("_", "-", SampleID),
         new_donor = gsub("_", "-", donor),
         cell_type_group = if_else(twoPassAnnotation_clean == "Per", "pericyte", "nonPericyte"))

  # Set the rownames to the original ones
row.names(joined_meta) <- joined_meta$orig_rownames

joined_meta <- joined_meta %>% select(-c(orig_rownames))



	# Save the joined metadata
save(joined_meta, file = "joined_meta.RData")
         
         
  # Filter the data to only include donors with at least 50 cells
joined_meta <- joined_meta %>%
  filter(freq_donor >= 50)
         
  # Save the filtered joined metadata
save(joined_meta, file = "joined_meta_filtered.RData")
         