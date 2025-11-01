# Load packages
library(Seurat)
library(dplyr)
library(ggplot2)
library(tidyr)


# Load original meta data extracted from the Seurat object
load("meta.RData")


# To ensure the final row names match the original ones
meta <- meta %>% mutate(orig_rownames = rownames(meta))


# Load marker-based cell types and join with original
marker_annot <- read.csv("desired_columns.csv",
                         sep = ",",
                         header = T)

meta <- left_join(meta, marker_annot)


# To switch to using another cell type, change this
meta$cell_type_annot <- meta$twoPassAnnotation_clean


### Before filtering


# Add column for pool, i.e. leave day out and unify the Sanger pools
meta <- meta %>%
  mutate(pool = sub("\\_d..", "", SampleID)) %>%
  mutate(pool = case_when(pool == "OT_G311310130744" |
                          pool == "OT_G31139991327" ~ "OT_PM_SangerPools",
                          .default = pool))


# Add the number of cells per donor, pool, and day
meta <- meta %>% group_by(SampleID, donor) %>% mutate(freq_donor = n()) %>% ungroup()

meta <- meta %>%
  group_by(pool) %>%
  mutate(first_day = min(Day_fixed),
         donors_per_pool = n_distinct(donor),
         proportion_d0 = 1/donors_per_pool) %>%
  ungroup() %>%
  group_by(SampleID) %>%
  mutate(cells_per_SampleID = sum(freq_donor)) %>%
  ungroup() %>%
  mutate(proportion_day_fixed = freq_donor/cells_per_SampleID) %>%
  mutate(R = proportion_day_fixed / proportion_d0)
  



## Adding donor meta

donor_meta <- read.table("hipsci.qc1_sample_info.20170927.tsv", sep="\t", header=T)



# Making donor meta easier to join with the original meta

colnames(donor_meta)[which(names(donor_meta) == "donor")] <- "donor_short"
colnames(donor_meta)[which(names(donor_meta) == "name")] <- "donor"

meta <- left_join(meta, donor_meta)


# Remove completely empty columns
meta <- meta[sapply(meta, function(x) sum(is.na(x)) != length(x))]


## Add gender/sex of donor from karyotypes when missing

# Check which donors are missing from donor_meta
#unique(meta$donor[is.na(match(meta$donor, donor_meta$donor))])

males <- c("HEL_312.3", "HEL_313.5", "HEL_314.1", "HEL_315.7", "HEL_318.3", "HEL_319.1",  "HEL_344.2",  "HPSI0516i-oadp_5_G6C3")
females <- c("HEL_316.5", "HEL_342.2", "HEL_317.1")

meta <- meta %>%
  mutate(gender = case_when(donor %in% males ~ "male",
                            donor %in% females ~ "female",
                            .default = gender))



## Proportion of patient lines per pool

# add the proportion of donors with a known disease in the first known day
# to compare pools based on the original setting

patient_prop_orig <- meta %>%
  filter(disease != "Normal" & !is.na(disease) & Day_fixed == first_day) %>%
  group_by(pool) %>%
  mutate(n_patient=n_distinct(donor)) %>%
  mutate(prop_patient_orig = n_patient/donors_per_pool, .keep="used") %>%
  distinct() %>%
  select(pool, prop_patient_orig)

meta <- left_join(meta, patient_prop_orig)


# Add the proportion of donors with a known disease for each day
# to see if the proportion of donors with a disease grows

patient_prop <- meta %>%
  filter(disease != "Normal" & !is.na(disease)) %>%
  group_by(SampleID) %>%
  mutate(n_patient = n_distinct(donor)) %>%
  mutate(prop_patient = n_patient/donors_per_pool, .keep="used") %>%
  distinct() %>%
  select(SampleID, prop_patient)

meta <- left_join(meta, patient_prop)


# Discretize the number of sendai reads
meta <- meta %>%
  mutate(sendai = case_when(
    rnaseq.sendai_reads == 0 ~ 0,
    rnaseq.sendai_reads > 0 ~ 1),
    sendai_approx = case_when(
      rnaseq.sendai_reads < 10 ~ 0,
      rnaseq.sendai_reads >= 10 ~ 1)
    )


#Defining top cell types

top_types <- meta %>%
  group_by(cell_type_annot) %>%
  summarize(n = n()) %>%
  arrange(desc(n)) %>%
  slice_head(n=6)


# Calculating cell type proportions for the top 6 types

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



joined_meta <- meta

row.names(joined_meta) <- joined_meta$orig_rownames

joined_meta <- joined_meta %>% select(-c(orig_rownames))



# Save the joined metadata

# Unfiltered
save(joined_meta, file = "joined_meta.RData")

# Filtered
save(joined_meta, file = "joined_meta_filtered.RData")



## Add meta to the pseudobulk objects


# Choose:
# 1) "", 2) "_filtered", 3) "_per", 4) "_others", or 5) "_filtered_per_nonper"
# for object pseudobulked:
# 1) the default way, unfiltered
# 2) the default way, filtered
# 3) the default way but contains only pericytes, filtered
# 4) the default way but contains only non-pericytes, filtered
# 5) by cell type category (pericyte/non-pericyte), filtered

pseudobulked <- ""
                    
joining <- function(data = "") {
              load(paste0("/scratch/project_2010414/Thesis_2.1/Thesis_2.1.2/Pseudobulk/pseudobulk_obj", data, ".RData"))
              
              pseudobulk_obj@meta.data <- pseudobulk_obj@meta.data %>%
                                            mutate(SampleID = gsub("g", "", SampleID))
                
              
              if(data != "") {
                joined_meta <- joined_meta %>% filter(freq_donor >= 50)
              }
  
              cols <- c("SampleID", "Day_fixed", "donorSample", 
                "pool", "geneMutatedOrCorrected", "proportion_d0",
                "proportion_day_fixed", "first_day", "donors_per_pool",
                "prop_patient_orig", "prop_patient", "R",
                "sendai", "sendai_approx", "freq_donor",
                "diff_method", "z_R", "z_pluri_novelty")
              
              donor_meta_joined_cols <- intersect(colnames(joined_meta), colnames(donor_meta))
  

              if(data == "" | data == "_filtered") {
                cols <- append(cols, c("proportion_oRG", "proportion_vRG",
                                       "proportion_panRG_O", "proportion_Per", 
                                       "proportion_PgS", "proportion_PgG2M",
                                       "z_proportion_Per", "z_proportion_vRG",
                                       "z_proportion_panRG_O", "z_proportion_PgS",
                                       "z_proportion_oRG", "z_proportion_PgG2M"))
              }
  
              
              pseudo_joined_meta <- joined_meta %>%
                                      select(all_of(donor_meta_joined_cols),
                                             all_of(cols)) %>%
                                      mutate(SampleID = gsub("_", "-", SampleID),
                                             donor = gsub("_", "-", donor)) %>%
                                      distinct()
              
              pseudobulk_obj@meta.data <- left_join(pseudobulk_obj@meta.data, pseudo_joined_meta,
                                                    by=join_by(SampleID, donor))
              
              return(pseudobulk_obj)
}


pseudobulk_obj <- joining(data = pseudobulked)

                    
# Add cell cycle scores

s.genes <- cc.genes$s.genes  
g2m.genes <- cc.genes$g2m.genes 

pseudobulk_counts <- pseudobulk_obj@assays$RNA$counts 

norm_pseudobulk_counts <- pseudobulk_counts %*% diag (mean(colSums(pseudobulk_counts)) / colSums(pseudobulk_counts))  

d_score <- cbind.data.frame(   
            orig.ident = colnames(pseudobulk_counts),  
            s = colSums(log10(norm_pseudobulk_counts+1)[s.genes[s.genes %in% rownames(norm_pseudobulk_counts)],]),  
            g2m = colSums(log10(norm_pseudobulk_counts+1)[g2m.genes[g2m.genes %in% rownames(norm_pseudobulk_counts)],]))

pseudobulk_obj@meta.data <- left_join(pseudobulk_obj@meta.data, d_score, by="orig.ident")

pseudobulk_obj$z_s <- (d_score$s - mean(d_score$s) ) / sd(d_score$s)  
pseudobulk_obj$z_g2m <- (d_score$g2m - mean(d_score$g2m) ) / sd(d_score$g2m)

row.names(pseudobulk_obj@meta.data) <- pseudobulk_obj$orig.ident


save(pseudobulk_obj, file = paste0("pseudobulk_obj_whole_meta", pseudobulked, ".RData")
