  # Load packages
library(dplyr)
library(ggplot2)
library(stringr)

  # Choose donor and chromosome of interest
don <- "riwg-2"; chr_of_interest <- "chr1"
#don <- "riwg-2"; chr_of_interest <- "chr15"
#don <- "HEL-318.3"; chr_of_interest <- "chr1"


  # Get the paths for relevant InferCNV and SCEVAN result folders:
inf_files <- list.files(path = "CNV_calling/InferCNV/Results",
                        pattern = paste0("*", don, "*"),
                        full.names = TRUE,
                        recursive = FALSE)

sce_files <- list.files(path = "CNV_calling/SCEVAN/output",
                        pattern=paste0("*", don, "*"),
                        full.names=TRUE,
                        recursive=FALSE)


  # Prepare a dataframe for the barplot data
barplot_data <- data.frame(CNV = numeric(),
                           chr = character(),
                           interval_start = numeric(),
                           interval_end = numeric(),
                           cell_type = character(),
                           donor = character(),
                           day = numeric(),
                           pool = character(),
                           annotation = character(),
                           method = character()
)


  # Add the InferCNV results to the barplot data
for (path in inf_files) {
  
    # Load data where cells are grouped by copy number state, the affected chromosome, and the location within
  data <- read.delim(paste0(path, "/HMM_CNV_predictions.HMMi6.leiden.hmm_mode-subclusters.Pnorm_0.5.pred_cnv_regions.dat"), header=T, sep="\t")
  
    # Add version, annotation type, pool, timepoint, donor, and cell type information
  dataname <- sapply(strsplit(path, "/"), tail, 1)
  dataname_split <- unlist(strsplit(dataname, "_"))
  version <- dataname_split[1]
  annotation <- dataname_split[2]
  pool <- dataname_split[3]
  day <- dataname_split[4]
  day <- as.numeric(sub("d", "", day))
  donor <- dataname_split[5]
  cell_type <- dataname_split[6]
  cell_type <- ifelse(annotation == "marker", gsub("-", " ", cell_type), cell_type)
  
  annot_title <- ifelse(annotation == "marker", "Marker-based\n annotation", "Two-pass annotation")
  
  data$cell_type <- cell_type
  data$donor <- donor
  data$day <- day
  data$pool <- pool
  data$annotation <- annotation
  data$method <- "InferCNV"
  data$version <- version
  data$annot_title <- annot_title
  
    # Remove reference cells from the data
  data <- data %>%
    filter(!grepl("reference",cell_group_name))
  
    # Change from the states to copy numbers
  data$state[data$state==1] <- 0
  data$state[data$state==2] <- 1
  data$state[data$state==3] <- 2
  data$state[data$state==4] <- 3
  data$state[data$state==5] <- 4
  data$state[data$state==6] <- 5
  colnames(data)[3] <- "CNV"
  
    # Load data with the group assignment of each cell
  obs_groups <- read.delim(paste0(path, "/infercnv.observation_groupings.txt"), sep=" ", header=T)
  
    # Compute the group sizes, i.e. how many cells per group
  obs_groups <- obs_groups %>%
    group_by(Dendrogram.Group) %>%
    count()
  
    # Change cell group names to match the data and add to the data
  obs_groups$cell_group_name <- paste0("all_observations.", obs_groups$Dendrogram.Group)
  
  obs_groups <- obs_groups %>%
    arrange(cell_group_name)
  
  data <- merge(data, obs_groups, by="cell_group_name", all=T)
  
    # Keep only the data for the chromosome of interest
  data <- data %>% filter(chr==chr_of_interest)
  
    # Create a sorted list with all sites where any CNV starts or ends
  cut_coordinates <- unique(c(1,data$start, data$end))
  cut_coordinates <- cut_coordinates[!is.na(cut_coordinates)]
  cut_coordinates <- sort(cut_coordinates)
  
  
  interval_data <- data.frame(matrix(ncol=3, nrow=0))
  colnames(interval_data) <- c("interval_start", "interval_end", "CNV")
  
  for (i in (seq(1,length(cut_coordinates)-1))) {
    groups <- data %>%
        # Get all CNVs in the interval between the current and the following coordinate
      filter((cut_coordinates[i] >= start & cut_coordinates[i] < end) | (cut_coordinates[i+1] > start & cut_coordinates[i+1] <= end)) %>%
        # Compute the mean of the CNVs of the interval weighed by the CNV size
      mutate(interval_start = cut_coordinates[i], interval_end = cut_coordinates[i+1],
             CNV = weighted.mean(CNV, n)) %>%
      distinct(interval_start, interval_end, CNV, chr, cell_type, donor, day, pool, annotation, method, version)
    
    
    interval_data <- bind_rows(interval_data,unique(groups))
  }
  
  interval_data <- interval_data %>% arrange(interval_start)
  
  
  barplot_data <- rbind(barplot_data, interval_data)
}

  # Only include InferCNV results from one chromosome at a time
  # and only include cell types with results available in all pools considered
inf_barplot_data <- barplot_data %>%
  filter(method == "InferCNV") %>%
  group_by(cell_type) %>%
  mutate(n = n_distinct(pool)) %>%
  ungroup() %>%
  mutate(max_n = max(n)) %>%
  filter(n == max_n)

  # Set the bar positions based on the number of cell types included in the plot
y_coords <- inf_barplot_data %>%
  group_by(annotation) %>%
  distinct(cell_type) %>%
  mutate(ymin = row_number()-0.3,
         ymax = row_number()+0.3)

inf_barplot_data <- left_join(inf_barplot_data, y_coords)

  # Plot the barplot for InferCNV results
ggplot(inf_barplot_data, aes(y = cell_type, fill = CNV)) + 
  geom_rect(aes(xmin=interval_start, xmax=interval_end, ymin = ymin, ymax = ymax, fill = CNV),
            size=1,
            stat="identity") +
  guides(alpha=F) +
  scale_fill_gradientn(
    colours = c("#1600A3", "slateblue1", "white", "hotpink", "#D10069"),
    limits = c(0,4)) +
  scale_x_continuous(
    name = "Genomic position (Mb)",
    labels = function(x) x / 1e6) +
  theme(plot.title = element_text(size = 12),
        legend.background = element_rect("grey90"),
        strip.text = element_text(size=8)) +
  ggtitle(paste0("InferCNV results for ", donor, " in chromosome ", sub("chr", "", chr_of_interest))) +
  ylab("Cell type") +
  facet_grid(vars(pool, method), rows=vars(annotation), space="free_y", scales="free_y")



    ## Add SCE results


for (path in sce_files) {

  if (file.exists(path)) {
      # Get the data
    data <- read.table(path, header=T, sep = "\t", stringsAsFactors = FALSE)

      # Select the necessary columns and make the column names cohesive with the previous ones
    data <- data %>%
      select(Chr, Pos, End, CN, segm.mean) %>%
      mutate(Chr = paste0("chr", Chr))
    
    colnames(data) <- c("chr", "start", "end", "CNV", "segment_mean")

      # Get the relevant information from the file name
    dataname <- sapply(strsplit(path, "/"), tail, 1)
    dataname_split <- unlist(strsplit(dataname, "_"))
    version <- dataname_split[1]
    annotation <- dataname_split[2]
    pool <- dataname_split[3]
    day <- dataname_split[4]
    day <- as.numeric(sub("d", "", day))
    donor <- dataname_split[5]
    cell_type <- dataname_split[6]
    cell_type <- ifelse(annotation == "marker", gsub("-", " ", cell_type), cell_type)

      # Needed for the plot text to ensure clarity
    annot_title <- ifelse(annotation == "marker", "Marker-based\n annotation", "Two-pass annotation")

      # Add the information to the plotting data
    data$cell_type <- cell_type
    data$donor <- donor
    data$day <- day
    data$pool <- pool
    data$annotation <- annotation
    data$method <- "SCEVAN"
    data$version <- version
    data$annot_title <- annot_title
    
      # Add to the data used to plot
    barplot_data <- bind_rows(barplot_data, data)
  }
  else {
    print(dataname)
    next
  }
}

    # Only include SCEVAN results from one chromosome at a time
    # and only include cell types with results available in all pools considered
sce_barplot_data <- barplot_data %>%
  filter(chr==chr_of_interest, method == "SCEVAN") %>%
  group_by(cell_type) %>%
  mutate(n = n_distinct(pool)) %>%
  ungroup() %>%
  mutate(max_n = max(n)) %>%
  filter(n == max_n)

    # Plot the barplot for SCEVAN
  ggplot(sce_barplot_data) + 
    geom_segment(aes(x=start, xend=end, y=cell_type, color = as.factor(CNV)), size=10) +
    guides(alpha=F) +
    scale_color_manual(values = c("0" = "#1600A3", "1"="slateblue1", "2"="white", "3"="hotpink", "4"="#D10069")) +
    scale_x_continuous(
      name = "Genomic position (Mb)",
      labels = function(x) x / 1e6) +
    theme(plot.title = element_text(size = 12),
          legend.background = element_rect("grey90"),
          strip.text = element_text(size=8)) +
    ggtitle(paste0("SCEVAN results for ", donor, " in chromosome ", sub("chr", "", chr_of_interest))) +
    ylab("Cell type") +
    labs(color = "CNV") +
    facet_grid(vars(pool, method), rows=vars(annot_title), space="free_y", scales="free_y")
  
