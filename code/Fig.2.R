pacman::p_load(dplyr, RColorBrewer, reshape, ggpubr, tidyverse, pheatmap, GenomicRanges, ggVennDiagram)

# Create output directory
output_dir <- "~/Desktop/mosaic_work"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

reads <- read.table("~/Google Drive/My Drive/edinburgh/vsg/cobalt_2025/cdhit-94-94/summary/read_counts/reads.txt", check.names = FALSE, header = TRUE)
tpm <- read.csv("~/Google Drive/My Drive/edinburgh/vsg/cobalt_2025/cdhit-94-94/summary/tpm/tpm.csv", check.names = FALSE, header = TRUE)
meta <- read.csv("~/Google Drive/My Drive/edinburgh/vsg/vsgseq2/data/meta.csv", check.names = FALSE, header = TRUE)

tpm <- tpm %>%
  mutate(Names = if_any(where(is.numeric), ~ . > 100000, vsg_id))

dominant <- tpm[tpm$Names != 'FALSE', ]

# Create the 4 dataframes
cow <- tpm %>%
  filter(str_starts(vsg_id, "9"))

wt_m <- tpm %>%
  filter(str_starts(vsg_id, "WT"))

brca_m <- tpm %>%
  filter(str_starts(vsg_id, "B"))

rad_m <- tpm %>%
  filter(str_starts(vsg_id, "R"))

blast <- read.table("~/Google Drive/My Drive/edinburgh/vsg/cobalt_2025/analysis/blast/start_blast/out/vsgs_blast.txt", check.names = FALSE, header = FALSE)

blast <- blast %>%
  filter(str_starts(V2, "haplotype2"))

# Define datasets and their configurations
datasets <- list(
  cow = list(
    data = cow, 
    prefix = "9",
    col_pattern = "^9",
    has_dpi = TRUE
  ),
  wt_m = list(
    data = wt_m, 
    prefix = "WT",
    col_pattern = "^WT",
    has_dpi = TRUE
  ),
  brca_m = list(
    data = brca_m, 
    prefix = "B",
    col_pattern = "^B",
    has_dpi = TRUE
  ),
  rad_m = list(
    data = rad_m, 
    prefix = "R",
    col_pattern = "^R",
    has_dpi = FALSE
  )
)

# Initialize results storage
all_results <- list()
summary_stats <- list()

# Redirect output to file
sink(file.path(output_dir, "mosaic_analysis_report.txt"))

# Print overall statistics first
cat("=== OVERALL STATISTICS ===\n")
cat("Total assembled transcripts across all experiments:", nrow(tpm), "\n")
cat("Transcripts from cow samples:", nrow(cow), "\n\n")

# Print metadata column names for debugging
cat("Metadata columns:", paste(names(meta), collapse = ", "), "\n\n")

# Loop through each dataset
for (dataset_name in names(datasets)) {
  cat("=== PROCESSING DATASET:", toupper(dataset_name), "===\n")
  
  current_config <- datasets[[dataset_name]]
  current_data <- current_config$data
  current_prefix <- current_config$prefix
  
  # Dataset-specific transcript count
  cat("Assembled transcripts in", dataset_name, ":", nrow(current_data), "\n")
  
  # Filter blast for current dataset
  blast_filtered <- blast %>%
    filter(str_starts(V1, current_prefix))
  
  # Filter reads for current dataset and calculate depth
  reads_filtered <- reads %>%
    select(vsg_id, length, starts_with(current_prefix)) %>%
    filter(str_starts(vsg_id, current_prefix))
  
  numeric_cols <- names(reads_filtered)[!names(reads_filtered) %in% c("vsg_id", "length")]
  
  if (length(numeric_cols) > 0) {
    depth_filtered <- reads_filtered
    depth_filtered[numeric_cols] <- (depth_filtered[numeric_cols] * 150) / depth_filtered$length
    
    # Apply 100x coverage and 2+ sample threshold
    core <- depth_filtered %>%
      filter(rowSums(across(all_of(numeric_cols), ~ .x >= 100), na.rm = TRUE) >= 2)
  } else {
    core <- reads_filtered
  }
  
  gold_standard_count <- length(unique(core$vsg_id))
  cat("'Gold standard' transcripts (100x coverage, 2+ samples):", gold_standard_count, "\n")
  
  # Blast analysis
  blast_core <- blast_filtered %>%
    filter(V1 %in% core$vsg_id)
  
  blast_hits_count <- length(unique(blast_core$V1))
  cat("Transcripts with blast hits:", blast_hits_count, "/", gold_standard_count)
  if (gold_standard_count > 0) {
    cat(" (", round((blast_hits_count/gold_standard_count)*100, 1), "%)\n", sep="")
  } else {
    cat("\n")
  }
  
  # Add length information for hit percentage calculation
  blast_core <- blast_core %>%
    left_join(core %>% select(vsg_id, length), by = c("V1" = "vsg_id"))
  
  blast_core$hit_perc <- (blast_core$V4 / blast_core$length) * 100
  
  # Full hits (99% identity over 99% of transcript)
  full_hits <- blast_core %>%
    filter(hit_perc >= 99 & V3 >= 99)
  
  full_hits_count <- length(unique(full_hits$V1))
  cat("Full length, complete transcripts (99% ID, 99% length):", full_hits_count)
  if (gold_standard_count > 0) {
    cat(" (", round((full_hits_count/gold_standard_count)*100, 1), "%)\n", sep="")
  } else {
    cat("\n")
  }
  
  # Good hits (99% identity, no length filter)
  good_hits <- blast_core %>%
    filter(V3 >= 99)
  
  good_hits_count <- length(unique(good_hits$V1))
  cat("Transcripts with 99% identity hits (no length filter):", good_hits_count, "\n")
  
  # Non-overlapping hits analysis
  if (nrow(good_hits) > 0) {
    non_overlapping_counts <- good_hits %>%
      group_by(V1) %>%
      summarise(
        total_hits = n(),
        non_overlapping_hits = {
          if (n() == 0) {
            0
          } else {
            ranges <- GRanges(
              seqnames = V2,
              ranges = IRanges(start = pmin(V9, V10), 
                               end = pmax(V9, V10))
            )
            reduced_ranges <- reduce(ranges)
            length(reduced_ranges)
          }
        },
        .groups = 'drop'
      )
    
    putative <- non_overlapping_counts %>%
      filter(non_overlapping_hits > 1) %>%
      filter(non_overlapping_hits < 200)
    
    putative_count <- length(unique(putative$V1))
    cat("Putative mosaics (non-overlapping hits > 1):", putative_count)
    if (gold_standard_count > 0) {
      cat(" (", round((putative_count/gold_standard_count)*100, 1), "%)\n", sep="")
    } else {
      cat("\n")
    }
  } else {
    putative_count <- 0
    cat("Putative mosaics (non-overlapping hits > 1): 0 (0%)\n")
    putative <- data.frame()
  }
  
  # Store summary statistics
  summary_stats[[dataset_name]] <- list(
    total_transcripts = nrow(current_data),
    gold_standard = gold_standard_count,
    with_blast_hits = blast_hits_count,
    full_hits = full_hits_count,
    good_hits_99 = good_hits_count,
    putative_mosaics = putative_count
  )
  
  # Save dataset-specific results
  if (exists("core") && nrow(core) > 0) {
    write.csv(core, file.path(output_dir, paste0(dataset_name, "_gold_transcripts.csv")), row.names = FALSE)
  }
  
  if (exists("full_hits") && nrow(full_hits) > 0) {
    write.csv(full_hits, file.path(output_dir, paste0(dataset_name, "_full_hits.csv")), row.names = FALSE)
  }
  
  # Create plots if we have putative mosaics
  if (nrow(putative) > 0) {
    # Save putative data first
    write.csv(putative, file.path(output_dir, paste0(dataset_name, "_putative_mosaics.csv")), row.names = FALSE)
    
    putative <- putative %>%
      separate(V1, into = c("col1", "col2", "col3"), sep = "_", remove = FALSE)
    
    cat("Sample col2 values:", head(unique(putative$col2)), "\n")
    
    # Merge putative with metadata using simple merge approach
    putative_with_meta <- merge(x = meta, y = putative, by.x = "isolate", by.y = "col2", all.y = TRUE)
    
    # Check if merge was successful and x_axis exists
    if ("x_axis" %in% names(putative_with_meta) && sum(!is.na(putative_with_meta$x_axis)) > 0) {
      cat("Successfully merged with metadata using isolate column\n")
      plot_x_var <- "x_axis"
      x_label <- "X Axis (from metadata)"
      putative_plot_data <- putative_with_meta %>% filter(!is.na(x_axis))
    } else {
      cat("Metadata merge failed or no x_axis column, using original col2 values\n")
      # Use original approach
      if (current_config$has_dpi) {
        putative_plot_data <- putative %>%
          mutate(plot_x = as.numeric(str_remove_all(col2, "d"))) %>%
          filter(!is.na(plot_x))
        plot_x_var <- "plot_x"
        x_label <- "DPI"
      } else {
        putative_plot_data <- putative
        plot_x_var <- "col2"
        x_label <- "Sample"
      }
    }
    
    # Create plots if we have data
    if (nrow(putative_plot_data) > 0) {
      # Create scatter plot
      if (current_config$has_dpi || plot_x_var == "x_axis") {
        scatter_plot <- ggscatter(putative_plot_data, x = plot_x_var, y = "non_overlapping_hits",
                                  xlab = x_label, 
                                  ylab = "Non-overlapping hits",
                                  add = "reg.line", 
                                  conf.int = FALSE,
                                  title = paste("Putative Mosaics -", toupper(dataset_name)))
        
        ggsave(file.path(output_dir, paste0(dataset_name, "_putative_mosaics_scatter.png")), 
               plot = scatter_plot, width = 6, height = 6, dpi = 300)
      }
      
      # Create box plot
      box_plot <- ggboxplot(putative_plot_data, x = plot_x_var, y = "non_overlapping_hits",
                            xlab = x_label, 
                            ylab = "Non-overlapping hits",
                            title = paste("Putative Mosaics -", toupper(dataset_name))) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
      
      ggsave(file.path(output_dir, paste0(dataset_name, "_putative_mosaics_boxplot.png")), 
             plot = box_plot, width = 6, height = 6, dpi = 300)
      
    } else {
      cat("No valid data for plotting after processing\n")
    }
  }
  
  # Store all results for this dataset
  all_results[[dataset_name]] <- list(
    data = current_data,
    core = if (exists("core")) core else data.frame(),
    blast_core = if (exists("blast_core")) blast_core else data.frame(),
    full_hits = if (exists("full_hits")) full_hits else data.frame(),
    good_hits = if (exists("good_hits")) good_hits else data.frame(),
    putative = putative,
    matches = if (exists("matches")) matches else character(),
    unique_to_dominant = if (exists("unique_to_dominant")) unique_to_dominant else character()
  )
  
  cat("\n")
}

# Print comprehensive summary table
cat("=== SUMMARY TABLE ===\n")
cat("Dataset\t\tTotal\tGold Std\tBlast Hits\tFull Hits\t99% ID Hits\tPutative Mosaics\n")
cat("-------\t\t-----\t--------\t----------\t---------\t-----------\t----------------\n")

for (dataset_name in names(datasets)) {
  stats <- summary_stats[[dataset_name]]
  cat(sprintf("%-10s\t%d\t%d\t\t%d\t\t%d\t\t%d\t\t%d\n", 
              toupper(dataset_name),
              stats$total_transcripts,
              stats$gold_standard,
              stats$with_blast_hits,
              stats$full_hits,
              stats$good_hits_99,
              stats$putative_mosaics))
}

# Calculate percentages for each dataset
cat("\n=== PERCENTAGE SUMMARY ===\n")
cat("Dataset\t\tFull Hits %\tPutative Mosaics %\n")
cat("-------\t\t-----------\t------------------\n")

for (dataset_name in names(datasets)) {
  stats <- summary_stats[[dataset_name]]
  if (stats$gold_standard > 0) {
    full_hits_pct <- round((stats$full_hits / stats$gold_standard) * 100, 1)
    mosaics_pct <- round((stats$putative_mosaics / stats$gold_standard) * 100, 1)
  } else {
    full_hits_pct <- 0
    mosaics_pct <- 0
  }
  cat(sprintf("%-10s\t%.1f%%\t\t%.1f%%\n", 
              toupper(dataset_name), full_hits_pct, mosaics_pct))
}

# Close the sink
sink()

# Save summary statistics as CSV
summary_df <- do.call(rbind, lapply(names(summary_stats), function(x) {
  data.frame(
    dataset = x,
    total_transcripts = summary_stats[[x]]$total_transcripts,
    gold_standard = summary_stats[[x]]$gold_standard,
    with_blast_hits = summary_stats[[x]]$with_blast_hits,
    full_hits = summary_stats[[x]]$full_hits,
    good_hits_99 = summary_stats[[x]]$good_hits_99,
    putative_mosaics = summary_stats[[x]]$putative_mosaics
  )
}))

write.csv(summary_df, file.path(output_dir, "summary_statistics.csv"), row.names = FALSE)

# Save all results as RData
save(all_results, summary_stats, file = file.path(output_dir, "mosaic_analysis_results.RData"))

cat("Analysis complete! All results saved to:", output_dir, "\n")