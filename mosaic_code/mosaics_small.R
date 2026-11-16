pacman::p_load(dplyr, RColorBrewer, reshape, ggpubr, tidyverse, pheatmap, GenomicRanges, ggVennDiagram)
data_dir <- "path/to/data"
out_dir <- "path/to/results"

# Create output directory
output_dir <- file.path(out_dir, "mosaic_work/small_blast")
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Function to calculate query coverage for blast hits
calculate_query_coverage <- function(transcript_id, blast_data, transcript_length) {
  # Get all hits for this transcript
  hits <- blast_data[blast_data$V1 == transcript_id, ]
  
  if (nrow(hits) == 0) {
    return(0)
  }
  
  # Extract query start and end positions (V7, V8)
  # Handle both forward and reverse orientations
  start_positions <- pmin(hits$V7, hits$V8)
  end_positions <- pmax(hits$V7, hits$V8)
  
  # Create a coverage vector
  covered <- rep(FALSE, transcript_length)
  
  # Mark covered positions
  for (i in 1:length(start_positions)) {
    start_pos <- max(1, start_positions[i])  # Ensure within bounds
    end_pos <- min(transcript_length, end_positions[i])  # Ensure within bounds
    
    if (start_pos <= transcript_length && end_pos >= 1) {
      covered[start_pos:end_pos] <- TRUE
    }
  }
  
  # Calculate coverage percentage
  coverage_bp <- sum(covered)
  coverage_percent <- (coverage_bp / transcript_length) * 100
  
  return(coverage_percent)
}

reads <- read.table(file.path(data_dir, "cobalt_2025/cdhit-94-94/summary/read_counts/reads.txt"), check.names = FALSE, header = TRUE)
tpm <- read.csv(file.path(data_dir, "cobalt_2025/cdhit-94-94/summary/tpm/tpm.csv"), check.names = FALSE, header = TRUE)
meta <- read.csv(file.path(data_dir, "cobalt_2025/analysis/mosaic_work/meta.csv"), check.names = FALSE, header = TRUE)

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

blast <- read.table(file.path(data_dir, "cobalt_2025/analysis/blast/start_blast/out/vsgs_blast_short_filtered.txt"), check.names = FALSE, header = FALSE)

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
      filter(non_overlapping_hits > 1) 
    
    putative_count <- length(unique(putative$V1))
    cat("Putative mosaics (non-overlapping hits > 1):", putative_count)
    if (gold_standard_count > 0) {
      cat(" (", round((putative_count/gold_standard_count)*100, 1), "%)\n", sep="")
    } else {
      cat("\n")
    }
    
    # Calculate query coverage for putative mosaics
    if (nrow(putative) > 0) {
      putative_with_coverage <- putative %>%
        rowwise() %>%
        mutate(
          transcript_length = core$length[core$vsg_id == V1][1],  # Get length from core
          query_coverage_percent = calculate_query_coverage(V1, good_hits, transcript_length)
        ) %>%
        ungroup()
      
      # Filter for high coverage putative mosaics (>= 100%)
      high_coverage_mosaics <- putative_with_coverage %>%
        filter(query_coverage_percent >= 100)
      
      high_coverage_count <- nrow(high_coverage_mosaics)
      
      cat("Putative mosaics with >= 95% query coverage:", high_coverage_count)
      if (putative_count > 0) {
        cat(" (", round((high_coverage_count/putative_count)*100, 1), "% of putative)\n", sep="")
      } else {
        cat("\n")
      }
      
      # Save coverage details
      if (nrow(putative_with_coverage) > 0) {
        write.csv(putative_with_coverage %>% 
                    select(V1, total_hits, non_overlapping_hits, transcript_length, query_coverage_percent), 
                  file.path(output_dir, paste0(dataset_name, "_putative_coverage_details.csv")), 
                  row.names = FALSE)
      }
      
      putative <- putative_with_coverage  # Keep the version with coverage info
    } else {
      high_coverage_count <- 0
    }
  } else {
    putative_count <- 0
    high_coverage_count <- 0
    cat("Putative mosaics (non-overlapping hits > 1): 0 (0%)\n")
    cat("Putative mosaics with >= 95% query coverage: 0\n")
    putative <- data.frame()
  }
  
  # Store summary statistics
  summary_stats[[dataset_name]] <- list(
    total_transcripts = nrow(current_data),
    gold_standard = gold_standard_count,
    with_blast_hits = blast_hits_count,
    full_hits = full_hits_count,
    good_hits_99 = good_hits_count,
    putative_mosaics = putative_count,
    high_coverage_mosaics = high_coverage_count
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
    
    # Extract sample name from vsg_id (e.g., "WT_L5" from "WT_L5_trinity.10")
    putative <- putative %>%
      mutate(sample_name = str_extract(V1, "^[^_]+_[^_]+"))
    
    cat("Sample names extracted:", head(unique(putative$sample_name)), "\n")
    
    # Merge putative with metadata using isolate column and sample_name
    putative_with_meta <- merge(x = meta, y = putative, by.x = "isolate", by.y = "sample_name", all.y = TRUE)
    
    # Check if merge was successful and axis column exists
    if ("axis" %in% names(putative_with_meta) && sum(!is.na(putative_with_meta$axis)) > 0) {
      cat("Successfully merged with metadata using isolate column\n")
      
      # Apply dataset-specific sorting
      if (dataset_name == "cow") {
        # For cow dataset: sort numerically by extracting numeric values from axis
        putative_plot_data <- putative_with_meta %>% 
          filter(!is.na(axis)) %>%
          mutate(axis_numeric = as.numeric(str_extract(axis, "\\d+"))) %>%
          arrange(axis_numeric) %>%
          mutate(axis = factor(axis, levels = unique(axis)))
        
        plot_x_var <- "axis"
        x_label <- "Days Post Infection"
      } else {
        # For other datasets: sort early then late
        putative_plot_data <- putative_with_meta %>% 
          filter(!is.na(axis)) %>%
          mutate(axis = factor(axis, levels = c("Early", "Late")))
        
        plot_x_var <- "axis"
        x_label <- "Time Point"
      }
    } else {
      cat("Metadata merge failed or no axis column, using original values\n")
      # Extract DPI or sample info from sample_name for fallback
      putative_plot_data <- putative %>%
        mutate(plot_x = str_extract(sample_name, "[^_]+$")) %>%  # Extract the part after the last underscore
        mutate(plot_x = if_else(str_detect(plot_x, "^[Ld]\\d+"), 
                                as.numeric(str_remove(plot_x, "[Ld]")), 
                                as.numeric(plot_x))) %>%
        filter(!is.na(plot_x))
      
      if (nrow(putative_plot_data) > 0) {
        plot_x_var <- "plot_x"
        x_label <- "DPI/Sample"
      } else {
        putative_plot_data <- putative
        plot_x_var <- "sample_name"
        x_label <- "Sample"
      }
    }
    
    # Create plots if we have data
    if (nrow(putative_plot_data) > 0) {
      # Create scatter plot (only for cow dataset with numeric data)
      if (dataset_name == "cow" && plot_x_var == "axis") {
        # Use axis_numeric for scatter plot to get proper regression line
        scatter_plot <- ggscatter(putative_plot_data, x = "axis_numeric", y = "non_overlapping_hits",
                                  xlab = x_label, 
                                  ylab = "Non-overlapping hits",
                                  add = "reg.line", 
                                  conf.int = FALSE,
                                  title = paste("Putative Mosaics -", toupper(dataset_name)))
        
        ggsave(file.path(output_dir, paste0(dataset_name, "_putative_mosaics_scatter.png")), 
               plot = scatter_plot, width = 6, height = 6, dpi = 300)
      } else if (current_config$has_dpi && plot_x_var != "axis") {
        # For other numeric fallback data
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
      
      # Save the updated putative data with metadata
      putative <- putative_plot_data
      
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
cat("Dataset\t\tTotal\tGold Std\tBlast Hits\tFull Hits\t99% ID Hits\tPutative Mosaics\tHigh Coverage Mosaics\n")
cat("-------\t\t-----\t--------\t----------\t---------\t-----------\t----------------\t---------------------\n")

for (dataset_name in names(datasets)) {
  stats <- summary_stats[[dataset_name]]
  cat(sprintf("%-10s\t%d\t%d\t\t%d\t\t%d\t\t%d\t\t%d\t\t\t%d\n", 
              toupper(dataset_name),
              stats$total_transcripts,
              stats$gold_standard,
              stats$with_blast_hits,
              stats$full_hits,
              stats$good_hits_99,
              stats$putative_mosaics,
              stats$high_coverage_mosaics))
}

# Calculate percentages for each dataset
cat("\n=== PERCENTAGE SUMMARY ===\n")
cat("Dataset\t\tFull Hits %\tPutative Mosaics %\tHigh Coverage Mosaics %\n")
cat("-------\t\t-----------\t------------------\t-----------------------\n")

for (dataset_name in names(datasets)) {
  stats <- summary_stats[[dataset_name]]
  if (stats$gold_standard > 0) {
    full_hits_pct <- round((stats$full_hits / stats$gold_standard) * 100, 1)
    mosaics_pct <- round((stats$putative_mosaics / stats$gold_standard) * 100, 1)
    high_cov_pct <- round((stats$high_coverage_mosaics / stats$gold_standard) * 100, 1)
  } else {
    full_hits_pct <- 0
    mosaics_pct <- 0
    high_cov_pct <- 0
  }
  cat(sprintf("%-10s\t%.1f%%\t\t%.1f%%\t\t\t%.1f%%\n", 
              toupper(dataset_name), full_hits_pct, mosaics_pct, high_cov_pct))
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
    putative_mosaics = summary_stats[[x]]$putative_mosaics,
    high_coverage_mosaics = summary_stats[[x]]$high_coverage_mosaics
  )
}))

write.csv(summary_df, file.path(output_dir, "summary_statistics.csv"), row.names = FALSE)

# Save all results as RData
save(all_results, summary_stats, file = file.path(output_dir, "mosaic_analysis_results.RData"))

cat("Analysis complete! All results saved to:", output_dir, "\n")