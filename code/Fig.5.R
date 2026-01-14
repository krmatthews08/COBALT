# First load required libraries
pacman::p_load(dplyr, RColorBrewer, reshape, ggpubr, tidyverse, circlize)
vsg <- read.table('~/Desktop/cobalt_assemblies/verkko_noref/blast/out/vsgs_blast.txt')
meta <- read.csv("~/Google Drive/My Drive/edinburgh/vsg/vsgseq2/data/meta.csv")
synteny <- read.csv('~/Desktop/cobalt_assemblies/verkko_noref/assembly/renamed/associations.tsv', sep = "\t")
tel_blast <- read.table('~/Desktop/cobalt_assemblies/verkko_noref/blast/out/tels_blast.txt')

vsg <- vsg %>%
  mutate(
    temp = ifelse(V9 > V10, V9, V10),
    V9 = ifelse(V9 > V10, V10, V9),
    V10 = temp
  ) %>%
  select(-temp)

names(vsg)[1] <- "vsg_id"

vsg <-vsg %>%
  mutate(sample = sub("_trinity.*", "", vsg_id))

# Read in data
tpm <- read.csv("~/Google Drive/My Drive/edinburgh/vsg/cobalt_2025/cdhit-94-94/summary/cluster/filtered_tpm_clusters.csv", check.names = FALSE)

remove_list <- c("R1_9", "R1_6", "R1_3", "R1_4REP", "R1_4", "R1_5", "R1_1", "R1_7", "R1_5REP") 

meta <- subset(meta, !(isolate %in% remove_list))

keep_tpm <- tpm %>%
  mutate(Names = if_any(where(is.numeric) & !matches("cluster"), ~ . > 100))

keep_tpm <- keep_tpm[keep_tpm$Names != 'FALSE', ]

merged_df <- merge(
  x = vsg,
  y = keep_tpm[, c("vsg_id", "cluster")],
  by = "vsg_id",
  all.x = TRUE
)

merged_df <- na.omit(merged_df)

result <- merged_df %>%
  group_by(cluster) %>%
  slice_max(V12, n = 1) %>%
  ungroup()

merged <- merge(
  x = merged_df,
  y = meta,
  by.x = "sample",
  by = "isolate",
  all.x = TRUE
)

dom_cow <- c("1003","1024","1051","1062","1084","1087","1088","1134","1151","116","117","119","1211","1230","1280","129","1318","1338","1388","1399","1402","1414","1440","1454","1459","1512","1615","1631","1635","1682","1691","1718","1724","1747","1792","1805","1871","1874","1888","1890","1943","1972","2093","2125","2169","2172","2222","223","2284","2315","2325","2334","2353","238","2418","253","2543","2550","2638","2696","2764","2787","2951","3023","330","378","483","486","542","546","569","570","574","583","6","609","612","652","68","683","718","775","820","821","836","854","865","879","88","907","917","95","956","979")
dom_mouse <- c("1003","1024","1051","1062","1084","1087","1088","1134","1151","116","117","119","1211","1230","1280","129","1318","1338","1388","1399","1402","1414","1440","1454","1459","1512","1615","1631","1635","1682","1691","1718","1724","1747","1792","1805","1871","1874","1888","1890","1943","1972","2093","2125","2169","2172","2222","223","2284","2315","2325","2334","2353","238","2418","253","2543","2550","2638","2696","2764","2787","2951","3023","330","378","483","486","542","546","569","570","574","583","6","609","612","652","68","683","718","775","820","821","836","854","865","879","88","907","917","95","956","979")

WT_early <- c("WT_E1", "WT_E2", "WT_E3", "WT_E4", "WT_E5", "WT_E6", "WT_E7", "WT_E8")
WT_late <- c("WT_L1", "WT_L2", "WT_L3", "WT_L4", "WT_L5", "WT_L6", "WT_L7", "WT_L8")
BRAC2_early <- c("B1_2","B1_4","B1_6","B1_8")
BRAC2_late <- c("B1_1","B1_3","B1_5","B1_7","B1_9")
RAD_early <- c("R2_4", "R2_5", "R2_6", "R2_7")
RAD_late <- c("R2_1", "R2_2REP", "R2_2", "R2_3REP", "R2_3", "R2_8", "R2_9")
cow_early <- c("908_5", "914_5", "908_6", "914_6", "900_10", "908_10", "924_10", "929_10", "908_d11", "924_d11", "900_13", "914_13", "929_13", "900_14", "914_14", "900_17", "914_d17", "924_17", "924_18", "924_19", "900_20", "908_20", "914_20", "924_20", "929_20", "900_21", "929_21", "900_24", "908_24", "914_24", "924_24", "929_24")
cow_late <- c("900_28","908_28","914_28","924_28","900_31","908_31","914_31","924_31","929_31")
cow_900 <- c("900_d4", "900_10", "900_13", "900_14", "900_17", "900_20", "900_21", "900_24", "900_26", "900_28", "900_31", "900_d33", "900_42", "900_d47", "900_54", "900_d70")
cow_908 <- c("908_d4","908_5","908_6","908_10","908_d11","908_20","908_24","908_26","908_28","908_31","908_35","908_38","908_40","908_42","908_d45","908_52","908_61","908_d70")
cow_914 <- c("914_d4","914_5","914_6","914_13","914_14","914_d17","914_20","914_24","914_26","914_28","914_31","914_33","914_35","914_42","914_d45","914_52","914_61","914_d68")
cow_924 <- c("924_d4","924_10","924_d11","924_17","924_18","924_19","924_20","924_24","924_26","924_28","924_31","924_33","924_35","924_38","924_40","924_d45","924_52","924_54","924_d66")
cow_929 <-  c("929_d4","929_10","929_13","929_20","929_21","929_24","929_26","929_31","929_d33","929_d45","929_52","929_d66")
cow <- c("900_d4","908_d4","914_d4","924_d4","929_d4","908_5","914_5","908_6","914_6","900_10","908_10","924_10","929_10","908_d11","924_d11","900_13","914_13","929_13","900_14","914_14","900_17","914_d17","924_17","924_18","924_19","900_20","908_20","914_20","924_20","929_20","900_21","929_21","900_24","908_24","914_24","924_24","929_24","900_26","908_26","914_26","924_26","929_26","900_28","908_28","914_28","924_28","900_31","908_31","914_31","924_31","929_31","900_d33","914_33","924_33","929_d33","908_35","914_35","924_35","908_38","924_38","908_40","924_40","900_42","908_42","914_42","908_d45","914_d45","924_d45","929_d45","900_d47","908_52","914_52","924_52","929_52","900_54","924_54","908_61","914_61","924_d66","929_d66","914_d68","900_d70","908_d70")
mouse <- c("B1_2", "B1_4", "B1_6", "B1_8", "R1_6", "R1_9", "R2_4", "R2_5", "R2_6", "R2_7", "WT_E1", "WT_E2", "WT_E3", "WT_E4", "WT_E5", "WT_E6", "WT_E7", "WT_E8", "B1_1", "B1_3", "B1_5", "B1_7", "B1_9", "R1_1", "R1_3", "R1_4REP", "R1_4", "R1_5REP", "R1_5", "R1_7", "R2_1", "R2_2REP", "R2_2", "R2_3REP", "R2_3", "R2_8", "R2_9", "WT_L1", "WT_L2", "WT_L3", "WT_L4", "WT_L5", "WT_L6", "WT_L7", "WT_L8")

mouse_list <- list(
  WT_early = WT_early,
  WT_late = WT_late
)

cows_list <- list(
  cow_900 = cow_900,
  cow_908 = cow_908,
  cow_914 = cow_914,
  cow_924 = cow_924,
  cow_929 = cow_929
)

cow_split <- list(
  cow_early = cow_early,
  cow_late = cow_late
)

mouse_cow_list <- list(
  cow = cow,
  mouse = mouse
)

dom_list <- list(
  dom_cow = dom_cow
)

mouse_list <- lapply(mouse_list, function(samples) {
  vsg[vsg$sample %in% samples, c("V2", "V9", "V10")]
})

cow_list <- lapply(cow_list, function(samples) {
  merged[merged$sample %in% samples, c("V2", "V9", "V10", "x_axis")]
})

cows_list <- lapply(cows_list, function(samples) {
  merged[merged$sample %in% samples, c("V2", "V9", "V10", "x_axis")]
})

mouse_cow_list <- lapply(mouse_cow_list, function(samples) {
  vsg[vsg$sample %in% samples, c("V2", "V9", "V10")]
})

filtered <- merged[merged$V3 == 100, ]
filtered <- filtered[filtered$V4 >= 800, ]

cow <- filtered[filtered$host == 'cow', c("V2", "V9", "V10", "x_axis", "type")]
mouse <- filtered[filtered$host == 'mouse', c("V2", "V9", "V10", "old_stage", "type")]

dom_list <- lapply(dom_list, function(samples){
  filtered[filtered$cluster %in% samples, c("V2", "V9", "V10",  "x_axis")]
  })

# Read and format chromosome data
chr <- read.table('~/Desktop/cobalt_assemblies/verkko_noref/assembly/chr.txt')
h1_chr <- chr[grep("haplotype1", chr$V1), ]
h2_chr <- chr[grep("haplotype2", chr$V1), ]
unassigned_chr <- chr[grep("unassigned", chr$V1), ]

h2_order <- c(
  "haplotype2-Tb11", "haplotype2-Tb10", "haplotype2-Tb9",
  "haplotype2-Tb8", "haplotype2-Tb7", "haplotype2-Tb6",
  "haplotype2-Tb5", "haplotype2-Tb4", "haplotype2-Tb3",
  "haplotype2-Tb2", "haplotype2-Tb1",
  "haplotype2-0000026", "haplotype2-0000025",
  "haplotype2-0000024", "haplotype2-0000023"
)

# Reorder haplotype2 chromosomes
h2_chr <- h2_chr[match(h2_order, h2_chr$V1), ]

# Combine back into one dataframe
chr <- rbind(h1_chr, h2_chr, unassigned_chr)

chr <- chr %>%
  mutate(V2 = 0) %>%
  select(V1, V2, V3, V5)

#chr <- chr %>%
#  separate(V1, 
#           into = c("haplotype", "chromosome"),
#           sep = "-",
#           remove = FALSE) %>%
#  select(V1, V2, V3, haplotype, chromosome)

telomere <- tel_blast %>%
  group_by(V1, V2) %>%
  dplyr::slice(which.max(V12)) %>%
  arrange(desc(V12)) %>%
  group_by(V1) %>%
  dplyr::slice(which.max(V12)) %>%
  ungroup() %>%
  mutate(
    temp = ifelse(V9 > V10, V9, V10),
    V9 = ifelse(V9 > V10, V10, V9),
    V10 = temp
  ) %>%
  select(-temp)

# Create telomere label dataframe with proper columns
telomere_label <- telomere %>%
  select(V10, V2, V9, V1)

circlize_plot <- function() {
  # Reset circos parameters
  circos.clear()
  
  # Set default small gaps
  gap_size <- c(rep(0.5, nrow(chr)-1), 1)  
  
  # Find the indices where each group starts
  h2_start <- which(grepl("haplotype2", chr$V1))[1]
  unassigned_start <- which(grepl("unassigned", chr$V1))[1]
  
  # Find the last unassigned index
  last_unassigned <- max(which(grepl("unassigned", chr$V1)))
  
  # Add larger gaps between groups
  gap_size[h2_start-1] <- 10      # Gap between haplotype1 and haplotype2
  gap_size[unassigned_start-1] <- 10  # Gap between haplotype2 and unassigned
  gap_size[last_unassigned] <- 10     # Gap between unassigned and haplotype1
  
  # Set the gaps using circos.par
  circos.par(gap.after = gap_size)
  
  # Initialize circos plot with chromosome data
  circos.genomicInitialize(chr, labels.cex = 1, sector.names = chr$V5)
  
  # Plot telomere track
  circos.genomicTrack(
    data.frame(chr = telomere$V2, start = telomere$V9, end = telomere$V10),
    ylim = c(0, 1),
    track.height = 0.1,
    panel.fun = function(region, value, ...) {
      circos.genomicRect(region, value, col = "white", border = "black")
    }
  )
  
  # Create region data frames for links
  region1 = data.frame(
    chr = synteny$Query,
    start = synteny$Q.start,
    end = synteny$Q.stop
  )
  
  region2 = data.frame(
    chr = synteny$Target,
    start = synteny$T.start,
    end = synteny$T.stop
  )
  
  # Add synteny links
  circos.genomicLink(
    region1 = region1,
    region2 = region2,
    col = rgb(0, 0, 1, 0.2),  # Blue with transparency
    border = NA
  )
  
  circos.clear()
}

png("/Users/goldriev/Google Drive/My Drive/edinburgh/vsg/cobalt_2025/figs/Fig.5/Fig.5.1.png", units="in", width=8, height=8, res=300)
circlize_plot()
dev.off()

circlize_plot <- function() {
  # Reset circos parameters
  circos.clear()
  
  # Set default small gaps
  gap_size <- c(rep(0.5, nrow(chr)-1), 1)  
  
  # Find the indices where each group starts
  h1_start <- which(grepl("haplotype1", chr$V1))[1]
  h2_start <- which(grepl("haplotype2", chr$V1))[1]
  unassigned_start <- which(grepl("unassigned", chr$V1))[1]
  
  # Find the last unassigned index
  last_unassigned <- max(which(grepl("unassigned", chr$V1)))
  
  # Add larger gaps between groups
  gap_size[h2_start-1] <- 10      # Gap between haplotype1 and haplotype2
  gap_size[unassigned_start-1] <- 10  # Gap between haplotype2 and unassigned
  gap_size[last_unassigned] <- 10     # Gap between unassigned and haplotype1
  
  # Set the gaps using circos.par
  circos.par(gap.after = gap_size)
  
  # Initialize circos plot with chromosome data
  circos.genomicInitialize(chr, labels.cex = 1, sector.names = chr$V5)
  
  circos.genomicTrack(
    data.frame(chr = telomere$V2, start = telomere$V9, end = telomere$V10),
    ylim = c(0, 1),
    track.height = 0.1,
    panel.fun = function(region, value, ...) {
      circos.genomicRect(region, value, col = "white", border = "black")
    }
  )
  
  x_axis_num <- as.numeric(unique(cow$x_axis))
  col_palette <- colorRampPalette(c("blue", "gold", "red"))(length(x_axis_num))
  names(col_palette) <- sort(x_axis_num)
  
  # Plot each type as a separate track
  track_height <- 0.1
  for(i in seq_along(cows_list)) {
    df <- cows_list[[i]]
    type_name <- names(cows_list)[i]
    
    circos.genomicTrack(
      data.frame(
        chr = df$V2,
        start = df$V9,
        end = df$V10,
        col = col_palette[as.character(df$x_axis)]
      ),
      ylim = c(0, 5),
      track.height = track_height,
      panel.fun = function(region, value, ...) {
        circos.genomicRect(
          region,
          value,
          col = NA,
          border = value$col
        )
      }
    )
  }
  legend(
    "bottomright",             # Position
    title = "Time point",      # Legend title
    legend = names(col_palette), # Labels
    fill = col_palette,        # Colors
    border = col_palette,      # Border colors
    cex = 0.8,                # Text size
    box.col = "white",        # Legend box color
    ncol = 1                  # Number of columns
  )
}

png("/Users/goldriev/Google Drive/My Drive/edinburgh/vsg/cobalt_2025/figs/Fig.5_100/Fig.5.2.png", units="in", width=20, height=15, res=300)
circlize_plot()
dev.off()

circlize_plot <- function() {
  # Reset circos parameters
  circos.clear()
  
  # Split mouse dataframe by type
  mouse_list <- split(mouse, mouse$type)
  
  # Set default small gaps
  gap_size <- c(rep(0.5, nrow(chr)-1), 1)  
  
  # Find the indices for group gaps
  h1_start <- which(grepl("haplotype1", chr$V1))[1]
  h2_start <- which(grepl("haplotype2", chr$V1))[1]
  unassigned_start <- which(grepl("unassigned", chr$V1))[1]
  last_unassigned <- max(which(grepl("unassigned", chr$V1)))
  
  # Set gaps
  gap_size[h2_start-1] <- 10
  gap_size[unassigned_start-1] <- 10
  gap_size[last_unassigned] <- 10
  
  circos.par(gap.after = gap_size)
  
  # Initialize circos plot with chromosome data
  circos.genomicInitialize(chr, labels.cex = 1, sector.names = chr$V5)
  
  circos.genomicTrack(
    data.frame(chr = telomere$V2, start = telomere$V9, end = telomere$V10),
    ylim = c(0, 1),
    track.height = 0.1,
    panel.fun = function(region, value, ...) {
      circos.genomicRect(region, value, col = "white", border = "black")
    }
  )
  
  # Create color mapping for stages
  unique_stages <- unique(mouse$old_stage)
  col_palette <- colorRampPalette(c("blue", "red"))(length(unique_stages))
  names(col_palette) <- sort(unique_stages)
  
  # Plot each type as a separate track
  track_height <- 0.1
  for(i in seq_along(mouse_list)) {
    df <- mouse_list[[i]]
    type_name <- names(mouse_list)[i]
    
    circos.genomicTrack(
      data.frame(
        chr = df$V2,
        start = df$V9,
        end = df$V10,
        col = col_palette[as.character(df$old_stage)]
      ),
      ylim = c(0, 5),
      track.height = track_height,
      panel.fun = function(region, value, ...) {
        circos.genomicRect(
          region,
          value,
          col = NA,
          border = value$col
        )
      }
    )
  }
  
  # Add legend
  legend(
    "bottomright",
    title = "Stage",
    legend = names(col_palette),
    fill = col_palette,
    border = col_palette,
    cex = 0.8,
    box.col = "white",
    ncol = 1
  )
}

png("/Users/goldriev/Google Drive/My Drive/edinburgh/vsg/cobalt_2025/figs/Fig.5_100/Fig.5.3.png", units="in", width=20, height=15, res=300)
circlize_plot()
dev.off()

circlize_plot <- function() {
  # Reset circos parameters
  circos.clear()
  
  # Split mouse dataframe by type
  mouse_list <- split(mouse, mouse$type)
  
  # Set default small gaps
  gap_size <- c(rep(0.5, nrow(chr)-1), 1)  
  
  # Find the indices for group gaps
  h1_start <- which(grepl("haplotype1", chr$V1))[1]
  h2_start <- which(grepl("haplotype2", chr$V1))[1]
  unassigned_start <- which(grepl("unassigned", chr$V1))[1]
  last_unassigned <- max(which(grepl("unassigned", chr$V1)))
  
  # Set gaps
  gap_size[h2_start-1] <- 10
  gap_size[unassigned_start-1] <- 10
  gap_size[last_unassigned] <- 10
  
  circos.par(gap.after = gap_size)
  
  # Initialize circos plot with chromosome data
  circos.genomicInitialize(chr, labels.cex = 1, sector.names = chr$V5)
  
  circos.genomicTrack(
    data.frame(chr = telomere$V2, start = telomere$V9, end = telomere$V10),
    ylim = c(0, 1),
    track.height = 0.1,
    panel.fun = function(region, value, ...) {
      circos.genomicRect(region, value, col = "white", border = "black")
    }
  )
  
  circos.genomicDensity(dom_list,
                        baseline = 0, 
                        area = FALSE, 
                        track.height = 0.2, 
                        col = c("black"), 
                        window.size = 10000,
                        count_by = "percent")

}


png("/Users/goldriev/Google Drive/My Drive/edinburgh/vsg/cobalt_2025/figs/Fig.5_100/Fig.5.4.png", units="in", width=8, height=8, res=300)
circlize_plot()
dev.off()
