# First load required libraries
pacman::p_load(dplyr, RColorBrewer, reshape, ggpubr, tidyverse, circlize)
data_dir <- "path/to/data"
out_dir <- "path/to/results"
fig_dir <- file.path(out_dir, "figs")
for (d in c("Fig.6")) dir.create(file.path(fig_dir, d), recursive = TRUE, showWarnings = FALSE)
vsg <- read.table(file.path(data_dir, "cobalt_2025/analysis/blast/start_blast/out/vsgs_blast_filtered.txt"))
meta <- read.csv(file.path(data_dir, "vsgseq2/data/meta.csv"))
tpm <- read.csv(file.path(data_dir, "cobalt_2025/cdhit-94-94/summary/cluster/filtered_tpm_clusters.csv"), check.names = FALSE)

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
merged$x_axis <- as.numeric(merged$x_axis)

cow <- merged[merged$host == 'cow', ]
mouse <- merged[merged$host == 'mouse', ]

png(file.path(fig_dir, "Fig.6/Fig.6.1.png"), units="in", width=20, height=15, res=300)
ggboxplot(cow, x = "x_axis", y = "V3", legend = "right") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  facet_wrap(~ V13, scales = 'free') +
  ylab("Percentage identity of expressed VSG vs genome assembly") +
  xlab("DPI")
dev.off()

png(file.path(fig_dir, "Fig.6/Fig.6.2.png"), units="in", width=20, height=15, res=300)
ggscatter(cow, x = "x_axis", y = "V3", color = "V13",
          add = "loess", 
          alpha = 0.01,
          conf.int = TRUE) +
  ylab("Percentage identity of expressed VSG vs genome assembly") +
  xlab("DPI")
dev.off()

png(file.path(fig_dir, "Fig.6/Fig.6.3.png"), units="in", width=20, height=15, res=300)
ggscatter(cow, x = "x_axis", y = "V3", color = "V13",
          add = "loess", 
          alpha = 0.01,
          conf.int = TRUE) +
  facet_wrap(~ V13, scales = 'free') +
  ylab("Percentage identity of expressed VSG vs genome assembly") +
  xlab("DPI")
dev.off()