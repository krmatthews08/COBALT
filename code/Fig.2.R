# Load packages
pacman::p_load(dplyr, RColorBrewer, reshape, ggpubr, tidyverse)
data_dir <- "path/to/data"
out_dir <- "path/to/results"
fig_dir <- file.path(out_dir, "figs")
for (d in c("Fig.2", "Fig.2/cow_individual_vsg", "Fig.2/mouse_individual_vsg")) dir.create(file.path(fig_dir, d), recursive = TRUE, showWarnings = FALSE)

# Read in data
tpm <- read.csv(file.path(data_dir, "cobalt_2025/cdhit-94-94/summary/tpm/cluster_tpm.csv"), check.names = FALSE)
meta <- read.csv(file.path(data_dir, "vsgseq2/data/meta.csv"))

remove_list <- c("R1_9", "R1_6", "R1_3", "R1_4REP", "R1_4", "R1_5", "R1_1", "R1_7", "R1_5REP") 
meta <- subset(meta, !(isolate %in% remove_list))
#count <- read.csv(file.path(data_dir, "vsgseq2/results/new_analyse/summary/vsg_count.csv"))

# Filter rows in tpm where any numeric column has a value greater than 100000
keep_tpm <- tpm %>%
  mutate(Names = if_any(where(is.numeric), ~ . > 100000, cluster))
keep_tpm$vsg <- keep_tpm$cluster 

keep_tpm$cluster <- as.character(keep_tpm$cluster)

# Replace 'vsg_id' with 'other' where 'Names' is not TRUE
keep_tpm <- keep_tpm %>% mutate(cluster = if_else(Names == "TRUE", cluster, "other"))

# Remove 'Names' column
keep_tpm <- subset(keep_tpm, select = -c(Names) )

# Reorder factor levels of 'vsg_id' to ensure 'other' comes last
keep_tpm$cluster <- factor(keep_tpm$cluster, levels = c(unique(keep_tpm$cluster[keep_tpm$cluster != "other"]), "other"))

# Reshape data from wide to long format
long <- melt(keep_tpm, id.vars = c("cluster"), variable.name = "isolate", value.name = "P")

# Merge 'long' and 'meta' data frames
long <- merge(x=meta, y=long, by.x="isolate", by.y="variable")[]
#long <- merge(x = long, y = orthogroups_long, by.x = "vsg", by.y = "value", all.x = TRUE)

# Define color palette
cow_long <- long[long$host == 'cow', ]
mouse_long <- long[long$host == 'mouse', ]

cow_long$x_axis <- as.numeric(cow_long$x_axis)
cow_long$cluster <- as.character(cow_long$cluster)
mouse_long$cluster <- as.character(mouse_long$cluster)
cow_long$value <- as.numeric(cow_long$value)
mouse_long$value <- as.numeric(mouse_long$value)

coul <- brewer.pal(11,"Set3") 
coul <- colorRampPalette(coul)(-1 + length(unique(cow_long$cluster)))
coul <- append(coul, "black")

png(file.path(fig_dir, "Fig.2/Fig.2.1.png"), units="in", width=13, height=12, res=300)
ggbarplot(cow_long, x = "x_axis", y = "value",
          add.params = list(size = 2), color = "cluster", fill = "cluster") +
  facet_wrap(~ old_stage, scales = 'free_x', ) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), legend.position = "right") +
  scale_fill_manual(values = coul) +
  scale_color_manual(values = coul) +
  ylab("TPM") +
  xlab("DPI") +
  guides(fill = guide_legend(nrow = 48), color = guide_legend(nrow = 6)) 
dev.off()

coul <- brewer.pal(11,"Set3") 
coul <- colorRampPalette(coul)(-1 + length(unique(mouse_long$cluster)))
coul <- append(coul, "black")

mouse_long <- mouse_long %>%
  mutate(type = factor(type, levels = c("WT", "BRCA2", "RAD51")))

png(file.path(fig_dir, "Fig.2/Fig.2.2.png"), units="in", width=13, height=12, res=300)
ggbarplot(mouse_long, x = "x_axis", y = "value",
          add.params = list(size = 2), color = "cluster", fill = "cluster") +
  facet_wrap(old_stage ~ type, scales = 'free_x') +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), legend.position = "right") +
  scale_fill_manual(values = coul) +
  scale_color_manual(values = coul) +
  ylab("TPM") +
  xlab("DPI") +
  guides(fill = guide_legend(nrow = 48),color = guide_legend(nrow = 6)) 
dev.off()

cow_filt <- cow_long[cow_long$cluster != 'other', ]
cow_filt$x_axis <- as.numeric(cow_filt$x_axis)
cow_filt$value <- as.numeric(cow_filt$value)

png(file.path(fig_dir, "Fig.2/Fig.2.3.png"), units="in", width=15, height=12, res=300)
ggline(cow_filt, x = "x_axis", y = "value",
       add.params = list(size = 2), color = "type", fill = "type") +
  facet_wrap( ~ cluster, scales = 'free_y') +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ylab("TPM") +
  xlab("DPI")
dev.off()

clusters <- unique(cow_filt$cluster)

for (cl in clusters) {
  df_sub <- subset(cow_filt, cluster == cl)
  p <- ggline(df_sub, x = "x_axis", y = "value",
              add.params = list(size = 2), color = "type", fill = "type") +
    ggtitle(cl) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    ylab("TPM") +
    xlab("DPI")
  png(sprintf(file.path(fig_dir, "Fig.2/cow_individual_vsg/%s.png"), cl),
      units="in", width=8, height=6, res=300)
  print(p)
  dev.off()
}

mouse_filt <- mouse_long[mouse_long$cluster != 'other', ]
mouse_filt$value <- as.numeric(mouse_filt$value)
mouse_filt <- mouse_filt %>%
  mutate(old_stage = factor(old_stage, levels = c("Early", "Late")))

png(file.path(fig_dir, "Fig.2/Fig.2.4.png"), units="in", width=20, height=15, res=300)
ggboxplot(mouse_filt, x = "old_stage", y = "value",
          add.params = list(size = 2), color = "type", fill = "type") +
  facet_wrap( ~ cluster, scales = 'free') +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ylab("TPM") +
  xlab("Stage")
dev.off()

clusters <- unique(mouse_filt$cluster)

for (cl in clusters) {
  df_sub <- subset(mouse_filt, cluster == cl)
  df_sub$old_stage <- factor(df_sub$old_stage, levels = c("Early", "Late"))
  df_sub$type <- factor(df_sub$type, levels = c("WT", "BRACA", "RAD51"))
  p <- ggboxplot(df_sub, x = "old_stage", y = "value",
                 add.params = list(size = 2), color = "type", fill = "type") +
    facet_wrap( ~ cluster, scales = 'free') +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    ylab("TPM") +
    xlab("Infection stage")
  png(sprintf(file.path(fig_dir, "Fig.2/mouse_individual_vsg/%s.png"), cl),
      units="in", width=8, height=6, res=300)
  print(p)
  dev.off()
}

count <- read.csv(file.path(data_dir, "cobalt_2025/cdhit-94-94/summary/vsgs/vsg_count.csv"), check.names = FALSE)

count <- count %>%
  mutate(sample = str_remove(sample, "_ORF"))

count <- count %>%
  mutate(sample = str_replace_all(sample, "_d", "_"))

count <- count %>%
  left_join(meta, by = c("sample" = "test"))

cow_count <- count[count$host == 'cow', ]
mouse_count <- count[count$host == 'mouse', ]
cow_count$x_axis <- as.numeric(cow_count$x_axis )

png(file.path(fig_dir, "Fig.2/Fig.2.5.png"), units="in", width=15, height=10, res=300)
ggbarplot(cow_count, x = "x_axis", y = "VSGs",
          fill = "type", color = "black") +
  facet_wrap(~ type, scales = 'free_x') +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  ylab("Number of assembled VSG") +
  xlab("DPI") 
dev.off()

mouse_count <- mouse_count %>%
  mutate(type = factor(type, levels = c("WT", "BRCA2", "RAD51")))

png(file.path(fig_dir, "Fig.2/Fig.2.6.png"), units="in", width=15, height=10, res=300)
ggbarplot(mouse_count, x = "x_axis", y = "VSGs",
          fill = "type", color = "black") +
  facet_wrap(key ~ type, scales = 'free_x') +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  ylab("Number of assembled VSG") +
  xlab("DPI") 
dev.off()
