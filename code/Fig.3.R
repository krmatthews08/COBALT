pacman::p_load(dplyr, RColorBrewer, reshape, ggpubr, tidyverse, pheatmap, UpSetR,ggVennDiagram)

#Persistance 

#extract clusters function
get_clusters <- function(sets, binary_data) {
  # Sets that must be present
  present <- sets
  # Sets that must be absent
  absent <- setdiff(set_names, sets)
  mask <- rep(TRUE, nrow(binary_data))
  for (set in present) {
    mask <- mask & (binary_data[[set]] == 1)
  }
  for (set in absent) {
    mask <- mask & (binary_data[[set]] == 0)
  }
  rownames(binary_data)[mask]
}

# Read in data

length <- read.csv("~/Google Drive/My Drive/edinburgh/vsg/cobalt_2025/cdhit-94-94/summary/cluster/champion_length.csv", check.names = FALSE)

# Read in data
tpm <- read.csv("~/Google Drive/My Drive/edinburgh/vsg/cobalt_2025/cdhit-94-94/summary/tpm/cluster_tpm.csv", check.names = FALSE)
meta <- read.csv("~/Google Drive/My Drive/edinburgh/vsg/vsgseq2/data/meta.csv")

remove_list <- c("R1_9", "R1_6", "R1_3", "R1_4REP", "R1_4", "R1_5", "R1_1", "R1_7", "R1_5REP") 
meta <- subset(meta, !(isolate %in% remove_list))

tpm <- tpm %>%
  left_join(length %>% select(cluster, length), by = "cluster")

# Filter rows in tpm where any numeric column has a value greater than 100
keep_tpm <- tpm %>%
  mutate(Names = if_any(where(is.numeric), ~ . > 100, cluster))

keep_tpm$vsg <- keep_tpm$cluster 

keep_tpm$cluster <- as.character(keep_tpm$cluster)

# Replace 'vsg_id' with 'other' where 'Names' is not TRUE
keep_tpm <- keep_tpm %>% mutate(cluster = if_else(Names == "TRUE", cluster, "other"))

# Remove 'Names' column
keep_tpm <- subset(keep_tpm, select = -c(Names) )

# Reorder factor levels of 'vsg_id' to ensure 'other' comes last
keep_tpm$cluster <- factor(keep_tpm$cluster, levels = c(unique(keep_tpm$cluster[keep_tpm$cluster != "other"]), "other"))

# Reshape data from wide to long format
long <- melt(tpm, id.vars = c("cluster", "length"), variable.name = "isolate", value.name = "P")

# Merge 'long' and 'meta' data frames
long <- merge(x=meta, y=long, by.x="isolate", by.y="variable")[]

# Define color palette
cow_long <- long[long$host == 'cow', ]
mouse_long <- long[long$host == 'mouse', ]

cow_long$x_axis <- as.numeric(cow_long$x_axis)
cow_long$cluster <- as.character(cow_long$cluster)
cow_long$value <- as.numeric(cow_long$value)

cow_long <- cow_long %>%
  mutate(value_adjusted = ifelse(value < 100, 0, value))

summed_data <- cow_long %>%
  group_by(cluster, x_axis, type) %>%
  summarise(value = sum(value_adjusted, na.rm = TRUE),     
            mean_length = mean(length, na.rm = TRUE),
            .groups = "drop")

# Create a single data frame with three entries for each cluster and type
cluster_type_summary <- summed_data %>%
  filter(value > 0) %>%  # Keep only rows where value > 0
  group_by(cluster, type) %>%
  summarise(
    # Row with the highest value
    mean_length = mean_length[which.max(value)],
    highest_x_axis = x_axis[which.max(value)],
    highest_value = max(value, na.rm = TRUE),
    
    # Row with the lowest x_axis
    earliest_x_axis = min(x_axis, na.rm = TRUE),
    earliest_value = value[which.min(x_axis)],
    
    # Row with the highest x_axis
    latest_x_axis = max(x_axis, na.rm = TRUE),
    latest_value = value[which.max(x_axis)],
    .groups = "drop"
  )
cluster_type_summary

cluster_type_summary <- cluster_type_summary %>%
  mutate(last_axis = latest_x_axis - highest_x_axis )

cluster_type_summary <- cluster_type_summary %>%
  mutate(first_axis = earliest_x_axis - highest_x_axis )

cluster_type_summary$highest_0_axis <- 0

cluster_type_summary <- cluster_type_summary %>%
  mutate(persistance = last_axis - first_axis)

cluster_type_summary$highest_x_axis <- as.numeric(cluster_type_summary$highest_x_axis)
cluster_type_summary$persistance_log <- log10(1+cluster_type_summary$persistance)

cluster_type_filtered <- cluster_type_summary %>% 
  filter(persistance > 1)

png("/Users/goldriev/Google Drive/My Drive/edinburgh/vsg/cobalt_2025/figs/Fig.3/Fig.3.1.png", units="in", width=6, height=6, res=300)
ggscatter(cluster_type_filtered, x = "persistance", y = "mean_length", color = 'type',
          add = "loess", 
          conf.int = T,
          add.params = list(size = 2), alpha = 0.5, 
          position = position_jitter( height = 0)) +  # Vertical dodging
  theme(axis.text.x = element_text(angle = 90, vjust = 0.1, hjust = 1)) +
  xlab("Persistance") +
  ylab("Length")
dev.off()

png("/Users/goldriev/Google Drive/My Drive/edinburgh/vsg/cobalt_2025/figs/Fig.3/Fig.3.2.png", units="in", width=6, height=6, res=300)
ggscatter(cluster_type_filtered, x = "highest_x_axis", y = "mean_length", color = 'type',
          add = "loess", 
          conf.int = T,
          add.params = list(size = 2), alpha = 0.5, 
          position = position_jitter( height = 0)) +  # Vertical dodging
  theme(axis.text.x = element_text(angle = 90, vjust = 0.1, hjust = 1)) +
  xlab("DPI") +
  ylab("Length")
dev.off()

png("/Users/goldriev/Google Drive/My Drive/edinburgh/vsg/cobalt_2025/figs/Fig.3/Fig.3.3.png", units="in", width=6, height=6, res=300)
ggscatter(cluster_type_filtered, x = "highest_x_axis", y = "persistance", color = 'type',
          add = "loess", 
          conf.int = T,
          add.params = list(size = 2), alpha = 0.5, 
          position = position_jitter( height = 0)) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.1, hjust = 1)) +
  xlab("DPI") +
  ylab("Persistance")
dev.off()


diverse_dpi <- cluster_type_filtered[cluster_type_filtered$highest_x_axis > 20 & cluster_type_filtered$highest_x_axis < 50, ]

png("/Users/goldriev/Google Drive/My Drive/edinburgh/vsg/cobalt_2025/figs/Fig.3/Fig.3.4.png", units="in", width=6, height=6, res=300)
ggscatter(diverse_dpi, x = "highest_x_axis", y = "persistance", color = 'type',
          add = "loess", 
          conf.int = T,
          add.params = list(size = 2), alpha = 0.5, 
          position = position_jitter( height = 0)) +  # Vertical dodging
  theme(axis.text.x = element_text(angle = 90, vjust = 0.1, hjust = 1)) +
  xlab("DPI") +
  ylab("Persistance")
dev.off()


png("/Users/goldriev/Google Drive/My Drive/edinburgh/vsg/cobalt_2025/figs/Fig.3/Fig.3.5.png", units="in", width=6, height=6, res=300)
ggscatter(cluster_type_summary, x = "highest_value", y = "persistance",color = 'type',
          add = "loess", 
          conf.int = T,
          add.params = list(size = 2), alpha = 0.5,
          position = position_jitter( height = 0)) +
  scale_x_log10() +
  xlab("TPM (log10)") +
  ylab("Persistance")
dev.off()

cluster_type_summary <- cluster_type_summary %>%
  mutate(
    DPI = case_when(
      highest_x_axis >= 0 & highest_x_axis <= 10 ~ "0-10",
      highest_x_axis > 10 & highest_x_axis <= 20 ~ "10-20",
      highest_x_axis > 20 & highest_x_axis <= 30 ~ "20-30",
      highest_x_axis > 30 & highest_x_axis <= 40 ~ "30-40",
      highest_x_axis > 40 & highest_x_axis <= 50 ~ "40-50",
      highest_x_axis > 50 & highest_x_axis <= 60 ~ "50-60",
      highest_x_axis > 60 & highest_x_axis <= 70 ~ "60-70",
      TRUE ~ NA_character_  # Handle cases outside the specified ranges
    )
  )

cluster_type_summary <- cluster_type_summary %>%
  mutate(
    expression = case_when(
      highest_value >= 0 & highest_value <= 100000 ~ "0-10",
      highest_value > 100000 & highest_value <= 200000 ~ "10-20",
      highest_value > 200000  & highest_value <= 300000  ~ "20-30",
      highest_value > 300000  & highest_value <= 400000  ~ "30-40",
      highest_value > 400000  & highest_value <= 500000  ~ "40-50",
      highest_value > 500000  & highest_value <= 600000  ~ "50-60",
      highest_value > 600000  & highest_value <= 700000  ~ "60-70",
      highest_value > 700000  & highest_value <= 800000  ~ "70-80",
      highest_value > 900000  & highest_value <= 1000000  ~ "90-100",
      TRUE ~ NA_character_  # Handle cases outside the specified ranges
    )
  )

long_form_data <- cluster_type_summary %>%
  pivot_longer(
    cols = c(first_axis, last_axis, highest_0_axis),
    names_to = "axis_type",
    values_to = "value"
  ) %>%
  select(cluster, persistance, type, DPI, expression, axis_type, value)

long_form_filtered <- long_form_data %>% 
  filter(value > 0 | value < 0)

png("/Users/goldriev/Google Drive/My Drive/edinburgh/vsg/cobalt_2025/figs/persistance/persistance_density.png", units="in", width=10, height=10, res=300)
ggdensity(long_form_filtered, x = "value", rug = F,
          color = "DPI") +
  ylab("Density") +
  xlab("Persistance") 
dev.off()

png("/Users/goldriev/Google Drive/My Drive/edinburgh/vsg/cobalt_2025/figs/persistance/persistance_density_facet.png", units="in", width=10, height=10, res=300)
ggdensity(long_form_filtered, x = "value", rug = F,
          color = "DPI") +
  ylab("Density") +
  xlab("Persistance") +
  facet_wrap(~ type, scales = 'free')
dev.off()

png("/Users/goldriev/Google Drive/My Drive/edinburgh/vsg/cobalt_2025/figs/Fig.3/Fig.3.6.png", units="in", width=6, height=6, res=300)
ggboxplot(long_form_filtered, x = "DPI", y = "persistance",
          color = "DPI", order = c("0-10", "10-20", "20-30", "30-40", "40-50", "50-60", "60-70")) +
  ylab("Persistance") +
  xlab("DPI") +
  facet_wrap(~ type, scales = 'free') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

png("/Users/goldriev/Google Drive/My Drive/edinburgh/vsg/cobalt_2025/figs/Fig.3/Fig.3.7.png", units="in", width=6, height=6, res=300)
ggboxplot(long_form_filtered, x = "DPI", y = "persistance",
          color = "DPI", order = c("0-10", "10-20", "20-30", "30-40", "40-50", "50-60", "60-70")) +
  ylab("Persistance") +
  xlab("DPI")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

# expression

long_form_filtered <- na.omit(long_form_filtered)

png("/Users/goldriev/Google Drive/My Drive/edinburgh/vsg/cobalt_2025/figs/persistance/expression_density.png", units="in", width=10, height=10, res=300)
ggdensity(long_form_filtered, x = "value", rug = F,
          color = "expression") +
  ylab("Density") +
  xlab("Persistance") 
dev.off()

png("/Users/goldriev/Google Drive/My Drive/edinburgh/vsg/cobalt_2025/figs/persistance/expression_density_facet.png", units="in", width=10, height=10, res=300)
ggdensity(long_form_filtered, x = "value", rug = F,
          color = "expression") +
  ylab("Density") +
  xlab("Persistance") +
  facet_wrap(~ type, scales = 'free')
dev.off()

png("/Users/goldriev/Google Drive/My Drive/edinburgh/vsg/cobalt_2025/figs/Fig.3/Fig.3.8.png", units="in", width=6, height=6, res=300)
ggboxplot(long_form_filtered, x = "expression", y = "persistance",
          color = "expression", order = c("0-10", "10-20", "20-30", "30-40", "40-50", "50-60", "60-70", "70-80", "80-90", "90-100")) +
  ylab("Persistance") +
  xlab("VSG percentage of population") +
  facet_wrap(~ type, scales = 'free') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

png("/Users/goldriev/Google Drive/My Drive/edinburgh/vsg/cobalt_2025/figs/Fig.3/Fig.3.9.png", units="in", width=6, height=6, res=300)
ggboxplot(long_form_filtered, x = "expression", y = "persistance",
          color = "expression", order = c("0-10", "10-20", "20-30", "30-40", "40-50", "50-60", "60-70", "70-80", "80-90", "90-100")) +
  ylab("Persistance") +
  xlab("VSG percentage of population") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

tpm_reduced <- tpm %>%
  group_by(cluster) %>%
  summarize(across(where(is.numeric), sum, na.rm = TRUE), .groups = "drop")

tpm_reduced <- tpm_reduced[, !names(tpm_reduced) %in% "length"]

rownames(tpm_reduced) <- tpm_reduced$cluster

tpm_reduced <- tpm_reduced %>%
  select(-cluster)

tpm_reduced <- tpm_reduced %>%
  select(!any_of(remove_list))

meta_subset <- meta %>%
  select(isolate, host, key, type)

cow_samples <- meta %>%
  filter(host == "cow") %>%
  pull(isolate) 

mouse_samples <- meta %>%
  filter(host == "mouse") %>%
  pull(isolate) 

cow_samples <- meta %>%
  filter(host == "cow") %>%
  pull(isolate) 

wt_samples <- meta %>%
  filter(host == "cow") %>%
  pull(isolate) 

tpm_reduced_cow <- tpm_reduced %>%
  select(any_of(cow_samples))

tpm_reduced_mouse <- tpm_reduced %>%
  select(any_of(mouse_samples))

rownames(meta_subset) <- meta_subset$isolate
meta_subset$isolate <- NULL #

cow_order <- c("900_d4", "908_d4", "914_d4", "924_d4", "929_d4",
               "908_5", "914_5",
               "908_6", "914_6",
               "900_10", "908_10", "924_10", "929_10",
               "908_d11", "924_d11",
               "900_13", "914_13", "929_13",
               "900_14", "914_14",
               "914_d17", "924_17",
               "924_18", "924_19",
               "900_20", "908_20", "914_20", "924_20", "929_20",
               "900_21", "929_21",
               "900_24", "908_24", "914_24", "924_24", "929_24",
               "900_26", "908_26", "914_26", "924_26", "929_26",
               "900_28", "908_28", "914_28", "924_28",
               "900_31", "908_31", "914_31", "924_31", "929_31",
               "900_d33", "914_33", "924_33", "929_d33",
               "908_35", "914_35", "924_35",
               "924_38",
               "908_40", "924_40",
               "900_42", "908_42", "914_42",
               "908_d45", "914_d45", "924_d45", "929_d45",
               "900_d47",
               "908_52", "914_52", "924_52", "929_52",
               "900_54", "924_54",
               "908_61", "914_61",
               "914_d68")

mouse_order <- c("WT_E1", "WT_E2", "WT_E3", "WT_E4", "WT_E5", "WT_E6", "WT_E7", "WT_E8", "WT_L1", "WT_L2", "WT_L3", "WT_L4", "WT_L5", "WT_L6", "WT_L7", "WT_L8",
                 "B1_2", "B1_4", "B1_6", "B1_8", "B1_1", "B1_3", "B1_5", "B1_7", "B1_9",
                 "R2_4", "R2_5", "R2_6", "R2_7", "R2_1", "R2_2", "R2_3", "R2_8", "R2_9", "R2_2REP", "R2_3REP")


order <- c("900_d4", "908_d4", "914_d4", "924_d4", "929_d4",
           "908_5", "914_5",
           "908_6", "914_6",
           "900_10", "908_10", "924_10", "929_10",
           "908_d11", "924_d11",
           "900_13", "914_13", "929_13",
           "900_14", "914_14",
           "914_d17", "924_17",
           "924_18", "924_19",
           "900_20", "908_20", "914_20", "924_20", "929_20",
           "900_21", "929_21",
           "900_24", "908_24", "914_24", "924_24", "929_24",
           "900_26", "908_26", "914_26", "924_26", "929_26",
           "900_28", "908_28", "914_28", "924_28",
           "900_31", "908_31", "914_31", "924_31", "929_31",
           "900_d33", "914_33", "924_33", "929_d33",
           "908_35", "914_35", "924_35",
           "924_38",
           "908_40", "924_40",
           "900_42", "908_42", "914_42",
           "908_d45", "914_d45", "924_d45", "929_d45",
           "900_d47",
           "908_52", "914_52", "924_52", "929_52",
           "900_54", "924_54",
           "908_61", "914_61",
           "914_d68",
           "WT_E1", "WT_E2", "WT_E3", "WT_E4", "WT_E5", "WT_E6", "WT_E7", "WT_E8", "WT_L1", "WT_L2", "WT_L3", "WT_L4", "WT_L5", "WT_L6", "WT_L7", "WT_L8",
           "B1_2", "B1_4", "B1_6", "B1_8", "B1_1", "B1_3", "B1_5", "B1_7", "B1_9",
           "R2_4", "R2_5", "R2_6", "R2_7", "R2_1", "R2_2", "R2_3", "R2_8", "R2_9", "R2_2REP", "R2_3REP")

wt_order <- c("WT_E1", "WT_E2", "WT_E3", "WT_E4", "WT_E5", "WT_E6", "WT_E7", "WT_E8", "WT_L1", "WT_L2", "WT_L3", "WT_L4", "WT_L5", "WT_L6", "WT_L7", "WT_L8",
              "900_d4", "908_d4", "914_d4", "924_d4", "929_d4",
              "908_5", "914_5",
              "908_6", "914_6",
              "900_10", "908_10", "924_10", "929_10",
              "908_d11", "924_d11",
              "900_13", "914_13", "929_13",
              "900_14", "914_14",
              "914_d17", "924_17",
              "924_18", "924_19",
              "900_20", "908_20", "914_20", "924_20", "929_20",
              "900_21", "929_21",
              "900_24", "908_24", "914_24", "924_24", "929_24",
              "900_26", "908_26", "914_26", "924_26", "929_26",
              "900_28", "908_28", "914_28", "924_28",
              "900_31", "908_31", "914_31", "924_31", "929_31",
              "900_d33", "914_33", "924_33", "929_d33",
              "908_35", "914_35", "924_35",
              "924_38",
              "908_40", "924_40",
              "900_42", "908_42", "914_42",
              "908_d45", "914_d45", "924_d45", "929_d45",
              "900_d47",
              "908_52", "914_52", "924_52", "929_52",
              "900_54", "924_54",
              "908_61", "914_61",
              "914_d68")

cow_day_order <- c(
  # 900_* values
  "900_d4", "900_10", "900_13", "900_14", "900_20", "900_21", 
  "900_24", "900_26", "900_28", "900_31", "900_d33", "900_42", 
  "900_d47", "900_54",
  
  # 908_* values  
  "908_d4", "908_5", "908_6", "908_10", "908_d11", "908_20", 
  "908_24", "908_26", "908_28", "908_31", "908_35", "908_40", 
  "908_42", "908_d45", "908_52", "908_61",
  
  # 914_* values
  "914_d4", "914_5", "914_6", "914_13", "914_14", "914_d17", 
  "914_20", "914_24", "914_26", "914_28", "914_31", "914_33", 
  "914_35", "914_42", "914_d45", "914_52", "914_61", "914_d68",
  
  # 924_* values
  "924_d4", "924_10", "924_d11", "924_17", "924_18", "924_19", 
  "924_20", "924_24", "924_26", "924_28", "924_31", "924_33", 
  "924_35", "924_38", "924_40", "924_d45", "924_52", "924_54",
  
  # 929_* values
  "929_d4", "929_10", "929_13", "929_20", "929_21", "929_24", 
  "929_26", "929_31", "929_d33", "929_d45", "929_52"
)

d4_vs_mice <- c("900_d4", "908_d4", "914_d4", "924_d4", "929_d4",
                "WT_E1", "WT_E2", "WT_E3", "WT_E4", "WT_E5", "WT_E6", "WT_E7", "WT_E8", "WT_L1", "WT_L2", "WT_L3", "WT_L4", "WT_L5", "WT_L6", "WT_L7", "WT_L8",
                "B1_2", "B1_4", "B1_6", "B1_8", "B1_1", "B1_3", "B1_5", "B1_7", "B1_9",
                "R2_4", "R2_5", "R2_6", "R2_7", "R2_1", "R2_2", "R2_3", "R2_8", "R2_9", "R2_2REP", "R2_3REP")

d4_10_vs_mice <- c("900_d4", "908_d4", "914_d4", "924_d4", "929_d4",
                   "900_10", "908_10", "924_10", "929_10",
                   "914_14",
                   "WT_E1", "WT_E2", "WT_E3", "WT_E4", "WT_E5", "WT_E6", "WT_E7", "WT_E8", "WT_L1", "WT_L2", "WT_L3", "WT_L4", "WT_L5", "WT_L6", "WT_L7", "WT_L8",
                   "B1_2", "B1_4", "B1_6", "B1_8", "B1_1", "B1_3", "B1_5", "B1_7", "B1_9",
                   "R2_4", "R2_5", "R2_6", "R2_7", "R2_1", "R2_2", "R2_3", "R2_8", "R2_9", "R2_2REP", "R2_3REP")

d4_vs_early_mice <- c("900_d4", "908_d4", "914_d4", "924_d4", "929_d4",
                      "WT_E1", "WT_E2", "WT_E3", "WT_E4", "WT_E5", "WT_E6", "WT_E7", "WT_E8",
                      "B1_2", "B1_4", "B1_6", "B1_8", "B1_1", "B1_3", "B1_5", "B1_7", "B1_9",
                      "R2_4", "R2_5", "R2_6", "R2_7", "R2_1", "R2_2", "R2_3", "R2_8", "R2_9", "R2_2REP", "R2_3REP")

d4_10_vs_early_mice <- c("900_d4", "908_d4", "914_d4", "924_d4", "929_d4",
                         "900_10", "908_10", "924_10", "929_10",
                         "914_14",
                         "WT_E1", "WT_E2", "WT_E3", "WT_E4", "WT_E5", "WT_E6", "WT_E7", "WT_E8",
                         "B1_2", "B1_4", "B1_6", "B1_8", "B1_1", "B1_3", "B1_5", "B1_7", "B1_9",
                         "R2_4", "R2_5", "R2_6", "R2_7", "R2_1", "R2_2", "R2_3", "R2_8", "R2_9", "R2_2REP", "R2_3REP")



tpm_reduced_cow <- tpm_reduced_cow[, cow_order]
tpm_reduced_cow_day <- tpm_reduced_cow[, cow_day_order]
tpm_reduced_mouse <- tpm_reduced_mouse[, mouse_order]
tpm_reduced <- tpm_reduced[, order]
tpm_reduced_wt <- tpm_reduced[, wt_order]
tpm_reduced_d4 <- tpm_reduced[, d4_vs_mice]
tpm_reduced_d4_10 <- tpm_reduced[, d4_10_vs_mice]
tpm_reduced_d4_early <- tpm_reduced[, d4_vs_early_mice]
tpm_reduced_d4_10_early <- tpm_reduced[, d4_10_vs_early_mice]

tpm_reduced_cow_day_filtered <- tpm_reduced_cow_day %>%
  filter(rowSums(across(everything())) != 0)

tpm_reduced_cow_filtered <- tpm_reduced_cow %>%
  filter(rowSums(across(everything())) != 0)

tpm_reduced_mouse_filtered <- tpm_reduced_mouse %>%
  filter(rowSums(across(everything())) != 0)

tpm_reduced_d4_filtered <- tpm_reduced_d4 %>%
  filter(rowSums(across(everything())) != 0)

tpm_reduced_d4_10_filtered <- tpm_reduced_d4_10 %>%
  filter(rowSums(across(everything())) != 0)

tpm_reduced_d4_early_filtered <- tpm_reduced_d4_early %>%
  filter(rowSums(across(everything())) != 0)

tpm_reduced_d4_10_early_filtered <- tpm_reduced_d4_10_early %>%
  filter(rowSums(across(everything())) != 0)

tpm_reduced_wt_filtered <- tpm_reduced_wt %>%
  filter(rowSums(across(everything())) != 0)

png("/Users/goldriev/Google Drive/My Drive/edinburgh/vsg/cobalt_2025/figs/Fig.3/Fig.3.10.png", units="in", width=11, height=6, res=300)
pheatmap(log (1 + tpm_reduced),
         col = rev(brewer.pal(8, 'PuOr')),
         cluster_rows = T, cluster_cols = F,
         clustering_distance_cols = 'euclidean',
         clustering_distance_rows = 'euclidean',
         clustering_method = 'ward.D',
         annotation_col = meta_subset,
         annotation_names_row = F,
         annotation_names_col = T,
         fontsize_row = 10,          # row label font size
         fontsize_col = 7,          # column label font size
         angle_col = 90, # sample names at an angle
         show_colnames = T, show_rownames = F, scale = "row")
dev.off()

png("/Users/goldriev/Google Drive/My Drive/edinburgh/vsg/cobalt_2025/figs/Fig.3/Fig.3.11.png", units="in", width=6, height=6, res=300)
pheatmap(log (1 + tpm_reduced_mouse_filtered),
         col = rev(brewer.pal(8, 'PuOr')),
         cluster_rows = T, cluster_cols = F,
         clustering_distance_cols = 'euclidean',
         clustering_distance_rows = 'euclidean',
         clustering_method = 'ward.D',
         annotation_col = meta_subset,
         annotation_names_row = F,
         annotation_names_col = T,
         fontsize_row = 10,          # row label font size
         fontsize_col = 7,          # column label font size
         angle_col = 45, # sample names at an angle
         show_colnames = T, 
         show_rownames = F, 
         scale = "row")
dev.off()

png("/Users/goldriev/Google Drive/My Drive/edinburgh/vsg/cobalt_2025/figs/Fig.3/Fig.3.12.png", units="in", width=8, height=6, res=300)
pheatmap(log (1 + tpm_reduced_cow_filtered),
         col = rev(brewer.pal(8, 'PuOr')),
         cluster_rows = T, cluster_cols = F,
         clustering_distance_cols = 'euclidean',
         clustering_distance_rows = 'euclidean',
         clustering_method = 'ward.D',
         annotation_col = meta_subset,
         annotation_names_row = F,
         annotation_names_col = T,
         fontsize_row = 10,          # row label font size
         fontsize_col = 7,          # column label font size
         angle_col = 90, # sample names at an angle
         show_colnames = T, show_rownames = F, scale = "row")
dev.off()

png("/Users/goldriev/Google Drive/My Drive/edinburgh/vsg/cobalt_2025/figs/Fig.3/Fig.3.13.png", units="in", width=11, height=6, res=300)
pheatmap(log (1 + tpm_reduced_wt),
         col = rev(brewer.pal(8, 'PuOr')),
         cluster_rows = T, cluster_cols = F,
         clustering_distance_cols = 'euclidean',
         clustering_distance_rows = 'euclidean',
         clustering_method = 'ward.D',
         annotation_col = meta_subset,
         annotation_names_row = F,
         annotation_names_col = T,
         fontsize_row = 10,          # row label font size
         fontsize_col = 7,          # column label font size
         angle_col = 45, # sample names at an angle
         show_colnames = T, show_rownames = F, scale = "row")
dev.off()

png("/Users/goldriev/Google Drive/My Drive/edinburgh/vsg/cobalt_2025/figs/Fig.3/Fig.3.14.png", units="in", width=11, height=6, res=300)
pheatmap(log (1 + tpm_reduced),
         col = rev(brewer.pal(8, 'PuOr')),
         cluster_rows = T, cluster_cols = T,
         clustering_distance_cols = 'euclidean',
         clustering_distance_rows = 'euclidean',
         clustering_method = 'ward.D',
         annotation_col = meta_subset,
         annotation_names_row = F,
         annotation_names_col = T,
         fontsize_row = 10,          # row label font size
         fontsize_col = 7,          # column label font size
         angle_col = 90, # sample names at an angle
         show_colnames = T, show_rownames = F, scale = "row")
dev.off()

png("/Users/goldriev/Google Drive/My Drive/edinburgh/vsg/cobalt_2025/figs/Fig.3/Fig.3.15.png", units="in", width=11, height=6, res=300)
pheatmap(log (1 + tpm_reduced_cow_filtered),
         col = rev(brewer.pal(8, 'PuOr')),
         cluster_rows = T, cluster_cols = T,
         clustering_distance_cols = 'euclidean',
         clustering_distance_rows = 'euclidean',
         clustering_method = 'ward.D',
         annotation_col = meta_subset,
         annotation_names_row = F,
         annotation_names_col = T,
         fontsize_row = 10,          # row label font size
         fontsize_col = 7,          # column label font size
         angle_col = 45, # sample names at an angle
         show_colnames = T, show_rownames = F, scale = "row")
dev.off()

png("/Users/goldriev/Google Drive/My Drive/edinburgh/vsg/cobalt_2025/figs/Fig.3/Fig.3.16.png", units="in", width=11, height=6, res=300)
pheatmap(log (1 + tpm_reduced_mouse_filtered),
         col = rev(brewer.pal(8, 'PuOr')),
         cluster_rows = T, cluster_cols = T,
         clustering_distance_cols = 'euclidean',
         clustering_distance_rows = 'euclidean',
         clustering_method = 'ward.D',
         annotation_col = meta_subset,
         annotation_names_row = F,
         annotation_names_col = T,
         fontsize_row = 10,          # row label font size
         fontsize_col = 7,          # column label font size
         angle_col = 45, # sample names at an angle
         show_colnames = T, show_rownames = F, scale = "row")
dev.off()

# Personalised heatmaps

png("/Users/goldriev/Google Drive/My Drive/edinburgh/vsg/cobalt_2025/figs/Fig.3/Fig.3.17.png", units="in", width=11, height=6, res=300)
pheatmap(log (1 + tpm_reduced_d4_filtered),
         col = rev(brewer.pal(8, 'PuOr')),
         cluster_rows = T, cluster_cols = T,
         clustering_distance_cols = 'euclidean',
         clustering_distance_rows = 'euclidean',
         clustering_method = 'ward.D',
         annotation_col = meta_subset,
         annotation_names_row = F,
         annotation_names_col = T,
         fontsize_row = 10,          # row label font size
         fontsize_col = 7,          # column label font size
         angle_col = 45, # sample names at an angle
         show_colnames = T, show_rownames = F, scale = "row")
dev.off()

png("/Users/goldriev/Google Drive/My Drive/edinburgh/vsg/cobalt_2025/figs/Fig.3/Fig.3.18.png", units="in", width=11, height=6, res=300)
pheatmap(log (1 + tpm_reduced_d4_10_filtered),
         col = rev(brewer.pal(8, 'PuOr')),
         cluster_rows = T, cluster_cols = T,
         clustering_distance_cols = 'euclidean',
         clustering_distance_rows = 'euclidean',
         clustering_method = 'ward.D',
         annotation_col = meta_subset,
         annotation_names_row = F,
         annotation_names_col = T,
         fontsize_row = 10,          # row label font size
         fontsize_col = 7,          # column label font size
         angle_col = 45, # sample names at an angle
         show_colnames = T, show_rownames = F, scale = "row")
dev.off()

png("/Users/goldriev/Google Drive/My Drive/edinburgh/vsg/cobalt_2025/figs/Fig.3/Fig.3.19.png", units="in", width=11, height=6, res=300)
pheatmap(log (1 + tpm_reduced_d4_early_filtered),
         col = rev(brewer.pal(8, 'PuOr')),
         cluster_rows = T, cluster_cols = T,
         clustering_distance_cols = 'euclidean',
         clustering_distance_rows = 'euclidean',
         clustering_method = 'ward.D',
         annotation_col = meta_subset,
         annotation_names_row = F,
         annotation_names_col = T,
         fontsize_row = 10,          # row label font size
         fontsize_col = 7,          # column label font size
         angle_col = 45, # sample names at an angle
         show_colnames = T, show_rownames = F, scale = "row")
dev.off()

png("/Users/goldriev/Google Drive/My Drive/edinburgh/vsg/cobalt_2025/figs/Fig.3/Fig.3.20.png", units="in", width=11, height=6, res=300)
pheatmap(log (1 + tpm_reduced_d4_10_early_filtered),
         col = rev(brewer.pal(8, 'PuOr')),
         cluster_rows = T, cluster_cols = T,
         clustering_distance_cols = 'euclidean',
         clustering_distance_rows = 'euclidean',
         clustering_method = 'ward.D',
         annotation_col = meta_subset,
         annotation_names_row = F,
         annotation_names_col = T,
         fontsize_row = 10,          # row label font size
         fontsize_col = 7,          # column label font size
         angle_col = 45, # sample names at an angle
         show_colnames = T, show_rownames = F, scale = "row")
dev.off()

png("/Users/goldriev/Google Drive/My Drive/edinburgh/vsg/cobalt_2025/figs/Fig.3/wt_lines.png", units="in", width=11, height=6, res=300)
pheatmap(log (1 + tpm_reduced_wt_filtered),
         col = rev(brewer.pal(8, 'PuOr')),
         cluster_rows = T, cluster_cols = T,
         clustering_distance_cols = 'euclidean',
         clustering_distance_rows = 'euclidean',
         clustering_method = 'ward.D',
         annotation_col = meta_subset,
         annotation_names_row = F,
         annotation_names_col = T,
         fontsize_row = 10,          # row label font size
         fontsize_col = 7,          # column label font size
         angle_col = 45, # sample names at an angle
         show_colnames = T, show_rownames = F, scale = "row")
dev.off()

png("/Users/goldriev/Google Drive/My Drive/edinburgh/vsg/cobalt_2025/figs/Fig.3/cow_day.png", units="in", width=11, height=6, res=300)
pheatmap(log (1 + tpm_reduced_cow_day_filtered),
         col = rev(brewer.pal(8, 'PuOr')),
         cluster_rows = T, cluster_cols = F,
         clustering_distance_cols = 'euclidean',
         clustering_distance_rows = 'euclidean',
         clustering_method = 'ward.D',
         annotation_col = meta_subset,
         annotation_names_row = F,
         annotation_names_col = T,
         fontsize_row = 10,          # row label font size
         fontsize_col = 7,          # column label font size
         angle_col = 45, # sample names at an angle
         show_colnames = T, show_rownames = F, scale = "row")
dev.off()


# Upset plots
tpm_cow_grouped <- tpm %>%
  mutate(
    `900` = apply(select(., starts_with("900")), 1, max),
    `908` = apply(select(., starts_with("908")), 1, max),
    `914` = apply(select(., starts_with("914")), 1, max),
    `924` = apply(select(., starts_with("924")), 1, max),
    `929` = apply(select(., starts_with("929")), 1, max)
  ) %>%
  select(cluster, `900`, `908`, `914`, `924`, `929`)  # Keep only the cluster and grouped columns

binary_data <- tpm_cow_grouped %>%
  column_to_rownames("cluster") %>%
  mutate(across(everything(), ~ ifelse(. > 10000, 1, 0)))  # Convert values > 10000 to 1, others to 0

png("/Users/goldriev/Google Drive/My Drive/edinburgh/vsg/cobalt_2025/figs/Fig.3/Fig.3.21.png", units="in", width=8, height=6, res=300)
upset(binary_data, sets = colnames(binary_data), order.by = "freq")
dev.off()

#extract cluster names 
set_names <- colnames(binary_data)
all_combos <- unlist(
  lapply(1:length(set_names), function(n) combn(set_names, n, simplify = FALSE)),
  recursive = FALSE
)

intersection_clusters <- set_names(all_combos, sapply(all_combos, paste, collapse = "_")) %>%
  map(~ get_clusters(.x, binary_data))

max_len <- max(sapply(intersection_clusters, length))

intersection_matrix <- sapply(intersection_clusters, function(x) {
  length(x) <- max_len
  x
})

intersection_df <- as.data.frame(intersection_matrix, stringsAsFactors = FALSE)
write.csv(intersection_df, "/Users/goldriev/Google Drive/My Drive/edinburgh/vsg/cobalt_2025/figs/Fig.3/cow_intersections.csv", row.names = FALSE, na = "")

tpm_mouse_grouped <- tpm %>%
  mutate(
    `BRCA2` = apply(select(., starts_with("B")), 1, max),
    `RAD51` = apply(select(., starts_with("R")), 1, max),
    `WT` = apply(select(., starts_with("W")), 1, max),
  ) %>%
  select(cluster, `WT`, `BRCA2`, `RAD51`)  # Keep only the cluster and grouped columns

binary_data <- tpm_mouse_grouped %>%
  column_to_rownames("cluster") %>%
  mutate(across(everything(), ~ ifelse(. > 10000, 1, 0)))  # Convert values > 10000 to 1, others to 0

png("/Users/goldriev/Google Drive/My Drive/edinburgh/vsg/cobalt_2025/figs/Fig.3/Fig.3.22.png", units="in", width=8, height=6, res=300)
upset(binary_data, sets = colnames(binary_data), order.by = "freq")
dev.off()

#extract cluster names 
set_names <- colnames(binary_data)
all_combos <- unlist(
  lapply(1:length(set_names), function(n) combn(set_names, n, simplify = FALSE)),
  recursive = FALSE
)

intersection_clusters <- set_names(all_combos, sapply(all_combos, paste, collapse = "_")) %>%
  map(~ get_clusters(.x, binary_data))

max_len <- max(sapply(intersection_clusters, length))

intersection_matrix <- sapply(intersection_clusters, function(x) {
  length(x) <- max_len
  x
})

intersection_df <- as.data.frame(intersection_matrix, stringsAsFactors = FALSE)
write.csv(intersection_df, "/Users/goldriev/Google Drive/My Drive/edinburgh/vsg/cobalt_2025/figs/Fig.3/mouse_intersections.csv", row.names = FALSE, na = "")

#All
tpm_grouped <- tpm %>%
  mutate(
    `900` = apply(select(., starts_with("900")), 1, max),
    `908` = apply(select(., starts_with("908")), 1, max),
    `914` = apply(select(., starts_with("914")), 1, max),
    `924` = apply(select(., starts_with("924")), 1, max),
    `929` = apply(select(., starts_with("929")), 1, max),
    `BRCA2` = apply(select(., starts_with("B")), 1, max),
    `RAD51` = apply(select(., starts_with("R")), 1, max),
    `WT` = apply(select(., starts_with("W")), 1, max),
  ) %>%
  select(cluster, `900`, `908`, `914`, `924`, `929`, `WT`, `BRCA2`, `RAD51`)  # Keep only the cluster and grouped columns

binary_data <- tpm_grouped %>%
  column_to_rownames("cluster") %>%
  mutate(across(everything(), ~ ifelse(. > 10000, 1, 0)))  

png("/Users/goldriev/Google Drive/My Drive/edinburgh/vsg/cobalt_2025/figs/Fig.3/Fig.3.23.png", units="in", width=8, height=6, res=300)
upset(binary_data, sets = colnames(binary_data), order.by = "freq")
dev.off()

#extract cluster names 
set_names <- colnames(binary_data)
all_combos <- unlist(
  lapply(1:length(set_names), function(n) combn(set_names, n, simplify = FALSE)),
  recursive = FALSE
)

intersection_clusters <- set_names(all_combos, sapply(all_combos, paste, collapse = "_")) %>%
  map(~ get_clusters(.x, binary_data))

max_len <- max(sapply(intersection_clusters, length))

intersection_matrix <- sapply(intersection_clusters, function(x) {
  length(x) <- max_len
  x
})

intersection_df <- as.data.frame(intersection_matrix, stringsAsFactors = FALSE)
write.csv(intersection_df, "/Users/goldriev/Google Drive/My Drive/edinburgh/vsg/cobalt_2025/figs/Fig.3/all_intersections.csv", row.names = FALSE, na = "")

#Split by DPI for upset

dpi_10_samples <- c("900_d4", "908_d4", "914_d4", "924_d4", "929_d4", "908_5", "914_5", "908_6", "914_6", "900_10", "908_10", "924_10", "929_10")
dpi_20_samples <- c("908_d11", "924_d11", "900_13", "914_13", "929_13", "900_14", "914_14", "914_d17", "924_17", "924_18", "924_19", "900_20", "908_20", "914_20", "924_20", "929_20")
dpi_30_samples <- c("900_21", "929_21", "900_24", "908_24", "914_24", "924_24", "929_24", "900_26", "908_26", "914_26", "924_26", "929_26", "900_28", "908_28", "914_28", "924_28")
dpi_40_samples <- c("900_31", "908_31", "914_31", "924_31", "929_31", "900_d33", "914_33", "924_33", "929_d33", "908_35", "914_35", "924_35", "924_38", "908_40", "924_40")
dpi_50_samples <- c("900_42", "908_42", "914_42", "908_d45", "914_d45", "924_d45", "929_d45", "900_d47")
dpi_60_samples <- c("908_52", "914_52", "924_52", "929_52", "900_54", "924_54")
dpi_70_samples <- c("908_61", "914_61", "914_d68")
wt_early_samples <- c("WT_E1", "WT_E2", "WT_E3", "WT_E4", "WT_E5", "WT_E6", "WT_E7", "WT_E8")
wt_late_samples <- c("WT_L1", "WT_L2", "WT_L3", "WT_L4", "WT_L5", "WT_L6", "WT_L7", "WT_L8")
brca_early_samples <- c("B1_2", "B1_4", "B1_6", "B1_8")
brca_late_samples <- c("B1_1", "B1_3", "B1_5", "B1_7", "B1_9")
rad_early_samples <- c("R2_4", "R2_5", "R2_6", "R2_7")
rad_late_samples <- c("R2_1", "R2_2", "R2_3", "R2_8", "R2_9", "R2_2REP", "R2_3REP")

tpm_grouped <- tpm %>%
  mutate(
    `dpi_10` = apply(select(., all_of(dpi_10_samples)), 1, max),
    `dpi_20` = apply(select(., all_of(dpi_20_samples)), 1, max),
    `dpi_30` = apply(select(., all_of(dpi_30_samples)), 1, max),
    `dpi_40` = apply(select(., all_of(dpi_40_samples)), 1, max),
    `dpi_50` = apply(select(., all_of(dpi_50_samples)), 1, max),
    `dpi_60` = apply(select(., all_of(dpi_60_samples)), 1, max),
    `dpi_70` = apply(select(., all_of(dpi_70_samples)), 1, max),
    `wt_early` = apply(select(., all_of(wt_early_samples)), 1, max),
    `wt_late` = apply(select(., all_of(wt_late_samples)), 1, max),
    `brca_early` = apply(select(., all_of(brca_early_samples)), 1, max),
    `brca_late` = apply(select(., all_of(brca_late_samples)), 1, max),
    `rad_early` = apply(select(., all_of(rad_early_samples)), 1, max),
    `rad_late` = apply(select(., all_of(rad_late_samples)), 1, max),
  ) %>%
  select(cluster, `dpi_10`, `dpi_20`, `dpi_30`, `dpi_40`, `dpi_50`, `dpi_60`, `dpi_70`, `wt_early`, `wt_late`,`brca_early`, `brca_late`, `rad_early`, `rad_late` )  # Keep only the cluster and grouped columns

binary_data <- tpm_grouped %>%
  column_to_rownames("cluster") %>%
  mutate(across(everything(), ~ ifelse(. > 10000, 1, 0)))  

png("/Users/goldriev/Google Drive/My Drive/edinburgh/vsg/cobalt_2025/figs/Fig.3/Fig.3.24.png", units="in", width=8, height=6, res=300)
upset(binary_data, sets = colnames(binary_data), order.by = "freq")
dev.off()

#extract cluster names 
set_names <- colnames(binary_data)
all_combos <- unlist(
  lapply(1:length(set_names), function(n) combn(set_names, n, simplify = FALSE)),
  recursive = FALSE
)

intersection_clusters <- set_names(all_combos, sapply(all_combos, paste, collapse = "_")) %>%
  map(~ get_clusters(.x, binary_data))

max_len <- max(sapply(intersection_clusters, length))

intersection_matrix <- sapply(intersection_clusters, function(x) {
  length(x) <- max_len
  x
})

intersection_df <- as.data.frame(intersection_matrix, stringsAsFactors = FALSE)
write.csv(intersection_df, "/Users/goldriev/Google Drive/My Drive/edinburgh/vsg/cobalt_2025/figs/Fig.3/all_dpi_intersections.csv", row.names = FALSE, na = "")

#Split by DPI for upset WT

tpm_grouped <- tpm %>%
  mutate(
    `dpi_10` = apply(select(., all_of(dpi_10_samples)), 1, max),
    `dpi_20` = apply(select(., all_of(dpi_20_samples)), 1, max),
    `dpi_30` = apply(select(., all_of(dpi_30_samples)), 1, max),
    `dpi_40` = apply(select(., all_of(dpi_40_samples)), 1, max),
    `dpi_50` = apply(select(., all_of(dpi_50_samples)), 1, max),
    `dpi_60` = apply(select(., all_of(dpi_60_samples)), 1, max),
    `dpi_70` = apply(select(., all_of(dpi_70_samples)), 1, max),
    `wt_early` = apply(select(., all_of(wt_early_samples)), 1, max),
    `wt_late` = apply(select(., all_of(wt_late_samples)), 1, max),
  ) %>%
  select(cluster, `dpi_10`, `dpi_20`, `dpi_30`, `dpi_40`, `dpi_50`, `dpi_60`, `dpi_70`, `wt_early`, `wt_late`)  # Keep only the cluster and grouped columns

binary_data <- tpm_grouped %>%
  column_to_rownames("cluster") %>%
  mutate(across(everything(), ~ ifelse(. > 10000, 1, 0)))  

png("/Users/goldriev/Google Drive/My Drive/edinburgh/vsg/cobalt_2025/figs/Fig.3/Fig.3.25.png", units="in", width=8, height=6, res=300)
upset(binary_data, sets = colnames(binary_data), order.by = "freq")
dev.off()

#extract cluster names 
set_names <- colnames(binary_data)
all_combos <- unlist(
  lapply(1:length(set_names), function(n) combn(set_names, n, simplify = FALSE)),
  recursive = FALSE
)

intersection_clusters <- set_names(all_combos, sapply(all_combos, paste, collapse = "_")) %>%
  map(~ get_clusters(.x, binary_data))

max_len <- max(sapply(intersection_clusters, length))

intersection_matrix <- sapply(intersection_clusters, function(x) {
  length(x) <- max_len
  x
})

intersection_df <- as.data.frame(intersection_matrix, stringsAsFactors = FALSE)
write.csv(intersection_df, "/Users/goldriev/Google Drive/My Drive/edinburgh/vsg/cobalt_2025/figs/Fig.3/wt_dpi_intersections.csv", row.names = FALSE, na = "")

#Actual venn diagram

tpm_grouped <- tpm %>%
  mutate(
    `cow` = apply(select(., starts_with("9")), 1, max),
    `WT` = apply(select(., starts_with("W")), 1, max),
    `BRCA2` = apply(select(., starts_with("B")), 1, max),
    `RAD51` = apply(select(., starts_with("R")), 1, max)
  ) %>%
  select(cluster, `cow`, `WT`, `BRCA2`, `RAD51`,) 

binary_data <- tpm_grouped %>%
  column_to_rownames("cluster") %>%
  mutate(across(everything(), ~ ifelse(. > 10000, 1, 0))) 

venn_sets <- lapply(colnames(binary_data), function(col) {
  rownames(binary_data[binary_data[[col]] == 1, , drop = FALSE])
})

names(venn_sets) <- colnames(binary_data)

png("/Users/goldriev/Google Drive/My Drive/edinburgh/vsg/cobalt_2025/figs/Fig.3/Fig.3.26.png", units="in", width=6, height=6, res=300)
ggVennDiagram(venn_sets, label_alpha = 0, label = "count") +
  scale_fill_gradient(low = "white", high = "dodgerblue") +
  theme(legend.position = "none") +
  ggtitle("")
dev.off()

#cow vs mouse
tpm_grouped <- tpm %>%
  mutate(
    `Cow` = apply(select(., starts_with("9")), 1, max),
    `Mouse` = apply(select(., starts_with("W")), 1, max),
  ) %>%
  select(cluster, `Cow`, `Mouse`) 

binary_data <- tpm_grouped %>%
  column_to_rownames("cluster") %>%
  mutate(across(everything(), ~ ifelse(. > 10000, 1, 0))) 

venn_sets <- lapply(colnames(binary_data), function(col) {
  rownames(binary_data[binary_data[[col]] == 1, , drop = FALSE])
})

names(venn_sets) <- colnames(binary_data)

png("/Users/goldriev/Google Drive/My Drive/edinburgh/vsg/cobalt_2025/figs/Fig.3/Fig.3.27.png", units="in", width=8, height=6, res=300)
ggVennDiagram(venn_sets, label_alpha = 0, label = "count") +
  scale_fill_gradient(low = "white", high = "dodgerblue") +
  theme(legend.position = "none") +
  ggtitle("")
dev.off()

tpm_grouped <- tpm %>%
  mutate(
    `900` = apply(select(., starts_with("900")), 1, max),
    `908` = apply(select(., starts_with("908")), 1, max),
    `914` = apply(select(., starts_with("914")), 1, max),
    `924` = apply(select(., starts_with("924")), 1, max),
    `929` = apply(select(., starts_with("929")), 1, max)
  ) %>%
  select(cluster, `900`, `908`, `914`, `924`, `929`) 

binary_data <- tpm_grouped %>%
  column_to_rownames("cluster") %>%
  mutate(across(everything(), ~ ifelse(. > 10000, 1, 0))) 

venn_sets <- lapply(colnames(binary_data), function(col) {
  rownames(binary_data[binary_data[[col]] == 1, , drop = FALSE])
})

names(venn_sets) <- colnames(binary_data)

png("/Users/goldriev/Google Drive/My Drive/edinburgh/vsg/cobalt_2025/figs/Fig.3/Fig.3.28.png", units="in", width=6, height=6, res=300)
ggVennDiagram(venn_sets, label_alpha = 0, label = "count") +
  scale_fill_gradient(low = "white", high = "dodgerblue") +
  theme(legend.position = "none") +
  ggtitle("")
dev.off()
