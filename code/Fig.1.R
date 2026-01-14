# Load packages
pacman::p_load(dplyr, RColorBrewer, reshape, ggpubr, scales, tibble, pheatmap, ggplotify, magick, pdftools)

# Fig.2

mapping <- read.csv("~/Google Drive/My Drive/edinburgh/vsg/vsgseq2_analysis/WOR/analysis/bench/mapping.csv")
bench <- read.csv("~/Google Drive/My Drive/edinburgh/vsg/vsgseq2_analysis/WOR/analysis/bench/cd-94-94/comparison.csv")

bench <- bench %>%
  group_by(file, tool) %>%
  mutate(total_expression = sum(expression),
         expression_pct = (expression / total_expression) * 100) %>%
  ungroup()

count_data <- bench %>%
  filter(expression_pct >= 0.01) %>%
  group_by(file, tool, library) %>%
  summarise(count = n(), .groups = 'drop')

count_data <- subset(count_data, file != "Expected")

count_data %>%
  group_by(tool, library) %>%
  summarise(mean_count = mean(count), .groups = 'drop')

vsg_order <- c("Tb427VSG-3", "Tb427VSG-2", "Tb427VSG-17", "Tb427VSG-417", "Tb427VSG-629", "Tb427VSG-11", "Tb427VSG-9", "Other")
bench$cross_id <- factor(bench$cross_id, levels = vsg_order)

comparisons <- list(
  c("VSGSeq_individual", "VSGSeq_combined"),
  c("VSGSeq_individual", "vsgseq2"),
  c("VSGSeq_combined", "vsgseq2")
)

fig.2a <- ggboxplot(count_data, x = "tool", y = "count",
                    xlab = "", 
                    ylab = "VSG count", fill= "tool", colour = "tool") +
  geom_hline(yintercept = 7, linetype = "dashed", color = "black") +
  facet_wrap(~ library, scales = "free_x") +
  stat_compare_means(comparisons = comparisons, method = "t.test", label = "p.signif")  +
  scale_fill_manual(values = brewer.pal(4, "Set3")) +
  theme(axis.text.x = element_blank()) + 
  theme(legend.position="none")

mapping$percent_aligned <- (mapping$aligned / mapping$total) * 100

mapping %>%
  group_by(tool, library, method) %>%
  summarise(mean_count = mean(percent_aligned), .groups = 'drop')

fig.2b <- ggboxplot(mapping, x = "tool", y = "percent_aligned", fill = "tool", legend = "right") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  facet_wrap(~ library, scales = "free_x") +
  ylab("Mapped reads (%)") +
  xlab("") +
  stat_compare_means(comparisons = comparisons, method = "t.test", label = "p.signif")  +
  scale_fill_manual(values = brewer.pal(4,"Set3") ) + 
  theme(axis.text.x = element_blank()) +
  theme(legend.position="none")


coul <- brewer.pal(7,"Set2") 
coul <- append(coul, "black")

lib_a <- subset(bench, library == "A")
lib_b <- subset(bench, library == "B")

fig.2c <- ggbarplot(lib_a, x = "file", y = "expression_pct",
          add.params = list(size = 2), color = "cross_id", fill = "cross_id") +
  facet_wrap(~ tool, scales = "free_x") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  scale_fill_manual(values = coul) +
  scale_color_manual(values = coul) +
  ylab("Library A population expression (%)") +
  xlab("") +
  theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust=1)) +
  geom_bar(
    data = subset(lib_a, file == "Expected_A"),
    aes(x = file, y = expression_pct),
    stat = "identity",
    fill = NA, color = "black", size = 1, width = 0.8) 

fig.2d <- ggbarplot(lib_b, x = "file", y = "expression_pct",
                    add.params = list(size = 2), color = "cross_id", fill = "cross_id") +
  facet_wrap(~ tool) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  scale_fill_manual(values = coul) +
  scale_color_manual(values = coul) +
  ylab("Library B population expression (%)") +
  xlab("") +
  theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust=1)) +
  geom_bar(
  data = subset(lib_b, file == "Expected_B"),
  aes(x = file, y = expression_pct),
  stat = "identity",
  fill = NA, color = "black", size = 1, width = 0.8) 

png("~/Google Drive/My Drive/edinburgh/vsg/vsgseq2/method_paper/figures/Fig.2.png", units="in", width=10, height=12, res=300)
row1 <- ggarrange(
  fig.2a + theme(legend.position = "none"),
  fig.2b + theme(legend.position = "right"),
  ncol = 2,
  labels = c("A", "B"),
  common.legend = FALSE,
  align = "hv"
)
ggarrange(
  row1,
  fig.2c + theme(legend.position = "right"),
  fig.2d + theme(legend.position = "right"),
  ncol = 1,
  heights = c(1, 1.2, 1.2),
  labels = c("", "C", "D"),
  font.label = list(
    size = 14, 
    color = "black", 
    face = "bold"
  ),
  common.legend = FALSE,
  align = "hv"
)
dev.off()

#Fig. 3

mapping <- read.csv("~/Google Drive/My Drive/edinburgh/vsg/vsgseq2_analysis/WOR/analysis/mouse/mapping.csv")
expression <- read.csv("~/Google Drive/My Drive/edinburgh/vsg/vsgseq2_analysis/WOR/analysis/mouse/combined_expression.csv")
speed <- read.csv("~/Google Drive/My Drive/edinburgh/vsg/vsgseq2_analysis/WOR/analysis/mouse/speed_test.csv")

mapping$percent_aligned <- (mapping$aligned / mapping$total) * 100

mapping %>%
  group_by(tool) %>%
  summarise(mean_count = mean(percent_aligned), .groups = 'drop')

fig.3b <- ggboxplot(mapping, x = "tool", y = "percent_aligned", fill = "tool", legend = "right") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  ylab("Mapped reads (%)") +
  xlab("") +
  stat_compare_means(comparisons = comparisons, method = "t.test", label = "p.signif")  +
  scale_fill_manual(values = brewer.pal(4,"Set3") ) +
  theme(legend.position="none")

mouse <- expression %>%
  group_by(file, tool) %>%
  mutate(total_expression = sum(expression),
         expression_pct = (expression / total_expression) * 100) %>%
  ungroup()

count_data <- mouse %>%
  filter(expression_pct >= 0.01) %>%
  group_by(file, tool) %>%
  summarise(count = n(), .groups = 'drop')

count_data %>%
  group_by(tool) %>%
  summarise(mean_count = mean(count), .groups = 'drop')

fig.3a <- ggboxplot(count_data, x = "tool", y = "count",
          xlab = "", 
          ylab = "VSG count", fill= "tool", colour = "tool") +
  stat_compare_means(comparisons = comparisons, method = "t.test", label = "p.signif")  +
  scale_fill_manual(values = brewer.pal(4, "Set3")) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(legend.position="none")

order <- c("VSGSeq", "vsgseq2") 
speed$tool <- factor(speed$tool, levels = order)

fig.4a <- ggboxplot(speed, x = "tool", y = "CO2e", fill = "tool", legend = "right") +
  ylab("CO2e (mg)") +
  xlab("") +
  stat_compare_means(method = "t.test", label = "p.signif") +
  scale_fill_manual(values = brewer.pal(4, "Set3")) 

fig.4b <- ggboxplot(speed, x = "tool", y = "total", fill = "tool", legend = "right") +
  ylab("Runtime (seconds)") +
  xlab("") +
  stat_compare_means(method = "t.test", label = "p.signif") +
  scale_fill_manual(values = brewer.pal(4, "Set3")) 

mouse$cluster <- as.character(mouse$cluster)

mouse <- mouse %>%
  mutate(Names = if_else(expression_pct > 10, cluster, "Other"))

coul <- brewer.pal(8,"Set3") 
coul <- colorRampPalette(coul)(-1 + length(unique(mouse$Names)))
coul <- append(coul, "black")

fig.3c <- ggbarplot(mouse, x = "tool", y = "expression_pct",
          add.params = list(size = 2), color = "Names", fill = "Names") +
  facet_wrap(mouse ~ day) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  scale_fill_manual(values = coul) +
  scale_color_manual(values = coul) +
  ylab("Population expression (%)") +
  xlab("") +
  theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust=1)) +
  theme(legend.position="none")

png("~/Google Drive/My Drive/edinburgh/vsg/vsgseq2/method_paper/figures/Fig.3.png", units="in", width=8, height=13, res=300)
ggarrange(
  ggarrange(fig.3a, fig.3b, ncol = 2, labels = c("A", "B"), font.label = list(size = 14, face = "bold")),
  fig.3c,
  ncol = 1, 
  heights = c(1, 3),  # fig.3c is 3x taller than the top row
  labels = c(NA, "C"),
  font.label = list(size = 14, face = "bold"),
  align = "hv",
  common.legend = TRUE,
  legend = "none"
)

dev.off()

png("~/Google Drive/My Drive/edinburgh/vsg/vsgseq2/method_paper/figures/Fig.4.png", units="in", width=5, height=3, res=300)
ggarrange(
  fig.4b, fig.4a,
  heights = c(1, 1), 
  common.legend = T,
  legend = "none",
  align = "hv", 
  labels = "auto",
  font.label = list(
    size = 14, 
    color = "black", 
    face = "bold", 
    family = NULL
  )
)
dev.off()

#Fig. 5

data <- read.table('~/Google Drive/My Drive/edinburgh/vsg/vsgseq2/method_paper/mosaics/mosaics.txt', header = TRUE, sep = "\t")

# Ensure unique row names by creating a unique identifier for internal use
data <- data %>%
  mutate(unique_id = paste(Donor, row_number(), sep = "_"))

# Prepare the data for the heatmap
# Select relevant columns and convert to matrix
heatmap_data <- data %>%
  select(unique_id, VSGSeq, vsgseq2) %>%
  column_to_rownames(var = "unique_id") %>%
  as.matrix()

# Ensure the cluster column is a factor
data$Expressed <- factor(data$Expressed, levels = unique(data$Expressed))

# Create annotation dataframe
annotation <- data.frame(Expressed = data$Expressed)
rownames(annotation) <- data$unique_id

# Create a black and white color palette
bw_palette <- colorRampPalette(c("gray90", "#BEBADA"))(100)

annotation_colors <- list(
  Expressed = c(
    "EF174487" = "#8DD3C7", "EF174488" = "#BEBADA", "EF174489" = "#FB8072",
    "EF174490" = "#80B1D3", "EF174491" = "#FDB462", "EF174492" = "#B3DE69",
    "EF174493" = "#FCCDE5", "EF560198" = "#D9D9D9", "EF560199" = "#BC80BD",
    "EF560200" = "#CCEBC5", "EF560201" = "#FFED6F", "EF560203" = "#1F78B4",
    "EF560204" = "#33A02C", "EF568372" = "#E31A1C"
  )
)

# Generate the heatmap
fig.5a <- as.ggplot(pheatmap(t(heatmap_data), cluster_rows = FALSE, cluster_cols = FALSE, 
                             color = bw_palette,
                             annotation_colors = annotation_colors,
                             show_rownames = TRUE, show_colnames = TRUE, 
                             annotation_col = annotation,
                             labels_col = data$Donor,
                             gaps_col = c(4, 8, 11, 15, 17, 20, 25, 27, 29, 31, 33, 35, 37)))

fig.5a <- fig.5a + theme(plot.margin = unit(c(0, 0, 0, 1), "cm"))

png_to_ggplot <- function(path, density = 600, pages = 1) {
  image <- image_read(path, density = density)
  raster <- as.raster(image)
  ggplot() +
    annotation_raster(raster, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf) +
    theme_void()
}

# Paths to the PDF files
path <- list(
  EF174492 = "~/Google Drive/My Drive/edinburgh/vsg/vsgseq2/method_paper/figures/Fig.5/EF174492.png",
  EF560198 = "~/Google Drive/My Drive/edinburgh/vsg/vsgseq2/method_paper/figures/Fig.5/EF560198.png",
  EF560204 = "~/Google Drive/My Drive/edinburgh/vsg/vsgseq2/method_paper/figures/Fig.5/EF560204.png"
)

# Create ggplot objects
fig.5b <- png_to_ggplot(path$EF174492)
fig.5c <- png_to_ggplot(path$EF560198)
fig.5d <- png_to_ggplot(path$EF560204)

png("~/Google Drive/My Drive/edinburgh/vsg/vsgseq2/method_paper/figures/Fig.5.png", units="in", width=10, height=15, res=300)
ggarrange(
  fig.5a,               
  ggarrange(fig.5b, fig.5c, ncol = 2, labels = c("b", "c")), 
  # Third row with dot plot and blank space
  ggarrange(fig.5d, ncol = 2, labels = c("d", "")), 
  nrow = 3, 
  labels = c("a", "", ""), 
  common.legend = TRUE, 
  legend = "right", 
  align = "hv", 
  font.label = list(size = 14, color = "black", face = "bold", family = NULL)
)
dev.off()


# Create box plot

comparisons <- list(
  c("VSGSeq_single", "VSGSeq_combined"),
  c("VSGSeq_single", "vsgseq2_single"),
  c("VSGSeq_single", "vsgseq2_combined"),
  c("VSGSeq_combined", "vsgseq2_single"),
  c("VSGSeq_combined", "vsgseq2_combined"),
  c("vsgseq2_single", "vsgseq2_combined")
)

mapping$percent_aligned <- (mapping$aligned / mapping$total) * 100

a <- ggboxplot(mapping, x = "tool", y = "percent_aligned", fill = "tool", legend = "right") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  ylab("Mapped reads (%)") +
  stat_compare_means(comparisons = comparisons, method = "t.test", label = "p.signif") +
  scale_fill_manual(values = brewer.pal(4,"Set3") ) 

b <- ggbarplot(mapping, x = "tool", y = "percent_aligned", fill = "tool", legend = "right") +
  theme(axis.text.x = element_blank(), # Hide x-axis text
        axis.title.x = element_blank(), # Hide x-axis title
        axis.ticks.x = element_blank()) + # Hide x-axis ticks
  facet_wrap(~ sample, scales = 'free_x') +
  ylab("Mapped reads (%)") +
  stat_compare_means(comparisons = comparisons, method = "t.test", label = "p.signif") +
  scale_fill_manual(values = brewer.pal(4, "Set3"))

comparisons <- list(
  c("VSGSeq", "vsgseq2")
)

speed %>%
  group_by(tool) %>%
  summarise(
    mean_total = mean(total, na.rm = TRUE),
    sum_total = sum(total, na.rm = TRUE),
    mean_CO2e = mean(CO2e, na.rm = TRUE),
    sum_CO2e = sum(CO2e, na.rm = TRUE)
  )

a <- a + theme(legend.position = "none")
b <- b + theme(legend.position = "none")
c <- c + theme(legend.position = "none")

png("~/Google Drive/My Drive/edinburgh/vsg/vsgseq2/method_paper/figures/Fig.3.png", units="in", width=10, height=12, res=300)
ggarrange(
  b,          
  ggarrange(a, c, ncol = 2, labels = c("b", "c")), 
  nrow = 2, 
  labels = "a", 
  common.legend = TRUE, 
  legend = "top", 
  align = "hv", 
  heights = c(2, 1),  # Adjust heights as needed
  font.label = list(size = 14, color = "black", face = "bold", family = NULL)
)
dev.off()

