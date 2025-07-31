setwd("/mnt/d/scTriH2/SEURAT/extracted/1.cellRangerCount/")

library(devtools)
install_github("jokergoo/ComplexHeatmap")
library(ComplexHeatmap)
library(circlize)     # for color functions
library(dplyr)
library(readr)

# Read the CSV file
file_path <- "/mnt/d/scTriH2/SEURAT/extracted/1.cellRangerCount/ranked_norm_gene_counts_high_only.csv"
gene_data <- read_csv(file_path)

# Clean column names by removing 'object_' prefix
colnames(gene_data) <- gsub("^object_", "", colnames(gene_data))

# Set gene names as row names and remove 'Gene' and 'Sum' columns
rownames(gene_data) <- gene_data$Gene
gene_data_matrix <- gene_data[1:(nrow(gene_data)), 2:(ncol(gene_data) - 1)]  # Remove first (Gene) column and last column (total)

# Convert to matrix
gene_matrix <- as.matrix(gene_data_matrix)

# Define color scale
col_fun <- colorRamp2(c(min(gene_matrix), max(gene_matrix)), c("white", "red"))


# Plot heatmap
Heatmap(gene_matrix,
        name = "Gene Count",
        col = col_fun,
        cluster_rows = TRUE,
        cluster_columns = TRUE,
        show_row_names = TRUE,
        show_column_names = TRUE,
        row_names_gp = gpar(fontsize = 6),
        column_names_gp = gpar(fontsize = 10),
        heatmap_legend_param = list(title = "Normalized Count"),
        column_title = "Library",
        row_title = "Bacterial Genes")

## Log transformation with norm_ranked_gene_counts

# Set gene names as row names and remove 'Gene' and 'Sum' columns
rownames(gene_data) <- gene_data$Gene
gene_data_matrix <- log(gene_data[1:120, 2:(ncol(gene_data) - 1)])  # Remove first (Gene) column and last column (total)
# also adjust the number of displayed genes



### Do it also for the absolute gene counts

abs_count_file_path <- "/mnt/d/scTriH2/SEURAT/extracted/1.cellRangerCount/gene_counts_table.csv"
abs_gene_data <- read_csv(abs_count_file_path)

# exclude the 4th column with very low count
abs_gene_data_filtered <- abs_gene_data[, -4]
abs_sample_cols <- 2:(ncol(abs_gene_data_filtered))
abs_gene_data_filtered$Sum <- rowSums(abs_gene_data_filtered[, abs_sample_cols])
abs_gene_data_filtered <- abs_gene_data_filtered %>% arrange(desc(Sum))

write_csv(abs_gene_data_filtered, "/mnt/d/scTriH2/SEURAT/extracted/1.cellRangerCount/gene_counts_high_only.csv")

# Read the CSV file
file_path <- "/mnt/d/scTriH2/SEURAT/extracted/1.cellRangerCount/gene_counts_high_only.csv"
gene_data <- read_csv(file_path)

# Clean column names by removing 'object_' prefix
colnames(gene_data) <- gsub("^object_", "", colnames(gene_data))

# Set gene names as row names and remove 'Gene' and 'Sum' columns
rownames(gene_data) <- gene_data$Gene
gene_data_matrix <- gene_data[1:100, 2:(ncol(gene_data) - 1)]  # Remove first (Gene) column and last column (total)

# Optional: Log transform the data, also adjust the number of displayed genes
gene_data_matrix <- log(gene_data[1:100, 2:(ncol(gene_data) - 1)])  # Remove first (Gene) column and last column (total)


# Convert to matrix
gene_matrix <- as.matrix(gene_data_matrix)
# optional: log-transform gene counts (add 1 to avoid log(0))
gene_matrix_log <- log(gene_matrix + 1)

# Color(without log transformation)
col_fun <- colorRamp2(c(min(gene_matrix), max(gene_matrix)), c("white", "red"))

# Define color scale for log values and color (with log transformation)
min_log <- min(gene_matrix_log, na.rm = TRUE)
max_log <- max(gene_matrix_log, na.rm = TRUE)

col_fun_log <- colorRamp2(c(min_log, max_log), c("white", "red"))

# Plot heatmap
png("gene_count_per_lib_36.png", width = 1200, height = 1000, res = 150)
Heatmap(gene_matrix,
        name = "Gene Count",
        col = col_fun,
        cluster_rows = TRUE,
        cluster_columns = TRUE,
        show_row_names = TRUE,
        show_column_names = TRUE,
        row_names_gp = gpar(fontsize = 6),
        column_names_gp = gpar(fontsize = 10),
        heatmap_legend_param = list(title = "Gene Count per Library"),
        column_title = "Library",
        row_title = "Bacterial Genes")
dev.off()

# Save heatmap as PNG with log-transformed values
png("gene_count_log10_per_lib.png", width = 1200, height = 1000, res = 150)
Heatmap(gene_matrix_log,
        name = "log(Count + 1)",
        col = col_fun_log,
        cluster_rows = TRUE,
        cluster_columns = TRUE,
        show_row_names = TRUE,
        show_column_names = TRUE,
        row_names_gp = gpar(fontsize = 6),
        column_names_gp = gpar(fontsize = 10),
        heatmap_legend_param = list(title = "log(Count + 1)"),
        column_title = "Library",
        row_title = "Bacterial Genes")
dev.off()



## Calculate fold change relative to row mean (normalized, ranked, DGE)
# Read the CSV file
file_path <- "/mnt/d/scTriH2/SEURAT/extracted/1.cellRangerCount/ranked_norm_gene_counts_high_only.csv"
gene_data <- read_csv(file_path)

# Clean column names by removing 'object_' prefix
colnames(gene_data) <- gsub("^object_", "", colnames(gene_data))

# Set gene names as row names and remove 'Gene' and 'Sum' columns
rownames(gene_data) <- gene_data$Gene
gene_data_matrix <- gene_data[1:(nrow(gene_data)), 2:(ncol(gene_data) - 1)]  # Remove first (Gene) column and last column (total)

# Convert to matrix
gene_matrix <- as.matrix(gene_data_matrix)

# calculated expression fold change
row_means <- rowMeans(gene_matrix)
fold_change_matrix <- sweep(gene_matrix, 1, row_means, "/")

# Filter genes with fold change > 1.5 in at least one sample
keep_genes <- apply(fold_change_matrix, 1, function(x) any(x > 1.5))
filtered_matrix <- fold_change_matrix[keep_genes, ]

col_fun <- colorRamp2(c(1, max(filtered_matrix)), c("white", "red"))
# Save heatmap to PNG
png("gene_fold_change_heatmap.png", width = 1200, height = 1000, res = 150)
Heatmap(filtered_matrix,
        name = "Expression Fold Change",
        col = col_fun,
        cluster_rows = TRUE,
        cluster_columns = TRUE,
        show_row_names = TRUE,
        show_column_names = TRUE,
        row_names_gp = gpar(fontsize = 6),
        column_names_gp = gpar(fontsize = 10),
        heatmap_legend_param = list(title = "Fold Change"),
        column_title = "Library",
        row_title = "Bacterial Genes")
dev.off()


##### Gene ranks in percentage

# Read the CSV file
file_path <- "/mnt/d/scTriH2/SEURAT/extracted/1.cellRangerCount/gene_counts_high_only.csv"
abs_gene_data <- read_csv(file_path)

# calculate the relative abundance of the genes
abs_gene_data_perc <- abs_gene_data[, 2:(ncol(abs_gene_data) - 1)]
col_sums <- colSums(abs_gene_data_perc)
abs_gene_data_perc <- sweep(abs_gene_data_perc, 2, col_sums, FUN = "/") * 100
abs_gene_data_perc <- cbind(Gene = abs_gene_data$Gene, abs_gene_data_perc, Sum = abs_gene_data$Sum)

write_csv(abs_gene_data_perc, "/mnt/d/scTriH2/SEURAT/extracted/1.cellRangerCount/gene_perc_high_only.csv")

abs_sample_cols <- 1:(ncol(abs_gene_data_perc))
abs_gene_data_perc <- abs_gene_data_perc[, abs_sample_cols]
abs_gene_data_perc <- abs_gene_data_perc[,2:8]/col_sums
# abs_gene_data_perc$Sum <- abs_gene_data_perc
# abs_gene_data_perc <- abs_gene_data_perc %>% arrange(desc())


library(ggplot2)
library(reshape2)
file_path <- "/mnt/d/scTriH2/SEURAT/extracted/1.cellRangerCount/gene_perc_high_only.csv"
gene_data <- read_csv(file_path)

# Clean column names by removing 'object_' prefix
colnames(gene_data) <- gsub("^object_", "", colnames(gene_data))
# Set gene names as row names 
rownames(gene_data) <- gene_data$Gene
gene_data_matrix <- gene_data[1:200, 2:(ncol(gene_data) - 1)]  # Remove first (Gene) column and last column (total)

# Convert to matrix
gene_matrix <- as.matrix(gene_data_matrix)
rownames(gene_matrix) <- gene_data$Gene[1:100]
# gene_data[1:(nrow(gene_data)), 2:(ncol(gene_data) - 1)]
# this does not work


# Take only the top 100 genes
gene_data_top <- gene_data[1:200, ]
# Reshape from wide to long format
gene_data_long <- melt(gene_data_top[, 1:(ncol(gene_data_top) - 1)], 
                       id.vars = "Gene",
                       variable.name = "Library",
                       value.name = "Percentage")
# compute top 30 genes for the legend
top_genes <- gene_data_long %>%
  group_by(Gene) %>%
  summarise(mean_perc = mean(Percentage, na.rm = TRUE)) %>%
  arrange(desc(mean_perc)) %>%
  slice(1:30) %>%
  pull(Gene)

# Truncate gene names to 20 characters
gene_data_long <- gene_data_long %>%
  mutate(Gene_trunc = ifelse(nchar(Gene) > 20,
                             paste0(substr(Gene, 1, 17), "..."),
                             Gene))


# Draw the stacked bar plot
png("gene_perc_stacked_bar.png", width = 1500, height = 1500, res = 150)
ggplot(gene_data_long, aes(x = Library, y = Percentage, fill = Gene)) +
  geom_bar(stat = "identity", position = "stack") +
  theme_minimal() +
  labs(title = "Relative Gene Percentage per Library",
       x = "Library",
       y = "Relative gene abundance per library (%)",
       fill = "Bacterial Genes") +
  theme(axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 10),
        legend.position = "right",
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8),
        legend.key.size = unit(6,"mm"), 
        panel.background = element_rect(fill="white")) +
  scale_fill_discrete(breaks = c(gene_data_long$Gene[1:30])) + 
  guides(fill = guide_legend(ncol = 1, byrow = TRUE))  # You can also use nrow instead
#ggsave("gene_perc_stacked_bar.png",plot = last_plot(),bg = NULL, units="px", width = 1500, height = 1500, dpi = 1000, device="png")
dev.off()


