suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(pheatmap))
suppressPackageStartupMessages(library(RColorBrewer))

# Perform log transformation on our count data
dds <- readRDS("dds.rds")
rld <- rlog(dds)

# Generate a PCA plot with DESeq2's plotPCA function
plotPCA(rld, intgroup = "condition")

# Calculate distances between samples in our log-transformed data
sample_dists <- dist(t(assay(rld)))

# Convert the output to a matrix
sample_dist_matrix <- as.matrix(sample_dists)

# Remove the column names of our matrix
colnames(sample_dist_matrix) <- NULL

# Set the colour palette for our heatmap
colours <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)

### Generate a heatmap using the pheatmap package for infected Reps and control Reps
labels_row <- c("Control Rep. 1", "Control Rep. 2", "Infected Rep. 1", "Infected Rep. 2")

pheatmap(sample_dist_matrix,
         clustering_distance_rows = sample_dists,
         clustering_distance_cols = sample_dists,
         col = colours,
         labels_row = labels_row,
         show_colnames = FALSE)


# Get the PCA data as a dataframe with returnData = TRUE
pcaData <- plotPCA(rld, intgroup=c("condition"), returnData = TRUE)
glimpse(pcaData)
# Using the attr() function, we will extract the percent variation explained of each axis of the pcaData object
percentVar <- round(100 * attr(pcaData, "percentVar"))
glimpse(percentVar)
  
### Dot plot for control and reps
pcaData %>% 
  ggplot(aes(x = PC1, y = PC2, colour = condition)) +
  geom_point(size = 3) +
  scale_y_continuous(limits = c(-10,10)) +
  scale_x_continuous(limits = c(-50,50)) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) + 
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  theme_test() +
  theme(
    axis.text = element_text(colour = "black", size = 11),
    axis.title = element_text(face = "bold"),
    legend.title = element_text(face = "bold"),
    panel.border = element_rect(linewidth = 2)) +
  labs(colour = "Condition")

pcaData$condition <- factor(
  pcaData$condition,
  levels = c("unstimulated_control", "TB_infected"),
  labels = c("Control", "Infected")
)

pcaData %>%
  ggplot(aes(x = PC1, y = PC2, colour = condition)) +
  geom_point(size = 3) +
  scale_y_continuous(limits = c(-10, 10)) +
  scale_x_continuous(limits = c(-50, 50)) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  theme_test() +
  theme(
    axis.text = element_text(colour = "black", size = 11),
    axis.title = element_text(face = "bold"),
    legend.title = element_text(face = "bold"),
    panel.border = element_rect(linewidth = 2)
  ) +
  labs(colour = "Condition")


# Summary of diferentially expressed genes
dat <- read_csv("res_no_NA_results.csv")
dat


### Volcano plot of differentially expressed genes fitting all data
# Cap -log10(padj) at 300 if padj == 0 for plotting
dat <- dat %>%
  mutate(log_padj = ifelse(padj == 0, 300, -log10(padj))) # cap infinite values

# Assign categories for legend coloring
labelled_dat <- dat %>%
  mutate(Genes = case_when(
    padj < 0.05 & log2FoldChange > 1 ~ "Upregulated",
    padj < 0.05 & log2FoldChange < -1 ~ "Downregulated",
    TRUE ~ "Not Significant"
  ))

# Custom color palette (optional)
volcano_colors <- c(
  "Downregulated" = "red",
  "Not Significant" = "green",
  "Upregulated" = "blue"
)

labelled_dat %>%
  ggplot(aes(x = log2FoldChange, y = log_padj, color = Genes)) +
  geom_point(size = 1) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", linewidth = 0.8) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", linewidth = 0.8) +
  scale_color_manual(values = volcano_colors) +
  labs(
    x = expression(log[2]("Fold Change")),
    y = expression(-log[10]("p-adjusted")),
    color = "Genes"
  ) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", size = 1.5)
  ) +
  coord_cartesian(xlim = c(-10, 10))


# Other plots 

# Take the 'gene_id' column in order and store it as a character vector
order <- top10Genes %>% 
  dplyr::select(gene_id) %>% 
  pull()

# Plot again but arrange the bars based on the order of the character vector
top10Genes %>% 
  ggplot(aes(x = gene_id, y = log2FoldChange)) +
  geom_col() +
  scale_x_discrete(limits = order) # use the scale_x_discrete() function with the limits= argument to order the plot
top10Genes %>% 
  ggplot(aes(x = gene_id, y = log2FoldChange)) +
  geom_col() +
  scale_x_discrete(limits = order) +
  coord_flip()
top10Genes %>% 
  ggplot(aes(x = gene_id, y = log2FoldChange)) +
  geom_errorbar(aes(ymin= log2FoldChange - lfcSE, ymax = log2FoldChange + lfcSE), width = 0.4) +
  geom_col() +
  scale_x_discrete(limits = order) +
  coord_flip()
# Bind our two top 10 genes dataframes using the bind_rows() column
joined_top10 <- bind_rows(top10Genes, bot10Genes) %>%
  mutate(Genes = if_else(log2FoldChange > 0, "UP", "DOWN")) %>%
  arrange(log2FoldChange)

# To organize our bar plot, we can also extract the order of our organized genes
order_joined <- joined_top10 %>% 
  dplyr::select(gene_id) %>% 
  pull()

joined_top10 %>% 
  ggplot(aes(x = gene_id, y = log2FoldChange, fill = Genes)) +
  geom_errorbar(aes(ymin= log2FoldChange - lfcSE, ymax =log2FoldChange + lfcSE), width = 0.4) +
  geom_col(color = "black", width = 0.8) +
  scale_x_discrete(limits = order_joined) +
  coord_flip()+ 
  labs(y = "log2(FC)", x= "Gene ID", fill = "Regulation") +
  geom_hline(yintercept = 0) +
  theme_light() +
  scale_y_continuous(limits = c(-17,17), breaks = seq(-15, 15, 5), expand = c(0,0)) +
  scale_fill_manual(values = c("red", "forestgreen")) +
  theme(
    panel.border = element_rect(color = "black"),
    panel.grid = element_line(colour = "white"),
    axis.text = element_text(colour = "black"))

#heat maps z score

vsd <- assay(vst(dds))

Z <- t(scale(t(vsd)))

# Check the Z-score matrix
head(Z)

# Remove rows with NA values - if not we will get errors!
Z_no_NA <- na.omit(Z)

# Now lets plot using pheatmap() 
pheatmap(Z_no_NA, 
         main = "Gene Expression in Control vs. TB Infected Mice")

# Extract the row indices (akak the location in the matrix) of the rows with the highest variation
topVarGenes <- rowVars(vsd) %>% 
  order(decreasing = TRUE) %>% 
  head(20) 

# Subset the Z-score matrix for only the row we're interested in (stored in topVarGenes)
Z_topGenes <- Z[topVarGenes,] 

pheatmap(Z_topGenes,
         main = "Top 20 Genes")

