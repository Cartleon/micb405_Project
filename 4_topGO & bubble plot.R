suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(pheatmap)) 
suppressPackageStartupMessages(library(topGO))
suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(stringr))

# Load the mapping file 
geneID2GO <- readMappings("geneID2GO.tsv")

# Create a character vector that contains the names of the all the gene IDs in the mapping file
geneUniverse <- names(geneID2GO)

# Load the differential expression data
data <- read_csv("entrez_results.csv")

# Filter for statistically significant upregulated genes
up_genes <- data %>% 
  filter(padj <= 0.05 & log2FoldChange >= 1)

# Filter for statistically signficant downregulated genes
down_genes <- data %>% 
  filter(padj <= 0.05 & log2FoldChange <= -1)

upregulated_genes <- as.character(up_genes$gene_id)
downregulated_genes <- as.character(down_genes$gene_id)

# Get binary values depending on if a gene is upregulated or not (or downregulated or not)
up_gene_list <- factor(as.integer(geneUniverse %in% upregulated_genes))
down_gene_list <- factor(as.integer(geneUniverse %in% downregulated_genes))

# Set names for the gene list.
names(up_gene_list) <- geneUniverse
names(down_gene_list) <- geneUniverse

# Build the GOdata object in topGO for upregulated
up_GO_data <- new("topGOdata", 
                  description = "Control_Infected", 
                  ontology = "BP", 
                  allGenes = up_gene_list,
                  annot = annFUN.gene2GO,
                  gene2GO = geneID2GO)

# Build the GOdata object in topGO for downregulated
down_GO_data <- new("topGOdata",
                    description = "Control_Infected",
                    ontology = "BP",
                    allGenes = down_gene_list,
                    annot = annFUN.gene2GO,
                    gene2GO = geneID2GO)

# Perform stats for upregulated data
up_result <- runTest(up_GO_data,
                     algorithm = "weight01",
                     statistic = "fisher")

# Perform stats for downregulated data
down_result <- runTest(down_GO_data,
                       algorithm = "weight01",
                       statistic = "fisher")

# Extract a summary of upregulated results

up_GO <- GenTable(up_GO_data,
                  weight01 = up_result,
                  orderBy = "weight01",
                  topNodes = 50,
                  numChar = 5000)

write.csv(up_GO, file = "upregulated_GO_results.csv", row.names = FALSE)

# Extract a summary of downregulated results
down_GO <- GenTable(down_GO_data,
                    weight01 = down_result,
                    orderBy = "down_result",
                    ranksOf = "down_result",
                    topNodes = 50,
                    numChar = 5000)

## Visualizing TopGo Analyses

# Filter out any non-significant data and calculate the gene ratio
up_GO_filtered <- up_GO %>%
  mutate(GeneRatio = Significant/Annotated, weight01 = as.numeric(weight01)) %>%
  filter(weight01 <= 0.05) %>%
  head(n = 20)

down_GO_filtered <- down_GO %>%
  mutate(GeneRatio = Significant/Annotated, weight01 = as.numeric(weight01)) %>%
  filter(weight01 <= 0.05) %>%
  head(n = 20)

# First, let's arrange the data based on the enrichment ratio. 

up_GO_filtered_arranged <- up_GO_filtered %>% 
  arrange(GeneRatio) %>%
  mutate(Term = factor(Term, levels = Term))  # Use Term as levels in arranged order

down_GO_filtered_arranged <- down_GO_filtered %>% 
  arrange(GeneRatio) %>%
  mutate(Term = factor(Term, levels = Term))

# Combine both datasets for determining common scales
combined_GO <- bind_rows(
  up_GO_filtered_arranged %>% mutate(Direction = "Upregulated"),
  down_GO_filtered_arranged %>% mutate(Direction = "Downregulated")
)

# Determine common ranges for color scale
pvalue_range <- range(combined_GO$weight01, na.rm = TRUE)

# Create custom breaks for size scale - showing 100, 200, and 300
size_breaks <- c(100, 200, 300)
size_limits <- c(min(c(min(combined_GO$Significant), 100)), 
                 max(c(max(combined_GO$Significant), 300)))

# Create the upregulated plot with fixed scales
up_plot <- ggplot(up_GO_filtered_arranged, 
                  aes(x = Term, y = GeneRatio, color = weight01, size = Significant)) +
  geom_point() +
  coord_flip() +
  scale_color_gradient(
    low = "red", 
    high = "blue", 
    name = "P-value",
    limits = pvalue_range,
    breaks = pretty(pvalue_range, n = 5)
  ) +
  scale_size_continuous(
    name = "Number of\nSignificant Genes",
    limits = size_limits,  # Set common limits
    range = c(3, 10),  # Adjust point size range as needed
    breaks = size_breaks  # Explicitly set breaks to show 100, 200, 300
  ) +
  theme_light() +
  labs(
    title = "Upregulated",
    x = "GO Term Description", 
    y = "Enrichment Ratio"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.y = element_text(size = 10),
    legend.position = "right",
    legend.title = element_text(size = 10)
  )

print(up_plot)

# Create the downregulated plot with the SAME scales
down_plot <- ggplot(down_GO_filtered_arranged, 
                    aes(x = Term, y = GeneRatio, color = weight01, size = Significant)) +
  geom_point() +
  coord_flip() +
  scale_color_gradient(
    low = "red", 
    high = "blue", 
    name = "P-value",
    limits = pvalue_range,
    breaks = pretty(pvalue_range, n = 5)
  ) +
  scale_size_continuous(
    name = "Number of\nSignificant Genes",
    limits = size_limits,  # Same limits as upregulated
    range = c(3, 10),  # Same point size range
    breaks = size_breaks  # Same breaks: 100, 200, 300
  ) +
  theme_light() +
  labs(
    title = "Downregulated",
    x = "GO Term Description", 
    y = "Enrichment Ratio"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.y = element_text(size = 10),
    legend.position = "right",
    legend.title = element_text(size = 10)
  )

print(down_plot)

# Save upregulated plot
ggsave("upregulated_GO_plot.png", 
       plot = up_plot,
       width = 10, 
       height = 8,
       dpi = 600)

# Save downregulated plot
ggsave("downregulated_GO_plot.png", 
       plot = down_plot,
       width = 10, 
       height = 8,
       dpi = 600)

## OLD
# Now let's extract the order of the term column
up_order_term <- up_GO_filtered_arranged %>% 
  pull(Term) # pull() extracts a column as a vector

down_order_term <- down_GO_filtered_arranged %>% 
  pull(Term) # pull() extracts a column as a vector

joined_GO_filtered_arranged <- bind_rows(
  up_GO_filtered_arranged %>% mutate(up_down = "UP"),
  down_GO_filtered_arranged %>% mutate(up_down = "DOWN")
)
joined_GO_filtered_arranged$up_down <- factor(joined_GO_filtered_arranged$up_down, levels = c("UP", "DOWN"))

all_terms <- unique(c(as.character(up_GO_filtered_arranged$Term), as.character(down_GO_filtered_arranged$Term)))
joined_GO_filtered_arranged$Term <- factor(joined_GO_filtered_arranged$Term, levels = all_terms)

ggplot(joined_GO_filtered_arranged, aes(x = Term, y = GeneRatio, color = weight01, size = Significant)) +
  geom_point() +
  coord_flip() +
  scale_color_gradient(low = "red", high = "blue") +
  facet_grid(. ~ up_down) +
  theme_light() +
  labs(x = "GO Term Description", y = "Enrichment Ratio", color = "P-value", size = "Number of Significant Genes")




