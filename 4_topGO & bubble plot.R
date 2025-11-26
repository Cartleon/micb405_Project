suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(pheatmap)) 
suppressPackageStartupMessages(library(topGO))
suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(ggplot2))

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
                  orderBy = "up_result",
                  ranksOf = "up_result",
                  topNodes = 50)

# Extract a summary of downregulated results
down_GO <- GenTable(down_GO_data,
                    weight01 = down_result,
                    orderBy = "down_result",
                    ranksOf = "down_result",
                    topNodes = 50)

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
  mutate(Term = factor(Term))

down_GO_filtered_arranged <- down_GO_filtered %>% 
  arrange(GeneRatio) %>%
  mutate(Term = factor(Term))

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
