suppressPackageStartupMessages(library(tidyverse))
library(dplyr)
library(stringr)

# Load GTF file
gtf <- read.delim("GCF_000001635.27_GRCm39_genomic.gtf", comment.char = "#", header = FALSE, stringsAsFactors = FALSE)

# Filter rows where 3rd column (feature type) is "gene"
# Extract gene_name and entrez IDs from 9th column (attributes) using regex
gtf_genes <- gtf[gtf$V3 == "gene", ] %>%
  mutate(
    gene_name = str_extract(V9, 'gene ([^;]+);'),
    entrez = str_extract(V9, 'GeneID:([0-9]+)')
  ) %>%
  mutate(
    gene_name = str_trim(str_replace_all(gene_name, 'gene ', '')),
    gene_name = str_replace_all(gene_name, "^;+|;+$", ""),  # Remove leading/trailing semicolons
    entrez = str_replace_all(entrez, 'GeneID:', '')
  ) %>%
  select(gene_name, entrez) %>%
  distinct()
view(gtf_genes)

# Load your gene list CSV (adjust your filename)
data <- read.csv("filtered_results.csv", stringsAsFactors = FALSE)

# Join your gene list with the GTF mapping to get Entrez IDs
data_mapped <- data %>%
  left_join(gtf_genes, by = c("gene_id" = "gene_name"))

# Replace gene_id with entrez where available, keep original where not
data_mapped$gene_id <- ifelse(!is.na(data_mapped$entrez), data_mapped$entrez, data_mapped$gene_id)

# Drop the extra 'entrez' column if desired
data_mapped <- data_mapped %>%
  select(-entrez)

# Save or use your data
write.csv(data_mapped, "entrez_results.csv", row.names = FALSE)
view(data_mapped)

## Entrez Results Unfiltered
# For GO term Volcano

# Load your gene list CSV (adjust your filename)
data_unfiltered <- read.csv("res_no_NA_results_geneid.csv", stringsAsFactors = FALSE)

# Join your gene list with the GTF mapping to get Entrez IDs
data_mapped_unfilered <- data_unfiltered %>%
  left_join(gtf_genes, by = c("gene_id" = "gene_name"))

# Replace gene_id with entrez where available, keep original where not
data_mapped_unfilered$gene_id <- ifelse(!is.na(data_mapped_unfilered$entrez), data_mapped_unfilered$entrez, data_mapped_unfilered$gene_id)

# Drop the extra 'entrez' column if desired
data_mapped_unfilered <- data_mapped_unfilered %>%
  select(-entrez)
