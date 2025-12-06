library(tidyverse)
library(DESeq2)

# Read the geneID2GO file
gene2go <- read_tsv("geneID2GO.tsv", col_names = c("gene_id", "go_terms"))

# GO terms of interes
go_terms_of_interest <- c(
  "GO:0006954" = "Inflammatory response",
  "GO:0050729" = "Positive regulation of inflammatory response",
  "GO:0043123" = "Positive regulation of canonical NF-kappaB signal transduction",
  "GO:0032729" = "Positive regulation of type II interferon production",
  "GO:0032728" = "Positive regulation of interferon-beta production",
  "GO:0071346" = "Cellular response to type II interferon",
  "GO:0035458" = "Cellular response to interferon-beta",
  "GO:0045087" = "Innate immune response"
)

## List to store genes for each GO term

# inflammatory response
inflammatory_genes <- gene2go %>%
  filter(str_detect(go_terms, "GO:0006954")) %>%
  pull(gene_id) %>%
  unique()

# positive regulation of inflammatory response
inflammatory_Reg_genes <- gene2go %>%
  filter(str_detect(go_terms, "GO:0050729")) %>%
  pull(gene_id) %>%
  unique()

# positive regulation of canonical NF-kappaB signal transduction
NFkB_Reg_genes <- gene2go %>%
  filter(str_detect(go_terms, "GO:0043123")) %>%
  pull(gene_id) %>%
  unique()

# positive regulation of type II interferon production
IFNG_Reg_genes <- gene2go %>%
  filter(str_detect(go_terms, "GO:0032729")) %>%
  pull(gene_id) %>%
  unique()

# positive regulation of interferon-beta production
IFNB1_Reg_genes <- gene2go %>%
  filter(str_detect(go_terms, "GO:0032728")) %>%
  pull(gene_id) %>%
  unique()

# cellular response to type II interferon
IFNG_genes <- gene2go %>%
  filter(str_detect(go_terms, "GO:0071346")) %>%
  pull(gene_id) %>%
  unique()

# cellular response to interferon-beta
IFNB1_genes <- gene2go %>%
  filter(str_detect(go_terms, "GO:0035458")) %>%
  pull(gene_id) %>%
  unique()

# innate immune response
Innate_genes <- gene2go %>%
  filter(str_detect(go_terms, "GO:0045087")) %>%
  pull(gene_id) %>%
  unique()

# load in data
dat <- read_csv("entrez_results_unfiltered.csv")
dat <- dat %>%
  mutate(log_padj = ifelse(padj == 0, 300, -log10(padj))) # cap infinite values

labelled_dat <- dat %>%
  mutate(
    Genes = case_when(
      padj < 0.05 & log2FoldChange > 1 ~ "Upregulated",
      padj < 0.05 & log2FoldChange < -1 ~ "Downregulated",
      TRUE ~ "Not Significant"
    ),
    Category = case_when(
      gene_id %in% inflammatory_genes ~ "Inflammatory response",
      gene_id %in% inflammatory_Reg_genes ~ "Positive regulation of inflammatory response",
      gene_id %in% NFkB_Reg_genes ~ "Positive regulation of canonical NF-kappaB signal transduction",
      gene_id %in% IFNG_Reg_genes ~ "Positive regulation of type II interferon production",
      gene_id %in% IFNB1_Reg_genes ~ "Positive regulation of interferon-beta production",
      gene_id %in% IFNG_genes ~ "Cellular response to type II interferon",
      gene_id %in% IFNB1_genes ~ "Cellular response to interferon-beta",
      gene_id %in% Innate_genes ~ "Innate immune response",
      TRUE ~ "Other"
    )
  )

# color scheme
volcano_colors <- c(
  "Inflammatory response" = "red",
  "Positive regulation of inflammatory response" = "purple",
  "Positive regulation of canonical NF-kappaB signal transduction" = "blue",
  "Positive regulation of type II interferon production" = "green",
  "Positive regulation of interferon-beta production" = "#ff7f0e",
  "Cellular response to type II interferon" = "#F0E442",
  "Cellular response to interferon-beta" = "brown",
  "Innate immune response" = "#e377c2",
  "Other" = "gray70"
)

labelled_dat %>%
  ggplot(aes(x = log2FoldChange, y = log_padj, color = Category)) +
  geom_point(aes(alpha = Category), size = 1) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", linewidth = 0.8) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", linewidth = 0.8) +
  scale_color_manual(
    values = volcano_colors,
    breaks = c(
      "Inflammatory response",
      "Positive regulation of inflammatory response",
      "Positive regulation of canonical NF-kappaB signal transduction",
      "Positive regulation of type II interferon production",
      "Positive regulation of interferon-beta production",
      "Cellular response to type II interferon",
      "Cellular response to interferon-beta",
      "Innate immune response"
    )
  ) +
  scale_alpha_manual(
    values = c(
      "Inflammatory response" = 1,
      "Positive regulation of inflammatory response" = 1,
      "Positive regulation of canonical NF-kappaB signal transduction" = 1,
      "Positive regulation of type II interferon production" = 1,
      "Positive regulation of interferon-beta production" = 1,
      "Cellular response to type II interferon" = 1,
      "Cellular response to interferon-beta" = 1,
      "Innate immune response" = 1,
      "Other" = 0.1  # Lower opacity for "Other"
    ),
    guide = "none"  # Hide alpha from legend
  ) +
  labs(
    x = expression(log[2]("Fold Change")),
    y = expression(-log[10]("p-adjusted")),
    color = "Gene Category"
  ) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", linewidth = 1.5)
  ) +
  coord_cartesian(xlim = c(-10, 10))


## Summary
total_upregulated <- labelled_dat %>%
  filter(Genes == "Upregulated") %>%
  nrow()
cat("Total upregulated genes:", total_upregulated, "\n")

total_downregulated <- labelled_dat %>%
  filter(Genes == "Downregulated") %>%
  nrow()
cat("Total downregulated genes:", total_upregulated, "\n")


upregulated_detailed <- labelled_dat %>%
  filter(Genes == "Upregulated") %>%
  group_by(Category) %>%
  summarize(
    Gene_Count = n(),
    Percentage = round(n() / total_upregulated * 100, 2)
  ) %>%
  arrange(desc(Gene_Count))
print(upregulated_detailed)

downregulated_detailed <- labelled_dat %>%
  filter(Genes == "Downregulated") %>%
  group_by(Category) %>%
  summarize(
    Gene_Count = n(),
    Percentage = round(n() / total_upregulated * 100, 2)
  ) %>%
  arrange(desc(Gene_Count))
print(downregulated_detailed)
