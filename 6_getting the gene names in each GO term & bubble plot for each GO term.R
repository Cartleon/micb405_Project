library(tidyverse)
library(stringr)
library(dplyr)

# Step 1: Load and process the GTF to extract Entrez-to-gene_name mapping
gtf <- read.delim("GCF_000001635.27_GRCm39_genomic.gtf", comment.char = "#", header = FALSE, stringsAsFactors = FALSE)

gtf_genes <- gtf[gtf$V3 == "gene", ] %>%
  mutate(
    gene_name = str_extract(V9, 'gene ([^;]+);'),
    entrez = str_extract(V9, 'GeneID:([0-9]+)')
  ) %>%
  mutate(
    gene_name = str_trim(str_replace_all(gene_name, 'gene ', '')),
    gene_name = str_replace_all(gene_name, "^;+|;+$", ""),
    entrez = str_replace_all(entrez, 'GeneID:', '')
  ) %>%
  dplyr::select(gene_name, entrez) %>%
  distinct()

# Step 2: Unnest go_to_gene as before
go_long <- go_to_gene %>%
  unnest(cols = GeneIDs) %>%
  mutate(entrez = as.character(GeneIDs)) %>%
  dplyr::select(GO, Term, entrez)

# Step 3: Map Entrez IDs to gene names in go_long
go_named <- go_long %>%
  left_join(gtf_genes, by = "entrez")

# Step 4: Summarize common gene names by GO
go_gene_names <- go_named %>%
  group_by(GO, Term) %>%
  summarise(
    gene_names = paste(na.omit(unique(gene_name)), collapse = ";"),
    .groups = "drop"
  )

# Step 5: View and/or save the table
view(go_gene_names)
write.table(go_gene_names, file = "go_common_gene_names.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

# Make a simple table: Each row = GO, Term, gene_name (one gene per row)
go_gene_name_table <- go_named %>%
  filter(!is.na(gene_name)) %>%
  arrange(GO, Term, gene_name)

go_gene_name_table <- go_gene_name_table[, c("GO", "Term", "gene_name")]

# View or export as needed
View(go_gene_name_table)
write.table(go_gene_name_table, file = "go_gene_name_table.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

# Getting stats from DSeq for genes
data$gene_id <- as.character(data$gene_id)
de_with_names <- merge(data, gtf_genes, by.x = "gene_id", by.y = "entrez", all.x = TRUE)

# --- Step 1: Prepare a vector of gene names associated to GO terms ---
genes_with_go <- unique(go_gene_name_table$gene_name)

# --- Step 2: Filter DE stats to genes ONLY in at least one GO term ---
de_with_names_in_go <- de_with_names %>%
  filter(gene_name %in% genes_with_go)

# --- Step 3: Join to GO associations for each gene (so every row is GO x gene_name x stats) ---
go_stats <- merge(go_gene_name_table, de_with_names_in_go, by = "gene_name", all.x = TRUE)

go_stats_table <- go_stats[, c("GO", "Term", "gene_name", "log2FoldChange", "lfcSE", "padj")]

# --- View or export result ---
View(go_stats_table)

# Pick the GO ID you want (from list: "GO:0006954", "GO:0050729", "GO:0043123", "GO:0032729", "GO:0032728", "GO:0071346", "GO:0035458", "GO:0045087")
go_id <- "GO:0045087"    # <----- Change this to the desired GO ID

# Create the sub_table for this GO term only
sub_table <- go_stats_table %>%
  filter(GO == go_id) %>%
  arrange(desc(log2FoldChange)) # Sort genes by highest FC

# You can print or View(sub_table) to inspect
view(sub_table)

# Get the GO term name from this sub_table
term_name  <- unique(sub_table$Term)[1]

# Make term name safe for filenames
safe_term  <- gsub("[^A-Za-z0-9_]", "_", term_name)

# Build filename in your existing directory (e.g. GO_sub_tables)
filename_table <- paste0("GO_sub_tables/", safe_term, "_sub_table.tsv")

# Save this GO-specific table
write.table(
  sub_table,
  file      = filename_table,
  sep       = "\t",
  row.names = FALSE,
  quote     = FALSE
)

  
## Plot

# Find min and max values for y-axis
y_min <- floor(min(sub_table$log2FoldChange, na.rm = TRUE))
y_max <- ceiling(max(sub_table$log2FoldChange, na.rm = TRUE))

# Make the term text safe for filename
safe_term <- gsub("[^A-Za-z0-9_]", "_", unique(sub_table$Term)[1])
filename <- paste0("Images/", safe_term, "_genes_plot.png")

p <- ggplot(sub_table, aes(
  x = reorder(gene_name, log2FoldChange),
  y = log2FoldChange,
  color = -log10(padj),
  size = lfcSE
)) +
  geom_point() +
  scale_y_continuous(breaks = seq(y_min, y_max, by = 1)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 9)) +
  labs(
    x = "Gene",
    y = expression(log[2]("Fold Change")),
    color = expression(-log[10]("p-adjusted")),
    size = "lfcSE",
    title = unique(sub_table$Term)[1]
  ) +
  scale_color_gradient(low = "blue", high = "orange")

p

#ggsave(filename, plot = p, width = 9, height = 6, dpi = 300)


## For GO:0045087 and GO:0006954, that have too many genes to be shown in one panel (change N_panels for number of panels)
go_id <- "GO:0006954"
sub_table <- go_stats_table %>%
  filter(GO == go_id) %>%
  arrange(desc(log2FoldChange))   # Sort genes by highest FC (optional)

# Step 1: Split into n (N_panels) panels with cut()
n_genes <- nrow(sub_table)
N_panels <- 2
sub_table <- sub_table %>%
  mutate(panel = cut(
    1:n_genes,
    breaks = N_panels,
    labels = paste0("Panel ", 1:N_panels),
    include.lowest = TRUE
  ))

# Step 2: Plot with facet_wrap (panel strips hidden)
p <- ggplot(sub_table, aes(
  x = reorder(gene_name, log2FoldChange),
  y = log2FoldChange,
  color = -log10(padj),
  size = lfcSE
)) +
  geom_point() +
  facet_wrap(~ panel, ncol = 1, scales = "free_x") +   # vertical panels
  scale_y_continuous(
    breaks = seq(
      floor(min(sub_table$log2FoldChange, na.rm = TRUE)),
      ceiling(max(sub_table$log2FoldChange, na.rm = TRUE)),
      by = 1
    )
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 7),
    strip.text = element_blank()   # <-- this hides panel names
  ) +
  labs(
    x = "Gene",
    y = expression(log[2]("Fold Change")),
    color = expression(-log[10]("p-adjusted")),
    size = "lfcSE",
    title = unique(sub_table$Term)[1]
  ) +
  scale_color_gradient(low = "blue", high = "orange")

# Step 3: Save to PNG (filename uses GO term safely)
safe_term <- gsub("[^A-Za-z0-9_]", "_", unique(sub_table$Term)[1])
filename <- paste0("Images/", safe_term, "_genes_faceted.png")
ggsave(filename, plot = p, width = 10, height = 6, dpi = 300)

