suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(pheatmap)) 
suppressPackageStartupMessages(library(topGO))
suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(purrr))
go_terms <- c(
  "GO:0006954", "GO:0050729", "GO:0071222", "GO:0032760", "GO:0043123", 
  "GO:0032729", "GO:0032728", "GO:0071346", "GO:0035458", "GO:0045071", 
  "GO:0140374", "GO:0045087"
)

geneID2GO <- readMappings("geneID2GO.tsv")
joined_GO_filtered_arranged <- read.table("joined_GO_filtered_arranged.tsv", header=TRUE, sep="\t", stringsAsFactors=FALSE)
data <- read_csv("entrez_results.csv")

geneUniverse <- names(geneID2GO)
up_genes <- data %>% filter(padj <= 0.05 & log2FoldChange >= 1) %>% pull(gene_id) %>% as.character()
up_gene_list <- factor(as.integer(geneUniverse %in% up_genes))
names(up_gene_list) <- geneUniverse

get_sig_genes <- function(go_terms, gene_list) {
  data.frame(GO = as.character(go_terms), stringsAsFactors = FALSE) %>%
    mutate(
      GeneIDs = lapply(GO, function(goid) {
        annotated <- names(geneID2GO)[sapply(geneID2GO, function(x) goid %in% x)]
        intersect(annotated, names(gene_list)[gene_list == 1])
      })
    )
}

go_to_gene <- get_sig_genes(go_terms, up_gene_list)
go_to_gene <- merge(
  go_to_gene,
  joined_GO_filtered_arranged[, c("GO.ID", "Term")],
  by.x = "GO", by.y = "GO.ID",
  all.x = TRUE
)
go_to_gene <- go_to_gene[, c("GO", "Term", "GeneIDs")]

# --- Replace dplyr::select with base R subsetting for go_long (no select error) ---
go_long <- unnest(go_to_gene, cols = GeneIDs)
go_long$gene_id <- as.character(go_long$GeneIDs)
go_long <- go_long[, c("GO", "Term", "gene_id")]

# (Continue your pipeline as before)
data$gene_id <- as.character(data$gene_id)
merged <- left_join(go_long, data, by = "gene_id")

go_means <- merged %>%
  group_by(GO, Term) %>%
  summarise(
    mean_log2FoldChange = mean(log2FoldChange, na.rm = TRUE),
    stderror_mean_log2FoldChange = sd(log2FoldChange, na.rm = TRUE) / sqrt(sum(!is.na(log2FoldChange))),
    mean_lfcSE = mean(lfcSE, na.rm = TRUE),
    .groups = "drop"
  )

view(go_means)
write.table(go_means, file = "go_mean_stats.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

## Plot

# merged <- left_join(go_long, data, by = "gene_id")

# Make a box plot of log2FoldChange grouped by GO term
ggplot(merged, aes(x = Term, y = log2FoldChange)) +
  geom_boxplot(fill = "lightblue") +
  xlab("GO Term") +
  ylab(expression(log[2]("Fold Change"))) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(size = 14),
    plot.margin = margin(t=10, r=10, b=10, l=40)  # Increase left side)
     )
