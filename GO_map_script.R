#Script to get GO mapping file
# packages
library(readr); library(dplyr); library(stringr)

# 1) Download & read NCBI gene2go (all species)
#   Columns: tax_id, GeneID, GO_ID, Evidence, Qualifier, GO_term, PubMed, Category
gene2go_url <- "https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2go.gz"
tmp <- tempfile(fileext = ".gz")
download.file(gene2go_url, tmp, mode = "wb")

go <- read_tsv("gene2go.gz", show_col_types = FALSE, comment = "#",
               col_types = "iicccccci",
               col_names = c("tax_id","GeneID","GO_ID","Evidence","Qualifier","GO_term","PubMed","Category"))

# 2) Keep mouse only (NCBI taxid 10090) and (optionally) all three namespaces
mouse_go <- go %>%
  dplyr::filter(tax_id == 10090) %>%
  dplyr::select(GeneID, GO_ID) %>%
  dplyr::distinct()

# 3) Collapse to one line per gene, comma-separated GO IDs
geneID2GO_tbl <- mouse_go %>%
  summarize(GO = paste(unique(GO_ID), collapse = ","), .by = GeneID) %>%
  arrange(GeneID)

# 4) Write in topGO mapping format
outfile <- "mouse_geneID2GO_from_NCBI_gene2go.txt"
write.table(geneID2GO_tbl, file = outfile, sep = "\t",
            row.names = FALSE, col.names = FALSE, quote = FALSE)

# 5) Use in topGO
library(topGO)
geneID2GO <- readMappings(outfile)   # returns a named list: Entrez GeneID -> vector of GO IDs
