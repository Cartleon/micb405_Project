suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(pheatmap))
suppressPackageStartupMessages(library(RColorBrewer))

# Set path to a directory with files and load files into R, can do setwd("YourPath/")
# Control replicate 1
ctrl_rep1 <- read_tsv("aligned_ctrl_rep1_ReadsPerGene.out.tab",
                      col_names = c("gene_id", "total","antisense", "sense"),
                      skip = 4)
# Control replicate 2
ctrl_rep2 <- read_tsv("aligned_ctrl_rep2_ReadsPerGene.out.tab",
                           col_names = c("gene_id", "total","antisense", "sense"),
                           skip = 4)
# Treatment replicate 1
sample_rep1 <- read_tsv("aligned_sample_rep1_ReadsPerGene.out.tab",
                            col_names = c("gene_id", "total","antisense", "sense"),
                            skip = 4)

# Treatment replicate 2
sample_rep2 <- read_tsv("aligned_sample_rep2_ReadsPerGene.out.tab",
                            col_names = c("gene_id", "total","antisense", "sense"),
                            skip = 4)

# Note the column names assigned - these will be important when we set our metadata file 
dat <- data.frame(row.names = ctrl_rep1$gene_id,
                  ctrl_rep1 = ctrl_rep1$sense,
                  ctrl_rep2 = ctrl_rep2$sense,
                  sample_rep1 = sample_rep1$sense,
                  sample_rep2 = sample_rep2$sense)

# Let's transform dat into a matrix
dat_matrix<- as.matrix(dat) 

# Look at the first 10 rows of the matrix
view(dat_matrix)

# Check the type of object that dat_matrix is
class(dat_matrix) 

# Make our Metadata file that contains our column information for our matrix
metadata <- data.frame(row.names = colnames(dat_matrix), 
                       condition = c("ctrl", "ctrl", "infected", "infected")
)

colnames(dat_matrix) == rownames(metadata)

## Running DESeq

# Create our DESeq2 object
dds_matrix <- DESeqDataSetFromMatrix(countData = dat_matrix, #matrix 
                                     colData = metadata, #metadata file
                                     design = ~condition)

# Set control condition using the relevel function
dds_matrix$condition <- relevel(dds_matrix$condition, ref = "ctrl")

# Check the levels of dds_matrix$condition
levels(dds_matrix$condition)

dds <- DESeq(dds_matrix)

saveRDS(dds, "dds.rds")

## Sanity checks

# Calculate distances between samples in our log-transformed data
sample_dists <- dist(t(assay(rld)))

# Convert the output to a matrix
sample_dist_matrix <- as.matrix(sample_dists)

# Remove the column names of our matrix
colnames(sample_dist_matrix) <- NULL

# Set the colour palette for our heatmap
colours <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)

# Generate a heatmap using the pheatmap package
pheatmap(sample_dist_matrix,
         clustering_distance_rows = sample_dists,
         clustering_distance_cols = sample_dists, 
         col = colours)

## Extracting data

# Names of the results that DESeq2 calculated
resultsNames(dds)

# Now we will extract the results for our comparison between the 12h timepoint and the 1h timepoint
res <- results(dds, name = "condition_infected_vs_ctrl") %>% as.data.frame() # we save it as a dataframe for easy manipulation with dplyr 
head(res)

## Data manipulation using dplyr
glimpse(res)

res_no_NA <- res %>% 
  drop_na()

# How many rows did we filter out!?
view(res_no_NA) 
write_csv(res_no_NA, "res_no_NA_results.csv") #no gene_id

# Results with an adjusted p-value < 0.05 for diferentially expressed genes
res_filtered <- res_no_NA %>% 
  filter(padj <= 0.05)

# How many rows did we filter out!?
glimpse(res_filtered)

# Biologically relevant genes
res_filtered_final <- res_filtered %>% 
  filter(log2FoldChange <= -1 | log2FoldChange >= 1) %>% # the '|' stand for OR here!
  rownames_to_column("gene_id") # Convert the rownames into a column so they can be saved in your CSV file

head(res_filtered_final)
write_csv(res_filtered_final, "filtered_results.csv")
