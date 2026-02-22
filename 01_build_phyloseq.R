# Build phyloseq object from example Mothur + metadata files
# Adapted from one of my older analysis scripts, cleaned up for yall's use.

library(NetCoMi)
library(phyloseq)

dir.create("results", showWarnings = FALSE)

shared_file   <- "data/example.shared"
taxonomy_file <- "data/example.taxonomy"
metadata_file <- "data/example_metadata.txt"

group_var    <- "Timepoint_Group"
group1_label <- "T1"
group2_label <- "T2"

build_phyloseq_from_mothur <- function(shared, taxonomy, metadata_path) {
  otu       <- import_mothur(mothur_shared_file = shared)
  tax       <- import_mothur(mothur_constaxonomy_file = taxonomy)
  
  sample_md <- read.delim(
    metadata_path,
    header = TRUE,
    sep = "\t",
    stringsAsFactors = FALSE,
    row.names = 1
  )
  
  phyloseq(otu_table(otu),
           tax_table(tax),
           sample_data(sample_md))
}

phylo_obj <- build_phyloseq_from_mothur(shared_file, taxonomy_file, metadata_file)

# Making sure it's valid
cat("Phyloseq object created:\n")
cat("  Samples:", nsamples(phylo_obj), "\n")
cat("  OTUs:", ntaxa(phylo_obj), "\n")
cat("  Sample groups:\n")
print(table(sample_data(phylo_obj)[[group_var]]))

save(phylo_obj, file = "results/phyloseq_object.RData")
