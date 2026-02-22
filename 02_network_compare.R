# Construct and compare SparCC networks - SIMPLIFIED VERSION

library(NetCoMi)
library(phyloseq)

load("results/phyloseq_object.RData")

group_var    <- "Timepoint_Group"
group1_label <- "T1"
group2_label <- "T2"
out_prefix   <- file.path("results", "timepoint_T1_vs_T2")

cat("Using 1 core (sequential processing for cleaner output)\n")

subset_two_groups <- function(phy, group_var, g1, g2) {
  sdata <- as(sample_data(phy), "data.frame")
  g1_samples <- sdata[[group_var]] == g1
  g2_samples <- sdata[[group_var]] == g2
  g1_phy <- prune_samples(g1_samples, phy)
  g2_phy <- prune_samples(g2_samples, phy)
  
  n1 <- nsamples(g1_phy)
  n2 <- nsamples(g2_phy)
  cat(sprintf("Group %s: %d samples\n", g1, n1))
  cat(sprintf("Group %s: %d samples\n", g2, n2))
  
  if (n1 < 10 || n2 < 10) {
    warning("Small sample size detected. Results may not be reliable.")
  }
  
  list(group1 = g1_phy, group2 = g2_phy)
}

collapse_to_genus <- function(phy, rank = "Rank6") {
  phy_genus <- phyloseq::tax_glom(phy, taxrank = rank, NArm = FALSE)
  taxtab <- phy_genus@tax_table@.Data
  genus_names <- make.unique(as.character(taxtab[, rank]))
  rownames(phy_genus@otu_table@.Data) <- genus_names
  phy_genus
}

groups       <- subset_two_groups(phylo_obj, group_var, group1_label, group2_label)
phy_g1_genus <- collapse_to_genus(groups$group1)
phy_g2_genus <- collapse_to_genus(groups$group2)

cat("\n=== Constructing networks ===\n")

# Simplified network construction with safer parameters
net <- netConstruct(
  data        = phy_g1_genus,
  data2       = phy_g2_genus,
  measure     = "sparcc",
  measurePar  = list(iter = 20, inner_iter = 10),  # Faster for small data
  normMethod  = "none",
  zeroMethod  = "multRepl",
  sparsMethod = "none",        # NO sparsification - keep all edges
  verbose     = 2,
  seed        = 123456
)

cat("\n=== Analyzing networks ===\n")

# Simplified network analysis
props <- netAnalyze(
  net,
  centrLCC   = TRUE,
  clustMethod = "cluster_fast_greedy",
  hubPar     = "degree",
  hubQuant   = 0.95,
  normDeg    = TRUE,
  normBetw   = FALSE,
  normClose  = FALSE,
  normEigen  = FALSE
)

cat("\n=== Comparing networks ===\n")

# Simplified comparison - no permutation test
comp <- netCompare(
  props,
  permTest = FALSE,
  verbose  = TRUE
)

res <- list(net = net, props = props, comp = comp)
save(res, file = paste0(out_prefix, "_results_simplified.RData"))

cat("\n=== Creating summary ===\n")

# Summary without permutation
summary_file <- paste0(out_prefix, "_summary_simplified.txt")
sink(file = summary_file)
cat("NETWORK COMPARISON SUMMARY (No Permutation Test)\n")
summary(comp, groupNames = c(group1_label, group2_label))
sink()

cat("\n=== Creating plot ===\n")

pdf_file <- paste0(out_prefix, "_network_plot_simplified.pdf")

# Checking with error because it wasn't working
tryCatch({
  pdf(file = pdf_file, width = 14, height = 7)
  
  plot(
    props,
    sameLayout    = TRUE,
    layoutGroup   = "union",
    repulsion     = 0.9,
    shortenLabels = "simple",
    labelLength   = 20,
    nodeSize      = "degree",
    nodeSizeSpread = 3,
    nodeColor     = "cluster",
    hubBorderCol  = "black",
    hubBorderWidth = 2,
    cexNodes      = 1.5,
    cexLabels     = 1.0,
    cexTitle      = 1.8,
    cexHubLabels  = 1.2,
    groupNames    = c(group1_label, group2_label),
    mar           = c(2, 2, 4, 2)
  )
  
  dev.off()
  cat("✓ Plot created successfully!\n")
  
}, error = function(e) {
  if (dev.cur() > 1) dev.off()
  cat("✗ Plot failed:", conditionMessage(e), "\n")
})

cat("\nSimplified analysis complete!\n")
cat("Results saved to:", paste0(out_prefix, "_results_simplified.RData"), "\n")
cat("Summary:", summary_file, "\n")
if (file.exists(pdf_file)) {
  cat("Plot:", pdf_file, "\n")
}