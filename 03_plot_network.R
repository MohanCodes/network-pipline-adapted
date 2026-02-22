# Plot network and write permutation summary

library(NetCoMi)

group1_label <- "T1"
group2_label <- "T2"
out_prefix   <- file.path("results", "timepoint_T1_vs_T2")

load(paste0(out_prefix, "_results.RData"))

# Write summary
summary_file <- paste0(out_prefix, "_summary_permutation.txt")
cat("Writing summary to:", summary_file, "\n")
sink(file = summary_file)
summary(res$comp, groupNames = c(group1_label, group2_label))
sink()

cat("\nChecking network structure...\n")
cat("Group 1 edges:", res$net$edgelist1 %>% nrow(), "\n")
cat("Group 2 edges:", res$net$edgelist2 %>% nrow(), "\n")

pdf_file <- paste0(out_prefix, "_network_plot.pdf")
cat("Generating network plot:", pdf_file, "\n")

tryCatch({
  pdf(file = pdf_file, width = 12, height = 8)
  
  plot(
    res$props,
    sameLayout    = TRUE,
    layoutGroup   = "union",
    shortenLabels = "simple",
    labelLength   = 18,
    charToRm      = "_unclassified",
    labelScale    = FALSE,
    nodeFilter    = "none",
    rmSingles     = "none",
    nodeSize      = "fix",
    nodeSizeSpread= 3,
    nodeColor     = "cluster",
    hubBorderCol  = "darkgray",
    cexNodes      = 1.2,
    cexLabels     = 0.9,
    cexTitle      = 1.5,
    groupNames    = c(group1_label, group2_label),
    mar           = c(2, 2, 4, 2)
  )
  
  dev.off()
  cat("\nPlot saved successfully!\n")
  
}, error = function(e) {
  dev.off()
  cat("\nError creating plot:", conditionMessage(e), "\n")
  cat("Network summary and permutation results are still saved in:\n")
  cat("  ", summary_file, "\n")
})

cat("\nAnalysis complete!\n")
cat("  Summary:", summary_file, "\n")
if (file.exists(pdf_file)) {
  cat("  Plot:", pdf_file, "\n")
}