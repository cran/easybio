litedown::reactor(warning = FALSE) # vignette setting

library(easybio)
library(Seurat)
library(data.table)

# The pbmc.markers dataset is included in easybio
head(pbmc.markers)

marker_matched <- matchCellMarker2(marker = pbmc.markers, n = 50, spc = "Human")

# Let's look at the top 2 potential cell types for each cluster
marker_matched[, head(.SD, 2), by = cluster]

cl2cell_auto <- marker_matched[, head(.SD, 1), by = .(cluster)]
cl2cell_auto <- setNames(cl2cell_auto[["cell_name"]], cl2cell_auto[["cluster"]])
print("Initial automated annotation:")
cl2cell_auto

plotPossibleCell(marker_matched[, head(.SD), by = .(cluster)], min.uniqueN = 2)

# Let's investigate clusters 1, 5, and 7
local_evidence <- check_marker(marker_matched, cl = c(1, 5, 7), topcellN = 2, cis = TRUE)
print(local_evidence)

canonical_markers <- check_marker(marker_matched, cl = c(1, 5, 7), topcellN = 2, cis = FALSE)
print(canonical_markers)

# For this example to be runnable, we need a Seurat object.
# We'll create a minimal one. In your real workflow, you would use your own srt object.
marker_genes <- unique(pbmc.markers$gene)
counts <- matrix(
  abs(rnorm(length(marker_genes) * 50, mean = 1, sd = 2)),
  nrow = length(marker_genes),
  ncol = 50
)
rownames(counts) <- marker_genes
colnames(counts) <- paste0("cell_", 1:50)
srt <- CreateSeuratObject(counts = counts)
# Assign clusters that match the pbmc.markers data
srt$seurat_clusters <- sample(0:8, 50, replace = TRUE)
Idents(srt) <- "seurat_clusters"


# Now, let's plot the evidence for clusters 1, 5, and 7
matchCellMarker2(marker = pbmc.markers, n = 50, spc = "Human") |>
  check_marker(cl = c(1, 5, 7), topcellN = 2, cis = TRUE) |>
  plotSeuratDot(srt = srt)

# Based on our exploration, we finalize the annotations
cl2cell_final <- finsert(
  list(
    c(3) ~ "B cell",
    c(8) ~ "Megakaryocyte",
    c(7) ~ "DC",
    c(1, 5) ~ "Monocyte",
    c(0, 2, 4) ~ "Naive CD8+ T cell",
    c(6) ~ "Natural killer cell"
  ),
  len = 9 # Ensure vector length covers all clusters (0-8)
)
print("Final curated annotation:")
cl2cell_final

custom_ref_list <- list(
  "T-cell" = c("CD3D", "CD3E", "CD3G"),
  "B-cell" = c("CD79A", "MS4A1"),
  "Myeloid" = c("LYZ", "CST3", "AIF1")
)
print(custom_ref_list)

custom_ref_df <- list2dt(custom_ref_list, col_names = c("cell_name", "marker"))
head(custom_ref_df)

marker_custom <- matchCellMarker2(
  marker = pbmc.markers,
  n = 50,
  ref = custom_ref_df
)
# Note that the cell_name column now contains our custom cell types
marker_custom[, head(.SD, 2), by = cluster]

get_marker(spc = "Human", cell = c("Monocyte", "Neutrophil"), number = 5, min.count = 1)

plotMarkerDistribution(mkr = "CD68")

