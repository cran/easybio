## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## -----------------------------------------------------------------------------
library(easybio)
head(pbmc.markers)

## -----------------------------------------------------------------------------
(marker <- matchCellMarker2(marker = pbmc.markers, n = 50, spc = 'Human')[, head(.SD, 2), by=cluster])

## ----fig.width=7.5, fig.height=6----------------------------------------------
plotPossibleCell(marker)

## -----------------------------------------------------------------------------
cl2cell <- marker[, head(.SD, 1), by = .(cluster)]
cl2cell <- setNames(cl2cell[["cell_name"]], cl2cell[["cluster"]])
cl2cell

## ----eval=FALSE---------------------------------------------------------------
#  cls <- list(
#    c(1, 5, 7),
#    c(8),
#    c(3),
#    c(0,2, 4, 6)
#  )
#  dotplotList <- plotSeuratDot(seuratObject, cls, marker = pbmc.markers, n = 50, spc = 'Human', topcellN = 2)

## -----------------------------------------------------------------------------
cl2cell <- finsert(
  expression(
  c(1, 5) == "Monocyte",
  c(7) == "DC",
  c(8) == "megakaryocyte",
  c(3) == "B.cell",
  c(0, 2) == "Naive.CD8.T.cell",
  c(4) == "Cytotoxic.T.Cell",
  c(6) == "Natural.killer.cell",
), len = 9)
cl2cell

## -----------------------------------------------------------------------------
get_marker(spc = 'Human', cell = c('Monocyte', 'Neutrophil'), number = 5, min.count = 1)

## ----fig.width=7.5, fig.height=7----------------------------------------------
plotMarkerDistribution(mkr = "CD68")

