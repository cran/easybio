litedown::reactor(warning = FALSE) # vignette setting

library(easybio)

data(CHOL_DEGs)
plotVolcano(
  data = CHOL_DEGs,
  x = logFC,
  y = -log10(adj.P.Val),
  color = tumor_vs_normal
)

library(fgsea)
data(examplePathways)
data(exampleRanks)

fgseaRes <- fgsea(
  pathways = examplePathways,
  stats = exampleRanks,
  minSize = 15,
  maxSize = 500
)
plotGSEA(
  fgseaRes,
  pathways = examplePathways,
  pwayname = "5991130_Programmed_Cell_Death",
  stats = exampleRanks,
  save = FALSE
)

