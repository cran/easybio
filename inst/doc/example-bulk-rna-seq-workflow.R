## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.align = 'center'
)

## -----------------------------------------------------------------------------
library(easybio)
# x <- prepare_geo('gseid')

## ----echo=FALSE---------------------------------------------------------------
knitr::include_graphics('limmaFilter.png')

## ----echo=FALSE---------------------------------------------------------------
knitr::include_graphics("limmaNorm.png")

## ----fig.dim=c(5, 3)----------------------------------------------------------
data(CHOL_DEGs)
plotVolcano(
  data = CHOL_DEGs,
  x = logFC,
  y = -log10(adj.P.Val),
  color = tumor_vs_normal
)

## ----eval=FALSE---------------------------------------------------------------
#  # devtools::install_github('person-c/easybio')
#  library(data.table)
#  # Over Presentation Analysis(ORA)
#  pathwayGO <- r4msigdb::query("Hs", pathway = "^GO(BP|CC|MF)_")
#  pathwayGO <- setNames(pathwayGO$symbol, pathwayGO$standard_name)
#  
#  oraRes <- fgsea::fora(
#    pathways = pathwayGO,
#    genes = CHOL_DEGs[.("Up"), gene_name, on = .(tumor_vs_normal)],
#    universe = unique(CHOL_DEGs$gene_name)
#  )
#  oraRes[, let(
#    category = fcase(
#      pathway %like% "GOBP", "BP",
#      pathway %like% "GOMF", "MF",
#      pathway %like% "GOCC", "CC"
#    )
#  )]
#  
#  oraRes <- oraRes[, .SD[order(padj)], by = .(category)]
#  oraRes[, let(pathwayGO = factor(pathway, levels = rev(pathway)))]
#  plotORA(
#    data = oraRes[, head(.SD, 5), by = category],
#    x = -log10(padj),
#    y = pathwayGO,
#    size = log10(overlap),
#    fill = category
#  )

## ----echo=FALSE---------------------------------------------------------------
knitr::include_graphics('ora.png')

## -----------------------------------------------------------------------------
library(fgsea)
data(examplePathways)
data(exampleRanks)

## ----fig.width=4, fig.height=2.7----------------------------------------------
fgseaRes <- fgsea(pathways = examplePathways, 
                  stats    = exampleRanks,
                  minSize  = 15,
                  maxSize  = 500)
plotGSEA(
  fgseaRes, 
  pathways = examplePathways, 
  pwayname = "5991130_Programmed_Cell_Death", 
  stats = exampleRanks, 
  save = FALSE
)

