## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## -----------------------------------------------------------------------------
library(easybio)
# x <- prepare_geo('gseid')

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

## -----------------------------------------------------------------------------
foraRes <- fora(examplePathways, genes=tail(names(exampleRanks), 200), universe=names(exampleRanks))

## ----fig.width=8, fig.height=2------------------------------------------------
# Adjust the pathway position on the y-axis based on the adjusted p-value (padj)
foraRes[, pathway := factor(pathway, levels = rev(pathway))]

plotORA(data = foraRes[1:8], x = -log10(padj), y = pathway, size = overlap, fill = 'constant')

