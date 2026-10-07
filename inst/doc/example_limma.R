litedown::reactor(warning = FALSE) # vignette setting

library(easybio)

data(CHOL_DEGs)
plot_volcano(
  data = CHOL_DEGs,
  x = logFC,
  y = -log10(adj.P.Val),
  color = tumor_vs_normal
)

library(fgsea)
data(examplePathways)
data(exampleRanks)

fgsea_res <- fgsea(
  pathways = examplePathways,
  stats = exampleRanks,
  minSize = 15,
  maxSize = 500,
  # Run on one thread. fgsea otherwise hands the multilevel step to a
  # BiocParallel worker, which intermittently fails to find fgsea's own
  # compiled function and aborts the vignette build ("could not find function
  # fgseaMultilevelCpp"). nproc = 1 does not help: setUpBPPARAM() assigns
  # SerialParam() for it and then overwrites it with MulticoreParam(workers = 1),
  # so the backend has to be passed in. The example is small enough that the
  # parallelism buys nothing.
  BPPARAM = BiocParallel::SerialParam()
)
plot_gsea(
  fgsea_res,
  pathways = examplePathways,
  pwayname = "5991130_Programmed_Cell_Death",
  stats = exampleRanks,
  save = FALSE
)

