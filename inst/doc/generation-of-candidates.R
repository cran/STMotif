## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
source("../R/subFunction.R")

## ---- echo=TRUE----------------------------------------------------------
head(STMotif::example_dataset[,1:10])
head(round(STSNormalization(vector = as.matrix(STMotif::example_dataset)),digits = 2)[,1:10])

## ---- echo=FALSE, fig.cap="SAX Encoding with 3 letters", out.width = '100%'----
knitr::include_graphics("saxencode.png")

## ---- echo=FALSE, fig.cap="Blocks creation", out.width = '100%'----------
knitr::include_graphics("partitioningintoblocks.png")

## ---- echo=FALSE, fig.cap="Combine the spatial-time series into each block", out.width = '100%'----
knitr::include_graphics("combineseries.png")

## ---- echo=FALSE, fig.cap="Application of motif discovery algorithm", out.width = '100%'----
knitr::include_graphics("motifDiscoveryAlgorithm.png")

