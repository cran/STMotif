## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  warning = FALSE,
  comment = "#>"
)

## ----load-package-------------------------------------------------------------
library(STMotif)

## ----eval = FALSE-------------------------------------------------------------
# install.packages("STMotif")

## ----normsax------------------------------------------------------------------
# Load the example dataset
dim(D <- STMotif::example_dataset)

# Normalization and SAX encoding
DS <- NormSAX(D = STMotif::example_dataset, a = 5)

# Preview the SAX-encoded dataset
head(NormSAX(D = STMotif::example_dataset, a = 5)[, 1:10])

## ----search-------------------------------------------------------------------
# Discover motifs
stmotifs <- SearchSTMotifs(D, DS, 4, 5, 4, 10, 2, 2)
stmotifs[[1]]

## ----rank---------------------------------------------------------------------
# Rank the discovered motifs
rstmotifs <- RankSTMotifs(stmotifs)
rstmotifs[[1]]

## ----csa----------------------------------------------------------------------
# Full CSA workflow in one call
rstmotifs <- CSAMiningProcess(D, DS, 4, 5, 4, 10, 2, 2)
rstmotifs[[1]]

## ----fig-heatmap, fig.height = 4, fig.width = 5, fig.align = "center"---------
display_motifsDataset(
  dataset = STMotif::example_dataset,
  rstmotifs[c(1:4)],
  5
)

## ----fig-series, fig.height = 4, fig.width = 5, fig.align = "center"----------
display_motifsSTSeries(
  dataset = STMotif::example_dataset,
  rstmotifs[c(1:4)],
  space = c(1:4, 10:12)
)

