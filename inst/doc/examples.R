## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ---- echo=FALSE---------------------------------------------------------
source(file = "../R/mainFunction.R")
source(file = "../R/subFunction.R")
source(file = "../R/visualization.R")
library(foreach)
library(stats)
library(ggplot2)
library(reshape2)

## ---- echo=TRUE----------------------------------------------------------

# The process is launched on the provided example dataset
dim(STMotif::example_dataset)


# Generation of candidates
candidates <- STSIdentifyCandidateSTMotifs(data = STMotif::example_dataset, tslice = 25, sslice = 10, alpha = 7, window_size = 3)


# Information of the block 1
# The candidates built 
head(candidates$motifs[[1]]$Subs)[1:3,]

# The candidates built with SAX
head(candidates$motifs[[1]]$Subs.SAX)[1:3,]

# Grouping identical candidates 
head(candidates$motifs[[1]]$Motif.raw)[1:3]

# Grouping identical motifs with SAX
head(candidates$motifs[[1]]$Motif.SAX)[1:3]

# Grouping the starting position of identical candidates in the combined series
head(candidates$motifs[[1]]$Indices)[1:3]


# Number of blocks in horizontal and vertical
c(candidates$nrows,candidates$ncols)


# Position of the block 1 into the original dataset
candidates$rectangles[[1]]

# Values of block 1
head(candidates$blocks$datasets[[1]])[1:3,]


# Encoded values of block 1
head(candidates$saxblocks$datasets[[1]])[1:3,]

## ---- echo=TRUE----------------------------------------------------------
stmotifs <- STSIdentifySTMotifs(candidates = candidates, ka = 1, si = 1)

# Output of the first selected motif
stmotifs[[1]]

## ---- echo=TRUE----------------------------------------------------------
sttightmotifs <- STSIdentifyTightSTMotifs(stmotifs = stmotifs, rectangles = candidates$rectangles)

# Output of the first motif after removing the isolated
sttightmotifs[[1]]

## ---- echo=TRUE----------------------------------------------------------
ranksttightmotifs <- STSRankTightSTMotifs(sttightmotifs = sttightmotifs)

# Motif with the better quality
ranksttightmotifs[[1]]

## ----fig, fig.height = 4, fig.width = 6, fig.align = "center"------------
# Plot the intensity of the dataset and highlight one selected motif
intensityDataset(dataset = STMotif::example_dataset,rankList = ranksttightmotifs,position = 1,alpha = 7)

## ----fig1, fig.height = 4, fig.width = 6, fig.align = "center"-----------
# Plot five specific spatial-series which some of them contain the best motif
displayPlotSeries(dataset = STMotif::example_dataset, rankList = ranksttightmotifs ,position = 1 ,space = c(1,2,5:7))

