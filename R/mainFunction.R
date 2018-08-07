#'  Generate candidates
#'
#' Find frequent motifs in a dataset
#' @param data Dataset containing numeric values
#' @param tslice "Time slice" Number of row in each block
#' @param sslice "Space slice" Number of column in each block
#' @param alpha Number of letter to do the encode
#' @param window_size Length of motifs
#' @return List [Motifs, Nrows, Ncols, Rectangles, Blocks, Saxblocks] cantining all information about the blocks : number of blocks in rows and columns, values of each block (numeric and encoded values) and information about frequent motifs in blocks.
#' @return Motifs : Contain for each block the Subs : data frame containing all the subsequences in original data frame, Subs.SAX : data frame containing all the subsequences in SAX representations, Motif.raw : a list showing the motifs discovered in the original data frame, Morif.SAX : a list showing the motifs discovered in SAX representations and Indices : a list showing the starting positions of subsequences for each motif discovered.
#' @return Nrows : Numeric value giving the number of row which delimited the blocks created
#' @return Ncols : Numeric value giving the number of column which delimited the blocks created
#' @return Rectangles : For each block gives information about the position of the block into the orginal dataset
#' @return Blocks : Contains the values of each blocks and the nrows and ncols
#' @return Saxblocks : Contains the values encoded bay SAX of each blocks and the nrows and ncols
#' @importFrom stats quantile sd na.omit
#' @import foreach
#' @note To see all the stages of the generation of candidates : \href{../inst/doc/generation-of-candidates.html}{Generation of candidates}
#' @examples
#' #Generation of candidates
#' candidates <- STSIdentifyCandidateSTMotifs(STMotif::example_dataset, 10, 10, 7, 3)
#' @export
STSIdentifyCandidateSTMotifs <- function(data, tslice, sslice, alpha = 7, window_size = 3) {
  dataset <- STSAdjustBlock(data, tslice, sslice)
  vector <- as.matrix(dataset)
  vector <- as.vector(vector)

  vectorNorm <- STSNormalization(vector)

  saxdataset <- STSSaxEncode(dataset, vectorNorm, alpha)
  saxblocks <- STSComputeBlocks(saxdataset, tslice, sslice)
  saxblocks$rectangles <- NULL

  blocks <- STSComputeBlocks(dataset, tslice, sslice)
  nrows = blocks$nrows
  ncols = blocks$ncols
  rectangles = blocks$rectangles
  blocks$rectangles <- NULL

  motifs<-list()
  size=length(blocks$datasets)
  for (i in 1:size) {
    block = blocks$datasets[[i]]
    saxblock = saxblocks$datasets[[i]]
    block = as.vector(as.matrix(block))
    saxblock = as.vector(as.matrix(saxblock))

    motifs[[i]] <- identifyMotifsInBlock(ts = block, tss = saxblock, sslice = sslice ,window.size = window_size, a = alpha)

    message(i, "/", size," - ", i/size*100,"%")
  }
  candidates = list(motifs = motifs, nrows = nrows, ncols = ncols, rectangles = rectangles, blocks = blocks, saxblocks = saxblocks)
  return (candidates)
}

#' Check of candidates
#'
#' First step in controlling the quality of candidates
#' @param candidates List which contains the discovered motifs
#' @param ka Support for Spatial Occurrence (SO)
#' @param si Support of Global Occurrence (GO)
#' @return Return a list of motifs [Isaxcode, Recmatrix, Vectst] which past the tresholds sigma and kappa, each motif contains :
#' @return Isaxcode : Motif sequences in character format
#' @return Recmatrix : Matrix giving as information the blocks containing this motif
#' @return Vectst : Coordinate of the start positions of the motif in the original dataset
#' @note To see how candidates are validated (see part 1) : \href{../inst/doc/check-candidates.html}{Validation of candidates}
#' @examples
#'#Generation of candidates
#'candidates <- STSIdentifyCandidateSTMotifs(STMotif::example_dataset, 10, 10, 7, 3)
#'
#'#Treatment of candidates
#'stmotifs <- STSIdentifySTMotifs(candidates, 1, 1)
#' @export
STSIdentifySTMotifs <- function(candidates, ka = 1, si = 1) {
  stmotifs <- list()
  blockOutput <- candidates$motifs
  nrows = candidates$nrows
  ncols = candidates$ncols
  rectangles <- candidates$rectangles
  blocks <- candidates$blocks
  saxblocks <- candidates$saxblocks
  for (i in 1:length(blockOutput)){
    stmotifs <- STSIdentifySTMotif(stmotifs, blockOutput[[i]], nrows, ncols, rectangles[[i]], ka = ka, si = si)
  }
  for (i in 1:length(stmotifs)) {
    stmotif = stmotifs[[i]]
    s = stmotif$vecs
    t = stmotif$vect
    stmotif$vecst = data.frame(s, t)
    stmotif$vecs <- NULL
    stmotif$vect <- NULL
    stmotifs[[i]] = stmotif
  }

  for (i in length(stmotifs):1){
    stmotif = stmotifs[[i]]
    if(length(unique(stmotif$vecst$s)) <= ka || length(stmotif$vecst$t) <= si){
      stmotifs[[i]] <- NULL
    }
  }
  return(stmotifs)
}


#' Filter the motifs
#'
#' Using the cluster to remove isolated motifs for each seleted motifs and after that recheck the global and spatial occurences.
#' @param stmotifs List of motifs which past the previous check
#' @param rectangles Information about the position of all blocks into the orginal dataset
#' @return Return a list of motifs [Isaxcode, Recmatrix, Vectst, SOtight, GOtight, SO, GO] which past the fliter, each motif contrains :
#' @return Isaxcode : Motif sequences in character format
#' @return Recmatrix : Matrix giving as information the blocks containing this motif
#' @return Vectst : Coordinate of the start positions of the motif in the original dataset
#' @return SOtight : Spatial Occurrences after removing outliers
#' @return GOtight : Global Occurrences after removing outliers
#' @return SO : Spatial Occurrences before removing outliers
#' @return GO : Global Occurrences before removing outliers
#' @note To see how motifs are filtered (see part 2) : \href{../inst/doc/check-candidates.html}{Motifs filtering}
#' @examples
#'#Generation of candidates
#'candidates <- STSIdentifyCandidateSTMotifs(STMotif::example_dataset, 10, 10, 7, 3)
#'
#'#Filter the selected motifs by their occurences
#'stmotifs <- STSIdentifySTMotifs(candidates, 1, 1)
#'
#'#Remove the outliers
#'sttightmotifs <- STSIdentifyTightSTMotifs(stmotifs, candidates$rectangles)
#' @export
STSIdentifyTightSTMotifs <- function(stmotifs, rectangles) {
  sttightmotifs <- list()
  sttightmotifsTemp  <- list()
  for(stmotif in (stmotifs)) {
    sttightmotifsTemp <- STSIdentifyTightSTMotif(stmotif, rectangles)
    for (item in (sttightmotifsTemp)) {
      pos = length(sttightmotifs)+1
      sttightmotifs[[pos]] <- item
      names(sttightmotifs)[pos] = item$isaxcod
    }
  }
  return (sttightmotifs)
}


#' Rank the motifs
#'
#' Rank motifs by their quality
#' @param sttightmotifs List of motifs selected motifs
#' @return The ranked version of the list provided as input.
#' @examples
#'#Generation of candidates
#'candidates <- STSIdentifyCandidateSTMotifs(STMotif::example_dataset, 10, 10, 7, 3)
#'
#'#Filter the selected motifs by their occurences
#'stmotifs <- STSIdentifySTMotifs(candidates, 1, 1)
#'
#'#Remove the outliers
#'sttightmotifs <- STSIdentifyTightSTMotifs(stmotifs, candidates$rectangles)
#'
#'#Ranking of the result
#'rank <- STSRankTightSTMotifs(sttightmotifs)
#' @export
STSRankTightSTMotifs <- function(sttightmotifs) {
  SOtight = rep(0, length(sttightmotifs))
  SO = rep(0, length(sttightmotifs))
  GOtight = rep(0, length(sttightmotifs))
  GO = rep(0, length(sttightmotifs))
  objfunc = rep(0, length(sttightmotifs))

  #Gathering of value for each TightSTmotifs
  ranktbl = data.frame(SOtight, SO, GOtight, GO)
  for (i in (1:length(sttightmotifs))) {
    motif = sttightmotifs[[i]]
    motif$objfunc = motif$SOtight/length(unique(sttightmotifs[[i]]$vecst[,2]))
    ranktbl$SOtight[i] = motif$SOtight
    ranktbl$SO[i] = motif$SO
    ranktbl$GOtight[i] = motif$GOtight
    ranktbl$GO[i] = motif$GO
    ranktbl$objfunc[i] = motif$objfunc
    sttightmotifs[[i]] = motif
  }
  #Order the id of each TightSTMotif by the SOtight and GOtight
  myrankorder = with(ranktbl, order(-SOtight,-GOtight))

  #Update of the stTightMotif list with new attributes and classify
  myranklist <- sttightmotifs[myrankorder]
  return(myranklist)
}
