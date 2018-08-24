#'  Generate candidates
#'
#' Find frequent candidates in a dataset
#' @param data Dataset containing numeric values
#' @param tslice "Time slice" Number of rows in each block
#' @param sslice "Space slice" Number of columns in each block
#' @param alpha Number of letters to do the encode
#' @param window_size Length of motifs
#' @return List [Motifs, Nrows, Ncols, Rectangles, Blocks, Saxblocks] containing all information about the blocks: number of blocks in horizontal and vertical, values of each block (numeric and SAX encoded) and information about frequent candidates in blocks
#' @return Motifs: Contain for each block the Subs: dataframe containing all the subsequences in the original dataframe; Subs.SAX: a dataframe containing all the subsequences in SAX representations; Motif.raw: a list showing the candidates discovered in the original dataframe; Motif.SAX: a list showing the candidates discovered in SAX representations and Indices: a list showing the starting positions of subsequences for each candidate discovered.
#' @return Nrows and Ncols: Number of blocks in horizontal and vertical
#' @return Rectangles: For each block gives information about the position of the block into the original dataset
#' @return Blocks: Contains the values of each block
#' @return Saxblocks: Contains the values encoded by SAX of each block
#' @importFrom stats quantile sd na.omit
#' @import foreach
#' @note To see all the stages of the generation of candidates: \href{../inst/doc/generation-of-candidates.html}{Generation of candidates}
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
#' First step in controlling the frequency of candidates
#' @param candidates List which contains the discovered candidates
#' @param ka Support for Spatial Occurrence (SO)
#' @param si Support of Global Occurrence (GO)
#' @return Return a list of motifs [Isaxcode, Recmatrix, Vectst] which past the tresholds sigma and kappa. Each motif contains:
#' @return Isaxcode: Motif sequences in character format
#' @return Recmatrix: Matrix giving as information the blocks containing this motif
#' @return Vectst: Coordinate of the start positions of the motif in the original dataset
#' @note To see how candidates are validated (part 1): \href{../inst/doc/validate-candidates.html}{Validation of candidates}
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
#' Using the cluster to remove isolated motifs for each seleted motifs and after that, recheck the global and spatial occurences.
#' @param stmotifs List of motifs which past the previous check
#' @param rectangles Information about the position of all blocks into the orginal dataset
#' @return Return a list of motifs [Isaxcode, Recmatrix, Vectst, SOtight, GOtight, SO, GO] which past the fliter, each motif contrains:
#' @return Isaxcode: Motif sequences in character format
#' @return Recmatrix: Matrix giving as information the blocks containing this motif
#' @return Vectst: Coordinate of the start positions of the motif in the original dataset
#' @return SO: Spatial Occurrences before removing isolated motifs
#' @return GO: Global Occurrences before removing isolated motifs
#' @return SOtight: Spatial Occurrences after removing isolated motifs
#' @return GOtight: Global Occurrences after removing isolated motifs
#' @note To see how motifs are filtered (part 2): \href{../inst/doc/validate-candidates.html}{Filter motifs}
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
#' @param sttightmotifs List of motifs after removing isolated motifs
#' @return The ranked version of the list provided as input.
#' @examples
#'#Generation of candidates
#'candidates <- STSIdentifyCandidateSTMotifs(STMotif::example_dataset, 10, 10, 7, 3)
#'
#'#Filter the selected motifs by their occurences
#'stmotifs <- STSIdentifySTMotifs(candidates, 1, 1)
#'
#'#Remove the isolated motifs
#'sttightmotifs <- STSIdentifyTightSTMotifs(stmotifs, candidates$rectangles)
#'
#'#Order the result list
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
