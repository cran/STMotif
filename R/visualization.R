#' Plot spatial series by selected columns
#'
#' Plot a selected motif and range of columns of the dataset
#' @param dataset Dataset containing numeric values
#' @param rankList List of ranked motifs
#' @param position Select by an integer a motif with his position
#' @param space Select a range of columns to plot the corresponding spatial series
#' @return Plot the spatial series
#' @import ggplot2
#' @import reshape2
#' @examples
#'#Launch all the workflow
#'candidates <- STSIdentifyCandidateSTMotifs(STMotif::example_dataset, 10, 10, 7, 3)
#'stmotifs <- STSIdentifySTMotifs(candidates, 1, 1)
#'sttightmotifs <- STSIdentifyTightSTMotifs(stmotifs, candidates$rectangles)
#'rankResult <- STSRankTightSTMotifs(sttightmotifs)
#'
#'#Plot the result
#'displayPlotSeries(dataset = STMotif::example_dataset, rankResult, 1, 2:11)
#' @export
displayPlotSeries <- function(dataset, rankList, position, space){
  size_motif <- length(rankList[[1]]$saxcod)
  namesCol <- paste("X",colnames(dataset),sep = "")
  data <- as.data.frame(dataset[,space])
  colnames(data) <- paste("X",colnames(dataset)[space], sep = "")
  data <- data.frame(x = 1:nrow(data),data)
  data <- reshape2::melt(data,id.vars = 1)
  data <- data.frame(data, color = "black")
  levels(data$color) <- c("black", "red")
  for(i in 1:length(rankList[[position]]$vecst$s)){
    if(rankList[[position]]$vecst$s[i]%in%space){
      data[data$variable==namesCol[rankList[[position]]$vecst$s[i]],][(rankList[[position]]$vecst$t[i]):(rankList[[position]]$vecst$t[i]+(size_motif-1)),4] <- "red"
    }
  }
  plot.series(data[1:nrow(data),])
}


#' Plot the intensity of values
#'
#' Display the intensity of values and higthlight one motif
#' @param dataset Dataset containing numeric values
#' @param rankList List of ranked motifs
#' @param position Select by an integer a motif with his position
#' @param alpha Number of letter used to do the encode
#' @return Pixelated dataset
#' @import ggplot2
#' @import reshape2
#' @import scales
#' @importFrom grDevices grey.colors
#' @examples
#'#Launch all the workflow
#'candidates <- STSIdentifyCandidateSTMotifs(STMotif::example_dataset, 10, 10, 7, 3)
#'stmotifs <- STSIdentifySTMotifs(candidates, 1, 1)
#'sttightmotifs <- STSIdentifyTightSTMotifs(stmotifs, candidates$rectangles)
#'rankResult <- STSRankTightSTMotifs(sttightmotifs)
#'
#'#Plot the result
#'intensityDataset(dataset = STMotif::example_dataset, rankResult, 1, 7)
#' @export
intensityDataset <- function(dataset,rankList,position,alpha){
  colorEncode <- 1:alpha
  datasetColor <- as.matrix(dataset)
  datasetColor <- as.vector(datasetColor)
  datasetColor <- STSNormalization(datasetColor)
  mybin <- binning(datasetColor, alpha)
  datasetColor <- colorEncode[mybin$bins_factor]
  datasetColor <- t(matrix(datasetColor, nrow = nrow(dataset), ncol = ncol(dataset)))
  datasetColor <- melt(datasetColor)
  datasetColor$motif <- FALSE
  k <- 1
  line <- NULL
  while(k<length(rankList[[position]]$vecst$t)){
    line[k] <- (rankList[[position]]$vecst$t[k]*length(dataset[1,])-length(dataset[1,]))+(rankList[[position]]$vecst$s[k])
    k <- k + 1
  }
  datasetColor[line,]$motif <- TRUE
  ggplot(data=datasetColor, aes(x=datasetColor$Var1, y=datasetColor$Var2, fill=datasetColor$value)) +
    geom_raster() +
    scale_fill_gradientn(colours = c("white","dimgrey"), values = scales::rescale(1:alpha), limits=c(1,alpha))+
    geom_point(aes(size=ifelse(datasetColor$motif, "dot", "no_dot"),colour=ifelse(datasetColor$motif, "dot", "no_dot")), na.rm=T) +
    scale_size_manual(values=c("dot"=1, "no_dot"=NA), guide="none") +
    scale_colour_manual(name="", values = c("dot"="red", "no_dot"=NA),labels = c(paste("Motif",position),"")) +
    theme_bw() +
    ggtitle("Intensity of values in the dataset") + xlab("Space") + ylab("Time") +
    guides(fill=guide_legend(title="SAX Encode")) + scale_y_reverse()
}


#' Plot the 5 motifs
#'
#' Display the intensity of values and higthlight the top five motifs
#' @param dataset Dataset containing numeric values
#' @param rankList List of ranked motifs
#' @param alpha Number of letters used to do the encode
#' @return Pixelated dataset
#' @import ggplot2
#' @import reshape2
#' @import scales
#' @importFrom grDevices grey.colors
#' @examples
#'#Launch all the workflow
#'candidates <- STSIdentifyCandidateSTMotifs(STMotif::example_dataset, 10, 10, 7, 3)
#'stmotifs <- STSIdentifySTMotifs(candidates, 1, 1)
#'sttightmotifs <- STSIdentifyTightSTMotifs(stmotifs, candidates$rectangles)
#'rankResult <- STSRankTightSTMotifs(sttightmotifs)
#'
#'#Plot the result
#'top5motifs(dataset = STMotif::example_dataset, rankList = rankResult, alpha = 7)
#' @export
top5motifs <- function(dataset,rankList,alpha){
  colorEncode <- 1:alpha
  datasetColor <- as.matrix(dataset)
  datasetColor <- as.vector(datasetColor)
  datasetColor <- STSNormalization(datasetColor)
  mybin <- binning(datasetColor, alpha)
  datasetColor <- colorEncode[mybin$bins_factor]
  datasetColor <- t(matrix(datasetColor, nrow = nrow(dataset), ncol = ncol(dataset)))
  datasetColor <- melt(datasetColor)
  rank <- rankList[1:5]
  datasetColor$motif <- 0
  for(j in 1:length(rank)){
    l <- 1
    line <- NULL
    while(l<length(rank[[j]]$vecst$t)){
      line[l] <- (rank[[j]]$vecst$t[l]*length(dataset[1,])-length(dataset[1,]))+(rank[[j]]$vecst$s[l])
      l <- l + 1
    }
    datasetColor[line,]$motif <- j
  }
  ggplot(data=datasetColor, aes(x=datasetColor$Var1, y=datasetColor$Var2, fill=datasetColor$value)) +
    geom_raster() +
    scale_fill_gradientn(colours = c("white","dimgrey"), values = scales::rescale(1:alpha), limits=c(1,alpha))+
    geom_point(aes(size=ifelse(datasetColor$motif == 1, "dot1", "no_dot1"),colour=ifelse(datasetColor$motif == 1, "dot1", "no_dot1")), na.rm=T) +
    geom_point(aes(size=ifelse(datasetColor$motif == 2, "dot2", "no_dot2"),colour=ifelse(datasetColor$motif == 2, "dot2", "no_dot2")), na.rm=T) +
    geom_point(aes(size=ifelse(datasetColor$motif == 3, "dot3", "no_dot3"),colour=ifelse(datasetColor$motif == 3, "dot3", "no_dot3")), na.rm=T) +
    geom_point(aes(size=ifelse(datasetColor$motif == 4, "dot4", "no_dot4"),colour=ifelse(datasetColor$motif == 4, "dot4", "no_dot4")), na.rm=T) +
    geom_point(aes(size=ifelse(datasetColor$motif == 5, "dot5", "no_dot5"),colour=ifelse(datasetColor$motif == 5, "dot5", "no_dot5")), na.rm=T) +
    scale_colour_manual(name="", values = c("dot1"="yellow", "dot2"="red", "dot3"="orange", "dot4"="blue", "dot5"="green", "no_dot1"=NA,"no_dot2"=NA,"no_dot3"=NA,"no_dot4"=NA,"no_dot5"=NA),labels = c('Motif1','Motif2','Motif3','Motif4','Motif5','','','','',''))+
    scale_size_manual(values=c("dot1"=1, "dot2"=1, "dot3"=1, "dot4"=1, "dot5"=1, "no_dot1"=NA,"no_dot2"=NA,"no_dot3"=NA,"no_dot4"=NA,"no_dot5"=NA), guide="none") +
    theme_bw() +
    ggtitle("Intensity of values in the dataset") + xlab("Space") + ylab("Time") +
    guides(fill=guide_legend(title="SAX Encode")) + scale_y_reverse()
}


#' Interactive visualization
#'
#' Launch a process to have an interactive visualization
#' @param dataset Dataset containing numeric values
#' @param rankList List of ranked motifs
#' @param alpha Number of letters used to do the encode
#' @import shiny
#' @import utils
#' @examples
#'\dontrun{
#'#Launch all the workflow
#'candidates <- STSIdentifyCandidateSTMotifs(STMotif::example_dataset, 10, 10, 7, 3)
#'stmotifs <- STSIdentifySTMotifs(candidates, 1, 1)
#'sttightmotifs <- STSIdentifyTightSTMotifs(stmotifs, candidates$rectangles)
#'rankResult <- STSRankTightSTMotifs(sttightmotifs)
#'
#'#Launch the process
#'runVisualization(dataset = STMotif::example_dataset, rankResult, 7)
#'}
#' @export
runVisualization <- function(dataset, rankList, alpha){
  pos <- 1
  envir = as.environment(pos)
  assign("datasetShiny", dataset, envir = envir)
  assign("listMotifShiny", rankList, envir = envir)
  assign("alphaShiny", alpha, envir = envir)
  appDir <- system.file("shiny-visualization", package = "STMotif")
  if (appDir == "") {
    stop("Could not find directory. Try re-installing `STMotif`.", call. = FALSE)
  }

  shiny::runApp(appDir, display.mode = "normal")
}
