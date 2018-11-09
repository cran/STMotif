#' Plot spatial series by selected columns
#'
#' Plot a selected motif and range of columns of the dataset
#' @param dataset Dataset containing numeric values
#' @param rmotifs List of ranked motifs
#' @param position Select by an integer a motif with his position
#' @param space Select a range of columns to plot the corresponding spatial series
#' @return Plot the spatial series
#' @import ggplot2
#' @import reshape2
#' @examples
#' #Launch all the workflow
#' #Plot the result
#' D  <- STMotif::example_dataset
#' DS <- NormSAX(STMotif::example_dataset,7)
#' stmotifs <- SearchSTMotifs(D,DS,3,7,10,10,3,7)
#' rstmotifs <- RankSTMotifs(stmotifs)
#' displayPlotSeries(dataset = D, rmotifs = rstmotifs ,position = 1 ,space = c(1,2,5:7))
#' @export
displayPlotSeries <- function(dataset, rmotifs, position, space){
  size_motif <- nchar(rmotifs[[1]]$isaxcod)
  namesCol <- paste("X",colnames(dataset),sep = "")
  data <- as.data.frame(dataset[,space])
  colnames(data) <- paste("X",colnames(dataset)[space], sep = "")
  data <- data.frame(x = 1:nrow(data),data)
  data <- reshape2::melt(data,id.vars = 1)
  data <- data.frame(data, color = "black")
  levels(data$color) <- c("black", "red")
  for(i in 1:length(rmotifs[[position]]$vecst$s)){
    if(rmotifs[[position]]$vecst$s[i]%in%space){
      data[data$variable==namesCol[rmotifs[[position]]$vecst$s[i]],][(rmotifs[[position]]$vecst$t[i]):(rmotifs[[position]]$vecst$t[i]+(size_motif-1)),4] <- "red"
    }
  }
  plot.series(data[1:nrow(data),])
}




#' Plot the intensity of values
#'
#' Display the intensity of values and higthlight one motif
#' @param dataset Dataset containing numeric values
#' @param rankList List of ranked motifs
#' @param alpha Number of letter used to do the encode
#' @return Pixelated dataset
#' @import ggplot2
#' @import reshape2
#' @import scales
#' @import RColorBrewer
#' @importFrom grDevices grey.colors
#' @examples
#' #Launch all the workflow
#' #Plot the result
#' D  <- STMotif::example_dataset
#' DS <- NormSAX(STMotif::example_dataset,7)
#' stmotifs <- SearchSTMotifs(D,DS,3,7,10,10,3,7)
#' rstmotifs <- RankSTMotifs(stmotifs)
#' intensityDataset(dataset = STMotif::example_dataset, rankList = rstmotifs,  7)
#' @export
intensityDataset <- function(dataset,rankList,alpha){
  colorEncode <- 1:alpha
  datasetColor.Org <- as.matrix(dataset)
  datasetColor.Org <- as.vector(datasetColor.Org)
  datasetColor.Org <- STSNormalization(datasetColor.Org)
  mybin <- binning(datasetColor.Org, alpha)
  datasetColor.Org <- colorEncode[mybin$bins_factor]
  datasetColor.Org <- t(matrix(datasetColor.Org, nrow = nrow(dataset), ncol = ncol(dataset)))
  datasetColor.Org <- melt(datasetColor.Org)
  datasetColor.Org$motif <- FALSE

  palhetaCores <- brewer.pal(length(rankList), 'Spectral')
  tam <- length(rankList)
  motifs.plot <-data.frame("s"=NULL, "t"=NULL, "g"= NULL)
  for (pos in 1:length(rankList)){
    motifs.plot<- rbind(motifs.plot ,data.frame("s"=rankList[[pos]]$vecst$s, "t"=rankList[[pos]]$vecst$t, "g"= pos, "color"=palhetaCores[pos]))
  }

  datasetColor <- merge(datasetColor.Org, motifs.plot, by.x=c('Var1', 'Var2'), by.y=c('s', 't'), all.x = TRUE)
  datasetColor$motif[!is.na(datasetColor$g)] <- TRUE
  datasetColor$g <- NULL
  datasetColor$color <- as.character(datasetColor$color)


  ggplot(data=datasetColor, aes(x=datasetColor$Var1, y=datasetColor$Var2, fill=datasetColor$value, color=datasetColor$color))   + geom_raster() +
    scale_fill_gradientn(colours = c("white","dimgrey"), values = scales::rescale(1:alpha), limits=c(1,alpha)) +
    theme_bw() + xlab("Space") + ylab("Time") + scale_y_reverse() +
    guides(fill=FALSE, color=FALSE) +
    geom_point(colour = ifelse(datasetColor$motif,datasetColor$color,NA), size = 1, show.legend = FALSE)

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
#' #Launch all the workflow
#' D  <- STMotif::example_dataset
#' DS <- NormSAX(STMotif::example_dataset,7)
#' stmotifs <- SearchSTMotifs(D,DS,3,7,10,10,3,7)
#' rstmotifs <- RankSTMotifs(stmotifs)
#' #Plot the result
#' top5motifs(dataset = STMotif::example_dataset, rankList = rstmotifs, alpha = 7)
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
    scale_size_manual(values=c("dot1"=1, "dot2"=1, "dot3"=1, "dot4"=1, "dot5"=1, "no_dot1"=NA,"no_dot2"=NA,"no_dot3"=NA,"no_dot4"=NA,"no_dot5"=NA), guide="none") +
    theme_bw() +
    ggtitle("") + xlab("") + ylab("") +
    guides(fill=guide_legend(title="")) + scale_y_reverse()
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
#' #Launch all the workflow
#' D  <- STMotif::example_dataset
#' DS <- NormSAX(STMotif::example_dataset,7)
#' stmotifs <- SearchSTMotifs(D,DS,3,7,10,10,3,7)
#' rstmotifs <- RankSTMotifs(stmotifs)
#' #Launch the process
#' runVisualization(dataset = STMotif::example_dataset, rstmotifs, 7)
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
