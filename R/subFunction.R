# Adjust a block
# Adjust the dimensions of a dataset to build the blocks
STSAdjustBlock  <- function(dataset, tslice, sslice) {
  c = ncol(dataset)
  r = nrow(dataset)
  ec = c %% sslice
  er = r %% tslice
  dataset = dataset[1:(r-er), 1:(c-ec)]
  return (dataset)
}

# Normalize the data
# Normalize the data using z-score
STSNormalization <- function (vector){
  return ((vector-mean(vector, na.rm = T))/sd(vector, na.rm = T))
}

# binning the dataset
# Build an encode for the values
binning <- function(v, alpha) {
  p <- seq(from = 0, to = 1, by = 1/alpha)
  q <- quantile(v, p)
  qf <- matrix(c(q[1:(length(q)-1)],q[2:(length(q))]), ncol=2)
  vp <- cut(v, q, FALSE, include.lowest=TRUE)
  m <- tapply(v, vp, mean)
  vm <- m[vp]
  mse <- mean( (v - vm)^2, na.rm = TRUE)
  return (list(binning=m, bins_factor=vp, q=q, qf=qf, bins=vm, mse=mse))
}

# Encode values
# Encode numeric values from a vector
STSSaxEncode <- function(dataset, vector, alpha) {
  mybin <- binning(vector, alpha)
  myletters <- letters[1:alpha]
  saxvector <- myletters[mybin$bins_factor]
  saxvector = matrix(saxvector, nrow = nrow(dataset), ncol = ncol(dataset))
  saxvector = data.frame(saxvector)
  colnames(saxvector) =  colnames(dataset)
  return(saxvector)
}

# Compute blocks
# Create blocks from the original dataset
STSComputeBlocks <- function(dataset, tslice, sslice) {
  datasets <- list()
  rectangles <- list()

  c = ncol(dataset)
  r = nrow(dataset)
  nc = c / sslice
  nr = r / tslice
  i = 1
  j = 1
  n = 1
  for (i in 1:nc) {
    sc = (i-1)*sslice + 1
    ec = (i)*sslice
    for (j in 1:nr) {
      sr = (j-1)*tslice + 1
      er = (j)*tslice
      ds = dataset[sr:er, sc:ec]
      datasets[[n]] = ds
      rect = c(sS = sc, eS = ec, sT = sr, eT = er, nr = j, nc = i)
      rectangles[[n]] = rect
      n = n + 1
    }
  }
  blocks = list(datasets = datasets, nrows = nr, ncols = nc, rectangles = rectangles)
  return (blocks)
}


# Generate candidates from a block
# Take a block and discover the frequent candidates
identifyMotifsInBlock <- function(ts, tss, window.size, sslice , a, overlap = 0) {
  #Normalization
  ts.nor <- STSNormalization(ts)

  #Calculate the increases value to the index (If overlap == 0): window.size
  b <- ifelse(overlap==1, yes = 1, no = round((1-overlap)*window.size, digits = 0))

  #Create a matrix with all the subsequences
  ts.subs <- foreach::foreach(i=seq(from = 1, to = length(ts), by = b), .combine = rbind) %do% {
    if(trunc(i/sslice)==trunc((i+window.size-1)/sslice)){ #Check if it's a fake motif
      c(i,subs.temp <- ts.nor[i:(i+window.size-1)])
    }
  }
  ts.subs <- na.omit(ts.subs)

  #Create a matrix with all the SAX subsequences
  ts.sax <- foreach::foreach(i=seq(from = 1, to = length(tss), by = b), .combine = rbind) %do% {
    if(trunc(i/sslice)==trunc((i+window.size-1)/sslice)){ #Check if it's a fake motif
      c(i,sax.temp <- tss[i:(i+window.size-1)])
    }
  }
  ts.sax <- na.omit(ts.sax)
  ts.sax <- as.data.frame(ts.sax, stringsAsFactors = FALSE)

  #give a name for each column. TD nmaybe give a name is not necessary
  colnames(ts.sax) <- c("StartPosition", 1:window.size)
  ts.sax$StartPosition <- as.numeric(ts.sax$StartPosition)

  #Creating a candidate list with a list of starpPosition of the same motifs
  i = j <- 1
  indices <- list()
  for (i in 1:nrow(ts.sax)){
    saxCandidate <- paste(ts.sax[i,-1], collapse = "")
    indices[[saxCandidate]] <- c(indices[[saxCandidate]],ts.sax[i,1])
  }
  while (j <= length(indices)){ #removing the candidate with just 1 or less candidate
    if(length(indices[[j]])<=1){indices[[j]]<-NULL}else{j<-j+1}
  }

  if(length(indices)>0){ #Check if there is repeated candidate found into the block
    #Create the output motif.sax
    #Each identical sequence is grouping to create a sub matrix of ts.sax
    motif.sax <- foreach::foreach(i = 1:length(indices)) %do% {
      ts.sax[which(ts.sax[,1] %in% indices[[i]]),]
    }
    #Create the output motif.raw
    #Each identical sequence is grouping to create a sub matrix of ts.sub
    motif.raw <- foreach::foreach(i = 1:length(indices)) %do% {
      ts.subs[which(ts.subs[,1] %in% indices[[i]]),]
    }
  }else{
    motif.sax <- NULL
    motif.raw <- NULL
  }

  return(list(Subs=ts.subs, Subs.SAX=ts.sax, Motif.raw=motif.raw, Motif.SAX=motif.sax, Indices=indices))
}


# Handle candidates
#
# Handle candidates from one block
STSIdentifySTMotif <- function(stmotifs, motif, nrows, ncols, rectangle, ka, si) {
  k <- length(stmotifs)

  #Get propreties of the block handled
  sS = rectangle["sS"]
  eS = rectangle["eS"]
  sT = rectangle["sT"]
  eT = rectangle["eT"]
  nr = rectangle["nr"]
  nc = rectangle["nc"]

  recMatrix = matrix(rep(0, nrows*ncols), nrow = nrows, ncol = ncols)
  tslice <- eT - sT + 1
  sslice <- eS - sS + 1
  #for candidate motif discoverd inside the block
  if(length(motif$Indices)>0){ #Check if there is repeated candidate found into the block
    for(a in 1:length(motif$Indices)){
      #vectorize the indices of the candidate
      vec <- motif$Indices[[a]]

      #check if the GO(motif[a] >= sigma)
      if(length(vec) >= si) {
        #scount: vector of 0 with the slice columns
        scount <- rep(0, sslice)

        #for each occurence of the candidate
        for(z in 1: length(vec)) {
          #mark each column wich contains the candidate
          i <- as.integer(vec[z] / tslice) + 1
          scount[i] <- 1
        }
        #count SO
        count <- sum(scount)


        #check if the SO(motif[a] >= kapa)
        if(count >= ka) {
          #take the SAX of the candidate
          saxcod <- motif$Motif.SAX[[a]][1,2:(length(motif$Subs.SAX))]
          isaxcod <- paste(saxcod, collapse = "")

          vect <- as.integer(vec %% tslice) + sT

          vecs <- as.integer(vec / tslice) + sS
          i <- match(isaxcod, names(stmotifs))
          if (is.na(i)) {
            k = k + 1
            stmotifs[[k]] <- list(saxcod=saxcod, isaxcod=isaxcod, vecs=vecs, vect=vect, recmatrix=recMatrix)
            stmotifs[[k]]$recmatrix[nr, nc] = 1
            names(stmotifs)[k] = isaxcod
          }
          else {
            list <- stmotifs[[i]]
            list$recmatrix[nr, nc] = max(list$recmatrix)+1
            list$vect <- c(list$vect, vect)
            list$vecs <- c(list$vecs, vecs)
            stmotifs[[i]] <- list
          }
        }
      }
    }
  }
  return (stmotifs)
}

# Handle one motif
# Remove the isolated motifs
STSIdentifyTightSTMotif <- function(stmotif, rectangles) {
  #We selected one motif with its information
  tight <- list()
  mat <- stmotif$recmatrix #Get the recmatrix of one motif
  vecst <- stmotif$vecst #Get start position of the motif
  #For each block
  for (i in 1:(nrow(mat)-1)) {
    for (j in 1:(ncol(mat)-1)) {
      #Checking blocks neighbor if there is a presence of this motif
      if (mat[i,j] != 0) {
        iP <- i + 1
        jP <- j + 1
        if ((iP <= nrow(mat)) && (mat[iP,j] != 0)) {
          k <- min(mat[iP,j], mat[i,j])
          mat[mat == mat[iP,j] | mat == mat[i,j]] = k
        }
        if ((jP <= ncol(mat)) && (mat[i,jP] != 0)) {
          k <- min(mat[i,jP], mat[i,j])
          mat[mat == mat[i,jP] | mat == mat[i,j]] = k
        }
        if ((iP <= nrow(mat)) && (mat[iP,j] != 0) && (jP <= ncol(mat)) && (mat[i,jP] != 0)) {
          k <- min(mat[iP,jP], mat[i,j])
          mat[mat == mat[iP,jP] | mat == mat[i,j]] = k
        }
      }
    }
  }
  vec <- as.vector(mat)
  vec <- vec[vec > 0]
  vec <- unique(vec)
  k <- 1
  SO = length(unique(vecst$s))
  GO = nrow(vecst)
  for (i in (vec)) {
    stmotif$recmatrix[mat != i] <- 0
    stmotif$recmatrix[mat == i] <- k
    vecrects <- as.vector(stmotif$recmatrix)
    #Get position of each blocks whitch contains this motif
    rects <- rectangles[vecrects>0]
    stmotif$vecst <- vecst
    conds = rep(FALSE, nrow(stmotif$vecst))
    for (rect in (rects)) {
      sS = rect["sS"]
      eS = rect["eS"]
      sT = rect["sT"]
      eT = rect["eT"]
      conds = conds | (stmotif$vecst$s >= sS & stmotif$vecst$s <= eS & stmotif$vecst$t >= sT & stmotif$vecst$t <= eT)
    }
    stmotif$vecst <- stmotif$vecst[conds,]
    stmotif$SOtight = length(unique(stmotif$vecst$s))
    stmotif$GOtight = nrow(stmotif$vecst)
    stmotif$SO = SO
    stmotif$GO = GO
    GO = nrow(vecst)

    tight[[k]] <- stmotif
    k <- k + 1
  }
  return(tight)
}


# Function to plot spatial series
plot.series <- function(series, label_series = "", label_x = "", label_y = "") {
  grf <- ggplot(data=series, ggplot2::aes(x = series$x, y = series$value, colour = series$color, group = 1))
  grf <- grf + scale_colour_identity(series$color) + geom_line() + geom_point(data=series, aes(x = series$x, y = series$value), size=0.5) + facet_grid(variable ~ .)
  grf <- grf + xlab(label_x)
  grf <- grf + ylab(label_y)
  grf <- grf + ggtitle("Motifs in spatial-time series")
  grf <- grf + theme_bw(base_size = 10)
  grf <- grf + theme(panel.grid.major = element_blank()) + theme(panel.grid.minor = element_blank())
  grf <- grf + theme(axis.text.x = element_blank(), axis.text.y = element_blank(), axis.ticks = element_blank())
  return(grf)
}
