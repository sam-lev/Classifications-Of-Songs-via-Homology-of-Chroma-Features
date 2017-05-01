" This work was created by Samuel Leventhal and Zahra Fahimar and 
licensed under a Creative Commons Attribution-NonCommercial-NoDerivatives
4.0 International License. The following uses chroma features of musical
compositions to extract homological features through persistence allowing
song classifaction through topological metrics such as bottleneck distance 
and p-Wasserstein via clustering and multi-dimensional scaling"



musicTopology <- function() {
  
  # must download h5 package for audio files
  if (!require(package = "h5")) {
    install.packages("h5")
  }
  library("h5")
  
  if (!require(package = "TDA")) {
    install.packages("TDA")
  }
  library("TDA")
  
  if (!require(package = "signal")) {
    install.packages("signal")
  }
  library("signal")
  
  if(!require(package = "stats")){
    install.packages("stats")
  }
  library("stats")
  
  if(!require(package = "rPython")){
    install.packages("rPython")
  }
  library("rPython")
  # connect to the sqlite file
  song = h5file("./MillionSongSubset/data/B/B/B/TRBBBFO128F931535D.h5", mode = "a")
  song2 = h5file("./MillionSongSubset/data/A/R/R/TRARRWA128F42A0195.h5", mode = "a")
  song3 = h5file("./MillionSongSubset/data/A/R/R/TRARRYC128F428CCDA.h5", mode = "a")

  # get the chroma feature as a data.frame
  chroma = song["/analysis/segments_pitches"] 
  sr = song["/analysis/sample_rate"]
  
  # second example song features 
  chroma2 = song2["/analysis/segments_pitches"] 
  sr2 = song2["/analysis/sample_rate"]
  
  # third example song features 
  chroma3 = song3["/analysis/segments_pitches"] 
  sr3 = song3["/analysis/sample_rate"]
  
  # Attributes held in songs .H5 file.
  attr <- list.datasets(song)
  print(attr)
  print("######")
  print("all attributes from song file")
  print(attr)
  print("#######")
  print("#######")
  print("Chroma features")
  print(chroma)
  
  # Calculate the chroma matrix.  Use a long FFT to discriminate
  # spectral lines as well as possible (2048 is the default value)
  #
  # Establish the sample rate used for current song
  # This is not the tue sample rate.
  # I can not find the field for sampling rate
  sr = 2# 45*1000#sample rate
  cfftlen = 2048
  hz = chroma[,1]
  tone = chroma[1,]
  ntone = ncol(chroma)
  nhz = nrow(chroma)
  
  # same parameters as above for song 2
  hz2 = chroma2[,1]
  tone2 = chroma2[1,]
  ntone2 = ncol(chroma2)
  nhz2 = nrow(chroma2)
  
  # same parameters as above for song 2
  hz3 = chroma3[,1]
  tone3 = chroma3[1,]
  ntone3 = ncol(chroma3)
  nhz3 = nrow(chroma3)
  
  # The frame advance is always one quarter of the FFT length.  Thus,
  # the columns  of chroma are at timebase of fftlen/4/sr
  t = 1:nhz*cfftlen/4/sr 
  t2 = 1:nhz2*cfftlen/4/sr
  t3 = 1:nhz3*(cfftlen/4/sr)
  
  # Plot the chroma features with resoect to time and Hz value
  chromMatrix = matrix(hz, t)#20*log10(hz))
  chromMatrix2 = matrix(hz2, t2)
  chromMatrix3 = matrix(hz3, t3)
  #plot(chromMatrix)

  
  # Plot the spectragram of the chroma features. Assuming sampling rate of 2
  # Note: you must comment out other plots to see this.
  specgram(hz, n = min(256, nhz), Fs = 2)

  
  # Compute the persistent homology of the observed chromatic features for
  # a song with respect to the observed tone's time signiture. Persistence 
  # in this case is calculated using the nesting of Vietoris-Rips complexes.
  persChroma <- ripsDiag(X = chromMatrix, maxdimension = 1, maxscale = 1, library="GUDHI", location = TRUE, printProgress=FALSE)$diagram
  
  pdf(paste(getwd(),"/plots/persChroma.pdf",sep=""),"pdf")
  plot(persChroma)
  dev.off()
  
  persChroma2 <- ripsDiag(X = chromMatrix2, maxdimension = 1, maxscale = 1, library="GUDHI", location = TRUE, printProgress=FALSE)$diagram
  
  pdf(paste(getwd(),"/plots/persChroma2.pdf",sep=""),"pdf")
  plot(persChroma2)  
  dev.off()
  
  persChroma3 <- ripsDiag(X = chromMatrix3, maxdimension = 1, maxscale = 1, library="GUDHI", location = TRUE, printProgress=FALSE)$diagram
  
  pdf(paste(getwd(),"/plots/perschroma3.pdf",sep=""),"pdf")
  plot(persChroma3)
  dev.off()


  # Calculate the q-Wasserstein distance between persistence diagrams
  # found for homologies of chomatic feature / time signature data. 
  # To determine the optimal degree q for computing the  q-Wasserstein
  # we iterate over q and plot the set of Wasserstein distances
  qWassersteiDistances = 1:10
  for(q in 1:10){
    qW = wasserstein(persChroma, persChroma2, p=q, dimension=1)
    qWassersteiDistances[q] = qW
  }
  # Here we compare two songs which are more similar (I assume-they were
  # taken from the same folder: MillionSongSubset/data/R/R
  qWassersteiDistances_similar = 1:10
  for(q in 1:10){
    qW = wasserstein(persChroma2, persChroma3, p=q, dimension=1)
    qWassersteiDistances_similar[q] = qW
  }
  print("")
  print("Maximal and Minial q-Wasserstein distance of two songs from iterating over q")
  print(max(qWassersteiDistances))
  print(min(qWassersteiDistances))
  print("")
  print("Now for the similar songs")
  print("Maximal and Minial q-Wasserstein distance of two songs from iterating over q")
  print(max(qWassersteiDistances_similar))
  print(min(qWassersteiDistances_similar))
  
 h5close(song)
 h5close(song2)
 h5close(song3)
 
 
 # Here we list all h5 files in a folder, calculate the persistant homology for each
 # song and add the corresponding diagram to 'persChromaSongs' and each songs chromatic
 # features to 'chromaFeatures'
  " Write persistence diagrams of all song's h5 files in
    data folders to csv files"
  #writeH5Persistence()
  
  
  
  ""
  " form list of all persistence diagrams."
  ""
  persChromaSongs <- readPersistenceCSV(paste(getwd(),"/output/persDiag/",sep=""))
  


   # Obtain the index of the persistance diagram affording the 
   # minimal bottleneck distance other than song itself. The index
   # correlating to minimal bottleneck distance also correlates to 
   # the song most similar in persistant homology and by ouy hypothsis
 # similar in musical composition.
 minBtl = min(bottleneckSongs[(2:length(bottleneckSongs))])
 # index of song nearest in persistence to first song
 simSong = match(minBtl, bottleneckSongs)
 # plot the persistence diagram of the song most similar to 
 # persitent homology of the first son
 pdf(paste(getwd(),"/plots/persexample.pdf",sep=""),"pdf")
 plot(persChromaSongs[[simSong]])
 dev.off()
 # Need to figure out how to find the attribute giving sonf title
 # I want to hear it!
 #print(h5attr(songData[[simSong]],"TITLE"))
 
 
} #End musicTopology()


pairwiseBottleneck <- function(fromDirPath, matrixData, persDiags, write){
  " Function computes the pairwise bottleneck distance between 
  persistence diagrams. If 'fromDirPath is present this function uses 
  readPersistenceCSV(fromDirPath) function's output, i.e. reads from csv 
  file 'pairwisePersistence.csv' to construct the matrix of pairwise 
  persistence diagrams between songs 
  If write is present it must be a string consisting of the desired 
  file name to write the matrix of pairwise bottleneck distances to.

  e.g. pairwiseBottleneck(fromDirPath = paste(getwd(), '/output/persDiag/A/A/A/',sep='') )
  "
  if (!require(package = "TDA")) {
    install.packages("TDA")
  }
  library("TDA")
  # computes pairwise bottleneck distance between
  # all diagram elements stored to file whose
  # pathnames are in the list 'fromDirPath'
  if(!missing(fromDirPath)){
    # Fill list of persistence diagrams from csv files
    # under path 'fromDirPath'
    persistenceData = readPersistenceCSV(dirPath = fromDirPath, songNames= TRUE)
    songNames = persistenceData$songNames
    persistenceDiagrams = persistenceData$diagrams
    # Compute all pairwise persistence diagrams and store
    # to matrix.
    bottleneckMatrix= matrix(data=NA, nrow=length(persistenceDiagrams), ncol=length(persistenceDiagrams))
    for(q in 1:length(persistenceDiagrams)){
      for(p in 1:length(persistenceDiagrams)){
        bdist = bottleneck(persistenceDiagrams[[q]], persistenceDiagrams[[p]], dimension = 1)
        bottleneckMatrix[q,p] = bdist
      }
    }
    # Label the bottleneck matrix X and Y axis to be that 
    # of the .h5 files i.e. (foo, bar) corresponds to 
    # the index of bottleneck distances from persistence diagrams
    # of foo.h5 and bar.h5
    rownames(bottleneckMatrix) <- paste(songNames)
    colnames(bottleneckMatrix) <- paste(songNames)
    
    if(!missing(write)){
      write.csv(bottleneckMatrix, row.names = TRUE, file = paste(getwd(),"/output/",write,".csv",sep=""))
    }
    return(bottleneckMatrix)
  }
  # computes pairwise bottleneck distance between
  # all diagram elements in diag
  if(!missing(persDiags)){
    bottleneckMatrix= matrix(data=NA, nrow=length(persDiags), ncol=length(persDiags))
    for(q in 1:length(persDiags)){
      for(p in 1:length(persDiags)){
        bottleneckMatrix[q,p] = bottleneck(persDiags[[q]], persDiags[[p]], dimension = 1)
      }
    }
    if(!missing(write)){
      write.csv(bottleneckMatrix, row.names = TRUE, file = paste(getwd(),"/output/",write,".csv",sep=""))
    }
    return(bottleneckMatrix)
  }
  # computes pairwise bottleneck distance between
  # all matrix data elements in matrixData, i.e.
  # matrixData assumed to be point cloud.
  # matrixData is of the form [ [dataset 1], [dataset 2], ...]
  if(!missing(matrixData)){
    persistenceDiagrams= matrix(data=NA, nrow=length(matrixData), ncol=length(matrixData))
    for( pc in 1:length(matrixData)){
      persistenceDiagrams[[pc]] = ripsDiag(matrixData[[pc]], maxdimension = 1, maxscale = 1, library="GUDHI", location = TRUE, printProgress=FALSE)$diagram
    }
    bottleneckMatrix= matrix(data=NA, nrow=length(persistenceDiagrams), ncol=length(persistenceDiagrams))
    for(q in 1:length(persistenceDiagrams)){
      for(p in 1:length(persistenceDiagrams)){
        bottleneckMatrix[q,p] = bottleneck(persistenceDiagrams[[q]], persistenceDiagrams[[p]], dimension = 1)
      }
    }
    if(!missing(write)){
      write.csv(bottleneckMatrix, row.names = FALSE, file = paste(getwd(),"/output/",write,".csv",sep=""))
    }
    return(bottleneckMatrix)
  }
  
}# end pairwiseBottleneck()



writeH5Persistence <- function(  ){
  "Computes persistence 
  diagrams of from all h5songs chroma features in MillionSongSubset 
  and writes persistence diagrams to csv files under same folder
  architecture in /output/perDiag/ as seen for h5 files. i.e.
  if foobar.h5 is in MillionSongSubset/data/A/G/H/foobar.h5 then
  function writes the persistence diagram of foobar.h5's chroma 
  features to /output/persDiag/A/G/H/foobar.csv
  "
  if (!require(package = "h5")) {
    install.packages("h5")
  }
  library("h5")
  
  if (!require(package = "TDA")) {
    install.packages("TDA")
  }
  library("TDA")
  
  dirName <- "./MillionSongSubset/data/"
  h5Files = list.files(dirName , pattern = ".h5",recursive=TRUE, all.files = TRUE)
  chromaFeatures = vector("list", length(h5Files)) #rep(NA, length(h5Files))
  persChromaSongs = vector("list", length(h5Files)) #rep(NA, length(h5Files))
  songData = vector("list", length(h5Files))
  for( s in 1:length(h5Files)){
    if(!file.exists(paste(getwd(),"/output/persDiag/",gsub('.{3}$', '', h5Files[s]),".csv",sep=""))){
      song = h5file(paste( dirName ,h5Files[s],sep=""), mode = "a")
      songData[[s]] = song
      chromaFeatures[[s]] = song["/analysis/segments_pitches"] 
      hz = chromaFeatures[[s]][,1]
      tone = chromaFeatures[[s]][1,]
      ntone = ncol(chromaFeatures[[s]])
      nhz = nrow(chromaFeatures[[s]])
      sr = 2 #unsure about how to find sr found out it 2205 or something
      cfftlen = 2048
      t = 1:nhz*cfftlen/4/sr 
      chromMatrix = matrix(hz, t)
      pers <- ripsDiag(X = chromMatrix, maxdimension = 1, maxscale = 1, library="GUDHI", location = TRUE, printProgress=FALSE)$diagram
      # Write to file. Create directory if does not 
      # exist then write file
      pathTo <- dirname(h5Files[s])
      #print(paste(getwd(),"/output/persDiag/",gsub('.{3}$', '', h5Files[s]),".csv",sep=""))
      #paste(getwd(),"/output/persDiag/",gsub('.{3}$', '', h5Files[s]),".csv",sep="")
      if(!file.exists(pathTo)){
        dir.create(paste(getwd(),"/output/persDiag/",pathTo,sep=""), showWarnings =FALSE, recursive = TRUE)
      }
      #print( paste(getwd(),"/output/persDiag/",gsub('.{3}$', '', h5Files[s]),".csv",sep="") )
      write.csv(pers, row.names = FALSE, file = paste(getwd(),"/output/persDiag/",gsub('.{3}$', '', h5Files[s]),".csv",sep=""))
      
      persChromaSongs[[s]] = pers
    }
  }
}# end writeH5Persistence()



readPersistenceCSV <- function( filePath, dirPath, songNames = FALSE){
  " given the filePath in 'filePath' to the assumed csv file
     the file is read in and returned as a matrix. 
     Note: To plot the persistence diagram one must use  
     plot.diagram() from the TDA package. Personal note, 
     use 'band'= some_range_of_lifespan to explain plot
     pink band along diaginal.

    e.g. sonHousePreachin <- readPersistenceCSV(filePath=paste(getwd(),
          '/output/persDiag/sonHouse/TRYUOOZ128F427E839.csv',sep=''))
        plot.diagram(sonHousePreachin)
   "
  if(!missing(filePath)){
    return(as.matrix(read.csv(filePath)))
  }
  if(!missing(dirPath)){
    csvFiles <- list.files( dirPath, pattern = ".csv",recursive=TRUE, all.files = TRUE )
    persistenceDiagrams <- vector("list", length(csvFiles))

    for( f in 1:length(csvFiles)){
      persistenceDiagrams[[f]] <- as.matrix(read.csv( paste(dirPath,csvFiles[f],sep="") ))
      csvFiles[f] = gsub('.{4}$', '', basename(csvFiles[f])) #csvFiles[f])
    }
    if(!songNames){
      return(persistenceDiagrams)
    }else{
      return(list(diagrams = persistenceDiagrams, songNames = csvFiles))
    }
  }
}

topomdscale <- function(bottleneck, pwasserstein, p, fromFile){
  " Performs Topological Multi-Dimensional Scaling between 
     the distance structure of pairwise bottleneck distances
     as provided by 'bottleneck' or the distance structure 'pwasserstein'
     using p-wasserstein distances with p given by argument 'p'.
     If 'fromFile' is provided csv file specified by fromFile is 
     read as the distance structure to be used.
    Takes the pairwise distance matrix as and attempts to embed 
    it isometrically in Euclidean space with minimum distortion 
    of the metric. 
  
  The representation is only determined up to location (cmdscale
  takes the column means of the configuration to be at the origin),
  rotations and reflections. The configuration returned is given in
  principal-component axes,
  
  pairwiseBottleneckDist <- as.data.frame(read.csv( pathToFile ))
  mdsBdist <- topomdscale(bottleneck = pairwiseBottleneckDist)
  plot(mdsBdistpoints[,1], mdsBdist$points[,2])
  "
  if (!require(package = "stats")) {
    install.packages("stats")
  }
  library("stats")
  
  if(!missing(fromFile)){
    distMat <- read.csv(fromFile)
    mds <- cmdscale(distMat, eig=TRUE, k=2)
    x <- mds$points[, 1]
    y <- mds$points[, 2]
    return('x'=x,'y'=y)
  }
  if(!missing(bottleneck)){
    pairwiseBdist <- as.matrix(bottleneck)
    mdsBdist <- cmdscale(pairwiseBdist, eig=TRUE, k=2)
    x <- mdsBdist$points[,1]
    y <- mdsBdist$points[,2]
    return(list('x' = x, 'y'= y))
  }
}

twelveBarBluesComposition <- function(pathToFile, X, mds = FALSE, cluster, K, plot=FALSE){
  "The function uses the persistence diagram of a 12-bar blues
   musical composition which is unique to a small subset of songs 
   to compare to other songs in order to identify which songs, as 
   found in comparison between homological features of the 12-bar 
   structure, also contain a 12-bar blues composition.
  "
  "
    As the base song to find other 12-Bar blues we use 
    Death Letter Blues by Son House 'TRYZZOJ128F1466CB7.h5'

    e.g. twelveBarBluesComposition(pathToFile = paste(getwd(),'/output/pairWiseBottleneckAASonhouse.csv',sep=''
                                                 ),mds = TRUE, plot = TRUE, cluster = kmeans)
  
     or
  
        twelveBarBluesComposition(pathToFile = paste(getwd(),'/output/pairWiseBottleneckAASonhouse.csv',sep=''
                                                ), mds = FALSE, plot = TRUE, cluster = 'kmeans', K=15)"
  
  if (!require(package = "TDA")) {
    install.packages("TDA")
  }
  library("TDA")
  
  if (!require(package = "ggplot2")) {
    install.packages("ggplot2")
  }
  library("ggplot2")

  pairwiseBottleneckDist <- as.data.frame(read.csv( pathToFile ))
  rownames(pairwiseBottleneckDist) = pairwiseBottleneckDist[,1]
  pairwiseBottleneckDist <- pairwiseBottleneckDist[,2:length(pairwiseBottleneckDist)]

  if(!missing(cluster)){
    
    if(cluster == 'kmeans'){
      if(missing(K)){ K = 10}
      if(mds){
        for(q in 1:nrow(pairwiseBottleneckDist)){
          for(p in 1:ncol(pairwiseBottleneckDist)){
            if( pairwiseBottleneckDist[q,p] > 0.4){
              pairwiseBottleneckDist[q,p] = 0
            }
          }
        }
        mdsBdist <- topomdscale(bottleneck = pairwiseBottleneckDist)
        pairwiseBdistClusters <- kmeans(pairwiseBottleneckDist, K)

        pairwiseBdistClusters$cluster <- as.factor(pairwiseBdistClusters$cluster)
        print(pairwiseBdistClusters)
        if(plot){
        plot(mdsBdist, col = pairwiseBdistClusters$cluster, xlab = "multiDimensional Scaling with Bottleneck Distance", ylab = "")
        }
      }else{
      pairwiseBdistClusters <- kmeans(pairwiseBottleneckDist, K)
      print(pairwiseBdistClusters)
      pairwiseBdistClusters$cluster <- as.factor(pairwiseBdistClusters$cluster)
      clusterPoints <- pairwiseBdistClusters$centers

      "
      Now of the K clusters, the five in each cluster with minimal bottleneck distance
      are found.
      "
      centersCopy <- pairwiseBdistClusters$centers
      fiveNearest <- pairwiseBdistClusters$centers[,1:6]
      fiveMinBdist <- pairwiseBdistClusters$centers[,1:6]
      colnames(fiveNearest) <- c(1:6)
      colnames(fiveMinBdist) <- c(1:6)
      print(centersCopy[2,134])
      for(i in 1:nrow(pairwiseBdistClusters$centers)){
        for( j in 1:6){
          minBCenterIndex <- which.min(centersCopy[i,])
          fiveNearest[i,j] <- colnames(pairwiseBdistClusters$centers)[minBCenterIndex]
          fiveMinBdist[i,j] <- pairwiseBdistClusters$centers[i,minBCenterIndex]
          centersCopy[i,minBCenterIndex] = 2
        }
      }
      print(fiveNearest)
      print(fiveMinBdist)
      if(plot){
        for(q in 1:nrow(pairwiseBdistClusters$centers)){
          for(p in 1:ncol(pairwiseBdistClusters$centers)){
            if( pairwiseBdistClusters$centers[q,p] > 0.4){
              pairwiseBdistClusters$centers[q,p] = 0
            }
          }
        }
        png('kmeansCenters.png')
        plot(pairwiseBdistClusters$centers, col = c(1:K), ylab = "",xlab = "k-Means Clustering Centers")
        dev.off()
      }
      "
      mydata <- dat
      wss <- (nrow(mydata)-1)*sum(apply(mydata,2,var))
      for (i in 2:15) wss[i] <- sum(kmeans(mydata,
                                       centers=i)$withinss)
      plot(1:15, wss, type='b', xlab='Number of Clusters',
     ylab='Within groups sum of squares',
     main='Assessing the Optimal Number of Clusters with the Elbow Method',
     pch=20, cex=2)"
      #op <- par(oma=c(5,7,1,1))
      #par(op)
      #ggplot(pairwiseBottleneckDist, aes( pairwiseBottleneckDist , pairwiseBottleneckDist, color = pairwiseBdistClusters$cluster)) + geom_point()

      }
    }
    if(cluster == 'knearest'){
      sonHouse <- c("TRTMHJA12903CD8F4D","TRTMBEM128F427E66C","TRTPAZU128F92FFC5A","TRTPKNW128F427E83C",
                    "TRTOAOZ128F42760C0","TRTYTVH12903CE7403","TRYCIVS12903C95AB5","TRYZZOJ128F1466CB7",
                    "TRYIOMM128F9331B2D","TRYAVAL128F930CA78","TRYUOOZ128F427E839")
      deathLetterBlues <-"TRYZZOJ128F1466CB7"
      if(!missing(K)){K=20}
      minBdist <- 2.0
      kCount <- 0
      kNearest <-  vector("list", length(colnames(pairwiseBottleneckDist)))
      colCopy <- pairwiseBottleneckDist$'TRYCIVS12903C95AB5'
      for(i in 1:(K+1)){
        minBIndex <- which.min(colCopy)
        kNearest[i] <- colnames(pairwiseBottleneckDist)[minBIndex]
        colCopy[minBIndex] <- 2.0

      }
      return(kNearest[1:K+1])
    }
    ""
    
  }
}


kmeansPersistence <- function( persChromaSong){
  
  
  # SimpleKMeans au
  library(RWeka)
  
  nClusters <- 3
  
  d= as.data.frame(persChromaSongs)
  
  # SimpleKmeans with random initial cluster assignment
  
  clusters <- SimpleKMeans(bottleneckSongs[-1], Weka_control(N=nClusters, init = 0, V=TRUE))
  clusters
  cbind(bottleneckSongs[-1],predict(clusters))
  
  # SimpleKMeans with Kmeans++ initial cluster assignment
  
  clusters <- SimpleKMeans(bottleneckSongs, Weka_control(N=nClusters, init=1,  V=TRUE))
  clusters
  predict(clusters)
  
  # SimpleKmeans with Kmeans++ initial cluster assignment and "weka.core.Manhattandistance"
  
  clusters <- SimpleKMeans(bottleneckSongs, Weka_control(N=nClusters, init=1, A="weka.core.ManhattanDistance", V=TRUE))
  clusters
  predict(clusters)
  
  
  bottleneckMatrix= matrix(data=NA, nrow=length(persChromaSongs), ncol=length(persChromaSongs))
  for(q in 1:length(persChromaSongs)){
    for(p in 1:length(persChromaSongs)){
      bottleneckMatrix[q,p] = bottleneck(persChromaSongs[[q]], persChromaSongs[[p]], dimension = 1)
    }}
  
  
  clusters <- SimpleKMeans(bottleneckMatrix, Weka_control(N=nClusters, init=1, A="weka.core.ManhattanDistance", V=TRUE))
  clusters
  predict(clusters)
  
  
  ###K medoid clustering
  clus <- cluster::pam(bottleneckMatrix, 3)
  plot(bottleneckMatrix, xlab='', ylab='', axes=FALSE, xpd=NA,
       cex=4, pch=as.character( clus$cluster ))
  clus$clustering
  box()
  
}



cleanPersCSV <- function(dirPath){
  problemExample <- "/Users/multivax/Documents/PhD/2.spring.17/Classifications-Of-Songs-via-Homology-of-Chroma-Features/project/code/output/persDiag/"
  csvFiles <- list.files( paste(problemExample,dirPath, sep=""), pattern = ".csv",recursive=TRUE, all.files = TRUE )
  persistenceDiagrams <- vector("list", length(csvFiles))
  
  for( f in 1:length(csvFiles)){
    problemCSV <- readPersistenceCSV(paste(problemExample, csvFiles[f],sep=""))
    if(ncol(problemCSV) == 4){
      #reformat
      fileName <- basename(csvFiles[f])
      fixedCSV <- problemCSV[,-1]
      write.csv(fixedCSV,file=paste(problemExample,fileName,sep=""),row.names=FALSE, append=FALSE)
    }
  }

}
