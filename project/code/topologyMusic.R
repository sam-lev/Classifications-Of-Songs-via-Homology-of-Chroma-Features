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
 


 

   
 
 ################          
 ###           ##
 "    Zahra    "
 ###           ##
 ################
 
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
 
 
 
} #End musicTopology()


pairwiseBottleneck <- function(fromDirPath, matrixData, persDiags, write){
  " Function computes the pairwise bottleneck distance between 
  persistence diagrams. If 'fromDirPath is present this function uses 
  readPersistenceCSV(fromDirPath) function's output, i.e. reads from csv 
  file 'pairwisePersistence.csv' to construct the matrix of pairwise 
  persistence diagrams between songs 
  If write is present it must be a string consisting of the desired 
  file name to write the matrix of pairwise bottleneck distances to.
  "
  # computes pairwise bottleneck distance between
  # all diagram elements stored to file whose
  # pathnames are in the list 'fromDirPath'
  if(!missing(fromDirPath)){
    # Fill list of persistence diagrams from csv files
    # under path 'fromDirPath'
    persistenceDiagrams = readPersistenceCSV(dirPath = fromDirPath)
    
    # Compute all pairwise persistence diagrams and store
    # to matrix.
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
      write.csv(bottleneckMatrix, row.names = FALSE, file = paste(getwd(),"/output/",write,".csv",sep=""))
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
  " Could be generalised later. Currently computes persistence 
  diagrams of from all h5songs chroma features in MillionSongSubset 
  and writes persistence diagrams to csv files under same folder
  architecture in /output/perDiag/ as seen for h5 files. i.e.
  if foobar.h5 is in MillionSongSubset/data/A/G/H/foobar.h5 then
  function writes the persistence diagram of foobar.h5's chroma 
  features to /output/persDiag/A/G/H/foobar.csv
  "
  
  h5Files = list.files("./MillionSongSubset/data/" , pattern = ".h5",recursive=TRUE, all.files = TRUE)
  chromaFeatures = vector("list", length(h5Files)) #rep(NA, length(h5Files))
  persChromaSongs = vector("list", length(h5Files)) #rep(NA, length(h5Files))
  songData = vector("list", length(h5Files))
  for( s in 1:length(h5Files)){
    if(!file.exists(paste(getwd(),"/output/persDiag/",gsub('.{3}$', '', h5Files[s]),".csv",sep=""))){
      song = h5file(paste("./MillionSongSubset/data/",h5Files[s],sep=""), mode = "a")
      songData[[s]] = song
      chromaFeatures[[s]] = song["/analysis/segments_pitches"] 
      hz = chromaFeatures[[s]][,1]
      tone = chromaFeatures[[s]][1,]
      ntone = ncol(chromaFeatures[[s]])
      nhz = nrow(chromaFeatures[[s]])
      sr = 2 #unsure about how to find sr 
      t = 1:nhz*cfftlen/4/sr 
      chromMatrix = matrix(hz, t)
      pers <- ripsDiag(X = chromMatrix, maxdimension = 1, maxscale = 1, library="GUDHI", location = TRUE, printProgress=FALSE)$diagram
      # Write to file. Create directory if does not 
      # exist then write file
      pathTo <- dirname(h5Files[s])
      print(paste(getwd(),"/output/persDiag/",gsub('.{3}$', '', h5Files[s]),".csv",sep=""))
      #paste(getwd(),"/output/persDiag/",gsub('.{3}$', '', h5Files[s]),".csv",sep="")
      if(!file.exists(pathTo)){
        dir.create(paste(getwd(),"/output/persDiag/",pathTo,sep=""), showWarnings =FALSE, recursive = TRUE)
      }
      print( paste(getwd(),"/output/persDiag/",gsub('.{3}$', '', h5Files[s]),".csv",sep="") )
      write.csv(pers, row.names = FALSE, file = paste(getwd(),"/output/persDiag/",gsub('.{3}$', '', h5Files[s]),".csv",sep=""))
      
      persChromaSongs[[s]] = pers
    }
  }
}# end writeH5Persistence()



readPersistenceCSV <- function( filePath, dirPath){
  " given the filePath in 'filePath' to the assumed csv file
     the file is read in and returned as a matrix. 
     Note: To plot the persistence diagram one must use  
     plot.diagram() from the TDA package. Personal note, 
     use 'band'= some_range_of_lifespan to explain plot
     pink band along diaginal.
   "
  if(!missing(filePath)){
    return(as.matrix(read.csv(filePath)))
  }
  if(!missing(dirPath)){
    csvFiles <- list.files(dirPath)
    persistenceDiagrams <- vector("list", length(csvFiles))
    for( f in 1:length(csvFiles)){
      persistenceDiagrams[[f]] <- as.matrix(read.csv(paste(dirPath,csvFiles[f],sep="")))
    }
    return(persistenceDiagrams)
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
    of the metric. "
  if (!require(package = "stats")) {
    install.packages("stats")
  }
  library("stats")
  
  if(!missing(fromFile)){
    distMat <- read.csv(fromFile)
    mds <- cmdscale(distMat, eig=TRUE, k=2)
    x <- mds$points[, 1]
    y <- mds$points[, 2]
    
    plot(x,y)
  }
}
