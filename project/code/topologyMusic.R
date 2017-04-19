musicTopology <- function() {
  
  setwd( getwd() )
  
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
  
  
 # Here we list all h5 files in a folder, calculate the persistant homology for each
 # song and add the corresponding diagram to 'persChromaSongs' and each songs chromatic
 # features to 'chromaFeatures'

 h5Files = list.files("./MillionSongSubset/data/A/R/R/" , pattern = ".h5", all.files = TRUE)
 chromaFeatures = vector("list", length(h5Files)) #rep(NA, length(h5Files))
 persChromaSongs = vector("list", length(h5Files)) #rep(NA, length(h5Files))
 songData = vector("list", length(h5Files))
 song1 = h5file(paste("./MillionSongSubset/data/A/R/R/",h5Files[1],sep=""), mode = "a")
 for( s in 1:length(h5Files)){
   song = h5file(paste("./MillionSongSubset/data/A/R/R/", h5Files[s],sep=""), mode = "a")
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
   persChromaSongs[[s]] = pers
 }
 
 bottleneckSongs = matrix(1:length(persChromaSongs))
 firstSongPers = persChromaSongs[[1]]
 # Compute the bottleneck distance between persistant diagrams of first song
 # in folder and all other songs. Eventually need to compute all bottleneck
 # distance between persistence diagrams, i.e. if n songs in folder 
 # then compute n! persistence diagrams to determine nearly topologically equivalent 
 # songs. 
 for(p in 1:length(persChromaSongs)){
   bottleneckSongs[p] = bottleneck(firstSongPers, persChromaSongs[[p]], dimension = 1)
 }
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
 
 pairwisePersistence < function(args , write , fromFolders){
   " Returns the pairwise persistence diagrams between 
    all elements in 'args' as a list if 'write' is missing.
    Otherwise, if 'fromFolders' is missing and 'write' is TRUE 
    then computes the pairwise persistence diagrams of all 
    elements in args and writes it to 'pairwisePersistence.csv' 
    by default or as 'fileName.csv' if file name specified as 
    string in 'fileName' "
   if(!missing(write) & !missing(args) | write == TRUE & !missing(args)){
     #write to csv file pairwise distance of all elements in 'args'
     
     if(!missing(fileName)){
       # write to specified filename
       
       #write.csv(dataFrame, row.names = FALSE, file = paste("output/",fileName,".csv",sep=""))
     }
     else{
       #write.csv(dataFrame, row.names = FALSE, file = "output/pairwisePersistence.csv")
     }
   }
   if(!missing(write) & !missing(fromFolders) | write == TRUE & !missing(fromFolders)){
     # write to csv file pairwise distance between all data files
     # recursively for all data files in folders specified 
     # by 'fromFolders'
   }
   else{
     pairwisePersistence <- c()
     
     return(pairwisePersistence)
   }
 }
 pairwiseBottleneck <- function(read, matrixData, diags){
 " Function computes the pairwise bottleneck distance between 
   persistence diagrams. If read = True this function uses pairwisePersistence() 
   function's output, i.e. reads from csv file 'pairwisePersistence.csv' 
   to construct the matrix of pairwise persistence diagrams between songs
 "
   # computes pairwise bottleneck distance between
   # all diagram elements stored to file 'read'
   if(!missing(read)){
     
   }
   # computes pairwise bottleneck distance between
   # all diagram elements in diag
   if(!missing(diags)){
     
   }
   # computes pairwise bottleneck distance between
   # all matrix data elements in natrixData
   if(!missing(matrixData)){
     
   }
 
 }
 
}
