import os
import sys

cwd = os.getcwd()
sys.path.append(cwd + '/hd5Python/PythonSrc/')
import hdf5_getters as h5get
import csv

def geth5Attributes(filPath):
    h5 = h5get.open_h5_file_read(filePath)
    duration = h5get.get_duration(h5)
    sampleRate = h5get.get_analysis_sample_rate(h5)
    artist = h5get.get_artist_name(h5)
    title = h5get.get_title(h5)
    trackID = h5get.get_track_id(h5)
    chromaFeatures = h5get.get_segments_pitches(h5)
    timbre = h5get.get_segments_timbre(h5)
    similarArtists = h5get.get_similar_artists(h5)
    md5 = h5get.get_audio_md5(h5)
    h5.close()
    # features not added to csv
    #chromaFeatures, timbre
    return trackID, duration, sampleRate, artist, title, md5


        
def writeCSV(songAttributes):
    with open(os.getcwd()+'/'+'attributes.csv','a') as csvfile:
        attributeFile = csv.writer(csvfile, delimiter=' ')
        attributeFile.writerow(songAttributes)

lenFullPath = len(os.getcwd()+'/MillionSongSubset/data/B/A/K')
allSongDir = [directory[0] for directory in os.walk(os.getcwd()+'/MillionSongSubset/data/') if len(directory[0]) == lenFullPath ]

for songDir in allSongDir:
    for h5File in os.listdir( songDir ):
        filePath = songDir +"/"+ h5File
        # Could now print csv file that in order of
        # appearance seen in folder each h5 files are
        # printed per row in the following attribute order
        # songH5FileName, duration, sampleRate, artist, title, trackID, md5AudoHash, similarArtists
        #duration, sampleRate, artist, title, trackID, chromaFeatures, timbre, similarArtists, md5 = geth5Attributes(filePath)

        "Current issue:
        Unable to open/create file '/Users/multivax/Documents/PhD/2.spring.17/Classifications-Of-Songs-via-Homology-of-Chroma-Features/project/code/MillionSongSubset/data/A/R/R/.DS_Store'
        "
        writeCSV(geth5Attributes(filePath))
      

"""
print("dir list", os.listdir("/Users/multivax/Documents/PhD/2.spring.17/Classifications-Of-Songs-via-Homology-of-Chroma-Features/project/code/MillionSongSubset/data/B/B/B/"))
print("duration",duration)
print("sr",sampleRate)
print("name",artist)
print("title",title)
print("chromaFeatures", chromaFeatures)
print("similar artists", similarArtists)
print("md5Audio", md5)
"""
