import os
import sys

cwd = os.getcwd()
sys.path.append(cwd + '/hd5Python/PythonSrc/')
import hdf5_getters as h5get
import csv

def geth5Attributes(filePath, wantSimilar = False):
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
    if wantSimilar == True:
        return trackID, [similarArtists]
    else:
        return trackID, duration, sampleRate, artist, title, md5


        
def writeCSV(songAttributes, fileName):
    with open(os.getcwd()+'/'+fileName+'.csv','a') as csvfile:
        attributeFile = csv.writer(csvfile, delimiter=' ')
        attributeFile.writerow(songAttributes)

lenFullPath = len(os.getcwd()+'/MillionSongSubset/data/B/A/K')
allSongDir = [directory[0] for directory in os.walk(os.getcwd()+'/MillionSongSubset/data/') if len(directory[0]) == lenFullPath ]

def writeAttributeCSV(fileName, wantSimilar = False):
    for songDir in allSongDir:
        for h5File in [h5 for h5 in os.listdir( songDir ) if h5.endswith(".h5")]:
            filePath = songDir +"/"+ h5File
            # Could now print csv file that in order of
            # appearance seen in folder each h5 files are
            # printed per row in the following attribute order
            # songH5FileName, duration, sampleRate, artist, title, trackID, md5AudoHash, similarArtists
            #duration, sampleRate, artist, title, trackID, chromaFeatures, timbre, similarArtists, md5 = geth5Attributes(filePath)
            writeCSV(geth5Attributes(filePath, wantSimilar), fileName)

def printAllSongAttributes():
    for songDir in allSongDir:
        for h5File in [h5 for h5 in os.listdir( songDir ) if h5.endswith(".h5")]:
            filePath = songDir +"/"+ h5File
            print(geth5Attributes(filePath, wantSimilar = True))

def findArtistMillionSongs(artist = None, title = None):
    h5 = h5get.open_h5_file_read('msd_summary_file.h5')
    artistSongs = []
    artistSong = None
    for k in range(1000000):
        a_name = h5get.get_artist_name(h5,k)
        a_title = h5get.get_title(h5,k)
        if a_name == artist:
            artistSongs.append(h5get.get_track_id(h5,k))
        if a_title == title:
            print(h5get.get_track_id(h5,k))
            artistSong = h5get.get_track_id(h5,k)
    h5.close()
    if artist != None and title != None:
        return artistSongs, artistSong
    if artist != None and title == None:
        return artistSongs
    if artist == None and title != None:
        return artistSong


def findInSubset(songFile = None, artist = None, title = None):
    allH5 = []
    for songDir in allSongDir:
        for h5File in [h5 for h5 in os.listdir( songDir ) if h5.endswith(".h5")]:
            filePath = songDir +"/"+ h5File
            allH5.append(filePath)
    if songFile is not None:
        matches  = [songh5 for songh5 in allH5 if songh5 == songFile]
    if artist is not None or title is not None:
        artistSongs = []
        artistSong = None
        for h5file in allH5:
            h5 = h5get.open_h5_file_read(h5file)
            a_name = h5get.get_artist_name(h5)
            a_title = h5get.get_title(h5)
            if a_name == artist:
                artistSongs.append(h5file)
            if a_title == title:
                print(h5get.get_track_id(h5,k))
                artistSong = h5file
            h5.close()
        if artist is not None and title is not None:
            return artistSongs, artistSong
        if artist is not None and title is None:
            return artistSongs
        if artist is None and title is not None:
            return artistSong
   

#sonHouse == 'TRGCDIW128EF33E1E9.h5' or sonHouse == 'TRHXYYB128F42795C6.h5' or sonHouse == 'TRRWDSF128F426CCFD.h5' or sonHouse == 'TRRRFYZ128F92EAB74.h5']

def readAttributeCSV(fileName):
    with open(fileName, 'r') as csvFile:
        reader = csv.reader(csvFile)
        print(len(next(reader)))


#findArtist("Son House", "Death Letter Blues",)
#print(DeathLetterBlues)
#writeAttributeCSV('similarArtists', wantSimilar = True)
#readAttributeCSV('similarArtists.csv')
print( findInSubset(artist = "son house") )
