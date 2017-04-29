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

def writeAttributestoCSV():
    for songDir in allSongDir:
        for h5File in [h5 for h5 in os.listdir( songDir ) if h5.endswith(".h5")]:
            filePath = songDir +"/"+ h5File
            # Could now print csv file that in order of
            # appearance seen in folder each h5 files are
            # printed per row in the following attribute order
            # songH5FileName, duration, sampleRate, artist, title, trackID, md5AudoHash, similarArtists
            #duration, sampleRate, artist, title, trackID, chromaFeatures, timbre, similarArtists, md5 = geth5Attributes(filePath)
            writeCSV(geth5Attributes(filePath))

def findArtist(artist, title):
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

#findArtist("Son House", "Death Letter Blues")
allH5 = []
for songDir in allSongDir:
    for h5File in [h5 for h5 in os.listdir( songDir ) if h5.endswith(".h5")]:
        allH5.append(h5File)
DeathLetterBlues = [sonHouse for sonHouse in allH5 if sonHouse.startswith("TRGCDIW128EF33E1E9") or sonHouse.startsWith(" TRHXYYB128F42795C6") or sonHouse.startsWith("TRRWDSF128F426CCFD") or sonHouse.startsWith("TRRRFYZ128F92EAB74")]

print(DeathLetterBlues)
