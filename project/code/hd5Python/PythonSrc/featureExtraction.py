import os
import sys

cwd = os.getcwd()
sys.path.append(cwd + '/hd5Python/PythonSrc/')
import hdf5_getters as h5get



filePath = "/Users/multivax/Documents/PhD/2.spring.17/Classifications-Of-Songs-via-Homology-of-Chroma-Features/project/code/MillionSongSubset/data/B/B/B/TRBBBEA128F93391BA.h5"
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

print("dir list", os.listdir("/Users/multivax/Documents/PhD/2.spring.17/Classifications-Of-Songs-via-Homology-of-Chroma-Features/project/code/MillionSongSubset/data/B/B/B/"))
print("duration",duration)
print("sr",sampleRate)
print("name",artist)
print("title",title)
print("chromaFeatures", chromaFeatures)
print("similar artists", similarArtists)
print("md5Audio", md5)
