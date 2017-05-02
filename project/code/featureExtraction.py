import os
import sys

cwd = os.getcwd()
sys.path.append(cwd + '/hd5Python/PythonSrc/')
import hdf5_getters as h5get
import csv

def geth5Attributes(filePath = None, wantSimilar = False):
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
        return trackID, " ".join(map(str, similarArtists))
    else:
        return trackID, duration, sampleRate, artist, title, md5


        
def writeCSV(songAttributes, fileName):
    with open(os.getcwd()+'/'+fileName+'.csv','a') as csvfile:
        attributeFile = csv.writer(csvfile, delimiter=' ')
        attributeFile.writerow(songAttributes)

def writeAttributeCSV(fileName, wantSimilar = False):
    lenFullPath = len(os.getcwd()+'/MillionSongSubset/data/B/A/K')
    allSongDir = [directory[0] for directory in os.walk(os.getcwd()+'/MillionSongSubset/data/') if len(directory[0]) == lenFullPath ]
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
    lenFullPath = len(os.getcwd()+'/MillionSongSubset/data/B/A/K')
    allSongDir = [directory[0] for directory in os.walk(os.getcwd()+'/MillionSongSubset/data/') if len(directory[0]) == lenFullPath ]
    for songDir in allSongDir:
        for h5File in [h5 for h5 in os.listdir( songDir ) if h5.endswith(".h5")]:
            filePath = songDir +"/"+ h5File
            print(geth5Attributes(filePath, wantSimilar = True))

def findArtistMillionSongs(artist = None, title = None):
    # for this you must use msd_summary_file.h5 which
    # can be found on the million song database
    h5 = h5get.open_h5_file_read('/Users/multivax/Documents/PhD/2.spring.17/msd_summary_file.h5')
    artistSongs = []
    artistSong = None
    for k in range(1000000):
        a_name = h5get.get_artist_name(h5,k)
        a_title = h5get.get_title(h5,k)
        if artist.lower() in a_name.lower():
            artistSongs.append(h5get.get_track_id(h5,k))
        if title.lower() in a_title.lower():
            artistSong = h5get.get_track_id(h5,k)
    h5.close()
    if artist is not None and title is not None:
        return artistSongs, artistSong
    if artist is not None and title is None:
        return artistSongs
    if artist is None and title is not None:
        return artistSong


def findInSubset(songFiles = None, artist = None, title = None):
    # you may need to change the allSongDir path to one that
    # is appropriate to you
    allSongDir = [directory[0] for directory in os.walk(os.getcwd()+'/MillionSongSubset/data/')]
    allH5 = []
    for songDir in allSongDir:
        for h5File in [h5 for h5 in os.listdir( songDir ) if h5.endswith(".h5")]:
            filePath = songDir +"/"+ h5File
            allH5.append(filePath)
    if songFiles is not None:
        songAttr = []
        for song in songFiles:
            matches  = [songh5 for songh5 in allH5 if song in songh5]
            songAttr.append(geth5Attributes(filePath = matches[0]))
        return songAttr
    if artist is not None or title is not None:
        matchSongFiles = []
        matchSongAttr = []
        artistSong = None
        for h5file in allH5:
            h5 = h5get.open_h5_file_read(h5file)
            a_name = h5get.get_artist_name(h5)
            a_title = h5get.get_title(h5)
            print(a_title)
            if artist is not None and artist.lower() in a_name.lower():
                matchSongFiles.append(h5file)
                matchSongAttr.append([a_name, a_title])
            if title is not None and a_title.lower() in  title.lower():
                print(h5get.get_track_id(h5,k))
                artistSong = h5file
            h5.close()
        if artist is not None and title is not None:
            return matchSongFiles, matchSongAttr, artistSong
        if artist is not None and title is None:
            return matchSongFiles, matchSongAttr
        if artist is None and title is not None:
            return artistSong
   

#sonHouse == 'TRGCDIW128EF33E1E9.h5' or sonHouse == 'TRHXYYB128F42795C6.h5' or sonHouse == 'TRRWDSF128F426CCFD.h5' or sonHouse == 'TRRRFYZ128F92EAB74.h5']

def readAttributeCSV(fileName):
    with open(fileName, 'r') as csvFile:
        reader = csv.reader(csvFile)
        print(len(next(reader)))


#print(findArtistMillionSongs(artist = "Son House", title = "Death Letter Blues"))
#print(DeathLetterBlues)
#writeAttributeCSV('similarArtists', wantSimilar = True)
#readAttributeCSV('similarArtists.csv')
#print( findInSubset(artist = "house", title = 'death letter blues') )
#print(geth5Attributes("./MillionSongSubset/data/A/A/E/TRAAEEH128E0795DFE.h5"))
def sonHouseMatch():
    sonHouseSongs = ['TRMGMNY12903CFCB14', 'TRMRUMM128F425158D', 'TRMAPTX128F429AF11', 'TRMUWJC128F1466CB8', 'TRWMBUJ128F426F03B', 'TRWBCXF128F4251576', 'TRWIWPV128F426FE06', 'TRWDQVZ128F4246845', 'TRWDTZU128F426BB05', 'TRWDERQ128F1466CB4', 'TRWOBUG128F4251578', 'TRWXVPW12903CE7400', 'TRGFBOA128F427E83A', 'TRGNUGG128F42760C9', 'TRGTORI12903CE7402', 'TRGLCBT128F9331B19', 'TRGJOWM12903D0034C', 'TRHXYYB128F42795C6', 'TRCCCFC128F4251560', 'TRCCRQR128F427E66E', 'TRCJCNQ128F427C360', 'TRRMKGH128F4246864', 'TRRWDSF128F426CCFD', 'TRRABAP128F92FFC4B', 'TRRORGN128F4251573', 'TRRONOR128F42597F2', 'TRBDTEM128F4251570', 'TRFNCIB128F1466CB1', 'TRFPYHU128F4246855', 'TRFYPPX128F4251587', 'TRQGDSI128F1466CBA', 'TRQFAUP128F4251580', 'TRQZTHI128F427E668', 'TRQSOPR128F425156D', 'TRZHPKM128F426C2BF', 'TRZCSDN128F425158B', 'TRZJGIK128F427E838', 'TRIUVLB128F427E83B', 'TRIOIEI128F92FA589', 'TRIOKJB128F425159C', 'TRAPRRF128F4251581', 'TRATNMT128F1466CB2', 'TRAOPNA128F427E83F', 'TRNBHKS128F427E66B', 'TRNBPIV128F1466CB6', 'TRNQGFK12903C95AB9', 'TRNNDGO128F425157B', 'TRNTDSB12903C95AA5', 'TRNLSZY128F4251577', 'TRNVCIO128F427E841', 'TRPMNOC128F426D66A', 'TRPBINQ12903CE73F9', 'TRTMHJA12903CD8F4D', 'TRTMBEM128F427E66C', 'TRTPAZU128F92FFC5A', 'TRTPKNW128F427E83C', 'TRTOAOZ128F42760C0', 'TRTYTVH12903CE7403', 'TRUMUKO128F4251564', 'TRUCOLD128F427E83D', 'TRUAYIE128F426CEAE', 'TRUXFNE12903CE73F3', 'TRLWXCL128F426FE01', 'TRLRLEN128F4251595', 'TRLOUIA128F4251585', 'TREWYAF128F92FFC53', 'TREQKDX128F427E83E', 'TREZLXL12903CEBF5E', 'TREPRWR128F1466CB5', 'TRETOIV128F426C12E', 'TREVJXS12903C95AA4', 'TREYFZS128F427E665', 'TREYLGE128F425157E', 'TREYSGX128F426E6AE', 'TRJSJGW12903C95AA8', 'TRSBSPH128F9331AE6', 'TRSBSFM128F427E837', 'TRSIGWP128F4251599', 'TRSJXGK128F427E840', 'TRSXRFS12903CE6EE4', 'TRSXXMI128F9331B3F', 'TRVCIMN128F427C433', 'TRDMWUL128F4270591', 'TRDWBVQ128F426CE1D', 'TRDBDNM128F426BB0B', 'TRDNTWI12903CE7401', 'TRDJADO128F426BB09', 'TRDDNWC128F9311E1C', 'TRDOZBS128F14A1B5A', 'TRDYWHG128F4251583', 'TROBDLX128F425157C', 'TRXIQQO128F425159A', 'TRXLKKL128F4251572', 'TRXEQKP12903C95AAE', 'TRXSWTF128F4251582', 'TRKMGOC128F426C12C', 'TRKIETH128F425156F', 'TRKEJGG128F4251597', 'TRKJKZT128F427E664', 'TRKVGOM128F4251589', 'TRKDNFN128F4251579', 'TRYCIVS12903C95AB5', 'TRYZZOJ128F1466CB7', 'TRYIOMM128F9331B2D', 'TRYAVAL128F930CA78', 'TRYUOOZ128F427E839']
    deathLetterBlues = 'TRYZZOJ128F1466CB7'
    for song in sonHouseSongs:
        if  'TRY' in song or 'TRT' in song:
            print song
#sonHouseMatch()
"""
TRTMHJA12903CD8F4D
TRTMBEM128F427E66C
TRTPAZU128F92FFC5A
TRTPKNW128F427E83C
TRTOAOZ128F42760C0
TRTYTVH12903CE7403
TRYCIVS12903C95AB5
TRYZZOJ128F1466CB7
TRYIOMM128F9331B2D
TRYAVAL128F930CA78
TRYUOOZ128F427E839
"""
nearPreachin= ["TRYUOOZ128F427E839","TRAAIRG128F93265E8","TRAAQPS128F429161D","TRAAQIH128F428BDEA","TRAAVIT128F92E657C","TRAAWGY128F4298A09","TRYZZOJ128F1466CB7","TRAAIAE128F42AC53D","TRAASYO128F4263883","TRAAGKF128F932D5A2","TRAAVAH128F4284D7C","TRAAAFD128F92F423A","TRTOAOZ128F42760C0","TRAABNV128F425CEE1","TRAATZQ128F425E5F1","TRAAERU128F930674F","TRAAQLJ128F428E870","TRAAJJG128F4284B27","TRAAARJ128F9320760","TRAAGHM128EF35CF8E","TRAAGVM128E0784D95"]

#print(findInSubset(songFiles = nearPreachin))# geth5Attributes(songList = nearPreachin))
nearMCB = ["TRYCIVS12903C95AB5","TRAAEEH128E0795DFE","TRAAXPA128F92FC706"
           ,"TRAASSO12903CDD2FF"
           ,"TRTMBEM128F427E66C"
           ,"TRAABDL12903CAABBA"
           ,"TRAAPLS12903CB2D40",
           "TRAACPE128F421C1B9",
           "TRAAUMZ12903CA97DC"]
nearMCB2 = ["TRAAKTK128F9343DCC",
           "TRAAVIT128F92E657C",
           "TRAAQYN128F92ED77E",
           "TRAAJJG128F4284B27",
           "TRAAIAX128F930531B",
           "TRAAZPJ128F4292DB5",
           "TRAAZQF128F4285BC3",
           "TRAAAAW128F429D538",
           "TRAAQLJ128F428E870",
           "TRAAERU128F930674F",
           "TRAAWVQ128F9313D31",
           "TRAAQPS128F429161D"]
#print(findInSubset(songFiles = nearPreachin))

#
# from kmeans clustering centers 5 minimum with 15 clusters
#
k1 = ["TRAAZPJ128F4292DB5" ,"TRAATMY128F149562F" ,"TRAAEEH128E0795DFE", "TRAAJJG128F4284B27" ,"TRYCIVS12903C95AB5", "TRAACPE128F421C1B9"]
k2 = ["TRAAKSB128F934D73A", "TRAAMGJ12903CD315C", "TRAAUSW128F426646E", "TRAAJVJ128F42570AD", "TRAAIYR12903CD290D", "TRAARHS128F931689D"]
k3=  ["TRAAGXN128F92DE52C", "TRAALHD12903D0BD6A" ,"TRAAIAN12903CFF16D", "TRAAIFG12903CB0BCB", "TRAACHN128F1489601", "TRAANLO128F422AE37"]
k4 =[ "TRAAHEH128F427FCEF", "TRAAVOU128F422F292" ,"TRAAHEQ128F93057B4", "TRAAPSQ12903CAE86A", "TRAAHEG128E07861C3", "TRAAPCC128F42519A7"]
k5 =[ "TRAAHYV128F932A425", "TRAARUN128F932B1A0" ,"TRAAQFC128F425ED83", "TRAAYTI128F933663B", "TRAAPPQ128F14961F5", "TRAACZK128F4243829"]
k6 =[ "TRAASFP128F42980BD", "TRAACZK128F4243829" ,"TRAAWUM128F426009A", "TRAAJMH128E0785407", "TRAAOTT128F9342E04", "TRAASZE128F93499AE"]
k7 =[ "TRAAJJG128F4284B27", "TRAAQPS128F429161D" ,"TRAACPE128F421C1B9", "TRAASSO12903CDD2FF", "TRAAVIT128F92E657C", "TRAAGMJ12903CAD7D4"]
k8  =["TRAAFON128F42926FB", "TRAAGKS128F9314884" ,"TRAAMKL128F9347D8A", "TRAAQZF12903CEDCB8", "TRAARSF12903CD1B12", "TRTMHJA12903CD8F4D"]
k9  =["TRAAJNM12903CF4750" ,"TRAAUGZ128F9307F78" ,"TRAAEDJ12903CCE390" ,"TRAAXUY128F426BA55", "TRAAIGD12903CA2F3C", "TRAAJJV128F42A2F99"]
k10 =["TRAAUSW128F426646E", "TRAASDM12903CD474E" ,"TRAABCL128F4286650", "TRAALCM128F4264370", "TRAAKDG128F42A0ECB", "TRAAUIH128F4254C9D"]
k11 =["TRAABVM128F92CA9DC", "TRAAXKI128F92DAC1B" ,"TRAAKAV128F931BE25" ,"TRAAEOU128F9309DE0", "TRAAPVV128F425D64C", "TRAANMF12903CEAFEF"]
k12 =["TRAATGM12903C9CD7C", "TRAATXO128F92D97A6" ,"TRAAENP128F147BF32" ,"TRAALUX12903CFD50B", "TRAADIV12903CB15C1", "TRYAVAL128F930CA78"]
k13 =["TRAAWUM128F426009A", "TRAASFP128F42980BD" ,"TRAAPIV128F1493132" ,"TRAACZK128F4243829", "TRAAVIH12903CFAC90" ,"TRAANMF12903CEAFEF"]
k14 =["TRAAAEF128F4273421", "TRAASDM12903CD474E", "TRAAKDG128F42A0ECB", "TRAAMXP128F4264F6A", "TRAAIBB12903CECB35" ,"TRAAUMN128F4273561"]
k15= ["TRAAERU128F930674F", "TRAAHDH128F93094B6", "TRAAQPS128F429161D" ,"TRAAIRG128F93265E8", "TRAAHTH12903CCFFE1", "TRAAVAH128F4284D7C"]

#print(findInSubset(songFiles = k15))
kmed1 =["TRAAGVM128E0784D95", "TRAAGVM128E0784D95" ,"TRAAZHO128F427C]EC", "TRAAHDH128F93094B6", "TRAAERU128F930674F", "TRAAHZP12903CA25F4" ,"TRAAKDG128F42A0ECB"]
kmed2 = ["TRAABYN12903CFD305", "TRAABYN12903CFD305", "TRAAKDF128F92FA302", "TRAAXUY128F426BA55", "TRAAGUA128F93236E8" ,"TRAAOTD128F42322FC", "TRAADLH12903CA70EE"]
kmed3 = ["TRAACCG128F92E8A55" ,"TRAACCG128F92E8A55", "TRAAZGC12903CB5DFB", "TRAAHZP12903CA25F4", "TRAAHTH12903CCFFE1", "TRAAMXP128F4264F6A", "TRAAPQQ128F427C750"]
kmed4 = ["TRAAIHL128F92E6DDA", "TRAAIHL128F92E6DDA", "TRAAMZR128F9315DCC", "TRAALDU128F930D56B", "TRAAUIH128F4254C9D", "TRAAPXK128F934DF25", "TRAADYI128E078FB38"]
kmed5 =[ "TRAAEUB12903CDD3EA", "TRAAEUB12903CDD3EA", "TRAAPSQ12903CAE86A", "TRAAVOU128F422F292" ,"TRAAPCC128F42519A7" ,"TRAAHEH128F427FCEF", "TRAAHEQ128F93057B4"]
kmed6 = ["TRAAENC128F1451DE9", "TRAAENC128F1451DE9", "TRAAPWD128F931ACD4", "TRAACZK128F4243829", "TRAAWUM128F426009A", "TRAAAVG12903CFA543" ,"TRAAHYV128F932A425"]
kmed7 = ["TRAAODT12903D0D16E", "TRAAODT12903D0D16E", "TRAAPPQ128F14961F5", "TRAAZBD128F427A97D" ,"TRTPKNW128F427E83C" ,"TRAARUN128F932B1A0", "TRAASZE128F93499AE"]
kmed8 = ["TRAAUMJ128F933BDB6", "TRAAUMJ128F933BDB6", "TRAAUEO128E0780CB9", "TRAAIII128F427D5D8", "TRAASYO128F4263883" ,"TRAAIAE128F42AC53D", "TRAABNV128F425CEE1"]
kmed9 = ["TRAADQX128F422B4CF", "TRAADQX128F422B4CF", "TRAABJV128F1460C49", "TRAADYB128F92D7E73", "TRAAIYR12903CD290D", "TRAAJVJ128F42570AD", "TRAAMGJ12903CD315C"]
kmed10 = ["TRAAWUH128F427F12C" ,"TRAAWUH128F427F12C", "TRAANMF12903CEAFEF", "TRAANTD128F9323B23", "TRAABLR128F423B7E3", "TRAAPIV128F1493132", "TRAACQT128F9331780"]
kmed11 = ["TRAABRB128F9306DD5" ,"TRAABRB128F9306DD5", "TRAAYPG128F933C941", "TRAALEO128F42A7AEA", "TRAAKVD12903CE8474", "TRAAORY128F934984F", "TRAADHS12903CE70A9"]
kmed12 = ["TRAAZKV12903CDDFE0", "TRAAZKV12903CDDFE0", "TRAAKVD12903CE8474", "TRAAKAV128F931BE25", "TRAAYPG128F933C941", "TRAAPVV128F425D64C" ,"TRAAORY128F934984F"]
kmed13 = ["TRAAGXN128F92DE52C", "TRAAGXN128F92DE52C", "TRAALHD12903D0BD6A", "TRAAIFG12903CB0BCB", "TRAAOIN128F9303790", "TRAAJEK128F422F55A" ,"TRAAFJW128F428A424"]
kmed14 = ["TRAADIV12903CB15C1" ,"TRAADIV12903CB15C1", "TRAASOD128F428154C", "TRAATGM12903C9CD7C", "TRYAVAL128F930CA78", "TRAAENP128F147BF32", "TRAATXO128F92D97A6"]
kmed15 = ["TRAAVJL128F14ACADD", "TRAAFON128F42926FB", "TRAAGKS128F9314884", "TRAAMKL128F9347D8A", "TRAAQZF12903CEDCB8", "TRAARSF12903CD1B12" ,"TRTMHJA12903CD8F4D"]


print(findInSubset(songFiles = kmed2))
#7, 15
