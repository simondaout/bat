
import os
import sys
import logging
import fnmatch

# add local directories to import path                      #hs+

sys.path.append ('../tools/')                     
sys.path.append ('../Common/')
                                                            #hs-
#

from   obspy.core.utcdatetime import UTCDateTime
from   obspy.core import  *
from   obspy.core import  read

from   obspy.taup.taup import getTravelTimes
from   obspy.taup.taup import locations2degrees

from   obspy.signal.trigger   import classicSTALTA,triggerOnset,recSTALTA,plotTrigger
from   obspy.signal.trigger   import pkBaer
import obspy.signal.cross_correlation

#

from   ConfigParser import SafeConfigParser
import numpy 

#       Import from common

import  Basic
import  Globals
import  Logfile
import  Debug
from    ConfigFile import ConfigObj 

#       Import from Tools

#       Import from Process


# --------------------------------------------------------------------------------------------------

logger = logging.getLogger('xcorr')
logger.setLevel(logging.DEBUG)

formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')

ch = logging.StreamHandler()
ch.setLevel(logging.DEBUG)
ch.setFormatter(formatter)

logger.addHandler(ch)

# --------------------------------------------------------------------------------------------------

class Event(object):
    def __init__(self, lat, lon, depth, time='', region='', strike=0, dip=0, rake=0):
        self.lat    = lat
        self.lon    = lon
        self.depth  = depth
        self.region = region
        self.time   = time
        self.strike = strike
        self.dip    = dip
        self.rake   = rake
# --------------------------------------------------------------------------------------------------
        
class Corr(object):
    def __init__(self, index, value, sampleindex):
        self.index       = index
        self.value       = value
        self.sampleindex = sampleindex
# --------------------------------------------------------------------------------------------------

class Station(object):
    def __init__(self, net, sta, loc, comp, lat=0, lon=0, ele=0, dip=0, azi=0, gain=0, takeoff=0, backazi=0, sazi=0):
        self.net  = net
        self.sta  = sta
        self.loc  = loc
        self.comp = comp
        self.lat  = lat
        self.lon  = lon
        self.ele  = ele
        self.dip  = dip
        self.azi  = azi
        self.gain = gain
        self.takeoff = takeoff
        self.backazi = backazi
        self.sazi = sazi

    def getName(self) :
        return self.net + '.' + self.sta + '.' + self.loc + '.' + self.comp
# --------------------------------------------------------------------------------------------------

def readMetaInfoFile (fname):

        MetaL = []
        Logfile.red ('Parsing MetaInfoFile')

        try:
                    evfile = fname
                    print 'METAFILE ',evfile
                    fobj = open(evfile, 'r')

                    for i in fobj:
                        i = i.split()
                        net = i[0]
                        sta = i[1]
                        loc = i[2]
                        comp = i[3]
                        lat = i[4]
                        lon = i[5]
                        ele = i[6]
                        dip = i[7]
                        azi = i[8]
                        gain = i[9]
                        #net, sta, loc, comp, lat, lon, ele, dip, azi, gain, inst = i.split()

                        if fnmatch.fnmatch(comp, '*HZ'):
                            MetaL.append(Station(net, sta, loc, comp, lat, lon, ele, dip, azi, gain))

                    Logfile.red ('%d ENTRIES IN METAFILE FOUND' % (len(MetaL)) )
        except:
            Logfile.exception ('METAFILE NOT READABLE')

        return MetaL
# --------------------------------------------------------------------------------------------------

def parseConfig (ofile):

        logger.info('\033[31m Parsing %s File \033[0m \n' % (ofile))
        cDict = {}
        parser = SafeConfigParser()
        parser.read(ofile)
        
        for section_name in parser.sections():
            for name, value in parser.items(section_name):
                cDict[name] = value

        return cDict
# --------------------------------------------------------------------------------------------------

def filterStations(StationList, Config, Origin, network):
    
    F = []
    minDist = int   (Config['mindist'])
    maxDist = int   (Config['maxdist'])
    o_lat   = float (Origin['lat'])
    o_lon   = float (Origin['lon'])
    
    Logfile.red ('Filter stations with configured parameters')
    
    for i in StationList:
        if i.loc == '--':
            i.loc = ''

        streamID = i.net + '.' + i.sta + '.' + i.loc + '.' + i.comp

        for j in network:
            if fnmatch.fnmatch(streamID, j):
                sdelta = locations2degrees(o_lat, o_lon, float(i.lat), float(i.lon))
                print streamID, sdelta,' degree'

                if sdelta > minDist and sdelta < maxDist:
                    #if i.net != 'GB':
                        F.append(Station(i.net, i.sta, i.loc, i.comp, i.lat, i.lon, i.ele, i.dip, i.azi, i.gain))


    Logfile.red  ('%d STATIONS LEFT IN LIST' % len(F))
    return F
# --------------------------------------------------------------------------------------------------

def calculateTimeWindows(mint, Origin):

    tw          = {}
    st          = str(Origin.time)[:-1]
    tw['start'] = UTCDateTime(UTCDateTime(st) + (mint - 100))
    tw['end']   = tw['start'] + 200
    
    tw['xcorrstart'] = UTCDateTime(UTCDateTime(st) + (mint - 8))
    tw['xcorrend']   = tw['xcorrstart'] + 20

    logger.info (' ORIGIN TIME %s  MINT %f' % (UTCDateTime(Origin.time),mint))
    logger.info (' OVERALL TIME WINDOW: %s - %s' % (tw['start'], tw['end']))
    logger.info (' XCROSS TIME WINDOW: %s - %s' % (tw['xcorrstart'], tw['xcorrend']))

    return tw
# --------------------------------------------------------------------------------------------------

def resampleWaveform(Waveform,end_frequence):

    new_frequence = end_frequence

    if Waveform.stats.sampling_rate == end_frequence:
        return Waveform
        
    if int(Waveform.stats.sampling_rate) == 20:
        oldsr = Waveform.stats.sampling_rate
        factor = int(Waveform.stats.sampling_rate / new_frequence + 0.5)
        logger.info('downsampling %s from %d Hz to %d Hz with factor %d' % (Waveform, oldsr, new_frequence, factor))
        Waveform.decimate(int(factor), strict_length=False, no_filter=False)
        return Waveform
    
    if Waveform.stats.sampling_rate == 25:
        #upsampling with factor = 4
        data = numpy.zeros((Waveform.stats.npts*2))
        test = numpy.fft.rfft(Waveform.data)
        zero = numpy.fft.rfft(data)
        t = numpy.concatenate((test,zero))
        tbdata = numpy.fft.irfft(t)
        
        tm = Trace(tbdata)
        tm.stats.network = Waveform.stats.network
        tm.stats.station = Waveform.stats.station
        tm.stats.sampling_rate = Waveform.stats.sampling_rate*4
        tm.stats.channel = Waveform.stats.channel
        tm.stats.starttime = utcdatetime.UTCDateTime(Waveform.stats.starttime)
        tm.stats._format = 'MSEED'
        
        oldsr = tm.stats.sampling_rate
        factor = int(tm.stats.sampling_rate / new_frequence + 0.5)
        logger.info('downsampling %s from %d Hz to %d Hz with factor after upsampling %d' % (tm, oldsr, new_frequence, factor))
        tm.decimate(int(factor), strict_length=False, no_filter=False)
        return tm
    
    if Waveform.stats.sampling_rate == 40:
        oldsr = Waveform.stats.sampling_rate
        factor = int(Waveform.stats.sampling_rate / new_frequence + 0.5)
        logger.info('downsampling %s from %d Hz to %d Hz with factor %d' % (Waveform, oldsr, new_frequence, factor))
        Waveform.decimate(int(factor), strict_length=False, no_filter=False)

        return Waveform
    
    if Waveform.stats.sampling_rate == 50:
        #upsampling with factor = 2
        data = numpy.zeros((Waveform.stats.npts))
        test = numpy.fft.rfft(Waveform.data)
        zero = numpy.fft.rfft(data)
        t = numpy.concatenate((test,zero))
        tbdata = numpy.fft.irfft(t)
        
        tm = Trace(tbdata)
        tm.stats.network = Waveform.stats.network
        tm.stats.station = Waveform.stats.station
        tm.stats.sampling_rate = Waveform.stats.sampling_rate*2
        tm.stats.channel = Waveform.stats.channel
        tm.stats.starttime = utcdatetime.UTCDateTime(Waveform.stats.starttime)
        tm.stats._format = 'MSEED'
        
        oldsr = tm.stats.sampling_rate
        factor = int(tm.stats.sampling_rate / new_frequence + 0.5)
        logger.info('downsampling %s from %d Hz to %d Hz with factor after upsampling %d' % (tm, oldsr, new_frequence, factor))
        tm.decimate(int(factor), strict_length=False, no_filter=False)
        return tm

    
    if Waveform.stats.sampling_rate == 80:
        oldsr = Waveform.stats.sampling_rate
        factor = int(Waveform.stats.sampling_rate / new_frequence + 0.5)
        logger.info('downsampling %s from %d Hz to %d Hz with factor %d' % (Waveform, oldsr, new_frequence, factor))
        Waveform.decimate(int(factor), strict_length=False, no_filter=False)
        
    if Waveform.stats.sampling_rate == 100:
        oldsr = Waveform.stats.sampling_rate
        factor = int(Waveform.stats.sampling_rate / new_frequence + 0.5)
        logger.info('downsampling %s from %d Hz to %d Hz with factor %d' % (Waveform, oldsr, new_frequence, factor))

        Waveform.decimate(int(factor), strict_length=False, no_filter=False)

        return Waveform
    
    else:
        oldsr  = Waveform.stats.sampling_rate
        factor = int(Waveform.stats.sampling_rate / new_frequence + 0.5)
        logger.info('downsampling %s from %d Hz to %d Hz with factor %d' % (Waveform, oldsr, new_frequence, factor))
        Waveform.decimate(int(factor), strict_length=False, no_filter=False)
        return Waveform

# --------------------------------------------------------------------------------------------------

def filterWaveform(Waveform):
    logger.info('\033[31m Filter Waveform: \033[0m' % ())
    
    new_frequence = 20
    
    st = obspy.core.stream.Stream()
    
    for i in Waveform:
        
        j = resampleWaveform(i,new_frequence)
        j.detrend("simple")
        j.filter("bandpass", freqmin=1, freqmax=4, corners=3, zerophase=False)
 
        st.append(j)
    
    return st
# --------------------------------------------------------------------------------------------------

def signoise(Waveform, ttime, Origin, path):
    
    st = str(Origin.time)[:-1]
    
    ponset = UTCDateTime(st) + ttime
    
    winnoise_start = Waveform.stats.starttime
    winnoise_end = ponset - 10
    winsig_start = ponset - 2
    winsig_end = ponset + 5
    
    winnoise = read(path, format="MSEED", starttime=winnoise_start, endtime=winnoise_end, nearest_sample=True)
    winsig   = read(path, format="MSEED", starttime=winsig_start, endtime=winsig_end, nearest_sample=True)
    
    psignal = abs(winsig.max()[0])
    pnoise = abs(winnoise.max()[0])

    signoise = float(psignal) / float(pnoise)
    print psignal, pnoise, signoise
    
    return signoise
# --------------------------------------------------------------------------------------------------

def readWaveformsCross(station, tw, Origin, ttime):
    
    time = Origin.time
    ts = time.split('T')
    
    datet = ts[0]
    datet = datet.split('-')
    year = datet[0].strip()
    month = datet[1]
    day = datet[2]
    #timep = ts[1][:-1]
     
    #print time,ts,year,month,day
    julday  = UTCDateTime(int(year), int(month), int(day)).julday
    julday  = "%03d" % julday
    sdspath = os.path.join(os.getcwd(), year)

    if station.loc == '--':
        station.loc = ''

    streamData = station.net + '.' + station.sta + '.' + station.loc + '.' + station.comp + '.D.' + str(year) + '.' + str(julday)
    
    entry = os.path.join(sdspath, station.net, station.sta, station.comp + '.D', streamData)
    st    = read(entry, format="MSEED", starttime=tw['start'], endtime=tw['end'], nearest_sample=True)

    if len(st.getGaps()) > 0:
        st.merge(method=0, fill_value='interpolate', interpolation_samples=0)
    
    snr = signoise(st[0], ttime, Origin, entry)
    
    stream = filterWaveform(st)
    
    print 'OVERALL: ', stream
    
    stream.write(streamData,format='MSEED')
    
    stream.trim(tw['xcorrstart'], tw['xcorrend'])
    print 'CROSSWINDOW: ', stream
    stream[0].stats.starttime = UTCDateTime(1000)
    #print 'AFTERFILTER: ',stream

    return stream, snr
# --------------------------------------------------------------------------------------------------

def readWaveformsPicker(station, tw, Origin, ttime):
    
    time = Origin.time
    ts   = time.split('T')
    
    datet = ts[0]
    datet = datet.split('-')
    year = datet[0].strip()
    month = datet[1]
    day = datet[2]
    julday = UTCDateTime(int(year), int(month), int(day)).julday
    julday = "%03d" % julday
    sdspath = os.path.join(os.getcwd(), year)

    if station.loc == '--':
        station.loc = ''
    streamData = station.net + '.' + station.sta + '.' + station.loc + '.' + station.comp + '.D.' + str(year) + '.' + str(julday)
    
    entry = os.path.join(sdspath, station.net, station.sta, station.comp + '.D', streamData)
    st = read(entry, format="MSEED", starttime=tw['start'], endtime=tw['end'], nearest_sample=True)
    if len(st.getGaps()) > 0:
        st.merge(method=0, fill_value='interpolate', interpolation_samples=0)
    
    stream = filterWaveform(st)
    print 'OVERALL: ', stream
    return stream
# --------------------------------------------------------------------------------------------------


def traveltimes(MetaDict, Event):

    logger.info('\033[31m Enter AUTOMATIC FILTER \033[0m')
    T = []
    Wdict = {}
    SNR = {}
    for i in MetaDict:
        
        de = locations2degrees(float(Event.lat), float(Event.lon), float(i.lat), float(i.lon))
        tt = getTravelTimes(delta=de, depth=float(Event.depth), model='ak135')
        if tt[0]['phase_name'] == 'P':
                        ptime = tt[0]['time']
                        T.append(ptime)

        logger.info('\033[31m \n\n+++++++++++++++++++++++++++++++++++++++++++++++++++ \033[0m') 
        print i.getName(), i.lat, i.lon, ptime
        ttime = ptime
        tw = calculateTimeWindows(ptime, Event)
        w, snr = readWaveformsCross(i, tw, Event, ttime)
        Wdict[i.getName()] = w
        SNR[i.getName()] = snr
    
    logger.info('\033[31m Exit AUTOMATIC FILTER \033[0m')
    
    return Wdict, SNR
# --------------------------------------------------------------------------------------------------

def writeStream(StreamDict):
    
    for s in StreamDict.iterkeys():
        fname = s + '.old'
        StreamDict[s].write(fname, format='MSEED')
        
def writeStreamFull(StreamDict):
    
    for s in StreamDict.iterkeys():
        fname = s + '.full'
        StreamDict[s].write(fname, format='MSEED')

def writeStreamREF(StreamDict):
    
    for s in StreamDict.iterkeys():
        fname = s + '.ref'
        StreamDict[s].write(fname, format='MSEED')

def f6(d1):
    return max(d1, key=d1.get)
# --------------------------------------------------------------------------------------------------


def doXcorr(StreamDict, SNRDict):
    
    t = f6(SNRDict)
    print 'REFERENCE: ',t
    for i in SNRDict.iterkeys():
        print 'STREAM: ',i,' SNR: ', SNRDict[i]
    
    corrDict = {}
    ref = StreamDict[t][0].data
    
    logger.info('\033[31m Enter Xcorr Procedure \033[0m')
    
    for stream in StreamDict.iterkeys():
        
        a, b = obspy.signal.cross_correlation.xcorr(ref, StreamDict[stream][0], 120)
        shift = a / StreamDict[stream][0].stats.sampling_rate
        corrDict[stream] = Corr(shift, b, a) 
        print 'Index: ', a, ' Value: ', b, ' ----> ', stream, StreamDict[stream][0].stats.sampling_rate 
        
    return corrDict,StreamDict[t]
# --------------------------------------------------------------------------------------------------

def writeShift(ShiftList):

    fobj = open('shift.dat', 'w')
    fobj.write('\n'.join(ShiftList))
    fobj.close()
# --------------------------------------------------------------------------------------------------    
    

def shiftSeismograms(StreamDict, XcorrDict,pickerShift):
    
    L = []
    S = []
    dsfactor=0.55
    for stream in StreamDict.iterkeys():
        for shift in XcorrDict.iterkeys():
                if stream == shift:
                    StreamDict[stream][0].stats.starttime = StreamDict[stream][0].stats.starttime + XcorrDict[shift].index+pickerShift
                    if XcorrDict[shift].value < 0:
                    	StreamDict[stream][0].data = -1 * StreamDict[stream][0].data
                    if abs(XcorrDict[shift].value) > dsfactor:
                        fname = stream + '.new'
                        StreamDict[stream].write(fname, format='MSEED')
                        print stream, XcorrDict[shift].value, XcorrDict[shift].sampleindex, XcorrDict[shift].index
                        if XcorrDict[shift].value < 0:
                           t = -1
                        else:
                           t = 1
                        info = '{0:10} {1:10} {2:4}'.format(stream, XcorrDict[shift].index,t)
                        L.append(info)
                        S.append(stream)
                    if abs(XcorrDict[shift].value) < dsfactor:
                        print 'OUT: ', stream, XcorrDict[shift].value, XcorrDict[shift].sampleindex, XcorrDict[shift].index
    
    return L, S
# --------------------------------------------------------------------------------------------------

def shiftReferenceSeismograms(StreamDict, XcorrDict):
    
    L = []
    S = []
    dsfactor=0.5

    for stream in StreamDict.iterkeys():
        for shift in XcorrDict.iterkeys():
                if stream == shift:
                    StreamDict[stream][0].stats.starttime = StreamDict[stream][0].stats.starttime + XcorrDict[shift].index

                    if XcorrDict[shift].value < 0:
                        StreamDict[stream][0].data = -1 * StreamDict[stream][0].data

                    if abs(XcorrDict[shift].value) > dsfactor:
                        fname = stream + '.new'
                        StreamDict[stream].write(fname, format='MSEED')
                        print stream, XcorrDict[shift].value, XcorrDict[shift].sampleindex, XcorrDict[shift].index
                        if XcorrDict[shift].value < 0:
                           t = -1
                        else:
                           t = 1
                        info = '{0:10} {1:10} {2:4}'.format(stream, XcorrDict[shift].index,t)
                        L.append(info)
                        S.append(stream)

                    if abs(XcorrDict[shift].value) < dsfactor:
                        print 'OUT: ', stream, XcorrDict[shift].value, XcorrDict[shift].sampleindex, XcorrDict[shift].index
    
    return L, S
# --------------------------------------------------------------------------------------------------

def getSCstationList(StreamList,net):
    
    t = ''

    for i in StreamList:
        t += i + '|'
        
    print t,len(StreamList),' of ',len(net)
# --------------------------------------------------------------------------------------------------

def searchMeta(sname,Metalist):
    
    for i in Metalist:
            mname = ('%s.%s.%s.%s')%(i.net,i.sta,i.loc,i.comp)

            if sname == i.getName():
                  return i
# --------------------------------------------------------------------------------------------------

def refTrigger(Waveform,Event,Meta):
    
    print Event
    name = ('%s.%s.%s.%s')%(Waveform[0].stats.network,Waveform[0].stats.station,Waveform[0].stats.location,Waveform[0].stats.channel)
    
    i = searchMeta(name,Meta)
    print i
    
    de = locations2degrees(float(Event.lat), float(Event.lon), float(i.lat), float(i.lon))
    tt = getTravelTimes(delta=de, depth=float(Event.depth), model='ak135')
    ptime = 0

    if tt[0]['phase_name'] == 'P':
       ptime = tt[0]['time']
    
    tw  = calculateTimeWindows (ptime, Event)
    stP = readWaveformsPicker  (i, tw, Event, ptime)
    trP = stP[0]
    
    cft = recSTALTA(trP.data, int(1 * trP.stats.sampling_rate), int(10 * trP.stats.sampling_rate))
    t = triggerOnset(cft,6,1.5)
    print len(trP),t,type(t)
    onset = t[0][0]/trP.stats.sampling_rate
    
    print 'TRIGGER ',trP.stats.starttime+onset
    print 'THEORETICAL: ',UTCDateTime(Event.time)+ptime
    tdiff = (trP.stats.starttime+onset)-(UTCDateTime(Event.time)+ptime) 
    #plotTrigger(trP,cft,6,1.5)
    print tdiff
    
    return tdiff
# --------------------------------------------------------------------------------------------------

if __name__ == '__main__':

    assert False                                       #hs

    fname = 'metainfo-2008-279.meta'
    orname = 'KYRGYZSTAN_2008-10-05_15-52-52.000000.origin'
    arrayconfig = 'KYRGYZSTAN_2008-10-05_15-52-52.000000.config'

    META = readMetaInfoFile(fname)
    Origin = parseConfig(orname)
    Config = parseConfig(arrayconfig)
    ev = Event(Origin['lat'], Origin['lon'], Origin['depth'], Origin['time'])
    
    L = []
    #REF = []
    ne = Config['networks'].split(',')
    
    for i in ne:
        network = Config[i].split('|')
        
        A = filterStations(META, Config, Origin, network)
        WD, SNR = traveltimes(A, ev)
        writeStream(WD)
        CD,ref = doXcorr(WD, SNR)
        onset = refTrigger(ref,ev,META)
        
        #REF.append(ref)
        
        C, Streams = shiftSeismograms(WD, CD,onset)
        L.extend(C)
    
    writeShift(L)
    getSCstationList(Streams,network)
    '''
    print '###########################################REF##########################################################'
    #start shifting of reference traces
     
    L = []#list for saving reference trace per array
    for i in REF:
        print 'REF '
        name = ('%s.%s.%s.%s')%(i[0].stats.network,i[0].stats.station,i[0].stats.location,i[0].stats.channel)
        t = searchMeta(name,META)
        L.append(t)

    WD, SNR = traveltimes(L, ev)
    writeStreamREF(WD)
    CD,ref = doXcorr(WD, SNR)
    C, Streams = shiftReferenceSeismograms(WD, CD)
    '''
