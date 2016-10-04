
import os
import sys

# add local directories to import path                      #hs+

sys.path.append ('../tools/')                     
sys.path.append ('../Common/')
                                                            #hs-
import logging
import fnmatch

from ConfigParser import SafeConfigParser

from   obspy.core import read
from   obspy.taup.taup import getTravelTimes
#from   obspy.taup.taup import locations2degrees       #hs
#from  obspy.core.util import locations2degrees
from   obspy.core.utcdatetime import UTCDateTime
import obspy.signal.cross_correlation
from   obspy.signal.trigger import pkBaer
from   obspy.core import *

from config import Station,Event

#    Import from NumPy

import numpy as np          

#    Import from common

import  Basic
import  Globals
import  Logfile
from    ObspyFkt   import loc2degrees 

#    Import from Process

#

logger = logging.getLogger('ARRAY-MP')


class Corr(object):
    def __init__(self, shift, value, sampleindex):
        self.shift = shift
        self.value = value
        self.sampleindex = sampleindex


class Xcorr(object):
    def __init__(self,Origin,StationMeta,EventPath,Config):
        self.Origin = Origin
        self.StationMeta = StationMeta
        self.EventPath = EventPath
        self.Config = Config


    def calculateTimeWindows(self,mint):

        tw = {}
        st = str(self.Origin.time)[:-1]

        tw['start'] = UTCDateTime(UTCDateTime(st) + (mint - 100))
        tw['end']   = tw['start'] + 200
    
        tw['xcorrstart'] = UTCDateTime(UTCDateTime(st) + (mint - 8))
        tw['xcorrend']   = tw['xcorrstart'] + 20

        Logfile.add (' ORIGIN TIME %s' % UTCDateTime (self.Origin.time))
        Logfile.add (' OVERALL TIME WINDOW: %s - %s' % (tw['start'], tw['end']))
        Logfile.add (' XCROSS TIME WINDOW: %s - %s'  % (tw['xcorrstart'], tw['xcorrend']))

        return tw

    def resampleWaveform(self,Waveform,end_frequence):

        new_frequence = end_frequence

        if Waveform.stats.sampling_rate == end_frequence:
            return Waveform
        
        if int(Waveform.stats.sampling_rate) == 20:
            oldsr  = Waveform.stats.sampling_rate
            factor = int(Waveform.stats.sampling_rate / new_frequence + 0.5)

            Logfile.add ('downsampling %s from %d Hz to %d Hz with factor %d' % (Waveform, oldsr, new_frequence, factor))
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
        tm.stats.network       = Waveform.stats.network
        tm.stats.station       = Waveform.stats.station
        tm.stats.sampling_rate = Waveform.stats.sampling_rate*4
        tm.stats.channel       = Waveform.stats.channel
        tm.stats.starttime     = utcdatetime.UTCDateTime(Waveform.stats.starttime)
        tm.stats._format       = 'MSEED'
        
        oldsr  = tm.stats.sampling_rate
        factor = int(tm.stats.sampling_rate / new_frequence + 0.5)

        Logfile.add ('downsampling %s from %d Hz to %d Hz with factor after upsampling %d' % (tm, oldsr, new_frequence, factor))
        tm.decimate(int(factor), strict_length=False, no_filter=False)
        return tm
    
    if Waveform.stats.sampling_rate == 40:
        oldsr  = Waveform.stats.sampling_rate
        factor = int(Waveform.stats.sampling_rate / new_frequence + 0.5)

        Logfile.add ('downsampling %s from %d Hz to %d Hz with factor %d' % (Waveform, oldsr, new_frequence, factor))
        Waveform.decimate(int(factor), strict_length=False, no_filter=False)

        return Waveform
    
    if Waveform.stats.sampling_rate == 50:
        #upsampling with factor = 2
        data   = numpy.zeros((Waveform.stats.npts))
        test   = numpy.fft.rfft(Waveform.data)
        zero   = numpy.fft.rfft(data)
        t      = numpy.concatenate((test,zero))
        tbdata = numpy.fft.irfft(t)
        
        tm = Trace(tbdata)
        tm.stats.network       = Waveform.stats.network
        tm.stats.station       = Waveform.stats.station
        tm.stats.sampling_rate = Waveform.stats.sampling_rate*2
        tm.stats.channel       = Waveform.stats.channel
        tm.stats.starttime     = utcdatetime.UTCDateTime(Waveform.stats.starttime)
        tm.stats._format       = 'MSEED'
        
        oldsr  = tm.stats.sampling_rate
        factor = int(tm.stats.sampling_rate / new_frequence + 0.5)

        Logfile.add ('downsampling %s from %d Hz to %d Hz with factor after upsampling %d' % (tm, oldsr, new_frequence, factor))
        tm.decimate(int(factor), strict_length=False, no_filter=False)
        return tm

    
    if Waveform.stats.sampling_rate == 80:
        oldsr  = Waveform.stats.sampling_rate
        factor = int(Waveform.stats.sampling_rate / new_frequence + 0.5)

        Logfile.add ('downsampling %s from %d Hz to %d Hz with factor %d' % (Waveform, oldsr, new_frequence, factor))
        Waveform.decimate(int(factor), strict_length=False, no_filter=False)
        
    if Waveform.stats.sampling_rate == 100:
        oldsr  = Waveform.stats.sampling_rate
        factor = int(Waveform.stats.sampling_rate / new_frequence + 0.5)

        Logfile.add ('downsampling %s from %d Hz to %d Hz with factor %d' % (Waveform, oldsr, new_frequence, factor))
        Waveform.decimate(int(factor), strict_length=False, no_filter=False)

        return Waveform
    
    else:
        oldsr  = Waveform.stats.sampling_rate
        factor = int(Waveform.stats.sampling_rate / new_frequence + 0.5)

        Logfile.add ('downsampling %s from %d Hz to %d Hz with factor %d' % (Waveform, oldsr, new_frequence, factor))
        Waveform.decimate(int(factor), strict_length=False, no_filter=False)
        return Waveform

    # ---------------------------------------------------------------------------------------------

    def filterWaveform(self,Waveform):
        logger.info('\033[31m Filter Waveform: \033[0m' % ())
    
        new_frequence = self.Config['new_frequence']
    
        st = Stream()
    
        for i in Waveform:
        
            j = resampleWaveform(i,new_frequence)
            j.detrend("simple")
            j.filter("bandpass", freqmin=0.4, freqmax=3, corners=3, zerophase=False)
            
            st.append(j)
    
        return st

    # ---------------------------------------------------------------------------------------------

    def signoise(self,Waveform, ttime, path):
    
        st = str(self.Origin.time)[:-1]
    
        ponset = UTCDateTime(st) + ttime
    
        winnoise_start = Waveform.stats.starttime
        winnoise_end   = ponset - 10
        winsig_start   = ponset - 2
        winsig_end     = ponset + 5

        try:
            winnoise = read(path, format="MSEED", starttime=winnoise_start, endtime=winnoise_end, nearest_sample=True)
            winsig   = read(path, format="MSEED", starttime=winsig_start, endtime=winsig_end, nearest_sample=True)
        except Exception, e:
            print e
        

        psignal = abs(winsig.max()[0])
        pnoise = abs(winnoise.max()[0])

        signoise = float(psignal) / float(pnoise)
        print psignal, pnoise, signoise
    
        return signoise

    # ---------------------------------------------------------------------------------------------

    def maxAmplitude(self,Waveform, ttime, Origin, path):
        st = str(Origin.time)[:-1]
    
        ponset = UTCDateTime(st) + ttime
        winamp_start = ponset
        winamp_end = ponset + 10

        winamp = read(path, format="MSEED", starttime=winamp_start, endtime=winamp_end, nearest_sample=True)
    
        id = winamp[0].stats.network+'.'+winamp[0].stats.station+'.'+winamp[0].stats.location+'.'+winamp[0].stats.channel
    
        print 'MAXAMP: ',winamp.max(),id
    
        return abs(winamp.max()[0])
    
    # ---------------------------------------------------------------------------------------------

    def readWaveformsCross(self,station, tw, ttime):
    
        time = self.Origin.time
        ts = time.split('T')
    
        datet = ts[0]
        datet = datet.split('-')
        year = datet[0].strip()
        month = datet[1]
        day = datet[2]
        #timep = ts[1][:-1]
     
        #print time,ts,year,month,day
        julday = UTCDateTime(int(year), int(month), int(day)).julday
        julday = "%03d" % julday

        sdspath = os.path.join(self.EventPath,'data', year)
        stream = ''
        snr = ''

        if station.loc == '--':
            station.loc = ''

        streamData = station.net + '.' + station.sta + '.' + station.loc + '.' + station.comp + '.D.' + str(year) + '.' + str(julday)
        
        entry = os.path.join(sdspath, station.net, station.sta, station.comp + '.D', streamData)
        #print entry

        st = read(entry, format="MSEED", starttime=tw['start'], endtime=tw['end'], nearest_sample=True)
        #print st,st.getGaps()

        if len(st.getGaps()) > 0:
            st.merge(method=0, fill_value='interpolate', interpolation_samples=0)

    
        snr = self.signoise(st[0], ttime, entry)
        #amp = self.maxAmplitude(st[0], ttime, Origin, entry)
    
        stream = self.filterWaveform(st)
        print 'OVERALL: ', st
        st.trim(tw['xcorrstart'], tw['xcorrend'])
        print 'CROSSWINDOW: ', st
        st[0].stats.starttime = UTCDateTime(1000)
    
        return stream, snr

    # ---------------------------------------------------------------------------------------------

    def traveltimes(self):

        logger.info('\033[31m Enter AUTOMATIC FILTER \033[0m')
        T = []
        Wdict = {}
        SNR = {}

        for i in self.StationMeta:
        
           #de = locations2degrees (float(self.Origin.lat), float(self.Origin.lon), float(i.lat), float(i.lon))
            de = loc2degrees       (self.Origin, i)
            tt = getTravelTimes(delta=de, depth=float(self.Origin.depth), model='ak135')

            if tt[0]['phase_name'] == 'P':
                        ptime = tt[0]['time']
                        T.append(ptime)

            logger.info('\033[31m \n\n+++++++++++++++++++++++++++++++++++++++++++++++++++ \033[0m') 
            print i.getName(), i.lat, i.lon, ptime
            ttime = ptime

            tw = self.calculateTimeWindows(ptime)
            try:
                w, snr = self.readWaveformsCross(i, tw, ttime)
                Wdict[i.getName()] = w
                SNR[i.getName()] = snr
            except:
                continue
    
            logger.info('\033[31m Exit AUTOMATIC FILTER \033[0m')

        return Wdict, SNR

    def f6(self,d1):
        return max(d1, key=d1.get)

    # ---------------------------------------------------------------------------------------------

    def doXcorr(self):
        
        StreamDict, SNRDict = self.traveltimes()
        
        t = self.f6(SNRDict)
        #print t
        Logfile.red ('Reference Station for Xcorr Procedure %s' % (t))

        #for i in SNRDict.iterkeys():
        #    print 'STREAM: ',i,' SNR: ', SNRDict[i]
    
        corrDict = {}
        ref      = StreamDict[t][0].data
    
        logger.info('\033[31m Enter Xcorr Procedure \033[0m')
    
        for stream in StreamDict.iterkeys():
            #print stream, StreamDict[stream][0]
           a, b  = obspy.signal.cross_correlation.xcorr(ref, StreamDict[stream][0], 50)
           shift = a / StreamDict[stream][0].stats.sampling_rate
           corrDict[stream] = Corr(shift, b, a) 
           xxx
           print 'Index: ', a, ' Value: ', b, ' ----> ', stream, StreamDict[stream][0].stats.sampling_rate , ' SHIFT IN TIME: ',shift
        #endfor

        logger.info ('\033[31m Finish Xcorr Procedure \033[0m')
        return corrDict

    # ---------------------------------------------------------------------------------------------

    def shiftSeismograms(self,StreamDict, XcorrDict):
    
        L = []
        S = []
        dsfactor=0.6

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

# -------------------------------------------------------------------------------------------------

def cmpFilterMetavsXCORR(XcorrMeta, StationMetaList):
    
    FilterList = []
 
    for i in StationMetaList:
        for j in XcorrMeta.iterkeys():
            if i.getName() == j:
                FilterList.append(i) 
 
    logger.info('\033[31m Xcorr Procedure finished %d of %d stations left for processing \033[0m'%(len(FilterList),len(StationMetaList)))
    return FilterList
