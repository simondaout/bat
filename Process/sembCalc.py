
import os
import sys
import platform

WINDOWS = (platform.system() == 'Windows')

import math
from   math     import radians, cos, sin, atan2
import logging
import time

# add local directories to import path

sys.path.append ('../tools/')                     
sys.path.append ('../Common/')
 
#
from  obspy.core.utcdatetime import UTCDateTime

#    Import from NumPy
import numpy as np          

#    Import from common

import  Basic
import  Globals
import  Logfile
import  DataTypes
from    DataTypes  import Location
from    ObspyFkt   import loc2degrees 
from    ConfigFile import ConfigObj, OriginCfg 

#    Import from Process

import  trigger

#import Csemblance                     # C-Code - unused
#import Ctrigger
#import sembPar                        # C-Code - unused

#USE_C_CODE =  True                    #16.12.2015
USE_C_CODE  =  False                   #16.12.2015

if USE_C_CODE : 
   import  CTrig                       # C-Code
   import  Cm                          # C-Code

else : 
   from semp import otest      

# -------------------------------------------------------------------------------------------------

logger = logging.getLogger('ARRAY-MP')

class SembMax (object):
    '''
    class to store sembmax object for each grid point
    '''
    def __init__(self,lat,lon,semb):

        self.lat  = lat
        self.lon  = lon
        self.semb = semb

# -------------------------------------------------------------------------------------------------

class FileSembMax(object):
    '''
    class to strore sembmax object for the sembmaxvalue file
    '''
    def __init__(self,istep,sembmaxX,sembmaxY,sembmax,usedarrays,delta,azi,deltakm):

        self.istep      = istep
        self.sembmaxX   = sembmaxX
        self.sembmaxY   = sembmaxY
        self.sembmax    = sembmax
        self.usedarrays = usedarrays
        self.delta      = delta
        self.azi        = azi
        self.deltakm    = deltakm

    def get(self):
        return ('%d %.2f %.2f %f %d %03f %f %03f\n' % (self.istep,self.sembmaxX,self.sembmaxY,
                                                       self.sembmax,self.usedarrays,self.delta,
                                                       self.azi,self.delta*119.19))

# -------------------------------------------------------------------------------------------------

def toAzimuth (latevent,lonevent,latsource,lonsource):
        '''
        method to calculate azimuth between two points
        '''
        # Convert to radians.
        lat1 = radians (latsource);
        lon1 = radians (lonsource);
        lat2 = radians (latevent);
        lon2 = radians (lonevent);

       # Compute the angle.
        x     =  sin(lon1-lon2 ) * cos(lat2);
        y     =  cos(lat1) * sin(lat2) - sin(lat1) * cos(lat2) * cos(lon1-lon2);
        angle = -atan2 (x,y);

        if angle < 0.0 :
         angle  += math.pi * 2.0;

       #And convert result to degrees.
        angle = math.degrees (angle)
        angle = '%02f'%angle

        return angle;
    
# -------------------------------------------------------------------------------------------------

def writeSembMatricesSingleArray (SembList,Config,Origin,arrayfolder,ntimes):
    '''
    method to write semblance matrizes from one processes to file for each timestep
    '''
    logger.info ('start write semblance matrices')
    
    cfg    = ConfigObj (dict=Config)
    origin = OriginCfg (Origin)

    dimX   = cfg.dimX()         # ('dimx')
    dimY   = cfg.dimY()         # ('dimy')
    winlen = cfg.winlen ()      # ('winlen')
    step   = cfg.step()         # ('step')

    latv   = []
    lonv   = []

    gridspacing = cfg.Float ('gridspacing')
    migpoints   = dimX * dimY

    o_lat       = origin.lat()         # float (Origin['lat'])
    o_lon       = origin.lon()         # float (Origin['lon'])
    oLatul      = 0
    oLonul      = 0
    
    z=0

    for i in xrange(dimX):
         oLatul = o_lat - ((dimX-1)/2) * gridspacing + i*gridspacing

         if z == 0 and i == 0:
             Latul = oLatul
         o=0    

         for j in xrange (dimY):
               oLonul = o_lon - ((dimY-1)/2) * gridspacing + j * gridspacing
               
               if o==0 and j==0:  Lonul = oLonul
               
               latv.append (oLatul)
               lonv.append (oLonul)
    #endfor

    rc  = UTCDateTime (Origin['time'])
    rcs = '%s-%s-%s_%02d:%02d:%02d'% (rc.day,rc.month,rc.year, rc.hour,rc.minute,rc.second)
    d   = rc.timestamp
    
    for a, i in enumerate(SembList):
        #logger.info('timestep %d' % a)

        fobj = open (os.path.join (arrayfolder,'%s_%03d.ASC' % (Origin['depth'],a)),'w')
        fobj.write ('# %s , %s\n' % (d,rcs))
        fobj.write ('# step %ds| ntimes %d| winlen: %ds\n' % (step,ntimes,winlen))
        fobj.write ('# \n')
        fobj.write ('# southwestlat: %.2f dlat: %f nlat: %f \n'%(Latul,gridspacing,dimX))
        fobj.write ('# southwestlon: %.2f dlon: %f nlon: %f \n'%(Lonul,gridspacing,dimY))
        fobj.write ('# ddepth: 0 ndepth: 1 \n')
        
        for j in range (migpoints):
            x    = latv[j]
            y    = lonv[j]
            z    = origin.depth()         # float(Origin['depth'])
            semb = i[j]

            fobj.write ('%.2f %.2f %.2f %.20f\n' % (x,y,z,semb))
        #endfor
            
        fobj.close()
    #endfor
         
# -------------------------------------------------------------------------------------------------

def collectSemb (SembList,Config,Origin,Folder,ntimes,arrays):
    '''
    method to collect semblance matrizes from all processes and write them to file for each timestep
    '''
    Logfile.add ('start collect in collectSemb')
    
    cfg    = ConfigObj (dict=Config)
    origin = ConfigObj (dict=Origin)

    dimX   = cfg.dimX()         # ('dimx')
    dimY   = cfg.dimY()         # ('dimy')
    winlen = cfg.winlen ()      # ('winlen')
    step   = cfg.step()         # ('step')

    latv        = []
    lonv        = []

    gridspacing = cfg.Float ('gridspacing')
    migpoints   = dimX * dimY
    o_lat       = origin.lat()         # float (Origin['lat'])
    o_lon       = origin.lon()         # float (Origin['lon'])
    oLatul      = 0
    oLonul      = 0
    
    z=0

    for i in xrange(dimX):
         oLatul = o_lat - ((dimX-1)/2) * gridspacing + i*gridspacing

         if z == 0 and i == 0 :
             Latul = oLatul
         o=0    

         for j in xrange (dimY):
               oLonul = o_lon - ((dimY-1)/2) * gridspacing + j*gridspacing
               
               if o==0 and j==0:
                    Lonul = oLonul
               
               latv.append (oLatul)
               lonv.append (oLonul)
    #endfor  
  
    #print 'SL: ',SembList, type(SembList),type(SembList[0]),SembList[0].ndim

    tmp=1
    for a in SembList:
        tmp *= a

    #sys.exit()
    
    sembmaxvaluev = np.ndarray (ntimes,dtype=float)
    sembmaxlatv   = np.ndarray (ntimes,dtype=float)
    sembmaxlonv   = np.ndarray (ntimes,dtype=float)
    
    rc         = UTCDateTime(Origin['time'])
    rcs        = '%s-%s-%s_%02d:%02d:%02d'% (rc.day,rc.month,rc.year, rc.hour,rc.minute,rc.second)
    d          = rc.timestamp
    usedarrays = 5
    
    folder      = Folder['semb']
    fobjsembmax = open (os.path.join (folder,'sembmax.txt'),'w')
    
    for a, i in enumerate(tmp):
        logger.info('timestep %d' % a)
        
        fobj  = open (os.path.join (folder, '%s_%03d.ASC' % (Origin['depth'],a)),'w')
        #fobj = open (os.path.join (folder, '%03d.ASC'    % a),'w')

        fobj.write ('# %s , %s\n' % (d,rcs))
        fobj.write ('# step %ds| ntimes %d| winlen: %ds\n' % (step,ntimes,winlen))
        fobj.write ('# \n')
        fobj.write ('# southwestlat: %.2f dlat: %f nlat: %f \n'%(Latul,gridspacing,dimX))
        fobj.write ('# southwestlon: %.2f dlon: %f nlon: %f \n'%(Lonul,gridspacing,dimY))
        fobj.write ('# ddepth: 0 ndepth: 1 \n')
        
        
        sembmax  = 0
        sembmaxX = 0
        sembmaxY = 0
        
        origin = DataTypes.dictToLocation (Origin)

        for j in range(migpoints):
            x    = latv[j]
            y    = lonv[j]
            semb = i[j]

            fobj.write ('%.2f %.2f %.20f\n' % (x,y,semb))
            
            if  semb > sembmax:
                sembmax  = semb;# search for maximum and position of maximum on semblance grid for given time step         
                sembmaxX = x;
                sembmaxY = y;
        #endfor

        #print sembmax,sembmaxX,sembmaxY
        
        delta = loc2degrees (Location (sembmaxX, sembmaxY), origin)
        azi   = toAzimuth   (float(Origin['lat']), float(Origin['lon']),float(sembmaxX), float(sembmaxY))

        sembmaxvaluev[a] = sembmax
        sembmaxlatv[a]   = sembmaxX
        sembmaxlonv[a]   = sembmaxY
        
        fobjsembmax.write ('%d %.2f %.2f %.20f %d %03f %f %03f\n' % (a*step,sembmaxX,sembmaxY,sembmax,usedarrays,delta,float(azi),delta*119.19))
        fobj.close()
    #endfor

    fobjsembmax.close()

    durationpath  = os.path.join (folder, "duration.txt")
    #print 'DD: ',durationpath
    trigger.writeSembMaxValue (sembmaxvaluev,sembmaxlatv,sembmaxlonv,ntimes,Config,Folder)
    print 'DD2: ',durationpath
    #print sembmaxvaluev,sembmaxlatv,sembmaxlonv,len(sembmaxvaluev),len(sembmaxlatv),len(sembmaxlonv)
    
    #print 'enter new trigger'
    #trigger.semblancestalta(sembmaxvaluev,sembmaxlatv,sembmaxlonv)
    #sys.exit()
    
    #Ctrigger.staltatriggering (ntimes,step,sembmaxlatv,sembmaxlonv,sembmaxvaluev,durationpath)

   #CTrig.staltatriggering  (ntimes,step,sembmaxlatv,sembmaxlonv,sembmaxvaluev,durationpath) #hs
    trigger.semblancestalta (sembmaxvaluev,sembmaxlatv,sembmaxlonv)                          #hs
# -------------------------------------------------------------------------------------------------

#def doCalc (flag,Config,WaveformDict,FilterMetaData,Gmint,Gmaxt,TTTGridMap,Folder,Origin, ntimesstart,ntimesend):
def  doCalc (flag,Config,WaveformDict,FilterMetaData,Gmint,Gmaxt,TTTGridMap,Folder,Origin, ntimes):
    '''
    method for calculating semblance of one station array
    '''
    Logfile.add ('PROCESS %d %s' % (flag,' Enters Semblance Calculation') )
    Logfile.add ('MINT  : %f  MAXT: %f Traveltime' % (Gmint,Gmaxt))

    cfg = ConfigObj (dict=Config)

    dimX   = cfg.dimX()         # ('dimx')
    dimY   = cfg.dimY()         # ('dimy')
    winlen = cfg.winlen ()      # ('winlen')
    step   = cfg.step()         # ('step')

    new_frequence   = cfg.newFrequency()          #  ('new_frequence')
    forerun         = cfg.Int  ('forerun')
    duration        = cfg.Int  ('duration')
    gridspacing     = cfg.Float('gridspacing')

    nostat          = len (WaveformDict)
    traveltimes     = {}
    recordstarttime = ''
    minSampleCount  = 999999999

    ntimes = int ((forerun + duration)/step)
    nsamp  = int (winlen * new_frequence)
    nstep  = int (step   * new_frequence)
    
    #for i in WaveformDict.iterkeys():
    #    print i,WaveformDict[i]
    
    ############################################################################
    calcStreamMap = WaveformDict

    for trace in calcStreamMap.iterkeys():
        recordstarttime = calcStreamMap[trace].stats.starttime
        d = calcStreamMap[trace].stats.starttime
        d = d.timestamp

        if calcStreamMap[trace].stats.npts < minSampleCount:
            minSampleCount = calcStreamMap[trace].stats.npts

    ############################################################################
    traces     = np.ndarray (shape=(len(calcStreamMap), minSampleCount), dtype=float)
    traveltime = np.ndarray (shape=(len(calcStreamMap), dimX*dimY), dtype=float)
    latv       = np.ndarray (dimX*dimY, dtype=float)
    lonv       = np.ndarray (dimX*dimY, dtype=float)
    ############################################################################

    #traces      = np.ndarray (nostat*minSampleCount,dtype=float)
    #traveltimes = np.ndarray (nostat*dimX*dimY,dtype=float)
    #latv        = np.ndarray (dimX*dimY,dtype=float)
    #lonv        = np.ndarray (dimX*dimY,dtype=float)
    #print 'minSC: ',minSampleCount,' LCSM: ',len(calcStreamMap)
    c=0
    streamCounter = 0
    
    for key in calcStreamMap.iterkeys():
        streamID = key
        c2       = 0
        #print streamID, len(calcStreamMap[key]),minSampleCount

        for o in calcStreamMap[key]:
            if c2 < minSampleCount:
                traces[c][c2] = o
                #print 'C: ',c,' C2: ',c2,' TRACES:',traces[c][c2]
                c2 += 1
        #endfor

        for key in TTTGridMap.iterkeys():
            
            if streamID == key:
                traveltimes[streamCounter] = TTTGridMap[key]
            else:
                "NEIN", streamID, key
        #endfor

        if not streamCounter in traveltimes : 
           continue                              #hs : thread crashed before

        g     = traveltimes[streamCounter]
        dimZ  = g.dimZ
        mint  = g.mint
        maxt  = g.maxt
        Latul = g.Latul
        Lonul = g.Lonul
        Lator = g.Lator
        Lonor = g.Lonor
        
        gridElem = g.GridArray
        
        for x in range(dimX):
            for y in range(dimY):
                elem = gridElem[x, y]
                
                traveltime [c][x * dimY + y] = elem.tt
                latv [x * dimY + y] = elem.lat
                lonv [x * dimY + y] = elem.lon
        #endfor

        c += 1
        streamCounter += 1
    #endfor

    ############################## CALCULATE PARAMETER FOR SEMBLANCE CALCULATION ##################
    nsamp     = winlen * new_frequence
    nstep     = int (step*new_frequence)
    migpoints = dimX * dimY

    dimZ = 0
    new_frequence = cfg.newFrequency ()              # ['new_frequence']
    maxp = int (Config['ncore'])
   #maxp = 20                                        #hs

    Logfile.add ('PROCESS %d  NTIMES: %d' % (flag,ntimes))
  
    #k = Csemblance.semb(flag,nostat,nsamp,ntimes,nstep,Gmint,Gmaxt,Lonul,Latul,minSampleCount,
    #                    dimZ,dimX,dimY,new_frequence,ntimesstart,ntimesend,winlen,step,gridspacing,
    #                    latv,lonv,traveltime,traces,backveclen)
    #k = sembPar.semb   (flag,nostat,nsamp,ntimes,nstep,Gmint,Gmaxt,Lonul,Latul,minSampleCount,dimZ,
    #                    dimX,dimY,new_frequence,ntimesstart,ntimesend,winlen,step,gridspacing,latv,
    #                    lonv,traveltime,traces,backveclen)
    
    if False :
       print ('nostat ',nostat,type(nostat))
       print ('nsamp ',nsamp,type(nsamp))
       print ('ntimes ',ntimes,type(ntimes))
       print ('nstep ',nstep,type(nstep))
       print ('dimX ',dimX,type(dimX))
       print ('dimY ',dimY,type(dimY))
       print ('mint ',Gmint,type(mint))
       print ('new_freq ',new_frequence,type(new_frequence))
       print ('minSampleCount ',minSampleCount,type(minSampleCount))
       print ('latv ',latv,type(latv))
       print ('traces',traces,type(traces))
       print ('traveltime',traveltime,type(traveltime))

    traveltime = traveltime.reshape (1,nostat*dimX*dimY)
    traces     = traces.reshape     (1,nostat*minSampleCount)
    #print 'traveltime2',traveltime,type(traveltime)
    
    t1 = time.time()

    if USE_C_CODE :
       k  = Cm.otest (maxp,nostat,nsamp,ntimes,nstep,dimX,dimY,Gmint,new_frequence,
                      minSampleCount,latv,lonv,traveltime,traces)
    else :
       k = otest (maxp,nostat,nsamp,ntimes,nstep,dimX,dimY,Gmint,new_frequence,
                  minSampleCount,latv,lonv,traveltime,traces)                       #hs

    t2 = time.time()
    
    Logfile.add ('%s took %0.3f s' % ('CALC:', (t2-t1)))
    #print 'K',k,len(k),' MUST ',ntimes*dimX*dimY,' RES ',k[1]
 
    partSemb = k

    #partSemb = partSemb.reshape (1,migpoints)
    partSemb  = partSemb.reshape (ntimes,migpoints)

    #print 'PARTSEMB FLAG: ',partSemb,type(partSemb),partSemb.ndim
    
    return partSemb
