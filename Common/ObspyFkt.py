
import platform

WINDOWS = (platform.system() == 'Windows')

import  obspy.core
#from   obspy.taup.taup            import locations2degrees
from    obspy.core.util            import locations2degrees
#from   obspy.core.util.geodetics  import gps2DistAzimuth    

from  obspy.taup.taup              import getTravelTimes


import  Basic

def loc2degrees (a, b) :

    if type(a) is dict : a1 = Basic.dictToLocation (a)
    else :               a1 = a

    if type(b) is dict : b1 = Basic.dictToLocation (b)
    else :               b1 = b

    delta = locations2degrees (float (a1.lat), float(a1.lon), float (b1.lat), float (b1.lon))
    return delta

# -------------------------------------------------------------------------------------------------

def obs_TravelTimes (delta1, depth1) :

   #Basic.checkFileExists ('ak135.model', isAbort=True)                          #17.12.2015
    return getTravelTimes (delta=delta1, depth = float (depth1), model='ak135')

def obs_kilometer2degrees (km) :

    return kilometer2degrees (float (km))