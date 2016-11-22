
import sys, tempfile, os, glob, re
from pyrocko import iris_ws, util, catalog, model, weeding, pile, io, orthodrome, trace, gui_util, evalresp_ext
import logging
import cPickle as pickle
import numpy as num
from subprocess import call




def get_stations(site, lat, lon, rmin, rmax, tmin, tmax, channel_pattern='BH?'):
    from pyrocko.fdsn import ws

    extra = {}
    if site == 'iris':
        extra.update(matchtimeseries=True)

    sx = ws.station(site=site, latitude=lat, longitude=lon, minradius=rmin, 
                    maxradius=rmax, startbefore=tmin, endafter=tmax,
                    channel=channel_pattern, format='text', level='channel',
                    includerestricted=False, **extra)

    return sx.get_pyrocko_stations()


def iris_get_data(stations, tmin, tmax, fn_template, neach=20):
    from pyrocko.fdsn import ws
    fns = []
    i = 0
    stations = sorted(stations, key=lambda s: (s.network, s.station))
    while i < len(stations):
        stations_now = stations[i:i+neach]
        selection = ws.make_data_selection( stations_now, tmin, tmax )

        f = tempfile.NamedTemporaryFile()
        try:
            data = ws.dataselect(site='iris', selection=selection)
            while True:
                buf = data.read(1024)
                if not buf:
                    break
                f.write(buf)

            f.flush()

            fns.extend(io.save(io.iload(f.name), fn_template))

        except ws.EmptyResult:
            pass

        f.close()
        i += neach

    return fns


def iris_get_responses(stations, time, fn_template):
    for sta in stations:
        fn = fn_template % dict(network=sta.network, station=sta.station, location=sta.location)
        if not os.path.isfile(fn):
            try:
                fi = iris_ws.ws_resp(sta.network, sta.station, sta.location, '*', time=time)

                util.ensuredirs(fn)
                fo = open(fn, 'w')
                while True:
                    data = fi.read(1024)
                    if not data:
                        break
                    fo.write(data)

                fo.close()
                fi.close()
            except iris_ws.NotFound:
                pass
        
        for cha in sta.get_channels():
            class DummyTrace:
                pass

            tr = DummyTrace()
            tr.tmin = time
            tr.tmax = time
            tr.nslc_id = (sta.network, sta.station, sta.location, cha.name)
            if os.path.exists(fn):
                cha.inv_response = trace.InverseEvalresp(fn, tr)
            else:
                cha.inv_response = None