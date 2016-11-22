import sys, tempfile, os, glob, re
from pyrocko import iris_ws, util, catalog, model, weeding, pile, io, orthodrome, trace, gui_util, evalresp_ext
import logging
import cPickle as pickle
import numpy as num
from subprocess import call

logger = logging.getLogger('iris_pull')

__autokiwi_usage__ = {
        'catalog_search':       'autokiwi catalog_search',
        'iris_pull':            'autokiwi iris_pull ( all | eventnames ... )',
        'manual_qc':            'autokiwi manual_qc ( all | eventnames ... )',
        'kiwi_setup':           'autokiwi kiwi_setup ( all | eventnames ... )',
        'plot_responses':       'autokiwi plot_responses ( all | eventnames ... )',
        'xlist':                'autokiwi xlist',
        'go':                   'autokiwi go', 
}

__autokiwi_commands__ = __autokiwi_usage__.keys()

def plural_s(i):
    if i == 1:
        return ''
    else:
        return 's'

def make_event_name(event):
    return util.time_to_str(event.time, format='%Y-%m-%d_%H-%M-%S')


def dump_event_infos(conf, event_name, event):
    p = conf.path('event_info_path', additional=dict(event_name=event_name) )
    util.ensuredirs(p)
    return event.olddump(p)

def __autokiwi_main__(command, options, conf, event_names):
    dispatch = dict(catalog_search=catalog_search, iris_pull=iris_pull, 
            manual_qc=manual_qc, plot_responses=plot_responses,
            kiwi_setup=kiwi_setup, xlist=xlist, go=go)

    dispatch[command](options, conf, event_names)

def event_names_from_pathes(base_config):
    # not nice; reverse lookup event names from pathes
    ev_dirs = glob.glob(base_config.path('event_dir', additional={ 'event_name': '*' }))
    event_names = []
    for ev_dir in ev_dirs:
        if not os.path.isdir(ev_dir): continue
        pat = base_config.path('event_dir') % {'event_name': 'xxxEVENTNAMExxx'}
        pat = re.escape(pat).replace('xxxEVENTNAMExxx', r'([^/]+)')
        m = re.match(pat, ev_dir)
        if m:
            ev_name = m.group(1)
            event_names.append(ev_name)

    return event_names

def inventory(options, config, event_names):
    base_conf = config['base_config']
    pull_conf = config['iris_pull_config']
    kiwi_conf = config['kiwi_config']

    event_names = sorted(event_names_from_pathes(base_conf))

    out = []

    for event_name in event_names:
        a = dict(event_name=event_name)
        event_dir = base_conf.path('event_dir', additional=a)
        data_dir = pull_conf.path('data_dir', additional=a)
        kiwi_dir = kiwi_conf.path('main_dir', additional=a)
        work_dir = kiwi_conf.path('work_dir', additional=a)
        report_dir = kiwi_conf.path('report_dir', additional=a)
        
        i = 0
        for d in (event_dir, data_dir, kiwi_dir, work_dir, report_dir):
            if not os.path.isdir(d):
                break

            if i == 1 and not os.path.exists(os.path.join(data_dir, 'stations.txt')):
                break
            
            i += 1

        fail = False
        fail_path = base_conf.path_or_none('fail_filename', additional=a)
        if fail_path and os.path.exists(fail_path):
            fail = True

        out.append((event_name, i, fail))

    return out

def xlist(options, config, event_names):
    for event_name, i, fail in inventory(options, config, event_names):
        print '%-6s %-4s %s' % ('#'*i, ['','fail'][fail], event_name)

def go(options, config, event_names):
    for event_name, i, fail in reversed(inventory(options, config, event_names)):
        if fail:
            logger.warn('Skipping event %s because of previous failure.' % event_name)
            continue

        if i < 5:
            commands = ','.join('iris_pull kiwi_setup process report post'.split()[i-1:])
            command_line = [ 'autokiwi', commands, event_name ]
            print ' '.join(command_line)
            call(command_line)
            break


def catalog_search(options, conf, event_names):
    conf = conf['catalog_search_config'] 
    
    gcmt = catalog. Geofon()
    for event in gcmt.iter_events(
            time_range=conf.timerange('time_range'),
            magmin=conf.minimum_magnitude):

        if conf.has('maximum_magnitude'):
            if conf.maximum_magnitude < event.magnitude:
                continue

        if conf.has('event_filter'):
            if not conf.event_filter(event):
                continue

        dump_event_infos(conf, make_event_name(event), event)


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


def combi_get_stations(*args, **kwargs):
    sites = ['geofon', 'iris']
    stations = {}
    for site in sites:
        stations_this = get_stations(site, *args, **kwargs)
        for station in stations_this:
            nsl = station.nsl()
            if nsl not in stations:
                stations[nsl] = station
                station.datacenters = [site]
            else:
                stations[nsl].datacenters.append(site)

    stations_list = [stations[nsl] for nsl in sorted(stations.keys())]
    return stations_list

def iris_get_vnets(codes, tmin, tmax):
    vnet = []
    for code in codes:
        data = iris_ws.ws_virtualnetwork(code=code, output='XML', timewindow=(tmin, tmax))
        vnet.extend(iris_ws.grok_virtualnet_xml(data))
    return vnet


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


def combi_get_data(stations, tmin, tmax, fn_template, neach=20):
    from pyrocko.fdsn import ws

    have = set()
    def remember(trs):
        for tr in trs:
            have.add(tr.nslc_id)
            yield tr

    fns = []
    for site in ['geofon', 'iris']:
        stations_site = [station for station in stations if site in station.datacenters]
        i = 0
        while i < len(stations_site):
            stations_now = stations[i:i+neach]
            selection = ws.make_data_selection(stations_now, tmin, tmax)
            selection = [sel for sel in selection if sel[:4] not in have]
            f = tempfile.NamedTemporaryFile()
            try:
                data = ws.dataselect(site=site, selection=selection)
                while True:
                    buf = data.read(1024)
                    if not buf:
                        break
                    f.write(buf)

                f.flush()

                fns.extend(io.save(remember(io.iload(f.name)), fn_template))

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


def combi_get_responses(stations, time, fn_template):
    from pyrocko.fdsn import ws
    from pyrocko.fdsn import station as fdsnstation

    def fn(net, sta, loc, cha):
        return fn_template % dict(
            network=net, station=sta, location=loc, channel=cha)

    def iter_nslcs(site=None, ignore=None):
        for station in stations:
            if site is not None and site not in station.datacenters:
                continue

            for channel in station.get_channels():
                nslc = station.nsl() + (channel.name,)
                if ignore is None or nslc not in ignore:
                    yield nslc

    responses = {}
    for nslc in iter_nslcs():
        if os.path.exists(fn(*nslc)):
            responses[nslc] = pload(fn(*nslc))

    for site in ['geofon', 'iris']:
        selection = []
        for nslc in iter_nslcs(site=site, ignore=responses):
            selection.append(nslc + (time, time+1.0))

        if selection:
            logger.info('downloading response information (%s)' % site)
            sxs= ws.station(
                site=site, level='response', selection=selection)
            
            for nslc_tspan in selection:
                nslc = nslc_tspan[:4]
                timespan = nslc_tspan[4:]
                try:
                    response = sxs.get_pyrocko_response(
                        nslc, 
                        timespan=timespan,
                        fake_input_units='M')

                    util.ensuredirs(fn(*nslc))
                    pdump(response, fn(*nslc))
                    responses[nslc] = response

                except (fdsnstation.NoResponseInformation,
                        fdsnstation.MultipleResponseInformation):
                    pass

    for station in stations:
        for channel in station.get_channels():
            nslc = station.nsl() + (channel.name,)
            if nslc in responses:
                channel.response = responses[nslc]
            else:
                channel.response = None


def select_stations(all_stations, nwanted, preferred_n=[], preferred_ns=[], preferred_nsl=[], blacklist_nsl=[]):
    badnesses = {} 
    stations = []
    for sta in all_stations:
        n = sta.network
        ns = sta.network, sta.station
        nsl = sta.network, sta.station, sta.location
        if nsl in blacklist_nsl: 
            continue

        stations.append(sta)

        badnesses[nsl] = 4.0
        if n in preferred_n:
            badnesses[nsl] = 1.0

        if ns in preferred_ns:
            badnesses[nsl] = 2.0

    for nsl in preferred_nsl:
        badnesses[nsl] =  1.0

    stations, _, _ = weeding.weed_stations(stations, nwanted=nwanted, badnesses_nsl=badnesses)
    
    stations_ns = {}
    for sta in stations:
        ns = sta.network, sta.station
        stations_ns.setdefault(ns,[]).append(sta)

    for k in stations_ns.keys():
        stations_ns[k].sort(key=lambda sta: sta.location)
        stations_ns[k] = stations_ns[k][0]

    return stations_ns.values()

def get_nsl(x):
    return x.network, x.station, x.location

class Preparator(object):

    def __init__(self, begin_phase, end_phase):
        self._begin_phase = begin_phase
        self._end_phase = end_phase

    def distance(self, event, station):
        return orthodrome.distance_accurate50m(event, station)

    def get_tpad(self, event, station):
        return 0.0

    def get_tmin(self, event, station):
        return event.time + self._begin_phase(self.distance(event,station), event.depth)

    def get_tmax(self, event, station):
        return event.time + self._end_phase(self.distance(event,station), event.depth)

    def get_tmin_limit(self, event, stations):
        return min(self.get_tmin(event, s) - self.get_tpad(event, s) for s in stations)
    
    def get_tmax_limit(self, event, stations):
        return max(self.get_tmax(event, s) + self.get_tpad(event, s) for s in stations)

    def prepare_station_traces(self, traces, event, station):
        return traces

    def iter_prepare(self, pile, event, stations): 
        stations_by_nsl = dict((get_nsl(s), s) for s in stations)
        nsl_avail = pile.gather_keys(gather=get_nsl)
        for nsl in nsl_avail:
            station = stations_by_nsl[nsl]
            tmin = self.get_tmin(event, station)
            tmax = self.get_tmax(event, station)
            tpad = self.get_tpad(event, station)
            
            for traces in pile.chopper(tmin=tmin, tmax=tmax, tpad=tpad,
                    trace_selector=lambda tr: get_nsl(tr) == nsl,
                    want_incomplete=False):

                if not traces:
                    continue

                station = stations_by_nsl[nsl]
                prepared_traces = self.prepare_station_traces(traces, event, station)

                yield prepared_traces

class InvResponsePreparator(Preparator):

    def __init__(self, begin_phase, end_phase, frequencyband, tfade_factor=2.0):
        Preparator.__init__(self, begin_phase, end_phase)
        self._frequencyband = frequencyband
        self._tfade_factor = tfade_factor

    def get_tpad(self, event, station):
        return self.get_tfade()
    
    def get_tfade(self):
        return self._tfade_factor/self._frequencyband[1]

    def prepare_station_traces(self, traces, event, station):
        out_traces = []
        for tr in traces:
            tr = tr.copy()
            cha = station.get_channel(tr.channel)
            try:
                if not cha.response:
                    logger.warn('No response for channel %s.%s.%s.%s' % tr.nslc_id)
                    raise

                tr2 = tr.transfer(
                    tfade=self.get_tfade(), 
                    freqlimits=self._frequencyband, 
                    transfer_function=cha.response,
                    invert=True)

                if not num.all(num.isfinite(tr2.get_ydata())):
                    logger.warn( 'Trace %s.%s.%s.%s has NaNs or Infs' % tr2.nslc_id )
                    raise 

                out_traces.append(tr2)
            except:
                logger.warn('restitution failed for trace %s.%s.%s.%s' % tr.nslc_id )

        
        return out_traces

def str_nsl_selection(nsl):
    by_n = {}
    for (n,s,l) in nsl:
        by_n.setdefault(n, []).append((s,l))

    snet = []
    for n in sorted(by_n.keys()):
        snet.append('  %s: ' % n + ', '.join( (sl[0], '%s.%s' % sl)[bool(sl[1])] for sl in sorted(by_n[n]) ))

    return '\n'.join(snet)

def pload(fn):
    f = open(fn, 'r')
    data = pickle.load(f)
    f.close()
    return data

def pdump(data, fn):
    f = open(fn,'w')
    pickle.dump(data, f)
    f.close()

def iris_pull(options, conf, event_names):
    conf = conf['iris_pull_config']

    if not event_names:
        sys.exit('need event name')

    preparator = InvResponsePreparator(conf.begin_phase, conf.end_phase, conf.inv_response_frequencyband)

    for event_name in event_names:
        conf.event_name = event_name
        event = _get_event_infos(conf)
        tevent = event.time
        
        station_query_save_path = conf.path('station_query_save_path')
        if os.path.exists(station_query_save_path):
            logger.info('Using stored station query.')
            all_stations = pload(station_query_save_path)
        else:
            logger.info('Querying for stations...')
            all_stations = combi_get_stations(lat=event.lat, lon=event.lon, 
                    rmin=conf.query_rmin, rmax=conf.query_rmax,
                    tmin=tevent, tmax=tevent+3600., channel_pattern=conf.query_channel_pattern)

            util.ensuredirs(station_query_save_path)
            pdump(all_stations, station_query_save_path)

        nstations = len(set( (sta.network, sta.station) for sta in all_stations ))
        logger.info('Station query returned %i station%s' % (nstations, plural_s(nstations)))
        preferred_ns = set(iris_get_vnets(conf.preferred_virtual_networks, tmin=tevent, tmax=tevent+3600.))
        preferred_n = set(conf.preferred_networks)

        for station in all_stations:
            station.set_event_relative_data(event)

        raw_trace_path = conf.path('raw_trace_path')

        nsl_all = set( get_nsl(s) for s in all_stations)
        
        state_save_path = conf.path('state_save_path')
        if os.path.exists(state_save_path):
            nsl_ok, nsl_blacklist, nsl_use = pload(state_save_path)
        else:
            nsl_ok = set()
            nsl_blacklist = set() 
            nsl_use = set()

        manual_blacklist_path = conf.path('manual_blacklist_path')
        nsl_blacklist.update(read_manual_blacklist(manual_blacklist_path))
        
        nsl_selected = set()
        nwanted = conf.get_or_none('nwanted')
        assert len(nwanted) == 2

        while True:
            if nwanted:
                selected = select_stations(all_stations, nwanted[1], preferred_n=preferred_n, preferred_ns=preferred_ns, preferred_nsl=nsl_ok, blacklist_nsl=nsl_blacklist)
            else:
                selected = all_stations

            nsl_selected = set( get_nsl(s) for s in selected )

            download = [ s for s in selected if get_nsl(s) not in nsl_ok ]
            nsl_download = set( get_nsl(s) for s in download )
             
            combi_get_responses(download, tevent, conf.path('resp_path'))

            tmin = preparator.get_tmin_limit(event, download)
            tmax = preparator.get_tmax_limit(event, download)

            logger.info('Loading data for event %s:\n%s' % (event_name, str_nsl_selection(nsl_download)))
            try: 
                fns = combi_get_data(download, tmin, tmax, raw_trace_path, neach=conf.query_nstations_per_datarequest)
                p = pile.make_pile(fns, show_progress=False)
                prepared_trace_path = conf.path('prepared_trace_path')
                for traces in preparator.iter_prepare(p, event, download):
                    for tr in traces:
                        nsl_ok.add(get_nsl(tr))
                    io.save(traces, prepared_trace_path)

            except iris_ws.NotFound:
                pass

            logger.info('Blacklisting:\n%s' % str_nsl_selection(nsl_download - nsl_ok))
            
            nsl_blacklist.update(nsl_download - nsl_ok)
            preferred_ns.update( set( nsl[:2] for nsl in nsl_ok ) )

            nsl_use = nsl_ok & nsl_selected

            logger.info('Have %i stations with data.' % len(nsl_use))

            if not nwanted:
                break

            else:
                if len(nsl_selected) == len(nsl_all):
                    break

                if len(nsl_use) > nwanted[0]:
                    break

        pdump( (nsl_ok, nsl_blacklist, nsl_use), state_save_path )

        if nwanted and len(nsl_use) >= nwanted[1]:
            selected = select_stations(selected, nwanted[1], preferred_n=preferred_n, preferred_ns=preferred_ns, blacklist_nsl=nsl_blacklist)
            nsl_selected = set( get_nsl(s) for s in selected )
            nsl_use = nsl_ok & nsl_selected

        stations = [ s for s in all_stations if get_nsl(s) in nsl_use ]
        
        model.dump_stations(stations, conf.path('stations_ok_path'))
        model.dump_stations(all_stations, conf.path('stations_all_path'))
            
        if nsl_use:
            logger.info('Data available for event %s:\n%s' % (event_name, str_nsl_selection(nsl_use)))
        else:
            logger.info('No data availabe for event %s' % event_name) 


class Stars:
    def __getitem__(self, k):
        return '*'

def read_manual_blacklist(fn):
    nsl_blacklist = set()
    if not os.path.exists(fn):
        return set()

    f = open(fn, 'r')
    for line in f:
        nsl = line.strip().split('.')
        assert len(nsl) == 3
        nsl_blacklist.add(tuple(nsl))

    f.close()
    return nsl_blacklist

def write_manual_blacklist(nsl_blacklist, fn):
    f= open(fn, 'w')
    for nsl in sorted(list(nsl_blacklist)):
        f.write('.'.join(nsl))
        f.write('\n')

    f.close()

def update_manual_blacklist(nsl_blacklist, fn):
    b = set(nsl_blacklist)
    b.update(read_manual_blacklist(fn))    
    write_manual_blacklist(b, fn)

def _get_event_infos(conf):
    return model.Event(load=conf.path('event_info_path', additional=dict(event_name=conf.event_name) ))

def _get_stations(conf):
    stations = model.load_stations(conf.path('stations_ok_path'))
    return stations

def _get_pile(conf, pathconf):
    fn_pattern = conf.path(pathconf)
    fns = glob.glob(fn_pattern % Stars())
    return pile.make_pile(fns, show_progress=False)

def _get_traces(conf, stations, pathconf):
    p = _get_pile(conf, pathconf)

    fn = conf.path('manual_blacklist_path')
    nsl_blacklist = read_manual_blacklist(fn)

    stations_by_nsl = dict((get_nsl(s), s) for s in stations)
    traces = p.all(trace_selector=lambda tr: get_nsl(tr) in stations_by_nsl and get_nsl(tr) not in nsl_blacklist)
    return traces

def _get_raw_traces(conf, stations):
    return _get_traces(conf, stations, 'raw_trace_path')

def _get_prepared_traces(conf, stations):
    return _get_traces(conf, stations, 'prepared_trace_path')

def manual_qc(options, conf, event_names):
    conf = conf['iris_pull_config']

    if not event_names:
        sys.exit('need event name')

    for event_name in event_names:
        conf.event_name = event_name
        event = _get_event_infos(conf)
        stations = _get_stations(conf)
        traces = _get_prepared_traces(conf, stations)
        
        fn = conf.path('manual_blacklist_path')
        nsl_blacklist = read_manual_blacklist(fn)

        retval, markers = trace.snuffle(traces, events=[event], stations=stations, want_markers=True)
        
        for m in markers:
            if type(m) is gui_util.Marker:
                try:
                    nslc = m.one_nslc()
                    nsl_blacklist.add(nslc[:3])

                except gui_util.MarkerOneNSLCRequired:
                    pass
    

        update_manual_blacklist(nsl_blacklist, fn)

class DummyAcc:
    def __init__(self, raw_traces):
        self.raw_traces = raw_traces

    def iter_traces(self, trace_selector):
        yield  [ tr for tr in self.raw_traces if trace_selector is None or trace_selector(tr) ] 

def kiwi_setup(options, config, event_names):
    from tunguska import prepare, gfdb

    conf = config['iris_pull_config']
    kiwi_conf = config['kiwi_config']

    if kiwi_conf.has('gfdb_path'):
        db = gfdb.Gfdb(kiwi_conf.path('gfdb_path'))
        deltat = db.dt
    else:
        if kiwi_conf.has('deltat'):
            deltat = kiwi_conf.deltat
        else:
            deltat = None
        db = None
    
    if not event_names:
        sys.exit('need event name')

    for event_name in event_names:
        conf.event_name = event_name
        kiwi_conf.event_name = event_name

        event = _get_event_infos(conf)
        stations = _get_stations(conf)
        for station in stations:
            station.set_event_relative_data(event)

        traces = _get_prepared_traces(conf, stations)
        raw_traces = _get_raw_traces(conf, stations)
         
        p = pile.Pile()
        buf = pile.MemTracesFile(None,traces)
        p.add_file(buf)

        processed = []
        for station in stations:
            
            tt1 = kiwi_conf.cut_span[0](station.dist_m, event.depth)
            tt2 = kiwi_conf.cut_span[1](station.dist_m, event.depth)
            if None in (tt1, tt2):
                continue

            tmin = tt1 + event.time
            tmax = tt2 + event.time

            traces = p.all( tmin=tmin, tmax=tmax, want_incomplete=False, 
                    trace_selector=lambda tr: get_nsl(tr) == get_nsl(station))

            for proj, in_channels, out_channels in station.guess_projections_to_rtu(out_channels=('R', 'T', 'Z')):
                proc =  trace.project(traces, proj, in_channels, out_channels) 
                processed.extend(proc)
                for tr in proc:
                    for ch in out_channels:
                        if ch.name == tr.channel:
                            station.add_channel(ch)

        for tr in processed:
            if deltat is not None:
                try:
                    tr.downsample_to(deltat, snap=True, allow_upsample_max=5)
                except util.UnavailableDecimation, e:
                    logger.warn( 'Cannot downsample %s.%s.%s.%s: %s' % (tr.nslc_id + (e,)))
                    continue
            
        stations_by_nsl = dict((get_nsl(s), s) for s in stations)

        acc = DummyAcc(raw_traces)
        prepare.save_kiwi_dataset(acc, stations_by_nsl, processed, event, kiwi_conf)


def plot_responses(options, config, event_names):
    import pylab as lab
    
    conf = config['iris_pull_config']
    if not event_names:
        sys.exit('need event name')

    for event_name in event_names:
        conf.event_name = event_name
        stations = _get_stations(conf)
        event = _get_event_infos(conf)
        fband = conf.inv_response_frequencyband
        

        combi_get_responses(stations, event.time, conf.path('resp_path'))

        for station in stations:

            for cha in station.get_channels():
                resp = cha.inv_response
                
                fmin, fmax = fband[0], fband[3]
                freqs = num.exp(num.linspace(num.log(fmin), num.log(fmax), 400))
                amps = num.abs(resp.evaluate(freqs))

                lab.loglog(freqs, amps)


    lab.show()
