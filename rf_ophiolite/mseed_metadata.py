"""Metadata that a miniSEED copy of the event data needs.

SAC files carry the event and station information in their headers; miniSEED files do
not. This command reads the SAC headers, the miniSEED headers and the StationXML (all
read only) and writes, to the folder --out:

    events_catalog.csv, events_catalog.xml   one line (QuakeML: one event) per event folder:
                                             origin time, location, depth, magnitude, and
                                             the time window of the files
    station_channels.csv                     coordinates, orientation, sample rate, sensor,
                                             datalogger and sensitivity of every channel
    mseed_index.csv                          one line per miniSEED file: stream, start and
                                             end time, samples, gaps, encoding, size

    python rfo.py mseed-metadata --sac data/raw_sac --mseed data/raw_mseed \
                                 --inventory stations.xml --out data/mseed_metadata
"""

import csv
import glob
import os
import sys
import time


from .common import is_inside


def add_arguments(p):
    p.add_argument('--sac', required=True, help='raw SAC event folders (read only)')
    p.add_argument('--mseed', required=True, help='raw miniSEED event folders (read only)')
    p.add_argument('--inventory', default=None, help='StationXML of the network (for station_channels.csv)')
    p.add_argument('--out', required=True, help='output folder, outside the two data folders')
    p.add_argument('--model', default='ak135', help='travel-time model for the window relative to P')


def _write_csv(path, rows):
    with open(path, 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    print('  wrote', path)


def event_catalog(sac_dir, mseed_dir, out, model_name):
    from obspy import read, UTCDateTime
    from obspy.core.event import Catalog, Comment, Event, Magnitude, Origin
    from obspy.taup import TauPyModel

    model = TauPyModel(model_name)
    rows = []
    for d in sorted(glob.glob(os.path.join(sac_dir, 'Event_*'))):
        files = sorted(glob.glob(os.path.join(d, '*.sac')))
        if not files:
            continue
        h = read(files[0], headonly=True)[0].stats.sac
        origin = UTCDateTime(
            year=int(h.nzyear),
            julday=int(h.nzjday),
            hour=int(h.nzhour),
            minute=int(h.nzmin),
            second=int(h.nzsec),
            microsecond=int(h.nzmsec) * 1000,
        ) + float(h.o)
        depth_km = float(h.evdp) / 1000.0 if float(h.evdp) > 1000 else float(h.evdp)  # SOD writes meters
        # SOD cuts every station the same way relative to its predicted P, so the window of the
        # first station describes the event
        row = dict(
            event_id=os.path.basename(d),
            origin_time_utc=origin.isoformat(),
            latitude=round(float(h.evla), 4),
            longitude=round(float(h.evlo), 4),
            depth_km=round(depth_km, 1),
            magnitude=round(float(h.mag), 1),
            n_stations=len(set(os.path.basename(f).split('.')[1] for f in files)),
            n_sac_files=len(files),
            n_mseed_files=len(glob.glob(os.path.join(mseed_dir, os.path.basename(d), '*.mseed'))),
            window_start_s_after_origin=round(float(h.b), 1),
            window_end_s_after_origin=round(float(h.e), 1),
            window_start_rel_P_s='',
            window_end_rel_P_s='',
        )
        try:
            tp = model.get_travel_times(
                source_depth_in_km=max(depth_km, 0.0),
                distance_in_degree=float(h.gcarc),
                phase_list=['P', 'Pdiff'],
            )[0].time
            row['window_start_rel_P_s'] = round(float(h.b) - tp, 1)
            row['window_end_rel_P_s'] = round(float(h.e) - tp, 1)
        except Exception:
            pass
        rows.append(row)
    if not rows:
        sys.exit('no Event_* folders with SAC files in %s' % sac_dir)
    _write_csv(os.path.join(out, 'events_catalog.csv'), rows)
    cat = Catalog()
    for r in rows:
        o = Origin(
            time=UTCDateTime(r['origin_time_utc']),
            latitude=r['latitude'],
            longitude=r['longitude'],
            depth=r['depth_km'] * 1000.0,
        )
        m = Magnitude(mag=r['magnitude'])
        ev = Event(
            origins=[o], magnitudes=[m], comments=[Comment(text='SOD event folder %s' % r['event_id'])]
        )
        ev.preferred_origin_id = o.resource_id
        ev.preferred_magnitude_id = m.resource_id
        cat.append(ev)
    path = os.path.join(out, 'events_catalog.xml')
    cat.write(path, format='QUAKEML')
    print('  wrote', path)
    return len(rows)


def station_channels(inventory, out):
    from obspy import read_inventory

    rows = []
    for net in read_inventory(inventory):
        for sta in net:
            for ch in sta:
                sens = ch.response.instrument_sensitivity if ch.response else None
                rows.append(
                    dict(
                        network=net.code,
                        station=sta.code,
                        location=ch.location_code,
                        channel=ch.code,
                        latitude=ch.latitude,
                        longitude=ch.longitude,
                        elevation_m=ch.elevation,
                        depth_m=ch.depth,
                        azimuth_deg=ch.azimuth,
                        dip_deg=ch.dip,
                        sample_rate_hz=ch.sample_rate,
                        start_date=str(ch.start_date)[:19] if ch.start_date else '',
                        end_date=str(ch.end_date)[:19] if ch.end_date else '',
                        sensor=(ch.sensor.description or ch.sensor.type) if ch.sensor else '',
                        datalogger=(
                            (ch.data_logger.description or ch.data_logger.type) if ch.data_logger else ''
                        ),
                        sensitivity=sens.value if sens else '',
                        sensitivity_frequency_hz=sens.frequency if sens else '',
                        input_units=sens.input_units if sens else '',
                    )
                )
    _write_csv(os.path.join(out, 'station_channels.csv'), rows)
    return len(rows)


def mseed_index(mseed_dir, out):
    """One line per miniSEED file. The index is appended event by event, so an interrupted run
    continues where it stopped."""
    from obspy import read

    path = os.path.join(out, 'mseed_index.csv')
    fields = [
        'event_id',
        'file',
        'network',
        'station',
        'location',
        'channel',
        'start_time_utc',
        'end_time_utc',
        'sample_rate_hz',
        'npts',
        'n_traces',
        'gaps',
        'encoding',
        'bytes',
    ]
    done = set()
    if os.path.exists(path):
        with open(path) as fh:
            done = set(r['event_id'] for r in csv.DictReader(fh))
    folders = [
        d for d in sorted(glob.glob(os.path.join(mseed_dir, 'Event_*'))) if os.path.basename(d) not in done
    ]
    t0 = time.time()
    with open(path, 'a', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        if not done:
            w.writeheader()
        for k, d in enumerate(folders, 1):
            for f in sorted(glob.glob(os.path.join(d, '*.mseed'))):
                row = dict.fromkeys(fields, '')
                row.update(event_id=os.path.basename(d), file=os.path.basename(f), bytes=os.path.getsize(f))
                try:
                    st = read(f, headonly=True)
                    tr = st[0]
                    row.update(
                        network=tr.stats.network,
                        station=tr.stats.station,
                        location=tr.stats.location,
                        channel=tr.stats.channel,
                        start_time_utc=str(min(t.stats.starttime for t in st))[:23],
                        end_time_utc=str(max(t.stats.endtime for t in st))[:23],
                        sample_rate_hz=tr.stats.sampling_rate,
                        npts=sum(t.stats.npts for t in st),
                        n_traces=len(st),
                        gaps=len(st.get_gaps()),
                        encoding=tr.stats.mseed.encoding,
                    )
                except Exception as e:
                    row['gaps'] = 'read error: %s' % str(e)[:60]
                w.writerow(row)
            fh.flush()
            if k % 50 == 0:
                print('  %d / %d event folders indexed (%.0f s)' % (k, len(folders), time.time() - t0))
    print('  wrote', path, '(%d event folders added, %d were already there)' % (len(folders), len(done)))


def run(a):
    for d in (a.sac, a.mseed):
        if not os.path.isdir(d):
            sys.exit('ERROR: folder not found: %s' % d)
        if is_inside(a.out, d):
            sys.exit('ERROR: --out (%s) must be outside the data folder %s' % (a.out, d))
    if not os.path.isdir(a.out):
        os.makedirs(a.out)
    n = event_catalog(a.sac, a.mseed, a.out, a.model)
    print('%d events' % n)
    if a.inventory:
        print('%d channels' % station_channels(a.inventory, a.out))
    mseed_index(a.mseed, a.out)
