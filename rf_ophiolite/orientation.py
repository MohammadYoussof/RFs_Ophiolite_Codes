"""Orientation of the horizontal components from P-wave particle motion.

Commands:
    orientation        one measurement per station and event, and the median offset of
                       every station with its 95 % interval
    orientation-apply  decide which sensors to correct and write a corrected StationXML
"""

import math
import os
import sys

import numpy as np

from .common import (
    CHANNEL_BAND,
    output_dir,
    read_tables,
    pyplot,
    save_figure,
    save_csv,
    wrap180,
    read_channel,
)


def _orient_event(args):
    folder, rows, band, win, inv_path = args
    from obspy import read_inventory, UTCDateTime

    inv = read_inventory(inv_path) if inv_path else None
    res = []
    for r in rows:
        try:
            tr = dict((c, read_channel(folder, r['station'], c)) for c in 'ZNE')
            if any(v is None for v in tr.values()):
                continue
            meta = {}
            for c, t in tr.items():
                if inv is not None:
                    ch = inv.select(station=r['station'], channel=CHANNEL_BAND + c, time=t.stats.starttime)
                    ch = ch[0][0][0] if len(ch) and len(ch[0]) and len(ch[0][0]) else None
                    if ch is not None:
                        meta[c] = ch
                        sens = (
                            ch.response.instrument_sensitivity.value
                            if ch.response and ch.response.instrument_sensitivity
                            else 1.0
                        )
                        t.data = t.data / sens
                t.detrend('demean')
                t.detrend('linear')
                t.taper(0.05, max_length=30)
                t.filter('bandpass', freqmin=band[0], freqmax=band[1], corners=4, zerophase=True)
            t0 = UTCDateTime(r['t_ref'])
            seg = {}
            for c, t in tr.items():
                x = t.slice(t0 + win[0], t0 + win[1]).data
                seg[c] = np.asarray(x, float)
            m = min(len(v) for v in seg.values())
            if m < 10:
                continue
            Z, N, E = seg['Z'][:m], seg['N'][:m], seg['E'][:m]
            if 'Z' in meta and meta['Z'].dip is not None and meta['Z'].dip > 0:
                Z = -Z  # vertical positive down in metadata
            aN = math.radians(meta['N'].azimuth if 'N' in meta and meta['N'].azimuth is not None else 0.0)
            aE = math.radians(meta['E'].azimuth if 'E' in meta and meta['E'].azimuth is not None else 90.0)
            Ng = N * math.cos(aN) + E * math.cos(aE)  # geographic north/east from metadata azimuths
            Eg = N * math.sin(aN) + E * math.sin(aE)
            # horizontal direction that moves with the upward vertical = away from the source
            phi = math.degrees(math.atan2(np.sum(Z * Eg), np.sum(Z * Ng)))
            H = Ng * math.cos(math.radians(phi)) + Eg * math.sin(math.radians(phi))
            cc = float(np.sum(Z * H) / math.sqrt(np.sum(Z * Z) * np.sum(H * H) + 1e-30))
            baz_obs = (phi + 180.0) % 360.0
            res.append(
                dict(
                    event=r['event'],
                    station=r['station'],
                    baz=r['baz'],
                    baz_obs=round(baz_obs, 2),
                    delta_deg=round(float(wrap180(baz_obs - r['baz'])), 2),
                    cc_ZR=round(cc, 3),
                    snr_P_bb=r['snr'],
                )
            )
        except Exception as e:
            print('  skipped %s %s: %s' % (r['event'], r['station'], e))
    return res


def cmd_orientation(a):
    import pandas as pd
    from concurrent.futures import ProcessPoolExecutor

    plt = pyplot()
    out = output_dir(a)
    ev, sq = read_tables(a.qc)
    acc = set(ev[ev.status == 'accepted'].event)
    pk = pd.read_csv(os.path.join(a.qc, 'tables', 'picks_all.csv'))
    pk = pk[pk.phase == 'P'].set_index(['event', 'station']).pick_time.to_dict()
    d = sq[sq.event.isin(acc) & (sq.snr_P_bb >= a.min_snr) & sq.usable_P.fillna(False).astype(bool)]
    jobs = []
    for evn, g in d.groupby('event'):
        rows = [
            dict(
                event=evn,
                station=r.station,
                baz=r.baz,
                snr=round(r.snr_P_bb, 1),
                t_ref=pk.get((evn, r.station), r.tP_pred),
            )
            for r in g.itertuples()
        ]
        jobs.append((os.path.join(a.data, evn), rows, (a.fmin, a.fmax), (a.t1, a.t2), a.inventory))
    print('%d events, %d station-events with P SNR (0.05-1 Hz) >= %g' % (len(jobs), len(d), a.min_snr))
    res = []
    if a.workers <= 1:
        for j in jobs:
            res += _orient_event(j)
    else:
        with ProcessPoolExecutor(max_workers=a.workers) as ex:
            for r in ex.map(_orient_event, jobs):
                res += r
    m = pd.DataFrame(res)
    if not len(m):
        sys.exit('no orientation measurements')
    save_csv(m, out, 'orientation_measurements.csv')
    good = m[m.cc_ZR >= a.min_cc]
    rows = []
    rng = np.random.RandomState(0)
    for sta, g in good.groupby('station'):
        x = g.delta_deg.values
        med = float(np.median(x))
        bs = [np.median(rng.choice(x, len(x))) for _ in range(500)] if len(x) >= 3 else [np.nan]
        lo, hi = (np.nanpercentile(bs, 2.5), np.nanpercentile(bs, 97.5)) if len(x) >= 3 else (np.nan, np.nan)
        rows.append(
            dict(
                station=sta,
                n=len(x),
                median_delta_deg=round(med, 1),
                ci95_low=round(lo, 1),
                ci95_high=round(hi, 1),
                mad_deg=round(float(np.median(np.abs(x - med))), 1),
                north_axis_azimuth_deg=round(-med, 1),
                flag=bool(len(x) >= a.min_n and abs(med) > a.flag_deg and not (lo <= 0 <= hi)),
            )
        )
    s = pd.DataFrame(rows).sort_values('station')
    save_csv(s, out, 'orientation_by_station.csv')
    fig, ax = plt.subplots(figsize=(11, 4))
    xi = dict((st, i) for i, st in enumerate(s.station))
    ax.scatter(good.station.map(xi) + rng.uniform(-0.15, 0.15, len(good)), good.delta_deg, s=5, c='0.6')
    ax.errorbar(
        range(len(s)),
        s.median_delta_deg,
        yerr=[s.median_delta_deg - s.ci95_low, s.ci95_high - s.median_delta_deg],
        fmt='o',
        color='k',
        ms=4,
        capsize=2,
    )
    for i, r in enumerate(s.itertuples()):
        if r.flag:
            ax.annotate(
                'check',
                (i, r.median_delta_deg),
                textcoords='offset points',
                xytext=(0, 8),
                ha='center',
                fontsize=7,
            )
    ax.axhline(0, color='k', lw=0.6)
    ax.axhspan(-a.flag_deg, a.flag_deg, color='0.92', zorder=0)
    ax.set_xticks(range(len(s)))
    ax.set_xticklabels(s.station, rotation=90, fontsize=6)
    ax.set_ylim(-60, 60)
    ax.set_ylabel('observed minus expected back-azimuth (deg)')
    ax.set_title(
        'Horizontal orientation from P particle motion (gray: single events, black: median and 95%% interval, '
        'cc >= %.1f)' % a.min_cc
    )
    save_figure(fig, out, 'orientation.png')
    plt.close(fig)
    print(
        'stations flagged (|median| > %g deg, >= %d events): %s'
        % (a.flag_deg, a.min_n, ', '.join(s[s.flag].station) or 'none')
    )


# -----------------------------------------------------------------------------
# orientation-apply: decide which stations to correct and write a corrected StationXML
# -----------------------------------------------------------------------------
def cmd_orientation_apply(a):
    """A correction is applied only where the offset is a property of the sensor, not of the
    structure below it: enough measurements, a 95 % interval that excludes zero and is narrow,
    |median| >= min_deg, and the same offset (within 'consistency' degrees) in every
    back-azimuth bin that holds enough measurements. Direction-dependent offsets are
    reported but not corrected."""
    import pandas as pd
    from obspy import read_inventory

    out = output_dir(a)
    m = pd.read_csv(os.path.join(out, 'orientation_measurements.csv'))
    s = pd.read_csv(os.path.join(out, 'orientation_by_station.csv'))
    m = m[m.cc_ZR >= a.min_cc].copy()
    m['bz'] = (m.baz // a.bin_deg * a.bin_deg).astype(int)
    rows = []
    for r in s.itertuples():
        g = m[m.station == r.station]
        bins = g.groupby('bz').delta_deg.agg(['median', 'size'])
        bins = bins[bins['size'] >= a.bin_min]
        spread = (
            float(
                wrap180(bins['median'].values - r.median_delta_deg).max()
                - wrap180(bins['median'].values - r.median_delta_deg).min()
            )
            if len(bins) >= 2
            else np.nan
        )
        sig = np.isfinite(r.ci95_low) and not (r.ci95_low <= 0 <= r.ci95_high)
        half = (r.ci95_high - r.ci95_low) / 2.0 if np.isfinite(r.ci95_low) else np.nan
        if r.n < a.min_n:
            cls = 'insufficient'
        elif not sig or abs(r.median_delta_deg) < a.min_deg:
            cls = 'nominal'
        elif not np.isfinite(half) or half > a.max_ci or r.mad_deg > a.max_mad:
            cls = 'uncertain'
        elif len(bins) >= 2 and spread > a.consistency:
            cls = 'direction-dependent'
        else:
            cls = 'correct'
        rows.append(
            dict(
                station=r.station,
                n=r.n,
                median_delta_deg=r.median_delta_deg,
                ci95_low=r.ci95_low,
                ci95_high=r.ci95_high,
                mad_deg=r.mad_deg,
                n_bins=len(bins),
                spread_between_bins_deg=round(spread, 1) if np.isfinite(spread) else np.nan,
                decision=cls,
                north_axis_azimuth_deg=(
                    round(float(-r.median_delta_deg) % 360.0, 1) if cls == 'correct' else 0.0
                ),
            )
        )
    t = pd.DataFrame(rows)
    save_csv(t, out, 'orientation_corrections.csv')
    print(t.decision.value_counts().to_string())
    print(
        'corrected:',
        ', '.join(
            '%s %+.0f' % (r.station, r.median_delta_deg) for r in t[t.decision == 'correct'].itertuples()
        ),
    )
    print('direction-dependent (not corrected):', ', '.join(t[t.decision == 'direction-dependent'].station))
    print('uncertain (not corrected):', ', '.join(t[t.decision == 'uncertain'].station))
    inv = read_inventory(a.inventory)
    corr = dict((r.station, r.north_axis_azimuth_deg) for r in t[t.decision == 'correct'].itertuples())
    n = 0
    for net in inv:
        for sta in net:
            if sta.code in corr:
                for ch in sta:
                    if ch.code.endswith('N'):
                        ch.azimuth = float(corr[sta.code])
                        n += 1
                    elif ch.code.endswith('E'):
                        ch.azimuth = float((corr[sta.code] + 90.0) % 360.0)
                        n += 1
    out_xml = a.out_xml or os.path.splitext(a.inventory)[0] + '_oriented.xml'
    inv.write(out_xml, format='STATIONXML')
    print('wrote %s with %d horizontal channels re-oriented' % (out_xml, n))
