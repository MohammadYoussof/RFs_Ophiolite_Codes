"""Overview of the data set after the quality control.

Commands:
    inventory    event coverage (distance, back-azimuth), event selection, station
                 availability and integrity problems per station
    snr          P and S signal-to-noise ratio against magnitude and distance
    response     amplitude and phase response of each sensor type
    noise-hour   pre-event noise (1 to 5 Hz) against local time of day
"""

import glob
import math
import os
import sys

import numpy as np

from .common import PSD_FMIN, PSD_FMAX, output_dir, read_tables, pyplot, save_figure, save_csv


def cmd_inventory(a):
    import pandas as pd

    plt = pyplot()
    out = output_dir(a)
    ev, sq = read_tables(a.qc)
    geo = sq.groupby('event').agg(dist=('dist', 'median'), baz=('baz', 'median')).reset_index()
    ev = ev.merge(geo, on='event', how='left')
    # --- coverage: back-azimuth (angle) and distance (radius) ---
    fig = plt.figure(figsize=(15.5, 4.8))
    ax = fig.add_subplot(1, 3, 1, projection='polar')
    ax.set_theta_zero_location('N')
    ax.set_theta_direction(-1)
    for status, col in (('accepted', 'k'), ('rejected', '0.6'), ('error', 'r')):
        e = ev[ev.status == status]
        if len(e):
            ax.scatter(
                np.radians(e.baz),
                e.dist,
                s=4 + (e.mag - 5.0).clip(0) * 25,
                c=col,
                alpha=0.7,
                label='%s (%d)' % (status, len(e)),
                edgecolor='none',
            )
    ax.set_rlim(0, 100)
    ax.set_title(
        '(a) All events: back-azimuth (angle), distance (radius, deg); size = magnitude', pad=14, fontsize=9
    )
    ax.legend(loc='lower left', bbox_to_anchor=(-0.25, -0.12), fontsize=8, frameon=False)
    # (b) accepted events by the number of class A-C P picks; the high-quality set highlighted
    axb = fig.add_subplot(1, 3, 2, projection='polar')
    axb.set_theta_zero_location('N')
    axb.set_theta_direction(-1)
    acc = ev[ev.status == 'accepted'].copy()
    npk = acc.n_P_picks_ABC.fillna(0) if 'n_P_picks_ABC' in acc else pd.Series(0, index=acc.index)
    hq = npk >= a.hq_picks
    axb.scatter(
        np.radians(acc.baz[~hq]),
        acc.dist[~hq],
        s=3 + npk[~hq],
        c='0.65',
        alpha=0.7,
        edgecolor='none',
        label='fewer than %d picks (%d)' % (a.hq_picks, int((~hq).sum())),
    )
    axb.scatter(
        np.radians(acc.baz[hq]),
        acc.dist[hq],
        s=3 + npk[hq],
        c='k',
        alpha=0.8,
        edgecolor='none',
        label='%d or more class A-C P picks (%d)' % (a.hq_picks, int(hq.sum())),
    )
    axb.set_rlim(0, 100)
    axb.set_title('(b) Accepted events; size = number of P picks', pad=14, fontsize=9)
    axb.legend(loc='lower left', bbox_to_anchor=(-0.25, -0.12), fontsize=8, frameon=False)
    ax2 = fig.add_subplot(1, 3, 3)
    bins = np.arange(5.0, 9.0, 0.25)
    for status, col in (('accepted', 'k'), ('rejected', '0.6')):
        e = ev[ev.status == status]
        ax2.hist(e.mag.dropna(), bins=bins, color=col, alpha=0.8, label=status)
    ax2.set_xlabel('magnitude')
    ax2.set_ylabel('events')
    ax2.legend(frameon=False)
    ax2.set_title('(c) Magnitude of accepted and rejected events', fontsize=9)
    plt.tight_layout()
    save_figure(fig, out, 'event_coverage.png')
    plt.close(fig)
    # --- funnel ---
    rows = []
    if a.listed:
        rows.append(dict(stage='events listed by the download', events=a.listed))
    rows.append(dict(stage='event folders with data', events=len(ev)))
    for r_, n_ in (
        ev[ev.status == 'rejected']
        .reason.fillna('')
        .astype(str)
        .str.split('(')
        .str[0]
        .str.strip()
        .value_counts()
        .items()
    ):
        rows.append(dict(stage='rejected: ' + r_, events=int(n_)))
    if (ev.status == 'error').any():
        rows.append(dict(stage='processing error', events=int((ev.status == 'error').sum())))
    acc = ev[ev.status == 'accepted']
    rows.append(dict(stage='accepted: 4 or more class A-C P picks', events=int(len(acc))))
    if 'n_P_picks_A' in ev:
        rows.append(
            dict(
                stage='accepted with 10 or more class A P picks',
                events=int((acc.n_P_picks_A.fillna(0) >= 10).sum()),
            )
        )
    if 'n_S_picks_ABC' in ev:
        rows.append(
            dict(stage='accepted with S picks (A-C)', events=int((acc.n_S_picks_ABC.fillna(0) > 0).sum()))
        )
    fn = pd.DataFrame(rows)
    save_csv(fn, out, 'event_funnel.csv')
    fig, ax = plt.subplots(figsize=(8, 0.45 * len(fn) + 1))
    ax.barh(range(len(fn)), fn.events, color='0.45')
    ax.set_yticks(range(len(fn)))
    ax.set_yticklabels(fn.stage)
    ax.invert_yaxis()
    for i, v in enumerate(fn.events):
        ax.text(v, i, ' %d' % v, va='center')
    ax.set_xlabel('events')
    ax.set_title('From the download to the accepted events')
    save_figure(fig, out, 'event_funnel.png')
    plt.close(fig)
    # --- availability: station x time ---
    sq2 = sq.merge(ev[['event', 'origin', 'status']], on='event', how='left')
    sq2['t'] = pd.to_datetime(sq2.origin.str.replace('Z', ''))
    stas = sorted(sq2.station.unique())
    fig, ax = plt.subplots(figsize=(11, 0.18 * len(stas) + 1.5))
    yi = dict((s, i) for i, s in enumerate(stas))
    for label, sel, col in (
        ('recorded', np.ones(len(sq2), bool), '0.8'),
        ('usable P', sq2.usable_P.fillna(False).astype(bool).values, 'k'),
    ):
        d = sq2[sel]
        ax.scatter(d.t, d.station.map(yi), s=6, c=col, marker='|', label=label)
    ax.set_yticks(range(len(stas)))
    ax.set_yticklabels(stas, fontsize=6)
    ax.invert_yaxis()
    ax.set_title('Station availability per event (gray: recorded, black: usable P)')
    ax.legend(loc='upper right', fontsize=7, markerscale=2, frameon=False)
    save_figure(fig, out, 'station_availability.png')
    plt.close(fig)

    # --- integrity per station ---
    def cnt(col):
        return (
            sq.groupby('station')[col].apply(lambda x: int(x.fillna(False).astype(bool).sum()))
            if col in sq
            else 0
        )

    integ = sq.groupby('station').agg(events=('event', 'nunique')).reset_index()
    for c in 'ZNE':
        for k in ('dead', 'clipped', 'flat'):
            integ['%s_%s' % (k, c)] = (
                integ.station.map(cnt('%s_%s' % (k, c))) if ('%s_%s' % (k, c)) in sq else 0
            )
    gcols = [c for c in ('gap_frac_Z', 'gap_frac_N', 'gap_frac_E') if c in sq]
    pcols = [c for c in ('pieces_Z', 'pieces_N', 'pieces_E') if c in sq]
    integ['events_with_gaps'] = integ.station.map(
        sq.assign(g=(sq[gcols].fillna(0) > 0).any(axis=1)).groupby('station').g.sum()
    )
    integ['events_split_files'] = integ.station.map(
        sq.assign(g=(sq[pcols].fillna(1) > 1).any(axis=1)).groupby('station').g.sum()
    )
    integ['spikes'] = integ.station.map(
        sq[['station'] + [c for c in ('spikes_Z', 'spikes_N', 'spikes_E') if c in sq]]
        .fillna(0)
        .groupby('station')
        .sum()
        .sum(axis=1)
    )
    integ['usable_P'] = integ.station.map(
        sq.groupby('station').usable_P.apply(lambda x: int(x.fillna(False).astype(bool).sum()))
    )
    integ['usable_P_pct'] = (100.0 * integ.usable_P / integ.events).round(1)
    save_csv(integ, out, 'integrity_by_station.csv')
    show = [c for c in integ.columns if c not in ('station', 'events', 'usable_P', 'usable_P_pct')]
    M = integ[show].values.astype(float)
    fig, ax = plt.subplots(figsize=(0.55 * len(show) + 2.5, 0.2 * len(integ) + 1.5))
    ax.imshow(np.where(M > 0, 1, 0), cmap='Greys', vmin=0, vmax=2.5, aspect='auto')
    for i in range(M.shape[0]):
        for j in range(M.shape[1]):
            if M[i, j] > 0:
                ax.text(j, i, '%d' % M[i, j], ha='center', va='center', fontsize=6)
    ax.set_xticks(range(len(show)))
    ax.set_xticklabels(show, rotation=60, ha='right', fontsize=7)
    ax.set_yticks(range(len(integ)))
    ax.set_yticklabels(integ.station, fontsize=6)
    ax.set_title('Integrity problems per station (number of events)')
    save_figure(fig, out, 'integrity_by_station.png')
    plt.close(fig)


# -----------------------------------------------------------------------------
# snr: phase visibility against magnitude and distance
# -----------------------------------------------------------------------------
def cmd_snr(a):
    import pandas as pd

    plt = pyplot()
    out = output_dir(a)
    ev, sq = read_tables(a.qc)
    d = sq.merge(ev[['event', 'mag', 'status']], on='event', how='left')
    fig, axs = plt.subplots(1, 4, figsize=(17, 3.8))
    for ax, col, lab in (
        (axs[0], 'snr_P_sp', 'P SNR (vertical, 0.5-2 Hz)'),
        (axs[1], 'snr_S', 'S SNR (transverse, 0.03-0.2 Hz)'),
    ):
        x = d[[col, 'mag', 'dist']].dropna()
        sc = ax.scatter(
            x.mag + np.random.RandomState(0).uniform(-0.03, 0.03, len(x)),
            x[col],
            c=x.dist,
            s=5,
            cmap='viridis',
            alpha=0.6,
        )
        ax.set_yscale('log')
        ax.axhline(3, color='k', lw=0.8, ls='--')
        ax.set_xlabel('magnitude')
        ax.set_ylabel(lab)
        fig.colorbar(sc, ax=ax, label='distance (deg)')
    mb = np.arange(5.5, 9.01, 0.25)
    d['mbin'] = pd.cut(d.mag, mb, right=False)
    vis = (
        d.groupby('mbin')
        .agg(
            stations=('station', 'size'),
            P_detected_pct=('snr_P_sp', lambda x: 100.0 * np.mean(x >= 3)),
            S_detected_pct=('snr_S', lambda x: 100.0 * np.mean(x >= 3)),
        )
        .reset_index()
    )
    vis = vis[vis.stations > 0]
    axs[2].plot([i.mid for i in vis.mbin], vis.P_detected_pct, 'k-o', ms=3, label='P')
    axs[2].plot([i.mid for i in vis.mbin], vis.S_detected_pct, '-o', color='0.55', ms=3, label='S')
    axs[2].set_xlabel('magnitude')
    axs[2].set_ylabel('station-events with SNR >= 3 (%)')
    axs[2].legend(frameon=False)
    db = np.arange(30, 100.01, 10)
    d['dbin'] = pd.cut(d.dist, db, right=False)
    visd = (
        d.groupby('dbin')
        .agg(
            stations=('station', 'size'),
            P_detected_pct=('snr_P_sp', lambda x: 100.0 * np.mean(x >= 3)),
            S_detected_pct=('snr_S', lambda x: 100.0 * np.mean(x >= 3)),
        )
        .reset_index()
    )
    visd = visd[visd.stations > 0]
    axs[3].plot([i.mid for i in visd.dbin], visd.P_detected_pct, 'k-o', ms=3, label='P')
    axs[3].plot([i.mid for i in visd.dbin], visd.S_detected_pct, '-o', color='0.55', ms=3, label='S')
    for i_, r_ in visd.iterrows():
        axs[3].text(r_.dbin.mid, 2, '%d' % r_.stations, ha='center', fontsize=6, color='0.4')
    axs[3].set_xlabel('distance (deg)')
    axs[3].set_ylabel('station-events with SNR >= 3 (%)')
    axs[3].set_ylim(0, 100)
    axs[3].legend(frameon=False)
    axs[0].set_title('(a) P visibility', fontsize=9)
    axs[1].set_title('(b) S visibility', fontsize=9)
    axs[2].set_title('(c) Detection rate by magnitude', fontsize=9)
    axs[3].set_title('(d) Detection rate by distance (numbers: station-events)', fontsize=9)
    plt.tight_layout()
    save_figure(fig, out, 'phase_visibility.png')
    plt.close(fig)
    vis['mbin'] = vis.mbin.astype(str)
    save_csv(vis.round(1), out, 'visibility_by_magnitude.csv')
    visd['dbin'] = visd.dbin.astype(str)
    save_csv(visd.round(1), out, 'visibility_by_distance.csv')


# -----------------------------------------------------------------------------
# response: sensor types compared
# -----------------------------------------------------------------------------
def cmd_response(a):
    import pandas as pd
    from obspy import read_inventory

    plt = pyplot()
    out = output_dir(a)
    inv = read_inventory(a.inventory)
    seen = {}
    for net in inv:
        for sta in net:
            for ch in sta:
                if ch.code.endswith('Z') and ch.sensor is not None and ch.sensor.model not in seen:
                    seen[ch.sensor.model] = ch
    f = np.logspace(math.log10(0.005), math.log10(45.0), 400)
    fig, axs = plt.subplots(1, 2, figsize=(11, 3.8))
    rows, ref = [], None
    for k, (model, ch) in enumerate(sorted(seen.items())):
        r = ch.response
        if not r.response_stages:
            print('  %s: sensitivity only, no response stages to plot' % model)
            continue
        h = r.get_evalresp_response_for_frequencies(f, output='VEL') / r.instrument_sensitivity.value
        ph = np.unwrap(np.angle(h))
        ls = ('-', '--', ':', '-.')[k % 4]
        axs[0].semilogx(f, 20 * np.log10(np.abs(h)), 'k', ls=ls, label=model)
        axs[1].semilogx(f, np.degrees(ph), 'k', ls=ls, label=model)
        row = dict(sensor=model, channel_example=ch.code)
        for fx in (0.05, 0.1, 0.5, 1.0, 2.0):
            row['phase_delay_ms_%gHz' % fx] = float(np.interp(fx, f, -ph / (2 * np.pi * f)) * 1000)
        for fx in (10.0, 20.0, 30.0):
            row['amp_dB_%gHz' % fx] = float(np.interp(fx, f, 20 * np.log10(np.abs(h))))
        rows.append(row)
    if not rows:
        print('No response stages in %s: use a StationXML requested at response level.' % a.inventory)
        return
    axs[0].set_ylim(-20, 10)
    axs[0].set_xlabel('frequency (Hz)')
    axs[0].set_ylabel('amplitude relative to sensitivity (dB)')
    axs[1].set_xlabel('frequency (Hz)')
    axs[1].set_ylabel('phase (deg)')
    for ax in axs:
        ax.grid(alpha=0.3, which='both')
        ax.axvspan(0.5, 2.0, color='0.9', zorder=0)
        ax.legend(frameon=False, fontsize=8)
    axs[0].set_title('Velocity response (gray band: P picking band)')
    axs[1].set_title('Phase response')
    plt.tight_layout()
    save_figure(fig, out, 'sensor_responses.png')
    plt.close(fig)
    t = pd.DataFrame(rows)
    ref = t.iloc[0]
    for fx in (0.5, 1.0, 2.0):
        c = 'phase_delay_ms_%gHz' % fx
        t['diff_vs_%s_ms_%gHz' % (ref.sensor.split('/')[0], fx)] = t[c] - ref[c]
    save_csv(t.round(2), out, 'sensor_response_comparison.csv')


# -----------------------------------------------------------------------------
# noise-hour: noise against local time of day
# -----------------------------------------------------------------------------
def cmd_noise_hour(a):
    import pandas as pd
    import pickle

    plt = pyplot()
    out = output_dir(a)
    ev, sq = read_tables(a.qc)
    n = int(round(math.log(PSD_FMAX / PSD_FMIN, 2) * 8)) + 1
    fc = PSD_FMIN * 2.0 ** (np.arange(n) / 8.0)
    sel = (fc >= 1.0) & (fc <= 5.0)
    tp = dict(((r.event, r.station), r.tP_pred) for r in sq[['event', 'station', 'tP_pred']].itertuples())
    rows = []
    for f in sorted(glob.glob(os.path.join(a.qc, 'work', '*.pkl'))):
        with open(f, 'rb') as fh:
            r = pickle.load(fh)
        for q in r.get('psd', []):
            if q['comp'] != 'Z':
                continue
            if len(q['psd']) != n:
                sys.exit(
                    'PSD grid differs from the quality-control settings (PSD_FMIN and PSD_FMAX in common.py)'
                )
            t = tp.get((q['event'], q['station']))
            if not isinstance(t, str):
                continue
            mid = pd.Timestamp(t.replace('Z', '')) - pd.Timedelta(seconds=58)  # center of the 100 s window
            loc = mid + pd.Timedelta(hours=a.utc_offset)
            rows.append(
                dict(
                    event=q['event'],
                    station=q['station'],
                    local_hour=loc.hour + loc.minute / 60.0,
                    noise_1_5Hz_dB=float(np.nanmean(np.asarray(q['psd'], float)[sel])),
                )
            )
    if not rows:
        sys.exit('no noise spectra found in %s/work' % a.qc)
    d = pd.DataFrame(rows)
    save_csv(d.round(2), out, 'noise_by_window.csv')
    day = (d.local_hour >= 8) & (d.local_hour < 18)
    night = (d.local_hour >= 22) | (d.local_hour < 4)
    s = (
        d.groupby('station')
        .apply(
            lambda g: pd.Series(
                dict(
                    n_day=int(((g.local_hour >= 8) & (g.local_hour < 18)).sum()),
                    n_night=int(((g.local_hour >= 22) | (g.local_hour < 4)).sum()),
                    median_day_dB=g.noise_1_5Hz_dB[(g.local_hour >= 8) & (g.local_hour < 18)].median(),
                    median_night_dB=g.noise_1_5Hz_dB[(g.local_hour >= 22) | (g.local_hour < 4)].median(),
                )
            )
        )
        .reset_index()
    )
    s['day_minus_night_dB'] = s.median_day_dB - s.median_night_dB
    save_csv(s.round(2), out, 'noise_day_night.csv')
    fig, axs = plt.subplots(1, 2, figsize=(12, 4))
    axs[0].scatter(d.local_hour, d.noise_1_5Hz_dB, s=3, c='0.6')
    hb = d.assign(h=np.floor(d.local_hour)).groupby('h').noise_1_5Hz_dB.median()
    axs[0].plot(hb.index + 0.5, hb.values, 'k-o', ms=3, label='network median per hour')
    axs[0].set_xlabel('local time (UTC%+g h)' % a.utc_offset)
    axs[0].set_ylabel('vertical noise, 1-5 Hz (dB rel. 1 (m/s^2)^2/Hz)')
    axs[0].set_xlim(0, 24)
    axs[0].legend(frameon=False)
    axs[0].set_title('Noise against time of day (%d day / %d night windows)' % (day.sum(), night.sum()))
    ok = s.dropna(subset=['day_minus_night_dB']).sort_values('day_minus_night_dB')
    axs[1].barh(range(len(ok)), ok.day_minus_night_dB, color='0.45')
    axs[1].set_yticks(range(len(ok)))
    axs[1].set_yticklabels(ok.station, fontsize=6)
    axs[1].set_xlabel('day minus night (dB)')
    axs[1].set_title('Day-night noise difference per station')
    plt.tight_layout()
    save_figure(fig, out, 'noise_vs_local_hour.png')
    plt.close(fig)
