"""Quality control, station noise, event screening and teleseismic P and S picking.

Input: event folders written by SOD (Event_YYYY_MM_DD_hh_mm_ss/*.sac). They are only
read, never written; the output folder must be outside the data folder, and the run
confirms at the end that no input file changed.

Per event (events run in parallel):
  1. integrity checks on raw counts: dead, clipped, flat segments, gaps, spikes
  2. counts to ground velocity (full response if the StationXML has response stages,
     otherwise the overall sensitivity, which gives flat-band velocity)
  3. pre-event noise spectrum of every channel (McNamara and Buland 2004)
  4. P (vertical) and S (transverse) signal-to-noise ratios against ak135
  5. screening: no usable data, no P detected, P not coherent across the array
  6. P and S arrival times: cross-correlation and stacking of an onset function,
     then of the waveforms (ICCS), multi-channel cross-correlation (MCCC; VanDecar and
     Crosson 1990), causal phase-weighted stack and AIC onset on the stack for the
     absolute time. Windows end before later phases with a different slowness
     (PcP for P; ScS, SKS and SKKS for S).
  7. processed waveforms: velocity, Z/R/T, band-passed, resampled, normalized, the same
     first sample on all components, time zero at the earthquake origin, so the SAC
     headers T0 (P pick) and T1 (S pick) are travel times
  8. one figure per event (record sections with predictions and picks)

After all events: event and station tables, noise probability density plots per
station, a network noise overview, a station performance table, a pick summary
figure, and the two pick files
    picks/picks_absolute.csv   absolute travel times (relocation, absolute tomography)
    picks/picks_relative.csv   relative delays (relative travel-time tomography)

Usage
-----
    python rfo.py qc --data data/raw_sac --inventory stations.xml --out qc_out --workers 6

    --events "Event_2015_*"   process a subset of the event folders
    --overwrite               process again events that are already done
    --no-figures, --no-waveforms
    --aggregate-only          rebuild tables and plots from finished events
    An interrupted run resumes by itself: finished events are skipped.

All thresholds are in CONFIG below; docs/METHODS_QC.md describes the method.
Requires Python 3.9 or later, ObsPy 1.3 or later, NumPy, SciPy, pandas and matplotlib.
"""

import argparse
import glob
import math
import os
import pickle
import shutil
import sys
import time
import traceback
import warnings
from concurrent.futures import ProcessPoolExecutor, as_completed

import numpy as np
from scipy.signal import hilbert, welch
import matplotlib

matplotlib.use('Agg')

warnings.filterwarnings('ignore')

# =============================================================================
# CONFIGURATION
# =============================================================================
CONFIG = {
    # --- data / metadata ---------------------------------------------------------
    'network_rename': {'5H2013': '5H'},    # SOD writes temporary-network ids with the start year
    'evdp_in_meters': True,                # SOD writes SAC EVDP in meters
    'model': 'ak135',                      # reference Earth model for predictions and residuals
    'figure_phases': ['P', 'pP', 'sP', 'PcP', 'PP', 'S', 'sS', 'ScS', 'SKS', 'SS'],

    # --- integrity (raw counts) --------------------------------------------------
    'full_scale_counts': 2 ** 23 - 1,      # 24-bit digitizer
    'clip_fraction': 0.95,
    'flat_run_s': 1.0,                     # identical consecutive samples longer than this = flat segment
    'spike_mad': 50.0,                     # isolated one-sample spikes (x MAD of the 2nd difference)

    # --- signal-to-noise (windows in s relative to the ak135 prediction) ---------
    'P_band': (0.5, 2.0),                  # short-period P
    'P_band_bb': (0.05, 1.0),              # broadband P (receiver-function band)
    'S_band': (0.03, 0.2),                 # S on the transverse component
    'snr_sig_P': (-1.0, 9.0), 'snr_noise_P': (-65.0, -5.0),
    'snr_sig_S': (-2.0, 18.0), 'snr_noise_S': (-65.0, -5.0),

    # --- screening ---------------------------------------------------------------
    'min_snr_trace_P': 2.0,                # a trace enters P picking only above this SNR
    'min_snr_trace_S': 2.0,
    'min_snr_detect': 3.0,                 # "phase detected" at a station
    'min_stations_alive': 3,               # fewer live vertical channels -> event rejected (no usable data)
    'min_stations_P': 3,                   # fewer stations with detected P -> event rejected (no main phase)
    'min_stations_S': 3,                   # fewer stations with detected S -> S not picked
    'min_picks': 4,                        # fewer A-C P picks -> event rejected (P not coherent)
    'min_picks_S': 3,                      # fewer A-C S picks -> all S picks of the event set to D
    'noise_contamination_ratio': 3.0,      # aggregation: pre-event noise this many times the station norm -> flag

    # --- picking -----------------------------------------------------------------
    #   onset_band: CAUSAL filter for the absolute onset (zero-phase filters create
    #               acausal precursors that bias onsets early by up to ~1 s)
    #   cf: coarse ICCS on an onset characteristic function (smoothed envelope of the
    #       causal onset band) - step-like, no oscillations, so it cannot skip a cycle
    #   wf: picking-band waveform ICCS; lag limit below half a dominant period
    #   interfere: later phases with a different slowness. Every correlation window
    #       ends 'margin' s before the earliest of them (never shorter than min_end):
    #       PcP arrives only 4-10 s after P at 76-82 deg
    #   skip_tol: max onset-function lag [s] against the leave-one-out stack
    #   classes: (min CC with the leave-one-out stack, max MCCC sigma [s], min SNR) for A, B, C
    #   max_dist: beyond ~80 deg S overlaps SKS/SKKS/ScS, so S is not picked there
    #   bootstrap: resamples of the stations to estimate the uncertainty of the common onset
    #       (sigma_onset); absolute uncertainty = sqrt(sigma_mccc^2 + sigma_onset^2)
    'pick': {
        'P': dict(comp='Z', band=(0.5, 2.0), rate=100.0, onset_band=(0.3, 8.0),
                  cf=dict(win=(-8.0, 6.0), lag=4.0, smooth=2.0),
                  wf=dict(win=(-2.0, 6.0), lag=0.35),
                  mccc_lag=0.25, mccc_pair_cc=0.5, iters=6, cc_stack=0.6,
                  aic_win=(-6.0, 8.0), aic_lookback=3.0, single_aic=(-1.5, 0.75),
                  skip_tol=0.3, interfere=('PcP',), margin=1.0, min_end=dict(cf=3.0, wf=2.5),
                  max_dist=180.0, bootstrap=200,
                  classes=((0.85, 0.05, 5.0), (0.75, 0.10, 3.0), (0.65, 0.20, 2.0))),
        'S': dict(comp='T', band=(0.03, 0.2), rate=10.0, onset_band=(0.02, 0.5),
                  cf=dict(win=(-25.0, 15.0), lag=10.0, smooth=0.2),
                  wf=dict(win=(-5.0, 10.0), lag=1.5),
                  mccc_lag=1.0, mccc_pair_cc=0.5, iters=6, cc_stack=0.6,
                  aic_win=(-20.0, 12.0), aic_lookback=10.0, single_aic=(-6.0, 3.0),
                  skip_tol=1.5, interfere=('ScS', 'SKS', 'SKKS'), margin=3.0, min_end=dict(cf=8.0, wf=6.0),
                  max_dist=80.0, bootstrap=200,
                  classes=((0.90, 0.25, 6.0), (0.85, 0.50, 4.0), (0.75, 1.00, 3.0))),
    },
    'vp_surface': 5.8, 'vs_surface': 3.36,  # km/s, station-elevation correction of residuals

    # --- noise PSD ---------------------------------------------------------------
    'psd_window_s': 100.0,                 # pre-P noise window length
    'psd_gap_before_P_s': 8.0,             # the window ends this long before predicted P
    'psd_nperseg_fraction': 0.5,           # Welch segment = window x fraction, 75 % overlap
    'psd_fmin': 0.05, 'psd_fmax': 30.0,    # limited by ~120 s of pre-event data per event file
    'psd_db_bins': (-200.0, -40.0, 1.0),

    # --- processed waveforms -----------------------------------------------------
    'write_waveforms': True,
    'write_only_usable': True,             # only stations with usable P (SNR >= min_snr_trace_P, no gap)
    'out_components': 'ZRT',               # any of Z, N, E, R, T
    'out_band': (0.02, 4.0),               # Hz, zero-phase Butterworth band-pass
    'out_rate': 20.0,                      # Hz; None keeps 100 Hz
    'normalize': 'station',                # 'station' (one factor for all components), 'trace', or None
    'time_reference': 'origin',            # t = 0 at the earthquake origin ('start' or 'P' also possible)

    # --- figures -----------------------------------------------------------------
    'figures': True,
    'plot_rejected': True,
}

VERSION = '1.1'
CLASS_NAMES = ('A', 'B', 'C')
_G = {}  # per-process globals: inventory, travel-time model


# =============================================================================
# SMALL HELPERS
# =============================================================================
def _rms(x):
    x = np.asarray(x, float)
    return float(np.sqrt(np.mean(x * x))) if x.size else np.nan


def _cut(tr, t0, n):
    """n samples of tr starting at UTCDateTime t0 (nearest sample), or None if outside."""
    i0 = int(round((t0 - tr.stats.starttime) * tr.stats.sampling_rate))
    if n <= 0 or i0 < 0 or i0 + n > tr.stats.npts:
        return None
    return np.asarray(tr.data[i0 : i0 + n], float)


def _win(tr, t0, t1):
    return _cut(tr, t0, int(round((t1 - t0) * tr.stats.sampling_rate)))


def _rmsw(tr, t0, t1):
    x = _win(tr, t0, t1)
    return _rms(x) if x is not None and len(x) else np.nan


def _resample(t, rate, fmax):
    if rate and rate < t.stats.sampling_rate:
        fac = int(round(t.stats.sampling_rate / rate))
        if fac > 1 and abs(t.stats.sampling_rate / fac - rate) < 1e-6 and fmax < 0.4 * rate:
            t.decimate(fac, no_filter=True)  # data already band-limited
        else:
            t.resample(rate)
    return t


def _bandpass(tr, f1, f2, rate=None):
    """Zero-phase Butterworth band-pass (alignment, SNR, display)."""
    t = tr.copy()
    t.data = np.asarray(t.data, float)
    t.detrend('demean')
    t.detrend('linear')
    t.taper(0.02, type='hann', max_length=10.0)
    t.filter('bandpass', freqmin=f1, freqmax=min(f2, 0.45 * t.stats.sampling_rate), corners=4, zerophase=True)
    return _resample(t, rate, f2)


def _causal(tr, f1, f2, rate=None):
    """Causal (minimum-phase) band-pass, used only for onset times."""
    t = tr.copy()
    t.data = np.asarray(t.data, float)
    t.detrend('demean')
    t.detrend('linear')
    t.taper(0.02, type='hann', max_length=10.0)
    t.filter(
        'bandpass', freqmin=f1, freqmax=min(f2, 0.45 * t.stats.sampling_rate), corners=2, zerophase=False
    )
    return _resample(t, rate, f2)


def _charfun(tr, smooth_hz):
    """Onset characteristic function: smoothed envelope of a causally filtered trace."""
    e = tr.copy()
    e.data = np.abs(hilbert(e.data))
    e.filter('lowpass', freq=min(smooth_hz, 0.4 * e.stats.sampling_rate), corners=2, zerophase=True)
    return e


def _ncc_valid(long_x, short_x):
    """Normalized cross-correlation of short_x slid along long_x (all 'valid' lags)."""
    n = len(short_x)
    s = short_x - short_x.mean()
    s = s / (np.linalg.norm(s) + 1e-30)
    c = np.correlate(long_x, s, mode='valid')
    cs = np.concatenate([[0.0], np.cumsum(long_x * long_x)])
    cm = np.concatenate([[0.0], np.cumsum(long_x)])
    e2 = (cs[n:] - cs[:-n]) - (cm[n:] - cm[:-n]) ** 2 / n
    return c / np.sqrt(np.maximum(e2, 1e-30))


def _peak(c):
    """Sub-sample (parabolic) position and value of the maximum of c."""
    k = int(np.argmax(c))
    if 0 < k < len(c) - 1:
        y0, y1, y2 = c[k - 1], c[k], c[k + 1]
        den = y0 - 2 * y1 + y2
        if den < 0:
            d = 0.5 * (y0 - y2) / den
            return k + d, float(y1 - 0.25 * (y0 - y2) * d)
    return float(k), float(c[k])


def _aic(x):
    """Maeda (1985) AIC function; its minimum marks the onset."""
    x = np.asarray(x, float)
    n = len(x)
    a = np.full(n, np.inf)
    if n < 6:
        return a
    k = np.arange(1, n - 1)
    c1, c2 = np.cumsum(x), np.cumsum(x * x)
    m1 = c1[k - 1] / k
    v1 = c2[k - 1] / k - m1**2
    m2 = (c1[-1] - c1[k - 1]) / (n - k)
    v2 = (c2[-1] - c2[k - 1]) / (n - k) - m2**2
    a[k] = k * np.log(np.maximum(v1, 1e-30)) + (n - k - 1) * np.log(np.maximum(v2, 1e-30))
    a[:2] = np.inf
    a[-2:] = np.inf
    return a


def _pws(segs, nu=2.0):
    """Phase-weighted stack (Schimmel & Paulssen 1997); returns (pws, linear stack)."""
    segs = np.asarray(segs, float)
    lin = segs.mean(axis=0)
    an = hilbert(segs, axis=1)
    ph = an / (np.abs(an) + 1e-30)
    return lin * np.abs(ph.mean(axis=0)) ** nu, lin


def _classify_snr(x):
    if x is None or not np.isfinite(x):
        return 'n/a'
    return 'clear' if x >= 5 else 'visible' if x >= 3 else 'weak' if x >= 2 else 'none'


def _gap_in(g, t0, t1, sr):
    """True if the zero-filled gap mask g = (fraction, start, mask) has samples between t0 and t1."""
    if not g:
        return False
    i0 = max(0, int(round((t0 - g[1]) * sr)))
    i1 = max(0, int(round((t1 - g[1]) * sr)))
    return bool(g[2][i0:i1].any())


def _utc(s):
    from obspy import UTCDateTime

    return UTCDateTime(s)


# =============================================================================
# METADATA AND TRAVEL TIMES
# =============================================================================
def inventory_info(inv):
    """Flatten an Inventory into {seed_id: [epoch dicts]} (picklable)."""
    info = {}
    for net in inv:
        for sta in net:
            for ch in sta.channels:
                sid = '%s.%s.%s.%s' % (net.code, sta.code, ch.location_code, ch.code)
                r = ch.response
                sens = (
                    r.instrument_sensitivity.value if (r is not None and r.instrument_sensitivity) else None
                )
                info.setdefault(sid, []).append(
                    dict(
                        lat=ch.latitude if ch.latitude is not None else sta.latitude,
                        lon=ch.longitude if ch.longitude is not None else sta.longitude,
                        elev=ch.elevation if ch.elevation is not None else sta.elevation,
                        azimuth=ch.azimuth,
                        dip=ch.dip,
                        sens=sens,
                        has_stages=bool(r is not None and r.response_stages),
                        sensor=(ch.sensor.model if ch.sensor and ch.sensor.model else ''),
                        start=ch.start_date,
                        end=ch.end_date,
                        site=(sta.site.name if sta.site and sta.site.name else ''),
                    )
                )
    return info


def channel_meta(seed_id, t):
    lst = _G['inv_info'].get(seed_id, [])
    for e in lst:
        if (e['start'] is None or e['start'] <= t) and (e['end'] is None or t <= e['end']):
            return e
    return lst[0] if lst else None


class TravelTimes(object):
    """First-arrival times and ray parameters, interpolated over the event's distance span."""

    def __init__(self, model, depth_km, dists, phases):
        lo, hi = max(0.5, min(dists) - 0.1), min(179.5, max(dists) + 0.1)
        n = max(3, int(math.ceil((hi - lo) / 0.25)) + 1)
        self.grid = np.linspace(lo, hi, n)
        phases = sorted(set(phases))
        self.t = dict((ph, np.full(n, np.nan)) for ph in phases)
        self.p = dict((ph, np.full(n, np.nan)) for ph in phases)
        for i, d in enumerate(self.grid):
            for a in model.get_travel_times(
                source_depth_in_km=depth_km, distance_in_degree=float(d), phase_list=phases
            ):
                if a.name in self.t and np.isnan(self.t[a.name][i]):
                    self.t[a.name][i] = a.time
                    self.p[a.name][i] = a.ray_param_sec_degree

    def _interp(self, y, d):
        if y is None or not np.isfinite(y).any():
            return np.nan
        ok = np.isfinite(y)
        if d < self.grid[ok][0] - 1e-6 or d > self.grid[ok][-1] + 1e-6:
            return np.nan
        return float(np.interp(d, self.grid[ok], y[ok]))

    def time(self, ph, d):
        return self._interp(self.t.get(ph), d)

    def rayp_km(self, ph, d):
        return self._interp(self.p.get(ph), d) / 111.195


def _init_worker(inv_path, cfg):
    from obspy import read_inventory
    from obspy.taup import TauPyModel

    _G['cfg'] = cfg
    _G['inv'] = read_inventory(inv_path) if inv_path else None
    _G['inv_info'] = inventory_info(_G['inv']) if _G['inv'] is not None else {}
    _G['model'] = TauPyModel(cfg['model'])


# =============================================================================
# READING AND PER-TRACE QC
# =============================================================================
def read_event(folder, cfg):
    from obspy import read, UTCDateTime

    st = read(os.path.join(folder, '*.sac'))
    pieces = {}
    for tr in st:
        tr.stats.network = cfg['network_rename'].get(tr.stats.network, tr.stats.network)
        tr.stats.location = ''
        pieces[tr.id] = pieces.get(tr.id, 0) + 1
    s = st[0].stats.sac
    origin = st[0].stats.starttime - (float(s.b) - float(s.get('o', 0.0)))
    origin = UTCDateTime(round(origin.timestamp, 3))  # SAC stores the origin to the millisecond
    depth = float(s.evdp) / 1000.0 if cfg['evdp_in_meters'] else float(s.evdp)
    ev = dict(
        name=os.path.basename(folder.rstrip('/')),
        origin=origin,
        lat=float(s.evla),
        lon=float(s.evlo),
        depth_km=depth,
        mag=float(s.get('mag', np.nan)),
    )
    hdr = {}
    for tr in st:
        h = tr.stats.sac
        hdr.setdefault(
            tr.stats.station,
            dict(
                gcarc=float(h.gcarc),
                baz=float(h.baz),
                az=float(h.az),
                stla=float(h.stla),
                stlo=float(h.stlo),
                stel=float(h.get('stel', 0.0)),
            ),
        )
    gaps = {}
    st.merge(method=1)  # joins *.1/*.2 pieces of a split channel
    for tr in st:
        if np.ma.isMaskedArray(tr.data):
            m = np.ma.getmaskarray(tr.data)
            gaps[tr.id] = (float(m.mean()), tr.stats.starttime, m)
            tr.data = tr.data.filled(0)
        tr.data = np.asarray(tr.data, np.float64)
    return st, ev, hdr, pieces, gaps


def integrity(tr, cfg, gapmask=None):
    x = np.asarray(tr.data, float)
    out = {
        'dead': bool(np.std(x) == 0.0),
        'clipped': bool(np.max(np.abs(x)) >= cfg['clip_fraction'] * cfg['full_scale_counts']),
    }
    xd = x.copy()
    if gapmask is not None:
        xd[gapmask] = np.nan  # zero-filled gaps are not flat segments
    same = np.concatenate([[0], (np.diff(xd) == 0).astype(np.int8), [0]])
    e = np.diff(same)
    starts, ends = np.where(e == 1)[0], np.where(e == -1)[0]
    run = int((ends - starts).max()) + 1 if starts.size else 1
    out['max_flat_s'] = round(run * tr.stats.delta, 3)
    out['flat'] = bool(out['max_flat_s'] > cfg['flat_run_s'] and not out['dead'])
    if not out['dead']:
        d2 = np.diff(x, 2)
        mad = np.median(np.abs(d2 - np.median(d2))) + 1e-12
        big = np.abs(d2) > cfg['spike_mad'] * mad
        out['spikes'] = int(np.sum(big[1:-1] & ~big[:-2] & ~big[2:]))
    else:
        out['spikes'] = 0
    out['dc_counts'] = float(np.mean(x))
    return out


def to_velocity(tr, meta):
    """Counts -> m/s. Full response if available, else flat-band sensitivity, else counts."""
    if meta is None or not meta.get('sens'):
        return tr.copy(), 'counts'
    t = tr.copy()
    if meta.get('has_stages') and _G.get('inv') is not None:
        try:
            t.detrend('demean')
            t.taper(0.02, max_length=10.0)
            sr = t.stats.sampling_rate
            t.remove_response(inventory=_G['inv'], output='VEL', pre_filt=(0.005, 0.01, 0.4 * sr, 0.45 * sr))
            return t, 'response'
        except Exception:
            t = tr.copy()
    t.data = t.data / float(meta['sens'])
    return t, 'sensitivity'


def rotate_components(z, n, e, metas, baz):
    """Common window and first sample for all components; metadata orientation; NE -> RT."""
    alive = [x for x in (z, n, e) if x is not None]
    t0 = max(x.stats.starttime for x in alive)
    t1 = min(x.stats.endtime for x in alive)
    for x in alive:
        x.trim(t0, t1, nearest_sample=True)
    m = min(x.stats.npts for x in alive)
    for x in alive:
        x.data = x.data[:m]
        if abs(x.stats.starttime - alive[0].stats.starttime) < 0.5 * x.stats.delta:
            x.stats.starttime = alive[0].stats.starttime
    R = T = None
    if n is not None and e is not None:
        mN, mE = metas.get('N') or {}, metas.get('E') or {}
        aN = math.radians(mN.get('azimuth') if mN.get('azimuth') is not None else 0.0)
        aE = math.radians(mE.get('azimuth') if mE.get('azimuth') is not None else 90.0)
        if abs(aN) > 1e-6 or abs(aE - math.pi / 2) > 1e-6:  # non-nominal orientation in metadata
            nn = n.data * math.cos(aN) + e.data * math.cos(aE)
            ee = n.data * math.sin(aN) + e.data * math.sin(aE)
            n.data, e.data = nn, ee
        b = math.radians(baz)
        R = n.copy()
        R.data = -n.data * math.cos(b) - e.data * math.sin(b)
        R.stats.channel = n.stats.channel[:2] + 'R'
        T = n.copy()
        T.data = n.data * math.sin(b) - e.data * math.cos(b)
        T.stats.channel = n.stats.channel[:2] + 'T'
    return R, T


# =============================================================================
# NOISE PSD (McNamara & Buland 2004 style, acceleration, dB rel. 1 (m/s^2)^2/Hz)
# =============================================================================
def psd_freqs(cfg):
    n = int(round(math.log(cfg['psd_fmax'] / cfg['psd_fmin'], 2) * 8)) + 1
    return cfg['psd_fmin'] * 2.0 ** (np.arange(n) / 8.0)


def noise_psd(vel_tr, t_end, cfg, fc):
    x = _win(vel_tr, t_end - cfg['psd_window_s'], t_end)
    if x is None or np.std(x) == 0:
        return None
    fs = vel_tr.stats.sampling_rate
    nper = int(len(x) * cfg['psd_nperseg_fraction'])
    f, p = welch(
        x - x.mean(),
        fs=fs,
        window='hann',
        nperseg=nper,
        noverlap=int(0.75 * nper),
        detrend='linear',
        scaling='density',
    )
    pa = p * (2 * np.pi * f) ** 2
    out = np.full(len(fc), np.nan)
    for i, c in enumerate(fc):
        sel = (f >= c / math.sqrt(2)) & (f <= c * math.sqrt(2)) & (f > 0)
        if sel.any():
            out[i] = 10 * np.log10(np.mean(pa[sel]) + 1e-40)
    return out


# =============================================================================
# PICKING
# =============================================================================
def _iccs(trs, t0, win, maxlag, iters, cc_min):
    """Iterative cross-correlation and stacking (Lou et al. 2013 style).
    The median shift is removed every iteration so the stack cannot drift; a trace whose
    shift would exceed maxlag keeps its previous shift and is excluded from the stack."""
    keys = list(trs)
    fs = trs[keys[0]].stats.sampling_rate
    n = int(round((win[1] - win[0]) * fs))
    m = int(round(maxlag * fs))
    shift = dict((k, 0.0) for k in keys)
    cc = dict((k, np.nan) for k in keys)
    weights = dict((k, 1.0) for k in keys)
    stack = None
    for _ in range(iters):
        segs = {}
        for k in keys:
            x = _cut(trs[k], t0[k] + shift[k] + win[0], n)
            if x is not None and np.std(x) > 0:
                segs[k] = (x - x.mean()) / np.std(x)
        use = [k for k in segs if weights.get(k, 0) > 0]
        if len(use) < 2:
            return shift, cc, None
        stack = np.average([segs[k] for k in use], axis=0, weights=[weights[k] for k in use])
        new = {}
        for k in segs:
            L = _cut(trs[k], t0[k] + shift[k] + win[0] - m / fs, n + 2 * m)
            if L is None:
                continue
            j, v = _peak(_ncc_valid(L, stack))
            new[k] = (shift[k] + (j - m) / fs, v)
        if not new:
            break
        ref = [new[k][0] for k in new if weights.get(k, 0) > 0] or [new[k][0] for k in new]
        med = float(np.median(ref))
        move = 0.0
        for k, (sh, v) in new.items():
            sh -= med
            if abs(sh) > maxlag:
                sh = shift[k]
                v = min(v, cc_min - 1e-3)
            move = max(move, abs(sh - shift[k]))
            shift[k] = sh
            cc[k] = v
        weights = dict((k, max(cc[k], 0.0) if np.isfinite(cc[k]) and cc[k] >= cc_min else 0.0) for k in keys)
        if sum(1 for x in weights.values() if x > 0) < 2:
            weights = dict((k, 1.0) for k in keys)
        if move < 0.5 / fs:
            break
    return shift, cc, stack


def _mccc(trs, t, win, maxlag, pair_cc):
    """Multi-channel cross-correlation (VanDecar & Crosson 1990): relative delays + sigma per trace."""
    keys = list(trs)
    fs = trs[keys[0]].stats.sampling_rate
    n = int(round((win[1] - win[0]) * fs))
    m = int(round(maxlag * fs))
    short, long_ = {}, {}
    for k in keys:
        a = _cut(trs[k], t[k] + win[0], n)
        b = _cut(trs[k], t[k] + win[0] - m / fs, n + 2 * m)
        if a is not None and b is not None and np.std(a) > 0:
            short[k], long_[k] = a, b
    keys = [k for k in keys if k in short]
    N = len(keys)
    if N < 3:
        return dict((k, 0.0) for k in keys), dict((k, np.nan) for k in keys), keys
    rows, rhs, w, pairs = [], [], [], []
    for i in range(N):
        for j in range(i + 1, N):
            p, v = _peak(_ncc_valid(long_[keys[j]], short[keys[i]]))
            if v < pair_cc:
                continue
            r = np.zeros(N)
            r[i], r[j] = -1.0, 1.0
            rows.append(r)
            rhs.append((p - m) / fs)
            w.append(v)
            pairs.append((i, j))
    if len(rows) < N:
        return dict((k, 0.0) for k in keys), dict((k, np.nan) for k in keys), keys

    def solve(rows, rhs, w):
        A = np.vstack(rows + [np.ones(N)])
        b = np.array(rhs + [0.0])
        W = np.sqrt(np.array(w + [10.0]))
        x = np.linalg.lstsq(A * W[:, None], b * W, rcond=None)[0]
        return x, b[:-1] - A[:-1].dot(x)

    x, res = solve(rows, rhs, w)
    keep = np.abs(res) <= max(3.0 * np.std(res), 2.0 / fs)
    if keep.sum() >= N and (~keep).any():
        rows = [r for r, kk in zip(rows, keep) if kk]
        rhs = [r for r, kk in zip(rhs, keep) if kk]
        w = [r for r, kk in zip(w, keep) if kk]
        pairs = [r for r, kk in zip(pairs, keep) if kk]
        x, res = solve(rows, rhs, w)
    sig = {}
    for idx, k in enumerate(keys):
        rr = [res[q] for q, (i, j) in enumerate(pairs) if i == idx or j == idx]
        s = float(math.sqrt(np.sum(np.square(rr)) / (len(rr) - 1))) if len(rr) >= 2 else np.nan
        sig[k] = max(s, 0.5 / fs) if np.isfinite(s) else s
    return dict((k, float(x[i])) for i, k in enumerate(keys)), sig, keys


def _limit(pc, end_limit):
    """Copy of the phase config with every correlation window ending before interfering phases."""
    q = dict(pc)
    for k in ('cf', 'wf'):
        w = pc[k]['win']
        q[k] = dict(pc[k], win=(w[0], max(pc['min_end'][k], min(w[1], end_limit))))
    return q


def _skip_lag(cf, t, k, sk, pc):
    """Onset-function lag of trace k against the leave-one-out stack (cycle-skip / bad-onset test)."""
    fs = cf[k].stats.sampling_rate
    w = pc['cf']['win']
    n = int(round((w[1] - w[0]) * fs))
    m = int(round(pc['cf']['lag'] * fs))
    ref = []
    for j in sk:
        if j == k:
            continue
        x = _cut(cf[j], t[j] + w[0], n)
        if x is not None and np.std(x) > 0:
            ref.append((x - x.mean()) / np.std(x))
    L = _cut(cf[k], t[k] + w[0] - m / fs, n + 2 * m)
    if not ref or L is None:
        return np.nan
    p, _ = _peak(_ncc_valid(L, np.mean(ref, axis=0)))
    return (p - m) / fs


def pick_phase(trs, raw, tpred, snr, pc, end_limit=1e9):
    """trs: {sta: picking-band Trace}, raw: {sta: velocity Trace}, tpred: {sta: UTCDateTime}.
    Returns ({sta: pick dict}, stack info)."""
    pc = _limit(pc, end_limit)
    keys = sorted(trs)
    if len(keys) < 2:
        return {}, None
    # (1) coarse ICCS on onset characteristic functions (cannot skip a cycle)
    cz = dict((k, _causal(raw[k], pc['onset_band'][0], pc['onset_band'][1], pc['rate'])) for k in keys)
    cf = dict((k, _charfun(cz[k], pc['cf']['smooth'])) for k in keys)
    s1, _, _ = _iccs(cf, tpred, pc['cf']['win'], pc['cf']['lag'], pc['iters'], 0.5)
    t1 = dict((k, tpred[k] + s1[k]) for k in keys)
    # (2) picking-band waveform ICCS, lag limit below half a dominant period
    s2, cc2, _ = _iccs(trs, t1, pc['wf']['win'], pc['wf']['lag'], pc['iters'], pc['cc_stack'])
    t2 = dict((k, t1[k] + s2[k]) for k in keys)
    good = [k for k in keys if np.isfinite(cc2.get(k, np.nan)) and cc2[k] >= pc['classes'][-1][0]]
    # (3) MCCC refinement of the relative times + uncertainty per trace
    if len(good) >= 3:
        d, sig, used = _mccc(
            dict((k, trs[k]) for k in good),
            dict((k, t2[k]) for k in good),
            pc['wf']['win'],
            pc['mccc_lag'],
            pc['mccc_pair_cc'],
        )
    else:
        d, sig, used = {}, {}, []
    t3 = dict((k, t2[k] + d.get(k, 0.0)) for k in keys)
    # (4) stack of the aligned, CAUSALLY filtered traces; AIC onset = absolute time
    fs = trs[keys[0]].stats.sampling_rate
    aw = pc['aic_win']
    sk = list(used if len(used) >= 3 else good)
    if not sk:
        return {}, None
    fso = cz[sk[0]].stats.sampling_rate
    n = int(round((aw[1] - aw[0]) * fso))
    segs, sk2 = [], []
    for k in sk:
        x = _cut(cz[k], t3[k] + aw[0], n)
        if x is not None and np.std(x) > 0:
            segs.append(x / (np.max(np.abs(x)) + 1e-30))
            sk2.append(k)
    if len(segs) < 2:
        return {}, None
    sk = sk2
    segs = np.asarray(segs)
    onset, pws, lin = _stack_onset(segs, aw[0], fso, pc['aic_lookback'])
    tt = aw[0] + np.arange(n) / fso
    # uncertainty of the common onset: bootstrap over stations (same procedure, fixed seed)
    nb = int(pc.get('bootstrap', 0) or 0)
    sig_onset = np.nan
    if nb >= 10 and len(segs) >= 3:
        rng = np.random.RandomState(0)
        ob = [
            _stack_onset(segs[rng.randint(0, len(segs), len(segs))], aw[0], fso, pc['aic_lookback'])[0]
            for _ in range(nb)
        ]
        sig_onset = max(float(np.std(ob)), 1.0 / fso)
    # (5) final picks: leave-one-out stack correlation, single-trace AIC check, cycle-skip test, class
    wf = pc['wf']['win']
    nw = int(round((wf[1] - wf[0]) * fs))
    W = {}
    for k in keys:
        x = _cut(trs[k], t3[k] + wf[0], nw)
        if x is not None and np.std(x) > 0:
            W[k] = (x - x.mean()) / np.std(x)
    in_ref = [k for k in sk if k in W]
    ref_sum = np.sum([W[k] for k in in_ref], axis=0) if in_ref else None
    picks = {}
    for k in keys:
        ccf = np.nan
        if k in W and ref_sum is not None:
            nr = len(in_ref) - (1 if k in in_ref else 0)
            if nr >= 1:
                ccf = float(np.corrcoef(W[k], ref_sum - W[k] if k in in_ref else ref_sum)[0, 1])
        tp = t3[k] + onset
        y = _win(cz[k], tp + pc['single_aic'][0], tp + pc['single_aic'][1])
        d_single = np.nan
        if y is not None and len(y) > 10:
            ay = _aic(y)
            if np.isfinite(ay).any():
                d_single = pc['single_aic'][0] + int(np.argmin(ay)) / cz[k].stats.sampling_rate
        sg = sig.get(k, np.nan)
        skip = _skip_lag(cf, t3, k, sk, pc)
        cls = 'D'
        for name, (cmin, smax, nmin) in zip(CLASS_NAMES, pc['classes']):
            if (
                np.isfinite(ccf)
                and ccf >= cmin
                and np.isfinite(sg)
                and sg <= smax
                and snr.get(k, 0) >= nmin
                and np.isfinite(skip)
                and abs(skip) <= pc['skip_tol']
            ):
                cls = name
                break
        picks[k] = dict(
            t_pick=tp,
            cc=ccf,
            sigma=sg,
            snr=snr.get(k, np.nan),
            cls=cls,
            in_mccc=k in used,
            d_single_aic=d_single,
            skip_lag=skip,
            shift_from_pred=float(tp - tpred[k]),
            sigma_onset=sig_onset,
            sigma_abs=float(math.sqrt(sg**2 + sig_onset**2)),
        )
    return picks, dict(pws=pws, lin=lin, t=tt, onset=onset, n=len(segs), wf_end=wf[1], sigma_onset=sig_onset)


def _stack_onset(segs, t_first, fso, lookback):
    """Absolute onset of a stack of aligned, causally filtered segments (rows of segs).
    The phase-weighted stack's envelope locates the arrival (first time it reaches half its
    maximum); the AIC minimum of the LINEAR stack within 'lookback' s before that point is the
    onset (PWS weights are acausal, the linear stack of causal traces is not)."""
    pws, lin = _pws(segs)
    n = len(lin)
    edge = max(2, int(0.06 * n))
    tap = np.ones(n)
    tap[:edge] = np.hanning(2 * edge)[:edge]
    tap[-edge:] = np.hanning(2 * edge)[edge:]
    env_s = np.abs(hilbert(pws * tap))
    env_s[:edge] = 0.0
    env_s[-edge:] = 0.0
    i_first = int(np.argmax(env_s >= 0.5 * env_s.max()))
    i0 = max(0, i_first - int(round(lookback * fso)))
    i1 = min(n, i_first + int(round(0.5 * fso)) + 1)
    a = _aic(lin[i0:i1])
    onset = t_first + ((i0 + int(np.argmin(a))) if np.isfinite(a).any() else i_first) / fso
    return onset, pws, lin


# =============================================================================
# PER-EVENT PROCESSING (runs in a worker process)
# =============================================================================
def process_event(folder, out_dir, cfg):
    t_start = time.time()
    st, ev, hdr, pieces, gaps = read_event(folder, cfg)
    name = ev['name']
    stations = sorted(set(tr.stats.station for tr in st))
    phases = set(
        ['P', 'S']
        + list(cfg['figure_phases'])
        + list(cfg['pick']['P']['interfere'])
        + list(cfg['pick']['S']['interfere'])
    )
    tt = TravelTimes(_G['model'], ev['depth_km'], [hdr[s]['gcarc'] for s in stations], phases)
    fc = psd_freqs(cfg)
    rows, psd_rows, data = [], [], {}
    for sta in stations:
        h = hdr[sta]
        dist, baz = h['gcarc'], h['baz']
        tPr, tSr = tt.time('P', dist), tt.time('S', dist)
        if not (np.isfinite(tPr) and np.isfinite(tSr)):
            rows.append(dict(event=name, station=sta, dist=round(dist, 3), note='no ak135 P/S prediction'))
            continue
        tP, tS = ev['origin'] + tPr, ev['origin'] + tSr
        row = dict(
            event=name, station=sta, dist=round(dist, 3), baz=round(baz, 2), tP_pred=str(tP), tS_pred=str(tS)
        )
        comps, metas, vel = {}, {}, {}
        for tr in st.select(station=sta):
            c = tr.stats.channel[-1]
            comps[c] = tr
            g = gaps.get(tr.id)
            ig = integrity(tr, cfg, g[2] if g else None)
            for k, v in ig.items():
                row['%s_%s' % (k, c)] = v
            row['gap_frac_' + c] = round(g[0], 4) if g else 0.0
            row['pieces_' + c] = pieces.get(tr.id, 1)
            metas[c] = channel_meta(tr.id, tr.stats.starttime)
            if not ig['dead'] and not ig['clipped']:
                v, how = to_velocity(tr, metas[c])
                vel[c] = v
                row['units'] = 'counts' if how == 'counts' else 'm/s (%s)' % how
        mz = metas.get('Z')
        if mz:
            row.update(sensor=mz['sensor'], elev=mz['elev'], lat=mz['lat'], lon=mz['lon'], site=mz['site'])
        else:
            row.update(sensor='', elev=h['stel'], lat=h['stla'], lon=h['stlo'], site='')
        row['components'] = ''.join(sorted(comps))
        # noise PSD of every live channel (raw orientation, before trimming/rotation)
        if row.get('units', '').startswith('m/s'):
            t_end = tP - cfg['psd_gap_before_P_s']
            for c, v in vel.items():
                if _gap_in(gaps.get(v.id), t_end - cfg['psd_window_s'], t_end, v.stats.sampling_rate):
                    continue
                p = noise_psd(v, t_end, cfg, fc)
                if p is not None:
                    psd_rows.append(
                        dict(event=name, seed_id=v.id, station=sta, comp=c, psd=p.astype(np.float32))
                    )
        if 'Z' not in vel:
            row.update(usable_P=False, usable_S=False, P_visibility='n/a', S_visibility='n/a')
            rows.append(row)
            continue
        R, T = rotate_components(vel.get('Z'), vel.get('N'), vel.get('E'), metas, baz)
        z = vel['Z']
        for tag, band in (('P_sp', cfg['P_band']), ('P_bb', cfg['P_band_bb'])):
            zf = _bandpass(z, *band)
            row['snr_' + tag] = _rmsw(zf, tP + cfg['snr_sig_P'][0], tP + cfg['snr_sig_P'][1]) / (
                _rmsw(zf, tP + cfg['snr_noise_P'][0], tP + cfg['snr_noise_P'][1]) + 1e-30
            )
        row['noiseZ_sp'] = _rmsw(_bandpass(z, *cfg['P_band']), tP - 110, tP - 5)
        row['snr_S'] = np.nan
        if T is not None:
            Tf = _bandpass(T, *cfg['S_band'])
            row['snr_S'] = _rmsw(Tf, tS + cfg['snr_sig_S'][0], tS + cfg['snr_sig_S'][1]) / (
                _rmsw(Tf, tS + cfg['snr_noise_S'][0], tS + cfg['snr_noise_S'][1]) + 1e-30
            )
        for k in ('snr_P_sp', 'snr_P_bb', 'snr_S'):
            if not np.isfinite(row.get(k, np.nan)) or row[k] > 1e6:
                row[k] = np.nan
        sr = z.stats.sampling_rate
        gap_P = _gap_in(gaps.get(z.id), tP - 70, tP + 30, sr)
        gap_S = any(
            _gap_in(gaps.get(v.id), tS - 70, tS + 30, sr)
            for v in (vel.get('N'), vel.get('E'))
            if v is not None
        )
        row['gap_in_P_window'], row['gap_in_S_window'] = gap_P, bool(gap_S)
        row['usable_P'] = bool(
            np.isfinite(row['snr_P_sp'])
            and row['snr_P_sp'] >= cfg['min_snr_trace_P']
            and not gap_P
            and not row.get('flat_Z', False)
        )
        row['usable_S'] = bool(
            T is not None
            and np.isfinite(row['snr_S'])
            and row['snr_S'] >= cfg['min_snr_trace_S']
            and not gap_S
            and not row.get('flat_N', False)
            and not row.get('flat_E', False)
        )
        row['P_visibility'] = _classify_snr(row['snr_P_sp'])
        row['S_visibility'] = _classify_snr(row['snr_S'])
        data[sta] = dict(
            Z=z,
            N=vel.get('N'),
            E=vel.get('E'),
            R=R,
            T=T,
            tP=tP,
            tS=tS,
            dist=dist,
            baz=baz,
            az=h['az'],
            lat=row['lat'],
            lon=row['lon'],
            elev=row['elev'] or 0.0,
        )
        rows.append(row)
    rowd = dict((r['station'], r) for r in rows)

    # ---------------- first-order event screening ----------------
    n_alive = len(data)
    n_P = sum(
        1 for r in rows if np.isfinite(r.get('snr_P_sp', np.nan)) and r['snr_P_sp'] >= cfg['min_snr_detect']
    )
    n_S = sum(1 for r in rows if np.isfinite(r.get('snr_S', np.nan)) and r['snr_S'] >= cfg['min_snr_detect'])
    evrow = dict(
        event=name,
        origin=str(ev['origin']),
        lat=ev['lat'],
        lon=ev['lon'],
        depth_km=ev['depth_km'],
        mag=ev['mag'],
        n_stations=len(stations),
        n_alive_Z=n_alive,
        n_P_detected=n_P,
        n_S_detected=n_S,
    )
    status, reason = 'accepted', ''
    if n_alive < cfg['min_stations_alive']:
        status, reason = 'rejected', 'no usable data (%d live vertical channels)' % n_alive
    elif n_P < cfg['min_stations_P']:
        status, reason = 'rejected', 'no P detected (SNR>=%.0f at %d stations)' % (cfg['min_snr_detect'], n_P)

    # ---------------- picking ----------------
    picks, stacks = [], {}
    if status == 'accepted':
        for phase in ('P', 'S'):
            pc = cfg['pick'][phase]
            if phase == 'S' and n_S < cfg['min_stations_S']:
                evrow['S_note'] = 'S not picked (detected at %d stations)' % n_S
                continue
            use = [
                s
                for s in data
                if data[s][pc['comp']] is not None
                and rowd[s]['usable_' + phase]
                and data[s]['dist'] <= pc['max_dist']
            ]
            if phase == 'S' and len(use) < cfg['min_stations_S']:
                evrow['S_note'] = 'S not picked (%d usable stations within %.0f deg)' % (
                    len(use),
                    pc['max_dist'],
                )
                continue
            if len(use) < 2:
                continue
            trs = dict(
                (s, _bandpass(data[s][pc['comp']], pc['band'][0], pc['band'][1], pc['rate'])) for s in use
            )
            raw = dict((s, data[s][pc['comp']]) for s in use)
            tpred = dict((s, data[s]['t' + phase]) for s in use)
            snr = dict((s, rowd[s]['snr_P_sp' if phase == 'P' else 'snr_S']) for s in use)
            lags = [
                tt.time(ip, data[s]['dist']) - tt.time(phase, data[s]['dist'])
                for s in use
                for ip in pc['interfere']
            ]
            lags = [g for g in lags if np.isfinite(g) and g > -2.0]
            end_limit = (min(lags) - pc['margin']) if lags else 1e9
            evrow['%s_window_end_s' % phase] = round(_limit(pc, end_limit)['wf']['win'][1], 2)
            pk, info = pick_phase(trs, raw, tpred, snr, pc, end_limit)
            stacks[phase] = info
            vsurf = cfg['vp_surface'] if phase == 'P' else cfg['vs_surface']
            new = []
            for s, p in pk.items():
                dist = data[s]['dist']
                tpred_rel = tt.time(phase, dist)
                pkm = tt.rayp_km(phase, dist)
                elev_km = max(0.0, (data[s]['elev'] or 0.0) / 1000.0)
                t_elev = elev_km * math.sqrt(max(1.0 / vsurf**2 - pkm**2, 0.0)) if np.isfinite(pkm) else 0.0
                obs = p['t_pick'] - ev['origin']
                fin = lambda v, nd: round(v, nd) if np.isfinite(v) else np.nan
                new.append(
                    dict(
                        event=name,
                        origin=str(ev['origin']),
                        ev_lat=ev['lat'],
                        ev_lon=ev['lon'],
                        ev_depth_km=ev['depth_km'],
                        mag=ev['mag'],
                        network=data[s]['Z'].stats.network,
                        station=s,
                        sta_lat=data[s]['lat'],
                        sta_lon=data[s]['lon'],
                        sta_elev_m=data[s]['elev'],
                        phase=phase,
                        component=pc['comp'],
                        band='%g-%g Hz' % pc['band'],
                        dist_deg=round(dist, 4),
                        baz=round(data[s]['baz'], 2),
                        pick_time=str(p['t_pick']),
                        tt_obs=round(obs, 4),
                        tt_pred=round(tpred_rel, 4),
                        model=cfg['model'],
                        res_abs=round(obs - tpred_rel, 4),
                        elev_corr=round(t_elev, 4),
                        res_abs_elev=round(obs - tpred_rel - t_elev, 4),
                        sigma_mccc=fin(p['sigma'], 4),
                        sigma_onset=fin(p['sigma_onset'], 4),
                        sigma_abs=fin(p['sigma_abs'], 4),
                        cc_stack=fin(p['cc'], 3),
                        snr=fin(p['snr'], 2),
                        quality=p['cls'],
                        single_aic_minus_pick=fin(p['d_single_aic'], 3),
                        onset_lag=fin(p['skip_lag'], 3),
                        rayp_s_per_deg=fin(pkm * 111.195, 4),
                    )
                )
            good = [q for q in new if q['quality'] in CLASS_NAMES]
            if phase == 'S' and len(good) < cfg['min_picks_S']:
                for q in good:
                    q['quality'] = 'D'
                good = []
                evrow['S_note'] = 'S not coherent (fewer than %d A-C picks)' % cfg['min_picks_S']
            mu = float(np.mean([q['res_abs_elev'] for q in good])) if good else np.nan
            for q in new:
                q['res_rel'] = (
                    round(q['res_abs_elev'] - mu, 4)
                    if (np.isfinite(mu) and q['quality'] in CLASS_NAMES)
                    else np.nan
                )
            picks += new
            evrow['n_%s_picks_ABC' % phase] = len(good)
            evrow['n_%s_picks_A' % phase] = sum(1 for q in good if q['quality'] == 'A')
            evrow['mean_%s_res_abs' % phase] = round(mu, 3) if np.isfinite(mu) else np.nan
        if evrow.get('n_P_picks_ABC', 0) < cfg['min_picks']:
            status, reason = 'rejected', 'P not coherent across the array (%d A-C picks)' % evrow.get(
                'n_P_picks_ABC', 0
            )
    evrow['status'], evrow['reason'] = status, reason

    # ---------------- outputs ----------------
    odir = os.path.join(out_dir, 'processed', name)
    if os.path.isdir(odir):  # results of an earlier run of this event (always inside --out)
        shutil.rmtree(odir)
    if status == 'accepted' and cfg['write_waveforms']:
        write_processed(out_dir, ev, data, rowd, picks, cfg)
    if cfg['figures'] and (status == 'accepted' or cfg['plot_rejected']):
        try:
            plot_event(
                os.path.join(out_dir, 'figures', 'qc_%s.png' % name),
                ev,
                evrow,
                data,
                rows,
                picks,
                stacks,
                tt,
                cfg,
            )
        except Exception:
            evrow['figure_error'] = traceback.format_exc(limit=1)
    evrow['seconds'] = round(time.time() - t_start, 1)
    return dict(event=evrow, stations=rows, picks=picks, psd=psd_rows)


# =============================================================================
# PROCESSED WAVEFORMS (time zero = earthquake origin by default)
# =============================================================================
def write_processed(out_dir, ev, data, rowd, picks, cfg):
    from obspy.core.util import AttribDict

    odir = os.path.join(out_dir, 'processed', ev['name'])
    if not os.path.isdir(odir):
        os.makedirs(odir)
    pk = dict(((q['station'], q['phase']), q) for q in picks)
    for sta, d in data.items():
        if cfg['write_only_usable'] and not rowd[sta].get('usable_P'):
            continue
        proc = [
            (c, _bandpass(d[c], cfg['out_band'][0], cfg['out_band'][1], cfg['out_rate']))
            for c in cfg['out_components']
            if d.get(c) is not None
        ]
        if not proc:
            continue
        fac = 1.0
        if cfg['normalize'] == 'station':
            z = dict(proc).get('Z', proc[0][1])
            w = _win(z, d['tP'] - 5, d['tP'] + 60)
            fac = float(np.max(np.abs(w))) if w is not None and len(w) else float(np.max(np.abs(z.data)))
        qP, qS = pk.get((sta, 'P')), pk.get((sta, 'S'))
        for c, t in proc:
            if cfg['normalize'] == 'trace':
                fac = float(np.max(np.abs(t.data)))
            if cfg['normalize'] in ('station', 'trace') and fac > 0:
                t.data = t.data / fac
            if cfg['time_reference'] == 'start':
                ref = t.stats.starttime
            elif cfg['time_reference'] == 'P':
                ref = _utc(qP['pick_time']) if (qP and qP['quality'] in CLASS_NAMES) else d['tP']
            else:
                ref = ev['origin']
            hdr = dict(
                evla=ev['lat'],
                evlo=ev['lon'],
                evdp=ev['depth_km'],
                mag=ev['mag'],
                stla=d['lat'],
                stlo=d['lon'],
                stel=d['elev'],
                gcarc=d['dist'],
                baz=d['baz'],
                az=d['az'],
                o=float(ev['origin'] - ref),
                iztype=11 if cfg['time_reference'] == 'origin' else 9,
                t8=float(d['tP'] - ref),
                kt8='P_' + cfg['model'][:5],
                t9=float(d['tS'] - ref),
                kt9='S_' + cfg['model'][:5],
                user0=fac if cfg['normalize'] else 1.0,
                kuser0='NORMFAC',
                idep=7,
                lcalda=0,
            )
            if qP and qP['quality'] in CLASS_NAMES:
                hdr.update(
                    t0=float(_utc(qP['pick_time']) - ref),
                    kt0='P_' + qP['quality'],
                    user1=qP['cc_stack'],
                    user2=qP['sigma_mccc'],
                )
            if qS and qS['quality'] in CLASS_NAMES:
                hdr.update(
                    t1=float(_utc(qS['pick_time']) - ref),
                    kt1='S_' + qS['quality'],
                    user3=qS['cc_stack'],
                    user4=qS['sigma_mccc'],
                )
            hdr['cmpaz'], hdr['cmpinc'] = {
                'Z': (0.0, 0.0),
                'N': (0.0, 90.0),
                'E': (90.0, 90.0),
                'R': ((d['baz'] + 180.0) % 360.0, 90.0),
                'T': ((d['baz'] - 90.0) % 360.0, 90.0),
            }[c]
            for k, v in list(hdr.items()):
                if isinstance(v, float) and not np.isfinite(v):
                    hdr[k] = -12345.0
            hdr.update(
                nzyear=ref.year,
                nzjday=ref.julday,
                nzhour=ref.hour,
                nzmin=ref.minute,
                nzsec=ref.second,
                nzmsec=int(round(ref.microsecond / 1000.0)) if ref.microsecond < 999500 else 999,
            )
            t.stats.sac = AttribDict(hdr)
            t.stats.channel = t.stats.channel[:2] + c
            t.data = t.data.astype(np.float32)
            t.write(
                os.path.join(odir, '%s.%s..%s.sac' % (t.stats.network, sta, t.stats.channel)), format='SAC'
            )


# =============================================================================
# EVENT FIGURE
# =============================================================================
QCOL = {'A': '#1a7f37', 'B': '#2E7D7A', 'C': '#C26A1B', 'D': '#B42318'}
VCOL = {'clear': '#1a7f37', 'visible': '#2E7D7A', 'weak': '#C26A1B', 'none': '#B42318', 'n/a': '#999999'}


def plot_event(fname, ev, evrow, data, rows, picks, stacks, tt, cfg):
    import matplotlib.pyplot as plt

    rowd = dict((r['station'], r) for r in rows)
    order = sorted(data, key=lambda s: data[s]['dist'])
    pk = dict(((q['station'], q['phase']), q) for q in picks)
    nst = max(len(order), 1)
    fig, axs = plt.subplots(
        1, 5, figsize=(24, max(6.0, 0.24 * nst + 2.6)), gridspec_kw=dict(width_ratios=[1.4, 1, 1, 1, 1])
    )
    yl = (-0.9, nst - 0.1)
    # (a) full record, time after origin
    ax = axs[0]
    for i, s in enumerate(order):
        tr = _bandpass(data[s]['Z'], 0.02, 1.0, 20.0)
        t = tr.times() + (tr.stats.starttime - ev['origin'])
        ax.plot(t, tr.data / (np.max(np.abs(tr.data)) + 1e-30) * 0.45 + i, color='k', lw=0.3)
        for ph, col in (('P', 'C0'), ('S', 'C3')):
            q = pk.get((s, ph))
            if q and q['quality'] in CLASS_NAMES:
                ax.plot([q['tt_obs']] * 2, [i - 0.42, i + 0.42], color=col, lw=1.4)
    dd = np.array([data[s]['dist'] for s in order])
    phs = cfg['figure_phases'] if ev['depth_km'] > 70 else [p for p in cfg['figure_phases'] if p[0] in 'PS']
    for j, ph in enumerate(phs):
        y = np.array([tt.time(ph, x) for x in dd])
        if not np.isfinite(y).any():
            continue
        col = 'C%d' % ((j + 2) % 10)
        ax.plot(y, np.arange(len(dd)), lw=0.9, alpha=0.55, color=col)
        ax.text(
            np.nanmean(y),
            nst - 0.05 + 0.6 * (j % 2),
            ph,
            fontsize=7,
            ha='center',
            va='bottom',
            color=col,
            clip_on=False,
        )
    ax.set_yticks(range(len(order)))
    ax.set_yticklabels(['%s %.1f°' % (s, data[s]['dist']) for s in order], fontsize=6)
    ax.set_ylim(*yl)
    ax.set_xlabel('time after origin (s)')
    ax.set_title('vertical, 0.02-1 Hz, %s phases; picks: blue P, red S' % cfg['model'], fontsize=9, pad=24)
    # (b) P and (d) S around the prediction
    for ax, ph, comp, band, w in (
        (axs[1], 'P', 'Z', cfg['P_band'], (-15.0, 30.0)),
        (axs[3], 'S', 'T', cfg['S_band'], (-60.0, 90.0)),
    ):
        for i, s in enumerate(order):
            tr0 = data[s][comp]
            if tr0 is None:
                ax.text(0, i, 'no horizontals', fontsize=6, ha='center', va='center', color=VCOL['none'])
                continue
            tr = _bandpass(tr0, band[0], band[1], cfg['pick'][ph]['rate'])
            t0 = data[s]['t' + ph]
            x = _win(tr, t0 + w[0], t0 + w[1])
            if x is None:
                continue
            q = pk.get((s, ph))
            col = QCOL[q['quality']] if q else '#9A9A9A'
            ax.plot(
                w[0] + np.arange(len(x)) / tr.stats.sampling_rate,
                x / (np.max(np.abs(x)) + 1e-30) * 0.45 + i,
                color=col,
                lw=0.5,
            )
            if q:
                tq = q['tt_obs'] - (t0 - ev['origin'])
                ax.plot([tq, tq], [i - 0.45, i + 0.45], color='k', lw=1.2)
            v = rowd[s].get('snr_P_sp' if ph == 'P' else 'snr_S', np.nan)
            vis = rowd[s].get('P_visibility' if ph == 'P' else 'S_visibility', 'n/a')
            ax.text(
                w[1] + 0.02 * (w[1] - w[0]),
                i,
                ('%.1f' % v) if np.isfinite(v) else '-',
                fontsize=6,
                va='center',
                color=VCOL.get(vis, '#999999'),
            )
        ax.axvline(0, color='r', lw=0.8, ls='--')
        ax.set_xlim(w[0], w[1] + 0.12 * (w[1] - w[0]))
        ax.set_ylim(*yl)
        ax.set_yticks([])
        ax.set_xlabel('time relative to %s %s (s)' % (cfg['model'], ph))
        ax.set_title(
            '%s on %s, %.2g-%.2g Hz: color = class, | = pick, number = SNR'
            % (ph, 'vertical' if comp == 'Z' else 'transverse', band[0], band[1]),
            fontsize=8,
        )
    # (c) P and (e) S aligned on the picks, causal onset band, with the stack
    for ax, ph, w in ((axs[2], 'P', (-5.0, 12.0)), (axs[4], 'S', (-30.0, 60.0))):
        pc = cfg['pick'][ph]
        for i, s in enumerate(order):
            q = pk.get((s, ph))
            tr0 = data[s][pc['comp']]
            if not q or tr0 is None:
                continue
            tr = _causal(tr0, pc['onset_band'][0], pc['onset_band'][1], pc['rate'])
            x = _win(tr, _utc(q['pick_time']) + w[0], _utc(q['pick_time']) + w[1])
            if x is None:
                continue
            ax.plot(
                w[0] + np.arange(len(x)) / tr.stats.sampling_rate,
                x / (np.max(np.abs(x)) + 1e-30) * 0.45 + i,
                color=QCOL[q['quality']],
                lw=0.6,
            )
            lab = ('%s %.2f' % (q['quality'], q['cc_stack'])) if np.isfinite(q['cc_stack']) else q['quality']
            ax.text(w[1] + 0.01 * (w[1] - w[0]), i, lab, fontsize=6, va='center')
        info = stacks.get(ph)
        if info:
            x = info['lin'] / (np.max(np.abs(info['lin'])) + 1e-30)
            ax.plot(info['t'] - info['onset'], x * 0.9 + nst + 0.6, color='k', lw=1.2)
            ax.text(
                w[0] + 0.02 * (w[1] - w[0]), nst + 1.4, 'stack of %d aligned traces' % info['n'], fontsize=7
            )
        ax.axvline(0, color='k', lw=0.8)
        ax.set_xlim(w[0], w[1] + 0.14 * (w[1] - w[0]))
        ax.set_ylim(yl[0], nst + 1.9)
        ax.set_yticks([])
        ax.set_xlabel('time relative to %s pick (s)' % ph)
        ax.set_title(
            '%s aligned on picks, causal %.2g-%.2g Hz (class, CC)'
            % (ph, pc['onset_band'][0], pc['onset_band'][1]),
            fontsize=8,
        )
    stt = evrow['status'].upper() + ('' if evrow['status'] == 'accepted' else ': ' + evrow['reason'])
    fig.suptitle(
        '%s   %s UTC   M%.1f   %.0f km   (%.2f, %.2f)   P picks A-C: %s   S picks A-C: %s   [%s]'
        % (
            ev['name'],
            ev['origin'].strftime('%Y-%m-%d %H:%M:%S'),
            ev['mag'],
            ev['depth_km'],
            ev['lat'],
            ev['lon'],
            evrow.get('n_P_picks_ABC', 0),
            evrow.get('n_S_picks_ABC', 0),
            stt,
        ),
        fontsize=10,
    )
    fig.text(
        0.5,
        0.005,
        'Pick classes: A green, B teal, C orange (used); D red (rejected); gray = not picked. '
        'SNR colors: >=5 clear, 3-5 visible, 2-3 weak, <2 none.',
        ha='center',
        fontsize=8,
    )
    plt.tight_layout(rect=(0, 0.015, 1, 0.97))
    d = os.path.dirname(fname)
    if not os.path.isdir(d):
        os.makedirs(d)
    plt.savefig(fname, dpi=110)
    plt.close(fig)


# =============================================================================
# AGGREGATION: tables, noise PDFs, station performance, pick summary
# =============================================================================
# columns of the two pick files: (name in picks_all.csv, name in the pick file)
_ID_COLS = [
    ('event', 'event'),
    ('origin', 'origin_time'),
    ('ev_lat', 'ev_lat'),
    ('ev_lon', 'ev_lon'),
    ('ev_depth_km', 'ev_depth_km'),
    ('mag', 'mag'),
    ('network', 'network'),
    ('station', 'station'),
    ('sta_lat', 'sta_lat'),
    ('sta_lon', 'sta_lon'),
    ('sta_elev_m', 'sta_elev_m'),
    ('phase', 'phase'),
    ('quality', 'quality'),
]
ABS_COLS = _ID_COLS + [
    ('pick_time', 'pick_time'),
    ('tt_obs', 'tt_obs_s'),
    ('tt_pred', 'tt_pred_s'),
    ('res_abs', 'res_s'),
    ('elev_corr', 'elev_corr_s'),
    ('res_abs_elev', 'res_elev_s'),
    ('sigma_abs', 'sigma_abs_s'),
    ('sigma_onset', 'sigma_onset_s'),
    ('sigma_mccc', 'sigma_mccc_s'),
    ('cc_stack', 'cc'),
    ('snr', 'snr'),
    ('dist_deg', 'dist_deg'),
    ('baz', 'baz_deg'),
    ('rayp_s_per_deg', 'rayp_s_per_deg'),
    ('model', 'model'),
    ('component', 'component'),
    ('band', 'band'),
]
REL_COLS = _ID_COLS + [
    ('dt_rel', 'dt_rel_s'),
    ('dt_rel_noelev', 'dt_rel_noelev_s'),
    ('sigma_mccc', 'sigma_mccc_s'),
    ('cc_stack', 'cc'),
    ('snr', 'snr'),
    ('n_event', 'n_picks_event'),
    ('dist_deg', 'dist_deg'),
    ('baz', 'baz_deg'),
    ('rayp_s_per_deg', 'rayp_s_per_deg'),
    ('tt_pred', 'tt_pred_s'),
    ('elev_corr', 'elev_corr_s'),
    ('model', 'model'),
    ('component', 'component'),
    ('band', 'band'),
]


def write_pick_files(pdir, good):
    """Two pick files, accepted events and classes A-C only, one row per station and phase:
      picks_absolute.csv  absolute travel times and residuals  -> relocation, absolute tomography
      picks_relative.csv  relative delays, event mean removed  -> relative travel-time tomography
    dt_rel_s is demeaned over exactly the rows of the same event and phase in this file."""
    if not os.path.isdir(pdir):
        os.makedirs(pdir)
    g = good.copy()
    for c in ('network', 'sigma_onset', 'sigma_abs'):
        if c not in g:  # work files written by v1.0
            g[c] = '' if c == 'network' else np.nan
    g['_ph'] = g['phase'].map({'P': 0, 'S': 1})
    g = g.sort_values(['origin', '_ph', 'station'])
    grp = g.groupby(['event', 'phase'])
    g['n_event'] = grp['station'].transform('size')
    g['dt_rel'] = (g['res_abs_elev'] - grp['res_abs_elev'].transform('mean')).round(4)
    g['dt_rel_noelev'] = (g['res_abs'] - grp['res_abs'].transform('mean')).round(4)
    g = g.round(
        {'ev_lat': 5, 'ev_lon': 5, 'sta_lat': 5, 'sta_lon': 5, 'ev_depth_km': 3, 'mag': 2}
    )  # SAC float32
    for fn, cols in (('picks_absolute.csv', ABS_COLS), ('picks_relative.csv', REL_COLS)):
        out = g[[a for a, _ in cols]].copy()
        out.columns = [b for _, b in cols]
        out.to_csv(os.path.join(pdir, fn), index=False)
    return g


def aggregate(out_dir, cfg, inv_info):
    import pandas as pd
    from obspy.signal.spectral_estimation import get_nlnm, get_nhnm

    work = sorted(glob.glob(os.path.join(out_dir, 'work', '*.pkl')))
    if not work:
        print('nothing to aggregate')
        return
    E, S, P, Q = [], [], [], []
    old = 0
    for f in work:
        r = _load_work(f)
        if r is None:
            E.append(
                dict(
                    event=os.path.basename(f)[:-4],
                    status='error',
                    reason='unreadable work file (run again without --aggregate-only)',
                )
            )
            continue
        old += r.get('version') != VERSION
        if 'error' in r:
            E.append(dict(event=r['event'], status='error', reason=r['error'].strip().splitlines()[-1]))
            continue
        E.append(r['event'])
        S += r['stations']
        P += r['picks']
        Q += r['psd']
    if old:
        print(
            'WARNING: %d events were processed by an older version; run again without --aggregate-only '
            '(and without --events) to update them' % old
        )
    obsolete = os.path.join(out_dir, 'tables', 'picks_ABC.csv')  # v1.0 output, replaced by picks/
    if os.path.exists(obsolete):
        os.remove(obsolete)
    ev, sq, pk = pd.DataFrame(E), pd.DataFrame(S), pd.DataFrame(P)
    for c in ('snr_P_sp', 'snr_S', 'noiseZ_sp'):
        if c not in sq:
            sq[c] = np.nan
    tdir = os.path.join(out_dir, 'tables')
    if not os.path.isdir(tdir):
        os.makedirs(tdir)
    # second-pass screening: pre-event noise far above the station's own norm (overlapping earthquake)
    med = sq.groupby('station')['noiseZ_sp'].median()
    sq['noise_ratio'] = sq['noiseZ_sp'] / sq['station'].map(med)
    ev['pre_event_noise_ratio'] = ev['event'].map(sq.groupby('event')['noise_ratio'].median()).round(2)
    ev['flag_noisy_pre_event'] = ev['pre_event_noise_ratio'] > cfg['noise_contamination_ratio']
    ev.to_csv(os.path.join(tdir, 'event_summary.csv'), index=False)
    sq.to_csv(os.path.join(tdir, 'station_event_qc.csv'), index=False)
    good = pd.DataFrame()
    if len(pk):
        pk = pk.merge(ev[['event', 'status']], on='event', how='left')
        pk.to_csv(os.path.join(tdir, 'picks_all.csv'), index=False)  # diagnostic: every pick, A-D
        good = pk[(pk.status == 'accepted') & (pk.quality.isin(CLASS_NAMES))].copy()
        write_pick_files(os.path.join(out_dir, 'picks'), good)
    # --- noise PSD probability density functions ---
    fc = psd_freqs(cfg)
    lo, hi, step = cfg['psd_db_bins']
    edges = np.arange(lo, hi + step, step)
    ndir = os.path.join(out_dir, 'noise')
    if not os.path.isdir(ndir):
        os.makedirs(ndir)
    nlT, nl = get_nlnm()
    nhT, nh = get_nhnm()
    lf = np.log10(1.0 / fc)
    nlnm = np.interp(lf, np.log10(nlT), nl)
    nhnm = np.interp(lf, np.log10(nhT), nh)
    by = {}
    for q in Q:
        by.setdefault(q['seed_id'], []).append(q['psd'])
    stats_rows, modes = [], {}
    bands = {
        'HF 5-20 Hz': (5, 20),
        'SP 1-5 Hz': (1, 5),
        'MS 0.1-0.4 Hz': (0.1, 0.4),
        'LP 0.05-0.1 Hz': (0.05, 0.1),
    }
    for sid, arr in sorted(by.items()):
        A = np.vstack(arr).astype(float)
        H = np.zeros((len(fc), len(edges) - 1))
        for i in range(len(fc)):
            col = A[:, i]
            col = col[np.isfinite(col)]
            if col.size:
                H[i] = np.histogram(col, bins=edges)[0] / float(col.size) * 100.0
        mode = edges[:-1][np.argmax(H, axis=1)] + step / 2.0
        mode[H.sum(axis=1) == 0] = np.nan
        medp = np.nanmedian(A, axis=0)
        modes[sid] = dict(
            H=H,
            mode=mode,
            median=medp,
            p10=np.nanpercentile(A, 10, axis=0),
            p90=np.nanpercentile(A, 90, axis=0),
            n=len(A),
        )
        row = dict(seed_id=sid, station=sid.split('.')[1], channel=sid.split('.')[3], n_windows=len(A))
        for bn, (f1, f2) in bands.items():
            sel = (fc >= f1) & (fc <= f2)
            row['median_dB ' + bn] = round(float(np.nanmean(medp[sel])), 1)
            row['above_NLNM_dB ' + bn] = round(float(np.nanmean(medp[sel] - nlnm[sel])), 1)
        sel = (fc >= 1) & (fc <= 10)
        row['pct_above_NHNM 1-10 Hz'] = round(float(np.nanmean(A[:, sel] > nhnm[sel]) * 100.0), 1)
        m = (inv_info.get(sid) or [{}])[0]
        row.update(sensor=m.get('sensor', ''), elev_m=m.get('elev', np.nan), site=m.get('site', ''))
        stats_rows.append(row)
    ns = pd.DataFrame(stats_rows)
    if len(ns):
        ns.to_csv(os.path.join(ndir, 'channel_noise_summary.csv'), index=False)
        _plot_station_pdfs(ndir, modes, fc, edges, nlT, nl, nhT, nh)
        _plot_network_noise(ndir, modes, ns, fc, nlT, nl, nhT, nh)
    # --- station performance table ---
    thr = cfg['min_snr_detect']
    perf = (
        sq.groupby('station')
        .agg(
            n_events=('event', 'nunique'),
            P_detected=('snr_P_sp', lambda x: int((x >= thr).sum())),
            S_detected=('snr_S', lambda x: int((x >= thr).sum())),
            median_snr_P=('snr_P_sp', 'median'),
            median_snr_S=('snr_S', 'median'),
        )
        .reset_index()
    )
    for c in ('sensor', 'elev', 'site'):
        if c in sq:
            perf[c] = perf['station'].map(sq.groupby('station')[c].first())
    if len(good):
        g = good.groupby(['station', 'phase']).size().unstack(fill_value=0)
        for ph in ('P', 'S'):
            perf['%s_picks_ABC' % ph] = perf['station'].map(g[ph] if ph in g else {}).fillna(0).astype(int)
        gz = good[good.phase == 'P'].groupby('station')['res_rel'].agg(['mean', 'std'])
        perf['P_station_term_s'] = perf['station'].map(gz['mean']).round(3)
        perf['P_res_rel_std_s'] = perf['station'].map(gz['std']).round(3)
        perf['P_pick_yield_pct'] = (100.0 * perf['P_picks_ABC'] / perf['n_events']).round(1)
    for c in (
        'dead_Z',
        'dead_N',
        'dead_E',
        'clipped_Z',
        'clipped_N',
        'clipped_E',
        'flat_Z',
        'flat_N',
        'flat_E',
    ):
        if c in sq:
            perf['n_' + c] = (
                perf['station']
                .map(sq.groupby('station')[c].apply(lambda x: int(x.fillna(False).astype(bool).sum())))
                .fillna(0)
                .astype(int)
            )
    if len(ns):
        z = ns[ns.channel.str.endswith('Z')].set_index('station')
        perf['Z_SP_above_NLNM_dB'] = perf['station'].map(z['above_NLNM_dB SP 1-5 Hz'])
        perf['Z_MS_above_NLNM_dB'] = perf['station'].map(z['above_NLNM_dB MS 0.1-0.4 Hz'])
        h = ns[~ns.channel.str.endswith('Z')].groupby('station')['above_NLNM_dB LP 0.05-0.1 Hz'].mean()
        perf['H_LP_above_NLNM_dB'] = perf['station'].map(h).round(1)
    perf.round(3).to_csv(os.path.join(tdir, 'station_performance.csv'), index=False)
    if len(good):
        _plot_pick_summary(os.path.join(out_dir, 'picks_summary.png'), good)
    n_acc = int((ev.status == 'accepted').sum())
    print(
        'events: %d accepted, %d rejected, %d errors'
        % (n_acc, int((ev.status == 'rejected').sum()), int((ev.status == 'error').sum()))
    )
    if (ev.status == 'rejected').any():
        print('rejection reasons:')
        for r_, n_ in ev[ev.status == 'rejected'].reason.str.split('(').str[0].value_counts().items():
            print('   %4d  %s' % (n_, r_.strip()))
    if len(good):
        print(
            'picks A-C: P %d, S %d  ->  %s/picks_absolute.csv and picks_relative.csv'
            % ((good.phase == 'P').sum(), (good.phase == 'S').sum(), os.path.join(out_dir, 'picks'))
        )
    print('tables in', tdir)


def _edges_log(fc):
    lf = np.log10(fc)
    d = np.diff(lf).mean()
    return 10 ** np.concatenate([[lf[0] - d / 2], (lf[:-1] + lf[1:]) / 2, [lf[-1] + d / 2]])


def _plot_station_pdfs(ndir, modes, fc, edges, nlT, nl, nhT, nh):
    import matplotlib.pyplot as plt

    try:
        from obspy.imaging.cm import pqlx as cmap
    except Exception:
        cmap = plt.get_cmap('jet')
    lab = {'Z': 'UD', 'E': 'EW', 'N': 'NS'}
    xe = _edges_log(fc)
    for sta in sorted(set(s.split('.')[1] for s in modes)):
        fig, axs = plt.subplots(3, 1, figsize=(7.2, 10.5), sharex=True)
        im = None
        for ax, c in zip(axs, 'ZEN'):
            sid = [s for s in modes if s.split('.')[1] == sta and s.endswith(c)]
            ax.plot(1.0 / nlT, nl, 'k', lw=1.0)
            ax.plot(1.0 / nhT, nh, 'k', lw=1.0)
            ax.text(0.06, -96, 'NHNM', fontsize=8)
            ax.text(0.5, -178, 'NLNM', fontsize=8)
            if sid:
                M = modes[sid[0]]
                Hm = np.ma.masked_where(M['H'].T <= 0, M['H'].T)
                im = ax.pcolormesh(xe, edges, Hm, cmap=cmap, vmin=0, vmax=20)
                ax.plot(fc, M['mode'], color='r', lw=1.4, label='Mode')
                ax.plot(fc, M['p10'], color='0.35', lw=0.7, ls='--')
                ax.plot(fc, M['p90'], color='0.35', lw=0.7, ls='--', label='10th / 90th percentile')
                ax.set_title(
                    'Seismic noise (%s)   %s   n = %d windows' % (lab[c], sid[0], M['n']), fontsize=9
                )
                ax.legend(loc='lower left', fontsize=7)
            else:
                ax.set_title('Seismic noise (%s)   %s   no data' % (lab[c], sta), fontsize=9)
            ax.set_xscale('log')
            ax.set_xlim(fc[0], fc[-1])
            ax.set_ylim(-200, -60)
            ax.set_ylabel('PSD (rel. to m/s²) in dB')
            ax.grid(alpha=0.2, which='both')
        axs[-1].set_xlabel('Frequency [Hz]')
        if im is not None:
            cb = fig.colorbar(im, ax=list(axs), shrink=0.35, pad=0.02)
            cb.set_label('Percent')
        fig.savefig(os.path.join(ndir, 'PSD_PDF_%s.png' % sta), dpi=110)
        plt.close(fig)


def _plot_network_noise(ndir, modes, ns, fc, nlT, nl, nhT, nh):
    import matplotlib.pyplot as plt

    fig, axs = plt.subplots(1, 4, figsize=(19, 4.8), gridspec_kw=dict(width_ratios=[1, 1, 1, 1.25]))
    lab = {'Z': 'vertical', 'E': 'east', 'N': 'north'}
    for ax, c in zip(axs[:3], 'ZEN'):
        sids = sorted(s for s in modes if s.endswith(c))
        lv = np.array([np.nanmean(modes[s]['median'][(fc >= 1) & (fc <= 5)]) for s in sids])
        if len(lv):
            norm = plt.Normalize(np.nanmin(lv), np.nanmax(lv))
            for s, v in zip(sids, lv):
                ax.plot(fc, modes[s]['mode'], color=plt.cm.viridis(norm(v)), lw=0.8, alpha=0.85)
        ax.plot(1.0 / nlT, nl, 'k', lw=1.2)
        ax.plot(1.0 / nhT, nh, 'k', lw=1.2)
        ax.text(0.06, -96, 'NHNM', fontsize=8)
        ax.text(0.5, -178, 'NLNM', fontsize=8)
        ax.set_xscale('log')
        ax.set_xlim(fc[0], fc[-1])
        ax.set_ylim(-200, -60)
        ax.grid(alpha=0.25, which='both')
        ax.set_title('Mode PSD per station, %s (yellow = noisiest at 1-5 Hz)' % lab[c], fontsize=9)
        ax.set_xlabel('Frequency [Hz]')
    axs[0].set_ylabel('PSD (rel. to m/s²) in dB')
    z = ns[ns.channel.str.endswith('Z')].sort_values('above_NLNM_dB SP 1-5 Hz')
    ax = axs[3]
    v = z['above_NLNM_dB SP 1-5 Hz']
    ax.barh(range(len(z)), v, color=plt.cm.RdYlGn_r(plt.Normalize(v.min(), v.max())(v)))
    ax.set_yticks(range(len(z)))
    ax.set_yticklabels(['%s (%s)' % (s, (sen or '')[:8]) for s, sen in zip(z.station, z.sensor)], fontsize=6)
    ax.set_xlabel('median Z noise above NLNM, 1-5 Hz (dB)')
    ax.set_title('Station ranking (quietest at top)', fontsize=9)
    ax.invert_yaxis()
    ax.grid(axis='x', alpha=0.3)
    plt.tight_layout()
    fig.savefig(os.path.join(ndir, 'network_noise_overview.png'), dpi=110)
    plt.close(fig)


def _plot_pick_summary(fname, good):
    import matplotlib.pyplot as plt

    fig, axs = plt.subplots(1, 4, figsize=(19, 4.6))
    P_ = good[good.phase == 'P']
    S_ = good[good.phase == 'S']
    ax = axs[0]
    st = P_.groupby('station').agg(
        lat=('sta_lat', 'first'), lon=('sta_lon', 'first'), term=('res_rel', 'mean'), n=('res_rel', 'size')
    )
    if len(st):
        v = float(np.nanmax(np.abs(st.term))) or 0.1
        sc = ax.scatter(
            st.lon,
            st.lat,
            c=st.term,
            cmap='RdBu_r',
            vmin=-v,
            vmax=v,
            s=30 + 4 * np.sqrt(st.n),
            marker='^',
            edgecolor='k',
            lw=0.4,
        )
        for s, r in st.iterrows():
            ax.text(r.lon + 0.03, r.lat + 0.03, s[-2:], fontsize=6)
        cb = fig.colorbar(sc, ax=ax)
        cb.set_label('mean relative P residual (s), + = slow')
    ax.set_title('P station terms (elevation-corrected)', fontsize=9)
    ax.set_xlabel('longitude')
    ax.set_ylabel('latitude')
    ax = axs[1]
    for ph, df, col in (('P', P_, 'C0'), ('S', S_, 'C3')):
        x = df.sigma_mccc.dropna()
        if len(x):
            ax.hist(
                x,
                bins=np.linspace(0, x.quantile(0.98) * 1.1 + 1e-3, 40),
                alpha=0.6,
                color=col,
                label='%s (median %.3f s)' % (ph, x.median()),
            )
    ax.set_xlabel('MCCC uncertainty (s)')
    ax.set_ylabel('picks')
    ax.legend(fontsize=7)
    ax.set_title('Relative-time uncertainty', fontsize=9)
    ax = axs[2]
    for ph, df, col in (('P', P_, 'C0'), ('S', S_, 'C3')):
        if len(df):
            ax.scatter(df.baz, df.res_rel, s=4, alpha=0.4, color=col, label=ph)
    ax.axhline(0, color='k', lw=0.6)
    ax.set_xlabel('back-azimuth (°)')
    ax.set_ylabel('relative residual (s)')
    ax.set_title('Relative residuals vs back-azimuth', fontsize=9)
    ax.legend(fontsize=7)
    ax = axs[3]
    q = good.groupby(['phase', 'quality']).size().unstack(fill_value=0)
    q.plot(kind='bar', stacked=True, ax=ax, color=[QCOL.get(c, 'gray') for c in q.columns], rot=0)
    ax.set_title('Picks by class (accepted events)', fontsize=9)
    ax.set_ylabel('picks')
    plt.tight_layout()
    fig.savefig(fname, dpi=110)
    plt.close(fig)


# =============================================================================
# DRIVER
# =============================================================================
def _load_work(pk):
    """Per-event result file, or None if missing or unreadable (e.g. a run killed while writing)."""
    try:
        with open(pk, 'rb') as fh:
            return pickle.load(fh)
    except Exception:
        return None


def _run_one(args):
    folder, out_dir, overwrite = args
    name = os.path.basename(folder.rstrip('/'))
    pk = os.path.join(out_dir, 'work', name + '.pkl')
    if not overwrite:
        r = _load_work(pk)
        # skip only complete results of this version; errors, unreadable or older results are redone
        if r is not None and r.get('version') == VERSION and 'error' not in r:
            return name, 'done earlier (skipped)'
    try:
        res = process_event(folder, out_dir, _G['cfg'])
        e = res['event']
        msg = '%s  P A-C %s  S A-C %s  (%.0f s)%s' % (
            e['status'],
            e.get('n_P_picks_ABC', 0),
            e.get('n_S_picks_ABC', 0),
            e['seconds'],
            ('  ' + e['reason']) if e['reason'] else '',
        )
    except Exception:
        res = dict(event=name, error=traceback.format_exc())
        msg = 'ERROR ' + res['error'].strip().splitlines()[-1]
    res['version'] = VERSION
    tmp = pk + '.part'
    with open(tmp, 'wb') as fh:  # write completely, then rename: never a half-written result
        pickle.dump(res, fh, protocol=pickle.HIGHEST_PROTOCOL)
    os.replace(tmp, pk)
    return name, msg


def _inside(path, folder):
    """True if path is folder itself or anything below it. Compares the folders themselves
    (device and inode) for path and each of its parents, so symbolic links, '..', and differences
    in upper/lower case on case-insensitive disks (macOS) cannot hide the overlap."""
    p, f = os.path.realpath(path), os.path.realpath(folder)
    if p == f or p.startswith(f.rstrip(os.sep) + os.sep):
        return True
    try:
        fst = os.stat(folder)
    except OSError:
        return False
    q = os.path.abspath(path)
    while True:
        if os.path.exists(q) and os.path.samestat(os.stat(q), fst):
            return True
        parent = os.path.dirname(q)
        if parent == q:
            return False
        q = parent


def _protected_folders(data):
    """The data folder and its sibling folders that also hold event folders (e.g. raw_mseed)."""
    out = [data]
    parent = os.path.dirname(os.path.abspath(data))
    for s in sorted(os.listdir(parent)):
        d = os.path.join(parent, s)
        if os.path.isdir(d) and not os.path.samefile(d, data) and glob.glob(os.path.join(d, 'Event_*')):
            out.append(d)
    return out


def _input_state(folders):
    """Size and modification time of every file in the event folders (the originals)."""
    state = {}
    for f in folders:
        for p in sorted(glob.glob(os.path.join(f, '*'))):
            s = os.stat(p)
            state[p] = (s.st_size, s.st_mtime_ns)
    return state


def add_arguments(ap):
    ap.add_argument(
        '--data', required=True, help='folder containing the Event_* directories, e.g. data/raw_sac'
    )
    ap.add_argument('--inventory', default=None, help='StationXML (overall sensitivity or full response)')
    ap.add_argument('--out', required=True, help='output folder')
    ap.add_argument('--events', default='Event_*', help='glob pattern for event folders')
    ap.add_argument('--workers', type=int, default=max(1, (os.cpu_count() or 2) - 1))
    ap.add_argument('--overwrite', action='store_true', help='re-process events already done')
    ap.add_argument('--no-figures', action='store_true', help='skip the per-event figures')
    ap.add_argument('--no-waveforms', action='store_true', help='skip the processed waveforms')
    ap.add_argument(
        '--aggregate-only', action='store_true', help='rebuild tables and plots from finished events'
    )


def run(a):
    cfg = dict(CONFIG)
    if a.no_figures:
        cfg['figures'] = False
    if a.no_waveforms:
        cfg['write_waveforms'] = False
    # The original data are only ever read. As a further safeguard the output folder must lie
    # completely outside the data folder (and the other way round).
    if not os.path.isdir(a.data):
        sys.exit('ERROR: data folder not found: %s' % a.data)
    for d in _protected_folders(a.data):
        if _inside(a.out, d) or _inside(d, a.out):
            sys.exit(
                'ERROR: --out (%s) must be a separate folder outside the data folders (%s), e.g. qc_out '
                'next to the data folder. Nothing was written.' % (a.out, d)
            )
    folders = []
    if not a.aggregate_only:
        folders = sorted(f for f in glob.glob(os.path.join(a.data, a.events)) if os.path.isdir(f))
        if not folders:
            sys.exit('ERROR: no event folders matching "%s" in %s' % (a.events, a.data))
    for d in ('work', 'figures', 'tables', 'noise', 'processed', 'picks'):
        p = os.path.join(a.out, d)
        if not os.path.isdir(p):
            os.makedirs(p)
    if a.inventory is None:
        print('WARNING: no inventory - data stay in counts, no noise PSDs')
    _init_worker(a.inventory, cfg)
    if a.inventory:
        full = any(e['has_stages'] for v in _G['inv_info'].values() for e in v)
        print(
            'inventory: %d channels, %s'
            % (
                len(_G['inv_info']),
                'full response' if full else 'overall sensitivity only -> flat-band velocity',
            )
        )
    if not a.aggregate_only:
        print(
            '%d event folders, %d workers, time zero = %s' % (len(folders), a.workers, cfg['time_reference'])
        )
        before = _input_state(folders)
        print('originals: %d files recorded (read only, never written)' % len(before))
        t0 = time.time()
        jobs = [(f, a.out, a.overwrite) for f in folders]
        if a.workers <= 1:
            for k, j in enumerate(jobs, 1):
                name, msg = _run_one(j)
                print('[%d/%d] %s: %s' % (k, len(jobs), name, msg))
                sys.stdout.flush()
        else:
            with ProcessPoolExecutor(
                max_workers=a.workers, initializer=_init_worker, initargs=(a.inventory, cfg)
            ) as ex:
                futs = [ex.submit(_run_one, j) for j in jobs]
                for k, fu in enumerate(as_completed(futs), 1):
                    name, msg = fu.result()
                    print('[%d/%d] %s: %s' % (k, len(jobs), name, msg))
                    sys.stdout.flush()
        print('processing time: %.1f min' % ((time.time() - t0) / 60.0))
        after = _input_state(folders)
        changed = sorted(p for p in set(before) | set(after) if before.get(p) != after.get(p))
        if changed:
            print(
                'WARNING: %d original files changed size or time stamp during the run, e.g. %s'
                % (len(changed), changed[0])
            )
        else:
            print(
                'originals check: all %d files in %d event folders unchanged (size and time stamp)'
                % (len(after), len(folders))
            )
    aggregate(a.out, cfg, _G['inv_info'])


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    add_arguments(ap)
    run(ap.parse_args(argv))


if __name__ == '__main__':
    main()
