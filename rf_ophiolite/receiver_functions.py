"""Receiver functions, H-k stacking and the Ps-only Moho depth.

The radial component is deconvolved by the vertical with the iterative time-domain
method of Ligorria and Ammon (1999). Single functions pass a set of acceptance rules
and a per-station screen; the accepted set of each station is stacked in H-k
(Zhu and Kanamori 2000) with bootstrap errors, and the Ps time is also converted to
depth with a fixed Vp/Vs for stations where the multiples are weak.
"""

import math
import os
import sys
import time

import numpy as np

from .common import processed_file, pyplot, save_figure, save_csv


def gauss_filter(x, dt, a):
    nfft = 2 ** int(np.ceil(np.log2(2 * len(x))))
    f = np.fft.rfftfreq(nfft, dt)
    return np.fft.irfft(np.fft.rfft(x, nfft) * np.exp(-((2 * np.pi * f) ** 2) / (4 * a**2)), nfft)[: len(x)]


def itdecon(R, Z, dt, a=2.0, shift=10.0, niter=200, minderr=0.001):
    """Iterative time-domain deconvolution (Ligorria & Ammon 1999). Returns the receiver function
    (amplitude relative to a unit spike, so the direct P is the R/Z amplitude ratio), the time axis
    (0 = direct P) and the fit of the R component in percent."""
    from scipy.signal import fftconvolve

    n = len(Z)
    Zg, Rg = gauss_filter(Z, dt, a), gauss_filter(R, dt, a)
    nsh = int(round(shift / dt))
    spikes = np.zeros(n)
    resid = Rg.copy()
    pw = float(np.sum(Rg**2)) + 1e-30
    zz = float(np.sum(Zg * Zg)) + 1e-30
    last = 1e9
    fit = 0.0
    for _ in range(niter):
        cc = fftconvolve(resid, Zg[::-1], mode='full')[n - 1 - nsh : 2 * n - 1 - nsh]
        k = int(np.argmax(np.abs(cc)))
        spikes[k] += cc[k] / zz
        resid = Rg - fftconvolve(spikes, Zg, mode='full')[nsh : nsh + n]
        fit = 100.0 * (1.0 - float(np.sum(resid**2)) / pw)
        if abs(last - fit) < minderr:
            break
        last = fit
    unit = np.zeros(n)
    unit[nsh] = 1.0
    rf = gauss_filter(spikes, dt, a) / gauss_filter(unit, dt, a).max()
    return rf, np.arange(n) * dt - shift, fit


def rf_quality(rf, t, fit, cfg):
    """Simple acceptance rules for a single receiver function."""
    direct = float(rf[(t > -1.0) & (t < 1.5)].max())
    pre = float(np.max(np.abs(rf[(t > -8.0) & (t < -1.5)])))
    post = float(np.max(np.abs(rf[(t > 1.5) & (t < 30.0)])))
    ok = (
        fit >= cfg['min_fit']
        and cfg['direct_min'] <= direct <= cfg['direct_max']
        and pre <= cfg['pre_ratio'] * direct
        and post <= cfg['post_ratio'] * direct
    )
    return ok, direct, pre, post


def _rf_event(args):
    """All receiver functions of one event folder (runs in a worker)."""
    folder, picks, cfg = args
    from obspy import read

    out = []
    for sta, q in picks.items():
        fz = processed_file(folder, sta, 'Z')
        fr = processed_file(folder, sta, 'R')
        if fz is None or fr is None:
            continue
        try:
            z = read(fz)[0]
            r = read(fr)[0]
            tP = float(z.stats.sac.t0)
            dt = z.stats.delta
            for x in (z, r):
                x.filter(
                    'bandpass', freqmin=cfg['band'][0], freqmax=cfg['band'][1], corners=2, zerophase=True
                )
            i0 = int(round((tP + cfg['win'][0] - z.stats.sac.b) / dt))
            n = int(round((cfg['win'][1] - cfg['win'][0]) / dt))
            if i0 < 0 or i0 + n > min(z.stats.npts, r.stats.npts):
                continue
            Z = np.asarray(z.data[i0 : i0 + n], float)
            R = np.asarray(r.data[i0 : i0 + n], float)
            m = int(2.0 / dt)
            w = np.hanning(2 * m)
            for x in (Z, R):
                x[:m] *= w[:m]
                x[-m:] *= w[m:]
            rf, t, fit = itdecon(R, Z, dt, a=cfg['gauss'], shift=-cfg['win'][0], niter=cfg['niter'])
            ok, direct, pre, post = rf_quality(rf, t, fit, cfg)
            keep = (t >= -5.0) & (t <= cfg['t_end'])
            out.append(
                dict(
                    event=os.path.basename(folder),
                    station=sta,
                    ok=ok,
                    fit=round(fit, 1),
                    direct=round(direct, 3),
                    pre=round(pre, 3),
                    post=round(post, 3),
                    baz=q['baz'],
                    rayp=q['rayp'],
                    dist=q['dist'],
                    quality=q['quality'],
                    rf=rf[keep].astype(np.float32),
                    t0=float(t[keep][0]),
                    dt=dt,
                )
            )
        except Exception as e:
            out.append(dict(event=os.path.basename(folder), station=sta, ok=False, error=str(e)[:80]))
    return out


def hk_stack_each(rfs, t0, dt, rayp_s_km, vp, H, K, weights=(0.6, 0.3, 0.1)):
    """Zhu & Kanamori (2000) H-k stack, one grid (len(H), len(K)) per receiver function, so that the
    stack of any subset (bootstrap) is a mean over these. rfs: (n, m) array, time axis t0 + i dt
    (0 = direct P); rayp in s/km."""
    n, m = rfs.shape
    tt = t0 + np.arange(m) * dt
    Hc = np.asarray(H, float)[:, None]
    b_of_k = np.asarray(K, float)[None, :]
    out = np.zeros((n, len(H), len(K)), dtype=np.float32)
    for i, p in enumerate(rayp_s_km):
        a = math.sqrt(max(1.0 / vp**2 - p**2, 1e-12))
        vs = vp / b_of_k
        b = np.sqrt(np.maximum(1.0 / vs**2 - p**2, 1e-12))
        r = rfs[i]
        g = (
            weights[0] * np.interp((Hc * (b - a)).ravel(), tt, r)
            + weights[1] * np.interp((Hc * (b + a)).ravel(), tt, r)
            - weights[2] * np.interp((2.0 * Hc * b).ravel(), tt, r)
        )
        out[i] = g.reshape(len(H), len(K))
    return out


def hk_stack(rfs, t0, dt, rayp_s_km, vp, H, K, weights=(0.6, 0.3, 0.1)):
    """Mean H-k stack of all functions (grid len(H) x len(K))."""
    return hk_stack_each(rfs, t0, dt, rayp_s_km, vp, H, K, weights).mean(axis=0).astype(float)


def hk_best(S, H, K):
    i, j = np.unravel_index(int(np.argmax(S)), S.shape)
    return float(H[i]), float(K[j]), float(S[i, j])


def cmd_rf(a):
    import pandas as pd
    import pickle
    from concurrent.futures import ProcessPoolExecutor

    out = os.path.join(a.qc, a.out)
    if not os.path.isdir(out):
        os.makedirs(out)
    cfg = dict(
        band=(a.fmin, a.fmax),
        win=(-10.0, 90.0),
        gauss=a.gauss,
        niter=200,
        t_end=60.0,
        min_fit=a.min_fit,
        direct_min=0.08,
        direct_max=1.2,
        pre_ratio=0.5,
        post_ratio=1.5,
    )
    pk = pd.read_csv(os.path.join(a.qc, 'picks', 'picks_absolute.csv'))
    pk = pk[(pk.phase == 'P') & pk.quality.isin(['A', 'B', 'C'][: ['A', 'B', 'C'].index(a.min_quality) + 1])]
    bad = set(a.exclude.split(',')) if a.exclude else set()
    pk = pk[~pk.station.isin(bad)]
    jobs = []
    for ev, g in pk.groupby('event'):
        folder = os.path.join(a.qc, 'processed', ev)
        if not os.path.isdir(folder):
            continue
        picks = dict(
            (
                r.station,
                dict(baz=r.baz_deg, rayp=r.rayp_s_per_deg / 111.195, dist=r.dist_deg, quality=r.quality),
            )
            for r in g.itertuples()
        )
        jobs.append((folder, picks, cfg))
    cache = os.path.join(out, 'rf_all.pkl')
    if os.path.exists(cache) and not a.recompute:
        with open(cache, 'rb') as fh:
            res = pickle.load(fh)
        print('loaded %d receiver functions from %s' % (len(res), cache))
    else:
        print(
            '%d events, %d station-events to deconvolve, %d workers'
            % (len(jobs), sum(len(j[1]) for j in jobs), a.workers)
        )
        res = []
        t0 = time.time()
        if a.workers <= 1:
            for k, j in enumerate(jobs, 1):
                res += _rf_event(j)
        else:
            with ProcessPoolExecutor(max_workers=a.workers) as ex:
                for k, r in enumerate(ex.map(_rf_event, jobs, chunksize=2), 1):
                    res += r
                    if k % 50 == 0:
                        print('  %d / %d events (%.0f s)' % (k, len(jobs), time.time() - t0))
                        sys.stdout.flush()
        with open(cache, 'wb') as fh:
            pickle.dump(res, fh, protocol=4)
        print('deconvolution done: %d receiver functions in %.0f s' % (len(res), time.time() - t0))
    good = [r for r in res if r.get('ok')]
    n_global = len(good)
    print(
        'accepted receiver functions: %d of %d (%.0f%%)'
        % (len(good), len(res), 100.0 * len(good) / max(len(res), 1))
    )
    # ---- per-station coherence screen: the few functions that do not resemble their station stack ----
    # correlation of each function with the leave-one-out stack of its station (-1 to 25 s) and its
    # RMS after the direct pulse (1.5 to 30 s) relative to the station median; a function is removed
    # when the correlation is below --station-cc or the RMS ratio above --station-amp (0 disables)
    for r in res:
        r['ok_global'] = bool(r.get('ok'))
        r['station_screen'] = 'not accepted' if not r.get('ok') else 'kept'
    if a.station_cc > 0 or a.station_amp > 0:
        by0 = {}
        for r in good:
            by0.setdefault(r['station'], []).append(r)
        keep = []
        removed = {}
        for sta in sorted(by0):
            lst = by0[sta]
            if len(lst) < 5:
                keep += lst
                continue
            m = min(len(r['rf']) for r in lst)
            R = np.vstack([r['rf'][:m] for r in lst]).astype(float)
            tt = lst[0]['t0'] + np.arange(m) * lst[0]['dt']
            w = (tt >= -1.0) & (tt <= 25.0)
            post = (tt > 1.5) & (tt < 30.0)
            X = R[:, w]
            Ssum = X.sum(axis=0)
            cc = np.array([np.corrcoef(X[i], (Ssum - X[i]) / (len(X) - 1))[0, 1] for i in range(len(X))])
            rms = np.sqrt((R[:, post] ** 2).mean(axis=1))
            ratio = rms / max(float(np.median(rms)), 1e-9)
            for i, r in enumerate(lst):
                r['station_cc'] = round(float(cc[i]), 3)
                r['amp_ratio'] = round(float(ratio[i]), 2)
                bad = (a.station_cc > 0 and cc[i] < a.station_cc) or (
                    a.station_amp > 0 and ratio[i] > a.station_amp
                )
                if bad:
                    r['ok'] = False
                    r['station_screen'] = 'removed: cc %.2f, amp ratio %.1f' % (cc[i], ratio[i])
                    removed[sta] = removed.get(sta, 0) + 1
                else:
                    keep.append(r)
        good = keep
        print(
            'station screen (cc < %.2f or amp ratio > %.1f): %d of %d accepted functions removed at %d stations'
            % (a.station_cc, a.station_amp, n_global - len(good), n_global, len(removed))
        )
        if removed:
            print('  ' + ', '.join('%s %d' % (k, v) for k, v in sorted(removed.items())))
    summ = pd.DataFrame([dict((k, v) for k, v in r.items() if k not in ('rf',)) for r in res])
    save_csv(summ.round(3), out, 'rf_summary.csv')
    # ---- H-k per station ----
    Hs = np.arange(a.h_min, a.h_max + 1e-9, 0.25)
    Ks = np.arange(a.k_min, a.k_max + 1e-9, 0.005)
    sta_info = pk.groupby('station').agg(
        lat=('sta_lat', 'first'), lon=('sta_lon', 'first'), elev=('sta_elev_m', 'first')
    )
    rows = []
    cand_rows = []
    rng = np.random.RandomState(0)
    vp_table = {}
    if getattr(a, 'vp_table', None):
        vt = pd.read_csv(a.vp_table)
        vp_table = dict(zip(vt.station, vt.vp.astype(float)))
        print(
            'station-specific Vp from %s for %d stations (others use %.2f km/s)'
            % (a.vp_table, len(vp_table), a.vp)
        )
    by = {}
    for r in good:
        by.setdefault(r['station'], []).append(r)
    for sta in sorted(by):
        lst = by[sta]
        if len(lst) < a.min_rf:
            rows.append(dict(station=sta, n_rf=len(lst), note='fewer than %d receiver functions' % a.min_rf))
            continue
        vp = float(vp_table.get(sta, a.vp))  # the Vp used for this station
        m = min(len(r['rf']) for r in lst)
        R = np.vstack([r['rf'][:m] for r in lst]).astype(float)
        t0r, dtr = lst[0]['t0'], lst[0]['dt']
        p = np.array([r['rayp'] for r in lst])
        G = hk_stack_each(R, t0r, dtr, p, vp, Hs, Ks)  # one grid per function
        S = G.mean(axis=0).astype(float)
        H, K, smax = hk_best(S, Hs, Ks)
        B = []
        Gf = G.reshape(len(lst), -1)
        for _ in range(a.bootstrap):
            idx = rng.randint(0, len(lst), len(lst))
            cnt = np.bincount(idx, minlength=len(lst)).astype(np.float32)
            Sb = (cnt @ Gf).reshape(S.shape) / len(lst)  # the same as stacking the resampled set
            B.append(hk_best(Sb, Hs, Ks)[:2])
        B = np.array(B)
        del G, Gf
        # Vp sensitivity: H changes with the assumed Vp
        Hv = [hk_best(hk_stack(R, t0r, dtr, p, v, Hs, Ks), Hs, Ks)[0] for v in (vp - 0.2, vp + 0.2)]
        # how clear is the maximum: fraction of the grid within 90 % of the maximum
        clear = float(np.mean(S >= 0.9 * smax)) * 100.0
        edge = bool(K <= Ks[0] + 0.011 or K >= Ks[-1] - 0.011 or H <= Hs[0] + 0.3 or H >= Hs[-1] - 0.3)
        # a maximum within 0.05 of the Vp/Vs limits is nearly as poorly constrained as one on the edge
        near = bool(edge or K <= Ks[0] + 0.051 or K >= Ks[-1] - 0.051)
        # Ps-only alternative: moveout-corrected Ps time with a fixed Vp/Vs (robust when multiples are weak)
        tps, Hps, Hps_err = ps_depth(R, t0r, dtr, p, vp, a.k_fixed, a.ps_hmin, a.ps_hmax)
        cands = ps_candidates(R, t0r, dtr, p, vp, a.k_fixed)
        for c in cands:
            cand_rows.append(dict(station=sta, depth_km=c[0], t_ps_s=c[1], amplitude=c[2], n_rf=len(lst)))
        row = dict(
            station=sta,
            lat=sta_info.lat.get(sta, np.nan),
            lon=sta_info.lon.get(sta, np.nan),
            elev_m=sta_info.elev.get(sta, np.nan),
            n_rf=len(lst),
            H_km=H,
            H_err_km=round(float(np.std(B[:, 0])), 2),
            kappa=round(K, 3),
            kappa_err=round(float(np.std(B[:, 1])), 3),
            vp_km_s=vp,
            H_vp_minus=Hv[0],
            H_vp_plus=Hv[1],
            stack_max=round(smax, 3),
            grid_pct_within_90=round(clear, 1),
            kappa_at_edge=edge,
            kappa_near_edge=near,
            t_ps_s=round(tps, 2),
            H_ps_km=round(Hps, 1),
            H_ps_err_km=round(Hps_err, 1),
            kappa_fixed=a.k_fixed,
            ps_candidates_km=' '.join('%.1f(%.2f)' % (c[0], c[2]) for c in cands),
            baz_mean_deg=round(
                float(
                    np.degrees(
                        np.arctan2(
                            np.mean(np.sin(np.radians([r['baz'] for r in lst]))),
                            np.mean(np.cos(np.radians([r['baz'] for r in lst]))),
                        )
                    )
                )
                % 360.0,
                1,
            ),
            quality=(
                'C'
                if near
                else (
                    'A'
                    if (np.std(B[:, 0]) <= 1.5 and clear <= 2.0)
                    else 'B' if (np.std(B[:, 0]) <= 3.0 and clear <= 5.0) else 'C'
                )
            ),
        )
        rows.append(row)
        np.savez_compressed(
            os.path.join(out, 'rf_%s.npz' % sta),
            rf=R.astype(np.float32),
            t0=t0r,
            dt=dtr,
            rayp=p,
            baz=np.array([r['baz'] for r in lst]),
            event=np.array([r['event'] for r in lst]),
            S=S.astype(np.float32),
            H=Hs,
            K=Ks,
        )
        _plot_rf_station(
            out,
            sta,
            R,
            t0r,
            dtr,
            p,
            np.array([r['baz'] for r in lst]),
            S,
            Hs,
            Ks,
            H,
            K,
            B,
            vp,
            row,
            gauss=a.gauss,
        )
    hk = pd.DataFrame(rows)
    save_csv(hk.round(3), out, 'hk_results.csv')
    save_csv(pd.DataFrame(cand_rows), out, 'ps_candidates.csv')
    ok = hk[hk.get('H_km', pd.Series(dtype=float)).notna()] if 'H_km' in hk else hk.iloc[0:0]
    print('H-k results for %d stations (%d with >= %d receiver functions)' % (len(hk), len(ok), a.min_rf))
    if len(ok):
        print(
            ok[
                [
                    'station',
                    'n_rf',
                    'H_km',
                    'H_err_km',
                    'kappa',
                    'kappa_err',
                    'kappa_at_edge',
                    'quality',
                    't_ps_s',
                    'H_ps_km',
                    'H_ps_err_km',
                ]
            ].to_string(index=False)
        )
        _plot_hk_map(out, ok)


def ps_depth(R, t0, dt, p, vp, kappa, hmin=24.0, hmax=55.0):
    """Moveout-corrected Ps stack: for each trial depth the Ps time of every receiver function is
    predicted from its ray parameter; the depth with the largest stacked amplitude is returned with
    the error from Vp/Vs +/- 0.05 and the time of the Ps peak for a typical ray parameter.

    The search range (hmin, hmax) should bracket the expected Moho. The defaults suit eastern
    Arabia, where the Moho is 25 to 50 km deep: the search starts at 24 km (Ps later than about
    3 s) so that the direct pulse and shallow conversions from the base of sediments or ophiolite
    are not taken for the Moho, and ends at 55 km because later phases are usually multiples of
    shallow interfaces. Always check the station figure."""
    tt = t0 + np.arange(R.shape[1]) * dt
    vs = vp / kappa
    a_ = np.sqrt(np.maximum(1.0 / vp**2 - p**2, 1e-12))
    b_ = np.sqrt(np.maximum(1.0 / vs**2 - p**2, 1e-12))
    Hgrid = np.arange(hmin, hmax + 1e-9, 0.25)
    amp = np.array(
        [np.mean([np.interp(h * (b_[i] - a_[i]), tt, R[i]) for i in range(len(R))]) for h in Hgrid]
    )
    Hbest = float(Hgrid[int(np.argmax(amp))])
    pm = float(np.median(p))
    tps = Hbest * (math.sqrt(max(1.0 / vs**2 - pm**2, 1e-12)) - math.sqrt(max(1.0 / vp**2 - pm**2, 1e-12)))
    err = []
    for kk in (kappa - 0.05, kappa + 0.05):
        vs2 = vp / kk
        b2 = math.sqrt(max(1.0 / vs2**2 - pm**2, 1e-12))
        err.append(abs(tps / (b2 - math.sqrt(max(1.0 / vp**2 - pm**2, 1e-12))) - Hbest))
    return tps, Hbest, max(err)


def ps_candidates(R, t0, dt, p, vp, kappa, hmin=12.0, hmax=70.0, rel=0.3):
    """All local maxima of the moveout-corrected Ps stack between hmin and hmax (km): candidate
    converters (base of sediments or ophiolite, intracrustal, Moho). Returns a list of
    (depth_km, time_s, amplitude) sorted by depth; amplitudes at least rel times the largest."""
    tt = t0 + np.arange(R.shape[1]) * dt
    vs = vp / kappa
    a_ = np.sqrt(np.maximum(1.0 / vp**2 - p**2, 1e-12))
    b_ = np.sqrt(np.maximum(1.0 / vs**2 - p**2, 1e-12))
    Hgrid = np.arange(hmin, hmax + 1e-9, 0.25)
    amp = np.array(
        [np.mean([np.interp(h * (b_[i] - a_[i]), tt, R[i]) for i in range(len(R))]) for h in Hgrid]
    )
    pm = float(np.median(p))
    ba = math.sqrt(max(1.0 / vs**2 - pm**2, 1e-12)) - math.sqrt(max(1.0 / vp**2 - pm**2, 1e-12))
    out = []
    top = float(amp.max())
    for k in range(1, len(Hgrid) - 1):
        if amp[k] > 0 and amp[k] >= amp[k - 1] and amp[k] > amp[k + 1] and amp[k] >= rel * top:
            out.append((float(Hgrid[k]), round(float(Hgrid[k] * ba), 2), round(float(amp[k]), 3)))
    return out


def _plot_rf_station(out, sta, R, t0, dt, p, baz, S, Hs, Ks, H, K, B, vp, row, gauss=2.0):
    """Station figure: the linear stack of all accepted functions on top, the functions sorted by
    back-azimuth below (with the PdS, PpPs and PpSs+PsPs times predicted by the H-k solution for the
    ray parameter of every trace), and the H-k grid on the right."""
    plt = pyplot()
    fig = plt.figure(figsize=(12, 7.5))
    ax_st = fig.add_axes([0.07, 0.75, 0.50, 0.13])
    ax = fig.add_axes([0.07, 0.08, 0.50, 0.64], sharex=ax_st)
    order = np.argsort(baz)
    tt = t0 + np.arange(R.shape[1]) * dt
    # ---- stack on top ----
    stack = R.mean(axis=0)
    dmax = float(np.max(np.abs(stack[(tt > -1.0) & (tt < 1.5)])) + 1e-9)
    stn = stack / dmax
    ax_st.fill_between(tt, 0, stn, where=stn > 0, color='0.25', lw=0)
    ax_st.fill_between(tt, 0, stn, where=stn < 0, color='0.8', lw=0)
    ax_st.plot(tt, stn, 'k', lw=0.6)
    late = tt >= 2.0
    ax_st.plot(
        tt[late], 4.0 * stn[late], color='0.45', lw=0.7
    )  # the same stack, four times, after the direct pulse
    ax_st.axhline(0, color='0.6', lw=0.4)
    pm = float(np.median(p))
    a_m = np.sqrt(max(1.0 / vp**2 - pm**2, 1e-12))
    b_m = np.sqrt(max((K / vp) ** 2 - pm**2, 1e-12))
    last_x = -99.0
    row_y = 1.03
    for lab, tx, c in (
        ('PdS', H * (b_m - a_m), '0.0'),
        ('PpPs', H * (b_m + a_m), '0.4'),
        ('PpSs+PsPs', 2 * H * b_m, '0.6'),
    ):
        ax_st.axvline(tx, color=c, lw=1.0, ls='--')
        row_y = 1.17 if (tx - last_x) < 2.5 and row_y == 1.03 else 1.03  # stagger labels that would touch
        ax_st.text(
            tx,
            row_y,
            lab,
            fontsize=7,
            ha='center',
            va='bottom',
            color=c,
            transform=ax_st.get_xaxis_transform(),
        )
        last_x = tx
    ax_st.set_ylim(-0.7, 1.1)
    ax_st.set_yticks([0, 1])
    ax_st.tick_params(labelsize=7, labelbottom=False)
    ax_st.set_ylabel('stack', fontsize=8)
    ax_st.text(
        0.5,
        1.33,
        'linear stack of %d functions, normalized to the direct pulse (gray line: x4 after 2 s); dashed: predicted times at the median ray parameter'
        % len(order),
        fontsize=6.5,
        ha='center',
        va='bottom',
        transform=ax_st.transAxes,
        color='0.3',
    )
    ax_st.set_title(
        '%s: %d receiver functions, Gaussian width %.1f' % (sta, len(order), gauss), fontsize=9, pad=34
    )
    # ---- the functions ----
    sc = 1.0 / (np.percentile(np.abs(R), 99) + 1e-9)
    for k, i in enumerate(order):
        y = R[i] * sc * 0.8 + k
        ax.plot(tt, y, 'k', lw=0.3)
        ax.fill_between(tt, k, y, where=y > k, color='0.3', lw=0)
    a_ = np.sqrt(np.maximum(1.0 / vp**2 - p[order] ** 2, 1e-12))
    b_ = np.sqrt(np.maximum((K / vp) ** 2 - p[order] ** 2, 1e-12))
    for lab, tx, c in (
        ('PdS', H * (b_ - a_), '0.0'),
        ('PpPs', H * (b_ + a_), '0.4'),
        ('PpSs+PsPs', 2 * H * b_, '0.6'),
    ):
        ax.plot(tx, np.arange(len(order)), color=c, lw=1.0, ls='--')
    ax.set_xlim(-2, 30)
    ax.set_ylim(-1, len(order) + 1)
    ax.set_xlabel('time after P (s)')
    ax.set_ylabel(
        'receiver functions sorted by back-azimuth (%.0f to %.0f deg)' % (baz[order][0], baz[order][-1])
    )
    ax.tick_params(labelsize=8)
    # ---- H-k grid ----
    ax2 = fig.add_axes([0.65, 0.08, 0.29, 0.76])
    im = ax2.pcolormesh(Ks, Hs, S, cmap='Greys', shading='auto')
    ax2.contour(Ks, Hs, S, levels=[0.9 * S.max()], colors='k', linewidths=0.8)
    ax2.scatter(B[:, 1], B[:, 0], s=4, c='0.5', alpha=0.5)
    ax2.plot(K, H, 'k+', ms=12, mew=1.5)
    ax2.errorbar(
        row['kappa_fixed'],
        row['H_ps_km'],
        yerr=row['H_ps_err_km'],
        fmt='o',
        mfc='white',
        mec='k',
        ecolor='k',
        ms=5,
        lw=0.8,
        capsize=2,
    )
    ax2.set_xlabel('Vp/Vs')
    ax2.set_ylabel('Moho depth H (km), Vp = %.2f km/s' % vp)
    ax2.set_title(
        'H-k: %.1f +/- %.1f km, Vp/Vs %.2f +/- %.2f, quality %s%s\nPdS only (Vp/Vs %.2f, circle): %.1f +/- %.1f km, PdS at %.1f s'
        % (
            H,
            row['H_err_km'],
            K,
            row['kappa_err'],
            row['quality'],
            ' (edge)' if row['kappa_at_edge'] else '',
            row['kappa_fixed'],
            row['H_ps_km'],
            row['H_ps_err_km'],
            row['t_ps_s'],
        ),
        fontsize=8,
    )
    ax2.tick_params(labelsize=8)
    fig.colorbar(im, ax=ax2, shrink=0.5, pad=0.03, label='stack amplitude')
    save_figure(fig, out, 'RF_%s.png' % sta)
    plt.close(fig)


def _plot_hk_map(out, ok):
    plt = pyplot()
    fig, axs = plt.subplots(1, 2, figsize=(12, 5))
    for ax, col, lab, cmap in (
        (axs[0], 'H_km', 'Moho depth (km)', 'viridis_r'),
        (axs[1], 'kappa', 'Vp/Vs', 'viridis'),
    ):
        sc = ax.scatter(ok.lon, ok.lat, c=ok[col], cmap=cmap, s=80, marker='^', edgecolor='k', lw=0.5)
        for r in ok.itertuples():
            ax.annotate(
                '%s\n%s' % (r.station[-2:], ('%.0f' % r.H_km) if col == 'H_km' else ('%.2f' % r.kappa)),
                (r.lon, r.lat),
                textcoords='offset points',
                xytext=(6, 2),
                fontsize=6,
            )
        fig.colorbar(sc, ax=ax, label=lab)
        ax.set_xlabel('longitude')
        ax.set_ylabel('latitude')
        ax.set_aspect(1.0 / math.cos(math.radians(float(ok.lat.mean()))))
    axs[0].set_title('Crustal thickness from H-k stacking', fontsize=9)
    axs[1].set_title('Vp/Vs from H-k stacking', fontsize=9)
    plt.tight_layout()
    save_figure(fig, out, 'hk_map.png')
    plt.close(fig)
