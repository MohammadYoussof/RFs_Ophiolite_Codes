"""Station and event terms of the travel-time residuals.

Model: residual(station, event) = station term + event term, solved by least squares
with the station terms summing to zero. Errors come from a bootstrap over events.
The command also maps the direction dependence of the delays (back-azimuth sectors)
and compares the S terms with the P terms.
"""

import math
import os

import numpy as np

from .common import CLASS_ORDER, output_dir, pyplot, save_figure, save_csv


def _lsq_terms(sta_idx, evt_idx, r, n_sta, n_evt, w_evt=None):
    """r_k = s_i + e_j, sum(s) = 0. Normal equations; w_evt = multiplicity per event (bootstrap)."""
    from scipy import sparse

    n = len(r)
    w = np.ones(n) if w_evt is None else w_evt[evt_idx]
    A = sparse.csr_matrix(
        (np.ones(2 * n), (np.r_[np.arange(n), np.arange(n)], np.r_[sta_idx, n_sta + evt_idx])),
        shape=(n, n_sta + n_evt),
    )
    W = sparse.diags(w)
    N = (A.T @ W @ A).toarray()
    b = A.T @ (w * r)
    c = np.r_[np.ones(n_sta), np.zeros(n_evt)]
    N += 1e3 * np.outer(c, c)  # constraint: station terms sum to zero
    N[np.diag_indices_from(N)] += 1e-9
    x = np.linalg.lstsq(N, b, rcond=None)[0]
    return x[:n_sta], x[n_sta:]


def _terms(df, nboot, seed=0, timing_mad=0.5, outlier_s=0.6, outlier_k=4.0):
    """Station and event terms with two robustness steps:
    1. stations whose residual spread (median absolute deviation) exceeds timing_mad are treated as
       timing suspects (clock faults give offsets that change with time); their picks are excluded
       from the solution and their terms reported separately;
    2. single picks whose residual exceeds max(outlier_s, outlier_k x MAD of all residuals) are dropped
       (cycle skips, mispicks). The solution is repeated after each step."""
    stas = sorted(df.station.unique())
    evts = sorted(df.event.unique())
    si = df.station.map(dict((s, i) for i, s in enumerate(stas))).values
    ei = df.event.map(dict((e, i) for i, e in enumerate(evts))).values
    r = df.res_elev_s.values.astype(float)
    n_sta, n_evt = len(stas), len(evts)
    use = np.ones(len(r), bool)
    timing = np.zeros(n_sta, bool)
    mad_sta = np.full(n_sta, np.nan)

    def solve(mask):
        s_, e_ = _lsq_terms_masked(si, ei, r, n_sta, n_evt, mask)
        return s_, e_, r - s_[si] - e_[ei]

    for it in range(3):
        s, e, resid = solve(use)
        for k in range(n_sta):  # spread over ALL picks of the station (a constant offset does not matter)
            rr = resid[si == k]
            mad_sta[k] = float(np.median(np.abs(rr - np.median(rr)))) if len(rr) >= 5 else np.nan
        new_timing = np.isfinite(mad_sta) & (mad_sta > timing_mad)
        mad_all = float(np.median(np.abs(resid[use] - np.median(resid[use]))))
        thr = max(outlier_s, outlier_k * mad_all)
        new_use = use & ~new_timing[si] & (np.abs(resid) <= thr)
        if np.array_equal(new_use, use) and np.array_equal(new_timing, timing):
            break
        use, timing = new_use, new_timing
    s, e, resid = solve(use)
    s[timing] = np.nan  # a timing suspect has no meaningful term
    rng = np.random.RandomState(seed)
    B = []
    for _ in range(nboot):
        w = np.bincount(rng.randint(0, n_evt, n_evt), minlength=n_evt).astype(float)
        sb, _ = _lsq_terms(si[use], ei[use], r[use], n_sta, n_evt, w)
        has = np.bincount(si[use], weights=w[ei[use]], minlength=n_sta) > 0
        sb[~has] = np.nan
        sb[timing] = np.nan
        B.append(sb - np.nanmean(sb))
    err = np.nanstd(np.array(B), axis=0) if nboot else np.full(n_sta, np.nan)
    err[timing] = np.nan
    return dict(
        stas=stas,
        evts=evts,
        s=s,
        e=e,
        err=err,
        si=si,
        ei=ei,
        resid=resid,
        use=use,
        timing=timing,
        mad_sta=mad_sta,
        n_rejected=int((~use & ~timing[si]).sum()),
    )


def _lsq_terms_masked(sta_idx, evt_idx, r, n_sta, n_evt, mask):
    """Least squares on the picks with mask True; stations/events without picks get 0."""
    s = np.zeros(n_sta)
    e = np.zeros(n_evt)
    sub_s, sub_e = _lsq_terms(sta_idx[mask], evt_idx[mask], r[mask], n_sta, n_evt)
    s[:] = sub_s
    e[:] = sub_e
    return s, e


def _components(si, ei, n_sta, n_evt):
    """Connected groups of stations (linked through shared events). One group = comparable terms."""
    from scipy import sparse
    from scipy.sparse.csgraph import connected_components

    g = sparse.coo_matrix((np.ones(len(si)), (si, n_sta + ei)), shape=(n_sta + n_evt, n_sta + n_evt))
    _, lab = connected_components(g, directed=False)
    return lab[:n_sta]


def _sector_edges(d, spec):
    """Sector edges in degrees from 'start'. spec = 'N' (N sectors with equal numbers of events, cut at the
    largest gap in back-azimuth) or explicit edges '0,60,90,180,360'."""
    if ',' in str(spec):
        e = np.array([float(x) for x in str(spec).split(',')])
        return 0.0, e
    n = max(1, int(spec))
    b = np.sort(d.groupby('event').baz_deg.median().values % 360.0)
    if len(b) < n + 1:
        return 0.0, np.linspace(0, 360, n + 1)
    gaps = np.diff(np.r_[b, b[0] + 360.0])
    k = int(np.argmax(gaps))
    start = (b[k] + gaps[k] / 2.0) % 360.0  # middle of the largest empty gap
    rel = np.sort((b - start) % 360.0)
    inner = [
        0.5 * (rel[int(round(i * len(rel) / n)) - 1] + rel[int(round(i * len(rel) / n))]) for i in range(1, n)
    ]
    return start, np.r_[0.0, inner, 360.0]


def _plot_timing(out, rp, t, phase, timing_mad):
    """Residual against time for the stations with the largest spread: a clock fault shows as offsets
    that change with time, a structural term as a constant."""
    import pandas as pd

    plt = pyplot()
    rp = rp.copy()
    rp['t'] = pd.to_datetime(rp.origin_time.str.replace('Z', ''))
    order = t.sort_values('resid_mad_s', ascending=False)
    show = list(order[order.timing_suspect].station) + list(order[~order.timing_suspect].station)
    show = show[:6]
    fig, axs = plt.subplots(len(show), 1, figsize=(11, 1.9 * len(show) + 0.8), sharex=True)
    axs = np.atleast_1d(axs)
    for ax, sta in zip(axs, show):
        g = rp[rp.station == sta]
        ax.scatter(g.t[g.used], g.resid_s[g.used], s=6, c='k', label='used')
        ax.scatter(g.t[~g.used], g.resid_s[~g.used], s=10, c='none', edgecolor='k', lw=0.6, label='excluded')
        ax.axhline(0, color='0.5', lw=0.6)
        row = t[t.station == sta].iloc[0]
        ax.set_ylabel('%s\nresidual (s)' % sta)
        ax.set_title(
            '%s: residual MAD %.2f s%s'
            % (sta, row.resid_mad_s, '  TIMING SUSPECT' if row.timing_suspect else ''),
            fontsize=8,
            loc='left',
        )
    axs[0].legend(fontsize=7, frameon=False, loc='upper right')
    axs[-1].set_xlabel('event date')
    plt.tight_layout()
    save_figure(fig, out, 'timing_check_%s.png' % phase)
    plt.close(fig)


def cmd_station_terms(a):
    import pandas as pd

    plt = pyplot()
    out = output_dir(a)
    pk = pd.read_csv(os.path.join(a.qc, 'picks', 'picks_absolute.csv'))
    keep = CLASS_ORDER[: CLASS_ORDER.index(a.min_quality) + 1]
    pk = pk[pk.quality.isin(keep)]
    if a.exclude:
        bad = [x.strip() for x in a.exclude.split(',') if x.strip()]
        pk = pk[~pk.station.isin(bad)]
        print('excluded by request: %s' % ', '.join(bad))
    terms = {}
    for phase in ('P', 'S'):
        d = pk[pk.phase == phase].copy()
        if d.station.nunique() < 3 or d.event.nunique() < 2:
            print('%s: not enough picks for station terms (%d picks)' % (phase, len(d)))
            continue
        T = _terms(d, a.bootstrap, timing_mad=a.timing_mad)
        stas, evts, s, e, err, si, ei, resid, use = (
            T[k] for k in ('stas', 'evts', 's', 'e', 'err', 'si', 'ei', 'resid', 'use')
        )
        comp = _components(si[use], ei[use], len(stas), len(evts))
        info = d.groupby('station').agg(
            lat=('sta_lat', 'first'),
            lon=('sta_lon', 'first'),
            elev_m=('sta_elev_m', 'first'),
            n_picks=('event', 'size'),
            rayp_s_per_deg=('rayp_s_per_deg', 'mean'),
        )
        t = pd.DataFrame(dict(station=stas, term_s=s, err_s=err)).merge(info.reset_index(), on='station')
        # simple alternative for comparison: mean of the event-demeaned residuals
        rel = d.res_elev_s - d.groupby('event').res_elev_s.transform('mean')
        t['mean_dt_rel_s'] = t.station.map(rel.groupby(d.station).mean())
        main = np.bincount(comp).argmax()
        t['connected'] = comp == main
        t['timing_suspect'] = T['timing']
        t['resid_mad_s'] = np.round(T['mad_sta'], 3)
        t['n_used'] = np.bincount(si[use], minlength=len(stas))
        t['reliable'] = (t.n_used >= a.min_picks) & t.connected & ~t.timing_suspect
        if T['timing'].any():
            print(
                'WARNING %s: timing suspects (residual MAD > %.2f s), excluded from the solution: %s'
                % (phase, a.timing_mad, ', '.join(t.station[t.timing_suspect]))
            )
        print('%s: %d single picks rejected as outliers' % (phase, T['n_rejected']))
        # per-pick residuals after the solution, for inspection and for the timing figure
        rp = d[['event', 'origin_time', 'station', 'quality', 'baz_deg', 'dist_deg']].copy()
        rp['resid_s'] = np.round(resid, 4)
        rp['used'] = use
        save_csv(rp, out, 'pick_residuals_%s.csv' % phase)
        _plot_timing(out, rp, t, phase, a.timing_mad)
        iso = t[~t.connected & ~t.timing_suspect & (t.n_used > 0)]
        if len(iso):
            print(
                'WARNING %s: stations %s share no events with the rest of the network; their terms are not '
                'comparable and are marked unreliable' % (phase, ', '.join(iso.station))
            )
        # back-azimuth coverage per station: a term averages over the directions that station saw
        rad = np.radians(d.baz_deg.values)
        cs = (
            pd.DataFrame(dict(station=d.station.values, c=np.cos(rad), s=np.sin(rad)))
            .groupby('station')
            .mean()
        )
        t['baz_mean_deg'] = t.station.map((np.degrees(np.arctan2(cs.s, cs.c)) % 360.0).round(1))
        t['baz_concentration'] = t.station.map(np.hypot(cs.s, cs.c).round(2))  # 1 = all from one direction
        t = t.round(4)
        save_csv(t, out, 'station_terms_%s.csv' % phase)
        terms[phase] = t
        ev_t = pd.DataFrame(dict(event=evts, event_term_s=e)).round(4)
        save_csv(ev_t, out, 'event_terms_%s.csv' % phase)
        if len(evts) < 10:
            print('%s: only %d events - bootstrap errors are not meaningful yet' % (phase, len(evts)))
        ok = t.reliable.values
        print(
            '%s: %d picks (%d used), %d stations (%d reliable), %d events; term range %.2f to %.2f s, median bootstrap '
            'error %.3f s, rms misfit of used picks %.3f s'
            % (
                phase,
                len(d),
                int(use.sum()),
                len(stas),
                int(ok.sum()),
                len(evts),
                np.nanmin(s[ok]) if ok.any() else np.nan,
                np.nanmax(s[ok]) if ok.any() else np.nan,
                np.nanmedian(err),
                float(np.sqrt(np.mean(resid[use] ** 2))),
            )
        )
        # --- map of the terms ---
        tt = t[t.reliable]
        if len(tt):
            fig, ax = plt.subplots(figsize=(7.5, 6))
            v = max(float(np.nanmax(np.abs(tt.term_s))), 0.05)
            sc = ax.scatter(
                tt.lon,
                tt.lat,
                c=tt.term_s,
                cmap='RdBu_r',
                vmin=-v,
                vmax=v,
                s=90,
                marker='^',
                edgecolor='k',
                lw=0.5,
            )
            for r in tt.itertuples():
                ax.annotate(
                    '%s\n%+.2f' % (r.station[-2:], r.term_s),
                    (r.lon, r.lat),
                    textcoords='offset points',
                    xytext=(6, 2),
                    fontsize=6,
                )
            fig.colorbar(sc, ax=ax, label='%s station term (s), positive = late (slow)' % phase)
            ax.set_xlabel('longitude')
            ax.set_ylabel('latitude')
            ax.set_aspect(1.0 / math.cos(math.radians(float(tt.lat.mean()))))
            ax.set_title(
                '%s station terms (elevation-corrected), stations with >= %d picks' % (phase, a.min_picks)
            )
            save_figure(fig, out, 'station_terms_%s.png' % phase)
            plt.close(fig)
        # --- direction dependence ---
        start, edges = _sector_edges(d, a.sectors)
        nsec = len(edges) - 1
        lab = ['%d-%d' % ((start + edges[k]) % 360, (start + edges[k + 1]) % 360 or 360) for k in range(nsec)]
        d['dev'] = resid
        d['sector'] = np.digitize((d.baz_deg.values - start) % 360.0, edges[1:-1])
        print(
            '%s back-azimuth sectors (deg): %s; events per sector: %s'
            % (
                phase,
                ', '.join(lab),
                ', '.join(str(int(d[d.sector == k].event.nunique())) for k in range(nsec)),
            )
        )
        dd = d[use & d.station.isin(t.station[t.reliable]).values]
        sec = dd.groupby(['station', 'sector']).agg(n=('dev', 'size'), dev_s=('dev', 'mean')).reset_index()
        sec['sector_deg'] = [lab[k] for k in sec.sector]
        sec = sec.merge(t[['station', 'lat', 'lon', 'term_s']], on='station')
        sec['term_plus_dev_s'] = sec.term_s + sec.dev_s
        save_csv(sec.round(4), out, 'station_terms_by_sector_%s.csv' % phase)
        fig, axs = plt.subplots(1, nsec, figsize=(3.6 * nsec, 3.8), sharex=True, sharey=True)
        axs = np.atleast_1d(axs)
        v = max(
            float(np.nanmax(np.abs(sec[sec.n >= 3].term_plus_dev_s))) if (sec.n >= 3).any() else 0.1, 0.05
        )
        sc = None
        for k, ax in enumerate(axs):
            g = sec[(sec.sector == k) & (sec.n >= 3)]
            ne = int(d[d.sector == k].event.nunique())
            if len(g):
                sc = ax.scatter(
                    g.lon,
                    g.lat,
                    c=g.term_plus_dev_s,
                    cmap='RdBu_r',
                    vmin=-v,
                    vmax=v,
                    s=50,
                    marker='^',
                    edgecolor='k',
                    lw=0.4,
                )
            ax.set_title('back-azimuth %s deg (%d events)' % (lab[k], ne), fontsize=8)
            ax.set_xlabel('longitude')
        axs[0].set_ylabel('latitude')
        if sc is not None:
            fig.colorbar(sc, ax=list(axs), label='%s delay for this direction (s)' % phase, shrink=0.8)
        save_figure(fig, out, 'station_terms_by_sector_%s.png' % phase)
        plt.close(fig)
    if 'P' in terms and 'S' in terms:
        m = terms['P'].merge(terms['S'], on='station', suffixes=('_P', '_S'))
        m = m[m.reliable_P & m.reliable_S]
        if len(m) >= 5:
            slope = float(np.polyfit(m.term_s_P, m.term_s_S, 1)[0])
            fig, ax = plt.subplots(figsize=(4.8, 4.5))
            ax.errorbar(
                m.term_s_P, m.term_s_S, xerr=m.err_s_P, yerr=m.err_s_S, fmt='o', color='k', ms=3, lw=0.6
            )
            xx = np.linspace(m.term_s_P.min(), m.term_s_P.max(), 10)
            ax.plot(xx, slope * xx, 'k--', lw=0.8, label='slope %.1f' % slope)
            ax.set_xlabel('P station term (s)')
            ax.set_ylabel('S station term (s)')
            ax.legend(frameon=False)
            ax.set_title('S against P station terms (%d stations)' % len(m))
            save_figure(fig, out, 's_vs_p_terms.png')
            plt.close(fig)
            print('S/P station-term slope: %.2f (%d stations)' % (slope, len(m)))
        else:
            print('S/P comparison skipped: fewer than 5 stations with reliable P and S terms')
