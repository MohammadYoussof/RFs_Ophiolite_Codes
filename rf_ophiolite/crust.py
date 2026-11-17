"""Crustal correction of the station terms.

Commands:
    crust-template  empty crust table to fill in by hand
    crust-from-rf   crust table from the receiver-function results, published Moho
                    depths and, optionally, a 3-D shear-velocity model
    crust           predicted crustal delay of every station and the residual
                    (mantle) delay
"""

import math
import os
import sys

import numpy as np

from .common import ANALYSIS_DIR, output_dir, pyplot, save_figure, save_csv


def cmd_crust_template(a):
    import pandas as pd

    out = output_dir(a)
    pk = pd.read_csv(os.path.join(a.qc, 'picks', 'picks_absolute.csv'))
    t = (
        pk.groupby('station')
        .agg(lat=('sta_lat', 'first'), lon=('sta_lon', 'first'), elev_m=('sta_elev_m', 'first'))
        .reset_index()
    )
    t['moho_km'] = ''  # Moho depth below the station (km), from receiver functions or a published map
    t['moho_err_km'] = ''  # its uncertainty (km); empty = --moho-err (default 3 km)
    t['ophiolite_km'] = ''  # ophiolite thickness below the station (km); 0 where there is none
    t['sediment_km'] = ''  # sediment thickness below the station (km); 0 or empty where there is none
    t['source'] = ''  # reference for each value
    p = os.path.join(out, 'crust_model_template.csv')
    if os.path.exists(p) and not a.force:
        sys.exit('%s exists (use --force to overwrite it)' % p)
    t.to_csv(p, index=False)
    print(
        '  wrote',
        p,
        '- fill moho_km, moho_err_km, ophiolite_km, sediment_km and source, then save it as crust_model.csv',
    )


def _col_time(layers, p):
    """Vertical travel time through layers [(thickness_km, vp_km_s)], ray parameter p (s/km)."""
    return sum(h * math.sqrt(max(1.0 / v**2 - p**2, 0.0)) for h, v in layers if h > 0)


def _orla_profiles(path, lat, lon, zmax):
    """Vs(z) at each station from a 3-D shear-velocity model in netCDF (variables depth, latitude,
    longitude, vs), in 0.5 km steps down to zmax. Written for ORLA2022 (Weidle et al. 2022)."""
    import netCDF4 as nc
    from scipy.interpolate import RegularGridInterpolator

    d = nc.Dataset(path)
    z, la, lo, V = d['depth'][:], d['latitude'][:], d['longitude'][:], np.array(d['vs'][:])
    f = RegularGridInterpolator((z, la, lo), V, bounds_error=False, fill_value=np.nan)
    zz = np.arange(0.0, zmax + 1e-9, 0.5)
    return zz, [f(np.c_[zz, np.full_like(zz, y), np.full_like(zz, x)]) for y, x in zip(lat, lon)]


def cmd_crust_from_rf(a):
    """Build crust_model.csv from the receiver-function results, with published Moho depths and a
    3-D shear-velocity model for comparison where they are given."""
    import pandas as pd

    plt = pyplot()
    out = output_dir(a)
    hk_all = pd.read_csv(os.path.join(a.qc, a.rf, 'hk_results.csv'))
    hk = (
        hk_all[hk_all.get('H_km', pd.Series(dtype=float)).notna()].copy()
        if 'H_km' in hk_all
        else hk_all.iloc[0:0]
    )
    rows = []
    for r in hk.itertuples():
        near = bool(getattr(r, 'kappa_near_edge', r.kappa_at_edge))
        use_hk = (r.quality in ('A', 'B')) and not near
        moho = r.H_km if use_hk else r.H_ps_km
        err = r.H_err_km if use_hk else max(r.H_ps_err_km, a.ps_err)
        rows.append(
            dict(
                station=r.station,
                lat=r.lat,
                lon=r.lon,
                elev_m=r.elev_m,
                n_rf=r.n_rf,
                moho_km=round(float(moho), 1),
                moho_err_km=round(max(float(err), a.min_err), 1),
                method=('H-k' if use_hk else 'Ps only, Vp/Vs %.2f (uncertain)' % r.kappa_fixed),
                H_hk_km=r.H_km,
                H_hk_err_km=r.H_err_km,
                kappa=r.kappa,
                kappa_at_edge=r.kappa_at_edge,
                kappa_near_edge=near,
                H_ps_km=r.H_ps_km,
                H_ps_err_km=r.H_ps_err_km,
                rf_quality=r.quality,
                ps_candidates_km=getattr(r, 'ps_candidates_km', ''),
                ophiolite_km=0.0,
                sediment_km=0.0,
                source='this study, receiver functions',
            )
        )
    m = pd.DataFrame(rows)
    if a.published:
        pub = pd.read_csv(a.published)
        pub = pub.rename(
            columns=dict(moho_km='moho_pub_km', moho_err_km='moho_pub_err_km', source='source_pub')
        )
        pub = pub[[c for c in ('station', 'moho_pub_km', 'moho_pub_err_km', 'source_pub') if c in pub]]
        # stations with a published value but no receiver-function result get a line too (prefer auto or published)
        if a.prefer != 'rf':
            pk = pd.read_csv(os.path.join(a.qc, 'picks', 'picks_absolute.csv'))
            info = pk.groupby('station').agg(
                lat=('sta_lat', 'first'), lon=('sta_lon', 'first'), elev_m=('sta_elev_m', 'first')
            )
            extra = [
                dict(
                    station=st,
                    lat=info.lat.get(st, np.nan),
                    lon=info.lon.get(st, np.nan),
                    elev_m=info.elev_m.get(st, np.nan),
                    n_rf=int(hk_all.set_index('station').n_rf.get(st, 0)) if len(hk_all) else 0,
                    moho_km=np.nan,
                    moho_err_km=np.nan,
                    method='none',
                    rf_quality='',
                    ophiolite_km=0.0,
                    sediment_km=0.0,
                    source='',
                )
                for st in pub.station
                if st not in set(m.station)
            ]
            if extra:
                m = (
                    pd.concat([m, pd.DataFrame(extra)], ignore_index=True)
                    .sort_values('station')
                    .reset_index(drop=True)
                )
        m = m.merge(pub, on='station', how='left')
        if a.prefer == 'auto':
            # the published value replaces a Ps-only (uncertain) depth, never an A or B H-k result
            k = m.moho_pub_km.notna() & (m.method != 'H-k')
            m.loc[k, 'moho_km'] = m.loc[k, 'moho_pub_km']
            m.loc[k, 'moho_err_km'] = m.loc[k, 'moho_pub_err_km'].fillna(a.min_err)
            m.loc[k, 'method'] = 'published'
            m.loc[k, 'source'] = m.loc[k, 'source_pub'].fillna('published')
        both = m[m.moho_pub_km.notna() & m.get('H_ps_km', pd.Series(np.nan, index=m.index)).notna()].copy()
        both['moho_rf_km'] = np.where(both.rf_quality.isin(['A', 'B']), both.H_hk_km, both.H_ps_km)
        both['moho_rf_err_km'] = np.where(
            both.rf_quality.isin(['A', 'B']), both.H_hk_err_km, both.H_ps_err_km
        )
        if len(both):
            dif = both.moho_rf_km - both.moho_pub_km
            print(
                '%d stations with a published Moho: this study minus published, mean %.1f km, MAD %.1f km'
                % (len(both), dif.mean(), (dif - dif.median()).abs().median())
            )
            fig, ax = plt.subplots(figsize=(5.2, 4.8))
            hk_ = both[both.rf_quality.isin(['A', 'B'])]
            hc_ = both[~both.rf_quality.isin(['A', 'B'])]
            ax.errorbar(
                hk_.moho_pub_km,
                hk_.H_hk_km,
                xerr=hk_.moho_pub_err_km.fillna(0.0),
                yerr=hk_.H_hk_err_km,
                fmt='s',
                color='k',
                ms=5,
                lw=0.6,
                label='H-k, quality A or B',
            )
            ax.errorbar(
                hc_.moho_pub_km,
                hc_.H_hk_km,
                xerr=hc_.moho_pub_err_km.fillna(0.0),
                yerr=hc_.H_hk_err_km,
                fmt='s',
                mfc='0.6',
                mec='0.6',
                color='0.6',
                ms=5,
                lw=0.6,
                label='H-k, quality C (Vp/Vs poorly constrained)',
            )
            ax.errorbar(
                both.moho_pub_km,
                both.H_ps_km,
                xerr=0,
                yerr=both.H_ps_err_km,
                fmt='o',
                mfc='white',
                mec='k',
                color='k',
                ms=4,
                lw=0.6,
                label='Ps only, Vp/Vs 1.75',
            )
            for r in both.itertuples():
                ax.annotate(
                    r.station[-2:],
                    (r.moho_pub_km, r.H_ps_km),
                    textcoords='offset points',
                    xytext=(4, -8),
                    fontsize=7,
                )
            ax.legend(frameon=False, fontsize=7, loc='upper left')
            lim = [
                min(both.moho_pub_km.min(), both.H_ps_km.min(), both.H_hk_km.min()) - 3,
                max(both.moho_pub_km.max(), both.H_ps_km.max(), both.H_hk_km.max()) + 3,
            ]
            ax.plot(lim, lim, 'k--', lw=0.8)
            ax.set_xlim(lim)
            ax.set_ylim(lim)
            ax.set_xlabel('published Moho depth (km)')
            ax.set_ylabel('this study (km)')
            ax.set_title('Moho depth: this study against published values', fontsize=9)
            save_figure(fig, out, 'moho_vs_published.png')
            plt.close(fig)
        if a.prefer == 'published':
            k = m.moho_pub_km.notna()
            m.loc[k, 'moho_km'] = m.loc[k, 'moho_pub_km']
            m.loc[k, 'moho_err_km'] = m.loc[k, 'moho_pub_err_km'].fillna(a.min_err)
            m.loc[k, 'source'] = m.loc[k, 'source_pub'].fillna('published')
            m.loc[k, 'method'] = 'published'
    if a.orla:
        zz, profs = _orla_profiles(a.orla, m.lat.values, m.lon.values, a.orla_zmax)
        vs05, vsmean, tcol = [], [], []
        for prof in profs:
            ok = ~np.isnan(prof)
            vs05.append(float(np.nanmean(prof[zz <= 5.0])))
            vsmean.append(float(np.nanmean(prof[ok])))
            vp = a.orla_kappa * prof
            # vertical P time through the ORLA column (NaN below the model bottom: filled with the last value)
            vp = pd.Series(vp).ffill().bfill().values
            tcol.append(float(np.sum(0.5 / vp[1:])))
        m['orla_vs_0_5_km_s'] = np.round(vs05, 2)
        m['orla_vs_mean_km_s'] = np.round(vsmean, 2)
        m['orla_tP_0_%g_s' % a.orla_zmax] = np.round(tcol, 3)
        m['orla_kappa'] = a.orla_kappa
        m['orla_zmax_km'] = a.orla_zmax
        print(
            'ORLA2022 column 0 to %g km: vertical P time %.2f to %.2f s (Vp = %.2f Vs), mean Vs %.2f to %.2f km/s'
            % (a.orla_zmax, min(tcol), max(tcol), a.orla_kappa, min(vsmean), max(vsmean))
        )
    save_csv(m, out, 'crust_model.csv')
    print(
        '%d stations in crust_model.csv: %d from H-k (quality A or B), %d published, %d from the Ps time only (uncertain)'
        % (
            len(m),
            int((m.method == 'H-k').sum()),
            int((m.method == 'published').sum()),
            int(m.method.str.startswith('Ps').sum()),
        )
    )


def cmd_crust(a):
    import pandas as pd

    plt = pyplot()
    out = output_dir(a)
    mod = pd.read_csv(a.model)
    for c in ('station', 'moho_km', 'ophiolite_km'):
        if c not in mod:
            sys.exit('%s needs a column %s' % (a.model, c))
    t = pd.read_csv(a.terms or os.path.join(a.qc, ANALYSIS_DIR, 'station_terms_P.csv'))
    cols = [c for c in ('station', 'moho_km', 'moho_err_km', 'ophiolite_km', 'sediment_km') if c in mod]
    d = t.merge(mod[cols], on='station', how='inner')
    d = d[d.reliable & pd.to_numeric(d.moho_km, errors='coerce').notna()].copy()
    d['moho_km'] = d.moho_km.astype(float)
    for c, default in (('ophiolite_km', 0.0), ('sediment_km', 0.0), ('moho_err_km', a.moho_err)):
        d[c] = pd.to_numeric(d[c], errors='coerce').fillna(default) if c in d else default
    if len(d) < 3:
        sys.exit('fewer than 3 stations with a Moho depth and a reliable station term')
    elev = (d.elev_m.fillna(0.0).clip(lower=0) / 1000.0).values
    z_ref = float((d.moho_km.values + d.moho_err_km.values - elev).max()) + 5.0

    def delay(r, el, moho):
        """Crustal column (from the station down, converted to below sea level, because the residuals are
        already elevation-corrected to sea level) minus the ak135 column, both down to z_ref."""
        p = r.rayp_s_per_deg / 111.195
        tops = np.array([0.0, r.sediment_km, r.sediment_km + r.ophiolite_km, moho])  # below the station
        z = np.clip(np.minimum(tops, moho) - el, 0.0, None)  # below sea level
        model = [
            (z[1] - z[0], a.vp_sediment),
            (z[2] - z[1], a.vp_ophiolite),
            (z[3] - z[2], a.vp_crust),
            (z_ref - z[3], a.vp_mantle),
        ]
        ref = [(min(20.0, z_ref), 5.8), (min(35.0, z_ref) - 20.0, 6.5), (z_ref - 35.0, 8.04)]  # ak135
        return _col_time(model, p) - _col_time(ref, p)

    pred, perr = [], []
    for r, el in zip(d.itertuples(), elev):
        pred.append(delay(r, el, r.moho_km))
        perr.append(
            abs(delay(r, el, r.moho_km + r.moho_err_km) - delay(r, el, r.moho_km - r.moho_err_km)) / 2.0
        )
    d['crust_delay_s'] = pred
    d['crust_err_s'] = perr  # from the Moho uncertainty
    d['crust_term_s'] = d.crust_delay_s - d.crust_delay_s.mean()  # relative, like the station terms
    d['obs_term_s'] = d.term_s - d.term_s.mean()
    d['residual_s'] = d.obs_term_s - d.crust_term_s
    d['residual_err_s'] = np.sqrt(d.err_s.fillna(0.0) ** 2 + d.crust_err_s**2)
    # travel time along the ray through 'mantle_km' of mantle: h / (v cos i)
    pv = d.rayp_s_per_deg.values / 111.195
    tpath = a.mantle_km / (a.vp_mantle * np.sqrt(np.clip(1.0 - (pv * a.vp_mantle) ** 2, 1e-6, None)))
    d['mantle_dVp_pct'] = -100.0 * d.residual_s / tpath
    if a.orla:
        # alternative crust: the ORLA2022 Vs model (Vp = kappa * Vs) from sea level down to orla_zmax, the uniform
        # crust below that to the Moho, then mantle; the same ak135 reference
        zz, profs = _orla_profiles(a.orla, d.lat.values, d.lon.values, a.orla_zmax)

        def delay_orla(r, el, moho, prof):
            p = r.rayp_s_per_deg / 111.195
            vp = pd.Series(a.orla_kappa * prof).ffill().bfill().values
            zb = max(moho - el, 0.0)  # Moho below sea level
            layers = [(0.5, vp[k]) for k in range(1, len(zz)) if zz[k] <= zb]
            top = min(zb, float(zz[-1]))
            layers += [(zb - top, a.vp_crust), (z_ref - zb, a.vp_mantle)]
            ref = [(min(20.0, z_ref), 5.8), (min(35.0, z_ref) - 20.0, 6.5), (z_ref - 35.0, 8.04)]
            return _col_time(layers, p) - _col_time(ref, p)

        d['crust_delay_orla_s'] = [
            delay_orla(r, el, r.moho_km, prof) for r, el, prof in zip(d.itertuples(), elev, profs)
        ]
        d['crust_term_orla_s'] = d.crust_delay_orla_s - d.crust_delay_orla_s.mean()
        d['residual_orla_s'] = d.obs_term_s - d.crust_term_orla_s
        d['mantle_dVp_orla_pct'] = -100.0 * d.residual_orla_s / tpath
    save_csv(d.round(4), out, 'crust_correction.csv')
    cc = float(np.corrcoef(d.obs_term_s, d.crust_term_s)[0, 1])
    print(
        '%d stations; observed term spread %.2f s, predicted crustal spread %.2f s, residual spread %.2f s; '
        'correlation observed-crust %.2f'
        % (len(d), d.obs_term_s.std(), d.crust_term_s.std(), d.residual_s.std(), cc)
    )
    print(
        'median uncertainty of the crustal correction %.2f s; of the residual %.2f s'
        % (d.crust_err_s.median(), d.residual_err_s.median())
    )
    print(
        'residual as a mantle anomaly over %g km: %.1f to %.1f %% Vp'
        % (a.mantle_km, d.mantle_dVp_pct.min(), d.mantle_dVp_pct.max())
    )
    if a.orla:
        cco = float(np.corrcoef(d.obs_term_s, d.crust_term_orla_s)[0, 1])
        ccr = float(np.corrcoef(d.residual_s, d.residual_orla_s)[0, 1])
        print(
            'with the ORLA2022 crust (0 to %g km, Vp = %.2f Vs): crustal spread %.2f s, residual spread %.2f s, '
            'correlation observed-crust %.2f; the two residuals correlate at %.2f, their difference has a spread of %.2f s'
            % (
                a.orla_zmax,
                a.orla_kappa,
                d.crust_term_orla_s.std(),
                d.residual_orla_s.std(),
                cco,
                ccr,
                (d.residual_s - d.residual_orla_s).std(),
            )
        )
        fig, axs = plt.subplots(1, 3, figsize=(13, 4.2))
        v = max(
            float(np.nanmax(np.abs(d[['obs_term_s', 'crust_term_orla_s', 'residual_orla_s']].values))), 0.05
        )
        for ax, col, lab in (
            (axs[0], 'crust_term_orla_s', 'predicted, ORLA2022 crust'),
            (axs[1], 'residual_orla_s', 'residual, ORLA2022 crust'),
        ):
            sc = ax.scatter(
                d.lon,
                d.lat,
                c=d[col],
                cmap='RdBu_r',
                vmin=-v,
                vmax=v,
                s=70,
                marker='^',
                edgecolor='k',
                lw=0.4,
            )
            ax.set_title('P station term: %s' % lab, fontsize=9)
            ax.set_xlabel('longitude')
            ax.set_aspect(1.0 / math.cos(math.radians(float(d.lat.mean()))))
        axs[0].set_ylabel('latitude')
        fig.colorbar(sc, ax=list(axs[:2]), label='s (positive = late)', shrink=0.8)
        axs[2].plot([-v, v], [-v, v], 'k--', lw=0.8)
        axs[2].errorbar(
            d.residual_s, d.residual_orla_s, xerr=d.residual_err_s, fmt='o', color='k', ms=3, lw=0.6
        )
        axs[2].set_xlabel('residual, uniform crust + Moho (s)')
        axs[2].set_ylabel('residual, ORLA2022 crust + Moho (s)')
        axs[2].set_title('Two crustal models (r = %.2f)' % ccr, fontsize=9)
        save_figure(fig, out, 'crust_correction_orla.png')
        plt.close(fig)
    fig, axs = plt.subplots(1, 4, figsize=(17, 4.2))
    v = max(float(np.nanmax(np.abs(d[['obs_term_s', 'crust_term_s', 'residual_s']].values))), 0.05)
    for ax, col, lab in (
        (axs[0], 'obs_term_s', 'observed'),
        (axs[1], 'crust_term_s', 'predicted from the crust'),
        (axs[2], 'residual_s', 'residual = observed - crust'),
    ):
        sc = ax.scatter(
            d.lon, d.lat, c=d[col], cmap='RdBu_r', vmin=-v, vmax=v, s=70, marker='^', edgecolor='k', lw=0.4
        )
        ax.set_title('P station term: %s' % lab, fontsize=9)
        ax.set_xlabel('longitude')
        ax.set_aspect(1.0 / math.cos(math.radians(float(d.lat.mean()))))
    axs[0].set_ylabel('latitude')
    fig.colorbar(sc, ax=list(axs[:3]), label='s (positive = late)', shrink=0.8)
    axs[3].errorbar(
        d.crust_term_s, d.obs_term_s, xerr=d.crust_err_s, yerr=d.err_s, fmt='o', color='k', ms=3, lw=0.6
    )
    lim = [-v, v]
    axs[3].plot(lim, lim, 'k--', lw=0.8)
    axs[3].set_xlabel('predicted crustal term (s)')
    axs[3].set_ylabel('observed station term (s)')
    axs[3].set_title('Observed against crust (r = %.2f)' % cc, fontsize=9)
    save_figure(fig, out, 'crust_correction.png')
    plt.close(fig)
    # east-west and north-south profiles of the residual
    fig, axs = plt.subplots(1, 2, figsize=(10, 3.8))
    for ax, x, lab in ((axs[0], d.lon, 'longitude'), (axs[1], d.lat, 'latitude')):
        ax.errorbar(x, d.residual_s, yerr=d.residual_err_s, fmt='o', color='k', ms=3, lw=0.6)
        ax.errorbar(x, d.obs_term_s, fmt='o', mfc='none', mec='0.6', ms=3, lw=0)
        for r, xx in zip(d.itertuples(), x):
            ax.annotate(
                r.station[-2:],
                (xx, r.residual_s),
                textcoords='offset points',
                xytext=(3, 3),
                fontsize=6,
                color='0.3',
            )
        ax.axhline(0, color='0.7', lw=0.6)
        ax.set_xlabel(lab)
        ax.set_ylabel('s (positive = late)')
    axs[0].set_title('Residual (black) and observed term (gray) against longitude', fontsize=9)
    axs[1].set_title('... against latitude', fontsize=9)
    save_figure(fig, out, 'residual_profiles.png')
    plt.close(fig)
    med = float(d.lon.median())
    w, e = d[d.lon < med], d[d.lon >= med]
    print(
        'residual west of %.2f E: mean %+.2f s (%d stations); east: %+.2f s (%d); correlation with longitude %.2f, with latitude %.2f'
        % (
            med,
            w.residual_s.mean(),
            len(w),
            e.residual_s.mean(),
            len(e),
            float(np.corrcoef(d.lon, d.residual_s)[0, 1]),
            float(np.corrcoef(d.lat, d.residual_s)[0, 1]),
        )
    )
