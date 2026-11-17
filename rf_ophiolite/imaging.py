"""Images from the accepted receiver functions of all stations.

Commands:
    section  one moveout-corrected stack per station, sorted along the network, with the
             Ps time of the Moho marked
    ccp      common-conversion-point stack: every function is mapped from time to depth
             along its own ray in a 1-D model, placed at the conversion points and
             binned on a longitude, latitude, depth grid; profiles and depth slices

Both read the station files rf_<station>.npz written by the rf command.
"""

import glob
import math
import os
import sys

import numpy as np

from .common import ANALYSIS_DIR, output_dir, pyplot, save_csv, save_figure

KM_PER_DEG = 111.195


# -----------------------------------------------------------------------------
# shared helpers
# -----------------------------------------------------------------------------
def load_station_rfs(rf_dir, min_rf=20):
    """Accepted receiver functions per station: {station: dict(rf, t0, dt, p, baz)}."""
    out = {}
    for f in sorted(glob.glob(os.path.join(rf_dir, 'rf_*.npz'))):
        sta = os.path.basename(f)[3:-4]
        d = np.load(f, allow_pickle=True)
        if d['rf'].shape[0] >= min_rf:
            out[sta] = dict(
                rf=d['rf'].astype(float),
                t0=float(d['t0']),
                dt=float(d['dt']),
                p=d['rayp'].astype(float),
                baz=d['baz'].astype(float),
            )
    return out


def vertical_slowness(p, vp, kappa):
    """Vertical slowness of P and S (s/km) for a ray parameter p (s/km)."""
    vs = vp / kappa
    return math.sqrt(max(1.0 / vp**2 - p**2, 1e-12)), math.sqrt(max(1.0 / vs**2 - p**2, 1e-12))


def moveout_stack(R, t0, dt, p, vp=6.3, kappa=1.75, pref=0.06, tmin=-2.0, tmax=30.0):
    """Stack of the functions R after stretching the time axis of each so that the Ps phase of
    any depth arrives at the time it would have for the reference ray parameter pref."""
    tt = t0 + np.arange(R.shape[1]) * dt
    tnew = np.arange(tmin, tmax + 1e-9, dt)
    aref, bref = vertical_slowness(pref, vp, kappa)
    st = np.zeros_like(tnew)
    for i in range(len(R)):
        a, b = vertical_slowness(p[i], vp, kappa)
        st += np.interp(tnew * (b - a) / (bref - aref), tt, R[i])
    return tnew, st / len(R)


def coastline(ax, lon0, lon1, lat0, lat1, color='0.45', lw=0.7):
    """Coastline from the land mask of the optional package global_land_mask."""
    try:
        from global_land_mask import globe
    except ImportError:
        return
    lat = np.arange(lat0, lat1 + 1e-9, 0.01)
    lon = np.arange(lon0, lon1 + 1e-9, 0.01)
    LO, LA = np.meshgrid(lon, lat)
    ax.contour(
        LO, LA, globe.is_land(LA, LO).astype(float), levels=[0.5], colors=[color], linewidths=lw, zorder=1
    )


def _station_table(a):
    """Coordinates and H-k results of the stations, and the crust table if there is one."""
    import pandas as pd

    hk = pd.read_csv(os.path.join(a.qc, a.rf, 'hk_results.csv')).set_index('station')
    path = a.crust or os.path.join(a.qc, ANALYSIS_DIR, 'crust_model.csv')
    crust = pd.read_csv(path).set_index('station') if os.path.exists(path) else None
    return hk, crust


# -----------------------------------------------------------------------------
# section: station stacks sorted along the network
# -----------------------------------------------------------------------------
def add_section_arguments(p):
    p.add_argument('--rf', default='rf', help='receiver-function folder below --qc (default rf)')
    p.add_argument(
        '--crust',
        default=None,
        help='crust table for the Moho markers (default <qc>/analysis/crust_model.csv)',
    )
    p.add_argument(
        '--sort', default='lon', choices=['lon', 'lat'], help='sort the stations by longitude or latitude'
    )
    p.add_argument('--min-rf', type=int, default=20, help='functions a station needs')
    p.add_argument('--tmax', type=float, default=10.0, help='end of the time axis (s)')
    p.add_argument('--vp', type=float, default=6.3, help='crustal Vp of the moveout correction (km/s)')
    p.add_argument('--kappa', type=float, default=1.75, help='Vp/Vs of the moveout correction')
    p.add_argument('--pref', type=float, default=0.06, help='reference ray parameter (s/km)')


def cmd_section(a):
    import pandas as pd

    plt = pyplot()
    out = output_dir(a)
    hk, crust = _station_table(a)
    rfs = load_station_rfs(os.path.join(a.qc, a.rf), a.min_rf)
    if not rfs:
        sys.exit('no station with %d or more receiver functions in %s' % (a.min_rf, os.path.join(a.qc, a.rf)))
    stas = sorted(rfs, key=lambda s: hk.loc[s, a.sort])
    aref, bref = vertical_slowness(a.pref, a.vp, a.kappa)
    stacks = {}
    fig, ax = plt.subplots(figsize=(max(7.0, 0.27 * len(stas) + 1.5), 5.6))
    gain = 0.9
    for k, sta in enumerate(stas):
        d = rfs[sta]
        t, st = moveout_stack(
            d['rf'], d['t0'], d['dt'], d['p'], a.vp, a.kappa, a.pref, tmin=-1.0, tmax=a.tmax
        )
        st = st / np.max(np.abs(st[(t > -1.0) & (t < 1.0)]))  # relative to the direct pulse
        stacks[sta] = st
        ax.fill_betweenx(t, k, k + gain * st, where=st > 0, color='0.15', lw=0)
        ax.fill_betweenx(t, k, k + gain * st, where=st < 0, color='0.75', lw=0)
        ax.plot(k + gain * st, t, color='k', lw=0.3)
        # Ps time of the Moho: of the depth in the crust table if there is one, else the Ps-only reading
        if crust is not None and sta in crust.index and np.isfinite(crust.loc[sta, 'moho_km']):
            firm = str(crust.loc[sta, 'method']) in ('H-k', 'published')
            ax.plot(
                k - 0.32,
                crust.loc[sta, 'moho_km'] * (bref - aref),
                marker='>',
                ms=5,
                mec='k',
                mfc='k' if firm else 'white',
                lw=0,
            )
        elif np.isfinite(hk.loc[sta, 't_ps_s']):
            ax.plot(k - 0.32, hk.loc[sta, 't_ps_s'], marker='>', ms=5, mec='k', mfc='white', lw=0)
    ax.set_xlim(-1.0, len(stas))
    ax.set_ylim(a.tmax, -1.0)
    ax.set_xticks(range(len(stas)))
    ax.set_xticklabels(stas, rotation=90, fontsize=7)
    ax.set_ylabel('time after P (s), corrected to a ray parameter of %.2f s/km' % a.pref)
    ax.set_title(
        'Station stacks sorted by %s (%s first); triangle: PdS time of the Moho, filled where it is from '
        'H-k or published'
        % ('longitude' if a.sort == 'lon' else 'latitude', 'west' if a.sort == 'lon' else 'south'),
        fontsize=8.5,
        loc='left',
    )
    for y in range(1, int(a.tmax) + 1):
        ax.axhline(y, color='0.9', lw=0.4, zorder=0)
    save_figure(fig, out, 'rf_section.png')
    plt.close(fig)
    tab = pd.DataFrame(stacks)
    tab.insert(0, 'time_s', np.round(t, 3))
    save_csv(tab.round(4), out, 'rf_section_stacks.csv')
    print('%d stations, %d functions' % (len(stas), sum(len(rfs[s]['rf']) for s in stas)))


# -----------------------------------------------------------------------------
# ccp: common-conversion-point stack
# -----------------------------------------------------------------------------
def add_ccp_arguments(p):
    p.add_argument('--rf', default='rf', help='receiver-function folder below --qc (default rf)')
    p.add_argument(
        '--crust',
        default=None,
        help='crust table for the Moho markers (default <qc>/analysis/crust_model.csv)',
    )
    p.add_argument('--min-rf', type=int, default=20, help='functions a station needs')
    p.add_argument(
        '--region',
        type=float,
        nargs=4,
        metavar=('LON0', 'LON1', 'LAT0', 'LAT1'),
        default=None,
        help='grid limits (default: the stations plus 0.4 degrees)',
    )
    p.add_argument('--bin-deg', type=float, default=0.1, help='horizontal bin size (degrees)')
    p.add_argument('--dz', type=float, default=0.5, help='depth step (km)')
    p.add_argument('--zmax', type=float, default=80.0, help='largest depth (km)')
    p.add_argument('--vp-crust', type=float, default=6.3)
    p.add_argument('--kappa-crust', type=float, default=1.75)
    p.add_argument(
        '--moho-model', type=float, default=40.0, help='Moho depth of the 1-D migration model (km)'
    )
    p.add_argument('--vp-mantle', type=float, default=8.04)
    p.add_argument('--kappa-mantle', type=float, default=1.79)
    p.add_argument(
        '--mute-km', type=float, default=3.0, help='depth above which the direct pulse is muted (km)'
    )
    p.add_argument(
        '--profile',
        type=float,
        nargs=4,
        action='append',
        metavar=('LAT0', 'LON0', 'LAT1', 'LON1'),
        default=None,
        help='profile end points; repeat for several (default: west-east and south-north through the network)',
    )
    p.add_argument('--half-width', type=float, default=45.0, help='half width of the profile swath (km)')
    p.add_argument('--step', type=float, default=10.0, help='bin size along a profile (km)')
    p.add_argument(
        '--slices', type=float, nargs='*', default=[25.0, 35.0, 45.0], help='depths of the map slices (km)'
    )
    p.add_argument(
        '--clip', type=float, default=0.10, help='color scale limit (amplitude relative to the direct P)'
    )


def ray_tables(p, z, vp_z, vs_z):
    """Ps delay time (s) and horizontal offset of the S leg (km) at the depths z for a ray
    parameter p (s/km) in the layered model vp_z, vs_z."""
    dz = z[1] - z[0]
    qp = np.sqrt(np.maximum(1.0 / vp_z**2 - p**2, 1e-12))
    qs = np.sqrt(np.maximum(1.0 / vs_z**2 - p**2, 1e-12))
    tps = np.concatenate([[0.0], np.cumsum((qs - qp)[1:] * dz)])
    xs = np.concatenate([[0.0], np.cumsum((p / qs)[1:] * dz)])
    return tps, xs


def ccp_volume(rfs, coords, lon_edges, lat_edges, z, vp_z, vs_z, mute_km=3.0):
    """Sum of the amplitudes and number of traces in every bin of the grid.

    rfs: output of load_station_rfs; coords: {station: (lat, lon)}. The conversion points lie
    on the source side of the station, at the offset of the S leg along the back-azimuth."""
    nlon, nlat, nz = len(lon_edges) - 1, len(lat_edges) - 1, len(z)
    dlon, dlat = lon_edges[1] - lon_edges[0], lat_edges[1] - lat_edges[0]
    vol = np.zeros((nlon, nlat, nz))
    cnt = np.zeros((nlon, nlat, nz))
    iz = np.arange(nz)
    n = 0
    for sta, d in rfs.items():
        slat, slon = coords[sta]
        tt = d['t0'] + np.arange(d['rf'].shape[1]) * d['dt']
        for i in range(len(d['rf'])):
            tps, xs = ray_tables(d['p'][i], z, vp_z, vs_z)
            amp = np.interp(tps, tt, d['rf'][i])
            amp[z < mute_km] = 0.0
            baz = np.radians(d['baz'][i])
            la = slat + xs * np.cos(baz) / KM_PER_DEG
            lo = slon + xs * np.sin(baz) / (KM_PER_DEG * np.cos(np.radians(slat)))
            ilo = np.floor((lo - lon_edges[0]) / dlon).astype(int)
            ila = np.floor((la - lat_edges[0]) / dlat).astype(int)
            ok = (ilo >= 0) & (ilo < nlon) & (ila >= 0) & (ila < nlat)
            np.add.at(vol, (ilo[ok], ila[ok], iz[ok]), amp[ok])
            np.add.at(cnt, (ilo[ok], ila[ok], iz[ok]), 1.0)
            n += 1
    return vol, cnt, n


def profile_section(
    vol, cnt, lon_c, lat_c, z, p0, p1, half_width_km=45.0, step_km=10.0, min_traces=3, mute_km=8.0
):
    """Project the bins within half_width_km of the line p0 to p1 (lat, lon) onto the line.
    Returns the distance of the profile bins (km), the section (distance, depth) and the
    profile length."""
    lat0, lon0 = p0
    lat1, lon1 = p1
    kx = KM_PER_DEG * np.cos(np.radians(0.5 * (lat0 + lat1)))
    ux, uy = (lon1 - lon0) * kx, (lat1 - lat0) * KM_PER_DEG
    length = float(np.hypot(ux, uy))
    ux, uy = ux / length, uy / length
    LO, LA = np.meshgrid(lon_c, lat_c, indexing='ij')
    X = (LO - lon0) * kx
    Y = (LA - lat0) * KM_PER_DEG
    along = X * ux + Y * uy
    across = -X * uy + Y * ux
    sel = (np.abs(across) <= half_width_km) & (along >= -20) & (along <= length + 20)
    s_edges = np.arange(-20, length + 20 + 1e-9, step_km)
    sec = np.full((len(s_edges) - 1, len(z)), np.nan)
    for k in range(len(s_edges) - 1):
        m = sel & (along >= s_edges[k]) & (along < s_edges[k + 1])
        if m.sum() == 0:
            continue
        num = vol[m].sum(axis=0)
        den = cnt[m].sum(axis=0)
        sec[k] = np.where(den >= min_traces, num / np.maximum(den, 1e-9), np.nan)
    # running mean over three bins along the profile; the shallowest part is masked
    smooth = np.full_like(sec, np.nan)
    for k in range(sec.shape[0]):
        blk = sec[[j for j in (k - 1, k, k + 1) if 0 <= j < sec.shape[0]]]
        good = ~np.isnan(blk)
        with np.errstate(invalid='ignore'):
            smooth[k] = np.where(
                good.sum(axis=0) > 0, np.nansum(blk, axis=0) / np.maximum(good.sum(axis=0), 1), np.nan
            )
    smooth[:, z < mute_km] = np.nan
    return 0.5 * (s_edges[:-1] + s_edges[1:]), smooth, length


def _stations_on_profile(coords, p0, p1, half_width_km):
    lat0, lon0 = p0
    lat1, lon1 = p1
    kx = KM_PER_DEG * np.cos(np.radians(0.5 * (lat0 + lat1)))
    ux, uy = (lon1 - lon0) * kx, (lat1 - lat0) * KM_PER_DEG
    length = float(np.hypot(ux, uy))
    ux, uy = ux / length, uy / length
    out = []
    for sta, (lat, lon) in coords.items():
        x, y = (lon - lon0) * kx, (lat - lat0) * KM_PER_DEG
        along, across = x * ux + y * uy, -x * uy + y * ux
        if abs(across) <= half_width_km and -20 <= along <= length + 20:
            out.append((sta, along))
    return out


def cmd_ccp(a):
    from scipy.ndimage import uniform_filter

    plt = pyplot()
    out = output_dir(a)
    hk, crust = _station_table(a)
    rfs = load_station_rfs(os.path.join(a.qc, a.rf), a.min_rf)
    if not rfs:
        sys.exit('no station with %d or more receiver functions in %s' % (a.min_rf, os.path.join(a.qc, a.rf)))
    coords = dict((s, (float(hk.loc[s, 'lat']), float(hk.loc[s, 'lon']))) for s in rfs)
    lats = np.array([c[0] for c in coords.values()])
    lons = np.array([c[1] for c in coords.values()])
    if a.region:
        lon0, lon1, lat0, lat1 = a.region
    else:
        lon0, lon1 = math.floor((lons.min() - 0.4) * 10) / 10, math.ceil((lons.max() + 0.4) * 10) / 10
        lat0, lat1 = math.floor((lats.min() - 0.4) * 10) / 10, math.ceil((lats.max() + 0.4) * 10) / 10
    lon_edges = np.arange(lon0, lon1 + 1e-9, a.bin_deg)
    lat_edges = np.arange(lat0, lat1 + 1e-9, a.bin_deg)
    lon_c = 0.5 * (lon_edges[:-1] + lon_edges[1:])
    lat_c = 0.5 * (lat_edges[:-1] + lat_edges[1:])
    z = np.arange(0.0, a.zmax + 1e-9, a.dz)
    vp_z = np.where(z < a.moho_model, a.vp_crust, a.vp_mantle)
    vs_z = np.where(z < a.moho_model, a.vp_crust / a.kappa_crust, a.vp_mantle / a.kappa_mantle)
    vol, cnt, n = ccp_volume(rfs, coords, lon_edges, lat_edges, z, vp_z, vs_z, a.mute_km)
    print('%d receiver functions of %d stations migrated' % (n, len(rfs)))
    raw = np.where(cnt > 0, vol / np.maximum(cnt, 1), np.nan)
    # light horizontal smoothing over 3 x 3 bins, weighted by the number of traces; a column
    # of 3 x 3 bins needs at least 5 traces
    num = uniform_filter(np.nan_to_num(vol), size=(3, 3, 1), mode='constant')
    den = uniform_filter(cnt, size=(3, 3, 1), mode='constant')
    stack = np.where(den * 9 >= 5, num / np.maximum(den, 1e-9), np.nan)
    path = os.path.join(out, 'ccp_volume.npz')
    np.savez_compressed(path, stack=stack, raw=raw, count=cnt, lon=lon_c, lat=lat_c, z=z)
    print('  wrote', path)

    profiles = a.profile or [
        [float(np.median(lats)), lon0, float(np.median(lats)), lon1],
        [lat0, float(np.median(lons)), lat1, float(np.median(lons))],
    ]
    names = 'ABCDEFGH'
    fig, axs = plt.subplots(len(profiles), 1, figsize=(7.5, 3.1 * len(profiles)), squeeze=False)
    for k, pr in enumerate(profiles):
        p0, p1 = (pr[0], pr[1]), (pr[2], pr[3])
        s_c, sec, length = profile_section(vol, cnt, lon_c, lat_c, z, p0, p1, a.half_width, a.step)
        ax = axs[k, 0]
        im = ax.pcolormesh(s_c, z, sec.T, cmap='RdBu_r', vmin=-a.clip, vmax=a.clip, shading='nearest')
        for sta, along in _stations_on_profile(coords, p0, p1, a.half_width):
            ax.plot(along, -1.5, marker='v', color='k', ms=4, clip_on=False)
            ax.text(along, -4.5, sta[-2:], ha='center', va='bottom', fontsize=5.5, rotation=90)
            if crust is not None and sta in crust.index and np.isfinite(crust.loc[sta, 'moho_km']):
                meth = str(crust.loc[sta, 'method'])
                ax.errorbar(
                    along,
                    crust.loc[sta, 'moho_km'],
                    yerr=crust.loc[sta, 'moho_err_km'],
                    fmt='s' if meth in ('H-k', 'published') else 'o',
                    ms=3.5,
                    mfc='k' if meth == 'H-k' else 'white',
                    mec='k',
                    ecolor='k',
                    lw=0.6,
                    capsize=1.5,
                    zorder=5,
                )
        ax.set_ylim(a.zmax - 5.0, 0)
        ax.set_xlim(-20, length + 20)
        ax.set_ylabel('depth (km)')
        ax.set_xlabel('distance along profile %s (km)' % names[k])
        ax.text(
            0.005,
            0.03,
            '%s: %.2fN %.2fE to %.2fN %.2fE' % (names[k], pr[0], pr[1], pr[2], pr[3]),
            transform=ax.transAxes,
            fontsize=7,
            va='bottom',
            bbox=dict(facecolor='white', edgecolor='none', alpha=0.8, pad=1),
        )
    fig.tight_layout()
    fig.colorbar(
        im, ax=list(axs[:, 0]), shrink=0.6, pad=0.02, label='stacked amplitude (relative to the direct P)'
    )
    save_figure(fig, out, 'ccp_profiles.png')
    plt.close(fig)

    if a.slices:
        fig, axs = plt.subplots(1, len(a.slices), figsize=(3.6 * len(a.slices), 3.4), squeeze=False)
        half = max(int(round(2.0 / a.dz)), 1)
        for ax, z0 in zip(axs[0], a.slices):
            k = int(np.argmin(np.abs(z - z0)))
            with np.errstate(invalid='ignore'):
                band = np.nanmean(stack[:, :, max(k - half, 0) : k + half + 1], axis=2)
            pm = ax.pcolormesh(
                lon_edges,
                lat_edges,
                band.T,
                cmap='RdBu_r',
                vmin=-1.2 * a.clip,
                vmax=1.2 * a.clip,
                shading='flat',
            )
            coastline(ax, lon0, lon1, lat0, lat1)
            ax.plot(lons, lats, '^', color='k', ms=3, mfc='none', mew=0.5)
            for j, pr in enumerate(profiles):
                ax.plot([pr[1], pr[3]], [pr[0], pr[2]], color='k', lw=0.7, ls='--')
                ax.text(pr[1], pr[0], names[j], fontsize=7, fontweight='bold', ha='right')
            ax.set_xlim(lon0, lon1)
            ax.set_ylim(lat0, lat1)
            ax.set_aspect(1.0 / math.cos(math.radians(0.5 * (lat0 + lat1))))
            ax.set_title('%.0f km (+/- 2 km)' % z0, fontsize=8)
            ax.tick_params(labelsize=7)
        fig.colorbar(pm, ax=list(axs[0]), shrink=0.8, label='stacked amplitude')
        save_figure(fig, out, 'ccp_slices.png')
        plt.close(fig)
