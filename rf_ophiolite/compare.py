"""Two receiver-function runs of the same data compared.

Typical use: the same data deconvolved with two Gaussian widths, or with and without
the station screen. Both runs are folders written by the rf command below --qc.

    python rfo.py compare-rf --qc qc_out --a rf --b rf_g25 --label-a "width 2.0" --label-b "width 2.5"

Written to <qc>/analysis:
    rf_runs_counts.csv       functions computed, accepted by the rules, removed by the
                             station screen, kept; stations with enough functions;
                             H-k quality classes
    rf_runs_by_station.csv   Ps time, H-k depth, Vp/Vs and quality class of every station
                             in the two runs, and the differences
    rf_runs_compare.png      station stacks of the two runs, Ps time and H-k depth of one
                             run against the other
"""

import os
import sys

import numpy as np

from .common import output_dir, pyplot, save_csv, save_figure
from .imaging import moveout_stack


def add_arguments(p):
    p.add_argument('--a', required=True, help='first receiver-function folder below --qc, e.g. rf')
    p.add_argument('--b', required=True, help='second receiver-function folder below --qc, e.g. rf_g25')
    p.add_argument('--label-a', default=None, help='name of the first run in tables and figure')
    p.add_argument('--label-b', default=None, help='name of the second run')
    p.add_argument('--min-rf', type=int, default=20, help='functions a station needs')
    p.add_argument(
        '--stations',
        default='',
        help='stations whose stacks are drawn, comma separated '
        '(default: the three with the most functions)',
    )


def _run_counts(folder, label, min_rf):
    import pandas as pd

    s = pd.read_csv(os.path.join(folder, 'rf_summary.csv'))
    ok = s[s.ok.astype(str) == 'True']
    rules = s[s.ok_global.astype(str) == 'True'] if 'ok_global' in s else ok
    per_sta = ok.groupby('station').size()
    hk = pd.read_csv(os.path.join(folder, 'hk_results.csv'))
    hk = hk[hk.n_rf >= min_rf] if 'H_km' in hk else hk.iloc[0:0]
    row = dict(
        run=label,
        computed=len(s),
        accepted_by_rules=len(rules),
        removed_by_screen=len(rules) - len(ok),
        kept=len(ok),
        kept_pct=round(100.0 * len(ok) / max(len(s), 1), 1),
        stations_with_enough=int((per_sta >= min_rf).sum()),
        median_per_station=(
            float(per_sta[per_sta >= min_rf].median()) if (per_sta >= min_rf).any() else np.nan
        ),
        quality_A=int((hk.quality == 'A').sum()),
        quality_B=int((hk.quality == 'B').sum()),
        quality_C=int((hk.quality == 'C').sum()),
    )
    return row, hk.set_index('station')


def cmd_compare_rf(a):
    import pandas as pd

    plt = pyplot()
    out = output_dir(a)
    la, lb = a.label_a or a.a, a.label_b or a.b
    fa, fb = os.path.join(a.qc, a.a), os.path.join(a.qc, a.b)
    for f in (fa, fb):
        if not os.path.exists(os.path.join(f, 'hk_results.csv')):
            sys.exit(
                '%s has no hk_results.csv (run the rf command with --out %s first)' % (f, os.path.basename(f))
            )
    ra, ha = _run_counts(fa, la, a.min_rf)
    rb, hb = _run_counts(fb, lb, a.min_rf)
    counts = pd.DataFrame([ra, rb])
    save_csv(counts, out, 'rf_runs_counts.csv')
    print(counts.to_string(index=False))

    both = ha.index.intersection(hb.index)
    cols = ['n_rf', 't_ps_s', 'H_km', 'H_err_km', 'kappa', 'quality']
    t = ha.loc[both, cols].add_suffix('_a').join(hb.loc[both, cols].add_suffix('_b'))
    t['d_t_ps_s'] = (t.t_ps_s_b - t.t_ps_s_a).round(2)
    t['d_H_km'] = (t.H_km_b - t.H_km_a).round(2)
    t['d_kappa'] = (t.kappa_b - t.kappa_a).round(3)
    t['same_quality'] = t.quality_a == t.quality_b
    save_csv(t.reset_index(), out, 'rf_runs_by_station.csv')
    far = list(t.index[t.d_t_ps_s.abs() > 1.0])
    print('%d stations in both runs (a = %s, b = %s)' % (len(both), la, lb))
    print(
        'Ps time: median difference %.2f s, 90 %% of the stations within %.2f s; more than 1 s at: %s'
        % (t.d_t_ps_s.abs().median(), t.d_t_ps_s.abs().quantile(0.9), ', '.join(far) or 'none')
    )
    print(
        'H-k depth: median difference %.2f km; Vp/Vs: %.3f; same quality class at %d of %d stations'
        % (t.d_H_km.abs().median(), t.d_kappa.abs().median(), int(t.same_quality.sum()), len(t))
    )

    show = [s.strip() for s in a.stations.split(',') if s.strip()] or list(
        t.sort_values('n_rf_a', ascending=False).index[:3]
    )
    fig, axs = plt.subplots(1, 3, figsize=(12.5, 4.2), gridspec_kw=dict(width_ratios=[1.35, 1, 1]))
    ax = axs[0]
    for k, sta in enumerate(show):
        for folder, col, lw in ((fa, 'k', 1.5), (fb, '#b2182b', 1.0)):
            f = os.path.join(folder, 'rf_%s.npz' % sta)
            if not os.path.exists(f):
                continue
            d = np.load(f, allow_pickle=True)
            tt, st = moveout_stack(
                d['rf'].astype(float),
                float(d['t0']),
                float(d['dt']),
                d['rayp'].astype(float),
                tmin=-2.0,
                tmax=12.0,
            )
            st = st / np.max(np.abs(st[(tt > -1.0) & (tt < 1.0)]))
            ax.plot(tt, 0.85 * st + (len(show) - 1 - k) * 1.25, color=col, lw=lw)
        ax.text(11.8, (len(show) - 1 - k) * 1.25 + 0.3, sta, ha='right', fontsize=9)
        ax.axhline((len(show) - 1 - k) * 1.25, color='0.8', lw=0.5, zorder=0)
    ax.plot([], [], color='k', lw=1.5, label=la)
    ax.plot([], [], color='#b2182b', lw=1.0, label=lb)
    ax.legend(frameon=False, fontsize=8, loc='upper right', ncol=2)
    ax.set_xlim(-2, 12)
    ax.set_ylim(-0.75, len(show) * 1.25 + 0.5)
    ax.set_yticks([])
    ax.set_xlabel('time after P (s)')
    ax.set_title('(a) station stacks', loc='left', fontsize=9)
    for ax, col, unit, title in (
        (axs[1], 't_ps_s', 's', '(b) PdS time'),
        (axs[2], 'H_km', 'km', '(c) H-k depth (black: quality A or B)'),
    ):
        x, y = t[col + '_a'], t[col + '_b']
        lo, hi = min(x.min(), y.min()), max(x.max(), y.max())
        pad = 0.06 * (hi - lo)
        ax.plot([lo - pad, hi + pad], [lo - pad, hi + pad], color='0.5', lw=0.8, ls='--')
        good = t.quality_a.isin(['A', 'B']) if col == 'H_km' else np.zeros(len(t), bool)
        ax.plot(x[~good], y[~good], 'o', mfc='white', mec='k' if col == 't_ps_s' else '0.45', ms=6, mew=1.0)
        ax.plot(x[good], y[good], 's', mfc='k', mec='k', ms=6.5)
        if col == 't_ps_s':
            for sta in far:
                ax.annotate(sta, (x[sta], y[sta]), textcoords='offset points', xytext=(6, 4), fontsize=7)
        ax.set_xlim(lo - pad, hi + pad)
        ax.set_ylim(lo - pad, hi + pad)
        ax.set_aspect('equal')
        ax.set_xlabel('%s (%s)' % (la, unit))
        ax.set_ylabel('%s (%s)' % (lb, unit))
        ax.set_title('%s, %d stations' % (title, len(t)), loc='left', fontsize=9)
    fig.tight_layout()
    save_figure(fig, out, 'rf_runs_compare.png')
    plt.close(fig)
