"""Manual P picking and its comparison with the automatic picks.

Commands:
    handpick       interactive, blind picking (the automatic pick is never shown)
    compare-picks  manual against automatic picks by quality class, absolute and relative
"""

import glob
import os
import sys

import numpy as np

from .common import (
    CLASS_ORDER,
    output_dir,
    read_tables,
    pyplot,
    save_figure,
    save_csv,
    read_channel,
)


class HandPicker(object):
    """One station at a time: raw vertical, causal band-pass, time relative to the ak135
    prediction. The automatic pick is never shown (blind test)."""

    KEYS = (
        'click = pick | 1 2 3 = quality (clear / fair / unsure) | n or enter = save and next\n'
        'x = cannot pick | b = back | f = filter on/off | q = save and quit | toolbar zoom works too'
    )

    def __init__(self, items, csv_path, picker, band=(0.3, 8.0), win=(-15.0, 20.0), fig=None):
        self.items, self.csv_path, self.picker, self.band, self.win = items, csv_path, picker, band, win
        self.k, self.pick, self.quality, self.filt = 0, None, 1, True
        self._shown = None
        self.done = self._load()
        self.fig = fig
        self._skip_done()

    def _load(self):
        import pandas as pd

        if os.path.exists(self.csv_path):
            return pd.read_csv(self.csv_path).to_dict('records')
        return []

    def _key(self, it):
        return (it['event'], it['station'])

    def _skip_done(self):
        have = set((r['event'], r['station']) for r in self.done)
        while self.k < len(self.items) and self._key(self.items[self.k]) in have:
            self.k += 1

    def save(self):
        import pandas as pd

        pd.DataFrame(
            self.done,
            columns=['event', 'station', 'phase', 'manual_pick_time', 'quality', 'picker', 'filter'],
        ).to_csv(self.csv_path, index=False)

    def record(self, cannot=False):
        from obspy import UTCDateTime

        if self.k >= len(self.items):
            return
        it = self.items[self.k]
        t = '' if (cannot or self.pick is None) else str(UTCDateTime(it['t_pred']) + self.pick)
        self.done = [r for r in self.done if (r['event'], r['station']) != self._key(it)]
        self.done.append(
            dict(
                event=it['event'],
                station=it['station'],
                phase='P',
                manual_pick_time=t,
                quality=(0 if not t else self.quality),
                picker=self.picker,
                filter=('causal %g-%g Hz' % self.band) if self.filt else 'none',
            )
        )
        self.save()
        self.k += 1
        self._skip_done()
        self.pick, self.quality = None, 1

    def on_click(self, x):
        if x is not None:
            self.pick = float(x)

    def on_key(self, key):
        if key in ('1', '2', '3'):
            self.quality = int(key)
        elif key in ('n', 'enter'):
            if self.pick is not None:
                self.record()
        elif key == 'x':
            self.record(cannot=True)
        elif key == 'b':
            self.k = max(0, self.k - 1)
            self.pick = None
        elif key == 'f':
            self.filt = not self.filt
        elif key == 'q':
            self.save()
            return 'quit'
        return None

    def trace(self):
        from obspy import UTCDateTime

        it = self.items[self.k]
        tr = read_channel(it['folder'], it['station'], 'Z')
        if tr is None:
            return None, None
        t0 = UTCDateTime(it['t_pred'])
        tr = tr.slice(t0 + self.win[0] - 60, t0 + self.win[1] + 20)
        tr.detrend('demean')
        tr.detrend('linear')
        tr.taper(0.05)
        if self.filt:
            tr.filter('bandpass', freqmin=self.band[0], freqmax=self.band[1], corners=2, zerophase=False)
        tt = tr.times() + float(tr.stats.starttime - t0)
        keep = (tt >= self.win[0]) & (tt <= self.win[1])
        y = tr.data[keep]
        return tt[keep], y / (np.max(np.abs(y)) + 1e-30)

    def draw(self):
        ax = self.fig.axes[0]
        same = self._shown == (self.k, self.filt)  # same trace: keep the current zoom
        lims = (ax.get_xlim(), ax.get_ylim()) if same else None
        self._shown = (self.k, self.filt)
        ax.clear()
        if self.k >= len(self.items):
            ax.text(0.5, 0.5, 'All done - press q', ha='center', transform=ax.transAxes)
            self.fig.canvas.draw_idle()
            return
        it = self.items[self.k]
        t, y = self.trace()
        if t is not None:
            ax.plot(t, y, 'k', lw=0.7)
        if self.pick is not None:
            ax.axvline(self.pick, color='r', lw=1)
        ax.set_xlim(self.win)
        if lims is not None:
            ax.set_xlim(lims[0])
            ax.set_ylim(lims[1])
        ax.set_xlabel('time relative to the ak135 P prediction (s)')
        ax.set_title(
            '%s   %s   %.1f deg   [%d / %d]   quality %d   filter %s\n%s'
            % (
                it['event'],
                it['station'],
                it['dist'],
                self.k + 1,
                len(self.items),
                self.quality,
                'on' if self.filt else 'off',
                self.KEYS,
            ),
            fontsize=8,
        )
        self.fig.canvas.draw_idle()


def _handpick_items(a):
    ev, sq = read_tables(a.qc)
    acc = ev[ev.status == 'accepted']
    if a.events:
        names = sorted(
            set(acc.event) & set(os.path.basename(f) for f in glob.glob(os.path.join(a.data, a.events)))
        )
    else:
        rich = sq[sq.usable_P.fillna(False).astype(bool)].groupby('event').size()
        cand = sorted(set(acc.event) & set(rich[rich >= 10].index))
        names = (
            sorted(
                str(x)
                for x in np.random.RandomState(a.seed).choice(cand, min(a.n_events, len(cand)), replace=False)
            )
            if cand
            else []
        )
    items = []
    for evn in names:
        g = sq[(sq.event == evn) & sq.usable_P.fillna(False).astype(bool)].sort_values('dist')
        items += [
            dict(
                event=evn, station=r.station, dist=r.dist, t_pred=r.tP_pred, folder=os.path.join(a.data, evn)
            )
            for r in g.itertuples()
        ]
    return names, items


def cmd_handpick(a):
    import matplotlib.pyplot as plt

    out = output_dir(a)
    names, items = _handpick_items(a)
    if not items:
        sys.exit('no events to pick')
    csv_path = os.path.join(out, 'manual_picks_%s.csv' % a.picker)
    print('%d events, %d stations to pick -> %s' % (len(names), len(items), csv_path))
    print('events: ' + ', '.join(names))
    print('Blind test: do not open the QC figures of these events before the picking is finished.')
    for k in ('keymap.fullscreen', 'keymap.quit'):  # f and q are keys of the picker
        if k in plt.rcParams:
            plt.rcParams[k] = []
    fig = plt.figure(figsize=(12, 4.5))
    fig.add_subplot(1, 1, 1)
    hp = HandPicker(items, csv_path, a.picker, fig=fig)

    def click(e):
        tb = getattr(fig.canvas, 'toolbar', None)
        if e.inaxes is not None and e.button == 1 and (tb is None or not tb.mode):
            hp.on_click(e.xdata)
            hp.draw()

    def key(e):
        if hp.on_key(e.key) == 'quit':
            plt.close(fig)
            return
        hp.draw()

    fig.canvas.mpl_connect('button_press_event', click)
    fig.canvas.mpl_connect('key_press_event', key)
    hp.draw()
    plt.show()
    hp.save()
    print('saved %d picks to %s' % (len(hp.done), csv_path))


# -----------------------------------------------------------------------------
# compare-picks: manual against automatic
# -----------------------------------------------------------------------------
def cmd_compare_picks(a):
    import pandas as pd
    from obspy import UTCDateTime

    plt = pyplot()
    out = output_dir(a)
    files = a.manual or sorted(glob.glob(os.path.join(out, 'manual_picks_*.csv')))
    if not files:
        sys.exit('no manual pick files (run handpick first)')
    man = pd.concat([pd.read_csv(f) for f in files], ignore_index=True)
    man = man[man.manual_pick_time.notna() & (man.manual_pick_time.astype(str) != '')]
    auto = pd.read_csv(os.path.join(a.qc, 'tables', 'picks_all.csv'))
    auto = auto[auto.phase == 'P'][['event', 'station', 'pick_time', 'quality']]
    m = man.merge(auto, on=['event', 'station'], how='inner')
    if not len(m):
        sys.exit('no station-events with both a manual and an automatic pick')
    m['diff_s'] = [UTCDateTime(x) - UTCDateTime(y) for x, y in zip(m.manual_pick_time, m.pick_time)]
    # relative comparison: remove each event's mean difference (what relative tomography sees)
    m['diff_rel_s'] = m.diff_s - m.groupby(['picker', 'event']).diff_s.transform('mean')
    save_csv(m.round(4), out, 'pick_comparison.csv')
    rows = []
    for cls, g in [('all A-C', m[m.quality_y.isin(['A', 'B', 'C'])])] + [
        (c, m[m.quality_y == c]) for c in CLASS_ORDER
    ]:
        if len(g):
            rows.append(
                dict(
                    auto_class=cls,
                    n=len(g),
                    mean_diff_s=g.diff_s.mean(),
                    median_diff_s=g.diff_s.median(),
                    std_diff_s=g.diff_s.std(),
                    std_relative_s=g.diff_rel_s.std(),
                    within_0p1s_pct=100.0 * np.mean(np.abs(g.diff_rel_s) <= 0.1),
                )
            )
    s = pd.DataFrame(rows).round(3)
    save_csv(s, out, 'pick_comparison_summary.csv')
    print(s.to_string(index=False))
    fig, axs = plt.subplots(1, 2, figsize=(11, 3.8))
    bins = np.linspace(-1, 1, 41)
    for c, col in zip(CLASS_ORDER, ('k', '0.35', '0.6', '0.85')):
        g = m[m.quality_y == c]
        if len(g):
            axs[0].hist(
                g.diff_s.clip(-1, 1), bins=bins, color=col, alpha=0.9, label='class %s (%d)' % (c, len(g))
            )
            axs[1].hist(g.diff_rel_s.clip(-1, 1), bins=bins, color=col, alpha=0.9, label='class %s' % c)
    axs[0].set_xlabel('manual minus automatic, absolute (s)')
    axs[1].set_xlabel('manual minus automatic, event mean removed (s)')
    for ax in axs:
        ax.axvline(0, color='k', lw=0.6)
        ax.legend(frameon=False, fontsize=7)
        ax.set_ylabel('picks')
    axs[0].set_title('Absolute times')
    axs[1].set_title('Relative times')
    plt.tight_layout()
    save_figure(fig, out, 'pick_comparison.png')
    plt.close(fig)
