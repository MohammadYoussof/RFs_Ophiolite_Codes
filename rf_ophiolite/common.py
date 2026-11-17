"""Helpers shared by the analysis commands.

Every analysis command reads the output folder of the quality control (--qc) and
writes its tables and figures to <qc>/analysis (or <qc>/analysis/<subdir> for a
robustness run). The original data folders are only read.
"""

import glob
import os
import sys

import numpy as np

# Noise spectra are stored on the frequency grid of the quality control
# (psd_fmin and psd_fmax in the CONFIG block of qc_pipeline.py).
PSD_FMIN, PSD_FMAX = 0.05, 30.0

# Pick classes, best first. Class D picks are kept in the tables but never used.
CLASS_ORDER = ('A', 'B', 'C', 'D')

# Folder below the quality-control output that receives the analysis results.
ANALYSIS_DIR = 'analysis'

# Band and instrument code of the channels in the raw event folders (HHZ, HHN, HHE).
CHANNEL_BAND = 'HH'


def is_inside(path, folder):
    """True if path is folder itself or below it (device and inode comparison)."""
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


def output_dir(a):
    """Folder for the results of an analysis command; created if missing."""
    out = os.path.join(a.qc, ANALYSIS_DIR, getattr(a, 'subdir', '') or '')
    if getattr(a, 'data', None) and not glob.glob(os.path.join(a.data, 'Event_*')):
        sys.exit('ERROR: no Event_* folders in %s (give the folder with the raw SAC event folders)' % a.data)
    if getattr(a, 'data', None) and (is_inside(out, a.data) or is_inside(a.data, out)):
        sys.exit('ERROR: the output folder %s must be outside the data folder %s' % (out, a.data))
    if not os.path.isdir(out):
        os.makedirs(out)
    return out


def read_tables(qc):
    """Event table and station-event table of the quality control."""
    import pandas as pd

    ev = pd.read_csv(os.path.join(qc, 'tables', 'event_summary.csv'))
    sq = pd.read_csv(os.path.join(qc, 'tables', 'station_event_qc.csv'))
    return ev, sq


def pyplot():
    import matplotlib.pyplot as plt

    plt.rcParams.update({'font.size': 9, 'axes.titlesize': 10})
    return plt


def save_figure(fig, out, name, dpi=150):
    p = os.path.join(out, name)
    fig.savefig(p, dpi=dpi, bbox_inches='tight')
    print('  wrote', p)


def save_csv(df, out, name):
    p = os.path.join(out, name)
    df.to_csv(p, index=False)
    print('  wrote', p)


def wrap180(x):
    """Angle or array of angles wrapped to -180 ... 180 degrees."""
    return (np.asarray(x, float) + 180.0) % 360.0 - 180.0


def read_channel(folder, sta, comp):
    """One channel of one station from a raw event folder, pieces merged. Read only."""
    from obspy import read, Stream

    st = Stream()
    for f in sorted(glob.glob(os.path.join(folder, '*.%s.*.%s%s*.sac' % (sta, CHANNEL_BAND, comp)))):
        st += read(f)
    if not len(st):
        return None
    st.merge(method=1, fill_value=0)
    tr = st[0]
    tr.data = np.asarray(tr.data, np.float64)
    return tr


def processed_file(folder, sta, comp):
    """Path of a processed trace (component Z, R or T) of a station, or None."""
    hits = sorted(glob.glob(os.path.join(folder, '*.%s.*.??%s.sac' % (sta, comp))))
    return hits[0] if hits else None
