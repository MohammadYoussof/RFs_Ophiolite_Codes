"""Draw the workflow chart of the package (docs/workflow.png, .svg and .pdf).

    python docs/draw_workflow.py

Edit the texts below and run it again to change the chart. Needs matplotlib only.
"""

import os

import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch, Polygon, Rectangle

HERE = os.path.dirname(os.path.abspath(__file__))
W, H = 13.0, 21.55  # size of the chart in inches

INK = '#1f1f1f'
SOFT = '#5c5c5c'
LINE = '#3a3a3a'
STAGES = {  # fill of the boxes, color of the stage band
    'data': ('#e6ebe8', '#55645d'),
    'qc': ('#e8ecda', '#66763f'),
    'tt': ('#f2ead7', '#8a7440'),
    'rf': ('#dfe7ea', '#4a6675'),
    'crust': ('#f1e2da', '#8c5744'),
}
FILE_FILL = '#fbfaf6'
MONO = 'DejaVu Sans Mono'

plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 10})
fig = plt.figure(figsize=(W, H))
ax = fig.add_axes([0, 0, 1, 1])
ax.set_xlim(0, W)
ax.set_ylim(0, H)
ax.axis('off')


def step(x, y, w, h, stage, number, title, commands=(), lines=(), writes=None, dashed=False):
    """A process box: number badge, title, commands in monospace, short lines, files written."""
    fill, edge = STAGES[stage]
    ax.add_patch(
        FancyBboxPatch(
            (x, y),
            w,
            h,
            boxstyle='round,pad=0,rounding_size=0.12',
            fc=fill,
            ec=edge,
            lw=1.3,
            ls='--' if dashed else '-',
            zorder=2,
        )
    )
    tx = x + 0.2
    if number is not None:
        ax.add_patch(plt.Circle((x + 0.36, y + h - 0.36), 0.22, fc=edge, ec='none', zorder=3))
        ax.text(
            x + 0.36,
            y + h - 0.37,
            str(number),
            color='white',
            fontsize=11,
            fontweight='bold',
            ha='center',
            va='center',
            zorder=4,
        )
        tx = x + 0.72
    ax.text(tx, y + h - 0.36, title, fontsize=11.5, fontweight='bold', color=INK, va='center', zorder=4)
    yy = y + h - 0.8
    for cmd in ([commands] if isinstance(commands, str) else commands):
        ax.text(x + 0.2, yy, cmd, fontsize=9.3, family=MONO, color=INK, va='center', zorder=4)
        yy -= 0.3
    yy -= 0.04
    for ln in lines:
        ax.text(x + 0.2, yy, ln, fontsize=9.2, color=INK, va='center', zorder=4)
        yy -= 0.29
    if writes:
        ax.text(
            x + 0.2,
            y + 0.21,
            'writes  ' + writes,
            fontsize=8.4,
            family=MONO,
            color=SOFT,
            va='center',
            zorder=4,
        )


def files(x, y, w, h, title, lines=()):
    """A data box: folders and files."""
    ax.add_patch(Rectangle((x, y), w, h, fc=FILE_FILL, ec='#8a8a8a', lw=1.0, zorder=2))
    ax.add_patch(
        Polygon(
            [[x + w - 0.22, y + h], [x + w, y + h - 0.22], [x + w - 0.22, y + h - 0.22]],
            closed=True,
            fc='#d9d7cf',
            ec='#8a8a8a',
            lw=0.8,
            zorder=3,
        )
    )
    ax.text(x + 0.15, y + h - 0.28, title, fontsize=9.6, fontweight='bold', color=INK, va='center', zorder=4)
    yy = y + h - 0.62
    for ln in lines:
        ax.text(x + 0.15, yy, ln, fontsize=8.5, family=MONO, color=SOFT, va='center', zorder=4)
        yy -= 0.27


def decision(cx, cy, w, h, text):
    """A decision diamond."""
    ax.add_patch(
        Polygon(
            [[cx - w / 2, cy], [cx, cy + h / 2], [cx + w / 2, cy], [cx, cy - h / 2]],
            closed=True,
            fc='white',
            ec=LINE,
            lw=1.3,
            zorder=2,
        )
    )
    ax.text(cx, cy, text, fontsize=9.2, color=INK, ha='center', va='center', zorder=4, linespacing=1.25)


def chain(x, y, w, h, texts, stage):
    """Sub-steps inside a box, in a row from x to x + w, joined by small arrows."""
    gap = 0.22
    pw = (w - gap * (len(texts) - 1)) / len(texts)
    for k, t in enumerate(texts):
        x0 = x + k * (pw + gap)
        ax.add_patch(
            FancyBboxPatch(
                (x0, y),
                pw,
                h,
                boxstyle='round,pad=0,rounding_size=0.18',
                fc='white',
                ec=STAGES[stage][1],
                lw=0.9,
                zorder=3,
            )
        )
        ax.text(
            x0 + pw / 2,
            y + h / 2,
            t,
            fontsize=8.6,
            color=INK,
            ha='center',
            va='center',
            zorder=4,
            linespacing=1.15,
        )
        if k:
            arrow([(x0 - gap, y + h / 2), (x0, y + h / 2)], lw=1.0, head=9)


def arrow(
    points, label=None, label_at=None, color=LINE, lw=1.4, dashed=False, ha='center', head=13, bg='white'
):
    """An arrow through the given points (the head is at the last one)."""
    xs, ys = [p[0] for p in points], [p[1] for p in points]
    if len(points) > 2:
        ax.plot(
            xs[:-1], ys[:-1], color=color, lw=lw, ls='--' if dashed else '-', solid_capstyle='butt', zorder=1
        )
    ax.add_patch(
        FancyArrowPatch(
            points[-2],
            points[-1],
            arrowstyle='-|>',
            mutation_scale=head,
            color=color,
            lw=lw,
            ls='--' if dashed else '-',
            shrinkA=0,
            shrinkB=0,
            zorder=5,
        )
    )
    if label:
        ax.text(
            label_at[0],
            label_at[1],
            label,
            fontsize=8.8,
            color=SOFT,
            ha=ha,
            va='center',
            style='italic',
            zorder=6,
            bbox=dict(fc=bg, ec='none', pad=1.5),
        )


def band(y0, y1, stage, text):
    """The colored band on the left that names a stage."""
    ax.add_patch(Rectangle((0.25, y0), 0.42, y1 - y0, fc=STAGES[stage][1], ec='none', zorder=2))
    ax.text(
        0.46,
        0.5 * (y0 + y1),
        text,
        color='white',
        fontsize=10.5,
        fontweight='bold',
        rotation=90,
        ha='center',
        va='center',
        zorder=3,
    )


# ----------------------------------------------------------------------------- title and legend
ax.text(
    0.25, H - 0.42, 'RFs_Ophiolite_Codes: workflow', fontsize=17, fontweight='bold', color=INK, va='center'
)
ax.text(
    0.25,
    H - 0.86,
    'From the data center to receiver-function images and crust-corrected delays. '
    'Every numbered step is one command:  python rfo.py <command>',
    fontsize=10,
    color=SOFT,
    va='center',
)
lx, ly = 7.55, H - 0.42
ax.add_patch(
    FancyBboxPatch(
        (lx, ly - 0.12),
        0.42,
        0.24,
        boxstyle='round,pad=0,rounding_size=0.06',
        fc=STAGES['qc'][0],
        ec=STAGES['qc'][1],
        lw=1.0,
    )
)
ax.text(lx + 0.52, ly, 'step', fontsize=8.8, color=SOFT, va='center')
ax.add_patch(
    Polygon(
        [[lx + 1.05, ly], [lx + 1.27, ly + 0.14], [lx + 1.49, ly], [lx + 1.27, ly - 0.14]],
        closed=True,
        fc='white',
        ec=LINE,
        lw=1.0,
    )
)
ax.text(lx + 1.59, ly, 'decision', fontsize=8.8, color=SOFT, va='center')
ax.add_patch(Rectangle((lx + 2.4, ly - 0.12), 0.36, 0.24, fc=FILE_FILL, ec='#8a8a8a', lw=0.9))
ax.text(lx + 2.86, ly, 'files', fontsize=8.8, color=SOFT, va='center')
ax.plot([lx + 3.4, lx + 3.82], [ly, ly], color=LINE, lw=1.2, ls='--')
ax.text(lx + 3.92, ly, 'optional or repeated', fontsize=8.8, color=SOFT, va='center')

# ----------------------------------------------------------------------------- stage 1: data
band(17.2, 20.4, 'data', 'DATA')
files(0.95, 19.05, 2.55, 1.05, 'Data center', ['FDSN web services:', 'events, stations, data'])
step(
    4.1,
    18.75,
    4.3,
    1.65,
    'data',
    1,
    'Download with SOD',
    'cd data;  sod -f ../sod/cool_5H_full.xml',
    ['events by magnitude and distance,', 'one time window per station'],
)
files(
    9.05,
    18.75,
    3.6,
    1.65,
    'Original event folders',
    ['data/raw_sac/Event_*/*.sac', 'data/raw_mseed/Event_*/*.mseed', 'stations.xml (StationXML)'],
)
arrow([(3.5, 19.57), (4.1, 19.57)])
arrow([(8.4, 19.57), (9.05, 19.57)])
step(
    4.1,
    17.2,
    4.3,
    1.3,
    'data',
    2,
    'Protect the originals',
    'rfo.py originals make',
    writes='originals_sha256.txt',
)
arrow(
    [(10.85, 18.75), (10.85, 17.85), (8.4, 17.85)],
    'read only from here on;\nfingerprint of every file',
    (10.85, 18.28),
)

# ----------------------------------------------------------------------------- stage 2: quality control
band(10.3, 16.8, 'qc', 'QUALITY CONTROL')
arrow([(6.25, 17.2), (6.25, 16.8)])
step(
    0.95,
    14.2,
    8.3,
    2.6,
    'qc',
    3,
    'Quality control and picking',
    'rfo.py qc --data data/raw_sac --inventory stations.xml --out qc_out',
    writes='tables/  picks/  processed/  noise/  figures/',
)
chain(
    1.15,
    14.98,
    7.9,
    0.58,
    [
        'integrity\nchecks',
        'instrument\nresponse',
        'noise\nspectra',
        'signal-to-\nnoise ratio',
        'P and S\npicks',
        'processed\nZ, R, T',
    ],
    'qc',
)
ax.text(
    1.15,
    14.7,
    'picks: cross-correlation and stacking (ICCS), then multi-channel cross-correlation (MCCC); classes A to D',
    fontsize=8.4,
    color=SOFT,
    va='center',
    zorder=4,
)
decision(10.95, 15.35, 2.7, 1.5, '4 or more\ngood P picks?')
arrow([(9.25, 15.35), (9.6, 15.35)])
arrow([(10.95, 16.1), (10.95, 16.38)], color=SOFT, lw=1.1, head=10)
ax.text(
    10.95,
    16.55,
    'no:  event rejected, kept in the tables',
    fontsize=8.8,
    color=SOFT,
    ha='center',
    va='center',
    style='italic',
)
arrow([(10.95, 14.6), (10.95, 13.85)], 'yes:  accepted', (11.75, 14.25))

step(
    0.95,
    12.55,
    4.0,
    1.3,
    'qc',
    None,
    'Check the originals',
    'rfo.py originals verify',
    ['must report: unchanged'],
)
step(
    5.25,
    12.55,
    4.0,
    1.3,
    'qc',
    4,
    'Overview of the data set',
    'inventory  snr  response  noise-hour',
    writes='analysis/*.png  *.csv',
)
files(
    9.6,
    12.55,
    3.05,
    1.3,
    'qc_out/',
    ['picks/picks_absolute.csv', 'picks/picks_relative.csv', 'processed/<event>/*.sac'],
)
arrow([(2.95, 14.2), (2.95, 13.85)])
arrow([(7.25, 14.2), (7.25, 13.85)])

arrow([(2.95, 12.55), (2.95, 12.25)])
step(
    0.95,
    10.3,
    4.0,
    1.95,
    'qc',
    5,
    'Sensor orientation',
    'rfo.py orientation',
    ['P particle motion against', 'the expected back-azimuth'],
    writes='orientation_by_station.csv',
)
decision(6.75, 11.25, 2.9, 1.3, 'stable offset of\n5 degrees or more?')
arrow([(4.95, 11.25), (5.3, 11.25)])
step(
    9.0,
    10.6,
    3.65,
    1.3,
    'qc',
    None,
    'Correct the metadata',
    'rfo.py orientation-apply',
    writes='stations_oriented.xml',
)
arrow([(8.2, 11.25), (9.0, 11.25)], 'yes', (8.6, 11.43))
arrow(
    [(12.65, 11.25), (12.88, 11.25), (12.88, 17.0), (8.6, 17.0), (8.6, 16.8)],
    'run step 3 again with the corrected StationXML',
    (10.75, 17.0),
    dashed=True,
)
arrow([(6.75, 10.6), (6.75, 10.0)], 'no', (7.0, 10.3))

# ----------------------------------------------------------------------------- stage 3: analysis
band(5.0, 9.8, 'rf', 'ANALYSIS')
ax.plot([2.95, 8.95], [10.0, 10.0], color=LINE, lw=1.4, zorder=1)  # from 'no' to the two branches
arrow([(2.95, 10.0), (2.95, 9.8)])
arrow([(8.95, 10.0), (8.95, 9.8)])
step(
    0.95,
    7.55,
    4.0,
    2.25,
    'tt',
    6,
    'Station terms',
    'rfo.py station-terms',
    [
        'residual = station term + event term',
        'bootstrap errors, timing suspects',
        'delays by back-azimuth sector',
    ],
    writes='station_terms_P.csv',
)
step(
    5.25,
    7.55,
    7.4,
    2.25,
    'rf',
    7,
    'Receiver functions',
    'rfo.py rf --gauss 2.0 --out rf',
    writes='rf/hk_results.csv  rf/RF_<station>.png  rf/rf_<station>.npz',
)
chain(
    5.45,
    8.05,
    7.0,
    0.62,
    ['iterative\ndeconvolution', 'acceptance\nrules', 'station\nscreen', 'H-k stack and\nPs-only depth'],
    'rf',
)
step(
    5.25,
    5.0,
    3.55,
    2.2,
    'rf',
    8,
    'Images',
    ['rfo.py section', 'rfo.py ccp'],
    ['station stacks along the network,', 'common-conversion-point stack'],
    writes='rf_section.png  ccp_*',
)
step(
    9.1,
    5.0,
    3.55,
    2.2,
    'rf',
    None,
    'Compare two runs',
    'rfo.py compare-rf --a rf --b rf_g25',
    ['for example Gaussian width', '2.0 against 2.5'],
    writes='rf_runs_*.csv  .png',
    dashed=True,
)
arrow([(7.0, 7.55), (7.0, 7.2)])
arrow([(10.9, 7.55), (10.9, 7.2)], dashed=True)

# ----------------------------------------------------------------------------- stage 4: crust
band(1.95, 4.7, 'crust', 'CRUST')
arrow([(8.95, 7.55), (8.95, 4.7)])
step(
    5.25,
    1.95,
    7.4,
    2.75,
    'crust',
    9,
    'Crust table',
    'rfo.py crust-from-rf --published data/moho_published.csv',
    writes='analysis/crust_model.csv',
)
decision(6.6, 3.15, 2.2, 1.1, 'H-k quality\nA or B?')
decision(9.35, 3.15, 2.2, 1.1, 'published\ndepth?')
ax.text(
    11.62,
    3.15,
    'Ps-only depth,\nuncertainty\n5 km or more',
    fontsize=8.6,
    color=INK,
    ha='center',
    va='center',
    zorder=4,
    linespacing=1.2,
)
arrow([(7.7, 3.15), (8.25, 3.15)], 'no', (7.97, 3.35), bg=STAGES['crust'][0])
arrow([(10.45, 3.15), (10.95, 3.15)], 'no', (10.7, 3.35), bg=STAGES['crust'][0])
ax.text(
    6.6, 2.44, 'yes:  H-k depth', fontsize=8.6, color=SOFT, ha='center', va='center', style='italic', zorder=4
)
ax.text(
    9.35,
    2.44,
    'yes:  published depth',
    fontsize=8.6,
    color=SOFT,
    ha='center',
    va='center',
    style='italic',
    zorder=4,
)
step(
    0.95,
    1.95,
    4.0,
    2.75,
    'crust',
    10,
    'Crustal correction',
    'rfo.py crust --model crust_model.csv',
    [
        'delay predicted for the crust',
        'below every station;',
        'station term minus crust',
        '= mantle residual',
    ],
    writes='crust_correction.csv',
)
arrow([(2.95, 7.55), (2.95, 4.7)], 'station terms', (2.95, 6.1))
arrow([(5.25, 3.3), (4.95, 3.3)])

# ----------------------------------------------------------------------------- products
ax.text(0.25, 1.55, 'PRODUCTS', fontsize=10.5, fontweight='bold', color=INK, va='center')
prod = [
    ('P and S travel times\nwith uncertainties', 'steps 3 and 6'),
    ('Moho depth and Vp/Vs\nat every station', 'steps 7 and 9'),
    ('Section and CCP image\nof the crust', 'step 8'),
    ('Crust-corrected\nmantle delays', 'step 10'),
]
pw = (W - 0.5 - 3 * 0.25) / 4
for k, (t, sub) in enumerate(prod):
    x = 0.25 + k * (pw + 0.25)
    ax.add_patch(
        FancyBboxPatch(
            (x, 0.3), pw, 1.0, boxstyle='round,pad=0,rounding_size=0.12', fc='white', ec=INK, lw=1.3, zorder=2
        )
    )
    ax.text(
        x + pw / 2,
        0.9,
        t,
        fontsize=10,
        fontweight='bold',
        color=INK,
        ha='center',
        va='center',
        linespacing=1.2,
        zorder=4,
    )
    ax.text(x + pw / 2, 0.47, sub, fontsize=8.6, color=SOFT, ha='center', va='center', zorder=4)

for ext in ('png', 'svg', 'pdf'):
    fig.savefig(os.path.join(HERE, 'workflow.' + ext), dpi=170, facecolor='white')
print('wrote workflow.png, workflow.svg, workflow.pdf in', HERE)
