"""Command line of the package: python rfo.py <command> [options].

The commands are listed in the order of the workflow (docs/WORKFLOW.md).
"""

import argparse
import os
import sys

from . import __version__

STEPS = '''
commands, in the order of the workflow

  data protection
    originals          record (make) or check (verify) SHA-256 fingerprints of the original files

  quality control
    qc                 integrity, noise, event screening, P and S picks, processed waveforms

  overview of the data set
    inventory          event coverage, event selection, station availability, integrity
    snr                P and S signal-to-noise ratio against magnitude and distance
    response           sensor responses compared
    noise-hour         noise against local time of day

  sensor orientation
    orientation        horizontal orientation from P particle motion
    orientation-apply  decide the corrections and write a corrected StationXML

  manual picks (optional)
    handpick           interactive blind picking
    compare-picks      manual against automatic picks

  travel times
    station-terms      station and event terms, direction dependence, S against P

  receiver functions
    rf                 deconvolution, acceptance rules, station screen, H-k stacking
    section            station stacks sorted along the network
    ccp                common-conversion-point stack, profiles and depth slices
    compare-rf         two receiver-function runs compared (for example two Gaussian widths)

  crust
    crust-template     empty crust table
    crust-from-rf      crust table from the receiver functions and published depths
    crust              crustal correction of the station terms, mantle residual

  other
    mseed-metadata     event catalog, station table and file index for the miniSEED copy
    demo               synthetic test of the deconvolution and the H-k stack (no data needed)

python rfo.py <command> --help lists the options of a command.
'''


def build_parser():
    from . import (
        compare,
        crust,
        imaging,
        manual_picks,
        mseed_metadata,
        orientation,
        originals,
        overview,
        qc_pipeline,
        receiver_functions,
        station_terms,
        synthetic,
    )

    ap = argparse.ArgumentParser(
        prog='rfo.py',
        description='Teleseismic quality control, travel times and receiver '
        'functions for SOD event folders (version %s).' % __version__,
        epilog=STEPS,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    ap.add_argument('--version', action='version', version=__version__)
    sub = ap.add_subparsers(dest='cmd', metavar='<command>')
    sub.required = True
    workers = max(1, (os.cpu_count() or 2) - 1)

    # ---- commands with their own argument lists
    p = sub.add_parser(
        'originals', description=originals.__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    originals.add_arguments(p)
    p.set_defaults(fn=originals.run)
    p = sub.add_parser(
        'qc', description=qc_pipeline.__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    qc_pipeline.add_arguments(p)
    p.set_defaults(fn=qc_pipeline.run)
    p = sub.add_parser(
        'mseed-metadata',
        description=mseed_metadata.__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    mseed_metadata.add_arguments(p)
    p.set_defaults(fn=mseed_metadata.run)
    p = sub.add_parser(
        'demo', description=synthetic.__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    synthetic.add_arguments(p)
    p.set_defaults(fn=synthetic.run)

    # ---- analysis commands: all read the quality-control output folder
    def add(name, fn, data=False, inv=False):
        p = sub.add_parser(
            name,
            description=(fn.__doc__ or sys.modules[fn.__module__].__doc__),
            formatter_class=argparse.RawDescriptionHelpFormatter,
        )
        p.add_argument('--qc', required=True, help='output folder of the qc command, e.g. qc_out')
        p.add_argument('--subdir', default='', help='write to <qc>/analysis/<subdir> (robustness runs)')
        if data:
            p.add_argument(
                '--data', required=True, help='raw SAC event folders (read only), e.g. data/raw_sac'
            )
        if inv:
            p.add_argument('--inventory', required=True, help='StationXML of the network')
        p.set_defaults(fn=fn, needs_qc=True)
        return p

    p = add('inventory', overview.cmd_inventory)
    p.add_argument('--listed', type=int, default=None, help='number of events listed by the download')
    p.add_argument(
        '--hq-picks', type=int, default=20, help='class A to C P picks that make a high-quality event'
    )
    add('snr', overview.cmd_snr)
    add('response', overview.cmd_response, inv=True)
    p = add('noise-hour', overview.cmd_noise_hour)
    p.add_argument('--utc-offset', type=float, default=4.0, help='local time minus UTC in hours (Oman: 4)')

    p = add('orientation', orientation.cmd_orientation, data=True, inv=True)
    p.add_argument('--fmin', type=float, default=0.05, help='band, low corner (Hz)')
    p.add_argument('--fmax', type=float, default=0.2, help='band, high corner (Hz)')
    p.add_argument('--t1', type=float, default=-1.0, help='window start relative to the P pick (s)')
    p.add_argument('--t2', type=float, default=6.0, help='window end relative to the P pick (s)')
    p.add_argument('--min-snr', type=float, default=5.0, help='P signal-to-noise ratio at 0.05 to 1 Hz')
    p.add_argument('--min-cc', type=float, default=0.7, help='minimum Z-R correlation of a measurement')
    p.add_argument('--min-n', type=int, default=5, help='measurements needed before a station is flagged')
    p.add_argument(
        '--flag-deg', type=float, default=10.0, help='offset (deg) above which a station is flagged'
    )
    p.add_argument('--workers', type=int, default=workers)
    p = add('orientation-apply', orientation.cmd_orientation_apply, inv=True)
    p.add_argument('--out-xml', default=None, help='corrected StationXML (default: <inventory>_oriented.xml)')
    p.add_argument('--min-cc', type=float, default=0.7, help='minimum Z-R correlation of a measurement')
    p.add_argument('--min-n', type=int, default=10, help='measurements a station needs')
    p.add_argument('--min-deg', type=float, default=5.0, help='smallest offset (deg) that is corrected')
    p.add_argument('--bin-deg', type=float, default=30.0, help='width of the back-azimuth bins (deg)')
    p.add_argument(
        '--consistency',
        type=float,
        default=12.0,
        help='largest spread of the bin medians (deg) for a sensor offset',
    )
    p.add_argument('--bin-min', type=int, default=8, help='measurements a back-azimuth bin needs to count')
    p.add_argument(
        '--max-ci', type=float, default=10.0, help='largest half-width of the 95%% interval (deg) to correct'
    )
    p.add_argument(
        '--max-mad',
        type=float,
        default=15.0,
        help='largest scatter (MAD, deg) of the single-event values to correct',
    )

    p = add('handpick', manual_picks.cmd_handpick, data=True)
    p.add_argument(
        '--picker', required=True, help='name or initials of the analyst (goes into the file name)'
    )
    p.add_argument('--n-events', type=int, default=10)
    p.add_argument('--seed', type=int, default=1, help='random choice of events (same seed = same events)')
    p.add_argument(
        '--events', default=None, help='glob pattern instead of a random choice, e.g. "Event_2015_*"'
    )
    p = add('compare-picks', manual_picks.cmd_compare_picks)
    p.add_argument(
        '--manual', nargs='*', default=None, help='manual pick files (default: all in <qc>/analysis)'
    )

    p = add('station-terms', station_terms.cmd_station_terms)
    p.add_argument('--min-quality', default='C', choices=['A', 'B', 'C'], help='lowest pick class used')
    p.add_argument(
        '--sectors',
        default='3',
        help='number of back-azimuth sectors with equal numbers of events '
        '(default 3), or explicit edges such as 0,60,90,180,360',
    )
    p.add_argument('--bootstrap', type=int, default=200)
    p.add_argument('--min-picks', type=int, default=10, help='picks needed for a reliable station term')
    p.add_argument(
        '--timing-mad',
        type=float,
        default=0.5,
        help='residual MAD (s) above which a station is a timing suspect',
    )
    p.add_argument('--exclude', default='', help='stations left out of the solution, comma separated')

    p = add('rf', receiver_functions.cmd_rf)
    p.add_argument('--workers', type=int, default=workers)
    p.add_argument('--fmin', type=float, default=0.05, help='band, low corner (Hz)')
    p.add_argument('--fmax', type=float, default=2.0, help='band, high corner (Hz)')
    p.add_argument(
        '--gauss',
        type=float,
        default=2.0,
        help='Gaussian width parameter a (2.0 is about 1 Hz, 2.5 about 1.2 Hz)',
    )
    p.add_argument(
        '--out',
        default='rf',
        help='output folder below --qc (default rf; e.g. rf_g25 for another Gaussian width)',
    )
    p.add_argument('--min-fit', type=float, default=70.0, help='minimum deconvolution fit (%%)')
    p.add_argument('--min-quality', default='C', choices=['A', 'B', 'C'], help='lowest P pick class used')
    p.add_argument(
        '--exclude', default='', help='stations to skip, comma separated (for example timing suspects)'
    )
    p.add_argument('--min-rf', type=int, default=20, help='receiver functions needed for an H-k result')
    p.add_argument('--vp', type=float, default=6.3, help='average crustal Vp (km/s)')
    p.add_argument(
        '--vp-table', default='', help='CSV with the columns station,vp: station-specific crustal Vp (km/s)'
    )
    p.add_argument(
        '--station-cc',
        type=float,
        default=0.3,
        help='remove an accepted function whose correlation with the '
        'leave-one-out station stack (-1 to 25 s) is below this (0 disables)',
    )
    p.add_argument(
        '--station-amp',
        type=float,
        default=2.0,
        help='remove an accepted function whose RMS after the direct '
        'pulse exceeds this multiple of the station median (0 disables)',
    )
    p.add_argument('--h-min', type=float, default=20.0, help='H-k search, smallest depth (km)')
    p.add_argument('--h-max', type=float, default=60.0, help='H-k search, largest depth (km)')
    p.add_argument('--k-min', type=float, default=1.60, help='H-k search, smallest Vp/Vs')
    p.add_argument('--k-max', type=float, default=2.00, help='H-k search, largest Vp/Vs')
    p.add_argument('--bootstrap', type=int, default=200)
    p.add_argument('--k-fixed', type=float, default=1.75, help='Vp/Vs for the Ps-only depth')
    p.add_argument('--ps-hmin', type=float, default=24.0, help='Ps-only depth search, smallest depth (km)')
    p.add_argument('--ps-hmax', type=float, default=55.0, help='Ps-only depth search, largest depth (km)')
    p.add_argument(
        '--recompute', action='store_true', help='redo the deconvolution even if rf_all.pkl exists'
    )

    p = add('section', imaging.cmd_section)
    imaging.add_section_arguments(p)
    p = add('ccp', imaging.cmd_ccp)
    imaging.add_ccp_arguments(p)
    p = add('compare-rf', compare.cmd_compare_rf)
    compare.add_arguments(p)

    p = add('crust-template', crust.cmd_crust_template)
    p.add_argument('--force', action='store_true', help='overwrite an existing template')
    p = add('crust-from-rf', crust.cmd_crust_from_rf)
    p.add_argument('--rf', default='rf', help='receiver-function folder below --qc (default rf)')
    p.add_argument(
        '--published', default=None, help='CSV with station, moho_km, moho_err_km, source (published values)'
    )
    p.add_argument(
        '--prefer',
        default='auto',
        choices=['auto', 'rf', 'published'],
        help='auto: H-k where quality A or B, else the published value, else the Ps-only depth; '
        'rf: never the published value; published: the published value wherever it exists',
    )
    p.add_argument(
        '--min-err', type=float, default=2.0, help='smallest Moho uncertainty (km) written to the table'
    )
    p.add_argument('--ps-err', type=float, default=5.0, help='smallest uncertainty (km) of a Ps-only depth')
    p.add_argument(
        '--orla',
        default=None,
        help='3-D Vs model in netCDF (ORLA2022 of Weidle et al. 2022) for the Vs columns',
    )
    p.add_argument('--orla-kappa', type=float, default=1.75, help='Vp/Vs that turns the model Vs into Vp')
    p.add_argument('--orla-zmax', type=float, default=33.0, help='bottom of the model column (km)')
    p = add('crust', crust.cmd_crust)
    p.add_argument('--model', required=True, help='filled crust table (from crust-template or crust-from-rf)')
    p.add_argument(
        '--terms', default=None, help='station terms to correct (default <qc>/analysis/station_terms_P.csv)'
    )
    p.add_argument('--vp-sediment', type=float, default=4.0)
    p.add_argument('--vp-ophiolite', type=float, default=7.0)
    p.add_argument('--vp-crust', type=float, default=6.3)
    p.add_argument('--vp-mantle', type=float, default=8.04)
    p.add_argument(
        '--mantle-km', type=float, default=150.0, help='depth range below the Moho for the %% estimate'
    )
    p.add_argument(
        '--moho-err', type=float, default=3.0, help='Moho uncertainty (km) where the table gives none'
    )
    p.add_argument(
        '--orla',
        default=None,
        help='3-D Vs model in netCDF: also correct with it as the crust above --orla-zmax',
    )
    p.add_argument('--orla-kappa', type=float, default=1.75)
    p.add_argument('--orla-zmax', type=float, default=33.0)
    return ap


def main(argv=None):
    a = build_parser().parse_args(argv)
    if a.cmd != 'handpick':
        import matplotlib

        matplotlib.use('Agg')
    if getattr(a, 'needs_qc', False) and not os.path.isdir(os.path.join(a.qc, 'tables')):
        sys.exit('ERROR: %s is not an output folder of the qc command (no tables/)' % a.qc)
    a.fn(a)


if __name__ == '__main__':
    main()
