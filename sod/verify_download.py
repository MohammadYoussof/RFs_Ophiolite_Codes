"""Compare a SOD download with the same data fetched directly from the data center.

Run it in the folder that holds raw_mseed/ and raw_sac/ (for example after the test
recipe). For every miniSEED file the same time window is requested from the FDSN
web service and compared sample by sample with the miniSEED and the SAC copy.
Every line should end with three times True.

    python verify_download.py                       GEOFON, network 5H
    python verify_download.py --base-url https://service.iris.edu --network IU

Needs ObsPy and an internet connection. The local files are only read.
"""

import argparse
import glob
import warnings

import numpy as np
from obspy import read
from obspy.clients.fdsn import Client

warnings.filterwarnings('ignore')


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--base-url', default='https://geofon.gfz.de', help='FDSN web service of the data center')
    ap.add_argument('--network', default='5H', help='network code at the data center')
    ap.add_argument('--max-files', type=int, default=30, help='largest number of files to check')
    a = ap.parse_args()
    files = sorted(glob.glob('raw_mseed/Event_*/*.mseed'))[: a.max_files]
    if not files:
        raise SystemExit('No miniSEED files found: run this in the folder that holds raw_mseed and raw_sac.')
    client = Client(a.base_url)
    all_ok = True
    for f in files:
        m = read(f)[0]
        s = read(f.replace('raw_mseed', 'raw_sac').replace('.mseed', '.sac'))[0]
        g = client.get_waveforms(
            a.network, m.stats.station, m.stats.location, m.stats.channel, m.stats.starttime, m.stats.endtime
        ).merge()[0]
        g.trim(m.stats.starttime, m.stats.endtime)
        n = min(len(m.data), len(s.data), len(g.data))
        same_start = g.stats.starttime == m.stats.starttime
        mseed_ok = bool(np.array_equal(m.data[:n], g.data[:n]))
        sac_ok = bool(np.allclose(s.data[:n], g.data[:n]))
        all_ok = all_ok and same_start and mseed_ok and sac_ok
        print(
            '%-22s same start: %s | miniSEED equals data center: %s | SAC equals data center: %s'
            % (m.id, same_start, mseed_ok, sac_ok)
        )
    print('\nRESULT:', 'all files match the data center' if all_ok else 'MISMATCH in at least one file')


if __name__ == '__main__':
    main()
