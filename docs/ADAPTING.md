# Using the codes on another network

The codes read event folders written by SOD and a StationXML file. Nothing else is
specific to the COOL network; the few settings that are sit in the places below.

## 1. Download

Copy `sod/cool_5H_full.xml` and change the network, the time range, the magnitude and
the distance range. The table in [../sod/README.md](../sod/README.md) says where each
one is. Keep the two writers (`sacWriter`, `mseedWriter`): the quality control reads the
SAC headers that SOD fills in (origin, location, depth, magnitude, distance,
back-azimuth).

## 2. Settings in the code

| Setting | Where | COOL value | Meaning |
|---|---|---|---|
| `network_rename` | `CONFIG` in `rf_ophiolite/qc_pipeline.py` | `{'5H2013': '5H'}` | SOD names a temporary network with its start year; this maps the name in the files to the network code of the StationXML. Use `{}` for a permanent network. |
| `CHANNEL_BAND` | `rf_ophiolite/common.py` | `'HH'` | first two letters of the channel codes in the raw files (`'BH'` for 20 or 40 Hz broadband channels) |
| `full_scale_counts` | `CONFIG` | `2**23 - 1` | full scale of the digitizer, for the clipping test |
| `evdp_in_meters` | `CONFIG` | `True` | SOD writes the event depth in meters into the SAC header |
| `pick` | `CONFIG` | see the file | bands, windows, lag limits and class thresholds of the P and S picks |
| `out_band`, `out_rate` | `CONFIG` | 0.02 to 4 Hz, 20 Hz | band and sampling rate of the processed waveforms |

The noise spectra need about 2 minutes of data before P in every event file (the window
of the SOD recipe starts 120 s before P). The traces are resampled to the rate given in
`pick` before picking, so other sampling rates work; the causal onset band of P reaches
8 Hz, so lower its upper corner for data sampled below 20 Hz.

## 3. Options to review for a new region

| Command | Option | COOL value | Why |
|---|---|---|---|
| `noise-hour` | `--utc-offset` | 4 | local time minus UTC |
| `station-terms` | `--timing-mad` | 0.5 s | scatter above which a station is treated as a clock problem |
| `rf` | `--vp` or `--vp-table` | 6.3 km/s | average crustal Vp of the H-k stack |
| `rf` | `--h-min`, `--h-max` | 20, 60 km | depth range of the H-k search |
| `rf` | `--ps-hmin`, `--ps-hmax` | 24, 55 km | depth range in which the Moho conversion is read; it should bracket the expected Moho and exclude shallow conversions |
| `rf` | `--gauss` | 2.0 | Gaussian width; compare two values with `compare-rf` |
| `ccp` | `--moho-model`, `--vp-crust`, `--kappa-crust` | 40 km, 6.3, 1.75 | 1-D model of the migration |
| `ccp` | `--region`, `--profile`, `--half-width` | network extent, two default profiles, 45 km | grid and profiles |
| `crust-from-rf` | `--published` | `data/moho_published.csv` | published depths of the region, or leave it out |
| `crust` | `--vp-crust`, `--vp-sediment`, `--vp-ophiolite` | 6.3, 4.0, 7.0 km/s | velocities of the crustal column |

## 4. Array size

The picks rest on the similarity of the waveforms across the array, which holds over an
aperture of a few degrees. For a much larger array, download and process sub-arrays
separately, with one SOD recipe and one output folder per sub-array.

## 5. A first run on new data

1. Download one event at a few stations and check it with `sod/verify_download.py`.
2. Run `rfo.py qc` on that event with `--events` and look at the figure in
   `qc_out/figures`.
3. Run the synthetic test (`rfo.py demo`) to see that the installation is complete.
4. Then run the steps in order, as in [WORKFLOW.md](WORKFLOW.md).
