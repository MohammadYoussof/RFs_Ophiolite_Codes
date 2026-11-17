# Output files

Everything is written below the output folder of the quality control (`qc_out` in the
examples). Nothing is written into the data folders.

```
qc_out/
    tables/        step 3   event and station tables
    picks/         step 3   the two pick files
    processed/     step 3   processed waveforms, one folder per event
    noise/         step 3   noise spectra per station
    figures/       step 3   one figure per event
    work/          step 3   intermediate results per event (used to continue a run)
    analysis/      steps 4 to 6, 8 to 10: tables and figures
    rf/            step 7   receiver functions (the folder name is the option --out)
```

## Step 3: `qc`

| File | Content |
|---|---|
| `tables/event_summary.csv` | per event: status, reason of rejection, counts, window limits, mean residual, noise flag |
| `tables/station_event_qc.csv` | per station and event: integrity, signal-to-noise ratios, visibility, usability |
| `tables/picks_all.csv` | every pick, classes A to D, rejected events included (diagnostic) |
| `tables/station_performance.csv` | per station: events, detections, pick yield, P station term, noise level, dead and flat counts |
| `picks/picks_absolute.csv` | absolute travel times and residuals (relocation, absolute tomography) |
| `picks/picks_relative.csv` | relative delays, event mean removed (relative travel-time tomography) |
| `processed/<event>/*.sac` | Z, R, T in velocity, band-passed, 20 Hz, time zero at the origin, picks in the headers |
| `noise/PSD_PDF_<station>.png` | probability density of the noise spectra, three components |
| `noise/network_noise_overview.png` | mode spectrum of all stations and a station ranking |
| `noise/channel_noise_summary.csv` | band-averaged noise levels per channel |
| `figures/qc_<event>.png` | full record, P and S around the predictions, P and S aligned on the picks |
| `picks_summary.png` | station terms, uncertainties, residuals against back-azimuth, classes |

The columns of the two pick files and the SAC headers of the processed files are
described in [METHODS_QC.md](METHODS_QC.md).

## Step 4: overview (`analysis/`)

| Command | Files |
|---|---|
| `inventory` | `event_funnel.csv`, `event_funnel.png`, `event_coverage.png`, `station_availability.png`, `integrity_by_station.csv`, `integrity_by_station.png` |
| `snr` | `phase_visibility.png`, `visibility_by_magnitude.csv`, `visibility_by_distance.csv` |
| `response` | `sensor_responses.png`, `sensor_response_comparison.csv` |
| `noise-hour` | `noise_vs_local_hour.png`, `noise_by_window.csv`, `noise_day_night.csv` |

## Step 5: sensor orientation (`analysis/`)

| Command | Files |
|---|---|
| `orientation` | `orientation_measurements.csv` (one line per station and event), `orientation_by_station.csv`, `orientation.png` |
| `orientation-apply` | `orientation_corrections.csv` (decision per station), the corrected StationXML (`--out-xml`) |

## Manual picks (`analysis/`)

| Command | Files |
|---|---|
| `handpick` | `manual_picks_<picker>.csv` |
| `compare-picks` | `pick_comparison.csv`, `pick_comparison_summary.csv`, `pick_comparison.png` |

## Step 6: `station-terms` (`analysis/`)

| File | Content |
|---|---|
| `station_terms_P.csv`, `station_terms_S.csv` | term, bootstrap error, picks used, flags `connected`, `timing_suspect`, `reliable` |
| `event_terms_P.csv`, `event_terms_S.csv` | event terms |
| `pick_residuals_P.csv`, `pick_residuals_S.csv` | residual of every pick after the solution |
| `station_terms_by_sector_P.csv`, `..._S.csv` | terms per back-azimuth sector |
| `station_terms_P.png`, `station_terms_S.png` | maps of the terms |
| `station_terms_by_sector_P.png`, `..._S.png` | maps per sector |
| `timing_check_P.png`, `timing_check_S.png` | residual against date for the most scattered stations |
| `s_vs_p_terms.png` | S terms against P terms |

## Step 7: `rf` (`rf/`)

| File | Content |
|---|---|
| `rf_summary.csv` | every function: fit, direct pulse, largest amplitude before and after, accepted or not, result of the station screen |
| `hk_results.csv` | per station: number of functions, H and Vp/Vs with bootstrap errors, sensitivity to Vp, sharpness of the maximum, Ps time, Ps-only depth, quality class |
| `ps_candidates.csv` | local maxima of the moveout-corrected stack (depth, Ps time, amplitude) |
| `RF_<station>.png` | stack of the accepted functions, the functions sorted by back-azimuth with the predicted times of the conversion and its multiples, and the H-k grid |
| `rf_<station>.npz` | accepted functions of the station: `rf`, `t0`, `dt`, `rayp`, `baz`, `event`, and the H-k grid `S`, `H`, `K` |
| `hk_map.png` | maps of the crustal thickness and Vp/Vs from the H-k stacks |
| `rf_all.pkl` | all deconvolved functions; a second run with other H-k settings reads it instead of deconvolving again (`--recompute` ignores it) |

## Step 8: images (`analysis/`)

| Command | Files |
|---|---|
| `section` | `rf_section.png`, `rf_section_stacks.csv` |
| `ccp` | `ccp_volume.npz` (`stack`, `raw`, `count`, `lon`, `lat`, `z`), `ccp_profiles.png`, `ccp_slices.png` |
| `compare-rf` | `rf_runs_counts.csv`, `rf_runs_by_station.csv`, `rf_runs_compare.png` |

## Steps 9 and 10: crust (`analysis/`)

| Command | Files |
|---|---|
| `crust-template` | `crust_model_template.csv` |
| `crust-from-rf` | `crust_model.csv`, `moho_vs_published.png` |
| `crust` | `crust_correction.csv`, `crust_correction.png`, `residual_profiles.png`, and `crust_correction_orla.png` with `--orla` |

Main columns of `crust_model.csv`: `station`, `moho_km`, `moho_err_km`, `method` (H-k,
published, Ps-only), the H-k and Ps-only values behind it, `ophiolite_km`,
`sediment_km`, and the published value for comparison.

Main columns of `crust_correction.csv`: the station term (`obs_term_s`), the predicted
crustal term (`crust_term_s`), the mantle residual (`residual_s`) with its error
(`residual_err_s`), and the residual as a velocity anomaly (`mantle_dVp_pct`).

## Other

| Command | Files |
|---|---|
| `originals make` | the manifest given with `-m` (SHA-256 and relative path per file) |
| `mseed-metadata` | `events_catalog.csv`, `station_channels.csv`, `mseed_index.csv` in `--out` |
| `demo` | `synthetic_test.png` in `--out` |
