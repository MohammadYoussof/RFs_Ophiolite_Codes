# Quality control and picking: methods

Command: `python rfo.py qc` (module `rf_ophiolite/qc_pipeline.py`).

Input is the SAC folder written by SOD. The SAC headers carry the event origin, location,
depth and magnitude and the station distance and back-azimuth. The miniSEED folder holds
the same samples and is not read by this step.

```
python rfo.py qc --data data/raw_sac --inventory stations.xml --out qc_out --workers 6
```

| Option | Effect |
|---|---|
| `--events "Event_2015_*"` | process a subset of the event folders |
| `--workers N` | number of events processed in parallel (default: CPU cores minus 1) |
| `--overwrite` | process again events that are already done |
| `--no-figures`, `--no-waveforms` | skip the per-event figures or the processed SAC files |
| `--aggregate-only` | rebuild the tables and summary plots from finished events |

An interrupted run continues by itself: completed events are skipped, while events that
ended with an error, half-written results and results of an older version are processed
again. Expect about 50 s per event per worker with the full instrument response, and
about 10 s with sensitivities only. Every threshold below is in the `CONFIG` block at the
top of the module.

## 1. Protection of the original data

The originals are protected in four ways:

1. the step only reads them;
2. it refuses an output folder inside the data folder or its sibling `raw_mseed` (and
   the other way round), comparing the folders themselves so that links or upper and
   lower case cannot hide an overlap;
3. it compares the size and the time stamp of every input file before and after the run;
4. `rfo.py originals` compares SHA-256 checksums of every file before and after, which
   proves byte-for-byte identity. The manifest can also be checked without Python:
   `shasum -a 256 -c --quiet originals_sha256.txt` prints nothing when all files match.

`chmod -R a-w` on the data folders adds a fifth barrier. It is not enforced for an
administrator using `sudo`, nor on exFAT or FAT external disks.

## 2. Processing steps

### 2.1 Integrity checks on raw counts, per channel

- Dead: constant value.
- Clipped: absolute counts at 95 % or more of the 24-bit full scale.
- Flat segments: identical samples for more than 1 s. Zero-filled gaps are excluded from
  this test.
- Gaps: a channel delivered in two or more files (`*.1.sac`, `*.2.sac`) is merged and the
  gap is recorded.
- Isolated one-sample spikes and DC offset.

Dead or clipped channels are not used.

### 2.2 Instrument correction

- If the StationXML has response stages, the full response is removed to velocity
  (ObsPy `remove_response`, pre-filter 0.005 to 0.01 Hz and 40 to 45 Hz, water level
  60 dB).
- If it has overall sensitivities only, the counts are divided by the sensitivity, which
  gives flat-band velocity.
- Without an inventory the data stay in counts and no noise spectra are computed.

Use the StationXML at response level. In a test on the COOL data, sensitivities alone
left a sensor-dependent bias of about 10 ms in the relative P times, shifted absolute S
times by up to 0.13 s, and overestimated the noise of one sensor type above 5 Hz by
about 3 dB.

### 2.3 Station noise

For every channel and event:

- a 100 s window before P, ending 8 s before the predicted arrival;
- Welch spectrum (50 s Hann segments, 75 % overlap), converted to acceleration and
  averaged over one-octave bands centered every 1/8 octave between 0.05 and 30 Hz
  (McNamara and Buland 2004);
- after all events: a probability density function of the spectra per channel (mode,
  10th and 90th percentiles) against the low and high noise models (Peterson 1993).

The lower frequency limit is set by the roughly 2 minutes of data before the event in
each event file.

### 2.4 Phase visibility

Signal-to-noise ratio (RMS of the signal window over RMS of the noise window), measured
against the ak135 predictions:

| Phase | Component and band | Signal window | Noise window |
|---|---|---|---|
| P | vertical, 0.5 to 2 Hz (also reported at 0.05 to 1 Hz) | P -1 to +9 s | P -65 to -5 s |
| S | transverse, 0.03 to 0.2 Hz | S -2 to +18 s | S -65 to -5 s |

Visibility classes: clear (5 or more), visible (3 to 5), weak (2 to 3), none (below 2).

### 2.5 Event screening

An event is rejected when

- fewer than 3 vertical channels are alive (no usable data);
- fewer than 3 stations have a P signal-to-noise ratio of 3 or more (no main phase);
- fewer than 4 P picks reach class A to C (P not coherent across the array).

After all events are processed, an event is also flagged when its noise before the event
is more than 3 times the median of the stations, which usually means that another
earthquake overlaps the window. S is picked only if at least 3 stations show S with a
ratio of 3 or more, and the S picks of an event are kept only if at least 3 reach class
A to C.

### 2.6 Travel-time picks

The stations are aligned together per event and phase, using traces with a
signal-to-noise ratio of 2 or more:

1. Coarse alignment by iterative cross-correlation and stacking (ICCS; Lou et al.
   2013) of an onset characteristic function, the smoothed envelope of the causally
   filtered trace. It has no oscillations, so this step cannot skip a cycle. The stack is
   re-centered in every iteration so that it cannot drift.
2. Fine alignment by ICCS on the waveforms in the picking band. The lag limit is
   below half a dominant period: 0.35 s for P, 1.5 s for S.
3. Multi-channel cross-correlation (MCCC; VanDecar and Crosson 1990). All station
   pairs are correlated and the relative delays are solved by weighted least squares
   with outlier rejection. This gives a relative time and a standard error per trace.
4. Absolute onset. The aligned traces are filtered causally (P: 0.3 to 8 Hz, S: 0.02
   to 0.5 Hz) and stacked. The first arrival is found where the envelope of the
   phase-weighted stack (Schimmel and Paulssen 1997) first reaches half its maximum; the
   onset is the AIC minimum (Maeda 1985) of the linear stack just before that point. The
   same onset offset applies to every station, so absolute and relative times stay
   consistent (compare Rawlinson and Kennett 2004).
5. Quality class of each pick:

| Phase | Class | Minimum correlation with the stack (leave-one-out) | Maximum standard error | Minimum signal-to-noise ratio |
|---|---|---|---|---|
| P | A | 0.85 | 0.05 s | 5 |
| P | B | 0.75 | 0.10 s | 3 |
| P | C | 0.65 | 0.20 s | 2 |
| S | A | 0.90 | 0.25 s | 6 |
| S | B | 0.85 | 0.50 s | 4 |
| S | C | 0.75 | 1.00 s | 3 |

   Every class also requires the onset function of the trace to agree with the
   leave-one-out stack within 0.3 s (P) or 1.5 s (S), which guards against cycle skips.
   Class D picks are written to `picks_all.csv` but never used.
6. Uncertainty of the absolute time. The onset of step 4 is common to all stations of
   the event, so its error is shared. It is estimated by repeating step 4 on 200
   bootstrap resamples of the stations (fixed random seed): `sigma_onset` is the
   standard deviation of the resampled onsets, at least one sample. The absolute
   uncertainty of each pick is `sqrt(sigma_mccc^2 + sigma_onset^2)`. `sigma_onset` does
   not affect the relative times.

Three design choices matter for the accuracy of the travel times:

- Causal filters for onsets. Zero-phase filters create an acausal precursor that
  moved onsets about 1 s early in testing.
- Windows end before interfering phases. PcP arrives only 4 to 10 s after P at 76 to
  82 degrees. Every correlation window ends 1 s before PcP (for P) or 3 s before ScS, SKS
  and SKKS (for S). Without this, some stations skipped one cycle (about 1.4 s).
- S only up to 80 degrees. Beyond that, S overlaps SKS, SKKS and ScS.

### 2.7 Residuals

- Absolute residual: observed minus ak135 (Kennett et al. 1995).
- Elevation correction: `h * sqrt(1/v^2 - p^2)`, with v = 5.8 km/s (P) or 3.36 km/s (S)
  and p the ak135 ray parameter.
- Relative residual: the event mean of the class A to C picks of the same phase removed
  (for teleseismic tomography), given with and without the elevation correction.

### 2.8 Processed waveforms

For accepted events and stations with usable P:

- components Z, R and T in velocity, with a common start and end;
- zero-phase band-pass 0.02 to 4 Hz, resampled to 20 Hz;
- normalized by the peak Z amplitude from P -5 to +60 s. The same factor is used for all
  three components, so relative amplitudes are preserved;
- time zero is the earthquake origin, so picks are travel times.

## 3. SAC headers of the processed files

| Header | Content |
|---|---|
| reference time (NZ*), O | earthquake origin; O = 0, IZTYPE = IO |
| B, E | time of the first and the last sample after the origin (s) |
| T0, KT0 | P pick (travel time, s) and class, for example `P_A` |
| T1, KT1 | S pick and class, if picked |
| T8, T9 | ak135 predictions of P and S |
| USER0 | normalization factor (m/s per unit) |
| USER1, USER2 | P: correlation with the stack, MCCC standard error (s) |
| USER3, USER4 | S: correlation with the stack, MCCC standard error (s) |
| EVDP | km (the raw SOD files use meters) |
| CMPAZ | R = back-azimuth + 180 degrees, T = back-azimuth - 90 degrees |

## 4. The two pick files

Both hold accepted events and classes A to C only, one row per station and phase, the
same rows in the same order. P is picked on the vertical component (0.5 to 2 Hz), S on
the transverse component (0.03 to 0.2 Hz, up to 80 degrees).

Common columns: `event`, `origin_time` (catalog origin, UTC), `ev_lat`, `ev_lon`,
`ev_depth_km`, `mag`, `network`, `station`, `sta_lat`, `sta_lon`, `sta_elev_m`, `phase`,
`quality` (A, B, C), `dist_deg`, `baz_deg`, `rayp_s_per_deg` (ak135), `model`,
`component`, `band`.

`picks_absolute.csv` (relocation, absolute travel-time tomography)

| Column | Meaning |
|---|---|
| `pick_time` | absolute arrival time (UTC) |
| `tt_obs_s` | observed travel time = pick minus origin |
| `tt_pred_s` | ak135 travel time |
| `res_s` | tt_obs minus tt_pred |
| `elev_corr_s` | delay of the station elevation |
| `res_elev_s` | res_s minus elev_corr_s |
| `sigma_abs_s` | absolute uncertainty |
| `sigma_onset_s`, `sigma_mccc_s` | the two parts of `sigma_abs_s` |
| `cc`, `snr` | correlation with the leave-one-out stack, signal-to-noise ratio |

`picks_relative.csv` (relative travel-time tomography)

| Column | Meaning |
|---|---|
| `dt_rel_s` | relative delay: res_elev_s minus its mean over the same event and phase |
| `dt_rel_noelev_s` | the same without the elevation correction (for models that include topography) |
| `sigma_mccc_s` | MCCC uncertainty of the relative time |
| `cc`, `snr` | as above |
| `n_picks_event` | picks of this event and phase (the mean is taken over them) |
| `tt_pred_s`, `elev_corr_s` | for reference |

Absolute residuals include the errors of the catalog origin time and hypocenter; relative
delays cancel the part common to the array. If picks are dropped from
`picks_relative.csv` (for example to keep class A only), subtract the event mean again
over the remaining picks of each event and phase.

## 5. Limitations

- Sensor responses are nominal. Calibrations of individual sensors are not part of
  the metadata.
- Noise spectra cover 0.05 to 30 Hz only, and their statistics improve with the
  number of events per station.
- Orientation is taken from the metadata. Check it with step 5 of the workflow, which
  needs many events per station.
- Correlation assumes similar waveforms across the array, which holds over an
  aperture of a few degrees. For much larger arrays, pick in sub-arrays.

## 6. References

- Kennett, B. L. N., Engdahl, E. R. and Buland, R. (1995). Constraints on seismic velocities in the Earth from traveltimes. Geophysical Journal International 122, 108-124.
- Lou, X., van der Lee, S. and Lloyd, S. (2013). AIMBAT: a Python/Matplotlib tool for measuring teleseismic arrival times. Seismological Research Letters 84, 85-93.
- Maeda, N. (1985). A method for reading and checking phase times in autoprocessing system of seismic wave data. Zisin 38, 365-379.
- McNamara, D. E. and Buland, R. P. (2004). Ambient noise levels in the continental United States. Bulletin of the Seismological Society of America 94, 1517-1527.
- Peterson, J. (1993). Observations and modeling of seismic background noise. U.S. Geological Survey Open-File Report 93-322.
- Rawlinson, N. and Kennett, B. L. N. (2004). Rapid estimation of relative and absolute delay times across a network by adaptive stacking. Geophysical Journal International 157, 332-340.
- Schimmel, M. and Paulssen, H. (1997). Noise reduction and detection of weak, coherent signals through phase-weighted stacks. Geophysical Journal International 130, 497-505.
- VanDecar, J. C. and Crosson, R. S. (1990). Determination of teleseismic relative phase arrival times using multi-channel cross-correlation and least squares. Bulletin of the Seismological Society of America 80, 150-169.
