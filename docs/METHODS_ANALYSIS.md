# Analysis: methods

Every command below reads the output of the quality control (`qc_out/tables`, `picks`,
`noise`, `processed`) and writes to `qc_out/analysis/`, or to
`qc_out/analysis/<subdir>` with `--subdir` for robustness runs. The receiver functions
are written to their own folder, `qc_out/rf` by default. The original SAC and miniSEED
files are only read, by `orientation` and `handpick`; nothing writes into the data
folders.

## 1. Overview of the data set (`inventory`, `snr`, `response`, `noise-hour`)

- Event selection. Events listed by the download, event folders, rejected (no P, not
  coherent, no data), accepted (4 or more class A to C P picks), accepted with 10 or
  more class A picks, accepted with S picks.
- Coverage. Back-azimuth and distance of all events; accepted events with the symbol
  size by the number of P picks; events with `--hq-picks` (20) or more picks marked as
  high quality.
- Availability and integrity. Stations that recorded each event and had usable P;
  per station the dead or constant traces, gaps and clipped samples.
- Phase visibility. P and S signal-to-noise ratio against magnitude and distance;
  detection rate per 0.25 magnitude bin and per 10 degree distance bin.
- Noise against the time of day. The noise windows of the quality control (vertical
  component, 1 to 5 Hz) against the local hour (`--utc-offset`), the network median per
  hour and the difference between day and night per station.
- Sensor responses. Amplitude response of every channel from the StationXML,
  compared with the network median at 0.01 to 50 Hz; a station whose sensitivity or
  corner differs is listed.

## 2. Sensor orientation (`orientation`, `orientation-apply`)

- For every class A to C P pick: Z, N, E at 0.05 to 0.2 Hz, window -1 to +6 s around the
  pick, P signal-to-noise ratio of 5 or more. The direction of the horizontal particle
  motion is the eigenvector of the horizontal covariance; its sign is fixed by the
  correlation of the radial component with the vertical. A measurement counts when the
  Z-R correlation is at least 0.7.
- Per station: the median of (observed minus expected back-azimuth), its 95 % interval
  (bootstrap), the scatter (MAD) and the medians in 30 degree back-azimuth bins.
- `orientation-apply` sorts the stations:

| Decision | Rule |
|---|---|
| insufficient | fewer than 10 measurements |
| nominal | offset below 5 degrees |
| uncertain | half-width of the 95 % interval above 10 degrees, or MAD above 15 degrees |
| direction-dependent | bin medians spread by more than 12 degrees among bins with at least 8 measurements; the deviation belongs to the structure, not to the sensor |
| correct | the rest |

  Sensors to correct get a new azimuth in a copy of the StationXML (north component =
  minus the offset, east component = north + 90 degrees). The quality control run again
  with this file changes only the horizontal components: S picks, radial and transverse
  traces, receiver functions.

## 3. Manual picks (`handpick`, `compare-picks`)

`handpick` shows the vertical traces of a random set of accepted events, one station at
a time, without the automatic pick, and stores the onset marked by the analyst.
`compare-picks` gives the difference between manual and automatic picks by quality
class: mean, median and standard deviation of the absolute times, and, with the mean
difference of each event removed (what a relative travel-time study sees), the standard
deviation and the fraction within 0.1 s. It is an independent check of the automatic
times and of the class thresholds.

## 4. Station terms and direction dependence (`station-terms`)

- Input. The absolute picks (`picks_absolute.csv`), classes down to `--min-quality`
  (default C), with the elevation-corrected residual `res_elev_s` (observed minus ak135,
  elevation reduced to sea level).
- Model. residual(i, j) = station term(i) + event term(j), solved by least squares
  with the condition that the station terms sum to zero. The event term absorbs the
  errors of origin time and hypocenter and the common deep path, so the station terms
  are relative delays within the network.
- Timing suspects. A station whose residual scatter (MAD over all its picks, around
  the network solution) exceeds `--timing-mad` (0.5 s) is left out of the solution and
  flagged. A clock fault looks like this and cannot be repaired by a constant term.
- Outliers. Picks with an absolute residual above max(0.6 s, 4 MAD) are dropped and
  the solution is repeated (3 times).
- Errors. Bootstrap over events (200 draws). A term is `reliable` with `--min-picks`
  (10) or more picks, when the station shares events with the main group of the network
  and is no timing suspect.
- Direction dependence. The events are split into back-azimuth sectors with equal
  numbers of events (`--sectors 3`, or explicit edges); per station and sector the mean
  residual after the solution is added to the term.
- S against P. The S terms are solved the same way and compared with the P terms
  (slope of S against P).

## 5. Receiver functions and crustal thickness (`rf`)

- Input. Processed Z and R traces of every accepted P pick (class A to C); stations
  in `--exclude` are skipped. Band 0.05 to 2 Hz, window -10 to +90 s around the P pick,
  2 s cosine taper.
- Deconvolution. Iterative time-domain deconvolution (Ligorria and Ammon 1999),
  radial by vertical, with the Gaussian width `--gauss` (2.0 is about 1 Hz), 200
  iterations. The function is normalized to a unit spike, so the amplitude of the direct
  P is the R/Z ratio.
- Acceptance rules for a single function (`rf_summary.csv` lists every function with
  these values):

| Quantity | Rule |
|---|---|
| fit of the deconvolution | at least 70 % |
| amplitude of the direct P | between 0.08 and 1.2 |
| before -1.5 s | nothing above half the direct pulse |
| between 1.5 and 30 s | nothing above 1.5 times the direct pulse |

- Station screen. The rules above look at one function at a time. The screen then
  compares every accepted function with the others of its station: its correlation with
  the leave-one-out station stack (-1 to 25 s) and its RMS after the direct pulse (1.5 to
  30 s) relative to the station median. A function is removed when the correlation is
  below `--station-cc` (0.3) or the RMS ratio is above `--station-amp` (2.0). Stations
  with fewer than 5 functions are not screened; 0 switches a test off.
- H-k stacking (Zhu and Kanamori 2000). Weights 0.6 (Ps), 0.3 (PpPs), -0.1
  (PpSs+PsPs); Vp 6.3 km/s (`--vp`, or one value per station with `--vp-table`); H 20 to
  60 km in 0.25 km; Vp/Vs 1.60 to 2.00 in 0.005; the ray parameter of each function; at
  least 20 functions per station. Errors from a bootstrap over functions (200 draws); the
  sensitivity to Vp from Vp plus and minus 0.2 km/s. `grid_pct_within_90` is the part of
  the grid within 90 % of the maximum (how sharp the maximum is), and `kappa_at_edge`
  flags a maximum at the edge of the Vp/Vs range: the multiples are too weak to constrain
  Vp/Vs, which is common under ophiolite and sediments.
- Ps-only depth. The Ps arrival time is read from the moveout-corrected stack within
  a depth range (`--ps-hmin 24`, `--ps-hmax 55` km; later phases are usually multiples of
  shallow interfaces) and converted with a fixed Vp/Vs (`--k-fixed 1.75`); its error from
  Vp/Vs plus and minus 0.05. `ps_candidates.csv` lists every local maximum of the stack
  between 12 and 70 km (depth, Ps time, amplitude), so that another phase can be chosen
  by hand where two conversions exist, for example the base of the ophiolite and the
  Moho.
- Quality. A: bootstrap error of H at most 1.5 km and a sharp maximum. B: at most
  3 km. C: otherwise, or Vp/Vs within 0.05 of the search limits.

The conversion at the Moho is written Ps in this text and in the table columns
(`t_ps_s`, `H_ps_km`). The figures label it PdS: the P-to-S conversion at a discontinuity
at depth d.

The Gaussian width sets the frequency content. A wider pulse (2.0) is more stable where
shallow layers ring; a narrower one (2.5) separates conversions that are close in time.
`compare-rf` shows what the choice changes (section 7).

## 6. Images (`section`, `ccp`)

- Section. For every station with `--min-rf` (20) or more functions the time axis of
  each function is stretched so that the Ps phase of any depth arrives at the time it
  would have for the reference ray parameter (`--pref 0.06` s/km) in the model `--vp`,
  `--kappa`. The functions are stacked and drawn side by side, sorted by longitude or
  latitude, with the Ps time of the Moho from the crust table marked.
- Common-conversion-point stack. Every function is mapped from time to depth along
  its own ray in a 1-D model (crust 6.3 km/s and Vp/Vs 1.75 down to `--moho-model`
  40 km, mantle 8.04 km/s and 1.79 below). At every depth the amplitude is placed at the
  conversion point, which lies on the source side of the station at the horizontal offset
  of the S leg. Amplitudes are summed in bins of `--bin-deg` (0.1 degrees) by `--dz`
  (0.5 km) and divided by the number of traces; the top `--mute-km` (3 km) is muted
  because it holds the direct pulse. The stored volume is smoothed over 3 by 3 bins
  horizontally, weighted by the number of traces, and a column of 3 by 3 bins needs at
  least 5 traces.
- Profiles and slices. Bins within `--half-width` (45 km) of a profile are projected
  onto it in steps of `--step` (10 km); a profile bin needs 3 or more traces, and a
  running mean over three bins is applied along the profile. Depth slices are means
  over plus and minus 2 km around `--slices` (25, 35, 45 km), drawn with the stations
  and, if the optional package is installed, the coastline.

Positive amplitudes are interfaces where the velocity increases downward (the Moho);
negative ones are decreases. Multiples of shallow interfaces appear as later arrivals
and are mapped too deep; they are the main source of misreading and are best recognized
by comparing with the H-k result of the stations.

## 7. Two runs compared (`compare-rf`)

For two receiver-function folders of the same data: the functions computed, accepted by
the rules, removed by the station screen and kept; the stations with enough functions;
the quality classes; and per station the Ps time, H-k depth, Vp/Vs and class in both
runs with their differences. The figure shows the station stacks of the two runs and the
Ps time and H-k depth of one run against the other.

## 8. Crust table and crustal correction (`crust-from-rf`, `crust`)

- Crust table. One line per station with the Moho depth, its uncertainty and the
  `method` behind it. Default rule (`--prefer auto`): the H-k depth where the maximum is
  well constrained (quality A or B); otherwise the published value where one exists
  (`--published`, a CSV with station, moho_km, moho_err_km, source); otherwise the
  Ps-only depth, flagged uncertain with an error of at least `--ps-err` (5 km).
  `--prefer published` takes the published value wherever it exists (also for manual
  readings written into such a file), `--prefer rf` never uses it. The receiver-function
  values are compared with the published ones in `moho_vs_published.png`. With `--orla`
  the table also carries the mean Vs of the top 5 km, the mean Vs of the column and the
  vertical P time through 0 to 33 km from a 3-D shear-velocity model (ORLA2022 of Weidle
  et al. 2022, with Vp = 1.75 Vs). `ophiolite_km` and `sediment_km` are 0 unless edited.
- Crustal correction. For every reliable station: the vertical travel time through
  the column below sea level (sediments at 4.0 km/s, ophiolite at 7.0, crust at 6.3 down
  to the Moho, mantle at 8.04 down to a common reference depth) minus the same column in
  ak135, at the mean ray parameter of the station. The crustal term is this delay minus
  its network mean; the mantle residual is the observed term minus the crustal term. Its
  error combines the bootstrap error of the term with the uncertainty of the Moho (3 km
  is about 0.1 s). The residual is also expressed as a Vp anomaly over `--mantle-km`
  (150 km) below the Moho.
- 3-D crust. With `--orla` the uniform crust above 33 km is replaced by the Vs column
  of the model (Vp = kappa Vs); the two residuals and their difference measure how much
  the result depends on the crustal model.
- Profiles. The residual is plotted against longitude and latitude, and its means
  west and east of the median longitude are printed.
- Robustness. Class A and B picks only (`station-terms --min-quality B --subdir ...`,
  then `crust --terms`), `--prefer published` for the Moho, `--vp-crust 6.1` and `6.5`,
  each in its own `--subdir`.

## 9. Synthetic test (`demo`)

A one-layer crust over a half-space. For a set of ray parameters the vertical component
is a band-limited P pulse with noise, and the radial component is built from delayed
copies of it: the direct P, the Ps conversion at the Moho and the two crustal multiples.
The receiver functions are recovered with the iterative deconvolution and stacked in
H-k; the test passes when the depth and the Vp/Vs that were put in are recovered. It
checks the numerical core of step 7 without any data.

## References

- Ismaiel, M., Rishal, M. and Patel, R. (2026). Moho topography and flexural response of Semail ophiolite in the Southern Oman Mountains: new constraints from teleseismic receiver functions and gravity anomalies. Geophysical Journal International 245, ggag025. https://doi.org/10.1093/gji/ggag025
- Kennett, B. L. N., Engdahl, E. R. and Buland, R. (1995). Constraints on seismic velocities in the Earth from traveltimes. Geophysical Journal International 122, 108-124.
- Ligorria, J. P. and Ammon, C. J. (1999). Iterative deconvolution and receiver-function estimation. Bulletin of the Seismological Society of America 89, 1395-1400.
- Weidle, C., Wiesenberg, L., El-Sharkawy, A., Krüger, F., Scharf, A., Agard, P. and Meier, T. (2022). A 3-D crustal shear wave velocity model and Moho map below the Semail Ophiolite, eastern Arabia. Geophysical Journal International 231, 817-834.
- Zhu, L. and Kanamori, H. (2000). Moho depth variation in southern California from teleseismic receiver functions. Journal of Geophysical Research 105, 2969-2980.
