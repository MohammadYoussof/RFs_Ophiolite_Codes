# Data

| File | Content |
|---|---|
| `moho_published.csv` | Moho depths at 12 COOL stations from H-k stacking (Vp 6.5 km/s) by Ismaiel et al. (2026), and the value of COO46 for the co-located station COO07; used by `rfo.py crust-from-rf --published` for comparison and as a fallback |

Columns of `moho_published.csv`: `station`, `moho_km`, `moho_err_km`, `source`. Any table
with these four columns can be given to `--published`, also one with manual readings.

Reference: Ismaiel, M., Rishal, M. and Patel, R. (2026). Moho topography and flexural
response of Semail ophiolite in the Southern Oman Mountains: new constraints from
teleseismic receiver functions and gravity anomalies. Geophysical Journal International
245, ggag025. https://doi.org/10.1093/gji/ggag025

## Not stored here

Waveforms. Download them with SOD into this folder (see [../sod/README.md](../sod/README.md)):

```
data/raw_sac/Event_YYYY_MM_DD_hh_mm_ss/*.sac
data/raw_mseed/Event_YYYY_MM_DD_hh_mm_ss/*.mseed
```

Both folders are listed in `.gitignore`. The COOL network (5H, 2013 to 2016) is open at
GEOFON: Weidle, C., Agard, P., Meier, T., Ducassou, C. and El-Hussain, I. (2013). COOL
seismic network. GEOFON Data Archive. https://doi.org/10.7914/SN/5H_2013

Station metadata. StationXML at response level:

```
curl -o stations.xml "https://geofon.gfz.de/fdsnws/station/1/query?network=5H&starttime=2013-10-01&endtime=2016-03-01&level=response&format=xml"
```

3-D crustal model (optional). `rfo.py crust-from-rf --orla` and `rfo.py crust --orla`
read the shear-velocity model ORLA2022 of Weidle et al. (2022) in netCDF, available from
the EarthScope (IRIS) Earth Model Collaboration:
https://ds.iris.edu/ds/products/emc-orla2022
