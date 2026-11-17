# Step 1: download with SOD

The event data are downloaded with SOD (Standing Order for Data, University of South
Carolina, https://www.seis.sc.edu/sod/). SOD works through three "arms": it reads an
event catalog (event arm), finds the stations and channels (network arm), and requests a
time window for every event and station that pass the distance test (waveform arm).
The recipes here were run with SOD 4.0.0rc1.

| File | Purpose |
|---|---|
| `cool_5H_test.xml` | one event at two stations; run this first (a few minutes) |
| `cool_5H_full.xml` | the whole COOL deployment, magnitude 5.5 and larger, 30 to 95 degrees |
| `verify_download.py` | compares a download, sample by sample, with the data center |

## Run

SOD writes into the folder it is started in, so start it in the folder that is to
hold the data (from the top of the repository):

```
sod -V -f sod/cool_5H_test.xml            # check the recipe

mkdir -p data_test && cd data_test
sod -f ../sod/cool_5H_test.xml            # run the test
python ../sod/verify_download.py          # every line should end with three times True
cd ..

mkdir -p data && cd data
sod -f ../sod/cool_5H_full.xml            # the full download (hours; SOD resumes if interrupted)
cd ..
```

SOD writes two copies of every seismogram, one folder per event:

```
raw_sac/Event_2015_05_30_11_23_02/5H2013.COO16..HHZ.sac
raw_mseed/Event_2015_05_30_11_23_02/5H2013.COO16..HHZ.mseed
```

The SAC files carry the event origin, location, depth, magnitude and the station
distance and back-azimuth in their headers, so the quality control reads the SAC
folder. The miniSEED folder holds the same samples; `python rfo.py mseed-metadata`
writes the event catalog and station table that a miniSEED user needs.

For a temporary network SOD puts the start year of the deployment into the file names
(`5H2013`), because temporary network codes are reused. The quality control maps this
back to the network code (`network_rename` in the `CONFIG` block of
`rf_ophiolite/qc_pipeline.py`).

## Station metadata

The quality control removes the instrument response, so it needs a StationXML file at
response level. From GEOFON:

```
curl -o stations.xml "https://geofon.gfz.de/fdsnws/station/1/query?network=5H&starttime=2013-10-01&endtime=2016-03-01&level=response&format=xml"
```

## Adapting the recipe to another network

| What | Where in the recipe |
|---|---|
| time range of the events | `eventArm` / `originTimeRange` |
| smallest magnitude | `eventArm` / `magnitudeRange` |
| data center | `host` in `fdsnStation` (network arm) and in `fdsnDataSelect` (waveform arm) |
| network, stations, channels | `networkCode`, `stationCode`, `channelCode` in `fdsnStation` |
| distance range | `waveformArm` / `distanceRange` |
| time window | `waveformArm` / `phaseRequest` (begin and end phase with offsets) |
| output folders | `mseedWriter` and `sacWriter` |

Always check a changed recipe with `sod -V -f recipe.xml` and try it on one event and
two stations before the full run.
