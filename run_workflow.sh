#!/bin/bash
# All steps of the workflow in order. Use it as a template: set the names below, then run
# the steps one at a time the first time and look at the results before going on.
# docs/WORKFLOW.md explains every step.
set -e

SAC=data/raw_sac                 # original SAC event folders (read only)
MSEED=data/raw_mseed             # original miniSEED event folders (read only)
XML=stations.xml                 # StationXML at response level
QC=qc_out                        # output of the quality control
WORKERS=6                        # events processed in parallel
MANIFEST=originals_sha256.txt

# ---- step 1: download (done once; see sod/README.md)
# mkdir -p data && cd data && sod -f ../sod/cool_5H_full.xml && cd ..

# ---- step 2: protect the originals
chmod -R a-w "$SAC" "$MSEED"
[ -f "$MANIFEST" ] || python rfo.py originals make "$SAC" "$MSEED" -m "$MANIFEST"

# ---- step 3: quality control and picking
python rfo.py qc --data "$SAC" --inventory "$XML" --out "$QC" --workers "$WORKERS"
python rfo.py originals verify -m "$MANIFEST"

# ---- step 4: overview of the data set
python rfo.py inventory  --qc "$QC"
python rfo.py snr        --qc "$QC"
python rfo.py response   --qc "$QC" --inventory "$XML"
python rfo.py noise-hour --qc "$QC" --utc-offset 4

# ---- step 5: sensor orientation
python rfo.py orientation       --qc "$QC" --data "$SAC" --inventory "$XML" --workers "$WORKERS"
python rfo.py orientation-apply --qc "$QC" --inventory "$XML" --out-xml stations_oriented.xml
# If sensors were corrected, run step 3 again with the corrected file into a new folder
# and continue with that folder:
#   python rfo.py qc --data "$SAC" --inventory stations_oriented.xml --out qc_out2 --workers "$WORKERS"
#   QC=qc_out2

# ---- step 6: station terms
#      add --exclude STATION1,STATION2 for stations judged unusable in the quality control
python rfo.py station-terms --qc "$QC"

# ---- step 7: receiver functions
python rfo.py rf --qc "$QC" --gauss 2.0 --workers "$WORKERS"
# another Gaussian width, and the comparison of the two runs:
#   python rfo.py rf --qc "$QC" --gauss 2.5 --out rf_g25 --workers "$WORKERS"
#   python rfo.py compare-rf --qc "$QC" --a rf --b rf_g25 --label-a "width 2.0" --label-b "width 2.5"

# ---- step 9: crust table (before step 8, so that the images carry the Moho markers)
python rfo.py crust-from-rf --qc "$QC" --published data/moho_published.csv

# ---- step 8: images
python rfo.py section --qc "$QC"
python rfo.py ccp     --qc "$QC"

# ---- step 10: crustal correction
python rfo.py crust --qc "$QC" --model "$QC/analysis/crust_model.csv"

python rfo.py originals verify -m "$MANIFEST"
