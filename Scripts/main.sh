#!/bin/bash
chmod +x Scripts/chim_cov.sh
mkdir -p DerivedData
mkdir -p DerivedData/chimers
mkdir -p DerivedData/coverages
./Scripts/chim_cov.sh
mkdir -p Results/
mkdir -p Results/my_rearrangments
mkdir -p Results/mutations
mkdir -p Results/imlicit_restr_sites
python3 ./Scripts/rearrangments_search.py 14-41
