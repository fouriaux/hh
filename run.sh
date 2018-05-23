#!/usr/bin/env bash

i_range=(
'0.000000005'
'0.00000005'
'0.0000005'
)

for i in ${i_range[@]}; do
./bin/fe $i > ./error_analysis/fe_$i.csv
./bin/cn $i > ./error_analysis/cn_$i.csv
done
