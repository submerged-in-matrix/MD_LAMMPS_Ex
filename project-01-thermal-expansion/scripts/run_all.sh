#!/bin/bash
for T in 100 200 300 400 500 600 700 800 900; do
    echo "Running T = ${T} K..."
    lmp -in in.cu_thermal -var T $T > ../outputs/log_${T}K.txt 2>&1
    echo "Done T = ${T} K"
done
echo "All temperatures complete."
