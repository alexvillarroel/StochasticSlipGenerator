#!/bin/bash

# Rango de magnitudes
cd /home/alex/StochasticSlipGenerator/util_functions
start_magnitude=8.1
end_magnitude=9.3
increment=0.1

# Otros parámetros
region="-77 -69 -37 -29"
northlat=-29
southlat=-36
n_subfaults=500
n_slip=1000

# Iterar sobre el rango de magnitudes
current_magnitude=$start_magnitude
while (( $(echo "$current_magnitude <= $end_magnitude" | bc -l) )); do
    echo "Ejecutando simulaciones para Mw = $current_magnitude"
    python Slipgenerator_coupling.py -r $region -nlat $northlat -slat $southlat -ns $n_subfaults -n $n_slip -m $current_magnitude
    current_magnitude=$(echo "$current_magnitude + $increment" | bc)
done

