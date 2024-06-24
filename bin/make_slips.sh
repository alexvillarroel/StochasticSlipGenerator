#!/bin/bash

# Rango de magnitudes
cd /home/alex/Desktop/StochasticSlipGenerator/util_functions
start_magnitude=8.1
end_magnitude=8.6
increment=0.1

# Otros parámetros
region="-77 -69 -37 -29"
northlat=-29
southlat=-36
dx=10
dy=10
n_slip=1000

# Iterar sobre el rango de magnitudes
current_magnitude=$start_magnitude
while (( $(echo "$current_magnitude <= $end_magnitude" | bc -l) )); do
    echo "Ejecutando simulaciones para Mw = $current_magnitude"
    python Slip_prueba2.py -r $region -nlat $northlat -slat $southlat -dx $dx -dy $dy -n $n_slip -m $current_magnitude
    current_magnitude=$(echo "$current_magnitude + $increment" | bc)
done

