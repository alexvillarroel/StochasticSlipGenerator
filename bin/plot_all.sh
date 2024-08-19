#!/bin/bash
cd python_codes
# Define el rango de simulaciones
start=8.2
end=9.3
step=0.1

# Establecer la configuración regional a una que use punto como separador decimal
export LC_NUMERIC="en_US.UTF-8"

# Ejecutar el script de Python para cada número en el rango
for sim_number in $(seq $start $step $end); do
    python plot.py $sim_number
done

# Restaurar la configuración regional original
unset LC_NUMERIC

