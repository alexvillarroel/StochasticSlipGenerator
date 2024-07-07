import os
import time
import numpy as np
from scipy.io import loadmat
from geostochpy import modfilters
import matplotlib.pyplot as plt
import csv

# Carpeta donde está el archivo simulations.mat
folder = '/home/alex/StochasticSlipGenerator/Output_data/Simulation_9.3_1000_coupling/'

# Cargar el archivo simulations.mat
mat_contents = loadmat(os.path.join(folder, 'simulations.mat'))

# Obtener la variable que contiene las simulaciones (asumiendo que se llama 'simulations')
simulations = mat_contents['simulations']
simulations = np.squeeze(simulations)
# Inicializar listas para estadísticas y resultados filtrados
maximous = []
filter_trues = []

# Iterar sobre cada simulación en simulations
for idx, sim in enumerate(simulations):
    # Extraer datos de la simulación
    Slip = sim['slip'].flatten()
    Y_grid = sim['lon'].flatten()
    depth = sim['depth'].flatten()

    # Calcular estadísticas de slip
    maximous.append(np.max(Slip))
    mean = np.mean(Slip)
    std_dev = np.std(Slip)
    min_val = np.min(Slip)
    max_val = np.max(Slip)
    median = np.median(Slip)

    # Escribir estadísticas al archivo
    with open(os.path.join(folder, f'stats_{idx}.csv'), 'w', newline='') as stats_file:
        writer = csv.writer(stats_file)
        writer.writerow(["Min", "Max", "Mean", "Median", "Standard Deviation"])
        writer.writerow([min_val, max_val, mean, median, std_dev])

    # Aplicar filtro físico y estadístico
    filt = modfilters.physical_filter(Slip, Y_grid, depth, 10000, 20000, -33.5, -31)
    if filt and np.max(Slip) > (np.mean(maximous) + 3 * np.std(maximous)) and np.min(Y_grid) <= -33 and np.max(Y_grid) >= -31:
        filter_trues.append(idx)

# Escribir resultados filtrados a archivo
filter_trues = [str(idx).zfill(5) for idx in filter_trues]
with open(os.path.join(folder, 'filter_trues.txt'), 'w') as f:
    for item in filter_trues:
        f.write(f"sim_{item}.mat\n")

print(filter_trues)
