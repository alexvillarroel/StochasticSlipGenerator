import sys
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
matplotlib.use('Agg')
import geostochpy
from scipy.io import loadmat
from multiprocessing import Pool
from tqdm import tqdm

if len(sys.argv) != 2:
    print("Uso: python file.py <número_de_simulación>")
    sys.exit(1)

# Capturar el número de simulación desde los argumentos de la línea de comandos
sim_number = float(sys.argv[1])

# Formatear las variables filename y folder usando el número ingresado
filename = f'/home/alex/StochasticSlipGenerator/Output_data/Simulation_{sim_number}_1000_coupling/simulations.mat'
folder = f'/home/alex/StochasticSlipGenerator/Output_data/Simulation_{sim_number}_1000_coupling/img/'

data = loadmat(filename)['simulations']
data = np.squeeze(data)

route_trench = geostochpy.get_data('trench-chile.txt')
lonsfosa, latsfosa, strikefosa = geostochpy.load_trench(route_trench)

def process_cell(cell_data):
    X_grid = cell_data['lon'].item()
    Y_grid = cell_data['lat'].item()
    Slip = cell_data['slip'].item()
    k = cell_data['k']
    geostochpy.plot_slip(X_grid, Y_grid, lonsfosa, latsfosa, Slip, folder + f'simulation_{k:04d}.png', show=False)

if __name__ == '__main__':
    with Pool() as pool:
        # Crear una barra de progreso para el bucle usando tqdm
        results = list(tqdm(pool.imap(process_cell, [{'lon': cell['lon'], 'lat': cell['lat'], 'slip': cell['slip'], 'k': i+1} for i, cell in enumerate(data)]), total=len(data)))
