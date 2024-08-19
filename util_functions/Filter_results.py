import os
import time
from scipy.io import loadmat
from geostochpy import modfilters
import geostochpy
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt
import csv

# El filtrado se basa en la física de la subducción,
# la profundidad de los deslizamientos y los límites geográficos.
# Se comienza filtrando los datos por estadística, con valores máximos de slip entre 
# 2 y 3 std

#
#
# folders=['/mnt/c/Users/axlph/OneDrive - Universidad de Concepción/magister/Proyecto de Tesis/Slips/Simulation_8.8/','/mnt/c/Users/axlph/OneDrive - Universidad de Concepción/magister/Proyecto de Tesis/Slips/Simulation_8.9/','/mnt/c/Users/axlph/OneDrive - Universidad de Concepción/magister/Proyecto de Tesis/Slips/Simulation_9.0/','/mnt/c/Users/axlph/OneDrive - Universidad de Concepción/magister/Proyecto de Tesis/Slips/Simulation_9.1/','/mnt/c/Users/axlph/OneDrive - Universidad de Concepción/magister/Proyecto de Tesis/Slips/Simulation_9.2/','/mnt/c/Users/axlph/OneDrive - Universidad de Concepción/magister/Proyecto de Tesis/Slips/Simulation_9.3/']
folders=['/mnt/c/Users/axlph/OneDrive - Universidad de Concepción/magister/Proyecto de Tesis/Slip_generation/Output_data/Simulation_9.0/']
# Get a list of all .mat files in the specified folder
for folder in folders:
    mat_files = [f for f in os.listdir(folder) if f.endswith('.mat')]
    start_time = time.time()  # Start timing
    # Load all .mat files
    mat_data = [loadmat(os.path.join(folder, f)) for f in mat_files]
    slip_data = [data['slip'] for data in mat_data if 'slip' in data]
    maximous = [np.max(slip) for slip in slip_data]
    slip_data = np.concatenate(slip_data)
    mean = np.mean(slip_data)
    std_dev = np.std(slip_data)
    min_val = np.min(slip_data)
    max_val = np.max(slip_data)
    median = np.median(slip_data)
    std_dev = np.std(slip_data)

    # Write statistics to file
    with open(folder+'stats.csv', 'w', newline='') as stats_file:
        writer = csv.writer(stats_file)
        writer.writerow(["Min", "Max", "Mean", "Median", "Standard Deviation"])
        writer.writerow([min_val, max_val, mean, median, std_dev])
    filter_trues=[]
    idx=0
    for data in mat_data:
        idx+=1
        X_grid=data['lat'].flatten()
        Y_grid=data['lon'].flatten()
        Slip=data['slip'].flatten()
        depth=data['depth'].flatten()
        filt=modfilters.physical_filter(Slip,Y_grid,depth,10000,20000,-33.5,-31)
        if filt==True and np.max(Slip) > (np.mean(maximous)+3*np.std(maximous))  and np.min(Y_grid)<=-33 and np.max(Y_grid)>=-31:
            filter_trues.append(idx)
    filter_trues = [str(idx).zfill(5) for idx in filter_trues]
    print(filter_trues)
    with open(folder+'filter_trues.txt', 'w') as f:
        for item in filter_trues:
            f.write("sim_%s.mat\n" % item)
    end_time = time.time()  # End timing
    print(f"Time taken to load all .mat files: {end_time - start_time} seconds")
    import os
import time
from scipy.io import loadmat
from geostochpy import modfilters
import geostochpy
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt
import csv

# El filtrado se basa en la física de la subducción,
# la profundidad de los deslizamientos y los límites geográficos.
# Se comienza filtrando los datos por estadística, con valores máximos de slip entre 
# 2 y 3 std

#
#
# folders=['/mnt/c/Users/axlph/OneDrive - Universidad de Concepción/magister/Proyecto de Tesis/Slips/Simulation_8.8/','/mnt/c/Users/axlph/OneDrive - Universidad de Concepción/magister/Proyecto de Tesis/Slips/Simulation_8.9/','/mnt/c/Users/axlph/OneDrive - Universidad de Concepción/magister/Proyecto de Tesis/Slips/Simulation_9.0/','/mnt/c/Users/axlph/OneDrive - Universidad de Concepción/magister/Proyecto de Tesis/Slips/Simulation_9.1/','/mnt/c/Users/axlph/OneDrive - Universidad de Concepción/magister/Proyecto de Tesis/Slips/Simulation_9.2/','/mnt/c/Users/axlph/OneDrive - Universidad de Concepción/magister/Proyecto de Tesis/Slips/Simulation_9.3/']
folders=['/home/alex/StochasticSlipGenerator/Output_data/Simulation_9.3_1000_coupling/']
filename='simulations.mat'
# Get a list of all .mat files in the specified folder
for folder in folders:
    data = loadmat(folders+filename)['simulations']
    data = np.squeeze(data)
    mat_files = [f for f in os.listdir(folder) if f.endswith('.mat')]
    start_time = time.time()  # Start timing
    # Load all .mat files
    mat_data = [loadmat(os.path.join(folder, f)) for f in mat_files]
    slip_data = [data['slip'] for data in mat_data if 'slip' in data]
    maximous = [np.max(slip) for slip in slip_data]
    slip_data = np.concatenate(slip_data)
    mean = np.mean(slip_data)
    std_dev = np.std(slip_data)
    min_val = np.min(slip_data)
    max_val = np.max(slip_data)
    median = np.median(slip_data)
    std_dev = np.std(slip_data)

    # Write statistics to file
    with open(folder+'stats.csv', 'w', newline='') as stats_file:
        writer = csv.writer(stats_file)
        writer.writerow(["Min", "Max", "Mean", "Median", "Standard Deviation"])
        writer.writerow([min_val, max_val, mean, median, std_dev])
    filter_trues=[]
    idx=0
    for data in mat_data:
        idx+=1
        X_grid=data['lat'].flatten()
        Y_grid=data['lon'].flatten()
        Slip=data['slip'].flatten()
        depth=data['depth'].flatten()
        filt=modfilters.physical_filter(Slip,Y_grid,depth,10000,20000,-33.5,-31)
        if filt==True and np.max(Slip) > (np.mean(maximous)+3*np.std(maximous))  and np.min(Y_grid)<=-33 and np.max(Y_grid)>=-31:
            filter_trues.append(idx)
    filter_trues = [str(idx).zfill(5) for idx in filter_trues]
    print(filter_trues)
    with open(folder+'filter_trues.txt', 'w') as f:
        for item in filter_trues:
            f.write("sim_%s.mat\n" % item)
    end_time = time.time()  # End timing
    print(f"Time taken to load all .mat files: {end_time - start_time} seconds")