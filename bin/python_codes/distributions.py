import numpy as np
import scipy.io
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
from sklearn.linear_model import LinearRegression
# Definir rutas de los archivos
fullpaths = [
    '/home/alex/StochasticSlipGenerator/Output_data/Simulation_8.1_1000_coupling/simulations.mat',
    '/home/alex/StochasticSlipGenerator/Output_data/Simulation_8.2_1000_coupling/simulations.mat',
    '/home/alex/StochasticSlipGenerator/Output_data/Simulation_8.3_1000_coupling/simulations.mat',
    '/home/alex/StochasticSlipGenerator/Output_data/Simulation_8.4_1000_coupling/simulations.mat',
    '/home/alex/StochasticSlipGenerator/Output_data/Simulation_8.5_1000_coupling/simulations.mat',
    '/home/alex/StochasticSlipGenerator/Output_data/Simulation_8.6_1000_coupling/simulations.mat',
    '/home/alex/StochasticSlipGenerator/Output_data/Simulation_8.7_1000_coupling/simulations.mat',
    '/home/alex/StochasticSlipGenerator/Output_data/Simulation_8.8_1000_coupling/simulations.mat',
    '/home/alex/StochasticSlipGenerator/Output_data/Simulation_8.9_1000_coupling/simulations.mat',
    '/home/alex/StochasticSlipGenerator/Output_data/Simulation_9.0_1000_coupling/simulations.mat',
    '/home/alex/StochasticSlipGenerator/Output_data/Simulation_9.1_1000_coupling/simulations.mat',
    '/home/alex/StochasticSlipGenerator/Output_data/Simulation_9.2_1000_coupling/simulations.mat',
    '/home/alex/StochasticSlipGenerator/Output_data/Simulation_9.3_1000_coupling/simulations.mat'
]

# Número de simulaciones por archivo
num_simulations = 1000

# Inicializar variables
total_simulations = num_simulations * len(fullpaths)
lengths = np.zeros(total_simulations)
widths = np.zeros(total_simulations)
mean_slips = np.zeros(total_simulations)
magnitudes = np.zeros(total_simulations)

# Asignar valores a 'magnitudes' para cada grupo de simulaciones
for i in range(len(fullpaths)):
    magnitudes[i*num_simulations:(i+1)*num_simulations] = 8.0 + 0.1 * (i + 1)
# Procesar cada archivo de simulación
for i, path in enumerate(fullpaths):
    # Cargar datos de simulación
    data = scipy.io.loadmat(path)
    simulations = np.squeeze(data['simulations'])  # Ajustar según la estructura correcta
    # Indices para la ubicación en las matrices finales
    idx_start = i * num_simulations
    idx_end = (i + 1) * num_simulations
    # Extraer y vectorizar los datos necesarios
    lengths[idx_start:idx_end] = np.array([sim['length'][0][0][0][0] * sim['lat'][0][0].shape[0] for sim in simulations])
    widths[idx_start:idx_end] = np.array([sim['width'][0][0][0][0] * sim['lat'][0][0].shape[1] for sim in simulations])
    mean_slips[idx_start:idx_end] = np.array([np.mean(sim['slip'][0][0]) for sim in simulations])

# Crear el DataFrame
unique_magnitudes = np.unique(magnitudes)
df = pd.DataFrame({'Magnitudes': magnitudes, 'Lengths': np.log10(lengths), 'Widths': np.log10(widths)})

# Colors for violins
palette = sns.color_palette("Set2", len(unique_magnitudes))
sns.set_style("darkgrid", {"grid.color": ".6", "grid.linestyle": ":"})
fig, (ax1, ax2) = plt.subplots(ncols=1, nrows=2, figsize=(10, 7))
# Primer gráfico
sns.violinplot(x='Magnitudes', y='Lengths', data=df, palette=palette, alpha=0.4, hue='Magnitudes', legend=False,
               native_scale=True, ax=ax1)

# Ajustar la regresión lineal para Lengths
X = df['Magnitudes'].values.reshape(-1, 1)
y = df['Lengths'].values
reg = LinearRegression().fit(X, y)
slope, intercept = reg.coef_[0], reg.intercept_
equation = f'L = {slope:.2f}MW + {intercept:.2f}'
sns.regplot(x='Magnitudes', y='Lengths', data=df, scatter_kws={"alpha": 0.0}, line_kws={"label": equation}, ax=ax1)
ax1.legend(prop={'size': 15}, loc='upper left')
ax1.set_title('Violin Plot of Lengths vs Magnitude')
ax1.set_xticks(unique_magnitudes)  # make sure each x-position is labeled

# Segundo gráfico
sns.violinplot(x='Magnitudes', y='Widths', data=df, palette=palette, alpha=0.5, hue='Magnitudes', legend=False,
               native_scale=True, ax=ax2)

# Ajustar la regresión lineal para Widths
df2 = pd.DataFrame({'Magnitudes': magnitudes[0:7999], 'Lengths': np.log10(lengths[0:7999]), 'Widths': np.log10(widths[0:7999])})
X = df2['Magnitudes'].values.reshape(-1, 1)
y = df2['Widths'].values
reg = LinearRegression().fit(X, y)
slope, intercept = reg.coef_[0], reg.intercept_
equation = f'W = {slope:.2f}Mw + {intercept:.2f}'

sns.regplot(x='Magnitudes', y='Widths', data=df2, scatter_kws={"alpha": 0.0}, line_kws={"label": equation}, ax=ax2)
ax2.legend(prop={'size': 15}, loc='upper left')
ax2.set_title('Violin Plot of Widths vs Magnitude')
ax2.set_xticks(unique_magnitudes)
plt.show()