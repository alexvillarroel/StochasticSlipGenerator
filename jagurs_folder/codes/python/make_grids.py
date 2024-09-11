import pygmt
import pandas as pd
import numpy as np

# Leer el archivo Excel asegurándose de que las columnas se lean como números
df_excel = pd.read_excel('/home/alex/StochasticSlipGenerator/jagurs_folder/grids/BahiaValparaiso.xlsx', header=None, names=['latitud', 'longitud', 'profundidad'], dtype=float, skiprows=8)
df_xyz = pd.read_csv('/home/alex/StochasticSlipGenerator/jagurs_folder/grids/valparaiso_1s.xyz', sep='\s+', header=None, names=['longitud', 'latitud', 'profundidad'])

# Restringir los valores de latitud y longitud a 6 decimales
df_excel['latitud'] = df_excel['latitud'].round(6)
df_excel['longitud'] = df_excel['longitud'].round(6)
df_xyz['latitud'] = df_xyz['latitud'].round(6)
df_xyz['longitud'] = df_xyz['longitud'].round(6)

# Obtener los límites de la región de df_excel
region_excel = [df_excel['longitud'].min(), df_excel['longitud'].max(), df_excel['latitud'].min(), df_excel['latitud'].max()]

# Obtener los límites de la región de df_xyz
region_xyz = [df_xyz['longitud'].min(), df_xyz['longitud'].max(), df_xyz['latitud'].min(), df_xyz['latitud'].max()]
# Filtrar las filas con valores de 0 en la profundidad
df_xyz = df_xyz[df_xyz['profundidad'] != 0]
df_xyz=pd.concat([df_xyz,df_excel],ignore_index=True)
# Definir el espaciado
spacing = 2/3  # 2/3 segundos de arco

# Ajustar los límites de la región para que sean múltiplos del espaciado
region_xyz[0] = np.floor(region_xyz[0] / spacing) * spacing
region_xyz[1] = np.ceil(region_xyz[1] / spacing) * spacing
region_xyz[2] = np.floor(region_xyz[2] / spacing) * spacing
region_xyz[3] = np.ceil(region_xyz[3] / spacing) * spacing

# Crear un grid usando blockmean
grid = pygmt.blockmean(data=df_xyz, region=region_xyz,spacing=f'{spacing}s')

# Crear una superficie basada en el grid calculado
valpo_grid = pygmt.surface(data=grid, region=region_xyz, spacing=f'{spacing}s',outgrid='valpo_grid.grd')


