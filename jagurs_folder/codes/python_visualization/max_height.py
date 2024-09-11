import pygmt
import netCDF4 as nc
import numpy as np
import xarray as xr

# Abrir el archivo NetCDF
dataset = nc.Dataset('/home/alex/Tsunamis/template/SD01.nc')

# Obtener las variables
lon = dataset.variables['lon'][:]
lat = dataset.variables['lat'][:]
wave_height = dataset.variables['wave_height'][:, :, :]  # wave height en todas las dimensiones (lon, lat, time)

# Seleccionar el primer momento en el tiempo para graficar
wave_height_at_time = wave_height[:, :, 0]  # Por ejemplo, el primer índice de tiempo

# Convertir a matriz de NumPy (en caso de que sea una MaskedArray)
wave_height_at_time = np.array(wave_height_at_time)

# Asegurarse de que lon y lat tengan las mismas dimensiones que los datos de wave_height
lon, lat = np.meshgrid(lon, lat)

# Crear un DataArray de xarray para que PyGMT pueda manejarlo como una cuadrícula
grid = xr.DataArray(wave_height_at_time, coords=[('lat', lat[:,0]), ('lon', lon[0,:])])

# Crear la figura
fig = pygmt.Figure()

# Definir la región y proyección
region = [lon.min(), lon.max(), lat.min(), lat.max()]
projection = 'M6i'

# Crear la paleta de colores
pygmt.makecpt(cmap='viridis', series=[wave_height_at_time.min(), wave_height_at_time.max()])

# Graficar los datos
fig.basemap(region=region, projection=projection, frame=True)
fig.grdimage(grid=grid, region=region, projection=projection, cmap='viridis')
fig.colorbar(frame='af+l"Wave height (m)"')

# Mostrar la figura
fig.show()

# Cerrar el archivo NetCDF
dataset.close()

