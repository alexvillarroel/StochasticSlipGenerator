import os
import re
import pygmt
import numpy as np
import geostochpy
import rockhound as rh
from scipy.io import loadmat
from scipy.interpolate import griddata

def get_simulation_number(filename):
    match = re.search(r'(\d+)', filename)
    return int(match.group(1)) - 1 if match else None

def plot_simulation(simulation, sim_num, lonfosa, latfosa, region, imgfolder):
    x = np.array(simulation[sim_num].lon).flatten()
    y = np.array(simulation[sim_num].lat).flatten()
    z = np.array(simulation[sim_num].slip).flatten()
    region_slip = [np.min(x), np.max(x), np.min(y), np.max(y)]

    grid=pygmt.nearneighbor(x=x, y=y, z=z,region=region_slip,spacing="1k",search_radius="40k")
    earth_grid = pygmt.datasets.load_earth_relief(resolution="30s", region=region)
    fig = pygmt.Figure()
    fig.basemap(region=region, projection='M12c', frame=['WSne', 'y1+laatitude(°)', 'x2+laongitude(°)', 'g'])
    fig.grdimage(grid=earth_grid, cmap='geo', shading=True, dpi=300)
    fig.colorbar(cmap=True, frame=["x+lElevation", "y+lm"],)
    cmap = pygmt.makecpt(cmap='jet', series=[0, np.max(np.round(z)+2), np.max(np.round(z)+2)/10], continuous=True)

    # Agregar eventos
    events = [
        {"x": -77.6, "y_range": [-29, -36], "text_x": -78, "text_y": -32.5, "label": "Valparaíso 1730"},
        {"x": -77, "y_range": [-34, -38.1], "text_x": -77.3, "text_y": -36, "label": "Concepción 1751"},
        {"x": -76, "y_range": [-32.5, -34.5], "text_x": -76.5, "text_y": -33.5, "label": "Valparaíso 1822,1906,1985"},
        {"x": -75, "y_range": [-30, -31.6], "text_x": -75.5, "text_y": -31, "label": "Illapel 1880,1943,2015"},
        {"x": -74, "y_range": [-34.3, -35.5], "text_x": -74.5, "text_y": -35, "label": "Talca 1928"},
        {"x": -74, "y_range": [-28, -30], "text_x": -74.5, "text_y": -29, "label": "Atacama 1922"},
        {"x": -75, "y_range": [-34, -38.1], "text_x": -75.5, "text_y": -36, "label": "Maule 2010"},
        {"x": -76, "y_range": [-35.2, -37.8], "text_x": -76.5, "text_y": -37, "label": "Concepción 1835"}
    ]

    for event in events:
        fig.plot(x=[event["x"], event["x"]], y=event["y_range"], fill='red', pen='2,red')
        fig.text(x=event["text_x"], y=event["text_y"], text=event["label"], fill='white', font="10p,Helvetica-Bold,black", angle=90)

    fig.grdimage(grid=grid, cmap=cmap, nan_transparent=True, dpi=300)

    depth_grid = rh.fetch_slab2('south_america').depth / -1000
    fig.grdcontour(grid=depth_grid, region=region, levels=20, annotation='40+e+f10p+gwhite')
    fig.coast(shorelines="1p,black", borders=["1/0.5p,black", "2/0.5p,gray", "3/0.5p,blue"])

    fig.plot(x=lonfosa, y=latfosa,
             projection='M12c',
             region=region,
             pen="1p",
             fill="white",
             style="f0.5i/0.1i+r+t+o1")
    fig.colorbar(
        cmap=cmap,
        position="g-79.8/-39.8+w6c/0.5c+v",
        box='+ggray+pblack',
        frame=["x+lCoseismic Slip (m)"],
    )

    fig.savefig(os.path.join(imgfolder, f'simulation_{sim_num+1:04d}.png'), dpi=300)

# Leer el archivo de la fosa
file_trench = geostochpy.get_data('trench-chile.txt')
trench = np.genfromtxt(file_trench, delimiter=" ")
lonfosa = trench[:, 0]
latfosa = trench[:, 1]
region = [-80, -69, -40, -27]

# Cargar los datos de la simulación
for i in np.arange(8.1,9.3,0.1):
    version = np.round(i,decimals=1)
    folder = f'/home/alex/StochasticSlipGenerator/Output_data/Simulation_{version}_1000_coupling/'
    imgfolder = f'/home/alex/StochasticSlipGenerator/Output_data/Simulation_{version}_1000_coupling/img_filtered'
    simulation = loadmat(folder + 'simulations.mat', struct_as_record=False, squeeze_me=True)
    simulation = simulation['simulations']

    # Listar los archivos en la carpeta de imágenes
    img_files = [f for f in os.listdir(imgfolder) if re.match(r'simulation_\d+\.png', f)]
    print(img_files)
    # Iterar sobre cada archivo de imagen y generar la figura
    for img_file in img_files:
        sim_num = get_simulation_number(img_file)
        if sim_num is not None:
            plot_simulation(simulation, sim_num, lonfosa, latfosa, region, imgfolder)
