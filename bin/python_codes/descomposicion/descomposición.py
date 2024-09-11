import os
import time
import numpy as np
from scipy.io import loadmat
from scipy.interpolate import griddata
from geostochpy import modfilters
import matplotlib.pyplot as plt
import csv
import geostochpy
import pygmt
import rockhound as rh
def plot_pygmt(X,Y,Z,latfosa,lonfosa,output):

    region_slip=[np.min(X),np.max(X),np.min(Y),np.max(Y)]
    grid=pygmt.nearneighbor(x=X,y=Y,z=Z,region=region_slip,spacing='1k',search_radius='40k')
    region=[-80,-69,-40,-27]
    region2=[-76,-68,-36,-28]
    earth_grid = pygmt.datasets.load_earth_relief(resolution="30s", region=region)
    #
    fig=pygmt.Figure()
    fig.basemap(region=region,projection='M12c',frame=['WSne','y1+laatitude(°)','x2+laongitude(°)','g'])
    fig.grdimage(grid=earth_grid,cmap='geo',shading=True,dpi=300)
    fig.colorbar(cmap=True,frame=["x+lElevation","y+lm"])
    #cmap=pygmt.makecpt(cmap='seis',series=[0, np.max(np.round(Z)), 3],continuous=True,reverse=True)
    if np.max(np.round(Z))>1:
        cmap=pygmt.makecpt(cmap='jet',series=[0, np.max(np.round(Z)), 3],continuous=True)
    else:
        cmap=pygmt.makecpt(cmap='matlab/hot',series=[0, 1],continuous=False,reverse=True)
    # fig.coast(shorelines=True, area_thresh=5000,land="gray")
    # valparaiso 1730
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

    ##
    fig.grdimage(grid=grid,cmap=cmap,nan_transparent=True,dpi=300,shading=False)
    # fig.grdimage(grid=file2,cmap=cmap,nan_transparent=True)
    #
    depth_grid=rh.fetch_slab2('south_america').depth/-1000
    fig.grdcontour(grid=depth_grid,region=region,levels=20,annotation='40+e+f8p+')
    fig.coast(shorelines="1p,black",borders=["1/0.5p,black", "2/0.5p,gray", "3/0.5p,blue"])

    fig.plot(x=lonfosa,y=latfosa,
            projection='M12c',
            region=region,
            pen="1p",
            fill="white",
            style="f0.5i/0.1i+r+t+o1")
    if np.max(np.round(Z))>1:
        fig.colorbar(
            cmap=cmap,
            # Colorbar positioned at map coordinates (g) longitude/latitude 0.3/8.7,
            # with a length/width (+w) of 4 cm by 0.5 cm, and plotted horizontally (+h)
            position="g-79.8/-39.8+w6c/0.5c+v",
            box='+ggray+pblack',
            frame=["x+lCoseismic Slip (m)"],
        )
    else:
        fig.colorbar(
            cmap=cmap,
            # Colorbar positioned at map coordinates (g) longitude/latitude 0.3/8.7,
            # with a length/width (+w) of 4 cm by 0.5 cm, and plotted horizontally (+h)
            position="g-79.8/-39.8+w6c/0.5c+v",
            box='+ggray+pblack',
            frame=["x+lTaper(coupling degree)"],)

    cities = {
        "        Valparaíso": {"lon": -71.63, "lat": -33.03},
        "        La Serena": {"lon": -71.4, "lat": -30.03},
        "Santiago": {"lon": -70.6, "lat": -33.45},
        "  Talca": {"lon": -71.8, "lat": -35.43},
        "  Concepción": {"lon": -73.039, "lat": -36.812}
    }
    for city, coords in cities.items():
        fig.plot(x=[coords["lon"]], y=[coords["lat"]], style="kcity/0.2c", pen="0.2p,black",fill='black')
        fig.text(x=coords["lon"], y=coords["lat"], text=city, font="8p,Helvetica-Bold,white",fill='black', justify="ML")
    fig.savefig(output,anti_alias=True,dpi=300)
    return
route_trench = geostochpy.get_data('trench-chile.txt')
lonsfosa, latsfosa, strikefosa = geostochpy.load_trench(route_trench)
slabdep, slabdip, slabstrike, slabrake = geostochpy.load_files_slab2(zone='south_america', rake=True)
route_file=geostochpy.get_data('median_lock_Herrera2023.txt')
route_file_mesh=geostochpy.get_data('mesh__Herrera2023.npy')
mesh=np.load(route_file_mesh,allow_pickle=True)
median_lock=np.loadtxt(route_file)
mesh1=mesh[0]
x=mesh1[:,0]
y=mesh1[:,1]
z=median_lock
slabcoupling=np.column_stack((x, y, z))
# datos
#Carpeta donde está el archivo simulations.mat
folder = '/home/alex/StochasticSlipGenerator/Output_data/Simulation_8.8_5000_coupling/'

simulation=loadmat(folder+'simulations.mat',struct_as_record=False, squeeze_me=True)
simulation=simulation['simulations']
#
n_sim=4173
X_grid=simulation[n_sim].lon
Y_grid=simulation[n_sim].lat
Slip=simulation[n_sim].slip
taper_coupling= griddata((slabcoupling[:,0], slabcoupling[:,1]), slabcoupling[:,2], (X_grid, Y_grid), method='linear', fill_value=0, rescale=True)
slip_sin_taper = np.where(taper_coupling != 0, Slip / taper_coupling, 0)
plot_pygmt(X_grid.flatten(),Y_grid.flatten(),taper_coupling.flatten(),latsfosa,lonsfosa,'taper.png')
plot_pygmt(X_grid.flatten(),Y_grid.flatten(),Slip.flatten(),latsfosa,lonsfosa,'slip.png')
plot_pygmt(X_grid.flatten(),Y_grid.flatten(),slip_sin_taper.flatten(),latsfosa,lonsfosa,'slip_sin_taper.png')
geostochpy.plot_slip(X_grid,Y_grid,lonsfosa,latsfosa,taper_coupling,'taper2.png',show=False,cmap='rainbow')
geostochpy.plot_slip(X_grid,Y_grid,lonsfosa,latsfosa,Slip,'slip2.png',show=False,cmap='rainbow')
geostochpy.plot_slip(X_grid,Y_grid,lonsfosa,latsfosa,slip_sin_taper,'slip_sin_taper2.png',show=False,cmap='rainbow')
