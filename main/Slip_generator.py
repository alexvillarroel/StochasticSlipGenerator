import argparse
import numpy as np 
import pygmt
import time
import modokada as mo 
import modfallas as mf 
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d,interp2d,LinearNDInterpolator
from geographiclib.geodesic import Geodesic
import json
import gc 
import os
import time
# change from Ms to Mw
def Mw_kauselramirez(ms):
    logm0kauselramirez = 1.5*ms+16.30
    Mw_calcs = 2.0/3.0*logm0kauselramirez-10.7 
    Mw_calcs = np.round(Mw_calcs,1) # 1 dec
    return Mw_calcs
# km2deg
def km2deg(km):
    return float(km/111.1)
#
def deg2km(deg):
    return float(deg*111.1)
#
def coincidencias_filas(A,B):
    coincidencias  =  [i for i in range(B.shape[0]) if np.any(np.all(A==B[i],axis=1))]
    if len(coincidencias)==0:
        return B[coincidencias]
    return np.unique(B[coincidencias],axis=0)
# main
def main(args):
    """What is executed upon running the program. 

        Runs for make a .xyz grid with NX x NY points  
        """
    # set command line variables (after python makegrid.py)
    # -29.500513696035412, -72.28001907930886
    plot_grid=args.plotgrid
    Mw=args.mw
    region=args.region
    width=args.width
    northlat=args.northlat
    southlat=args.southlat
    nx=args.nx
    ny=args.ny
    dx=width/nx
    n_slip=args.nslip
    dy=deg2km((northlat-southlat))/ny
    timestr = time.strftime("%Y%m%d-%H%M%S")
    # 
    np.disp('The width of your soon grid is: '+str(deg2km((northlat-southlat)))+' km ')
    points=nx*ny
    data_trench=np.loadtxt("../trench-chile.txt")
    # make grid
    grid=open('../Output_data/'+timestr+'_'+str(northlat)+str(southlat)+'.xyz','w')
    lats=[]
    lons=[]
    for i in range(nx):
        for j in range(ny):
            lat=northlat-km2deg(dy*(2*j+1)/2)
            lats.append(northlat-km2deg(dy*j))
            lon=np.interp(lat,data_trench[:,1],data_trench[:,0])+km2deg(dx*(2*i+1)/2)
            lons.append(np.interp(lat,data_trench[:,1],data_trench[:,0])+km2deg(dx*i))
            grid.write(str(lon)+','+str(lat)+','+'0'+'\n')
    grid.close()
    # lat abajo, lon al lado
    # ruta del archivo de la fosa
    ruta_fosa = "../Slab/"
    # archivo fosa ( primera columna: longitudes, segunda columna: latitudes)
    arc_fosa  = ruta_fosa + "trench-chile.txt"
    # carga de fosa usando funcion del modulo modfallas
    lonfosa, latfosa  = mf.carga_fosa(arc_fosa)

    directorio = "../Slab/"
    files=["sam_slab2_dep_02.23.18.xyz","sam_slab2_dip_02.23.18.xyz","sam_slab2_str_02.23.18.xyz","sam_rake.xyz"]
    # archivo de profundidad
    slabdep   = np.sort(np.genfromtxt(directorio+files[0], delimiter = ","),axis=0) # se lee el archivo a un array
    # archivo de dip
    slabdip    = np.sort(np.genfromtxt(directorio+files[0], delimiter = ","),axis=0) # se lee el archivo a un array
    # archivo de strike
    slabstrike = np.sort(np.genfromtxt(directorio+files[0], delimiter = ","),axis=0) # se lee el archivo a un array
    ## archivo de rake
    slabrake   = np.sort(np.genfromtxt(directorio+files[0], delimiter = ","),axis=0)
    # las longitudes estan en formato 0 - 360, se cambian a -180 - 180
    slabdep[:,0],slabdip[:,0],slabstrike[:,0]   = slabdep[:,0] - 360, slabdip[:,0] - 360,slabstrike[:,0] - 360
    # longitudes y latitudes del modelo
    lonmod,latmod = slabdep[:,0], slabdep[:,1]
    # longitudes y latitudes unicas del modelo
    lonunimod,latunimod = np.unique(lonmod),np.unique(latmod)
    X,Y=np.meshgrid(lonunimod,latunimod)
    minlat,maxlat=np.min(lats),np.max(lats)
    minlon,maxlon=np.min(lons),np.max(lons)
    X_grid,Y_grid=np.reshape(lons,(ny,nx)),np.reshape(np.sort(lats)[::-1].T,(ny,nx))
    #
    interp_dep=LinearNDInterpolator(list(zip(slabdep[:,0],slabdep[:,1])),slabdep[:,-1],fill_value=0)
    interp_dip=LinearNDInterpolator(list(zip(slabdip[:,0],slabdip[:,1])),slabdip[:,-1],fill_value=0)
    interp_strike=LinearNDInterpolator(list(zip(slabstrike[:,0],slabstrike[:,1])),slabstrike[:,-1],fill_value=0)
    interp_rake=LinearNDInterpolator(list(zip(slabrake[:,0],slabrake[:,1])),slabrake[:,-1],fill_value=0)
    dep=interp_dep(X_grid,Y_grid)
    dip=interp_dip(X_grid,Y_grid)
    strike=interp_strike(X_grid,Y_grid)
    rake=interp_rake(X_grid,Y_grid)

    ## Creation slip models
    # mean matrix
    mu   = mf.matriz_medias(11, dep)
    # matriz de covarianza
    C    = mf.matriz_covarianza(dip, dep, X_grid, Y_grid)
    # creacion de mapas de slip    
    Slip=mf.distribucion_slip(C, mu, 20)
    ventana = mf.ventana_taper_slip_fosa(Slip, lons,lats,2) # ventana de taper
    Slip    = mf.taper_slip_fosa(Slip,ventana)
    Slip    = mf.escalar_magnitud_momento(Mw+0.1, Slip, dep, X_grid, Y_grid) # se escala el Slip a la magnitud deseada <--------- Slip final
    # descomentar solo en fase de pruebas
    fig = plt.figure()
    ax = fig.add_subplot(111)
    im = ax.imshow(Slip)
    fig.colorbar(im)
    plt.show()

    fig = pygmt.Figure()
    # Chilean trench
    # Load sample grid (3 arc-minutes global relief) in target area
    grid = pygmt.datasets.load_earth_relief(resolution="30s", region=region)
    data_grid=np.genfromtxt('../Output_data/'+timestr+'_'+str(northlat)+str(southlat)+'.xyz',delimiter=',')
    # Plot original grid



    if plot_grid==True:
        fig.basemap(region=region, projection="M0/0/12c", frame=True)
        fig.grdimage(grid=grid, cmap="oleron")
        fig.plot(
            data=data_trench,
            region=region,
            pen="0.2p",
            fill="white",
            style="f0.5i/0.1i+r+t+o1",
        )
        fig.plot(data=data_grid, style="p0.2c", pen="2p,red",fill="red")
        fig.colorbar(frame=["x+lTopography", "y+lm"])
        flag=input('Do you want save figure? [y/n]')
        if flag=='y':
            fig.savefig('./Output_images/'+timestr+'_'+str(northlat)+str(southlat)+'.png')
        fig.show()

    return


if __name__ == "__main__":
    desc = """
        Expected region to segmentation "lonmin" "lonmax" "latmin" "latmax" (-r), north corner fault "lat" (-nlat),
        ,south corner fault "lat" (slat), number of points in x(-nx), number of points in y (lat) (-ny), width of fault (-w)
        FORMAT:  LAT from -90 to 90
                 LON from -180 to 180
        Writes file to ${-r}.xyz with 0 values in z

        Example 1 input:
            python trench_grid.py -r -77 -69 -37 -29 -nlat -29.5 -slat -36.5 -nx 10 -ny 30 -w 150

        Example 1 output (Output_data and Output_image folders respectly):
            -29.5-36.5.xyz
            -29.5-36.5.png
        """
    parser = argparse.ArgumentParser(
        description=desc, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "-n",
        "--nslip",
        dest="nslip",
        type=int,
        help="number of iterations for slip models.",
        required=True,
    )
    parser.add_argument(
        "-m",
        "--Mw",
        dest="mw",
        type=float,
        help="Moment magnitude of slip models.",
        required=True,
    )
    parser.add_argument(
        "-p",
        "--plotgrid",
        dest="plotgrid",
        type=bool,
        help="True if you want plot grid, False if not or omit this argument.",
        required=False,
    )
    parser.add_argument(
        "-r",
        "--region",
        metavar=("lonmin", "lonmax", "latmin", "latmax"),
        dest="region",
        type=float,
        nargs=4,
        help="limits of map [lonmin lonmax latmin latmax] ",
        required=True,
    )
    parser.add_argument(
        "-w",
        "--width",
        dest="width",
        type=float,
        required=True,
        help="width of fault (km)",
    )
    parser.add_argument(
        "-nx",
        "--numofx",
        dest="nx",
        type=int,
        required=True,
        help="number of x points",
    )
    parser.add_argument(
        "-ny",
        "--numberofy",
        dest="ny",
        type=int,
        required=True,
        help="number of y points",
    )
    parser.add_argument(
        "-nlat",
        "--northlat",
        dest="northlat",
        type=float,
        help="north lat of geometry ",
        required=True)
    parser.add_argument(
        "-slat",
        "--southlat",
        dest="southlat",
        type=float,
        help="south lat of geometry ",
        required=True)
    pargs = parser.parse_args()

    main(pargs)