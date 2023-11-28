import argparse
import numpy as np
import sys 
# appending a path 
sys.path.append('../') 
import matplotlib.cm as cm
import matplotlib as mpl
import random
import matplotlib.pyplot as plt
from scipy.io import savemat 
import main as slgen
import os
import time
from tqdm import trange
# main
def main(args):
    """What is executed upon running the program. 

        Runs for make a .xyz grid with NX x NY points  
        """
    # set command line variables (after python makegrid.py)
    # -29.500513696035412, -72.28001907930886
    region=args.region
    Mw=args.mw
    width=args.width
    northlat=args.northlat
    okada='No'
    southlat=args.southlat
    nx=args.nx
    length=args.length
    ny=args.ny
    dx=width/nx
    n_slip=args.nslip
    dy=length/ny
    timestr = time.strftime("%Y%m%d-%H%M%S")
    # 
    print('Starting Grid generation')
    route_trench = "../Slab/trench-chile.txt"
    lons_fosa, lats_fosa  = slgen.load_trench(route_trench)
    # load slab files
    slabdep,slabdip,slabstrike,slabrake=slgen.load_files_slab2(rake=True)
    dir='Simulation_'+timestr
    os.chdir('../Output_data/')
    os.mkdir(dir)
    os.chdir(dir)
    print('Do you wanna make multifault.ctl file?')
    #answer=input('[yes/no]')
    answer="no"
    if answer=='yes':
        lat_griddomain=input('Insert latitude of Domain 01 layer')
        lon_griddomain=input('Insert longitude of Domain 01 layer')
    #
    # make grid
    for i in trange(1,n_slip+1,desc=f'Generation process: '):
        random_north = northlat - random.random() * (np.abs(northlat-southlat)-slgen.km2deg(length))
        lons,lons_ep,lats,lats_ep=slgen.make_fault_alongtrench(lons_fosa,lats_fosa,random_north, nx,ny,width,length)
        # Interpolate Slab data to the new gridded fault
        [X_grid,Y_grid,dep,dip,strike,rake]=slgen.interp_slabtofault(lons,lats,nx,ny,slabdep,slabdip,slabstrike,slabrake)
        ## Creation slip models
        # mean matrix
        mu   = slgen.matriz_medias(10, dep)
        # matriz de covarianza
        C    = slgen.matriz_covarianza(dip, dep, X_grid, Y_grid)

        # for comcot simulation
        Slip=slgen.distribucion_slip(C, mu, 10)
        ventana = slgen.ventana_taper_slip_fosa(Slip, X_grid,Y_grid,2) # ventana de taper
        Slip    = slgen.taper_slip_fosa(Slip,ventana)
        Slip    = slgen.escalar_magnitud_momento(Mw, Slip, dep, X_grid, Y_grid) # se escala el Slip a la magnitud deseada <--------- Slip final
        # Hypocenter=slgen.hypocenter(X_grid,Y_grid,dep,length,width)
        slip=Slip.flat
        archivo_salida='sim_'+str(i).zfill(len(str(n_slip)))
        file_multifault='fault_multi_'+str(i).zfill(len(str(n_slip)))+'.ctl'
        mdic_multifault={'depth':dep.flat,'length':dy*np.ones((1,nx*ny)),'width':dx*np.ones((1,nx*ny)),
              'slip':slip,'strike':strike.flat,'dip':dip.flat,'rake':rake.flat
              ,'lat':lats_ep,'lon':lons_ep,'time':np.zeros((1,nx*ny))}
        mdic={'depth':dep,'length':dy*np.ones((1,nx*ny)),'width':dx*np.ones((1,nx*ny)),
              'slip':Slip,'strike':strike,'dip':dip,'rake':rake
              ,'lat':Y_grid,'lon':X_grid,'time':np.zeros((1,nx*ny))}
        savemat(archivo_salida+'.mat',mdic) # se guarda los arrays de la falla

        if answer=='yes':
            print('ok')
            # multifault.make_multifault(archivo_salida+'.mat',file_multifault,float(lat_griddomain),float(lon_griddomain))
        slgen.plot_slip(X_grid,Y_grid,lons_fosa,lats_fosa,Slip)
        #slgen.plot_slip_gmt(region,X_grid,Y_grid,lons_fosa,lats_fosa,Slip,archivo_salida,dx,dy)

    return

if __name__ == "__main__":
    desc = """
        Expected region to segmentation "lonmin" "lonmax" "latmin" "latmax" (-r), north corner fault "lat" (-nlat),
        ,south corner fault "lat" (slat), number of points in x(-nx), number of points in y (lat) (-ny), width of fault (-w),
        number of simulations (-n), magnitude(Mw, -m) and -p True or False if you want to show the plot
        
        FORMAT:  LAT from -90 to 90
                 LON from -180 to 180

        Example 1 input:
            python Slip_generator.py -r -77 -69 -37 -29 -nlat -29.5 -slat -37.5 -nx 12 -ny 48 -w 150 -l 450 -n 5 -m 9.0

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
        "-l",
        "--length",
        dest="length",
        type=float,
        required=True,
        help="length of fault (km)",
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
        help="north lat of possible segment",
        required=True)
    parser.add_argument(
        "-slat",
        "--southlat",
        dest="southlat",
        type=float,
        help="south lat of possible segment ",
        required=True)
    pargs = parser.parse_args()

    main(pargs)