import argparse
import numpy as np
import sys 
# appending a path 
sys.path.append('../') 
import matplotlib.cm as cm
import matplotlib as mpl
import random
import matplotlib.pyplot as plt
from scipy.signal.windows import general_hamming
from scipy.io import savemat 
import geostochpy
import os
import time
from tqdm import trange
from clawpack.geoclaw import dtopotools, topotools
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
    route_trench = geostochpy.get_data('trench-chile.txt')
    lonsfosa, latsfosa,strikefosa  = geostochpy.load_trench(route_trench)
    slabdep,slabdip,slabstrike,slabrake=geostochpy.load_files_slab2(zone='south_america',rake=True)
    dir='Simulation_'+timestr
    os.chdir('../Output_data/')
    os.mkdir(dir)
    os.chdir(dir)
    os.mkdir('img')
    #print('Do you wanna make multifault.ctl file?')
    #answer=input('[yes/no]')
    answer="no"
    if answer=='yes':
        lat_griddomain=input('Insert latitude of Domain 01 layer')
        lon_griddomain=input('Insert longitude of Domain 01 layer')
    #
    # make grid
    for i in trange(1,n_slip+1,desc=f'Generation process: '):
        random_north = northlat - random.random() * (np.abs(northlat-southlat)-geostochpy.km2deg(length))
        lon,lat,lon_flat,lat_flat=geostochpy.make_fault_alongstriketrench(lonsfosa, latsfosa,strikefosa,random_north, nx, ny, width, length)
        X_grid,Y_grid,dep,dip,strike,rake=geostochpy.interp_slabtofault(lon_flat,lat_flat,nx,ny,slabdep,slabdip,slabstrike,slabrake)
        ## Creation slip models
        # mean matrix
        media,rigidez=geostochpy.media_slip(Mw,dx*1000,dy*1000,dep)
        leveque_taper=geostochpy.taper_LeVeque(dep,55000)
        # leveque_taper=leveque_taper/np.max(leveque_taper)
        villarroel_taper=geostochpy.taper_except_trench_tukey(dep,alpha_dip=0.3,alpha_strike=0.3)
        taper=leveque_taper*villarroel_taper
        # taper=geostochpy.taper_except_trench_tukey(dep,alpha_dip=0.6,alpha_strike=0.4,dip_taperfunc=geostochpy.taper_LeVeque,strike_taperfunc=geostochpy.tukey_window_equal)

        mu = geostochpy.matriz_medias_villarroel(media,taper)
        # matriz de covarianza
        C    = geostochpy.matriz_covarianza_optimized(dip, dep, X_grid, Y_grid,length*1000,width*1000)
        # for comcot simulation
        Slip=geostochpy.distribucion_slip(C, mu, 20)
        Slip,rigidez,Mo_original,Mo_deseado=geostochpy.escalar_magnitud_momento(Mw, Slip, dep, dy*1000, dx*1000,prem=True) # se escala el Slip a la magnitud deseada <--------- Slip final
        # Hypocenter=geostochpy.hypocenter(X_grid,Y_grid,dep,length,width) se tiene en cuenta la rigidez con el modelo PREM incluido @fetched with Rockhound
        archivo_salida='sim_'+str(i).zfill(len(str(n_slip)))
        # file_multifault='fault_multi_'+str(i).zfill(len(str(n_slip)))+'.ctl'
        mdic_multifault={'depth':dep,'length':dy*np.ones((1,nx*ny)),'width':dx*np.ones((1,nx*ny)),
              'slip':Slip,'strike':strike,'dip':dip,'rake':rake
              ,'lat':Y_grid,'lon':X_grid,'time':np.zeros((1,nx*ny))}
        # mdic={'depth':dep,'length':dy*np.ones((1,nx*ny)),'width':dx*np.ones((1,nx*ny)),
        #       'slip':Slip,'strike':strike,'dip':dip,'rake':rake
        #       ,'lat':Y_grid,'lon':X_grid,'time':np.zeros((1,nx*ny))}
        savemat(archivo_salida+'.mat',mdic_multifault) # se guarda los arrays de la falla
        # for make deformation file
        # header='Longitude,Latitude,Depth(km),Length,Width,Strike,Dip,Rake,Slip,UnitSrc'
        # array = np.column_stack((lons_ep.flat, lats_ep.flat, dep.flat, (dy*np.ones((1,nx*ny))).flat, (dx*np.ones((1,nx*ny))).flat, strike.flat, dip.flat, rake.flat, slip, np.ones((1,nx*ny)).flat))
        # np.savetxt(archivo_salida+'.csv', array, delimiter=",", header=header, comments='')
        # subfault_fname = archivo_salida + '.csv'
        # input_units = {"length":"km", "width":"km", "depth":"m", "slip":"m"}
        # fault = dtopotools.CSVFault()
        # fault.read(subfault_fname, input_units=input_units, coordinate_specification="noaa sift")
        # resolucion=1/30.
        # tamano_buffer=2.
        # x, y = fault.create_dtopo_xy(dx=resolucion, buffer_size=tamano_buffer)
        # dtopo = fault.create_dtopography(x, y, times=[1.], verbose=True)
        # dtopofile=archivo_salida+'.tt3'
        # dtopo.write(dtopofile, dtopo_type=3)
        #

        # filename='./img/'+archivo_salida+'.png'
        # if answer=='yes':
        #     print('ok')
            # multifault.make_multifault(archivo_salida+'.mat',file_multifault,float(lat_griddomain),float(lon_griddomain))
        # geostochpy.plot_slip_gmt([-78,-68,-38,-28],X_grid,Y_grid,lons_fosa,lats_fosa,Slip,10,10,filename)

        #geostochpy.plot_slip(X_grid,Y_grid,lons_fosa,lats_fosa,Slip,filename)
        # geostochpy.plot_slip_css(region,lons,lats,lons_fosa,lats_fosa,Slip)
        #geostochpy.plot_slip_gmt(region,X_grid,Y_grid,lons_fosa,lats_fosa,Slip,dx,dy,filename)
    return

if __name__ == "__main__":
    desc = """
        Expected region to segmentation "lonmin" "lonmax" "latmin" "latmax" (-r), north corner fault "lat" (-nlat),
        ,south corner fault "lat" (slat), number of points in x(-nx), number of points in y (lat) (-ny), width of fault (-w),
        number of simulations (-n), magnitude(Mw, -m) and -p True or False if you want to show the plot
        
        FORMAT:  LAT from -90 to 90
                 LON from -180 to 180

        Example 1 input:
            python Slip_generator.py -r -77 -69 -37 -29 -nlat -28 -slat -36 -nx 15 -ny 45 -w 150 -l 450 -n 10 -m 9.0

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