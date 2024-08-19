import numpy as np
import argparse
import os
import scipy.io
import sys
sys.path.append('../') 
import geostochpy
from geostochpy import modfilters
import json
from skimage.feature import peak_local_max
import cv2
def detect_main_patches(matrix, min_distance=5, threshold_abs=None):
    """
    Detecta parches principales en una matriz.

    Parameters:
    - matrix: Matriz de entrada.
    - min_distance: Distancia mínima entre los parches detectados.
    - threshold_abs: Valor absoluto para establecer un umbral de detección.

    Returns:
    - Coordenadas de los parches principales.
    """

    # Encuentra los máximos locales en la matriz
    coordinates = peak_local_max(matrix, min_distance=min_distance, threshold_abs=threshold_abs)

    return coordinates

def main(args):
    """What is executed upon running the program. 
    Filter slip program, according physics (coupling), depth, and
    by limits(limitations in borders) 
    """
    # Load args
    couplingfilename=args.couplingfile
    simulationfolder=args.simulationfolder
    # names of files
    files=os.listdir('../Output_data/'+simulationfolder)
    route_trench = "../Slab/trench-chile.txt"
    # load trench
    lons_fosa, lats_fosa  = geostochpy.load_trench(route_trench)
    # load files
    dicc={}
    max_slips=[]
    for filename in files:
        if filename.endswith('.mat'):
            fmat=scipy.io.loadmat('../Output_data/'+simulationfolder+'/'+filename)
            X_grid=fmat['lon']
            Y_grid=fmat['lat']
            Slip=fmat['slip']
            depth=fmat['depth']
            # filt by depth slip max
            np.append(max_slips,np.max(Slip))
            Slip_filter_physic_flag=modfilters.physical_filter(Slip,Y_grid,depth,10000,20000,-33.5,-31)
            if np.max(Slip)<10 or np.max(Slip) > 14 or detect_main_patches(Slip,min_distance=10,threshold_abs=0.8).size == 0 or np.min(Y_grid)>-33 or np.max(Y_grid)<-31:
                Slip_filter_physic_flag=False
            # filt for location of slip max. we need that its in 
            #
            dicc.update({filename:Slip_filter_physic_flag})
    with open('../Output_data/'+simulationfolder+'/Slip_filters.json', 'w') as archivo:
        json.dump(dicc, archivo)
    filtered_true=list(dict(filter(lambda item: item[1], dicc.items())).keys())
    print(len(filtered_true))
    print(np.nanmean(max_slips))
    print(filtered_true)
    # Now this list we filt by physical limitations
    # dicc_physics={}
    # for element in filtered_true:
    #     fmat=scipy.io.loadmat('../Output_data/'+simulationfolder+'/'+element)
    #     Slip=fmat['slip']
    #     Depth=fmat['depth']
    #     dicc_physics.update({element:Slip_filter_physic_flag})
    # we made a list with items that meet the conditions of depth and physical conditions
    # filtered_by_depth_and_physics=list(dict(filter(lambda item: item[1], dicc_physics.items())).keys())
    # print(filtered_by_depth_and_physics)

if __name__ == "__main__":
    desc = """
        Expected region to segmentation "lonmin" "lonmax" "latmin" "latmax" (-r), north corner fault "lat" (-nlat),
        ,south corner fault "lat" (slat), number of points in x(-nx), number of points in y (lat) (-ny), width of fault (-w),
        number of simulations (-n), magnitude(Mw, -m) and -p True or False if you want to show the plot
        
        FORMAT:  LAT from -90 to 90
                 LON from -180 to 180

        Example 1 input:
            python Filter_slips.py -c coupling.xyz

        """
    parser = argparse.ArgumentParser(
        description=desc, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "-c",
        "--couplingfile",
        dest="couplingfile",
        type=str,
        help="Name of coupling file, MUST BE IN auxiliar folder",
        required=False)
    parser.add_argument(
        "-f",
        "--simulationfolder",
        dest="simulationfolder",
        type=str,
        help="Name of simulation folder, MUST BE inside of Output_data folder",
        required=True)

    pargs = parser.parse_args()

    main(pargs)