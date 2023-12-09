import numpy as np
import argparse
import os
import scipy.io
import sys
sys.path.append('../') 
import main as slgen
from main import modfilters
import json
import shutil

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
    lons_fosa, lats_fosa  = slgen.load_trench(route_trench)
    # load files
    dicc={}
    for filename in files:
        if filename.endswith('.mat'):
            fmat=scipy.io.loadmat('../Output_data/'+simulationfolder+'/'+filename)
            X_grid=fmat['lon']
            Y_grid=fmat['lat']
            Slip=fmat['slip']
            depth=fmat['depth']
            Slip_filter_flag=modfilters.depthfilter(X_grid,Y_grid,Slip,depth)
            dicc.update({filename:Slip_filter_flag})
    with open('../Output_data/'+simulationfolder+'/Slip_filters.json', 'w') as archivo:
        json.dump(dicc, archivo)
    filtered_true=list(dict(filter(lambda item: item[1], dicc.items())).keys())
    # Now this list we filt by physical limitations
    dicc_physics={}
    for element in filtered_true:
        fmat=scipy.io.loadmat('../Output_data/'+simulationfolder+'/'+element)
        Slip=fmat['slip']
        Slip_filter_physic_flag=modfilters.physical_filter(Slip)
        dicc_physics.update({element:Slip_filter_physic_flag})
    # we made a list with items that meet the conditions of depth and physical conditions
    filtered_by_depth_and_physics=list(dict(filter(lambda item: item[1], dicc_physics.items())).keys())
    print(len(filtered_by_depth_and_physics))
    print(filtered_by_depth_and_physics)

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