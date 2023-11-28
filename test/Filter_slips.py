import numpy as np
import argparse
import scipy.io
import sys
sys.path.append('../') 
import main as slgen
from main import modfilters

def main(args):
    """What is executed upon running the program. 
    Filter slip program, according physics (coupling), depth, and
    by limits(limitations in borders) 
    """
    # Load args
    couplingfilename=args.couplingfile
    simulationfolder=args.simulationfolder
    # names of files
    file='sim_1.mat'
    route_trench = "../Slab/trench-chile.txt"
    # load trench
    lons_fosa, lats_fosa  = slgen.load_trench(route_trench)
    # load files
    fmat=scipy.io.loadmat('../Output_data/'+simulationfolder+'/'+file)
    X_grid=fmat['lon']
    Y_grid=fmat['lat']
    Slip=fmat['slip']
    modfilters.couplingfilter(X_grid,Y_grid,Slip,couplingfilename,lons_fosa,lats_fosa)
    return

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