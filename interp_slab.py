import argparse
import math
import os
import numpy as np 
import pygmt
from shapely.geometry import Point, MultiPoint
from shapely.ops import nearest_points
def main(args):
    file=args.file[0]
    #output=args.output
    date=args.date_slab[0]
    suffix=args.suffix_nameslab[0]
    outgrid=args.outgrid[0]
    dir='./Output_data/'
    dir_slab='./Slab2/'
    list=['_unc_','_thk_','_str_','_dip_','_dep_']
    list2 = [suffix+cadena+date+'.xyz' for cadena in list]
    mesh=np.genfromtxt(dir+file,delimiter='\t')
    mesh[:,0]=mesh[:,0]+360
    os.chdir('./Slab2')
    #print(list)
    flag=True
    for index,file2 in enumerate(filter(os.path.isfile,os.listdir())):
        if file2 in list2:
            temp=np.genfromtxt(file2,delimiter=',')
            if flag==True:
                destinations=[]
                for j in range(np.shape(temp)[0]):
                    destinations.append(Point(temp[j,0],temp[j,1]))
                destinations2=MultiPoint(destinations)
                flag=False
                #
            filename=open(file2.replace(suffix,outgrid),'w+')
            for i in range(np.shape(mesh)[0]):
                orig=Point(mesh[i,0],mesh[i,1])
                nearest_point=nearest_points(orig,destinations2)
                idx=destinations.index(nearest_point[1])
                filename.write(str(mesh[i,0])+'\t'+str(mesh[i,1])+'\t'+str(temp[idx,2])+'\n')
            filename.close()
    return

if __name__ == "__main__":
    desc = """
        This program interpolates and make files for uncertainties,
        thickness, strike, dip and depth , with an mesh with .xyz extension like input and 
        the geometry model for a region with Slab2.0 
        (https://www.sciencebase.gov/catalog/item/5aa1b00ee4b0b1c392e86467).

        FORMAT:  
        Writes file to ${-r}.xyz with z column full of 0
        Expected region to segmentation "lonmin" "lonmax" "latmin" "latmax" (-r), -f for files to insert(only 2)
        and -o for output name file
        Write files with output name inserted (-o) and the extension _unc, _thk, _str, _dip, _dep

        Example 1 input:
            python interp_slab.py -s sam_slab2 -d 02.23.18 -f concatenated_grid.xyz

        Example 1 output:
            concatenated_grid.xyz
            concatenated_grid.png

        """
    parser = argparse.ArgumentParser(
        description=desc, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "-o",
        "--outgrid",
        dest="outgrid",
        type=str,
        nargs=1,
        required=True,
        help="output suffix for the grids")
    parser.add_argument(
        "-f",
        "--file",
        dest="file",
        type=str,
        nargs=1,
        required=True,
        help="file to interpolate his geometry from Slab2.0")
    parser.add_argument(
        "-s",
        "--suffix_slabfiles",
        dest="suffix_nameslab",
        type=str,
        nargs=1,
        required=True,
        help="suffix of Slab2.0 filenames (ex. From  sam_slab2_unc_02.23.18.xyz, -n sam_slab2)")
    parser.add_argument(
        "-d",
        "--date_slabfiles",
        dest="date_slab",
        type=str,
        nargs=1,
        required=True,
        help="suffix date of Slab2.0 files(ex. From  sam_slab2_unc_02.23.18.xyz, -d 02.23.18)")
    pargs = parser.parse_args()
    main(pargs)