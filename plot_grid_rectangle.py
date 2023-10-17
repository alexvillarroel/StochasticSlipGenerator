import argparse
import numpy as np 
import pygmt
import time
import os
import scipy
# km2deg
def km2deg(km):
    return float(km/111.1)
def deg2km(deg):
    return float(deg*111.1)
# main
def main(args):
    """What is executed upon running the program. 

        Runs for make a .xyz grid with NX x NY points  
        """
    # set command line variables (after python makegrid.py)
    # -29.500513696035412, -72.28001907930886
    region=args.region
    width=args.width
    northlat=args.northlat
    southlat=args.southlat
    nx=args.nx
    ny=args.ny
    dx=width/nx
    dy=deg2km((northlat-southlat))/ny
    timestr = time.strftime("%Y%m%d-%H%M%S")
    # 
    np.disp('The width of your soon grid is: '+str(deg2km((northlat-southlat)))+' km ')
    points=nx*ny
    data_trench=np.loadtxt("trench-chile.txt")
    # make grid
    grid=open('temp.xyz','w')

    for i in range(nx):
        for j in range(ny):
            lat=northlat-km2deg(dy*j)
            lon=np.interp(lat,data_trench[:,1],data_trench[:,0])+km2deg(dx*i)
            grid.write(str(lon)+'\t'+str(lat)+'\t'+'0'+'\n')
    grid.close()

    #  

    fig = pygmt.Figure()
    # Chilean trench
    # Load sample grid (3 arc-minutes global relief) in target area
    grid = pygmt.datasets.load_earth_relief(resolution="30s", region=region)
    data_grid=np.loadtxt('temp.xyz')
    os.remove('temp.xyz')
    data_sort_lat=np.unique(np.sort(data_grid,axis=1)[:,1])
    print(data_sort_lat.shape)
    data_sort_lon=np.unique(np.sort(data_grid,axis=0)[:,0])
    print(data_sort_lon.shape)
    # Plot original grid
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
            python trench_grid.py -r -77 -69 -37 -29 -nlat -29.5 -s lat -36.5 -nx 10 -ny 30 -w 150

        Example 1 output (Output_data and Output_image folders respectly):
            -29.5-36.5.xyz
            -29.5-36.5.png
        """
    parser = argparse.ArgumentParser(
        description=desc, formatter_class=argparse.RawDescriptionHelpFormatter
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