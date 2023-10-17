import argparse
import numpy as np 
import pygmt
import time
# km2deg
def km2deg(km):
    return float(km/111.1)
# main
def main(args):
    """What is executed upon running the program. 

        Runs for make a .xyz grid with NX x NY points  
        """
    # set command line variables (after python makegrid.py)
    # -29.500513696035412, -72.28001907930886
    region=args.region
    strike=np.deg2rad(args.strike)
    width=args.width
    length=args.length
    corner=args.corner
    nx=args.nx
    ny=args.ny
    dx=width/nx
    dy=length/ny
    timestr = time.strftime("%Y%m%d-%H%M%S")
    # 
    np.disp('The south-west corner of your soon grid is: '+str(corner[0])+' '+str(corner[1]))
    north_west_lat=corner[0]
    north_west_lon=corner[1]
    points=nx*ny
    # make grid
    grid=open('./Output_data/'+timestr+'_'+str(corner[0])+str(corner[1])+'.xyz','w')
    for i in range(nx):
        for j in range(ny):
            norm=np.sqrt((dx*(2*i+1)/2)**2+(dy*(2*j+1)/2)**2)
            deg=np.arctan(((dx*(2*i+1)/2))/((dy*(2*j+1)/2)))+strike
            lat=float(north_west_lat)+(km2deg(norm*np.cos(deg)))
            lon=float(north_west_lon)+(km2deg(norm*np.sin(deg)))
            grid.write(str(lon)+'\t'+str(lat)+'\t'+str(strike)+'\n')
    grid.close()

    #  

    fig = pygmt.Figure()
    # Chilean trench
    data_trench=np.loadtxt("trench-chile.txt")
    # Load sample grid (3 arc-minutes global relief) in target area
    grid = pygmt.datasets.load_earth_relief(resolution="30s", region=region)
    data_grid=np.loadtxt('./Output_data/'+timestr+'_'+str(corner[0])+str(corner[1])+'.xyz')
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
    fig.savefig('./Output_images/'+timestr+'_'+str(corner[0])+str(corner[1])+'.png')
    fig.show()

    return


if __name__ == "__main__":
    desc = """
        Expected region to segmentation "lonmin" "lonmax" "latmin" "latmax" (-r), strike in degrees (-s), southwestern corner fault "lat" "lon" (-c),
        number of points in x(-nx), number of points in y (lat) (-ny), width of fault (-w) and length(-l)
        FORMAT:  LAT from -90 to 90
                 LON from -180 to 180
        Writes file to ${-r}.xyz with z column like strike values

        Example 1 input:
            python makegrid.py -r -74 -70 -37 -29 -c -33.0 -72.72589 -s 5 -nx 10 -ny 20 -w 200 -l 400

        Example 1 output (Output_data and Output_image folders respectly):
            -33.0-72.72589.xyz
            -33.0-72.72589.png

        Example 2 input:
            python makegrid.py -r -75 -69 -38 -29 -c -36.5 -74.237303 -s 23 -nx 10 -ny 20 -w 200 -l 430

        Example 2 output (Output_data and Output_image folders respectly):
            -36.5-74.237303.xyz
            -36.5-74.237303.png
        
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
        "-s",
        "--strike",
        dest="strike",
        type=float,
        required=True,
        help="strike of fault(degree respect North)",
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
        "-c",
        "--corner",
        metavar=("lat","lon"),
        dest="corner",
        type=float,
        nargs=2,
        help="southwestern corner fault ",
        required=True)
    pargs = parser.parse_args()

    main(pargs)