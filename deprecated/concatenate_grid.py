import argparse
import math
import os
import numpy as np 
import pygmt
def main(args):
    files=args.files
    output=args.output[0]
    region=args.region
    # concatenate files
    with open('./Output_data/'+str(output), "w") as new_file:
        for name in files:
            with open('./Output_data/'+name) as f:
                for line in f:
                    new_file.write(line)
                
                new_file.write("\n")
    #
    fig = pygmt.Figure()
    # Chilean trench
    data_trench=np.loadtxt("trench-chile.txt")
    data_grid=np.loadtxt('./Output_data/'+str(output))
    # Load sample grid (3 arc-minutes global relief) in target area
    grid = pygmt.datasets.load_earth_relief(resolution="30s", region=region)
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
    fig.colorbar(frame=["x+lTopography", "y+lm"])
    fig.savefig('grid.png')
    fig.show()
    return

if __name__ == "__main__":
    desc = """
        This program concatenates different grids and generates a graph with the concatenated grid.

        FORMAT:  
        Writes file to ${-r}.xyz with z column full of 0
        Expected region to segmentation "lonmin" "lonmax" "latmin" "latmax" (-r), -f for files to insert(only 2)
        and -o for output name file

        Example 1 input:
            python concatenate_grid.py -r -74 -70 -37 -29 -f 20230930-111539_-36.5-74.237303.xyz 20230930-111433_-33.0-72.72589.xyz -o concatenated_grid.xyz

        Example 1 output:
            concatenated_grid.xyz
            concatenated_grid.png

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
        "-o",
        "--outputfile",
        dest="output",
        type=str,
        nargs=1,
        help="output name for the new grid",
        required=True,
    )
    parser.add_argument(
        "-f",
        "--files",
        dest="files",
        metavar=("file1","file2"),
        type=str,
        nargs=2,
        required=True,
        help="name of files to concatenate")
    pargs = parser.parse_args()
    main(pargs)