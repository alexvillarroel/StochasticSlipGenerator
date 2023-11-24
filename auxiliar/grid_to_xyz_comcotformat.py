import pygmt
import pandas
import numpy as np

def grid_to_xyz(file,outputfile):
    """What is executed upon running the program. 

        Run for make a .xyz with comcot format with a .grd file.  
        """
    if '.grd' or '.nc' in file:
        xyz=pygmt.grd2xyz(file,output_type='pandas')
        xyz=xyz.sort_values(['y','x'],ascending=True)
        # bathymetry in positive values, topography in negative
        xyz['z'][~np.isnan(xyz['z'])]*=-1
        xyz.to_csv(outputfile,index=False,na_rep='nan',header=False)
    return
# use the function
#example
grid_to_xyz('centralchile.grd','centralchile.xyz')