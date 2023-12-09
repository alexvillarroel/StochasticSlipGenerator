import numpy as np
import sys 
# appending a path 
sys.path.append('../')
import main as slgen
from scipy.io import savemat,loadmat 
route='../Output_data/Simulation_8.8/'
file='sim_0178.mat'
#route trench
route_trench = "../Slab/trench-chile.txt"
fmat=loadmat(route+file)
X_grid=fmat['lon']
Y_grid=fmat['lat']
Slip=fmat['slip']
depth=fmat['depth']
region=[-77,-69,-37,-29]
lons_fosa, lats_fosa  = slgen.load_trench(route_trench)
slgen.plot_slip_gmt(region,X_grid,Y_grid,lons_fosa,lats_fosa,Slip,250/15,650/25,'Plot_filtered')