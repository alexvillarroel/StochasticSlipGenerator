import geostochpy
import numpy as np
from scipy.io import loadmat
import os
file='/home/alex/StochasticSlipGenerator/Output_data/Simulation_9.3_5000_coupling/simulations.mat'

def make_faultfile_mat(matfile,n_sim,filename):
    NEWLINE_SIZE_IN_BYTES = 1  # -2 on Windows?
    DATA=loadmat(matfile,struct_as_record=False,squeeze_me=True)
    DATA=DATA['simulations'][n_sim-1]
    lat=DATA.lat.flatten()
    lon=DATA.lon.flatten()
    depth=DATA.depth.flatten()
    depth=depth/1000
    length=DATA.length.flatten()
    width=DATA.width.flatten()
    dip=DATA.dip.flatten()
    strike=DATA.strike.flatten()
    rake=DATA.rake.flatten()
    slip=DATA.slip.flatten()
    DATA = np.column_stack((lat, lon, depth, length, width, dip, strike, rake, slip))
    with open(filename, 'w') as file:
        file.write("! lat[deg] lon[deg] depth[km] length[km] width[km] dip[deg] strike[deg] rake [deg] slip_amp[m]\n")
        np.savetxt(file, DATA, fmt='%6.4f', delimiter=' ')
        file.seek(0, os.SEEK_END) # Go to the end of the file.
        # Go backwards one byte from the end of the file.
        file.seek(file.tell() - NEWLINE_SIZE_IN_BYTES, os.SEEK_SET)
        file.truncate() # Truncate the file to this point.
    return
make_faultfile_mat(file,50,'faults.txt')