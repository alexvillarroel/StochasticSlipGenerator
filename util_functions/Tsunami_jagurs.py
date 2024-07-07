from geostochpy import modjagurs
from scipy.io import loadmat
import numpy as np
import os
folder='/home/alex/StochasticSlipGenerator/Output_data/Simulation_9.0_1000_coupling/'
filename='simulations.mat'
index_fault=0

def make_faultfile(matfile,index_fault,filename_output):
    NEWLINE_SIZE_IN_BYTES = 1  # -2 on Windows?
    simulation=loadmat(folder+filename,struct_as_record=False, squeeze_me=True)
    simulation=simulation['simulations']
    lat = np.array(simulation[index_fault-1].lon).flatten()
    lon = np.array(simulation[index_fault-1].lat).flatten()
    slip = np.array(simulation[index_fault-1].slip).flatten()
    depth = np.array(simulation[index_fault-1].depth).flatten()
    length = np.array(simulation[index_fault-1].length).flatten()
    width = np.array(simulation[index_fault-1].width).flatten()
    dip = np.array(simulation[index_fault-1].dip).flatten()
    strike = np.array(simulation[index_fault-1].strike).flatten()
    rake = np.array(simulation[index_fault-1].rake).flatten()
    DATA = np.column_stack((lat, lon, depth, length, width, dip, strike, rake, slip))
    with open(filename_output, 'w') as file:
        file.write("! lat[deg] lon[deg] depth[km] length[km] width[km] dip[deg] strike[deg] rake [deg] slip_amp[m]\n")
        np.savetxt(file, DATA, fmt='%6.4f', delimiter=' ')
        file.seek(0, os.SEEK_END)
        file.seek(file.tell() - NEWLINE_SIZE_IN_BYTES, os.SEEK_SET)

make_faultfile(folder+filename,1,'faults.txt')