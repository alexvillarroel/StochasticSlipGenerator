"""
Slip and deformation figure of a Stochastic generation
============================================================================================================
    This example makes a plot of Slip distribution and of deformation in meters.

"""
import scipy.io
import geostochpy
import matplotlib.pyplot as plt
import numpy as np
from geostochpy import modokada as mo
from clawpack.clawutil.data import get_remote_file
# shorelines
shorelines_file = geostochpy.get_data('pacific_shorelines_east_4min.npy')
shore = np.load(shorelines_file)
#
route_trench = geostochpy.get_data('trench-chile.txt')
# load trench
lonfosa, latfosa  = geostochpy.load_trench(route_trench)
files=['sim_0293.mat', 'sim_2053.mat', 'sim_3036.mat', 'sim_3076.mat', 'sim_3384.mat', 'sim_3677.mat', 'sim_4115.mat', 'sim_4128.mat', 'sim_4445.mat']
folder='../Output_data/Simulation_20240207-085309/'
# folder='../Output_data/Simulation_9_150_450/'
# files=['sim_0582.mat', 'sim_1595.mat', 'sim_1666.mat', 'sim_2004.mat', 'sim_2148.mat', 'sim_2877.mat', 'sim_2978.mat', 'sim_3528.mat', 'sim_3930.mat']
region=[-78,-68,-38,-28]

#
for file in files:
    filename=folder+'/img/'+file.replace(".mat", ".png")
    filename_def=folder+'/img/'+file.replace(".mat", "def_.png")
    filename_def2=folder+'/img/'+file.replace(".mat", "def_cp_.png")
    fmat=scipy.io.loadmat(folder+file)
    X_grid=fmat['lon']
    Y_grid=fmat['lat']
    Slip=fmat['slip']
    Strike=fmat['strike']
    Dip=fmat['dip']
    Rake=fmat['rake']
    depth=fmat['depth']
    # plt.hist(Slip.flatten())
    # plt.axvline(np.percentile(Slip,85))
    # plt.show()
    geostochpy.plot_slip(X_grid,Y_grid,lonfosa,latfosa,Slip,filename)
    dtopo = mo.okada_solucion_optimized( X_grid, Y_grid, 550/180, Strike, Dip, depth, Rake, Slip, 550000, resolucion = 1/30., tamano_buffer = 1., verbose = False ) # calculo deformacion
    deformation=dtopo.dZ_at_t(0)
    X_deformation=dtopo.X
    Y_deformation=dtopo.Y
    geostochpy.plot_deformation(X_deformation,Y_deformation,lonfosa,latfosa,deformation.reshape(X_deformation.shape),filename_def)
    dtopo.plot_dZ_colors(t=0,dZ_interval=1)
    plt.plot(shore[:,0]-360, shore[:,1], 'g')
    plt.plot(lonfosa,latfosa,'darkgreen')
    plt.axis(region)
    plt.grid(visible=True,axis='both')
    plt.plot(-71.63, -33.03,'o',color='gold')
    plt.text(-71.63, -33.03,'Valparaíso',ha='left')
    plt.savefig(filename_def2,dpi=1200)
    # dtopo.write(filename_def,dtopo_type=3)
    # geostochpy.plot_slip_gmt(region,X_grid,Y_grid,lonfosa,latfosa,Slip,dx,dy,filename)