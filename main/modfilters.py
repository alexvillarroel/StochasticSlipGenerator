import scipy.interpolate
import scipy.spatial.distance
import numpy as np
import random
import main as slgen
def corr2_coeff(A, B):
    # Rowwise mean of input arrays & subtract from input arrays themeselves
    A_mA = A - A.mean(1)[:, None]
    B_mB = B - B.mean(1)[:, None]

    # Sum of squares across rows
    ssA = (A_mA**2).sum(1)
    ssB = (B_mB**2).sum(1)

    # Finally get corr coeff
    return np.dot(A_mA, B_mB.T) / np.sqrt(np.dot(ssA[:, None],ssB[None]))
def depthfilter(X_grid,Y_grid,Slip,depth):
    
    return Slip
def couplingfilter(X_grid,Y_grid,Slip,couplingfilename,lonfosa,latfosa):
    coupling=np.genfromtxt('../auxiliar/'+couplingfilename,delimiter='\t')
    # erasing nans for interpolation
    mask = np.isnan(coupling).any(axis=1)
    coupling_sorted=coupling[~mask]
    #
    [X,Y]=np.meshgrid(np.unique(coupling_sorted[:,0]),np.unique(coupling_sorted[:,1]))
    interp_coupling=scipy.interpolate.LinearNDInterpolator(list(zip(coupling_sorted[:,0],coupling_sorted[:,1])),coupling_sorted[:,-1],fill_value=0)
    Coupling_interpolated=interp_coupling(X_grid,Y_grid)
    dist=scipy.spatial.distance.euclidean(Slip.flat/np.max(Slip.flat),Coupling_interpolated.flat)
    print(dist)
    #slgen.plot_slip(X_grid,Y_grid,lonfosa,latfosa,corr2D,cmap="hot_r")
    return Slip