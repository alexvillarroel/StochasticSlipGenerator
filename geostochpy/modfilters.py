import scipy.interpolate
import scipy.spatial.distance
import numpy as np
import random
import geostochpy as slgen
def obtain_border(matrix):
    border = []
    
    # Obtener la primera fila
    border.extend(matrix[0])

    # Obtener la última columna (sin incluir la primera fila y la última fila)
    border.extend(fila[-1] for fila in matrix[1:-1])

    # Obtener la última fila en orden inverso
    border.extend(reversed(matrix[-1]))

    # Obtener la primera columna (sin incluir la primera y la última fila)
    border.extend(fila[0] for fila in reversed(matrix[1:-1]))

    return border
def corr2_coeff(A, B):
    # Rowwise mean of input arrays & subtract from input arrays themeselves
    A_mA = A - A.mean(1)[:, None]
    B_mB = B - B.mean(1)[:, None]

    # Sum of squares across rows
    ssA = (A_mA**2).sum(1)
    ssB = (B_mB**2).sum(1)

    # Finally get corr coeff
    return np.dot(A_mA, B_mB.T) / np.sqrt(np.dot(ssA[:, None],ssB[None]))

def depthfilter(Slip,depth,prctile=10):
    flag=True
    q=np.percentile(Slip.flat,prctile)
    idx=np.argwhere(depth>=55000)
    if len(idx)>0:
        for i,j in idx:
            if Slip[i,j]>=q:
                flag=False
                return flag
    return flag

def couplingfilter(X_grid,Y_grid,Slip,couplingfilename,lonfosa,latfosa):
    coupling=np.genfromtxt('../auxiliar/'+couplingfilename,delimiter='\t')
    # erasing nans for interpolation
    mask = np.isnan(coupling).any(axis=1)
    coupling_sorted=coupling[~mask]
    #
    [X,Y]=np.meshgrid(np.unique(coupling_sorted[:,0]),np.unique(coupling_sorted[:,1]))
    interp_coupling=scipy.interpolate.LinearNDInterpolator(list(zip(coupling_sorted[:,0],coupling_sorted[:,1])),coupling_sorted[:,-1],fill_value=0)
    Coupling_interpolated=interp_coupling(X_grid,Y_grid)
    #slgen.plot_slip(X_grid,Y_grid,lonfosa,latfosa,corr2D,cmap="hot_r")
    return Slip

def physical_filter(Slip,lats,depth,profmin,profmax,latmin,latmax):
    flag=True
    idx=np.argmax(Slip)
    lat=lats.flatten()[idx]
    prof=depth.flatten()[idx]
    if (profmin<=prof<=profmax) and (latmin<=lat<=latmax):
        return flag
    else:
        return False

def depth_max_slip(Slip,depth,prof_range):
    flag=True

    return flag

