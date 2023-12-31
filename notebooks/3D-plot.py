"""
    Visualize a 3D grid of points using the 3D grid .

    :return: [description]
    :rtype: [type]
    """
import scipy.io
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import cm
from matplotlib.ticker import LinearLocator
from matplotlib.colors import ListedColormap, Normalize

file='../Output_data/Simulation_9.0/sim_0004.mat'


filemat=scipy.io.loadmat(file)
Slip=filemat['slip']
X_array=filemat['lon']
Y_array=filemat['lat']
depth=filemat['depth']/-1000
region=[-76,-70,-36,-28]

def plot_3d(X_array,Y_array,depth,Slip,filename=None):
    """
    Plot a 3D plot of Slip with a depth of the surface .

    :param X_array: [description]
    :type X_array: [type]
    :param Y_array: [description]
    :type Y_array: [type]
    :param depth: [description]
    :type depth: [type]
    :param Slip: [description]
    :type Slip: [type]
    :param filename: [description], defaults to None
    :type filename: [type], optional
    :return: [description]
    :rtype: [type]
    """
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    # Plot the surface with face colors taken from the array we made.
    surf = ax.plot_surface(X_array,Y_array, depth, facecolors=cm.rainbow(Slip/np.max(Slip.flatten())), linewidth=0)
    # Customize the z axis.
    # Agregar barra de colores con el rango correcto
    norm = plt.Normalize(np.min(Slip.flat), np.max(Slip.flat))
    sm = plt.cm.ScalarMappable(cmap=plt.cm.rainbow, norm=norm)
    sm.set_array([])  # este paso es necesario para que la barra de colores refleje correctamente los valores
    plt.title('Slip_0004_Simulation')
    # Configurar etiquetas
    ax.set_xlabel('Longitude')
    ax.set_ylabel('Latitude')
    ax.set_zlabel('Depth')
    cbar = fig.colorbar(sm, ax=ax, pad=0.1)
    cbar.set_label('Slip[m]')
    # Configurar posición de la barra de colores
    cbar.ax.yaxis.set_label_position('right')  # Puedes ajustar la posición según tus necesidades
    if filename!=None:
        fig.savefig(filename)
    return plt.show()


plot_3d(X_array,Y_array,depth,Slip)
