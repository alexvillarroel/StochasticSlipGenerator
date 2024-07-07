import argparse
import numpy as np
import sys
import random
from scipy.io import savemat
from scipy.interpolate import griddata
import geostochpy
import os
import time
from tqdm import trange
from joblib import Parallel, delayed  # Para paralelizar



def generate_grid(i, n_slip,n_subfaults, northlat, southlat, Mw, lonsfosa, latsfosa, strikefosa, slabdep, slabdip, slabstrike, slabrake,slabcoupling, scenarios):
    # Seleccionar aleatoriamente un escenario (L, W) de la lista de escenarios
    length, width = random.choice(scenarios)
    length = int(np.round(length))
    width = int(np.round(width))

    # Encontrar los mejores valores de nx, dx, ny y dy
    nx, ny,dx,dy = find_best_factors(length, width,n_subfaults)
    # Calcular dx y dy basados en length y width

    # Estimar el norte aleatoriamente dentro del rango permitido
    random_north = northlat - random.random() * (np.abs(northlat - southlat) - geostochpy.km2deg(length))
    lon, lat, lon_flat, lat_flat = geostochpy.make_fault_alongstriketrench(lonsfosa, latsfosa, strikefosa, random_north, nx, ny, width, length)
    X_grid, Y_grid, dep, dip, strike, rake = geostochpy.interp_slabtofault(lon_flat, lat_flat, nx, ny, slabdep, slabdip, slabstrike, slabrake)
    
    # Creación de modelos de deslizamiento
    media, rigidez = geostochpy.media_slip(Mw, dx * 1000, dy * 1000, dep)
    taper_coupling= griddata((slabcoupling[:,0], slabcoupling[:,1]), slabcoupling[:,2], (X_grid, Y_grid), method='linear', fill_value=0, rescale=True)
    mu = geostochpy.matriz_medias_villarroel(media, taper_coupling)
    C = geostochpy.matriz_covarianza_von_karman(dip, dep, X_grid, Y_grid, length, width)
    Slip = geostochpy.distribucion_slip_optimizada(C, mu, nx*ny-1)
    Slip, rigidez, Mo_original, Mo_deseado = geostochpy.escalar_magnitud_momento(Mw, Slip, dep, dy * 1000, dx * 1000, prem=True)
    
    return {
        'depth': dep,
        'length': dy * np.ones((ny, nx)),
        'width': dx * np.ones((ny, nx)),
        'slip': Slip,
        'strike': strike,
        'dip': dip,
        'rake': rake,
        'lat': Y_grid,
        'lon': X_grid,
        'time': np.zeros((ny, nx))
    }

def main(args):
    """Se ejecuta al correr el programa para generar una cuadrícula .xyz con NX x NY puntos"""
    region = args.region
    Mw = args.mw
    northlat = args.northlat
    southlat = args.southlat
    n_slip = args.nslip
    n_subfaults=args.ns
    
    # Generar escenarios de longitud y ancho para la magnitud dada
    scenarios = generate_scenarios(Mw, n_slip)

    print('Starting Grid generation')
    route_trench = geostochpy.get_data('trench-chile.txt')
    lonsfosa, latsfosa, strikefosa = geostochpy.load_trench(route_trench)
    slabdep, slabdip, slabstrike, slabrake = geostochpy.load_files_slab2(zone='south_america', rake=True)
    route_file=geostochpy.get_data('median_lock_Herrera2023.txt')
    route_file_mesh=geostochpy.get_data('mesh__Herrera2023.npy')
    mesh=np.load(route_file_mesh,allow_pickle=True)
    median_lock=np.loadtxt(route_file)
    mesh1=mesh[0]
    x=mesh1[:,0]
    y=mesh1[:,1]
    z=median_lock
    slabcoupling=np.column_stack((x, y, z))
    dir = f'Simulation_{Mw}_{n_slip}_coupling'
    os.makedirs(os.path.join('../Output_data/', dir, 'img'), exist_ok=True)
    
    # Paralelización del proceso de generación de modelos de deslizamiento
    results = Parallel(n_jobs=-1)(
        delayed(generate_grid)(i, n_slip,n_subfaults, northlat, southlat, Mw, lonsfosa, latsfosa, strikefosa, slabdep, slabdip, slabstrike, slabrake,slabcoupling, scenarios)
        for i in trange(1, n_slip + 1, desc='Generation process')
    )

    # Guardar todos los resultados en un único archivo .mat
    savemat(os.path.join('../Output_data/', dir, 'simulations.mat'), {'simulations': results})

def generate_scenarios(M, n_slips, sigma_L=0.18**2, alpha_W=0.17**2):
    """
    Genera escenarios de longitud (L) y ancho (W) para un evento sísmico dado.

    Parámetros:
    M (float): Magnitud del evento sísmico.
    num_scenarios (int): Número de escenarios a generar.
    sigma_L (float): Desviación estándar para log10(L).
    alpha_W (float): Desviación estándar para log10(W).
    Paper:"Scaling Relations of Earthquake Source Parameter Estimates
        with Special Focus on Subduction Environment
        by Lilian Blaser, Frank Krüger, Matthias Ohrnberger, and Frank Scherbaum"
    Retorna:
    list: Lista de tuplas (L, W) con los escenarios generados.
    """
    # Ecuaciones para log10(L) y log10(W)
    mu_L = -2.37 + 0.57 * M
    mu_W = -1.86 + 0.46 * M

    # Generar valores de log10(L) y log10(W)
    log10_L = np.random.normal(mu_L, sigma_L, n_slips)
    log10_W = np.random.normal(mu_W, alpha_W, n_slips)

    # Convertir log10(L) y log10(W) a L y W
    L = 10 ** log10_L
    W = 10 ** log10_W
    W[W > 180] = 180
    # Combinar L y W en una lista de escenarios
    scenarios = list(zip(L, W))
    return scenarios

def find_best_factors(length, width, total_subfallas):
    # Encontrar factores de total_subfallas
    factores_subfallas = np.arange(1, total_subfallas + 1)[total_subfallas % np.arange(1, total_subfallas + 1) == 0]

    # Calcular diferencias mínimas y encontrar la mejor proporción
    dx = width / factores_subfallas
    dy = length / (total_subfallas // factores_subfallas)
    diferencias = np.abs(dx - dy)
    indice_min_diferencia = np.argmin(diferencias)

    mejor_nx = factores_subfallas[indice_min_diferencia]
    mejor_ny = total_subfallas // mejor_nx
    return mejor_nx, mejor_ny, width / mejor_nx, length / mejor_ny

if __name__ == "__main__":
    desc = """
        Expected region to segmentation "lonmin" "lonmax" "latmin" "latmax" (-r), north corner fault "lat" (-nlat),
        south corner fault "lat" (slat), number of points in x(-nx), number of points in y (lat) (-ny), width of fault (-w),
        number of simulations (-n), magnitude(Mw, -m) and -p True or False if you want to show the plot
        
        FORMAT:  LAT from -90 to 90
                 LON from -180 to 180

        Example 1 input:
            python Slipgenerator_coupling.py -r -77 -69 -37 -29 -nlat -28 -slat -36 -ns 600 -n 10 -m 9.0
    """
    parser = argparse.ArgumentParser(
        description=desc, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("-n", "--nslip", dest="nslip", type=int, help="number of iterations for slip models.", required=True)
    parser.add_argument("-m", "--Mw", dest="mw", type=float, help="Moment magnitude of slip models.", required=True)
    parser.add_argument("-r", "--region", metavar=("lonmin", "lonmax", "latmin", "latmax"), dest="region", type=float, nargs=4, help="limits of map [lonmin lonmax latmin latmax] ", required=True)
    parser.add_argument("-ns", "--numbersubfaults", dest="ns", type=int, required=True, help="number of subfaults")
    parser.add_argument("-nlat", "--northlat", dest="northlat", type=float, help="north lat of possible segment", required=True)
    parser.add_argument("-slat", "--southlat", dest="southlat", type=float, help="south lat of possible segment", required=True)
    pargs = parser.parse_args()
    main(pargs)
