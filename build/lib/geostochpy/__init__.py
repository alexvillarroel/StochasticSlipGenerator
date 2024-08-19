"""
Module with functions necessary for the creation of faults and the calculation of slip distributions
    Author/Github: Alex Villarroel Carrasco / @alexvillarroel
    University of Concepcion,Chile.

    SOME FUNCTIONS ARE COPYRIGHTED BY RODRIGO CIFUENTES LOBOS
    in his github repository: @RCifuentesLobos : PaleoTsunami
    
    /


"""

import numpy as np
import warnings
from numpy import linalg as la
import matplotlib.pyplot as plt 
from mpl_toolkits.basemap import Basemap, cm
import cartopy.crs as ccr
import cartopy.feature as cf
import collections as col
from geopy import distance
import geographiclib as geo 
from geographiclib.geodesic import Geodesic
from scipy.interpolate import LinearNDInterpolator, interp1d,RegularGridInterpolator,NearestNDInterpolator
from scipy import signal
from scipy.signal import lfilter,windows
from scipy.special import kv
import scipy.ndimage.filters as filters
from clawpack.geoclaw import dtopotools, topotools
from multiprocessing import Pool, Process, cpu_count
import rockhound as rh
from scipy import linalg as la
from rockhound.slab2 import ZONES
import pygmt
import os
_ROOT = os.path.abspath(os.path.dirname(__file__))
# Ms to Mw kausel Ramirez function
def Mw_kauselramirez(ms):
    """
    calculate the Mw with a ms

    :param ms: [description]
    :type ms: [type]
    :return: [description]
    :rtype: [type]
    """
    # Ms to Mw kausel Ramirez function
    logm0kauselramirez = 1.5*ms+16.30
    Mw_calcs = 2.0/3.0*logm0kauselramirez-10.7 
    Mw_calcs = np.round(Mw_calcs,1) # 1 dec
    return Mw_calcs
# funcion para cargar datos de modelo slab2.0
def load_files_slab2(zone='south_america',rake=False):
    """

    Load data from the hayes 2018 slab2.0 model.

    Inputs:
    zone : str
    Subduction zone to fech the model.
    Available zones:
    - ``alaska``: Alaska
    - ``calabria``: Calabria
    - ``caribbean``: Caribbean
    - ``cascadia``: Cascadia
    - ``central_america``: Central America
    - ``cotabalo``: Cotabalo
    - ``halmahera``: Halmahera
    - ``hellenic``: Hellenic Arc
    - ``himalaya``: Himalaya
    - ``hindu_kush``: Hindu Kush
    - ``izu_bonin``: Izu-Bonin
    - ``kamchatka``: Kamchatka-Kuril Islands-Japan
    - ``kermadec``: Kermadec
    - ``makran``: Makran
    - ``manila_trench``: Manila Trench
    - ``muertos_trough``: Muertos Trough
    - ``new_guinea``: New Guinea
    - ``pamir``: Pamir
    - ``philippines``: Philippines
    - ``puysegur``: Puysegur
    - ``ryukyu``: Ryukyu
    - ``scotia_sea``: Scotia Sea
    - ``solomon_islands``: Solomon Islands
    - ``south_america``: South America
    - ``sulawesi``: Sulawesi
    - ``sumatra_java``: Sumatra-Java
    - ``vanuatu``: Vanuatu
    If you have and wanna load rake file,
    invoke function like load_file_slab2(zone,True)

    """
    
    # fetch slab2 model for zone with RockHound
    slab_zone=rh.fetch_slab2(zone)
    slabdep=pygmt.grd2xyz(slab_zone.depth).to_numpy()
    slabdip=pygmt.grd2xyz(slab_zone.dip).to_numpy()
    slabstrike=pygmt.grd2xyz(slab_zone.strike).to_numpy()
    #
    directory = get_data('sam_rake.xyz')
    slabdep[:,0],slabdip[:,0],slabstrike[:,0]   = slabdep[:,0] - 360, slabdip[:,0] - 360,slabstrike[:,0] - 360
    # Check if depth its in meters
    slabdep[:,2]*=-1
    ## rake
    if rake==True:
        slabrake   = np.genfromtxt(directory, delimiter = ",")
        return slabdep,slabdip,slabstrike,slabrake

    return slabdep,slabdip,slabstrike
# carga posicion de la fosa

def load_trench( trench_file ):
    """
    upload file with the location of the south american plate trench
    Input: trench file
    """

    trench = np.genfromtxt( trench_file, delimiter = " " )
    lons = trench[:,0]
    lats = trench[:,1]
    strike = trench[:,2]
    return lons, lats,strike

# carga datos del modelo PREM

def load_PREM():
    """
    carga los datos de velocidad de onda de corte del modelo prem
    directorio: directorio donde esta el archivo del modelo slab2
    """
    PREM=rh.fetch_prem(load=True)
    prem_prof = PREM.depth.to_numpy()*1000 # en metros
    prem_vs = PREM.Vsv.to_numpy()*1000 # velocidad Vsv en metros/segundos

    return prem_vs, prem_prof


# Define la falla, ancho, largo y tamano de subfallas

def crea_falla( lats, lons, prof, dip, strike, latini, latfin, area_sf, profundidad, razon_aspecto ):

    """
    Crea una falla a partir de datos del modelo Slab2.0 de Hayes, 2018
    Entradas:
    lats: latitudes de Slab2.0
    lons: longitudes de Slab2.0
    prof: profundidades de Slab2.0
    dip: dip de Slab2.0
    strike: strike de Slab2.0
    latini: latitud inicial (norte)
    latfin: latitud final (sur)
    area_sf: Area deseada para cada subfalla
    profundidad: profundidad media de la falla, sobre la cual se calculara el punto medio de la distribucion
    razon_aspecto : Razon de aspecto de la falla (largo/ancho) se calcula el ancho de la falla a partir del largo, para respetar la relacion L/A = const. de Kanamori Anderson, 1975
    """
    
    # se pasa los arrays de lats y lons a arrays unidimensionales que contienen las coordenadas sin repeticion

    # longitudes
    vector_lon_input = lons[0,:] # primera fila de matriz de lons, columnas se repiten
    # se chequea si son crecientes monotonos, util para interpolacion 
    if all( x < y for x, y in zip( vector_lon_input, vector_lon_input[1:] ) ):
        vector_lon_input = vector_lon_input
    else:
        vector_lon_input = vector_lon_input[::-1]

    # latitudes
    vector_lat_input = lats[:,0] # primera columna de matriz de lats, filas se repiten
    # se chequea si son crecientes monotonos, util para interpolacion 
    if all( x < y for x, y in zip( vector_lat_input, vector_lat_input[1:] ) ):
        vector_lat_input = vector_lat_input
    else:
        vector_lat_input = vector_lat_input[::-1]


    lim_norte = latini  # nuevo limite superior
    dif_lim_norte = np.abs( lats-lim_norte ) # diferencias entre array de latitudes y valor del limite superior
    idx_lim_norte = ( np.where( dif_lim_norte == dif_lim_norte.min() )[0][0], np.where( dif_lim_norte == dif_lim_norte.min() )[1][0] )# indice del valor de Slab2.0 que mas se aproxima 

    lim_sur = latfin  # nuevo limite inferior
    dif_lim_sur = np.abs( lats-lim_sur ) # diferencias entre array de latitudes y valor del limite inferior
    idx_lim_sur = ( np.where( dif_lim_sur == dif_lim_sur.min() )[0][0], np.where( dif_lim_sur == dif_lim_sur.min() )[1][0] )# indice del valor de Slab2.0 que mas se aproxima 

    # se calcula la distancia entre los limites (largo de la falla) en metros
    largo_falla = Geodesic.WGS84.Inverse(lats[idx_lim_norte], lons[idx_lim_norte], lats[idx_lim_sur], lons[idx_lim_sur] )[ "s12" ]
    largo_subfalla = np.sqrt( area_sf ) # subfallas cuadradas
    n_fallas_filas = np.floor_divide( largo_falla, largo_subfalla ) # cantidad de fallas en sentido norte - sur  
    # a partir del numero de fallas en el sentido norte sur (ctdad de latitudes) se crea un vector de latitudes equidistantes
    lats_fallas = np.reshape( np.linspace( lim_norte, lim_sur, int( n_fallas_filas ) ),( int( n_fallas_filas ),1 ) )
    
    # se busca la latitud del medio para referenciarla a la profundidad deseada
    if len(lats_fallas)%2 != 0:
        lat_mediana = lats_fallas[ np.floor_divide( len( lats_fallas ), 2) ]
    else:
        lat_mediana = lats_fallas[ np.floor_divide( len( lats_fallas ), 2) - 1 ]

    # busca indice de la latitud del medio
    dif_lat_mediana = np.abs( lats - lat_mediana )
    # primer indice, muestra la linea de profundidades para esta latitud
    idx_lat_mediana = np.where( dif_lat_mediana == dif_lat_mediana.min() )[0][0] 
    # se busca indice de la profundidad en la linea de la latitud media
    dif_profundidad = np.abs( profundidad - prof[ idx_lat_mediana, ] )
    idx_profundidad = np.where( dif_profundidad == dif_profundidad.min() )[0][0]
    
    # indice elemento central de la falla creada, a partir de la latitud central y la profundidad
    idx_subfalla_central = ( idx_lat_mediana, idx_profundidad )

    # longitud de la subfalla central
    lon_subfalla_central = lons[ idx_subfalla_central ]#[0][0]
    # profundidad de la subfalla central (punto con la profundidad mas cercana a la ingresada)
    prof_subfalla_central = prof[ idx_subfalla_central ]#[0][0]

    # se busca los indices de los elementos mas cercanos a las latitudes de las fallas creadas por el linespace
    dif_lats = np.ones( (len( lats_fallas ), ) + np.shape( lats ) ) # inicializacion de array para diferencias de latitudes
    for i in range( len( lats_fallas ) ):
        dif_lats[i] = np.abs( lats - lats_fallas[i] )
    
    idx_fallas = np.ones( (len( lats_fallas ), ) + ( 1,2 ) ) # inicializacion de array con los indices de las latitudes 
    for j in range( len( lats_fallas ) ):
        idx_fallas[j] = ( np.where( dif_lats[j] == dif_lats[j].min() )[0][0], np.where( dif_lats[j] == dif_lats[j].min() )[1][0] )
    
    # ancho de la falla
    ancho_falla = largo_falla/razon_aspecto
    n_fallas_columnas = np.floor_divide( ancho_falla, largo_subfalla ) # numero de subfallas en el sentido este-oeste
    # completar array de latitudes con el nuevo ancho
    #matriz_latitudes = np.reshape(np.tile(lats_fallas, int(n_fallas_columnas)),(int(n_fallas_columnas),(len(lats_fallas))))
    matriz_latitudes = np.tile( lats_fallas, int( n_fallas_columnas ) )
    # creacion de array con  longitudes a usarse
    # calculo de longitudes de los centros de las subfallas a partir del ancho de la falla
    # es necesario saber si la cantidad es par o impar
    if n_fallas_columnas%2 != 0:
        mitad_ancho = ancho_falla / 2 # en metros
        n_fallas_xlado = int( n_fallas_columnas ) // 2 # cantidad de subfallas a ambos lados de falla central
        lon_limite_oeste = Geodesic.WGS84.Direct( lat_mediana, lon_subfalla_central, 270, mitad_ancho )[ "lon2" ]
        lon_limite_este = Geodesic.WGS84.Direct( lat_mediana, lon_subfalla_central, 90, mitad_ancho )[ "lon2" ]
        lons_subfallas_oeste = np.linspace( lon_limite_oeste, lon_subfalla_central, ( n_fallas_xlado + 1 ) )
        lons_subfallas_este = np.linspace( lon_subfalla_central, lon_limite_este, ( n_fallas_xlado + 1 ) )
        lons_subfallas = np.append( lons_subfallas_oeste[:-1], lons_subfallas_este ) # vector con las longitudes de las subfallas
        lons_subfallas = np.reshape( lons_subfallas, ( 1, int( n_fallas_columnas ) ) )
    else:
        mitad_ancho = ancho_falla / 2 
        n_fallas_oeste = int( n_fallas_columnas ) / 2 - 1 # -1 para no contar 2 veces la subfalla del medio
        n_fallas_este = int( n_fallas_columnas ) / 2
        lon_limite_oeste = Geodesic.WGS84.Direct( lat_mediana, lon_subfalla_central, 270, ( mitad_ancho - largo_subfalla ) )[ "lon2" ]
        lon_limite_este = Geodesic.WGS84.Direct( lat_mediana, lon_subfalla_central, 90, mitad_ancho )[ "lon2" ]
        lons_subfallas_oeste = np.linspace( lon_limite_oeste, lon_subfalla_central, ( int( n_fallas_oeste ) + 1 ) )
        lons_subfallas_este = np.linspace( lon_subfalla_central, lon_limite_este, ( int( n_fallas_este ) + 1 ) )
        lons_subfallas = np.append( lons_subfallas_oeste[:-1], lons_subfallas_este ) # vector con las longitudes de las subfallas
        lons_subfallas = np.reshape( lons_subfallas, ( 1, int( n_fallas_columnas ) ) )

    # creacion de matriz de longitudes
    matriz_longitudes = np.tile( lons_subfallas, ( int( n_fallas_filas ), 1 ) ) # matriz con longitudes de las subfallas

    # se debe encontrar las profundidades, dips y strikes correspondientes a estas latitudes y longitudes de cada subfalla
    # profundidades correspondientes a cada subfalla:
    # se interpolara para encontrar los valores de profundidad correspondientes a cada subfalla
    
    vec_lons_subfallas_todas = np.reshape( matriz_longitudes, 
        ( int( n_fallas_filas * n_fallas_columnas ),  ) ) # vector con todos los elementos de la matriz de longitudes de las subfallas creadas
    vec_lats_subfallas_todas = np.reshape( matriz_latitudes, 
        ( int( n_fallas_filas * n_fallas_columnas ),  ) ) # vector con todos los elementos de la matriz de latitudes de las subfallas creadas


    # objeto de interpolacion de profundidades
    profs_int = RegularGridInterpolator( ( vector_lat_input, vector_lon_input ), prof )
    # inicializacion array de valores interpolados de profundidades
    prof_subfallas = np.ones( ( int( n_fallas_columnas * n_fallas_filas ), 1) )
    for p in range( int( n_fallas_columnas*n_fallas_filas ) ):
        prof_subfallas[p] = profs_int( ( vec_lats_subfallas_todas[p], vec_lons_subfallas_todas[p] ) )
    prof_subfallas = np.reshape( prof_subfallas, ( int( n_fallas_filas ), int( n_fallas_columnas ) ) )
    
    # dips correspondientes a cada subfalla:
    # se interpolara para encontrar los valores de dip correspondientes a cada subfalla

    # objeto de interpolacion de dips
    dips_int = RegularGridInterpolator( ( vector_lat_input, vector_lon_input ), dip )
    # inicializacion array de valores interpolados de dip
    dip_subfallas = np.ones( ( int( n_fallas_columnas * n_fallas_filas ), 1) )
    for d in range( int( n_fallas_columnas * n_fallas_filas ) ):
        dip_subfallas[d] = dips_int( ( vec_lats_subfallas_todas[d], vec_lons_subfallas_todas[d] ) )
    dip_subfallas = np.reshape( dip_subfallas, (int( n_fallas_filas ), int( n_fallas_columnas ) ) )
    
    # strike correspondiente a cada subfalla:
    # se interpolara para encontrar los valores de strike correspondientes a cada subfalla

    # objeto de interpolacion de strikes
    strikes_int = RegularGridInterpolator( ( vector_lat_input, vector_lon_input ), strike )
    # inicializacion array de valores interpolados de strike
    strike_subfallas = np.ones( ( int( n_fallas_columnas*n_fallas_filas ), 1) )
    for s in range( int( n_fallas_columnas*n_fallas_filas ) ):
        strike_subfallas[s] = strikes_int( ( vec_lats_subfallas_todas[s], vec_lons_subfallas_todas[s] ) )
    strike_subfallas = np.reshape( strike_subfallas, ( int( n_fallas_filas ), int( n_fallas_columnas ) ) )
    # revisar, quiza sea necesario invertir los valores de la latitud




    return largo_falla, matriz_longitudes, matriz_latitudes, prof_subfallas, dip_subfallas, strike_subfallas

# crea la falla acomodandose a la fosa
def make_fault_alongtrench(lons_trench,lats_trench,northlat, nx,ny,width,length):
    """
    Function that creates a fault by calculating the latitudes and longitudes of various subfaults, considering the trench as a boundary.
    Inputs:
    """
    dx=width/nx
    dy=length/ny
    lats=[]
    lons=[]
    lats_ep=[] # _ep : center of each subfault 
    lons_ep=[] #
    for j in range(ny):
        for i in range(nx):
            lat=northlat-km2deg(dy*(2*j+1)/2)
            lats.append(northlat-km2deg(dy*j))
            lats_ep.append(lat)
            lon=np.interp(lat,lats_trench,lons_trench)+km2deg(dx*(2*i+1)/2)
            lons.append(np.interp(lat,lats_trench,lons_trench)+km2deg(dx*i))
            lons_ep.append(lon)
    return lons,lons_ep,lats,lats_ep
# LATS SUBFALLAS
def make_fault_alongstriketrench(lons_trench, lats_trench, strike_trench, northlat, nx, ny, width, length):
    dy = length / ny
    j = np.arange(ny)
    lat_trench = northlat - km2deg(dy * j)
    lon_trench = np.interp(lat_trench, lats_trench, lons_trench)

    # Interpolar los datos
    strike = np.interp(lat_trench, lats_trench, strike_trench)

    # Crear matrices de índices
    I, J = np.meshgrid(np.arange(nx), j)

    # Repetir el array strike en las columnas nx veces
    strike_repeated = np.tile(strike[:, np.newaxis], (1, nx))

    # Calcular latitud y longitud sin usar bucles
    dx = width / nx
    lat_offset = km2deg(np.sin(np.deg2rad(strike_repeated)) * (dx * (I + 1) - 1 / 2))
    lon_offset = km2deg(np.cos(np.deg2rad(strike_repeated)) * (dx * (I + 1) - 1 / 2))

    lat = lat_trench[J] - lat_offset
    lon = lon_trench[J] + lon_offset
    # posicion aleatorea de la falla a lo largo del dip
    desplazamiento_aleatorio = np.random.uniform(0,1.71-km2deg(width))
    lon += desplazamiento_aleatorio
    lat_flat = lat.flatten()
    lon_flat = lon.flatten()
    return lon,lat,lon_flat,lat_flat

def make_fault_alongtrench_optimized(lons_trench, lats_trench, northlat, nx, ny, width, length):
    dx = width / nx
    dy = length / ny
    #
    j = np.arange(ny)
    i = np.arange(nx)
    #
    lats = northlat - km2deg(dy * j[:, None])
    lats_ep = northlat - km2deg(dy * (2*j[:, None] + 1) / 2)
    #
    lons = np.interp(lats, lats_trench, lons_trench) + km2deg(dx * i[None, :])
    lons_ep = np.interp(lats_ep, lats_trench, lons_trench) + km2deg(dx * (2 * i[None, :] + 1) / 2)
    return lons, lons_ep, lats.repeat(nx, axis=1).ravel(), lats_ep.repeat(nx, axis=1).ravel()
def interp_slabtofault(lons,lats,nx,ny,slabdep,slabdip,slabstrike,slabrake):
    # # limits stochastic grid
    lonwest=np.min(lons)
    lonest=np.max(lons)
    latnorth=np.max(lats)
    latsouth=np.min(lats)

    ind_sup_lons=np.argwhere(slabdep[:,0]<lonwest-1) # +-1 for best interpolation
    ind_inf_lons=np.argwhere(slabdep[:,0]>lonest+1)
    index_lons=np.concatenate((ind_sup_lons,ind_inf_lons))
    # delete 
    slabdep=np.delete(slabdep,index_lons,0)
    slabdip=np.delete(slabdip,index_lons,0)
    slabstrike=np.delete(slabstrike,index_lons,0)
    #
    ind_sup_lats=np.argwhere(slabdep[:,1]<(latsouth-1)) # +-1 for best interpolation
    ind_inf_lats=np.argwhere(slabdep[:,1]>(latnorth+1))
    index_lats=np.concatenate((ind_sup_lats,ind_inf_lats))
    # delete
    slabdep=np.delete(slabdep,index_lats,0)
    slabdip=np.delete(slabdip,index_lats,0)
    slabstrike=np.delete(slabstrike,index_lats,0)
    #
    mask = np.isnan(slabdep).any(axis=1)
    slabdep=slabdep[~mask]
    slabdip=slabdip[~mask]
    slabstrike=slabstrike[~mask]
    #
    X_grid,Y_grid=np.reshape(lons,(ny,nx)),np.reshape(np.sort(lats)[::-1].T,(ny,nx))
    interp_dep=LinearNDInterpolator(list(zip(slabdep[:,0],slabdep[:,1])),slabdep[:,-1])
    interp_dip=LinearNDInterpolator(list(zip(slabdip[:,0],slabdip[:,1])),slabdip[:,-1])
    interp_strike=LinearNDInterpolator(list(zip(slabdip[:,0],slabstrike[:,1])),slabstrike[:,-1])
    interp_rake=LinearNDInterpolator(list(zip(slabrake[:,0],slabrake[:,1])),slabrake[:,-1])

    dep=interp_dep(X_grid,Y_grid)
    dip=interp_dip(X_grid,Y_grid)
    strike=interp_strike(X_grid,Y_grid)
    rake=interp_rake(X_grid,Y_grid)
    return X_grid,Y_grid,dep,dip,strike,rake

# Estima el tiempo de ruptura (tau) de la falla
def interp_slab(lons,lats,slabdep,slabdip,slabstrike,slabrake):
        # # limits stochastic grid
    lonwest=np.min(lons)
    lonest=np.max(lons)
    latnorth=np.max(lats)
    latsouth=np.min(lats)

    ind_sup_lons=np.argwhere(slabdep[:,0]<lonwest) # +-1 for best interpolation
    ind_inf_lons=np.argwhere(slabdep[:,0]>lonest)
    index_lons=np.concatenate((ind_sup_lons,ind_inf_lons))
    # delete 
    slabdep=np.delete(slabdep,index_lons,0)
    slabdip=np.delete(slabdip,index_lons,0)
    slabstrike=np.delete(slabstrike,index_lons,0)
    #
    ind_sup_lats=np.argwhere(slabdep[:,1]<(latsouth)) # +-1 for best interpolation
    ind_inf_lats=np.argwhere(slabdep[:,1]>(latnorth))
    index_lats=np.concatenate((ind_sup_lats,ind_inf_lats))
    # delete
    slabdep=np.delete(slabdep,index_lats,0)
    slabdip=np.delete(slabdip,index_lats,0)
    slabstrike=np.delete(slabstrike,index_lats,0)
    #
    mask = np.isnan(slabdep).any(axis=1)
    slabdep=slabdep[~mask]
    slabdip=slabdip[~mask]
    slabstrike=slabstrike[~mask]
    interp_dep=LinearNDInterpolator(list(zip(slabdep[:,0],slabdep[:,1])),slabdep[:,2])
    interp_dip=LinearNDInterpolator(list(zip(slabdip[:,0],slabdip[:,1])),slabdip[:,2])
    interp_strike=LinearNDInterpolator(list(zip(slabdip[:,0],slabstrike[:,1])),slabstrike[:,2])
    interp_rake=LinearNDInterpolator(list(zip(slabrake[:,0],slabrake[:,1])),slabrake[:,2])
    dep=interp_dep(lons,lats)
    dip=interp_dip(lons,lats)
    strike=interp_strike(lons,lats)
    rake=interp_rake(lons,lats)
    return dep,dip,strike,rake
def tau_ruptura_BilekModel( largo_subfalla, beta = 3500,prem=False,profs=None):

    """ 
    Para modelo sin PREM, es decir, velocidad de corte constante
    #
    Calcula el tiempo de ruptura de la falla dada a partir del largo de falla, considerando
    la velocidad de ruptura como 0.8 de la velocidad de propagacion de onda S y 
    una velocidad de cizalle de 3.5 km/s utilizando la aproximacion dada en Bilek y Lay, 1999
    tau approx L/(0.8beta) 
    Entradas: 
    largo_falla: largo de la falla en metros dado por la funcion que crea la falla
    beta: velocidad de la onda de cizalle en metros/segundo
    #
    Para modelo prem, es decir, velocidad de corte variable 
    Calcula el tiempo de ruptura de la falla dada a partir del largo de falla, considerando
    la velocidad de ruptura como 0.8 de la velocidad de propagacion de onda S y 
    una velocidad de cizalle dependiente de la profundidad media de la falla
    y la velocidad vs asociada a esta en el modelo prem, utilizando la aproximacion dada en Bilek y Lay, 1999
    tau approx L/(0.8beta) 
    Entradas: 
    largo_falla: largo de la falla en metros dado por la funcion que crea la falla
    prof_media: profundidad media de la falla
    prof_prem: profundidades del modelo prem
    vs: modelo de velocidades del modelo prem asociadas a las profundidades de prof_prem
    """

    # se chequean unidades de medida de beta, tiene que estar en metros/s, no km/s
    if prem==False:
        if beta < 1000:
            beta = beta * 1000
        else: 
            beta = beta

        tau = largo_subfalla/( 0.8 * beta )
    else:
        PREM=rh.fetch_prem(load=True)
        prem_prof = PREM.depth.to_numpy()*1000 # en metros
        prem_vsv = PREM.Vsv.to_numpy()*1000 # velocidad Vsv en metros/segundos
        vs_prof_media = np.interp(profs,prem_prof,prem_vsv)

        # puede que la profundidad media de la falla se encuentre justo entre 2 valores del modelo
        # si esto sucede, se calcula el promedio de las velocidades de ambos
        tau = largo_subfalla/( 0.8 * vs_prof_media )


    return tau

# Version alternativa de estimacion de tiempo de ruptura (tau) de la falla, utilizando la profundidad de la falla y 
# la de velocidad de onda de cizalle del modelo PREM
# suele sobreestimar los valores de rigidez usualmente presentes en la zona de subduccion 
# Estima rigidez a partir de la ecuacion presentada en Bilek y Lay, 1999

def estima_rigidez_BilekModel( largo_subfalla, tau,prem=None,profs=None ):

    """
    Estima la rigidez de la falla dado el tiempo de ruptura y el largo de la falla
    segun la relacion dada en Bilek y Lay, 1999 mu = densidad*largo_falla**2/(0.8**2*tau**2).
    entradas: 
    largo_falla: largo de la falla en metros
    tau: tiempo de duracion de la ruptura
    """
    PREM=rh.fetch_prem(load=True)
    prem_prof = PREM.depth.to_numpy()*1000 # en metros
    prem_rho = PREM.density.to_numpy() * 1000 # densidad en kg/m3
    interp_rho = interp1d( prem_prof, prem_rho)
    rho = interp_rho(profs.flatten()).reshape((profs.shape))
    #
    mu = ( rho * ( largo_subfalla**2 ) )/( ( 0.8**2 ) * ( tau**2 ) )
    return mu

# estima rigidez a partir del modelo PREM

def estima_rigidez_prem(profs ):
    """
    Estima la rigidez de cada subfalla en funcion de su profunidad a partir de los 
    datos de velocidad de onda S y densdidad del modelo PREM, dada Vs=sqrt(mu/rho)
    Entradas:
    profs: profundidades de las subfallas
    """

    """
    carga los datos de velocidad de onda de corte del modelo prem
    """
    PREM=rh.fetch_prem(load=True)
    prem_prof = PREM.depth.to_numpy()*1000 # en metros
    prem_vsv = PREM.Vsv.to_numpy()*1000 # velocidad Vsv en metros/segundos
    prem_vsh = PREM.Vsh.to_numpy()*1000 # velocidad Vsh en metros/segundos
    prem_vs=np.sqrt((prem_vsv**2+prem_vsh**2)/2) # promedio de las velocidades de onda de corte
    prem_rho = PREM.density.to_numpy() * 1000 # densidad en m3/kg
    rigidez = ( prem_vsv * prem_vsv ) * prem_rho # rigidez en Pa

    # se crea el interpolador 
    interp_rigidez = interp1d( prem_prof, rigidez)

    # se inicializa el array que contendra la rigidez de cada subfalla

    rigidez_sf=interp_rigidez(profs.flatten()).reshape((profs.shape))
    return rigidez_sf


# calcula las distancias entre subfallas
# sin uso, utilizar version alternativa
def dist_haversine(lon1, lon2, lat1, lat2):
    # Convertir grados a radianes
    lon1, lon2, lat1, lat2 = np.radians(lon1), np.radians(lon2), np.radians(lat1), np.radians(lat2)

    # Diferencias en longitudes y latitudes
    dlon = lon2 - lon1
    dlat = lat2 - lat1

    # Fórmula de Haversine
    a = np.sin(dlat/2)**2 + np.cos(lat1) * np.cos(lat2) * np.sin(dlon/2)**2
    c = 2 * np.arcsin(np.sqrt(a))

    # Radio de la Tierra en metros (aproximación para latitudes medias)
    radius = 6371000

    # Distancia en metros
    distance = radius * c

    return distance
def dist_sf_vectorized(lon1, lon2, lat1, lat2):
    """
    Compute the distance between two points on the sphere .

    :param lon1: [description]
    :type lon1: [type]
    :param lon2: [description]
    :type lon2: [type]
    :param lat1: [description]
    :type lat1: [type]
    :param lat2: [description]
    :type lat2: [type]
    :return: [description]
    :rtype: [type]
    """
    # Vectorizar la función de distancia geodésica
    vectorized_dist = np.vectorize(lambda lon1, lon2, lat1, lat2: Geodesic.WGS84.Inverse(lat1, lon1, lat2, lon2)[ "s12" ])

    # Calcular las distancias para cada par de coordenadas
    distances = vectorized_dist(lon1, lon2, lat1, lat2)

    return distances
def dist_sf( lon1, lon2, lat1, lat2 ):

    """
    Calcula la distancia en metros entre dos subfallas (sf) i y j utilizando el metodo de karney, 2013 y el elipsoide wgs84
    Entradas:
    lon1: longitud de subfalla i-esima
    lon2: longitud de subfalla j-esima
    lat1: latitud de subfalla i-esima
    lat2: latitud de subfalla j-esima
    """
    #dx=deg2km(np.abs(lon1-lon2))*1000
    #dy=deg2km(np.abs(lat1-lat2))*1000
    #distancia=np.sqrt(dx**2+dy**2)
    subfalla_i = (lon1, lat1)
    subfalla_j = (lon2, lat2)
    distancia = distance.distance( subfalla_i, subfalla_j ).meters
    return distancia
# metodo alternativo de calculo de distancias entre subfallas
def dist_sf_alt( lon1, lon2, lat1, lat2 ):

    """
    Calcula la distancia en metros entre dos subfallas (sf) i y j resolviendo el problema inverso 
    de la geodesia. Utiliza el elipsoide WGS84. 
    Entradas:
    lon1: longitud de subfalla i-esima
    lon2: longitud de subfalla j-esima
    lat1: latitud de subfalla i-esima
    lat2: latitud de subfalla j-esima
    """

    dist = Geodesic.WGS84.Inverse( lat1, lon1, lat2, lon2 )[ "s12" ]
    
    return dist

# taper para medias

def taper_LeVeque( d,dmax):
    
    """
    realiza el taper propuesto por leveque et al., 2016
    t(d)=1-exp(-20|d-dmax|/dmax) 
    donde dmax es la profunidad maxima y de la falla y d es la profundidad
    de cada subfalla
    Entrada:
    d: matriz con profundidades
    """
    d[d>=dmax]=dmax
    taper = ( 1 - np.exp( -20 * abs( d - dmax )/dmax ) )
    taper[taper==0]=np.min(taper[taper!=0])
    return taper
# calculo de medias
def matriz_medias( media, prof ):
    """
    calcula las medias mu_i para cada subfalla segun mu_i = log(mu*tau_i)-1/2log(alpha**2+1)
    Entradas:
    media: valor promedio alrededor del que se desea que se centre el slip
    prof: matriz con profundidades de cada subfalla, necesaria para el taper
    """
    tau = taper_LeVeque( prof,55000) # taper
    alpha = 0.5 # valor sugerido por LeVeque et al, 2016
    mu = np.log( media*tau ) - 1/2 * np.log( alpha**2+1 )
    return mu
def matriz_medias_villarroel(media,taper,alpha=0.5):
    #alpha = 0.75 valor sugerido por LeVeque et al, 2016
    taper[taper==0]=np.min(taper[taper!=0])
    mu = np.log( media*taper) - 1/2 * np.log( alpha**2+1 )
    return mu
# calculo matriz de covarianza

# Matriz de covarianza_optimized
def matriz_covarianza_optimized( dip, prof, lons, lats,largo_falla,ancho_falla,alpha=0.5 ):

    """
    calcula la matriz de correlacion:
    Cij = exp(-(dstrike(i,j)/rstrike)-(ddip(i,j)/rdip)
    donde:
    dstrike y ddip son estimados de las distancias entre las subfallas i & j
    ddip = dprof/sin(dip)
    d es la distancia euclidiana entre dos puntos aproximada segun sus lats y lons
    rstrike y rdip son los largos de correlacion en cada direccion 
    Entradas: 
    dip: matriz con los dips de cada punto
    prof: matriz con las profundidades de cada punto
    lons: matriz con longitud de cada subfalla, necesaria para el calculo de las distancias entre estas
    lats: matriz con latitud de cada subfalla, necesaria para el calculo de las distancias entre estas
    Observacion: en esta primera instancia se considera la aproximacion  1 grado approx 111.12 km (generalizar para fallas mas grandes)
    Observacion 2: aproximacion para calculo de distancia corregida con calculo geodesico de distancia (aumenta el costo computacional bastante)
    Obeservacion 3: el ancho y largo de la falla debe ser en metros
    """

    rdip = 0.4*ancho_falla # largo de correlacion en direccion de dip
    rstrike = 0.4*largo_falla # largo de correlacion en direccion de strike
    n_filas = np.shape( prof )[0] # dimension 0
    n_columnas = np.shape( prof )[1] # dimension 1
    # Vectorizar cálculos
    vector_dip, vector_prof, vector_lon, vector_lat = dip.flatten(), prof.flatten(), lons.flatten(), lats.flatten()
    ddip = (vector_prof[:, np.newaxis] - vector_prof) / np.sin(np.deg2rad((vector_dip[:, np.newaxis] + vector_dip) / 2))
    # Reemplazar infinitos y NaN por 0
    ddip = np.nan_to_num(ddip, nan=0.0, posinf=0.0, neginf=0.0)

    # d
    d = dist_haversine( vector_lon, vector_lon[:,np.newaxis], vector_lat, vector_lat[:,np.newaxis])

    dstrike = np.sqrt(np.abs(d**2 - ddip**2))
    # Calcular Cij de manera vectorizada
    C = np.exp(-((dstrike / rstrike) + (ddip / rdip)))
    
    # Matriz de covarianza
    mat_cova = np.log(alpha**2 * C + 1)

    return mat_cova

# distribucion de slip
# Función de correlación von Karman
def von_karman_correlation(r, H):
    r = np.nan_to_num(r, nan=0.0, posinf=1e10, neginf=-1e10)
    H = np.nan_to_num(H, nan=0.0, posinf=1e10, neginf=-1e10)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        return (r ** H) * kv(H, r)

# Calculo matriz de covarianza usando VK-ACF
def matriz_covarianza_von_karman(dip, prof, lons, lats,length_fault,width_fault, H=0.7):
    """
    Calcula la matriz de correlación utilizando la función von Karman.
    length_fault and width_fault en km, prof in meters
    """
    # Vectorizar cálculos
    vector_dip, vector_prof, vector_lon, vector_lat = dip.flatten(), prof.flatten(), lons.flatten(), lats.flatten()
    
    # Calculo de ddip
    ddip = (vector_prof[:, np.newaxis] - vector_prof) / np.sin(np.deg2rad((vector_dip[:, np.newaxis] + vector_dip) / 2))
    ddip = np.abs(np.nan_to_num(ddip, nan=0.0, posinf=0.0, neginf=0.0))
    
    # Calculo de d
    d = dist_haversine(vector_lon, vector_lon[:, np.newaxis], vector_lat, vector_lat[:, np.newaxis])
    dstrike = np.sqrt(np.abs(d**2 - ddip**2))
    # Calculo de correlation length . 
    a_s=17.7+0.35*length_fault
    a_d=6.7+0.41*width_fault
    # Calculo de r_ij

    r_ij = np.sqrt((dstrike / a_s) ** 2 + (ddip / a_d) ** 2)/1000 # en km
    # Calculo de Cij usando la función von Karman
    C = np.where(r_ij!=0,von_karman_correlation(r_ij, H) / von_karman_correlation(1e-10,H),1)
    # Matriz de covarianza
    mat_cova = np.log(C + 1)
    mat_cova[np.isnan(mat_cova)]=0
    return mat_cova

def distribucion_slip( C, mu, N ):
    """
    Calcula la distribucion de slip con la expansion de karhunen-loeve log-normal
    exp(mu)exp(sum(zk*sqrt(lambdak)*vk))
    Entradas: 
    C: matriz de covarianza, se le calculan los eigenvalues y vectors
    mu: matriz con medias
    N: numero de modos a contar para la sumatoria
    """

    n_cols_cova = np.shape( C )[0] # C es una matriz cuadrada de dim x dim 

    dim_mu = np.shape( mu ) # dimensiones de mu, para volver a rearmar matriz de slip
    mu = np.reshape( mu, ( n_cols_cova, ) )
    # se calculan los valores y vectores propios de la matriz de covarianza
    eig_val, eig_vecs = la.eigh( C )
    eig_val = eig_val.real # valores propios (lambda_k)
    eig_vecs = eig_vecs.real # vectores propios (v_k) (columnas de eig_vecs)

    z = np.random.normal( 0, 1, n_cols_cova ) # distribucion gaussiana aleatoria z~N(0,1)

    # iniciacion array de slip
    S = np.ones( ( n_cols_cova, N+1 ) )
    for i in np.arange(1,N+1):
      S[:,i] = z[i]*np.sqrt( np.abs(eig_val[i]) )*np.real(eig_vecs[:,i])
    S = np.multiply( np.exp(mu), np.exp( np.sum( S, axis = 1 ) ) )
    S = np.reshape( S, dim_mu )
    return S



def distribucion_slip_optimizada(C, mu, N):
    """
    Calcula la distribucion de slip con la expansion de karhunen-loeve log-normal
    exp(mu)exp(sum(zk*sqrt(lambdak)*vk))
    Entradas: 
    C: matriz de covarianza, se le calculan los eigenvalues y vectors
    mu: matriz con medias
    N: numero de modos a contar para la sumatoria
    """

    n_cols_cova = np.shape(C)[0]  # C es una matriz cuadrada de dim x dim 
    dim_mu = np.shape(mu)  # dimensiones de mu, para volver a rearmar matriz de slip
    mu = np.reshape(mu, (n_cols_cova,))  # reshaping mu

    # se calculan los valores y vectores propios de la matriz de covarianza
    eig_val, eig_vecs = la.eigh(C)
    eig_val = eig_val.real  # valores propios (lambda_k)
    eig_vecs = eig_vecs.real  # vectores propios (v_k) (columnas de eig_vecs)

    z = np.random.normal(0, 1, n_cols_cova)  # distribucion gaussiana aleatoria z~N(0,1)

    # Cálculo vectorizado del array de slip
    sqrt_eig_val = np.sqrt(np.abs(eig_val[:N]))
    S = np.exp(mu) * np.exp(np.sum(z[:N] * sqrt_eig_val * eig_vecs[:, :N], axis=1))

    # Reshape final para volver a la dimensión original de mu
    S = np.reshape(S, dim_mu)
    return S
# calculo de magnitud

def magnitud_momento(slip, prof, largo_subfalla,ancho_subfalla,prem=False):
    """
    Calcula la magnitud de momento Mw de la distribucion de slip con la ecuacion
    Mw = 2/3log(Mo)-6.05
    donde Mo = mu * area de ruptura * slip
    Entradas: 
    slip: corresponde a la distribucion de slip 
    prof: matriz con profundidades
    lons: matriz con longitudes
    lats: matriz con latitudes
    """

    #mu = 3*10**10 # rigidez en N/m
    # se estima el tiempo que tarda la ruptura, tau, para luego estimar la rigidez
    tau = tau_ruptura_BilekModel( largo_subfalla,beta=2500,prem=prem,profs=prof)
    # se estima la rigidez de la interfaz dependiendo de la profundidad
    rigidez = estima_rigidez_BilekModel( largo_subfalla, tau,prem=prem,profs=prof )
    # calculo Mo
    #
    area_subfalla = largo_subfalla*ancho_subfalla
    Mo = rigidez*area_subfalla*slip
    Mo = np.sum(Mo)
    Mw = 2.0/3.0*np.log10(Mo)-2/3*9.1

    return Mw, area_subfalla,rigidez

# estima la magnitud de momento utilizando la rigidez obtenida del modelo PREM

# calcula la media de slip
def media_slip(Mw,largo_subfalla,ancho_subfalla,prof):
    """
    Largo, ancho y profundidad deben estar en metros
    """
    tau=tau_ruptura_BilekModel( largo_subfalla,beta=2500,prem=True,profs=prof)
    rigidez=estima_rigidez_BilekModel(largo_subfalla,tau,prem=True,profs=prof)
    Mo=10**(3/2*Mw+9.1)
    matriz_area=np.ones_like(prof)*largo_subfalla*ancho_subfalla
    media_slip=Mo/np.sum(matriz_area*rigidez)
    return media_slip/np.size(matriz_area),rigidez
# escala la magnitud
def escalar_magnitud_momento(Mw, slip, prof, largo_subfalla,ancho_subfalla,prem=False):
    """
    Escala la distribucion de slip a la magnitud de momento deseada de entrada.
    crea una matriz de slip uniforme con magnitud de momento igual a la primitiva
    para normalizar el slip. Despues escala con matriz uniforme con magnitud de momento
    deseada
    Entradas:
    Mw: magnitud de momento deseada
    slip: distribucion "primitiva" de slip
    prof: matriz de profundidades de la interfaz
    lons: longitudes de las subfallas
    lats: latitudes de las subfallas
    observacion: los 4 ultimos argumentos son para calcular la magnitud de momento 
    con la funcion "magnitud_momento" definida anteriormente
    """
    # dimensiones array de slip
    dims                = np.shape( slip )
    if prem==False:    
        Mw_prim, areasf,rigidez = magnitud_momento( slip, prof, largo_subfalla,ancho_subfalla,prem=False ) # magnitud de momento de dist slip primitiva
    else:
        Mw_prim, areasf,rigidez = magnitud_momento( slip, prof, largo_subfalla,ancho_subfalla,prem=True ) # magnitud de momento de dist slip primitiva
    Mo_original     = 10**(3./2. *  Mw_prim + 9.1 ) # momento de la distribucion original
    Mo_deseado      = 10**(3./2. *  Mw + 9.1 )  # momento deseado para la magnitud a escalar
    razon_Mo        = Mo_deseado/Mo_original # razones entre los momentos, se cancelan los otros factores que no sean slip (rigidez y area)
    slip_escalado   = slip*razon_Mo
    #D_original = np.divide( ( Mo_original/n_subfallas ) * slip_prim_homogeneo, (rigidez*areasf)*np.ones( dims )   )
    #slip_prim_homogeneo = np.exp( ( 3.0/2.0 ) * ( Mw_prim+6.05 ) )/area_falla # array de slip homogeneo con Mw primitiva
    #slip_norm = np.divide( slip, slip_prim_homogeneo ) # array slip normalizado por slip con Mw primitiva
    #slip_deseado_homogeneo = np.exp( ( 3.0/2.0 ) * ( Mw + 6.05 ) )/area_falla # array de slip homogeneo con Mw deseado
    #slip_escalado = np.multiply(slip_norm, slip_deseado_homogeneo)
    return slip_escalado,rigidez,Mo_original,Mo_deseado


def corr2d_fourier(slip, test, lonslip = [], latslip = [], lontest = [], lattest = [], normalizar = True):
    """
    Calcula la correlacion cruzada entre 2 arrays 2d utilizando el metodo de la 
    transformada rapida de fourier para obtener indices para filtrar
    Entradas:
    slip: array que contiene el slip generado
    test: arrray de prueba con el cual se busca realzar el filtro (e.g.: anom. de gravedad)
    lonslip, latslip: longitud y latitud de las subfallas que componen el slip
    lontest, lattest: longitud y latitud de las subfallas que componen el test
    """


    if normalizar: 
        xcorr2d = signal.correlate(slip/np.size(slip), test/np.size(test), mode='full', method='fft')
    else:
        xcorr2d = signal.correlate(slip, test, mode='full', method='fft')
    
    
    valor_corr = np.sum(xcorr2d)


    return valor_corr

def sigmoid(x):
    """
    Funcion auxiliar para crear ventana sigmoidal para el taper
    """
    S = np.ones((np.size(x),))
    for i in range(len(x)):
        S[i] = 1/(1+np.exp(-x[i]))
    return S
# calculo laplaciano para estudiar la suavidad

def calcula_laplaciano(tamano_sf,F):

    """
    Entradas 
    tamano_sf: tamano subfalla
    F: Campo escalar a calcularle el laplaciano
    """

    step = np.sqrt(tamano_sf)

    laplaciano = filters.laplace(F,mode="wrap")
    sum_lap    = np.sum(laplaciano)

    return sum_lap

# calculo de la cantidad maxima de valores y vectores propios a utilizarse

def N_max_matriz_covarianza(C):
    """
    Busca la esquina de la curva l de los valores propios
    Entrada:
    C : matriz de covarianza de la falla
    """
    # valores auxiliares
    n_filas_cova = np.shape( C )[0]
    n_cols_cova  = np.shape( C )[1]
    # valores y vectores propios
    eig_val, eig_vec = la.eig( C )
    eig_vals         = eig_val.real # valores propios (lambda_k)
    eig_vecs         = eig_vec.real # vectores propios (v_k) (columnas de eig_vecs)
    # vector con indices de valores propios
    idx_vec = np.arange((np.size(eig_vals)))
    # transformamos los valores propios y los indices a un espacio log-log
    x = idx_vec+1
    y = np.abs(eig_vals)
    # Triangular/circumscribed circle simple approximation to curvature 
    # (after Roger Stafford)

    # the series of points used for the triangle/circle
    x1 = x[:-2]
    x2 = x[1:-1]
    x3 = x[2:]
    y1 = y[:-2]
    y2 = y[1:-1]
    y3 = y[2:]

    # the side lengths for each triangle
    a = np.sqrt(np.square(x3-x2)+np.square(y3-y2))
    b = np.sqrt(np.square(x1-x3)+np.square(y1-y3))
    c = np.sqrt(np.square(x2-x1)+np.square(y2-y1))
    # semi perimetro
    s = (a+b+c)/2.
    # radio de cada circulo
    R = (a*b*c)/(4*np.sqrt((s*(s-a)*(s-b)*(s-c))))
    # The curvature for each estimate for each value which is
    # the reciprocal of its circumscribed radius. Since there aren't circles for 
    # the end points they have no curvature
    kappa       = np.ones((n_filas_cova))
    kappa[0]    = 0.
    kappa[-1]   = 0.
    kappa[1:-1] = np.reciprocal(R)
    idx_max = np.where(kappa == np.max(kappa))[0][0] - 1
    return idx_max

def km2deg(km):
    return km/111.1
#
def deg2km(deg):
    return float(deg*111.1)
# based in Melgar and Hayes 2017. along strike and along-dip values
def nearest_value(value, matrix):
    return matrix.flat[np.abs(matrix - value).argmin()]

def hypocenter(lons,lats,depth,length,width):
    """
    Generate a random hypocenter from a given grid .

    :param lons: [description]
    :type lons: [type]
    :param lats: [description]
    :type lats: [type]
    :param depth: [description]
    :type depth: [type]
    :param length: [description]
    :type length: [type]
    :param width: [description]
    :type width: [type]
    :return: [description]
    :rtype: [type]
    """
    # along dip normal pdf parameters
    mu=0.046
    sigma=0.208
    # along strike exponential pdf parameters
    lamda=5.435
    #
    along_dip=np.random.normal(loc=mu,scale=sigma)
    along_strike=np.random.exponential(scale=1/lamda)
    while (along_dip < -0.5) or (along_dip > 0.5):
        along_dip=np.random.normal(loc=mu,scale=sigma)
    while (along_strike > 0.5) or (along_strike < 0):
        along_strike=np.random.exponential(scale=1/lamda)
    # once computed a random value for each component, its neccesary 
    # transform to lat and lon hypocenter.

    L=along_strike*length # linear regression
    W=along_dip*width+width/2 # linear regression
    flag=np.random.randint(0,1)
    # if its 0, the hypocenter will to north of fault
    # else, to south
    middle= (np.max(lats[:])+np.min(lats[:]))/2
    if flag == 0:
        y=middle+km2deg(L)
        nearest_y_pos=np.argwhere(lats==nearest_value(y,lats[:]))
        x_min=lons[nearest_y_pos[0][0]][0] # catch the first column
        x=x_min+km2deg(W)
        nearest_x_pos=np.argwhere(lons==nearest_value(x,lons[:]))
        z=depth[nearest_x_pos[0][0]][nearest_x_pos[0][1]]
    else:
        y=middle-km2deg(L)
        nearest_y_pos=np.argwhere(lats==nearest_value(y,lats[:]))
        x_min=lons[nearest_y_pos[0][0]][0]
        x=x_min+km2deg(W)
        nearest_x_pos=np.argwhere(lons==nearest_value(x,lons[:]))
        z=depth[nearest_x_pos[0][0]][nearest_x_pos[0][1]]
    hypocenter=[x,y,z]
    return hypocenter
################ FILTER FUNCTIONS
# crea ventana para hacerle un taper al patron de slip
def ventana_taper_slip_fosa( Slip, lon, lat, ventana_flag ):
    """
    Crea una ventana para aplicar un taper al slip en cercanias de la fosa
    Entradas:
    Slip: array con los valores de slip por subfalla
    lon: array de longitudes de las subfallas
    lat: array con las latitudes de las subfallas
    """

    # valores auxiliares

    n_subfallas = np.size(Slip)
    n_lons      = np.shape(Slip)[1]
    n_lats      = np.shape(Slip)[0]
    ventana     = []

    # seleccionar la ventana
    # 1: tukey
    # 2: gauss
    # 3: hanning
    # 4: hamming

    if   ventana_flag == 1:
        ventana     = signal.tukey(n_lons, 0.02)
    elif ventana_flag == 2:
        ventana     = signal.gaussian(n_lons,std=7)
    elif ventana_flag == 3:
        ventana     = np.hanning(n_lons)
    elif ventana_flag == 4:
        ventana     = np.hamming(n_lons)
    elif ventana_flag == 5:
        n = np.linspace(1,100,n_lons)
        ventana     = sigmoid(n)

        
    ventana[int(n_lons/2):n_lons] = 1
    ventana = np.tile(ventana,(n_lats,1))

    return ventana


# crea el taper al patron de slip a partir de la ventana creada con la funcion ventana_taper_slip_fosa()
def taper_borders(Slip,taperfunc=np.hanning):
    filas, columnas = Slip.shape
    ventana_filas = taperfunc(filas)
    ventana_columnas = taperfunc(columnas)
    taper_2d = np.outer(ventana_filas, ventana_columnas)
    Slip_taper = Slip * taper_2d

    return Slip_taper
def taper_except_trench(Slip,taperfunc=np.hanning):
    filas, columnas = Slip.shape
    ventana_filas = taperfunc(filas)
    ventana_columnas = taperfunc(columnas)
    # ventana_columnas[0:int(columnas)//2] = np.linspace(0.2,1,np.size(ventana_columnas[0:int(columnas)//2]))
    ventana_columnas[0:int(columnas)//2] = 1
    taper_2d = np.outer(ventana_filas, ventana_columnas)
    Slip_taper = Slip * taper_2d

    return Slip_taper
def taper_slip_fosa(Slip, ventana):
    """
    Aplica un taper para atenuar el slip en las cercanias de la fosa
    Entradas:
    Slip: array con los valores de slip por subfalla
    ventana: ventana a usarse para aplicar el taper
    """

    # valores auxiliares
    n_subfallas    = np.size(Slip)
    n_lons         = np.shape(Slip)[1]
    n_lats         = np.shape(Slip)[0]
    tamano_ventana = np.shape(ventana)
    tamano_slip    = np.shape(Slip)

    # se asegura que la ventana tenga el mismo tamagno que el slip
    assert tamano_ventana == tamano_slip, "ventana mal creada"

    tapered_slip = np.multiply(Slip,ventana)
    
    return tapered_slip
# FILTER ALONG DIP IN SUBDUCTION ZONE ()
def tukey_window(N, alpha=0.5):
    """
    Genera una ventana Tukey.Se ha
    ajustado los valores de la parte de Hann para que no alcancen completamente cero al comienzo

    Parámetros:
    - N: Longitud de la ventana.
    - alpha: Parámetro de apertura (0 para una ventana rectangular, 1 para una ventana Hann).

    Retorna:
    - Ventana Tukey.
    """
    if alpha <= 0:
        return np.ones(N)
    elif alpha >= 1:
        return np.hanning(N)
    else:
        
        x = np.linspace(0, 1, N, endpoint=False)
        w = np.ones_like(x)

        # Aplica la parte de la ventana Hann
        first_condition = x < alpha / 2
        last_condition = x >= 1 - alpha / 2

        w[first_condition] = 0.5 * (1 + np.cos(2 * np.pi / alpha * (x[first_condition] - alpha / 2)))
        w[last_condition] = 1

        # Ajusta los valores para que no alcancen completamente cero al comienzo y al final
        w[first_condition] = 0.5 + 0.5 * w[first_condition]
        w[last_condition] = 0 + 1 * w[last_condition]

        return w
def tukey_window_equal(N, alpha=0.5):
    """
    Genera una ventana Tukey.Se ha
    ajustado los valores de la parte de Hann para que no alcancen completamente cero al comienzo

    Parámetros:
    - N: Longitud de la ventana.
    - alpha: Parámetro de apertura (0 para una ventana rectangular, 1 para una ventana Hann).

    Retorna:
    - Ventana Tukey.
    """
    if alpha <= 0:
        return np.ones(N)
    elif alpha >= 1:
        return np.hanning(N)
    else:
        
        x = np.linspace(0, 1, N, endpoint=False)
        w = np.ones_like(x)

        # Aplica la parte de la ventana Hann
        first_condition = x < alpha / 2
        last_condition = x >= 1 - alpha / 2

        w[first_condition] = 0.5 * (1 + np.cos(2 * np.pi / alpha * (x[first_condition] - alpha / 2)))
        w[last_condition] = 0.5 * (1 + np.cos(2 * np.pi / alpha * (x[last_condition] - 1 + alpha / 2)))

        # Ajusta los valores para que no alcancen completamente cero al comienzo y al final
        w[first_condition] = 0.5 + 0.5 * w[first_condition]
        w[last_condition] = 0.5 + 0.5 * w[last_condition]

        return w
def taper_except_trench_tukey(depth,dip_taperfunc=tukey_window,strike_taperfunc=tukey_window_equal,alpha_dip=0.5,alpha_strike=0.5):
    """
    2-D Taper to a Slip. Tukey windows along the strike and along the dip.

    :param Slip: [description]
    :type Slip: [type]
    :param dip_taperfunc: [description], defaults to tukey_window
    :type dip_taperfunc: [type], optional
    :param strike_taperfunc: [description], defaults to windows.tukey
    :type strike_taperfunc: [type], optional
    :param alpha_dip: [description], defaults to 0.5
    :type alpha_dip: float, optional
    :param alpha_strike: [description], defaults to 0.5
    :type alpha_strike: float, optional
    :return: [description]
    :rtype: [type]
    """
    filas, columnas = depth.shape
    ventana_filas = strike_taperfunc(filas,alpha_strike)
    ventana_columnas = dip_taperfunc(columnas,alpha_dip)
    # ventana_columnas[0:int(columnas)//2] = np.linspace(0.2,1,np.size(ventana_columnas[0:int(columnas)//2]))
    # ventana_columnas[0:int(columnas)//2] = 1
    taper_2d = np.outer(ventana_filas, ventana_columnas)

    return taper_2d #Slip_taper, taper_2d
################ PLOT FUNCTIONS
def plot_grid(region,lons_trench,lats_trench,lons,lats):
    """
    Plot a grid of regions in a region of land .

    :param region: [description]
    :type region: [type]
    :param data_trench: [description]
    :type data_trench: [type]
    :param lons: [description]
    :type lons: [type]
    :param lats: [description]
    :type lats: [type]
    """
    # plotting grid
    fig = pygmt.Figure()
    # Chilean trench
    # Load sample grid (3 arc-minutes global relief) in target area
    grid = pygmt.datasets.load_earth_relief(resolution="30s", region=region)
    # Plot original grid
    fig.basemap(region=region, projection="M0/0/12c", frame=True)
    fig.grdimage(grid=grid, cmap="oleron")
    fig.plot(
        x=lons_trench,
        y=lats_trench,
        region=region,
        pen="0.2p",
        fill="white",
        style="f0.5i/0.1i+r+t+o1",
    )
    fig.plot(x=lons,y=lats, style="p0.2c", pen="2p,red",fill="red")
    fig.colorbar(
    # Place the colorbar inside the plot (lower-case "j") with justification
    # Bottom Right and an offset ("+o") of 0.7 centimeters and
    # 0.3 centimeters in x or y directions, respectively
    # Move the x label above the horizontal colorbar ("+ml")
    position="jBR+o0.7c/0.8c+h+w5c/0.3c+ml",
    # Add a box around the colobar with a fill ("+g") in "white" color and
    # a transparency ("@") of 30 % and with a 0.8-points thick black
    # outline ("+p")
    box="+gwhite@30+p0.8p,black",
    # Add x and y labels ("+l")
    frame=["x+lElevation", "y+lm"])
    fig.savefig('grid.png')
    fig.show()
    return
def plot_slip(X_grid,Y_grid,lonfosa,latfosa,Slip,filename,show=False,cmap='rainbow'):
        fig = plt.figure()
        # iniciliazar mapa
        m = Basemap(projection='merc', lat_0=35, lon_0=210,
            resolution = 'h',
            llcrnrlon=-78, llcrnrlat=-38,
            urcrnrlon=-68, urcrnrlat=-28)
        # transformar coordenadas geograficas a coord de mapa
        mlons, mlats         = m(X_grid,Y_grid)
        mfosalons, mfosalats = m(lonfosa, latfosa)
        # anexos
        m.plot(mfosalons, mfosalats, marker=None, color='k')
        m.drawcoastlines(color='black')
        m.drawcountries(linewidth=0.25)
        m.fillcontinents(color='gray',lake_color='aqua')
        m.drawmapboundary(fill_color='white')
        # m.drawmeridians(np.arange(-180,180,2),labels=[1,1,0,1])
        # m.drawparallels(np.arange(-50,-30,2),labels=[1,1,0,1])
        m.drawparallels(np.arange(-40, -25, 2), labels=[1,0,0,0])
        m.drawmeridians(np.arange(-80, -68, 2), labels=[0,0,0,1])
        # promedio
        my_cmap = plt.get_cmap(cmap)
        csave = m.pcolormesh(mlons, mlats, Slip,cmap=my_cmap,shading='gouraud')
        cbar     = m.colorbar(csave,location='bottom',pad="5%")
        cbar.set_label('Slip [m]')
        if show==False:
            plt.savefig(filename)
            plt.close()
            return fig
        else:
            return fig
def plot_slip_gmt(region,X_grid,Y_grid,lonfosa,latfosa,Slip,dx,dy,filename=False):
    """
    Plot a Slip grid of data on a grid

    :param region: [description]
    :type region: [type]
    :param X_grid: [description]
    :type X_grid: [type]
    :param Y_grid: [description]
    :type Y_grid: [type]
    :param lonfosa: [description]
    :type lonfosa: [type]
    :param latfosa: [description]
    :type latfosa: [type]
    :param Slip: [description]
    :type Slip: [type]
    :param dx: [description]
    :type dx: [type]
    :param dy: [description]
    :type dy: [type]
    :param archivo_salida: [description]
    :type archivo_salida: [type]
    """
    import numpy as np
    # define parameters for plotting
    latmax=np.max(Y_grid.flat)
    latmin=np.min(Y_grid.flat)
    lonmax=np.max(X_grid.flat)
    lonmin=np.min(X_grid.flat)
    reg2=(lonmin,lonmax,latmin,latmax)
    #
    fig=pygmt.Figure()
    spacing=f"{dx}k/{dy}k"
    cmap2use = 'rainbow'
    cmap=pygmt.makecpt(cmap = 
              cmap2use, 
              series = [np.min(Slip.flatten())-1, np.max(Slip.flatten())+3,5], 
              background='i',
              continuous=True)
    grid=pygmt.xyz2grd(x=X_grid.flatten(),y=Y_grid.flatten(),z=Slip.flatten(),region = reg2,spacing=spacing)
    grid=pygmt.grdsample(grid,spacing=0.001,region=reg2,verbose='q')
    cmap=pygmt.grd2cpt(grid=grid,cmap='rainbow',nlevels=True,continuous=True)
    fig.basemap(region=region, projection="Q12c", frame='ag')
    fig.coast(land="darkgray",projection='Q12c',frame=True,water="lightblue",borders=["1/0.5p,black", "2/0.5p,gray", "3/0.5p,blue"])
    fig.grdimage(grid, 
            cmap=cmap,
            projection='Q12c',
            nan_transparent=True,
            interpolation='l+a+t1.0'
            )
    # fig.grdcontour(grid=grid,projection='Q12c')
    # fig.grdcontour(grid=grid,interval=3,annotation=0)
    fig.plot(x=lonfosa,y=latfosa,
        projection='Q12c',
        region=region,
        pen="1p",
        fill="white",
        style="f0.5i/0.1i+r+t+o1")
    fig.colorbar(frame=['a+3',"x+lSlip[m]"],cmap=cmap,projection='Q12c')
    with fig.inset(position="jTL+w3.5c+o0.2c", margin=0, box="+p1.5p,gold"):
        # Create a figure in the inset using coast. This example uses the azimuthal
        # orthogonal projection centered at 47E, 20S. The land color is set to
        # "gray" and Madagascar is highlighted in "red3".
        fig.coast(
            region="g",
            projection="G-70/-33/?",
            land="gray",
            water="royalblue",
            dcw="CL+gred3",
        )
    fig.shift_origin(yshift="0c",xshift="13c")
    #
    data_lat=np.ones((np.unique(Y_grid.flatten()).size,2))
    data_lat[:,1]=np.unique(Y_grid.flatten())[::-1]
    data_lat[:,0]=np.mean(Slip,axis=1)
    fig.plot(
        projection="X2c/12c",
        # Note that the y-axis annotation "Counts" is shown in x-axis direction
        # due to the rotation caused by horizontal=True
        frame=["ag", "WSne","xaf+lSlip mean[m]"],
        region=[0,np.max(np.mean(Slip,axis=1))+1,region[2],region[3]],
        x=data_lat[:,0],
        y=data_lat[:,1],
        pen="2p,black",
        )
    fig.shift_origin(yshift='13c',xshift='-13c')
    fig.plot(
    projection="X12c/3c",
    # Note that the y-axis annotation "Counts" is shown in x-axis direction
    # due to the rotation caused by horizontal=True
    frame=["ag", "WSne","yaf+lSlip mean[m]"],
    region=[region[0],region[1],0,np.max(np.mean(Slip,axis=0))+1],
    x=np.mean(X_grid,axis=0),
    y=np.mean(Slip,axis=0),
    pen="2p,black",
    )
    np.mean(X_grid,axis=0),np.mean(Slip,axis=0)
    # Shift the plot origin and add right margin histogram
    if filename != False:
        fig.savefig(filename) 
    else:
        fig.show()
    return

def plot_slip_gmt_relief(region,X_grid,Y_grid,lonfosa,latfosa,Slip,dx,dy,filename=False):
    """
    Plot a Slip grid of data on a grid

    :param region: [description]
    :type region: [type]
    :param X_grid: [description]
    :type X_grid: [type]
    :param Y_grid: [description]
    :type Y_grid: [type]
    :param lonfosa: [description]
    :type lonfosa: [type]
    :param latfosa: [description]
    :type latfosa: [type]
    :param Slip: [description]
    :type Slip: [type]
    :param dx: [description]
    :type dx: [type]
    :param dy: [description]
    :type dy: [type]
    :param archivo_salida: [description]
    :type archivo_salida: [type]
    """
    import numpy as np
    # define parameters for plotting
    latmax=np.max(Y_grid.flat)
    latmin=np.min(Y_grid.flat)
    lonmax=np.max(X_grid.flat)
    lonmin=np.min(X_grid.flat)
    reg2=(lonmin,lonmax,latmin,latmax)
    #
    fig=pygmt.Figure()
    spacing=f"{dx}k/{dy}k"
    cmap2use = 'rainbow'
    relief=pygmt.datasets.load_earth_relief(resolution='05m')
    cmap=pygmt.makecpt(cmap = 
              cmap2use, 
              series = [np.min(Slip.flatten())-1, np.max(Slip.flatten())+3,5], 
              background='i',
              continuous=True)
    grid=pygmt.xyz2grd(x=X_grid.flatten(),y=Y_grid.flatten(),z=Slip.flatten(),region = reg2,spacing=spacing)
    # grid=pygmt.grdsample(grid,spacing=spacing2,region=reg2,verbose='q')
    cmap=pygmt.grd2cpt(grid=grid,cmap='rainbow',nlevels=True,continuous=True)
    fig.basemap(region=region, projection="Q12c", frame='ag')
    fig.coast(land="darkgray",projection='Q12c',frame=True,water="lightblue",borders=["1/0.5p,black", "2/0.5p,gray", "3/0.5p,blue"])
    fig.grdimage(relief,projection='Q12c')
    fig.grdimage(grid, 
            cmap=cmap,
            projection='Q12c',
            nan_transparent=True,
            interpolation='l+a+t1.0'
            )
    fig.grdcontour(grid=grid,projection='Q12c')
    # fig.grdcontour(grid=grid,interval=3,annotation=0)
    fig.plot(x=lonfosa,y=latfosa,
        projection='Q12c',
        region=region,
        pen="1p",
        fill="white",
        style="f0.5i/0.1i+r+t+o1")
    fig.colorbar(frame=['a+3',"x+lSlip[m]"],cmap=cmap,projection='Q12c')
    with fig.inset(position="jTL+w3.5c+o0.2c", margin=0, box="+p1.5p,gold"):
        # Create a figure in the inset using coast. This example uses the azimuthal
        # orthogonal projection centered at 47E, 20S. The land color is set to
        # "gray" and Madagascar is highlighted in "red3".
        fig.coast(
            region="g",
            projection="G-70/-33/?",
            land="gray",
            water="royalblue",
            dcw="CL+gred3",
        )
    fig.shift_origin(yshift="0c",xshift="13c")
    #
    data_lat=np.ones((np.unique(Y_grid.flatten()).size,2))
    data_lat[:,1]=np.unique(Y_grid.flatten())[::-1]
    data_lat[:,0]=np.mean(Slip,axis=1)
    fig.plot(
        projection="X2c/12c",
        # Note that the y-axis annotation "Counts" is shown in x-axis direction
        # due to the rotation caused by horizontal=True
        frame=["ag", "WSne","xaf+lSlip mean[m]"],
        region=[0,np.max(np.mean(Slip,axis=1))+1,region[2],region[3]],
        x=data_lat[:,0],
        y=data_lat[:,1],
        pen="2p,black",
        )
    fig.shift_origin(yshift='13c',xshift='-13c')
    fig.plot(
    projection="X12c/3c",
    # Note that the y-axis annotation "Counts" is shown in x-axis direction
    # due to the rotation caused by horizontal=True
    frame=["ag", "WSne","yaf+lSlip mean[m]"],
    region=[region[0],region[1],0,np.max(np.mean(Slip,axis=0))+1],
    x=np.mean(X_grid,axis=0),
    y=np.mean(Slip,axis=0),
    pen="2p,black",
    )
    np.mean(X_grid,axis=0),np.mean(Slip,axis=0)
    # Shift the plot origin and add right margin histogram
    if filename != False:
        fig.savefig(filename) 
    else:
        fig.show()
    return
def plot_slip_css(region,lons, lats, lonfosa, latfosa, Slip, cmap='rainbow'):
    """
    Plots the CSS projection of a SlipLine

    :param region: [description]
    :type region: [type]
    :param lons: [description]
    :type lons: [type]
    :param lats: [description]
    :type lats: [type]
    :param lonfosa: [description]
    :type lonfosa: [type]
    :param latfosa: [description]
    :type latfosa: [type]
    :param Slip: [description]
    :type Slip: [type]
    :param cmap: [description], defaults to 'rainbow'
    :type cmap: str, optional
    """
    ax = plt.axes(projection=ccrs.Mercator())
    ax.set_extent(region)
    ax.gridlines(draw_labels=True, dms=True, x_inline=False, y_inline=False)
    # Transformar coordenadas geográficas a coordenadas del mapa
    # mlons, mlats = ax.projection.transform_points(ccrs.Mercator(), lons, lats)
    # mfosalons, mfosalats = ax.projection.transform_point(ccrs.Mercator(),lonfosa, latfosa)
    # Crear la visualización de los datos de Slip
    my_cmap = plt.get_cmap(cmap)
    csave = ax.pcolormesh(np.unique(lons), np.unique(lats), Slip, cmap=my_cmap, shading='gouraud')

    # Añadir barra de color
    cbar = plt.colorbar(csave, ax=ax, orientation='horizontal', pad=0.05)
    cbar.set_label('m')

    # Añadir límites de la fosa
    # ax.plot(mfosalons, mfosalats, marker=None, color='k')

    # Añadir características geográficas
    ax.add_feature(cf.COASTLINE)
    ax.add_feature(cf.LAND)
    ax.add_feature(cf.RIVERS)
    ax.add_feature(cf.BORDERS)
    ax.gridlines(draw_labels=True, linestyle='--', linewidth=0.5)
    plt.show()
    plt.close()
    return

def plot_3d(X_array,Y_array,depth,Slip,filename=None):
    """
    Plot a 3D plot of a 3D Slip grid with a depth of the surface .

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
def plot_deformation_gmt(region,X_grid,Y_grid,dx,dy,lonfosa,latfosa,Deformation,filename):
    """
    Plot a Slip grid of data on a grid

    :param region: [description]
    :type region: [type]
    :param X_grid: [description]
    :type X_grid: [type]
    :param Y_grid: [description]
    :type Y_grid: [type]
    :param lonfosa: [description]
    :type lonfosa: [type]
    :param latfosa: [description]
    :type latfosa: [type]
    :param Slip: [description]
    :type Slip: [type]
    :param dx: [description]
    :type dx: [type]
    :param dy: [description]
    :type dy: [type]
    :param archivo_salida: [description]
    :type archivo_salida: [type]
    """
    import numpy as np
    # define parameters for plotting
    latmax=np.max(Y_grid.flat)
    latmin=np.min(Y_grid.flat)
    lonmax=np.max(X_grid.flat)
    lonmin=np.min(X_grid.flat)
    reg2=(lonmin,lonmax,latmin,latmax)
    #
    fig=pygmt.Figure()
    spacing=f"{dx}k/{dy}k"
    cmap2use = 'polar'
    cmap=pygmt.makecpt(cmap = 
              cmap2use, 
              series = [-np.max(Deformation.flatten()), np.max(Deformation.flatten()),2], 
              background=False,
              continuous=True)
    grid=pygmt.xyz2grd(x=X_grid.flatten(),y=Y_grid.flatten(),z=Deformation.flatten(),region = reg2,spacing=spacing)
    grid=pygmt.grdsample(grid,spacing=0.001,region=reg2,verbose='q')
    cmap=pygmt.grd2cpt(grid=grid,cmap='polar',nlevels=True,continuous=True)
    fig.basemap(region=region, projection="Q12c", frame='ag')
    fig.coast(land="darkgray",projection='Q12c',frame=True,water="lightblue",borders=["1/0.5p,black", "2/0.5p,gray", "3/0.5p,blue"])
    fig.grdimage(grid, 
            cmap=cmap,
            projection='Q12c',
            nan_transparent='+z0',
            interpolation='l+a+t1.0'
            )
    # fig.grdcontour(grid=grid,projection='Q12c')
    # fig.grdcontour(grid=grid,interval=3,annotation=0)
    fig.plot(x=lonfosa,y=latfosa,
        projection='Q12c',
        region=region,
        pen="1p",
        fill="white",
        style="f0.5i/0.1i+r+t+o1")
    fig.colorbar(frame=['a+3',"x+lVertical Displacement[m]"],cmap=cmap,projection='Q12c')
    #cities
    #valparaiso
    fig.text(x=-71.3,y= -33.03,text='Valparaíso',fill='white',justify='ML',pen="1p,black")
    fig.plot(x=-71.63,y=-33.03,style='c0.2c',fill='white',pen="1p,black")
    # la serena
    fig.text(x=-71.2,y= -30.03,text='La Serena',fill='white',justify='ML',pen="1p,black")
    fig.plot(x=-71.4,y=-30.03,style='c0.2c',fill='white',pen="1p,black")
    #santiago
    fig.text(x=-70.4,y= -33.45,text='Santiago',fill='white',justify='ML',pen="1p,black")
    fig.plot(x=-70.6,y=-33.45,style='c0.2c',fill='white',pen="1p,black")
    # talca
    fig.text(x=-71.6,y= -35.43,text='Talca',fill='white',justify='ML',pen="1p,black")
    fig.plot(x=-71.8,y=-35.43,style='c0.2c',fill='white',pen="1p,black")
    # concepcion
    fig.text(x=-72.8,y= -36.812,text='Concepción',fill='white',justify='ML',pen="1p,black")
    fig.plot(x=-73.039,y=-36.812,style='c0.2c',fill='white',pen="1p,black")
        # Shift the plot origin and add right margin histogram
    if filename != False:
        fig.savefig(filename) 
    else:
        fig.show()
    return
def plot_deformation(X_grid,Y_grid,lonfosa,latfosa,Deformation,filename,show=False):
    import matplotlib.colors as colors
    from matplotlib import cm
    fig = plt.figure()
    # iniciliazar mapa
    m = Basemap(projection='merc', lat_0=35, lon_0=210,
        resolution = 'h',
        llcrnrlon=-78, llcrnrlat=-38,
        urcrnrlon=-68, urcrnrlat=-28)
    # transformar coordenadas geograficas a coord de mapa
    mlons, mlats         = m(X_grid,Y_grid)
    mfosalons, mfosalats = m(lonfosa, latfosa)
    # anexos
    m.plot(mfosalons, mfosalats, marker=None, color='k')
    m.drawcoastlines(color='black')
    m.drawcountries(linewidth=0.25)
    m.fillcontinents(color='gray',lake_color='aqua')
    m.drawmapboundary(fill_color='white')
    # m.drawmeridians(np.arange(-180,180,2),labels=[1,1,0,1])
    # m.drawparallels(np.arange(-50,-30,2),labels=[1,1,0,1])
    m.drawparallels(np.arange(-40, -25, 2), labels=[1,0,0,0])
    m.drawmeridians(np.arange(-80, -68, 2), labels=[0,0,0,1])
    # promedio
    csave = m.pcolormesh(mlons, mlats, Deformation,norm=colors.CenteredNorm(),cmap='RdBu_r',shading='gouraud')
    cbar     = m.colorbar(csave,location='bottom',pad="5%")
    cbar.set_label('Seafloor deformation (dZ) [m]')
    if show==False:
        plt.savefig(filename,dpi=1200)
        plt.close()
        return
    else:
        return plt.show()

#### FUNCTION FOR DATABASE
def get_data(file):
    return os.path.join(_ROOT, 'data', file)