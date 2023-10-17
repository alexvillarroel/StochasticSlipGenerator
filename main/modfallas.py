"""
Modulo con funciones necesarias para la creacion de fallas y el calculo de
distibuciones de slip
"""

"""
VERSION ORIGINAL CON PRUEBAS
"""

import numpy as np 
from numpy import linalg as la
import matplotlib.pyplot as plt 
import collections as col
from geopy import distance
import geographiclib as geo 
from geographiclib.geodesic import Geodesic
from scipy.interpolate import RegularGridInterpolator, interp1d
from scipy import signal
from scipy.signal import lfilter
import scipy.ndimage.filters as filters
from clawpack.geoclaw import dtopotools, topotools
from multiprocessing import Pool, Process, cpu_count



# funcion para cargar datos de modelo slab2.0
def carga_datos_slab2( directorio, latini, latfin, cambiarnan = True ):

    """
    Carga los datos del modelo slab2.0 de hayes 2018
    Entradas:
    directorio: directorio donde esta el archivo del modelo slab2
    cambiarnan: flag para reemplazar nan por 0, por defecto es falsa
    """

    # se chequea si el input directorio es o no string
    if not isinstance( directorio, basestring ) :
        directorio = str(directorio)
    else:
        directorio = directorio

    # se chequea el formateo del string (se desea que no termine con /)
    if not directorio.endswith("/"):
        directorio = directorio + "/"


    # archivo de profundidad
    proffile = directorio + "sam_slab2_dep_02.23.18.xyz" # nombre del archivo de prof
    slabprof = np.genfromtxt(proffile, delimiter = ",") # se lee el archivo a un array
    # archivo de dip
    dipfile = directorio + "sam_slab2_dip_02.23.18.xyz" # nombre del archivo de dip
    slabdip = np.genfromtxt(dipfile, delimiter = ",") # se lee el archivo a un array
    # archivo de strike
    strfile = directorio + "sam_slab2_str_02.23.18.xyz"
    slabstrike = np.genfromtxt(strfile, delimiter = ",") # se lee el archivo a un array

    # las longitudes estan en formato 0 - 360, se cambian a -180 - 180
    slabprof[:,0] = slabprof[:,0] - 360
    slabdip[:,0] = slabdip[:,0] - 360
    slabstrike[:,0] = slabstrike[:,0] - 360

    # se cambia dimensiones de los array para graficar
    repslat = col.Counter( slabprof[:,1] ).values( )[0] # (n cols) formato xyz repite valores de latitud para cada longitud, se obtiene cuantas veces se repite este valor para reshape
    repslon = len( slabprof )/repslat # (n filas)

    lon = np.reshape( slabprof[:,0], ( repslon, repslat ) )
    lat = np.reshape( slabprof[:,1], ( repslon, repslat ) )
    prof = np.reshape( slabprof[:,2], ( repslon, repslat ) ) * -1
    dip = np.reshape( slabdip[:,2], ( repslon, repslat ) )
    strike = np.reshape( slabstrike[:,2], ( repslon, repslat ) )

    idx = ( lat <= latini ) & ( lat >= latfin ) # indices de las latitudes dentro del area de interes
    # numero de columnas se mantiene (repslat), disminuye solo numero de filas (repslon)
    lon_adi = lon[idx] # adi: area de interes
    lat_adi = lat[idx]
    prof_adi = prof[idx]
    dip_adi = dip[idx]
    strike_adi = strike[idx]

    #idx_lonini = np.where(lon[0,] == lonini)[0][0]
    #idx_lonfin = np.where(lon[0,] == lonfin)[0][0]

    # redimensionar arrays
    filas = len(lat_adi)/repslat # cantidad de filas en array cortado nuevo
    lon = np.reshape( lon_adi, ( filas, repslat ) )
    lat = np.reshape( lat_adi, ( filas, repslat ) )
    prof = np.reshape( prof_adi, ( filas, repslat ) ) * 1000
    dip = np.reshape( dip_adi, ( filas, repslat ) )
    strike = np.reshape( strike_adi, ( filas, repslat ) )

    # si se desea se puede cambiar los valores nan por 0
    if cambiarnan:
        prof[ np.isnan( prof ) ] = 0
        dip[ np.isnan( dip ) ] = 0
        strike[ np.isnan( strike ) ] = 0

    # se debe revisar que la profundidad este en metros
    if prof.max( ) < 1000 :
        prof *= 1000

    return lon, lat, prof, dip, strike, repslat, repslon


# carga posicion de la fosa

def carga_fosa( archivo_fosa ):
    """
    carga archivo con la ubicacion de la fosa de la placa sudamericana
    Entrada: archivo fosa
    """

    fosa = np.genfromtxt( archivo_fosa, delimiter = " " )
    lons = fosa[:,0]
    lats = fosa[:,1]

    return lons, lats




# carga datos del modelo PREM

def carga_datos_PREM( directorio ):
    """
    carga los datos de velocidad de onda de corte del modelo prem
    directorio: directorio donde esta el archivo del modelo slab2
    """

    # se chequea si el input directorio es o no string
    if not isinstance( directorio, basestring):
        directorio = str( directorio )
    else:
        directorio = directorio

    # se chequea el formateo del string (se desea que no termine con /)
    if not directorio.endswith( "/" ):
        directorio = directorio + "/"

    premfile =  directorio + "PREM_1s.csv" # radio, profundidad, densidad, Vpv, Vph, Vsv, Vsh, eta, Q-mu, Q-kappa
    PREM = np.genfromtxt( premfile, delimiter = "," )
    prem_prof = PREM[:,1]*1000 # en metros
    prem_vs = PREM[:,5]*1000 # velocidad Vsv en metros/segundos

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
def crea_falla_fosa(lats, lons, prof, dip, strike, latini, latfin, area_sf, razon_aspecto):
    """
    Funcion que crea una falla calculando las latitudes y longitudes de diversas subfallas, considerando la fosa como un limite
    Entradas: 
    """
    # ruta del archivo de la fosa
    ruta_fosa = "../Slab/"
    # archivo fosa ( primera columna: longitudes, segunda columna: latitudes)
    arc_fosa  = ruta_fosa + "SAM2.txt"
    lonfosa, latfosa = carga_fosa(arc_fosa)
    # largo de falla
    largo_falla_completa = Geodesic.WGS84.Inverse(limnorte, -72, limsur, -72)[ "s12" ]
    # tamano subfallas (subfallas cuadradas)
    tamano_subfalla      = np.sqrt(area_sf)
    # cantidad de fallas
    cant_subfallas_lats  = largo_falla_completa // tamano_subfalla
    cant_subfallas_lons  = largo_falla_completa // razon_aspecto
    # latitudes subfallas
    lats_subfallas       = np.flip(np.linspace(limsur, limnorte, int(cant_subfallas_lats))) 

    # longitudes subfallas
    ancho_falla          = tamano_subfalla * cant_subfallas_lons
    # interpolacion de longitud de la fosa para las latitudes limites
    interpolador_lons_fosa   = interp1d(latfosa, lonfosa)
    interpolador_lats_fosa   = interp1d(lonfosa, latfosa)
    # longitud de la fosa a la latitud de cada subfalla
    lons_fosa_para_subfallas = interpolador_lons_fosa(lats_subfallas)
    # teniendo las longitudes de la fosa para las latitudes de las subfallas se tiene
    # el limite oeste de la falla, falta encontrar el limite este. Como se conoce el ancho de la falla, 
    # basta encontrar la longitud de este ancho para cada latitud
    lons_limite_este = np.ones(np.shape(lons_fosa_para_subfallas))

    for ilon in range(len(lons_limite_este)):
        lons_limite_este[ilon] = Geodesic.WGS84.Direct(lats_subfallas[ilon], lons_fosa_para_subfallas[ilon], 90, ancho_falla)[ "lon2" ]

    # teniendo los limites este y oeste, basta encontrar las longitudes de las subfallas intermedias

    array_lons = np.ones((int(n_subfallas_lats),int(n_subfallas_lons)))  # LONS SUBFALLAS
    for jlat in range(int(n_subfallas_lats)):
        array_lons[jlat,:] = np.linspace(lons_fosa_para_subfallas[jlat],lons_limite_este[jlat],int(n_subfallas_lons))

    array_lats = np.tile(np.reshape(lats_subfallas,(int(n_subfallas_lats),1)),int(n_subfallas_lons)) # LATS SUBFALLAS



# Estima el tiempo de ruptura (tau) de la falla

def tau_ruptura( largo_falla, beta = 2500 ):

    """
    Calcula el tiempo de ruptura de la falla dada a partir del largo de falla, considerando
    la velocidad de ruptura como 0.8 de la velocidad de propagacion de onda S y 
    una velocidad de cizalle de 2.5 km/s utilizando la aproximacion dada en Bilek y Lay, 1999
    tau approx L/(0.8beta) 
    Entradas: 
    largo_falla: largo de la falla en metros dado por la funcion que crea la falla
    beta: velocidad de la onda de cizalle en metros/segundo
    """

    # se chequean unidades de medida de beta, tiene que estar en metros/s, no km/s
    if beta < 1000:
        beta = beta * 1000
    else: 
        beta = beta

    tau = largo_falla/( 0.8 * beta )

    return tau

# Version alternativa de estimacion de tiempo de ruptura (tau) de la falla, utilizando la profundidad de la falla y 
# la de velocidad de onda de cizalle del modelo PREM
# suele sobreestimar los valores de rigidez usualmente presentes en la zona de subduccion 

def tau_ruptura_prem( largo_falla, prof_media, prof_prem, vs_prem ):

    """
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

    dif_profs = np.abs( prof_prem - prof_media )
    idx_prof_media = np.where( dif_profs == dif_profs.min() )
    vs_prof_media = vs_prem[idx_prof_media]

    # puede que la profundidad media de la falla se encuentre justo entre 2 valores del modelo
    # si esto sucede, se calcula el promedio de las velocidades de ambos
    if len(vs_prof_media) > 1:
        vs_prof_media = vs_prof_media.mean()
    else:
        vs_prof_media = vs_prof_media

    tau = largo_falla/( 0.8 * vs_prof_media )

    return tau



# Estima rigidez a partir de la ecuacion presentada en Bilek y Lay, 1999

def estima_rigidez( largo_falla, tau ):

    """
    Estima la rigidez de la falla dado el tiempo de ruptura y el largo de la falla
    segun la relacion dada en Bilek y Lay, 1999 mu = densidad*largo_falla**2/(0.8**2*tau**2).
    entradas: 
    largo_falla: largo de la falla en metros
    tau: tiempo de duracion de la ruptura
    """

    rho = 2700. # kg/m3
    mu = ( rho * ( largo_falla**2 ) )/( ( 0.8**2 ) * ( tau**2 ) )

    return mu

# estima rigidez a partir del modelo PREM

def estima_rigidez_prem( directorio, profs ):
    """
    Estima la rigidez de cada subfalla en funcion de su profunidad a partir de los 
    datos de velocidad de onda S y densdidad del modelo PREM, dada Vs=sqrt(mu/rho)
    Entradas:
    directorio: directorio donde se encuentra el archivo con el modelo PREM
    profs: profundidades de las subfallas
    """
        # se chequea si el input directorio es o no string
    if not isinstance( directorio, basestring):
        directorio = str( directorio )
    else:
        directorio = directorio

    # se chequea el formateo del string (se desea que no termine con /)
    if not directorio.endswith( "/" ):
        directorio = directorio + "/"

    premfile =  directorio + "PREM_1s.csv"
    PREM = np.genfromtxt( premfile, delimiter = "," )
    prem_vs = PREM[:,5] * 1000 # velocidad Vsv en metros/segundos
    prem_rho = PREM[:,2] * 1000 # densidad en m3/kg
    prem_prof = PREM[:,1] * 1000 # profundidad en metros
    rigidez = ( prem_vs * prem_vs ) * prem_rho # rigidez en Pa

    # se crea el interpolador 
    interp_rigidez = interp1d( prem_prof, rigidez )

    # se inicializa el array que contendra las rigideces de cada subfalla
    rigidez_sf = np.ones( np.shape( profs ) )
    for prow in range( np.shape( profs )[0] ):
        for pcol in range( np.shape( profs )[1] ):
            rigidez_sf[ prow ][ pcol ] = interp_rigidez( profs[ prow ][ pcol ] )

    return rigidez_sf


# calcula las distancias entre subfallas
# sin uso, utilizar version alternativa


def dist_sf( lon1, lon2, lat1, lat2 ):

    """
    Calcula la distancia en metros entre dos subfallas (sf) i y j utilizando el metodo de karney, 2013 y el elipsoide wgs84
    Entradas:
    lon1: longitud de subfalla i-esima
    lon2: longitud de subfalla j-esima
    lat1: latitud de subfalla i-esima
    lat2: latitud de subfalla j-esima
    """

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

def taper_LeVeque( d ):
    
    """
    realiza el taper propuesto por leveque et al., 2016
    t(d)=1-exp(-20|d-dmax|/dmax) 
    donde dmax es la profunidad maxima y de la falla y d es la profundidad
    de cada subfalla
    Entrada:
    d: matriz con profundidades
    """

    dmax = np.max( d )+500
    taper = ( 1 - np.exp( -20 * abs( d - dmax )/dmax ) )
    return taper

# calculo de medias

def matriz_medias( media, prof ):
    """
    calcula las medias mu_i para cada subfalla segun mu_i = log(mu*tau_i)-1/2log(alpha**2+1)
    Entradas:
    media: valor promedio alrededor del que se desea que se centre el slip
    prof: matriz con profundidades de cada subfalla, necesaria para el taper
    """
    tau = taper_LeVeque( prof ) # taper
    alpha = 0.5 # valor sugerido por LeVeque et al, 2016
    mu = np.log( media*tau ) - 1/2 * np.log( alpha**2+1 )
    return mu

# calculo matriz de covarianza

def matriz_covarianza( dip, prof, lons, lats ):

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
    """

    maxlat = np.max(lats) # latitud maxima
    lon_maxlat = lons[np.where( lats == maxlat )[0][0], np.where( lats == maxlat )[1][0]] # longitud correspondiente a la latitud maxima
    minlat = np.min(lats) # latitud minima
    lon_minlat = lons[np.where( lats == minlat )[0][0], np.where( lats == minlat )[1][0]] # longitud correspondiente a la latitud minima
    maxlon = np.max(lons) # longitud maxima
    lat_maxlon = lats[np.where( lons == maxlon )[0][0], np.where( lons == maxlon )[1][0]] # latitud correspondiente a la longitud maxima
    minlon = np.min(lons) # longitud minima
    lat_minlon = lats[np.where( lons == minlon )[0][0], np.where( lons == minlon )[1][0]]
    maxprof = np.max(prof) # profundidad maxima
    minprof = np.min(prof) # profundidad minima
    
    largo_falla = dist_sf_alt( lon_minlat, lon_maxlat, minlat, maxlat ) # largo de la falla en metros
    #ancho_falla = np.sqrt((np.max(prof)-np.min(prof))**2+(abs(np.max(lons)-np.min(lons))*111000.12)**2)
    ancho_falla = np.sqrt( ( np.max( prof ) - np.min( prof ) )**2 + dist_sf_alt( maxlon, minlon, lat_maxlon, lat_minlon )**2 )
    rdip = 0.4*ancho_falla # largo de correlacion en direccion de dip
    rstrike = 0.4*largo_falla # largo de correlacion en direccion de strike
    n_filas = np.shape( prof )[0] # dimension 0
    n_columnas = np.shape( prof )[1] # dimension 1
    vector_dip = np.reshape( dip, ( n_filas * n_columnas, 1 ) )
    vector_prof = np.reshape( prof, ( n_filas * n_columnas, 1 ) )
    vector_lon = np.reshape( lons, ( n_filas * n_columnas, 1 ) )
    vector_lat = np.reshape( lats, ( n_filas * n_columnas, 1 ) )
    alpha = 0.5 # LeVeque et al 2016

    # calculo de ddip
    ddip = np.ones( ( len( vector_dip ), len( vector_dip ) ) ) # iniciacion matriz de distancia a lo largo de dip
    for i in range( len( vector_dip ) ):
        for j in range( len( vector_dip ) ):
            if np.sin( np.deg2rad( ( vector_dip[i] + vector_dip[j] )/2 ) ) != 0:
                ddip[i][j] = ( vector_prof[i] - vector_prof[j] )/np.sin( np.deg2rad( ( vector_dip[i] + vector_dip[j] )/2 ) )
            else:
                ddip[i][j] = 0
    
    # calculo de dstrike
    dstrike = np.ones( ( len( vector_prof ), len( vector_prof ) ) ) # iniciacion matriz de distancia a lo largo de strike
    d = np.ones( ( len( vector_prof ), len( vector_prof ) ) ) # iniciacion matriz de distancia entre fallas
    for k in range( len( vector_prof ) ):
        for l in range( len( vector_prof ) ):
            #d[k][l] = dist_sf(vector_lat[k],vector_lat[l],vector_lon[k],vector_lon[l])
            d[k][l] = dist_sf_alt( vector_lon[k], vector_lon[l], vector_lat[k], vector_lat[l] )
            #d[k][l] = np.sqrt((abs(vector_lon[k]-vector_lon[l])*111000.12)**2+(abs(vector_lat[k]-vector_lat[l])*111000.12)**2)
            dstrike[k][l] = np.sqrt( abs(d[k][l]**2 - ( ddip[k][l] )**2 ) )
    
    # calculo de Cij
    C = np.ones( ( len( vector_dip ), len( vector_dip ) ) ) # matriz de correlacion
    for m in range( len( vector_prof ) ):
        for n in range( len( vector_prof ) ):
            C[m][n] = np.exp( -( ( dstrike[m][n] )/rstrike )-( ( ddip[m][n] )/rdip ) )


    mat_cova = np.log( alpha**2*C+1 ) # matriz de covarianza

    return mat_cova

# distribucion de slip

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
    eig_val, eig_vecs = la.eig( C )
    eig_val = eig_val.real # valores propios (lambda_k)
    eig_vecs = eig_vecs.real # vectores propios (v_k) (columnas de eig_vecs)

    z = np.random.normal( 0, 1, n_cols_cova ) # distribucion gaussiana aleatoria z~N(0,1)

    # iniciacion array de slip
    S = np.ones( ( n_cols_cova, N ) )
    for i in range(N):
      S[:,i] = z[i]*np.sqrt( np.abs(eig_val[i]) )*np.real(eig_vecs[:,i])
    S = np.multiply( mu, np.exp( np.sum( S, axis = 1 ) ) )
    S = np.reshape( S, dim_mu )
    return S


# calculo de magnitud

def magnitud_momento(slip, prof, lons, lats):
    """
    Calcula la magnitud de momento Mw de la distribucion de slip con la ecuacion
    Mw = 2/3log(Mo)-6.06
    donde Mo = mu * area de ruptura * slip
    Entradas: 
    slip: corresponde a la distribucion de slip 
    prof: matriz con profundidades
    lons: matriz con longitudes
    lats: matriz con latitudes
    """

    #mu = 3*10**10 # rigidez en N/m
    
    # calculos previos 

    maxlat     = np.max(lats) # latitud maxima
    lon_maxlat = lons[np.where(lats == maxlat)[0][0], np.where(lats == maxlat)[1][0]] # longitud correspondiente a la latitud maxima
    minlat     = np.min(lats) # latitud minima
    lon_minlat = lons[np.where(lats == minlat)[0][0], np.where(lats == minlat)[1][0]] # longitud correspondiente a la latitud minima
    maxlon     = np.max(lons) # longitud maxima
    lat_maxlon = lats[np.where(lons == maxlon)[0][0], np.where(lons == maxlon)[1][0]] # latitud correspondiente a la longitud maxima
    minlon     = np.min(lons) # longitud minima
    lat_minlon = lats[np.where(lons == minlon)[0][0], np.where(lons == minlon)[1][0]]
    maxprof    = np.max(prof) # profundidad maxima
    minprof    = np.min(prof) # profundidad minima
    
    largo_falla = dist_sf_alt(lon_minlat, lon_maxlat, minlat, maxlat) # largo de la falla en metros
    ancho_falla = np.sqrt((np.max(prof)-np.min(prof))**2 + dist_sf_alt( maxlon, minlon, lat_maxlon, lat_minlon)**2) # ancho de la falla

    n_fils = np.shape(lons)[0] # numero de filas
    n_cols = np.shape(lons)[1] # numero de columnas

    # tamano subfallas
    largo_subfalla = largo_falla/n_fils 
    ancho_subfalla = ancho_falla/n_cols 
    area_subfalla  = largo_subfalla*ancho_subfalla

    # se estima el tiempo que tarda la ruptura, tau, para luego estimar la rigidez
    # se utiliza beta = 2500 m/s
    tau = tau_ruptura( largo_falla )
    # se estima la rigidez de la interfaz dependiendo de la profundidad
    rigidez = estima_rigidez( largo_falla, tau )
    
    # calculo Mo

    Mo = rigidez*area_subfalla*slip
    Mo = np.sum(Mo)
    Mw = 2.0/3.0*np.log10(Mo)-6.06

    return Mw, area_subfalla

# estima la magnitud de momento utilizando la rigidez obtenida del modelo PREM

def magnitud_momento_PREM(slip, prof, lons, lats, directorio):
    """
    Calcula la magnitud de momento Mw de la distribucion de slip con la ecuacion
    Mw = 2/3log(Mo)-6.06. Calcula la rigidez individual de cada subfalla dependiendo de
    su profundidad, a partir de los datos de velocidad de onda de corte y densidad del 
    modelo PREM
    donde Mo = mu * area de ruptura * slip
    Entradas: 
    slip: corresponde a la distribucion de slip 
    prof: matriz con profundidades
    lons: matriz con longitudes
    lats: matriz con latitudes
    """

    # calculos previos 

    maxlat     = np.max(lats) # latitud maxima
    lon_maxlat = lons[np.where(lats == maxlat)[0][0], np.where(lats == maxlat)[1][0]] # longitud correspondiente a la latitud maxima
    minlat     = np.min(lats) # latitud minima
    lon_minlat = lons[np.where(lats == minlat)[0][0], np.where(lats == minlat)[1][0]] # longitud correspondiente a la latitud minima
    maxlon     = np.max(lons) # longitud maxima
    lat_maxlon = lats[np.where(lons == maxlon)[0][0], np.where(lons == maxlon)[1][0]] # latitud correspondiente a la longitud maxima
    minlon     = np.min(lons) # longitud minima
    lat_minlon = lats[np.where(lons == minlon)[0][0], np.where(lons == minlon)[1][0]]
    maxprof    = np.max(prof) # profundidad maxima
    minprof    = np.min(prof) # profundidad minima
    
    largo_falla = dist_sf_alt(lon_minlat, lon_maxlat, minlat, maxlat) # largo de la falla en metros
    ancho_falla = np.sqrt((np.max(prof)-np.min(prof))**2 + dist_sf_alt( maxlon, minlon, lat_maxlon, lat_minlon)**2) # ancho de la falla

    n_fils = np.shape(lons)[0] # numero de filas
    n_cols = np.shape(lons)[1] # numero de columnas

    # tamano subfallas
    largo_subfalla = largo_falla/n_fils 
    ancho_subfalla = ancho_falla/n_cols 
    area_subfalla  = largo_subfalla*ancho_subfalla

    rigidez = estima_rigidez_prem( directorio, prof )

    # calculo Mo y Mw

    Mo = np.multiply( rigidez, slip ) * area_subfalla
    Mo = np.sum(Mo)
    Mw = 2.0/3.0*np.log10(Mo)-6.06

    return Mw, area_subfalla





# escala la magnitud
def escalar_magnitud_momento(Mw, slip, prof, lons, lats):
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
    n_subfallas         = np.size( slip )
    slip_prim_homogeneo = np.ones( dims )

    # calculos previos 

    maxlat      = np.max(lats) # latitud maxima
    lon_maxlat  = lons[np.where(lats == maxlat)[0][0], np.where(lats == maxlat)[1][0]] # longitud correspondiente a la latitud maxima
    minlat      = np.min(lats) # latitud minima
    lon_minlat  = lons[np.where(lats == minlat)[0][0], np.where(lats == minlat)[1][0]] # longitud correspondiente a la latitud minima
    maxlon      = np.max(lons) # longitud maxima
    lat_maxlon  = lats[np.where(lons == maxlon)[0][0], np.where(lons == maxlon)[1][0]] # latitud correspondiente a la longitud maxima
    minlon      = np.min(lons) # longitud minima
    lat_minlon  = lats[np.where(lons == minlon)[0][0], np.where(lons == minlon)[1][0]]
    maxprof     = np.max(prof) # profundidad maxima
    minprof     = np.min(prof) # profundidad minima
    
    largo_falla = dist_sf_alt(lon_minlat, lon_maxlat, minlat, maxlat) # largo de la falla en metros
    ancho_falla = np.sqrt((np.max(prof)-np.min(prof))**2 + dist_sf_alt( maxlon, minlon, lat_maxlon, lat_minlon)**2) # ancho de la falla

    area_falla  = largo_falla * ancho_falla    
    # se utiliza beta = 2500 m/s
    tau = tau_ruptura( largo_falla )
    # se estima la rigidez de la interfaz dependiendo de la profundidad
    rigidez     = estima_rigidez( largo_falla, tau )


    Mw_prim, areasf = magnitud_momento( slip, prof, lons, lats ) # magnitud de momento de dist slip primitiva
    Mo_original     = 10**( (3./2.) * ( Mw_prim + 6.06 ) ) # momento de la distribucion original
    Mo_deseado      = 10**( (3./2.) * ( Mw + 6.06 ) )  # momento deseado para la magnitud a escalar
    razon_Mo        = Mo_deseado/Mo_original # razones entre los momentos, se cancelan los otros factores que no sean slip (rigidez y area)
    slip_escalado   = slip*razon_Mo
    #D_original = np.divide( ( Mo_original/n_subfallas ) * slip_prim_homogeneo, (rigidez*areasf)*np.ones( dims )   )


    #slip_prim_homogeneo = np.exp( ( 3.0/2.0 ) * ( Mw_prim+6.06 ) )/area_falla # array de slip homogeneo con Mw primitiva
    #slip_norm = np.divide( slip, slip_prim_homogeneo ) # array slip normalizado por slip con Mw primitiva
    #slip_deseado_homogeneo = np.exp( ( 3.0/2.0 ) * ( Mw + 6.06 ) )/area_falla # array de slip homogeneo con Mw deseado
    #slip_escalado = np.multiply(slip_norm, slip_deseado_homogeneo)
    return slip_escalado


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
        ventana     = signal.hamming(n_lons)
    elif ventana_flag == 5:
        n = np.linspace(1,100,n_lons)
        ventana     = sigmoid(n)

        
    ventana[int(n_lons/2):n_lons] = 1
    ventana = np.tile(ventana,(n_lats,1))

    return ventana


# crea el taper al patron de slip a partir de la ventana creada con la funcion ventana_taper_slip_fosa()
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


# funcion para calcular sta/lta
def sta_lta(data, dt, min_period):
    """
    STA/LTA as used in FLEXWIN.

    :param data: The data array.
    :param dt: The sample interval of the data.
    :param min_period: The minimum period of the data.
    """
    Cs = 10 ** (-dt / min_period)
    Cl = 10 ** (-dt / (12 * min_period))
    TOL = 1e-9

    noise = data.max() / 1E5

    # 1000 samples should be more then enough to "warm up" the STA/LTA.
    extended_syn = np.zeros(len(data) + 1000, dtype=np.float64)
    # copy the original synthetic into the extended array, right justified
    # and add the noise level.
    extended_syn += noise
    extended_syn[-len(data):] += data

    # This piece of codes "abuses" SciPy a bit by "constructing" an IIR
    # filter that does the same as the decaying sum and thus avoids the need to
    # write the loop in Python. The result is a speedup of up to 2 orders of
    # magnitude in common cases without needing to write the loop in C which
    # would have a big impact in the ease of installation of this package.
    # Other than that its quite a cool little trick.
    a = [1.0, -Cs]
    b = [1.0]
    sta = lfilter(b, a, extended_syn)

    a = [1.0, -Cl]
    b = [1.0]
    lta = lfilter(b, a, extended_syn)

    # STA is now STA_LTA
    sta /= lta

    # Apply threshold to avoid division by very small values.
    sta[lta < TOL] = noise
    return sta[-len(data):]