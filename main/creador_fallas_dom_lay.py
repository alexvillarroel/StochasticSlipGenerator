# Fallas para estudio de resolucion del modelo
# datos de limites de latitud y longitud de Moreno et al 2009

import numpy as np
import modokada as mo 
import modfallas as mf 
import modrestricciones as mr
import matplotlib.pyplot as plt
import os
from scipy.interpolate import RegularGridInterpolator, interp1d, interp2d ,griddata
from geographiclib.geodesic import Geodesic
import json
import gc 

gc.collect()

"""
GUARDAR O NO DEFORMACIONES
"""
#okada = raw_input("Calcular deformaciones? [s/n] ")
okada = "s"

"""
RUTA DE SALIDA
"""
ruta_out = "../fallas3"

"""
DEFINICION FUNCIONES AUXILIARES
"""

# transfomar de Ms a Mw segun formula de kausel y ramirez 1997
def Mw_kauselramirez(ms):
    logm0kauselramirez = 1.5*ms+16.30
    Mw_calcs = 2.0/3.0*logm0kauselramirez-10.7 
    Mw_calcs = np.round(Mw_calcs,1) # un decimal
    return Mw_calcs


"""
DEFINICION DE INTERVALOS DEL ARBOL LOGICO
"""
N_arbol         = np.array([14,16,18,20])       # complejidad
Mw_arbol        = np.linspace(8.6,9.0,5)        # magnitud 
LN_arbol        = np.linspace(-39.0,-38.0,5)    # Limite norte
LS_arbol        = np.linspace(-44.5,-43.0,5)    # Limite sur
sflats_arbol    = np.array([80.])               # cantidad de subfallas en latitud (along strike) para setear razon de aspecto 
dist_fosa_arbol = np.array([0., 20., 40., 60.]) # distancias (en km) del limite oeste de la ruptura con respecto de la fosa
niterslips      = 100                           # cantidad de modelos de slip creados para cada combinacion de parametros

"""
CARGA DE DATOS GEOMETRICOS DEL SLAB
"""

# ruta del archivo de la fosa
ruta_fosa = "../Slab/"
# archivo fosa ( primera columna: longitudes, segunda columna: latitudes)
arc_fosa  = ruta_fosa + "SAM2.txt"
# carga de fosa usando funcion del modulo modfallas
lonfosa, latfosa  = mf.carga_fosa(arc_fosa)
# directorio con archivos del modelo SLAB2.0 Hayes, 2018
directorio = "../Slab/"
# archivo de profundidad
proffile   = directorio + "sam_slab2_dep_02.23.18.xyz" # nombre del archivo de prof
slabprof   = np.genfromtxt(proffile, delimiter = ",") # se lee el archivo a un array
# archivo de dip
dipfile    = directorio + "sam_slab2_dip_02.23.18.xyz" # nombre del archivo de dip
slabdip    = np.genfromtxt(dipfile, delimiter = ",") # se lee el archivo a un array
# archivo de strike
strfile    = directorio + "sam_slab2_str_02.23.18.xyz"
slabstrike = np.genfromtxt(strfile, delimiter = ",") # se lee el archivo a un array
# archivo de rakes
rakefile   = directorio + "sam_rake.xyz"
slabrake   = np.genfromtxt(rakefile, delimiter = ",")

"""
PREACONDICIONAMIENTO DE DATOS PARA CALCULOS POSTERIORES
"""

# las longitudes estan en formato 0 - 360, se cambian a -180 - 180
slabprof[:,0]   = slabprof[:,0]   - 360
slabdip[:,0]    = slabdip[:,0]    - 360
slabstrike[:,0] = slabstrike[:,0] - 360

# longitudes y latitudes del modelo
lonmod = slabprof[:,0]
latmod = slabprof[:,1]
# longitudes y latitudes unicas del modelo
lonunimod = np.unique(lonmod) # x
latunimod = np.unique(latmod) # y
# se grillan los arrays
lonmodgrid, latmodgrid = np.meshgrid(lonunimod, np.flip(latunimod)) # X, Y datos grillados para graficar modelo completo

# profundidades, dips y strikes del modelo
profmod   = slabprof[:,2]*-1000. # metros positivos hacia abajo
dipmod    = slabdip[:,2]
strikemod = slabstrike[:,2]
# se grillan los datos
profmod   = np.reshape(profmod,np.shape(lonmodgrid))
dipmod    = np.reshape(dipmod,np.shape(lonmodgrid))
strikemod = np.reshape(strikemod,np.shape(lonmodgrid))
# cambiar nan por 0
idx_nan_prof_mod   = np.isnan(profmod)
idx_nan_dip_mod    = np.isnan(dipmod)
idx_nan_strike_mod = np.isnan(strikemod)
profmod[idx_nan_prof_mod]     = 0 # profundidad grillada con 0 en lugar de nan
dipmod[idx_nan_dip_mod]       = 0 # dip grillado con 0 en lugar de nan
strikemod[idx_nan_strike_mod] = 0 # strike grillado con 0 en lugar de nan

"""
CICLO DE CREACION DE MAPAS DE SLIP
"""

# numero inicial de identificadores (por defecto debe ser 1, de ser requerido se puede aumentar para no sobreescribir archivos preexistentes)
nslips = 1
for ln in range(len(LN_arbol)):
    for ls in range(len(LS_arbol)):
        for sf in range(len(sflats_arbol)):
            for df in range(len(dist_fosa_arbol)):    
                # tamano subfalla
                # razon de aspecto de subfallas (largo/ancho)
                razon_aspecto        = 7./2.
                # cantidad de subfallas
                n_subfallas_lons      = 20. # numero de subfallas along dip, junto a sflats_arbol[sf] define la razon de aspecto
                n_subfallas_lats      =  sflats_arbol[sf] #np.floor(n_subfallas_lons*razon_aspecto)
                #np.floor(n_subfallas_lats/razon_aspecto)
                print(n_subfallas_lats, n_subfallas_lons)
                lats_subfallas      = np.flip(np.linspace(LS_arbol[ls], LN_arbol[ln], int(n_subfallas_lats)))

                # determinacion tamano subfallas
                largo_falla = mf.dist_sf_alt(-72., -72., LN_arbol[ln], LS_arbol[ls])
                largo_subfallas = largo_falla / n_subfallas_lats
                # subfallas cuadradas
                ancho_subfallas = largo_subfallas
                ancho_falla     = ancho_subfallas * n_subfallas_lons

                # interpolacion de longitud de la fosa para las latitudes limites
                interpolador_lons_fosa   = interp1d(latfosa, lonfosa)
                interpolador_lats_fosa   = interp1d(lonfosa, latfosa)
                # longitud de la fosa a la latitud de cada subfalla
                lons_fosa_para_subfallas = interpolador_lons_fosa(lats_subfallas) 
                # corrimiento de limite oeste de las subfallas (corrimiento hacia el este, alejamiento de la fosa) de acuerdo a la distancia con la fosa 
                for idxlon in range(len(lons_fosa_para_subfallas)):
                    lons_fosa_para_subfallas[idxlon] = Geodesic.WGS84.Direct(lats_subfallas[idxlon], lons_fosa_para_subfallas[idxlon], 90, dist_fosa_arbol[df]*1000)[ "lon2" ]

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
                #print("latitudes y longitudes de fallas creadas")
                #print(LS_arbol[ls], LN_arbol[ln])
                """
                ANADE ATRIBUTOS
                """
                

                """
                INTERPOLACION
                """

                # se crea un interpolador: se debe limitar el area total del modelo slab 2.0 al area de interes mas un buffer para ahorrar costos computacionales
                # se corta el modelo completo al area de interes

                latnorte = LN_arbol[ln]     + 1.
                latsur   = LS_arbol[ls]     - 1.
                loneste  = array_lons.max() + 1. 
                lonoeste = array_lons.min() - 1.

                idx_mascara_lons = np.argwhere( (lonmodgrid > lonoeste) & (lonmodgrid < loneste) )
                idx_mascara_lats = np.argwhere( (latmodgrid > latsur) & ( latmodgrid < latnorte) )


                def coincidencias_filas(A,B):
                    coincidencias  =  [i for i in range(B.shape[0]) if np.any(np.all(A==B[i],axis=1))]
                    if len(coincidencias)==0:
                        return B[coincidencias]
                    return np.unique(B[coincidencias],axis=0)

                mascara = coincidencias_filas(idx_mascara_lons, idx_mascara_lats)
                # se cuenta la cantidad de columnas y filas con las que cuenta el rectangulo recortado
                primer_elemento_mascara = mascara[0,0] # se cuenta cuants veces se repite este elemento (que corresponde al indice de la primera fila del rectangulo de interes) para saber cuantas columnas tiene el rectangulo
                n_columnas_area_interes = np.shape(mascara[mascara[:,0]==primer_elemento_mascara,:])[0]
                n_filas_area_interes    = np.shape(mascara)[0]/n_columnas_area_interes

                # dimensiones de la mascara
                dim_mascara          = np.shape(mascara)
                dim_mascara_filas    = dim_mascara[0]
                dim_mascara_columnas = dim_mascara[1]

                # longitud cortada al area de interes
                lonmodgrid_cortado = np.ones((dim_mascara_filas,1))
                for i in range(dim_mascara_filas):
                    lonmodgrid_cortado[i] = lonmodgrid[mascara[i,0],mascara[i,1]]

                # latitud cortadad al area de interes
                latmodgrid_cortado = np.ones((dim_mascara_filas,1))
                for j in range(dim_mascara_filas):
                    latmodgrid_cortado[j] = latmodgrid[mascara[j,0],mascara[j,1]]

                # se crea la grilla cortada
                lonmodgrid_cortado_grilla, latmodgrid_cortado_grilla = np.meshgrid(lonmodgrid_cortado, latmodgrid_cortado)

                # profundidad del modelo cortado al area de interes
                profmod_cortado = np.ones((dim_mascara_filas,1))
                for h in range(dim_mascara_filas):
                    profmod_cortado[h] = profmod[mascara[h,0],mascara[h,1]]

                # se amolda a las dimensiones correctas
                profmod_cortado = np.reshape(profmod_cortado,(n_filas_area_interes, n_columnas_area_interes))

                # dip del modelo cortado al area de interes
                dipmod_cortado = np.ones((dim_mascara_filas,1))
                for d in range(dim_mascara_filas):
                    dipmod_cortado[d] = dipmod[mascara[d,0],mascara[d,1]]

                # se amolda a las dimensiones correctas
                dipmod_cortado = np.reshape(dipmod_cortado,(n_filas_area_interes, n_columnas_area_interes))

                # strike del modelo cortado al area de interes
                strikemod_cortado = np.ones((dim_mascara_filas,1))
                for s in range(dim_mascara_filas):
                    strikemod_cortado[s] = strikemod[mascara[s,0],mascara[s,1]]

                # se amolda a las dimensiones correctas
                strikemod_cortado = np.reshape(strikemod_cortado, (n_filas_area_interes, n_columnas_area_interes))

                # grilla del area de interes
                lonmodgrid_cortado_area = np.reshape(lonmodgrid_cortado, (n_filas_area_interes, n_columnas_area_interes))
                latmodgrid_cortado_area = np.reshape(latmodgrid_cortado, (n_filas_area_interes, n_columnas_area_interes))


                """
                INTERPOLADORES
                """

                # interpolador de profundidad
                interpolador_prof   = interp2d(lonmodgrid_cortado_area[0,:], latmodgrid_cortado_area[:,0], profmod_cortado, kind='cubic')

                # interpolador de dip
                interpolador_dip    = interp2d(lonmodgrid_cortado_area[0,:], latmodgrid_cortado_area[:,0], dipmod_cortado, kind='cubic')

                # interpolador de strike
                interpolador_strike = interp2d(lonmodgrid_cortado_area[0,:], latmodgrid_cortado_area[:,0], strikemod_cortado, kind='cubic')

                # interpolador de rake
                interpolador_rake   = interp2d(slabrake[:,0], slabrake[:,1], slabrake[:,2], kind = 'cubic')

                ###################
                # interpolaciones #
                ###################

                # profundidad
                prof_subfallas_int  = np.ones(np.shape(array_lons))
                for i in range(np.shape(array_lons)[0]):
                    for j in range(np.shape(array_lons)[1]):
                        prof_subfallas_int[i,j]   = interpolador_prof(array_lons[i,j], array_lats[i,j])
                prof_subfallas_int = np.abs(prof_subfallas_int)
                # dip
                dip_subfallas_int   = np.ones(np.shape(array_lons))
                for i in range(np.shape(array_lons)[0]):
                    for j in range(np.shape(array_lons)[1]):
                        dip_subfallas_int[i,j]    = interpolador_dip(array_lons[i,j], array_lats[i,j])

                # strike
                strike_subfallas_int = np.ones(np.shape(array_lons))
                for i in range(np.shape(array_lons)[0]):
                    for j in range(np.shape(array_lons)[1]):
                        strike_subfallas_int[i,j] = interpolador_strike(array_lons[i,j], array_lats[i,j])

                # rake
                rake_subfallas_int   = np.ones(np.shape(array_lons))
                for i in range(np.shape(array_lons)[0]):
                    for j in range(np.shape(array_lons)[1]):
                        rake_subfallas_int[i,j] = interpolador_rake(array_lons[i,j], array_lats[i,j])

                assert not np.isnan(np.sum(prof_subfallas_int)), "Error de interpolacion"
                print("interpolacion completa")

                
                """
                CALCULO DE SLIP
                """
                # matriz de medias
                mu   = mf.matriz_medias(11, prof_subfallas_int)
                # matriz de covarianza
                C    = mf.matriz_covarianza(dip_subfallas_int, prof_subfallas_int, array_lons, array_lats)
                if np.isnan(C).any() or np.isinf(C).any(): # se chequea que este correcta C
                    continue
                # creacion de mapas de slip    
                for comp in range(len(N_arbol)):
                    for mag in range(len(Mw_arbol)):
                        for i in range(niterslips):        
                            print(nslips, "ln: ", LN_arbol[ln], "ls: ", LS_arbol[ls], "Mw: ", Mw_arbol[mag] , "Complejidad: ", N_arbol[comp], "Dist. Fosa: ", dist_fosa_arbol[df])
                            Slip    = mf.distribucion_slip(C, mu, int(N_arbol[comp])) # se crea la distribucion de slip
                            #Mw      = Mw_kauselramirez(Ms_arbol[mag])  #para fallas de la franche          # se calcula el Mw a partir de Ms con la formula de Kausel y Ramirez
                            # se crea un taper para modular slip en la fosa
                            Mw = Mw_arbol[mag]
                            ventana = mf.ventana_taper_slip_fosa(Slip, array_lons, array_lats,2) # ventana de taper
                            Slip    = mf.taper_slip_fosa(Slip,ventana)
                            Slip    = mf.escalar_magnitud_momento(Mw+0.1, Slip, prof_subfallas_int, array_lons, array_lats) # se escala el Slip a la magnitud deseada <--------- Slip final
                            # descomentar solo en fase de pruebas
                            #fig = plt.figure()
                            #ax = fig.add_subplot(111)
                            #im = ax.imshow(Slip)
                            #fig.colorbar(im)
                            #plt.show()
                            """
                            GUARDADO DE DATOS
                            """
                            # nombre archivo con datos de la falla
                            nombre_archivo = "falla_%d" %(nslips) 
                            archivo_salida = os.path.join(ruta_out, nombre_archivo)
                            # nombre diccionario con parametros
                            readme_nombre  = "readme_%d.json" %(nslips)
                            readme_salida  = os.path.join(ruta_out, readme_nombre)
                            # definicion diccionario con parametros del arbol
                            dict_ramas     = {"N": N_arbol[comp],
                                            "Mw": Mw,
                                            "LN": LN_arbol[ln],
                                            "LS": LS_arbol[ls],
                                            "AR": np.floor(n_subfallas_lats/n_subfallas_lons) } 
                            json_ramas   = json.dumps(dict_ramas)      # definicion diccionario json
                            json_archivo = open(readme_salida,'w')     # Creacion del archivo donde se guardara el diccionario
                            json_archivo.write(json_ramas)             # Escritura diccionario en archivo json 
                            json_archivo.close()                       # cierre archivo 
                            np.savez(archivo_salida,Slip=Slip,array_lons=array_lons,array_lats=array_lats,
                                                        dip_subfallas_int=dip_subfallas_int,strike_subfallas_int=strike_subfallas_int,
                                                        rake_subfallas_int=rake_subfallas_int, prof_subfallas_int=prof_subfallas_int,largo_falla=largo_falla) # se guarda los arrays de la falla
                            if okada == "s" or okada == "S":
                                falla = mo.crea_falla_dtopo( array_lons, array_lats, razon_aspecto, strike_subfallas_int, dip_subfallas_int, 
                                                            prof_subfallas_int, rake_subfallas_int, Slip, largo_falla)                  # creacion objeto falla
                                dtopo = mo.okada_solucion( array_lons, array_lats, razon_aspecto, strike_subfallas_int, dip_subfallas_int, prof_subfallas_int, rake_subfallas_int,
                                                            Slip, largo_falla, resolucion = 1/30., tamano_buffer = 1., verbose = True ) # calculo deformacion
                                #falla.plot_subfaults( slip_color = True )
                                # nombre archivo deformacion
                                if not np.isnan(dtopo.dZ_max()):
                                    nombre_archivo_def = "deformacion_%d.tt3" % (nslips)
                                    nombre_out_def = os.path.join(ruta_out, nombre_archivo_def)
                                    # guardado de archivo de deformacion
                                    dtopo.write(nombre_out_def, dtopo_type = 3)
                                else:
                                    continue        
                            nslips += 1
                    
plt.show()