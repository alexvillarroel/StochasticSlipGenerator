import numpy as np
from scipy.io import loadmat,savemat
import geostochpy
import shutil
def save_matformat(lons_ep,lats_ep,depth,dip,rake,strike,slip,dx,dy,output_file,time=[]):
    # Save grids in a .mat file
    # if you have a kinematic rupture, you can to insert a time array of nx*ny size
    # else, all subfaults broke in t=0
    if time!=[]:
        mdic_multifault={'depth':depth.flat,'length':dy*np.ones_like(depth.flat),'width':dx*np.ones_like(depth.flat),
                'slip':slip,'strike':strike.flat,'dip':dip.flat,'rake':rake.flat
                ,'lat':lats_ep,'lon':lons_ep,'time':np.zeros_like(depth.flat)}
        if output_file.endswith('.mat'):
            savemat(output_file,mdic_multifault)
        else:
            savemat(output_file+'.mat',mdic_multifault)
    else:
        mdic_multifault={'depth':depth.flat,'length':dy*np.ones_like(depth.flat),'width':dx*np.ones_like(depth.flat),
                'slip':slip,'strike':strike.flat,'dip':dip.flat,'rake':rake.flat
                ,'lat':lats_ep,'lon':lons_ep,'time':time}
        if output_file.endswith('.mat'):
            savemat(output_file,mdic_multifault)
        else:
            savemat(output_file+'.mat',mdic_multifault)

    return print(f'{output_file} saved :) \n')


def make_multifault(filemat,fname,fname_layerdomain):
    # Load the fault parameters from the MATLAB file
    layerdomain=np.loadtxt(fname_layerdomain)
    X_start=np.min(layerdomain[:,0])
    Y_start=np.min(layerdomain[:,1])
    data = loadmat(filemat)
    depth = data['depth'].flatten()
    length = data['length'].flatten()
    width = data['width'].flatten()
    slip = data['slip'].flatten()
    strike = data['strike'].flatten()
    dip = data['dip'].flatten()
    rake = data['rake'].flatten()
    lat = data['lat'].flatten()
    lon = data['lon'].flatten()
    time = data['time'].flatten()
    n_segment = depth.size
    # You need to remember that this file write since second subfault,
    # the first you need insert in COMCOT 
    # Open the COMCOT control file for writing
    with open(fname, 'w') as f:
        # Write the COMCOT control file header
        f.write('#################################################\n')
        f.write('#                                               #\n')
        f.write('# Control file for COMCOT program (v1.7)        #\n')
        f.write('#       MULTI-FAULT CONFIGURATION               #\n')
        f.write('#################################################\n')
        f.write('#--+-----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8\n')
        f.write('#===============================================:===============================\n')
        f.write('# Parameters for Multiple Fault Plane Setup     : Value Field                  |\n')
        f.write('#===============================================:===============================\n\n')
        # Write the fault parameters for each segment
        for i in range(1, n_segment):
            f.write('###############################################:===============================\n')
            f.write('# Parameters for Fault Segment %i               :Values                        |\n' % int(i+1))
            f.write('###############################################:===============================\n')
            f.write(' Fault Rupture Time (seconds)                   :   %f\n' % time[i])
            f.write(' Faulting Option (0: Model; 1- Data;)           :   0\n')
            f.write(' Focal Depth                             (meter):     %f\n' % (depth[i] ))
            f.write(' Length of source area                   (meter):    %f\n' % (length[i] * 1000))
            f.write(' Width of source area                    (meter):   %f\n' % (width[i] * 1000))
            f.write(' Dislocation of fault plate              (meter):     %f\n' % (slip[i] ))
            f.write(' Strike direction (theta)               (degree):   %f\n' % strike[i])
            f.write(' Dip  angle       (delta)               (degree):    %f\n' % dip[i])
            f.write(' Slip angle       (lambda)              (degree):   %f\n' % rake[i])
            f.write(' Origin of Comp. Domain (Layer 01) (Lat, degree):   %f\n' % Y_start)
            f.write(' Origin of Comp. Domain (Layer 01) (Lon, degree):   %f\n' % X_start)
            f.write(' Epicenter Location: Latitude           (degree):   %f\n' % lat[i])
            f.write(' Epicenter Location: Longitude          (degree):   %f\n' % (lon[i]))
            f.write(' File Name of Deformation Data                  : raukumara_gps_coupling.xyz\n')
            f.write(' Data Format Option (0-COMCOT; 1-MOST; 2-XYZ)   :     2\n\n')

    # Close the COMCOT control file
    f.close()
    list_segment1=[time[0],depth[0],length[0]*1000,width[0]*1000,slip[0],strike[0],dip[0],rake[0],Y_start,X_start,lat[0],lon[0]]
    return list_segment1
    # MAKE COMCOT.CTL FILE
def make_comcot_ctl(total_runtime,time_save_interval,num_subfaults,list_segment1,filename_layers):
    # Especifica el nombre del archivo
    route_comcot=geostochpy.get_data('comcot.ctl')

    # Abre el archivo en modo escritura
    with open(route_comcot, 'r') as f:
        # Escribe el encabezado del archivo de control COMCOT
        lines=f.readlines()
    ##### read filename_layers
    gridsizes=[]
    #
    X_start=[]
    X_end=[]
    Y_Start=[]
    Y_end=[]
    #
    for layer in filename_layers:
        file=np.loadtxt(layer)
        if gridsizes==[]:
            gridsizes.append(np.abs(file[1,0]-file[0,0])*60)
        else:
            gridsizes.append(gridsizes[0]/(np.abs(file[1,0]-file[0,0])*60))
        X_start.append(np.nanmin(file[:,0]))
        X_end.append(np.nanmax(file[:,0]))
        Y_Start.append(np.nanmin(file[:,1]))
        Y_end.append(np.nanmax(file[:,1]))
    ###### Start with edition
    lines[10]=f' Total run time (Wall clock, seconds)           : {total_runtime}\n'
    lines[11]=f' Time interval to Save Data    ( unit: sec )    : {time_save_interval}\n'
    #### editing fault parameters#####
    lines[25]=f' No. of FLT Planes (With fault_multi.ctl if >1) : {num_subfaults}\n'
    lines[26]=f' Fault Rupture Time (seconds)                   : {list_segment1[0]}\n'
    lines[28]=f' Focal Depth                             (meter): {list_segment1[1]}\n'
    lines[29]=f' Length of source area                   (meter): {list_segment1[2]}\n'
    lines[30]=f' Width of source area                    (meter): {list_segment1[3]}\n'
    lines[31]=f' Dislocation of fault plate              (meter): {list_segment1[4]}\n'
    lines[32]=f' Strike direction (theta)               (degree): {list_segment1[5]}\n'
    lines[33]=f' Dip  angle       (delta)               (degree): {list_segment1[6]}\n'
    lines[34]=f' Slip angle       (lamda)               (degree): {list_segment1[7]}\n'
    lines[35]=f' Origin of Comp. Domain (Layer 01) (Lat, degree): {list_segment1[8]}\n'
    lines[36]=f' Origin of Comp. Domain (Layer 01) (Lon, degree): {list_segment1[9]}\n'
    lines[37]=f' Epicenter: Latitude                    (degree): {list_segment1[10]}\n'
    lines[38]=f' Epicenter: Longitude                   (degree): {list_segment1[11]}\n'
    ###### editing grid 1 parameters
    lines[66]=f' Run This Layer ?       (0:Yes,       1:No     ):     0\n'
    lines[67]=f' Coordinate System    (0:spherical, 1:cartesian):     0\n'
    lines[68]=f' Governing Equations  (0:linear,    1:nonlinear):     0\n'
    lines[69]=f' Grid Size  (dx, sph:minute, cart:meter)        :     {gridsizes[0]}\n'
    lines[70]=f' Time step                            ( second ):     1.0\n'
    lines[71]=f' Bottom Friction Switch? (0:Yes,1:No,2:var. n ) :     1\n'
    lines[72]=f' Manning\'s Roughness Coef. (For fric.option=0)  :     0.013\n'
    lines[73]=f' Layer Ouput Option? (0:Z+Hu+Hv;1:Z Only;2:NONE):     1\n'
    lines[74]=f' X_start                                        :  {X_start[0]}\n'
    lines[75]=f' X_end                                          :  {X_end[0]}\n'
    lines[76]=f' Y_Start                                        :  {Y_Start[0]}\n'
    lines[77]=f' Y_end                                          :  {Y_end[0]}\n'
    lines[78]=f' File Name of Bathymetry Data                   : {filename_layers[0]}\n'
    # editin another grids
    for i in np.arange(1,len(gridsizes)):
        k=0
        lines[87+20*k]=f' Run This Layer ?       (0:Yes,       1:No     ):     0\n'
        lines[88+20*k]=f' Coordinate           (0:spherical, 1:cartesian):     0\n'
        lines[89+20*k]=f' Governing Eqn.       (0:linear,    1:nonlinear):     0\n'
        lines[90+20*k]=f' Bottom Friction Switch? (0:Yes,1:No,2:var. n ) :     1\n'
        lines[91+20*k]=f' Manning\'s Roughness Coef. (For fric.option=0)  :     0.013\n'
        lines[92+20*k]=f' Layer Ouput Option? (0:Z+Hu+Hv;1:Z Only;2:NONE):     1\n'
        lines[93+20*k]=f' GridSize Ratio of Parent layer to current layer:     {gridsizes[i]}\n'
        lines[94+20*k]=f' X_start                                        :  {X_start[i]}\n'
        lines[95+20*k]=f' X_end                                          :  {X_end[i]}\n'
        lines[96+20*k]=f' Y_Start                                        :  {Y_Start[i]}\n'
        lines[97+20*k]=f' Y_end                                          :  {Y_end[i]}\n'
        lines[98+20*k]=f' File Name of Bathymetry Data                   : {filename_layers[i]}\n'
        k+=1
    with open('comcot.ctl', 'w') as newfile:
        newfile.writelines(lines)
        return

