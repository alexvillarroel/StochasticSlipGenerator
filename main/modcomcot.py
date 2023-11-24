import numpy as np
from scipy.io import loadmat
def make_multifault(file,fname,lat_griddomain,lon_griddomain):
    # Load the fault parameters from the MATLAB file
    data = loadmat(file)
    depth = data['depth']
    length = data['length']
    width = data['width']
    slip = data['slip']
    strike = data['strike']
    dip = data['dip']
    rake = data['rake']
    lat = data['lat']
    lon = data['lon']
    time = data['time']
    n_segment = depth.size
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
            f.write(' Fault Rupture Time (seconds)                   :   %f\n' % time[0,i])
            f.write(' Faulting Option (0: Model; 1- Data;)           :   0\n')
            f.write(' Focal Depth                             (meter):     %f\n' % (depth[0,i] ))
            f.write(' Length of source area                   (meter):    %f\n' % (length[0,i] * 1000))
            f.write(' Width of source area                    (meter):   %f\n' % (width[0,i] * 1000))
            f.write(' Dislocation of fault plate              (meter):     %f\n' % (slip[0,i] ))
            f.write(' Strike direction (theta)               (degree):   %f\n' % strike[0,i])
            f.write(' Dip  angle       (delta)               (degree):    %f\n' % dip[0,i])
            f.write(' Slip angle       (lambda)              (degree):   %f\n' % rake[0,i])
            f.write(' Origin of Comp. Domain (Layer 01) (Lat, degree):   %f\n' % lat_griddomain)
            f.write(' Origin of Comp. Domain (Layer 01) (Lon, degree):   %f\n' % lon_griddomain)
            f.write(' Epicenter Location: Latitude           (degree):   %f\n' % lat[0,i])
            f.write(' Epicenter Location: Longitude          (degree):   %f\n' % (lon[0,i]))
            f.write(' File Name of Deformation Data                  : raukumara_gps_coupling.xyz\n')
            f.write(' Data Format Option (0-COMCOT; 1-MOST; 2-XYZ)   :     2\n\n')

    # Close the COMCOT control file
    f.close()
    return