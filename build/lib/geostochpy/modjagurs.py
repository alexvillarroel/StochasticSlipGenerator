import numpy as np
import netCDF4 as nc
import os
import pygmt
from scipy.io import loadmat
import subprocess
def cell2node(Z):
    Z = np.array(Z, dtype=float)
    Nz = Z.shape
    aux = np.zeros((Nz[0]+1, Nz[1]+1))

    aux[0, 0] = Z[0, 0]
    aux[0, -1] = Z[0, -1]
    aux[-1, 0] = Z[-1, 0]
    aux[-1, -1] = Z[-1, -1]

    for i in range(1, Nz[1]-1):
        aux[0, i] = 0.5 * (Z[0, i] + Z[0, i+1])
        aux[-1, i] = 0.5 * (Z[-1, i] + Z[-1, i+1])

    for i in range(1, Nz[0]-1):
        aux[i, :] = 0.5 * (Z[i, 0] + Z[i+1, 0])
        aux[i, -1] = 0.5 * (Z[i, -1] + Z[i+1, -1])

    aux[1:-1, 1:-1] = 0.25 * (Z[:-1, :-1] + Z[1:, 1:] + Z[:-1, 1:] + Z[1:, :-1])

    return aux
def node2cell(X):
    X = np.array(X, dtype=float)
    Nx = X.shape
    Y = np.zeros((Nx[0]-1, Nx[1]-1))

    Y = 0.25 * (X[1:, 1:] + X[:-1, :-1] + X[:-1, 1:] + X[1:, :-1])

    return Y
def read_ffm_usgs(file):
    fid = open(file, 'r')
    if fid != -1:
        line = fid.readline().strip()
    else:
        print('error al buscar archivo')
        return

    sep = line.find('=')
    fault_segment = round(float(line[sep+1:]))

    x_fault = []
    y_fault = []
    z_fault = []

    x_corner = []
    y_corner = []
    z_corner = []

    plane_corners = []

    x_centroid = []
    y_centroid = []
    z_centroid = []
    slip = []
    rake = []
    dip = []
    strike = []
    t_rup = []
    t_rise = []
    t_fal = []
    M0 = []

    new_line = 0
    nx = []
    ny = []
    Dx = []
    Dy = []

    for i in range(fault_segment):
        line = fid.readline().strip()

        sep_eqs = line.find('=')
        sep_Dx = line.find('Dx')
        sep_km = line.find('km')
        sep_Dy = line.find('Dy')

        nx.append(round(float(line[sep_eqs+1:sep_Dx-1])))
        Dx.append(float(line[sep_eqs+1:sep_km-1]))
        ny.append(round(float(line[sep_eqs+1:sep_Dy-1])))
        Dy.append(float(line[sep_eqs+1:sep_km-1]))

        line = fid.readline().strip()
        sepA = line.find(' ')
        sepB = line.find('.')
        sepC = line.find(',')
        sepD = line.find(':')

        nx0 = round(float(line[sepA+1:sepC-1]))
        ny0 = round(float(line[sepC+1:sepB-1]))

        Lon_epi = float(line[sepD+1:sepD+2])
        Lat_epi = float(line[sepD+2:])

        boundaries = []
        for _ in range(5):
            boundary = fid.readline().strip().split()
            boundaries.append([float(boundary[0]), float(boundary[1]), float(boundary[2])])
        plane_corners.append(boundaries)

        data = []
        for _ in range(nx[i]*ny[i]):
            values = fid.readline().strip().split()
            data.append([float(values[0]), float(values[1]), float(values[2]), float(values[3]), float(values[4]), float(values[5]), float(values[6]), float(values[7]), float(values[8]), float(values[9]), float(values[10])])

        nd = ny[i]
        ns = nx[i]

        x_fault_i = [[0.0] * (ns+1) for _ in range(nd+1)]
        y_fault_i = [[0.0] * (ns+1) for _ in range(nd+1)]
        z_fault_i = [[0.0] * (ns+1) for _ in range(nd+1)]

        x_corner_i = [[0.0] * ns for _ in range(nd)]
        y_corner_i = [[0.0] * ns for _ in range(nd)]
        z_corner_i = [[0.0] * ns for _ in range(nd)]

        x_centroid_i = [[0.0] * ns for _ in range(nd)]
        y_centroid_i = [[0.0] * ns for _ in range(nd)]
        z_centroid_i = [[0.0] * ns for _ in range(nd)]

        slip_i = [[0.0] * ns for _ in range(nd)]
        rake_i = [[0.0] * ns for _ in range(nd)]
        strike_i = [[0.0] * ns for _ in range(nd)]
        dip_i = [[0.0] * ns for _ in range(nd)]
        t_rup_i = [[0.0] * ns for _ in range(nd)]
        t_rise_i = [[0.0] * ns for _ in range(nd)]
        t_fal_i = [[0.0] * ns for _ in range(nd)]
        M0_i = [[0.0] * ns for _ in range(nd)]

        for m in range(nd):
            for n in range(ns):
                k = n + ns*(m-1)
                x_centroid_i[m][n] = data[k][1]
                y_centroid_i[m][n] = data[k][0]
                z_centroid_i[m][n] = data[k][2]
                slip_i[m][n] = data[k][3] / 100
                rake_i[m][n] = data[k][4]
                strike_i[m][n] = data[k][5]
                dip_i[m][n] = data[k][6]
                t_rup_i[m][n] = data[k][7]
                t_rise_i[m][n] = data[k][8]
                t_fal_i[m][n] = data[k][9]
                M0_i[m][n] = data[k][10]

        sd = 0
        if data[0][0] > data[1][0]:
            sd = 1
        else:
            sd = -1

        vs = [data[1][1] - data[0][1], data[1][0] - data[0][0], data[1][2] - data[0][2]]
        vd = [data[ns+1][1] - data[0][1], data[ns+1][0] - data[0][0], data[ns+1][2] - data[0][2]]

        for m in range(nd+1):
            for n in range(ns+1):
                lambda_ = n - 1.5
                mu = m - 1.5
                x_fault_i[m][n] = x_centroid_i[0][0] + lambda_*vs[0] + mu*vd[0]
                y_fault_i[m][n] = y_centroid_i[0][0] + lambda_*vs[1] + mu*vd[1]
                z_fault_i[m][n] = z_centroid_i[0][0] + lambda_*vs[2] + mu*vd[2]

        for m in range(nd):
            for n in range(ns):
                lambda_ = n - 0.5
                mu = m - 1.5
                x_corner_i[m][n] = x_centroid_i[0][0] + lambda_*vs[0] + mu*vd[0]
                y_corner_i[m][n] = y_centroid_i[0][0] + lambda_*vs[1] + mu*vd[1]
                z_corner_i[m][n] = z_centroid_i[0][0] + lambda_*vs[2] + mu*vd[2]

        x_fault.append(x_fault_i)
        y_fault.append(y_fault_i)
        z_fault.append(z_fault_i)

        x_corner.append(x_corner_i)
        y_corner.append(y_corner_i)
        z_corner.append(z_corner_i)

        x_centroid.append(x_centroid_i)
        y_centroid.append(y_centroid_i)
        z_centroid.append(z_centroid_i)

        slip.append(slip_i)
        rake.append(rake_i)
        strike.append(strike_i)
        dip.append(dip_i)
        t_rup.append(t_rup_i)
        t_rise.append(t_rise_i)
        t_fal.append(t_fal_i)
        M0.append(M0_i)

        new_line = 1

    fid.close()

    S = {
        'planes_number': fault_segment,
        'x_fault': x_fault,
        'y_fault': y_fault,
        'z_fault': z_fault,
        'x_corner': x_corner,
        'y_corner': y_corner,
        'z_corner': z_corner,
        'x_centroid': x_centroid,
        'y_centroid': y_centroid,
        'z_centroid': z_centroid,
        'slip': slip,
        'rake': rake,
        'strike': strike,
        'dip': dip,
        't_rup': t_rup,
        't_rise': t_rise,
        't_fal': t_fal,
        'M0': M0,
        'Epicenter': [Lon_epi, Lat_epi],
        'cell_epicenter': [nx0, ny0],
        'nx': nx,
        'ny': ny,
        'Dx': Dx,
        'Dy': Dy,
        'plane_corners': plane_corners
    }

    return S
def read_jagurs(file, type):
    dataset = nc.Dataset(file, 'r')
    
    x_range = dataset.variables['x_range'][:]
    y_range = dataset.variables['y_range'][:]
    z = dataset.variables['z'][:]
    N = dataset.variables['dimension'][:]
    
    z = np.flipud(np.reshape(z, (N[0], N[1])).T)
    
    if type == 0:
        z = -z
    elif type == 1:
        z = np.flipud(z)
    
    varargout = [x_range, y_range, z]
    
    if len(varargout) == 4:
        spacing = dataset.variables['spacing'][:]
        varargout.append(spacing)
    
    dataset.close()
    
    return varargout
def generate_slip_file(file_input, directory, kinematic, DT, usgs_format, file_out):
    def cell2node(Z):
        # TODO: Implement cell2node function
        pass

    def node2cell(X):
        # TODO: Implement node2cell function
        pass

    def read_ffm_usgs(file):
        # TODO: Implement read_ffm_usgs function
        pass

    def read_jagurs(file, type):

        # TODO: Implement read_jagurs function
        pass

    def source_time(t):
        # TODO: Implement source_time function
        pass

    def generate_fault_step(file_input, DT, tv, tr, U, slip, n):
        fault_step = f"{file_input[:-4]}.{(n-1)*DT}-{n*DT}sec.txt"
        slp_rate = slip * (source_time((n*DT - tv) / tr) - source_time(((n-1)*DT - tv) / tr))
        DATA = np.column_stack((U, slp_rate))

        with open(fault_step, 'w') as file:
            file.write("! lat[deg] lon[deg] depth[km] length[km] width[km] dip[deg] strike[deg] rake [deg] slip_amp[m]\n")
            np.savetxt(file, DATA, fmt='%6.4f', delimiter=' ')

        return fault_step

    def generate_fault_list(file_input, Nt):
        file_list = f"{file_input[:-4]}.list"
        with open(file_list, 'w') as file:
            for n in range(1, Nt+1):
                fault_step = generate_fault_step(file_input, DT, tv, tr, U, slip, n)
                file.write(f"{fault_step}\n" if n < Nt else fault_step)

        return file_list

    def generate_static_file(file_out):
        DATA = np.column_stack((U, slip))

        with open(file_out, 'w') as file:
            file.write("! lat[deg] lon[deg] depth[km] length[km] width[km] dip[deg] strike[deg] rake [deg] slip_amp[m]\n")
            np.savetxt(file, DATA, fmt='%6.4f', delimiter=' ')

    if usgs_format == 1:
        B = read_ffm_usgs(file_input)
        A = np.zeros((np.dot(B.nx, B.ny), 12))
        Np = B.planes_number
        i = 0
        for n in range(Np):
            no = B.nx[n] * B.ny[n]
            j = i + no
            A[i:j, 1] = B.t_rup[n][:]
            A[i:j, 2] = B.t_rise[n][:]
            A[i:j, 3] = B.slip[n][:]
            A[i:j, 4] = B.x_centroid[n][:]
            A[i:j, 5] = B.y_centroid[n][:]
            A[i:j, 6] = np.abs(B.z_centroid[n][:])
            A[i:j, 7] = A[i:j, 6] * 0 + B.Dx[n]
            A[i:j, 8] = A[i:j, 6] * 0 + B.Dy[n]
            A[i:j, 9] = B.strike[n][:]
            A[i:j, 10] = B.dip[n][:]
            A[i:j, 11] = B.rake[n][:]
            i = j + 1

        A = A[A[:, 1].argsort()]
        A[:, 0] = np.arange(1, np.dot(B.nx, B.ny) + 1)

    else:
        I = np.genfromtxt(file_input)
        A = I[:, 1:]

    U = A[:, [5, 4, 6, 7, 8, 10, 9, 11]]
    slip = A[:, 3]

    if kinematic:
        T = np.max(A[:, 1]) + np.max(A[:, 2])
        Nt = int(np.ceil(T / DT))
        T = Nt * DT

        tv = A[:, 1]
        tr = A[:, 2]

        file_list = generate_fault_list(file_input, Nt)
        type(file_list)

    else:
        generate_static_file(file_out)
    return
def convert_grd_levels(directory, level_name):
    load_path = os.path.join(directory, 'levels.mat')
    levels = np.load(load_path)
    N_levels = len(levels)

    for n in range(N_levels):
        file = f"{level_name}{n}.grd"
        filename = os.path.join(directory, file)

        if os.path.exists(filename):
            os.remove(filename)

        with nc.Dataset(filename, 'w', format='NETCDF4_CLASSIC') as dataset:
            dataset.createDimension('side', 2)

            x_range = dataset.createVariable('x_range', 'f4', ('side',))
            x_range[:] = levels[n]['x']
            x_range.units = 'x'

            y_range = dataset.createVariable('y_range', 'f4', ('side',))
            y_range[:] = levels[n]['y']
            y_range.units = 'latitude [degrees_north]'

            z_range = dataset.createVariable('z_range', 'f4', ('side',))
            z_range[:] = [np.min(-levels[n]['level']), np.max(-levels[n]['level'])]
            z_range.units = 'z'

            spacing = dataset.createVariable('spacing', 'f4', ('side',))
            spacing[:] = levels[n]['ds'] + [0, 0]

            dimension = dataset.createVariable('dimension', 'i4', ('side',))
            dimension[:] = [levels[n]['level'].shape[1], levels[n]['level'].shape[0]]

            z = dataset.createVariable('z', 'f4', ('xyside',))
            z[:] = -np.flipud(levels[n]['level']).flatten()
            z.scale_factor = 1
            z.add_offset = 0
            z._FillValue = np.nan
            z.node_offset = 0

            dataset.title = f"level {n}"
            dataset.source = 'geostochpy_Package'

        print(f"Level {n} successfully converted")

    g = 9.81
    for n in range(N_levels):
        h = -np.min(levels[n]['level'])
        dx = levels[n]['ds'] * 1.11e5
        dt_CFL = dx / np.sqrt(2 * g * h)
        print(f"Level {n} | CFL cond. | Max. dt = {dt_CFL:.4f} s")
    return
def make_faultfile(lat, lon, depth, length, width, dip, strike, rake, slip, file):
    """
    This function is responsible for creating a fault file. The fault file typically contains 
    information about the geological fault, such as its location, orientation, and slip-rate. 
    The format its for jagurs software.
    """
    NEWLINE_SIZE_IN_BYTES = -1  # -2 on Windows?
    DATA = np.column_stack((lat, lon, depth, length, width, dip, strike, rake, slip))
    with open(file, 'w') as file:
        file.write("! lat[deg] lon[deg] depth[km] length[km] width[km] dip[deg] strike[deg] rake [deg] slip_amp[m]\n")
        np.savetxt(file, DATA, fmt='%6.4f', delimiter=' ')
    file.seek(0, os.SEEK_END) # Go to the end of the file.
    # Go backwards one byte from the end of the file.
    file.seek(file.tell() - NEWLINE_SIZE_IN_BYTES, os.SEEK_SET)
    file.truncate() # Truncate the file to this point.
    return

def make_faultfile_mat(matfile,filename):
    NEWLINE_SIZE_IN_BYTES = 1  # -2 on Windows?
    DATA=loadmat(matfile)
    lat=DATA['lat'].flatten()
    lon=DATA['lon'].flatten()
    depth=DATA['depth'].flatten()
    depth=depth/1000
    length=DATA['length'].flatten()
    width=DATA['width'].flatten()
    dip=DATA['dip'].flatten()
    strike=DATA['strike'].flatten()
    rake=DATA['rake'].flatten()
    slip=DATA['slip'].flatten()
    DATA = np.column_stack((lat, lon, depth, length, width, dip, strike, rake, slip))
    with open(filename, 'w') as file:
        file.write("! lat[deg] lon[deg] depth[km] length[km] width[km] dip[deg] strike[deg] rake [deg] slip_amp[m]\n")
        np.savetxt(file, DATA, fmt='%6.4f', delimiter=' ')
        file.seek(0, os.SEEK_END) # Go to the end of the file.
        # Go backwards one byte from the end of the file.
        file.seek(file.tell() - NEWLINE_SIZE_IN_BYTES, os.SEEK_SET)
        file.truncate() # Truncate the file to this point.
    return
def plot_grids(file_list):
    # the file_list is a list of files to plot, must be have the same size of levels and
    # bath is positive and land is negative
    fig = pygmt.Figure()
    #!gmt grdmath -1 $i  MUL = $j
    # list=['bathy.SD01.grd','bathy.SD02.grd','bathy.SD03.grd','bathy.SD04.grd','bathy.SD05.grd']
    lvl=5
    with fig.subplot(nrows=2, ncols=3, subsize=("20c", "30c"), frame="lrtb",autolabel=True, margins=False):
        for i in file_list:
            j=i.replace('.grd','_topo.grd')
            subprocess.run(f'gmt grdmath -1 {i}  MUL = {j}',shell=True)
            info=pygmt.grdinfo(j,per_column=True)
            region=info.split(' ')[0:4]
            fig.basemap(region=region,projection='M15c',frame=True,panel=True)
            fig.grdimage(grid=j, projection='M15c', frame=True,cmap='geo')
            fig.colorbar()
            fig.text(x=region[1],y=region[2],text=f'Bath lvl {lvl}',font='10p,Helvetica-Bold,black',justify='L')
            lvl-=1
    fig.show()