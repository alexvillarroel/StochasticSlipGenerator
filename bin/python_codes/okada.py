import numpy as np
import matplotlib.pyplot as plt
from clawpack.geoclaw import dtopotools

def set_fault_vectorized(strike, dips, rakes, depths, slips, lengths, widths, longitudes, latitudes):
    if not (len(dips) == len(rakes) == len(depths) == len(slips) == len(lengths) == len(widths) == len(longitudes) == len(latitudes)):
        raise ValueError("All input arrays must have the same length")

    subfaults = []
    for dip, rake, depth, slip, length, width, longitude, latitude in zip(dips, rakes, depths, slips, lengths, widths, longitudes, latitudes):
        subfault = dtopotools.SubFault()
        subfault.strike = strike
        subfault.dip = dip
        subfault.rake = rake
        subfault.depth = depth
        subfault.slip = slip
        subfault.length = length
        subfault.width = width
        subfault.longitude = longitude
        subfault.latitude = latitude
        subfault.coordinate_specification = "top center"
        subfaults.append(subfault)

    fault = dtopotools.Fault()
    fault.subfaults = subfaults
    return fault

def plot_okada(strike, dips, rakes, depths, slips, lengths, widths, longitudes, latitudes, verbose=False):
    """
    Make 3 plots to illustrate the Okada solution.
    """
    fault = set_fault_vectorized(strike, dips, rakes, depths, slips, lengths, widths, longitudes, latitudes)
    
    fig, axs = plt.subplots(2, 2, figsize=(12, 10))
    ax1, ax2, ax3, ax4 = axs.flatten()

    # Subfault projection on surface on ax1:
    ax = fault.plot_subfaults(axes=ax1, plot_rake=True, xylim=[-.5, 1.5, -1, 1])
    for i, (dip, rake, depth) in enumerate(zip(dips, rakes, depths)):
        ax1.text(0.6, 0.8 - 0.02 * i, f"Subfault {i + 1}: Strike = {strike:.1f}, Dip = {dip:.1f}, Rake = {rake:.1f}, Depth = {depth / 1e3:.1f} km", fontsize=6)
    ax1.set_ylabel('latitude (degrees)')

    # Depth profile on ax3:
    for subfault in fault.subfaults:
        z_top = -subfault.centers[0][2] / 1e3  # convert to km
        z_bottom = -subfault.centers[2][2] / 1e3  # convert to km
        ax3.plot([0, np.cos(subfault.dip * np.pi / 180.) * subfault.width / 1.e3], [z_top, z_bottom])
    ax3.set_xlim(-50, 150)
    ax3.set_ylim(-55, 0)
    ax3.set_xlabel('distance orthogonal to strike')
    ax3.set_ylabel('depth (km)')
    ax3.set_title('Depth profile')

    # Grid to use for evaluating and plotting dz
    x = np.linspace(-0.5, 1., 101)
    y = np.linspace(-1., 1., 101)
    times = [1.]

    # color map of deformation dz on ax2:
    fault.create_dtopography(x, y, times, verbose=verbose)
    dtopo = fault.dtopo
    dtopo.plot_dZ_colors(t=1., axes=ax2)

    # transect of dz on ax4:
    dZ = dtopo.dZ[-1, 50, :]
    ax4.plot(x, dZ)
    ax4.set_ylim(-0.5, 0.5)
    ax4.set_title('Transect of dz along y=0')
    ax4.set_xlabel('Longitude (degrees)')
    ax4.set_ylabel('Seafloor deformation (m)')

    plt.tight_layout()
    plt.show()

# Ejemplo con 500 subfallas
num_subfaults = 500
strike = 300
dips = np.random.uniform(10, 30, num_subfaults)  # Valores aleatorios entre 10 y 30 grados
rakes = np.random.uniform(80, 100, num_subfaults)  # Valores aleatorios entre 80 y 100 grados
depths = np.random.uniform(5e3, 25e3, num_subfaults)  # Valores aleatorios entre 5 y 25 km
slips = np.random.uniform(0.5, 2.0, num_subfaults)  # Valores aleatorios entre 0.5 y 2 metros
lengths = np.random.uniform(80e3, 150e3, num_subfaults)  # Valores aleatorios entre 80 y 150 km
widths = np.random.uniform(40e3, 70e3, num_subfaults)  # Valores aleatorios entre 40 y 70 km
longitudes = np.random.uniform(-0.5, 0.5, num_subfaults)  # Valores aleatorios entre -0.5 y 0.5 grados
latitudes = np.random.uniform(-0.5, 0.5, num_subfaults)  # Valores aleatorios entre -0.5 y 0.5 grados

plot_okada(strike, dips, rakes, depths, slips, lengths, widths, longitudes, latitudes)
