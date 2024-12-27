import matplotlib.pyplot as plt 
import numpy as np
from astropy.stats import sigma_clipped_stats
from astropy.io import fits
import os
from astropy.visualization import SqrtStretch, simple_norm
from astropy.visualization.mpl_normalize import ImageNormalize
from photutils.aperture import CircularAperture
from photutils.detection import DAOStarFinder

directory = 'images/night2'
flux_values = []
mag_values = []

# sort directory in alphabetical order
file_list = sorted(os.listdir(directory))

for filename in file_list:
    # open file
    filepath = os.path.join(directory, filename)
    hdul = fits.open(filepath)
    data = hdul[0].data

    # compute stats
    mean, median, std = sigma_clipped_stats(data, sigma=3.0)
    #print(np.array((mean, median, std)))

    # find stars
    daofind = DAOStarFinder(fwhm=5.0, threshold=5.*std, brightest=1) # keeps only the brightest star
    sources = daofind(data - median)
    #sources.sort('flux', reverse=True)
    print(f'filename {filename} with flux {sources['flux'][0]}')
    flux_values.append(sources['flux'][0])
    mag_values.append(sources['mag'][0]) # calculated as -2.5 * log10(flux)

# UNCOMMENT TO VISUALIZE AND DEBUG INDIVIDAUL FILES
# filepath = os.path.join('images/night2/Light_Polaris_10.0s_IRCUT_20241127-202252.fit')
# hdul = fits.open(filepath)
# data = hdul[0].data

# # compute stats
# mean, median, std = sigma_clipped_stats(data, sigma=3.0)
# #print(np.array((mean, median, std)))

# # # find stars
# daofind = DAOStarFinder(fwhm=20.0, threshold=300.*std)
# sources = daofind(data - median)
# print(sources['flux'][0])
# flux_values.append(sources['flux'][0])

# # visualize and plot locations of detected sources
# positions = np.transpose((sources['xcentroid'], sources['ycentroid']))
# apertures = CircularAperture(positions, r=20.)
# norm = simple_norm(data, 'sqrt', percent=99.9)
# plt.imshow(data, cmap='Greys', origin='lower', norm=norm, interpolation='nearest')
# apertures.plot(color='blue', lw=1.5, alpha=0.5)
# plt.show()

# visualize the flux and magnitude values
plt.figure(figsize=(10, 5))

# Plot flux values
plt.subplot(1, 2, 1)
plt.plot(flux_values, label='Flux values')
plt.xlabel('File index')
plt.ylabel('Flux value')
plt.title('Flux values of detected stars')
# Add a trendline for flux values
z_flux = np.polyfit(range(len(flux_values)), flux_values, 2)
p_flux = np.poly1d(z_flux)
plt.plot(range(len(flux_values)), p_flux(range(len(flux_values))), "r--", label='Trendline')
plt.legend()

# Plot magnitude values
plt.subplot(1, 2, 2)
plt.plot(mag_values, label='Magnitude values', color='orange')
plt.xlabel('File index')
plt.ylabel('Magnitude value')
plt.title('Magnitude values of detected stars')
# Add a trendline for magnitude values
z_mag = np.polyfit(range(len(mag_values)), mag_values, 2)
p_mag = np.poly1d(z_mag)
plt.plot(range(len(mag_values)), p_mag(range(len(mag_values))), "r--", label='Trendline')
plt.legend()

plt.tight_layout()
plt.show()