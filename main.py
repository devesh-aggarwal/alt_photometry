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

    # # find stars
    daofind = DAOStarFinder(fwhm=30.0, threshold=400.*std)
    sources = daofind(data - median)
    print(f'filename {filename} with flux {sources['flux'][0]}')
    flux_values.append(sources['flux'][0])

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

# visualize the flux values
plt.plot(flux_values)
plt.xlabel('File index')
plt.ylabel('Flux value')
plt.title('Flux values of detected stars')
plt.show()