import matplotlib.pyplot as plt 
import numpy as np
import os
from astropy.stats import sigma_clipped_stats
from astropy.io import fits
from photutils.detection import DAOStarFinder
from astropy.stats import sigma_clipped_stats
from astropy.visualization import simple_norm
from photutils.aperture import CircularAperture, aperture_photometry, CircularAnnulus, ApertureStats


directory = 'images/night1'
flux_values = []
mag_values = []

# sort directory in alphabetical order
file_list = sorted(os.listdir(directory))

for filename in file_list:
    # open file
    filepath = os.path.join(directory, filename)
    try:
        with fits.open(filepath) as hdul:
            data = hdul[0].data
    except OSError: # occurs when system files such as .DS_Store are encountered
        print(f"ERROR: {filename}")
        continue
    
    # compute stats to remove background (median)
    mean, median, std = sigma_clipped_stats(data, sigma=3.0)
    print(np.array((mean, median, std)))

    # find stars
    daofind = DAOStarFinder(fwhm=13, threshold=2500.*std)
    sources = daofind(data - median)

    # calculate flux using primary aperture
    positions = np.transpose((sources['xcentroid'], sources['ycentroid']))
    apertures = CircularAperture(positions, r=19.)
    phot_table = aperture_photometry(data, apertures)

    # subtract local background
    annulus_apertures = CircularAnnulus(positions, r_in=35, r_out=45)
    background_stats = ApertureStats(data, annulus_apertures)
    bkg_mean = background_stats.mean
    apertures_area = apertures.area_overlap(data)
    bkg_sum = bkg_mean * apertures_area
    phot_table['bkg_sum'] = bkg_sum
    phot_table['aperture_sum_bkgsub'] = phot_table['aperture_sum'] - bkg_sum
    phot_table.sort('aperture_sum_bkgsub', reverse=True)

    
    print(f'filename {filename} with true flux {phot_table['aperture_sum_bkgsub'][0]} (total flux {phot_table['aperture_sum'][0]} and background {phot_table['bkg_sum'][0]}).')
    flux_values.append(phot_table['aperture_sum_bkgsub'][0])


# calculate magnitude from flux
def calculate_magnitude(flux_values):
    return -2.5 * np.log10(flux_values)

mag_values = calculate_magnitude(np.array(flux_values))
delta_mag_values = mag_values - mag_values[0]

delta_flux_values = flux_values - flux_values[0]

# visualize the flux and magnitude values
plt.figure(figsize=(10, 5))

# Plot flux values
plt.subplot(1, 2, 1)
plt.plot(delta_flux_values, label='Flux values')
plt.xlabel('File index')
plt.ylabel('Flux value')
plt.title('Flux values of detected stars')

# Plot magnitude values
plt.subplot(1, 2, 2)
plt.plot(delta_mag_values, label='Magnitude values', color='orange')
plt.xlabel('File index')
plt.ylabel('Magnitude value')
plt.title('Magnitude values of detected stars')

plt.tight_layout()
plt.show()