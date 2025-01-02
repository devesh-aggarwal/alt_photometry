import matplotlib
matplotlib.use('TkAgg') # better backend
import matplotlib.pyplot as plt 
import numpy as np
import os
from astropy.stats import sigma_clipped_stats
from astropy.io import fits
from photutils.detection import DAOStarFinder
from astropy.stats import sigma_clipped_stats
from astropy.visualization import simple_norm
from photutils.aperture import CircularAperture, aperture_photometry, CircularAnnulus, ApertureStats

filepath = os.path.join('images/night2_alt_separated/12005/Light_Polaris_10.0s_IRCUT_20241127-202314.fit')
hdul = fits.open(filepath)
data = hdul[0].data

# subtract background (median)
mean, median, std = sigma_clipped_stats(data, sigma=3.0)

# find stars
daofind = DAOStarFinder(fwhm=20, threshold=3000.*std)
sources = daofind(data - median)

# create apertures and perform photometry
positions = np.transpose((sources['xcentroid'], sources['ycentroid']))
apertures = CircularAperture(positions, r=30.)
phot_table = aperture_photometry(data, apertures)
print(phot_table)

# subtract local background
annulus_apertures = CircularAnnulus(positions, r_in=40, r_out=50)
background_stats = ApertureStats(data, annulus_apertures)
bkg_mean = background_stats.mean
print(bkg_mean)
apertures_area = apertures.area_overlap(data)
bkg_sum = bkg_mean * apertures_area
phot_table['bkg_sum'] = bkg_sum
phot_table['aperture_sum_bkgsub'] = phot_table['aperture_sum'] - bkg_sum
phot_table.sort('aperture_sum', reverse=True)
print(phot_table)

# visualize and plot locations of detected sources
norm = simple_norm(data, 'sqrt', percent=99.9)

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 6), dpi=150)

# plot the raw data
ax1.imshow(data, cmap='Greys', origin='lower', norm=norm, interpolation='nearest')
ax1.set_title('Raw Image (Stretched)')

# plot the zoomed-in data with apertures
ax2.imshow(data, cmap='Greys', origin='lower', norm=norm, interpolation='nearest')
ap_patches = apertures.plot(color='white', lw=2, label='Analysis aperture', ax=ax2)
ann_patches = annulus_apertures.plot(color='red', lw=2, label='Background annulus', ax=ax2)
handles = (ap_patches[0], ann_patches[0])
ax2.legend(loc='center', bbox_to_anchor=(0.5, 0.1), facecolor='#458989', labelcolor='white', handles=handles, prop={'weight': 'bold'})
ax2.set_title('Magnified Image of Polaris with Apertures')
plt.subplots_adjust(wspace=0.04)
plt.show()