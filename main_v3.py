import matplotlib.pyplot as plt 
import numpy as np
import os
from astropy.stats import sigma_clipped_stats
from astropy.io import fits
from photutils.detection import DAOStarFinder
from astropy.stats import sigma_clipped_stats
from photutils.aperture import CircularAperture, aperture_photometry, CircularAnnulus, ApertureStats

base_directory = 'images/night2_alt_separated'
flux_values_avg = []
flux_values_all = []
altitudes = []

# sort directory in numerical order, high to low
folders = sorted(os.listdir(base_directory), key=lambda x: int(''.join(filter(str.isdigit, x))), reverse=True)

for folder in folders:
    # open file
    folder_path = os.path.join(base_directory, folder)
    file_list = sorted(os.listdir(folder_path))

    flux_values = []
    for filename in file_list: 
        filepath = os.path.join(folder_path, filename)
        try:
            with fits.open(filepath) as hdul:
                data = hdul[0].data
        except OSError: # occurs when system files (such as .DS_Store) are encountered
            print(f"ERROR: {filename}")
            continue
        
        # subtract background (median)
        mean, median, std = sigma_clipped_stats(data, sigma=3.0)

        # find stars
        daofind = DAOStarFinder(fwhm=13, threshold=2500.*std)
        sources = daofind(data - median)

        # create apertures and perform photometry
        if sources is not None:
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
        else:
            print(f"No sources found in {filename}")
    
    # calculate average flux and magnitude values for each folder
    flux_values_avg.append(np.mean(flux_values))
    print(f'AVERAGE FLUX VALUE FOR FOLDER {folder}: {np.mean(flux_values)}')
    flux_values_all.extend(flux_values)

    # get altitude from folder name
    altitudes.append(int(folder))


# calculate magnitude from flux
def calculate_magnitude(flux_values):
    return -2.5 * np.log10(flux_values)

magnitude_all = calculate_magnitude(np.array(flux_values_all))
magnitude_avg = calculate_magnitude(np.array(flux_values_avg))

# calculate change in magnitude since no calibration point exists
delta_magnitude_all = magnitude_all - magnitude_all[0]
delta_magnitude_avg = magnitude_avg - magnitude_avg[0]

# visualize the flux and magnitude values for each file and folder
plt.figure(figsize=(15, 10))

# plot for flux_values_all
plt.subplot(2, 2, 1)
plt.plot(flux_values_all, label='Flux Values All')
plt.xlabel('File Index')
plt.ylabel('Flux Value')
plt.title('Flux Values All')
plt.legend()

# plot for flux_values_avg with altitudes as x-axis
plt.subplot(2, 2, 2)
plt.plot(altitudes, flux_values_avg, label='Flux Values Avg', color='orange')
plt.xlabel('Altitude')
plt.ylabel('Flux Value')
plt.title('Flux Values Avg by Altitude')
plt.legend()

# plot for delta_magnitude_all
plt.subplot(2, 2, 3)
plt.plot(delta_magnitude_all, label='Change in Magnitude All')
plt.xlabel('File Index')
plt.ylabel('Change in Magnitude')
plt.title('Change in Magnitude All')
plt.legend()

# plot for delta_magnitude_avg with altitudes as x-axis
plt.subplot(2, 2, 4)
plt.plot(altitudes, delta_magnitude_avg, label='Change in Magnitude Avg', color='orange')
plt.xlabel('Altitude')
plt.ylabel('Change in Magnitude')
plt.title('Change in Magnitude Avg by Altitude')
plt.legend()

plt.tight_layout()
plt.show()

# UNCOMMENT TO VISUALIZE AND DEBUG INDIVIDAUL FILES
# filepath = os.path.join('images/night2/Light_Polaris_10.0s_IRCUT_20241127-202337.fit')
# hdul = fits.open(filepath)
# data = hdul[0].data

# # subtract background (median)
# mean, median, std = sigma_clipped_stats(data, sigma=3.0)

# # find stars
# daofind = DAOStarFinder(fwhm=13, threshold=2500.*std)
# sources = daofind(data - median)

# # create apertures and perform photometry
# positions = np.transpose((sources['xcentroid'], sources['ycentroid']))
# apertures = CircularAperture(positions, r=19.)
# phot_table = aperture_photometry(data, apertures)
# print(phot_table)

# # subtract local background
# annulus_apertures = CircularAnnulus(positions, r_in=35, r_out=45)
# background_stats = ApertureStats(data, annulus_apertures)
# bkg_mean = background_stats.mean
# print(bkg_mean)
# apertures_area = apertures.area_overlap(data)
# bkg_sum = bkg_mean * apertures_area
# phot_table['bkg_sum'] = bkg_sum
# phot_table['aperture_sum_bkgsub'] = phot_table['aperture_sum'] - bkg_sum
# phot_table.sort('aperture_sum', reverse=True)
# print(phot_table)

# # visualize and plot locations of detected sources
# norm = simple_norm(data, 'sqrt', percent=99.9)
# plt.imshow(data, cmap='Greys', origin='lower', norm=norm, interpolation='nearest')
# ap_patches = apertures.plot(color='white', lw=2,
#                            label='Photometry aperture')
# ann_patches = annulus_apertures.plot(color='red', lw=2,
#                                     label='Background annulus')
# handles = (ap_patches[0], ann_patches[0])
# plt.legend(loc=(0.17, 0.05), facecolor='#458989', labelcolor='white',
#            handles=handles, prop={'weight': 'bold', 'size': 11})
# plt.show()