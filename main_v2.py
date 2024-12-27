import matplotlib.pyplot as plt 
import numpy as np
import os
from astropy.stats import sigma_clipped_stats
from astropy.io import fits
from photutils.detection import DAOStarFinder
from astropy.stats import sigma_clipped_stats

base_directory = 'images/night2_alt_separated'
flux_values_avg = []
mag_values_avg = []
flux_values_all = []
mag_values_all = []

# sort directory in alphabetical order
folders = sorted(os.listdir(base_directory))

for folder in folders:
    # open file
    folder_path = os.path.join(base_directory, folder)
    file_list = sorted(os.listdir(folder_path))

    flux_values = []
    mag_values = []
    for filename in file_list: 
        filepath = os.path.join(folder_path, filename)
        try:
            with fits.open(filepath) as hdul:
                data = hdul[0].data
        except OSError: # occurs when system files such as .DS_Store are encountered
            print(f"ERROR: {filename}")
            continue

        # compute stats
        mean, median, std = sigma_clipped_stats(data, sigma=3.0)
        #print(np.array((mean, median, std)))

        # find stars
        daofind = DAOStarFinder(fwhm=12.0, threshold=5.*std, brightest=1) # keeps only the brightest star
        sources = daofind(data - median)
        #sources.sort('flux', reverse=True)
        print(f'filename {filename} with flux {sources['flux'][0]}')
        flux_values.append(sources['flux'][0])
        mag_values.append(sources['mag'][0]) # calculated as -2.5 * log10(flux)
    
    # calculate average flux and magnitude values for each folder
    flux_values_avg.append(np.mean(flux_values))
    print(f'AVERAGE FLUX VALUE FOR FOLDER {folder}: {np.mean(flux_values)}')
    mag_values_avg.append(np.mean(mag_values))
    print(f'AVERAGE MAGNITUDE FOR FOLDER {folder}: {np.mean(mag_values)}')

    flux_values_all.extend(flux_values)
    mag_values_all.extend(mag_values)

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

# Plotting
fig, axs = plt.subplots(2, 2, figsize=(12, 10))

# Plot flux values
axs[0, 0].plot(flux_values_all, 'bo-', label='Flux Values')
z = np.polyfit(range(len(flux_values_all)), flux_values_all, 1)
p = np.poly1d(z)
axs[0, 0].plot(range(len(flux_values_all)), p(range(len(flux_values_all))), "r--", label='Trendline')
axs[0, 0].set_title('Flux Values')
axs[0, 0].set_xlabel('Index')
axs[0, 0].set_ylabel('Flux')
axs[0, 0].legend()

# Plot magnitude values
axs[0, 1].plot(mag_values_all, 'ro-', label='Magnitude Values')
z = np.polyfit(range(len(mag_values_all)), mag_values_all, 1)
p = np.poly1d(z)
axs[0, 1].plot(range(len(mag_values_all)), p(range(len(mag_values_all))), "b--", label='Trendline')
axs[0, 1].set_title('Magnitude Values')
axs[0, 1].set_xlabel('Index')
axs[0, 1].set_ylabel('Magnitude')
axs[0, 1].legend()

# Plot average flux values
axs[1, 0].plot(flux_values_avg, 'go-', label='Average Flux Values')
z = np.polyfit(range(len(flux_values_avg)), flux_values_avg, 1)
p = np.poly1d(z)
axs[1, 0].plot(range(len(flux_values_avg)), p(range(len(flux_values_avg))), "r--", label='Trendline')
axs[1, 0].set_title('Average Flux Values')
axs[1, 0].set_xlabel('Folder Index')
axs[1, 0].set_ylabel('Average Flux')
axs[1, 0].legend()

# Plot average magnitude values
axs[1, 1].plot(mag_values_avg, 'mo-', label='Average Magnitude Values')
z = np.polyfit(range(len(mag_values_avg)), mag_values_avg, 1)
p = np.poly1d(z)
axs[1, 1].plot(range(len(mag_values_avg)), p(range(len(mag_values_avg))), "b--", label='Trendline')
axs[1, 1].set_title('Average Magnitude Values')
axs[1, 1].set_xlabel('Folder Index')
axs[1, 1].set_ylabel('Average Magnitude')
axs[1, 1].legend()

plt.tight_layout()
plt.show()