#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Dec 28 20:40:42 2025

@author: mikako
"""
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Dec 27 21:14:55 2025

@author: mikako
"""
from ccdproc import CCDData, wcs_project, Combiner
from astropy.io import fits
from astropy.wcs import WCS
from astropy import units as u
import numpy as np
import os
import matplotlib.pyplot as plt
import glob 

### PATH ####
# to call montage from Spyder, this path setting is needed
# If calling python from terminal, this does not need.
os.environ['PATH'] += ':/Users/mikako/bin/montage/bin'
dir = '/Users/mikako/work/teaching/2025_ObsTech/LCO/data/E3/V/'

#file_names = glob.glob(dir+'lsc*.fits.fz')

# 50 s exposure - too short
#file_names = glob.glob(dir+'lsc1m004-fl04-20141027-014?-e90.fits')

# 150 s exposure reasonablly good, with slight tracking problem
file_names = glob.glob(dir+'lsc1m004-fl04-20141120-012?-e90.fits')
# reasonablly good, apart from bias levels are not perfect
#file_names = glob.glob(dir+'lsc1m004-fl04-20141121-01??-e90.fits')

# not the good set
#file_names = glob.glob(dir+'lsc1m004-fl04-20141124-005?-e90.fits')

header_layer = 0

hdu_list = {}
for i in range(len(file_names)):
    hdu_list[i] = fits.open(file_names[i])[header_layer]
    print(i, hdu_list[i].header['DATE-OBS'], hdu_list[i].header['EXPTIME'], hdu_list[i].header['FILTER'])

# If decided to chose the best images, make a list of files, here.
selected_indices = range(len(hdu_list))

reproject_list = {}
for i in range(len(selected_indices)):  
    index = selected_indices[i]
    ccd = CCDData(hdu_list[index].data, meta=hdu_list[index].header, unit=u.adu, wcs=WCS(hdu_list[index].header))
    reprojected_image = wcs_project(ccd, WCS(hdu_list[selected_indices[0]].header))
    reproject_list[i] = reprojected_image


hdu_list = {}
for i in range(len(selected_indices)):
    index  = selected_indices[i]
    hdu_list[i] = fits.open(file_names[index])[header_layer]
    print(index, file_names[index], hdu_list[i].header['DATE-OBS'], hdu_list[i].header['EXPTIME'], hdu_list[i].header['FILTER'])


combiner = Combiner( list(reproject_list.values()))
                     
combiner.sigma_clipping(low_thresh=2, high_thresh=10, func=np.ma.median)
final_image = combiner.average_combine()    


# Export the final combined image to FITS with header from HDU1
# Fill masked pixels with NaN if present before writing
data_to_write = final_image.data.filled(np.nan) if np.ma.isMaskedArray(final_image.data) else final_image.data
hdr = hdu_list[0].header.copy()
hdr['HISTORY'] = 'Sigma-clipped average combination from frames 1-'+str(len(hdu_list))+'.'
# Update exposure-related keywords to reflect 4 combined frames
for key in ['EXPTIME', 'EXP_TIME', 'EXPOSURE', 'ITIME']:
        if key in hdr:
                try:
                        hdr[key] = float(hdr[key]) * len(hdu_list)
                        hdr['HISTORY'] = f'{key} scaled by ' + str(len(hdu_list))+'x for combined exposure.'
                except Exception:
                        pass
hdr['NCOMBINE'] = len(hdu_list)
out_fits = os.path.join(dir, 'stacked_image.fits')
fits.PrimaryHDU(data=data_to_write, header=hdr).writeto(out_fits, overwrite=True)


# Compare the original frame (HDU1) with the combined image
fig = plt.figure()
ax1 = plt.subplot(1,2,1, projection=WCS(hdu_list[0].header))
ax1.imshow(hdu_list[0].data, origin='lower', vmin=-1, vmax=400)
ax1.coords['ra'].set_axislabel('Right Ascension')
ax1.coords['dec'].set_axislabel('Declination')
ax1.set_title('Original (HDU1)')

ax2 = plt.subplot(1,2,2, projection=WCS(hdu_list[0].header))
ax2.imshow(final_image.data, origin='lower', vmin=-1, vmax=400)
ax2.coords['ra'].set_axislabel('Right Ascension')
ax2.coords['dec'].set_axislabel('Declination')
ax2.set_title('Combined (sigma-clipped)')

plt.tight_layout()
plt.show()

