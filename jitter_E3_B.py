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
dir = '/Users/mikako/work/teaching/2025_ObsTech/LCO/data/E3/B/'

#file_names = glob.glob(dir+'lsc*.fits.fz')
file_names = glob.glob(dir+'lsc1m004-fl04-20141120-010?-e90.fits')
# satellite running across - not good
#file_names = glob.glob(dir+'lsc1m004-fl04-20141120-011?-e90.fits')
#file_names = glob.glob(dir+'lsc1m004-fl04-20141121-011?-e90.fits')
#not good image
#file_names = glob.glob(dir+'lsc1m004-fl04-20141124-0051-e90.fits')

header_layer = 0

hdu_list = {}
for i in range(len(file_names)):
    hdu_list[i] = fits.open(file_names[i])[header_layer]
    print(i, hdu_list[i].header['DATE-OBS'], hdu_list[i].header['EXPTIME'], hdu_list[i].header['FILTER'])


# Build CCDData objects so wcs_project can access WCS metadata
ccd0 = CCDData(hdu_list[0].data, meta=hdu_list[0].header, unit=u.adu, wcs=WCS(hdu_list[0].header))



#reproject_list = {}
#for i in range(len(file_names)):  
#    ccd = CCDData(hdu_list[i].data, meta=hdu_list[i].header, unit=u.adu, wcs=WCS(hdu_list[i].header))
#    reprojected_image = wcs_project(ccd, WCS(hdu_list[0].header))
#    reproject_list[i] = reprojected_image

# no jitter - simply add together
#reproject_list = {}
#for i in range(len(file_names)):  
#    ccd = CCDData(hdu_list[i].data, meta=hdu_list[i].header, unit=u.adu, wcs=WCS(hdu_list[i].header))
#    reproject_list[i] = ccd


# sigma clilpping
#combiner = Combiner( list(reproject_list.values()))
#combiner = Combiner([reproject_list[1], reproject_list[2], reproject_list[3]])


##This is 150 sec exposure files - no jitter 
#selected_indices = [3, 14, 15, 20, 17, 24, 29, 31]
# There is some satellinte running - image 3
#selected_indices = [3, 15, 29, 31]
# tracking problem 
#selected_indices = [14, 17, 20, 24]

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


#combiner = Combiner([reproject_list[i] for i in selected_indices])
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
ax1.imshow(hdu_list[0].data, origin='lower', vmin=-1, vmax=900)
ax1.coords['ra'].set_axislabel('Right Ascension')
ax1.coords['dec'].set_axislabel('Declination')
ax1.set_title('Original (HDU1)')

ax2 = plt.subplot(1,2,2, projection=WCS(hdu_list[0].header))
ax2.imshow(final_image.data, origin='lower', vmin=-1, vmax=900)
ax2.coords['ra'].set_axislabel('Right Ascension')
ax2.coords['dec'].set_axislabel('Declination')
ax2.set_title('Combined (sigma-clipped)')

plt.tight_layout()
plt.show()

