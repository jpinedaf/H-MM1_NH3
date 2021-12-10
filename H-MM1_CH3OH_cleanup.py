from astropy.io import fits
import numpy as np

data, hd = fits.getdata('Oph-H-MM1_ACA+12M_CH3OH_0_natural_mscale_v2_TdV.fits', header=True)
key_list = ['NAXIS3', 'NAXIS4', 'PC03_01', 'PC04_01', 
            'PC03_02', 'PC04_02', 
            'PC01_03', 'PC02_03', 'PC03_03', 'PC04_03',
            'PC01_04', 'PC02_04', 'PC03_04', 'PC04_04',
            'VELREF', 'CTYPE3', 'CRVAL3', 'CDELT3', 'CRPIX3', 
            'CUNIT3', 'CTYPE4', 'CRVAL4', 'CDELT4', 'CRPIX4', 'CUNIT4']

for key_i in key_list:
    hd.remove(key_i)

hd['NAXIS'] = 2

fits.writeto('H-MM1_CH3OH_TdV.fits', np.squeeze(data), hd, overwrite=True)
