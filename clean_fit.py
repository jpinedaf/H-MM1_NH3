import numpy as np
from config import \
    file_NH3_thick, file_NH3_thin, \
    file_NH3_Tex, file_NH3_eTex, file_NH3_Tkin, file_NH3_eTkin, \
    file_NH3_Ncol, file_NH3_eNcol, file_NH3_Vlsr, file_NH3_eVlsr, \
    file_NH3_dv,  file_NH3_edv

from astropy.io import fits

#
# NH3 combination
#
fit_NH3_thick, hd_NH3 = fits.getdata(file_NH3_thick,
                                     header=True)
fit_NH3_thin = fits.getdata(file_NH3_thin)
# prepare FITS header for 2D files
hd_NH3_2d = hd_NH3.copy()
key_list = ['NAXIS3', 'CRPIX3', 'CDELT3', 'CUNIT3', 'CTYPE3', 'CRVAL3',
            'OBSGEO-Z', 'OBSGEO-Y', 'OBSGEO-X', 'SPECSYS']
for key_i in key_list:
    hd_NH3_2d.remove(key_i)
hd_NH3_2d['NAXIS'] = 2
hd_NH3_2d['WCSAXES'] = 2
# Optically thick fit good when Tk +- <1K otherwise fix Tk
mask = (fit_NH3_thick[6, :, :] < 1.0) & \
       (fit_NH3_thick[6, :, :] != 0.0)
fit_NH3 = mask * fit_NH3_thick + (1 - mask) * fit_NH3_thin
#
NH3_Tkin = fit_NH3[0, :, :]
NH3_eTkin = fit_NH3[6, :, :]
NH3_Tex = fit_NH3[1, :, :]
NH3_eTex = fit_NH3[7, :, :]
NH3_Ncol = fit_NH3[2, :, :]
NH3_eNcol = fit_NH3[8, :, :]
NH3_dv = fit_NH3[3, :, :]
NH3_edv = fit_NH3[9, :, :]
NH3_Vlsr = fit_NH3[4, :, :]
NH3_eVlsr = fit_NH3[10, :, :]
# Vlsr +- 0.05 km/s
bad = (NH3_eVlsr > 0.02) | (NH3_eVlsr == 0.0)
NH3_Vlsr[bad] = np.nan
NH3_eVlsr[bad] = np.nan
# dv +- 0.05 km/s
bad = (NH3_edv > 0.015) | np.isnan(NH3_Vlsr)
NH3_dv[bad] = np.nan
NH3_edv[bad] = np.nan
# Tkin +- 1 K
bad = (NH3_eTkin > 1.0) | (NH3_eTkin == 0.0) | np.isnan(NH3_dv)
NH3_Tkin[bad] = np.nan
NH3_eTkin[bad] = np.nan
# Tex +- 1 K
bad = (NH3_eTex > 1.0) | np.isnan(NH3_eTkin)
NH3_Tex[bad] = np.nan
NH3_eTex[bad] = np.nan

NH3_Ncol[bad] = np.nan
NH3_eNcol[bad] = np.nan
#
hd_NH3_2d['BUNIT'] = 'K'
fits.writeto(file_NH3_Tex, NH3_Tex, hd_NH3_2d, overwrite=True)
fits.writeto(file_NH3_eTex, NH3_eTex, hd_NH3_2d, overwrite=True)
hd_NH3_2d['BUNIT'] = 'K'
fits.writeto(file_NH3_Tkin, NH3_Tkin, hd_NH3_2d, overwrite=True)
fits.writeto(file_NH3_eTkin, NH3_eTkin, hd_NH3_2d, overwrite=True)
hd_NH3_2d['BUNIT'] = 'log(cm-2)'
fits.writeto(file_NH3_Ncol, NH3_Ncol, hd_NH3_2d, overwrite=True)
fits.writeto(file_NH3_eNcol, NH3_eNcol, hd_NH3_2d, overwrite=True)
hd_NH3_2d['BUNIT'] = 'km/s'
fits.writeto(file_NH3_Vlsr, NH3_Vlsr, hd_NH3_2d, overwrite=True)
fits.writeto(file_NH3_eVlsr, NH3_eVlsr, hd_NH3_2d, overwrite=True)
hd_NH3_2d['BUNIT'] = 'km/s'
fits.writeto(file_NH3_dv, NH3_dv, hd_NH3_2d, overwrite=True)
fits.writeto(file_NH3_edv, NH3_edv, hd_NH3_2d, overwrite=True)
