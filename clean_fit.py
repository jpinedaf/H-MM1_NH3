import numpy as np
from config import file_oNH2D_thin, file_oNH2D_thick, \
    file_oNH2D_Tex, file_oNH2D_eTex, file_oNH2D_tau, file_oNH2D_etau, \
    file_oNH2D_Vlsr, file_oNH2D_eVlsr, file_oNH2D_dv, file_oNH2D_edv, \
    file_pNH2D_thin, file_pNH2D_thick, \
    file_pNH2D_Tex, file_pNH2D_eTex, file_pNH2D_tau, file_pNH2D_etau, \
    file_pNH2D_Vlsr, file_pNH2D_eVlsr, file_pNH2D_dv, file_pNH2D_edv, \
    file_NH3_thick, file_NH3_thin, \
    file_NH3_Tex, file_NH3_eTex, file_NH3_Tkin, file_NH3_eTkin, \
    file_NH3_Ncol, file_NH3_eNcol, file_NH3_Vlsr, file_NH3_eVlsr, \
    file_NH3_dv,  file_NH3_edv

from astropy.io import fits

do_NH3 = True
do_oNH2D = True
do_pNH2D = True

if do_NH3:
    # NH3 combination
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


if do_oNH2D:
    fit_oNH2D_thick, hd_oNH2D = fits.getdata(file_oNH2D_thick,
                                             header=True)
    fit_oNH2D_thin = fits.getdata(file_oNH2D_thin)
    hd_oNH2D_2d = hd_oNH2D.copy()
    key_list = ['NAXIS3', 'CRPIX3', 'CDELT3', 'CUNIT3', 'CTYPE3', 'CRVAL3',
                'OBSGEO-Z', 'OBSGEO-Y', 'OBSGEO-X']
    for key_i in key_list:
        hd_oNH2D_2d.remove(key_i)
    hd_oNH2D_2d['NAXIS'] = 2
    hd_oNH2D_2d['WCSAXES'] = 2

    snr_min_tau = 3.0
    # Optically thick if tau > 3* etau
    # Tex != 20 K   error_tau < 1
    oNH2D_mask = (fit_oNH2D_thick[1, :, :] > snr_min_tau * fit_oNH2D_thick[5, :, :]) & \
                 (fit_oNH2D_thick[1, :, :] != 50.0) & (fit_oNH2D_thick[0, :, :] != 20.0) & \
                 (fit_oNH2D_thick[5, :, :] < 1.0)
    fit_oNH2D = fit_oNH2D_thick * oNH2D_mask + (1 - oNH2D_mask) * fit_oNH2D_thin
    oNH2D_Tex = fit_oNH2D[0, :, :]
    oNH2D_eTex = fit_oNH2D[4, :, :]
    oNH2D_tau = fit_oNH2D[1, :, :]
    oNH2D_etau = fit_oNH2D[5, :, :]
    oNH2D_Vlsr = fit_oNH2D[2, :, :]
    oNH2D_eVlsr = fit_oNH2D[6, :, :]
    oNH2D_dv = fit_oNH2D[3, :, :]
    oNH2D_edv = fit_oNH2D[7, :, :]
    #
    hd_oNH2D_2d['BUNIT'] = 'K'
    op_thin = (oNH2D_tau == 0.1) | (oNH2D_etau == 0.0) | (oNH2D_Tex == 0.0)
    oNH2D_Tex[op_thin] = np.nan
    oNH2D_eTex[op_thin] = np.nan
    fits.writeto(file_oNH2D_Tex, oNH2D_Tex, hd_oNH2D_2d, overwrite=True)
    fits.writeto(file_oNH2D_eTex, oNH2D_eTex, hd_oNH2D_2d, overwrite=True)
    hd_oNH2D_2d['BUNIT'] = ''
    oNH2D_tau[op_thin] = np.nan
    oNH2D_etau[op_thin] = np.nan
    fits.writeto(file_oNH2D_tau, oNH2D_tau, hd_oNH2D_2d, overwrite=True)
    fits.writeto(file_oNH2D_etau, oNH2D_etau, hd_oNH2D_2d, overwrite=True)
    hd_oNH2D_2d['BUNIT'] = 'km/s'
    bad_v = (oNH2D_Vlsr == 0.0) | (oNH2D_eVlsr == 0.0)
    oNH2D_Vlsr[bad_v] = np.nan
    oNH2D_eVlsr[bad_v] = np.nan
    fits.writeto(file_oNH2D_Vlsr, oNH2D_Vlsr, hd_oNH2D_2d, overwrite=True)
    fits.writeto(file_oNH2D_eVlsr, oNH2D_eVlsr, hd_oNH2D_2d, overwrite=True)
    hd_oNH2D_2d['BUNIT'] = 'km/s'
    bad_dv = (oNH2D_edv == 0.0) | (oNH2D_edv == 0.0)
    oNH2D_dv[bad_dv] = np.nan
    oNH2D_edv[bad_dv] = np.nan
    fits.writeto(file_oNH2D_dv, oNH2D_dv, hd_oNH2D_2d, overwrite=True)
    fits.writeto(file_oNH2D_edv, oNH2D_edv, hd_oNH2D_2d, overwrite=True)

if do_pNH2D:
    fit_pNH2D_thick, hd_pNH2D = fits.getdata(file_pNH2D_thick,
                                             header=True)
    fit_pNH2D_thin = fits.getdata(file_pNH2D_thin)
    hd_pNH2D_2d = hd_pNH2D.copy()
    key_list = ['NAXIS3', 'CRPIX3', 'CDELT3', 'CUNIT3', 'CTYPE3', 'CRVAL3',
                'OBSGEO-Z', 'OBSGEO-Y', 'OBSGEO-X']
    for key_i in key_list:
        hd_pNH2D_2d.remove(key_i)
    hd_pNH2D_2d['NAXIS'] = 2
    hd_pNH2D_2d['WCSAXES'] = 2

    snr_min_tau = 3.0
    # Optically thick if tau > 3* etau
    # Tex != 20 K   error_tau < 1
    pNH2D_mask = (fit_pNH2D_thick[1, :, :] > snr_min_tau * fit_pNH2D_thick[5, :, :]) & \
                 (fit_pNH2D_thick[1, :, :] != 50.0) & (fit_pNH2D_thick[0, :, :] != 20.0) & \
                 (fit_pNH2D_thick[5, :, :] < 1.0)
    fit_pNH2D = fit_pNH2D_thick * pNH2D_mask + (1 - pNH2D_mask) * fit_pNH2D_thin

    pNH2D_Tex = fit_pNH2D[0, :, :]
    pNH2D_eTex = fit_pNH2D[4, :, :]
    pNH2D_tau = fit_pNH2D[1, :, :]
    pNH2D_etau = fit_pNH2D[5, :, :]
    pNH2D_Vlsr = fit_pNH2D[2, :, :]
    pNH2D_eVlsr = fit_pNH2D[6, :, :]
    pNH2D_dv = fit_pNH2D[3, :, :]
    pNH2D_edv = fit_pNH2D[7, :, :]
    #
    op_thin = (pNH2D_tau == 0.1) | (pNH2D_etau == 0.0) | (pNH2D_Tex == 0.0)
    pNH2D_Tex[op_thin] = np.nan
    pNH2D_eTex[op_thin] = np.nan
    hd_pNH2D_2d['BUNIT'] = 'K'
    fits.writeto(file_pNH2D_Tex, pNH2D_Tex, hd_pNH2D_2d, overwrite=True)
    fits.writeto(file_pNH2D_eTex, pNH2D_eTex, hd_pNH2D_2d, overwrite=True)
    hd_pNH2D_2d['BUNIT'] = ''
    pNH2D_tau[op_thin] = np.nan
    pNH2D_etau[op_thin] = np.nan
    fits.writeto(file_pNH2D_tau, pNH2D_tau, hd_pNH2D_2d, overwrite=True)
    fits.writeto(file_pNH2D_etau, pNH2D_etau, hd_pNH2D_2d, overwrite=True)
    hd_pNH2D_2d['BUNIT'] = 'km/s'
    bad_v = (pNH2D_Vlsr == 0.0) | (pNH2D_eVlsr == 0.0)
    pNH2D_Vlsr[bad_v] = np.nan
    pNH2D_eVlsr[bad_v] = np.nan
    fits.writeto(file_pNH2D_Vlsr, pNH2D_Vlsr, hd_pNH2D_2d, overwrite=True)
    fits.writeto(file_pNH2D_eVlsr, pNH2D_eVlsr, hd_pNH2D_2d, overwrite=True)
    hd_pNH2D_2d['BUNIT'] = 'km/s'
    bad_dv = (pNH2D_edv == 0.0) | (pNH2D_edv == 0.0)
    pNH2D_dv[bad_dv] = np.nan
    pNH2D_edv[bad_dv] = np.nan
    fits.writeto(file_pNH2D_dv, pNH2D_dv, hd_pNH2D_2d, overwrite=True)
    fits.writeto(file_pNH2D_edv, pNH2D_edv, hd_pNH2D_2d, overwrite=True)
