from spectral_cube import SpectralCube
import astropy.units as u
from astropy.io import fits

file_in = 'data/Oph-H-MM1_ACA+12M_oNH2D_natural_mscale_beam4.fits'
file_out = 'data/Oph-H-MM1_oNH2D_TdV.fits'
cube = SpectralCube.read(file_in).with_spectral_unit(u.km / u.s,
                                                     velocity_convention='radio')
subcube = cube[75:90, :, :]
mom0 = subcube.moment(order=0)
mom0.write(file_out, overwrite=True)

hdul = fits.open(file_out)  # open a FITS file
hdr = hdul[0].header
del hdr['OBSGEO-X']
del hdr['OBSGEO-Y']
del hdr['OBSGEO-Z']
hdul.writeto(file_out, overwrite=True)
hdul.close()

file_in = 'data/Oph-H-MM1_ACA+12M_pNH2D_natural_mscale_beam4.fits'
file_out = 'data/Oph-H-MM1_pNH2D_TdV.fits'
cube = SpectralCube.read(file_in).with_spectral_unit(u.km / u.s,
                                                     velocity_convention='radio')
subcube = cube[75:89, :, :]
mom0 = subcube.moment(order=0)
mom0.write(file_out, overwrite=True)

hdul = fits.open(file_out)  # open a FITS file
hdr = hdul[0].header
del hdr['OBSGEO-X']
del hdr['OBSGEO-Y']
del hdr['OBSGEO-Z']
hdul.writeto(file_out, overwrite=True)
hdul.close()
