from  astropy.utils.data import download_file
import os

# url with data cubes stored in Dataverse
url_list = ['https://dataverse.harvard.edu/api/access/datafile/6092603',
			'https://dataverse.harvard.edu/api/access/datafile/6092604']
file_list = ['data/H-MM1_NH3_11.fits', 'data/H-MM1_NH3_22.fits']

# download and rename files
for url_i, file_out_i in zip(url_list, file_list):
	file_i = download_file(url_i)
	os.rename(file_i, file_out_i)
