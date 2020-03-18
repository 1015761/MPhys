import os
import ast

if os.name == 'nt':

	import math
	import numpy as np
	from numpy import ones, arange
	import astropy.io.fits as pf
	from astropy import units as u 
	from astropy.coordinates import SkyCoord
	import matplotlib.pyplot as plt 
	import pandas as pd
	import csv
	from statistics import mean, median
	from random import randrange, random
	import scipy as sp
	from time import time, sleep
	import lmfit as lm
	#from progress.bar import ChargingBar



	#On Personal windows device
	dir ="C:/Users/blake"


elif os.name == 'posix':


	import math
	import numpy as np
	from numpy import ones, arange
	from numpy import nanmean as mean
	import astropy.io.fits as pf
	#import matplotlib.pyplot as plt 
	#import eleanor
	import pandas as pd
	import csv
	from random import randrange, random
	import scipy as sp
	from time import time
	import os
	import fnmatch
	import re
	import lmfit as lm

	dir ="/mnt/users/kepler/kepler2/TESS/planethunters"

bin_fac_default = 5


def unpack_params(params):
	params_dict = params.valuesdict()

	x0 = [params_dict['tmid']]
	x0.append(params_dict['dur'])
	x0.append(params_dict['v'])
	x0.append(params_dict['h'])
	x0.append(params_dict['p3'])
	x0.append(params_dict['p2'])
	x0.append(params_dict['p1'])
	x0.append(params_dict['p0'])
	return x0




def gen_model(x0, t):
    [tmid, dur, v, h] = x0[:4]

    d2 = dur/2.
    ingr = v*d2
    full = d2 - ingr

    tcen = np.array(t) - tmid
    in_full = abs(tcen) <full
    in_gr = (~in_full) & (abs(tcen)<d2)

    transit = np.ones_like(tcen)
    transit[in_full] -= h
    transit[in_gr] -= h*(d2 - abs(tcen[in_gr]))/ingr

    trend = np.polyval(x0[4:], tcen) + 1

    return transit * trend


def write_params(x0, poly = True):

	#Get params in lmfit.Parameters() instance

	if len(x0) <8:
		x0 = (*x0[:4], 0., 0., 0., 0.)

	params = lm.Parameters()
	params.add('tmid', value = float(x0[0]), min = -1.5, max = 1.5, vary = True)
	params.add('dur', value = float(x0[1]), min= 0., vary = True)
	params.add('v', value = float(x0[2]),min = 0., max = 1., vary = True)
	params.add('h', value = float(x0[3]),min = 0.,  max = 1., vary = True)
	params.add('p3', value = float(x0[4]), vary = poly )
	params.add('p2', value = float(x0[5]), vary = poly )
	params.add('p1', value = float(x0[6]), vary = poly )
	params.add('p0', value = float(x0[7]), vary = poly )

	return params


def resids(params, times, data, eps = None):
    #unpack parameters
        
    x0 = unpack_params(params)


    times = np.array(times)    
    model = gen_model(x0, times)
    if eps is None:
        return data - model
    return (data-model)/eps



def fit(x0, transit, transit_times, transit_errs, poly = True):

    '''
    Takes in:

    x0: list of parameters

    transit: light curve, centered on marked time from PHT

    transit_times: time relative to marked time (t=0 is centre)

    transit_errs: as above for errors


    Returns:

    Success: bool to see if fitting succeeds

    fit_values: Best fit parameters
    fit_stderrs: Best fit uncertainties (excluding correlations)

    '''

    params = write_params(x0, poly)
    soln = lm.minimize(resids, params, nan_policy = 'omit', method = 'leastsq', args = (transit_times, transit, transit_errs))
    #Unpack parameters to list
    #fit_names = []
    fit_values = []
    fit_stderrs = []
    #fit_correls = []
    for name in soln.params:
        par = soln.params[name]
        #fit_names.append(par.name)
        fit_values.append(par.value)
        fit_stderrs.append(par.stderr)
        #fit_correls.append(par.correl)

    #return soln.success, fit_names, fit_values, fit_stderrs , fit_correls
    return soln





def bin_index(index, bin_fac = bin_fac_default):

	#Takes in index of point in unbinned data and returns index of
	# corresponding point in binned data

	index_bin = math.ceil(index/bin_fac)
	return index_bin


def remove_nans(lc, err ,nan_value = 1.):
	
	#Replaces all nan values with a predetermined value
	
	df = pd.DataFrame(lc)
	df_error = pd.DataFrame(err)

	df.fillna(nan_value, inplace = True)
	df_error.fillna(np.inf, inplace = True)
	
	y=df[0].values.tolist()
	err = df_error[0].values.tolist()
	return y, err

def read_file(tic, sec, indir = dir, plot = False, index = False, bin = True):

	'''
	Reads a locally downloaded FITS file for a given TIC in a given sector.

	Takes in:
		tic: Int - TIC ID of star
		sec: Int - Observational sector
		indir: Str - location of files
		plot: Bool - whether to plot light curve from fits file (default False)
		bin: Bool - whether to bin datapoints to 10 min cadences (default True)

	Returns:
		alltimes: Array -array of times
		lc: Array - Time series flux (systematics corrected)
		err: Array - Error on flux measurements
		qual: Array - Bitwise flags for known systematics
		bkg: Array - Background flux
		x1, y1: Array - Centroid position in pixels
	'''

	if os.name == 'nt':

		fits_dir = "{}/{}".format(indir,"MPhys_Data/FITS files")
		sector =  "s{}".format(str(10000 + sec)[1:])
		filename = "{}_{}.fits".format(sector, tic)

		lchdu = pf.open("{}/{}".format(fits_dir, filename))




		lcdata = lchdu[1].data


		sapf = lcdata["SAP_FLUX"]
		err = lcdata["SAP_FLUX_ERR"]
		pdc_sapf = lcdata['PDCSAP_FLUX']
		qual = lcdata["QUALITY"]
		bkg = lcdata["SAP_BKG"]
		alltimes = np.arange(len(pdc_sapf)) if index else lcdata['TIME']

		x1 = lcdata['MOM_CENTR1']
		y1 = lcdata['MOM_CENTR2']


		
		
	elif os.name == 'posix':

		if sec <=9:
			fits_dir = '/mnt/zfsusers/kepler/kepler2/TESS/Sector{}/light_curves/two_min'.format(sec)
		elif sec ==11:
			fits_dir = '/mnt/zfsusers/kepler/kepler2/TESS/planethunters/Rel11/Sector11_skip/light_curves/two_min'
		elif sec == 16:
			fits_dir = '/mnt/zfsusers/kepler/kepler2/TESS/planethunters/Rel16/Sector16/light_curves/twp_min'
		else:
			fits_dir = '/mnt/zfsusers/kepler/kepler2/TESS/planethunters/Rel{}/Sector{}/light_curves/two_min'.format(sec, sec)

		



		sector =  "s{}".format(str(10000 + sec)[1:])
		ticname = str(int(1e16) + int(tic))[1:]
		main = "{}-{}".format(sector, ticname)
		file_format = r"^tess.+{}.+$".format(main)

		filename="None"

		for file in os.listdir(fits_dir):
			if re.match(file_format, file):


				filename = file


				lchdu = pf.open("{}/{}".format(fits_dir, filename))




				lcdata = lchdu[1].data


				sapf = lcdata["SAP_FLUX"]
				err = lcdata["PDCSAP_FLUX_ERR"]
				pdc_sapf = lcdata['PDCSAP_FLUX']
				qual = lcdata["QUALITY"]
				bkg = lcdata["SAP_BKG"]
				alltimes = np.arange(len(pdc_sapf)) if index else lcdata['TIME']
				x1 = lcdata['MOM_CENTR1']
				y1 = lcdata['MOM_CENTR2']
				break
		else:
			message = "File not found: {} {}  --- Obtained from FFI".format(tic, sec)

			raise ValueError(message)


		

	def binning(bin_fac=bin_fac_default):
		#Bins all quantities by a predetermined factor
		N = len(sapf)
		lc_bin = []
		t_bin = []
		err_bin=[]
		pdc_lc_bin = []
		bkg_bin = []
		x1_bin = []
		y1_bin = []

		for i in range(N//bin_fac):
			start = i*bin_fac

			end = min([N, (i+1)*bin_fac - 1 ])

			if np.all(np.isnan(pdc_sapf[start:end])):
				pdc_lc_bin.append(np.nan)
			else:
				pdc_lc_bin.append(np.nanmean(pdc_sapf[start:end]))

			t_bin.append(np.nanmean(alltimes[start:end]))

			if np.all(np.isnan(err[start:end])):
				err_bin.append(np.nan)
			else:
				err_bin.append(np.nanmean(err[start:end]))

			if np.all(np.isnan(bkg[start:end])):
				bkg_bin.append(np.nan)
			else:
				bkg_bin.append(np.nanmean(bkg[start:end]))
			x1_bin.append(np.nanmean(x1[start:end]))
			y1_bin.append(np.nanmean(y1[start:end]))




		return  pdc_lc_bin, t_bin, err_bin, bkg_bin, x1_bin, y1_bin



	if bin:
		lc, alltimes, err, bkg, x1, y1 = binning()

	else:
		lc = pdc_sapf

	scale = np.nanmedian(sapf)
	scale_pdc = np.nanmedian(pdc_sapf)

	lc_norm=[float(i)/scale for i in lc]
	err_norm = [float(i)/scale for i in err]



	return np.array(alltimes), np.array(lc_norm), np.array(err_norm), np.array(qual), np.array(bkg), np.array(x1), np.array(y1)




def time_to_index(time, start_time =0.008333,  cadence = 2):
	#Turns a BJD time into index in the data
	interval = cadence / 1440
	index=[]
	if not isinstance(time, list):
		time = [time]

	for t in time:


		index.append(int((t-start_time)/interval))


	return index

def get_tics(indir = dir):
	'''
	Reads list of TICs that have been passed through PHT

	Takes in:
		indir: Str - location of list of TICs
	Returns:
		times: Dict - dictionary of each marked transit for each TIC in each 	sector - reference as times[tic][sec] = [t1, t2, ...]
		star_radii: Dict - Stellar radius (in solar radii) for each TIC
		TIC_list: List - list of TICs to be analysed
		labels: Dict - List of classifications given to each light curve
		n_trans_dict: Dict - number of transits marked in each light curve
	'''

	if os.name == 'nt':
		dirname = "{}/{}".format(indir, "MPhys")
	elif os.name == 'posix':
		dirname = "/mnt/zfsusers/blakeland/Documents/MPhys"

	filename = "{}/TICs.csv".format(dirname)


	file = pd.read_csv(filename)

	TICs = list(file['TIC'])
	TIC_list = list(set(TICs)) #Contains each TIC no more than once
	sectors = list(file['Sector'])
	times_str = list(file['Times'])
	radius = list(file['Radius'])
	label = list(file['Label'])
	n_trans = list(file['N_transits'])


	assert len(times_str) == len(TICs)

	N_lc = len(TICs)

	times_list=[]
	for i in times_str:
		try:
			items = ast.literal_eval(i)
			if items != ['']:
				try:
					this_list = [float(j) for j in items]
				except:
					print("Failed on {}: type : str".format(items))
					break
			else:
				this_list = []
			times_list.append(this_list)
		except:
			if i != ['']:
				items = i[1:-1].split(' ') 
				this_list = [float(time) for time in items]
			else:
				this_list = []


	star_radii= {TICs[i]: radius[i] if radius[i] != -99 else 1 for i in range(N_lc)}
	n_trans_dict = {TICs[i]: {sectors[i]: n_trans[i]} for i in range(N_lc)}

	output = {}
	labels = {}
	n_trans_dict = {}

	for i in range(N_lc):
		tic = TICs[i]
		sec = sectors[i]
		keep = label[i]
		times = times_list[i]

		if tic in output:
			output[tic][sec] = times
		else:
			output[tic] = {sec: times}
			
		if tic in labels:
			labels[tic][sec] = keep
		else:
			labels[tic] = {sec: keep}

		if tic in n_trans_dict:
			n_trans_dict[tic][sec] = n_trans[i]
		else:
			n_trans_dict[tic] = {sec: n_trans[i]}


	return output, star_radii, TIC_list, labels, n_trans_dict
