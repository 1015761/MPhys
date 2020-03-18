from utils import *
from numpy import zeros
import matplotlib as mpl 
mpl.use('Agg')

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from operator import itemgetter
import ast
from math import sqrt
from wotan import flatten


tic_list, radii, list_of_tics, labels, n_trans= get_tics()



def mom_dump_dist(t0, quality, alltimes):
	'''
	Returns time to nearest momentum dump.

	Takes in:
		t0: float - marked time of transit
		quality: array - bitwise flags for systematics
		alltimes: array - time series

	Returns:
		dt: float - time to nearest momentum dump
	'''

	dump_times=[]
	mom_dump = np.bitwise_and(quality, 2**5)
	for i, qual in enumerate(mom_dump):
		if qual>1:
			dump_times.append(alltimes[i])
	
	t0_array = t0*np.ones(len(dump_times))
	dt = t0_array - dump_times
	i_min = min(enumerate(abs(dt)), key=itemgetter(1))[0] 

	return dt[i_min]

	
def safe_mode_dist(lc, t0, times, SM_dur = 10):
	'''
	Finds time to the nearest downlink safe mode

	Takes in:
	 	lc: array - light curve
	 	t0: float - marked time
	 	times: array - time series
	 	SM_dur: int - number of consecutive NaNs to define a safe mode (	   		default = 10)
	Returns:
		dt: float - time to nearest safe mode

	'''

	idx = np.isnan(np.array(lc))

	SM_dur_t = times[SM_dur + 1]  -times[0] 


	nan_times = times[~idx]

	diff_times = [nan_times[i] - nan_times[i-1] for i in range(1, len(nan_times))]	

	
	edge_times = []
	for i, t in enumerate(diff_times):
		if t >= SM_dur_t:
			#Gap of at least SM_dur nans

			#If diff_times[i] == 15, then nan_times[i] and nan_times[i+1] are 15 apart - so both are 'edge  times'
	
			edge_times.append(nan_times[i])
			edge_times.append(nan_times[i+1])


	
	t_array = t0*np.ones(len(edge_times))
	dt = min(abs(t_array - edge_times))
	#dt = 1


	 #This only gives distance to nans, not direction
	return dt


def background(bkg, t):

	'''
	Finds background jitter.

	Takes in:
		bkg: array - background flux
		t: float - marked time
	Returns:
		MSE: float - Metric to quantify background jitter (See section 4.4)
	
	'''


	#Mask out transit - find residuals - MSE of residuals

	t = np.array(t)


	mask = (abs(t)<0.25) & (~np.isnan(bkg))

	

	times_mask = t[mask]

	bkg_mask = bkg[mask]


	if len(times_mask) ==0:
		return np.nan


	deg = 3

	poly_coeffs = np.polyfit(times_mask, bkg_mask, deg)

	trend = np.polyval(poly_coeffs, t)




	jitter = np.array(bkg)/trend


	jitter_in = jitter[mask]

	jitter_out = jitter[~mask]

	not_nan = (~np.isnan(bkg))


	res2_in = (jitter_in-1)**2
	res2_out = (jitter_out -1)**2
	
	MSE = np.nanmean(res2_in)/np.nanmean(res2_out)
	
	return MSE





def SNR(alltimes, lc, dip, x0):
	
	'''
	Calculates SNR for event

	Takes in:
		alltimes: array - Time series
		lc: array - light curve
		dip: float - depth of transit
		x0: tuple - best fit parameters
	Returns:
		SNR: float - Signal to noise ratio

	'''

	(tmid, dur) = x0[:2]

	lc = np.array(lc)

	if np.all(np.isnan(lc)):
		return np.nan

	#To find SNR - take detrended curve to find overall stddev - wotan preserves SNR (I think ...)

	lc_flat = flatten(alltimes, lc, window_length = 0.5, method = 'biweight', return_trend = False)

	masked_lc = lc_flat
	masked_lc[int(tmid - dur*0.5):int(tmid+dur*0.5)]=np.nan

	sigma_oot = np.nanstd(masked_lc)




	return dip/sigma_oot

	#Strictly sqrt(n_transit)*dip/sigma_oot but only consider one transit


def centr_jitter(x_centr, y_centr, alltimes, t0):

	'''
	Finds the centroid jitter

	Takes in:
		x_centr, y_centr: array - Pixel location of flux weighted centroid
		alltimes: array - Time series
		t0: float - marked time

	Returns:
		sigma_r : float - Jitter metric (See section 4.6) 

	'''

	idx = abs(t) < 1.5


	x_trans = np.array(x_centr)[idx]
	y_trans = np.array(y_centr)[idx]
	alltimes = alltimes[idx]
	


	mask = (abs(t)>0.25)

	x_mask = x_trans[mask]
	y_mask = y_trans[mask]
	t_mask = alltimes[mask]

	
	x_coeffs = np.polyfit(t_mask, x_mask, 3)
	y_coeffs = np.polyfit(t_mask, y_mask, 3)

	x_trend = np.polyval(x_coeffs, alltimes)
	y_trend = np.polyval(y_coeffs, alltimes)

	x_detrend = x_trans/x_trend
	y_detrend = y_trans/y_trend

	sigma_x = np.nanstd(x_trans[~mask])/np.nanstd(x_mask)
	sigma_y = np.nanstd(y_trans[~mask])/np.nanstd(y_mask)
	sigma_r = sqrt(sigma_x**2 + sigma_y**2)

	return sigma_r 















	