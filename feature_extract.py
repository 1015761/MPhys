import matplotlib as mpl 
mpl.use('Agg')

from utils import *
import systematics as syst

import matplotlib.pyplot as plt
import numpy as np 
import os
import scipy.optimize as opt 
from wotan import flatten


def get_params(tic, sec, keep):

	"""
	For each light curve, extract features of interest.

	See section 4 of report.
	"""


	times, radii, list_of_tics, labels, n_trans = get_tics()
	#First check if the tic is a TOI

	is_TOI = False
	is_Rand = False

	if os.name == 'nt':
		file_dir = 'C://Users/blake/MPhys'
	elif os.name == 'posix':
		file_dir = '/mnt/zfsusers/blakeland/Documents/MPhys'
	file = '{}/CTC_output_{}.csv'.format(file_dir, sec)
	file_TOI = '{}/TOI_list.csv'.format(file_dir)

	if not os.path.exists(file):
		is_TOI = True
	else:
		CTC_file = pd.read_csv(file)
		TOI_file = pd.read_csv(file_TOI)
		CTC_list = list(CTC_file['TIC_ID'])
		TOI_list = list(TOI_file['TIC ID'])
		if tic in TOI_list:
			is_TOI = True
		elif not tic in CTC_list:
			#Not taken from CTC or TOI list
			is_Rand = True


	try:
		alltimes_bin, lc_bin, err_bin, qual, bkg_bin, x_bin, y_bin = read_file(tic, sec)
		alltimes, lc, err, qual, bkg, x, y = read_file(tic, sec, bin = False)
	except ValueError:
		#utils.get_tics() throws ValueError if the fits file exist
		return {'TIC': tic, 'Sector': sec, 'Times': times[tic][sec], 'Keep': -1, 'Params':np.nan , 'Uncs': np.nan}

	if is_TOI:
		#To check that TOI times are marked
		fig_dir = '{}/figures/all_sector'.format(file_dir)
		figname = '{}/{}_{}_allcurve.png'.format(fig_dir, tic, sec)
		plt.clf()
		plt.plot(alltimes_bin, lc_bin, '.')
		for t in times[tic][sec]:
			plt.axvline(0, ymin = 0, ymax = 0.2, color = 'r')
		plt.title("TIC {} in sector {} - Taken from TOI - binned".format(tic, sec))
		plt.savefig(figname)

	elif is_Rand:
		fig_dir = '{}/figures/all_sector'.format(file_dir)
		figname = '{}/{}_{}_allcurve.png'.format(fig_dir, tic, sec)
		plt.clf()
		plt.plot(alltimes_bin, lc_bin, '.')
		for t in times[tic][sec]:
			plt.axvline(0, ymin = 0, ymax = 0.2, color = 'r')
		plt.title("TIC {} in sector {} - Randomly selected - binned".format(tic, sec))
		plt.savefig(figname)



	this_tic={}
	this_tic_unc = {}

	n=1

	if len(times[tic][sec]) == 0 :
		return {'TIC': tic, 'Sector': sec, 'Times': times[tic][sec], 'Keep': 0, 'Params':np.nan , 'Uncs': np.nan}

	for i, t in enumerate(times[tic][sec]):
		# t is in days from start of sector

		[idx] = time_to_index(t)
		bindx = bin_index(idx)

		#Look at small segment around marked time - ensures segment doesn't exceed bounds of data or include next marked transit

		if i>0:
			[prev_idx] = time_to_index(times[tic][sec][i-1])
			prev_idx = bin_index(prev_idx)
		else:
			prev_idx = 0

		if i+1<len(times[tic][sec]):
			[next_idx] = time_to_index(times[tic][sec][i+1])
			next_idx = bin_index(next_idx)
		else:
			next_idx = len(lc_bin)




		all_t_rel = alltimes_bin - (alltimes[0] + t)

		mask = abs(all_t_rel) <1.5

		transit = np.array(lc_bin)[mask]
		transit_times = np.array(all_t_rel)[mask]

		plt.plot(transit_times, transit)
		plt.show()
		t_rel = all_t_rel[mask]
		transit_err = np.array(err_bin)[mask]
		transit_bkg = np.array(bkg_bin)[mask]

		around_mark = lc_bin[bindx - 16 : bindx + 16]

		if np.count_nonzero(~np.isnan(around_mark)) < 16:
			#Return similar but impossible values and asign 0
			continue  

		x0 = [0., 0.5 , 0.5, 0.05, 0., 0., 0., 0.]

	
		soln = fit(x0, transit, t_rel, transit_err) 
	

		for i in range(4):
			tmid = np.random.normal(0., 0.01)
			dur = np.random.normal(0.5, 0.05)
			v = random()
			h = 0.05 * random()
			x1 = [tmid, dur, v, h, 0., 0., 0., np.nanmedian(transit) - 1]

			new_soln = fit(x1, transit, t_rel, transit_err) 

			if new_soln.chisqr <= soln.chisqr and new_soln.success:
				soln = new_soln 


		if not soln.success:

			flat_lc = flatten(alltimes, lc, window_length = 0.5, method = 'biweight', return_trend = False)
			flat_lc_bin = flatten(alltimes_bin, lc_bin, window_length=0.5, method = 'biweight', return_trend=False)

			flat_transit = flat_lc_bin[mask]
			flat_x0 = (0., 0.5, 0.5, 0.05)
			flat_soln = fit(flat_x0, flat_transit, t_rel, transit_err, poly =False)

			if flat_soln.success:
				soln = flat_soln
			else:
				#Failed to fit
				this_tic[t] = (np.nan for i in range(7))
				this_tic_unc[t] = (np.nan for i in range(7))
				continue

		params_fit = unpack_params(soln.params)

		param_uncs= [soln.params[name].stderr for name in soln.params]
		#Note:
		#	param_uncs will return NoneType if the best fit exceeds the allowed range - in this case give all uncertainties as infinity (no constraint on value)


		(tmid, dur, v, h) = params_fit[:4]


		tmid = int(tmid)
		dur = int(dur) if int(dur) != 0 else 1

		t_fit = tmid + alltimes[0] + t


		md_dist = syst.mom_dump_dist(t_fit, qual, alltimes)
		nan_dist = syst.safe_mode_dist(lc, t_fit, alltimes)
		bkg_jitter = syst.background(transit_bkg, t_rel)


		if bkg_jitter == 0:

			figname = "{}/figures/bkg/{}_{}_{}_bkg.png".format(file_dir, tic, sec, n)
			plt.clf()
			plt.plot(transit_times, bkg_bin[tmid - dur//2:tmid+dur//2], '.')

			plt.title("Background {} for TIC {} in sector {}.".format(n, tic, sec))
			plt.savefig(figname)


		centroid = syst.centr_jitter(x_bin, y_bin, all_t_rel,  t)





		t0 = int(tmid - dur//2)
		t3 = int(tmid + dur//2)

		try:
			d_h = param_uncs[3]/h if not param_uncs[3] is None else np.inf
		
		except TypeError:
			print(param_uncs)

		if radii[tic]>0:
			r_p = math.sqrt(h)*radii[tic]
			d_r_p = 0.5 * d_h *r_p
		else:
			#If star radius is unknown - take a larger uncertainty and assume r_star = 1 r_sun. Give the uncertainty as 0.5

			r_p = math.sqrt(h)
			d_r_p = (0.5*d_h + 0.5)*r_p

		
		SNR = syst.SNR(alltimes_bin, lc_bin, h, params_fit)


		u = random()
		if r_p * radii[tic] >= 0.5:
			figname = "{}/figures/{}_{}_{}_{}.png".format(file_dir, tic, sec, n, radii[tic])
			plt.clf()
			plt.plot(transit_times, transit, '.')
			model = gen_model(params_fit, transit_times)
			plt.plot(transit_times, model)
			plt.axvline(x = 0., ymin = 0, ymax = 0.2, color = 'r')
			plt.title("Transit {} for TIC {} in sector {}. R_star = {} - {}".format(n, tic, sec, radii[tic], keep))
			plt.savefig(figname)


		dt = param_uncs[0] if not param_uncs[0] is None else np.inf
		dv = param_uncs[2] if not param_uncs[0] is None else np.inf
 
		n+=1
		this_tic[t] = (v, r_p, md_dist, nan_dist, bkg_jitter, SNR, centroid, n_trans[tic])


		if not param_uncs[0] is None:
			this_tic_unc[t] = (dv, d_r_p, dt, dur, 0.1*bkg_jitter, d_h*SNR, 0.1*centroid, 0)
		else:
			this_tic_unc[t] = (np.inf, np.inf, np.inf, np.inf, np.inf, np.inf, np.inf, np.inf)

	return {'TIC': tic, 'Sector': sec, 'Times': times[tic][sec], 'Keep': keep, 'Params': this_tic, 'Uncs': this_tic_unc}

