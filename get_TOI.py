import pandas as pd
import numpy as np
from astropy.time import Time
from random import randint
import csv
from utils import *
import os
import math


def read_TOIs():
	#Reads list of TOIs to put in same format as PHT
	if os.name == 'nt':
		dir ="C:/Users/blake"
		dirname = "{}/{}".format(dir, "MPhys")
		filename_TOIs = "{}/{}".format(dirname, "TOI_list.csv")
	elif os.name == 'posix':
		dir ="/mnt/users/kepler/kepler2/TESS/planethunters"
		dirname = "/mnt/zfsusers/blakeland/Documents/MPhys"
		filename_TOIs = "{}/{}".format(dirname, "TOI_list.csv")
	
	TOI_file = pd.read_csv(filename_TOIs)

	TICs = list(TOI_file['TIC ID'])
	sectors_list = list(TOI_file['Sectors'])
	epochs_list = list(TOI_file['Epoch (BJD)'])
	periods_list = list(TOI_file['Period (days)'])


	rads = TOI_file['Stellar Radius (R_Sun)']

	#sectors = pd.DataFrame(sectors_list, index = TICs).to_dict()[0]
	radii = {TICs[i]: rads[i] for i in range(len(TICs))}
	#print(radii)
	sectors = {}

	n_1=0

	for i, tic in enumerate(TICs):

		sec = sectors_list[i]
		if isinstance(sec, str):
			lst = [int(a) for a in sec.split(',')]
			n_1+=1

		sectors[tic] = lst 





	epochs = {TICs[i]:float(epochs_list[i]) for i in range(len(TICs))}
	periods = {TICs[i]:float(periods_list[i]) for i in range(len(TICs))}

	BJD_to_TESSdate = -2457000

	


	out_dict = {sec:{} for sec in range(1,50)} # Allow for many sectors

	n=0
	for tic in TICs:

		t0 = epochs[tic] + BJD_to_TESSdate
		p = periods[tic]

		if p == 0:
			continue
		
		if np.isnan(radii[tic]):
			radius = -99

		else:
			radius = radii[tic]

		for sec in sectors[tic]:
			try:
				alltimes, lc, err, qual, bkg, x, y = read_file(tic, sec, bin = False)
			except:
				continue

			start = alltimes[0]
			end = alltimes[len(alltimes)-1] 


			all_transits = np.array([t0 + p*i for i in range(2000)])


			mask = (all_transits >= start) & (all_transits <= end )
			

			dt = list(all_transits[mask] - start)
			lc = np.array(lc)

			dt_out=[]
			n_times = 0

			while n_times < 6 :
				if dt:
					u = randint(0, len(dt)-1)

					test_time = dt[u]

					mask = abs(np.array(alltimes) - test_time) <= 0.25

					if np.count_nonzero(np.isnan(lc[mask]))<0.5 * len(lc[mask]):
						dt_out.append(dt.pop(u)) 
					else:
						dt.pop(u)

					n_times+=1 
				else:
					n_times = 6 






			else:
				dt_out = dt		
			


			out_dict[sec][tic] = [tic, sec, dt_out, radius, 1 , 'TOI', len(dt)]



		
	
	print("TOIs done")


			
	return out_dict


