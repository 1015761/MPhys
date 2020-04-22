#This script is intended to get all the various inputs (TOIs, CTC etc. ) in the same format
import pandas as pd
import csv
import os
import numpy as np
from get_TOI import read_TOIs
from random import choice

def time_to_index(time, cadence = 2, min_time = 0.008333):
	#Turns a BJD time into index in the data
	interval = cadence / 1440
	index=[]
	if not isinstance(time, list):
		time = [time]

	for t in time:
		index.append(int((t-min_time)/interval))
	return index

def keep_cond(label):
	return True if int(label)==1 else False

def read_PHT(sec, mode = 'Time', with_laels = True):

	'''
	Reads PHT output for a given sector
	
	'''


	if os.name == 'nt':
		indir = 'C:/Users/blake/MPhys'
	elif os.name == 'posix':
		indir = '/mnt/zfsusers/blakeland/Documents/MPhys'

	filename_pht = '{}_Data/PHT_output_{}.csv'.format(indir, sec)
	filename_CTC = '{}/CTC_output_{}.csv'.format(indir, sec)

	PHT = pd.read_csv(filename_pht, error_bad_lines = False)
	CTC = pd.read_csv(filename_CTC, error_bad_lines = False)

	all_tics = PHT.TIC_ID
	times = PHT.db_peak
	
	df = pd.DataFrame(list(times), index = all_tics).to_dict()



	TICs = CTC.TIC_ID
	radii = CTC.rad
	TOIs = CTC.TOI
	assert len(TOIs) == len(TICs)
	keep_list = CTC.p_Keep
	#print(keep_list)

	#is_TOI = pd.DataFrame(TOIs, index = TICs)
	is_TOI = {}

	for i, tic in enumerate(TICs):
		is_TOI[tic] = True if not np.isnan(TOIs[i]) else False

	label = {TICs[i]: keep_list[i] if not np.isnan(keep_list[i]) else 0 for i in range(len(TICs))}
	#print(label)
	#	label = pd.DataFrame(keep_list, index = TICs)['p_Keep'].to_dict()

	times = []
	indices = []
	output={}

	for tic in TICs:
		if tic != -999:
			times_list = [float(j) for j in df[0][tic][1:-1].split(',')]
			times.append(times_list)
			indices_list = time_to_index(times_list, np.mean(PHT.min_time))
			indices.append(indices_list)
			if mode == 'index':
				output[tic] = indices_list
			else:
				output[tic] = times_list


	tic_radius={}

	for i, tic in enumerate(TICs):
		tic_radius[tic] = radii[i] 

	this_sec = {}

	for tic in TICs:
		if not is_TOI[tic]:
			if with_label:
				this_sec[tic] = [tic, sec, output[tic], tic_radius[tic], label[tic], "CTC", len(output[tic])]
			else:
				this_sec[tic] = [tic, sec, output[tic], tic_radius[tic], len(output[tic])]
	print("CTC sector {} done".format(sec))
	return this_sec






def write_CSV(sec_list, TOIs = True, balance = True, with_label = True):

	'''
	Writes TICs.csv - the list of TICS to be analysed
	'''

	#TODO: Look at sec_list for lcs to be predicted - how to input?

	if with_label:

		if os.name == 'nt':
			indir = 'C:/Users/blake/MPhys'
		elif os.name == 'posix':
			indir = '/mnt/zfsusers/blakeland/Documents/MPhys'

		filename = "{}/TICs.csv".format(indir)
		if os.name == 'nt':
			with open(filename, 'w', newline = '') as file:
				writer = csv.writer(file)
				writer.writerow(['ID', 'TIC', 'Sector', 'Times', 'Radius', 'Label', 'Mode', 'N_transits'])
		elif os.name=='posix':
			with open(filename, 'w') as file:
				writer = csv.writer(file)
				writer.writerow(['ID', 'TIC', 'Sector', 'Times', 'Radius', 'Label', 'Mode', 'N_transits'])
		



		sec_tics = {}
		if os.name=='nt':
			n_tic=0
			n_keep=0
			n_trans=0

			for sec in sec_list:
				sec_tics[sec] = read_PHT(sec)
			for sec in sec_tics:
				for tic in sec_tics[sec]:
					n_tic += 1

					#Number of marked times
					n_trans += len(sec_tics[sec][tic][2])

					if keep_cond(sec_tics[sec][tic][4]):
						n_keep += len(sec_tics[sec][tic][2])

					with open(filename, 'a', newline = '') as file:
						writer = csv.writer(file)
						writer.writerow([n_tic, *sec_tics[sec][tic]])




			if TOIs:
				TOI_tics = read_TOIs()
				for sec in range(1, 19):
					for tic in TOI_tics[sec]:

						n_tic+=1
						n_trans += len(TOI_tics[sec][tic][2])
						n_keep += len(TOI_tics[sec][tic][2])

						with open(filename, 'a', newline = '') as file:
							writer = csv.writer(file)
							writer.writerow([n_tic, *TOI_tics[sec][tic]])

		elif os.name == 'posix':
			n_tic = 0
			n_keep = 0
			n_trans = 0
			for sec in sec_list:
				sec_tics[sec] = read_PHT(sec)

				for sec in sec_tics:
					for tic in sec_tics[sec]:

						n_tic += 1
						n_trans += len(sec_tics[sec][tic][2])

						if keep_cond(sec_tics[sec][tic][4]):
							n_keep += len(sec_tics[sec][tic][2])

						with open(filename, 'a') as file:
							writer = csv.writer(file)
							writer.writerow([n_tic, *sec_tics[sec][tic]])

							#n_trans += len(sec_tics[sec][tic][2])


			if TOIs:
				TOI_tics = read_TOIs()
				#print(TOI_tics)
				for sec in range(1, 39):
					#try:
					for tic in TOI_tics[sec]:
						n_tic+=1
						n_trans += len(TOI_tics[sec][tic][2])
						n_keep += len(TOI_tics[sec][tic][2])

						with open(filename, 'a') as file:
							writer = csv.writer(file)
							writer.writerow([n_tic, *TOI_tics[sec][tic]])
	#				except KeyError:
	#					print("Couldn't find TOIs for sector {}".format(sec))

		if balance:
			#Balance data set with too many 'planets' by taking random transits
			n_rand =0
			#I want to first check that there are more positives than negatives

			print(n_keep)
			print(n_trans)


			while n_keep >= 0.5*n_trans:
				#If so, take random tic (that has marked transits)

				#Constrain to chosen sectors
				

				rand_sec = choice(sec_list)

				filename_PHT = '{}_Data/PHT_output_{}.csv'.format(indir, rand_sec)

				PHT_file = pd.read_csv(filename_PHT, error_bad_lines = False)
				

				TIC_list = np.array(PHT_file['TIC_ID'])

				times_list = np.array(PHT_file['db_peak'])

				#print(times_list)

				TIC_times = pd.DataFrame(times_list, index = TIC_list)

				this_sec = list(sec_tics[rand_sec].keys())



				if TOIs:
					all_TOIs = list(TOI_tics[rand_sec].keys())
					mask = (TIC_list>0) & (times_list != '[0]') & (~np.isin(TIC_list, this_sec)) & (~np.isin(TIC_list, all_TOIs)) # & (not np.isnan(times_list))
				else:
					mask = (TIC_list>0) & (times_list != '[0]') & (~np.isin(TIC_list, this_sec))# & (~np.isnan(times_list))


				good_tics = list(TIC_list[mask])
				#print(len(good_tics))
				#print(np.array(good_tics))

				if len(good_tics)==0:
					
					continue



				rand_tic = choice(good_tics)

				rand_times = [ float(i) for i in TIC_times[0][rand_tic][1:-1].split(',')]


				if os.name=='posix':
					with open(filename, 'a') as file:
						writer = csv.writer(file)
						writer.writerow([n_tic, rand_tic, rand_sec, rand_times, -99, 0, 'Random'])
				elif os.name=='nt':
					with open(filename, 'a', newline = '') as file:
						writer = csv.writer(file)
						writer.writerow([n_tic, rand_tic, rand_sec, rand_times, -99, 0, 'Random'])
						
				n_tic+=1
				n_rand += 1
				n_trans += len(rand_times)

				if n_rand%100==0:
					print("{} random TICs added".format(n_rand))



			print('Chosen {} random TICs from sectors {} to balance the training set'.format(n_rand, sec_list))



	else:

		if os.name == 'nt':
			indir = 'C:/Users/blake/MPhys'
		elif os.name == 'posix':
			indir = '/mnt/zfsusers/blakeland/Documents/MPhys'

		filename = "{}/TICs_no_label.csv".format(indir)
		if os.name == 'nt':
			with open(filename, 'w', newline = '') as file:
				writer = csv.writer(file)
				writer.writerow(['ID', 'TIC', 'Sector', 'Times', 'Radius', 'N_transits'])
		elif os.name=='posix':
			with open(filename, 'w') as file:
				writer = csv.writer(file)
				writer.writerow(['ID', 'TIC', 'Sector', 'Times', 'Radius', 'N_transits'])
		


		sec_tics = {}
		if os.name=='nt':
			n_tic=0
			n_keep=0
			n_trans=0

			for sec in sec_list:
				sec_tics[sec] = read_PHT(sec)
			for sec in sec_tics:
				for tic in sec_tics[sec]:
					n_tic += 1

					#Number of marked times
					n_trans += len(sec_tics[sec][tic][2])

					if keep_cond(sec_tics[sec][tic][4]):
						n_keep += len(sec_tics[sec][tic][2])

					with open(filename, 'a', newline = '') as file:
						writer = csv.writer(file)
						writer.writerow([n_tic, *sec_tics[sec][tic]])




			if TOIs:
				TOI_tics = read_TOIs()
				for sec in range(1, 19):
					for tic in TOI_tics[sec]:

						n_tic+=1
						n_trans += len(TOI_tics[sec][tic][2])
						n_keep += len(TOI_tics[sec][tic][2])

						with open(filename, 'a', newline = '') as file:
							writer = csv.writer(file)
							writer.writerow([n_tic, *TOI_tics[sec][tic]])

		elif os.name == 'posix':
			n_tic = 0
			n_keep = 0
			n_trans = 0
			for sec in sec_list:
				sec_tics[sec] = read_PHT(sec)

				for sec in sec_tics:
					for tic in sec_tics[sec]:

						n_tic += 1
						n_trans += len(sec_tics[sec][tic][2])

						if keep_cond(sec_tics[sec][tic][4]):
							n_keep += len(sec_tics[sec][tic][2])

						with open(filename, 'a') as file:
							writer = csv.writer(file)
							writer.writerow([n_tic, *sec_tics[sec][tic]])

							#n_trans += len(sec_tics[sec][tic][2])









write_CSV([12, 13, 14, 15], TOIs = True, balance = False)