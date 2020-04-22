import pandas as pd 
import math
import numpy as np

def str_to_dict(str):

	#splits "{a: (b,c), d: (e,f)}" to ['a', '(b,c)', 'd', '(e,f)']    

	#Either {<time> : -99, ...} or {<time>, (<params>), ...}
	str = str[1:-1].replace('),', '):')
	str = str.replace(': -99,', ': -99:')
	str = str.replace(': -101,', ': -101:')
	str = str.replace(' (', '(')
	items = str.split(':')



	assert(len(items)%2==0)
	out_dict = {}
	for i in range(0, len(items), 2):
		if items[i+1] == ' -99':
			out_dict[float(items[i])] = -99
		elif items[i+1] == ' -101':
			out_dict[float(items[i])] = -101
		else:
			out_dict[float(items[i])] = str_to_list(items[i+1])
	return out_dict


def str_to_list(str, int_mode = False):
	lst = []
	for i in str[1:-1].split(','):
		try:
			lst.append(float(i))
		except:
			lst.append(np.nan)
	df = pd.DataFrame(lst)
	#df.fillna(-1e10, inplace = True)
	lst = df[0].values.tolist()


	
	if int_mode:
		lst = [int(i) for i in lst]

	return lst


def read_manifest(PRF = False, n_trans = False):

	
	if not PRF:
		filename = 'C://Users/blake/MPhys/PRF_manifest.csv'
		file = pd.read_csv(filename, error_bad_lines = False)
		labels = file['Keep']
		params_df = file['Params']
		X=[]
		y=[]
		n99=0
		for i, item in enumerate(params_df):
			try:
				this_tic = str_to_dict(item)
				
				for t in this_tic:
					params = this_tic[t]
					if  isinstance(params, list):
					#if params[0]<=1 and params[0]>=0:
						if n_trans:
							X.append(params)
						else:
							X.append(params[:-1])
					
					#Change this to reflect threshold 
					#this_y = 1 if labels[i]>0.4 else 0
						this_y = labels[i]

						y.append(this_y)

					else:
						if params == -99:
							n99+=1
						X.append((np.nan for i in range(7)))

						y.append(0)
			except:
				pass
		return X, y


	else:
		filename = 'C://Users/blake/MPhys/PRF_manifest.csv'
		file = pd.read_csv(filename, error_bad_lines = False)
		labels = file['Keep']
		params_df = file['Params']
		uncs_df = file['Uncs']
		X=[]
		y=[]
		dX = []
		n99=0
		for i, item in enumerate(params_df):
			try:
				this_tic = str_to_dict(item)
				this_tic_unc = str_to_dict(uncs_df[i])
				
				for t in this_tic:
					params = this_tic[t]
					uncs = this_tic_unc[t]
					if  isinstance(params, list):
					#if params[0]<=1 and params[0]>=0:
						if n_trans:
							X.append(params)
							dX.append(uncs)
						else:

							X.append(params[:-1])
							dX.append(uncs[:-1])
					
					#Change this to reflect threshold 
					#this_y = 1 if labels[i]>0.4 else 0
						this_y = labels[i]

						y.append(this_y)

					else:
						if params == -99:
							n99+=1
						X.append((-0.5, -0.5, -100000, -1000000, -0.01, -0.01, -0.001))
						dX.append((0., 0., 0., 0., 0., 0., 0.))

						y.append(0)
			except:
				pass

		#print(n99)
	
		return X, y, dX

def read_unseen():
	#TODO: Adapt writing csv code
	#TODO: Add MPI code for prediction code

	filename = 'C://Users/blake/MPhys/PRF_manifest.csv'
	file = pd.read_csv(filename, error_bad_lines = False)
	labels = file['Keep']
	params_df = file['Params']
	sectors = file['Sector']
	TICs = file['TIC_ID']
	X={}
	n99=0
	for i, item in enumerate(params_df):

		sec = sectors[i]
		tic = TICs[i]

		if not sec in X:
			X[sec] = {}

		try:
			this_tic = str_to_dict(item)
			for t in this_tic:
				params = this_tic[t]
				if  isinstance(params, list):
					if n_trans:
						X[sec][tic] = params
					else:
						X[sec][tic] = params[:-1]

				else:
					if params == -99:
						n99+=1
					X[sec][tic] = (np.nan for i in range(7))
		except:
			pass
	return X

