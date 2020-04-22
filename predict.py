from utils import *
import pandas as pd 
import numpy as np 
from PRF import prf 
from read_manifest import read_manifest 
from time import time
from os.path import exists
import csv


'''

This script will use a PRF to predict whether to keep unseen light curves.

It requires the training set of features and labels to be stored in the file "manifest.csv" and the features of the light curves to be predicted in the file "unseen.csv".

The TICs, sector numbers, and whether to keep will be stored in "prediction.csv". NB if this file already exists, it will be overwritten.

'''


def train(indir = '', n_tree = 100):
	'''
	This will train the PRF from the test set in indir/manifest.csv.

	By default indir is assumed to be the current working directory and the PRF will have 100 decision trees.
	'''

	if indir != '':
		filename = '{}/manifest.csv'.format(indir)
	else:
		filename = 'manifest.csv'

	assert exists(filename), "The training set could not be found. Please make sure it is saved in the correct folder as 'manifest.csv'"

	X, y0, dX = read_manifest(PRF = True, n_trans = n_trans_bool)

	y = [int(i) for i in np.round(y0)]



	bag_size  = int(len(X) * 0.9)

	test_X = []
	test_dX = []
	test_y = []

	for i in range(len(X) - bag_size):
		idx = randint(0, len(X)-1)
		test_X.append(X.pop(idx))
		test_dX.append(dX.pop(idx))
		test_y.append(y.pop(idx))

	train_X =  np.array(X)
	train_y = np.array(y)
	#train_py = np.array(py)
	train_dX = np.array(dX)
	test_X = np.array(test_X)
	test_y = np.array(test_y)



	test_dX = np.array(test_dX)

	

	start = time()

	TN = 0
	TP = 0
	FN = 0
	FP = 0
	actual_neg = 0
	actual_pos = 0
	pred_neg = 0
	pred_pos = 0reg
	n = len(test_X)


	#dX_train = np.zeros(np.shape(train_X))

	print('No uncertainties')
	clf = prf(n_estimators = n_tree)
	clf.fit(X = train_X, y = train_y) 


	pred_y = clf.predict(X = test_X) #, dX = dX)



	for i, j in enumerate(test_y):
 
		if j == 0:
			actual_neg+=1
			if pred_y[i] == 0:
				pred_neg+=1
				TN +=1
			else:
				pred_pos +=1
				FP += 1
		else:
			actual_pos+=1
			if pred_y[i] == 0:
				pred_neg+=1
				FN +=1
			else:
				pred_pos +=1
				TP += 1




	print(len(test_X) + len(train_X))
	print(clf.score(test_X, test_y))

	print("PRF trained successfully on {} light curves. {} decision trees were used. Took {:.2f} seconds.".format(len(test_X), n_tree, time()-start))

			 
	matrix = [
		["n = {}".format(n), "Predicted positive = {}" . format(pred_pos), "Predicted negative = {}".format(pred_neg)],
		["Actual positive = {}".format(actual_pos), "TP = {}: {:.2f}%".format(TP, TP*100/actual_pos), "FN = {}: {:.2f}%".format(FN, FN*100/actual_pos)],
		["Actual negative = {}".format(actual_neg), "FP = {}: {:.2f}%".format(FP, FP*100/actual_neg), "TN = {}: {:.2f}%".format(TN, TN*100/actual_neg)]
	]

	s = [[str(e) for e in row] for row in matrix]
	lens = [max(map(len, col)) for col in zip(*s)]
	fmt = '\t'.join('{{:{}}}'.format(x) for x in lens)
	table = [fmt.format(*row) for row in s]
	print('\n'.join(table))


	print('\nReady to predict\n')

	return clf

def predict(clf, indir = ''):
	if indir != '':
		filename = '{}/unseen.csv'
	else:
		filename = 'unseen.csv'

	assert exists(filename), "Could not find features to be used for prediction. Please make sure they are saved as 'uneen.csv'"

	

	pred_y = {}

	for sec in X:
		if not sec in pred_y:
			pred_y[sec] = {}

			for tic in X[sec]:
				pred_y[sec][tic] = clf.predict(X = X[sec][tic][0])

	return pred_y


if __name__ == 'main':
	
	clf = train()
	pred_y = predict(clf)

	filename = 'prediction.csv'

	with open(filename, 'w', newline = '') as file:
		writer = csv.writer(file)
		writer.writerow(['TIC', 'Sector', 'Keep (1) or reject (0)?'])

	for sec in pred_y:
		for tic in pred_y[sec]:

			with open(filename, 'a', newline = '') as file:
				writer = csv.writer(file)
				writer.writerow([tic, sec, pred_y[sec][tic]])

	print("Predictions complete. Saved to 'prediction.csv'")

