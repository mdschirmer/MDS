#!/usr/bin/env python
"""
:Summary:
Program to create subsets of a data sets that are likely to be not statistically significantly different.

:Description:
Subsets are picked based on the psuedo gower distance, where the distance between subjects picked in each iteration is minimized.

As input the code requires 
	a) phenotype file (only numerical values, i.e. please convert information such as sex (M/F) into 0/1)
	b) phenotype list to be used for comparison (see example phenotype_list.csv file)
	c) Number of subsets
	
Output
CSV file with added column including subset IDs.

Execution
./subset_data.py INFILE OUTFILE NUM_GROUPS (PHENOTYPE_LIST)

Example
./subset_data.py phenotypes.csv phenotypes_split.csv 5 phenotype_list.csv
 

:Requires:
numpy
scipy

:TODO:
 - add statistical tests into this code
 - simplify the import statements

:AUTHOR: MDS
:ORGANIZATION: MGH/HMS
:CONTACT: software@markus-schirmer.com
:SINCE: 2019-05-15
:VERSION: 0.1
"""
#=============================================
# Metadata
#=============================================
__author__ = 'mds'
__contact__ = 'software@markus-schirmer.com'
__date__ = '2019-05'
__version__ = '0.1'

#=============================================
# Import statements
#=============================================
import sys
import numpy as np
import scipy.spatial.distance as ssd
import scipy.optimize
import csv
import scipy.optimize as so

import random

import pdb
#=============================================
# Helper functions
#=============================================
def intersection(lst1, lst2):
    return list(set(lst1) & set(lst2))

def get_pseudo_gower_distance(X):
	# calculate gower distance
	# X needs to be num_subject x num_features

	# need to normalize
	maximums = np.max(X,axis=0)
	minimums = np.min(X,axis=0)

	# This is to normalize the numeric values between 0 and 1.
	X_norm = np.divide(X - minimums[None,:] ,maximums,out=np.zeros_like(X), where=maximums!=0)

	# Calculate distances
	distance = ssd.squareform(ssd.pdist(X_norm,'minkowski', p=1))/X_norm.shape[1]

	return distance


def match_subjects(subject_list, comparison_list=None , keys=None, ngroups=2):

	assert isinstance(subject_list, list),'Subject list should be an instance of list.'
	
	# handle comparison between two lists
	if comparison_list is not None:
		assert isinstance(subject_list, list),'Subject list should be an instance of list.'
		two_list = True
		num_sub_in_subject_list = len(subject_list)
		for sub in comparison_list:
			subject_list.append(sub)
	else:
		two_list=False

	# find the list of keys to compare to
	if keys is None:
		# get all the keys in 
		keys = []
		for sub in subject_list:
			if not keys:
				keys = sub.keys()
			else:
				keys = intersection(keys, sub.keys())

	# get the data
	X = np.asarray([[sub[key] for key in keys] for sub in subject_list])

	assert (X.shape[0]==len(subject_list)),'Issue with data array. Wrong dimensions.'

	# calculate distance

	if two_list:
		# pdb.set_trace()
		cost = get_pseudo_gower_distance(X)[:num_sub_in_subject_list,num_sub_in_subject_list:]
	else:
		dist = get_pseudo_gower_distance(X)
		# set diagonal elements to inf as pseudo cost
		cost = dist + np.where(np.eye(dist.shape[0])>0, np.inf, 0)

	# gather output 
	output = []
	for ii in np.arange(ngroups):
		output.append([])
	keeping_count = []
	loop_idx = np.arange(cost.shape[0])
	np.random.shuffle(loop_idx)
	for ii in loop_idx:
		if ii not in keeping_count:
			output[0].append(ii)
			keeping_count.append(ii)

			possible_jj = np.argsort(cost[ii,:])
			# randomly permute groups
			groups = np.arange(1, ngroups)
			np.random.shuffle(groups)
			for group in groups:
				idx = 0
				while idx<cost.shape[1]:
					jj = possible_jj[idx]
					if jj not in keeping_count:
						output[group].append(jj)
						keeping_count.append(jj)
						idx = np.inf
					idx = idx+1

	return output


def load_csv(file):
	with open(file,'r') as fid:
		reader = csv.reader(fid)
		header = next(reader)

		phenos = {}
		for idx, row in enumerate(reader):
			phenos[idx] = {}
			for ii, item in enumerate(row):
				try: 
					phenos[idx][header[ii]] = float(item)
				except:
					phenos[idx][header[ii]] = item
	return phenos

def load_key_list(file):
	with open(file,'r') as fid:
		reader = csv.reader(fid)
		key_list = []
		for idx, row in enumerate(reader):
			for clm in row:
				key_list.append(clm)

	return key_list

def shuffle_forward(l):
    order = range(len(l)); random.shuffle(order)
    return list(np.array(l)[order]), order

def shuffle_backward(l, order):
    l_out = [0] * len(l)
    for i, j in enumerate(order):
        l_out[j] = l[i]
    return l_out

#=============================================
# Main method
#=============================================
def main(argv):

	if len(argv)==1:
		print_help()

	## catch input
	infile=argv[1] # phenotype information file (csv)
	outfile=argv[2] # output file
	num_groups = int(argv[3])
	if len(argv)<5:
		key_file = 'phenotype_list.csv'

	# load phenotype list and data
	key_list = load_key_list(key_file)
	data = load_csv(infile)

	# initialize lists
	subject_list = [data[ii] for ii in sorted(data.keys())]
	comparison_list=None

	# create subsets
	output = match_subjects(subject_list, comparison_list=comparison_list, keys=key_list, ngroups=num_groups)

	# write output
	with open(outfile,'w') as fid:
		fid.write('Group,')
		for key in data[list(data.keys())[0]].keys():
			fid.write('%s,' %key)
		fid.write('\n')
		for ii, group in enumerate(output):
			for sub in group:
				fid.write('Set_%i,' %ii)
				for key in data[sub].keys():
					fid.write('%s,' %(data[sub][key]))
				fid.write('\n')

if __name__ == "__main__":
	main(sys.argv)
