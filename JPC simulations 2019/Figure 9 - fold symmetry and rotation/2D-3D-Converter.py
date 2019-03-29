#-------------------------------------------------------------------------------
# Name:        2D-3D Converter
#
# Author:      Mark Tingey
#
# Created:     05/02/2019
#-----------------------------------------------------------------------------

import math
import random
import csv
import os
from numpy.random import choice
import numpy as np
from scipy.optimize import curve_fit
import time
import matplotlib.pyplot as plt
from scipy import stats
import pandas as pd
import pylab as pl

def make_hist(data, r, d_r):

			hist_values, bin_edges = np.histogram(data, bins = 2 * int(r/d_r), range = (-r, r))
			new_bin_edges = []
			for i in bin_edges:
				if i >= 0:
					new_bin_edges.append(i)

			new_hist_values = smoothdata(hist_values, r, d_r)

			return new_hist_values, new_bin_edges

def smoothdata(data, r, d_r):
			"""smoothds data with 3 moving window and takes abs value average"""
			smooth_data = []
			r += 1

			##comment out for smoothing
			smooth_data = []
			for i in range(len(data)):
				smooth_data.append(data[i])

			##adds + and - bins
			final_smooth_data = []
			for i in range(int(r/d_r)):
				final_smooth_data.append(smooth_data[i] + smooth_data[len(smooth_data)-1-i])

			return list(reversed(final_smooth_data))

def deconvolution(hv, be, r, d_r):
			"""hv = hist_values, be = bin_edges"""

			density = []
			matrix = gen_matrix(r, d_r)

			while len(hv) > len(matrix):
				hv.pop()
			while len(matrix) > len(hv):
				matrix.pop()

			rev_hv = list(reversed(hv))
			x = 0
			for i in range(len(rev_hv)):
				##calculate how much to subtract from bin
				density_sub = 0
				y = 0
				for j in range(x):
					density_sub += density[y] * matrix[j][i]
					y += 1

				##calculate final bin value
				density.append((rev_hv[i] - density_sub) / matrix[i][i])
				x += 1

			unrev_hv = list(reversed(density))

			smooth_data = []
			for i in range(len(unrev_hv)):
				if i == 0 or i == (len(unrev_hv) - 1):
					smooth_data.append(unrev_hv[i])
				else:
					smooth_data.append(np.average([unrev_hv[i-1], unrev_hv[i], unrev_hv[i+1]]))

			return unrev_hv, smooth_data, hv


def gen_matrix(r, d_r):

			##'be' is short for bin edges
			if r%d_r > 0:
				be = range(0, r+r%d_r, d_r)
			else:
				be = range(0, r+d_r, d_r)

			matrix = []
			for i in range(len(be)-1):
				matrix.append([])

			x = 0
			for i in range(len(matrix)):
				for j in range(x):
					matrix[i].append(0)
				x += 1

			##generate areas of sections closest to x axis
			for i in range(len(matrix)):
				theta = np.arccos(float(be[len(be)-2-i])/float(be[len(be)-1-i]))
				arc_area = (theta/(2*np.pi)) * np.pi * float(be[len(be)-1-i])**2
				tri_area = 0.5 * float(be[len(be)-2-i]) * (np.sin(theta) * float(be[len(be)-1-i]))
				matrix[i].append(4 * (arc_area - tri_area))

			##skipping factor
			x = 2
			##generate areas of layers going further out from x axis
			while len(matrix[0]) < len(matrix):
				for i in range(len(matrix) - len(matrix[0])):
					num = 0
					for j in range(len(matrix)):
						for k in range(len(matrix[i]) + 1):
							if j == i and k < len(matrix[i]):
								num += matrix[j][k]
							elif j > i:
								num += matrix[j][k]
					theta = np.arccos(float(be[len(be)-1-x-i])/float(be[len(be)-1-i]))
					arc_area = (theta/(2*np.pi)) * np.pi * float(be[len(be)-1-i])**2
					tri_area = 0.5 * float(be[len(be)-1-x-i]) * (np.sin(theta) * float(be[len(be)-1-i]))
					matrix[i].append(4 * (arc_area - tri_area) - num)
				x += 1

			return matrix

##global variables

cwd = os.getcwd()
nwd = cwd+"\\data set"
os.chdir(nwd)
binsize = 10
filename = "simulation0.csv"

##read data file into dataframe
xy = pd.read_csv(filename, header = None, usecols=[0,1])
xcoord = xy[0].tolist()
ycoord = xy[1].tolist()

##set radius
minrad = np.min(ycoord)
maxrad = np.max(ycoord)
radius = int(((maxrad-minrad)/2)+(0.5*((maxrad - minrad)/2)))

##run conversion
a,b = make_hist(ycoord, radius, binsize)
final = deconvolution(a,b,radius,binsize)
export = final[0]
neglist = [x for x in export]

##convert to list
export = [float(np_float) for np_float in export]
del export[-1]
neglist = [float(np_float) for np_float in neglist]
del neglist[-1]
neglist.reverse()

##make list of bins
midbin = (binsize/2)
posbins = pl.frange(midbin,radius,binsize)
posbins = posbins.tolist()
negbins = pl.frange(-midbin,-radius,-binsize)
negbins = np.flip(negbins)
negbins = negbins.tolist()
binlist = negbins + posbins
del binlist[-1]
del binlist[:1]
final_hist = neglist + export

##ensure bins are same length as histagram
while len(binlist) != len(final_hist):
    if len(binlist) > len(final_hist):
        del binlist[-1]
        del binlist[:1]
    if len(binlist) < len(final_hist):
        [0]+binlist
        binlist.append(0)



outfile = filename[:-4] + "-converted" + "-bin-" + str(binsize) + ".csv"

##write export histogram to csv
with open(outfile, 'w', newline='') as writefile:
    writer = csv.writer(writefile)
    writer.writerows(zip(binlist, final_hist))
writefile.close()


##Export final histogram