# -*- coding: utf-8 -*-
"""
Created on Thu Apr 20 10:52:38 2017

@author: Andrew Ruba

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.
    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.
    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""

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

## below arrays are for saving simulation data for statistical analysis
global gausslist
gausslist = []
global bimodallist
bimodallist = []
global bimodalmean
bimodalmean = []
global bimodalsd
bimodalsd = []
global bimodalheight
bimodalheight = []
global bimodalauc
bimodalauc = []

def sim(gui, PTNUM, RADIUS, PREC, ITER, BINSIZE, PERCERROR, A, ROTATION):
	def simulation(num_points, radius, dr, ss, mm, aa, rotation):
		
		def area_fn(X):
			X = float(X)
			A = -(dr**2)*np.pi
			B = dr*2*np.pi
			return X*B+A
			
		def gauss_fn(x, s, m):
			a = area_fn(m)
			x = float(x)
			s = float(s)
			m = float(m)
			return a*np.e**(-(x-m)**2.0/(2.0*s**2.0))
			
		def combine(x):
			s = ss
			m = mm
			return (area_fn(x) * gauss_fn(x, s, m))
		
		##starting with perfect x,y and adding error
		xydata = []
		mm = mm + 0.00001
		while len(xydata) < num_points:
			theta = np.radians((np.random.random()*360.0*aa)-0.5*(1.0-aa)*360.0)
			## precision distribution sampling
#             ss = choice([3,5,7,9], p=[0.1475,0.2775,0.3075,0.2675])
#             ss = choice([4.5,5.5,6.5,7.5,8.5,9.5], p=[0.02,0.05,0.07,0.11,0.2,0.55])
			y_prec = np.random.normal(0.0, ss)
			z_prec = np.random.normal(0.0, ss)
			yy = mm*np.cos(theta)
			zz = mm*np.sin(theta)
			yyy = yy*np.cos(np.radians(rotation))+zz*np.sin(np.radians(rotation))
			zzz = -yy*np.sin(np.radians(rotation))+zz*np.cos(np.radians(rotation))
			xydata.append((yyy+y_prec, zzz+z_prec))
	
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
			
		def make_hist(data, r, d_r):
			
			hist_values, bin_edges = np.histogram(data, bins = 2 * int(r/d_r), range = (-r, r))
			
			new_bin_edges = []    
			for i in bin_edges:
				if i >= 0:
					new_bin_edges.append(i)
					
			new_hist_values = smoothdata(hist_values, r, d_r)
			  
			return new_hist_values, new_bin_edges
			
		def csv_read(path):
			with open(path, 'rb') as csvfile:
				reader = csv.reader(csvfile, delimiter = ',')
				holdlist = []
				for row in reader:
					holdlist.append(float(row[1]))
			return holdlist
	
		jkl = []
		for y,z in xydata:
			jkl.append(y)
			
		radius = int(np.floor(radius/dr))*dr
		
		if num_points == PTNUM + 1:
			## decide the proper bin size
			minbinsize = 2
			binsizes = []
			binsizesdata = [[] for variable in range(1, int(PREC)+1)]
			gui.message.set('0% done calculating ideal bin size...')
			gui.update()
			for binoptimization in range(10):
				for binsize in range(1, int(PREC)+1):
					if binsize >= minbinsize:
						error = 0
						# print ('binsize ' + str(binsize))
						jkl = []
						mm = mm + 0.00001
						while len(jkl) < num_points-1:
							theta = np.random.random()*360.0
							## precision distribution sampling
#                             ss = choice([3,5,7,9], p=[0.1475,0.2775,0.3075,0.2675])
#                             ss = choice([4.5,5.5,6.5,7.5,8.5,9.5], p=[0.02,0.05,0.07,0.11,0.2,0.55])
							y_prec = np.random.normal(0.0, ss)
							jkl.append(mm*np.cos(theta)+y_prec)
						a,b = make_hist(jkl, radius, binsize)
						final_unsmooth, final_smooth, final_2d = deconvolution(a, b, radius, binsize)
						holdlist = []
						addZero = False
						for val in list(reversed(final_unsmooth)):
							if not addZero:
								if val >= 0.0:
									holdlist.append(val)
								else:
									addZero = True
									holdlist.append(0.0)
							else:
								holdlist.append(0.0)
						final_unsmooth = list(reversed(holdlist))
						##rescale ideal data
						matrix = gen_matrix(radius, binsize)
						newmatrix = []
						for i in matrix:
							newmatrix.append(list(reversed(i)))
						matrix = list(reversed(newmatrix))
						# print (a)
						# print (final_unsmooth)
						while len(a) > len(matrix):
							a.pop()
						while len(matrix) > len(a):
							matrix.pop()
						
						for ncol in range(len(matrix[0])):
							binsub = 0.0
							for mcol in range(len(matrix)):
								binsub += float(matrix[mcol][ncol]*final_unsmooth[mcol])
							try:
								if a[ncol] != 0.0:
									# print (binsub)
									error += np.square(a[ncol] - binsub) / a[ncol]
							except:
								pass
								
						popped = a.pop()
						while popped == 0:
							popped = a.pop()
						binsizesdata[binsize-1].append((error, len(a)+1,1-stats.chi2.cdf(error, len(a)+1),binsize))
					else:
						binsizesdata[binsize-1].append((1000000.0,1,0.0,binsize))
					gui.message.set(str((binoptimization*10) + 10) + ' % done calculating ideal bin size...')
					gui.update()
			
			finalbinsizes = []
			for bintrial in range(len(binsizesdata)):
				errhold = []
				dfhold = []
				pvalhold = []
				binhold = []
				for trial in range(len(binsizesdata[bintrial])):
					chisq, df, pval, binsize = binsizesdata[bintrial][trial]
					errhold.append(chisq)
					dfhold.append(df)
					pvalhold.append(pval)
					binhold.append(binsize)
				chisq = np.average(errhold)
				df = np.round(np.average(dfhold))
				pval = 1-stats.chi2.cdf(chisq,df)
				binsize = binhold[0]
				finalbinsizes.append((chisq,df,pval,binsize))
			# print (finalbinsizes)
			for binsizedata in finalbinsizes:
				chisq, df, pval, binsize = binsizedata
				if pval >= 0.95:
					dr = binsize
					break
				else:
					dr = int(PREC)

		a,b = make_hist(jkl, radius, dr)
		final = deconvolution(a,b,radius,dr)
		
		if num_points != PTNUM + 1:
			def gauss_fn(x, a, s, m):
				return a*np.e**(-(x-m)**2.0/(2.0*s**2.0))
			
			def bimodal(x,mu1,sigma1,A1,mu2,sigma2,A2):
				return gauss_fn(x, A1, sigma1, mu1)+gauss_fn(x, A2, sigma2, mu2)
			
			try:
				guess = [np.max(final[0]), ss, mm]
				tempbins = list(range(int(dr/2), radius+int(dr/2), dr))
				tempdensity = final[0]
				holdlist = []
				addZero = False
				for val in list(reversed(tempdensity)):
					if not addZero:
						if val >= 0.0:
							holdlist.append(val)
						else:
							addZero = True
							holdlist.append(0.0)
					else:
						holdlist.append(0.0)
				tempdensity = list(reversed(holdlist))
				while len(tempdensity) > len(tempbins):
					tempdensity.pop()
				while len(tempbins) > len(tempdensity):
					tempbins.pop()
				revtempbins = list(np.negative(list(reversed(tempbins))))
				revtempdensity = list(reversed(tempdensity))
				bins = revtempbins + tempbins
				density = revtempdensity + tempdensity
				params, var = curve_fit(gauss_fn, bins, density, p0 = guess)
				params_gauss = np.abs(params)
				## computes 1 SD errors
				var_gauss = np.sqrt(np.diag(var))
				
				def frange(beg, end, step):
					f_range = []
					while beg < end - (step/2.0):
						f_range.append(beg)
						beg += step
					return f_range
				
				guess = [-mm, ss, np.max(final[0]), mm, ss, np.max(final[0])]
				
				tempbins = frange(dr/2.0, radius, dr)
					
				tempdensity = final[0]
				holdlist = []
				addZero = False
				for val in list(reversed(tempdensity)):
					if not addZero:
						if val >= 0.0:
							holdlist.append(val)
						else:
							addZero = True
							holdlist.append(0.0)
					else:
						holdlist.append(0.0)
				tempdensity = list(reversed(holdlist))
				while len(tempdensity) > len(tempbins):
					tempdensity.pop()
				while len(tempbins) > len(tempdensity):
					tempbins.pop()
				revtempbins = list(np.negative(list(reversed(tempbins))))
				revtempdensity = list(reversed(tempdensity))
				bins = revtempbins + tempbins
				density = revtempdensity + tempdensity
				params, var = curve_fit(bimodal, bins, density, p0 = guess)
				params = np.abs(params)
				## computes 1 SD errors
				var = np.sqrt(np.diag(var))
				## average paramters
				stdev = np.average((params[1], params[4]))
				mean = np.average((params[0], params[3]))
				height = np.average((params[2], params[5]))
				stdev_e = np.average((var[1], var[4]))
				mean_e = np.average((var[0], var[3]))
				height_e = np.average((var[2], var[5]))
				params_bimodal = [height, stdev, mean]
				var_bimodal = [height_e, stdev_e, mean_e]
				
				## uncomment following for comparing central vs. peripheral peak fitting errors
#                bimodalmean.append(params_gauss[0])
				bimodalmean.append(mean)
#                bimodalmean.append(tempdensity)
				bimodalsd.append(stdev)
				bimodalheight.append(height)
				auc = 0.0
				step = mean - 5.0*stdev
				while step < mean + 5.0*stdev:
					auc+=0.01*gauss_fn(step,height,stdev,mean)
					step += 0.01
				bimodalauc.append(auc)
#                bimodallist.append(var_bimodal[1])
				gausslist.append(var_gauss[1])
		#        if np.sum(var_bimodal) < np.sum(var_gauss):
				params = params_bimodal
				var = var_bimodal
		#        else:
#                params = params_gauss
#                var = var_gauss
			except RuntimeError:
				params = []
				var = []
			hist_mids = []
			for i in range(len(b)-1):
				hist_mids.append(np.average((b[i],b[i+1])))
			norm_values = []
			for i in final[1]:
				norm_values.append(i/np.max(final[1]))
			return params, var, norm_values, hist_mids, dr
			
		else:
			return dr

	pt_min = PTNUM
	pt_max = PTNUM
	rt_min = RADIUS
	rt_max = RADIUS
	prec_min = PREC
	prec_max = PREC
	iterations = ITER
	PREC = float(PREC)
	one_diff = []
	perc_err = PERCERROR*0.01
	
	
	def roundup(x):
		val = int(math.ceil(x / 10.0)) * 10
		if val >= 30:
			return val
		else:
			return 30

	ptlist = range(pt_min, pt_max+100, 100)
	for pt in ptlist:
		for rt in range(rt_min, rt_max+1, 1):
			for prec in range(prec_min, prec_max+1, 1):
				prec = prec+0.000001
				xrng = roundup(rt + prec*5.0)
				# DR = simulation(pt+1, xrng, BINSIZE, float(prec), float(rt))
				## uncomment below to manually set bin size
				DR = PTNUM
				# print ('ideal bin size: '+ str(DR))
				p_ab, v, d, h_m, DR = simulation(1000000, xrng, DR, float(prec), float(rt), A, ROTATION)
				p_normal, v, d, h_m, DR = simulation(1000000, xrng, DR, float(prec), float(rt), 1.0, 0)
				# print (p)
				a, s, m_ab = p_ab
				a, s, m_normal = p_normal

	return '%.3f'%(m_ab-m_normal)