# -*- coding: utf-8 -*-
"""
Created on Mon Apr 17 11:06:15 2017

@author: tuf61935
"""


import tkinter as tk
import csv
import random
import os
import tkinter.filedialog as tkfd
import sys
import numpy as np
import math

class categorization(tk.Frame):
	
	def __init__(self, master = None):
		tk.Frame.__init__(self, master)
		self.message = tk.StringVar()
		self.pack()
		self.createwidgets()
		self.master = master
		
	def createwidgets(self):
		self.TOPCONTAINER = tk.Frame(self)
		self.TOPCONTAINER.pack({'side': 'top'})           
		
		self.LEFTTOPCONTAINER = tk.Frame(self.TOPCONTAINER)
		self.LEFTTOPCONTAINER.pack({'side': 'left'})   
		
		self.BOTTOMCONTAINER = tk.Frame(self)
		self.BOTTOMCONTAINER.pack({'side': 'bottom'})
		
		self.PERCERR = tk.Button(self.LEFTTOPCONTAINER)
		self.PERCERR['text'] = 'Categorize route'
		self.PERCERR['command'] = lambda: self.select_folder()
		self.PERCERR.pack({'side': 'top'})        

		self.QUIT = tk.Button(self.LEFTTOPCONTAINER)
		self.QUIT['text'] = 'Quit'
		self.QUIT['fg'] = red
		self.QUIT['command'] = lambda: self.quit()
		self.QUIT.pack({'side': 'top'})     

		self.BINSIZECONTAINER = tk.Frame(self.TOPCONTAINER)
		self.BINSIZECONTAINER.pack({'side': 'top'})        
		
		self.BINSIZELABEL = tk.Label(self.BINSIZECONTAINER)
		self.BINSIZELABEL['text'] = 'ideal bin size (nm): '
		self.BINSIZELABEL.pack({'side': 'left'})
		
		self.BIN_SIZE = tk.Entry(self.BINSIZECONTAINER)
		self.BIN_SIZE['width'] = 3
		self.BIN_SIZE.insert(0, '10')
		self.BIN_SIZE.pack({'side': 'left'})
		
		self.PTNUMCONTAINER = tk.Frame(self.TOPCONTAINER)
		self.PTNUMCONTAINER.pack({'side': 'top'})        
		
		self.PTNUMLABEL = tk.Label(self.PTNUMCONTAINER)
		self.PTNUMLABEL['text'] = 'ground truth point #: '
		self.PTNUMLABEL.pack({'side': 'left'})
		
		self.PTNUM = tk.Entry(self.PTNUMCONTAINER)
		self.PTNUM['width'] = 8
		self.PTNUM.insert(0, '10000')
		self.PTNUM.pack({'side': 'left'})
		
		self.MESSAGE = tk.Label(self.BOTTOMCONTAINER)
		self.MESSAGE['textvariable'] = self.message
		self.message.set('Click "categorize route" to begin')
		self.MESSAGE.pack({'side': 'right'})

	def categorization_algorithm(self, DATA, BINSIZE, PTNUM):

		def roundup(x):
			val = int(math.ceil(x / BINSIZE)) * BINSIZE
			if val >= 50:
				return val
			else:
				return 50
		global radius
		global peripheral_radius
		global precision
		global num_points
		## dist types: 'peripheral', 'uniform', 'central', 'bimodal'
		global dist_type
		global yrng
		yrng = roundup(np.max([np.abs(float(y)) for x,y in DATA]))
		global smooth
		## set to True for +/- 1 bin smoothing window
		smooth = False
		## set to True to remove high error values from 3D
		global allZero
		allZero = True
		## each bimodal distribution has equal density in central and outer peak
		global bimodal_factor
		npc_range = yrng
		sigma_range = yrng
		dr = BINSIZE

		def smoothdata(data, r, d_r):
			"""smoothds data with 3 moving window and takes abs value average, related to smooth variable"""
			smooth_data = []
			r += 1

			smooth_data = []
			for i in range(len(data)):
				smooth_data.append(data[i])

			## adds + and - bins
			final_smooth_data = []
			for i in range(int(r/d_r)):
				final_smooth_data.append(smooth_data[i] + smooth_data[len(smooth_data)-1-i])
			return list(reversed(final_smooth_data))
					
		def make_hist(data, r, d_r):
			r = int(r)
			hist_values, bin_edges = np.histogram(data, bins = int(2 * r/d_r), range = (-r, r))
			new_bin_edges = []    
			for i in bin_edges:
				if i >= 0:
					new_bin_edges.append(i)
			new_hist_values = smoothdata(hist_values, r, d_r)
			return new_hist_values, new_bin_edges
			
		def gen_matrix(r, d_r):
			
			r = int(r)
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

			## generate areas of sections closest to x axis
			for i in range(len(matrix)):
				theta = np.arccos(float(be[len(be)-2-i])/float(be[len(be)-1-i]))
				arc_area = (theta/(2*np.pi)) * np.pi * float(be[len(be)-1-i])**2
				tri_area = 0.5 * float(be[len(be)-2-i]) * (np.sin(theta) * float(be[len(be)-1-i]))
				matrix[i].append(4 * (arc_area - tri_area))
			
			## skipping factor
			x = 2
			## generate areas of layers going further out from x axis
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
				## calculate how much to subtract from bin
				density_sub = 0
				y = 0
				for j in range(x):
					density_sub += density[y] * matrix[j][i]
					y += 1
				
				## calculate final bin value
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

		## simulate points from all four ground truth models
		def gauss_dist():
			if dist_type != 'bimodal':
				yzdata = []
				while len(yzdata) < num_points:
					theta = np.random.random()*360.0
					y_prec = np.random.normal(0.0, precision)
					z_prec = np.random.normal(0.0, precision)
					yzdata.append((radius*np.cos(theta)+y_prec, radius*np.sin(theta)+z_prec))
			else:
				yzdata = []
				while len(yzdata) < num_points:
					theta = np.random.random()*360.0
					y_prec = np.random.normal(0.0, precision)
					z_prec = np.random.normal(0.0, precision)
					if random.random() < bimodal_factor:
						yzdata.append((0.0*np.cos(theta)+y_prec, 0.0*np.sin(theta)+z_prec))
					else:
						yzdata.append((peripheral_radius*np.cos(theta)+y_prec, peripheral_radius*np.sin(theta)+z_prec))
			return yzdata

		def uniform_dist():
			yzdata = []
			while len(yzdata) < num_points:
				rand_radius = random.random()
				if random.random() < rand_radius:
					rand_radius = rand_radius * peripheral_radius
					theta = np.random.random()*360.0
					y_prec = np.random.normal(0.0, precision)
					z_prec = np.random.normal(0.0, precision)
					yzdata.append((rand_radius*np.cos(theta)+y_prec, rand_radius*np.sin(theta)+z_prec))
			return yzdata

		## experimental simulation
		exp_trials = []
		# for i in range(num_trials):
		# if dist_type == 'uniform':
			# exp = uniform_dist()
		# else:
			# exp = gauss_dist()
		exp = [float(y) for x,y in DATA]
		a,b = make_hist(exp, yrng, dr)
		final_unsmooth_exp, final_smooth_exp, final_2d_exp = deconvolution(a, b, yrng, dr)
		temp = []
		addZero = False
		if smooth:
			final_unsmooth_exp = [i for i in final_smooth_exp]
		for i in list(reversed(final_unsmooth_exp)):
			if not addZero:
				if i >= 0.0:
					temp.append(i)
				else:
					# if allZero:
						# addZero = True
					temp.append(0.0)
			else:
				temp.append(0.0)
		final_unsmooth_exp = [i/np.max(temp) for i in list(reversed(temp))]
		for asdf in final_unsmooth_exp:
			print (asdf)
		## outlier detection
		bin_edges = []
		for i in range(len(b)-1):
			bin_edges.append((b[i]+b[i+1])/2)
		_mean = np.average(bin_edges, weights=final_unsmooth_exp)
		_stdev =  math.sqrt(np.average((bin_edges-_mean)**2, weights=final_unsmooth_exp))
		for i in range(len(bin_edges)):
			if bin_edges[i] < _mean-2*_stdev or bin_edges[i] > _mean+2*_stdev:
				final_unsmooth_exp[i] = 0.0
		print (bin_edges)
		print (final_unsmooth_exp)
		exp_trials.append(final_unsmooth_exp)

		all_ideals = {'peripheral':{}, 'central':{}, 'bimodal':{}, 'uniform':{}}
		for all in ['peripheral','bimodal','uniform']:
			for r in range(1,npc_range+1):
				all_ideals[all][str(r)] = {}
		all_ideals['central']['0'] = {}

		## calculate ground truth models for SAR
		## peripheral
		num_points = PTNUM
		# hold_radius = radius
		# radius = peripheral_radius
		# hold_type = dist_type
		dist_type = 'peripheral'
		for r in range(1,npc_range+1):
			for p in range(1, sigma_range+1):
				print (r,p, dist_type)
				radius = r
				precision = p
				ideal = gauss_dist()
				ideal = [float(y) for y,z in ideal]
				a,b = make_hist(ideal, yrng, dr)
				final_unsmooth_ideal, final_smooth_ideal, final_2d_ideal = deconvolution(a, b, yrng, dr)
				temp = []
				addZero = False
				if smooth:
					final_unsmooth_ideal = [i for i in final_smooth_ideal]
				for i in list(reversed(final_unsmooth_ideal)):
					if not addZero:
						if i >= 0.0:
							temp.append(i)
						else:
							temp.append(0.0)
					else:
						temp.append(0.0)
				final_unsmooth_ideal_peripheral = [i/np.max(temp) for i in list(reversed(temp))]
				all_ideals['peripheral'][str(r)][str(p)] = [i/np.max(temp) for i in list(reversed(temp))]

		# ## uniform
		for r in range(1,npc_range+1):
			for p in range(1, sigma_range+1):
				print (r, p, 'uniform')
				peripheral_radius = r
				precision = p
				ideal = uniform_dist()
				ideal = [float(y) for y,z in ideal]
				a,b = make_hist(ideal, yrng, dr)
				final_unsmooth_ideal, final_smooth_ideal, final_2d_ideal = deconvolution(a, b, yrng, dr)
				temp = []
				addZero = False
				if smooth:
					final_unsmooth_ideal = [i for i in final_smooth_ideal]
				for i in list(reversed(final_unsmooth_ideal)):
					if not addZero:
						if i >= 0.0:
							temp.append(i)
						else:
							temp.append(0.0)
					else:
						temp.append(0.0)
				final_unsmooth_ideal_uniform = [i/np.max(temp) for i in list(reversed(temp))]
				all_ideals['uniform'][str(r)][str(p)] = [i/np.max(temp) for i in list(reversed(temp))]

		## bimodal
		dist_type = 'bimodal'
		for r in range(1,npc_range+1):
			for p in range(1, sigma_range+1):
				print (r,p, dist_type)
				peripheral_radius = r
				precision = p
				##calculate bimodal_factor
				inner_area = np.pi*dr**2.0
				outer_area = np.pi*(((peripheral_radius+dr/2.0)**2.0)-((peripheral_radius-dr/2.0)**2.0))
				bimodal_factor = inner_area / (inner_area + outer_area)
				ideal = gauss_dist()
				ideal = [float(y) for y,z in ideal]
				a,b = make_hist(ideal, yrng, dr)
				final_unsmooth_ideal, final_smooth_ideal, final_2d_ideal = deconvolution(a, b, yrng, dr)
				temp = []
				addZero = False
				if smooth:
					final_unsmooth_ideal = [i for i in final_smooth_ideal]
				for i in list(reversed(final_unsmooth_ideal)):
					if not addZero:
						if i >= 0.0:
							temp.append(i)
						else:
							temp.append(0.0)
					else:
						temp.append(0.0)
				final_unsmooth_ideal_bimodal = [i/np.max(temp) for i in list(reversed(temp))]
				all_ideals['bimodal'][str(r)][str(p)] = [i/np.max(temp) for i in list(reversed(temp))]

		## central
		radius = 0.0000001
		dist_type = 'central'
		print (radius, dist_type)
		for p in range(1, sigma_range+1):
			print (p, dist_type)
			precision = p
			ideal = gauss_dist()
			ideal = [float(y) for y,z in ideal]
			a,b = make_hist(ideal, yrng, dr)
			final_unsmooth_ideal, final_smooth_ideal, final_2d_ideal = deconvolution(a, b, yrng, dr)
			temp = []
			addZero = False
			if smooth:
				final_unsmooth_ideal = [i for i in final_smooth_ideal]
			for i in list(reversed(final_unsmooth_ideal)):
				if not addZero:
					if i >= 0.0:
						temp.append(i)
					else:
						temp.append(0.0)
				else:
					temp.append(0.0)
			final_unsmooth_ideal_central = [i/np.max(temp) for i in list(reversed(temp))]
			all_ideals['central']['0'][str(p)] = [i/np.max(temp) for i in list(reversed(temp))]

		num_correct = {'peripheral':{}, 'central':{}, 'bimodal':{}, 'uniform':{}}
		min_sar = []
		lowest_sar = [999999.0,'','','']
		for trial in exp_trials:
			for asdf in trial:
				print (asdf)
			for i in all_ideals:
				for j in all_ideals[i]:
					for k in all_ideals[i][j]:
						sar = 0.0
						for l in range(len(trial)):
							exp_val = trial[l]
							sar += np.abs(exp_val - all_ideals[i][j][k][l])
						if not math.isnan(sar):
							min_sar.append([sar,i,j,k])

		min_sar.sort()
		with open('categorization_results.csv', 'w', newline='') as csvfile:
			spamwriter = csv.writer(csvfile, delimiter=',')
			spamwriter.writerow(['sar','type','mean','sigma'])
			for row in min_sar:
				spamwriter.writerow(row)
		self.message.set('Top 3 results:\n'+'type: '+str(min_sar[0][1])+', mean: '+str(min_sar[0][2])+', sigma: '+str(min_sar[0][3])+', sar: '+str(min_sar[0][0])+'\n'+'type: '+str(min_sar[1][1])+', mean: '+str(min_sar[1][2])+', sigma: '+str(min_sar[1][3])+', sar: '+str(min_sar[1][0])+'\n'+'type: '+str(min_sar[2][1])+', mean: '+str(min_sar[2][2])+', sigma: '+str(min_sar[2][3])+', sar: '+str(min_sar[2][0]))
		
	def selectfile(self):
		
		file_opt = options = {}
		options['defaultextension'] = '.csv'
		options['filetypes'] = [('CSV (Comma Separated Value)', '.csv'), ('all files', '.*')]
		options['initialdir'] = 'C:\\'
		options['parent'] = self
		options['title'] = 'Open file...'
		file = tkfd.askopenfilename(**file_opt)
		
		return file
		
	def selectfolder(self):
		"""
		Returns the filepath of selected folder. Requires self Tk frame as argument.
		"""

		folder = tkfd.askdirectory(initialdir = 'C:\\', parent = self, title = 'Select a folder', mustexist = True)
		
		return folder

	def extractdata(self, file_path):
		"""
		Returns data from a file as a list of lines. Works with csv, xls, and xlxs.
		"""    
		
		data = []
		file_ext = file_path.split(os.extsep).pop()
		
		if file_ext == 'csv':
			with open(file_path, 'r') as f:
				reader = csv.reader(f)
				for row in reader:
					x,y = row
					data.append(row)
		elif file_ext == 'xls':
			print ('must be csv')
		elif file_ext == 'xlsx':
			print ('must be csv')
		else:
			print ('must be csv')
		
		return data
	
	def select_folder(self):
		file = self.selectfile()
		data = self.extractdata(file)
		binsize = int(self.BIN_SIZE.get())
		ptnum = int(self.PTNUM.get())
		results = self.categorization_algorithm(data, binsize, ptnum)


root = tk.Tk()
root.title(string = 'Categorization algorithm')
root.geometry('400x150')

red = '#FF0000'

app = categorization(master = root)
app.mainloop()
root.destroy()