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
from scipy import stats

class bin_size_optimization(tk.Frame):
	
	def __init__(self, master = None):
		tk.Frame.__init__(self, master)
		self.message = tk.StringVar()
		self.folder = None
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
		
		self.SELECT = tk.Button(self.LEFTTOPCONTAINER)
		self.SELECT['text'] = 'Optimize bin size'
		self.SELECT['command'] = lambda: self.select_folder()
		self.SELECT.pack({'side': 'top'})        
	   
		self.QUIT = tk.Button(self.LEFTTOPCONTAINER)
		self.QUIT['text'] = 'Quit'
		self.QUIT['fg'] = red
		self.QUIT['command'] = lambda: self.quit()
		self.QUIT.pack({'side': 'top'})          

		self.MESSAGE = tk.Label(self.BOTTOMCONTAINER)
		self.MESSAGE['textvariable'] = self.message
		self.message.set('Click "optimize bin size" to begin')
		self.MESSAGE.pack({'side': 'right'})

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
		
	def selectfolder(self):
		"""
		Returns the filepath of selected folder. Requires self Tk frame as argument.
		"""

		folder = tkfd.askdirectory(initialdir = 'C:\\', parent = self, title = 'Select a folder', mustexist = True)
		
		return folder
		
	def selectfile(self):
	
		file_opt = options = {}
		options['defaultextension'] = '.csv'
		options['filetypes'] = [('CSV (Comma Separated Value)', '.csv'), ('all files', '.*')]
		options['initialdir'] = 'C:\\'
		options['parent'] = self
		options['title'] = 'Open file...'
		file = tkfd.askopenfilename(**file_opt)
		
		return file
		
	def sim(self, DATA):

		def simulation():
		
			def gen_matrix(r, d_r):
				
				##'be' is short for bin edges
				if r%d_r > 0:
					be = list(range(0, r+r%d_r, d_r))
				else:
					be = list(range(0, r+d_r, d_r))
			
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

			jkl = []
			for x,y in DATA:
				jkl.append(float(y))
			maxbinsize = int(np.max(jkl)+1)
			minbinsize = 1
			binsizes = []
			binsizesdata = [[] for variable in range(1, maxbinsize+1)]
			self.message.set('calculating ideal bin size...')
			for binoptimization in range(1):
				for binsize in range(1, maxbinsize+1):
					if binsize >= minbinsize:
						error = 0
						a,b = make_hist(jkl, maxbinsize, binsize)
						final_unsmooth, final_smooth, final_2d = deconvolution(a, b, maxbinsize, binsize)
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
						## outlier detection
						# print('asdfasdfasdf')
						# bin_edges = []
						# for i in range(len(b)-1):
							# bin_edges.append((b[i]+b[i+1])/2)
						# print (binsize)
						# print(len(bin_edges))
						# print(len(final_unsmooth))
						# while len(bin_edges) != len(final_unsmooth):
							# if len(bin_edges) > len(final_unsmooth):
								# bin_edges.pop()
							# else:
								# bin_edges.append(bin_edges[-1]+binsize)
						# _mean = np.average(bin_edges, weights=final_unsmooth)
						# _stdev =  math.sqrt(np.average((bin_edges-_mean)**2, weights=final_unsmooth))
						# for i in range(len(bin_edges)):
							# if bin_edges[i] < _mean-3*_stdev or bin_edges[i] > _mean+3*_stdev:
								# try:
									# final_unsmooth[i] = 0.0
								# except:
									# pass
						matrix = gen_matrix(maxbinsize, binsize)
						newmatrix = []
						for i in matrix:
							newmatrix.append(list(reversed(i)))
						matrix = list(reversed(newmatrix))
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
									error += np.square(a[ncol] - binsub) / a[ncol]
							except:
								pass
						popped = a.pop()
						while popped == 0:
							popped = a.pop()
						binsizesdata[binsize-1].append((error, len(a)+1,1-stats.chi2.cdf(error, len(a)+1),binsize))
					else:
						binsizesdata[binsize-1].append((1000000.0,1,0.0,binsize))         
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
				print (('binsize: '+str(binsize), 'p-value: '+str(pval)))
				if pval > 0.99:
					dr = binsize
					break
				else:
					dr = maxbinsize
			with open('binsize_optimization_results.csv', 'w', newline='') as csvfile:
				spamwriter = csv.writer(csvfile, delimiter=',')
				spamwriter.writerow(['chisq','df','pvalue','binsize'])
				for row in finalbinsizes:
					spamwriter.writerow(row)
			return dr

		DR = simulation()

		return 'ideal bin size: '+ str(DR)
	
	def select_folder(self):
		file = self.selectfile()
		data = self.extractdata(file)
		results = self.sim(data)     
		self.message.set(results)


root = tk.Tk()
root.title(string = 'Bin Size Optimization')
root.geometry('300x100')

red = '#FF0000'

app = bin_size_optimization(master = root)
app.mainloop()
root.destroy()