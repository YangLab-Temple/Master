# -*- coding: utf-8 -*-
"""
Created on Mon Apr 17 11:06:15 2017

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


import tkinter as tk
import easymode as ez
import csv
import random
import os
import tkinter.filedialog as tkfd
import sys
import easymode as ez
import numpy as np
from numpy.random import choice

class csvmerger(tk.Frame):
	
	def __init__(self, master = None):
		tk.Frame.__init__(self, master)
		
		self.message = tk.StringVar()

		self.binsize = 10
		
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
		
		self.PERCERR = tk.Button(self.LEFTTOPCONTAINER)
		self.PERCERR['text'] = 'Run deconvolution'
		self.PERCERR['command'] = lambda: self.can_percerror()
		self.PERCERR.pack({'side': 'top'})        

		self.QUIT = tk.Button(self.LEFTTOPCONTAINER)
		self.QUIT['text'] = 'Quit'
		self.QUIT['fg'] = red
		self.QUIT['command'] = lambda: self.quit()
		self.QUIT.pack({'side': 'top'})     
		
		self.BINSIZECONTAINER = tk.Frame(self.TOPCONTAINER)
		self.BINSIZECONTAINER.pack({'side': 'top'})        
		
		self.BINSIZELABEL = tk.Label(self.BINSIZECONTAINER)
		self.BINSIZELABEL['text'] = 'binsize (nm): '
		self.BINSIZELABEL.pack({'side': 'left'})
		
		self.BINSIZE = tk.Entry(self.BINSIZECONTAINER)
		self.BINSIZE['width'] = 3
		self.BINSIZE.insert(0, '10')
		self.BINSIZE.pack({'side': 'left'})
		
		self.MESSAGE = tk.Label(self.BOTTOMCONTAINER)
		self.MESSAGE['textvariable'] = self.message
		self.message.set('Set parameters and run simulation.')
		self.MESSAGE.pack({'side': 'right'})
		
	def gen_matrix(self, r, d_r):
		
		##'be' is short for bin edges
		print(r)
		print(d_r)
		if r%d_r > 0:
			be = range(0, int(r+r%d_r), int(d_r))
		else:
			be = range(0, int(r+d_r), int(d_r))
	
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
	
	def smoothdata(self, data, r, d_r):
		"""smoothds data with 3 moving window and takes abs value average"""
		smooth_data = []
		r += 1
		
		for i in range(len(data)):
			if i == 0 or i == (len(data) - 1):
				smooth_data.append(data[i])
			else:
				smooth_data.append(np.average([data[i-1], data[i], data[i+1]]))
		
		##comment out for smoothing
		# smooth_data = []
		# for i in range(len(data)):
			# smooth_data.append(data[i])
	
		##adds + and - bins
		final_smooth_data = []
		for i in range(int(r/d_r)):
			final_smooth_data.append(smooth_data[i] + smooth_data[len(smooth_data)-1-i])
	
		return list(reversed(final_smooth_data))
		
	def deconvolution(self, hv, be, r, d_r):
		"""hv = hist_values, be = bin_edges"""
		
		density = []
		matrix = self.gen_matrix(r, d_r)

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
		
	def make_hist(self, data, r, d_r):
		
		hist_values, bin_edges = np.histogram(data, bins = 2 * int(r/d_r), range = (-r, r))
		
		new_bin_edges = []    
		for i in bin_edges:
			if i >= 0:
				new_bin_edges.append(i)
				
		new_hist_values = self.smoothdata(hist_values, r, d_r)
		  
		return new_hist_values, new_bin_edges
		
	def csv_read(self, path):
		with open(path, 'r') as csvfile:
			reader = csv.reader(csvfile, delimiter = ',')
			holdlist = []
			for row in reader:
				holdlist.append(float(row[1]))
		return holdlist
			
	def can_percerror(self):
		
		if True:
			self.binsize = int(self.BINSIZE.get())
			self.select_folder_percerr()
		else:
			self.message.set('Insert a value for every field.')
	
	def select_folder_percerr(self):
		xy_data_file = ez.selectfile(self)
		xy_data = self.csv_read(xy_data_file)
		dr = int(self.BINSIZE.get())
		radius = np.round(np.max(np.abs(xy_data))+2*dr)
		a,b = self.make_hist(xy_data, radius, dr)
		final = self.deconvolution(a,b,radius,dr)
		tempbins = list(range(int(dr/2), int(radius+dr/2), dr))
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
		smooth_data = []
		for i in range(len(density)):
			if i == 0 or i == (len(density) - 1):
				smooth_data.append(density[i])
			else:
				smooth_data.append(np.average([density[i-1], density[i], density[i+1]]))
		smooth_data = smooth_data / np.max(smooth_data)
		tocsv = []
		for i in range(len(bins)):
			tocsv.append((bins[i],smooth_data[i]))
		ez.savelist(self, tocsv, 'csv')   
		self.message.set(str('Done'))


root = tk.Tk()
root.title(string = 'YangLab Deconvolution')

orange = '#FF7F22'
black = '#000000'
gray = '#D3D3D3'
white = '#FFFFFF'
red = '#FF0000'
green = '#008000'

app = csvmerger(master = root)
app.mainloop()
root.destroy()