# -*- coding: utf-8 -*-
"""
Created on Wed Feb 15 12:40:19 2017

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

import csv
import os
import tkinter.filedialog as tkfd
import sys

def filetype(file_path):
	
	split_ext = file_path.split(os.extsep).pop()
	
	return split_ext

def selectfile(root):
	
	file_opt = options = {}
	options['defaultextension'] = '.csv'
	options['filetypes'] = [('CSV (Comma Separated Value)', '.csv'), ('all files', '.*')]
	options['initialdir'] = 'C:\\'
	options['parent'] = root
	options['title'] = 'Open file...'
	file = tkfd.askopenfilename(**file_opt)
	
	return file

def listfiles(folder, file_ext):
	"""
	Returns list of absolute filepaths from folder of files with 
	selected file extension.
	
	file_ext should be in the format 'csv', not '.csv'
	"""

	folder_list = os.listdir(folder)

	file_list = []
	
	for file_path in folder_list:
		split_ext = file_path.split(os.extsep).pop()
		if file_ext == split_ext:
			file_list.append(os.path.normpath(os.path.join(folder, file_path)))
	
	return file_list
	
def extractdata(file_path):
	"""
	Returns data from a file as a list of lines. Works with csv, xls, and xlxs.
	"""    
	
	data = []
	file_ext = file_path.split(os.extsep).pop()
	
	if file_ext == 'csv':
		with open(file_path, 'rb') as f:
			reader = csv.reader(f)
			for row in reader:
				data.append(row)
	elif file_ext == 'xls':
		print ('file type not supported')
	elif file_ext == 'xlsx':
		print ('file type not supported')
	else:
		print ('file type not supported')
	
	return data

def selectfolder(root):
	"""
	Returns the filepath of selected folder. Requires root Tk frame as argument.
	"""

	folder = tkfd.askdirectory(initialdir = 'C:\\', parent = root, title = 'Select a folder', mustexist = True)
	
	return folder

def savelist(root, lst, file_ext):
	"""
	Takes a list and exports it as the desired filetype in the desired location.
	Requires root Tk frame as argument as well as list to be saved.
	file_ext should be in the format 'csv', not '.csv'
	"""
	
	file_opt = options = {}
	options['defaultextension'] = '.csv'
	options['filetypes'] = [('CSV (Comma Separated Value)', '.csv'), ('all files', '.*')]
	options['initialdir'] = 'C:\\'
	options['initialfile'] = 'New Document.csv'
	options['parent'] = root
	options['title'] = 'Save file as...'
	filename = tkfd.asksaveasfilename(**file_opt)
	
	if file_ext == 'csv':
		with open(filename, 'wb') as csvfile:
			writer = csv.writer(csvfile)
			for i in lst:
				writer.writerow(i)
	elif file_ext == 'xls':
		print ('file type not supported')
	elif file_ext == 'xlxs':
		print ('file type not supported')
		
def resourcepath(relative_path):
	""" Get absolute path to resource, works for dev and for PyInstaller """
	try:
		# PyInstaller creates a temp folder and stores path in _MEIPASS
		base_path = sys._MEIPASS
	except Exception:
		base_path = os.path.abspath(".")

	return os.path.join(base_path, relative_path)