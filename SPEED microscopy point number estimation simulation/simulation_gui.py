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


import Tkinter as tk
import easymode as ez
import csv
import random
import os
import tkFileDialog as tkfd
import sys
import easymode as ez
import simulation as sim
import numpy as np
from numpy.random import choice

class csvmerger(tk.Frame):
    
    def __init__(self, master = None):
        tk.Frame.__init__(self, master)
        
        self.message = tk.StringVar()
        
        self.radius = None
        self.numpoints = None
        self.prec = None
        self.canSim = False
        
        self.binsize = 10
        self.percerror = 10
        self.iter = None
        
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
        self.PERCERR['text'] = 'Simulate % error'
        self.PERCERR['command'] = lambda: self.can_percerror()
        self.PERCERR.pack({'side': 'top'})        
        
        self.SELECT = tk.Button(self.LEFTTOPCONTAINER)
        self.SELECT['text'] = 'Simulate data set'
        self.SELECT['command'] = lambda: self.can_simulate()
        self.SELECT.pack({'side': 'top'})
       
        self.QUIT = tk.Button(self.LEFTTOPCONTAINER)
        self.QUIT['text'] = 'Quit'
        self.QUIT['fg'] = red
        self.QUIT['command'] = lambda: self.quit()
        self.QUIT.pack({'side': 'top'})     
        
        self.NUMPOINTSCONTAINER = tk.Frame(self.TOPCONTAINER)
        self.NUMPOINTSCONTAINER.pack({'side': 'top'})        
        
        self.NUMPOINTSLABEL = tk.Label(self.NUMPOINTSCONTAINER)
        self.NUMPOINTSLABEL['text'] = '# points: '
        self.NUMPOINTSLABEL.pack({'side': 'left'})
        
        self.NUMPOINTS = tk.Entry(self.NUMPOINTSCONTAINER)
        self.NUMPOINTS['width'] = 6
        self.NUMPOINTS.insert(0, '100')
        self.NUMPOINTS.pack({'side': 'left'})

        self.RADIUSCONTAINER = tk.Frame(self.TOPCONTAINER)
        self.RADIUSCONTAINER.pack({'side': 'top'})        
        
        self.RADIUSLABEL = tk.Label(self.RADIUSCONTAINER)
        self.RADIUSLABEL['text'] = 'radius (nm): '
        self.RADIUSLABEL.pack({'side': 'left'})
        
        self.RADIUS = tk.Entry(self.RADIUSCONTAINER)
        self.RADIUS['width'] = 3
        self.RADIUS.insert(0, '0')
        self.RADIUS.pack({'side': 'left'})
        
        self.PRECCONTAINER = tk.Frame(self.TOPCONTAINER)
        self.PRECCONTAINER.pack({'side': 'top'})
        
        self.PRECLABEL = tk.Label(self.PRECCONTAINER)
        self.PRECLABEL['text'] = 'precision (nm): '
        self.PRECLABEL.pack({'side': 'left'})
        
        self.PREC = tk.Entry(self.PRECCONTAINER)
        self.PREC['width'] = 3
        self.PREC.insert(0, '10')
        self.PREC.pack({'side': 'left'})
        
        self.ITERCONTAINER = tk.Frame(self.TOPCONTAINER)
        self.ITERCONTAINER.pack({'side': 'top'})
        
        self.ITERLABEL = tk.Label(self.ITERCONTAINER)
        self.ITERLABEL['text'] = '# iterations: '
        self.ITERLABEL.pack({'side': 'left'})
        
        self.ITER = tk.Entry(self.ITERCONTAINER)
        self.ITER['width'] = 6
        self.ITER.insert(0, '10000')
        self.ITER.pack({'side': 'left'})
        
        self.MESSAGE = tk.Label(self.BOTTOMCONTAINER)
        self.MESSAGE['textvariable'] = self.message
        self.message.set('Set parameters and save simulation.')
        self.MESSAGE.pack({'side': 'right'})

    def can_simulate(self):
        
        try:
            self.radius = int(self.RADIUS.get())
            self.numpoints = int(self.NUMPOINTS.get())
            self.prec = int(self.PREC.get())
            self.message.set('Select a save folder.')
            self.select_folder_sim()
        except:
            self.message.set('Each parameter needs a value.')
    
    def can_percerror(self):
        
        if True:
            self.radius = int(self.RADIUS.get())
            self.numpoints = int(self.NUMPOINTS.get())
            self.prec = int(self.PREC.get())
            self.iter = int(self.ITER.get())
            self.select_folder_percerr()
        else:
            self.message.set('Bin size must be even number.')
            
    def select_folder_sim(self):
        
        self.folder = ez.selectfolder(self)
        self.message.set('Simulating and saving...')
        
        def area_fn(X):
            X = float(X)
            A = -(10.0**2)*np.pi
            B = 10.0*2*np.pi
            return X*B+A
        
        def gauss_fn(x, s, m):
            a = area_fn(m)
            x = float(x)
            s = float(s)
            m = float(m)
            return a*np.e**(-(x-m)**2.0/(2.0*s**2.0))
            
        def combine(x):
            s = int(self.PREC.get())
            m = int(self.RADIUS.get())
            return (area_fn(x) * gauss_fn(x, s, m))
                
        ##starting with perfect x,y and adding error
        xydata = []
        while len(xydata) < int(self.NUMPOINTS.get()):
            theta = np.random.random()*360.0
            ## modify here for precision distribution sampling
#            ss = choice([3,5,7,9], p=[0.1475,0.2775,0.3075,0.2675])
#            ss = choice([4.5,5.5,6.5,7.5,8.5,9.5], p=[0.02,0.05,0.07,0.11,0.2,0.55])
            ss = int(self.PREC.get())
            y_prec = np.random.normal(0.0, ss)
            z_prec = np.random.normal(0.0, ss)
            xydata.append((int(self.RADIUS.get())*np.cos(theta)+y_prec, int(self.RADIUS.get())*np.sin(theta)+z_prec))
        
        with open(self.folder + '/' + 'simulation.csv', 'wb') as z:
            writer = csv.writer(z)
            writer.writerow(('y (nm)', 'z (nm)'))
            for i in range(len(xydata)):
                writer.writerow(xydata[i])
        
        self.message.set('Success. Quit or simulate again.')
    
    def select_folder_percerr(self):
        
        results = sim.sim(self, self.numpoints, self.radius, self.prec, self.iter, self.binsize, self.percerror)
        for i in results:
            print i
        one, two, three, four, five, six, seven, eight, nine, ten  = results        
        self.message.set(str(one) + '% confidence +/- 0.1*precsion\n' + str(two) + '% confidence +/- 0.2*precsion\n' + str(three) + '% confidence +/- 0.3*precsion\n' + str(four) + '% confidence +/- 0.4*precsion\n' + str(five) + '% confidence +/- 0.5*precsion\n' + str(six) + '% confidence +/- 0.6*precsion\n' + str(seven) + '% confidence +/- 0.7*precsion\n' + str(eight) + '% confidence +/- 0.8*precsion\n' + str(nine) + '% confidence +/- 0.9*precsion\n' + str(ten) + '% confidence +/- 1.0*precsion\n')


root = tk.Tk()
root.title(string = 'YangLab Simulation')

orange = '#FF7F22'
black = '#000000'
gray = '#D3D3D3'
white = '#FFFFFF'
red = '#FF0000'
green = '#008000'

app = csvmerger(master = root)
app.mainloop()
root.destroy()