# -*- coding: utf-8 -*-
"""
Created on Wed Mar 14 20:37:04 2018

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
import numpy as np
from scipy import stats

def roundup(x):
    val = int(math.ceil(x / 10.0)) * 10
    if val >= 30:
        return val
    else:
        return 30

global radius
## don't edit the 0.0000001, just the integer value
## for a bimodal distribution, put outermost radius
radius = 23 + 0.000001
global peripheral_radius
peripheral_radius = 23 + 0.000001
global precision
precision = 10
global num_points
num_points = 100
global dist_type
## dist types: 'peripheral', 'uniform', 'central', 'bimodal'
dist_type = 'peripheral'
global num_trials
num_trials = 1000
global yrng
yrng = roundup(peripheral_radius + precision*5.0)
global smooth
## set to True for +/- 1 bin smoothing window
smooth = False
## set to True to remove high error values from 3D
global allZero
allZero = True
global binsizemult
## for increasing the range of the bin size optimization algorithm (precision * binsizemult)
binsizemult = 2
global bimodal_factor
## must be modified so that the ground truth peripheral and central peaks in the bimodal distribution in the 3D distribution are equal
bimodal_factor = 0.0

def smoothdata(data, r, d_r):
    """smoothds data with 3 moving window and takes abs value average, related to smooth variable"""
    smooth_data = []
    r += 1

    smooth_data = []
    for i in range(len(data)):
        smooth_data.append(data[i])

    ## adds + and - bins
    final_smooth_data = []
    for i in range(r/d_r):
        final_smooth_data.append(smooth_data[i] + smooth_data[len(smooth_data)-1-i])
    return list(reversed(final_smooth_data))
            
def make_hist(data, r, d_r):
    r = int(r)
    hist_values, bin_edges = np.histogram(data, bins = 2 * r/d_r, range = (-r, r))
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
            if random.random() < 1.0/(bimodal_factor*peripheral_radius):
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
    
def binsize_determination(type_of_dist = None):
    minbinsize = 3
    binsizesdata = [[] for variable in range(1, binsizemult*int(precision)+1)]
    for binoptimization in range(1000):
    
        for binsize in range(1, binsizemult*int(precision)+1):
            if binsize >= minbinsize:
                error = 0
                jkl = []
                if type_of_dist == 'uniform':
                    jkl = uniform_dist()
                else:
                    jkl = gauss_dist()
                a,b = make_hist(jkl, yrng, binsize)
                final_unsmooth, final_smooth, final_2d = deconvolution(a, b, yrng, binsize)
                holdlist = []
                addZero = False
                for val in list(reversed(final_unsmooth)):
                    if not addZero:
                        if val >= 0.0:
                            holdlist.append(val)
                        else:
                            if allZero:
                                addZero = True
                            holdlist.append(0.0)
                    else:
                        holdlist.append(0.0)
                final_unsmooth = list(reversed(holdlist))
                ## rescale ideal data
                matrix = gen_matrix(yrng, binsize)
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

    print 'binsize determination'
    for chisq, df, pval, binsize in finalbinsizes:
        print chisq
    for chisq, df, pval, binsize in finalbinsizes:
        print pval
    for chisq, df, pval, binsize in finalbinsizes:
        print binsize 

    for chisq, df, pval, binsize in finalbinsizes:
        if pval > 0.99:
            dr = binsize
            break
        else:
            dr = int(precision)
    return dr

dr = binsize_determination()
print 'dr: '+str(dr)

## experimental simulation
exp_trials = []
for i in range(num_trials):
    if dist_type == 'uniform':
        exp = uniform_dist()
    else:
        exp = gauss_dist()
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
                if allZero:
                    addZero = True
                temp.append(0.0)
        else:
            temp.append(0.0)
    final_unsmooth_exp = [i/np.max(temp) for i in list(reversed(temp))]
    exp_trials.append(final_unsmooth_exp)

## calculate ground truth models for SAR
## peripheral
num_points = 1000000
hold_radius = radius
radius = peripheral_radius
hold_type = dist_type
dist_type = 'peripheral'
ideal = gauss_dist()
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

## uniform
ideal = uniform_dist()
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

## bimodal
dist_type = 'bimodal'
ideal = gauss_dist()
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

## central
radius = 0.0000001
dist_type = 'central'
ideal = gauss_dist()
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
dist_type = hold_type

## determine SAR and categorize
num_correct = {'peripheral':0, 'central':0, 'bimodal':0, 'uniform':0}
for trial in exp_trials:
    peripheral_sar = 0.0
    central_sar = 0.0
    bimodal_sar = 0.0
    uniform_sar = 0.0
    for i in range(len(trial)):
        exp_val = trial[i]
        peripheral_sar += np.abs(exp_val - final_unsmooth_ideal_peripheral[i])
        central_sar += np.abs(exp_val - final_unsmooth_ideal_central[i])
        bimodal_sar += np.abs(exp_val - final_unsmooth_ideal_bimodal[i])
        uniform_sar += np.abs(exp_val - final_unsmooth_ideal_uniform[i])
    sars = [(peripheral_sar, 'peripheral'),(central_sar, 'central'),(bimodal_sar, 'bimodal'),(uniform_sar, 'uniform')]
    sars.sort()
    num_correct[sars[0][1]] += 1

print 'Categorization success: ' + str(100.0*num_correct[dist_type]/num_trials)