# -*- coding: utf-8 -*-
"""
Created on Tue Dec 20 12:16:16 2016

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
import numpy as np

def pt_on_sphere(radius = 1.0):
    """Generates a single point on the surface of a sphere with a given radius
    according to the Cook (1957) method. Cook, J. M. "Technical Notes and Short
    Papers: Rational Formulae for the Production of a Spherically Symmetric 
    Probability Distribution." Math. Tables Aids Comput. 11, 81-82, 1957."""
    
    ##Samples four numbers from a uniform distribution from (-1, 1) whose
    ##summed squares are greater than or equal to 1.
    nums = [2, 2, 2, 2]
    while nums[0]**2 + nums[1]**2 + nums[2]**2 + nums[3]**2 >= 1:
        nums = []
        for i in range(4):
            if np.random.random() > 0.5:
                nums.append(-1 * np.random.random())
            else:
                nums.append(np.random.random())
    
    ##The rules of quaternion transformation are used to generate (x, y, z)
    ##Cartesian coordinates.
    x = 2 * ((nums[1] * nums[3] + nums[0] * nums[2]) / 
            (nums[0]**2 + nums[1]**2 + nums[2]**2 + nums[3]**2))         
    y = 2 * ((nums[2] * nums[3] - nums[0] * nums[1]) / 
            (nums[0]**2 + nums[1]**2 + nums[2]**2 + nums[3]**2))      
    z = ((nums[0]**2 + nums[3]**2 - nums[1]**2 - nums[2]**2) / 
        (nums[0]**2 + nums[1]**2 + nums[2]**2 + nums[3]**2))
    
    ##Accounts for the desired radius.
    x = radius * x
    y = radius * y
    z = radius * z
    
    return (x, y, z)

def histogram_msd(histogram):
    """Compares the distance histogram of the simulated data with the distance
    histogram of the experimental data."""
    
    ##Normalizes the distance histogram of the simulated data.
    norm_histogram = []
    for val in histogram:    
        norm_histogram.append(float(val) / np.max(histogram))
    
    ##Experimentally determined distance histogram values.
    sstr3_rab8a_exp_values = [0.826, 1.000, 0.804, 0.500, 0.261, 0.130]
    
    ##Makes simulated histogram and experimental histogram bins same length.
    while len(norm_histogram) > len(sstr3_rab8a_exp_values):
        sstr3_rab8a_exp_values.append(0.0)
    
    ##Calculates mean square difference (MSD) between simulated and
    ##experimental histogram.
    msd = 0.0
    for i in range(len(norm_histogram)):
        msd += (norm_histogram[i] - sstr3_rab8a_exp_values[i])**2
    
    return msd
    
def random_copy_nums():
    """Samples a random copy number for SSTR3 and Rab8a on the surface of the 
    simulated vesicle."""
    
    ##Experimental copy number and probability.
    rab8a_copy_num_elements = [2, 3, 4, 5, 6, 7]
    rab8a_copy_num_weights = [0.07, 0.11, 0.23, 0.26, 0.24, 0.09]
    sstr3_copy_num_elememts = [3, 4, 5, 6, 7, 8, 9]
    sstr3_copy_num_weights = [0.05, 0.07, 0.26, 0.23, 0.21, 0.14, 0.04]
    
    ##Randomly samples a copy number for SSTR3 and Rab8a based on above
    ##probabilites.
    rab8a_sample = np.random.choice(rab8a_copy_num_elements, 
                                    p = rab8a_copy_num_weights)                      
    sstr3_sample = np.random.choice(sstr3_copy_num_elememts, 
                                    p = sstr3_copy_num_weights)

    return (rab8a_sample, sstr3_sample)    

##Output headers for final file.
simulation_msd = [('Vesicle Radius', 'MSD')]

##Number of distances simulated.
num_points = 1000

##Ciliary radius is 125 nm. Therefore, that is the limit for simualated
##vesicle size.
cilia_radius = 126

##Average and standard deviation of single molecule localization precision.
prec_avg, prec_std = (18.8, 6.7)

##Bins for the experimental distance data.
sstr3_rab8a_exp_bin_edges = [0, 20, 40, 60, 80, 100, 120]
sstr3_rab8a_exp_bin_mids = [10, 30, 50, 70, 90, 110]

for r in range(cilia_radius):
    print r
    temp_dist_data = []
    for n in range(num_points):
        
        ##Random copy numbers are sampled for each SSTR3 and Rab8a. Individual
        ##locations are simulated for each copy and averaged to simulated the 
        ##final point spread function of the sin
        rab8a_copy_num, sstr3_copy_num = random_copy_nums()
        rab8a_x = []
        rab8a_y = []
        sstr3_x = []
        sstr3_y = []
        for i in range(rab8a_copy_num):
            x, y, z = pt_on_sphere(radius = r)
            rab8a_x.append(x)
            rab8a_y.append(y)
        for i in range(sstr3_copy_num):
            x, y, z = pt_on_sphere(radius = r)
            sstr3_x.append(x)
            sstr3_y.append(y)
        
        ##Random localization precision sampled from a normal distribution
        ##of average prec_avg and standard deviation prec_std. Random
        ##localization precision is then added to simulated (x, y) location to
        ##model localization error.
        precision = np.max((0.0001, np.random.normal(prec_avg, prec_std)))
        rab8a_avg = (np.average(rab8a_x) + np.random.normal(0, precision), 
                     np.average(rab8a_y) + np.random.normal(0, precision))
                     
        precision = np.max((0.0001, np.random.normal(prec_avg, prec_std)))
        sstr3_avg = (np.average(sstr3_x) + np.random.normal(0, precision), 
                     np.average(sstr3_y) + np.random.normal(0, precision))
        
        ##Calculates final distance between two simulated SSTR3 and Rab8a
        ##point spread functions.             
        distance = np.sqrt((rab8a_avg[0] - sstr3_avg[0])**2 + 
                           (rab8a_avg[1] - sstr3_avg[1])**2)
                            
        temp_dist_data.append(distance)
        
    values, bin_edges = np.histogram(temp_dist_data, 
                                     bins = sstr3_rab8a_exp_bin_edges)
                                     
    simulation_msd.append((r, histogram_msd(values)))

with open('simulation.csv', 'wb') as csvfile:
    writer = csv.writer(csvfile, delimiter = ',')
    for row in simulation_msd:
        writer.writerow(row)