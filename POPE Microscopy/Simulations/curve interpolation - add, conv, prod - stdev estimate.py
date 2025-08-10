# -*- coding: utf-8 -*-
"""
Created on Mon Oct 17 10:12:28 2022

@author: tuh05145
"""
from scipy.signal import convolve
from scipy.stats import norm
import numpy as np
import matplotlib.pyplot as plt

def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array-value)).argmin()
    return array[idx]

#data
mu1 =  0
mu2 = 0
sigma1 = 321.118
sigma2 =330.2

##bounds
lower = mu1 - 7*sigma1
upper = mu1 + 7*sigma1
points = 500000

x = np.linspace(lower, upper, points)
convFWHM = 2.35482004503
##normalization factor for convolution
dx = x[1] - x[0]
##Gaussian PDF 1
y = norm.pdf(x, loc=mu1, scale=sigma1)
##Gaussian PDF 2
y1 = norm.pdf(x, loc=mu2, scale=sigma2)
##Addition of gPDF1 and 2
z = y+y1
##Confolution of gPDF1 and 2
z1 = convolve(y,y1, mode='same', method='fft')*dx

##Calculate Sigma 3 and Mu3
sig_sqr1 = sigma1**2
sig_sqr2 = sigma2**2
mu_sqr1 = mu1**2
mu_sqr2 = mu2**2
sigma3 = np.sqrt((sig_sqr1 * sig_sqr2)/(sig_sqr1+sig_sqr2))
mu3 = ((mu1*sig_sqr2)+(mu2*sig_sqr1))/(sig_sqr1+sig_sqr2)

##Scale the gaussian product
sum_of_sigma_sq = (sigma1**2) + (sigma2 **2)
diff_of_mean = (mu1 - mu2)
diff_of_mean_sq = diff_of_mean**2
scale_factor = (1/np.sqrt(2 * np.pi * sum_of_sigma_sq)) * np.exp(-1 * diff_of_mean_sq / (2 * sum_of_sigma_sq))

##Alternative scale factor of gaussian product
base = 1/(np.sqrt(2*np.pi*(sig_sqr1*sig_sqr2)/sigma3**2))
exp = -(1/2)*((mu_sqr1/sig_sqr1)+(mu_sqr2/sig_sqr2)-((mu3**2)/(sigma3**2)))
alt_scale_factor = base * np.exp(exp)

##Scaled product of gPDF1 and 2
prod = (1/alt_scale_factor)*y*y1

##Plot original gPDFs
plt.subplot(2, 1, 1)

plt.plot(x, y, lw=3, alpha=0.4)
plt.plot(x, y1, lw=3, alpha=0.4)
plt.xlim([lower, upper])

#find the maximum value
##max  value of addition of two PDFs
max = z.argmax()
MaxVal = z[max]

##max value of convolution of two PDFs
max1 = z1.argmax()
MaxVal1 = z1[max1]

##max value of product of two PDFs
max2 = prod.argmax()
MaxVal2 = prod[max2]

#half maximum value
hMax = MaxVal/2 ##half max value of addition
hMax1 = MaxVal1/2 ##half max value of convolution
hMax2 = MaxVal2/2 ##half max value of the product

#index location of half maximum
#loc of addition
loc = find_nearest(z,hMax)
locHMAX = np.where(z == loc)
#loc of convolution
loc1 = find_nearest(z1,hMax1)
locHMAX1 = np.where(z1 == loc1)
#loc of product
loc2 = find_nearest(prod,hMax2)
locHMAX2 = np.where(prod == loc2)

#half width at half maximum
##addition
hwhm_loc = locHMAX[0][0]
halfA = abs(x[hwhm_loc]-x[max])
fwhm = halfA*2
stdev = fwhm/convFWHM
##Convolution
hwhm_loc1 = locHMAX1[0][0]
halfC = abs(x[hwhm_loc1]-x[max1])
fwhm1 = halfC*2
stdev1 = fwhm1/convFWHM

##Product
hwhm_loc2 = locHMAX2[0][0]
halfP = abs(x[hwhm_loc2]-x[max2])
fwhm2 = halfP*2
stdev2 = fwhm2/convFWHM

#roll new curve to mu located between mu1 and mu2
space = abs(lower-upper) #provides the length of the array
perNM = points/space #converts array position to nm

##mean of Product
meanProd = ((mu1*(sigma2**2))+(mu2*(sigma1**2)))/((sigma1**2)+(sigma2**2))

#if statment for shifting convolved PSF if Mu1 and 2 are identical
if mu1 == mu2:
    mu3 = mu1
else:
    mu3 = np.mean([mu1,mu2]) #provides mu for convolved array

#don't shift if mu1==mu2
if mu1 == mu2:
    shift = 0
else:
    pos1 = int(mu3*perNM) #provides the number of positions array needs to shift per nm
    pos2 = int(mu2*perNM)
    shift = abs(pos1-pos2)
    shift = int(shift)

#shift mu3 if mu2 is greater than mu1
if mu2 >= mu1:
    shift = int((shift*-1))
else:
    shift = shift

convolved = np.roll(z1, shift)
print("----------Addition of two Gaussians----------")
print("FWHM added distribution", fwhm)
print("stdev added distribution", stdev)
print()
print("----------Convolution of two Gaussians----------")
print("FWHM convolved distribution", fwhm1)
print("stdev convolved distribution", stdev1)
print()
print("----------Product of two Gaussians----------")
print("FWHM product distribution", fwhm2)
print("stdev product distribution", stdev2)
print("the mu of the product:", meanProd)


plt.subplot(2, 1, 2)
plt.plot(x, y, lw=3, alpha=0.5)
plt.plot(x, y1, lw=3, alpha=0.5)
plt.plot(x, z, lw=3, alpha=1, label = "added PSF")
plt.plot(x, convolved, lw=3, alpha=1,label="convolved PSF")
plt.plot(x, prod, lw=3, alpha=1, label="product PSF")
plt.legend()
plt.xlim([lower, upper])

plt.tight_layout()
plt.show()