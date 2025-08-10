# -*- coding: utf-8 -*-
"""
Created on Thu May 12 13:00:58 2022

@author: tuh05145
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.colors import ListedColormap


# Our 2-dimensional distribution will be over variables X and Y
N = 500
X = np.linspace(145, 155, N)
Y = np.linspace(10, 20, N)
X, Y = np.meshgrid(X, Y)

# Mean vector and covariance matrix
mu = np.array([150.3959, 14.1567])
Sigma = np.array([[0.79711, 0], [0, 0.79711]])
mu2 = np.array([150.7162, 13.68461])
Sigma2 = np.array([[0.6469,0], [0,0.6469]])
mu3 = np.array([150.589, 13.8677])
Sigma3 = np.array([[0.64691,0],[0,0.64691]])

# Pack X and Y into a single 3-dimensional array mu1
pos = np.empty(X.shape + (2,))
pos[:, :, 0] = X
pos[:, :, 1] = Y


def multivariate_gaussian(pos, mu, Sigma):
    """Return the multivariate Gaussian distribution on array pos."""

    n = mu.shape[0]
    Sigma_det = np.linalg.det(Sigma)
    Sigma_inv = np.linalg.inv(Sigma)
    N = np.sqrt((2*np.pi)**n * Sigma_det)
    # This einsum call calculates (x-mu)T.Sigma-1.(x-mu) in a vectorized
    # way across all the input variables.
    fac = np.einsum('...k,kl,...l->...', pos-mu, Sigma_inv, pos-mu)

    return np.exp(-fac / 2) / N

# The distribution on the variables X, Y packed into pos.
Z = multivariate_gaussian(pos, mu, Sigma)
Z2 = multivariate_gaussian(pos, mu2, Sigma2)
Z3 = multivariate_gaussian(pos, mu3, Sigma3)

# plot using subplots
fig = plt.figure(dpi=900)
ax1 = fig.add_subplot(2,1,1,projection='3d')


##ax1.plot_surface(X, Y, Z, rstride=3, alpha=.4, cstride=3, linewidth=3, antialiased=True, cmap=cm.Blues)
##ax1.plot_surface(X, Y, Z2, rstride=3, alpha=.4, cstride=3, linewidth=3, antialiased=False, cmap=cm.Oranges)
c=ax1.contour(X, Y, Z3, alpha=1, zdir='z', offset=0, linewidths=.4, colors='#000000')
a=ax1.contour(X, Y, Z, alpha=1, zdir='z', offset=0, linewidths=.2, colors='#5B9BD5')
b=ax1.contour(X, Y, Z2, alpha=.6, zdir='z', offset=0, linewidths=.2, colors='#ED7D31')

##ax1.view_init(55,-70)
ax1.view_init(90, 270)
ax1.set_xticks([])
ax1.set_yticks([])
ax1.set_zticks([])
#ax1.set_xlabel(r'$x_1$')
#ax1.set_ylabel(r'$x_2$')

ax2 = fig.add_subplot(2,1,2,projection='3d')

a=ax2.contour(X, Y, Z, alpha=1, zdir='z', offset=0, linewidths=.2, colors='#5B9BD5')
b=ax2.contour(X, Y, Z2, alpha=.6, zdir='z', offset=0, linewidths=.2, colors='#ED7D31')
ax2.view_init(90, 270)

ax2.grid(False)
ax2.set_xticks([])
ax2.set_yticks([])
ax2.set_zticks([])
#ax2.set_xlabel(r'$x_1$')
#ax2.set_ylabel(r'$x_2$')

plt.show()