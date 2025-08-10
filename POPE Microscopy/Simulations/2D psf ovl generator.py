import numpy as np
import matplotlib.pyplot as plt
from statistics import NormalDist
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.colors import ListedColormap

def OVLscore (mu1, mu2, sigma1, sigma2):
    OVL = NormalDist(mu1, sigma1).overlap(NormalDist(mu2, sigma2))
    return(OVL)

def overlap(OVL,mu,sigma):
    """estimates the position of the POPE PSF so the OVL is equal to the required value"""

    mu1 = mu
    sigma1 = sigma
    mu2 = mu
    sigma2 = (sigma*1.15)
    step = 0.001
    count = 1

    result = OVLscore(mu1,mu2,sigma1,sigma2)
    while result >= OVL:
        mu2 = mu2 + step
        result = NormalDist(mu1, sigma1).overlap(NormalDist(mu2, sigma2))
        #print(result)
        #print(str(count))
        count = count+1
    return(mu2,sigma2)



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

# set the degree of overlap required as well as position and matrix for reference PSF
OVL = 0.75
muX = 0.0
muY = 0.0
sigmaTest = 0.7
Shape = 0.0

# Our 2-dimensional distribution will be over variables X and Y
N = 400
X = np.linspace(-5, 5, N)
Y = np.linspace(-5, 5, N)
X, Y = np.meshgrid(X, Y)

#collect x, Y, and Sigma for 2D gaussian
posX = np.array(overlap(OVL,muX,sigmaTest))
posY = np.array(overlap(OVL,muY,sigmaTest))
varSig = posX[1]

# Mean vector and covariance matrix
mu = np.array([muX, muY])
Sigma = np.array([[sigmaTest, 0], [0, sigmaTest]])
mu2 = np.array([posX[0],posY[0]])
Sigma2 = np.array([[varSig,Shape], [Shape,varSig]])

# Pack X and Y into a single 3-dimensional array mu1
pos = np.empty(X.shape + (2,))
pos[:, :, 0] = X
pos[:, :, 1] = Y

# The distribution on the variables X, Y packed into pos.
Z = multivariate_gaussian(pos, mu, Sigma)
Z2 = multivariate_gaussian(pos, mu2, Sigma2)

##print covariance matrix for simulated PSF
print("Mu for simulated PSF:")
print(mu2)
print("Shape of simulated PSF:")
print(Sigma2)

# plot figures
#plot 2D contour
fig = plt.figure(dpi=1200)

#plot 3D curve
fig2 = plt.figure(dpi=1200)

ax1 = fig2.add_subplot(1,1,1,projection='3d')


ax1.plot_surface(X, Y, Z, rstride=3, alpha=.4, cstride=3, linewidth=1, antialiased=True, cmap=cm.gist_gray)
ax1.plot_surface(X, Y, Z2, rstride=3, alpha=.6, cstride=3, linewidth=1, antialiased=True, cmap=cm.gist_gray)
ax1.view_init(55,-70)
ax1.set_xticks([])
ax1.set_yticks([])
ax1.set_zticks([])
ax1.view_init()


ax2 = fig.add_subplot()

a=ax2.contour(X, Y, Z, alpha=1, linewidths=0.5, colors='#5B9BD5')
b=ax2.contour(X, Y, Z2, alpha=1, linewidths=0.5,colors='#ED7D31')

fig.gca().set_aspect('equal', adjustable='box')
ax2.grid(False)
ax2.set_xticks([])
ax2.set_yticks([])
ax2.set_xlabel(r'$X$')
ax2.set_ylabel(r'$Y$')

plt.show()