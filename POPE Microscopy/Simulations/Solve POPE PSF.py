#-------------------------------------------------------------------------------
# Name:        module2
# Purpose:
#
# Author:      angel
#
# Created:     05/05/2022
# Copyright:   (c) angel 2022
# Licence:     <your licence>
#-------------------------------------------------------------------------------

import math
import cmath

def SolveSigma (Sigma1,Sigma2):
    sqF = Sigma2**2
    sqA = Sigma1**2
    step1 = sqA*sqF
    step2 = sqA-sqF
    step3 = step1/step2
    final = cmath.sqrt(step3)
    return final

def SolveMu (Mu1,Mu2,Sigma1,SigmaO):
    var1 = Sigma1**2
    var2 = SigmaO**2
    den = var1+var2
    prod1 = Mu1*var2
    left = den*Mu2
    Nleft = left-prod1
    final = Nleft/var1
    return final


#Single PSF not using POPE module
Sigma1 = 150
Mu1 = 1


#combined PSF using POPE module
Sigma2 = 172.723101
Mu2 = 31.4


SigmaO = SolveSigma(Sigma1,Sigma2)
Mu0 = SolveMu(Mu1,Mu2,Sigma1,SigmaO)

print('Sigma: ' + str(SigmaO))
print('Mu: ' + str(Mu0))