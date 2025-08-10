#-------------------------------------------------------------------------------
# Name:        module1
# Purpose:
#
# Author:      Mark Tingey
#
# Created:     12/10/2022
# Copyright:   (c) angel 2022
# Licence:     <your licence>
#-------------------------------------------------------------------------------

import math

def Precision(F,s,a,b,N):
    left = (16*(s**2+(a**2/12)))/(9*N)
    #print(left)
    right = ((8*math.pi*b**2)*(s**2+(a**2/12))**2)/(a**2*N**2)
    #print(right)
    precision = math.sqrt(F*(left+right))
    return(precision)

def Moving(s,D,t):
    MovingW = math.sqrt(s**2 +((1/3)*D*t))
    return(MovingW)

## Variables
F = 2
a10 = 1600
a60 = 267
N = 100000
N_POPE = N*2
b = 2
b_POPE = b*1
D = 0.5
t = 2

##PSF stdev for 60x WFFM
stdev_WFFM = 200
stdev_WFFM_POPE = 258

##PSF stdev for 10x WFFM
stdev_WFFM10 = 382.5
stdev_WFFM10_POPE = 270.46
##WFFM precision
WFFM60 = Precision(F,stdev_WFFM, a60,b,N)

WFFM10 = Precision(F,stdev_WFFM10, a10,b,N)


##WFFM precision
WFFM_POPE60 = Precision(F,stdev_WFFM_POPE, a60, b_POPE, N_POPE)

WFFM_POPE10 = Precision(F,stdev_WFFM10_POPE, a10, b_POPE, N_POPE)


print("---------------Stationary Particle---------------")
print("WFFM Localization Precision at 10x:")
print("WFFM Localization Precision 10x: ", WFFM10)
print("WFFM POPE Localization Precision 10x: ", WFFM_POPE10)
print("WFFM fold increase 10x: ", WFFM10/WFFM_POPE10)
print("")
print("WFFM Localization Precision at 60x:")
print("WFFM Localization Precision 60x: ", WFFM60)
print("WFFM POPE Localization Precision 60x: ", WFFM_POPE60)
print("WFFM fold increase 60x: ", WFFM60/WFFM_POPE60)
print("")

##Test for equality
CorrectionF = 0.001
test_WFFM_POPE10 = WFFM_POPE10
test_WFFM_POPE60 = WFFM_POPE60
count10 = 0
count60 = 0

PSF10 = stdev_WFFM_POPE
PSF60 = stdev_WFFM_POPE

if test_WFFM_POPE10 > WFFM10:
    count10 = count10 + 1

if test_WFFM_POPE60 > WFFM60:
    count60 = count60 + 1

while test_WFFM_POPE10 <= WFFM10:
    PSF10 = PSF10 + CorrectionF
    test_WFFM_POPE10 = Precision(F,PSF10,a10,b_POPE,N_POPE)

while test_WFFM_POPE60 <= WFFM60:
    PSF60 = PSF60 + CorrectionF
    test_WFFM_POPE60 = Precision(F,PSF60,a60,b_POPE,N_POPE)

##how much can the PSF increase before being equal to No-POPE and what is the Max PSF size to be equal to No-POPE
increase10 = PSF10-stdev_WFFM
Max_PSF10 = increase10 + stdev_WFFM
ratio10 = PSF10/stdev_WFFM
increase60 = PSF60-stdev_WFFM
Max_PSF60 = increase60 + stdev_WFFM
ratio60 = PSF60/stdev_WFFM

## how much can the PSF StDev increase before it is equal?
dif_PSF10 = Max_PSF10 - stdev_WFFM10_POPE
dif_PSF60 = Max_PSF60 - stdev_WFFM_POPE

if count10 == 1:
    print("-------Stationary particle 10X POPE localization precision is larger than No-POPE precision-------")
else:
    print("POPE and no-POPE in 10x are equal if the stdev increases by: ", ratio10*100,"%")
    print("10x maximum PSF stdev size for POPE is: ", Max_PSF10)
    print("10x StDev nm available increase: ", dif_PSF10,"nm" )
print("")
if count60 == 1:
    print("-------Stationary particle 60X POPE localization precision is larger than No-POPE precision-------")
else:
    print("POPE and no-POPE in 60x are equal if the stdev increases by: ", ratio60*100,"%")
    print("60x maximum PSF stdev size for POPE is: ", Max_PSF60)
    print("60x StDev nm available increase: ", dif_PSF60,"nm" )
print("")



##Moving Particles
moving = Moving(stdev_WFFM, D, t)
moving10 = Moving(stdev_WFFM10, D, t)

moving_POPE = Moving(stdev_WFFM_POPE, D, t)
moving_POPE10 = Moving(stdev_WFFM10_POPE, D, t)

Moving_WFFM10 = Precision(F,moving10, a10,b,N)
Moving_WFFM60 = Precision(F,moving, a60, b, N)
Moving_WFFM10_POPE = Precision(F,moving_POPE10,a10,b_POPE,N_POPE)
Moving_WFFM60_POPE = Precision(F,moving_POPE,a60,b_POPE,N_POPE)
Diff10 = Moving_WFFM10 - WFFM10
Diff60 = Moving_WFFM60 - WFFM60
Diff10_POPE = Moving_WFFM10 - Moving_WFFM10_POPE
Diff60_POPE = Moving_WFFM60 - Moving_WFFM60_POPE

##moving particle equality test
move_PSF10 = moving_POPE
move_PSF60 = moving_POPE
movetest_WFFM_POPE10 = Moving_WFFM10_POPE
movetest_WFFM_POPE60 = Moving_WFFM60_POPE
move_count10 = 0
move_count60 = 0

if movetest_WFFM_POPE10 > Moving_WFFM10:
    move_count10 = move_count10 + 1

if movetest_WFFM_POPE60 > Moving_WFFM60:
    move_count60 = move_count60 + 1

while movetest_WFFM_POPE10 <= Moving_WFFM10:
    move_PSF10 = move_PSF10 + CorrectionF
    #print(move_PSF10)
    movetest_WFFM_POPE10 = Precision(F,move_PSF10,a10,b_POPE,N_POPE)

while movetest_WFFM_POPE60 <= Moving_WFFM60:
    move_PSF60 = move_PSF60 + CorrectionF
    movetest_WFFM_POPE60 = Precision(F,move_PSF60,a60,b_POPE,N_POPE)

##how much can the PSF increase before being equal to No-POPE and what is the Max PSF size to be equal to No-POPE
move_increase10 = move_PSF10-moving
move_increase60 = move_PSF60-moving
moving_Max_PSF10 = move_increase10 + moving
moving_Max_PSF60 = move_increase60 + moving
move_ratio10 = move_PSF10/moving
move_ratio60 = move_PSF60/moving

## how much can the PSF StDev increase before it is equal?
dif_moving_PSF10 = moving_Max_PSF10 - stdev_WFFM10_POPE
dif_moving_PSF60 = moving_Max_PSF60 - stdev_WFFM_POPE

print("---------------Moving Particle---------------")
print("WFFM Localization Precision at 10x:")
print("WFFM Localization Precision 10x: ", Moving_WFFM10)
print("WFFM POPE Localization Precision 10x: ", Moving_WFFM10_POPE)
print("WFFM fold increase 10x: ", Moving_WFFM10/Moving_WFFM10_POPE)
print("")
print("WFFM Localization Precision at 60x:")
print("WFFM Localization Precision 60x: ", Moving_WFFM60)
print("WFFM POPE Localization Precision 60x: ", Moving_WFFM60_POPE)
print("WFFM fold increase 60x: ", Moving_WFFM60/Moving_WFFM60_POPE)
print("")
if move_count10 == 1:
    print("-------moving particle 10X POPE localization precision is larger than No-POPE precision-------")
else:
    print("Moving POPE and no-POPE in 10x are equal if the stdev increases by: ", move_ratio10*100,"%")
    print("10x maximum PSF stdev size for POPE is: ", moving_Max_PSF10)
    print("10x StDev nm available increase: ", dif_moving_PSF10,"nm" )
print("")
if move_count60 == 1:
    print("-------moving particle 60X POPE localization precision is larger than No-POPE precision-------")
else:
    print("Moving POPE and no-POPE in 60x are equal if the stdev increases by: ", move_ratio60*100,"%")
    print("60x maximum PSF stdev size for POPE is: ", moving_Max_PSF60)
    print("60x StDev nm available increase: ", dif_moving_PSF60,"nm" )
print("")
print("Difference between 10x moving and stationary: ", Diff10)
print("Difference between 60x moving and stationary: ", Diff60)
print("Difference between 10x moving, POPE and No-POPE: ", Diff10_POPE)
print("Difference between 60x moving, POPE and No-POPE: ", Diff60_POPE)
