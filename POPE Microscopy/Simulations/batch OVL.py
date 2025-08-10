#-------------------------------------------------------------------------------
# Name:        module1
# Purpose:
#
# Author:      angel
#
# Created:     11/04/2022
# Copyright:   (c) angel 2022
# Licence:     <your licence>
#-------------------------------------------------------------------------------

from statistics import NormalDist
import pandas as pd
import csv
import os



curves = pd.read_csv ('OVL.csv')

#determine X OVL

curves["X OVL"] = 0
curves["Y OVL"] = 0
runs = len(curves)
X_clock = 0
y_clock = 0

while X_clock < runs:
    curves.at[X_clock,"X OVL"] = NormalDist(curves.iat[X_clock,0], curves.iat[X_clock,2]).overlap(NormalDist(curves.iat[X_clock,4],curves.iat[X_clock,6]))
    X_clock = X_clock + 1
    print(X_clock)


while y_clock < runs:
    curves.at[y_clock,"Y OVL"] = NormalDist(curves.iat[y_clock,1], curves.iat[y_clock,3]).overlap(NormalDist(curves.iat[y_clock,5],curves.iat[y_clock,7]))
    y_clock = y_clock + 1
    print(y_clock)

curves.to_csv("OVL Calculated.csv", encoding = 'utf-8', index=False)
