#-------------------------------------------------------------------------------
# Name:        Point_Circle_Generator
#
# Author:      Mark Tingey
#
# Created:     23/01/2019
#-------------------------------------------------------------------------------
##import Statements
import random as rnd
import numpy as np
import math
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import csv
import os
import pandas as pd

##Define center point class
class Point:
    def __init__(self, x, y):
        self.x = x
        self.y = y

    def __str__(self):
        return str((self.x, self.y))

    def reset(self):
        self.__init__()

##define circle class
class Circle:
    def __init__(self, origin, radius):
        self.origin = origin
        self.radius = radius
    def reset(self):
        self.__init__()

##define x coordinate for circle center
def x_coordinate(radius, angle):
    x = radius * math.sin(math.radians(angle))
    return(x)

##define y coordinate for circle center
def y_coordinate(radius, angle):
    y = radius * math.cos(math.radians(angle))
    return(y)

##Update list of angles
def get_angle(angle_count):
    list_of_angles = [0,45,90,135,180,225,270,315]
    choose = np.random.choice(list_of_angles)
    rad = choose + angle_count
    return(rad)

##Global Variables
rnd.seed = 1234
np.random.seed(5)
print_angles = [0,45,90,135,180,225,270,315]
list_of_x = []
list_of_y = []
list_of_z = []
radius_of_pore = 35
t = 0
points = 1000
iteration = 0
csv_file = 0
angle_count = 0

##Create folder to hold data
path = os.getcwd()

data_folder = os.path.join(path, "Simulated_Dataset")
if os.path.isdir(data_folder)==False:
    os.mkdir("Simulated_Dataset")
os.chdir(data_folder)

##iteration for each angle position
while iteration != 45:

    ##perform simulation
    while t != points:
        ##define center point of the area

        angle = get_angle(angle_count)

        origin_x = x_coordinate(radius_of_pore, angle)
        origin_y = y_coordinate(radius_of_pore, angle)

        ##define circle
        origin = Point(origin_x,origin_y)
        radius = 13.394
        circle = Circle(origin, radius)

        ##generate point in circle
        p = rnd.random() * 2 * math.pi
        r = circle.radius * math.sqrt(rnd.random())
        x = math.cos(p) * r
        y = math.sin(p) * r

        x = origin_x + x
        y = origin_y + y
        z = 0

        ##apply gaussian error to each point
        error_x = np.random.normal(x, 10)
        error_y = np.random.normal(y, 10)
        error_z = np.random.normal(z, 10)

        ##write coordinates onto list
        list_of_x.append(error_x)
        list_of_y.append(error_y)
        list_of_z.append(error_z)

        ##counter
        t = t+1


    ##Write the coordinates to csv

    angle_count = angle_count + 1
    output = zip(list_of_x, list_of_y)

    outputname = ("simulation" + str(csv_file)+".csv")

    ##write file to csv
    with open(outputname, 'w', newline='') as writefile:
        writer = csv.writer(writefile)
        writer.writerows(zip(list_of_x, list_of_y))
    writefile.close()

    print(str(print_angles))

    ##reset class
    Point.reset(None)
    Circle.reset(None)


    ##change angle by 1 degree
    csv_file = csv_file + 1
    print_angles = [x+1 for x in print_angles]

    ##uncomment the below section to see a preview graph of generated points
##    plt.scatter(list_of_x,list_of_y, s = 1)
##    plt.show()

    iteration = iteration + 1
    t = 0
    list_of_x = []
    list_of_y = []
    list_of_z = []


print(str(print_angles))



