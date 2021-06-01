import math
import numpy as np
from scipy import linalg
from channel import channel
from I_beam import I_beam
from RHS import RHS
from stress_from_strain import stress_from_strain
from get_data import get_data
import time
import matplotlib.pyplot as plt
import pandas as pd 
from finite_strip_moment import finitestrip_shape
import os
import xlrd
from xlutils.copy import copy
from xlwt import Workbook

#cd c:/Users/saulg/Documents/year\ 4/IIB_Project/code/mycode


def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx

Material, Forming_process, Bending_type, Actual_Length, shape, Material_flat, Material_corner = get_data(100)

print(Material)
print(Forming_process)
print(Bending_type)
print(shape)
print(Actual_Length)
print(Material_flat)
print(Material_corner)


M, A, max_stress = finitestrip_shape(1700, shape, Material_flat, Material_corner, 0.00001)
print(M, A)


print(shape[2])
y=np.array([])
x1=np.array([])
x2=np.array([])
for i in range(101):
     y = np.append(y,i*shape[2]/200)
     s,E =stress_from_strain(A * (i*shape[2]/200), Material_flat)
     print(A*i*shape[2]/200)
     x1 = np.append(x1,s)
     x2 = np.append(x2,E)

y_m = -1*y[::-1]
x_m1 = -1*x1[::-1]
x_m2 = x2[::-1]

x1 = np.concatenate((x_m1, x1), axis = None)
x2 = np.concatenate((x_m2, x2), axis = None)
y = np.concatenate((y_m, y), axis = None)
fig,ax = plt.subplots(1)


print("y")
for n in y:
     print(n)
print("stress")
for n in x1:
     print(n)
print("EEEEEEEEEEEEEEEEEEEEEEEt")
for n in x2:
     print(n)


ax.plot(x2,y)

# tell matplotlib which yticks to plot 
plt.show()