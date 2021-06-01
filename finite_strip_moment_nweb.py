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

factor = 5

list_M =[]
list_C =[]
list_stress = []
list_stress_c = []
list_Et= []
list_Et_c= []
list_type = []

for data_point in range(104,105):
     Material, Forming_process, Bending_type, Actual_Length, shape, Material_flat, Material_corner = get_data(data_point)
     if Bending_type == 4:
          """print(Material)
          print(Forming_process)
          print(Bending_type)
          print(shape)
          print(Actual_Length)
          print(Material_flat)
          print(Material_corner)"""

          number = 100
          max_L = Actual_Length*5 
          min_L = shape[2]/5
          r = (max_L/min_L) ** (1/(number-1))
          Moment = []
          stress_list = []
          Moment_1 = []
          Curvature_1 = []
          Curvature = []
          Length = []
          Length_factor= True
          tic = time.perf_counter()
          A=0.00001
          for i in range(100,101):
               L = Actual_Length
               Length.append(L)
               A_initial = A
               

               M, A, max_stress = finitestrip_shape(L, shape, Material_flat, Material_corner, A_initial, n=i)
               C=A
               if M =="Fail":
                    Length.pop()
                    break
               print(i)
               Moment.append(M)
               Curvature.append(C)
               stress_list.append(max_stress)
          toc = time.perf_counter()
          

print("M")
for n in Moment:
     print( n)







"""
print(finitestrip_shape(L = 1000, shape = "channel", b = 138.60, d = 202.05, r = 5, t_flange = 6.11, t_web = 6.01, c = 23, Eel = 195000, spr = 520, n = 7.5, v = 0.3, k = -0.46))

x1=[]
y1=[]

for f in file:
     x1.append(f[0])
     y1.append(f[1])
plt.semilogx(x1, y1, linewidth = 0.4, color = "red", label='CUFSM')"""