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

Material, Forming_process, Bending_type, Actual_Length, shape, Material_flat, Material_corner = get_data(126)


number = 150
max_L = 10000
min_L = 5
r = (max_L/min_L) ** (1/(number-1))
Moment = []
stress_list = []
Moment_1 = []
Length = []
Length_factor= False
Length_list=[]
Moment_list=[]
tic = time.perf_counter()
for delta in [0,1/300,1/100,1/10]:
     Length=[]
     Moment=[]
     for i in range(number):
          L = min_L * r**i
          Length.append(L)
          if i == 0:
               A_initial = 0.00001
          elif i == 1:
               A_initial = A
               A_2 = A
          elif i == 2:
               A_initial = abs(2 * A - A_2)
               A_3 = A_2
               A_2 = A
          else:
               A_initial = abs(3 * A - 3 * A_2 + A_3)
               A_3 = A_2
               A_2 = A

          M, A, max_stress = finitestrip_shape(L, shape, Material_flat, Material_corner, A_initial,delta)

          if M =="Fail":
               Length.pop()
               break
          if Length_factor:
               Moment_1.append(M)
               for num in range(2,round(L/Length[0])+1):
                    M_2 = Moment[find_nearest(Length, L/num)]
                    if M > M_2:
                         M = M_2
                         max_stress = stress_list[find_nearest(Length, L/num)]

          Moment.append(M)
          stress_list.append(max_stress)
     Moment_list.append(Moment)
     Length_list.append(Length)

toc = time.perf_counter()
print(f"Done in {toc - tic:0.4f} seconds")
c=["black","red","blue","green","red","blue","green","red","blue","green"]
labels = ["No imperfections","D/300","D/100","D/10"]
for n in range(len(Moment_list)):

     plt.semilogx(Length_list[n], Moment_list[n], linewidth = 0.3, color = c.pop(0),label=labels.pop(0))

#plt.semilogx([1000], [151], 'rx',  label='Buckling')
#plt.semilogx([1000], [169], 'bx',  label='Failure')

plt.xlabel('L / mm')
plt.ylabel('Moment / KNm')
plt.grid(True,'both')
plt.legend()
plt.show()









"""
print(finitestrip_shape(L = 1000, shape = "channel", b = 138.60, d = 202.05, r = 5, t_flange = 6.11, t_web = 6.01, c = 23, Eel = 195000, spr = 520, n = 7.5, v = 0.3, k = -0.46))

x1=[]
y1=[]

for f in file:
     x1.append(f[0])
     y1.append(f[1])
plt.semilogx(x1, y1, linewidth = 0.4, color = "red", label='CUFSM')"""