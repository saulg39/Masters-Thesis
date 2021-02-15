import math
import numpy as np
from scipy import linalg
from channel import channel
from I_beam import I_beam
from RHS import RHS
from stress_from_strain import stress_from_strain
import time
import matplotlib.pyplot as plt
import pandas as pd 
from finite_strip_moment import finitestrip_shape
import os
import xlrd
from xlutils.copy import copy
from xlwt import Workbook

#C:\Users\saulg\Documents\year 4\IIB_Project\code

def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx

THIS_FOLDER = os.path.dirname(os.path.abspath(__file__))
my_file = os.path.join(THIS_FOLDER, 'Data','beam_column_data.xlsx')
book = xlrd.open_workbook(my_file)
# get the first worksheet
sheet = book.sheet_by_index(1)
   # read a row slice
   
row_number = 116

row = sheet.row_slice(rowx=row_number-1,
                    start_colx=0,
                    end_colx=40)


Forming_process = row[4].value
Bending_type = row[5].value
Actual_Length = row[6].value
shape = []
for i in range(7,14):
     shape.append(row[i].value)
Material_flat = []
for i in range(14,20):
     Material_flat.append(row[i].value)
Material_corner = [row[20].value,[]]
for i in range(21,27):
     Material_corner[1].append(row[i].value)

print(Forming_process)
print(Bending_type)
print(shape)
print(Actual_Length)
print(Material_flat)
print(Material_corner)

number = 500
max_L = 20000
min_L = 20
r = (max_L/min_L) ** (1/(number-1))
Moment = []
Moment_1 = []
Length = []
Length_factor= True
tic = time.perf_counter()

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

     M, A = finitestrip_shape(L, shape, Material_flat, Material_corner, A_initial)

     if M =="Fail":
          break
     if Length_factor:
          Moment_1.append(M)
          for num in range(2,round(L/Length[0])+1):
               M_2 = Moment[find_nearest(Length, L/num)]
               if M > M_2:
                    M = M_2

     Moment.append(M)
toc = time.perf_counter()
print(f"Done in {toc - tic:0.4f} seconds")
M, A = finitestrip_shape(Actual_Length, shape, Material_flat, Material_corner, A_initial)

if Length_factor and M !="Fail":
     for num in range(2,round(Actual_Length/Length[0])+1):
          M_2 = Moment[find_nearest(Length, Actual_Length/num)]
          if M > M_2:
               M = M_2
print(Actual_Length,M)
plt.semilogx([Actual_Length], [M],'rx')
if Length_factor:
     plt.semilogx(Length, Moment_1, linewidth = 0.4, color = "black", label='Singe Half Wavelength')
     plt.semilogx(Length, Moment, linewidth = 1.2, color = "black", label='Multiple Half Wavelengths')
else:
     plt.semilogx(Length, Moment, linewidth = 1.4, color = "black")

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