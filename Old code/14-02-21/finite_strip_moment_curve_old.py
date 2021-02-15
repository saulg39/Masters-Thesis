import math
import numpy as np
from scipy import linalg
from channel import channel
from I_beam import I_beam
from RHS import RHS
from stress_from_strain_old import stress_from_strain
import time
import matplotlib.pyplot as plt
import pandas as pd 
from finite_strip_moment_old import finitestrip_shape
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
   
row_number = 64
if row_number<73:
     row = sheet.row_slice(rowx=row_number-1,
                         start_colx=28,
                         end_colx=45)
     shape = "RHS"

     L, b, d, r, t_flange, t_web, c, spr, s_ult, n, Eel, s_pl, curve_pl = [row[0].value, row[2].value, row[1].value, row[4].value, row[3].value, row[3].value, row[3].value, row[9].value, row[10].value, row[13].value, row[12].value, row[15].value, row[16].value]
else:
     row = sheet.row_slice(rowx=row_number-1,
                         start_colx=10,
                         end_colx=33)
     shape = "I Beam"

     L, b, d, r, t_flange, t_web, c, spr, s_ult, n, Eel, s_pl, curve_pl = [row[0].value, row[2].value, row[1].value, row[5].value, row[4].value, row[3].value, row[3].value, row[16].value, row[17].value, row[18].value, row[19].value, row[21].value, row[22].value]

          
A_initial=0.0001
print(L,finitestrip_shape(L = L, shape = shape, b = b, d = d, r = r, t_flange = t_flange, t_web = t_web, c = c, Eel = Eel, spr = spr, n = n, v = 0.3, k = -0.46, A_initial = A_initial, s_ult = s_ult))
"""
number = 300
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

     M, A = finitestrip_shape(L = L, shape = shape, b = b, d = d, r = r, t_flange = t_flange, t_web = t_web, c = c, Eel = Eel, spr = spr, n = n, v = 0.3, k = -0.46, A_initial = A_initial, s_ult = s_ult)

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








"""
print(finitestrip_shape(L = 1000, shape = "channel", b = 138.60, d = 202.05, r = 5, t_flange = 6.11, t_web = 6.01, c = 23, Eel = 195000, spr = 520, n = 7.5, v = 0.3, k = -0.46))

x1=[]
y1=[]

for f in file:
     x1.append(f[0])
     y1.append(f[1])
plt.semilogx(x1, y1, linewidth = 0.4, color = "red", label='CUFSM')"""