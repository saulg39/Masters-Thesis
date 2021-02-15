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
   
Buckling_list=[]
for row_number in range(19,100):
     if sheet.cell(row_number,0).value == 1:
          if row_number<73:
               row = sheet.row_slice(rowx=row_number,
                                   start_colx=28,
                                   end_colx=45)
               shape = "RHS"

               L, b, d, r, t_flange, t_web, c, spr, s_ult, n, Eel, s_pl, curve_pl = [row[0].value, row[2].value, row[1].value, row[4].value, row[3].value, row[3].value, row[3].value, row[9].value, row[10].value, row[13].value, row[12].value, row[15].value, row[16].value]
          else:
               row = sheet.row_slice(rowx=row_number,
                                   start_colx=10,
                                   end_colx=33)
               shape = "I Beam"

               L, b, d, r, t_flange, t_web, c, spr, s_ult, n, Eel, s_pl, curve_pl = [row[0].value, row[2].value, row[1].value, row[5].value, row[4].value, row[3].value, row[3].value, row[16].value, row[17].value, row[18].value, row[19].value, row[21].value, row[22].value]

          A_initial = 0.00001
          L_current = L
          ratio = 1
          M_crit=10**20
          while L_current>d/10:
               L_current= L/ratio
               M, A = finitestrip_shape(L_current,[shape, b, d, r, t_web, t_flange , c], [Eel , spr , n , 0.3, -0.46, s_ult], ["N", [Eel , spr , n , 0.3, -0.46, s_ult]], A_initial = A_initial)
               if M<M_crit:
                    M_crit=M
               ratio = math. ceil(ratio*1.2)
          Buckling_list.append([row_number,M_crit])
print(Buckling_list)















  

"""n = 6
max_L = 100000
min_L = 1
r = (max_L/min_L) ** (1/(n-1))
Moment = []
Moment_1 = []
Length = []
Length_factor= True
tic = time.perf_counter()

for i in range(n):
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

     M, A = finitestrip_shape(L, shape = "RHS", b = 25, d = 60, r = 0.01, t_flange = 1, t_web = 1, c = 12.7, Eel = 203000, spr = 579, n = 7.6, v = 0.3, k = -0.46, A_initial = A_initial, s_ult = 700)

     if M =="Fail":
          break
     if Length_factor:
          Moment_1.append(M)
          for n in range(2,round(L/Length[0])+1):
               M_2 = Moment[find_nearest(Length, L/n)]
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
plt.show()"""









"""
print(finitestrip_shape(L = 1000, shape = "channel", b = 138.60, d = 202.05, r = 5, t_flange = 6.11, t_web = 6.01, c = 23, Eel = 195000, spr = 520, n = 7.5, v = 0.3, k = -0.46))

x1=[]
y1=[]

for f in file:
     x1.append(f[0])
     y1.append(f[1])
plt.semilogx(x1, y1, linewidth = 0.4, color = "red", label='CUFSM')"""