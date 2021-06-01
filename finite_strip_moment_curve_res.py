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
from finite_strip_moment_res import finitestrip_shape_res
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

for data_point in [196]:
     Material, Forming_process, Bending_type, Actual_Length, shape, Material_flat, Material_corner = get_data(data_point)
     if Bending_type == 4:
          """print(Material)
          print(Forming_process)
          print(Bending_type)
          print(shape)
          print(Actual_Length)
          print(Material_flat)
          print(Material_corner)"""

          number = 500
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

               M, A, max_stress = finitestrip_shape(L, shape, Material_flat, Material_corner, A_initial)
               C=A
               if M =="Fail":
                    Length.pop()
                    break
               if Length_factor:
                    Moment_1.append(M)
                    Curvature_1.append(M)
                    
                    for num in range(2,round(factor*L/(shape[2]))):
                         M_2 = Moment[find_nearest(Length, L/num)]
                         C_2 = Curvature[find_nearest(Length, L/num)]
                         if M > M_2:
                              M = M_2
                              max_stress = stress_list[find_nearest(Length, L/num)]
                         if C > C_2:
                              C = C_2
 
               Moment.append(M)#
               Curvature.append(C)
               stress_list.append(max_stress)
          toc = time.perf_counter()
          print(f"Done in {toc - tic:0.4f} seconds")
          typ = "Flexural"
          M, A, max_stress= finitestrip_shape(Actual_Length, shape, Material_flat, Material_corner, A_initial)
          C=A
          if Length_factor and M !="Fail":
               for num in range(2,round(factor*Actual_Length/(shape[2]))):

                    M_2 = Moment[find_nearest(Length, Actual_Length/num)]
                    C_2 = Curvature[find_nearest(Length, Actual_Length/num)]
                    if M_2<M:
                         M=M_2
                         typ = "Local"
                         max_stress = stress_list[find_nearest(Length, Actual_Length/num)]
                    if C > C_2:
                         C = C_2
          print(data_point, Actual_Length,M,C, max_stress)
          list_M.append(M)
          list_C.append(C)
          list_stress.append(max_stress[0][0])
          list_stress_c.append(max_stress[1][0])
          list_Et.append(max_stress[0][1])
          list_Et_c.append(max_stress[1][1])
          list_type.append(typ)

          plt.semilogx([Actual_Length], [M],'kx')
          if Length_factor:
               plt.semilogx(Length, Moment_1, linewidth = 0.4, color = "black", label='No Residual Stresses Singe Half Wavelength')
               plt.semilogx(Length, Moment, linewidth = 1.2, color = "black", label='No Residual Stresses Multiple Half Wavelengths')
          else:
               plt.semilogx(Length, Moment, linewidth = 1.4, color = "black")

          number = 500
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

               M, A, max_stress = finitestrip_shape_res(L, shape, Material_flat, Material_corner, A_initial)
               C=A
               if M =="Fail":
                    Length.pop()
                    break
               if Length_factor:
                    Moment_1.append(M)
                    Curvature_1.append(M)
                    
                    for num in range(2,round(factor*L/(shape[2]))):
                         M_2 = Moment[find_nearest(Length, L/num)]
                         C_2 = Curvature[find_nearest(Length, L/num)]
                         if M > M_2:
                              M = M_2
                              max_stress = stress_list[find_nearest(Length, L/num)]
                         if C > C_2:
                              C = C_2
 
               Moment.append(M)#
               Curvature.append(C)
               stress_list.append(max_stress)
          toc = time.perf_counter()
          print(f"Done in {toc - tic:0.4f} seconds")
          typ = "Flexural"
          M, A, max_stress= finitestrip_shape_res(Actual_Length, shape, Material_flat, Material_corner, A_initial)
          C=A
          if Length_factor and M !="Fail":
               for num in range(2,round(factor*Actual_Length/(shape[2]))):

                    M_2 = Moment[find_nearest(Length, Actual_Length/num)]
                    C_2 = Curvature[find_nearest(Length, Actual_Length/num)]
                    if M_2<M:
                         M=M_2
                         typ = "Local"
                         max_stress = stress_list[find_nearest(Length, Actual_Length/num)]
                    if C > C_2:
                         C = C_2
          print(data_point, Actual_Length,M,C, max_stress)
          list_M.append(M)
          list_C.append(C)
          list_stress.append(max_stress[0][0])
          list_stress_c.append(max_stress[1][0])
          list_Et.append(max_stress[0][1])
          list_Et_c.append(max_stress[1][1])
          list_type.append(typ)

          plt.semilogx([Actual_Length], [M],'rx')
          if Length_factor:
               plt.semilogx(Length, Moment_1, linewidth = 0.4, color = "red", label='Residual Stresses Singe Half Wavelength')
               plt.semilogx(Length, Moment, linewidth = 1.2, color = "red", label='Residual Stresses Multiple Half Wavelengths')
          else:
               plt.semilogx(Length, Moment, linewidth = 1.4, color = "black")

          plt.xlabel('L / mm')
          plt.ylabel('Buckling Moment / KNm')
          plt.grid(True,'both')
          plt.title(data_point)
          plt.legend()
          plt.show()
     else:
          list_M.append("")
          list_C.append("")
          list_stress.append("")
          list_stress_c.append("")
          list_Et.append("")
          list_Et_c.append("")
          list_type.append("")
print("M")
for n in list_M:
     print(n)
print("list_C")
for n in list_C:
     print(n)
print("list_stress")
for n in list_stress:
     print(n)
print("list_stress_c")
for n in list_stress_c:
     print(n)
print("list_Et")
for n in list_Et:
     print(n)
print("list_Et_c")
for n in list_Et_c:
     print(n)
print("type")
for n in list_type:
     print(n)






"""
print(finitestrip_shape(L = 1000, shape = "channel", b = 138.60, d = 202.05, r = 5, t_flange = 6.11, t_web = 6.01, c = 23, Eel = 195000, spr = 520, n = 7.5, v = 0.3, k = -0.46))

x1=[]
y1=[]

for f in file:
     x1.append(f[0])
     y1.append(f[1])
plt.semilogx(x1, y1, linewidth = 0.4, color = "red", label='CUFSM')"""