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
from finite_strip_Compression_new import finitestrip_shape
import os
import xlrd
from xlutils.copy import copy
from xlwt import Workbook

#C:\Users\saulg\Documents\year 4\IIB_Project\code

def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx



factor = 6
list_F =[]
list_strain =[]
list_stress = []
list_Et= []
list_type=[]

for data_point in range(237, 376):
     Material, Forming_process, Bending_type, Actual_Length, shape, Material_flat, Material_corner = get_data(data_point)
     if Bending_type == "C":

          number = 100
          max_L = Actual_Length*3
          min_L = shape[2]/10
          r = (max_L/min_L) ** (1/(number-1))
          Force = []
          stress_list = []
          Force_1 = []
          strain_list = []
          Et_list = []
          Length = []
          Length_factor= True
          tic = time.perf_counter()


          for i in range(number):
               L = min_L * r**i
               Length.append(L)
               if i == 0:
                    sg = Material_flat[1]
               elif i == 1:
                    sg = s1
                    s2 = s1
               else:
                    sg = abs(2 * s1 - s2)
                    s3 = s2
                    s2 = s1
               
               F, sg, strain, Et= finitestrip_shape(L, shape, Material_flat, Material_corner, sg)
               s1 = sg
               if Length_factor:
                    Force_1.append(F)
                    for n in range(2,round(factor*L/(shape[2]))):
                         F_2 = Force[find_nearest(Length, L/n)]
                         if F > F_2:
                              F = F_2
                              sg =stress_list[find_nearest(Length, L/n)]
                              strain = strain_list[find_nearest(Length, L/n)]
                              Et = Et_list[find_nearest(Length, L/n)]
               Force.append(F)
               stress_list.append(sg)
               strain_list.append(strain)
               Et_list.append(Et)
               

          typ = "Flexural"
          F, sg, strain, Et= finitestrip_shape(Actual_Length, shape, Material_flat, Material_corner, sg)

          if Length_factor and F !="Fail":
               for n in range(2,round(factor*Actual_Length/(shape[2]))):

                    F_2 = Force[find_nearest(Length, Actual_Length/n)]
                    if F > F_2:
                         F = F_2
                         typ = "Local"
                         sg =stress_list[find_nearest(Length, Actual_Length/n)]
                         strain = strain_list[find_nearest(Length, Actual_Length/n)]
                         Et = Et_list[find_nearest(Length, Actual_Length/n)]
          print(data_point, Actual_Length,F, sg, strain, Et)
          list_F.append(F)
          list_stress.append(sg)
          list_strain.append(strain)
          list_Et.append(Et)
          list_type.append(typ)

          plt.semilogx([Actual_Length], [F],'rx')
          if Length_factor:
               plt.semilogx(Length, Force_1, linewidth = 0.4, color = "black", label='Singe Half Wavelength')
               plt.semilogx(Length, Force, linewidth = 1.2, color = "black", label='Multiple Half Wavelengths')
          else:
               plt.semilogx(Length, Force, linewidth = 1.4, color = "black")

          #plt.semilogx([1000], [151], 'rx',  label='Buckling')
          #plt.semilogx([1000], [169], 'bx',  label='Failure')

          plt.xlabel('L / mm')
          plt.ylabel('Buckling Moment / KNm')
          plt.grid(True,'both')
          plt.title(data_point)
          plt.legend()
          #plt.show()
     else:
          list_F.append("")
          list_stress.append("")
          list_strain.append("")
          list_Et.append("")
print("Force")
for n in list_F:
     print(n)
print("Strain")
for n in list_strain:
     print(n)
print("list_stress")
for n in list_stress:
     print(n)
print("list_Et")
for n in list_Et:
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