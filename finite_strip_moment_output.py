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
   
for row_number in [100,101,102,103,108,109,110,111,116,117,118]:

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

     A_initial = 0.00001
     L_current = Actual_Length
     ratio = 1
     M_crit=10**20
     while L_current>shape[1]/20:
          L_current= Actual_Length/ratio
          M, A_initial, max_stress = finitestrip_shape(L_current, shape, Material_flat, Material_corner, A_initial)
          if M<M_crit:
               M_crit=M
               Stress_crit=max_stress
          ratio = math. ceil(ratio*1.1)

     print([row_number,M_crit,Stress_crit])
     print("")
     print("")






