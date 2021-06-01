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


for n in range(100,126):
     Material, Forming_process, Bending_type, Actual_Length, shape, Material_flat, Material_corner = get_data(n)
     if shape[0] == "I Beam":
          
          x, y, connections = I_beam(b = shape[1], d = shape[2], t_web = shape[4], t_flange = shape[5], r = shape[3])
          J = (2*shape[1]*shape[5]**3 + shape[2]*shape[4]**3)/3
          
     elif shape[0] == "channel":

          x, y, connections = channel(b = shape[1], c = shape[6], d = shape[2], t = shape[4], r = shape[3])

     elif shape[0] == "RHS":

          x, y, connections = RHS(b = shape[1], d = shape[2], t = shape[4], r = shape[3])
          J = 1.9*((shape[1]*shape[2])**2 * shape[4])/(shape[1]+shape[2])

     else:
          x, y, connections = I_beam(b = shape[1], d = shape[2], t_web = shape[4], t_flange = shape[5], r = shape[3])
     Iyy = 0
     a=0
     for con in connections:
          i,j,t,f = con

          area = math.sqrt((x[j]-x[i])**2 + (y[j]-y[i])**2) * t
          a+= area
          Iyy +=((x[i] +x[j])/2) **2 * area
     moment = 0

     for con in connections:
          i,j,t,f = con
          if Material_corner[0] == "Y" and f == True:
               if y[i] + y[j] > 0:
                    s1 = Material_corner[1][1]
               if y[i] + y[j] < 0:
                    s1 = -Material_corner[1][1]

          else:
               if y[i] + y[j] > 0:
                    s1 = Material_flat[1]
               if y[i] + y[j] < 0:
                    s1 = -Material_flat[1]

          

          area = math.sqrt((x[j]-x[i])**2 + (y[j]-y[i])**2) * t
          moment += (((y[i] +y[j]) * s1)/2) * area
     


     G = Material_flat[0]/(2*(1+Material_flat[3]))

     M_cr=math.sqrt(((math.pi/(Actual_Length)) * math.sqrt(Material_flat[0]*Iyy*G*J))**2+((math.pi/(Actual_Length))**2*Material_flat[0]*Iyy*shape[2]/2)**2)
     print(M_cr/1000000)









     