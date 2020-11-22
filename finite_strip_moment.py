import math
import numpy as np
from scipy import linalg
from channel import channel
from I_beam import I_beam
from stress_from_strain import stress_from_strain
import time
import matplotlib.pyplot as plt
import pandas as pd 
#C:\Users\saulg\Documents\year 4\IIB_Project\code
#function [scr] = finitestrip_shape(L)
def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx

def finitestrip_shape(L, shape, b, d, r, t_web, t_flange, c, Eel, spr, n, v, k):
     ##########################################################################
     #
     #  This program finds the critical buckling stress of a channel
     #  section using the finite strip method.
     #
     #  It requires an initial guess of the critical buckling stress
     #  to calculate the tangent modulus and then calculates a buckling stress
     #  from that value.
     #
     #  If the guess and the result are too far apart then the program changes
     #  the value of the guess to the average of the previous two and calculates
     #  a new critical stress.
     #
     #  The program repeats this iteration until the difference between 
     #  the guess and the resulting critical stress is less than 0.001 MPa.
     #
     #  Notation:
     #
     #  scr    critical buckling stress
     #  b      width of section
     #  c      length of lip
     #  d      depth of section
     #  t      thickness of section
     #  r      radius of corners
     #  Eel    initial modulus
     #  spr    0.2# proof stress
     #  n      defines roundness of stress strain curve
     #  k      defines ratio between shear stress and strain increments
     #  v      poisson's ratio
     #  L      half-wavelength of locally buckled section
     #  Et     tangent modulus
     #
     ##########################################################################

     ## Run channel.m

     if shape == "I Beam":
          
          x, y, connections = I_beam(b, d, t_web, t_flange, r)
          
     elif shape == "channel":

          x, y, connections = channel(b, c, d, t_web, r)

     else:
          x, y, connections = I_beam(b, d, t_web, t_flange, r)


     ## Initialising variables
     lamda = 2

     A=0.00001
     A_upper = False
     A_lower = False
     B=0
     stop = False
     count = 0

     while not stop:
          count +=1
          if A_upper == False or A_lower == False:
               A = A * lamda
          else:
               ratio = - (math.log(A_upper[1]) / math.log(A_lower[1]))
               ratio = ratio ** 0.8
               if ratio > 1:
                    A = (A_lower[0] * A_upper[0] ** (1/ratio))**(1/((1/ratio) + 1))
               else:
                    A = (A_upper[0] * A_lower[0] ** ratio)**(1/(ratio + 1))

          K = np.zeros((4*len(x),4*len(x)))

          G = np.zeros((4*len(x),4*len(x)))

          ## Calculating Et and phi
          stress_list = []
          for i in range(len(x)):
               stress_list.append(stress_from_strain(A * (y[i]-B), n, Eel, spr))

          ## Forming global K and G matrices loop

          for con in connections:
               ## assigning stresses
               i,j,t = con
               s1 = stress_list[i][0]

               s2 = stress_list[j][0]

               ## Calculating Et and phi

               Et = math.sqrt(stress_list[i][1]*stress_list[j][1]) 

               phi = Eel / ((1 + v * k) * Eel - (v + k) * v * Et)

               ## Calculating width of element, bel, and thickness, t

               bel = math.sqrt((x[j]-x[i])**2 + (y[j]-y[i])**2)
               

               ## Calculating K11

               S11 = phi * (t ** 3 / 12) * Et * (math.pi ** 4 * bel / L ** 3)
               S12 = phi * (t ** 3 / 12) * Et * v * (math.pi ** 2 / (bel * L))
               S21 = phi * (t ** 3 / 12) * ((k + v) * Et - k * Eel) * (math.pi ** 2 / (bel * L))
               S22 = phi * (t ** 3 / 12) * Eel * (L / bel ** 3)
               S33 = (t ** 3 / 12) * (Eel * Et / ((1 + k + 2 * v) * Et + (1 - k) * Eel)) * (math.pi ** 2 / (bel * L))

               K11 = np.array([[S11 / 2, S11 / 4, S11 / 6 - S12, S11 / 8 - 3 * S12 / 2],
               [S11 / 4, S11 / 6 + 2 * S33, S11 / 8 - S12 / 2 + 2 * S33, S11 / 10 - S12 + 2 * S33],
               [S11 / 6 - S21, S11 / 8 - S21 / 2 + 2 * S33, S11 / 10 - S21 / 3 - S12 / 3 + 2 * S22 + 8 * S33 / 3, S11 / 12 - S21 / 4 - 3 * S12 / 4 + 3 * S22 + 3 * S33],
               [S11 / 8 - 3 * S21 / 2, S11 / 10 - S21 + 2 * S33, S11 / 12 - 3 * S21 / 4 - S12 / 4 + 3 * S22 + 3 * S33, S11 / 14 - 3 * S21 / 5 - 3 * S12 / 5 + 6 * S22 + 18 * S33 / 5]])


               ## Calculating K22

               F11 = phi * Et * (math.pi ** 2 * bel * t / L)
               F12 = phi * v * Et * math.pi * t
               F21 = phi * ((v + k) * Et - k * Eel) * math.pi * t
               F22 = phi * Eel * (L * t / bel)
               F33 = (Eel * Et / ((1 - k) * Eel + (1 + k + 2 * v) * Et)) * (math.pi ** 2 * bel * t / L)

               K22 = np.array([[F33 / 2, F33 / 4, 0, F33 * L / (2 * bel * math.pi)],
               [F33 / 4, F22 / 2 + F33 / 6, -F21 / 2, F33 * L / (4 * bel * math.pi) - F21 / 4],
               [0, -F12 / 2, F11 / 2, F11 / 4],
               [F33 * L / (2 * bel * math.pi), F33 * L / (4 * bel * math.pi) - F12 / 4, F11 / 4, F11 / 6 + F33 * L ** 2 / (2 * bel ** 2 * math.pi ** 2)]])


               ## Calculating K element matrix

               C = np.array([[1, 0, 0, 0, 0, 0, 0, 0], 
                    [0, bel, 0, 0, 0, 0, 0, 0], 
                    [-3, -2*bel, 3, -bel, 0, 0, 0, 0], 
                    [2, bel, -2, bel, 0, 0, 0, 0], 
                    [0, 0, 0, 0, 1, 0, 0, 0], 
                    [0, 0, 0, 0, -1, 0, 1, 0], 
                    [0, 0, 0, 0, 0, 1, 0, 0], 
                    [0, 0, 0, 0, 0, -1, 0, 1]])

               hyp = math.sqrt((x[j]-x[i])**2 + (y[j]-y[i])**2)

               R = np.array([[0, -(y[j]-y[i])/hyp, (x[j]-x[i])/hyp, 0, 0, 0, 0, 0],
                    [0, 0, 0, 1, 0, 0, 0, 0],
                    [0, 0, 0, 0, 0, -(y[j]-y[i])/hyp, (x[j]-x[i])/hyp, 0],
                    [0, 0, 0, 0, 0, 0, 0, 1],
                    [0, (x[j]-x[i])/hyp, (y[j]-y[i])/hyp, 0, 0, 0, 0, 0],
                    [1, 0, 0, 0, 0, 0, 0, 0],
                    [0, 0, 0, 0, 0, (x[j]-x[i])/hyp, (y[j]-y[i])/hyp, 0],
                    [0, 0, 0, 0, 1, 0, 0, 0]])

               K1122 = np.zeros((8,8))

               K1122[:4,:4] = K11

               K1122[4:,4:] = K22

               Kel = np.transpose(R).dot(np.transpose(C).dot(K1122.dot(C.dot(R))))


               ## Calculating G element matrix

               H = (bel * t * math.pi ** 2 / (2 * L)) * np.array(            
                    [[s1/2 + s2/2, s1/6 + s2/3, s1/12 + s2/4, s1/20 + s2/5, 0, 0, 0, 0],
                    [s1/6 + s2/3, s1/12 + s2/4, s1/20 + s2/5, s1/30 +s2/6, 0, 0, 0, 0],
                    [s1/12 + s2/4, s1/20 + s2/5, s1/30 +s2/6, s1/42 + s2/7, 0, 0, 0, 0],
                    [s1/20 + s2/5, s1/30 +s2/6, s1/42 + s2/7, s1/56 +s2/8, 0, 0, 0, 0],
                    [0, 0, 0, 0, s1/2 + s2/2, s1/6 + s2/3, 0, 0],
                    [0, 0, 0, 0, s1/6 + s2/3, s1/12 + s2/4, 0, 0],
                    [0, 0, 0, 0, 0, 0, s1/2 + s2/2, s1/6 + s2/3],
                    [0, 0, 0, 0, 0, 0, s1/6 + s2/3, s1/12 + s2/4]])

               Gel = np.transpose(R).dot(np.transpose(C).dot(H.dot(C.dot(R))))


               ## Assembling K matrix

               K[(4*i):(4*i+4),(4*i):(4*i+4)] = K[(4*i):(4*i+4),(4*i):(4*i+4)] + Kel[0:4,0:4]

               K[(4*j):(4*j+4),(4*i):(4*i+4)] = K[(4*j):(4*j+4),(4*i):(4*i+4)] + Kel[4:8,0:4]

               K[(4*i):(4*i+4),(4*j):(4*j+4)] = K[(4*i):(4*i+4),(4*j):(4*j+4)] + Kel[0:4,4:8]

               K[(4*j):(4*j+4),(4*j):(4*j+4)] = K[(4*j):(4*j+4),(4*j):(4*j+4)] + Kel[4:8,4:8]


               ## Assembling G matrix

               G[(4*i):(4*i+4),(4*i):(4*i+4)] = G[(4*i):(4*i+4),(4*i):(4*i+4)] + Gel[0:4,0:4]

               G[(4*j):(4*j+4),(4*i):(4*i+4)] = G[(4*j):(4*j+4),(4*i):(4*i+4)] + Gel[4:8,0:4]

               G[(4*i):(4*i+4),(4*j):(4*j+4)] = G[(4*i):(4*i+4),(4*j):(4*j+4)] + Gel[0:4,4:8]

               G[(4*j):(4*j+4),(4*j):(4*j+4)] = G[(4*j):(4*j+4),(4*j):(4*j+4)] + Gel[4:8,4:8]


               
           
           
          ## Calculating buckling stress

          w = linalg.eigvals(K,G)
          w = w.real
          index = np.where(w > 0, w, np.inf).argmin()
          lamda = w[index]

          if lamda<1:
               if A_upper == False:
                    A_upper = [A,lamda]
               elif A_upper[0] > A:
                     A_upper = [A,lamda]
          elif lamda>1:
               if A_lower == False:
                    A_lower = [A,lamda]
               elif A_lower[0] < A:
                     A_lower = [A,lamda]
          if A_upper == False or A_lower == False:
               if abs(lamda-1)<0.001 or count > 60:
                    stop = True
          elif count > 60 or (A_upper[0]-A_lower[0])/A_lower[0] < 0.001:
               stop = True
          
          #[scr index] = min(eig(K,G)) 
     if count > 200:
          return "Fail"

     moment = 0
     for con in connections:
          area = math.sqrt((x[con[1]]-x[con[0]])**2 + (y[con[1]]-y[con[0]])**2) * con[2]
          s1 = stress_list[con[0]][0]
          s2 = stress_list[con[1]][0]
          moment += (((2 * y[con[0]] +y[con[1]]) * s1 + (2 * y[con[1]] +y[con[0]]) * s2)/6) * area
     return moment/10**6
  

n = 130
max_L = 4000
min_L = 20
r = (max_L/min_L) ** (1/(n-1))
Moment = []
Moment_1 = []
Length = []
Length_factor= False

for i in range(n):
     L = min_L * r**i
     Length.append(L)
     M = finitestrip_shape(L, shape = "channel", b = 50.8, d = 152.4, r = 0.01, t_flange = 1, t_web = 1, c = 12.7, Eel = 203000, spr = 579, n = 7.6, v = 0.3, k = -0.46)
     if Length_factor:
          Moment_1.append(M)
          for n in range(2,round(L/Length[0])+1):
               M_2 = Moment[find_nearest(Length, L/n)]
               if M > M_2:
                    M = M_2

     Moment.append(M)
     
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