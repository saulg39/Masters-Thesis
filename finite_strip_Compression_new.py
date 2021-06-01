import math
import numpy as np
from scipy import linalg
from channel import channel
from channel_nl import channel_nl
from channel_imp import channel_imp
from I_beam import I_beam
from RHS import RHS
from stress_from_strain import stress_from_strain
import time
import matplotlib.pyplot as plt
import pandas as pd 
# cd c:/Users/saulg/Documents/IIB_Project
#C:\Users\saulg\Documents\year 4\IIB_Project\code
#function [scr] = finitestrip_shape(L)

def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx

def finitestrip_shape(L, shape, Material_flat, Material_corner, sg,delta=0):
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
     #  sg     a guess at the critical buckling stress
     #  Et     tangent modulus
     #
     ##########################################################################

     
     ## Run channel.m

     if shape[0] == "I Beam":
          
          x, y, connections = I_beam(b = shape[1], d = shape[2], t_web = shape[4], t_flange = shape[5], r = shape[3])
          
     elif shape[0] == "channel":

          x, y, connections = channel(b = shape[1], c = shape[6], d = shape[2], t = shape[4], r = shape[3]) 
     elif shape[0] == "channel_imp":

          x, y, connections = channel_imp(b = shape[1], c = shape[6], d = shape[2], t = shape[4], r = shape[3],delta = delta)

     elif shape[0] == "RHS":

          x, y, connections = RHS(b = shape[1], d = shape[2], t = shape[4], r = shape[3])
     elif shape[0] == "C":

          x, y, connections = channel_nl(b = shape[1], d = shape[2], t = shape[4], r = shape[3])
     elif shape[0] == "L":

          x, y, connections = RHS(b = shape[1], d = shape[2], t = shape[4], r = shape[3])

     else:
          x, y, connections = I_beam(b = shape[1], d = shape[2], t_web = shape[4], t_flange = shape[5], r = shape[3])
     
     ## Initialising variables
     scr = 0
     lamda = 2

     sg=sg/2
     s_upper = False
     s_lower = False
     stop = False
     count = 0


     while not stop:
          count +=1
          if s_upper == False or s_lower == False:
               sg = sg * lamda
          else:
               ratio = - (math.log(s_upper[1]) / math.log(s_lower[1]))
               ratio = ratio ** 0.88
               if ratio > 1:
                    sg = (s_lower[0] * s_upper[0] ** (1/ratio))**(1/((1/ratio) + 1))
               else:
                    sg = (s_upper[0] * s_lower[0] ** ratio)**(1/(ratio + 1))

          K = np.zeros((4*len(x),4*len(x)))

          G = np.zeros((4*len(x),4*len(x)))


          if Material_corner[0] == "helol":
               strain = sg / Material_flat[0] + 0.002 * (sg / Material_flat[0]) ** Material_flat[2] 
               stress_corner = stress_from_strain(strain, Material_corner[1])[0]



          ## Forming global K and G matrices loop

          for con in connections:


               ## assigning stresses
               i,j,t,f = con

               if Material_corner[0] == "hello" and f == True:
                    Eel = Material_corner[1][0]
                    spr = Material_corner[1][1]
                    v = Material_corner[1][3]
                    k = Material_corner[1][4]
                    n = Material_corner[1][2]
                    stress = stress_corner

               else:
                    Eel = Material_flat[0]
                    spr = Material_flat[1]
                    v = Material_flat[3]
                    k = Material_flat[4]
                    n = Material_flat[2]
                    stress = sg
               ## Calculating Et and phi
               Et = (spr * Eel) / (spr + 0.002 * n * Eel * (stress / spr) ** (n - 1))

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

               H = (stress * bel * t * math.pi ** 2 / (2 * L)) * np.array(            
                    [[1/1, 1/2, 1/3, 1/4, 0, 0, 0, 0],
                    [1/2, 1/3, 1/4, 1/5, 0, 0, 0, 0],
                    [1/3, 1/4, 1/5, 1/6, 0, 0, 0, 0],
                    [1/4, 1/5, 1/6, 1/7, 0, 0, 0, 0],
                    [0, 0, 0, 0, 1/1, 1/2, 0, 0],
                    [0, 0, 0, 0, 1/2, 1/3, 0, 0],
                    [0, 0, 0, 0, 0, 0, 1/1, 1/2],
                    [0, 0, 0, 0, 0, 0, 1/2, 1/3]])

               Gel = np.transpose(R).dot(np.transpose(C).dot(H.dot(C.dot(R))))

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
               if s_upper == False:
                    s_upper = [sg,lamda]
               elif s_upper[0] > sg:
                     s_upper = [sg,lamda]
          elif lamda>1:
               if s_lower == False:
                    s_lower = [sg,lamda]
               elif s_lower[0] < sg:
                     s_lower = [sg,lamda]
          if s_upper == False or s_lower == False:
               if abs(lamda-1)<0.001 or count > 60:
                    stop = True
          elif count > 60 or (s_upper[0]-s_lower[0])/s_lower[0] < 0.001:
               stop = True

     Force = 0
     Area = 0
     for con in connections:
          area = math.sqrt((x[con[1]]-x[con[0]])**2 + (y[con[1]]-y[con[0]])**2) * con[2]
          if Material_corner[0] == "hello" and con[3] == True:
               stress = stress_corner

          else:
               stress = sg

          Force += stress * area
          Area += area
          strain = sg / Material_flat[0] + 0.002 * (sg / Material_flat[0]) ** Material_flat[2] 

     return Force/1000, sg, strain, Et
           
         
      
       
     

"""

n = 300
max_L = 10000
min_L = 20
r = (max_L/min_L) ** (1/(n-1))
Stress = []
Stress_1 = []
Length = []
Length_factor= True
S_2 = 300

tic = time.perf_counter()
for i in range(n):
     L = min_L * r**i
     Length.append(L)
     if i == 0:
          sg = 300
     elif i == 1:
          sg = s1
          s2 = s1
     else:
          sg = abs(2 * s1 - s2)
          s3 = s2
          s2 = s1
     
     S = finitestrip_shape(L, shape = "channel", b = 45, d = 125, r = 4, t_flange = 1.98, t_web = 1.98, c = 23, Eel = 208000, spr = 328, n = 7.5, v = 0.3, k = -0.46, sg = sg)
     s1 = S 
     if Length_factor:
          Stress_1.append(S)
          for n in range(2,round(L/Length[0])+1):
               S_2 = Stress[find_nearest(Length, L/n)]
               if S > S_2:
                    S = S_2
     Stress.append(S)
toc = time.perf_counter()
print(f"Done in {toc - tic:0.4f} seconds")


if Length_factor:
     plt.semilogx(Length, Stress_1, linewidth = 0.4, color = "black", label='Singe Half Wavelength')
     plt.semilogx(Length, Stress, linewidth = 1.2, color = "black", label='Multiple Half Wavelengths')
else:
     plt.semilogx(Length, Stress, linewidth = 0.8, color = "black", label='Inelastic buckling')
     

#plt.semilogx([1000], [151], 'rx',  label='Buckling')
#plt.semilogx([1000], [169], 'bx',  label='Failure')
plt.xlabel('L / mm')
plt.ylabel('Stress / MPa')
plt.grid(True,'both')
plt.legend()
plt.show()"""




























"""lambdas, modes = linalg.eig(K,G)

     mode = (modes[:,index])
     zz = []
     xx = []
     yy = []
     for i in range(len(x)):
          zz.append(mode[4*i-3]*100)
          xx.append(mode[4*i-2]*100)
          yy.append(mode[4*i-1]*100)
           

     coordinates = [np.transpose(xx), np.transpose(yy), np.transpose(zz)]

     pd.DataFrame(coordinates).to_csv("coordinates.cvs", header=None, index=None)"""