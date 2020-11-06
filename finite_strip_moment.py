import math
import numpy as np
from scipy import linalg
from channel import channel
from stress_from_strain import stress_from_strain
import time
import matplotlib.pyplot as plt
import pandas as pd 
#C:\Users\saulg\Documents\year 4\IIB_Project\code
#function [scr] = finitestrip_shape(L)
def finitestrip_shape(L):
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

     Eel = 208000
     spr = 328
     n = 7.5
     v = 0.3
     k = -0.46
     sg = 150
     t = 1.98

     ## Run channel.m

     x, y, t_list = channel(45,23,125,1.98,0.01)
     

     ## Initialising variables
     lamda = 2

     A=0.00001
     B=0
     
     stop = 0

     while abs(lamda-1)>0.001 and stop < 200:
          stop +=1
          A = A * lamda ** ((1/stop)**0.2)
          K = np.zeros((4*len(x),4*len(x)))

          G = np.zeros((4*len(x),4*len(x)))

          ## Calculating Et and phi
          stress_list = []
          for i in range(len(x)):
               stress_list.append(stress_from_strain(A * (y[i]-B), n, Eel, spr))

          ## Forming global K and G matrices loop

          for i in range(len(x)-1):
               ## assigning stresses
               s1 = stress_list[i][0]

               s2 = stress_list[i+1][0]

               ## Calculating Et and phi

               Et1 = (spr * Eel) / (spr + 0.002 * n * Eel * (abs(s1) / spr) ** (n - 1))

               Et2 = (spr * Eel) / (spr + 0.002 * n * Eel * (abs(s2) / spr) ** (n - 1))

               Et = math.sqrt(stress_list[i][1]*stress_list[i+1][1]) 

               phi = Eel / ((1 + v * k) * Eel - (v + k) * v * Et)

               ## Calculating width of element, bel, and thickness, t

               bel = math.sqrt((x[i+1]-x[i])**2 + (y[i+1]-y[i])**2)
               
               t = t_list[i]

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

               hyp = math.sqrt((x[i+1]-x[i])**2 + (y[i+1]-y[i])**2)

               R = np.array([[0, -(y[i+1]-y[i])/hyp, (x[i+1]-x[i])/hyp, 0, 0, 0, 0, 0],
                    [0, 0, 0, 1, 0, 0, 0, 0],
                    [0, 0, 0, 0, 0, -(y[i+1]-y[i])/hyp, (x[i+1]-x[i])/hyp, 0],
                    [0, 0, 0, 0, 0, 0, 0, 1],
                    [0, (x[i+1]-x[i])/hyp, (y[i+1]-y[i])/hyp, 0, 0, 0, 0, 0],
                    [1, 0, 0, 0, 0, 0, 0, 0],
                    [0, 0, 0, 0, 0, (x[i+1]-x[i])/hyp, (y[i+1]-y[i])/hyp, 0],
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

               K[(4*i):(4*i+8),(4*i):(4*i+8)] = K[(4*i):(4*i+8),(4*i):(4*i+8)] + Kel


               ## Assembling G matrix

               G[(4*i):(4*i+8),(4*i):(4*i+8)] = G[(4*i):(4*i+8),(4*i):(4*i+8)] + Gel

               
           
           
          ## Calculating buckling stress

          w = linalg.eigvals(K,G)
          w = w.real
          index = np.where(w > 0, w, np.inf).argmin()
          lamda = w[index]
          #print(A,lamda )
          
          
          #[scr index] = min(eig(K,G)) 
     if stop > 200:
          return "Fail"

     moment = 0
     for i in range(len(x)-1):
          area = math.sqrt((x[i+1]-x[i])**2 + (y[i+1]-y[i])**2) * t_list[i]
          s1 = stress_list[i][0]
          s2 = stress_list[i+1][0]
          moment += (((2 * y[i] +y[i + 1]) * s1 + (2 * y[i + 1] +y[i]) * s2)/6) * area

     return moment/10**6
  

#print(finitestrip_shape(3000))

x = []
y = []

for i in range(300):
     number = 10 * 1.016**i
     y.append(finitestrip_shape(number))
     x.append(number)

plt.semilogx(x,y)
plt.xlabel('L/mm')
plt.ylabel('Moment KNm')
plt.show()