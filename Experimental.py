import math
import numpy as np
from scipy import linalg
from channel import channel
from I_beam import I_beam
from stress_from_strain import stress_from_strain
import time
import matplotlib.pyplot as plt
import pandas as pd 
import xlrd
import os
THIS_FOLDER = os.path.dirname(os.path.abspath(__file__))
my_file = os.path.join(THIS_FOLDER, 'Data','beam_column_data.xlsx')

book = xlrd.open_workbook(my_file)



# get the first worksheet
first_sheet = book.sheet_by_index(1)



# read a row slice
row = first_sheet.row_slice(rowx=78,
                            start_colx=10,
                            end_colx=14)

print(row[2].value)























"""
names = ['Original', 'Current']
values = [4,0.02]

plt.figure(figsize=(9, 9))

plt.subplot(111)
plt.bar(names, values)
plt.suptitle('Time taken to calculate the buckling moment for a half-sine length of a cross section ')
plt.show()"""


"""plt.rcdefaults()
fig, ax = plt.subplots()

# Example data
people = ('Original', 'Current')
y_pos = np.arange(len(people))
performance = [4,0.02]
error = np.random.rand(len(people))

ax.bar(y_pos, performance, 0.3)
ax.set_xticks(y_pos)
ax.set_xticklabels(people)
#ax.invert_xaxis()  # labels read top-to-bottom
ax.set_ylabel('Time  /  s')
ax.set_title('Time taken to calculate the buckling moment\nfor a half-sine length of a cross section ')

plt.show()
"""












"""
n = 300
max_L = 4000
min_L = 20
r = (max_L/min_L) ** (1/(n-1))
Moment = []
Moment_1 = []
Length = []
Length_factor= False

for i in range(n):
     L = min_L * r**i
     print(L)




def value(B):
	H = 0.1 *2
	V = 0.4 
	su = 0.080
	M = 4 * H
	e = M/V
	A = B*(B-2*e)
	sc = 1.2
	Vult = A*sc*(2+math.pi)*su
	Hult = A*su
	return (H/Hult + (2*V/Vult - 1)**2)

for i in range(100):
	n = 0.01*i + 4.1
	print(value(n),n)"""










"""
def value(V):
	
	return 0.5*(math.pi + 2) *(1 + math.sqrt(1-0.25*V)) - V
for i in range(100):
	n = i*0.0001 +3.48
	print(value(n),n)


def value(H):
	D = 30
	V = 50 + 15
	su = 0.050
	M = 60 * H - 150
	e = M/V
	th = 2 * math.acos(2*e/D)
	A = D**2 * th/4 - 2*e*math.sqrt((D/2)**2-e**2)
	sc = 1 + 0.18 *(math.sqrt((D-2*e)/(D+2*e)))
	Vult = A*sc*(2+math.pi)*su
	Hult = A*su
	return (H/Hult + (2*V/Vult - 1)**2)

for i in range(100):
	n = 25 + 0.01*i
	print(value(n),n)


def value(H):
	return (0.233**2 + ((60*H/1420785)*(1-0.3*H/35343))**2 +(abs(H/35343))**3)-1

for i in range(100):
	print(value(23760+0.1*i),23760+0.1*i)"""