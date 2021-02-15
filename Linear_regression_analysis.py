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
from sklearn.linear_model import LinearRegression
import statsmodels.api as sm
THIS_FOLDER = os.path.dirname(os.path.abspath(__file__))
my_file = os.path.join(THIS_FOLDER, 'Data','data.xls')

book = xlrd.open_workbook(my_file)
# get the first worksheet
sheet = book.sheet_by_index(0)

th=[]
ex=[]

for i in range(sheet.nrows):
     th.append([sheet.cell(i,0).value])
     ex.append(sheet.cell(i,1).value)


x, y = np.array(th), np.array(ex)

x = sm.add_constant(x)

model = sm.OLS(y, x)

results = model.fit()

print(results.summary())
print('coefficient of determination:', results.rsquared)

model = LinearRegression(fit_intercept=True).fit(x, y)

r_sq = model.score(x, y)
print('coefficient of determination:', r_sq)

print('intercept:', model.intercept_)

print('slope:', model.coef_)

fig, ax = plt.subplots(figsize=(6,4))
plt.scatter(th, ex,c='k',marker='x', label='Data points')
plt.plot([0,1000],[0,1000], linewidth = 0.1, color = "black")
plt.plot([0,1000],[ model.intercept_, model.intercept_+1000*model.coef_[1]], linewidth = 0.6,linestyle = '--', color = "black", label='Line of best fit')
plt.ylabel('Experimental Buckling Moment / KNm')
plt.xlabel('Theoretical Buckling Moment / KNm')
ax.set_xlim([0,300])
ax.set_ylim([0,300])         
plt.grid(True,'both')
plt.legend()
plt.show()