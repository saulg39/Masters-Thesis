import numpy as np
import math
from channel import channel
import time

def stress_from_strain(strain, n, Eel, spr):
  sign = np.sign(strain)
  strain = abs(strain)
  if strain == 0:
    return 0
  stress_est = strain * Eel
  strain_est = stress_est / Eel + 0.002 * (stress_est / spr) ** n 
  stop = 0
  while  abs((strain - strain_est) / strain) > 0.0001 and stop < 10000:
    
    stress_est2 = strain / (1 / Eel + (0.002 / spr ** n) * stress_est ** (n-1))
    stress_est3 = (stress_est ** n * stress_est2) ** (1 / (n + 1))
    strain_est3 = stress_est3 / Eel + 0.002 * (stress_est3 / spr) ** n
    Et_est = (spr * Eel) / (spr + 0.002 * n * Eel * (stress_est3 / spr) ** (n - 1))
    stress_est = (strain - strain_est3) * Et_est + stress_est3
    strain_est = stress_est / Eel + 0.002 * (stress_est / spr) ** n
    stop+=1
  if stop < 10000:
    return sign * stress_est
  else:
    return "Time out Error"
"""tic = time.perf_counter()
for n in range(1000):
  for t in range(1000):
    if stress_from_strain(t**1.5*0.00001 + n**2*0.0000001 ,7.5,208000,328)=="Time out Error":
      print("Error stress is t*0.000001")

toc = time.perf_counter()
print(f"Done in {toc - tic:0.4f} seconds")"""