import numpy as np
import math
from channel import channel
import time
import matplotlib.pyplot as plt

def stress_from_strain(strain, n, Eel, spr, s_ult = False):
      s_ult = False
      sign = np.sign(strain)
      strain = abs(strain)
      if strain == 0:
            return [0, Eel]
      else:
            strain_det = spr / Eel + 0.002
            if strain < strain_det:
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
                  Et_est = (spr * Eel) / (spr + 0.002 * n * Eel * (stress_est / spr) ** (n - 1))

            else:
                  Ey = Eel / (1 + 0.002 * n * Eel / spr)
                  if s_ult==False:
                        s_ult = spr * (1 - 0.0375 * (n - 5)) / (0.2 + 185 * spr / Eel)
                  m = 1 + 3.5 * spr / s_ult
                  strain_ult = 1 - spr / s_ult
                  stress_est = spr + (strain - strain_det) * Ey
                  strain_est = strain_det + (stress_est - spr) / Ey + strain_ult * ((stress_est - spr) / (s_ult - spr)) ** m
                  stop = 0
                  while  abs((strain - strain_est) / strain) > 0.0001 and stop < 10000:
                        Et_est = ((s_ult - spr) * Ey) / ((s_ult - spr) + strain_ult * m * Ey * ((stress_est - spr) / (s_ult - spr)) ** (m - 1))
                        stress_est = (strain - strain_est) * Et_est + stress_est
                        strain_est = strain_det + (stress_est - spr) / Ey + strain_ult * ((stress_est - spr) / (s_ult - spr)) ** m
                        stop+=1
                  Et_est = ((s_ult - spr) * Ey) / ((s_ult - spr) + strain_ult * m * Ey * ((stress_est - spr) / (s_ult - spr)) ** (m - 1))

            if stop < 10000:
                  return [sign * stress_est, Et_est]
            else:
                  return "Time out Error"

"""
tic = time.perf_counter()
x = []
y1 = []
y2 = []
for i in range(10000):
      strain = i * 0.000002
      value1 = stress_from_strain(strain,7.5,208000,328,True)
      value2 = stress_from_strain(strain,7.5,208000,328,False)
      if value1=="Time out Error" or value2=="Time out Error":
            print("Error stress is ",strain)
      else:
            y1.append(value1)
            y2.append(value2)
            x.append(strain)
toc = time.perf_counter()
print(f"Done in {toc - tic:0.4f} seconds")
plt.plot(x,y1)
plt.plot(x,y2)
plt.show()"""