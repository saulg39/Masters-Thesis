import numpy as np
import math
from channel import channel
from stress_from_strain import stress_from_strain
import time

def mid_axes(A, n, Eel, spr, x, y, t_list):
  tot_area = 0
  tot_yarea = 0
  max_y = y[0]
  min_y = y[0]
  for i in range(len(x)-1):
    area = math.sqrt((x[i+1]-x[i])**2 + (y[i+1]-y[i])**2) * t_list[i]
    tot_area += area
    tot_yarea += area * 0.5 * (y[i] + y[i+1])
    max_y = max(max_y, y[i + 1])
    min_y = min(min_y, y[i + 1])
  B = tot_yarea / tot_area
  print(B, max_y, min_y)
  # get reference moment
  half_height = (max_y - min_y) / 2
  print(A * half_height)
  print(stress_from_strain(A * half_height, n, Eel, spr))
  ref_moment = stress_from_strain(A * half_height, n, Eel, spr) * tot_area * half_height
  print(ref_moment)





x, y, t_list = channel(45,23,125,1.98,4)
mid_axes(0.0001, 7.5,208000,328, x, y, t_list)
