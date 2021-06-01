import numpy as np
import math
import matplotlib.pyplot as plt
from channel import channel
from I_beam import I_beam
from RHS import RHS
from stress_from_strain import stress_from_strain
import os
import xlrd


def plastic_moment(shape, b, d, r, t_web, t_flange, c, Eel, spr, n, v, k):
    if shape == "I Beam":
          
          x, y, connections = I_beam(b, d, t_web, t_flange, r)
          
    elif shape == "channel":

          x, y, connections = channel(b, c, d, t_web, r)

    elif shape == "channel_imp":

          x, y, connections = channel_imp(b, c, d, t_web, r)

    elif shape == "RHS":

          x, y, connections = RHS(b, d, t_web, r)

    else:
          x, y, connections = I_beam(b, d, t_web, t_flange, r)
    max_abs_y = y[0]
    B=0
    moment_list = []
    Apl_list = []
    for num in y:
        if abs(num-B) > max_abs_y:
            max_abs_y = abs(num-B)
    
    A =spr / max_abs_y / Eel
    stress_list = []
    for i in range(len(x)):
        stress_list.append([A *(y[i]-B)*Eel,1])
    moment = 0
    for con in connections:
        area = math.sqrt((x[con[1]]-x[con[0]])**2 + (y[con[1]]-y[con[0]])**2) * con[2]
        s1 = stress_list[con[0]][0]
        s2 = stress_list[con[1]][0]
        moment += (((2 * y[con[0]] +y[con[1]]) * s1 + (2 * y[con[1]] +y[con[0]]) * s2)/6) * area
    print(moment/10**6,A)
    return (moment/10**6, A)
        





"""THIS_FOLDER = os.path.dirname(os.path.abspath(__file__))
my_file = os.path.join(THIS_FOLDER, 'Data','beam_column_data.xlsx')

book = xlrd.open_workbook(my_file)
# get the first worksheet
first_sheet = book.sheet_by_index(1)
row = first_sheet.row_slice(rowx=84,
                            start_colx=9,
                            end_colx=31)

b, d, r, t_flange, t_web, c, spr, s_ult, n, Eel, s_pl, curve_pl = [row[1].value, row[0].value, row[4].value, row[3].value, row[2].value, row[3].value, row[15].value, row[16].value, row[17].value, row[18].value, row[20].value, row[21].value]

plastic_moment(shape = "I Beam", b = b, d = d, r = r, t_flange = t_flange, t_web = t_web, c = c, Eel = Eel, spr = spr, n = n, v = 0.3, k = -0.46)"""
