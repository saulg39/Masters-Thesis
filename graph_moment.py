import numpy as np
import math
import matplotlib.pyplot as plt
from channel import channel
from I_beam import I_beam
from RHS import RHS
from stress_from_strain import stress_from_strain


def moment_graph(shape, b, d, r, t_web, t_flange, c, Eel, spr, n, v, k, last):
    if shape == "I Beam":
          
          x, y, connections = I_beam(b, d, t_web, t_flange, r)
          
    elif shape == "channel":

          x, y, connections = channel(b, c, d, t_web, r)

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
    
    #max_A = 20 * spr / max_abs_y / Eel
    A_list = np.linspace(0,last*0.0000226,500)
    for A in A_list:
        stress_list = []
        for i in range(len(x)):
            stress_list.append(stress_from_strain(A *(y[i]-B), n, Eel, spr))
        moment = 0
        for con in connections:
            area = math.sqrt((x[con[1]]-x[con[0]])**2 + (y[con[1]]-y[con[0]])**2) * con[2]
            s1 = stress_list[con[0]][0]
            s2 = stress_list[con[1]][0]
            moment += (((2 * y[con[0]] +y[con[1]]) * s1 + (2 * y[con[1]] +y[con[0]]) * s2)/6) * area
        moment_list.append(moment/10**6/185)
        Apl_list.append(A/0.0000226)
        
    return(moment_list, Apl_list)




"""moment, A = moment_graph(shape = "I Beam", b = 139, d = 199.27, r = 6, t_flange = 10.26, t_web = 7.99, c = 12.7, Eel = 203000, spr = 571, n = 7.5, v = 0.3, k = -0.46)
plt.figure(figsize=(11, 8))
plt.plot(A, moment, 'r-')
plt.ylabel('Moment / KNm')
plt.grid(True,'both')
plt.legend()
plt.show()"""