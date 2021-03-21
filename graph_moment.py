import numpy as np
import math
import matplotlib.pyplot as plt
from channel import channel
from I_beam import I_beam
from RHS import RHS
from stress_from_strain import stress_from_strain
from plastic_moment import plastic_moment

def moment_graph(shape, Material_flat, Material_corner, last, curve_pl=1, s_pl=1, s_ult = False):


    if shape[0] == "I Beam":
          
        x, y, connections = I_beam(b = shape[1], d = shape[2], t_web = shape[4], t_flange = shape[5], r = shape[3])
          
    elif shape[0] == "channel":

        x, y, connections = channel(b = shape[1], c = shape[6], d = shape[2], t = shape[4], r = shape[3])

    elif shape[0] == "RHS":

        x, y, connections = RHS(b = shape[1], d = shape[2], t = shape[4], r = shape[3])

    else:
        x, y, connections = I_beam(b = shape[1], d = shape[2], t_web = shape[4], t_flange = shape[5], r = shape[3])


    max_abs_y = y[0]
    B=0
    moment_list = []
    Apl_list = []
    for num in y:
        if abs(num-B) > max_abs_y:
            max_abs_y = abs(num-B)

    #max_A = 20 * spr / max_abs_y / Eel
    A_list = np.linspace(0,last*curve_pl,150)

    for A in A_list:
        stress_list = []
        for i in range(len(x)):
            stress_list.append(stress_from_strain(A * (y[i]-B), Material_flat))
        moment = 0
        for con in connections:
            i,j,t,f = con
            if Material_corner[0] == "Y" and f == True:
                stress_list_1 = stress_from_strain(A * (y[i]-B), Material_corner[1])
                s1 = stress_from_strain(A * (y[i]-B), Material_corner[1])

            else:
                stress_list_1 = stress_list[i]
                stress_list_2 = stress_list[j]
            
            s1 = stress_list_1[0]
            s2 = stress_list_2[0]

            area = math.sqrt((x[j]-x[i])**2 + (y[j]-y[i])**2) * t
            moment += (((2 * y[i] +y[j]) * s1 + (2 * y[i] +y[j]) * s2)/6) * area
        moment_list.append(moment/10**6/s_pl)
        Apl_list.append(A/curve_pl)
    return(moment_list, Apl_list)




"""moment, A = moment_graph(shape = "I Beam", b = 139, d = 199.27, r = 6, t_flange = 10.26, t_web = 7.99, c = 12.7, Eel = 203000, spr = 571, n = 7.5, v = 0.3, k = -0.46)
plt.figure(figsize=(11, 8))
plt.plot(A, moment, 'r-')
plt.ylabel('Moment / KNm')
plt.grid(True,'both')
plt.legend()
plt.show()"""