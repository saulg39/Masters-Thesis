import numpy as np
import math
def RHS(b,d,t,r):

     # returns the x and y coordinates of a lipped channel
     # 4 elements per corner
     # 10 elements for the web
     # 4 elements for flange and lip

     nweb = 10
     nflange = 4
     ncorner = 4
     

     bflat = b-t-2*r
     dflat = d-t-2*r

     y_flange = np.array([dflat/2 + r] * (math.ceil(nflange/2+1)))
     x_flange = np.linspace(0,bflat/2,math.ceil(nflange/2+1))

     x_web = np.array([bflat/2 + r] * (math.ceil(nweb/2+1)))
     y_web = np.linspace(dflat/2,0,math.ceil(nweb/2+1))

     x = np.concatenate((x_flange, x_web), axis = None)
     y = np.concatenate((y_flange, y_web), axis = None)

     x_m1 = x[-2::-1]
     y_m1 = -1*y[-2::-1]
     x_m2 = -1*x[1:]
     y_m2 = -1*y[1:]
     x_m3 = -1*x[-2:0:-1]
     y_m3 = y[-2:0:-1]

     x = np.concatenate((x, x_m1, x_m2, x_m3), axis = None)
     y = np.concatenate((y, y_m1, y_m2, y_m3), axis = None)

     connections = []

     for i in range(len(x)-1):
          connections.append([i, i+1, t])
     connections.append([i+1, 0, t])
     return x, y, connections
