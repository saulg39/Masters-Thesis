import numpy as np
import math
def channel_nl(b,d,t,r):

     # returns the x and y coordinates of a lipped channel
     # 4 elements per corner
     # 10 elements for the web
     # 4 elements for flange and lip

     nweb = 10
     nflange = 4
     ncorner = 4

     bflat = b-t/2-r
     dflat = d-t-2*r
  



     y_flange = np.array([dflat/2 + r] * (nflange+1))
     x_flange = np.linspace(bflat+r,r,nflange+1)

     y_corner2 = np.array([])
     x_corner2 = np.array([])

     for n in range(1,ncorner):
          x_corner2 = np.append(x_corner2, r-r*math.sin(n*math.pi/(2*ncorner)))
          y_corner2 = np.append(y_corner2, dflat/2+r*math.cos(n*math.pi/(2*ncorner)))

     x_web = np.array([0] * (math.ceil(nweb/2+1)))
     y_web = np.linspace(dflat/2,0,math.ceil(nweb/2+1))

     x = np.concatenate((x_flange, x_corner2, x_web), axis = None)
     y = np.concatenate((y_flange, y_corner2, y_web), axis = None)

     x_m = x[-2::-1]
     y_m = -1*y[-2::-1]

     x = np.concatenate((x, x_m), axis = None)
     y = np.concatenate((y, y_m), axis = None)

     connections = []


     for j in range(nflange):
          i =  j
          connections.append([i, i+1, t,False])
     for j in range(ncorner):
          i = + nflange + j
          connections.append([i, i+1, t,True])
     for j in range(nweb):
          i = ncorner + nflange + j
          connections.append([i, i+1, t,False])
     for j in range(ncorner):
          i = nflange + nweb + j
          connections.append([i, i+1, t,True])
     for j in range(nflange):
          i = 2 * ncorner + nflange + nweb + j
          connections.append([i, i+1, t,False])



     return x, y, connections

