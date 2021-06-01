import numpy as np
import math
def channel_imp(b,c,d,t,r,delta):

     # returns the x and y coordinates of a lipped channel
     # 4 elements per corner
     # 10 elements for the web
     # 4 elements for flange and lip


     nweb = 10
     nflange = 4
     nlip = 4
     ncorner = 4

     bflat = b-t-2*r
     dflat = d-t-2*r
     cflat = c-t/2-r


     x_lip = np.array([bflat + 2*r] * (nlip+1))
     y_lip = np.linspace(0,cflat,nlip+1)
     y_lip = y_lip + 0.5*dflat-cflat

     y_corner1 = np.array([])
     x_corner1 = np.array([])

     for n in range(1,ncorner):
          x_corner1 = np.append(x_corner1, bflat+r+r*math.cos(n*math.pi/(2*ncorner)))
          y_corner1 = np.append(y_corner1, dflat/2+r*math.sin(n*math.pi/(2*ncorner)))

     y_flange = np.array([dflat/2 + r] * (nflange+1))
     x_flange = np.linspace(bflat+r,r,nflange+1)

     y_corner2 = np.array([])
     x_corner2 = np.array([])

     for n in range(1,ncorner):
          x_corner2 = np.append(x_corner2, r-r*math.sin(n*math.pi/(2*ncorner)))
          y_corner2 = np.append(y_corner2, dflat/2+r*math.cos(n*math.pi/(2*ncorner)))

     
     y_web = np.linspace(dflat/2,0,math.ceil(nweb/2+1))
     x_web = np.array([dflat*delta*math.cos(math.pi*num/dflat)  for num in y_web])
     
     x = np.concatenate((x_lip, x_corner1, x_flange, x_corner2, x_web), axis = None)
     y = np.concatenate((y_lip, y_corner1, y_flange, y_corner2, y_web), axis = None)

     x_m = x[-2::-1]
     y_m = -1*y[-2::-1]

     x = np.concatenate((x, x_m), axis = None)
     y = np.concatenate((y, y_m), axis = None)

     connections = []

     for i in range(nlip):
          connections.append([i, i+1, t,False])
     for j in range(ncorner):
          i = nlip + j
          connections.append([i, i+1, t,True])
     for j in range(nflange):
          i = nlip + ncorner + j
          connections.append([i, i+1, t,False])
     for j in range(ncorner):
          i = nlip +ncorner + nflange + j
          connections.append([i, i+1, t,True])
     for j in range(nweb):
          i = nlip + 2 * ncorner + nflange + j
          connections.append([i, i+1, t,False])
     for j in range(ncorner):
          i = nlip + 2 * ncorner + nflange + nweb + j
          connections.append([i, i+1, t,True])
     for j in range(nflange):
          i = nlip + 3 * ncorner + nflange + nweb + j
          connections.append([i, i+1, t,False])
     for j in range(ncorner):
          i = nlip + 3 * ncorner + 2 * nflange + nweb + j
          connections.append([i, i+1, t,True])
     for j in range(nflange):
          i = nlip + 4 * ncorner + 2 * nflange + nweb + j
          connections.append([i, i+1, t,False])


     return x, y, connections

