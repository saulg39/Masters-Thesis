import numpy as np
import math
def channel(b,c,d,t,r):

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

     y_flange = np.array([dflat/2 + r] * (nflange))
     x_flange = np.linspace(bflat+r,r,nflange)

     x_web = np.array([0] * (math.ceil(nweb/2)))
     y_web = np.linspace(dflat/2,0,math.ceil(nweb/2))

     x = np.concatenate((x_lip, x_flange, x_web), axis = None)
     y = np.concatenate((y_lip, y_flange, y_web), axis = None)

     x_m = x[-2::-1]
     y_m = -1*y[-2::-1]

     x = np.concatenate((x, x_m), axis = None)
     y = np.concatenate((y, y_m), axis = None)

     connections = []

     for i in range(len(x)-1):
          connections.append([i, i+1, t])

     return x, y, connections