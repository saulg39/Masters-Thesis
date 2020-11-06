import numpy as np
import math
def I_beam(b,d,t,r):

     # returns the x and y coordinates of a lipped channel
     # 4 elements per corner
     # 10 elements for the web
     # 4 elements for flange and lip

     nweb = 8
     nflange = 4
     ncorner = 4

     bflat = (b-t-2*r)/2
     dflat = d-2*t-2*r


     y_flange = np.array([dflat/2 + r + t/2] * (nflange+1))
     t_flange = np.array([t] * (nflange+1))
     x_flange = np.linspace(bflat+r+ t/2,r+ t/2,nflange+1)

     x_web = np.array([0] * (math.ceil(nweb+1)))
     t_web = np.array([t] * (math.ceil(nweb+1)))
     y_web = np.linspace(dflat/2,-dflat/2,math.ceil(nweb+1))

     x = np.concatenate((x_flange, -1 * x_flange, x_flange, -1 * x_flange, x_web), axis = None)
     y = np.concatenate((y_flange, y_flange, -1 * y_flange, -1 * y_flange, y_web), axis = None)
     t = np.concatenate((t_flange, t_flange, t_flange, t_flange, t_web), axis = None)

     connections = []
     for i in range(4):
          for j in range(nflange):
               connections.append([i*(nflange+1)+j,i*(nflange+1)+j+1])
     for j in range(nweb):
          connections.append([4*(nflange+1)+j,4*(nflange+1)+j+1])
     connections.append([nflange,nflange*2+1])
     connections.append([nflange*3+2,nflange*4+3])
     connections.append([nflange,nflange*4+4])
     connections.append([nflange*4+4,nflange*2+1])
     connections.append([nflange*3+2,nflange*4+nweb+4])
     connections.append([nflange*4+nweb+4,nflange*4+3])

     return x, y, t, connections

