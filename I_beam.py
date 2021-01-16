import numpy as np
import math
def I_beam(b, d, t_web, t_flange, r):

     # returns the x and y coordinates of a lipped channel
     # 4 elements per corner t
     # 10 elements for the web
     # 4 elements for flange and lip

     nweb = 8
     nflange = 3
     ncorner = 4
     
     bflat = (b-t_web-2*r)/2
     dflat = d-2*t_flange-2*r


     y_flange = np.array([dflat/2 + r + t_flange/2] * (nflange+1))
     x_flange = np.linspace(bflat+r+ t_web/2,r+ t_web/2,nflange+1)

     x_web = np.array([0] * (math.ceil(nweb+1)))
     y_web = np.linspace(dflat/2,-dflat/2,math.ceil(nweb+1))

     x = np.concatenate((x_flange, -1 * x_flange, x_flange, -1 * x_flange, x_web), axis = None)
     y = np.concatenate((y_flange, y_flange, -1 * y_flange, -1 * y_flange, y_web), axis = None)

     connections = []
     for i in range(4):
          for j in range(nflange):
               connections.append([i*(nflange+1)+j,i*(nflange+1)+j+1, t_flange])
     for j in range(nweb):
          connections.append([4*(nflange+1)+j,4*(nflange+1)+j+1, t_web])
     connections.append([nflange,nflange*2+1, t_flange])
     connections.append([nflange*3+2,nflange*4+3, t_flange])
     connections.append([nflange,nflange*4+4, t_web])
     connections.append([nflange*4+4,nflange*2+1, t_web])
     connections.append([nflange*3+2,nflange*4+nweb+4, t_web])
     connections.append([nflange*4+nweb+4,nflange*4+3, t_web])

     return x, y, connections

