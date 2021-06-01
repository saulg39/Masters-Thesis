import numpy as np
import math
def I_beam(b, d, t_web, t_flange, r):

     # returns the x and y coordinates of a lipped channel
     # 4 elements per corner t
     # 10 elements for the web
     # 4 elements for flange and lip

     nweb = 10
     nflange = 4
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
               connections.append([i*(nflange+1)+j,i*(nflange+1)+j+1, t_flange, True,131.8])
     for j in range(nweb):
          connections.append([4*(nflange+1)+j,4*(nflange+1)+j+1, t_web, False,115.6])
     connections.append([nflange,nflange*2+1, t_flange*0.5, True])
     connections.append([nflange*3+2,nflange*4+3, t_flange*0.5, True])
     connections.append([nflange,nflange*4+4, t_web*0.5, True])
     connections.append([nflange*4+4,nflange*2+1, t_web*0.5, True])
     connections.append([nflange*3+2,nflange*4+nweb+4, t_web*0.5, True])
     connections.append([nflange*4+nweb+4,nflange*4+3, t_web*0.5, True])

     return x, y, connections

print(I_beam(140,220,6,6,0))
"""(array([ 70.  ,  53.25,  38.6 ,  19.3,   3.  , -70.  , -53.25, -38.6 ,
       -19.3,  -3.  ,  70.  ,  53.25,  38.6 ,  19.3,   3.  , -70.  ,
       -53.25, -38.6 , -19.3,  -3.  ,   0.  ,   0.  ,   0.  ,   0.  ,
         0.  ,   0.  ,   0.  ,   0.  ,   0.  ,   0.  ,   0.  ]), array([ 107. ,  107. ,  107. ,  107. ,  107. ,  107. ,  107. ,  107. ,
        107. ,  107. , -107. , -107. , -107. , -107. , -107. , -107. ,
       -107. , -107. , -107. , -107. ,  104. ,   83.2,   62.4,   41.6,
         20.8,    0. ,  -20.8,  -41.6,  -62.4,  -83.2, -104. ]), [[0, 1, 6, True], [1, 2, 6, True], [2, 3, 6, True], [3, 4, 6, True], [5, 6, 6, True], [6, 7, 6, True], [7, 8, 6, True], [8, 9, 6, True], [10, 11, 6, True], [11, 12, 6, True], [12, 13, 6, True], [13, 14, 6, True], [15, 16, 6, True], [16, 17, 6, True], [17, 18, 6, True], [18, 19, 6, True], [20, 21, 6, False], [21, 22, 6, False], [22, 23, 6, False], [23, 24, 6, False], [24, 25, 6, False], [25, 26, 6, False], [26, 27, 6, False], [27, 28, 6, False], [28, 29, 6, False], [29, 30, 6, False], [4, 9, 3.0, True], [14, 19, 3.0, True], [4, 20, 3.0, True], [20, 9, 3.0, True], [14, 30, 3.0, True], [30, 19, 3.0, True]])
"""