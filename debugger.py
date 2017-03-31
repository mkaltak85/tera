from mpl_toolkits.mplot3d import Axes3D
from matplotlib.collections import PolyCollection
from matplotlib.colors import colorConverter
import matplotlib.pyplot as plt
import numpy as np


# -----------------------------------------------------------------------------
# 
# Debugger for hull class
#
# -----------------------------------------------------------------------------
def QhullDebugger( sx ) :
   from scipy.spatial import ConvexHull
   from poscar import AttachArray
   from numpy import zeros, append, arange, random
   from sys import exit 

   #
   # set number of configuration points 
   #
   npoints=len(sx.db.keys())
   
   #
   # consistency check
   #
   if npoints <= 0: 
      sys.exit('QhullDebugger has no points to work with')

   sx.hull.npoints=npoints

   #
   # coordinates in configuration space 
   #
   coord = zeros((sx.hull.npoints,3))
   dist  = zeros(sx.hull.npoints)
   #
   # attach object to class
   #
   sx.hull.points=AttachArray(coord)

   #
   # loop over all entries and store coordinates 
   #
   i=0
   keys=[]
   for key in sx.db.keys():
      sx.hull.points[i,0]=sx.db[key]["cart-coord"][0]
      sx.hull.points[i,1]=sx.db[key]["cart-coord"][1]
      sx.hull.points[i,2]=sx.db[key]["formationE"]
      #
      # add index to db
      #
      sx.db[key]["hull_idx"]=i
      keys=append(keys,key)
      i+=1

   #
   # attach keys 
   #
   sx.hull.keys=AttachArray(keys)

   #
   # compute convex hull 
   #
   hull = ConvexHull(sx.hull.points,qhull_options='Qv')

   # 
   # plot figure 
   # 
   fig = plt.figure()
   ax = fig.gca(projection='3d')
   
   def cc(arg):
       return colorConverter.to_rgba(arg, alpha=0.6)
   
   x = sx.hull.points[:,0]
   y = sx.hull.points[:,1]
   z = sx.hull.points[:,2]
   ax.scatter(x,y,z)

   plt.show()

 #  exit('STOP')
 #  verts = []
 #  zs = [0.0, 1.0, 2.0, 3.0]
 #  for z in zs:
 #      ys = random.rand(len(xs))
 #      ys[0], ys[-1] = 0, 0
 #      verts.append(list(zip(xs, ys)))
 #  
 #  poly = PolyCollection(verts, facecolors=[cc('r'), cc('g'), cc('b'),
 #                                           cc('y')])
 #  poly.set_alpha(0.7)
 #  ax.add_collection3d(poly, zs=zs, zdir='y')
 #  
 #  ax.set_xlabel('X')
 #  ax.set_xlim3d(0, 10)
 #  ax.set_ylabel('Y')
 #  ax.set_ylim3d(-1, 4)
 #  ax.set_zlabel('Z')
 #  ax.set_zlim3d(0, 1)
 #  
 #  plt.show()
