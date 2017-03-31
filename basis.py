import numpy as np
from numpy import dot as np_dot
from numpy import sqrt as np_sqrt
from numpy import abs as np_abs
from numpy import sin as np_sin
from numpy import cos as np_cos
from numpy import arccos as np_arccos
from numpy import savetxt as np_savetxt
from numpy import array as np_array
from numpy import ndarray as np_ndarray
from numpy import linalg as np_linalg
from numpy import append as np_append
import sys 

# define some important constants
pi=3.14159265358979323846
sqrt3=1.7320508075688772

global directoryOUT
directoryOUT='./output'

global eps, epsh, epss
eps=1.E-6
epss=1.E-8
epsh=1.E-4  #

#coordinates of endpoints
Cx=0; Cy=0
Ax=0.5; Ay=0.5*sqrt3
Bx=1; By=0

#phase name dictionary  
System={"phase": ["Li", "Ag", "Mn_{8}O_{16}"] } 

#maximum concentrations
MaxX=8.0

#some HEX color codes 
global HEXred
global HEXblue
global HEXgreen
HEXred='#ae0000'
HEXblue='#0033cc'
HEXgreen='#009933'

#sample rate: 1/NConcentY0Points
NConcentY0Points=1000
# treshhold for neighborhood of line 
epsH=1./1000.

global ConcentY0

# Returns the directory the current script (or interpreter) is running in
def GetAppDir():
    from os import path as OsPath
     
    path = OsPath.realpath(sys.argv[0])
    if OsPath.isdir(path):
        return path
    else:
        return OsPath.dirname(path)

#------------------------------------------------------------------------------
# prints content to file 
#------------------------------------------------------------------------------
def PrintTo(handle,string):
   #write string to file if its open
   if not(handle.closed):
      handle.write(string+'\n')
   #write after opening
   else:
      sys.exit('Error: PrintTo reports undefined handle'+str(handle))

#------------------------------------------------------------------------------
# creates name A_x B_x C_x with xvec = [x,y,z]
#------------------------------------------------------------------------------
def CreateSystemName(xvec):
   string=System["phase"][0]+'_{'+xvec[0]+'}'
   string+=System["phase"][1]+'_{'+xvec[1]+'}'
   string+=System["phase"][2]
   return string

#------------------------------------------------------------------------------
# computes length between two points in 2D 
#------------------------------------------------------------------------------
def LenBetPts2D( x1,y1,x0,y0 ):
    v=np_sqrt((x1-x0)**2+(y1-y0)**2)
    return v

#------------------------------------------------------------------------------
# computes normvector of plane defined by three points 
#------------------------------------------------------------------------------
def norm_vector( simplex, allpoints ):
    "vector orthogonal to plane defined by three points i,j,k in points"
    "not normalized !"
    xi=allpoints[simplex[0]]
    xj=allpoints[simplex[1]]
    xk=allpoints[simplex[2]]
    n1=(xj[1]-xi[1])*(xk[2]-xi[2])-(xj[2]-xi[2])*(xk[1]-xi[1])
    n2=(xj[2]-xi[2])*(xk[0]-xi[0])-(xj[0]-xi[0])*(xk[2]-xi[2])
    n3=(xj[0]-xi[0])*(xk[1]-xi[1])-(xj[1]-xi[1])*(xk[0]-xi[0])
    norm=np_sqrt(n1*n1+n2*n2+n3*n3)
    n1=n1/norm
    n2=n2/norm
    n3=n3/norm
    if ( n3 < 0 ) :
       array=np_array([-n1,-n2,-n3])
    else:
       array=np_array([n1,n2,n3])
    return array

#------------------------------------------------------------------------------
# computes unit vector between to points
#------------------------------------------------------------------------------
def unit_vector(A,B):
   v=B[0]-A[0]
   v=np_append(v,B[1]-A[1])
   
   norm=np_sqrt(np_dot(v,v))

   print "vector and norm",norm,v

   if norm>eps:
      v=v/norm
   else:
      sys.exit("Unit_vector reports coincident points")
   
   return v 

#------------------------------------------------------------------------------
# computes the length of vector A,B
#------------------------------------------------------------------------------
def length_of_vector(A,B):
   v=B[0]-A[0]
   v=np_append(v,B[1]-A[1])
   
   norm=np_sqrt(np_dot(v,v))
   
   return norm


#------------------------------------------------------------------------------
# Get all barycentric coordinates from a point p=[x0,y0]
#------------------------------------------------------------------------------
def GetBaryCentricCoordinates(p):
   x0=p[0]
   y0=p[1]
   x1=Cx;y1=Cy
   x2=Ax;y2=Ay
   x3=Bx;y3=By

   b=[0,0,0]
   x=[0,0,0]
   # baricentric coordinates [b] are obtained from the solution 
   # 
   # x0 = b1 x1 + b2 x2 + b3 x3
   # y0 = b1 y1 + b2 y2 + b3 y3
   #  1 = b1 + b2 + b3
   #  
   # inserting b3= 1-b1-b2 into Eq. 1 and two and solving the system for b1,b2
   # yields:
    
   # determinant needs to be non-zero
   d=(x1-x3)*(y2-y3)+(y1-y3)*(x3-x2)
   if np_abs(d)<eps: 
      #in this case p is an endpoint find out which one
      if np_abs(x0-Cx)<eps and np_abs(y0-Cy) : 
         return [epss,epss,1]
      elif np_abs(x0-Ax)<eps and np_abs(y0-Ay) : 
         return [epss,1,epss]
      elif np_abs(x0-Bx)<eps and np_abs(y0-By) : 
         return [1,epss,epss]
      else:
         sys.exit('Error in GetConcentrations'+str(p))
   else:
      b[0]=((y2-y3)*(x0-x3)+(x3-x2)*(y0-y3)) /d   # phase C in %
      b[1]=((y3-y1)*(x0-x3)+(x1-x3)*(y0-y3)) /d   # pahse B in %
      b[2]=1-b[0]-b[1]                            # phase A in %
 
#      for i in [0,1,2]:
#         b[i]=round( b[i] , 12)
#         if np_abs(b[i])< 1.E-14:
#            b[i]=0

      # normalize stoichiometries with C_z==1, this works only if b[0]>0
      if np_abs(b[1]-1)<eps:
         return [epss,1,epss]
      elif np_abs(b[2]-1)<eps:
         return [1,epss,epss]
      elif np_abs(b[0]-1)<eps:
         return [epss,epss,1]
      else:
         return [b[2],b[1],b[0]]

#------------------------------------------------------------------------------
# Get all concentrations from carthesian coordinates for a point inside 
# the triangle 
#------------------------------------------------------------------------------
def GetConcentrations(p):
   x0=p[0]
   y0=p[1]
   x1=Cx;y1=Cy
   x2=Ax;y2=Ay
   x3=Bx;y3=By

   b=[0,0,0]
   x=[0,0,0]
   # baricentric coordinates [b] are obtained from the solution 
   # 
   # x0 = b1 x1 + b2 x2 + b3 x3
   # y0 = b1 y1 + b2 y2 + b3 y3
   #  1 = b1 + b2 + b3
   #  
   # inserting b3= 1-b1-b2 into Eq. 1 and two and solving the system for b1,b2
   # yields:
    
   # determinant needs to be non-zero
   d=(x1-x3)*(y2-y3)+(y1-y3)*(x3-x2)
   if np_abs(d)<eps: 
      #in this case p is an endpoint find out which one
      if np_abs(x0-Cx)<eps and np_abs(y0-Cy) : 
         x[2]=1
         return x
      elif np_abs(x0-Ax)<eps and np_abs(y0-Ay) : 
         x[1]=1
         return x
      elif np_abs(x0-Bx)<eps and np_abs(y0-By) : 
         x[0]=1
         return x
      else:
         sys.exit('Error in GetConcentrations'+str(p))
   else:
      b[0]=((y2-y3)*(x0-x3)+(x3-x2)*(y0-y3)) /d   # phase C in %
      b[1]=((y3-y1)*(x0-x3)+(x1-x3)*(y0-y3)) /d   # pahse B in %
      b[2]=1-b[0]-b[1]                            # phase A in %
   
      for i in [0,1,2]:
        if np_abs(b[i]) < eps:
           b[i]=0

      # normalize stoichiometries with C_z==1, this works only if b[0]>0
      if np_abs(b[1]-1)<eps:
         x[1]=1
         return x 
      elif np_abs(b[2]-1)<eps:
         x[0]=1
         return x 
      else:
         if 1-b[2]-b[1]>eps:
            x[0]=round(b[2]/(1-b[2]-b[1]),12)
            x[1]=round(b[1]/(1-b[2]-b[1]),12)
            if b[0]>eps: 
               x[2]=1
            else: 
               x[2]=0
         else:
            if b[1]>b[2]:
               x[0]=round(b[2]/(1-b[2]),12)
               x[1]=round(b[1]/(1-b[1]),12)
               x[0]=x[0]*x[1]
            else:
               x[0]=round(b[2]/(1-b[2]),12)
               x[1]=round(b[1]/(1-b[1]),12)
               x[1]=x[0]*x[1]
         return x 

#------------------------------------------------------------------------------
# Computes the concentration of Phase A from projection
#------------------------------------------------------------------------------
def ConcentA(Xc,Yc):
    " determine Li concentration via projection "
    " works as follows: "
    " determine intersection point of line from 100% Ag point (Ax,Ay) to (Xc,Yc)"
    if np_abs(Xc-Ax)<eps and np_abs(Yc-Ay)<eps:
       t=0
    elif np_abs(Xc-Bx)<eps and np_abs(Yc-By)<eps:
       x=1
       return x
    else:
       t=Ax+Ay*(Xc-Ax)/(Ay-Yc)
       if np_abs(t)<epsh:
          t=0
       if t>1: t=1
       if t<0: t=0
    if np_abs(t-1)<eps:
       x=0
    else:
       x=t/(1-t)
    return round(x,12)

#------------------------------------------------------------------------------
# Computes the concentration of Phase B from projection
#------------------------------------------------------------------------------
def ConcentB(Xc,Yc):
    " same as for Li but now for Ag" 
    #this is the Li point
    if np_abs(Xc-Bx)<eps and np_abs(Yc-By)<eps:
       t=0
    elif np_abs(Xc-Ax)<eps and np_abs(Yc-Ay)<eps:
       x=1
       return x
    #this is the Mn point
    elif np_abs(Xc)<eps and np_abs(Yc)<eps:
       t=0
    else:
       #general point 
       t=(By*(Xc-Bx)-Bx*(Yc-By))/(Ay*(Xc-Bx)-Ax*(Yc-By))
       if t>1: t=1
       if t<0: t=0
    if np_abs(t-1)<eps:
       x=0
    else:
       x=t/(1-t)
    return round(x,12)

#------------------------------------------------------------------------------
# Computes the concentration of Phase C from projection
#------------------------------------------------------------------------------
def ConcentC(Xc,Yc):
    #check if X,Y is on the line B-A, if so the concentration is 0 else 1
    k=(Ay-By)/(Ax-Bx)
    if ( np_abs(Yc-By-k*(Xc-Bx))>eps ):
       x=1.0
    else:
       x=0.0
    return x 

#------------------------------------------------------------------------------
# Computes binary reaction weights for 
#
#  A_x0 B_y0 C_z0 -> a[0] A_x1 B_y1 C_z1 + a[1] A_x2 B_y2 C_z2 +
#
#------------------------------------------------------------------------------
def GetBinReactionWeights(x0,x1,x2,typ):

   a = np_ndarray(shape=(2,2), dtype=float, order='F')

   #type defines on which concentrations one is interested
   if typ == 'BC': 
      #in this case x[0] is always zero 
      b=np_array( [ x0[1],x0[2] ])
      # reactants 
      a[0,0]=x1[1]
      a[1,0]=x1[2]
      a[0,1]=x2[1] 
      a[1,1]=x2[2] 
  
   elif typ == 'AC':
      #in this case x[0] is always zero 
      b=np_array( [ x0[0],x0[2] ])
      # reactants 
      a[0,0]=x1[0]
      a[1,0]=x1[2]
      a[0,1]=x2[0] 
      a[1,1]=x2[2] 

   elif typ == 'AB':
      #in this case x[0] is always zero 
      b=np_array( [ x0[0],x0[1] ])
      # reactants 
      a[0,0]=x1[0]
      a[1,0]=x1[1]
      a[0,1]=x2[0] 
      a[1,1]=x2[1] 

   else: 
      sys.exit('Error in GetBinReactionWeights:'+type+' does not make sense')

 
   lamb=np_linalg.solve(a,b)

   return lamb 

#------------------------------------------------------------------------------
# Computes reaction weights for 
#
#  A_x0 B_y0 C_z0 -> a[0] A_x1 B_y1 C_z1 + a[1] A_x2 B_y2 C_z2 +
#                    a[2] A_x3 B_y3 C_z3
#
# here the l.h.s. is given by (x',y') Carthesian coordinates of the triangle 
# as well as the r.h.s in form of a 3x3 matrix s in the same coordinate system  
# the (x_i,y_i)_i=1..3 are the first and second column of this matrix s
#------------------------------------------------------------------------------
def GetReactionWeights(p,s):
 
   #first convert every point (x'_i, y'_i)_i=0...3 to Carthesian coordinates
   #product 
   b = np_ndarray(3, dtype=float, order='F')
   b=GetConcentrations(p) 

   # reactants 
   a = np_ndarray(shape=(3,3), dtype=float, order='F')
   xvec=GetConcentrations([s[0,0],s[0,1]])
   a[0,0]=xvec[0]
   a[1,0]=xvec[1]
   a[2,0]=xvec[2]
   xvec=GetConcentrations([s[1,0],s[1,1]])
   a[0,1]=xvec[0] 
   a[1,1]=xvec[1] 
   a[2,1]=xvec[2] 
   xvec=GetConcentrations([s[2,0],s[2,1]])
   a[0,2]=xvec[0] 
   a[1,2]=xvec[1] 
   a[2,2]=xvec[2] 

   lamb=np_linalg.solve(a,b)

   return lamb 

#------------------------------------------------------------------------------
# Rotates cut through hull back to x-axis 
#------------------------------------------------------------------------------
def RotateHullCut(x,y):
    #endpoint   (1,0)
    #angle is 
    alpha=pi+np_arccos( (x[0]-1 )/np_sqrt( (x[0]-1)**2+y[0]**2)  )
    #rotate back
    x0=np_cos(alpha)*(x[0]-1)+np_sin(alpha)*y[0]+1
    y0=-np_sin(alpha)*(x[0]-1)+np_cos(alpha)*y[0]
    for i in range(0,len(x)):
       xp=np_cos(alpha)*(x[i]-1)+np_sin(alpha)*y[i]
       yp=-np_sin(alpha)*(x[i]-1)+np_cos(alpha)*y[i]
       x[i]=(xp+1-x0)*(2-np_sqrt(3)/2)
       y[i]=yp+0-y0

#------------------------------------------------------------------------------
# gives the ith point of the path with constant B-concentration defined by 
# ConcentY0
#------------------------------------------------------------------------------
def ConstY0Path(i,n,ConcentY0):
    x1=Bx
    y1=By
    x0=ConcentY0/(ConcentY0+1)*np_cos(pi/3)
    y0=ConcentY0/(ConcentY0+1)*np_sin(pi/3)
    #vector between points (x0,y0)->(x1,y1)
    vx=x1-x0
    vy=y1-y0
    norm=np_sqrt(vx*vx+vy*vy)
    vx=vx/norm
    vy=vy/norm
    scale=norm/float(n)
    pxi=x0+i*scale*vx
    pyi=y0+i*scale*vy
    array=np_array([pxi,pyi])
    return array 

#------------------------------------------------------------------------------
# Computes Centre of mass for a triangle in 3D defined by three points 
#------------------------------------------------------------------------------
def CoM(triangle, allpoints ):
    n1=allpoints[triangle[0],0]+allpoints[triangle[1],0]+\
       allpoints[triangle[2],0]
    n2=allpoints[triangle[0],1]+allpoints[triangle[1],1]+\
       allpoints[triangle[2],1]
    n3=allpoints[triangle[0],2]+allpoints[triangle[1],2]+\
       allpoints[triangle[2],2]
    n1=n1/3.
    n2=n2/3.
    n3=n3/3.
    com=np_array([n1,n2,n3])
    return com

#------------------------------------------------------------------------------
#checks if point (x[i],y[i],z[i]) is not (0,0,0), i.e. not a endpoint   
#------------------------------------------------------------------------------
def IsNot0(i,x,y,z):
    if np_abs(x[i])<epsh or np_abs(y[i])<epsh or np_abs(z[i])<epsh:
       return 0
    else:
       return 1

#------------------------------------------------------------------------------
#checks if point (x[i],y[i]) is on boundary line:
#------------------------------------------------------------------------------
def IsOnEdge(i,j,xp,yp,zp):
    " Checks if point in on boundary "
    " j=1: on line C-B, A_x=0"
    " j=2: on line C-A, B_y=0 "
    " j=3: on line A-B, C_z=0 "
    if ( np_abs(1-j)<eps ): 
       if np_abs(xp[i])<eps :
          return 1
       else:
          return 0
    elif ( np_abs(2-j)<eps ):
       if np_abs(yp[i])<eps :
          return 1
       else:
          return 0
    elif ( np_abs(3-j)<eps ):
       if np_abs(zp[i])<eps :
          return 1
       else:
          return 0
    else:
        sys.exit("Error in IsOnEdge: crazy j value"+str(j))

#------------------------------------------------------------------------------
# Voltage between two points i , j in the list xx,yy,xx having energy e
#------------------------------------------------------------------------------
def Voltage(i,j,xx,yy,zz,e,EndPointE ):
    "gives voltage between points"
    V = e[j]-e[i]
    dx = xx[j]-xx[i]
    dy = yy[j]-yy[i]
    dz = zz[j]-zz[i]
    d = np_sqrt(dx*dx+dy*dy+dz*dz)
    if d < 1.E-6:
       V = 0
    else:
       V =np_abs((V - dx*EndPointE[0]-dy*EndPointE[1]-dz*EndPointE[2])/d)
    return np_round(V,1)

#------------------------------------------------------------------------------
# computes distace to simplex j of point i
# currently not used 
#------------------------------------------------------------------------------
def distance( i , j, simplices,allpoints):
    #this is the normvector on the polygon plane
    nvec=norm_vector(simplices[j],allpoints)
    #this is the first edge point of polygon
    xi=allpoints[simplices[j][0]]
    #this is the distance of plane to origin
    d=-nvec[0]*xi[0]-nvec[1]*xi[1]-nvec[2]*xi[2]
    #this is the projected distance (along z-axis) to polygon plane 
    if np_abs(nvec[2]) > eps:
       r=-np_abs((d-nvec[0]*x[i]-nvec[1]*y[i]-nvec[2]*z[i])/nvec[2])
    else:
       r=0
    return r


#------------------------------------------------------------------------------
# computes distance to hull of given 2D coordinates (X,Y)
#------------------------------------------------------------------------------
def distance_to_hull( X , Y, simplices, allpoints, Sign=None):
    "computes distace to simplex of point v=(X,Y,0) "
    "distance to plane is n.v-n.x0 = n.v-d"
    nvec=norm_vector(simplices,allpoints)
    #this is the first edge point of polygon
    x0=allpoints[simplices[0]]
    #this is the distance of plane to origin
    d=nvec[0]*x0[0]+nvec[1]*x0[1]+nvec[2]*x0[2]
    #this is the projected distance (along z-axis) to polygon plane 
    if np_abs(nvec[2]) > eps:
       if Sign :
          r=(d-nvec[0]*X-nvec[1]*Y)/nvec[2]
       else:
          r=-np_abs((d-nvec[0]*X-nvec[1]*Y)/nvec[2])
    else:
       r=0
    return r

#------------------------------------------------------------------------------
# x, y -- x and y coordinates of point
# poly -- a list of tuples [(x, y), (x, y), ...]
# (inout algorithm)
#------------------------------------------------------------------------------
def isPointInPath(x, y, poly):
    num = len(poly)
    i = 0
    j = num - 1
    c = False
    for i in range(num):
            if  ((poly[i][1] > y) != (poly[j][1] > y)) and \
                    (x < (poly[j][0] - poly[i][0]) * (y - poly[i][1]) / \
                    (poly[j][1] - poly[i][1]) + poly[i][0]):
                c = not c
            j = i
    return c

#------------------------------------------------------------------------------
# checks if simplex is boundary
#------------------------------------------------------------------------------
def IsNotBoundary( i , simplex, allpoints ):
    n=norm_vector(simplex[i],allpoints)
    if np_abs(n[2])>epsh and np_abs(np_abs(n[2])-1)>epsh:
       return 1
    else:
       return 0

#------------------------------------------------------------------------------
# Calculates fromation energies 
#
# calculate the Formation energy DG(x',y',z'), where
# x'=x/(x+y+z), y'=y/(x+y+z), z'=z/(x+y+z) 
# DG(x',y',z')=E(x',y',z')-x'*EA-y'*EB-z'*EC
# 
# Here E(x',y',z')=E(x,y,z)/(x+y+z) and 
# E(x,y,z) is the energy of the compound A_x B_y C_z
#
#------------------------------------------------------------------------------
def GetFormEnAndLatt(ename,fname):
   #load formation energies from formen
   #read file 
   x,y,z,E,a,b,c,v =  np_loadtxt(ename, comments='#', unpack=True)
   npoints=len(x)
   
   #
   # obtain Endpoint energies EA,EB,EC
   #
   k=[0,0,0]
   for i in range(0,npoints):
      if abs(x[i]-1.)<eps and abs(y[i])<eps and abs(z[i])<eps:
         EB = E[i]
         k[0]=k[0]+1
      if abs(x[i])<eps and abs(y[i]-1.)<eps and abs(z[i])<eps:
         EC = E[i]
         k[1]=k[1]+1
      if abs(x[i])<eps and abs(y[i])<eps and abs(z[i]-1.)<eps:
         EA = E[i]
         k[2]=k[2]+1
   
   # Print Warnings if ename has more than 1 endpoint energy.
   if k[0]>1: 
      print "Warning: "+ename+" contains more than one energy for phase B!"+\
        " Will use E(B)=",EB
   if k[1]>1: 
      print "Warning: "+ename+" contains more than one energy for phase C!"+\
        " Will use E(C)=",EC
   if k[2]>1: 
      print "Warning: "+ename+" contains more than one energy for phase A!"+\
        " Will use E(A)=",EA
   
   DG = np_ndarray(shape=(npoints,10), dtype=float, order='F')
   Latt = np_ndarray(shape=(npoints,4), dtype=float, order='F')
   xp=[] ; yp=[] ; zp=[]
   for i in range(0,npoints):
      ntot=x[i]+y[i]+z[i]
      e=E[i]-x[i]*EB-y[i]*EC-z[i]*EA
      xp=np_append(xp,x[i]/ntot)
      yp=np_append(yp,y[i]/ntot)
      zp=np_append(zp,z[i]/ntot)
      #also Lattice info
      Latt[i,0]=a[i]
      Latt[i,1]=b[i]
      Latt[i,2]=c[i]
      Latt[i,3]=v[i]
      DG[i,0]=x[i]
      DG[i,1]=y[i]
      DG[i,2]=z[i]
      DG[i,3]=xp[i]
      DG[i,4]=yp[i]
      DG[i,5]=zp[i]
      #also calculate position on triangle from isometric point of view
      #here (0,0) cooincides with the A phase 
      DG[i,6]=xp[i]+yp[i]*np_cos(pi/3.)
      DG[i,7]=yp[i]*np_sin(pi/3.)
      DG[i,8]=e/ntot
      DG[i,9]=E[i]
    
   np_savetxt(fname,DG,header='  x            y            z            x\''+\
      '           y\'           z\'           x_plot       y_plot       '+\
      'z_plot       E',fmt='%12.6f')
   return Latt


#------------------------------------------------------------------------------
# finds boundary points
#------------------------------------------------------------------------------
def FindBoundary(xp,yp,zp,x,y,z,bound1,bound2,bound3,idx1,idx2,idx3):
   npoints = len(xp)
   for i in range(0,npoints):
       if IsNot0(i,xp,yp,zp) == 0:
          d=0.
          if IsOnEdge(i,1,xp,yp,zp) == 1 : 
             # choosing (y,z) pairs as relevant plane
             bound1=np_append(bound1,[[ y[i],z[i] ]], axis=0)
             idx1=np_append(idx1,i)
          if IsOnEdge(i,2,xp,yp,zp) == 1 : 
             # choosing (x,z) pairs as relevant plane
             bound2=np_append(bound2,[[ x[i],z[i] ]], axis=0)
             idx2=np_append(idx2,i)
          if IsOnEdge(i,3,xp,yp,zp) == 1 : 
             # choosing (x,z) pairs as relevant plane
             bound3=np_append(bound3,[[ x[i],z[i] ]], axis=0)
             idx3=np_append(idx3,i)

#------------------------------------------------------------------------------
# Total energy of point i of file ename
#------------------------------------------------------------------------------
def Energy( i, e ):
    "gives the energy for point i"
    return e[i];

#------------------------------------------------------------------------------
# determines total energy based on endpoint energies 
#------------------------------------------------------------------------------
def TotalEnergy(Xc,Yc,Zc,EndPointE):
    xvec=GetConcentrations([Xc,Yc])
    x=xvec[0]
    y=xvec[1]
    z=xvec[2]
    e=Zc*(x+y+z)+x*EndPointE[0]+y*EndPointE[1]+z*EndPointE[2]
    return e

#------------------------------------------------------------------------------
# Calculates first avergage of 2 lattice constants 
#------------------------------------------------------------------------------
def CalculateLattAvgFirst(latt,points):
   from numpy import minimum as np_minimum
   from numpy import maximum as np_maximum
   npoints=len(points)
   nlatt=len(latt)
   ymin=10000
   ymax=-10000
   if npoints == nlatt: 
      #find concentrations ymin < ConcentY0 < ymax
      for i in range(0,npoints):
         xvec=GetConcentrations([points[i,0],points[i,1]])
         x=xvec[0]
         y=xvec[1]
         z=xvec[2]
         if np_abs(z-1)<epsh and abs(x)<epsh: 
            if y>ConcentY0:
               ymin=np_minimum(y,ymin)
               l=i
            if y<ConcentY0:
               ymax=np_maximum(y,ymax)
               k=i
   lattl=latt[l,:]
   lattk=latt[k,:]
   lattl.sort()
   lattk.sort()
   l1=(lattl[0]+lattk[0])/2
   l2=(lattl[1]+lattk[1])/2
   l3=(lattl[2]+lattk[2])/2
   l4=(lattl[3]+lattk[3])/2
   r=[[l1,l2,l3,l4]]
   return r

#------------------------------------------------------------------------------
# Calculates the avergage of 2 lattice constants 
#------------------------------------------------------------------------------
def CalculateLattAvg2Only(latt,i,j):
    l1=latt[i]
    l2=latt[j]
    l1.sort()
    l2.sort()
    r=[]
    r=np_append(r,(l1[0]+l2[0])/2)
    r=np_append(r,(l1[1]+l2[1])/2)
    r=np_append(r,(l1[2]+l2[2])/2)
    r=np_append(r,(l1[3]+l2[3])/2)
    return r

#------------------------------------------------------------------------------
# Consistentcy check for Stoichiometries and reactions
#------------------------------------------------------------------------------
def IsConsistent(a,x0,x1,x2):
   # points are consistent if two non-zero entries of a conicide
   result=False
   if np_abs(a[0])>epsh and np_abs(a[1])>epsh and np_abs(a[0]-a[1])<epsh:
      result=True
   elif np_abs(a[0])>epsh and np_abs(a[2])>epsh and np_abs(a[0]-a[2])<epsh:
      result=True
   elif np_abs(a[1])>epsh and np_abs(a[2])>epsh and np_abs(a[1]-a[2])<epsh:
      result=True
   elif np_abs(a[0])>epsh and np_abs(a[1])<epsh and np_abs(a[2])<epsh:
      result=True
   elif np_abs(a[1])>epsh and np_abs(a[0])<epsh and np_abs(a[2])<epsh:
      result=True
   elif np_abs(a[2])>epsh and np_abs(a[0])<epsh and np_abs(a[1])<epsh:
      result=True
   else: 
      result=False
   return result

#------------------------------------------------------------------------------
# calculates average of lattice constants for a reaction 
#
#  A_x0 B_y0 C_z0 -> lam[0] A_x1 B_y1 C_z1 + lamb[1] A_x2 B_y2 C_z2 +
#                    lamb[2] A_x3 B_y3 C_z3
#
# the reactants on r.h.s. are represented by a simplex s pointing a the entries
# in Latt
#------------------------------------------------------------------------------
def CalculateLattAvg(lamb,s,Latt):
   npoints=len(Latt)
  
   if npoints < s[0] or npoints < s[1] or npoints < s[2] :
      sys.exit('Error in CalculateLattAvg: Index out of bound'+str(s)+\
         ' '+str(npoints))
   else:
      Latt1=Latt[s[0]]
      Latt2=Latt[s[1]]
      Latt3=Latt[s[2]]
      Latt1.sort()  
      Latt2.sort()  
      Latt3.sort()  

#      print lamb
#      print Latt1
#      print Latt2
#      print Latt3
  
      #calculate average 
      LattAvg=[0,0,0,0]
      if np_abs(sum(lamb)-1)>eps: 
         LattAvg[0]=(Latt1[0]*lamb[0]+Latt2[0]*lamb[1]+Latt3[0]*lamb[2])/sum(lamb)
         LattAvg[1]=(Latt1[1]*lamb[0]+Latt2[1]*lamb[1]+Latt3[1]*lamb[2])/sum(lamb)
         LattAvg[2]=(Latt1[2]*lamb[0]+Latt2[2]*lamb[1]+Latt3[2]*lamb[2])/sum(lamb)
         LattAvg[3]=(Latt1[3]*lamb[0]+Latt2[3]*lamb[1]+Latt3[3]*lamb[2])/sum(lamb)
      else:
         LattAvg[0]=Latt1[0]*lamb[0]+Latt2[0]*lamb[1]+Latt3[0]*lamb[2]
         LattAvg[1]=Latt1[1]*lamb[0]+Latt2[1]*lamb[1]+Latt3[1]*lamb[2]
         LattAvg[2]=Latt1[2]*lamb[0]+Latt2[2]*lamb[1]+Latt3[2]*lamb[2]
         LattAvg[3]=Latt1[3]*lamb[0]+Latt2[3]*lamb[1]+Latt3[3]*lamb[2]
       
      print ' Average ' 
      print LattAvg
      print ' '

   return LattAvg

#------------------------------------------------------------------------------
# calculates average of lattice constants for a reaction 
#
#  A_x0 B_y0 C_z0 -> lam[0] A_x1 B_y1 C_z1 + lamb[1] A_x2 B_y2 C_z2 +
#                    lamb[2] A_x3 B_y3 C_z3
#
# neglects contributions from enpoint B 
#------------------------------------------------------------------------------
def CalculateLattAvgChange(lamb,s,Latt):
   npoints=len(Latt)

   # id of endpoint B in ename list 
   I=1
  
   if npoints < s[0] or npoints < s[1] or npoints < s[2] :
      sys.exit('Error in CalculateLattAvg: Index out of bound'+str(s)+\
         ' '+str(npoints))
   else:
      
      Latt1=Latt[s[0]]
      Latt2=Latt[s[1]]
      Latt3=Latt[s[2]]
      Latt1.sort()  
      Latt2.sort()  
      Latt3.sort()  
#      # in case the endpoint B is involved don't take it into account 
#      if s[0]==I or s[1]==I or s[2] == I:
#         print ' ENDPOINT '
#         print s
#         print lamb,sum(lamb)
#         print Latt1
#         print Latt2
#         print Latt3
#      else:
#         print 'Sould be 1 ',lamb,sum(lamb)
        
      #calculate average 
      LattAvg=[0,0,0,0]
      if np_abs(sum(lamb)-1)>eps: 
         # in the case the weight do not add up to 1 average only over those 
         # constants which add up to 1 
         s1=lamb[0]+lamb[1]
         s2=lamb[0]+lamb[2]
         s3=lamb[1]+lamb[2]
           
         if np_abs(s1-1)<eps: 
            LattAvg[0]=Latt1[0]*lamb[0]+Latt2[0]*lamb[1]
            LattAvg[1]=Latt1[1]*lamb[0]+Latt2[1]*lamb[1]
            LattAvg[2]=Latt1[2]*lamb[0]+Latt2[2]*lamb[1]
            LattAvg[3]=Latt1[3]*lamb[0]+Latt2[3]*lamb[1]
#            print ' Avering over 1 and 2 ' 
         elif np_abs(s2-1)<eps: 
            LattAvg[0]=Latt1[0]*lamb[0]+Latt3[0]*lamb[2]
            LattAvg[1]=Latt1[1]*lamb[0]+Latt3[1]*lamb[2]
            LattAvg[2]=Latt1[2]*lamb[0]+Latt3[2]*lamb[2]
            LattAvg[3]=Latt1[3]*lamb[0]+Latt3[3]*lamb[2]
#            print ' Avering over 1 and 3 ' 
         elif np_abs(s3-1)<eps: 
            LattAvg[0]=Latt2[0]*lamb[1]+Latt3[0]*lamb[2]
            LattAvg[1]=Latt2[1]*lamb[1]+Latt3[1]*lamb[2]
            LattAvg[2]=Latt2[2]*lamb[1]+Latt3[2]*lamb[2]
            LattAvg[3]=Latt2[3]*lamb[1]+Latt3[3]*lamb[2]
#            print ' Avering over 2 and 3 ' 
      else:
         LattAvg[0]=Latt1[0]*lamb[0]+Latt2[0]*lamb[1]+Latt3[0]*lamb[2]
         LattAvg[1]=Latt1[1]*lamb[0]+Latt2[1]*lamb[1]+Latt3[1]*lamb[2]
         LattAvg[2]=Latt1[2]*lamb[0]+Latt2[2]*lamb[1]+Latt3[2]*lamb[2]
         LattAvg[3]=Latt1[3]*lamb[0]+Latt2[3]*lamb[1]+Latt3[3]*lamb[2]
       
#      print ' Average ' 
#      print LattAvg
#      print ' '

   return LattAvg

#------------------------------------------------------------------------------
# calculates average of lattice constants for a reaction 
#
#  A_x0 B_y0 C_z0 -> lam[0] A_x1 B_y1 C_z1 + lamb[1] A_x2 B_y2 C_z2 +
#                    lamb[2] A_x3 B_y3 C_z3
#
# neglects contributions from enpoint B 
#------------------------------------------------------------------------------
def LossOf(I,lamb,s):

   # in case the endpoint B is involved don't take it into account 
   loss=0
   if s[0]==I or s[1]==I or s[2] == I:
#      print ' ENDPOINT '
#      print s
#      print lamb,sum(lamb)

      # one end point is involved in reaction only if sum of weight is not 1
      if np_abs(sum(lamb)-1)>eps: 
         # constants which add up to 1 
         s1=lamb[0]+lamb[1]
         s2=lamb[0]+lamb[2]
         s3=lamb[1]+lamb[2]
           
         if np_abs(s1-1)<eps: 
#            print ' Loss of '+str(I)+':',lamb[2]
            loss=lamb[2]
         elif np_abs(s2-1)<eps: 
#            print ' Loss of '+str(I)+':',lamb[1]
            loss=lamb[1]
         elif np_abs(s3-1)<eps: 
#            print ' Loss of '+str(I)+':',lamb[0]
            loss=lamb[0]
      else:
         sys.exit('Error in LossOf')
       
   return loss
 
   
#------------------------------------------------------------------------------
# calcualte Voltages for arbitrary constant B concentrations 
# asuming that phase A is the intercalating material 
# this is done by a cut through the convex hull and checking which points 
# are on tie lines (representing stable mixtures)
#------------------------------------------------------------------------------
def CalculateVoltAndLatt(handle,vname,cname,simplices,allpoints,stoich,e,Latt,\
   L2D):
   #calculate endpoint energies 
   EndPointE=[]
   EndPointE=np_append(EndPointE,Energy(0,e))
   EndPointE=np_append(EndPointE,Energy(1,e))
   EndPointE=np_append(EndPointE,Energy(2,e))
   
   #initialize weights of reaction
   weights=[0,0,0]

   #initialize auxilary arrays 
   xrot=[];yrot=[];zrot=[]
   ConcentAx=[];ConcentBy=[];ConcentCz=[]
   Ten=[]
   Isimplex=[]
   LattAvg=CalculateLattAvgFirst(Latt,allpoints)
   k=0
   l=0
   nsimplices=len(simplices)
   loss=[]
   for i in range(1,NConcentY0Points): 
          X=ConstY0Path(i,NConcentY0Points,ConcentY0)[0]
          Y=ConstY0Path(i,NConcentY0Points,ConcentY0)[1]
          for j in range(0,nsimplices): 
              if IsNotBoundary(j,simplices,allpoints) == 1:
                 # get the simplex id of current point (X,Y)
                 if isPointInPath(X,Y,allpoints[simplices[j]]) : 
                    Z=distance_to_hull( X, Y, simplices[j], allpoints)
                    xrot=np_append(xrot,[X])
                    yrot=np_append(yrot,[Y])
                    zrot=np_append(zrot,[Z])
                    #determine Li, Ag, Mn concentrations and the corresponding 
                    #energy 
                    xvec=GetConcentrations([X,Y])
                    ConcentAx=np_append(ConcentAx,xvec[0])
                    ConcentBy=np_append(ConcentBy,xvec[1])
                    ConcentCz=np_append(ConcentCz,xvec[2])
                    #the total energy is obtained from a linear interpolation of	
                    #of the convex hull 
                    Ten=np_append(Ten,[TotalEnergy(X,Y,Z,EndPointE)])
                    #also store the simplex id for every point i of the cut 
                    Isimplex=np_append(Isimplex,[j])

                    #determine reactin weights for the simplex 
                    lamb=GetReactionWeights([X,Y],allpoints[simplices[j]])

                    #calculate loss of phase B ... has index 1 in list allpoints 
                    loss=np_append(loss,LossOf(1,lamb,simplices[j]))

                    #average lattice parameter based on these weights 
                    LattAvg=np_append(LattAvg,[CalculateLattAvgChange(lamb,simplices[j],Latt)],axis=0)
                    
                    k=k+1

   #RotateHullCut( xrot, yrot);
   npointscut=k
   
   #*************************************VOLTAGE*****************************************
   #these will be stable points of the convex hull line
   d=(zrot[1]-zrot[0])/LenBetPts2D(xrot[1],yrot[1],xrot[0],yrot[0])
   dx=0
   lminimum=[]
   #the first point is a minimum 
   lminimum=np_append(lminimum,True)
   k=0
   l=0
   for i in range(0,npointscut-1):
      lminimum=np_append(lminimum,False)
      dnew=(zrot[i+1]-zrot[i])/LenBetPts2D(xrot[i+1],yrot[i+1],xrot[i],yrot[i])
      #skip first point 
      #check slope first
      if abs(dnew-d)>epsh:
         #update Li step
         dx=ConcentAx[i]-ConcentAx[k]
         #update tagger
         lminimum[i]=True
         #update slope minimum counter 
         d=dnew 
         k=k+1
   
   lminimum=np_append(lminimum,True)
   
   #clean double and triple entries:
   l=0
   k=0
   #create indexarray
   idxmin=[]
   idxmin=np_append(idxmin,0)
   for i in range(1,npointscut):
      if lminimum[i]:
         if i-l>3: 
            l=i
            k=k+1
            idxmin=np_append(idxmin,i)
         else:
            lminimum[i]=False

   k=k+1    
   lminimum[npointscut-1]=True
   idxmin=np_append(idxmin,npointscut-1)
   nmin=k+1

   #store the information about cut into array
   CutB = np_ndarray(shape=(npointscut,12), dtype=float, order='F')
   count=0
   for i in range(0,npointscut):
      CutB[i,0]=xrot[i]
      CutB[i,1]=yrot[i]
      CutB[i,2]=zrot[i]
      CutB[i,3]=ConcentAx[i]
      CutB[i,4]=ConcentBy[i]
      CutB[i,5]=ConcentCz[i]
      CutB[i,6]=Ten[i]
      CutB[i,7]=LattAvg[i,0]
      CutB[i,8]=LattAvg[i,1]
      CutB[i,9]=LattAvg[i,2]
      CutB[i,10]=LattAvg[i,3]
      CutB[i,11]=loss[i]

   #plot a point for every minima and compute voltages
   volt=[]
   elec=[]
   stab=0
   bina=0
   for i in range(0,nmin-1):
      if L2D == 1:
          PrintTo(handle,' set label at '+str(xrot[idxmin[i]])+\
             ','+str(yrot[idxmin[i]])+\
              ' "" front point pt 6 ps 2.0 lc rgb "#000000" ')
      else:
          PrintTo(handle,' set label at '+str(xrot[idxmin[i]])+','+\
             str(yrot[idxmin[i]])+','+str(zrot[idxmin[i]])+\
             '"" front  point  pt 6 ps 2.0 lc rgb "#000000"')
      
      k=idxmin[i]
      l=idxmin[i+1]
      dLi=(ConcentAx[l]-ConcentAx[k])
      volt=np_append(volt,-(Ten[l]-Ten[k]-dLi*EndPointE[0])/dLi)
      elec=np_append(elec,ConcentAx[k]+dLi/float(2))
   
   VoltageOfHull = np_ndarray(shape=(nmin-2,7), dtype=float, order='F')
   for i in range(0,nmin-2):
       VoltageOfHull[i,0]=elec[i]
       VoltageOfHull[i,1]=ConcentY0
       VoltageOfHull[i,2]=volt[i]
       k=idxmin[i]
       l=idxmin[i+1]
       VoltageOfHull[i,3]=Ten[l]
       VoltageOfHull[i,4]=Ten[k]
       VoltageOfHull[i,5]=ConcentAx[l]
       VoltageOfHull[i,6]=ConcentAx[k]
   np_savetxt(vname,VoltageOfHull,header='  x          Ag         V'+\
      '        E1         E2           Li1        Li2 ',fmt='%10.4f')

   # determine tie lines for binary 

   np_savetxt(cname,CutB,fmt='%12.6f')

#------------------------------------------------------------------------------
# prints header of GNUplot file
#------------------------------------------------------------------------------
def PrintHeader(LPlotBCut,L2D,MaxDist,handle):
   PrintTo(handle,'set terminal wxt size 800,800 enhanced')
   PrintTo(handle,' #set terminal pngcairo size 600, 600 enhanced'+\
                       'font "Helvetica, 10" #small figure' )       
   PrintTo(handle,' #set terminal pngcairo size 1200, 1200 enhanced'+\
                       'font "Helvetica, 16" ')
   PrintTo(handle,' #set terminal postscript hull.eps size 5.0,3.9'+\
                       'enhanced color font "Helvetica,16" linewidth 2')
   if LPlotBCut:
      PrintTo(handle, ' #set output "LixAg'+str(ConcentY0)+'_hull.png" ')
      PrintTo(handle, ' #set output "LixAg'+str(ConcentY0)+'_hull.eps" ')
   else:
      PrintTo(handle, ' #set output "LixAgy_hull.png" ' )
      PrintTo(handle, ' #set output "LixAgy_hull.eps" ' )
   
   PrintTo(handle, ' set style line 1 pt 7 lc rgb "'+HEXblue+'" ps 1'  )
   PrintTo(handle, ' set style line 2 pt 7 lc rgb "'+HEXred+'"  ps 1.5' )
   PrintTo(handle, ' set style line 3 pt 7 lc rgb "'+HEXred+'"  lw 0.5 ')
   PrintTo(handle, ' set style line 4 pt 7 lc rgb "#000000"  ps 0.1 ')  
   PrintTo(handle, ' set style line 5 pt 7 lc rgb "'+HEXgreen+'"  ps 1.5 ') 
   PrintTo(handle, ' set style line 6 pt 7 lc rgb "'+HEXred+'"  ps 1.5 ')  
   PrintTo(handle, ' set view 0,0')
   PrintTo(handle, ' unset ytics')
   PrintTo(handle, ' unset ztics')
   PrintTo(handle, ' unset xtics')
   if L2D == 1:
      PrintTo(handle, ' set xrange [-0.2:1.2]')
      PrintTo(handle, ' set yrange [-0.1:1.1]')
   else:
      PrintTo(handle, ' set xrange [-0.1:1.1]')
      PrintTo(handle, ' set yrange [-0.1:1.1]')
      PrintTo(handle, ' set zrange [-3:0.1]')
   PrintTo(handle, ' unset border')
   PrintTo(handle, ' # creating the palette by specifying H,S,V ')
   PrintTo(handle, ' h2 = 0/360.0')
   PrintTo(handle, ' h1 = 240/360.0')
   PrintTo(handle, ' set palette model HSV functions (1-gray)*(h2-h1)+h1,1,0.68\
      ;set title "Distance to hull [eV]')
   PrintTo(handle,'#set palette rgb 7,5,15;set title "Distance to hull [eV]')
   PrintTo(handle,'#set palette rgb 3,11,6;set title "Distance to hull [eV]')
   PrintTo(handle,'#set palette rgb 23,28,3;set title "Distance to hull [eV]')
   PrintTo(handle,'#set palette rgb 21,22,23;set title "Distance to hull [eV]')
   PrintTo(handle,'#set palette rgb 30,31,32;set title "Distance to hull [eV]')
   PrintTo(handle,'#set palette rgb 30,13,10;set title "Distance to hull [eV]')
   PrintTo(handle,'#set palette rgb 34,35,36;set title "Distance to hull [eV]')
   PrintTo(handle, ' set cbrange [0:'+str(MaxDist)+']')

#------------------------------------------------------------------------------
# prints auxilary lines to GNUplot file
#------------------------------------------------------------------------------
def PrintAuxilary(L2D,handle):   
   #number of intermediate auxilary lines 
   ntot=10
   if L2D ==1 :
      PrintTo(handle, ' set label "Ag [%]" at '+str(.4*np_cos(pi/3)-0.1)+','+\
         str(.4*np_sin(pi/3)+0.1)+' rotate by 60 left')
      PrintTo(handle, ' set label "Li [%]" at '+str(.5*np_cos(pi/3)+0.6)+','+\
         str(.5*np_sin(pi/3)+0.1)+' rotate by -60 left')
      PrintTo(handle, ' set label "Mn_{8}O_{16} [%]" at 0.4,-0.1')
   else:
      PrintTo(handle, ' set label "Ag [%]" at '+str(.4*np_cos(pi/3)-0.1)+','+\
         str(.4*np_sin(pi/3)+0.1)+',0 rotate by 60 left')
      PrintTo(handle, ' set label "Li [%]" at '+str(.5*np_cos(pi/3)+0.6)+','+\
         str(.5*np_sin(pi/3)+0.1)+',0 rotate by -60 left')
      PrintTo(handle, ' set label "Mn_{8}O_{16} [%]" at 0.4,-0.1,0')
   
   for i in range(0,ntot+1):
       xi=float(i)/ntot
       y1=xi*np_cos(pi/3)
       y2=xi*np_sin(pi/3)
       y01=np_cos(pi/3)
       y02=np_sin(pi/3)
       if L2D == 1 :
          PrintTo(handle, ' set arrow from '+str(xi)+',0 to '+str(y1)+','+
             str(y2)+' nohead ls 3')
          PrintTo(handle, ' set arrow from '+str(y01-y1+xi)+','+str(y02-y2)+\
             ' to '+str(xi)+',0 nohead ls 3')
          PrintTo(handle, ' set arrow from '+str(y01-y1)+','+str(y02-y2)+' to '\
             +str(y01-y1+xi)+','+str(y02-y2)+' nohead ls 3')
       else:
          PrintTo(handle, ' set arrow from '+str(xi)+',0,0 to '+str(y1)+','\
             +str(y2)+',0 nohead ls 3')
          PrintTo(handle, ' set arrow from '+str(y01-y1+xi)+','+str(y02-y2)+\
             ',0 to '+str(xi)+',0,0 nohead ls 3')
          PrintTo(handle, ' set arrow from '+str(y01-y1)+','+str(y02-y2)+\
             ',0 to '+str(y01-y1+xi)+','+str(y02-y2)+',0 nohead ls 3')
   
       if i==0 :  
           PrintTo(handle, ' set label "0" at -0.05,0.00')
       else:
           PrintTo(handle,' set label "'+str(i*10)+'" at '+str(y1-0.1)+','+\
              str(y2+0.005))
           PrintTo(handle,' set label "'+str(i*10)+'" at '+str(1-(xi-0.01))+\
              ',-0.01 rotate by -60 left')
           PrintTo(handle,' set label "'+str(i*10)+'" at '+str(y01-y1+xi+0.01)+\
              ','+str(y02-y2)+' rotate by 60 left')
   
#------------------------------------------------------------------------------
# prints simplices to GNUplot file 
#------------------------------------------------------------------------------
def PrintSimplices(L2D,handle,simplices,points,bhull1,bhull2,bhull3,idx1,idx2,idx3):
   nsimplices=len(simplices)
   if L2D == 1 : 
      #polygons
      for i in range(0,nsimplices):
          if IsNotBoundary(i,simplices,points) == 1:
             PrintTo(handle,' set object '+str(i+1)+" polygon from \\")
             npoints = len(points[simplices[i],:])
             for j in range(0,npoints):
                 PrintTo(handle,'   '+str(points[simplices[i],:][j,0])+\
                    ','+str(points[simplices[i],:][j,1])+' to \\')
             PrintTo(handle,'   '+str(points[simplices[i],:][0,0])+','+\
                str(points[simplices[i],:][0,1]))
             PrintTo(handle,'   ')
   else:
      for i in range(0,nsimplices):
          if IsNotBoundary(i,simplices,points) == 1:
             PrintTo(handle,' set object '+str(i+1)+" polygon from \\")
             npoints = len(points[simplices[i],:])
             for j in range(0,npoints):
                 PrintTo(handle,'   '+str(points[simplices[i],:][j,0])+\
                    ','+str(points[simplices[i],:][j,1])+','+\
                    str(points[simplices[i],:][j,2])+' to \\')
             PrintTo(handle,'   '+str(points[simplices[i],:][0,0])+\
                ','+str(points[simplices[i],:][0,1])+','+\
                str(points[simplices[i],:][0,2]))
             PrintTo(handle,'   ')
    
      #boundary
      k=nsimplices+1
      npoints = len(bhull1.vertices)
      PrintTo(handle,' set object '+str(k)+" polygon back from \\" )
      for j in range(0,npoints):
         i=int(idx1[bhull1.vertices[j]])
         PrintTo(handle,' '+str(points[i,0])+','+str(points[i,1])+\
         ','+str(points[i,2])+' to \\')
      i=int(idx1[bhull1.vertices[0]])
      PrintTo(handle,' '+str(points[i,0])+','+str(points[i,1])+','+\
         str(points[i,2]))
      PrintTo(handle,'   ')

      k=nsimplices+2
      npoints = len(bhull2.vertices)
      PrintTo(handle,' set object '+str(k)+" polygon back from \\" )
      for j in range(0,npoints):
         i=int(idx2[bhull2.vertices[j]])
         PrintTo(handle,' '+str(points[i,0])+','+str(points[i,1])+\
         ','+str(points[i,2])+' to \\')
      i=int(idx2[bhull2.vertices[0]])
      PrintTo(handle,' '+str(points[i,0])+','+str(points[i,1])+','+\
         str(points[i,2]))
      PrintTo(handle,'   ')

      k=nsimplices+3
      npoints = len(bhull3.vertices)
      PrintTo(handle,' set object '+str(k)+" polygon back from \\" )
      for j in range(0,npoints):
         i=int(idx3[bhull3.vertices[j]])
         PrintTo(handle,' '+str(points[i,0])+','+str(points[i,1])+\
         ','+str(points[i,2])+' to \\')
      i=int(idx3[bhull3.vertices[0]])
      PrintTo(handle,' '+str(points[i,0])+','+str(points[i,1])+','+\
         str(points[i,2]))
      PrintTo(handle,'   ')

#------------------------------------------------------------------------------
# Set style of polygons 
#------------------------------------------------------------------------------
def PrintPolygonStyle(L2D,handle,simplices,points):
   nsimplices=len(simplices)
   for i in range(0,nsimplices):
       if IsNotBoundary(i,simplices,points) == 1:
          PrintTo(handle,' set object '+str(i+1)+' fc rgb "'+HEXred+\
             '" fs solid 0.3 border rgb"'+HEXred+'" lw 1.5 ')
   
   if L2D == 0 : 
      PrintTo(handle,' set object '+str(nsimplices+1)+'fc rgb "'+HEXred+\
         '" fs solid 0.3 border rgb "'+HEXred+'" lw 1.5' )
      PrintTo(handle,' set object '+str(nsimplices+2)+'fc rgb "'+HEXred+\
         '" fs solid 0.3 border rgb "'+HEXred+'" lw 1.5' )
      PrintTo(handle,' set object '+str(nsimplices+3)+'fc rgb "'+HEXred+\
         '" fs solid 0.3 border rgb "'+HEXred+'" lw 1.5' )

#------------------------------------------------------------------------------
#show normal vectors of polygons
#------------------------------------------------------------------------------
def PrintNormalVectors(L2D,simplices,points,e):
   #calculate endpoint energies 
   EndPointE=[]
   EndPointE=np_append(EndPointE,Energy(0,e))
   EndPointE=np_append(EndPointE,Energy(1,e))
   EndPointE=np_append(EndPointE,Energy(2,e))
   nsimplices=len(simplices)
   for j in range(0,nsimplices): 
       if IsNotBoundary(j,simplices,points) == 1:
          p=CoM(simplices[j],points)
          n=norm_vector(simplices[j],points)
          if L2D == 1 :
             PrintTo(handle,' set arrow from '+str(p[0])+','+str(p[1])+','+\
                str(p[2])+' to '+str(p[0]+n[0]/10)+','+str(p[1]+n[1]/10))
          else:
             PrintTo(handle,' set arrow from '+str(p[0])+','+str(p[1])+','+\
                str(p[2])+' to '+str(p[0]+n[0]/10)+','+str(p[1]+n[1]/10)+\
                ','+str(p[2]+n[2]/10)+' ls 1 ')

   #print voltages
   for i in range(0,nsimplices):
       p1x=points[simplices[i],:][0,0]
       p1y=points[simplices[i],:][0,1]
       p1z=points[simplices[i],:][0,2]
       p2x=points[simplices[i],:][1,0]
       p2y=points[simplices[i],:][1,1]
       p2z=points[simplices[i],:][1,2]
       p3x=points[simplices[i],:][2,0]
       p3y=points[simplices[i],:][2,1]
       p3z=points[simplices[i],:][2,2]
       px=p1x+(p2x-p1x)/2
       py=p1y+(p2y-p1y)/2
       pz=p1z+(p2z-p1z)/2
       ux=p1x+(p3x-p1x)/2
       uy=p1y+(p3y-p1y)/2
       uz=p1z+(p3z-p1z)/2
       vx=p2x+(p3x-p2x)/2
       vy=p2y+(p3y-p2y)/2
       vz=p2z+(p3z-p2z)/2
       i1=simplices[i][0]
       i2=simplices[i][1]
       i3=simplices[i][2]
       if L2D == 1 :
          PrintTo(handle,' set label "'+str(Voltage(i1,i2,xx,yy,zz,e,EndPointE\
             ))+'" at '+str(px)+','+str(py)+' font "Helvetica,8" tc rgb'+\
             ' "#009933" ')
          PrintTo(handle,' set label "'+str(Voltage(i1,i3,xx,yy,zz,e,EndPointE\
             ))+'" at '+str(ux)+','+str(uy)+' font "Helvetica,8" tc rgb'+\
             ' "#009933" ')
          PrintTo(handle,' set label "'+str(Voltage(i2,i3,xx,yy,zz,e,EndPointE\
             ))+'" at '+str(vx)+','+str(vy)+' font "Helvetica,8" tc rgb'+\
             ' "#009933" ')
       else:
          PrintTo(handle,' set label "'+str(Voltage(i1,i2,xx,yy,zz,e,EndPointE\
             ))+'" at '+str(px)+','+str(py)+','+str(pz)+' font "Helvetica,8"'+\
             ' tc rgb "#009933"')
          PrintTo(handle,' set label "'+str(Voltage(i1,i3,xx,yy,zz,e,EndPointE\
             ))+'" at '+str(ux)+','+str(uy)+','+str(uz)+' font "Helvetica,8"'+\
             ' tc rgb "#009933"')
          PrintTo(handle,' set label "'+str(Voltage(i2,i3,xx,yy,zz,e,EndPointE\
             ))+'" at '+str(vx)+','+str(vy)+','+str(vz)+' font "Helvetica,8"'+\
             ' tc rgb "#009933"')

#------------------------------------------------------------------------------
# print Tail of gnuplot file 
#------------------------------------------------------------------------------
def PrintTail(handle,hname,cname,LPlotBCut,L2D):
   
   if L2D == 1 :
      PrintTo(handle,'plot "'+hname+'" u ($1):($2):($4) w p ls 2 palette'+\
         ' notitle\\')
      if LPlotBCut:
         PrintTo(handle,'    ,"'+cname+'" u ($1):($2):($3) w p ls 4 '+\
         ' notitle\\')
   else:
      PrintTo(handle,'splot "'+hname+'" u ($1):($2):($3):($4) w p ps 1.5 pt 7'+\
         'lc palette notitle\\')
      if LPlotBCut:
         PrintTo(handle,'  ,"'+cname+'" u ($1):($2):($3):(0) w p ps 0.5 pt 7'+\
         ' lc palette  notitle\\')
   PrintTo(handle, '    ')

#------------------------------------------------------------------------------
# For visualizing the structural change
#------------------------------------------------------------------------------
def PrintStructureChange(handle,cname):
   PrintTo(handle,'set terminal wxt size 800,800 enhanced')
   PrintTo(handle,' #set terminal pngcairo size 600, 600 enhanced'+\
                       ' font "Helvetica, 10" #small figure' )       
   PrintTo(handle,' #set terminal pngcairo size 1200, 1200 enhanced'+\
                       ' font "Helvetica, 16" ')
   PrintTo(handle,' #set terminal postscript hull.eps size 5.0,3.9'+\
                       'enhanced color font "Helvetica,16" linewidth 2')
   PrintTo(handle,' #set output "StructuralChange_y'+str(ConcentY0)+'.png" ')
   PrintTo(handle,' #set output "StructuralChange_y'+str(ConcentY0)+'.eps" ')
   PrintTo(handle,' set style line 1 pt 7 lc rgb "'+HEXblue+'" lw 2.0 ps 1'  )
   PrintTo(handle,' set style line 2 pt 7 lc rgb "'+HEXred+'"  lw 2.0ps 1.5' )
   PrintTo(handle,' set style line 3 pt 7 lc rgb "'+HEXgreen+'"  lw 2.0 ')
   PrintTo(handle,' set style line 4 pt 7 lc rgb "#000000"  ps 0.1 ')  
   PrintTo(handle,' set style line 5 pt 7 lc rgb "'+HEXgreen+'"  ps 1.5 ') 
   PrintTo(handle,' set style line 6 pt 7 lc rgb "'+HEXred+'"  ps 1.5 ')  
   PrintTo(handle,' set xrange [0.0:'+str(MaxX)+']')
   PrintTo(handle,' set yrange [9.0:12.0]')
   PrintTo(handle,' set xlabel "'+CreateSystemName(["x",str(ConcentY0),""])+'"')
   PrintTo(handle,' set ylabel "lattice parameter [A]"')
   PrintTo(handle,' set y2label "unit cell volume [A]"')
   PrintTo(handle,' set xtics (1,2,3,4,5,6,7,8)')
   PrintTo(handle,' unset y2tics')
   PrintTo(handle,' set key top left box')
   PrintTo(handle,' set y2tics  ("150" 8.5, "200" 9, "250" 9.5, "300" 10, '+\
     '"350" 10.5,"400"  11, "450"  11.5, "500" 12, "550" 12.5)')
   PrintTo(handle,'plot "'+cname+'" u ($4):($9) w l ls 1 title "a"\\')
   PrintTo(handle,'    ,"'+cname+'" u ($4):($10) w l ls 2 title "b"\\')
   PrintTo(handle,'    ,"'+cname+'" u ($4):($11/100+7) w l ls 3 title "Vol."\\')
   PrintTo(handle,'    ')

#------------------------------------------------------------------------------
# For visualizing the Phase B loss
#------------------------------------------------------------------------------
def PrintLossOf(handle,cname):
   PrintTo(handle,'set terminal wxt size 800,600 enhanced font "Helvetica,16"')
   PrintTo(handle,' #set terminal pngcairo size 800, 600 enhanced'+\
                       ' font "Helvetica, 10" #small figure' )       
   PrintTo(handle,' #set terminal pngcairo size 1200, 768 enhanced'+\
                       ' font "Helvetica, 16" ')
   PrintTo(handle,' #set terminal postscript hull.eps size 5.0,3.9'+\
                       'enhanced color font "Helvetica,16" linewidth 2')
   PrintTo(handle,' #set output "LossInB_y'+str(ConcentY0)+'.png" ')
   PrintTo(handle,' #set output "LossInB_y'+str(ConcentY0)+'.eps" ')
   PrintTo(handle,' set style line 1 pt 7 lc rgb "'+HEXblue+'" lw 2.0 ps 1'  )
   PrintTo(handle,' set style line 2 pt 7 lc rgb "'+HEXred+'"  lw 2.0ps 1.5' )
   PrintTo(handle,' set style line 3 pt 7 lc rgb "'+HEXgreen+'"  lw 2.0 ')
   PrintTo(handle,' set style line 4 pt 7 lc rgb "#000000"  ps 0.1 ')  
   PrintTo(handle,' set style line 5 pt 7 lc rgb "'+HEXgreen+'"  ps 1.5 ') 
   PrintTo(handle,' set style line 6 pt 7 lc rgb "'+HEXred+'"  ps 1.5 ')  
   PrintTo(handle,' set xrange [0.0:'+str(MaxX)+']')
   PrintTo(handle,' set yrange [0:10]')
   PrintTo(handle,' set ytics '+\
     '("10"1,"20"2,"30"3,"40"4,"50"5,"60"6,"70"7,"80"8,"90"9,"100"10)')
   PrintTo(handle,' set y2tics '+\
     '("1.02"1,"1.04"2,"1.06"3,"1.08"4,"1.1"5,'+\
     '"1.12"6,"1.14"7,"1.16"8,"1.18"9,"1.20"10)')
   PrintTo(handle,' set xlabel "'+CreateSystemName(["x",str(ConcentY0),""])+'"')
   PrintTo(handle,' set ylabel "Loss in %"')
   PrintTo(handle,' set y2label "Lattice constant ratio a/b "')
   PrintTo(handle,' set xtics (1,2,3,4,5,6,7,8)')
   PrintTo(handle,' set key top left box')
   PrintTo(handle,'plot "'+cname+'" u ($4):(10*$12/'+str(ConcentY0)+')'+\
     ' w l ls 1 title "Ag-Loss"\\')
   PrintTo(handle,'    ,"'+cname+'" u ($4):(50*($10/$9-1)) w l ls 2 title "a/b"\\')
   PrintTo(handle,'    ')
  
  
