import wx
import icons35 as ico
from matplotlib.backends.backend_wx import NavigationToolbar2Wx


# -----------------------------------------------------------------------------
# prints reaction equation into logger
# -----------------------------------------------------------------------------
def PrintReactionInfo(Project,x,y,logger,key):
   from poscar import CalculateFormation
   from poscar import CalculateReaction

   text = CalculateFormation(Project,x,y,key)
   logger.Add2Log(text)
   text = CalculateReaction(Project,x,y)
   logger.Add2Log(text)

   # print list of composition with higher energies 
   if len( key ) > 0 : 
      l = len( Project.db[key]["SameBC"] ) 
      if l >= 1 : 
         text = '   {:} compositions with higher energy found:'.format( l )
	 f=Project.dbF[ key ]["maxForce"]

         text += '\n   Energy [eV]   maximum residual force on ions (x,y,z) [eV/A]   name'
	 text += '\n       {:.4f}              {:.4f}    {:.4f}   {:.4f}            {:}'.\
           	format( Project.dbF[ key ]["formationE"], f[0],f[1],f[2], key )
         logger.Add2Log(text)
         for l in Project.db[key]["SameBC"] : 
	    f=Project.dbF[ l ]["maxForce"]
            text = '       {:.4f}              {:.4f}    {:.4f}   {:.4f}            {:}'.\
            	format( Project.dbF[ l ]["formationE"], f[0],f[1],f[2], l )
            logger.Add2Log(text)
   

#------------------------------------------------------------------------------
#creates hexagon patches for plotting binding energy strength
#------------------------------------------------------------------------------
def CreateHexagonPatches(nsample, sx, status=None, \
   Reaction = None, Hill=None, Hull=None,NormIDX=None, function=None):
   from basis import distance_to_hull,GetBaryCentricCoordinates
   from matplotlib.patches import Polygon
   from numpy import array as np_array
   from numpy import zeros as np_zeros
   from poscar import FindSimplices
   from poscar import GetConcentrations2
   from numpy import append as np_append

   if nsample < 0 : 
      nsample =64 
  

   #small number 
   small=1.E-6

   #sample the plane with following lattice vectors 
   sq3=1.732050807568877
   sq32=0.8660254037844386
   v1x= 0.866025403784439;v1y= 0.500000000000000
   v2x= 0.000000000000000;v2y= 1.000000000000000
   v3x=-0.866025403784439;v3y= 0.500000000000000
   v4x=-0.866025403784439;v4y=-0.500000000000000
   v5x= 0.000000000000000;v5y=-1.000000000000000
   v6x= 0.866025403784439;v6y=-0.500000000000000

   #store simplices
   simplices2 = sx.hull.tri.simplices.copy() 
   simplices = sx.hull.inner_simplices.copy() 

   allpoints2=np_zeros( ( len(sx.hull.points),3) ) 
   allpoints2[:,0]=sx.hull.points[:,0]
   allpoints2[:,1]=sx.hull.points[:,1]
   allpoints2[:,2]=sx.hull.distance[:]

   x=[]
   y=[]
   z=[]
   X=[];Y=[]
   X=np_append(X,1.0000000000000000/float(nsample))
   X=np_append(X,0)
   Y=np_append(Y,1.0000000000000000/(2*nsample))
   Y=np_append(Y,sq32/nsample)
   a=(sq3/3)/nsample
   patches=[]
   colors=[]


   #loop over sample points 
   if status is not None: 
      status.SetRange(int(nsample*nsample))
      if Hull : 
         text = 'Updating formation energy contour'
      elif Hill:
         if NormIDX is None : 
            text = 'Updating stability contour'
         else:
            text = 'Updating reaction energy contour'
      status.SetStatusText(text)
      iter = 0 

   i=-1
   for J in range(0,nsample):
      for I in range(0,nsample):
         if status is not None: 
            status.SetProgress(iter)
            iter+=1
         xp=I*X[0]+J*Y[0]
         yp=I*X[1]+J*Y[1]

         # determine simplex of hull and simplex of hill 
         if Hull or Hill :
            idx_i, idx_u = FindSimplices( sx, xp, yp, Hill, Hull ) 
         # both must be positive when found 
         # print round(xp,4),round(yp,4),'simplex',idx_i,idx_u
         zp=0
         lappend=False
         if Hill:
            if idx_i >= 0:
               zp=-distance_to_hull( xp, yp, simplices2[idx_i], \
                allpoints2)
               lappend=True
         elif Hull:
            if idx_u >= 0 and idx_u < len(simplices):
               zp=zp+distance_to_hull( xp, yp, simplices[idx_u], \
                  sx.hull.points)
               lappend=True
         if lappend: 
            #scale 
            if NormIDX is not None:
               Xb=GetBaryCentricCoordinates( [xp,yp] ) 
               Xn, idx =  GetConcentrations2([xp,yp],NormIDX=NormIDX,ReturnIDX=True)
               alpha = Xn[ idx ] / Xb[ idx ] 
               zp  = zp*sum(Xn)
                    
            colors=np_append(colors,zp)
            z=np_append(z,zp)
            x=np_append(x,xp)
            y=np_append(y,yp)
            i+=1
   
            ox=xp
            oy=yp
            #endpoints
            if i==0:
               polygon = Polygon([np_array([0.5/nsample,0]),
                                  np_array([a*v1x,a*v1y]),
                                  np_array([0.25/nsample, 0.5/nsample*sq32]),
                                  np_array([0,0])], True )
               patches.append(polygon)
            #elif i < nsample : 
            elif J==0:

               polygon = Polygon([np_array([ox+0.5/nsample,0]),
                                  np_array([ox+a*v1x,a*v1y]),
                                  np_array([ox+a*v2x,oy+a*v2y]),
                                  np_array([ox+a*v3x,oy+a*v3y]),
                                  np_array([ox-0.5/nsample,0])], True )
               patches.append(polygon)
            #left border
            #elif abs(oy-sq3*ox)<small:
            elif I==0:
               polygon = Polygon([np_array([ox+0.5/nsample,oy]),
                                  np_array([ox+a*v1x,oy+a*v1y]),
                                  np_array([ox+0.25/nsample,oy+0.5/nsample*sq32]),
                                  np_array([ox-0.25/nsample,oy-0.5/nsample*sq32]),
                                  np_array([ox+a*v5x,oy+a*v5y]),
                                  np_array([ox+a*v6x,oy+a*v6y])
                                  ], True )
               patches.append(polygon)
            #right border 
            #elif abs(oy+sq3*ox-sq3)<small:
            elif I+J == nsample:
               polygon = Polygon([  
                                  np_array([ox+0.25/nsample,oy-0.5/nsample*sq32]),
                                  np_array([ox-0.25/nsample,oy+0.5/nsample*sq32]),
                                  np_array([ox+a*v3x,oy+a*v3y]),
                                  np_array([ox+a*v4x,oy+a*v4y]),
                                  np_array([ox+a*v5x,oy+a*v5y])
                                  ], True )
               patches.append(polygon)
            
            #endpoints 
            elif oy < sq3*ox and oy+small < -sq3*ox+sq3 and oy>0:
               polygon = Polygon([np_array([ox+a*v1x,oy+a*v1y]),
                                  np_array([ox+a*v2x,oy+a*v2y]),
                                  np_array([ox+a*v3x,oy+a*v3y]),
                                  np_array([ox+a*v4x,oy+a*v4y]),
                                  np_array([ox+a*v5x,oy+a*v5y]),
                                  np_array([ox+a*v6x,oy+a*v6y])], True )
               patches.append(polygon)

               
   # add right endpoint
   polygon = Polygon([  
                      np_array([1,0]),
                      np_array([1-0.25/nsample,0.5/nsample*sq32]),
                      np_array([1+a*v3x,a*v3y]),
                      np_array([1-0.5/nsample,0])
                      ], True )
   patches.append(polygon)
   colors=np_append(colors,0)
   # add endpoints here 
   x=np_append(x,1)
   y=np_append(y,0)
   z=np_append(z,0)

   # add top endpoint
   ox=0.5; oy=sq32
   polygon = Polygon([  
                      np_array([ox,oy]),
                      np_array([ox-0.25/nsample,oy-0.5/nsample*sq32]),
                      np_array([ox+a*v5x,oy+a*v5y]),
                      np_array([ox+0.25/nsample,oy-0.5/nsample*sq32])
                      ], True )
   patches.append(polygon)
   colors=np_append(colors,0)
   # add endpoints here 
   x=np_append(x,ox)
   y=np_append(y,oy)
   z=np_append(z,0)

   if status is not None: 
      status.prog.Hide()

   
   # add endpoints 
   if function : 
       return x,y,z

   else:
      return patches, colors

#------------------------------------------------------------------------------
# wrapper to calculate binding energy patches
#------------------------------------------------------------------------------
def GetBindingContour( settings , status = None ):
   import matplotlib.tri as tri
   #
   # calculate triangulation points
   #
   x,y,z = CreateHexagonPatches( settings.nsample, settings, \
           status = status,  Hull=True, function=True )
   T = tri.Triangulation(x,y)
   return T,x,y,z

#------------------------------------------------------------------------------
# wrapper to calculate stability patches
#------------------------------------------------------------------------------
def GetStabilityContour( settings , status = None ):
   import matplotlib.tri as tri
   #
   # calculate triangulation points
   #
   x,y,z = CreateHexagonPatches( settings.nsample, settings, \
           status = status,  Hill=True, function=True )
   T = tri.Triangulation(x,y)
   return T,x,y,z

#------------------------------------------------------------------------------
# calculates reaction energies
#------------------------------------------------------------------------------
def GetReactionEnergiesOfPoints(sx):
   from basis import distance_to_hull,GetBaryCentricCoordinates
   from numpy import array as np_array
   from numpy import zeros as np_zeros
   from poscar import FindSimplices
   from poscar import GetConcentrations2
   from numpy import append as np_append

   if sx.nsample < 0 : 
      nsample =64 
   else:
      nsample =64 
  

   #small number 
   small=1.E-6

   #store simplices
   simplices2 = sx.hull.tri.simplices.copy() 
   simplices = sx.hull.inner_simplices.copy() 

   npoints=len(sx.hull.points)
   allpoints2=np_zeros( ( npoints,3) ) 
   allpoints2[:,0]=sx.hull.points[:,0]
   allpoints2[:,1]=sx.hull.points[:,1]
   allpoints2[:,2]=sx.hull.distance[:]

   x=[]
   y=[]
   z=[]

   for i in range(0,npoints):
         xp=sx.hull.points[i,0]
         yp=sx.hull.points[i,1]

         # determine simplex of hull and simplex of hill 
         idx_i, idx_u = FindSimplices( sx, xp, yp, True, False ) 
         # both must be positive when found 
         zp=0
         if idx_i >= 0:
            zp=-distance_to_hull( xp, yp, simplices2[idx_i], \
             allpoints2)
            Xb=GetBaryCentricCoordinates( [xp,yp] ) 
            Xn, idx =  GetConcentrations2([xp,yp],NormIDX=sx.NormIDX,ReturnIDX=True)
            alpha = Xn[ idx ] / Xb[ idx ] 
            zp  = zp*sum(Xn)
            z = np_append(z,zp)
   if (len(z)!=npoints) :
      print 'ERROR not all points found'
   return z 

#------------------------------------------------------------------------------
# wrapper to calculate binding energy patches
#------------------------------------------------------------------------------
def GetReactionContour( settings , status = None ):
   import matplotlib.tri as tri
   #
   # calculate triangulation points
   #
   x,y,z = CreateHexagonPatches( settings.nsample, settings, \
           status = status,  Hill=True, NormIDX=settings.NormIDX, \
           function=True )
   T = tri.Triangulation(x,y)
   #
   # Store Reaction energy for every point
   #
   settings.hull.ReactionEnergies = GetReactionEnergiesOfPoints(settings)
   
   return T,x,y,z

#------------------------------------------------------------------------------
# creates coordinates for grid of ternary 
#------------------------------------------------------------------------------
def GetBorderCollection(settings,For3D=False ):
   from numpy import append as np_append
   from numpy import array as np_array
   from matplotlib.collections import LineCollection
   from mpl_toolkits.mplot3d.art3d import  Line3DCollection
  
   #coordinates of endpoints 
   A=[0.0000000000000000,0.0000000000000000]
   B=[1.0000000000000000,0.0000000000000000]
   C=[0.5000000000000000,0.8660254037844386]

   #set some default values: 
   if settings.lwBorder is not None: 
      LW = settings.lwBorder
   else :
      LW = '2'

   if settings.BorderColor is not None: 
      COLOR = settings.BorderColor
   else :
      COLOR = 'black'

   ZORDER = 8 

   if settings.ls is None: 
      LS = '-'
   else:
      LS = settings.ls

   y = []
   z = []
   lines=[]
   lc = []
   x = [ A[0],B[0],C[0],A[0] ]
   y = [ A[1],B[1],C[1],A[1] ]
   z = [0,0,0,0]
   if For3D :
      lines.append(list( zip(x,y,z) ))
   else:
      lines.append(list( zip(x,y) ))
   lc.append( COLOR )

   if For3D:
      line_segments = Line3DCollection( lines, linewidths = LW , colors=lc, \
         linestyles=LS ,zorder=ZORDER)
   else:
      line_segments = LineCollection( lines, linewidths = LW , colors=lc, \
         linestyles=LS ,zorder=ZORDER)

   return line_segments
    

#------------------------------------------------------------------------------
# creates coordinates for grid of ternary 
#------------------------------------------------------------------------------
def GetGridCollection(settings,For3D=False ):
   from numpy import append as np_append
   from numpy import array as np_array
   from matplotlib.collections import LineCollection
   from mpl_toolkits.mplot3d.art3d import  Line3DCollection
  
   #coordinates of endpoints 
   A=[0.0000000000000000,0.0000000000000000]
   B=[1.0000000000000000,0.0000000000000000]
   C=[0.5000000000000000,0.8660254037844386]

   #set some default values: 
   if settings.lwGrid is not None: 
      LW = settings.lwGrid
   else :
      if settings.UseUniform:
         LW = '0.1'
      else:
         LW = '0.1'
 

   if settings.GridColors is not None: 
      x,y,z =settings.GridColors
      if settings.UseUniform:
         COLOR = (x,y,z)
      else:
         COLOR = (y,x,z)
   else :
      COLOR = ('gray','gray','gray')

   ZORDER =8 

   if settings.ls is not None: 
      LS = settings.ls
   else :
      LS = 'dashed'
   
   x = []
   y = []
   z = []
   lines=[]
   ls = []
   lw = []
   lc = []
   seg=[]
   if settings.UseUniform: 
      #y-axis
      for i in range(0,10): 
          x = []
          y = []
          z = []
          seg=[]
          x.append(float(i)/10.*C[0])
          y.append(float(i)/10.*C[1])
          x.append(1.-float(i)/10.*C[0])
          y.append(float(i)/10.*C[1])
          z.append(0)
          z.append(0)
          if For3D :
             seg.append(list( zip(x,y,z) ))
          else:
             seg.append(list( zip(x,y) ))
          lines.append(seg[0])
          ls.append( LS ) 
          lw.append( LW ) 
          lc.append( COLOR[1] )

      #x-axis
      for i in range(0,10): 
          x = []
          y = []
          z = []
          seg=[]
          x.append(float(i)/10.)
          y.append(0.)
          x.append(1.-float(10-i)/10.*C[0])
          y.append(float(10-i)/10.*C[1])
          z.append( 0 ) 
          z.append( 0 ) 
          if For3D :
             seg.append(list( zip(x,y,z) ))
          else:
             seg.append(list( zip(x,y) ))
          lines.append(seg[0])
          lw.append( LW ) 
          ls.append( LS ) 
          lc.append( COLOR[0] )

      #z-axis
      for i in range(0,10): 
          x = []
          y = []
          z = []
          seg=[]
          x.append(float(i+1)/10.*C[0])
          y.append(float(i+1)/10.*C[1])
          x.append(float(i+1)/10.)
          y.append(0)
          z.append( 0 ) 
          z.append( 0 ) 
          if For3D :
             seg.append(list( zip(x,y,z) ))
          else:
             seg.append(list( zip(x,y) ))
          lines.append(seg[0])
          lw.append( LW ) 
          ls.append( LS ) 
          lc.append( COLOR[2] )

      #add ticks for horizontal lines
      tic=0.025
      for i in range(0,11): 
          x = []
          y = []
          z = []
          seg=[]
          x.append(float(i)/10.*C[0])
          y.append(float(i)/10.*C[1])
          x.append(float(i)/10.*C[0]-tic*C[0])
          y.append(float(i)/10.*C[1]+tic*C[1])
          z.append( 0 ) 
          z.append( 0 ) 
          if For3D :
             seg.append(list( zip(x,y,z) ))
          else:
             seg.append(list( zip(x,y) ))
          lines.append(seg[0])
          lw.append( LW ) 
          ls.append( 'solid' ) 
          lc.append( COLOR[2] )

      #add ticks in solid lines for AC line 
      for i in range(0,11): 
          x = []
          y = []
          z = []
          seg=[]
          x.append(float(i)/10.)
          y.append(0)
          x.append(-tic*C[0]+float(i)/10.)
          y.append(-tic*C[1])
          z.append( 0 ) 
          z.append( 0 ) 
          if For3D :
             seg.append(list( zip(x,y,z) ))
          else:
             seg.append(list( zip(x,y) ))
          lines.append(seg[0])
          lw.append( LW ) 
          ls.append( 'solid' ) 
          lc.append( COLOR[0] )
      #add ticks in solid lines for AC line 
      for i in range(0,11): 
          x = []
          y = []
          z = []
          seg=[]
          x.append(1.-float(10-i)/10.*C[0])
          y.append(float(10-i)/10.*C[1])
          x.append(1.-float(10-i)/10.*C[0]+tic)
          y.append(float(10-i)/10.*C[1])
          z.append( 0 ) 
          z.append( 0 ) 
          if For3D :
             seg.append(list( zip(x,y,z) ))
          else:
             seg.append(list( zip(x,y) ))
          lines.append(seg[0])
          lw.append( LW ) 
          ls.append( 'solid' ) 
          lc.append( COLOR[1] )

   else:
      if settings.xi is not None: 
         for Xi in settings.xi: 
             x.append( float(Xi)/float(1+Xi) ) 
             y.append( 0  )
             x.append( C[0] )
             y.append( C[1] )
             ls.append( LS ) 
             z.append( 0 ) 
             z.append( 0 ) 
         if For3D :
            seg.append(list( zip(x,y,z) ))
         else:
            seg.append(list( zip(x,y) ))
         lc.append( COLOR[0] )

      x = []
      y = []
      z = []
      if settings.yi is not None: 
         for Yi in settings.yi: 
             l = float(Yi)/float(1+Yi)
             x.append( l*C[0] ) 
             y.append( l*C[1] )
             x.append( B[0] )
             y.append( B[1] )
             z.append( 0 ) 
             z.append( 0 ) 
             ls.append( LS ) 
         if For3D :
            seg.append(list( zip(x,y,z) ))
         else:
            seg.append(list( zip(x,y) ))
         lc.append( COLOR[1] )

      x = []
      y = []
      z = []
      if settings.zi is not None: 
         for Zi in settings.zi: 
             l = float(Zi)/float(1+Zi)
             x.append( l*C[0]-l*B[0]+B[0] ) 
             y.append( l*C[1]-l*B[1]+B[1] )
             x.append( A[0] )
             y.append( A[1] )
             z.append( 0 ) 
             z.append( 0 ) 
             ls.append( LS ) 
         if For3D :
            seg.append(list( zip(x,y,z) ))
         else:
            seg.append(list( zip(x,y) ))
         lc.append( COLOR[2] )
      lines=seg

   if For3D:
      line_segments = Line3DCollection( lines, linewidths = LW , colors=lc, \
         linestyles=ls ,zorder=ZORDER)
   else:
      line_segments = LineCollection( lines, linewidths = LW , colors=lc, \
         linestyles=ls ,zorder=ZORDER)

   return line_segments

#------------------------------------------------------------------------------
# creates coordinates for grid of ternary 
#------------------------------------------------------------------------------
def GetTieLineCollection(sx,For3D=False,debug=None):
   from numpy import append as np_append
   from numpy import array as np_array
   from numpy.random import rand as np_random_rand
   from matplotlib.collections import LineCollection
   from matplotlib.collections import PatchCollection
   from matplotlib.patches import Polygon


   #set some default values: 
   if sx.lw is not None: 
      LW = sx.lw
   else :
      LW = '1'

   if sx.tl_color is not None: 
      COLOR = sx.tl_color
   else :
      COLOR = 'black'

   ZORDER = 8

   if sx.ls is not None: 
      LS = sx.ls
   else :
      LS = 'solid'

   lines=[]
   X=[]
   if debug : 
      patches = []  
      k=0
   for j in range(0,len(sx.hull.inner_simplices)):
      if debug: 
         k+=1
      x=[];y=[]
      for i in [0,1,2]: 
         xi=sx.hull.points[sx.hull.inner_simplices[j],0][i]
         yi=sx.hull.points[sx.hull.inner_simplices[j],1][i]
         x=np_append(x,xi)
         y=np_append(y,yi)
      #append last point
      xi=sx.hull.points[sx.hull.inner_simplices[j],0][0]
      yi=sx.hull.points[sx.hull.inner_simplices[j],1][0]
      x=np_append(x,xi)
      y=np_append(y,yi)

      lines.append(list(zip(x,y))) 

      if debug : 
         polygon =\
         Polygon([np_array([sx.hull.points[sx.hull.inner_simplices[j],0][0],\
                            sx.hull.points[sx.hull.inner_simplices[j],1][0]]),\
                  np_array([sx.hull.points[sx.hull.inner_simplices[j],0][1],\
                            sx.hull.points[sx.hull.inner_simplices[j],1][1]]),\
                  np_array([sx.hull.points[sx.hull.inner_simplices[j],0][2],\
                            sx.hull.points[sx.hull.inner_simplices[j],1][2]]),\
                  np_array([sx.hull.points[sx.hull.inner_simplices[j],0][0],\
                            sx.hull.points[sx.hull.inner_simplices[j],1][0]])],
         True )
         patches.append(polygon)
   if debug: 
      line_segments  = PatchCollection(patches, alpha=1.0, edgecolor = 'none')
      line_segments.set_array(np_random_rand(k))
   else:
      line_segments = LineCollection( lines, linewidths = LW , colors=COLOR, \
      linestyles=LS ,zorder=ZORDER)

   return line_segments 

###############################################################################
# subclassed from src in: /usr/local/lib/python2.7/dist-packages/matplotlib/
class CustomNavigationToolbar2Wx(NavigationToolbar2Wx):
    def __init__(self, canvas,parent,Project):
        from matplotlib.backend_bases import NavigationToolbar2
        wx.ToolBar.__init__(self, canvas.GetParent(), -1)
        NavigationToolbar2.__init__(self, canvas)
        self.canvas = canvas
        self.Project = Project
        self._idle = True
        self.statbar = None
        self.logger = self.GetTopLevelParent().C.PageLog
        self.statusbar = self.GetTopLevelParent().statusbar
        self.prevZoomRect = None
        self.LineToolActive = False
        self.Backgrounds = None
        # for now, use alternate zoom-rectangle drawing on all
        # Macs. N.B. In future versions of wx it may be possible to
        # detect Retina displays with window.GetContentScaleFactor()
        # and/or dc.GetContentScaleFactor()
        self.retinaFix = 'wxMac' in wx.PlatformInfo
        #
        # cids for Snap Tool 
        #
        self.Snapcids = [ ]
        #
        # collected coordinates cids for LineTool
        #
        self.LineToolcids = [ ]
        self.AddCoordinatesOn = False   # used to add coordinates when Line tool selected
        #
        # nbind coordinates info in statubar on init
        #
        self.Initcids = [ ] 
        self.Initcids.append(canvas.mpl_connect('motion_notify_event',\
           self.OnInitMotion))
        #
        # cursor point for Snap and Line tool 
        #
        self.CursorPoint = None
        self.LineMoving = None

    def _init_toolbar(self):
        self._parent = self.canvas.GetParent()
        #
        # store project data 
        #
        self.Project = self._parent.Project
        #
        # initialize Tool ID dictionary 
        #
        self.wx_ids = {}
        self.button_is_enabled = {}
        # list of toolitems to add to the toolbar, format is:
        # (
        #   text, # the text of the button (often not visible to users)
        #   tooltip_text, # the tooltip shown on hover (where possible)
        #   icon, # name of the image for the button (without the extension)
        #   name_of_method, # name of the method in NavigationToolbar2 to call
        # )
        self.items=(
        ('Home', 'Home', 'Reset original view', 'home'), \
        ('Back', 'Back', 'Back to previous view', 'back'), \
        ('Forward', 'Forward', 'Forward to next view', 'forward'), \
        ('Pan', 'Pan', 'Pan axes with left mouse, zoom with right', 'pan'), \
        ('Zoom', 'Zoom', 'Zoom in (out) to rectangle with left (right) mouse', 'zoom'), \
        ('Opt', 'Opt', 'Set optimum layout', 'OnOpt'), \
        (None, None, None, None), \
        ('Points', 'DataPoints', 'Show data points', 'OnPoints'), \
        ('TieLines', 'DataPoints2', 'Show convex hull', 'OnTieLines'), \
        ('Grid', 'Grid', 'Show grid', 'OnGrid'), \
        ('Binding', 'BindingEnergy', 'Show binding energy', 'OnBinding'), \
        ('Stability', 'Stability', 'Show stability', 'OnStability'), \
        ('Reaction', 'Reaction', 'Show reaction stability', 'OnReaction'), \
        (None, None, None, None), \
        ('LineTool', 'LineTool', 'Start line tool', 'OnLineTool'), \
        ('Snap', 'SnapCursor', 'Snap cursor to data points', 'OnSnap'), \
        ('Toggle3D', 'Toggle3D', 'Show diagram in 3D', 'OnToggle3D'), \
        )
        self.tools=[]
        for text, icon, tooltip_text, callback in self.items:
            if text is None:
                self.AddSeparator()
                continue
            self.tools.append( text )
            self.wx_ids[text] = wx.NewId()

            # buttons are enabled by default
            self.button_is_enabled[ text ] = True

            img=getattr(ico,icon).GetImage().ConvertToBitmap()
            self._AddTool(self, self.wx_ids, text,
                        img,
                        tooltip_text)
            if hasattr(self,callback):
               self.Bind(wx.EVT_TOOL, getattr(self, callback),
                      id=self.wx_ids[text])

        self.Realize()
        # 
        # after toolbar set up, plot project  
        # 
        self.InitPlot() 

    def _AddTool(self,parent, wx_ids, text, bmp, tooltip_text):
        if text in ['Pan', 'Zoom','Opt','Points','TieLines','Grid',\
            'Snap','Binding','Stability','Reaction','Snap','LineTool','Toggle3D']:
            parent.AddCheckTool(
                wx_ids[text],
                bmp,
                shortHelp=tooltip_text,
                longHelp=tooltip_text)
        else:
            parent.AddSimpleTool(wx_ids[text], bmp, tooltip_text, tooltip_text)
        #
        # store tool id here 
        #
        setattr(self,text+'ID',wx_ids[text])

    # -------------------------------------------------------------------------
    def OnPoints(self, *args):
       a = self.canvas.figure.get_axes()[0]
       #
       # switch status of visibility
       #
       self.Project.PointsScatterOn = not self.Project.PointsScatterOn 
       if self.Project.Toggle3DOn :
          # make sure collection is loaded 
          if self.Project.PointsScatter3D is None: 
             self.InitPoints(a)
          self.Project.PointsScatter3D.set_visible(self.Project.PointsScatterOn)
       else:
          # make sure collection is loaded 
          if self.Project.PointsScatter is None: 
             self.InitPoints(a)
          self.Project.PointsScatter.set_visible(self.Project.PointsScatterOn)
       #
       # Switch Color of data points depending on on or off
       #
       if self.Project.PointsScatterOn:
          obj = self.GetCurrentContourObject() 
          if len(obj)>0:
             self.ChangePointColorMap()
       #
       # replot 
       #
       self.canvas.draw()


    # -------------------------------------------------------------------------
    def OnTieLines(self, *args):
       a = self.canvas.figure.get_axes()[0]
       #
       # switch status of visibility
       #
       if self.Project.TieLinesCollectionOn:
          if self.Project.Toggle3DOn :
             for c in self.Project.TieLinesCollection3D:
                 c.remove()
          else:
             for c in self.Project.TieLinesCollection:
                 c.remove()
       else:
          if self.Project.Toggle3DOn:
             for c in self.Project.TieLinesCollection3D:
                 a.add_collection( c )
          else:
             for c in self.Project.TieLinesCollection:
                 a.add_collection( c )
       #
       # switch status of visibility
       #
       self.Project.TieLinesCollectionOn = not self.Project.TieLinesCollectionOn
       #
       # replot 
       #
       self.canvas.draw()
 
    # -------------------------------------------------------------------------
    def OnGrid(self, *args):
       a = self.canvas.figure.get_axes()[0]
       #
       # switch status of visibility
       #
       if self.Project.GridLinesCollectionOn:

          if self.Project.Toggle3DOn :
             for c in self.Project.GridLinesCollection3D:
                 c.remove()
          else:
             for c in self.Project.GridLinesCollection:
                 c.remove()

       else:
          if self.Project.Toggle3DOn :
             for c in self.Project.GridLinesCollection3D:
                 a.add_collection( c )
          else:
             for c in self.Project.GridLinesCollection:
                 a.add_collection( c )
       #
       # switch status of visibility
       #
       self.Project.GridLinesCollectionOn = \
          not self.Project.GridLinesCollectionOn
       #
       # replot 
       #
       self.canvas.draw()

    # -------------------------------------------------------------------------
    def OnBinding(self, *args):
       a = self.canvas.figure.get_axes()[0]
       #
       # first, get rid of all other Contoursets
       #
       if self.Project.StabilityContourSetOn:
          self.OnStability()
       if self.Project.ReactionContourSetOn:
          self.OnReaction()
       #
       # switch status of visibility
       #
       obj = 'BindingContourSet'
       if self.Project.Toggle3DOn :
          obj +='3D'
       if self.Project.BindingContourSetOn:
          self.RemoveContourSet(obj)
       else:
          self.AddContourSet(obj)
       #
       # switch status of visibility
       #
       self.Project.BindingContourSetOn = not self.Project.BindingContourSetOn
       #
       # Switch Color of data points depending on on or off
       #
       if self.Project.PointsScatterOn:
          self.ChangePointColorMap(uniform=not self.Project.BindingContourSetOn)
       #
       # instead of replotting only, set optimal view too
       #
       self.OnOpt()

    # -------------------------------------------------------------------------
    def OnStability(self, *args):
       a = self.canvas.figure.get_axes()[0]
       #
       # first, get rid of all other Contoursets
       #
       if self.Project.BindingContourSetOn:
          self.OnBinding()
       if self.Project.ReactionContourSetOn:
          self.OnReaction()
       #
       # switch status of visibility
       #
       obj = 'StabilityContourSet'
       if self.Project.Toggle3DOn :
          obj +='3D'

       if self.Project.StabilityContourSetOn:
          self.RemoveContourSet(obj)
       else:
          self.AddContourSet(obj)
       #
       # switch status of visibility
       #
       self.Project.StabilityContourSetOn = not self.Project.StabilityContourSetOn
       #
       # Switch Color of data points depending on on or off
       #
       if self.Project.PointsScatterOn:
          self.ChangePointColorMap(uniform=not self.Project.StabilityContourSetOn)
       #
       # instead of replotting only, set optimal view too
       #
       self.OnOpt()

    def OnReaction(self, *args):
       a = self.canvas.figure.get_axes()[0]
       #
       # first, get rid of all other Contoursets
       #
       if self.Project.BindingContourSetOn:
          self.OnBinding()
       if self.Project.StabilityContourSetOn:
          self.OnStability()
       #
       # switch status of visibility
       #
       obj = 'ReactionContourSet'
       if self.Project.Toggle3DOn :
          obj +='3D'
       if self.Project.ReactionContourSetOn:
          self.RemoveContourSet(obj)
       else:
          self.AddContourSet(obj)
       #
       # switch status of visibility
       #
       self.Project.ReactionContourSetOn = not self.Project.ReactionContourSetOn
       #
       # Switch Color of data points depending on on or off
       #
       if self.Project.PointsScatterOn:
          self.ChangePointColorMap(uniform=not self.Project.ReactionContourSetOn)
       #
       # instead of replotting only, set optimal view too
       #
       self.OnOpt()

    # -------------------------------------------------------------------------
    # tools 
    def OnSnap(self, *args):
       #
       # in case hook is present, disconnect
       #
       if self.Project.SnapOn: 
          for  cid in self.Snapcids:
             self.canvas.mpl_disconnect(cid)       
          self.Snapcids = [ ]
       else:
          #
          # capture figure background, so that we dont need to replot everything
          #
          self.SwitchOnVerboseMode(debug=True)
          #
          # bind functions
          #
          self.Snapcids.append(self.canvas.mpl_connect(
             'button_press_event', self.OnSnapClick) )
          self.Snapcids.append(self.canvas.mpl_connect(
             'motion_notify_event', self.OnSnapMove ))
          #
          # add a cursor to the axis 
          #
          a = self.canvas.figure.get_axes()[0]
          if self.CursorPoint is None : 
             self.CursorPoint=a.scatter(0,0,s=self.Project.ps, color='gold',
                                              zorder=11,
                                              alpha=1.0)
          else:
             self.CursorPoint.set_color('gold')
          self.CursorPoint.set_visible(False)

       self.Project.SnapOn = not self.Project.SnapOn

    def OnLineTool(self, *args):
       #
       # switch cursor, depending on status of tool
       #
       self.Project.LineToolOn = not self.Project.LineToolOn
       self.canvas.Bind(wx.EVT_ENTER_WINDOW, self.ChangeCursor)
       #
       # start collecting datapoints 
       #
       if self.Project.LineToolOn :
          self.LineToolcids = [ ] 
          text="Line tool: press ENTER or \"q\" when finished,"+\
                "ESC to delete last point"
          self.statusbar.SetStatusText(text,0)
          #
          # clean coordinates
          # 
          self.Project.PickedCoordinates = [ ] 
          self.LineToolcids.append(\
             self.canvas.mpl_connect('button_press_event', self.OnLineToolClick) )
          self.LineToolcids.append(\
             self.canvas.mpl_connect('motion_notify_event', self.OnLineToolMove) )
          self.LineToolcids.append(\
             self.canvas.mpl_connect('key_press_event', self.OnKeyPressed) )
          #
          # capture figure background, so that we dont need to replot everything
          #
          self.SwitchOnVerboseMode(debug=True)
          #
          # add a cursor to the axis 
          #
          a = self.canvas.figure.get_axes()[0]
          if self.CursorPoint is None : 
             self.CursorPoint=a.scatter(0,0,s=self.Project.ps, color='red',
                                           zorder=11,
                                           alpha=1.0)
          else:
             self.CursorPoint.set_color('red')
          self.CursorPoint.set_visible(False)
       else:
          l = len( self.Project.PickedCoordinates )
          for i in range( l-1, -1,-1):
             del self.Project.PickedCoordinates[i]
             self.Project.LineToolPointCollection[i].remove()
             del self.Project.LineToolPointCollection[i]
          for i in range( l-2, -1,-1):
             self.Project.LineToolLinesCollection[i].remove()
             del self.Project.LineToolLinesCollection[i]

          if len( self.Project.PickedCoordinates ) == 0 : 
             self.Project.PickedCoordinates = [ ] 
          for  cid in self.LineToolcids:
             self.canvas.mpl_disconnect(cid)       
          self.LineToolcids = [ ]
          self.canvas.draw()

    #--------------------------------------------------------------------------
    # this captures the currently plotted background, disables most of the
    # buttons which can change the background. 
    # It should be used when OnSnap or OnLineTool are chosen, i.e. preparation
    # to blit figure. 
    #--------------------------------------------------------------------------
    def SwitchOnVerboseMode(self, debug=None):
       # Let's capture the background of the figure
       axes = self.canvas.figure.get_axes()
       self.Backgrounds = [ self.canvas.copy_from_bbox(ax.bbox) for ax in axes ]

    def ChangeCursor(self,event):
       from matplotlib.backends.backend_wx import  wxc
       if self.Project.LineToolOn :
           self.canvas.SetCursor(wxc.StockCursor(wx.CURSOR_CROSS))
       else:
           self.canvas.SetCursor(wxc.StockCursor(wx.CURSOR_ARROW))

    def OnKeyPressed(self,event):
       from dialog import ShowStatusErrorMessage
       from MyTools import LineToolForm
       from MyTools2D import LineTool2DForm 
       l = len( self.Project.PickedCoordinates )
       if event.key == 'escape' : 
          if l>0 : 
             del self.Project.PickedCoordinates[l-1]
             self.Project.LineToolPointCollection[l-1].remove()
             del self.Project.LineToolPointCollection[l-1]
          if l>1:
             self.Project.LineToolLinesCollection[l-2].remove()
             del self.Project.LineToolLinesCollection[l-2]
          self.canvas.draw()
          #
          # capture figure background, 
          # so that we dont need to replot everything
          #
          self.SwitchOnVerboseMode(debug=True)

          if len( self.Project.PickedCoordinates ) == 0 : 
             self.Project.PickedCoordinates = [ ] 

       elif event.key == 'enter' or event.key == 'q': 
          if l > 1 : 
             if self.Load2DLineToolForm():
                dlg = LineTool2DForm(self.GetTopLevelParent(), 
                style=wx.DEFAULT_FRAME_STYLE)
             else:
                dlg = LineToolForm(self.GetTopLevelParent(), 
                style=wx.DEFAULT_FRAME_STYLE)

             dlg.ShowModal()
             dlg.Destroy() 
          else: 
             text = 'At least two points are need'
             ShowStatusErrorMessage( self.statusbar, text )
             
    def OnToggle3D(self, *args):
       from mpl_toolkits.mplot3d import Axes3D

       self.Project.Toggle3DOn = not self.Project.Toggle3DOn 
       #
       # clear figure and initialize plot 
       #
       self.ClearFigure()
       self.InitPlot()


    def OnOpt(self, *args):
       #
       # set opt z-lim and view in 3D
       #
       if self.Project.Toggle3DOn:
          a = self.canvas.figure.get_axes()[0]
          Z = self.GetCurrentZData()
          if Z is not None:
             a.set_zlim(min(Z)*1.1,max(Z)*1.1)
       # 
       # remove Color bar if any of the contours is set visible
       # 
       obj = self.GetCurrentContourObject() 
       if len(obj)>0: 
          self.RemoveContourSet(obj) 
          self.AddContourSet(obj) 
       else: 
          self.SetAxesToOpt()
          
       #
       # replot 
       #
       self.canvas.draw()
    # -------------------------------------------------------------------------
    # settings panel 
    def OnLabels(self, *args):
       a = self.canvas.figure.get_axes()[0]
       if self.Project.LabelsOn :
          self.Project.LabelsOn = False
          self.Project.LabelAHandle.remove()
          self.Project.LabelBHandle.remove()
          self.Project.LabelCHandle.remove()
          del self.Project.LabelAHandle
          del self.Project.LabelBHandle
          del self.Project.LabelCHandle
          self.Project.LabelAHandle = None
          self.Project.LabelBHandle = None
          self.Project.LabelCHandle = None
       else:
          self.Project.LabelsOn = True
          self.InitLabeling(a)
       #
       # replot 
       #
       self.canvas.draw()

    # --------------------------------------------------------------------------
    # Helper routines to add and remove colorbar 
    # -------------------------------------------------------------------------
    def InitPlot(self):
        """ initializes Plot, should be called only once from _init_toolbar  """
        from matplotlib.pyplot import imread 
        from io import BytesIO
        from mpl_toolkits.mplot3d import Axes3D
        import matplotlib.gridspec as gridspec
        # 
        # Plot an empty page
        # 
        if self.Project is None : 
            a=self.canvas.figure.add_subplot(111)
            #
            # display Logo 
            # 
            img=getattr(ico,'Acronym').GetImage()
            bio = BytesIO()
            bs = wx.OutputStream(bio)
            img.SaveStream(bs, wx.BITMAP_TYPE_PNG)
            bio.seek(0) # rewind stream
            #a.set_frame_on(False)
            self.canvas.figure.patch.set_facecolor('white')
            a.get_xaxis().tick_bottom()
            a.get_yaxis().tick_left()
            a.get_yaxis().tick_right()
            a.axes.get_yaxis().set_visible(False)
            a.axes.get_xaxis().set_visible(False)
            #
            # save instance for empty page 
            #
            self.Project.EmptyCollection = a.imshow(imread(bio))
            self.Project.EmptyCollectionOn = True 
        # 
        # Plot diagram and its collections
        # 
        else:
            #
            # check diagram type
            #
            if self.Project.Toggle3DOn: 
               #
               # add 3D subplot axes and set aspect ratio to 1 
               #
               a=self.canvas.figure.add_subplot(111,projection='3d')
               a.set_aspect('equal')
            else:
               #
               # add 2D subplot axes and set aspect ratio to 1 
               #
               a=self.canvas.figure.add_subplot(111)
               a.set_aspect('equal','datalim')
               #
               # change background to white and remove frame 
               #
               a.set_frame_on(False)
               a.get_xaxis().tick_bottom()
               a.get_yaxis().tick_left()
               a.get_yaxis().tick_right()
               a.axes.get_yaxis().set_visible(False)
               a.axes.get_xaxis().set_visible(False)
#DEBUG:
#               a.set_frame_on(True)
#               a.axhline(y=0.920025404, xmin=0, xmax=1,lw=2,c='red')
#               a.axhline(y=0.866025404, xmin=0, xmax=1)
#               a.axhline(y=0.839025404, xmin=0, xmax=1)
#               a.axhline(y=-0.050000000, xmin=0, xmax=1,lw=2,c='red')
#DEBUG:
#            self.canvas.figure.patch.set_facecolor(self.Project.facecolor)
            #
            # adjust subplot  
            #
            self.canvas.figure.subplots_adjust(bottom=0.0,left=0.0,top=1,right=0.80)
            #
            # Plot Border
            #
            self.InitBorderLines(a)
            #
            # Plot points 
            #
            self.InitPoints(a)
            #
            # Plot convex hull tie lines 
            #
            self.InitTieLines(a)
            #
            # Plot grid lines 
            #
            self.InitGridLines(a)
            #
            # Plot binding energy with colormap
            #
            self.InitBindingContour(a)
            #
            # Plot stability energy with colormap
            #
            self.InitStabilityContour(a)
            #
            # Plot binding energy with colormap
            #
            self.InitReactionContour(a)
            #
            # set the correct ContourSet to visible
            #
            if len(self.GetCurrentContourObject())>0:
               self.AddContourSet(self.GetCurrentContourObject())
            #
            # Reset all toolbar buttons to off
            #
            tb = self._parent.parent.toolbar
            color= wx.SystemSettings.GetColour(wx.SYS_COLOUR_BACKGROUND)
            for text,icon,tooltip_text,callback in self.items:
               if text is None:
                  continue
               t=getattr(tb,text)
               t.SetBackgroundColour( color )
               t.Refresh()
            #
            # initialize labeling
            #
            self.InitLabeling(a)

            #
            # replot points
            #
            self.OnPoints()
            self.OnPoints()

        self.canvas.draw()

    #--------------------------------------------------------------------------
    def InitBorderLines(self,a):
        from mpl_toolkits.mplot3d.art3d import Poly3DCollection        
        from matplotlib.colors import colorConverter
        from numpy import zeros 
        #
        # plot Borderlines, decide weather to use 3D or 2D style 
        #
        if self.Project.Toggle3DOn:
           #
           # convert tie line colors to rgb tuple
           #
           c =  colorConverter.to_rgb( self.Project.tl_color)
           c = c + (self.Project.alpha,)
           #
           # plot border lines as polygon edges 
           #
           self.Project.BorderLinesCollection3D = [ ] 
           t = GetBorderCollection( self.Project )
           a.add_collection3d( t )
           self.Project.BorderCollection3D.append( t ) 

        else:
           self.Project.BorderLinesCollection = [ ] 
           t = GetBorderCollection( self.Project )
           a.add_collection(t)
           self.Project.BorderLinesCollection.append( t ) 
        #
        # remove them if they should not be shown intially 
        #
        if not self.Project.BorderCollectionOn:
           if self.Project.Toggle3DOn :
              for c in self.Project.BorderCollection3D:
                  c.remove()
           else:
              for c in self.Project.BorderCollection:
                  c.remove()

    #--------------------------------------------------------------------------
    def InitPoints(self,a):
        """ Initializes Points object """
        #
        # set data points 
        #
        X = self.Project.hull.points[:,0]
        Y = self.Project.hull.points[:,1]
        #
        # plot points, decide weather to use 3D or 2D style 
        #
        if self.Project.Toggle3DOn:
           Z = self.Project.hull.points[:,2]
           Cmd = self.Project.Cmd
           self.Project.PointsScatter3D=a.scatter(X,Y,Z, c=Z,cmap=Cmd,\
              s=self.Project.ps, zorder=9, edgecolors=self.Project.pt_color )
              #s=self.Project.ps, zorder=9, edgecolors='black' )
           #
           # remove them if they should not be shown intially 
           #
           self.Project.PointsScatter3D.set_visible(self.Project.PointsScatterOn)
        else: 
           Cmd = None
           Z = self.Project.pt_color
           self.Project.PointsScatter=a.scatter(X,Y,c=Z, s=self.Project.ps, \
             cmap=Cmd, zorder=9, edgecolors=Z )
           #  cmap=Cmd, zorder=9, edgecolors='black' )
           #
           # remove them if they should not be shown intially 
           #
           self.Project.PointsScatter.set_visible(self.Project.PointsScatterOn)
        
      
    #--------------------------------------------------------------------------
    def InitTieLines(self,a):
        from mpl_toolkits.mplot3d.art3d import Poly3DCollection        
        from matplotlib.colors import colorConverter
        from numpy import zeros 
        #
        # plot Tielines, decide weather to use 3D or 2D style 
        #
        if self.Project.Toggle3DOn:
           #
           # convert tie line colors to rgb tuple
           #
           c =  colorConverter.to_rgb( self.Project.tl_color)
           c = c + (self.Project.alpha,)
           #
           # plot tie lines as polygon edges 
           #
           self.Project.TieLinesCollection3D = [ ] 
           for i in range(0,len(self.Project.hull.inner_simplices)):
              v=self.Project.hull.inner_simplices[ i ]
              simplex=self.Project.hull.points[ v ]
              vertices=[zip( simplex[:,0],simplex[:,1] ,simplex[:,2] )]
              #z = zeros(len(self.Project.hull.points))
              #vertices=[zip( simplex[:,0],simplex[:,1] ,z )]
              face = Poly3DCollection(vertices, 
                  alpha=self.Project.alpha, linewidth=self.Project.lw, edgecolors='black')
              face.set_facecolor(c)
              t = a.add_collection3d(face)
              self.Project.TieLinesCollection3D.append( face ) 
        else:
           self.Project.TieLinesCollection = [ ] 
           t = GetTieLineCollection( self.Project )
           a.add_collection(t)
           self.Project.TieLinesCollection.append( t ) 

        #
        # remove them if they should not be shown intially 
        #
        if not self.Project.TieLinesCollectionOn:
           if self.Project.Toggle3DOn :
              for c in self.Project.TieLinesCollection3D:
                  c.remove()
           else:
              for c in self.Project.TieLinesCollection:
                  c.remove()
                   
    #--------------------------------------------------------------------------
    def InitGridLines(self,a):
        from mpl_toolkits.mplot3d.art3d import Poly3DCollection        
        from matplotlib.colors import colorConverter
        #
        # plot gridlines, decide weather to use 3D or 2D style 
        #
        if self.Project.Toggle3DOn: 
           self.Project.GridLinesCollection3D = [ ] 
           g = GetGridCollection(self.Project,For3D=True) 
           a.add_collection( g )
           self.Project.GridLinesCollection3D.append( g ) 
        else:
           self.Project.GridLinesCollection = [ ] 
           g = GetGridCollection(self.Project) 
           a.add_collection( g )
           self.Project.GridLinesCollection.append( g ) 
        #
        # remove them if they should not be shown intially 
        #
        if not self.Project.GridLinesCollectionOn:

           if self.Project.Toggle3DOn :
              for c in self.Project.GridLinesCollection3D:
                  c.remove()
           else:
              for c in self.Project.GridLinesCollection:
                  c.remove()

    #--------------------------------------------------------------------------
    def InitBindingContour(self,a):
        from mpl_toolkits.mplot3d.art3d import Poly3DCollection        
        from matplotlib.colors import colorConverter
        from numpy import linspace
        #
        # use Trianglulation 
        #
        status=self.GetTopLevelParent().statusbar
        T, x,y,z = GetBindingContour(self.Project,status=status) 
        lev = linspace(z.min(),z.max(),self.Project.CmdIntersect)
        #
        # plot contour, decide weather to use 3D or 2D style 
        #
        if self.Project.Toggle3DOn:
           #self.Project.BindingContourSet3D = a.tricontourf(x,y,T.triangles,z,N=20, \
           #   levels=lev, zdir = 'z', offset = 1.2*z.min(),cmap=self.Project.Cmd)
           self.Project.BindingContourSet3D = a.plot_trisurf(x,y,T.triangles,z, \
              cmap=self.Project.Cmd, lw=0.1)
           b=self.Project.BindingContourSet3D
        else:
           self.Project.BindingContourSet = a.tricontourf(x,y,T.triangles,z,N=20, \
              levels=lev,cmap=self.Project.Cmd)
           b=self.Project.BindingContourSet

#DEBUG : plot triangles 
#           tri=self.Project.hull.tri
#           a.triplot(tri.points[:,0], tri.points[:,1], tri.simplices.copy())
#DEUUB

        #
        # unset visibility, since initialized
        #
        b.set_alpha(0.0)

    #--------------------------------------------------------------------------
    def InitStabilityContour(self,a):
        from mpl_toolkits.mplot3d.art3d import Poly3DCollection        
        from matplotlib.colors import colorConverter
        from numpy import linspace
        #
        # use Trianglulation 
        #
        status=self.GetTopLevelParent().statusbar
        T, x,y,z = GetStabilityContour(self.Project,status=status) 
        lev = linspace(z.min(),z.max(),self.Project.CmdIntersect)
        #
        # plot contour, decide weather to use 3D or 2D style 
        #
        if self.Project.Toggle3DOn:
           #self.Project.StabilityContourSet3D = a.tricontourf(x,y,T.triangles,z,N=20, \
           #   levels=lev, zdir = 'z', offset = 1.2*z.min(),cmap=self.Project.Cmd)
           self.Project.StabilityContourSet3D = a.plot_trisurf(x,y,T.triangles,z,\
               cmap=self.Project.Cmd, lw=0.1)
           b=self.Project.StabilityContourSet3D
        else:
           self.Project.StabilityContourSet = a.tricontourf(x,y,T.triangles,z,N=20, \
              levels=lev,cmap=self.Project.Cmd)
           b=self.Project.StabilityContourSet
        #
        # unset visibility, since initialized
        #
        b.set_alpha(0.0)

    #--------------------------------------------------------------------------
    def InitReactionContour(self,a):
        from mpl_toolkits.mplot3d.art3d import Poly3DCollection        
        from matplotlib.colors import colorConverter
        from numpy import linspace
        #
        # use Trianglulation 
        #
        status=self.GetTopLevelParent().statusbar
        T, x,y,z = GetReactionContour(self.Project,status=status) 
        lev = linspace(z.min(),z.max(),self.Project.CmdIntersect)
        #
        # plot contour, decide weather to use 3D or 2D style 
        #
        if self.Project.Toggle3DOn:
           #self.Project.ReactionContourSet3D = a.tricontourf(x,y,T.triangles,z,N=20, \
           #   levels=lev, zdir = 'z', offset = 1.2*z.min(),cmap=self.Project.Cmd)
           self.Project.ReactionContourSet3D = a.plot_trisurf(x,y,z ,\
              cmap=self.Project.Cmd,lw=0.1)
           b=self.Project.ReactionContourSet3D
        else:
           self.Project.ReactionContourSet = a.tricontourf(x,y,T.triangles,z,N=20, \
              levels=lev,cmap=self.Project.Cmd)
           b=self.Project.ReactionContourSet
        #
        # unset visibility, since initialized
        #
        b.set_alpha(0.0)

    #--------------------------------------------------------------------------
    def InitLabeling(self,a):
       """ This routine intializes labeling of endpoints and ticks """
       sqrt32 =0.866025404
       items = ( ('LabelA','C','LabelAOffset', (0,0), 'LabelAHandle') ,\
                 ('LabelB','A','LabelBOffset', (1,0), 'LabelBHandle'),\
                 ('LabelC','B','LabelCOffset', (0.5,sqrt32),'LabelCHandle'))
       for label , phase, offset, pos, handle in items:
          #
          # destroy handle if its already initialized
          #
          if getattr(self.Project,handle) is not None:
              l= getattr(self.Project,handle)
              l.remove()
              del l 
              setattr(self.Project,handle,None)

          #
          # set endpoint labels, their offsets and color 
          #
          l = getattr(self.Project,label)
          o = getattr(self.Project,offset)
          c = (pos[0] - o[0], pos[1]-o[1])
          color = self.Project.LabelColor[phase]
          #
          # initialize them to standard labels, if they are not set yet
          #
          if l is None:
             text = self.Project.pdb[phase]['LaTeXName']
          else:
             text = l 
          #
          # store  handle 
          #
          if self.Project.Toggle3DOn:
             l = a.text(c[0],c[1],0, text, color = color , \
                fontsize=self.Project.fontsize)
          else:
             l = a.text(c[0],c[1], text, color = color , \
                fontsize=self.Project.fontsize)
          setattr( self.Project,handle,l )
          #
          # switch off labels optionally
          #
          if not self.Project.LabelsOn:
             getattr(self.Project,handle).set_text(' ')

    # -------------------------------------------------------------------------
    def ClearFigure(self):
       axes = self.canvas.figure.get_axes()
       I = len(axes)
       for i in range( 0, I-1 ):
          self.canvas.figure.delaxes( axes[I-i-1] ) 
       self.canvas.figure.delaxes( axes[0] ) 
       # also clear figure, just to be save 
       #
       self.canvas.figure.clear()
       #
       # replot
       #
       #self.canvas.draw()

    def CreateColorBar(self,b, pos=None):
       """ create an axes on the right side of ax. The width of cax will be 5%
           of ax and the padding between cax and ax will be fixed at 0.05 inch. """
       import matplotlib.pyplot as plt
       #
       # TODO: calculate optimal size 
       #       l      b     w      h 
       if pos is not None:
           box = pos
       else:
           box = [0.80,  0.11, 0.03, 0.840]
       self.Project.CBAxes = self.canvas.figure.add_axes(box ) 
       self.Project.ColorBar = plt.colorbar(b, cax = self.Project.CBAxes,\
          format='%1.2f')  

       self.Project.CBAxes.set_ylabel(self.Project.CBLabel, \
            fontsize=self.Project.fontsize)
       self.Project.CBAxes.yaxis.tick_right()
       self.Project.CBAxes.yaxis.set_label_position("right")
       self.Project.CBAxes.tick_params(labelsize=self.Project.fontsize)

    def RemoveContourSet(self,obj):
       """ Call this to remove a contour plot object and its colorbar """
       getattr(self.Project,obj).set_alpha(0.0)
       if self.Project.ColorBar is not None :
          self.Project.ColorBar.remove()
          I = len(self.canvas.figure.axes) 
          for i in range( 0, I-1 ):
              self.canvas.figure.delaxes( self.canvas.figure.axes[I-i-1] ) 
          del self.Project.ColorBar
          self.Project.ColorBar = None
    
    def AddContourSet(self,obj):
       """ Call this to add a contour plot object and its colorbar """
       b =  getattr(self.Project,obj)
       b.set_alpha(1.0)
       self.SetAxesToOpt()
       boxsize = self.GetOptCBBoxSize()
       self.CreateColorBar( b, pos = boxsize ) 

    def ChangePointColorMap(self,uniform=False):
       from matplotlib.colors import Normalize
       from numpy import zeros 
       """ Call this function to switch from uniform to colormap for points and
           vice versa """
       #
       # check if its 3D
       #
       if self.Project.Toggle3DOn: 
          p =  self.Project.PointsScatter3D
       else:
          p =  self.Project.PointsScatter
       #
       # Switch Color of data points to uniform color of settings 
       #
       if uniform : 
          p.set_facecolor(self.Project.pt_color)
       else:   
          #
          # get current z data displayed
          #
          Z = self.GetCurrentZData()
          if Z is not None:
             if self.Project.Toggle3DOn:
                # 
                # here, we need to remove the handle first
                # 
                p.remove()
                X=self.Project.hull.points[:,0]
                Y=self.Project.hull.points[:,1]
                Cmd = self.Project.Cmd
                a = self.canvas.figure.get_axes()[0]
                self.Project.PointsScatter3D=a.scatter(X,Y,Z, c=Z,cmap=Cmd,\
                   s=self.Project.ps, zorder=9)
                p=self.Project.PointsScatter3D 
                a.set_zlim(min(Z)*1.1,max(Z)*1.1)
             p.set_array(Z)
             p.set_cmap(self.Project.Cmd)
             norm = Normalize(vmin=min(Z),vmax=max(Z))
             p.set_norm(norm)

    def GetCurrentContourObject(self):
       """ gets current displayed contour """
       sets = ('BindingContourSet','StabilityContourSet','ReactionContourSet')
       obj = ''
       for s in sets: 
          S = s+'On'
          if getattr(self.Project,S):
             obj = s
             if getattr(self.Project,'Toggle3DOn'):
                 obj += '3D' 
       
       return obj

    def GetCurrentZData(self,IncludeTieLines=False):
       """ gets current z data displayed """
       from numpy import append
       sets = (('BindingContourSetOn','points'),
               ('StabilityContourSetOn','distance'),
               ('ReactionContourSetOn','ReactionEnergies'))
       bind = self.Project.hull.points[:,2]
       dist = self.Project.hull.distance
       obj = ''
       for s, array in sets:
           if getattr(self.Project,s):
              obj = array

       dat = None
       if len(obj)>0:
          d= getattr(self.Project.hull,obj)
          if obj == 'points' : 
             dat = d[:,2]
          else: 
             dat = d

       if IncludeTieLines:
          if self.Project.TieLinesCollectionOn:
             if dat is None:
                dat = []
             
             dat = append(dat,min(self.Project.hull.points[:,2]))
             dat = append(dat,max(self.Project.hull.points[:,2]))
       return dat 

    def SetAxesToOpt(self):
       """ sets first axes to optimum """
       from numpy import abs as abs
       sqrt32 =0.866025404
       # 
       # Get main axes
       # 
       fig = self.canvas.figure
       a = fig.get_axes()[0]
       # 
       # in 3D just set top view and return
       # 
       if self.Project.Toggle3DOn:
          #
          # get current z data displayed
          #
          Z = self.GetCurrentZData(IncludeTieLines=True)
          if Z is not None:
             a.set_zlim(min(Z)*1.1,max(Z)*1.1)
          #
          # set view
          #
          a.view_init(45,-45)
          return 

       # 
       # this is the minimum distance to axes border in % for x and y axis 
       # 
       Dx0,Dy0 = 0.2,0.1
       # 
       # Get current axis distance in % limits
       # 
       xmin,xmax = a.get_xlim()
       ymin,ymax = a.get_ylim()
       # 
       # Get size of first axes (aw, ah) in pixel 
       # 
       bbox = a.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
       ratio = bbox.width/ bbox.height
       DyAct = (bbox.height/(ymax-ymin-sqrt32))/2.
       DxAct = (bbox.width/(xmax-xmin-1.))/2.
       #
       # in case width is larger than height
       #
       if ratio > 1 : 
          DxPer,DyPer = max(Dx0,DxAct), Dy0
       else:
          DxPer,DyPer = Dx0, max(Dy0,DyAct)
       #
       # calculate optimum x and y-range
       #
       xminOpt = -(xmax-xmin)*DxPer
       xmaxOpt = 1.-xminOpt
       yminOpt = -(ymax-ymin)*DyPer
       ymaxOpt = sqrt32-yminOpt
       #
       # set the data range
       #

       a.set_xlim( [xminOpt,xmaxOpt] )
       a.set_ylim( [yminOpt,ymaxOpt] )

       self.canvas.draw()

    def GetOptCBBoxSize(self):
       """ returns the size of the first axes """
       # 
       # default box size is: 
       # 
       l, b, w, h = 0.80,  0.11, 0.03, 0.720
       if self.Project.Toggle3DOn: 
          return  [l,b,w,h]

       sqrt32 =0.866025404
       # 
       # Get main axes
       # 
       fig = self.canvas.figure
       a = fig.get_axes()[0]
       # 
       # Get axis limits
       # 
       xmin,xmax = a.get_xlim()
       ymin,ymax = a.get_ylim()
#DEBUG:
#       print 'has range',ymin,ymax
#DEBUG:
       # 
       # Get size of first axes (aw, ah) in pixel 
       # 
#DEBUG:
       bbox = a.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
       aw, ah = bbox.width, bbox.height
#DEBUG
       #
       # optimum height
       #
       h = sqrt32 / (ymax - ymin)
       #
       # position at bottom 
       #
       b = (abs( ymax-ymin - sqrt32 ) / 2.)/(ymax-ymin)
#DEBUG:
#       print l, b, w, h 
#DEBUG:
       
       return [l, b, w, h ]

    # -------------------------------------------------------------------------
    # event bindings for mouse and Snap
    # -------------------------------------------------------------------------
    def OnSnapClick(self,e):
       a = e.inaxes
       if a is not None :
          fig = e.canvas.figure
          size = fig.get_size_inches()*fig.dpi
          #
          # pointer should be in right axes 
          #
          if e.x < size[0]*0.8 :
              #
              # also pointer should be inside triangle 
              #
              # if (ix,iy) is close to a point existing in the database
              Inside, ix,iy,key = IsInsideOrClose(e,self.Project)

              if Inside:
                  if e.button == 1 :
                     PrintReactionInfo(self.Project, ix, iy , self.logger, key) 
                  #
                  # in case of a double click open structure or edit file 
                  #
                  if e.dblclick :
                     if e.button == 1:
                         OpenStructure(self.Project,key,self.statusbar,\
                                       self.logger)
                     #
                     # right mouse button opens base file in editor 
                     #
                     elif e.button ==3 : 
                         EditStructureFiles(self.Project,key,\
                                            self.statusbar,self.logger)

    def OnSnapMove(self,e):
       from numpy import sin,pi
       a = e.inaxes
       if a is not None :
          fig = e.canvas.figure
          size = fig.get_size_inches()*fig.dpi
          #
          # pointer should be in right axes 
          #
          if e.x < size[0]*0.8 :
              #
              # also pointer should be inside triangle 
              #
              # if (ix,iy) is close to a point existing in the database
              Inside, ix,iy,key = IsInsideOrClose(e,self.Project)
   
              if Inside: 
                  #
                  # restore saved backgrounds
                  #
                  for background in self.Backgrounds:
                      self.canvas.restore_region(background)
                  #
                  # make visible and shift cursor
                  #
                  self.CursorPoint.set_visible(True)
                  self.CursorPoint.set_offsets((ix, iy))
                  #
                  # blit canvas, whatever this means. anyway its definitely
                  # faster compared to redrawing everthing on the canvas
                  #
                  k=0
                  for axes in self.canvas.figure.get_axes():
                     k+=1
                     if k==1:
                         #
                         # update point on first axis (cursor)
                         #
                         axes.draw_artist( self.CursorPoint )
                     self.canvas.blit(axes.bbox)
                  if e.button == 1 :
                      print 'pressedn in OnMove'
              else:
                  self.CursorPoint.set_visible(False)

              #self.canvas.draw()

    # -------------------------------------------------------------------------
    # event bindings for mouse and LineToolClick 
    # -------------------------------------------------------------------------
    def OnLineToolClick(self,e):
       from poscar import GetBarycentricCoord
       from matplotlib.collections import LineCollection
       from numpy import append as np_append 

       b = e.inaxes
       a = self.canvas.figure.get_axes()[0]
       if b is not None :
          fig = e.canvas.figure
          size = fig.get_size_inches()*fig.dpi
          #
          # pointer should be in right axes 
          #
          if e.x < size[0]*0.8 :
              #
              # also pointer should be inside triangle 
              #
              # if (ix,iy) is close to a point existing in the database
              Inside, ix,iy,key = IsInsideOrClose(e,self.Project)
              if Inside: 

                 # add two coordinate pairs, even without key pressed
                 if e.button == 1:
                    E=[ix,iy]
                    x0=GetBarycentricCoord( self.Project, ix,iy )
                    E.append(x0[0]);E.append(x0[1]);E.append(x0[2])
                    self.Project.PickedCoordinates.append( E )
                    self.Project.LineToolPointCollection.append(
                        a.scatter( ix, iy, c = 'red' , s = 100,zorder=9 ) )
                    l = len(self.Project.PickedCoordinates)
                    if l > 1 : 
                       x=[];y=[]
                       x=np_append(x,self.Project.PickedCoordinates[l-2][0])
                       x=np_append(x,self.Project.PickedCoordinates[l-1][0])
                       y=np_append(y,self.Project.PickedCoordinates[l-2][1])
                       y=np_append(y,self.Project.PickedCoordinates[l-1][1])
                       lines = []
                       lines.append( list(zip(x,y)))
                       line_segments = LineCollection( lines, color ='red',\
                          lw=2,zorder=9)

                       self.Project.LineToolLinesCollection.append(
                       a.add_collection(line_segments))
                       
                    self.canvas.draw()
                    #
                    # capture figure background, 
                    # so that we dont need to replot everything
                    #
                    self.SwitchOnVerboseMode(debug=True)

    def OnLineToolMove(self,e):
       from numpy import append as np_append 
       from matplotlib.collections import LineCollection
       b = e.inaxes
       a = self.canvas.figure.get_axes()[0]
       if self.LineMoving is not None: 
          self.LineMoving.remove()
          del self.LineMoving 
          self.LineMoving = None 
          
       if b is not None :
          fig = e.canvas.figure
          size = fig.get_size_inches()*fig.dpi
          #
          # pointer should be in right axes 
          #
          if e.x < size[0]*0.8 :
              #
              # also pointer should be inside triangle 
              #
              # if (ix,iy) is close to a point existing in the database
              Inside, ix,iy,key = IsInsideOrClose(e,self.Project)
   
              if Inside: 
                  #
                  # restore saved backgrounds
                  #
                  for background in self.Backgrounds:
                      self.canvas.restore_region(background)
                  #
                  # make visible and shift cursor
                  #
                  self.CursorPoint.set_visible(True)
                  self.CursorPoint.set_offsets((ix, iy))

                  l = len(self.Project.PickedCoordinates)
                  if l > 0 : 
                     x=[];y=[]
                     x=np_append(x,self.Project.PickedCoordinates[l-1][0])
                     x=np_append(x,ix)
                     y=np_append(y,self.Project.PickedCoordinates[l-1][1])
                     y=np_append(y,iy)
                     lines = []
                     lines.append( list(zip(x,y)))
                     line_segments = LineCollection( lines, color ='red',\
                        lw=2,zorder=9)
                     self.LineMoving = a.add_collection(line_segments)
                     #
                     # blit canvas, whatever this means. anyway its definitely
                     # faster compared to redrawing everthing on the canvas
                     #
                     k=0
                     for axes in self.canvas.figure.get_axes():
                        k+=1
                        if k==1:
                            #
                            # update point on first axis (cursor)
                            #
                            axes.draw_artist( self.LineMoving )
                            axes.draw_artist( self.CursorPoint )
                        self.canvas.blit(axes.bbox)
              else:
                 self.CursorPoint.set_visible(False)
                 if self.LineMoving is not None: 
                    del self.LineMoving 
                    self.LineMoving = None 
                    self.LineMoving.remove()
              #self.canvas.draw()

    # -------------------------------------------------------------------------
    # event bindings present when intializing 
    # -------------------------------------------------------------------------
    def OnInitMotion(self,e):
       from poscar import GetFormulaAll,PHASES, GetGrandPotential, kJmol
       
       a = e.inaxes
       if a is not None :
          fig = e.canvas.figure
          size = fig.get_size_inches()*fig.dpi
          #
          # pointer should be in right axes 
          #
          if e.x < size[0]*0.8 :
              #
              # also pointer should be inside triangle 
              #
              # if (ix,iy) is close to a point existing in the database
              Inside, ix,iy,key = IsInsideOrClose(e,self.Project)
    
              if Inside: 
                  g=GetGrandPotential(self.Project,ix,iy)
                  text = 'Formation energy: {:.4f} [eV] ({:.4f} [kJ/mol])'.\
                         format(g,g*kJmol)
                  self.statusbar.SetStatusText(text,0)

                  text = GetFormulaAll(self.Project, ix,iy)

                  for i in [0,1,2]:
                      p=PHASES[i]
                      self.statusbar.SetStatusText(text[i],i+1)
              else:
                  for i in [0,1,2]:
                    self.statusbar.SetStatusText('',i+1)

    # -------------------------------------------------------------------------
    # decide wether to load 2D or 1D linetool form
    # -------------------------------------------------------------------------
    def Load2DLineToolForm(self):
       from numpy import dot, array,subtract
       from poscar import small
       n = len( self.Project.PickedCoordinates )
       if n == 5 : 
           x0 = array( self.Project.PickedCoordinates[ 0 ] )
           x1 = array( self.Project.PickedCoordinates[ 4 ] )
           d = subtract( x0,x1)
           if dot( d,d ) < small:
              return True
           else:
              return False
       else:
          return False

# *****************************************************************************
# HELPER routines used by GUI app 
# *****************************************************************************
def ClosestPoint(x,y,sx):
    from numpy import argmin
    #find index in loaded data closest to 
    tmp = ((x-sx.hull.points[:,0])**2,(y-sx.hull.points[:,1])**2)
    idx = argmin(tmp[0]+tmp[1])
    tmp=sx.hull.points[idx,:]
    dist=(x-tmp[0])**2 + (y-tmp[1])**2
    return dist,tmp[0],tmp[1],sx.hull.keys[idx]
 
def IsInsideOrClose(e,sx):
    """ checks wheather point is in the triangle, if true and point is as close
    to one data point as tol the function returns the coordinates of the data
    point otherwise the input coordinates are returned """
    from poscar import small 
    tol=0.5E-4
    # squareroot of 3
    sq3=1.7320508075688772

    #x and y coordinates are:
    x=e.xdata; y=e.ydata

    if x is not None and x is not None:
    # 3 conditions have to be satisfied such that X is inside the equilateral
    # tirangle:  (i) y>0 
    #           (ii) y < sqrt3 x 
    #          (iii) y < sqrt3 (1-x)
    #
    # we add a small tolarance 0 -> small to make comparisions safe
       if y > small and (sq3*x-y) > small and (sq3*(1-x))-y > small:
          IsInside = True 
       else:
          IsInside = False

    #check if (x,y) is close to a point in sx.hull.points
       d,xi,yi,key=ClosestPoint(x,y,sx)
      
       if d>tol:
          xi=x
          yi=y
          key=''
       else: 
          IsInside = True
    else:
       IsInside=False
       xi=1000
       yi=1000
         
    return IsInside, xi,yi,key

def OpenStructure(Project,key, status,logger):
   """ opens structure with VESTA """
   from poscar import MakePOSCAR
   from poscar import MakeVESTA
   import subprocess
   from os import name as osname
   from os import path,getcwd
   from dialog import ShowStatusErrorMessage


   if len(key)==0:
      text = 'Entry not found'
      ShowStatusErrorMessage( status, text )
      return
   else:
      if not "coordinates" in Project.db[key]:
         text = 'No structure found for {:}'.format( key )
         ShowStatusErrorMessage( status, text )
         return


   if len( Project.pathVESTA )>0 : 
      if len( Project.OutDir ) > 0: 
         d = str( Project.OutDir  )
      else:
         text = 'No output directory given'
         ShowStatusWarningMessage( status, text )
         d = getcwd()
      #
      # create an POSCAR file 
      #
      Path=path.join(d,key)  #name without extension
      fname="{:}vesta".format(Path)
      #MakePOSCAR( key, Project, fname ) 
      MakeVESTA( key, Project, fname ) 
      if osname == 'nt':
         f=fname
	 fname='\"'+f+'\"'
      l = len( Project.db[key]["SameBC"] )
      if l >=  1: 
         for i in Project.db[key]["SameBC"] :
            Path=path.join(d,i)  #name without extension
            f='{:}vesta'.format(Path)
            #MakePOSCAR( i, Project, f ) 
            MakeVESTA( i, Project, f ) 
            if osname == 'nt':
	       fname+=' \"'+f+'\"'
            else:
	       fname+=' '+f
      #
      # set up command and exectute
      #
      if osname == 'nt':
         cmd="\"{:}\" {:}".format(Project.pathVESTA, fname)
      else:
         cmd="{:} {:}".format(Project.pathVESTA, fname)
      p = subprocess.Popen(cmd, shell=True)
      
      text="               {:} structures stored in {:}".format(l+1,d)
      logger.Add2Log(text)
   else:
      text = 'VESTA exectuable not found!'
      ShowStatusErrorMessage( status, text )

def EditStructureFiles(Project,key, status,logger):
   """ opens structure with Editor """
   import subprocess
   from os import name as osname
   from os import path,getcwd
   from dialog import ShowStatusErrorMessage


   if len(key)==0:
      text = 'Entry not found'
      ShowStatusErrorMessage( status, text )
      return
   else:
      if not "coordinates" in Project.db[key]:
         text = 'No structure found for {:}'.format( key )
         ShowStatusErrorMessage( status, text )
         return

   print Project.db[ key ]["path"]
   l = len( Project.db[key]["SameBC"] )
   if l >=  1: 
      for i in Project.db[key]["SameBC"] :
         print Project.dbF[ i ]["path"]
