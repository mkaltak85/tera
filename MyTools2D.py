import wx 
import icons25 as ico
from matplotlib.backends.backend_wx import NavigationToolbar2Wx
from poscar import Data

class PlotSettings(Data):
    def __init__(self):
       self.plot = 0 
       self.norm = 2 
       self.segregation = 0 
       self.xlabel = ''

###############################################################################
# subclassed from src in: /usr/local/lib/python2.7/dist-packages/matplotlib/
class CustomNavigationToolbar2Wx(NavigationToolbar2Wx):
    def __init__(self, canvas,parent,Project, Data):
        from matplotlib.backend_bases import NavigationToolbar2
        wx.ToolBar.__init__(self, canvas.GetParent(), -1)
        NavigationToolbar2.__init__(self, canvas)
        self.canvas = canvas
        self.Project = Project
        self.Settings = PlotSettings()
        self.Data = Data
        self._idle = True
        self.statbar = None
        self.prevZoomRect = None
        # for now, use alternate zoom-rectangle drawing on all
        # Macs. N.B. In future versions of wx it may be possible to
        # detect Retina displays with window.GetContentScaleFactor()
        # and/or dc.GetContentScaleFactor()
        self.retinaFix = 'wxMac' in wx.PlatformInfo
        
        self.Initcids = [ ] 
        self.Initcids.append( canvas.mpl_connect('motion_notify_event',\
           self.OnInitMotion) )

    def _init_toolbar(self):
        self._parent = self.canvas.GetParent()
        #
        # store project data 
        #
        self.Project = self._parent.Project
        self.Settings = PlotSettings()
        self.Data = self._parent.Data
        #
        # initialize Tool ID dictionary 
        #
        self.wx_ids = {}
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
        ('Save', 'Save', 'Save figure', 'OnSave'), \
        #('Save', 'Save', 'Save figure', 'save_figure'), \
        )
        for text, icon, tooltip_text, callback in self.items:
            if text is None:
                self.AddSeparator()
                continue
            self.wx_ids[text] = wx.NewId()
            img=getattr(ico,icon).GetImage().ConvertToBitmap()
            self._AddTool(self, self.wx_ids, text,
                        img,
                        tooltip_text)
            if hasattr(self,callback):
               self.Bind(wx.EVT_TOOL, getattr(self, callback),
                      id=self.wx_ids[text])
        self.Realize()

        # 
        #  plot data
        # 
        self.InitPlot() 

    def _AddTool(self,parent, wx_ids, text, bmp, tooltip_text):
        if text in ['Pan', 'Zoom']:
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

    def OnSave(self, *args):
       from os import getcwd
       from numpy import savetxt
       wildcard = \
           "Joint Photographic Experts Group (*.jpeg,*.jpg)|*.jpeg;*.jpg|"\
           "Tagged Image File Format (*.tiff,*.tif)|*.tiff;*.tif|"\
           "Portable Network Graphics (*.png)|*.png|"\
           "Encapsulated Postscript Vector (*.eps)|*.eps|"\
           "Portable Document Format (*.pdf)|*.pdf|"\
           "Scalable Vector Graphics (*.svg)|*.svg|"\
           "Raw data (*.dat,*.txt)|*.dat;*.txt"
       #store extentions
       EXT = ["jpg","tif","png","eps","svg","pdf","dat"]
       currentDirectory = getcwd()
       dlg = wx.FileDialog(
           self, message="Save figure as",
           wildcard=wildcard,
           defaultDir=currentDirectory, 
           defaultFile="",
           style=wx.SAVE|wx.FD_OVERWRITE_PROMPT
           )
       if dlg.ShowModal() == wx.ID_OK:
           path = dlg.GetPaths()[0]
           #
           # in case one file is chosen
           #
           if len( path ) > 0: 
              fname='{:}.{:}'.format(path,EXT[dlg.GetFilterIndex()])
              # in case of a text file, the raw data has to be obtained from
              # figure
              if dlg.GetFilterIndex() == 6: 
                 ax = self.canvas.figure.get_axes()[0]
                 xydata=ax.lines[0].get_xydata()
                 p =  self.GetTopLevelParent().statusbar
                 text = 'Data saved to: {:}'.format( fname )
                 savetxt( fname, xydata, fmt='%20.6f')
              else:
                 fig = self.canvas.figure
                 fig.savefig(fname)
                 text = 'Figure saved to: {:}'.format( fname )
              p =  self.GetTopLevelParent().statusbar
              p.SetStatusText(text)

       dlg.Destroy()

    def InitPlot(self,Select=None):
        """ Plot data """
        Iam = self.GetParent().Iam
        
        if Select is not None: 
           self.Settings.plot = Select[0]
           self.Settings.norm = Select[1]
           self.Settings.xlabel = Select[2]
           self.Settings.segregation = Select[3]

        # decide what to plot 
        if Iam == 'PageVolt' : 
           self.PlotChemPot() 

        # decide what to plot 
        if Iam == 'PageStab' : 
           self.PlotStability() 

        # decide what to plot 
        if Iam == 'PageSegr' : 
           self.PlotSegregation() 

    def PlotStability(self): 
        """ Plots Formation energy  """
        from numpy import zeros, abs, linspace
        from poscar import small
	from scipy.spatial import ConvexHull
        import matplotlib.pyplot as plt

        # decide which axis to plot 
        color = [ 'blue','red' ] 
        linewidth=[ 1.5,1.5 ] 

        self.canvas.figure.clear()

        a = self.canvas.figure.add_subplot(111)

        x,y,z,levels = self.GetStabData()

        cp=a.contourf( x, y,z,levels)
        self.canvas.figure.colorbar(cp)

        ax = self.canvas.figure.get_axes()[0]
        ax.set_xlim( [ 0, 8 ] )
        ax.set_ylim( [ 0, 2 ] )
        a.set_aspect('equal','datalim')


        ylabel = 'Formation Energy'
        ax.set_ylabel( ylabel ) 
        ax.set_xlabel( self.Settings.xlabel ) 


        self.canvas.draw()
	
    def PlotChemPot(self, debug=None): 
        """ Plot derivative of grand potential with a negative sign """
        from numpy import zeros , append
	from scipy.spatial import ConvexHull
        # decide which chemical potential to plot
        color = [ 'red', 'green','blue' ] 
        linewidth=[ 1.5, 1.5, 1.5 ] 
        self.canvas.figure.clear()
        a = self.canvas.figure.add_subplot(111)
        ax = self.canvas.figure.get_axes()[0]
        ylabel = 'Chemical potential [V]'
        ax.set_ylabel( ylabel ) 
        ax.set_xlabel( self.Settings.xlabel ) 

        if debug :
           print " ---------------------------------------------- "
           print " xscaled ", x
      
        self.canvas.draw()

    def PlotSegregation(self): 
        """ Plots segregation of one enpoint """
        from numpy import zeros 
        from poscar import PHASES
	from scipy.spatial import ConvexHull
        # decide which axis to plot 
        iseg=self.Settings.segregation
        color = [ 'blue','red' ] 
        linewidth=[ 1.5,1.5 ] 
        self.canvas.figure.clear()
        a = self.canvas.figure.add_subplot(111)
        ax = self.canvas.figure.get_axes()[0]
        ylabel =' Segregation of {:} in %'.format(\
            self.Project.pdb[ PHASES[ iseg ] ]["LaTeXName"])
        ax.set_ylabel( ylabel ) 
        ax.set_xlabel( self.Settings.xlabel ) 


        self.canvas.draw()


    def GetSegData(self,ScaleTo100=None ):
        """ Returns segregation data """
        from numpy import zeros, abs
        from poscar import small, CalculateSegregation
        x=[0]
        y=[0]
        return x, y

    def GetStabData(self):
        """ Returns stability data """
        from numpy import zeros, linspace
        from numpy import append as np_append
        from matplotlib.mlab import griddata

        isel = self.Settings.plot
        x=zeros( len( self.Data.X ) )
        y=zeros( len( self.Data.X ) )
        x = self.Data.X[:,isel  ]
        y = self.Data.X[:,isel+1]

        dim=self.Project.nsample
        z = self.Data.G[:,isel  ]
        levels=linspace( z.min(), z.max(),dim)
        Z=zeros( shape=( dim ,dim ) )
        X=zeros( shape=( dim ,dim ) )
        Y=zeros( shape=( dim ,dim ) )
        for i in range( 0, dim):
           Zp=[ ]
           Xp=[ ]
           Yp=[ ]
           for j in range( 0, dim) :
              Zp = np_append( Zp, z[ i*dim+ j ] )
              Xp = np_append( Xp, x[ i*dim+ j ] )
              Yp = np_append( Yp, y[ i*dim+ j ] )
           Z[i] = Zp
           X[i] = Xp
           Y[i] = Yp

        return X,Y,Z,levels

    def OnInitMotion(self,e):
        p =  self.GetTopLevelParent().statusbar
        x=e.xdata ; y=e.ydata
        if x is not None and y is not None:
           text = ' x = {:.4}  y = {:.4} '.format(x,y)
           p.SetStatusText(text)

class PanelSettings(wx.Panel):
    """ Settings Canvas: contains settings to select """
    #----------------------------------------------------------------------
    def __init__(self, parent, Project):
        wx.Panel.__init__(self, parent=parent, id=wx.ID_ANY)
        self.parent=parent
        self.Project = Project
        self.Build()

    def Build(self):
        # *********************************************************************
        sizer = wx.BoxSizer(wx.VERTICAL)
        # *********************************************************************

        # 
        # initial normalization selection 
        # 
        t1 = " {:} ".format(\
           self.Project.pdb['B']["name"])
        t2 = " {:} ".format(\
           self.Project.pdb['C']["name"])

        self.NormalizationEntries = [ t1,t2 ]
        self.SelectorBox = wx.RadioBox(self, 
                        label = 'Normalize to', 
                        choices = self.NormalizationEntries,
                        style = wx.RA_VERTICAL ) 
        tip = 'Select how plots should be normalized'
        text = wx.ToolTip( tip )
        self.SelectorBox.SetToolTip( text )

        # Chemical Potential selector 
        lblList=[self.Project.pdb['A']["name"],
                 self.Project.pdb['B']["name"],
                 self.Project.pdb['C']["name"] ]
        self.RadioBox = wx.RadioBox(self, 
                        label = 'Concentration', 
                        choices = lblList,
                        style = wx.RA_VERTICAL ) 
        self.RadioBox.Bind(wx.EVT_RADIOBOX,self.OnEndpointSelector)  
        tip = 'Select which concentration to plot on the x-axis'
        text = wx.ToolTip( tip )
        self.RadioBox.SetToolTip( text )

        self.SelectorBox.Bind(wx.EVT_RADIOBOX,self.OnNormSelector)  

        # Segregation selector
        self.SegregationBox = wx.RadioBox(self, 
                        label = 'Segregation', 
                        choices = lblList,
                        style = wx.RA_VERTICAL ) 
        self.SegregationBox.Bind(wx.EVT_RADIOBOX,self.OnEndpointSelector)  
        tip = 'Select which segregation to plot on the y-axis'
        text = wx.ToolTip( tip )
        self.SegregationBox.SetToolTip( text )

        # Chemical Potential selector 
        sizer.Add(self.RadioBox, 1, wx.EXPAND)

        # segregation box
        sizer.Add(self.SegregationBox, 1, wx.EXPAND)

        # Normalization 
        sizer.Add(self.SelectorBox, 1, wx.EXPAND)
        #
        # finish sizer
        #
        # *********************************************************************
        self.SetSizer(sizer)
        # *********************************************************************


    def SetNormSelector(self,*args):
       # 
       # select normalization
       # 
       t1 = " {:} ".format(\
          self.Project.pdb['B']["name"])
       t2 = " {:} ".format(\
          self.Project.pdb['C']["name"])

       t3 = " {:} ".format(\
          self.Project.pdb['A']["name"])
       t4 = " {:} ".format(\
          self.Project.pdb['C']["name"])

       t5 = " {:} ".format(\
          self.Project.pdb['A']["name"])
       t6 = " {:} ".format(\
          self.Project.pdb['B']["name"])

       plot=self.RadioBox.GetSelection() 
       if  plot == 0 :
          self.SelectorBox.SetString( 0, t1 )
          self.SelectorBox.SetString( 1, t2 )
       elif  plot == 1 :
          self.SelectorBox.SetString( 0, t3 )
          self.SelectorBox.SetString( 1, t4 )
       elif  plot == 2 :
          self.SelectorBox.SetString( 0, t5 )
          self.SelectorBox.SetString( 1, t6 )

    def GetNorm(self):
       # get Norm selector choice
       plot=self.RadioBox.GetSelection() 
       n=self.SelectorBox.GetSelection()
       if plot == 0 :
           norm = 2
           xlabel = 'Concentration {:}$_x$ per {:}'.\
                    format(\
                    self.Project.pdb[ 'A' ]['LaTeXName'],
                    self.Project.pdb[ 'C' ]['LaTeXName']
                    )
           if n == 0 :
              norm = 1
              xlabel = 'Concentration {:}$_x$ per {:}'.\
                    format(\
                    self.Project.pdb[ 'A' ]['LaTeXName'],
                    self.Project.pdb[ 'B' ]['LaTeXName']
                    )
       elif plot == 1 :
           norm = 2
           xlabel = 'Concentration {:}$_y$ per {:}'.\
                    format(\
                    self.Project.pdb[ 'B' ]['LaTeXName'],
                    self.Project.pdb[ 'C' ]['LaTeXName']
                    )
           if n == 0 :
              norm = 0
              xlabel = 'Concentration {:}$_y$ per {:}'.\
                    format(\
                    self.Project.pdb[ 'B' ]['LaTeXName'],
                    self.Project.pdb[ 'A' ]['LaTeXName']
                    )
       elif plot == 2 :
           norm = 1
           xlabel = 'Concentration {:}$_z$ per {:}'.\
                    format(\
                    self.Project.pdb[ 'C' ]['LaTeXName'],
                    self.Project.pdb[ 'B' ]['LaTeXName']
                    )
           if n == 0 :
              norm = 0
              xlabel = 'Concentration {:}$_z$ per {:}'.\
                    format(\
                    self.Project.pdb[ 'C' ]['LaTeXName'],
                    self.Project.pdb[ 'A' ]['LaTeXName']
                    )
       return norm, xlabel

    def OnEndpointSelector(self,*args):
       plot=self.RadioBox.GetSelection() 
       # set Norm selector
       self.SetNormSelector()
       # get Norm selector choice
       norm,xlabel = self.GetNorm()
       # get segregation
       segregation = self.SegregationBox.GetSelection()
       # and plot 
       nb=self.parent.GetParent().GetParent().A.nb
       current_page = nb.GetSelection()
       if current_page >= 0 : 
          tb=nb.GetPage(current_page).MPLToolBar
          tb.InitPlot( Select=[plot,norm,xlabel,segregation] )

    def OnNormSelector(self,e):
       plot=self.RadioBox.GetSelection() 
       # get segregation
       segregation = self.SegregationBox.GetSelection()
       # and plot 
       norm,xlabel = self.GetNorm()
       nb=self.parent.GetParent().GetParent().A.nb
       current_page = nb.GetSelection()
       if current_page >= 0 : 
          tb=nb.GetPage(current_page).MPLToolBar
          tb.InitPlot( Select=[plot,norm,xlabel,segregation] )


###############################################################################
class TabPanel(wx.Panel):
    """
    This describes the page layout. every page has its own matplotlib figure
    """
    #--------------------------------------------------------------------------
    def __init__(self, parent, name, Project =None ):
        from matplotlib.figure import Figure
        from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg as \
           FigureCanvas

        wx.Panel.__init__(self, parent=parent, id=wx.ID_ANY)
        self.parent=parent
        self.Iam = name 
        self.Project = Project
        # store reference to parent containing the requested data 
        self.Data = parent.GetParent()
        # *********************************************************************
        sizer = wx.BoxSizer(wx.VERTICAL)
        # *********************************************************************
        #
        # initialize MPL figure canvas
        #
        self.figure = Figure()
        self.figure.add_subplot(111)
        self.canvas = FigureCanvas(self, -1, self.figure)
        #
        # add customized matplotlib toolbar
        #
        self.MPLToolBar = CustomNavigationToolbar2Wx(self.canvas,parent,\
           self.Project, self.Data )
        sizer.Add(self.MPLToolBar, 0, wx.EXPAND)
        sizer.Add(self.canvas, 1, wx.EXPAND)
        #
        # finish sizer
        #
        # *********************************************************************
        self.SetSizer(sizer)
        # *********************************************************************
        

###############################################################################
class PanelCanvas(wx.Panel):
    """ Panel Canvas: contains Notebook with MPL canvas"""
    #----------------------------------------------------------------------
    def __init__(self, parent, Project):
        from matplotlib.figure import Figure
        from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg as \
           FigureCanvas
        wx.Panel.__init__(self, parent)
# DEBUG:
#        self.SetBackgroundColour('green')
# DEBUG:
        self.Project = Project
        # 
        # attach data to project 
        lcontinue = self.CalculateData(debug=False ) 

        self.nb = wx.Notebook(self,style=wx.BK_BOTTOM)
        items = ( 
                ('PageLatt','LattCanvas','LattTB','Lattice Parameter'),
                ('PageVolt','VoltCanvas','VoltTB','Chemical Potential'),
                ('PageSegr','SegrCanvas','SegrTB','Segregation'),
                ('PageStab','StabCanvas','StabTB','Formation Energy'), )
        
        sizer=wx.BoxSizer()
        for page, canvas, toolbar, title in items:
           tab=TabPanel(self.nb, page, Project = self.Project)
           setattr(self,page,tab)
           self.nb.AddPage(tab, title , select=True)

        sizer.Add(self.nb,1,wx.EXPAND|wx.ALL)

        self.nb.Bind(wx.EVT_NOTEBOOK_PAGE_CHANGED, self.OnPageChanged)
        #
        # add sizer 
        #  
        sizer=wx.BoxSizer()
        sizer.Add(self.nb, 1, wx.EXPAND|wx.ALL)
        self.SetSizer(sizer)
    
    def OnPageChanged(self,*args):
       t=self.GetTopLevelParent().B
       current_page = self.nb.GetSelection()
       if current_page >=0  : 
          if hasattr(t, 'OnEndpointSelector'):
             t.OnEndpointSelector()

    # -------------------------------------------------------------------------
    def CalculateData(self,debug=None):
       """ calculate data, like lattice deformation, voltage, Convex hull cut"""
       from numpy import array, zeros,copy, abs
       from numpy import append as np_append
       from poscar import Get2DFormationEnergies
       from poscar import small
       
       nsteps=len( self.Project.PickedCoordinates ) 
       self.Project.nsample=50
     
       x=zeros( shape=( 4, 3 ) ) 
       x[0]=copy( self.Project.PickedCoordinates[0][2:5])
       x[1]=copy( self.Project.PickedCoordinates[1][2:5])
       x[2]=copy( self.Project.PickedCoordinates[2][2:5])
       x[3]=copy( self.Project.PickedCoordinates[4][2:5])

       # make a list of indices for which the concentration is non-zero all four points,
       self.idx=[]
       lidx=[]
       for j in [0,1,2]:
          k=-1
          for i in [ 0,1,2,3]:
             if abs( x[i,j] ) > small:
                k+=1
          if k == 3 : 
             self.idx.append( j ) 
             if j==0 :
                 a=[1,2]
             elif j==1 :
                 a=[2,0]
             else:
                 a=[0,1]
             lidx.append( a ) 

       n=self.Project.nsample*self.Project.nsample
       # formation energy 
       self.G = zeros( shape= ( n , 2*len( self.idx ) ) ) 
       # data points 
       self.X = zeros( shape= ( n , 2*len( self.idx ) ) ) 

       ITER=-1
       for norm in self.idx:
          ITER+=1

          x0=[] ;  y0=[]
          l = lidx[ ITER ]
          x0=np_append(x0, x[0][ l[0] ] / x[0][ norm ] )
          x0=np_append(x0, x[1][ l[0] ] / x[1][ norm ] )
          y0=np_append(y0, x[3][ l[1] ] / x[3][ norm ] )
          y0=np_append(y0, x[2][ l[1] ] / x[2][ norm ] )
          
          if debug: 
             print "Norm:", norm, l
             print x0, y0
          k=-1
          # use uniform spacing for line segments
          # 
          for i in range( 0, self.Project.nsample) : 
             xi = x0[0]+ i*(x0[1]-x0[0])/float( self.Project.nsample - 1)
             for j in range( 0, self.Project.nsample) : 
                k+=1
                yi = y0[0]+ j*(y0[1]-y0[0])/float( self.Project.nsample - 1)
                self.X[k, ITER*(len(self.idx)-1)  ] = xi 
                self.X[k, ITER*(len(self.idx)-1)+1] = yi 
                
                xb = array( [ 0.,0.,0.] )
                ni=( xi+yi+1. )
                xb[ norm ] = 1./ni
                xb[ l[0] ] = xi/ni
                xb[ l[1] ] = yi/ni
              
                r=Get2DFormationEnergies( self.Project, xb )

                self.G[ k, ITER*(len(self.idx)-1)  ] = r[0]*ni
                self.G[ k, ITER*(len(self.idx)-1)+1] = r[1]*ni

                if debug:
                   print norm, l, xi,yi, xb, r[0],r[1], r[0]*ni, r[1]*ni

###############################################################################
class LineTool2DForm(wx.Dialog):
    
    def __init__(self, *args, **kw ):
        from poscar import Data
        super(LineTool2DForm, self).__init__(*args, **kw) 

        self.MainNB = self.GetParent().A.NoteBook
        current_page = self.MainNB.GetSelection()      
        if current_page >= 0 : 
           self.Project=self.MainNB.GetPage(current_page).Project
           nsteps=len( self.Project.PickedCoordinates ) 
        self.InitUI()
        self.SetTitle("2D Line Tool ")
        self.SetSize((800, 600))
        self.Bind(wx.EVT_CLOSE, self.OnClose)

    def InitUI(self):
        from widgets import ProportionalSplitter

        self.MainPanel = wx.Panel(self)           
        self.MainPanel.SetBackgroundColour('green') 
        # *********************************************************************
        vbox = wx.BoxSizer(wx.VERTICAL) # main box contains other horiz. boxes
        # *********************************************************************
        self.vSplitter = ProportionalSplitter(self.MainPanel, 
             id=-1, proportion = 0.78, style = wx.SP_BORDER) 
        self.A = PanelCanvas( self.vSplitter, self.Project )
        self.B = PanelSettings( self.vSplitter, self.Project )
        self.vSplitter.SplitVertically(self.A, self.B)
        self.vSplitter.SetSashGravity(0.5)
        vbox.Add(self.vSplitter,1, wx.EXPAND)
        self.statusbar = wx.StatusBar(self.MainPanel)
        vbox.Add(self.statusbar,0,wx.EXPAND)
        # *********************************************************************
        self.MainPanel.SetSizer(vbox) # finialize main vertical box 
        # *********************************************************************
              
    def OnClose(self, e):
        # 
        # disconnect matplotlib hook
        #
        items = [ 'PageLatt', 'PageVolt', 'PageSegr', 'PageStab' ]
        for page in items:
          p = getattr( self.A,page)
          tb = p.MPLToolBar
          for  cid in tb.Initcids:
              tb.canvas.mpl_disconnect(cid)       
              tb.Initcids = [ ]
        self.Destroy()

    def OnSave(self, e):
       from os import getcwd
       wildcard = \
           "Joint Photographic Experts Group (*.jpeg,*.jpg)|*.jpeg;*.jpg|"\
           "Tagged Image File Format (*.tiff,*.tif)|*.tiff;*.tif|"\
           "Portable Network Graphics (*.png)|*.png|"\
           "Encapsulated Postscript Vector (*.eps)|*.eps|"\
           "Portable Document Format (*.pdf)|*.pdf|"\
           "Scalable Vector Graphics (*.svg)|*.svg|"\
           "Raw data (*.dat)|*.dat"
       currentDirectory = getcwd()
       dlg = wx.FileDialog(
           self, message="Save figure as",
           wildcard=wildcard,
           defaultDir=currentDirectory, 
           defaultFile="",
           style=wx.SAVE|wx.FD_OVERWRITE_PROMPT
           )
       if dlg.ShowModal() == wx.ID_OK:
           path = dlg.GetPaths()[0]
           #
           # in case one file is chosen
           #
           if len( path ) > 0: 
               print path
       dlg.Destroy()
        

