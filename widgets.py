# -*- coding: utf-8 -*-
import wx
import icons25 as ico
import icons35 as ico25
import wx.lib.agw.aui as aui
from  wx.lib.scrolledpanel import ScrolledPanel
from numpy import arange, sin, pi

# matplotlib src found here /usr/local/lib/python2.7/dist-packages/matplotlib/
import matplotlib
matplotlib.use('WXAgg')

###############################################################################
class ProportionalSplitter(wx.SplitterWindow):
    def __init__(self,parent, id = -1, proportion=0.66, \
    size = wx.DefaultSize, **kwargs):
       wx.SplitterWindow.__init__(self,parent,id,wx.Point(0, 0),size, **kwargs)
       self.SetMinimumPaneSize(1) #the minimum size of a pane.
       self.proportion = proportion
       if not 0 < self.proportion < 1:
          raise ValueError, \
          "proportion value for ProportionalSplitter must be between 0 and 1."
       self.ResetSash()
       self.Bind(wx.EVT_SIZE, self.OnReSize)
       self.Bind(wx.EVT_SPLITTER_SASH_POS_CHANGED, self.OnSashChanged, id=id)
       ##hack to set sizes on first paint event
       self.Bind(wx.EVT_PAINT, self.OnPaint)
       self.firstpaint = True

    def SplitHorizontally(self, win1, win2):
       if self.GetParent() is None: return False
       i=int(round(self.GetParent().GetSize().GetHeight() * self.proportion))
       return wx.SplitterWindow.SplitHorizontally(self, win1, win2,i)

    def SplitVertically(self, win1, win2):
       if self.GetParent() is None: return False
       i=int(round(self.GetParent().GetSize().GetWidth() * self.proportion))
       return wx.SplitterWindow.SplitVertically(self, win1, win2, i)

    def GetExpectedSashPosition(self):
       size=self.GetParent().GetClientSize()
       if self.GetSplitMode() == wx.SPLIT_HORIZONTAL:
          tot = max(self.GetMinimumPaneSize(),size.height)
       else:
          tot = max(self.GetMinimumPaneSize(),size.width)
       return int(tot*self.proportion)

    def ResetSash(self):
       self.SetSashPosition(self.GetExpectedSashPosition())

    def OnReSize(self, event):
       """Window has been resized, so we need to adjust the sash based on
          self.proportion."""
       self.ResetSash()

    def OnSashChanged(self, event):
       """We'll change self.proportion now based on where user dragged the
          sash."""
       pos = float(self.GetSashPosition())
       if self.GetSplitMode() == wx.SPLIT_HORIZONTAL:
          tot = self.GetParent().GetClientSize().height
          self.proportion = pos/tot
       else:
          tot = \
          self.GetParent().GetClientSize().width
          self.proportion = pos/tot
       event.Skip()

    def OnPaint(self,event):
       if self.firstpaint:
          if self.GetSashPosition() != self.GetExpectedSashPosition():
             self.ResetSash()
          self.firstpaint = False
       event.Skip()

###############################################################################
class TestPanel(wx.Panel):
    """ Test class showing standard paths of application """
    def __init__(self, parent):
        wx.Panel.__init__(self, parent, -1)

        sizer = wx.FlexGridSizer(0, 3, 5, 5)
        sizer.AddGrowableCol(1)
        box = wx.BoxSizer(wx.VERTICAL)
        fs = self.GetFont().GetPointSize()
        bf = wx.Font(fs+4, wx.SWISS, wx.NORMAL, wx.BOLD)

        t = wx.StaticText(self, -1, "StandardPaths")
        t.SetFont(bf)
        box.Add(t, 0, wx.CENTER|wx.ALL, 4)
        box.Add(wx.StaticLine(self, -1), 0, wx.EXPAND)

        # get the global (singleton) instance of wx.StandardPaths
        sp = wx.StandardPaths.Get()

        # StandardPaths will use the value of wx.App().GetAppName()
        # for some of the stnadard path components.  Let's set it to
        # something that makes that obvious for the demo.  In your own
        # apps you'll set it in to something more meaningfull for your
        # app in your OnInit, (or just let it default.)
        wx.GetApp().SetAppName("tera")

        self.help = {}

        # Loop through all of the getters in wx.StandardPaths and make
        # a set of items in the sizer for each.
        def makeitem(name, *args):
            func = getattr(sp, name)
            sizer.Add(wx.StaticText(self, -1, "%s%s:" %(name, repr(args))),
                      0, wx.ALIGN_RIGHT|wx.ALIGN_CENTER_VERTICAL)
            sizer.Add(wx.TextCtrl(self, -1, func(*args),
                                  size=(275,-1), style=wx.TE_READONLY),
                      0, wx.EXPAND|wx.ALIGN_CENTER_VERTICAL)

            btn = wx.Button(self, wx.ID_HELP)
            sizer.Add(btn)
            self.help[btn] = func.__doc__

        for x in ['GetConfigDir',
                  'GetUserConfigDir',
                  'GetDataDir',
                  'GetLocalDataDir',
                  'GetUserDataDir',
                  'GetUserLocalDataDir',
                  'GetDocumentsDir',
                  'GetPluginsDir',
                  'GetInstallPrefix',
                  'GetResourcesDir',
                  'GetTempDir',
                  'GetExecutablePath',
                  ]:
            makeitem(x)

        # this one needs parameters
        makeitem('GetLocalizedResourcesDir', 'en',
                 wx.StandardPaths.ResourceCat_Messages )

        self.Bind(wx.EVT_BUTTON, self.OnShowDoc, id=wx.ID_HELP)

        box.Add(sizer, 0, wx.CENTER|wx.EXPAND|wx.ALL, 20)
        self.SetSizer(box)

    def OnShowDoc(self, evt):
        doc = self.help[evt.GetEventObject()]

        # trim the whitespace from each line
        lines = []
        for line in doc.split('\n'):
            lines.append(line.strip())
        doc = '\n'.join(lines)
        wx.TipWindow(self, doc, 500)    

###############################################################################
class PreferencesSettings(ScrolledPanel):
    """ Preference class showing settings which can be changed in preferences """
    def __init__(self, parent):
        ScrolledPanel.__init__(self, parent)
        self.SetupScrolling(scroll_y = False)
# DEBUG:
#        self.SetBackgroundColour('green')
# DEBUG:
        #
        # construct sizer 
        #  
        box=wx.BoxSizer(wx.VERTICAL)

        def makeItem(Label,tip,cols=2,bold=False):
            nm = wx.StaticBox(self, -1, Label+':') 
            fs = self.GetFont().GetPointSize()
            bf = wx.Font(fs+1, wx.SWISS, wx.NORMAL, wx.BOLD)  
            if bold:
               nm.SetFont(bf)
            CategoryBox=wx.StaticBoxSizer(nm,wx.VERTICAL)
            sizer = wx.FlexGridSizer(0, cols, 1, 1)
            sizer.AddGrowableCol(1)

            b = wx.Button(self, wx.ID_ANY,'Browse')
            t =wx.TextCtrl(self, -1, size=(200,b.GetDefaultSize()[1])) 
            text = wx.ToolTip( tip )
            text.SetDelay(0)
            t.SetToolTip(text )
            sizer.Add(t,0, wx.ALIGN_RIGHT)
            sizer.Add(b,wx.ALIGN_LEFT)
            obj=Label.replace(' ','')
            setattr(self,obj,t)
            b.Bind(wx.EVT_BUTTON,getattr(self,'On'+obj))
            CategoryBox.Add(sizer, 0, wx.CENTER|wx.EXPAND|wx.ALL, 10)
            box.Add(CategoryBox,0, wx.CENTER|wx.EXPAND|wx.ALL, 0)
        #
        # make vesta path item 
        #
        makeItem('VESTA path','Set path to VESTA executable')
        #
        # make vesta path item 
        #
        makeItem('Output directory','Set output directory')

        self.SetSizer(box)


    def OnVESTApath(self,e):
       """ opens file dialog """
       wildcard = "VESTA executable (*)|*" 
       dlg = wx.FileDialog(
           self, message="Choose VESTA executable",
           wildcard=wildcard,
           style=wx.OPEN | wx.CHANGE_DIR )
       if dlg.ShowModal() == wx.ID_OK:
           paths = dlg.GetPaths()[0]
           #
           # in case one file is chosen
           #
           if len( paths ) > 0: 
              self.VESTApath.SetValue(str(paths))
       dlg.Destroy()

    def OnOutputdirectory(self,e):
       """ opens directory dialog """
       from os import getcwd
       dlg = wx.DirDialog(self, "Choose directory with structures:",
                           style=wx.DD_DEFAULT_STYLE,
                           defaultPath = getcwd(), 
                           )
       if dlg.ShowModal() == wx.ID_OK:
          #
          # set directory control entry 
          #
          self.Outputdirectory.SetValue(dlg.GetPath())
       dlg.Destroy()

###############################################################################
#class PlotSettings(wx.Panel):
class PlotSettings(ScrolledPanel):
    """ Settings class showing settings which can be changed in plot """
    def __init__(self, parent):
        from wx.lib.buttons import GenButton
        from color import maps
        import wx.combo
        import icons20 as ico20
        ScrolledPanel.__init__(self, parent)
        self.cmaps = maps().AllMapNames
# DEBUG:
        #self.SetBackgroundColour('green')
# DEBUG:
        self.SetupScrolling() 
        box = wx.BoxSizer(wx.VERTICAL)
        #
        # helper to create entries 
        #
        def makeItem(x,entries,cols=2,bold=False):
            nm = wx.StaticBox(self, -1, x+':') 
            fs = self.GetFont().GetPointSize()
            bf = wx.Font(fs+1, wx.SWISS, wx.NORMAL, wx.BOLD)  
            if bold:
               nm.SetFont(bf)
            CategoryBox=wx.StaticBoxSizer(nm,wx.VERTICAL)
            sizer = wx.FlexGridSizer(0, cols, len(entries), len(entries))
            sizer.AddGrowableCol(1)
            for name in entries: 
               if name is None: 
                  continue
  
               if name =='Color x':
                  te= wx.StaticText(self, -1, 'Color y')
               elif name =='Color y':
                  te= wx.StaticText(self, -1, 'Color z')
               elif name =='Color z':
                  te= wx.StaticText(self, -1, 'Color x')
               elif name =='Ticks x':
                  te= wx.StaticText(self, -1, 'Ticks y')
               elif name =='Ticks y':
                  te= wx.StaticText(self, -1, 'Ticks z')
               elif name =='Ticks z':
                  te= wx.StaticText(self, -1, 'Ticks x')
               else:
                  te= wx.StaticText(self, -1, name)
               sizer.Add(te,0, wx.ALIGN_LEFT|wx.ALIGN_CENTER_VERTICAL)
               if name in ['Color','Color x','Color y','Color z']:
                  t = GenButton(self, label='', size=(60, 27))
                  t.SetBezelWidth(1)
                  t.SetUseFocusIndicator(False)
               elif name in ['Size','Resolution']:
                  imax=200
                  imin=1
                  if name == 'Resolution':
                     imin=50
                     imax=500
                  t = wx.Slider(self,value=imin,minValue=imin,maxValue=imax,size=(60,27))
                  text = wx.ToolTip(name+':'+str(t.GetValue()))
                  text.SetDelay(0)
                  t.SetToolTip(text)
               elif name in ['Width']:
                  t = wx.Slider(self,value=1,minValue=0.1,maxValue=5,size=(60,27))
                  text = wx.ToolTip(str(t.GetValue()))
                  text.SetDelay(0)
                  t.SetToolTip(text)
               elif name in ['Show','Barycentric']:
                  if x == 'Contour':
                     t = wx.combo.BitmapComboBox(self, size=(100,-1))
                     items=(\
                            ('BindingEnergy','Formation energy'),
                            ('Stability','Stability'),
                            ('ReactionA','Reaction energy per A'),
                            ('ReactionB','Reaction energy per B'),
                            ('ReactionC','Reaction energy per C'),
                            ('New','None'))
                     for ico, m in items:
                        img=getattr(ico20,ico).GetImage().ConvertToBitmap()
                        t.Append(m,img,m)
                  else:
                     t = wx.CheckBox(self,-1,'')
               elif name in ['Map']:
                 # t = wx.ComboBox(self,-1,'',size=(125,-1),choices=self.cmaps)
                  t = wx.combo.BitmapComboBox(self, size=(100,-1))
                  for m in self.cmaps:
                     if hasattr(ico20,'Map'+m):
                        img=getattr(ico20,'Map'+m).GetImage().ConvertToBitmap()
                     else:
                        img=getattr(ico20,'MapTera').GetImage().ConvertToBitmap()
                     t.Append(m,img,m)
               else:
                  t =wx.TextCtrl(self, -1, size=(100,-1)) 
                  tp = ''

                  if x=='Label':
                      if name in ['Offset A','Offset B','Offset C']:
                         tp = 'Label offset format: x"space"y'
                      elif name in ['Label A','Label B','Label C']:
                         tp = 'For subscripts use LaTeX format, e.g. MnO$_2$'
                  if x=='Diagram Title':
                     tp = 'Set title of diagram'
                  if x=='Grid Lines':
                     tp = 'Enter ticks subdivisions separated with "space"'

                  text = wx.ToolTip( tp )
                  text.SetDelay(0)
                  t.SetToolTip(text )
                  
               # 
               # bind here 
               # 
               if name in ['Color x','Color y', 'Color z']:
                  callback='OnColor'
               else:
                  callback='On'+name
               if hasattr(self,callback):
                  if name in ['Size','Width','Resolution']:
                     t.Bind(wx.EVT_SCROLL_CHANGED,getattr(self,callback))
                     t.Bind(wx.EVT_SCROLL,self.OnChangeToolTip)
                  elif name in ['Color', 'Color x','Color y', 'Color z']:
                     t.Bind(wx.EVT_BUTTON,getattr(self,callback))
                  elif name in ['Map','Normalization']:
                     self.Bind(wx.EVT_COMBOBOX, self.OnCombo, t)
               obj=x.replace(' ','')+name.replace(' ','')
#DEBUG
#               print name, obj
#DEBUG
               setattr(self,obj,t)
               sizer.Add(t,0, wx.ALIGN_RIGHT)
            CategoryBox.Add(sizer, 0, wx.CENTER|wx.EXPAND|wx.ALL, 10)
            box.Add(CategoryBox,0, wx.CENTER|wx.EXPAND|wx.ALL, 0)
        #
        # now create items
        #
        makeItem('Diagram Title',('Name',None))
        makeItem('Label',('Show','Label A','Label B','Label C',\
                 'Offset A','Offset B','Offset C'))
        makeItem('Points',('Show','Size','Color'))
        makeItem('Tie Lines',('Show','Width','Color'))
        makeItem('Grid Lines',('Show','Barycentric','Width',\
                  'Color x','Color y','Color z',\
                  'Ticks x','Ticks y','Ticks z'))
        makeItem('Contour',('Show','Resolution','Map','Label'))

        self.SetSizer(box)

    def OnCombo(self, evt):
        bcb = evt.GetEventObject()
        idx = evt.GetInt()
        st  = bcb.GetString(idx)
        cd  = bcb.GetClientData(idx)

    def OnWidth(self,e):
        b = e.GetEventObject()
        e.Skip()

    def OnSize(self,e):
        b = e.GetEventObject()
        e.Skip()

    def OnResolution(self,e):
        b = e.GetEventObject()
        e.Skip()

    def OnChangeToolTip(self,e):
        b = e.GetEventObject()
        b.SetToolTipString(str(b.GetValue()))
        e.Skip()

    def OnColor(self,e):
        """
        This is mostly from the wxPython Demo!
        """
        b = e.GetEventObject()
        dlg = wx.ColourDialog(self)
 
        # Ensure the full colour dialog is displayed, 
        # not the abbreviated version.
        dlg.GetColourData().SetChooseFull(True)
 
        if dlg.ShowModal() == wx.ID_OK:
            data = dlg.GetColourData()
            b.SetBackgroundColour(data.GetColour().Get())
            b.Refresh()
        dlg.Destroy()
    

########################################################################
class ProgressStatusBar( wx.StatusBar ):
    '''Custom StatusBar with a built-in progress bar'''
    def __init__( self, parent, id_=wx.ID_ANY,
                  style=wx.SB_FLAT, name='ProgressStatusBar' ):
        super( ProgressStatusBar, self ).__init__( parent, id_, style, name )

        self._changed = False
        self.busy = False
        self.IamBusy = False
        self.timer = wx.Timer( self )
        self.prog = wx.Gauge( self, style=wx.GA_HORIZONTAL )
#        self.prog.Hide()

        self.SetFieldsCount( 4 )
        self.SetStatusWidths( [-1, 150,150,150] )

        self.Bind( wx.EVT_IDLE, lambda evt: self.__Reposition() )
        self.Bind( wx.EVT_TIMER, self.OnTimer )
        self.Bind( wx.EVT_SIZE, self.OnSize )

    def __del__( self ):
        if self.timer.IsRunning():
            self.timer.Stop()

    def __Reposition( self ):
        '''Repositions the gauge as necessary'''
        if self._changed:
            lfield = self.GetFieldsCount() - 1
            rect = self.GetFieldRect( lfield )
            prog_pos = (rect.x + 2, rect.y + 2)
            self.prog.SetPosition( prog_pos )
            prog_size = (rect.width - 8, rect.height - 4)
            self.prog.SetSize( prog_size )
        self._changed = False

    def OnSize( self, evt ):
        self._changed = True
        self.__Reposition()
        evt.Skip()

    def OnTimer( self, evt ):
        if not self.prog.IsShown():
            self.timer.Stop()

        if self.busy:
            self.prog.Pulse()

    def Run( self, rate=100 ):
        if not self.timer.IsRunning():
            self.timer.Start( rate )

    def GetProgress( self ):
        return self.prog.GetValue()

    def SetProgress( self, val ):
        if not self.prog.IsShown():
            self.ShowProgress( True )

        if val == self.prog.GetRange():
            self.prog.SetValue( 0 )
            self.ShowProgress( False )
        else:
            self.prog.SetValue( val )

    def SetRange( self, val ):
        if val != self.prog.GetRange():
            self.prog.SetRange( val )

    def ShowProgress( self, show=True ):
        self.__Reposition()
        self.prog.Show( show )

    def StartBusy( self, rate=100 ):
        self.busy = True
        self.__Reposition()
        self.ShowProgress( True )
        if not self.timer.IsRunning():
            self.timer.Start( rate )

    def SetBusy( self ):
        self.IamBusy = True

    def UnsetBusy( self ):
        self.IamBusy = False


    def StopBusy( self ):
        self.timer.Stop()
        self.ShowProgress( False )
        self.prog.SetValue( 0 )
        self.busy = False

    def IsBusy( self ):
        return self.busy

########################################################################
class CustomToolbar(wx.Panel):
    def __init__(self, parent):
        from wx.lib import buttons
        from matplotlib.backends.backend_wx import  wxc
        """Constructor"""
        wx.Panel.__init__(self, parent)
        btn_size = (30, 30)
        self.ID = {}
        self.ID_ = {}
        self.items=(
        ('New', 'New', 'New diagram', 'OnNew'), \
        ('Open', 'Open', 'Open diagram', 'OnOpen'), \
        ('Save', 'Save', 'Save current diagram', 'OnSave'), \
        (None, None, None, None), \
        ('Home', 'Home', 'Reset original view', 'OnHome'), \
        ('Back', 'Back', 'Back to previous view', 'OnBack'), \
        ('Forward', 'Forward', 'Forward to next view', 'OnForward'), \
        ('Pan', 'Pan', 'Pan axes with left mouse, zoom with right', 'OnMove'), \
        ('Zoom', 'Zoom', 'Zoom in (out) to rectangle with left (right) mouse', 'OnZoom'), \
        ('Opt', 'Opt', 'Set optimum plot layout', 'OnOpt'), \
        (None, None, None, None), \
        ('Points', 'DataPoints', 'Show data points', 'OnPoints'), \
        ('TieLines', 'DataPoints2', 'Show convex hull', 'OnTieLines'), \
        ('Grid', 'Grid', 'Show grid', 'OnGrid'), \
        (None, None, None, None), \
        ('Binding', 'BindingEnergy', 'Show binding energy', 'OnBinding'), \
        ('Stability', 'Stability', 'Show stability', 'OnStability'), \
        ('Reaction', 'Reaction', 'Show reaction stability', 'OnReaction'), \
        (None, None, None, None), \
        ('Snap', 'SnapCursor', 'Inspect data points', 'OnSnap'), \
        ('LineTool', 'LineTool', 'Start line tool', 'OnLineTool'), \
        ('Toggle3D', 'Toggle3D', 'Show diagram in 3D', 'OnToggle3D'), \
        )
        
        top_sizer = wx.BoxSizer(wx.VERTICAL)
        tb_sizer = wx.BoxSizer(wx.HORIZONTAL)

        for name, icon, tooltip, callback in self.items: 
            if name is None: 
               t = wx.StaticLine(self, -1, style=wx.LI_VERTICAL)
               tb_sizer.Add( t , 1, wx.EXPAND, 1)
               continue
            else: 
               img=getattr(ico,icon).GetImage().ConvertToBitmap()
               self.ID[name] = []
               self.ID_[name] = wx.NewId()
               if name in ['New', 'Open','Save','Back','Forward','Home','Opt']:
                  t = buttons.GenBitmapButton(self,self.ID_[name], 
                     style=wx.NO_BORDER, size=btn_size, bitmap=img)
               else:
                  #t = buttons.GenBitmapToggleButton(self,self.ID_[name], 
                  t = buttons.GenBitmapButton(self,self.ID_[name], 
                     style=wx.NO_BORDER, size=btn_size, bitmap=img)
               # 
               # bind a hover event 
               # 
               wx.EVT_ENTER_WINDOW(t, self.OnWindow)
               wx.EVT_LEAVE_WINDOW(t, self.OffWindow)
               # 
               # set tool information here
               # 
               t.SetToolTip(wx.ToolTip(tooltip))
               setattr(self,name,t)
               tb_sizer.Add(getattr(self,name), 0, wx.ALL, 0)
               #
               # bind button 
               # 
               if hasattr( parent.GetParent(),callback ):
                  self.Bind(wx.EVT_BUTTON, getattr(parent.GetParent(),callback),
                            t)
               # 
               # initialize most of the buttons as disabled, because no diagram
               # is open yet
               # 
               if name not in ['New', 'Open','Save']:
                  t.Disable()
                  t.Refresh()
               # 
               # initialize ID dictionary for states of every page  
               # 
               self.ID[name] = []

        top_sizer.Add(tb_sizer)
        
        h_sizer = wx.BoxSizer(wx.HORIZONTAL)
        h_sizer.Add(wx.StaticLine(self), 1, wx.EXPAND)
        top_sizer.Add(h_sizer, 1, wx.EXPAND)
        
        self.SetSizer(top_sizer)

    def OnWindow(self,event):
        b=event.GetEventObject()
        #
        # use b.SetBezelWidth(1) to set frame 
        #b.SetBackgroundColour( color )
        #
        #color='Red'
        b.SetBezelWidth(1)
        b.Refresh()
        event.Skip()

    def OffWindow(self,event):
        b=event.GetEventObject()
        #color= wx.SystemSettings.GetColour(wx.SYS_COLOUR_BACKGROUND)
        #b.SetBackgroundColour( color )
        b.SetBezelWidth(0)
        b.Refresh()
        event.Skip()   


###############################################################################
class Logger(wx.Panel):
    """ Logger class conatins text field to enter logs """
    #--------------------------------------------------------------------------
    def __init__(self, parent, id=-1 ):
        wx.Panel.__init__(self, parent=parent, id=wx.ID_ANY)
        self.parent = parent
        #
        # add text control field with multilines to write text
        #
        self.text = wx.TextCtrl(self,style = wx.TE_MULTILINE)     
        self.text.SetEditable(False)
	self.text.SetFont(wx.Font(9, wx.MODERN, wx.NORMAL, wx.NORMAL, 0, "Courier"))
        self.SetInitInfo()
        #
        # add sizer 
        #  
        sizer=wx.BoxSizer()
        sizer.Add(self.text, 1, wx.EXPAND|wx.ALL)
        self.SetSizer(sizer)
 
    def Add2Log(self, record):
        self.text.SetEditable(True)
        text = record+'\n'
        self.text.WriteText(text)
        self.text.SetEditable(False)

    def SetInitInfo(self):
        version = self.parent.GetParent().GetParent().GetParent().GetParent().info.GetVersion()
        info = ' TERA version: '+str(version)
        self.Add2Log(info)

###############################################################################
class Commenter(wx.Panel):
    """ Comment class conatins text field to enter comments """
    #--------------------------------------------------------------------------
    def __init__(self, parent, id=-1 ):
        wx.Panel.__init__(self, parent=parent, id=wx.ID_ANY)
        #
        # add text control field with multilines to write text
        #
        self.text = wx.TextCtrl(self,style = wx.TE_MULTILINE)     
        #
        # add sizer 
        #  
        sizer=wx.BoxSizer()
        sizer.Add(self.text, 1, wx.EXPAND|wx.ALL)
        self.SetSizer(sizer)
        
    def Add(self, record):
        print 'adding'



###############################################################################
class TabPanel(wx.Panel):
    """
    This describes the page layout. every page has its own matplotlib figure
    """
    #--------------------------------------------------------------------------
    def __init__(self, parent, Project =None,  id=-1 ):
        from plotter import CustomNavigationToolbar2Wx
        from matplotlib.figure import Figure
        from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg as \
           FigureCanvas

        wx.Panel.__init__(self, parent=parent, id=wx.ID_ANY)
        self.parent=parent
        self.Project = Project
        # *********************************************************************
        sizer = wx.BoxSizer(wx.VERTICAL)
        # *********************************************************************
        #
        # initialize MPL figure canvas
        #
        self.figure = Figure()
        self.figure.add_subplot(111)
        self.canvas = FigureCanvas(self, -1, self.figure)
        sizer.Add(self.canvas, 1, wx.EXPAND)
        #
        # add customized matplotlib toolbar
        #
        self.MPLToolBar = CustomNavigationToolbar2Wx(self.canvas,parent,Project)
        sizer.Add(self.MPLToolBar, 0, wx.EXPAND)
        self.MPLToolBar.Hide()
        #
        #
        # finish sizer
        #
        # *********************************************************************
        self.SetSizer(sizer)
        # *********************************************************************


###############################################################################
class Notebook(aui.AuiNotebook):
    """
    Notebook class, handles the plots 
    """
    #--------------------------------------------------------------------------
    def __init__(self, parent,toolbar):
        """Constructor"""
        aui.AuiNotebook.__init__(self, parent, wx.ID_ANY, style=aui.AUI_NB_MIDDLE_CLICK_CLOSE )
        # 
        # store reference for toolbar 
        # 
        self.toolbar = toolbar 
        # 
        # bind event when closing tab with mouse 
        # 
        self.Bind(aui.EVT_AUINOTEBOOK_PAGE_CLOSE, self.OnCloseTab, self)
        # 
        # bind event when closing tab with mouse 
        # 
        self.Bind(aui.EVT_AUINOTEBOOK_PAGE_CLOSED, self.OnClosedTab, self)
        #
        # bind parent toolbar to MPLToolBar
        #
        self.Bind(aui.EVT_AUINOTEBOOK_PAGE_CHANGED, self.OnChanged, self)


    #--------------------------------------------------------------------------
    def AddTab(self,Project = None,title=None): 
        #
        # Open a tab and send data 
        #
        tab=TabPanel(self, Project = Project, id=self.GetPageCount())

        if Project is None: 
           self.AddPage(tab, 'New diagram '+str(self.GetPageCount()),\
               select=True)
        else: 
           self.AddPage(tab, Project.title, select=True)
           # 
           # eventually enable buttons, if there is at least one tab open
           # 
           pages = self.GetPageCount()
           if pages == 1 : 
              for name, icon, tooltip, callback in self.toolbar.items: 
                 if name in ['New','Open','Save',None]: 
                    continue
                 else:
                    b=getattr(self.toolbar,name)
                    b.Enable()
                    b.Refresh()
           #
           # initialize button states of toolbar, only if not done yet
           #
           for key in self.toolbar.ID.keys():
              if self.GetPageCount() > len(self.toolbar.ID[key]) : 
                 self.toolbar.ID[key].append(False)
           #
           # get mpl toolbar of selected tab and change status of every item in
           # toolbar 
           #
           current_page = pages-1
           SelectedMPLToolBar=self.GetPage(current_page).MPLToolBar
           items = ( ( 'Grid', 'GridLinesCollectionOn' ), \
                     ( 'Points', 'PointsScatterOn' ) ,\
                     ( 'TieLines', 'TieLinesCollectionOn' ) ,\
                     ( 'Binding', 'BindingContourSetOn' ) ,\
                     ( 'Stability', 'StabilityContourSetOn' ) ,\
                     ( 'Reaction', 'ReactionContourSetOn' ) ,\
                    )
           for text, switch  in items:

               if hasattr( Project, switch ) : 
                  if getattr( Project, switch ) : 
                     if text is None:
                         continue
                     t=getattr(self.toolbar,text)
                     self.toolbar.ID[text][current_page] =True 
                     if self.toolbar.ID[text][current_page] : 
                        color= wx.SystemSettings.GetColour(wx.SYS_COLOUR_BTNHILIGHT)
                        t.SetBackgroundColour( color )
                        t.Refresh()


    #--------------------------------------------------------------------------
    def RemoveTab(self): 
        current_page = self.GetSelection()
        if current_page >= 0 : 
           self.DeletePage(current_page)     
        # 
        # eventually disable buttons, if this was the last tab, which was closed
        # 
        pages = self.GetPageCount()
        if pages < 1 : 
           for name, icon, tooltip, callback in self.toolbar.items: 
              if name in ['New','Open','Save',None]: 
                 continue
              else:
                 b=getattr(self.toolbar,name)
                 b.Disable()
                 b.Refresh()

    #--------------------------------------------------------------------------
    def OnCloseTab(self,event): 
        current_page = event.GetSelection()
        MPLToolBar = self.GetPage(current_page).MPLToolBar
        
        for  cid in MPLToolBar.Initcids:
           MPLToolBar.canvas.mpl_disconnect(cid)       
        MPLToolBar.Initcids = [ ]

    #--------------------------------------------------------------------------
    def OnClosedTab(self,event): 
        current_page = event.GetSelection()
        # 
        # eventually disable buttons, if this was the last tab, which was closed
        # 
        pages = self.GetPageCount()
        if pages < 1 : 
           for name, icon, tooltip, callback in self.toolbar.items: 
              if name in ['New','Open','Save',None]: 
                 continue
              else:
                 b=getattr(self.toolbar,name)
                 b.Disable()
                 b.Refresh()

    # -------------------------------------------------------------------------
    def OnChanged(self,e):
        #
        # get current active page number 
        #
        current_page = self.GetSelection()
        #
        # initialize button states of toolbar, only if not done yet
        #
        for key in self.toolbar.ID.keys():
           if self.GetPageCount() > len(self.toolbar.ID[key]) : 
              self.toolbar.ID[key].append(False)
        #
        # write error if nothing open yet, that's probably never the case
        #
        if current_page < 0  : 
           print 'ERROR, no tabs open'
        else:
           #
           # get mpl toolbar of selected tab and change status of every item in
           # toolbar 
           #
           SelectedMPLToolBar=self.GetPage(current_page).MPLToolBar
           for text, icon, tooltip_text, callback in SelectedMPLToolBar.items:
               if text is None :
                   continue
               t=getattr(self.toolbar,text)
               if self.toolbar.ID[text][current_page] : 
                  color= wx.SystemSettings.GetColour(wx.SYS_COLOUR_BTNHILIGHT)
                  t.SetBackgroundColour( color )
                  t.Refresh()
               else:
                  color= wx.SystemSettings.GetColour(wx.SYS_COLOUR_BACKGROUND)
                  t.SetBackgroundColour( color )
                  t.Refresh()

           self.SetProjectFromSettings(current_page)


    def SetProjectFromSettings(self,current_page,default=False):
       from poscar import Data
       #
       # fill settings pane 
       #
       items=( ('Label','LabelA','LabelB','LabelC'),\
                ('Points','Show','Size','Color'),\
                ('TieLines','Show','Width','Color'),\
                ('Label','Show',None,None),\
                ('Label','OffsetA','OffsetB','OffsetC'),\
                ('GridLines','Show','Barycentric','Width'),\
                ('GridLines','Colorx','Colory','Colorz'),\
                ('GridLines','Ticksx','Ticksy','Ticksz'),\
              )
       if default:
           Project = Data()
       else:
           Project = self.GetPage(current_page).Project
       PagePlot=self.GetTopLevelParent().B.PagePlot
       for name, A, B, C in items: 
          #
          # get values  in settings pane
          #
          ObjA = getattr(PagePlot,name+A)
          if B is None:
             ObjB = None
          else:
             ObjB = getattr(PagePlot,name+B)
          if C is None:
             ObjC = None
          else:
             ObjC = getattr(PagePlot,name+C)
          if name == 'Label' : 
             if A == 'Show':
                ObjA.SetValue(Project.LabelsOn)
             elif A == 'OffsetA':
                oa = Project.LabelAOffset
                ob = Project.LabelBOffset
                oc = Project.LabelCOffset
                ObjA.SetValue(str(oa[0])+' '+str(oa[1]))
                ObjB.SetValue(str(ob[0])+' '+str(ob[1]))
                ObjC.SetValue(str(oc[0])+' '+str(oc[1]))
             else:
                if default:
                   ObjA.SetValue('')
                   ObjB.SetValue('')
                   ObjC.SetValue('')
                else:
                   ObjA.SetValue(Project.pdb['C']["name"])
                   ObjB.SetValue(Project.pdb['A']["name"])
                   ObjC.SetValue(Project.pdb['B']["name"])
          if name == 'Points' : 
             ObjA.SetValue(Project.PointsScatterOn)
             ObjB.SetValue(Project.ps)
             ObjC.SetBackgroundColour(Project.pt_color)
          if name == 'TieLines' : 
             ObjA.SetValue(Project.TieLinesCollectionOn)
             ObjB.SetValue(Project.lw)
             ObjC.SetBackgroundColour(Project.tl_color)
          if name == 'GridLines' : 
             if A == 'Show':
                ObjA.SetValue(Project.GridLinesCollectionOn)
                ObjB.SetValue(Project.UseUniform)
                ObjC.SetValue(Project.lwGrid)
             if A == 'Colorx':
                ObjA.SetBackgroundColour(Project.GridColors[0])
                ObjB.SetBackgroundColour(Project.GridColors[1])
                ObjC.SetBackgroundColour(Project.GridColors[2])
             if A == 'Ticksx':
                xi=''
                yi=''
                zi=''
                if (Project.zi is not None):
                   zi = " ".join(str(x) for x in Project.zi)
                if (Project.yi is not None):
                   yi = " ".join(str(x) for x in Project.yi)
                if (Project.xi is not None):
                   xi = " ".join(str(x) for x in Project.xi)
                ObjA.SetValue(xi)
                ObjB.SetValue(yi)
                ObjC.SetValue(zi)
                
       #
       # diagram title 
       #
       PagePlot.DiagramTitleName.SetValue(Project.title)
       #
       # contour
       #
       if Project.BindingContourSetOn:
          PagePlot.ContourShow.SetValue('Formation energy')
       elif Project.StabilityContourSetOn:
          PagePlot.ContourShow.SetValue('Stability')
       elif Project.ReactionContourSetOn:
          text='Reaction energy'
          if Project.NormIDX == 0:
             text+=' per B'
          elif Project.NormIDX == 1:
             text+=' per C'
          else:
             text+=' per A'
          PagePlot.ContourShow.SetValue(text)
       else:
          PagePlot.ContourShow.SetValue('None')
          
       PagePlot.ContourMap.SetValue(PagePlot.cmaps[Project.CmdIDX])
       PagePlot.ContourLabel.SetValue(Project.CBLabel)
       PagePlot.ContourResolution.SetValue(Project.nsample)

       #
       # preferences
       #
       PagePref=self.GetTopLevelParent().B.PagePref
       PagePref.VESTApath.SetValue( Project.pathVESTA )
       PagePref.Outputdirectory.SetValue( Project.OutDir )



###############################################################################
class PanelA(wx.Panel):
    """ Panel A: contains Notebook with MPL canvas"""
    #----------------------------------------------------------------------
    def __init__(self, parent, toolbar ):
        wx.Panel.__init__(self, parent)
# DEBUG:
        #self.SetBackgroundColour('blue')
# DEBUG:
        #
        # Add NoteBook here
        # 
        sizer=wx.BoxSizer()
        self.NoteBook = Notebook(self,toolbar)       
        sizer.Add(self.NoteBook, 1, wx.EXPAND|wx.ALL)
        self.SetSizer(sizer)

###############################################################################
class PanelB(wx.Panel):
    """ Panel B: contains settings notebook """
    #----------------------------------------------------------------------
    def __init__(self, parent):
        wx.Panel.__init__(self, parent)
# DEBUG:
        #self.SetBackgroundColour('green')
# DEBUG:
        nb = wx.Notebook(self)
        #
        # add plot tab
        #
        self.PagePlot = PlotSettings(nb)
        #
        # add preference tab
        #
        self.PagePref = PreferencesSettings(nb)
        #self.PagePref = TestPanel(nb)
        #
        # add those tabs to the notebook
        # 
        nb.AddPage(self.PagePlot,'Appearance')
        nb.AddPage(self.PagePref,'Preferences')
        #
        # add sizer 
        #  
        sizer=wx.BoxSizer(wx.VERTICAL)
        sizer.Add(nb, 1, wx.EXPAND|wx.ALL)

        #
        # add ok, reset and cancel button
        #
        hbox = wx.BoxSizer(wx.HORIZONTAL)
        self.IsDefaultOn = wx.CheckBox(self,-1,'Save as default')
        hbox.Add(self.IsDefaultOn,0,wx.ALIGN_LEFT)
        sizer.Add(hbox,0,wx.ALIGN_LEFT)
      
        hbox = wx.BoxSizer(wx.HORIZONTAL)
        o = wx.Button(self, wx.ID_OK, "OK")
        o.Bind(wx.EVT_BUTTON,self.OnOK)
        c = wx.Button(self, wx.ID_CANCEL, "Cancel")
        c.Bind(wx.EVT_BUTTON,self.OnCancel)
        hbox.Add(o,0,wx.CENTER)
        hbox.Add(c,0,wx.CENTER)
        sizer.Add(hbox,0,wx.CENTER)

        self.SetSizer(sizer)
        #      
        # bind event when on focus     
        #      

    def OnOK(self,e):
       from dialog import ShowStatusErrorMessage
       IsOK = self.Apply()
       #      
       # replot if everything is OK      
       #      
       self.SetCursor( wx.StockCursor(wx.CURSOR_WAIT) )
       self.GetTopLevelParent().MainPanel.SetCursor( wx.StockCursor(wx.CURSOR_WAIT) )

       if IsOK is None:
          nb = self.GetTopLevelParent().A.NoteBook
          current_page = nb.GetSelection()
          if current_page >=0:
             nb.GetPage(current_page).MPLToolBar.ClearFigure()
             nb.GetPage(current_page).MPLToolBar.InitPlot()
             if nb.GetPage(current_page).MPLToolBar.Project.SnapOn:
                nb.GetPage(current_page).MPLToolBar.OnSnap()
             if nb.GetPage(current_page).MPLToolBar.Project.LineToolOn:
                nb.GetPage(current_page).MPLToolBar.OnLineTool()
       else:
          #
          # in this case IsOK is the error id displayed in the statusbar
          #
          ShowStatusErrorMessage(self.GetTopLevelParent().statusbar,IsOK)

       self.IsDefaultOn.SetValue(False)
       self.SetCursor( wx.StockCursor(wx.CURSOR_ARROW) )
       self.GetTopLevelParent().MainPanel.SetCursor(
           wx.StockCursor(wx.CURSOR_ARROW) )

    def OnCancel(self,e):
       from dialog import ReadWriteCfgFile
       nb = self.GetTopLevelParent().A.NoteBook 
       current_page = nb.GetSelection()
       if current_page < 0 :
          return
       else:
          self.GetTopLevelParent().statusbar.SetStatusText('Defaults restored')
          nb.SetProjectFromSettings(current_page,default=True)
       self.IsDefaultOn.SetValue(False)

    def Apply(self):
       from dialog import IsFloat
       from color import maps
       from dialog import ReadWriteCfgFile, ShowStatusWarningMessage
       from os import path 

       def RGB2HEX(x): 
          new = ( max( 0, min(x[0],255)),\
                  max( 0, min(x[1],255)),\
                  max( 0, min(x[2],255)))
          return "#{0:02x}{1:02x}{2:02x}".format(new[0], new[1], new[2])
       #
       # fill settings pane 
       #
       items=( ('Label','LabelA','LabelB','LabelC'),\
                ('Label','Show',None,None),\
                ('Label','OffsetA','OffsetB','OffsetC'),\
                ('Points','Show','Size','Color'),\
                ('TieLines','Show','Width','Color'),\
                ('GridLines','Show','Barycentric','Width'),\
                ('GridLines','Colorx','Colory','Colorz'),\
                ('GridLines','Ticksx','Ticksy','Ticksz'),\
              )
       nb = self.GetTopLevelParent().A.NoteBook
       current_page = nb.GetSelection()
       if current_page <0:
          return 'No data imported'
       else:
          Project = nb.GetPage(current_page).Project
          for name, A, B, C in items: 
             #
             # get values  in settings pane
             #
             ObjA = getattr(self.PagePlot,name+A)
             if B is None:
                ObjB = None
             else:
                ObjB = getattr(self.PagePlot,name+B)
             if C is None:
                ObjC = None
             else:
                ObjC = getattr(self.PagePlot,name+C)
             if name == 'Label' : 
                if A == 'LabelA':
                   Project.LabelA = ObjA.GetValue()
                   Project.LabelB = ObjB.GetValue()
                   Project.LabelC = ObjC.GetValue()
                if A == 'Show':
                   Project.LabelsOn = ObjA.GetValue()
                elif A == 'OffsetA':
                   oa = str(ObjA.GetValue()).split()
                   ob = str(ObjB.GetValue()).split()
                   oc = str(ObjC.GetValue()).split()
                   if IsFloat(oa[0]) and IsFloat(oa[1]):
                      Project.LabelAOffset=(float(oa[0]),float(oa[1]))
                   else:
                      return  'Cannot read Offset A'

                   if IsFloat(ob[0]) and IsFloat(ob[1]):
                      Project.LabelBOffset=(float(ob[0]),float(ob[1]))
                   else:
                      return  'Cannot read Offset B'

                   if IsFloat(oc[0]) and IsFloat(oc[1]):
                      Project.LabelCOffset=(float(oc[0]),float(oc[1]))
                   else:
                      return  'Cannot read Offset C'

             if name == 'Points' : 
               Project.PointsScatterOn = ObjA.GetValue()
               Project.ps = ObjB.GetValue()
               a = RGB2HEX(ObjC.GetBackgroundColour().Get())
               Project.pt_color = a
               
             if name == 'TieLines' : 
               Project.TieLinesCollectionOn= ObjA.GetValue()
               Project.lw = ObjB.GetValue()
               a = RGB2HEX(ObjC.GetBackgroundColour().Get())
               Project.tl_color = a
             if name == 'GridLines' : 
                if A == 'Show':
                  Project.GridLinesCollectionOn = ObjA.GetValue()
                  Project.UseUniform = ObjB.GetValue()
                  Project.lwGrid = ObjC.GetValue()
                if A == 'Colorx':
                  a = ObjA.GetBackgroundColour().Get()
                  b = ObjB.GetBackgroundColour().Get()
                  c = ObjC.GetBackgroundColour().Get()
                  x,y,z = RGB2HEX(a), RGB2HEX(b), RGB2HEX(c) 
                  Project.GridColors = (x,y,z)
                  #
                  # set the same colors also for the label
                  #
                  Project.LabelColor['A']=x
                  Project.LabelColor['B']=y
                  Project.LabelColor['C']=z
                if A == 'Ticksx':
                   xi = str(ObjA.GetValue()).split()
                   yi = str(ObjB.GetValue()).split()
                   zi = str(ObjC.GetValue()).split()
               
                   a=[]
                   for i in xi: 
                      if IsFloat(i) : 
                         a.append(float(i))
                      else:
                         return 'Cannot read Ticks x'
                   if len(a)>0:
                      Project.xi=a
                   else:
                      Project.xi=None
                   a=[]
                   for i in yi: 
                      if IsFloat(i) : 
                         a.append(float(i))
                      else:
                         return 'Cannot read Ticks y'
                   if len(a)>0:
                      Project.yi=a
                   else:
                      Project.yi=None
                   a=[]
                   for i in zi: 
                      if IsFloat(i) : 
                         a.append(float(i))
                      else:
                         return 'Cannot read Ticks z'
                   if len(a)>0:
                      Project.zi=a
                   else:
                      Project.zi=None
          #
          # diagram title 
          #
          Project.title = self.PagePlot.DiagramTitleName.GetValue()
          nb.SetPageText(current_page,Project.title)
          #
          # contour selection
          #
          i = self.PagePlot.ContourShow.GetValue()
          if i=='Formation energy':
             Project.BindingContourSetOn = True
             Project.StabilityContourSetOn = False
             Project.ReactionContourSetOn = False
          elif i=='Stability':
             Project.StabilityContourSetOn = True
             Project.ReactionContourSetOn = False
             Project.BindingContourSetOn = False
          elif i.split()[0] == 'Reaction' :
             if i.split()[3]=='A':
                Project.NormIDX = 2
             elif i.split()[3]=='B':
                Project.NormIDX = 0
             else:
                Project.NormIDX = 1
             Project.ReactionContourSetOn = True
             Project.StabilityContourSetOn = False
             Project.BindingContourSetOn = False
          else:
             Project.ReactionContourSetOn = False
             Project.StabilityContourSetOn = False
             Project.BindingContourSetOn = False
          #
          # color map
          #
          y = self.PagePlot.ContourMap.GetValue()
          def GetPosInList(y,x):
             j=-1
             for i in x : 
               j+=1
               if i==y : 
                  return j
             return 0
          Project.CmdIDX = GetPosInList(y,self.PagePlot.cmaps)
          Project.Cmd = getattr(maps(),maps().AllMapNames[Project.CmdIDX])
          Project.CmdIntersect = self.PagePlot.ContourResolution.GetValue()
          Project.CBLabel = self.PagePlot.ContourLabel.GetValue()
          
          Project.nsample = Project.CmdIntersect/2
          #
          # restore  map 
          #
          nb.GetPage(current_page).Project = Project

          #
          # set preferences 
          #
          items=(('VESTApath' , 'pathVESTA', 'VESTA execturable'),\
                 ('Outputdirectory', 'OutDir','Output directory' ))

          for i,j,info in items:
            f = str( getattr(self.PagePref,i).GetValue() )
            if path.exists( f ) and len(f)>0:
               setattr(Project,j, f)
            else:
               text="{:} {:} not found!".format(info,f)
               ShowStatusWarningMessage(self.GetTopLevelParent().statusbar,\
                  text,delay=1)
               setattr(Project,j, '' )
          #
          # if this is new default save to file 
          #
          if self.IsDefaultOn.GetValue(): 
             ReadWriteCfgFile(write=True,p=Project)
             text=' New default settings saved successfully '
             self.GetTopLevelParent().C.PageLog.Add2Log( text )
             
          return None

###############################################################################
class PanelC(wx.Panel):
    """ Panel C: contains logger tab """
    #----------------------------------------------------------------------
    def __init__(self, parent):
        wx.Panel.__init__(self, parent)
# DEBUG:
        #self.SetBackgroundColour('green')
# DEBUG:
        nb = wx.Notebook(self,style=wx.BK_BOTTOM)
        #
        # add log tab
        #
        self.PageLog = Logger(nb)
        #
        # add comment tab
        #
        self.PageComment = Commenter(nb)
        #
        # add those tabs to the notebook
        # 
        nb.AddPage(self.PageLog,'Log')
        nb.AddPage(self.PageComment,'Comments')
        #
        # add sizer 
        #  
        sizer=wx.BoxSizer()
        sizer.Add(nb, 1, wx.EXPAND|wx.ALL)
        self.SetSizer(sizer)
