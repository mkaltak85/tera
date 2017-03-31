import wx 
import wx.html
from Documentation import documentation

class wxHTML(wx.html.HtmlWindow):
     def OnLinkClicked(self, link):
         parent = self.GetParent()
         leftP  = parent.GetWindow1()
         rightP = parent.GetWindow2()
         key = link.GetHref()
	 if not key in documentation :  
            text='Entry "{:}" not found'.format(key)
            dlg = wx.MessageDialog(self,text, 'Error',wx.OK|wx.ICON_ERROR)
            result = dlg.ShowModal()
            return
         else:
            rightP.SetPage( documentation[ key ] )
         # add entry into Browser History
	 History = parent.GetParent().GetParent().BrowserHistory
	 CurrentPage = parent.GetParent().GetParent().CurrentPage 

         if not History[ CurrentPage ] == key : 
             # delete all entries up to current page
	     for i in range( len( History)-1, CurrentPage,-1):
                 del History[ -1 ] 
	     CurrentPage = len( History ) - 1
	     CurrentPage += 1
	     History.append( '{:}'.format( key ) ) 
	     parent.GetParent().GetParent().BrowserHistory = History
	     parent.GetParent().GetParent().CurrentPage = CurrentPage

	 print History, parent.GetParent().GetParent().CurrentPage

class HelpDialog(wx.Dialog):
    def __init__(self, *args, **kw):
        from icons25 import Tera as ico_Tera
        super(HelpDialog, self).__init__(*args, **kw) 

        self.SetIcon(ico_Tera.GetIcon())
        panel = wx.Panel(self,-1)
            
        self.InitHelpWindow(panel)
        self.SetMinSize((400,345))

        self.BrowserHistory=[ self.home ]
        self.CurrentPage = 0

        
    #--------------------------------------------------------------------------
    # initialize window
    #--------------------------------------------------------------------------
    def InitHelpWindow(self,panel):
        from widgets import ProportionalSplitter
        mainsizer = wx.BoxSizer(wx.VERTICAL)
 
        self.InitToolBar(panel)
	mainsizer.Add(self.toolbar, 0,wx.ALL,0 )

        self.vSplitter = ProportionalSplitter(panel, 
             id=-1, proportion = 0.25, style = wx.SP_BORDER) 
        self.InitLeftPanel()
        self.InitRightPanel()
        self.vSplitter.SplitVertically(self.leftP, self.rightP)
        self.vSplitter.SetSashGravity(0.5)
        mainsizer.Add(self.vSplitter , 1, wx.EXPAND)

        panel.SetSizerAndFit(mainsizer)

    #--------------------------------------------------------------------------
    # Init toolbar
    #--------------------------------------------------------------------------
    def InitToolBar(self,panel):
        import icons25 as ico
        toolitems=(
        ('ID_HOME', 'Home', 'Go to home', 'OnHome'), \
        ('ID_UNDO', 'Back', 'Go back', 'OnBack'), \
        ('ID_REDO', 'Forward', 'Go forward', 'OnForward'), \
        )
        self.toolbar = wx.ToolBar(panel)
	for idname, icon, tip, c in toolitems: 
           img=getattr(ico,icon).GetImage().ConvertToBitmap()
	   id = getattr( wx, idname )
	   callb = getattr( self, c )
           t = self.toolbar.AddLabelTool(id, '', img,shortHelp=tip)

	   self.Bind(wx.EVT_TOOL, callb, t)
           self.toolbar.Realize()

    #--------------------------------------------------------------------------
    # display catalogue
    #--------------------------------------------------------------------------
    def InitRightPanel(self):
        mainsizer = wx.BoxSizer()
        self.rightP = wxHTML(self.vSplitter)
	#self.home="Introduction"
	self.home="Tutorial"
        self.rightP.SetPage( documentation[ self.home ]  )
        self.leftP.SetSizerAndFit(mainsizer)

    #--------------------------------------------------------------------------
    # displays site 
    #--------------------------------------------------------------------------
    def InitLeftPanel(self):
        mainsizer = wx.BoxSizer()
        # create catalouge 
        catalogue=''
        for key in documentation.keys():
            entry = key.replace("_", " ")
            catalogue+='\n<p><b><a href="{:}">{:}</a></b></p>'.format( key,entry )

        self.leftP = wxHTML(self.vSplitter)
        self.leftP.SetPage( catalogue )
        self.leftP.SetSizerAndFit(mainsizer)

    def OnHome(self,event):
        if not self.BrowserHistory[ self.CurrentPage ] == self.home : 
           # delete all entries up to current page
	   for i in range( len( self.BrowserHistory)-1, self.CurrentPage,-1):
               del self.BrowserHistory[ -1 ] 
	   self.CurrentPage = len( self.BrowserHistory ) - 1
           self.CurrentPage += 1 
           self.BrowserHistory.append( self.home )
           self.rightP.SetPage( documentation[ self.home ]  )
	print 'HOME', self.BrowserHistory, self.CurrentPage

    def OnBack(self,event):
	if self.CurrentPage > 0 : 
           self.CurrentPage -= 1 
	   goto = self.BrowserHistory[ self.CurrentPage ] 
           self.rightP.SetPage( documentation[ goto ]  )
	print 'BACK',self.BrowserHistory, self.CurrentPage


    def OnForward(self,event):
	if len( self.BrowserHistory )-1 > self.CurrentPage : 
           self.CurrentPage += 1 
	   goto = self.BrowserHistory[ self.CurrentPage ] 
           self.rightP.SetPage( documentation[ goto ]  )
	print 'FORWARD',self.BrowserHistory, self.CurrentPage

