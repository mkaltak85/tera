import wx 

#------------------------------------------------------------------------------
def IsFloat(value):
  try:
    float(value)
    return True
  except ValueError:
    return False

#------------------------------------------------------------------------------
def get_extension(filename):
    """ gets extension of filename """
    from os import path
    basename = path.basename(filename)  # os independent

#    ext = '.'.join(basename.split('.')[1:])
#    return ext if ext else None

    all=basename.split('.')
    if len( all ) > 1 : 
       ext = basename.split('.')[len(all)-1]
    else:
       ext = None
    return ext

#------------------------------------------------------------------------------
def GetRandomData():
    """ Writes a random table """
    from numpy.random import randint
    from numpy.random import random
    from numpy import zeros,argmin
    n=200
    nb=n/10
    a=-100
    b=-1
    I=10
    s = 2
    def g(x,i):
       from numpy import pi,exp,append,sum
       from poscar import small
       xb=[]
       a=0.5
       b=1.5
       x0=0.3
       for j in range(0,3):
          xb = append(xb, x[j]/sum( x ) )
          if abs(x[j])<small:
            x0=0.5

       f = exp( -(xb[i]-x0)**2 )*\
          ( (b-a)*random()+a)
       return s*f

    # energies
    e = (b-a)*random(3) + a
    data = zeros( ( n+3+3*nb, 4 ))

    data[0,0] = 1.
    data[0,1] = 0.
    data[0,2] = 0.
    data[0,3] = e[0]
    data[1,0] = 0.
    data[1,1] = 1.
    data[1,2] = 0.
    data[1,3] = e[1]
    data[2,0] = 0.
    data[2,1] = 0.
    data[2,2] = 1.
    data[2,3] = e[2]
    #
    # interior 
    #
    for i in range(0,n):
       for j in range(0,3):
          data[i+3,j] = float(randint(1,I))/float(randint(1,I))
       for j in range(0,3):
          r = g(data[i+3,:3],j)
          data[i+3,3] = data[i+3,3] + e[j]*data[i+3,j]*r

    # now add some points on the borders, but avoid endpoints
    for k in [ 0,1,2]:
       for i in range(n+3+k*nb,n+3+(k+1)*nb):
          for j in range(0,3):
             if k != j : 
                data[i,j] = float(randint(1,I))/float(randint(1,I))
          for j in range(0,3):
             if k != j : 
                r = g(data[i,:3],j)
                data[i,3] = data[i,3] + e[j]*data[i,j]*r
    return data 

#------------------------------------------------------------------------------
def ReadDataTable(f,species,nion,statusbar,IsRandom,log):
    """  loads endpoint information and raw txt data from file  """
    from poscar import Data,UpdateEndpoints,IntegerStoichiometry
    from numpy import loadtxt,array
    from poscar import PHASES,CMAX,is_number,GetLCM
    
    def argsort(seq):
       return sorted(range(len(seq)), key=seq.__getitem__)
   

    sx = UpdateEndpoints(species,nion,log=log)
    #
    # read data 
    #
    if IsRandom:
       data = GetRandomData()
    else:
       data = loadtxt(f,comments='#')

    a = array(sx.pdb["phase_matrix"])

    # import data into sx.db
    for i in range(0,len(data)):
       xb=[]
       for j in range(0,3):
          xb.append(data[i, j ] )
       si=IntegerStoichiometry(xb,sx)[0]
       # 
       # number of ions
       # 
       n = a.dot( array( si ) ).tolist()
       #
       SPECIES = []
       NION = [ ]
       for j in range(0,len(n) ) :
           if n[j] > 0:
              NION.append( int( n[j] ) ) 
              SPECIES.append(sx.pdb["EndpointAtoms"][j])

       NSPECIES = len(SPECIES)
       NIONS = int(sum(NION))

       if NIONS > 0:
          #build id of structure 
          IDKEY=str(NIONS)+'-'+str(NSPECIES)+':'
          #UNIQUE is True if it is a unique configuration for this stoichiometry
          UNIQUE=True
          for j in range(0,NSPECIES):
             IDKEY+=str(SPECIES[j])+'.'+str(NION[j])+'.'
          #look if database contains an entry with the same name 
          #if so, create a new IDKEY 
          if IDKEY in sx.db : 
             #and change uniqueness to false
             sx.db[IDKEY]["IsUnique"]=False
             sx.dbF[IDKEY]["IsUnique"]=False
             UNIQUE = False
             for j in range(1,CMAX):
                if IDKEY+'c-'+str(j)+'.' not in sx.db : 
                   IDKEY+='c-'+str(j)+'.'
                   break
             if j>CMAX-1:
                sys.exit('Error: maximum number of distinct configurations for one'+\
                   ' stoichiometry reached!')

          #
          # get multiplicity for energy based on integer stoichiometry and
          # barycentric coordinates
          #
          Xb = []
          Xin = []
          Xj = data[i,:3] 
          for j in [0,1,2] :
             Xb.append( float( Xj[ j] )/float(sum( data[i,:3] )))
             Xin.append( Xj[ j]  ) 

          alpha = GetLCM( Xin ) 

          if is_number( data[i,3] ):
             try:
                #store energy per pseudo-ion 
                EN=float(data[i,3])*alpha
             except ValueError:
                text='Endpoint B missing'
                ShowStatusErrorMessage(statusbar,text)
                return None

          #add dictionary with IDKEY-Key to database
          sx.db[IDKEY]={}
          #name of configuration
          sx.db[IDKEY]["name"]=IDKEY
          #add number of ions 
          sx.db[IDKEY]["nions"]=NIONS
          #add number of ions per species
          sx.db[IDKEY]["nion"]=NION
          #add number of species
          sx.db[IDKEY]["nspecies"]=NSPECIES
          #add names of species
          sx.db[IDKEY]["species"]=SPECIES
          #add information of uniquess of entry 
          sx.db[IDKEY]["IsUnique"]=UNIQUE
          #finally store energy 
          sx.db[IDKEY]["energy"] = EN
          sx.db[IDKEY]["maxForce"] = [ 0.,0.,0.]
          sx.db[IDKEY]["degeneracy"]=1
          sx.db[IDKEY]["path"]=''

          #add dictionary with IDKEY-Key to database
          sx.dbF[IDKEY]={}
          sx.dbF[IDKEY]["name"]=IDKEY
          sx.dbF[IDKEY]["nions"]=NIONS
          sx.dbF[IDKEY]["nion"]=NION
          sx.dbF[IDKEY]["nspecies"]=NSPECIES
          sx.dbF[IDKEY]["species"]=SPECIES
          sx.dbF[IDKEY]["IsUnique"]=UNIQUE
          sx.dbF[IDKEY]["energy"] = EN
          sx.dbF[IDKEY]["maxForce"] = [ 0.,0.,0.]
          sx.dbF[IDKEY]["degeneracy"]=1
          sx.dbF[IDKEY]["path"]=''
    return sx 

#------------------------------------------------------------------------------
def SaveRawDataFile(project,fname):
    """  loads database and endpoint information from file  """
    from numpy import savetxt, zeros


    npoints=len(project.db.keys())
    data = zeros( shape=( npoints, 4 ) )

    k=-1
    for key in project.db.keys():
        k+=1
        for i in [ 0, 1, 2 ]:
           data[k,i] = project.db[ key ]["bary-coord"][i]
        data[k,3] = project.db[ key ]["formationE"]

    a=project.pdb['A']["name"]
    b=project.pdb['B']["name"]
    c=project.pdb['C']["name"]
    Header='{:} \nEndpoints: {:}   {:}    {:}'.format(
            project.title,a,b,c, )
    savetxt( fname, data, fmt="%20.10f",header=Header )

#------------------------------------------------------------------------------
def ReadTeraFile(filename):
    """  loads database and endpoint information from file  """
    import pickle
    from poscar import Data
 
    sx = Data()

    with open(filename, 'r') as input:
        sx.SaveLoad = pickle.load(input)
        #
        # load full file 
        #
        for a in sx.SaveLoad : 
           setattr(sx,a, pickle.load(input))
    return sx 


#------------------------------------------------------------------------------
def SaveTeraFile(sx, filename ):
    """ saves database and endpoint information to a file """
    import pickle

    with open(filename, 'wb') as output:
        #
        # and the list of attributes
        #
        pickle.dump(sx.SaveLoad,output)
        #
        # save selected attributes 
        #
        for  a in sx.SaveLoad : 
           pickle.dump(getattr(sx,a), output)

#----------------------------------------------------------------------
def SaveProjectAs(frame, nb, logger,statusbar):
    """
    Create and show the Save FileDialog
    """
    from os import getcwd
    wildcard = "TERA diagram (*.tera,*.ter)|*.tera;*.ter|" \
        "Raw text format (*.dat,*.txt)|*.dat;*.txt|"\
        "Joint Photographic Experts Group (*.jpeg,*.jpg)|*.jpeg;*.jpg|"\
        "Tagged Image File Format (*.tiff,*.tif)|*.tiff;*.tif|"\
        "Portable Network Graphics (*.png)|*.png|"\
        "Encapsulated Postscript Vector (*.eps)|*.eps|"\
        "Portable Document Format (*.pdf)|*.pdf|"\
        "Scalable Vector Graphics (*.svg)|*.svg"

    # store for later
    EXT = ["tera","ter","dat","txt","jpeg","jpg","tiff","tif","png","eps","svg","pdf"]
    currentDirectory = getcwd()

    dlg = wx.FileDialog(
        frame, message="Save file as ...", 
        defaultDir=currentDirectory, 
        defaultFile="", wildcard=wildcard, 
        style=wx.SAVE|wx.FD_OVERWRITE_PROMPT
        )
    if dlg.ShowModal() == wx.ID_OK:
        path = dlg.GetPath()
        # 
        # find extension
        # 
        ext=get_extension(path)
        # 
        # add extension from current dialog choice 
        # 
        if ext is None: 
           path+='.'+EXT[ dlg.GetFilterIndex() ]
           ext =  EXT[ dlg.GetFilterIndex() ]
        # 
        # get data from currently active notebook pad 
        # 
        current_page = nb.GetSelection()
        text = ''
        if current_page >= 0 : 
           project = nb.GetPage(current_page).Project 
           # 
           # decide how to write data 
           # 
           if ext == 'tera' or ext =='ter': 
              SaveTeraFile( project, path)

           elif ext == 'jpeg' or ext =='jpg' or \
                ext == 'tiff' or ext =='tif' or\
                ext == 'png' or\
                ext == 'pdf' or\
                ext == 'eps' or\
                ext == 'svg' :
              fig = nb.GetPage(current_page).canvas.figure
              fig.savefig(path)
           else:
              SaveRawDataFile( project, path)

           text ='          Project saved as: {:}'.format(path)
           logger.Add2Log( text )
        else: 
           print 'ERROR: page not open'
    dlg.Destroy()

########################################################################
class ImportForm(wx.Dialog):
    def __init__(self, *args, **kw ):
        from poscar import Data
        import os 
	from sys import argv
        from icons25 import Tera as ico_Tera
        super(ImportForm, self).__init__(*args, **kw) 

        panel = wx.Panel(self, -1)
        self.SetIcon(ico_Tera.GetIcon())
        self.parent = self.GetParent()
        self.statusbar = self.parent.statusbar
        self.currentDirectory = os.getcwd()
        #
        # initialize empty list for imported projects
        #
        self.Projects = [  ] 
        #
        # store an instance for the logger
        #
        self.logger = self.parent.C.PageLog
        #
        # intialize UI
        #
        self.InitUI(panel)
        #
        # bind closing event to cancel callback
        #
        self.Bind(wx.EVT_CLOSE, self.OnCancel)
        
        
        self.SetMaxSize((400,375))
        self.SetMinSize((400,375))

        self.OpenArguments()

    # --------------------------------------------------------------------------
    def InitUI(self,panel):
        # ********************************************************************** 
        mainsizer = wx.BoxSizer(wx.VERTICAL)
        # ********************************************************************** 
        #
        # file control entry 
        #
        hbox = wx.BoxSizer(wx.HORIZONTAL)
        b = wx.Button(panel, label="Browse...")
        label=wx.StaticText(panel, -1, "File:" , size=(60,b.GetDefaultSize()[1]), 
           style=wx.ALIGN_RIGHT )
        hbox.Add(label,0, wx.ALL, 5)
        self.filectrl = wx.TextCtrl(panel, -1, size=b.GetDefaultSize())
        hbox.Add(self.filectrl, 1, wx.ALL, 5)
        hbox.Add(b, 0, wx.ALL, 5)
        mainsizer.Add(hbox, 0, wx.EXPAND, 1)
        b.Bind(wx.EVT_BUTTON, self.OnOpenFile)
        tip='path to TERA file(s)\nmultiple files are separated by "space"'
        self.filectrl.SetToolTip(wx.ToolTip(tip))
	tip='Browse for TERA file(s)'
        b.SetToolTip(wx.ToolTip(tip))
        #
        # directory control entry 
        #
        hbox = wx.BoxSizer(wx.HORIZONTAL)
        label=wx.StaticText(panel, -1, "Directory:" , size=(60,b.GetDefaultSize()[1]),
           style=wx.ALIGN_RIGHT ) 
        hbox.Add(label,0, wx.ALL, 5)
        b = wx.Button(panel, label="Browse...")
        self.dirctrl = wx.TextCtrl(panel, -1, size=b.GetDefaultSize())
        hbox.Add(self.dirctrl, 1, wx.ALL, 5)
        hbox.Add(b, 0, wx.ALL, 5)
        mainsizer.Add(hbox, 0, wx.EXPAND, 1)
        b.Bind(wx.EVT_BUTTON, self.OnOpenDir)
	tip='directory with POSCAR structure files\neach file must contain energy of structure\n"Energy : 0.000 unit"'
        self.dirctrl.SetToolTip(wx.ToolTip(tip))
	tip='Browse for directory'
        b.SetToolTip(wx.ToolTip(tip))
        #
        # Top row with entry C
        #
        hbox = wx.BoxSizer(wx.HORIZONTAL)
        label=wx.StaticText(panel, -1, "   " ) 
        hbox.Add(label,0, wx.ALL, 5)
        label=wx.StaticText(panel, -1, " Phase A " ) 
        hbox.Add(label, 1, wx.EXPAND, 1)
        label=wx.StaticText(panel, -1, " Phase B " ) 
        hbox.Add(label, 1, wx.EXPAND, 1)
        label=wx.StaticText(panel, -1, " Phase C " ) 
        hbox.Add(label, 1, wx.EXPAND, 1)
        mainsizer.Add(hbox, 0, wx.EXPAND, 1)

        hbox = wx.BoxSizer(wx.HORIZONTAL)
        label=wx.StaticText(panel, -1, "   " ) 
        hbox.Add(label,0, wx.ALL, 5)
        self.Cctrl = wx.TextCtrl(panel, -1, size=(60,b.GetDefaultSize()[1]),)
        hbox.Add(self.Cctrl, 1, wx.EXPAND, 1)
        self.Actrl = wx.TextCtrl(panel, -1, size=b.GetDefaultSize())
        hbox.Add(self.Actrl, 1, wx.EXPAND, 1)
        self.Bctrl = wx.TextCtrl(panel, -1, size=b.GetDefaultSize())
        hbox.Add(self.Bctrl, 1, wx.EXPAND, 1)
        label=wx.StaticText(panel, -1, "   " ) 
        hbox.Add(label,0, wx.ALL, 5)
        mainsizer.Add(hbox, 0, wx.EXPAND, 1)
        tip='Enter formula (e.g. "Li" )\nLeave empty to detect from file'
        self.Actrl.SetToolTip(wx.ToolTip(tip))
        tip='Enter formula (e.g. "Ag" )\nLeave empty to detect from file'
        self.Bctrl.SetToolTip(wx.ToolTip(tip))
        tip='Enter formula (e.g. "Mn8O16" )\nLeave empty to detect from file'
        self.Cctrl.SetToolTip(wx.ToolTip(tip))

        #
        # Force criterion
        #
        hbox = wx.BoxSizer(wx.HORIZONTAL)
        label=wx.StaticText(panel, -1, "   " ) 
        hbox.Add(label,0, wx.ALL, 5)
        label=wx.StaticText(panel, -1, " max ev/A in x " ) 
        hbox.Add(label, 1, wx.EXPAND, 1)
        label=wx.StaticText(panel, -1, " max eV/A in y " ) 
        hbox.Add(label, 1, wx.EXPAND, 1)
        label=wx.StaticText(panel, -1, " max eV/A in z " ) 
        hbox.Add(label, 1, wx.EXPAND, 1)
        mainsizer.Add(hbox, 0, wx.EXPAND, 1)
        #
        # force threshold
        #
        hbox = wx.BoxSizer(wx.HORIZONTAL)
        label=wx.StaticText(panel, -1, "   " ) 
        hbox.Add(label,0, wx.ALL, 5)
        self.FxCtrl = wx.TextCtrl(panel, -1, size=(60,b.GetDefaultSize()[1]),)
        hbox.Add(self.FxCtrl, 1, wx.EXPAND, 1)
        self.FyCtrl = wx.TextCtrl(panel, -1, size=b.GetDefaultSize())
        hbox.Add(self.FyCtrl, 1, wx.EXPAND, 1)
        self.FzCtrl = wx.TextCtrl(panel, -1, size=b.GetDefaultSize())
        hbox.Add(self.FzCtrl, 1, wx.EXPAND, 1)
        label=wx.StaticText(panel, -1, "    " ) 
        hbox.Add(label,0, wx.ALL, 5)
        mainsizer.Add(hbox, 0, wx.EXPAND, 1)
        tip='Allowed maximum residual force on ions in eV/A (e.g. 0.05).'+\
             '\nStructures with higher residual forces are neglected.'+\
             '\nThis entry is ingored when importing datafiles.'
        self.FxCtrl.SetToolTip(wx.ToolTip(tip))
        self.FyCtrl.SetToolTip(wx.ToolTip(tip))
        self.FzCtrl.SetToolTip(wx.ToolTip(tip))
        self.FxCtrl.SetValue('0.05')
        self.FyCtrl.SetValue('0.05')
        self.FzCtrl.SetValue('0.05')
        #
	# Disorder options:
        #
        hbox = wx.BoxSizer(wx.HORIZONTAL)
        label=wx.StaticText(panel, -1, "   " ) 
        hbox.Add(label,0, wx.ALL, 5)
        label=wx.StaticText(panel, -1, "Temperature ",size=(120,b.GetDefaultSize()[1]) ) 
        hbox.Add(label,0, wx.EXPAND, 1)
        self.TCtrl = wx.TextCtrl(panel, -1, size=(10,b.GetDefaultSize()[1]),)
        hbox.Add(self.TCtrl, 1, wx.EXPAND, 1)
        label=wx.StaticText(panel, -1, "    " ) 
        hbox.Add(label,0, wx.ALL, 5)
        mainsizer.Add(hbox, 0, wx.EXPAND, 1)
        tip='Temperature for averaing configurations with same stoichiometry,'+\
             '\ne.g. 300.0 [K] for room temperature'
        self.TCtrl.SetToolTip(wx.ToolTip(tip))
        self.TCtrl.SetValue('0.0000001')

        hbox = wx.BoxSizer(wx.HORIZONTAL)
        self.IsRandom = wx.CheckBox(panel, -1,"Random Diagram")
        hbox.Add(self.IsRandom, 1, wx.ALL, 5)
        hbox.Add(wx.StaticText(panel, -1,""), 0, wx.ALL, 5)
        self.IsRandom.SetValue(False)
        self.IsRandom.Bind(wx.EVT_CHECKBOX,self.OnChecked)
        mainsizer.Add(hbox, 0, wx.EXPAND, 1)
        tip='Select this to produce a random diagram'
        self.IsRandom.SetToolTip(wx.ToolTip(tip))
        #
        # OK and cancel buttons
        #
        buttonsizer = wx.BoxSizer(wx.HORIZONTAL)
        okbutton = wx.Button(panel, wx.ID_OK, "OK")
        cancelbutton = wx.Button(panel, wx.ID_CANCEL, "Cancel")
        buttonsizer.Add(okbutton, 1, wx.ALIGN_CENTER | wx.ALL, 5)
        buttonsizer.Add(cancelbutton, 1, wx.ALIGN_CENTER | wx.ALL, 5)
        okbutton.SetDefault()
        okbutton.Bind(wx.EVT_BUTTON, self.OnOK)
        cancelbutton.Bind(wx.EVT_BUTTON, self.OnCancel)
        mainsizer.Add(buttonsizer,0, wx.ALIGN_CENTER | wx.ALL, 1 )

        panel.SetSizer(mainsizer)
        panel.Layout()
   
    def OnChecked(self, event):
       if self.IsRandom.GetValue():
          self.Actrl.Enable(False)           
          self.Bctrl.Enable(False)           
          self.Cctrl.Enable(False)           
          self.Actrl.SetValue('B')
          self.Bctrl.SetValue('C')
          self.Cctrl.SetValue('A')
          self.filectrl.Enable(False)           
          self.dirctrl.Enable(False)           
          self.FxCtrl.Enable(False)           
          self.FyCtrl.Enable(False)           
          self.FzCtrl.Enable(False)           
       else:
          self.Actrl.Enable(True)           
          self.Bctrl.Enable(True)           
          self.Cctrl.Enable(True)           
          self.filectrl.Enable(True)           
          self.dirctrl.Enable(True)           
          self.FxCtrl.Enable(True)           
          self.FyCtrl.Enable(True)           
          self.FzCtrl.Enable(True)           

    #----------------------------------------------------------------------
    def OnOpenDir(self, event):
        """
        Open DirDialog and import POSCAR files 
        """
        from os import getcwd
        dlg = wx.DirDialog(self, "Choose directory with structures:",
                           style=wx.DD_DEFAULT_STYLE,
                           defaultPath = getcwd(), 
                           )
        if dlg.ShowModal() == wx.ID_OK:
           #
           # set directory control entry 
           #
           self.dirctrl.SetValue(dlg.GetPath())
           #
           # clean file control entry 
           #
           self.filectrl.SetValue('')
           #
           # clean endpoint text entries
           #
           self.Actrl.SetValue('')
           self.Bctrl.SetValue('')
           self.Cctrl.SetValue('')
        dlg.Destroy()
 
    #----------------------------------------------------------------------
    def OnOpenFile(self, event):
        """
        Create and show the Open FileDialog
        """
        wildcard = "TERA diagram (*.tera,*.ter)|*.tera;*.ter|" \
            "Data file (*.dat,*.txt)|*.dat;*.txt"

        dlg = wx.FileDialog(
            self, message="Choose a file",
            defaultDir=self.currentDirectory, 
            defaultFile="",
            wildcard=wildcard,
            style=wx.OPEN | wx.MULTIPLE | wx.CHANGE_DIR
            )
        if dlg.ShowModal() == wx.ID_OK:
            paths = dlg.GetPaths()

            s = ''
            for path in paths:
               s += path+' '
            self.filectrl.SetValue(s)
            self.dirctrl.SetValue('')

            #
            # in case one file is chosen
            #
            if len( paths ) == 1 : 
               f = paths[0]
               #
               # try reading it and edit endpoints 
               #
               ext = get_extension(f)
               if ext == 'tera' or ext == 'ter' : 
                  sx = ReadTeraFile(f)
                  A=self.Actrl.GetValue()
                  B=self.Bctrl.GetValue()
                  C=self.Cctrl.GetValue()
                  if len(A) == 0 : 
                     self.Actrl.SetValue(sx.pdb['A']['name'])
                  if len(B) == 0 : 
                     self.Bctrl.SetValue(sx.pdb['B']['name'])
                  if len(C) == 0 : 
                     self.Cctrl.SetValue(sx.pdb['C']['name'])
        dlg.Destroy()

    #----------------------------------------------------------------------
    def OnOK(self, event):
        from poscar import ImportStructures, UpdateEndpoints, \
                           UpdateFormationEnergies, GetConvexHull, \
                           GetDelaunay
        from debugger import QhullDebugger
        from os.path import exists
        try:
           #
           # in case directory was chosen, import all data 
           #
           d = self.dirctrl.GetValue()
           #                    
           # import structues with energies                    
           #                    
           T=self.GetTemperature()
           if T < 0 : 
              text='Cannot Read Temperature (not a positive number)'
              ShowStatusErrorMessage(self.statusbar,text)
              return
           fthresh=self.GetForceThreshold()
           if len(fthresh) == 1 : 
               if fthresh[0] == -1 :
                  text='Cannot Read force threshold for x direction'
                  ShowStatusErrorMessage(self.statusbar,text)
                  return
               elif fthresh[0] == -2 : 
                   text='Cannot Read force threshold for y direction'
                   ShowStatusErrorMessage(self.statusbar,text)
                   return
               elif fthresh[0] == -3 : 
                   text='Cannot Read force threshold for z direction'
                   ShowStatusErrorMessage(self.statusbar,text)
                   return

           if len( d ) > 0 : 
              if not exists(d):
                 text='Directory "{:}" not found'.format(d)
                 ShowStatusErrorMessage(self.statusbar,text)
                 return
              if len(self.Actrl.GetValue()) == 0 : 
                 text='Endpoint B missing'
                 ShowStatusErrorMessage(self.statusbar,text)
              elif len(self.Bctrl.GetValue()) == 0 : 
                 text='Endpoint C missing'
                 ShowStatusErrorMessage(self.statusbar,text)
              elif len(self.Cctrl.GetValue()) == 0 : 
                 text='Endpoint A missing'
                 ShowStatusErrorMessage(self.statusbar,text)
              else:
                 ax = ImportStructures( directory=d ,log =  self.parent.C.PageLog,
                         force=fthresh, Temperature=T ) 
                 #                 
                 # get endpoints                 
                 #                 
                 species,nion=self.GetEndpoints()
                 #                 
                 # update endpoint in database 
                 #                 
                 sx = UpdateEndpoints( species, nion, log=self.parent.C.PageLog, Temperature=T) 
                 #                 
                 # calculate formation energies 
                 #                 
                 failed=UpdateFormationEnergies(sx, ax,log=self.parent.C.PageLog, Temperature=T) 
                 #                 
                 # continue only if none has failed 
                 #                 
                 if failed is None: 
                    #                 
                    # compute convex hull 
                    #                 
                    GetConvexHull(sx)
                    #                 
                    # get lowest energetic structures
                    #                 
                    x = UpdateEndpoints( species, nion, DoLog=False , Temperature=T) 
                    #                 
                    # calculate lowest formation energies 
                    #                 
                    UpdateFormationEnergies(x, sx, lowest=True,\
                        verbosity=0, log=self.parent.C.PageLog, Temperature=T) 
                    #                 
                    # compute convex hull 
                    #                 
                    GetConvexHull(x,ldepth=True)
                    #                 
                    # compute delaunay triangulation 
                    #                 
                    GetDelaunay(x)
                    #                 
                    # store reference of x here                  
                    #              
                    self.Projects.append( x )  
                    #
                    # destroy if everything went well 
                    # 
                    self.Destroy()
                 else:
                    text='Endpoint '+failed+' not found'
                    ShowStatusErrorMessage(self.statusbar,text)
           #
           # in case file was chosen, load file 
           #
           F = self.filectrl.GetValue().split()
           if self.IsRandom.GetValue(): 
              F = [ 'random.dat' ] 
           if len(F)>0:  
              for f in F: 
                 #
                 # Get file extension
                 #                 
                 ext = get_extension(f)
                 #                 
                 # read file depending on format
                 #                 
                 #                 
                 # if file is in tera format, read it                 
                 #                 
                 lcontinue=False
                 if ext == 'tera' or ext == 'ter': 
                    ax = ReadTeraFile(f)
                    if len(self.Actrl.GetValue()) == 0 : 
                       self.Actrl.SetValue(ax.pdb['A']['name'])
                    if len(self.Bctrl.GetValue()) == 0 : 
                       self.Bctrl.SetValue(ax.pdb['B']['name'])
                    if len(self.Cctrl.GetValue()) == 0 : 
                       self.Cctrl.SetValue(ax.pdb['C']['name'])
                    lcontinue = True
                 if ext == 'dat' or ext =='txt': 
                    if len(self.Actrl.GetValue()) == 0 : 
                       self.Actrl.SetValue('B')
                    if len(self.Bctrl.GetValue()) == 0 : 
                       self.Bctrl.SetValue('C')
                    if len(self.Cctrl.GetValue()) == 0 : 
                       self.Cctrl.SetValue('A')
                    species,nion=self.GetEndpoints()
                    ax = ReadDataTable(f,species,nion, \
                       self.statusbar,self.IsRandom.GetValue(),log=self.parent.C.PageLog)
                    lcontinue = True
                    
                 if lcontinue:
                    #                 
                    # get endpoints                 
                    #                 
                    if ext == 'tera' or ext =='ter':
                       species,nion=self.GetEndpoints()
                       #                 
                       # update endpoint in database 
                       #                 
                    sx = UpdateEndpoints( species, nion, log=self.parent.C.PageLog, Temperature=T ) 
                    #                 
                    # calculate formation energies 
                    #                 
                    failed=UpdateFormationEnergies(sx, ax,log=self.parent.C.PageLog, lowest=True, Temperature=T) 
                    #                 
                    # continue only if none has failed 
                    #                 
                    if failed is None: 
                       #                 
                       # get lowest energetic structures
                       #                 
                       x = UpdateEndpoints( species, nion, DoLog=False, Temperature=T ) 
                       #                 
                       # compute convex hull 
                       #                 
                       GetConvexHull(sx)
                       #                 
                       # calculate lowest formation energies 
                       #                 
                       UpdateFormationEnergies(x, sx, lowest=True,\
                           log=self.parent.C.PageLog, Temperature=T) 
                       #                 
                       # compute convex hull 
                       #                 
                       GetConvexHull(x,ldepth=True)
                       #                 
                       # compute delaunay triangulation 
                       #                 
                       GetDelaunay(x)
                       #                 
                       # store reference of x here                  
                       #              
                       self.Projects.append( x )  
                       #
                       # destroy if everything went well 
                       # 
                       self.Destroy()
                    else:
                       self.logger.Add2Log( '\n ERROR: \"'+f+\
'\" COULD NOT BE READ!\n')
        finally:
            self.GetParent().Enable(True)

    #----------------------------------------------------------------------
    def OnCancel(self, event):
        self.Destroy()

    #----------------------------------------------------------------------
    def GetTemperature(self):
       " gets Temperature from entry box " 
       tmin=0.0
       try: 
          tmin=float( self.TCtrl.GetValue() ) 
          return tmin
       except ValueError: 
          tmin =-1000.0

    #----------------------------------------------------------------------
    def GetForceThreshold(self):
       " gets force threshold from entry boxes " 
          
       f=[]
       try: 
          fmin=float( self.FxCtrl.GetValue() ) 
          f.append(fmin)
       except ValueError: 
          f=[]
          f.append(-1)
          return f

       try: 
          fmin=float( self.FyCtrl.GetValue() ) 
          f.append(fmin)
       except ValueError: 
          f=[]
          f.append(-2)
          return f

       try: 
          fmin=float( self.FzCtrl.GetValue() ) 
          f.append(fmin)
       except ValueError: 
          f=[]
          f.append(-3)
          return f
       
       return f 
     
    #----------------------------------------------------------------------
    # get species and nions from control entries A,B and C
    def GetEndpoints(self):
       import re
       from numpy import argsort, asarray
       species=[]
       nion=[]

       a=self.Actrl.GetValue()
       b=self.Bctrl.GetValue()
       c=self.Cctrl.GetValue()
       A=re.findall(r'([A-Z][a-z]*)(\d*)', a)
       B=re.findall(r'([A-Z][a-z]*)(\d*)', b)
       C=re.findall(r'([A-Z][a-z]*)(\d*)', c)
        
       if len(A)==0 and len(B)==0 and len(C)==0: 
          return species, nion

       #phase C is largest
       Atmp=[]
       Btmp=[]
       Ctmp=[]
       Atmp2=[]
       Btmp2=[]
       Ctmp2=[]

       maxiter=max( len(C), len(A), len(B) )

       for i in range(0,maxiter):
           if i > len(A)-1: 
              Atmp.append( '' )
              Atmp2.append( 0 )
           else:
              Atmp.append( A[i][0] )
              if len(A[i][1]) == 0:
                 Atmp2.append( 1 )
              else:
                 Atmp2.append( int(A[i][1]) )

           if i > len(B)-1: 
              Btmp.append( '' )
              Btmp2.append( 0 )
           else:
              Btmp.append( B[i][0] )
              if len(B[i][1]) == 0:
                 Btmp2.append( 1 )
              else:
                 Btmp2.append( int(B[i][1]) )

           if i > len(C)-1: 
              Ctmp.append( '' )
              Ctmp2.append( 0 )
           else:
              Ctmp.append( C[i][0] )
              if len(C[i][1]) == 0:
                 Ctmp2.append( 1 )
              else:
                 Ctmp2.append( int(C[i][1]) )

       species=[ Atmp,Btmp,Ctmp ] 
       nion=[ Atmp2,Btmp2,Ctmp2 ] 

       return species, nion

    #----------------------------------------------------------------------
    # in case files are passed as arguments via terminal, open them 
    def OpenArguments(self):
       from os.path import exists,isdir,isfile,join,abspath
       from os import getcwd
       from sys import argv
       if len(argv[1:]) > 0 : 
          sf=''
          sd=''
          for a in argv[1:]: 
             #
             # check if file exists, if not auto-complete path
             #
             if exists(a):
                if ( isfile( a ) ):
                   sf+='{:} '.format(abspath(a))
                elif ( isdir( a ) ):
                   sd+='{:}'.format(abspath(a))
                   break
             else:
                cd = getcwd()
                afull=join(cd,a)
                if exists( afull ):
                   sf+='{:} '.format(afull)
                else:
                   if ( isdir( a ) ) : 
                      sd+='{:} '.format(a)
                   else:
                      print 'file or directory {:} not found'.format(a)
          if len( sf ) >0:
             self.filectrl.SetValue(sf)

          if len( sd ) >0: 
             self.dirctrl.SetValue(sd)
# DEBUG:
#          #self.dirctrl.SetValue("/home/mkaltak/github/TERA/wxpython/structures")
#          #self.filectrl.SetValue("/home/mkaltak/github/TERA/wxpython/diagrams/File1.tera /home/mkaltak/github/TERA/wxpython/diagrams/2D.tera ")
          self.Bctrl.SetValue("Ag")
          self.Cctrl.SetValue("Mn8O16")
          self.Actrl.SetValue("Li")
#          #self.filectrl.SetValue("/home/mkaltak/github/TERA/wxpython/data/energies.dat")
#          #self.filectrl.SetValue("/home/mkaltak/github/TERA/wxpython/data/random.dat")
#          self.filectrl.SetValue("/home/mkaltak/github/TERA/wxpython/diagrams/Newest.ter")
# DEBUG:



#------------------------------------------------------------------------------
def ShowStatusErrorMessage(status, text):
   """ 
     Shows error in status bar 
   """ 
   timer=status.GetTopLevelParent().timer
   status.SetBackgroundColour('RED')
   status.SetStatusText(text)
   status.Refresh()
   timer.Start(50)

#------------------------------------------------------------------------------
def ShowStatusWarningMessage(status, text, delay=0):
   from time import sleep
   """ 
     Shows Warning in status bar 
   """ 
   timer=status.GetTopLevelParent().timer
   status.SetBackgroundColour('YELLOW')
   status.SetStatusText(text)
   status.Refresh()
   timer.Start(10)
   if delay>0:
      sleep(delay)
   
#------------------------------------------------------------------------------
def ReadWriteCfgFile(write=False,p=None):
   """
     Writes a default config file to user's home directory
   """
   from os.path import exists,isdir,join
   from os import mkdir 
   from wx import StandardPaths
   import pickle
   from poscar import Data

   LOCAL_FN = "tera.cfg"
   ULDD = StandardPaths.Get().GetUserLocalDataDir()
   if not isdir(ULDD):
      mkdir(ULDD)
   path = join(ULDD, LOCAL_FN)

   if exists(path) and p is not None:
      #
      # in case file exists and settings given, write them to the file
      #
      if write: 
         with open(path, 'wb') as output:
             pickle.dump(p.CfgSaveLoad,output)
             #
             # save selected attributes 
             #
             for  a in p.CfgSaveLoad : 
                pickle.dump(getattr(p,a), output)
      else:
         with open(path, 'r') as input:
             p.CfgSaveLoad = pickle.load(input)
             #
             # load selected attributes 
             #
             for  a in p.CfgSaveLoad : 
                setattr(p, a, pickle.load(input) )

   elif not exists(path):
      sx=Data(cfg=True)
      with open(path, 'wb') as output:
          pickle.dump(sx.CfgSaveLoad,output)
          #
          # save selected attributes 
          #
          for  a in sx.CfgSaveLoad : 
             pickle.dump(getattr(sx,a), output)

