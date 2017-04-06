import sys
import numpy as np 
from scipy.spatial import Delaunay
# a small number defining zero
small=1.E-6

# very smal
very_small=1.E-12

# define some important constants
#number of distinct configurations for same stoichiometry
CMAX=100 


#path of structures 
DirDB='./structures'

#this defines the maximum size of the supercell
global SCMAX
SCMAX=100

#which end point is used as reference point 
global PHASES
PHASES=[ 'A', 'B', 'C' ] 

#barycentric coordinates
global PHASE_COORD
PHASE_COORD=np.zeros((len(PHASES),len(PHASES)-1))
PHASE_COORD[0,0] = 1.0000000000000000
PHASE_COORD[0,1] = 0.0000000000000000
PHASE_COORD[1,0] = 0.5000000000000000
PHASE_COORD[1,1] = 0.8660254037844386
PHASE_COORD[2,0] = 0.0000000000000000
PHASE_COORD[2,1] = 0.0000000000000000


#multiplicative converstion factor from eV to kJ/mol 
kJmol=96.48531
kBolt=0.000086173303  # [ eV / K ]

#******************************************************************************
# convex hull  class
#******************************************************************************
class hull_class: 
   def __init__(self):
      #number of configurations in db
      self.npoints=0
      #initialize database 
      self.points=0
      #vertices of boundary on lines AB, AC and BC
      self.AB_vertices=0
      self.AC_vertices=0
      self.BC_vertices=0
      #initialize inner simplices
      self.inner_simplices=0
      #store simplices also as paths
      self.path=[]
      #all vertices
      self.vertices=0
      #distance of every point to the hull
      self.distance=0
      #keys of points 
      self.keys=[]
      #Delaunay triangulation simplices
      self.tri=Delaunay
      self.tri_points=0
      self.ReactionEnergies = 0 
   
#******************************************************************************
# database class
#******************************************************************************
class Data: 
   def __init__(self,cfg=False):
      from os.path import exists
      from os import makedirs, getcwd
      from color import maps
      #
      # output directory
      #
      self.OutDir = getcwd()
      #
      # initialize database 
      #
      self.db={}
      #
      # Full datababse, contains points with high energy as well 
      #
      self.dbF={}
      #
      #database of phases
      self.pathVESTA=''
      self.MinBondLength=2.6
      self.pdb={}
      #hull 
      self.hull=hull_class()
      # stable phases
      self.stable_keys=[]
      #
      # Points and tie line settings 
      # 
      #
      self.fontsize=24               # fontsize 
      self.nsample= 15               # number of patches in one carthesian direction
      self.ps=75                     # point size of data points
      self.pt_color = '#000000'      # point color
      self.lw=2                      # line width of tie lines 
      self.tl_color = '#000000'      # line color 
      self.alpha = 0.25              # transperancy of 3D polygons
      #
      # Grid specific 
      #
      self.xi=[0,0.25, 0.5,  1., 2, 4, 8]     # grid default subdivision for x concent
      self.yi=[0,0.25, 0.5,  1., 2, 4, 8]     # grid default subdivision for y concent
      self.zi=[0,0.25, 0.5,  1., 2, 4, 8]     # grid default subdivision for z concent
      self.Resolution = None                  # resolution of grid subdivision
      self.nlines = None                      # number of subdivision
      self.GridColors = ('#009900','#cc6699','#668cff') # grid colors 
      self.BorderColor = '#000000'            # Border colors 
      self.zorder = 3                         # zorder of grid
      self.ls = None                          # line style 
      self.lwGrid = 3                         # line style 
      self.lwBorder = 3                         # line style 
      self.UseUniform = True                  # use uniform gird 
      #
      # normalization endpoint and name
      #
      self.NormIDX=2                 # normalization index 
      #
      # what is shown by default
      #
      #
      # empty page
      #
      self.EmptyCollection = None
      self.EmptyCollectionOn = False
      #
      # scatter object for data poins 
      #
      self.PointsScatter = None   
      self.PointsScatter3D = None   
      self.PointsScatterOn = True          # show points 
      #
      # Tie lines collection 
      #
      self.TieLinesCollection = None
      self.TieLinesCollection3D = None
      self.TieLinesCollectionOn = True             # show tie lines 
      #
      # Grid lines collection
      #
      self.GridLinesCollection = None
      self.GridLinesCollection3D = None
      self.GridLinesCollectionOn = True            # show grid 
      #
      # Binding energy patch collection
      #
      self.BindingContourSet = None
      self.BindingContourSet3D = None
      self.BindingCollection3D = None
      self.BindingContourSetOn = False         # show formation energy
      #
      # Stability patch collection
      #
      self.StabilityContourSet = None
      self.StabilityContourSet3D = None
      self.StabilityCollection3D = None
      self.StabilityContourSetOn = False      # show stability
      #
      # Reaction patch collection
      #
      self.ReactionContourSet = None
      self.ReactionContourSet3D = None
      self.ReactionCollection3D = None
      self.ReactionContourSetOn = False      # show stability

      self.BorderCollectionOn = True                   # show border
      self.BorderCollection = []             # show tie lines 
      self.BorderCollection3D = []             # show tie lines 

      self.LineToolOn = False        # linetool on?
      self.LineToolPointCollection = []             # show tie lines 
      self.LineToolLinesCollection = []             # show tie lines 
      #
      # Label
      #
      self.LabelsOn = True                               # show labels
      self.LabelAHandle = None                           # Handle bottom left 
      self.LabelBHandle = None                           # Handle bottom right
      self.LabelCHandle = None                           # Handle top 
      self.LabelA = None                                 # label bottom left 
      self.LabelB = None                                 # label bottom right
      self.LabelC = None                                 # label top 
      self.LabelAOffset = (.15,.075)                       # offset for A
      self.LabelBOffset = (-.025,.075)                    # offset for B
      self.LabelCOffset = (.025,-.025)                   # offset for C
      self.LabelColor = { "A":'#009900', \
                          "B" :'#cc6699',\
                          "C" :'#668cff' }               # colors for A, B, C

      self.SnapOn = False           # snap cursor
      self.BoundaryOn = False       # show Boundary
      self.Toggle3DOn = False      # 3D Mode switch
      self.PickedCoordinates = [ ] 
      # 
      # data title 
      # 
      self.title = 'New diagram'
      #
      # colorbar specific 
      #
      self.ColorBar = None 
      self.CBAxes = None
      self.CBAxes2 = None
      self.CmdIDX = 70
      self.Cmd = getattr(maps(),maps().AllMapNames[self.CmdIDX])
      self.CmdIntersect = 50
      self.CBLabel = 'Formation energy (eV)'
      #
      # face color
      #
      self.facecolor = 'white'
      #
      # what to save and load
      #
      self.SaveLoad = ['OutDir',\
                       'db',\
                       'dbF',\
                       'pdb',\
                       'pathVESTA',\
                       'MinBondLength',\
                       'hull',\
                       'stable_keys',\
                       'fontsize',\
                       'nsample',\
                       'ps',\
                       'pt_color',\
                       'lw',\
                       'tl_color',\
                       'alpha',\
                       'xi',\
                       'yi',\
                       'zi',\
                       'Resolution',\
                       'nlines',\
                       'GridColors',\
                       'BorderColor',\
                       'zorder',\
                       'lwGrid',\
                       'lwBorder',\
                       'ls',\
                       'UseUniform',\
                       'NormIDX',\
                       'EmptyCollectionOn',\
                       'PointsScatterOn',\
                       'TieLinesCollectionOn',\
                       'GridLinesCollectionOn',\
                       'BindingContourSetOn',\
                       'StabilityContourSetOn',\
                       'ReactionContourSetOn',\
                       'LabelsOn',\
                       'LabelA',\
                       'LabelB',\
                       'LabelC',\
                       'LabelAOffset',\
                       'LabelBOffset',\
                       'LabelCOffset',\
                       'LabelColor',\
                       'SnapOn',\
                       'BoundaryOn',\
                       'Toggle3DOn',\
                       'title',\
                       'Cmd',\
                       'CmdIDX',\
                       'CmdIntersect',\
                       'CBLabel',\
                       'facecolor' ]
      #
      # what to save and load
      #
      self.CfgSaveLoad = ['OutDir',\
                       'pathVESTA',\
                       'MinBondLength',\
                       'fontsize',\
                       'nsample',\
                       'ps',\
                       'pt_color',\
                       'lw',\
                       'tl_color',\
                       'alpha',\
                       'xi',\
                       'yi',\
                       'zi',\
                       'Resolution',\
                       'nlines',\
                       'GridColors',\
                       'BorderColor',\
                       'zorder',\
                       'lwGrid',\
                       'lwBorder',\
                       'ls',\
                       'UseUniform',\
                       'NormIDX',\
                       'EmptyCollectionOn',\
                       'PointsScatterOn',\
                       'TieLinesCollectionOn',\
                       'GridLinesCollectionOn',\
                       'BindingContourSetOn',\
                       'StabilityContourSetOn',\
                       'ReactionContourSetOn',\
                       'LabelsOn',\
                       'LabelA',\
                       'LabelB',\
                       'LabelC',\
                       'LabelAOffset',\
                       'LabelBOffset',\
                       'LabelCOffset',\
                       'LabelColor',\
                       'SnapOn',\
                       'BoundaryOn',\
                       'Toggle3DOn',\
                       'title',\
                       'Cmd',\
                       'CmdIDX',\
                       'CmdIntersect',\
                       'CBLabel',\
                       'facecolor' ]

      self.LoadFromFile=[ 'fontsize',\
                'nsample',\
                'ps',\
                'pt_color',\
                'lw',\
                'tl_color',\
                'alpha',\
                'xi',\
                'yi',\
                'zi',\
                'Resolution',\
                'nlines',\
                'GridColors',\
                'BorderColor',\
                'zorder',\
                'lwGrid',\
                'lwBorder',\
                'ls',\
                'UseUniform',\
                'NormIDX',\
                'EmptyCollectionOn',\
                'PointsScatterOn',\
                'TieLinesCollectionOn',\
                'GridLinesCollectionOn',\
                'BindingContourSetOn',\
                'StabilityContourSetOn',\
                'ReactionContourSetOn',\
                'LabelsOn',\
                'LabelA',\
                'LabelB',\
                'LabelC',\
                'LabelAOffset',\
                'LabelBOffset',\
                'LabelCOffset',\
                'LabelColor',\
                'SnapOn',\
                'BoundaryOn',\
                'Toggle3DOn',\
                'title',\
                'Cmd',\
                'CmdIDX',\
                'CmdIntersect',\
                'CBLabel',\
                'facecolor' ]
      #
      # set specific defaults, stored in config file 
      #
      if  not cfg : 
         self.SetDefaultsFromFile()
      
   def SetDefaultsFromFile(self):
      from dialog import ReadWriteCfgFile
      from color import maps
      ReadWriteCfgFile(p=self)
      self.Cmd = getattr(maps(),maps().AllMapNames[self.CmdIDX])
       
   #
   # copies some specific settings from x to it self
   #
   def CopyPlottingSettings(self, x):
      for entry in self.LoadFromFile: 
	  if hasattr( self, entry ) and hasattr( x, entry ):  
             setattr( self, entry, getattr( x,entry) )

#******************************************************************************
#Class "structure" containes information of one configuration
#******************************************************************************
class structure:
   def __init__(self):
      #small letters for same variables as in DB
      self.idkey=''
      self.title=''
      self.nions=0
      self.nion=[]
      self.species=[]
      self.nspecies=0
      self.coord=0
      self.bravais=0
      #reciprocal lattice 
      recipro=0

#******************************************************************************
#Class for attaching a np array to another class as object 
#******************************************************************************
class AttachArray(np.ndarray):
   def __new__(cls, a):
      obj = np.asarray(a).view(cls)
      return obj

   def Clean(self):
      self[:,:] = 0
      return self

#-------------------------------------------------------------------------------
# adds an entry to the database 
#
# Every entry in DB is characterized by 
# IDKEY ...... "NIONS-NSPECIES:SPECIES[0].NSPECIES[0]. ... "
# NIONS ...... total number of ions 
# NION ....... array of ions per species 
# NSPECIES ... array of number of ions per species 
# SPECIES .... array of species name 
# BRAVAIS .... primitive lattice vectors
# COORD ...... coordinates of ions 
# UNIQUE ..... True if its a unique entry in the database
#
#-------------------------------------------------------------------------------
def AddToDB(x,IDKEY,NIONS,NION,NSPECIES,SPECIES,BRAVAIS,COORD,EN,UNIQUE,\
   FORCE,DEGEN, PATH, fname= None, logger=None):
   #add dictionary with IDKEY-Key to database
   x.db[IDKEY]={}
   #name of configuration
   x.db[IDKEY]["name"]=IDKEY
   #add number of ions 
   x.db[IDKEY]["nions"]=NIONS
   #add number of ions per species
   x.db[IDKEY]["nion"]=NION
   #add number of species
   x.db[IDKEY]["nspecies"]=NSPECIES
   #add names of species
   x.db[IDKEY]["species"]=SPECIES
   #add Bravais lattice 
   x.db[IDKEY]["lattice"]=BRAVAIS
   #add coordinates of ions
   x.db[IDKEY]["coordinates"]=COORD
   #add energy 
   x.db[IDKEY]["energy"]=EN
   #add maximum force
   x.db[IDKEY]["maxForce"]=FORCE
   #add information of uniquess of entry 
   x.db[IDKEY]["IsUnique"]=UNIQUE

   #degeneracy 
   x.db[IDKEY]["degeneracy"]=DEGEN
   # file path
   x.db[IDKEY]["path"]=PATH

   x.dbF[IDKEY] = x.db[IDKEY]

   if logger is not None  and fname is not None: 
      text = '\"'+fname+'\" imported as: '+IDKEY
      logger.Add2Log(text)

#-------------------------------------------------------------------------------
# adds an entry to the subdatabase 
#-------------------------------------------------------------------------------
def AddToSubDB(s,DB,key,stoichiometry ):
   if key not in DB:
      sys.exit('Error in AddToSubDB: database DB does not contain entry with '+\
      key)

   #add dictionary with key to database
   s.db[key]={}
   #name of configuration
   s.db[key]["name"]=DB[key]["name"]
   #add number of ions 
   s.db[key]["nions"]=DB[key]["nions"]
   #add number of ions per species
   s.db[key]["nion"]=DB[key]["nion"]
   #add number of species
   s.db[key]["nspecies"]=DB[key]["nspecies"]
   #add names of species
   s.db[key]["species"]=DB[key]["species"]
 
   if 'lattice' in DB[key]: 
      #add Bravais lattice 
      s.db[key]["lattice"]=DB[key]["lattice"]
   if 'coordinates' in DB[key]:
      #add coordinates of ions
      s.db[key]["coordinates"]=DB[key]["coordinates"]
   #add energy 
   s.db[key]["energy"]=DB[key]["energy"]
   #add max force
   s.db[key]["maxForce"]=DB[key]["maxForce"]
   #add information of uniquess of entry 
   s.db[key]["IsUnique"]=DB[key]["IsUnique"]
   #add stoichiometry 
   s.db[key]["stoichiometry"]=stoichiometry
   # list of equivalent keys, i.e. same barycentric coordinates
   s.db[key]["SameBC"]=[]
   s.db[key]["SameBC_E"]=[]
   #degeneracy 
   s.db[key]["degeneracy"]=DB[key]["degeneracy"]
   #file path
   s.db[key]["path"]=DB[key]["path"]
 
#-------------------------------------------------------------------------------
# creates a POSCAR of the key IDKEY from the database 
#-------------------------------------------------------------------------------
def MakePOSCAR(IDKEY,x,fname):
   from numpy import savetxt
   import os 
   prec='%20.14f'
   #remove file if it exists
   try:
      os.remove(fname)
   except OSError:
      pass

   f = open(fname , 'a')
   #write string to file if its open
   f.write( x.dbF[IDKEY]["name"]+'\n' )
   f.write( '   1  \n' )
   savetxt(f,x.dbF[IDKEY]["lattice"],fmt=prec)
   for i in range(0,len(x.dbF[IDKEY]["species"])):
      f.write(' '+x.dbF[IDKEY]["species"][i]+' ')
   f.write('\n')
   for i in range(0,len(x.dbF[IDKEY]["species"])):
      f.write(' '+str(x.dbF[IDKEY]["nion"][i])+' ')
   f.write('\n')
   f.write('Direct \n')
   savetxt(f,x.dbF[IDKEY]["coordinates"],fmt=prec)

#-------------------------------------------------------------------------------
# creates a vesta file of the key IDKEY from the database 
#-------------------------------------------------------------------------------
def MakeVESTA(IDKEY,x,fname):
   from numpy import savetxt
   from numpy import dot, sqrt, arccos, degrees
   import os 
   prec='%20.14f'
   #remove file if it exists
   try:
      os.remove(fname)
   except OSError:
      pass

   f = open(fname , 'a')
   
   text = '#VESTA_FORMAT_VERSION 3.3.0\n#created by TERA\n\nCRYSTAL\n\nTITLE\n{:}\n\n'.\
           format(x.dbF[IDKEY]["name"])
   f.write( text )
   # calculate cell parameters
   avec=x.dbF[IDKEY]["lattice"][0]
   bvec=x.dbF[IDKEY]["lattice"][1]
   cvec=x.dbF[IDKEY]["lattice"][2]
   a=sqrt(dot(avec,avec))
   b=sqrt(dot(bvec,bvec))
   c=sqrt(dot(cvec,cvec))

   s=dot(avec,bvec)/(a*b)
   alpha=degrees( arccos(s) ) 

   s=dot(bvec,cvec)/(b*c)
   beta=degrees( arccos(s) ) 

   s=dot(cvec,avec)/(c*a)
   gamma=degrees( arccos(s) ) 

   text = 'CELLP\n  {:.6f}   {:.6f}   {:.6f}  {:.6f}  {:.6f}  {:.6f}\n  0.000000   0.000000   0.000000   0.000000   0.000000   0.000000\nSTRUC\n'.\
           format(a,b,c,beta,gamma,alpha)
   f.write( text )

   k=-1
   ks=-1
   for itype in x.dbF[IDKEY]["species"]:
      ks+=1
      for iion in range(0,x.dbF[IDKEY]["nion"][ks]):
            k+=1
            c=x.dbF[IDKEY]["coordinates"][k]
            sr='{:}{:}'.format(itype,iion+1)
            text = '{:3d} {:3s} {:10s} 1.0000   {:.6f}   {:.6f}   {:.6f}    1a   1\n                            0.000000   0.000000   0.000000  0.00\n'.\
            format(k+1,itype,sr,c[0],c[1],c[2])
            f.write( text )
   text='  0 0 0 0 0 0 0\nTHERI 0 \n'
   f.write( text )
   k=-1
   ks=-1
   for itype in x.dbF[IDKEY]["species"]:
      ks+=1
      for iion in range(0,x.dbF[IDKEY]["nion"][ks]):
            k+=1
            c=x.dbF[IDKEY]["coordinates"][k]
            text = '{:}  {:}{:}  1.0000\n'.format(k+1,itype,iion+1)
            f.write( text )
   WriteBondsToVesta(x,IDKEY,f)

   #savetxt(f,x.dbF[IDKEY]["lattice"],fmt=prec)
   #for i in range(0,len(x.dbF[IDKEY]["species"])):
   #   f.write(' '+x.dbF[IDKEY]["species"][i]+' ')
   #f.write('\n')
   #for i in range(0,len(x.dbF[IDKEY]["species"])):
   #   f.write(' '+str(x.dbF[IDKEY]["nion"][i])+' ')
   #f.write('\n')
   #f.write('Direct \n')
   #savetxt(f,x.dbF[IDKEY]["coordinates"],fmt=prec)

def WriteBondsToVesta(x,IDKEY,f,bondlength=None):
   TMLIST=['Sc','Ti','V','Cr','Mn','Fe','Co','Ni','Cu','Zn']
   TMLIST=TMLIST+\
   ['Y','Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd']
   
   isTMO=False
  
   if bondlength is None:
      bondlength = 2.6

   TMs=[]
   if 'O' in x.dbF[IDKEY]["species"]:
      for TM in TMLIST: 
         if TM in x.dbF[IDKEY]["species"]:
            TMs.append(TM)
      if len(TMs)>0 :
         text='   0 0 0\nSBOND\n'
         f.write( text )
         k=0
         for TM in TMs: 
            k+=1
            text='{:3d} {:5s}{:6s}    0.00000{:11.6f}  0  1  1  0  1  0.250  1.000 180 180 180\n'.\
               format(k,TM,'O',bondlength)
            f.write( text )
         text='  0 0 0 0\nSTYLE\nDISPF 147487\nMODEL   2  1  0\nSURFS   0  1  1\nSECTS  96  1\nFORMS   0  1\nATOMS   0  0  1\nBONDS   1\nPOLYS   1\nVECTS 1.000000'
         f.write(text)

#-------------------------------------------------------------------------------
# Determine endpoints of ternary diagram
# nion = 3xnions matrix of stoichiometries for every subsystem
# species is a 3xnions  array of strings, where
# species[0,:]=["Asub_1", "Asub_2", ..., "Asub_nions"] is the first sub system
# species[1,:]=["Bsub_1", "Bsub_2", ..., "Bsub_nions"] is the second sub system
# species[2,:]=["Csub_1", "Csub_2", ..., "Csub_nions"] is the third sub system
# 
# phase A = (Asub_1)_nion[0,0] (Asub_2)_nion[0,1] ... (Asub_nions)_nion[0,nions]
# phase B = (Bsub_1)_nion[1,0] (Bsub_2)_nion[1,1] ... (Bsub_nions)_nion[1,nions]
# phase C = (Csub_1)_nion[2,0] (Csub_2)_nion[2,1] ... (Csub_nions)_nion[2,nions]
#-------------------------------------------------------------------------------
def UpdateEndpoints(species,nion,log=None,DoLog=True,Temperature=0):
   from basis import GetConcentrations
   from numpy import append,zeros,array,argsort
   def argsort(seq):
       return sorted(range(len(seq)), key=seq.__getitem__)
   #initialize database 
   dat=Data()
  
   #initialize barycentric lattice vectors 
   barylatt=zeros( (2,2) ) 
   barylatt[0,0]=PHASE_COORD[0,0]-PHASE_COORD[2,0]
   barylatt[1,0]=PHASE_COORD[0,1]-PHASE_COORD[2,1]
   barylatt[0,1]=PHASE_COORD[1,0]-PHASE_COORD[2,0]
   barylatt[1,1]=PHASE_COORD[1,1]-PHASE_COORD[2,1]
   
   dat.pdb["BaryLatt"]=barylatt
   dat.pdb["muvector"]=[]

   for i in PHASES: 
      dat.pdb[ i ] = {} 
   
   #create the sub systems
   smax=len(nion[0][:])
   #distinct atoms list
   atoms=[]
   #define phase A
   i=0
   for phases in PHASES:
      si=[] ; ni=[]
      name=''
      if len(species[i][:]) != len(nion[i][:]) :
         sys.exit('DefineTernaryEndpoints reports inconsistent array sisze')
      #store subsystem definitions
      for j in range(0,smax): 
         if nion[i][j]>0:
            si=append(si,species[i][j])
            ni=append(ni,nion[i][j])
      #sort the arrays alphabetically in case user puts them in weired order 
      idx1=argsort(si)
      #build the name 
      n=[]
      s=[]
      for j in range(0,len(si)):
         s.append(si[idx1][j])
         n.append(int(ni[idx1][j]))  
         if ni[idx1][j]==1:
            name+=si[idx1][j]
         else:
            name+=si[idx1][j]+str(n[j] )

      dat.pdb[ phases ]["species"]=s
      dat.pdb[ phases ]["nion"]=n
      dat.pdb[ phases ]["name"]=name
      if DoLog: 
         if log is None : 
            print '     Phase '+phases+':', name
         else:
            log.Add2Log('     Phase '+phases+': '+name)

      for j in s:
         atoms.append( str(j) ) 
      i+=1

   #determine all distinct atoms
   tmp= list( set().union( atoms ))
   #sort them
   tmp.sort()
   dat.pdb["EndpointAtoms"]=tmp
 
   #determine phase coefficient matrix 
   A=zeros( ( len(tmp),len(PHASES) ))
   l=0
   for i in PHASES:
      k=[0]*len(tmp)
      for species_i in dat.pdb[ i ]["species"]:
         idx=dat.pdb["EndpointAtoms"].index(species_i)
         jdx=dat.pdb[ i ]["species"].index(species_i)
         A[idx,l]=dat.pdb[ i ]["nion"][jdx]
         k[idx]=dat.pdb[ i ]["nion"][jdx]
      dat.pdb[ i ]["nion_sorted"]=array(k)
      #create latexname
      x=GetConcentrations( [ PHASE_COORD[l,0] , PHASE_COORD[l,1] ] )
      dat.pdb[ i ]["LaTeXName"]=CreateName(dat.pdb,x)
      l=l+1

   dat.pdb["phase_matrix"]=A

   return dat
#-------------------------------------------------------------------------------
# Maps fractional set of stoichiometries s=(x',y',z') to si=(x0,y0,z0)
#-------------------------------------------------------------------------------
def IntegerStoichiometry(s,sx,All=None,debug=None):
   from fractions import Fraction,gcd
  
   na=sum(sx.pdb["A"]["nion"])
   nb=sum(sx.pdb["B"]["nion"])
   nc=sum(sx.pdb["C"]["nion"])
   maxconc=max([na*SCMAX, nb*SCMAX, nc*SCMAX])

   #get lowest common multiple
   S = GetLCMVector( s ) 
   
   if debug : 
      print '**Stoichiometry and its integer **',s, sint 

   sint = [ S ] 
   #number of returned stoichiometries
   if All:
      for i in range(2,maxconc): 
         if  S[0]*i <=maxconc and S[1]*i <=maxconc and S[2]*i <=maxconc:
            sint.append([ S[0]*i, i*S[1], i*S[2] ] )

   if debug : 
       print sint 
        
   return sint

#-------------------------------------------------------------------------------
# calculates the greatest common divisor of a list of integers
# and returns its smallest value + multiples 
#-------------------------------------------------------------------------------
def GetGCDVector(N): 
   from fractions import gcd
   from numpy import abs 
  
   for i in range(0,len(N)):
      if abs( N[i] ) > small :
         g = N[i]
         break


   for i in range(0,len(N)):
      if abs( N[i] ) > small :
         m = N[i] 
         #check w.r.t. to others
         for j in range( i , len(N) ) :
            if abs( N[j] ) > small :
               l = N[j]
               g = gcd(g,l)

   #get prime factors of g 
   primes = list(GetPrimeFactors( g ) )


   ND=[]
   ND.append( N )
   i_=1
   for i in primes:
      X=[]
      for j in range(0,len(N)):
         X.append(N[j]/(i*i_))
      i_=i
      ND.append( X )

   return ND  


#-------------------------------------------------------------------------------
# calculates the lowest common multiple of a float vector s = [s0,s1,s2]
#-------------------------------------------------------------------------------
def GetLCM(s): 
   from fractions import Fraction
   from numpy import abs 

   def LCM(a,b): 
      from fractions import gcd
      return  a*b/gcd(a,b)

   f1=Fraction.from_float(s[0]).limit_denominator(1000)
   f2=Fraction.from_float(s[1]).limit_denominator(1000)
   f3=Fraction.from_float(s[2]).limit_denominator(1000)
   p1=f1.numerator
   q1=f1.denominator
   p2=f2.numerator
   q2=f2.denominator
   p3=f3.numerator
   q3=f3.denominator
  
   lcm = 1 

   if abs( s[0] ) > small : 
      lcm  = LCM( q1, lcm )

   if abs( s[1] ) > small : 
      lcm  = LCM( q2, lcm )

   if abs( s[2] ) > small : 
      lcm  = LCM( q3, lcm )

   return lcm 

#-------------------------------------------------------------------------------
# calculates the lowest common multiple of a float vector s = [s0,s1,s2]
# returns the complete vector mulitplied by lowest common mulitple 
#-------------------------------------------------------------------------------
def GetLCMVector(s,prec=None, GetLCM = None): 
   from fractions import Fraction
   from numpy import abs 

   def LCM(a,b): 
      from fractions import gcd
      return  a*b/gcd(a,b)

   if prec is None:
      f1=Fraction.from_float(s[0]).limit_denominator(1000)
      f2=Fraction.from_float(s[1]).limit_denominator(1000)
      f3=Fraction.from_float(s[2]).limit_denominator(1000)
   else: 
      f1=Fraction.from_float(s[0]).limit_denominator(prec)
      f2=Fraction.from_float(s[1]).limit_denominator(prec)
      f3=Fraction.from_float(s[2]).limit_denominator(prec)

   p1=f1.numerator
   q1=f1.denominator
   p2=f2.numerator
   q2=f2.denominator
   p3=f3.numerator
   q3=f3.denominator
  
   lcm = 1 

   if abs( s[0] ) > small : 
      lcm  = LCM( q1, lcm )

   if abs( s[1] ) > small : 
      lcm  = LCM( q2, lcm )

   if abs( s[2] ) > small : 
      lcm  = LCM( q3, lcm )


   if GetLCM :
      return [ (p1*lcm)/q1 , (p2*lcm)/q2, (p3*lcm)/q3  ] , lcm 
   else:
      return [ (p1*lcm)/q1 , (p2*lcm)/q2, (p3*lcm)/q3  ]


#-------------------------------------------------------------------------------
# Looks if a specific stoichiometry s=[x,y,z] is in the database
#-------------------------------------------------------------------------------
def IsInDatabase(s,sx,debug=None):
   import basis as b

   if debug: 
      print '=================================================================='
      print sx.pdb['A']["species"], sx.pdb['B']["species"],sx.pdb['C']["species"]
      print sx.pdb['A']["nion"], sx.pdb['B']["nion"],sx.pdb['C']["nion"]
      print s
 
    
   #map fractional numbers to integers 
   si=IntegerStoichiometry(s,sx, All=True)

   if debug: 
      print ' Integer stoichiometry:',si[0:min(2,len(si))]
      l = min(2,len(si))
   else: 
      l = len(si)     
 
   if len(si)==0: 
     si=IntegerStoichiometry(s,sx, All=True) 
     sys.exit('Error in IsInDatabase, no integer stoichiometry found')

   keys=[]
   for i in range(0,l): 
      IdFromStoichiometry(si[i],sx,keys)
 
   if debug : 
      print 'these keys found',keys[0:min(3,len(keys))]
      print '=================================================================='
     
   if len(keys)==0:
      return keys
   else:
      return keys
   
#-------------------------------------------------------------------------------
# builds a proper key for searching the x.db based on the stoichiometry and
# s.pdb
#-------------------------------------------------------------------------------
def IDKEYFromStoichiometryAndDB_PHASES(s, pdb):

   NionA=pdb["A"]["nion"]
   NionB=pdb["B"]["nion"]
   NionC=pdb["C"]["nion"]
   sA=pdb["A"]["species"]
   sB=pdb["B"]["species"]
   sC=pdb["C"]["species"]
  
   #multiply number of ions of every subspecies
   nionA=[]
   nionB=[]
   nionC=[]
   for i in range(0,len(NionA)):
      nionA.append(NionA[i]*s[0])
   for i in range(0,len(NionB)):
      nionB.append(NionB[i]*s[1])
   for i in range(0,len(NionC)):
      nionC.append(NionC[i]*s[2])
 
   #in case two or more endpoints contain the same species, special care must be
   #taken 
   #construct array of distinct atoms:
   species=[]
   n=[]
   #the first one is easy 
   for i in range(0,len(sA)):
      if nionA[i]>0 : 
         species.append( sA[i] ) 
         n.append( nionA[i] ) 

   l=[-1]*len(sB)
   #ions from phase B need to be scaned and compared if they occur in phase A
   for i in range(0,len(species)): 
      for j in range(0,len(sB)): 
         #in case one cooincident species is found do not add the name but the 
         if species[i] == sB[j] : 
            #mark the position 
            l[j]=i
   #now add only new species 
   for j in range(0,len(sB)):
      if l[j]<0: 
         if nionB[j]>0 : 
            species.append( sB[j] ) 
            n.append( nionB[j] ) 
      else:
         #add number of ions for existing species 
         n[ l[j] ] = n[ l[j] ] + nionB[j] 

   l=[-1]*len(sC)
   #ions from phase C need to be scaned and compared if they occur in phase A+B
   for i in range(0,len(species)): 
      for j in range(0,len(sC)): 
         #in case one cooincident species is found do not add the name but the 
         if species[i] == sC[j] : 
            #mark the position 
            l[j]=i
   #now add only new species 
   for j in range(0,len(sC)):
      if l[j]<0: 
         if nionC[j]>0 : 
            species.append( sC[j] ) 
            n.append( nionC[j] ) 
      else:
         #add number of ions for existing species 
         n[ l[j] ] = n[ l[j] ] + nionC[j] 
 
   #sort arrays 
   idx=np.argsort(species)
   species.sort()

   Nion=[]
   for i in range(0,len(idx)):
      Nion.append( n[ idx[i] ] ) 

   nions=sum(Nion)
   nspecies=len(species)

   #construct key
   idkey=''
   #create idkey of configuration
   idkey=str(nions)+'-'+str(nspecies)+':'
   for i in range(0,nspecies):
      idkey+=str(species[i])+'.'+str(Nion[i])+'.'


   #add c-iadd if iadd>-1
   if iadd is not None:
      if iadd >-1:
         idkey+='c-'+str(iadd)+'.'

   return idkey 

#-------------------------------------------------------------------------------
# id from stoichiometry
#-------------------------------------------------------------------------------
def IdFromStoichiometry(s, x, keys,debug=None):

   if debug : 
      print ' looking for',s 

   NionA=x.pdb["A"]["nion"]
   NionB=x.pdb["B"]["nion"]
   NionC=x.pdb["C"]["nion"]
   sA=x.pdb["A"]["species"]
   sB=x.pdb["B"]["species"]
   sC=x.pdb["C"]["species"]
  
   #multiply number of ions of every subspecies
   nionA=[]
   nionB=[]
   nionC=[]
   for i in range(0,len(NionA)):
      nionA.append(NionA[i]*s[0])
   for i in range(0,len(NionB)):
      nionB.append(NionB[i]*s[1])
   for i in range(0,len(NionC)):
      nionC.append(NionC[i]*s[2])
 
   #in case two or more endpoints contain the same species, special care must be
   #taken 
   #construct array of distinct atoms:
   species=[]
   n=[]
   #the first one is easy 
   for i in range(0,len(sA)):
      if nionA[i]>0 : 
         species.append( sA[i] ) 
         n.append( nionA[i] ) 

   l=[-1]*len(sB)
   #ions from phase B need to be scaned and compared if they occur in phase A
   for i in range(0,len(species)): 
      for j in range(0,len(sB)): 
         #in case one cooincident species is found do not add the name but the 
         if species[i] == sB[j] : 
            #mark the position 
            l[j]=i
   #now add only new species 
   for j in range(0,len(sB)):
      if l[j]<0: 
         if nionB[j]>0 : 
            species.append( sB[j] ) 
            n.append( nionB[j] ) 
      else:
         #add number of ions for existing species 
         n[ l[j] ] = n[ l[j] ] + nionB[j] 

   l=[-1]*len(sC)
   #ions from phase C need to be scaned and compared if they occur in phase A+B
   for i in range(0,len(species)): 
      for j in range(0,len(sC)): 
         #in case one cooincident species is found do not add the name but the 
         if species[i] == sC[j] : 
            #mark the position 
            l[j]=i
   #now add only new species 
   for j in range(0,len(sC)):
      if l[j]<0: 
         if nionC[j]>0 : 
            species.append( sC[j] ) 
            n.append( nionC[j] ) 
      else:
         #add number of ions for existing species 
         n[ l[j] ] = n[ l[j] ] + nionC[j] 
 
   #sort arrays 
   idx=np.argsort(species)
   species.sort()

   Nion=[]
   for i in range(0,len(idx)):
      Nion.append( n[ idx[i] ] ) 

   nions=sum(Nion)
   nspecies=len(species)

   NionDivided=GetGCDVector(Nion) 

   if debug : 
      print 'Species:',species
      print 'nion-all',NionDivided

   #loop over all stoichiometries found 
   for k in range(0, len( NionDivided ) ): 
      nions=sum(NionDivided[k])
      #construct key
      key=''
      #create key of configuration
      key=str(nions)+'-'+str(nspecies)+':'
      for i in range(0,nspecies):
         key+=str(species[i])+'.'+str(NionDivided[k][i])+'.'

      if debug: 
         print 'constructed key:',key

      if key in x.db: 

         # do not add if the same key is already present in the list
         Ladd=True
         for old in keys: 
            if old == key : 
              Ladd = False
                  
         if Ladd: 
            keys.append(key)

            if debug: 
               print 'key added',keys
      
      #test other configurations 
      key_=key
      for j in range(0,CMAX):
         key=key_+'c-'+str(j)+'.'
         if key in x.db:
            # do not add if the same key is already present in the list
            Ladd=True
            for old in keys: 
               if old == key : 
                 Ladd = False
                     
            if Ladd: 
               keys.append(key) 
 
#-------------------------------------------------------------------------------
# checks if string is a number
#-------------------------------------------------------------------------------
def is_number(s):
   try:
       float(s)
       return True
   except ValueError:
       return False

#-------------------------------------------------------------------------------
# Checks if string describes direct coordinates in POSCAR 
#-------------------------------------------------------------------------------
def IsCarthesian(s): 
   C=-1;c=-1;K=-1;k=-1
   #look for the first three letters
   C=s.rfind('C',0 ,2)
   c=s.rfind('c',0 ,2)
   K=s.rfind('K',0 ,2)
   k=s.rfind('k',0 ,2)
   
   if C>=0 or c>=0 or K>=0 or k>=0 : 
      return True
   else:
      return False

#-------------------------------------------------------------------------------
# Checks if string describes selective dynamics in POSCAR
#-------------------------------------------------------------------------------
def IsSelective(s): 
   s1=-1; S1=-1
   s1=s.find('se')
   s1=s.find('SE')
   s1=s.find('sE')
   s1=s.find('Se')
   
   if s1>=0 or S1>=0: 
      return True
   else:
      return False

def ToReciprocal( fhandle, nions, bravais, logger=None, forcethresh=None ):
  from numpy import zeros, dot
  from numpy.linalg import inv
  binv=inv( bravais ) 
  c = zeros(3)
  coord=zeros( shape=( nions, 3 ) ) 
  for i in range( 0, nions ): 
      s=fhandle.readline().split()
      #check if there are 3 coordinates 
      if len(s)< 3:   
         sys.exit(fhandle.name+' contain less than 3 ionic coordinates for ion #'+str(\
           i))
      c[0] = eval( s[0] ) 
      c[1] = eval( s[1] ) 
      c[2] = eval( s[2] ) 
      
      cr=dot( c, binv ) 
      coord[i,0] = cr[0]
      coord[i,1] = cr[1]
      coord[i,2] = cr[2]

  return coord

#-------------------------------------------------------------------------------
# Imports POSCAR structure to database
#-------------------------------------------------------------------------------
def ImportPOSCAR(fname,x,Temp=None,logger=None,forcethresh=None,
        Temperature=0, debug=None ):
   from numpy import argsort, ndarray, append, log
   import os

   # constants
   kB=0.000086173
  
   if Temp is None:
      Temp = 300.
 
   f=open(fname)
   tmp=f.readline().splitlines()
   name=tmp[0]

   s=f.readline().splitlines()
   scale=eval(s[0])
   if scale<0: 
      sys.exit('Sorry: no POSCARs with volume definition can be understand yet')
 
   #read the bravais lattice   
   BRAVAIS = ndarray(shape=(3,3), dtype=float, order='F')
   BRAVAIS[:,:]=0
   for i in range(0,3):
      s = f.readline().split()
      b1=eval(s[0]);b2=eval(s[1]);b3=eval(s[2])
      BRAVAIS[i,0]=scale*b1
      BRAVAIS[i,1]=scale*b2
      BRAVAIS[i,2]=scale*b3

   #read SPECIES in POSCAR 
   s = f.readline().split()
   # in case ions cannot be identified, write waring
   if is_number(s[0]): 
      print 'WARNING: "'+str(fname)+'" contains no information about '+\
         'ion type, file not imported.'
      return 
   species=[]
   for i in range(0,len(s)):
      species=append(species,str(s[i]))

   idx=argsort(species)
   SPECIES = sorted(species)

   #
   # number of species
   #
   NSPECIES=len(SPECIES)
      
   #number of ions of each SPECIES 
   b = f.readline().split()
   nion = []
   for i in range(0,NSPECIES):
      if is_number(eval(b[i])):
         #sort the array as well
         nion=append(nion,eval(b[i]))
      else:
         sys.exit('Error: number of individual ions in "'+str(fname)+'"'+\
           ' are corrupt')

   NION=nion[idx]
   #total number of ions
   NIONS = int(sum(NION))
 
   NION_=[0]*len(SPECIES)
   #store integer only 
   for i in range(0,len(NION)):
     NION_[i]=int(NION[i])


   #build id of structure 
   IDKEY=str(NIONS)+'-'+str(NSPECIES)+':'
   #UNIQUE is True if it is a unique configuration for this stoichiometry
   UNIQUE=True
   for i in range(0,NSPECIES):
      IDKEY+=str(SPECIES[i])+'.'+str(NION_[i])+'.'


   #look if database contains an entry with the same name 
   #if so, create a new IDKEY 
   if IDKEY in x.db : 
      #and change uniqueness to false
      x.db[IDKEY]["IsUnique"]=False
      UNIQUE = False
      for i in range(1,CMAX):
         if IDKEY+'c-'+str(i)+'.' not in x.db : 
            IDKEY+='c-'+str(i)+'.'
            break
      if i>CMAX-1:
         sys.exit('Error: maximum number of distinct configurations for one'+\
            ' stoichiometry reached!')
   #read coordinates of ions, this needs to be done NION times 
   coord = ndarray(shape=(NIONS,3), dtype=float, order='F')
   COORD = ndarray(shape=(NIONS,3), dtype=float, order='F')
   coord[:,:]=0
   COORD[:,:]=0
   
   #check if coordinates are direct 
   s = f.readline().split()

   #maybe POSCAR contains selective dynamics
   if IsSelective(s[0]):
      # read the next line to determine type of coordinates and ignore selective
      # dynamics 
      s = f.readline().split()
   elif IsCarthesian(s[0]):
      coord = ToReciprocal( f, NIONS, BRAVAIS, logger=logger )  
   else:
      #store unsorted coordinates temporarily
      for i in range(0,NIONS):
         s = f.readline().split()
         #check if there are 3 coordinates 
         if len(s)<3:   
            sys.exit(fname+' contain less than 3 ionic coordinates for ion #'+str(\
              i))
         coord[i,0]=eval(s[0])
         coord[i,1]=eval(s[1])
         coord[i,2]=eval(s[2])

   #sort coordinates accourding index array idx and store them
   l=0
   for isp in range(0,NSPECIES):
      #determine position of storing point
      #this the number of ions shifted for this species 
      ni=int(nion[idx[isp]])
      #this is the index in unsorted list 
      ipos = int(sum(nion[0:idx[isp]]))
      for i in range(0, int(NION[isp]) ) : 
         j = ipos +i
         COORD[l,:]=coord[j,:]
         l=l+1
 
#   sys.exit('STOP')
   NION=[]
   for i in range(0,NSPECIES):
      NION.append( int(nion[idx[i]]) )
    
   # search for energy in file should be the last line
   s = f.readline().split()

   if len(s) > 0:
      # check if line has Energy tag
      isEnergy = False
      for i in s : 
         if i == 'energy' or i == 'Energy' : 
            isEnergy = True 

         if is_number( i ) and isEnergy : 
            EN=eval( i )
   else:
       sys.exit('Something went wrong when reading file '+fname)

   # check if force line is present, if so add it 
   FORCE=[]
   # search for energy in file should be the last line
   s = f.readline().split()

   if len(s) > 0:
      # check if line has Energy tag
      isForce = False
      l=-1
      for i in s : 
         l+=1
         if i == 'Force' or i == 'force' : 
            isForce = True 
            l_=l

         if is_number( i ) and isForce : 
            FORCE.append( eval(i) ) 
   
   if len( FORCE ) != 3  : 
       FORCE = [ 0.,0.,0.]

   else:
      # do not import if force criterion is not satisfied 
      if forcethresh is not None : 
         k=-1
         for fi in forcethresh : 
             k+=1
             if FORCE[k] > fi : 
                f.close()
                return None
 
   # check for degeneracy entry (for entropy effects)
   f.close()
   
   DEGEN=1
   if Temperature>small :
      pattern='DEGENERACY'
      with open(fname) as f:
         for line in f:
             if pattern in line:
                DEGEN=eval(line.split()[2])
                # correct energy 
                if debug:
                   print fname, DEGEN, EN, Temperature

   if Temperature>small :
      EN = AddVibrationalEnergy( fname, Temperature,EN, NION, SPECIES,
                                verbosity=0)

   # extract file path 
   PATH = os.path.basename(fname)

   #add this information to the database 
   AddToDB(x,IDKEY,NIONS,NION,NSPECIES,SPECIES,BRAVAIS,COORD,EN, 
      UNIQUE, FORCE, DEGEN, PATH,fname = fname, logger=logger)

   return IDKEY


#-------------------------------------------------------------------------------
# adds vibrational energy contribution 
#-------------------------------------------------------------------------------
def AddVibrationalEnergy( fname,T, EN, NION, SPECIES,verbosity=0):
    def Free( y) :
        from numpy import abs
        f0=0.

        # T=303.10K
        f0p5=-0.04335034
        f1p0=-0.11264543
        f2p0=-0.16645551
        # T=404.10K
#        f0p5=-0.07295442
#        f1p0=-0.17989218
#        f2p0=-0.28115652
        # T=302.10K from diff
        f0p5=-0.084139
        f1p0=-0.345939
        f2p0=-0.407234

        if abs(y) < small:
            return 0.0
        elif y>small and y <= 0.5 :
            ya =0.0
            yb =0.5
            fa = 0.0
            fb = f0p5
        elif y>0.5 and y <=1.0 :
            ya =0.5
            yb =1.0
            fa = f0p5
            fb = f1p0
        elif y>1.0 and y <=2.0 :
            ya =1.0
            yb =2.0
            fa = f1p0
            fb = f2p0
        return (fb-fa)/(yb-ya)*( y-ya )+fa

    if T < small:
        return EN

    F=0.
    if 'Ag' in SPECIES and 'Mn' in SPECIES:
        idxA=SPECIES.index('Ag')
        idxM=SPECIES.index('Mn')
        s=NION[ idxM ]/8.
        y=NION[ idxA ]/s
        F=Free(y)*s
        if verbosity >0:
           if 'Ag' in SPECIES and not ('Li' in SPECIES ):
              print T, NION, SPECIES, s,y, Free(y),F

    if 'Ag' in SPECIES and not ('Mn' in SPECIES or 'Li' in SPECIES):
        F=-0.059007
        if verbosity >0:
           print T, NION, SPECIES,F

    if T>small :
       EN = EN+F
    return EN


#-------------------------------------------------------------------------------
# Creates the database of structures 
#-------------------------------------------------------------------------------
def ImportStructures(directory=None,log=None,status=None,force=None, Temperature=0):
   import os.path
   from numpy import dot,sqrt

   #initialize full database
   x=Data()

   fout=0
   #make list of entries
   if directory is not None: 
      if os.path.exists(directory):
         l=0
         for file in os.listdir(directory):
             if file.endswith(".vasp"):
                l+=1

             if status is not None: 
                status.Initialize(l)

         for file in os.listdir(directory):
             if file.endswith(".vasp"):
               f = os.path.join(directory,file)
               key=ImportPOSCAR(f,x,logger=log,forcethresh=force, Temperature=Temperature)
               if key is None:
                    fout +=1
               if status is not None:
                  status.UpdateProgressbar(' Importing files...')
         if log is not None: 
            log.Add2Log(str(l-fout)+' structures imported ')
      else:
         sys.exit('No directory "'+str(directory)+'" found!')
   return x

#-------------------------------------------------------------------------------
# Compute formation energy and filters entries in x corresponding to endpoints 
# and stores the result in sx, also a list of equivalence classes is formed where
# each entry contains the keys with same barycentric coordinates but different 
# energy in ascending order 
#-------------------------------------------------------------------------------
def UpdateFormationEnergies(sx, x, lowest=None, Temperature=0, verbosity=0, log=None):
   import basis as b 
   from numpy import amax as np_amax
   from numpy import amin as np_amin
   from numpy import append as np_append
   from numpy import array, dot, argsort
   from numpy import concatenate
   from numpy import copy as np_copy
   from copy import copy 
   global SCMAX, CMAX

   # 
   # copy plotting settings from x to sx
   # 
   sx.CopyPlottingSettings( x ) 

   uc=[]
   #determine keys of all phases and their energies 
   for phase in PHASES: 

      species = sx.pdb[phase]["species"]
      nspecies = len(species)
      i=0
      lfound=False
      for key in x.db.keys():
         i+=1
         nspecies_x = len( x.db[key]["species"] )
         if  nspecies_x == nspecies : 
            #
            # now check if species are same
            #
            lfound=False
            k=0
            for j in range(0,len( species ) ):
               if species[j] == str(x.db[key]["species"][j]):
                  k+=1
            if k == len(species) : 
               lfound = True
            #
            # now check their abudance 
            #
            if lfound : 
               x0 = sx.pdb[phase]["nion"]
               x1 = x.db[key]["nion"] 
               alpha=GetScalingFactor(x0,x1)
               if alpha <= 0 : 
                  lfound = False
               elif alpha > 0 : 
                  break
      #      
      # return if no key found
      #
      if not lfound : 
         return phase
      #
      # otherwise continue
      #
      sx.pdb[phase]["key"]=key
      sx.pdb[phase]["energy"]=x.db[key]["energy"]/float(alpha)
      sx.pdb["muvector"] = np_append( sx.pdb["muvector"],sx.pdb[phase]["energy"])

      if log is not None: 
          if verbosity>0:
             txt='Endpoint {:} : {:} : {:.4f}'.format( phase,key,sx.pdb[phase]["energy"])
             log.Add2Log(txt)

   #get a list of compositions being in the sub-system spaned by sx.pdb
   for key in x.dbF.keys():
      lcheck,stoichiometry = PartOfSubDB( key, sx, x.dbF,uc  )
      if lcheck: 
         #determine largest supercell
         tmp=[SCMAX]
         tmp = np_append( tmp, stoichiometry) 
         SCMAX=int(max( SCMAX, np_amin( tmp ) ) )
         #add entry to subspace 
         AddToSubDB(sx, x.dbF, key, stoichiometry)
         #determine formation energy 
         DE=x.dbF[key]["energy"]
         l=0
         for phase in PHASES:
            DE = DE-stoichiometry[l]*sx.pdb[phase]["energy"]
            l+=1
         #scale formation energy 
         sx.db[key]["formationE_SC"]=DE
         DE=DE/sum(stoichiometry)
         sx.db[key]["formationE"]=DE
         #also calculate barycentric coordinates x'+y'+z'=1
         sx.db[key]["bary-coord"]=array(stoichiometry,dtype=float)/\
            sum(stoichiometry)
         #calculate carthesian coordinates 
         sx.db[key]["cart-coord"]=sx.pdb["BaryLatt"].dot(\
            sx.db[key]["bary-coord"][0:2])
         #store multiplication factor for every phase 
         sx.db[key]["stoichiometry"]=stoichiometry
         # 
         # keep a copy of entry for the full catalog
         # 
         sx.dbF[key] = sx.db[key]

         # in case one wants to keep only those entries with lowest formation
         # energy 
         if lowest:
            # loop over current keys in sub database and replace those entries             
            # with higher energies and create a list of keys with higher energy
            for key_sub in sx.db.keys():
               # compare current bary-centric coordinates 
               bc_diff= sx.db[key]["bary-coord"] - sx.db[key_sub]["bary-coord"]
               diff=dot(bc_diff,bc_diff)

               if diff < very_small and key_sub != key :

                  #delete the entry in sx with higher energy 
                  fe_tot = sx.db[key]["formationE"] 
                  fe_sub = sx.db[key_sub]["formationE"] 
                  # here key_sub in sx has to be removed
                  if fe_tot < fe_sub : 
                     #merge lists of stoichiometrically equivalent structres
		     sx.db[ key ][ "SameBC" ]=sx.db[ key ][ "SameBC" ]+sx.db[ key_sub ][ "SameBC" ]
                     dummy = concatenate((sx.db[ key ][ "SameBC_E"],sx.db[key_sub ][ "SameBC_E" ]),axis=0)
                     sx.db[ key ][ "SameBC_E" ]=np_copy(dummy)
                     del sx.db[key_sub]
		     sx.db[ key ][ "SameBC" ].append( key_sub ) 
		     sx.db[ key ][ "SameBC_E" ] = np_append( \
	                sx.db[ key ][ "SameBC_E" ], fe_sub ) 
                     break
                  # here key in sx has to be removed
                  else:
		     sx.db[ key_sub ][ "SameBC" ]=sx.db[ key ][ "SameBC" ]+sx.db[ key_sub ][ "SameBC" ]
                     del sx.db[key]
		     sx.db[ key_sub ][ "SameBC" ].append( key ) 
		     sx.db[ key_sub ][ "SameBC_E" ] = np_append( \
	                sx.db[ key_sub ][ "SameBC_E" ], fe_tot ) 
                     break

   # sort equivalent entries, such that energies are in increasing order 
   for key in sx.db.keys():
     idx = argsort( sx.db[ key ][ "SameBC_E" ] ) 
     l=copy(sx.db[ key ]["SameBC"]) 
     del sx.db[ key ]["SameBC"]
     del sx.db[ key ]["SameBC_E"]
     sx.db[ key ]["SameBC"] = [] 
     for i in idx : 
         sx.db[ key ][ "SameBC" ].append( l[ i ] ) 

     
   for i in PHASES:
      if not sx.pdb[i]['key'] in sx.db.keys(): 
         # check which key is in this database
         IDKEY = sx.pdb[i]['key']
         #print 'prolem',sx.pdb[i]['key'], sx.pdb[i]['key'] in sx.db.keys()
         for j in range( 1, CMAX) :
            newkey='{:}c-{:}.'.format( IDKEY , j )
            if newkey in sx.db : 
               IDKEY=newkey
               break
         # replace endpoint key by new one
         sx.pdb[i]['key']=IDKEY
         
   if lowest:   
      if log is not None: 
          k=0
          for key in sx.db.keys():
             k+=len( sx.db[key]["SameBC"] ) + 1  
          txt = '{:} stoichiometrically distinict structures found out of {:}'.\
              format(len(sx.db.keys()), k )
          log.Add2Log(txt)
 
          # this averages energies          
          AverageEnergies( sx, Temperature, verbosity=verbosity) 

          if verbosity>1:
             entry=[]
             for key in sx.db.keys():
                entry.append(key)
             entry.sort()
             for key in entry:
	       text='   {:}: {:.4f}'.format(sx.db[key]["name"],sx.db[key]["formationE"])
               log.Add2Log(text)
               if len( sx.db[ key ][ "SameBC" ] ) > 1 : 
	          for k in sx.db[ key ][ "SameBC" ] : 
	             text='     {:}: {:.4f}'.format(k,sx.dbF[k]["formationE"])
                     log.Add2Log(text)
      else:
         print len(sx.db.keys()), ' structures with lowest energies imported'
   else:
      if log is not None: 
         log.Add2Log('includes '+str(len(sx.db.keys()))+' structures in total')
      else:
         print len(sx.db.keys()), ' structures in total imported'

   return None

#-------------------------------------------------------------------------------
# averages energies via the partition sum < E >  = Sum_n=1^N E_n exp( -E_n/kb T)
# using degeneracy this is <E> = sum_m=1^M Omega_m E_m exp( -E_m/kb T) 
# 
# In general the idea is to write Z = Z' Z_0 where Z_0 = exp( -E_0/(kB T) )
# with E_0 being the lowest energy in the ensamble. 
# Z' refers to E'_n = E_n - E_0, so that 
# F = -kB T ln( Z ) = -kbT ln( Z') -kB T ln( Z_0 ) = 
#   = F' + E_0 
#-------------------------------------------------------------------------------
def AverageEnergies(sx, T, verbosity=0):
   from numpy import zeros, argsort, exp , dot, log

   beta = 1./( kBolt*T)

   verbosity=0


   for key in sx.db.keys():
#       if sx.db[ key ][ 'stoichiometry' ][ 0 ] == 0 : 
#           verbosity = -5
#       else:
#           verbosity =0
       k=0
       fenergiesT0=zeros( len( sx.db[key]["SameBC"] ) +1 )
       fenergiesT0[ k ] = sx.db[ key ][ "formationE_SC" ] 
       nconf = sx.db[key]["degeneracy"]
       Zp =  sx.db[key]["degeneracy"] # *1 
       s = sum( sx.db[ key ][ 'stoichiometry' ] )
       E0 = sx.db[ key ][ 'energy' ]

       for same_s in sx.db[key]["SameBC"]:
           k+=1
           alpha=1
           if not sum( sx.dbF[ same_s ][ 'stoichiometry' ]) == s :
              alpha=s/sum( sx.dbF[ same_s ][ 'stoichiometry' ])
           fenergiesT0[ k ] = alpha*sx.dbF[ same_s ][ "formationE_SC" ]
           nconf += sx.dbF[same_s]["degeneracy"]
           if T> small :
              Em = alpha*sx.dbF[ same_s ][ "formationE_SC" ]-\
                   sx.db[key]['formationE_SC']
              if verbosity > 0 :
                  print same_s, beta, Em, sx.dbF[ same_s ]["degeneracy"]
              Zp += sx.dbF[ same_s ]["degeneracy"] * exp( -beta* Em )

       # calculate Boltzmann Weights , first entry has lowest energy -> exp =1
       k=0
       BoltzmannWeight=zeros( len( sx.db[key]["SameBC"] ) +1 )
       BoltzmannWeight[ k ] =  sx.db[key]["degeneracy"]/Zp
       E = sx.db[ key ][ "formationE_SC" ]*BoltzmannWeight[ k ] 

       if verbosity > 0 :
           print '-------------------------------------------------------------'
           print '    Degen E[SC] Ef[SC]  E-E_0 [BC]  Boltz  key   Beta=',beta,'Zp=',Zp  
       if verbosity > 0 :
           txt='       {:2}  {:.4f}   {:.4f}  {:.4f} {:.4f} {:}'.\
              format( sx.db[key]["degeneracy"],\
              sx.db[key]['energy'],
              sx.db[key]["formationE_SC"],0. ,
              BoltzmannWeight[k ],
              sx.dbF[key]["name"])
           print txt 

       for same_s in sx.db[key]["SameBC"]:
           k+=1
           alpha=1
           if not sum( sx.dbF[ same_s ][ 'stoichiometry' ]) == s : 
              alpha=s/sum( sx.dbF[ same_s ][ 'stoichiometry' ])              
           Em = alpha*sx.dbF[ same_s ][ "formationE_SC" ]-\
                sx.db[key]['formationE_SC']
           BoltzmannWeight[ k ] = sx.dbF[same_s]["degeneracy"]*\
                                  exp( -Em*beta )/Zp 
           E += sx.dbF[ same_s ][ "formationE_SC" ]*BoltzmannWeight[ k ] 
           if verbosity ==-5 :
              txt='       {:2}  {:.4f}   {:.4f}  {:.4f} {:.4f} {:}'.\
                 format( sx.dbF[same_s]["degeneracy"],\
                 sx.dbF[same_s]['energy'],
                 sx.dbF[same_s]["formationE_SC"],Em,
                 BoltzmannWeight[ k],
                 sx.dbF[same_s]["name"])
              print txt

       if verbosity > 0 :
          print 'fenergies,weights:',fenergiesT0, BoltzmannWeight
          print 'Zp,<E>,E0,F= ',Zp, dot( fenergiesT0,BoltzmannWeight ),E0,-kBolt*T*log(Zp)+E0
       if verbosity ==-5 :
          print   sx.db[key]['stoichiometry'][1],-kBolt*T*log(Zp)+E0,-kBolt*T*log(Zp),key
   
       sx.db[ key ][ 'EnsembleFormationE' ] = dot( fenergiesT0,BoltzmannWeight )/s
   return


#-------------------------------------------------------------------------------
# Checks if a composition is in the sub-space defined by phase A, B and C
#-------------------------------------------------------------------------------
def PartOfSubDB(key,sx,DB,uc):
   from math import ceil
   from numpy import abs as np_abs
   from numpy import zeros as np_zeros
   #initialize return values
   lcheck=True         #is entry in database?
   scale=[0]*len(PHASES)    #multiplicity factor for every phase

   #first check if composition contains an ion not being in phase A, B or C
   #and set up system of linear equation, most of the time overdetermined but
   #should have a unique solution apart from trivial one x=(0,0,...0), 
   #where x = len(Phases), the least square solution takes this issue into account 
   b=np_zeros( len(sx.pdb["EndpointAtoms"]) ) 
   l=0
   for i in DB[key]["species"]: 
      try:
         idx=sx.pdb["EndpointAtoms"].index( i )
         b[idx]=DB[key]["nion"][l]
         l+=1
      except ValueError:
         lcheck=False
         break

   #solving the system is only necessary if test above passed
   xret=np_zeros(3)
   if lcheck: 
      x =  np.linalg.lstsq(sx.pdb["phase_matrix"],b)[0]

      if len(x) != len(scale):
         sys.exit('Fatal error: IsInDB_PHASES report inconsistent array sizes!')

      #get rid of small numbers
      for i in range(0,len(scale)):
         if x[i] < 0 and np_abs(x[i])>small : 
            x=np_zeros(3)
            break
         else:
            xret = x
      if sum(x) < small:
         lcheck=False
   

   return lcheck, xret


#===============================================================================
# Math helper routines 
#===============================================================================
#-------------------------------------------------------------------------------
# Calculates reciprocal lattice vectors from a given Bravais matrix 
# in units of 2*pi (same as VASP)
#-------------------------------------------------------------------------------
def ReciprocalLattice(bravais):
   result=np.ndarray(shape=(3,3),dtype=float, order='F')
   
   b=np.cross( bravais[1,:],bravais[2,:] )
   s=np.dot( bravais[0,:], b )
   result[0,0]=b[0]/s
   result[0,1]=b[1]/s
   result[0,2]=b[2]/s

   b=np.cross( bravais[2,:],bravais[0,:] )
   s=np.dot( bravais[1,:], b )
   result[1,0]=b[0]/s
   result[1,1]=b[1]/s
   result[1,2]=b[2]/s

   b=np.cross( bravais[0,:],bravais[1,:] )
   s=np.dot( bravais[2,:], b )
   result[2,0]=b[0]/s
   result[2,1]=b[1]/s
   result[2,2]=b[2]/s

   return result 

#-------------------------------------------------------------------------------
# calculates convex hull of subsystem spanned by endpoints 
#-------------------------------------------------------------------------------
def GetConvexHull(sx,ldepth=None):
   from scipy.spatial import ConvexHull
   #set number of configuration points 
   npoints=len(sx.db.keys())
   
   #consistency check
   if npoints <= 0: 
      sys.exit('GetConvexHull has no points to work with')

   sx.hull.npoints=npoints

   #coordinates in configuration space 
   coord = np.zeros((sx.hull.npoints,3))
   dist  = np.zeros(sx.hull.npoints)
   #attach object to class
   sx.hull.points=AttachArray(coord)

   #loop over all entries and store coordinates 
   i=0
   keys=[]
   for key in sx.db.keys():
      sx.hull.points[i,0]=sx.db[key]["cart-coord"][0]
      sx.hull.points[i,1]=sx.db[key]["cart-coord"][1]
      #do not take entries with positive formation energy into account
      #those configurations are unstable on the one hand and on the other 
      #will deteroriate the convex hull
      if 'EnsembleFormationE' in sx.db[ key ] : 
         if sx.db[key]["EnsembleFormationE"]>0:
            sx.hull.points[i,2]=-small
         else:
            sx.hull.points[i,2]=sx.db[key]["EnsembleFormationE"]

      else:
         if sx.db[key]["formationE"]>0:
            sx.hull.points[i,2]=-small
         else:
            sx.hull.points[i,2]=sx.db[key]["formationE"]
      #add index to db
      sx.db[key]["hull_idx"]=i
      keys=np.append(keys,key)
      i+=1

   #attach keys 
   sx.hull.keys=AttachArray(keys)

   #compute convex hull 
   hull = ConvexHull(sx.hull.points,qhull_options='Qv')

   #remove boundary simplices and add corresponding polygons
   #also store keys of stable phases 
   RegularizeBoundaries(sx,hull)

   #make up for wrong formation energies
   i=0
   for key in sx.db.keys():
      if 'EnsembleFormationE' in sx.db[ key ] : 
         sx.hull.points[i,2]=sx.db[key]["EnsembleFormationE"]
      else:
         sx.hull.points[i,2]=sx.db[key]["formationE"]
      i+=1

   #optionally, calculate depth of each point to the hull 
   if ldepth:
      CalculateHullDistance(sx)


#-------------------------------------------------------------------------------
# removes boundary simplices and regularizes 
#-------------------------------------------------------------------------------
def RegularizeBoundaries(sx,hull):
   from basis import IsNotBoundary
   from scipy.spatial import ConvexHull
   from matplotlib.path import Path

   #determine boundary simplices with 2D convex hull method. 
   CBadd=True
   ACadd=True
   ABadd=True
   for key in sx.db.keys(): 
      #for C-B line 
      if sx.db[key]["bary-coord"][0]<small:
         if CBadd:
            CBdat=[ np.array( [ 
                    sx.db[key]["cart-coord"][0], 
                    sx.db[key]["formationE"]
                   ] )] 
            CBadd=False
         else:
            CBdat=np.append( CBdat, \
                [ np.array([\
                            sx.db[key]["cart-coord"][0],
                            sx.db[key]["formationE"]
                           ])
                ], axis=0 )

      #for A-C line 
      if sx.db[key]["bary-coord"][1]<small:
         if ACadd:
            ACdat=[ np.array( [ 
                    sx.db[key]["cart-coord"][0], 
                    sx.db[key]["formationE"]
                   ] )] 
            ACadd=False
         else:
            ACdat=np.append( ACdat, \
                [ np.array([\
                            sx.db[key]["cart-coord"][0],
                            sx.db[key]["formationE"]
                           ])
                ], axis=0 )

      #for A-B line 
      if sx.db[key]["bary-coord"][2]<small:
         if ABadd:
            ABdat=[ np.array( [ 
                    sx.db[key]["cart-coord"][0], 
                    sx.db[key]["formationE"]
                   ] )] 
            ABadd=False
         else:
            ABdat=np.append( ABdat, \
                [ np.array([\
                            sx.db[key]["cart-coord"][0],
                            sx.db[key]["formationE"]
                           ])
                ], axis=0 )

   CBhull=ConvexHull( CBdat, qhull_options = 'Qv' ) 
   AChull=ConvexHull( ACdat, qhull_options = 'Qv' ) 
   ABhull=ConvexHull( ABdat, qhull_options = 'Qv' ) 
 
   #mark position of boundary vertices in the points list
   simplexCB=[]
   simplexAC=[]
   simplexAB=[]
   i=0
   for key in sx.db.keys():
      #for C-B line 
      if sx.db[key]["bary-coord"][0]<small:
         for j in CBhull.vertices:
            if np.abs(CBdat[j,0]-sx.hull.points[i,0])<small and\
               np.abs(CBdat[j,1]-sx.db[key]["formationE"])<small:
               #print 'CB',k1,i, CBdat[j,0], CBdat[j,1], sx.hull.points[i,:]
               simplexCB=np.append(simplexCB,i)

      #for A-C line 
      if sx.db[key]["bary-coord"][1]<small:
         for j in AChull.vertices:
            if np.abs(ACdat[j,0]-sx.hull.points[i,0])<small and\
               np.abs(ACdat[j,1]-sx.db[key]["formationE"])<small:
               #print 'AC',k2,i, ACdat[j,0], ACdat[j,1], sx.hull.points[i,:]
               simplexAC=np.append(simplexAC,i)

      #for A-B line 
      if sx.db[key]["bary-coord"][2]<small:
         for j in ABhull.vertices:
            if np.abs(ABdat[j,0]-sx.hull.points[i,0])<small*10 and\
               np.abs(ABdat[j,1]-sx.db[key]["formationE"])<small*10:
               #print 'AB',k3,i, ABdat[j,0], ABdat[j,1], sx.hull.points[i,:]
               simplexAB=np.append(simplexAB,i)
      i+=1

   #attach simplices 
   #calculate interior simplices only 
   int_simp=[]
   vert=[]
   for j in range(1, len(hull.simplices)) : 
      if IsNotBoundary(j,hull.simplices,sx.hull.points) == 1:
         if len(int_simp)>0:
            int_simp=np.append(int_simp, [hull.simplices[j]], axis = 0 )
            vert=np.append(vert, hull.simplices[j])
         else:
            int_simp=[ hull.simplices[j] ] 
            vert=hull.simplices[j]

   #add remaining vert from boundaries 
   vert=np.append(vert, simplexAB )
   vert=np.append(vert, simplexAC )
   vert=np.append(vert, simplexCB )
   vertices= np.unique(vert)

   j=0
   i=0
   stab=[]
   for key in sx.db.keys():
      if j < len(vertices):
         if i == vertices[j]: 
            stab.append(key)
            j+=1
      i+=1

   sCB = Sorter(sx, simplexCB )
   sAB = Sorter(sx, simplexAB )
   sAC = Sorter(sx, simplexAC )

   #store simplices 
   sx.hull.inner_simplices=AttachArray( int_simp ) 
   sx.hull.AB_vertices=AttachArray( sAB ) 
   sx.hull.AC_vertices=AttachArray( sAC ) 
   sx.hull.BC_vertices=AttachArray( sCB ) 


   sx.hull.vertices=AttachArray( vertices ) 
   sx.stable_keys = stab

   sx.hull.path = []
   for i in range(0,len(int_simp)) : 
      s = sx.hull.points[ int_simp[i],0:2 ]
      j = s.tolist()
      p = Path( j ) 
      sx.hull.path = np.append( sx.hull.path , p ) 


#-------------------------------------------------------------------------------
# helper 
#-------------------------------------------------------------------------------
def Sorter(sx,simplex):

   tmp = simplex 

   #sort along x-axis only
   coord=np.zeros(len(simplex))
   j=0
   for i in simplex: 
      i_=int(i)
      coord[j]=sx.hull.points[i_,0]
      j+=1

   idx = np.argsort(coord)
   return tmp[ idx[ : ] ] 


#-------------------------------------------------------------------------------
# calculates distance to hull of every point in sx for ternary diagrams
#-------------------------------------------------------------------------------
def CalculateHullDistance(sx,debug = None ):
   from basis import isPointInPath
   from numpy import append as np_append
   from numpy import abs as np_abs
   from numpy import minimum as np_minimum
   from basis import distance_to_hull
   from numpy import savetxt
   dist=[]
   l=0
   mb=0
   mv=0


   #for i in range(0,sx.hull.npoints):
   for key in sx.db.keys():
     
      i=sx.db[key]["hull_idx"]
      x=sx.hull.points[i,0]
      y=sx.hull.points[i,1]
      e=sx.hull.points[i,2]
      bc = sx.db[key]["bary-coord"] 

      if debug: 
         print '****VERBOSE****',key,'hull_idx',i
         
      #check if point is a vertex 
      m=0
      for k in sx.hull.vertices:
         if i==k: 
            m+=1
            l+=1
            mv+=1
            dist=np_append(dist, 0 ) 
            if debug :
               print 'VERTEX',bc, 0
            break
      
      if m == 0 :
         ms=0
         #loop over inner simplices 
         for j in range(0,len(sx.hull.inner_simplices)):
            simplex=sx.hull.inner_simplices[j]
            # get the simplex id of current point (X,Y)
            if isPointInPath(x,y,sx.hull.points[simplex]) : 
               z=distance_to_hull( x, y, simplex, sx.hull.points)
               dist=np_append(dist, np_abs(z-e) ) 
               l+=1
               if debug :
                  print 'INSIDE',bc, np_abs(z-e),key
               break
            #count down 
            else:
               #check if it is on the boundary 
               ms+=1

         mbi=0
         tie=0

         #in this case , point must be on the boundary line 
         if ms == len(sx.hull.inner_simplices):
            #this is the BC line 
            if bc[0] < small : 
               mb+=1
               mbi+=1
               #print key, 'should be on BC', bc, x,y,e,\
               #   sx.db[key]["energy"]
               d=10000.
               simplex=sx.hull.BC_vertices

               npts=len(simplex)
               for j in range(0,npts-1):
                  #rotate back on x-axis 
                  j_=int(simplex[j])
                  j1_=int(simplex[j+1])
                  x1=0.5000000000000000*sx.hull.points[j_,0]+\
                      0.8660254037844386*sx.hull.points[j_,1]
                  y1=sx.hull.points[j_,2]
                  x2=0.5000000000000000*sx.hull.points[j1_,0]+\
                      0.8660254037844386*sx.hull.points[j1_,1]
                  y2=sx.hull.points[j1_,2]
                  x0=0.5000000000000000*x+0.8660254037844386*y
                  k=(y1-y2)/(x1-x2)
                  D=(x1*y2-x2*y1)/(x1-x2)
                  z0=k*x0+D
                  d=np_minimum(d,np_abs(z0-e))
               dist=np_append(dist,d ) 
               if debug :
                  print 'BC',bc, d,key

            #this is the AC line 
            elif bc[1] < small : 
               mb+=1
               mbi+=1
               #print key, 'should be on AC', bc, e,\
               #   sx.db[key]["energy"]
               d=10000.
               simplex=sx.hull.AC_vertices
               npts=len(simplex)
               for j in range(0,npts-1):
                  j_=int(simplex[j])
                  j1_=int(simplex[j+1])
                  x1=sx.hull.points[j_,0]
                  y1=sx.hull.points[j_,2]
                  x2=sx.hull.points[j1_,0]
                  y2=sx.hull.points[j1_,2]
                  x0=x

                  k=(y1-y2)/(x1-x2)
                  D=(x1*y2-x2*y1)/(x1-x2)
                  z0=k*x0+D
                  d=np_minimum(d,np_abs(z0-e))
               dist=np_append(dist, d ) 
               if debug :
                  print 'AC',bc, d,key

            #this is the AB line 
            elif bc[2] < small : 
               mb+=1
               mbi+=1
               #print key, 'should be on AB', bc, e,\
               #   sx.db[key]["energy"]
               d=10000.
               simplex=sx.hull.AB_vertices
               npts=len(simplex)
               for j in range(0,npts-1):
                  j_=int(simplex[j])
                  j1_=int(simplex[j+1])
                  x1=1+0.5000000000000000*(sx.hull.points[j_,0]-1)-\
                      0.8660254037844386*sx.hull.points[j_,1]
                  y1=sx.hull.points[j_,2]
                  x2=1+0.5000000000000000*(sx.hull.points[j1_,0]-1)-\
                      0.8660254037844386*sx.hull.points[j1_,1]
                  y2=sx.hull.points[j1_,2]
                  x0=1+0.5000000000000000*(x-1)-0.8660254037844386*y
                  k=(y1-y2)/(x1-x2)
                  D=(x1*y2-x2*y1)/(x1-x2)
                  z0=k*x0+D
                  d=np_minimum(d,np_abs(z0-e))
               dist=np_append(dist,d ) 
               if debug :
                  print 'AB',bc, d,key,sx.db[key]["energy"]/sum(sx.db[key]["stoichiometry"]),\
                        sx.db[key]["energy"]
            else:
               #check tie lines of simplices 
               for j in range(1,len(sx.hull.inner_simplices)):
                  simplex=sx.hull.inner_simplices[j]
                                           
                  for k in [ 0, 1 ] : 
                     A=sx.hull.points[simplex[k],:]
                     B=sx.hull.points[simplex[k+1],:]
                     v1=B[0:2]-A[0:2]
                     v2=np.array( [x,y] ) - A[0:2]
                     norm1=np.sqrt(np.dot(v1,v1))
                     norm2=np.sqrt(np.dot(v2,v2))
                     prod=np.dot(v1,v2)/(norm1*norm2)
                     #in this case we have the correct line 
                     if np_abs( prod -1 ) <small: 
                        if np_abs(B[1]-A[1])>small:
                           k=(y-A[1])/(B[1]-A[1])
                           d=e-A[2]-(B[2]-A[2])*k
                           dist=np_append(dist, np_abs(d) ) 
                           tie+=1
                           break
                        elif np_abs(B[0]-A[0])>small:
                           k=(x-A[0])/(B[0]-A[0])
                           d=e-A[2]-(B[2]-A[2])*k
                           dist=np_append(dist, np_abs(d) ) 
                           tie+=1
                           break
                        else:
                           sys.exit('Error')  

                  # there are more than one simplex with same tie line,
                  # take care of this
                  if tie>0:
                     break

               if tie == 0 :
                  print bc, key
                  dist=np_append(dist, 1. ) 
                  #sys.exit('POINT not FOUND')
               else:
                  if debug :
                     print 'TIELINE',bc ,key
      #increase counter   
      i+=1
            
   if len(dist) != sx.hull.npoints : 
      print l , 'points in interior', mb ,'on the boundary',tie
      sys.exit('Error in GetConvexHull: '+str(len(dist))+' distances of '\
        +str(sx.hull.npoints)+' found!')

   sx.hull.distance=AttachArray(dist)

#-------------------------------------------------------------------------------
# Calculates the Delaunay triangulation of the data points
# To this end, the z-coordinate is the stability, i.e. the distance to the hull
# meaning that the Delaunay triangulation is a linear approximation to the 
# Free energy surface 
#-------------------------------------------------------------------------------
def GetDelaunay(sx):
   from scipy.spatial import Delaunay
   import numpy as np
   import matplotlib.pyplot as plt

   #store convex hull points to 2D array 
   x = sx.hull.points[:,0]
   y = sx.hull.points[:,1]
   points = np.zeros( (len(x),2) ) 

   points[:,0] = x[:]
   points[:,1] = y[:]

   #calculate trianglulation
   tri = Delaunay(points)
   #attach class to database
   sx.hull.tri = tri 

   #attach tri points as well 
   tri_points = np.zeros( (len(x),3) ) 
   tri_points[:,0] = x[:]
   tri_points[:,1] = y[:]
   tri_points[:,2] = sx.hull.distance[:]

   sx.hull.tri_points=AttachArray(tri_points)


def CalculateReaction_OLD(sx,X,ReturnWeights=None,LaTeX=True,NormIDX=None):
   from basis import GetConcentrations 
   from basis import isPointInPath
   from basis import GetBaryCentricCoordinates
   from basis import GetBinReactionWeights
   from numpy import append as np_append
   from numpy import array as np_array
   from numpy import dot as np_dot
   from numpy import sqrt as np_sqrt

   #triangle centre point
   C = [0.5000000000000000, 0.2886751345948129]

   #this is what is returned
   #intialize reaction list 
   r=[]
   #r is list of strings in the latex format specifying the reaction happening
   # the format here is [ name1,name2,...name_n, energy]
   # if the length of r is 
   #    2  ...  stable point or vertex chosen with energy =0              (i)
   #    3  ...  binary system of general point on tie line, energy = 0    (ii)
   #    4  ...  binary system of point point on tie line, energy != 0     (iii)
   #    5  ...  ternary system of general point, energy =0                (iv)
   #    5  ...  ternary system of point with energy !=0                   (v)

   #obtain stoichiometry of point X
   #x0=GetConcentrations([X[0],X[1]])
   x0=GetConcentrations2([X[0],X[1]],NormIDX=NormIDX, debug=False)
   name0=CreateName(sx.pdb, x0,LaTeX=LaTeX)
   #append name only if entry 
   r.append(name0)

   #collect keys if they can be found in database 
   keys=[]
   f = IsInDatabase(x0,sx,debug=False)
   if len(f)>1:
      #need to loop over double entries and use smallest energy 
      print f 
      sys.exit('Not implemented yet')
   elif len(f)==1:
      if f[0] in sx.db:  
         keys.append(f[0])

   #check if point is on a boundary 
   Xb=GetBaryCentricCoordinates( X ) 
   DeltaE=0

   #check if point is a vertex inside triangle 
   k=0
   #check if point is a vertex
   for v in sx.hull.vertices:
      p=sx.hull.points[int(v)]
      if np.abs(p[0]-X[0])<small and np.abs(p[1]-X[1])<small:
         k+=1
         break
   
   #only do something if it is not a vertex inside the triangle
   if k==0 : 
      #loop over all simplices
      ms=0
      for i in range(0,len(sx.hull.inner_simplices)):
         simplex=sx.hull.inner_simplices[i]
         points=sx.hull.points[simplex]
         # get the simplex id of current point (X,Y)
         if isPointInPath( X[0] , X[1] , points ): 
            ms+=1
            break
         else:
            #shift point towards centre of triangle 
            V=np_array( [ C[0]-X[0] , C[1]-X[1] ] )
            n=np_dot(V,V)
            if n>small: 
               v = [small*V[0]/np_sqrt( n ) ,small*V[1]/np_sqrt(n) ] 
               if isPointInPath( X[0]+v[0] , X[1]+v[1] , points ): 
                  ms+=1
                  break
            else:
               if isPointInPath( X[0]+small , X[1]+small , points ): 
                  ms+=1
                  break
               #else shake it even more
               elif isPointInPath( X[0]-small , X[1]+small , points ): 
                  ms+=1
                  break
               #else shake it even more
               elif isPointInPath( X[0]-small , X[1]-small , points ): 
                  ms+=1
                  break
               #else shake it even more
               elif isPointInPath( X[0]+small , X[1]-small , points ): 
                  ms+=1
                  break
               #else shake it even more
               elif isPointInPath( X[0] , X[1]-small , points ): 
                  ms+=1
                  break
               #else shake it even more
               elif isPointInPath( X[0] , X[1]+small , points ): 
                  ms+=1
                  break
               #else shake it even more
               elif isPointInPath( X[0]+small , X[1] , points ): 
                  ms+=1
                  break
               #else shake it even more
               elif isPointInPath( X[0]-small , X[1] , points ): 
                  ms+=1
                  break

      #Stop here if a propper simplex has not been found 
      if ms==0 : 
         print ' ERROR: no simplex found for ('+str(X[0])+','+str(X[1])+')'
         sys.exit('This is a fatal Error')

      #simplex has been found, calculate reaction weights and names of 
      #reactants, they must be in the database !
      else:
         #get weights 
         weights=GetReactionWeights2( X , points, NormIDX=NormIDX )

         keystmp=np.array(keys)
         for k in range(0,len(PHASES)):
            x=GetConcentrations2( [ points[k][0] , points[k][1] ],NormIDX=NormIDX )

            #append name for nonzero weights only
            if np.abs(weights[k])>small:
               name1=CreateName(sx.pdb,x,weights[k],LaTeX=LaTeX)
               r.append(name1)

            #look for a key
            f = IsInDatabase(x,sx,debug=False)
            
            if len(f)==1:
               if f[0] in sx.db:  
                  keystmp=np.append(keystmp,f[0] )
                  if np.abs(weights[k])>small:
                      keys.append(f[0])
            else:
               print x0 ,name0
               #sys.exit('Error: key not found!'+name0)

         #calculate reaction energy 
         DeltaE=ReactionEnergy3D(keystmp,sx,weights,Xb,X,\
            debug=False,NormIDX=NormIDX )

         #maybe weights need to be returned
         if ReturnWeights:
            Lamb=weights

   #append reaction energy      
   r.append(DeltaE)

   #maybe weights need to be returned
   if ReturnWeights:
      return keys,r,Lamb
   else:
      return keys,r
  

#-------------------------------------------------------------------------------
# Calculates reaction energy from keys 
#-------------------------------------------------------------------------------
def ReactionEnergy(keys,x0,xi,lam,sx,EApprox=None,debug=None):
   from numpy import sum as np_sum

   E0=0
   #there is a reaction only in the case of four keys
   if len(keys)==4 or len(keys)==3:
      #find supercell: 
      beta =  FindSuperCell( keys[0], sx, x0, debug = debug)
      
      #if no super cell was found, one has an arbitrary composition
      #in this case FindSuperCell returns -1 
      if beta > 0 : 
         #we set E0 to the scaled energy of the supercell 
         E0=sx.db[keys[0]]["energy"]/beta
 
         #loop over remaining three key entries 
         for i in range(1,len(keys)):
            alpha=lam[i-1]
            #keep non-zero weights only
            if np.abs(alpha)>small:
               #find super cell of configuration 
               beta = FindSuperCell( keys[i], sx, xi[i-1])
               #print an error if something went wrong here
               if beta <0 :
                  sys.exit(' Error in ReactionEnergy 1, beta is negative '+str(i))
               else:
                  E0 = E0 - alpha*sx.db[keys[i]]["energy"]/beta

      if EApprox is not None and len(keys)==3: 
       
         E0 = EApprox*np_sum(x0)
         k=0
         for i in PHASES : 
            E0 = E0 + x0[k]*sx.pdb[i]["energy"]
            k+=1

         #loop over remaining three key entries 
         for i in range(0,len(keys)):
            alpha=lam[i]
            #keep non-zero weights only
            if np.abs(alpha)>small:
               #find super cell of configuration 
               beta = FindSuperCell( keys[i], sx, xi[i],debug=debug)
               #print an error if something went wrong here
               if beta <0 :
                  print 'Error for key',keys[i],xi[i]
                  sys.exit(' Error in ReactionEnergy 2, beta is negative '+str(i))
               else:
                  E0 = E0 - alpha*sx.db[keys[i]]["energy"]/beta
         E0=-E0

   return round(E0,6)


def SWITCH( keys ) : 
#debug : 
   #K=['29-4:Ag.2.Li.3.Mn.8.O.16.c-1.','54-4:Ag.1.Li.5.Mn.16.O.32.', '28-3:Li.4.Mn.8.O.16.','1-1:Ag.1.']
   K=['52-4:Ag.3.Li.1.Mn.16.O.32.', '103-3:Ag.7.Mn.32.O.64.', '1-1:Ag.1.']
   #K=['55-4:Ag.3.Li.4.Mn.16.O.32.']
   debug = False
   if len(keys) == len(K) : 
      for i in range(0,len(K)):
         if keys[i] is not K[i] : 
            debug = True 

   return debug 


#-------------------------------------------------------------------------------
# Calculates reaction energy from keys 
#-------------------------------------------------------------------------------
def ReactionEnergy3D(keys,sx,lam,Xb,X,debug=None,NormIDX=None):
   from numpy import sum as np_sum
   from numpy import abs as np_abs
   
   E0=0
   # if four keys found we have a reaction of a configuration point, which is in
   # the database. The energy of the l.h.s. of the reaction is exact 
   if len(keys)==4:
     
      # energy of the l.h.s. 
      #we set E0 to the scaled energy of the supercell 
      E0=sx.db[keys[0]]["energy"]/sum(sx.db[keys[0]]["stoichiometry"]) 
      if NormIDX is not None:
         Xi = sx.db[keys[0]]["cart-coord"]
         Xbi = sx.db[keys[0]]["bary-coord"]
         Xn, idx =  GetConcentrations2([Xi[0],Xi[1]],NormIDX=NormIDX,ReturnIDX=True)
         alpha = Xn[ idx ] / Xbi[ idx ] 
         E0  = E0*sum(Xn)
 
      #loop over remaining three key entries 
      for i in range(1,4):
          #print an error if something went wrong here
          E0i=lam[i-1]*sx.db[keys[i]]["energy"]/sum(sx.db[keys[i]]["stoichiometry"])

          if NormIDX is not None:
             Xi = sx.db[keys[i]]["cart-coord"]
             Xbi = sx.db[keys[i]]["bary-coord"]
             Xn, idx =  GetConcentrations2([Xi[0],Xi[1]],NormIDX=NormIDX,ReturnIDX=True)
             alpha = Xn[ idx ] / Xbi[ idx ] 
             E0i = E0i*sum(Xn)

          E0 = E0 - E0i

   # if keys have length 3, an arbitrary coniguration is chosen 
   # for which the energy has to be approximated using the approximated Grand potential
   # energy
   if len(keys)==3:

      E0=0  
      k=0
      # sum individual contributions of each endpoint weighted with
      # bayecentric coordinates
      for i in PHASES : 
         E0 = E0 + Xb[k]*sx.pdb[i]["energy"]
         k+=1

      #add the Grand potential energy difference in bayecentric coordinates 
      E0=E0+GetGrandPotential(sx,X[0],X[1])

      #scale energy 
      if NormIDX is not None:
         Xn, idx =  GetConcentrations2([X[0],X[1]],NormIDX=NormIDX,ReturnIDX=True)
         alpha = Xn[ idx ] / Xb[ idx ] 
         E0  = E0*sum(Xn)

      #remove energy of remaining three entries 
      for i in range(0,3):

         E0i=lam[i]*sx.db[keys[i]]["energy"]/sum(sx.db[keys[i]]["stoichiometry"])

         if NormIDX is not None:
            Xi = sx.db[keys[i]]["cart-coord"]
            Xbi = sx.db[keys[i]]["bary-coord"]
            Xn, idx =  GetConcentrations2([Xi[0],Xi[1]],NormIDX=NormIDX,ReturnIDX=True)
            alpha = Xn[ idx ] / Xbi[ idx ] 
            E0i = E0i*sum(Xn)

         E0 = E0 - E0i


      E0=-E0

   return round(E0,6)

#-------------------------------------------------------------------------------
# Find supercell size, based on stoichiometry x0
# meaning multiplication factor of normalized stoichiometry x0, so that only
# integers x01, x02,x03 in A_x01 B_x02 C_x03 appear
#-------------------------------------------------------------------------------
def FindSuperCell(key,sx,x0,ReturnSuperCell=None,debug=None,Sdb=None):
   from fractions import Fraction
   from numpy import abs as np_abs
   from numpy import append as np_append
   #number of ions 
   nions0=sx.db[ key ]["nions"]

   # find integer stoichiometries
   s=IntegerStoichiometry( x0, sx, All=True)

   if debug: 
      print 'FindSuperCell:',s[0:min(3,len(s))],x0,nions0

   l=0
   for si in range(0,len(s)):
      l+=1
      SUM=0
      sc=s[si]
      k=0  

      for i in PHASES: 
         for m in range(0,len(sx.pdb[ i ]["species"])):
            SUM+=sc[k] * sx.pdb[ i ]["nion"][m]  
         k+=1

      if SUM == nions0 :  
         break
     
   if l>len(s) or l==0: 
      #return -1 if no supercell was found
      beta=-1
   else:
      k=0
      alpha=[]
      for i in PHASES:
         #compute multiplicity
         if x0[k]>0:
            alpha=np_append( alpha, sc[k]/x0[k])
         k+=1

     # if Sdb : 
     #    lcheck=False
     #    SUM =0
     #    for i in range(0,len(sc)):
     #       SUM+=np_abs(sc[i]-sx.db[ key ]["stoichiometry"][i])
     #    if SUM < small: 
     #       betaFloat=1

      betaFloat=alpha[0]
      for  i in range(0,len(alpha)):
         if np_abs(alpha[i]-betaFloat)>small:
             sys.exit('Error wrong stoichiometry ' ) 


      beta=int(betaFloat)
      beta=Fraction.from_float(betaFloat).limit_denominator(1).numerator

   if debug: 
      print 'multiplication factor',betaFloat,'->',beta

   if ReturnSuperCell:
      return int(beta),sc
   else:
      return int(beta)

#-------------------------------------------------------------------------------
# Checks if coordinates X=[ x, y ] are in triangle
#-------------------------------------------------------------------------------
def IsInsideTriangleOrClose(X,sx,tol):
   if len(X)==0 or len(X)==1 or len(PHASES)!=3:
      print len(X), X, len(PHASES), PHASES
      sys.exit('IsInsideTriangle reports probably wrong calling!')

   #squareroot of 3
   sq3=1.7320508075688772

   #x and y coordinates are:
   x=X[0]; y=X[1]
 
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
      d,xi,yi=ClosestPoint(x,y,sx)
     
      if d>tol:
         xi=x
         yi=y
      else: 
         IsInside = True
   else:
      IsInside=False
      xi=1000
      yi=1000
        
   return IsInside, xi,yi


#-------------------------------------------------------------------------------
# Checks if coordinates X=[ x, y ] are in triangle
#-------------------------------------------------------------------------------
def IsInsideTriangle(X,tol):
   if len(X)==0 or len(X)==1 or len(PHASES)!=3:
      print len(X), X, len(PHASES), PHASES
      sys.exit('IsInsideTriangle reports probably wrong calling!')

   #squareroot of 3
   sq3=1.7320508075688772

   #x and y coordinates are:
   x=X[0]; y=X[1]
 
   IsInside=False
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

   return IsInside

#-------------------------------------------------------------------------------
# closes point to sx
#-------------------------------------------------------------------------------
def ClosestPoint(x,y,sx):
  #find index in loaded data closest to 
  tmp = ((x-sx.hull.points[:,0])**2,(y-sx.hull.points[:,1])**2)
  idx = np.argmin(tmp[0]+tmp[1])
  tmp=sx.hull.points[idx,:]
  dist=(x-tmp[0])**2 + (y-tmp[1])**2

  return dist,tmp[0],tmp[1]

#-------------------------------------------------------------------------------
# returns a string in LaTeX format for the stoichiometry x0, based on the 
# endpoints defined in the dictionary pdb 
#-------------------------------------------------------------------------------
def CreateName(pdb,x0,lam=None,LaTeX=True,debug=None):
   from fractions import Fraction
   from numpy import abs as np_abs
   #number of significant digits
   ndigits=2
   #name  of compound to be build 
   name0=''
   #we have len(PHASES) total phases 
   phase=['']*len(PHASES)
   xi=[0]*len(PHASES)
   #counter to 
   for ic in [2,0,1]:
      i=PHASES[ic]
      xi[ic]=round(x0[ic],ndigits)
      if len(pdb[i]["species"])>1:
         k=0
         if np_abs(xi[ic]-1)>small and xi[ic]>small:
            phase[ic]+='['
         for spec in pdb[i]["species"]:
            if pdb[i]["nion"][k]==1:
               phase[ic]+=spec
            else:
               if LaTeX:
                  phase[ic]+=spec+'$_{'+str(pdb[i]["nion"][k])+'}$'
               else:
                  phase[ic]+=spec+str(pdb[i]["nion"][k])
            k+=1
         if np_abs(xi[ic]-1)>small and xi[ic]>small:
            phase[ic]+=']'
      elif len(pdb[i]["species"])==1:          
         k=0
         if pdb[i]["nion"][k]==1:
            phase[ic]+=pdb[i]["species"][0]
         else:
            if LaTeX:
               phase[ic]+=pdb[i]["species"][0]+'$_{'+str(pdb[i]["nion"][k])+\
              '}$'
            else:
               phase[ic]+=pdb[i]["species"][0]+str(pdb[i]["nion"][k])+']'

      if xi[ic]%1==0:
         if int(xi[ic])>1:
            if LaTeX:
               name0+=phase[ic]+'$_{'+str(int(xi[ic]))+'}$'
            else:
               name0+=phase[ic]+str(int(xi[ic]))
         else:
            if xi[ic]>small:
               name0+=phase[ic]
      else:
         if LaTeX:
            name0+=phase[ic]+'$_{'+str(xi[ic])+'}$'
         else:
            name0+=phase[ic]+str(xi[ic])
      ic+=1

   #in case we have a multiplication factor lam
   if lam is not None:
      # and it is not small
      if np.abs(lam)>0.00005 : 
         #in case it is 1
         if LaTeX:
            if np.abs(lam-1)>small:
               f=Fraction.from_float(lam).limit_denominator(100)
               p=f.numerator
               q=f.denominator
#               print lam, f, q 
               if q < 128 : 
                   name0='$\\frac{'+str(p)+'}{'+str(q)+'}\\times$ '+name0
               else:
                   name0='$'+str(round(lam,ndigits))+'\\times$ '+name0
         else:
            name0=str(round(lam,ndigits))+'x'+name0

   return name0

#------------------------------------------------------------------------------
# calculates cut coordinates (straight lines from one endpoint)
# kind =[ i, c,dir,npoints] defines the endpoint i and its start
#------------------------------------------------------------------------------
def GetCutCoordinates(pdb,kind,xmin,xmax):
   from numpy import append as np_append
   from numpy import abs as np_abs
   from numpy import dot as np_dot
   from numpy import sqrt as np_sqrt

   #this are the coordinates 
   result=None
   npoints=kind[3]

   #store concentration ratio
   c=kind[1]
   cc=c
   #make it safe
   if np_abs(c) <small:
      c=0.0005

   #scale back 
   c=c/(1+c) 

   # concentration 
   if kind[0] == 0 : 
      #calculate starting point 
      #direction is endpoint 0-> bottom right
      #these is the starting carthesian coordinates 
      v=c*pdb["BaryLatt"][:,1]
      #this is the final point 
      x_final = xmax/(xmax+kind[1]+1)
      #y_final = kind[1]/(xmax+kind[1]+1)
      y_final = cc/(xmax+cc+1)
      v_final = x_final*pdb["BaryLatt"][:,0]+y_final*pdb["BaryLatt"][:,1]
      #this is the direction vector of the cut
      n=v_final-v
      #normalize it to 1 
      norm = np_sqrt(np_dot(n,n))
      #shift point slightly  
      v=v+n/10000

   elif kind[0] == 1: 

      #calculate starting point 
      #direction is endpoint top->x-line
      #these is the starting carthesian coordinates 
      v=c*pdb["BaryLatt"][:,0]
      #this is the final point 
      x_final = xmax/(xmax+kind[1]+1)
      y_final = cc/(xmax+cc+1)
      v_final = x_final*pdb["BaryLatt"][:,1]+y_final*pdb["BaryLatt"][:,0]
      #this is the direction vector of the cut
      n=v_final-v
      norm = np_sqrt(np_dot(n,n))
      #shift point slightly  
      v=v+n/10000

   elif kind[0] == 2: 

      #calculate starting point 
      #direction is endpoint bottom left 
      #these is the starting carthesian coordinates 
      v=pdb["BaryLatt"][:,1]
      v=v+c*(pdb["BaryLatt"][:,0]-pdb["BaryLatt"][:,1])
      #this is the direction vector of the cut
      x_final = 1./(xmax+kind[1]+1)
      y_final = cc/(xmax+cc+1)
      v_final = y_final*pdb["BaryLatt"][:,1]+x_final*pdb["BaryLatt"][:,0]
      #this is the direction vector of the cut
      n=v_final-v
      #normalize it to 1 
      norm = np_sqrt(np_dot(n,n))
      #in case of an error return None
      #shift point slightly  
      v=v+n/10000

   #in case of an error return None
   if norm < small: 
      return result
   result=[v]

   #now loop over points 
   for i in range(1,npoints):
      if IsInsideTriangle(v+i*n/npoints,0.0001): 
         result = np_append(result, [v+i*n/npoints], axis=0)   
   
   result = np_append(result, [ v_final ] , axis  = 0 ) 
   

   return result


#------------------------------------------------------------------------------
# Calculates average of lattice constants along a set of carthesian coordinates
#------------------------------------------------------------------------------
def AverageLatticeConstants(coordinates,sx,NormIDX=None,status=None):
   from numpy import append as np_append
   from numpy import array as np_array
   from numpy import sum as np_sum
   from numpy import abs as np_abs
   from numpy import zeros as np_zeros
   from basis import GetConcentrations 

   npoints=len(coordinates)
   
   LattAvg=None
   #loop over coordinates 

   if status is not None: 
      status.Initialize(npoints)

   for i in range(0,npoints):
      #update status bar if wanted
      if status is not None: 
         status.UpdateProgressbar(' Calculating... ')

      x=coordinates[i,0]
      y=coordinates[i,1]
      #coordinates must be in the triangle 
      if IsInsideTriangle(coordinates[i,:],0.0001): 
         #get reaction
         keys,r,lam=CalculateReaction_OLD(sx,coordinates[i,:],\
             NormIDX = NormIDX, ReturnWeights=True)
         #for the following case a set of configurations and weights has been
         #found
         average=None
         #calculate average lattice constants 
         loss=np_zeros(len(lam))
         if len(keys)!=0 and len(keys)==len(lam) and len(keys)<4: 
            
            if NormIDX is not None:
               average,e = GetLatticeAverage(sx,keys,lam,NormIDX = NormIDX)
            else:
               average = GetLatticeAverage(sx,keys,lam)

            #add loss to average 
            nl=0
            for key in keys: 
               nc=0
               for j in PHASES: 
                  #check if endpoint is present
                  if key == sx.pdb[j]["key"]:
                     loss[nc] = lam[nl]
                  nc+=1
               nl+=1

         #contract lam to fit keys 
         elif len(keys) < len(lam): 

            lam_refined=[]
            #add loss to average 
            nc=0
            for i in range(0,len(lam)): 
               if np_abs(lam[i]) >small: 
                  lam_refined = np_append( lam_refined, lam[i] ) 
                  key = keys[nc]
                  #loop over phases, check if endpoint is present
                  for j in PHASES: 
                     #check if endpoint is present
                     if key == sx.pdb[j]["key"]:
                        loss[nc] = lam[i]
                  nc+=1
               
            #now check again 
            if len(keys)!=0 and len(keys)==len(lam_refined): 

               if NormIDX is not None:
                  average,e = GetLatticeAverage(sx,keys,lam_refined,NormIDX=NormIDX)
               else:
                  average = GetLatticeAverage(sx,keys,lam_refined)

         elif len(keys)==4 and len(lam)==3 :
            #calculate average lattice constants 
            if NormIDX is not None:
               average,e = GetLatticeAverage(sx,keys[1:],lam, NormIDX = NormIDX)
            else:
               average = GetLatticeAverage(sx,keys[1:],lam)

            #add loss to average 
            nl=0
            #calculate average lattice constants 
            for key in keys[1:]: 
               nc=0
               for j in PHASES: 
                  #check if endpoint is present
                  if key == sx.pdb[j]["key"]:
                     loss[nc] = lam[nl]
                  nc+=1
               nl+=1

         if average is not None:
           
          # if np_sum( average ) > small:
              x0=GetConcentrations2([x,y],NormIDX = NormIDX )
              x0=np_append(x0,average)
              x0=np_append(x0,e)
              x0=np_append(x0,loss)
              add = np_array(x0)
              if LattAvg is None:
                 LattAvg = [np_array(add)]
              else:
                 LattAvg = np_append( LattAvg, [add], axis=0)

         
            #  print add[0:2],add[8:]
   return LattAvg
#------------------------------------------------------------------------------
# averages lattice constans for structures keys in sx.db with weighting factors
# lam
#------------------------------------------------------------------------------
def GetLatticeAverage(sx,keys,lam,NormIDX=None):
   from numpy import sum as np_sum
   from numpy import abs as np_abs
   from numpy import dot as np_dot
   from numpy import cross as np_cross
   from numpy import sqrt as np_sqrt
   from numpy import array as np_array
   from basis import GetConcentrations 

   #intialize averate: a,b,c,V,formationE
   LattAvg=np_array( [0,0,0,0] ) 
   EAvg=0

   #check length of keys
   if len(lam)==3 and len(keys)==3:
      #exclude endpoints in the average, meaning use only weights suming to one
      loop = None
      #average over all three coordinates
      if np_abs(lam[0]+lam[1]+lam[2]-1)<small:
         loop=[0,1,2]
      #average over the first two
      elif np_abs(lam[0] + lam[1]-1)<small:
         loop=[0,1]

      #average over the first and last
      elif np_abs(lam[0] + lam[2]-1)<small:
         loop=[0,2]

      #average over the second and last
      elif np_abs(lam[1] + lam[2]-1)<small:
         loop=[1,2]
  
   elif len(lam)==2 and len(keys)==2:

      loop = None
      #average over all three coordinates
      if np_abs(lam[0]-1)<small:
         loop=[0]
      #average over the first two
      elif np_abs(lam[1] -1)<small:
         loop=[1]
  
   if loop is not None:
      #loop over correct indices 
      for i in loop :
         latt = sx.db[keys[i]]["lattice"]
         #this is the formation energy 
         z=sx.db[keys[i]]["formationE"]
         a=np_sqrt( np_dot( latt[0,:],latt[0,:] ) )
         b=np_sqrt( np_dot( latt[1,:],latt[1,:] ) )
         c=np_sqrt( np_dot( latt[2,:],latt[2,:] ) )

         #calculate volume of parallelepiped
         v = lam[i]*np_abs(np_dot( latt[0,:] , np_cross(latt[1,:],latt[2,:])))
         
         lsorted=[ a*lam[i], b*lam[i], c*lam[i] ]

         #stupid hack, don't know what should be done here
         if NormIDX is not None: 
            # get super cell 
            SC = SuperCell(keys[i],sx,NormIDX=NormIDX)
            alpha=1.

            for j in range(0,len(SC)): 
               if np_abs(SC[j])>0 : 
                  lsorted[j]=lsorted[j]/float(SC[j])
                  alpha=alpha*SC[j]
               else: 
                  lsorted[j]=0
                  e=0
            #scale Volume and Energy 
            v=v/alpha
            e=sx.db[keys[i]]["energy"]/alpha

            #calculate average energy 
            EAvg= EAvg + e*lam[i]
         lsorted.sort()
         lsorted.append(v)                 

         #add 
         LattAvg = LattAvg + np_array( lsorted ) 

   if NormIDX  is not None:
      return LattAvg,EAvg
   else:
      return LattAvg

#------------------------------------------------------------------------------
# determines hill and hull simplex 
# very fast version (uses matplotlib function path)
#------------------------------------------------------------------------------
def FindSimplices(sx, xp, yp,Hill=None, Hull=None):
   from numpy import array as np_array 
   from basis import GetBaryCentricCoordinates
   sq3=1.732050807568877

   #point has to be in the Triangle
   LF=False
   LG=False
   IsInside = False
   if yp<(sq3*xp+small): 
      LF = True
   if yp<(-sq3*xp+sq3+small): 
      LG = True

   if (LF is True) and (LG is True) : 
      IsInside = True 
   else:
      IsInside = False

   idx_i = -1    #index of hill simplex
   idx_u = -1    #index of hull simplex

   if Hill and IsInside: 
      # Delaunay simplex is found trivially with 
      p = np_array( (xp,yp) )
      idx_i = sx.hull.tri.find_simplex( p ) 
 

   # only if point can be in triangle at all 
   if Hull and IsInside: 
      # hull simplex is trickier, use matplotlibs path class for this 
      # and loop over all path
      LF=True
      for i in range(0, len(sx.hull.path) ) :
         p = sx.hull.path[i]
         if p.contains_point( [ xp,yp ] ): 
            LF=False
            idx_u = i 
            break 
 
      #shake point slightly    
      if LF:
         for k in range(0, len(sx.hull.path) ) :
            p = sx.hull.path[k]
            if p.contains_point( [ xp +small,yp ] ): 
               LF=False
               idx_u = k 
               break 

         if LF:
            for i in range(0, len(sx.hull.path) ) :
               p = sx.hull.path[i]
               if p.contains_point( [ xp-small,yp ] ): 
                  LF=False
                  idx_u = i 
                  break 

            #shake point slightly    
            if LF:
               for k in range(0, len(sx.hull.path) ) :
                  p = sx.hull.path[k]
                  if p.contains_point( [ xp,yp+small ] ): 
                     LF=False
                     idx_u = k 
                     break 

               if LF:
                  for i in range(0, len(sx.hull.path) ) :
                     p = sx.hull.path[i]
                     if p.contains_point( [ xp,yp-small ] ): 
                        LF=False
                        idx_u = i 
                        break 

                  #shake point slightly    
                  if LF:
                     for k in range(0, len(sx.hull.path) ) :
                        p = sx.hull.path[k]
                        if p.contains_point( [ xp+small,yp+small ] ): 
                           LF=False
                           idx_u = k 
                           break 

                     if LF:
                        for i in range(0, len(sx.hull.path) ) :
                           p = sx.hull.path[i]
                           if p.contains_point( [ xp+small,yp-small ] ): 
                              LF=False
                              idx_u = i 
                              break 
 
                        #shake point slightly    
                        if LF:
                           for k in range(0, len(sx.hull.path) ) :
                              p = sx.hull.path[k]
                              if p.contains_point( [ xp-small,yp+small ] ): 
                                 LF=False
                                 idx_u = k 
                                 break 

                           if LF:
                              for i in range(0, len(sx.hull.path) ) :
                                 p = sx.hull.path[i]
                                 if p.contains_point( [ xp-small,yp-small ] ): 
                                    LF=False
                                    idx_u = i 
                                    break 

                              if LF: 
                                 print 'NOT FOUND',xp, yp, GetBaryCentricCoordinates( [xp,yp] ) 
   return idx_i , idx_u 

  
#------------------------------------------------------------------------------
# determines hill and hull simplex 
# slow version (uses point in triangle check auxilary function)
#------------------------------------------------------------------------------
def FindSimplices_SLOW(sx, xp, yp,Hill=None, Hull=None):
   from numpy import array as np_array 
   from basis import GetBaryCentricCoordinates
   sq3=1.732050807568877

   def PointInsideTriangle2(pt,tri):
      '''checks if point pt(2) is inside triangle tri(3x2). @Developer'''
      a = 1/(-tri[1,1]*tri[2,0]+tri[0,1]*(-tri[1,0]+tri[2,0])+ \
          tri[0,0]*(tri[1,1]-tri[2,1])+tri[1,0]*tri[2,1])
      s = a*(tri[2,0]*tri[0,1]-tri[0,0]*tri[2,1]+(tri[2,1]-tri[0,1])*pt[0]+ \
          (tri[0,0]-tri[2,0])*pt[1])
      if s<0: return False
      else: t = a*(tri[0,0]*tri[1,1]-tri[1,0]*tri[0,1]+\
         (tri[0,1]-tri[1,1])*pt[0]+ (tri[1,0]-tri[0,0])*pt[1])
      return ((t>0) and (1-s-t>0))

   #point has to be in the Triangle
   LF=False
   LG=False
   IsInside = False
   if yp<(sq3*xp+small): 
      LF = True
   if yp<(-sq3*xp+sq3+small): 
      LG = True

   if (LF is True) and (LG is True) : 
      IsInside = True 
   else:
      IsInside = False

   idx_i = -1    #index of hill simplex
   idx_u = -1    #index of hull simplex

   if Hill and IsInside: 
      # Delaunay simplex is found trivially with 
      p = np_array( (xp,yp) )
      idx_i = sx.hull.tri.find_simplex( p ) 
 
   # only if point can be in triangle at all 
   if Hull and IsInside: 
      # hull simplex is trickier, use matplotlibs path class for this 
      # and loop over all path
      LF=False
      for i in range(0, len(sx.hull.inner_simplices) ) :
         tri = sx.hull.points[sx.hull.inner_simplices[i]]
         if PointInsideTriangle2( [xp,yp], tri ):
            idx_u = i 
            LF=True  #found simplex
            break 
         elif PointInsideTriangle2( [xp+small,yp], tri ):
            idx_u = i 
            LF=True  #found simplex
            break 
         elif PointInsideTriangle2( [xp-small,yp], tri ):
            idx_u = i 
            LF=True  #found simplex
            break 
         elif PointInsideTriangle2( [xp,yp+small], tri ):
            idx_u = i 
            LF=True  #found simplex
            break 
         elif PointInsideTriangle2( [xp,yp-small], tri ):
            idx_u = i 
            LF=True  #found simplex
            break 
         elif PointInsideTriangle2( [xp+small,yp+small], tri ):
            idx_u = i 
            LF=True  #found simplex
            break 
         elif PointInsideTriangle2( [xp-small,yp+small], tri ):
            idx_u = i 
            LF=True  #found simplex
            break 
         elif PointInsideTriangle2( [xp+small,yp-small], tri ):
            idx_u = i 
            LF=True  #found simplex
            break 
         elif PointInsideTriangle2( [xp-small,yp-small], tri ):
            idx_u = i 
            LF=True  #found simplex
            break 
 
      if not LF :    
         print ' NOT FOUND ', xp,yp,GetBaryCentricCoordinates( [xp,yp] )

   return idx_i , idx_u 


#------------------------------------------------------------------------------
# determines the Grand potential
#------------------------------------------------------------------------------
def GetGrandPotential(sx, xp, yp,NormIDX=None,WithDistance=None):
   from basis import distance_to_hull

   # determine simplex of hull and simplex of hill 
   idx_i, idx_u = FindSimplices( sx, xp, yp, Hill=True, Hull=True ) 

   zp=0
   # both must be positive when found 
   du=0
   di=0
   if idx_i >= 0 and idx_u >= 0 :
      du=distance_to_hull( xp, yp, sx.hull.tri.simplices[idx_i],\
         sx.hull.tri_points)

      di = distance_to_hull( xp, yp, sx.hull.inner_simplices[idx_u], \
         sx.hull.points,Sign=True)
      zp=di-du
 
   if NormIDX is not None : 
      Xn =  GetConcentrations2([xp,yp],NormIDX=NormIDX)
      zp = zp*sum( Xn ) 

   if WithDistance: 
      return zp,du,di
   else:
      return zp 

#------------------------------------------------------------------------------
# determines the energy per ion 
#------------------------------------------------------------------------------
def GetEnergyPerIons(sx, xp, yp, X, NormIDX=None):
   from basis import distance_to_hull
   from numpy import array, append, dot

   # determine simplex of hull and simplex of hill 
   idx_i, idx_u = FindSimplices( sx, xp, yp, Hill=True, Hull=True ) 

   zp=0
   # both must be positive when found 
   if idx_i >= 0 and idx_u >= 0 :
      zp=distance_to_hull( xp, yp, sx.hull.tri.simplices[idx_i],\
         sx.hull.tri_points)
      zp=zp+distance_to_hull( xp, yp, sx.hull.inner_simplices[idx_u], \
         sx.hull.points,Sign=True)

   if NormIDX is not None : 
      Xn =  GetConcentrations2([xp,yp],NormIDX=NormIDX)
      zp = zp*sum( Xn ) 

   mu=[]
   for k in [ 'A', 'B', 'C' ]: 
      mu = append( mu, sx.pdb[k]["energy"])
    
   E = zp + dot( X, mu )

   return  E

#------------------------------------------------------------------------------
# determine raction energy 
# if NormIDX = None the energy is given w.r.t. barycentric coordinates 
#------------------------------------------------------------------------------
def GetReactionEnergy(sx, xp, yp,NormIDX=None):
   from numpy import array as np_array

   X = np_array( [ xp, yp ] )
 
   keys, r = CalculateReaction_OLD(sx, X,NormIDX=NormIDX )

   E = r[len(r)-1]

   return E

#------------------------------------------------------------------------------
# determines all reaction energies energy 
#------------------------------------------------------------------------------
def GetAllReactionEnergies(sx,status=None):
   from numpy import append as np_append
 
   E=[]
   for i in range(0,len(sx.hull.points)):
      X = sx.hull.points[i,0:2]
      keys, r = CalculateReaction_OLD(sx, X,NormIDX=sx.NormIDX )
      E =np_append( E, r[len(r)-1] ) 
      if status is not None:
         k = int(float(i+1)/float(len(sx.hull.points))*100)
         status.SetStatusText(str(k)+'% done')

   status.SetStatusText('')
   return E 


#------------------------------------------------------------------------------
# determines the distance to hull 
#------------------------------------------------------------------------------
def GetStability(sx, xp, yp,NormIDX=None):
   from basis import distance_to_hull,GetBaryCentricCoordinates

   # determine simplex of hull and simplex of hill 
   idx_i, idx_u= FindSimplices( sx, xp, yp, Hill=True, Hull=None ) 

   zp=0
   # both must be positive when found 
   if idx_i >= 0 :
      zp=distance_to_hull( xp, yp, sx.hull.tri.simplices[idx_i],\
         sx.hull.tri_points)
 
   #normailze 
   if NormIDX is not None:
      Xb=GetBaryCentricCoordinates( [xp,yp] ) 
      Xn, idx =  GetConcentrations2([xp,yp],NormIDX=NormIDX,ReturnIDX=True)
      alpha = Xn[ idx ] / Xb[ idx ] 
      zp  = zp*sum(Xn)
  
   return zp 

#-------------------------------------------------------------------------------
# alternative for integer stoichiometries s = [x',y',z'] to si = [x0,y0,z0]
# where x',y',z' are floats and x0,y0,z0 are integers
#-------------------------------------------------------------------------------
def IntegerStoichiometry2(s,sx,All=None):
   from numpy import abs as np_abs 


   na=sum(sx.pdb["A"]["nion"])
   nb=sum(sx.pdb["B"]["nion"])
   nc=sum(sx.pdb["C"]["nion"])
   maxconc=max([na*SCMAX, nb*SCMAX, nc*SCMAX])

   alpha=[]
   sint=[]
   for i in range(0,len(PHASES)):
      I,R = divmod(s[i],1)
      alpha.append( 1 ) 

      #check only non-integer stoichios
      if abs(R)>small:
         
         #multiply rest by increasing number until 1 is almost reached: 
         for j in range(1,maxconc): 
            if np_abs(R*j-1)<0.01:  
               print 'fractional stoichio',j,s, round(R*j,1)
               alpha[i] = j 
               break
    
   for j in range(max(alpha),maxconc,max(alpha)):
      sint.append([int(round(j*s[0],1)),int(round(j*s[1],1)),int(round(j*s[2],1))])

   return sint

#------------------------------------------------------------------------------
# returns closest point in LoadedData[:,0:1] to (x,y) 
#------------------------------------------------------------------------------
def ClosestPoint(x,y,sx):
  from numpy import argmin as np_argmin
  tmp = ((x-sx.hull.points[:,0])**2,(y-sx.hull.points[:,1])**2)
  idx = np_argmin(tmp[0]+tmp[1])
  tmp=sx.hull.points[idx,:]

  return tmp

#------------------------------------------------------------------------------
# Get all concentrations from carthesian coordinates by first calculating 
# barycentric coordinates
#------------------------------------------------------------------------------
def GetConcentrations2(p,debug=None,NormIDX=None,ReturnIDX=None):
   from basis import GetBaryCentricCoordinates
   from numpy import abs as np_abs
   from numpy import argmax as np_argmax

   if debug : 
      print 'CARTESIAN',p 

   Xb = GetBaryCentricCoordinates(p) 
  
   if debug: 
      print 'BARYCENTRIC',Xb


   #get integer stoichiometry   
   S  = GetLCMVector( Xb, prec=1000000000 )   

   if debug: 
      print 'INTEGER',S
  
   #norm stoichiometry to specific index 
   if NormIDX is not None:
      d = S[ NormIDX ] 
      idx = NormIDX 
      #in case this coordinate is small use the maxium 
      if np_abs(Xb[ NormIDX ]) < 0.001 : 
         d = max( S[0],S[1],S[2] )   
         idx = int( np_argmax( S ) )

   # when no index given use barycentric stoichiometry
   else:
      S = Xb
      d = 1.
 
   if np_abs(d) > small: 
      X = [ float(S[0])/float(d), float(S[1])/float(d), float(S[2])/float(d) ]
   else:
      sys.exit('ERROR in GetConcentrations2 dviding through 0')

   if debug: 
      print 'RETURN',S,X,d

   if ReturnIDX and NormIDX is not None :
      return  X, idx 
   else:
      return  X
 

#------------------------------------------------------------------------------
# Computes reaction weights for 
#
#  A_x0 B_y0 C_z0 -> a[0] A_x1 B_y1 C_z1 + a[1] A_x2 B_y2 C_z2 +
#                    a[2] A_x3 B_y3 C_z3
#
# does everything in barycentric cooridinates
#------------------------------------------------------------------------------
def GetReactionWeights2(X,s,NormIDX=None):
   from numpy import ndarray as np_ndarray
   from numpy import linalg as np_linalg
   from basis import GetBaryCentricCoordinates

   b = np_ndarray(3, dtype=float, order='F')
   b=GetConcentrations2([X[0],X[1]],NormIDX=NormIDX)

   # reactants 
   a = np_ndarray(shape=(3,3), dtype=float, order='F')

   for j in range(0,3):
      xvec=GetConcentrations2([s[j,0],s[j,1]],NormIDX=NormIDX)
      for i in range(0,3):
         a[i,j]=xvec[i]

   lamb=np_linalg.solve(a,b)

   return lamb 

#------------------------------------------------------------------------------
# get list of prime factors for inter number n 
# use as p = list(GetPrimeFactors(n)) to get a list p of prime factors for n
#------------------------------------------------------------------------------
def GetPrimeFactors(n):
    for p in GetPrimes(n):
        while n%p==0:
            yield p
            n=n//p
            if n==1: return
def GetPrimes(n):
   if n < 2: return
   yield 2
   plist = [2]
   for i in range(3,n):
       test = True
       for j in plist:
           if j>n**0.5:
               break
           if i%j==0:
               test = False
               break
       if test:
           plist.append(i)
           yield i

#------------------------------------------------------------------------------
# calculates the super cell of key in database sx.db w.r.t. NormIDX 
#------------------------------------------------------------------------------
def SuperCell(key,sx,NormIDX=None,debug=None):
   from numpy import amax as np_amax
   from numpy import amin as np_amin
   from numpy import dot as np_dot
   from numpy import sqrt as np_sqrt
   nions=sum(sx.db[key]["nion"])
   
   if NormIDX is not None : 
      nkey=sx.pdb[PHASES[NormIDX]]["key"]

      #find length of vectors
      norm_uc=[0]*3
      norm_sc=[0]*3
      for i in range(0,3):
        v=sx.db[nkey]["lattice"][i,:]
        v_=sx.db[key]["lattice"][i,:]
        norm_uc[i]=np_sqrt(np_dot(v,v))
        norm_sc[i]=np_sqrt(np_dot(v_,v_))
      sc=[ int(round(norm_sc[0]/norm_uc[0])), \
           int(round(norm_sc[1]/norm_uc[1])), \
           int(round(norm_sc[2]/norm_uc[2])) 
         ]
   
      SC = sc
      if debug:
         print '*******************'
         print 'length of vectors (uc)',norm_uc
         print 'length of vectors (sc)',norm_sc
         print 'supercell :',SC
         print '*******************'

   #in case no norm index is given return supercell w.r.t. to every endpoint
   else: 
      SC=[]
      for j in range(0,len(PHASES)): 
         nkey=sx.pdb[PHASES[j]]["key"]
         norm_uc=[0]*3
         norm_sc=[0]*3
         for i in range(0,3):
           v=sx.db[nkey]["lattice"][i,:]
           v_=sx.db[key]["lattice"][i,:]
           norm_uc[i]=np_sqrt(np_dot(v,v))
           norm_sc[i]=np_sqrt(np_dot(v_,v_))
         sc=[ int(round(norm_sc[0]/norm_uc[0])), \
              int(round(norm_sc[1]/norm_uc[1])), \
              int(round(norm_sc[2]/norm_uc[2])) 
            ]
	 for k in range(0,3):
	    if sc[k] == 0 :
               sc=[0,0,0]

         SC.append( sc ) 

   if debug:
      print SC 
   return SC 

#------------------------------------------------------------------------------
# returns carthesian coordinates from reciprocal coordinates
#------------------------------------------------------------------------------
def ReciprocalToCarthesian(key,sx):
   from numpy import zeros as np_zeros
   coord = np_zeros((len(sx.db[key]["coordinates"]),3))
   for i in range(0,len(sx.db[key]["coordinates"])):
      for j in range(0,3): 
         for k in range(0,3): 
            coord[i,k] = coord[i,k]+sx.db[key]["coordinates"][i,j]*sx.db[key]["lattice"][k,j]
 
   return coord

#------------------------------------------------------------------------------
# returns True if x1 is a multiply of x0
#------------------------------------------------------------------------------
def GetScalingFactor(x0,x1):
   from numpy import abs
 
   if len(x0) != len(x1) : 
      return 0 
 
   if len(x0) == 1 : 
      l1 = float(x1[0])/float(x0[0])
   elif len(x0) >1 : 
      for i in range(0,len(x0)-1): 
         l1 = float(x1[i])/float(x0[i])
         l2 = float(x1[i+1])/float(x0[i+1])
         if abs( l1-l2 ) > small : 
            return -1

   if abs(l1) < small : 
        return -1
   else:
        return l1

# ------------------------------------------------------------------------------
# Calculates formation energy from point (x,y)
# ------------------------------------------------------------------------------
def GetFormulaAll(Project,x,y):

   x0 = GetConcentrations2([x,y])
   s=[]
   for i in [2,0,1]:
      p=PHASES[i]
      if abs(x0[i]>0.001):
         s.append( [ x0[0]/x0[i], x0[1]/x0[i], x0[2]/x0[i] ] )
      else:
         s.append( [ 0,0,0] )

   k=-1
   text = []
   for i in [2,0,1]:
      p=PHASES[i]
      k+=1
      if abs(s[k][i]>small):
        if abs(s[k][i])<100:
           name = CreateName( Project.pdb, s[k], LaTeX = False)
           text.append(name)
        else:
           text.append(' ')
      else:
        text.append(' ')

   return text 

# ------------------------------------------------------------------------------
# Calculates formation energy from point (x,y)
# ------------------------------------------------------------------------------
def CalculateFormation(Project,x,y,f):
   from numpy import abs

   x0 = GetConcentrations2([x,y])

   if len(f) > 0:
       text ="\n   Database entry {} found for: ".format(f)
   else:
      text ="\n   Approximated point:"

   name = CreateName(Project.pdb, x0,LaTeX=False) 

   text += " {:}".format(name)
 
   E = GetGrandPotential(Project,x,y)
   text += "\n     Formation energy/Ion: {:.4f} [eV]".format(E)
   s=[]
   for i in [2,0,1]:
      p=PHASES[i]
      E+=Project.pdb[p]["energy"]*x0[i]
      if abs(x0[i]>small):
         s.append( [ x0[0]/x0[i], x0[1]/x0[i], x0[2]/x0[i] ] )
      else:
         s.append( [ 0,0,0] )

   k=-1
   if len(f)>0:
      force=Project.dbF[f]["maxForce"]
      text += "\n     Structure has maximum residual force"+\
           " {:.4f} {:.4f} {:.4f} [eV/A]".format(force[0],force[1],force[2])

   for i in [2,0,1]:
      p=PHASES[i]
      k+=1
      if abs(s[k][i]>small):
         Ei = E/x0[i]
         name = CreateName( Project.pdb, s[k], LaTeX = False)
         text += "\n     Energy for {:}:   {:.6} eV".format(\
         name,Ei)
   
   return text

# ------------------------------------------------------------------------------
# Calculates energy in bary centric coordinates 
# ------------------------------------------------------------------------------
def GetBarycentricCoord(Project,x,y):
   x0 = GetConcentrations2([x,y])
   return x0

# ------------------------------------------------------------------------------
# Calculates segregation iseg from point (x,y)
# ------------------------------------------------------------------------------
def CalculateSegregation(Project,x,iseg, debug=None):
   from basis import distance_to_hull
   from numpy import array, linalg,abs
   from numpy import zeros as np_zeros

   #store simplices
   simplices2 = Project.hull.tri.simplices.copy() 
   simplices = Project.hull.inner_simplices.copy() 
   allpoints2=np_zeros( ( len(Project.hull.points),3) ) 
   allpoints2[:,0]=Project.hull.points[:,0]
   allpoints2[:,1]=Project.hull.points[:,1]
   allpoints2[:,2]=Project.hull.distance[:]

   lsimplices = True


   # imediatley return in case of endpoint
   if abs(x[ 2+iseg ]-1)<small :
      if debug: 
	  print " Endpoint segregated ",x,[2+iseg]
      return 1.0

   idx_i, idx_u = FindSimplices( Project, x[0], x[1], Hill=True, Hull=True ) 
   if idx_i >= 0: 
      E=-distance_to_hull( x[0], x[1], simplices2[idx_i], allpoints2)

   if abs(E)<small:
      if debug: 
	  print " stable Endpoint selected ", x, E
      return 0.0

   a = np_zeros(shape=(3,3))
   if idx_u >= 0 and idx_u < len(simplices):
      keys=[]
      k=-1
      lseg=-1
      for i in simplices[idx_u]:
         keys.append(Project.hull.keys[i])
         k+=1
         xi=Project.db[Project.hull.keys[i]]["bary-coord"]
         for j in range(0,3):
            a[j,k] = xi[j]

         # decide if its an endpoint
	 if abs( xi[ iseg ] - 1 ) < small : 
            lseg = k 
         
      # decide if endpoint included in reaction
      if lseg >= 0 : 
         lamb=linalg.solve(a,x[2:5])
         if debug: 
	    print " reaction weights and keys:", lamb, keys
      else:
         if debug: 
            print " no endpoint present",keys
         return 0 
     
   return lamb[lseg] 

# ------------------------------------------------------------------------------
# Calculates reaction from point (x,y)
# ------------------------------------------------------------------------------
def CalculateReaction(Project,x,y):
   from basis import distance_to_hull
   from numpy import array, linalg,abs
   from numpy import zeros as np_zeros

   #store simplices
   simplices2 = Project.hull.tri.simplices.copy() 
   simplices = Project.hull.inner_simplices.copy() 
   allpoints2=np_zeros( ( len(Project.hull.points),3) ) 
   allpoints2[:,0]=Project.hull.points[:,0]
   allpoints2[:,1]=Project.hull.points[:,1]
   allpoints2[:,2]=Project.hull.distance[:]


   x0 = GetConcentrations2([x,y])
   name = [ CreateName(Project.pdb, x0,LaTeX=False) ]
   lsimplices = True
   # skip in case of endpoint 
   for i in [2,0,1]:
      p = PHASES[i]
      if abs(x0[i]-1)<small :
         return "     {:} is stable".format(name[0])

   idx_i, idx_u = FindSimplices( Project, x, y, Hill=True, Hull=True ) 
   if idx_i >= 0: 
      E=-distance_to_hull( x, y, simplices2[idx_i], allpoints2)

   if abs(E)<small:
      return "     {:} is stable".format(name[0])

   text = "     {:} --> ".format(name[0])


   b=array( x0 )
   a = np_zeros(shape=(3,3))
   if idx_u >= 0 and idx_u < len(simplices):
      keys=[]
      k=-1
      for i in simplices[idx_u]:
         keys.append(Project.hull.keys[i])
         k+=1
         xi=Project.db[Project.hull.keys[i]]["bary-coord"]
         name.append( CreateName(Project.pdb, xi ,LaTeX=False) )
         for j in range(0,3):
            a[j,k] = xi[j]
      lamb=linalg.solve(a,b)

   i=-1
   for key in keys:
      i+=1
      if i==0:
         if abs( lamb[i] ) > small:
            text+=" {:.2} x {:}".format(lamb[i],name[i+1])
      else:
         if abs(lamb[i]) < small :
            continue
         text+="+ {:.2} x {:}".format(lamb[i],name[i+1])

   text += "\n     Reaction energy: {:.4f} eV ( {:.4} kJ/mol )".\
   format(E,E*kJmol )
     
   return text


#------------------------------------------------------------------------------
# Get carthesian coordinates  from barycentric ones
#------------------------------------------------------------------------------
def GetCarthesianCoordinates(x0):
   x=0
   y=0
   for i in [2,0,1]:
      x = x+ x0[i] * PHASE_COORD[i,0]
      y = y+ x0[i] * PHASE_COORD[i,1]

   return [ x,y ]

def GetGrandDifferential(Project, xi, xf, debug=None):
   """ computes differential of Grand potential from barycentric point x0 """
   from numpy import array, copy, append
   dx=0.001

   Xa = copy( xi[2:5] ) 
   xa=array( GetCarthesianCoordinates( Xa ), copy=True ) 
   Ga=GetGrandPotential( Project, xa[ 0 ], xa[ 1 ] )

   Xb = copy( xf[2:5] ) 
   xb=array( GetCarthesianCoordinates( Xb ), copy=True ) 
   Gb=GetGrandPotential( Project, xb[ 0 ], xb[ 1 ] )

   lcmbv, lcmb  = GetLCMVector( Xb, GetLCM=True ) 
   lcmav, lcma  = GetLCMVector( Xa, GetLCM=True ) 
   if debug : 
      print Ga, Gb, Xb, lcmbv, lcmb, lcmav, lcma  
   D=[]
   DX=copy( Xb ) 
   for k in [ 0, 1, 2] :
      #
      # normalize w.r.t. to k 
      #
      DiffG = Gb*lcmb*lcmav[k] - Ga*lcma*lcmbv[k]
      dx = lcmb*lcmav[k] - lcma*lcmbv[k] 

      DG=0
      if abs( dx ) > small : 
	 DG = DiffG/dx 
      D = append( D, DG ) 

      if debug : 
	 print k, Gb*lcmb*lcmav[k], Ga*lcma*lcmbv[k], dx , DG 

   return DX, D 


def GetGrandDifferentialOLD(Project, xj, dx=None, debug=None):
   """ computes differential of Grand potential from barycentric point x0 """
   from numpy import array, copy, append
   # displacement is small 
   if dx is None: 
      dx=0.01

   Xb = copy( xj[2:5] ) 
   x0=array( GetCarthesianCoordinates( Xb ), copy=True ) 
   G0=GetGrandPotential( Project, x0[ 0 ], x0[ 1 ] )
   D=[]
   DX=copy( Xb ) 
   for k in [ 0, 1, 2] :
      #
      # displacement in k-direction, calculated from stoichiometry 
      #
      Xk = copy( Xb ) 
      Xk[ k ] = Xk[ k ] + dx 
      Xb_k = copy( Xk ) 
      Xb_k[ 0 ] = Xk[ 0 ]/( 1.+dx ) 
      Xb_k[ 1 ] = Xk[ 1 ]/( 1.+dx ) 
      Xb_k[ 2 ] = Xk[ 2 ]/( 1.+dx ) 
      Xc_k=array( GetCarthesianCoordinates( Xb_k ), copy=True ) 
      Gk=( 1. + dx )*GetGrandPotential( Project, Xc_k[ 0 ], Xc_k[ 1 ] )
      Dk = (Gk - G0)/dx
      D = append(D,  Dk ) 
   if debug : 
      print xj, x0, G0, D
   return DX, D 

def Get2DFormationEnergies(Project, xj, debug=None):
   """ computes 1D formation energies """
   from numpy import copy, zeros
 
   result=zeros( 2 ) 
   
   # current point
   x0=GetCarthesianCoordinates( xj )
   Ej,du,di=GetGrandPotential( Project, x0[ 0 ], x0[ 1 ], WithDistance=True )

   result[0]= Ej
   result[1]= di
   if debug : 
      print Ej,du,di
   return result

def Get1DFormationEnergies(Project, xj, dx=None, debug=None):
   """ computes 1D formation energies """
   from numpy import copy, zeros
 
   result=zeros( 2 ) 

   # current point
   Xj = copy( xj[2:5] ) 
   x0=GetCarthesianCoordinates( Xj )
   Ej,du,di=GetGrandPotential( Project, x0[ 0 ], x0[ 1 ], WithDistance=True )

   result[0]= Ej
   result[1]= di
   if debug : 
      print Ej,du,di
   return result

def GetLatticeParameter(Project, x, debug = None):
   """ computes lattice parameter for a barycentric point from reaction
       weights """
   from basis import distance_to_hull
   from numpy import array, linalg,abs
   from numpy import zeros 
   from numpy import dot, sqrt, cross, cos, arccos
   
   def angle(a,b): 
     num= dot( a,b)
     den= sqrt( dot( a,a) )*sqrt( dot( b,b ))
     if abs( den ) > small:
        return arccos( num/den )
     else:
        return 0.

   #store simplices
   simplices2 = Project.hull.tri.simplices.copy() 
   simplices = Project.hull.inner_simplices.copy() 
   allpoints2=zeros( ( len(Project.hull.points),3) ) 
   allpoints2[:,0]=Project.hull.points[:,0]
   allpoints2[:,1]=Project.hull.points[:,1]
   allpoints2[:,2]=Project.hull.distance[:]
   Null = array( [0,0,0] )

   # phase permutation order
   phaselist=[0,1,2]
   # result
   L = zeros( shape=( 3,7 ) )

   x0 = GetConcentrations2(x)
   lsimplices = True
   idx_i, idx_u = FindSimplices( Project, x[0], x[1], Hill=True, Hull=True ) 
   if idx_i >= 0: 
      E=-distance_to_hull( x[0], x[1], simplices2[idx_i], allpoints2)

   b=array( x0 )
   a = zeros(shape=(3,3))
   if idx_u >= 0 and idx_u < len(simplices):
      keys=[]
      k=-1
      for i in simplices[idx_u]:
         keys.append(Project.hull.keys[i])
         k+=1
         xi=Project.db[Project.hull.keys[i]]["bary-coord"]
         for j in range(0,3):
            a[j,k] = xi[j]
      lamb=linalg.solve(a,b)

      k=-1
      l=[0]*12
      if debug :
         print " ************************* "  , keys, lamb 
      for key in keys: 
	 k+=1
         # normalize w.r.t. endpoints 
         if debug:
    	    print " PART:",key,lamb[k],sum(lamb)
	 latt=Project.db[key]["lattice"]
	 avec=latt[0,:]
	 bvec=latt[1,:]
	 cvec=latt[2,:]
	 SC= SuperCell( key, Project )     

         # normalize w.r.t. endpoints 
	 alpha=0.
	 beta=0.
	 gamma=0.
         a=0.
         b=0.
         c=0.
	 if debug: 
            print 'SuperCell',SC[j], avec, bvec, cvec
         for j in phaselist:
            p = PHASES[j]
	    if abs( SC[j][0])>small:
               Avec=avec/float(SC[j][0]) 
            else:
	       Avec=Null
	    if abs( SC[j][1])>small:
               Bvec=bvec/float(SC[j][1]) 
            else:
	       Bvec=Null
	    if abs( SC[j][2])>small:
               Cvec=cvec/float(SC[j][2]) 
            else:
               Cvec=Null


            a_ = sqrt( dot( Avec, Avec ) ) 
            b_ = sqrt( dot( Bvec, Bvec ) )
            c_ = sqrt( dot( Cvec, Cvec ) )
            alpha_ = angle( Avec, Bvec ) 
            beta_  = angle( Bvec, Cvec )
            gamma_ = angle( Cvec, Avec )
	    if debug: 
               print 'SuperCell',SC[j], a_,b_,c_,alpha_,beta_,gamma_

            a = sqrt( dot( Avec, Avec ) ) * lamb[ k ]  
            b = sqrt( dot( Bvec, Bvec ) ) * lamb[ k ]  
            c = sqrt( dot( Cvec, Cvec ) ) * lamb[ k ]  
            #v = abs(dot( Avec , cross( Bvec,Cvec)))

            # average angles 
            alpha = angle( Avec, Bvec ) * lamb[ k ] 
            beta  = angle( Bvec, Cvec ) * lamb[ k ] 
            gamma = angle( Cvec, Avec ) * lamb[ k ] 

 

            if debug :
               print " ",Project.pdb[p]["key"],SC[j],a_,b_,c_, alpha_, beta_, gamma_
               print "      ",a,b,c, alpha, beta, gamma

            L[j]=L[j]+array( [ a, b, c, alpha, beta, gamma, 0. ] )
       
      # volume of parallelepiped 
      for j in phaselist: 
         a=L[j][0]; b=L[j][1]; c=L[j][2]; al=L[j][3]; be=L[j][4]; ga=L[j][5];
         v = a*b*c*sqrt( 1.+2*cos(al)*cos(be)*cos(ga) - \
                  ( cos(al)*cos(al)+cos(be)*cos(be)+cos(ga)*cos(ga) ) )
         L[j][6]=v
      if debug:
         print "L",L
   return L

