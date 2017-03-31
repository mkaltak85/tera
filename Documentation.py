import wx
import  wx.lib.wxpTag


global documentation
documentation={}

documentation[ "About" ]=\
"""
<h2>About TERA</h2>
<p><b>
TERA is a tool for drawing and analyzing phase diagrams of <a
href="Introduction">ternary systems</a>.
</h3></p>
"""

documentation[ "RAW_File_Example" ] = \
"""
<h2>Data file example for
Li<sub>x</sub>Ag<sub>y</sub>[Mn<sub>8</sub>O<sub>16</sub>]<sub>z</sub></h2>
<p>
<FONT FACE="courier">
# data points and energies of Li_x Ag_y [Mn_8 O_16]_z <br>
# x      y      z         energy   <br> 
# Li-bcc<br>
  1.00   0.00   0.00      -3.6686  <br> 
# Ag-fcc<br>
  0.00   1.00   0.00      -0.6496  <br> 
# Li_0   Ag_y   Mn_8 O_16          <br> 
  0.00   0.00   1.00    -180.1014  <br> 
  0.00   0.25   1.00    -180.2740  <br> 
  0.00   0.50   1.00    -180.4880  <br> 
  0.00   0.75   1.00    -180.7380  <br> 
  0.00   1.00   1.00    -180.9910  <br> 
  0.00   1.25   1.00    -181.1750  <br> 
  0.00   1.50   1.00    -181.3540  <br> 
  0.00   1.75   1.00    -181.5230  <br> 
  0.00   2.00   1.00    -181.5410  <br> 
# Li_x   Ag_0   Mn_8  O_16         <br> 
  0.25   0.00   1.00    -181.7970  <br> 
  0.50   0.00   1.00    -183.5350  <br>  
  1.00   0.00   1.00    -186.7690  <br> 
  1.50   0.00   1.00    -190.4000  <br>      
  2.00   0.00   1.00    -193.9110  <br> 
  3.00   0.00   1.00    -200.4550  <br> 
  4.00   0.00   1.00    -207.3958  <br> 
  6.00   0.00   1.00    -219.8714  <br> 
  7.00   0.00   1.00    -225.2778  <br> 
  8.00   0.00   1.00    -232.6754  <br> 
# Lix    Ag_0.25 Mn_8 O_16 <br>
  0.25   0.25   1.00    -182.00430<br>
  0.75   0.25   1.00    -185.25    <br> 
  1.75   0.25   1.00    -192.5439  <br> 
  3.75   0.25   1.00    -205.7510  <br> 
  7.00   0.25   1.00    -225.3510  <br> 
# Lix    Ag_0.5 Mn_8 O_16<br>
  0.25   0.50   1.00    -182.2738  <br> 
  0.50   0.50   1.00    -183.9965  <br> 
  1.00   0.50   1.00    -187.2189  <br> 
  1.50   0.50   1.00    -191.0144  <br> 
  2.00   0.50   1.00    -194.2698  <br> 
  2.50   0.50   1.00    -197.8514  <br> 
  3.00   0.50   1.00    -200.6810  <br> 
  3.50   0.50   1.00    -204.2649  <br> 
  4.00   0.50   1.00    -207.3635  <br> 
  4.50   0.50   1.00    -210.3676  <br> 
  5.50   0.50   1.00    -216.6371  <br> 
  8.00   0.50   1.00    -231.9477  <br> 
# Lix    Ag_0.75 Mn_8 O_16<br>
  0.25   0.75   1.00    -182.4406  <br> 
# Li_x Ag Mn_8 O_16<br>
#  0.00  1.00  1.00    -180.6861   <br>
  0.50   1.00   1.00    -184.5830  <br> 
  1.00   1.00   1.00    -188.1094  <br> 
  1.50   1.00   1.00    -190.9553  <br> 
  2.00   1.00   1.00    -194.8214  <br> 
  3.00   1.00   1.00    -201.3244  <br>  
  4.00   1.00   1.00    -207.4398  <br> 
  5.00   1.00   1.00    -213.6364  <br> 
  6.00   1.00   1.00    -219.9136  <br> 
  7.00   1.00   1.00    -226.5611  <br> 
  8.00   1.00   1.00    -231.3060  <br> 
# Li_x   Ag_0.25 Mn_8 O_16<br>
  0.25   1.25   1.00    -182.8798  <br> 
  0.75   1.25   1.00    -186.1182  <br> 
  0.50   1.25   1.00    -184.4618  <br> 
  1.00   1.25   1.00    -187.8036  <br> 
#Li_x    Ag_1.5 Mn_8 O_16<br>
  0.25   1.50   1.00    -183.1601  <br> 
  0.50   1.50   1.00    -185.0125  <br> 
  1.00   1.50   1.00    -188.2271  <br> 
  1.50   1.50   1.00    -191.5372  <br> 
  2.00   1.50   1.00    -194.5702  <br> 
  3.00   1.50   1.00    -201.2488  <br> 
  4.00   1.50   1.00    -207.6336  <br>  
  4.50   1.50   1.00    -210.7589  <br> 
  5.50   1.50   1.00    -216.5106  <br> 
  8.00   1.50   1.00    -230.9297  <br> 
# Li_x  Ag_1.75 Mn_8 O_16<br>
  0.25   1.75   1.00    -183.0822  <br> 
# Li_x   Ag_2   Mn_8 O_16<br>
  1.00   2.00   1.00    -187.9647  <br> 
  2.00   2.00   1.00    -195.1177  <br> 
  3.00   2.00   1.00    -201.4613  <br> 
  4.00   2.00   1.00    -207.9336  <br> 
  5.00   2.00   1.00    -213.4962  <br> 
  6.00   2.00   1.00    -220.8288  <br> 
  7.00   2.00   1.00    -225.4183  <br>     
  8.00   2.00   1.00    -230.2316  <br> 
# these structures are taken from materialsproject.org <br>
# but relaxed with same parameters as used for calculations above <br>
# Li_x   Ag_y <br>
  1.00   1.00   0.00      -4.7877  <br>  
  1.00   3.00   0.00      -6.2509  <br>  
  3.00   1.00   0.00     -12.4126  <br>   
</FONT>
</p>
"""


documentation[ "Supported_File_Types" ]=\
"""
</p>
<p><b>*.ter|*.tera:</b>  TERAs own data format. 
</p>
<p><b>*.txt|*.dat:</b> Textfiles with four columns, where 
lines starting with # are ignored:
<FONT FACE="courier">
<br>
<br>
    # x    y     z       Energy <br>
    1.000  1.200 0.0000  -6.24  <br>
    0.000  0.000 1.0000  -3.73  <br>
    1.000  0.000 0.0000  -1.18  <br>
     ...
<br>
</FONT><br>
An example file is given <a href="RAW_File_Example">here</a>.
<p><b> *.vasp:</b> Directory containing files in POSCAR format of <a
href="VASP">VASP</a> with an additional entry for the corresponding energy of
the structure:  
<FONT FACE="courier">
<br>
<br>
Li1Ag1 <br>
 1.0 <br>
     3.143834489 -0.000000076 -0.000000086 <br>  
    -0.000000077  3.143834335  0.000000068 <br>  
    -0.000000087  0.000000068  3.143834290 <br>  
  Li Ag  <br>
 1 1    <br>
carthesian <br>
      1.57192      1.57192      1.57192    <br>  
      0.00000      0.00000      0.00000    <br>  
 Energy         :    -4.787713 [eV]  <br>
 Force          :     0.010202      0.001005      0.020005 [eV/A]  <br>
 DEGENERACY     :  1 
<br>
<br>
</FONT>
The lines with <FONT FACE="courier"> Force, DEGENERACY </FONT> are optional and can be left
out (see <a href="Residual_Force">Residual forces</a>). 
</p>
<p> In both cases the energies for the endpoints A, B and C must be present. 
"""


documentation[ "Introduction" ]=\
"""
<h2>Introduction</h2>
<p>
<a href="About">TERA</a> allows to visualize and analyze the phase space of
ternary systems. A ternary system contains three components A,B and C, where
each component stands for an element in the periodic table (e.g. H, Li,  ... )
or a chemical composition (e.g. H<sub>2</sub>, H<sub>2</sub>O, ... ). 
</h3></p>
<p> To analyze a ternary system, the configuration space has to be sampled
accurately, that is the energy of A<sub>x</sub>B<sub>y</sub>C<sub>z</sub> has to
be calculated for a sufficiently large set of distininct stoichiometries
(x,y,z), including the energy of the endpoints A, B and C. 
For more accuracy one should investigate the energies of iso-stoichiometric
configuration of different symmetry, see <a href="Disorder">disorder</a> for
more information.
</p>
<p> 
<b>Stable Phases</b>
</p>
<p> 
Stable compositions are found using the <a href="Convex_Hull_Method">convex hull
method</a>. First, the formation energy &#x394;F(x,y,z) of every configuration
A<sub>x</sub>B<sub>y</sub>C<sub>z</sub> is calculated as
</p>
<p> 
&#x394;F(x,y,z) = E(x,y,z) - x' &#x3BC<sub>x</sub>-y &#x3BC<sub>y</sub>-z
&#x3BC<sub>z</sub>,
</p>
<p> 
where E(x,y,z) is the energy of A<sub>x</sub>B<sub>y</sub>C<sub>z</sub> and
&#x3BC<sub>x</sub> indicates the chemical potential of A particles and so on.
Note that compounds with negative formation energy &#x394;F(x,y,z) &#60 0 are taken
into account and those with &#x394;F(x,y,z) &#62 0 are ignored. In the next
step, a radial projection 
</p>
<p> 
&#x3C6: (x,y,z)  &#8614; (x',y',z') = (x,y,z)/(x+y+z)
</p>
<p> 
is used to single out stoichiometrically distinct configurations mapping the
latter ones onto an equilateral triangle in the 111 plane of R<sup>3</sup>
(barycentric coordinates). This allows to draw the corresponding formation
energy &#x394;F(x',y',z') = &#x394;F(x,y,z)/(x+y+z) in the 111 direction and to
determine all minima of the function geometrically with the <a
href="Convex_Hull_Method">convex hull</a>. Due to the convexity of &#x3C6, the
minima of &#x394;F(x,y,z), i.e. stable points, will be on the convex hull of
&#x394;F(x',y',z'). 
<p> 
A small tutorial can be found <a href="Tutorial">here</a>.
</p>
"""

documentation[ "Tutorial" ]=\
"""
<h2>Tutorial</h2>
<p> Starting TERA will give following user interface:
<p> <wxp module="TTutorial" class="BlankGUI"> </wxp> </p>
</p>
<p><b>
Import Data
</b></h3></p>
<p> Click on <wxp module="TTutorial" class="Open"></wxp> to import data.
Alternatively one can use the standard keyboard shortcut "CRTL+O". This will
show the import form:
<p> <wxp module="TTutorial" class="ImportForm"> </wxp> </p>
<p>
The user has to specify either a file or a directory path. An example file
is given <a href="RAW_File_Example">here</a>. All supported file formats
are explained on this <a href="Supported_File_Types">page</a>.
Furthermore the endpoints of the ternary system "Phase A", "Phase B" and "Phase
C" have to be specified. The chemical formulae are not crucial for raw text
files, but are important if a directory with structure files is chosen. 
The force text fields are explained <a href="ImportFormEnhanced">elsewhere</a>.
</p>

<p>
<b>Convex Hull
</b></h3></p>
<p>
Click on <wxp module="TTutorial" class="DataPoints2"></wxp>, this should draw
the convex hull of the data points
</p>
<p>
<wxp module="TTutorial" class="ConvexHull"></wxp>
</p>
"""

