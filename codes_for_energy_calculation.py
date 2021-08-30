print "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@"
print "#############Lennard_Jonnes_potential_precalculation_on grid############"
print "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@"
print
# ATOMS_DICT_CONSISTING_OF_LENNARD-JONNES_CONSTANTS
ATOMS={'C-C':[{'e':0.150},{'C12':2516582.400},{'C6':1228.800000}], \
    'C-N':[{'e':0.155},{'C12':1198066.249},{'C6':861.634784}], \
    'C-O':[{'e':0.173},{'C12':820711.722},{'C6':754.059521}], \
    'C-S':[{'e':0.173},{'C12':2905899.052},{'C6':1418.896022}], \
    'C-H':[{'e':0.055},{'C12':29108.222},{'C6':79.857949}], \
    'N-C':[{'e':0.155},{'C12':1198066.249},{'C6':861.634784}], \
    'N-N':[{'e':0.160},{'C12':540675.281},{'C6':588.245000}], \
    'N-O':[{'e':0.179},{'C12':357365.541},{'C6':357365.541}], \
    'N-S':[{'e':0.179},{'C12':1383407.742},{'C6':994.930149}], \
    'N-H':[{'e':0.057},{'C12':10581.989},{'C6':48.932922}], \
    'O-C':[{'e':0.173},{'C12':820711.722},{'C6':754.059521}], \
    'O-N':[{'e':0.179},{'C12':357365.541},{'C6':505.677729}], \
    'O-O':[{'e':0.200},{'C12':230584.301},{'C6':429.496730}], \
    'O-S':[{'e':0.200},{'C12':947676.268},{'C6':870.712934}], \
    'O-H':[{'e':0.063},{'C12':6035.457},{'C6':39.075098}], \
    'S-C':[{'e':0.173},{'C12':2905899.052},{'C6':1418.896022}], \
    'S-N':[{'e':0.179},{'C12':1383407.742},{'C6':994.930149}], \
    'S-O':[{'e':0.200},{'C12':947676.268},{'C6':870.712934}], \
    'S-S':[{'e':0.200},{'C12':3355443.200},{'C6':1638.400000}], \
    'S-H':[{'e':0.063},{'C12':33611.280},{'C6':92.212017}], \
    'S-H':[{'e':0.063},{'C12':33611.280},{'C6':92.212017}], \
    'H-C':[{'e':0.055},{'C12':29108.222},{'C6':79.857949}], \
    'H-N':[{'e':0.057},{'C12':10581.989},{'C6':48.932922}], \
    'H-O':[{'e':0.063},{'C12':6035.457},{'C6':39.075098}], \
    'H-S':[{'e':0.063},{'C12':33611.280},{'C6':92.212017}], \
    'H-H':[{'e':0.020},{'C12':81.920},{'C6':2.560000}]}
#EMPTY_LIST_FOR_STORING_QUERY_PROTEIN_COORDINATES_AND_ATOMS        
empty=[]
x=[]
#FILE_CATING_BEGINS_HERE
with open('C:\\Users\\user\\Desktop\\ligand_eug.pdb') as f:
    for i in f:
        if i.startswith('HETATM'):
            columns=i.split()
            x=columns[5:8]
            x.append(columns[-1])
            empty.append(x)
#FINDING_THE_EXTREME_POINTS_OF_3D_COORDINATES
x_axis, y_axis, z_axis, atoms =zip(*empty)
xmax=max(map(float, x_axis))
xmin=min(map(float, x_axis))
ymax=max(map(float, y_axis))
ymin=min(map(float, y_axis))
zmax=max(map(float, z_axis))
zmin=min(map(float, z_axis))
#FINDING_THE_CENTER
x_center = int( xmax + xmin )/2
y_center = int( ymax + ymin )/2
z_center = int( zmax + zmin )/2
sp = 1.0
#FINDING_NPTS
x_npts=int((xmax-(xmin))/sp)
y_npts=int((ymax-(ymin))/sp)
z_npts=int((zmax-(zmin))/sp)
npts=[x_npts+1, y_npts+1, z_npts+1]
center=[x_center, y_center, z_center]
grid_points =[]
#CREATING_GRID_COORDINATES_FOR EACH_GRID_POINT
for z in xrange(-(npts[2]/2), npts[2]/2 + 1):
    for y in xrange(-(npts[1]/2), npts[1]/2 + 1):
        for x in xrange(-(npts[0]/2), npts[0]/2 + 1):
            grid_points.append([x*sp+center[0], y*sp+center[1], z*sp+center[2]])
            
distance=''
lig_atoms=['C']
Atom_grid_LJ_pot=[]
FE_vdW_coeff=0.1485
import math
print "<<<<<<<<<<<<<<CALCULATION_OF_LENNARD_JONNES_GRID_POTENTIAL>>>>>>>>>>>>>>"
print
for atom_type in lig_atoms:
    grid_LJ_pot=[]                
    for i in grid_points:
        Total_LJ_pot=0
        for j in empty:
            distance = math.sqrt((i[0] - float(j[0]))**2 + (i[1] - float(j[1]))**2 + (i[2] - float(j[2]))**2)
            Total_LJ_pot+=(((ATOMS[j[3]+ '-' +atom_type][1]['C12'])/float (distance)**12)-2.0*((ATOMS[j[3]+ '-' +atom_type][1]['C12'])/float(distance)**6))*FE_vdW_coeff
        grid_LJ_pot.append(Total_LJ_pot)
    Atom_grid_LJ_pot.append(grid_LJ_pot) 
print "@@@@@@@@@Display_OF_LJ_POTENTIAL_FOR_CORRESPONDING_GRID_POINTS@@@@@@@@@@" 
print   
for i in range(1):
	counter=i
	print   '       grid_points'+'                               '+ 'potential'
	for j in grid_points:
		print  '       '+str (j)+'           :::          '+str (Atom_grid_LJ_pot[0][counter])
		counter+=1
# GRID_COORDINATE_FILE_CREATION        
with open('C:\\Users\\user\\Desktop\\grid.pdb','w') as f:
    for i, a in enumerate(grid_points):
        x,y,z = a
        f.write("ATOM %6d  O   GRD A   1    %8.3f%8.3f%8.3f  1.00  0.00           O\n"%(i,x,y,z))     
print "#######################SUCCESSFUL_COMPLETION############################"        