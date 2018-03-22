from minirats.utils.py.cellutils import py_cell_utils as cu
from matplotlib import pyplot as plt
import numpy as np
import math

# read cells
lmax = 6 # read all levels 

idens,ipre,ixhii,ixheii,ixheiii = 1,6,7,8,9 
center = [0.5,0.5,0.5]
radius = 1.0

ramDir = '/Users/mauerhof/Documents/Ramses/stromgren/'
timestep = 11
ncells = cu.count_cells(ramDir,timestep,lmax,center,radius) 
cells,cell_pos,cell_l = cu.read_cells_hydro(ramDir,timestep,lmax,ncells,[idens,ipre,ixhii,ixheii,ixheiii,10,14,18],center,radius,True)
#print (ncells)

#print cells[10]


nx, ny = cu.get_map_nxny(lmax, 0, 1, 0, 1)
weights = np.ones_like(cells[:,2])
HII, wHII = cu.make_map(lmax, True, 0, 1, 0, 1, 0, 0.01, cells[:,2], cells[:,0], cell_pos[:,0], cell_pos[:,1], cell_pos[:,2], cell_l, nx, ny)
HeII, wHeII = cu.make_map(lmax, True, 0, 1, 0, 1, 0, 0.01, cells[:,3], cells[:,0], cell_pos[:,0], cell_pos[:,1], cell_pos[:,2], cell_l, nx, ny)
HeIII, wHeIII = cu.make_map(lmax, True, 0, 1, 0, 1, 0, 0.01, cells[:,4], cells[:,0], cell_pos[:,0], cell_pos[:,1], cell_pos[:,2], cell_l, nx, ny)

radius = np.empty([ncells])
for i in range(ncells):
	radius[i] = math.sqrt(cell_pos[i][0]**2 + cell_pos[i][1]**2 + cell_pos[i][2]**2)

radius = radius*3.086e21
print HII
#print radius.max()
#print TempMap.min()
#print TempMap.max()
#print cells[:,2].min()
#print cells[:,2].max()
plt.imshow(HII, origin='lower')
#plt.imshow(HeII, origin='lower')
#print TempMap[32][32:64:1]
plt.colorbar()
#plt.savefig('stromgren_HII'+str(timestep)+'.png')
plt.show()
#print cells[:,3].size
#print cells[:,3][900:1000:1] + cells[:,4][900:1000:1]
#print HeII[0,0]
#print HeIII[0,0]
#print cell_l
#print cells[0,4]
#print cell_pos.size
#print len(cell_pos)
#print cell_pos[:,0].size
#print cell_pos[0:98:1][:]
#print cell_pos.max()
#print cells[:,3].max()
#print TempMap.min()
#print TempMap.max()
