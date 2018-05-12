from minirats.utils.py.cellutils import py_cell_utils as cu
from minirats.utils.py.readwrite import ramses_info as ri
from matplotlib import pyplot as plt
import numpy as np


def plotOutput(ramDir, snapshot, indices, lmax, center, radius, xymin, xymax, name):
	
	#ncells = cu.count_cells(ramDir,snapshot,lmax,center,radius) 
	#cells,cell_pos,cell_l = cu.read_cells_hydro(ramDir,snapshot,lmax,ncells,indices,center,radius,True)
	
	#Tgas = cells[:,1]/cells[:,0]*1.67e-24/1.38062e-16*(unit_length/unit_time)**2
	ri.rd_info(snapshot, ramDir)
	
	



plotOutput('/Users/mauerhof/Documents/RamsesFiles/idealized/G8/lmax11Ngroups5/', 12, [1,5,7,8,9], 11, [0.5,0.5,0.5], 1, 0.47, 0.53, 'test')
