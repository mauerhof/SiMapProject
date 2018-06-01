import numpy as np
import math
from matplotlib import pyplot as plt



def order(orderVector, toBeOrdered):
	
	
	data = np.column_stack((orderVector, toBeOrdered))
	data = data[np.argsort(data[:,0])]
	
	return data[:,1]



def orderCell(cell_pos, toBeOrdered):
	
	orderVector = np.zeroes_like(cell_pos[:,0])
	for i in range(3):
		orderVector = orderVector + math.sqrt(i)*cell_pos[:,i] 
	
	return order(orderVector, toBeOrdered)





def ScatterPlot_equal(x,y,size,xlabel,ylabel,name):

	xy_line = (0, 1)
	plt.style.use('dark_background')
	fig, ax = plt.subplots(figsize=(8,8))
	ax.scatter(x, y, s=1, marker='.', facecolor='0.5', lw=0)
	plt.xlim(0,1)
	plt.ylim(0,1)
	ax.set_xlabel(xlabel, fontsize=19)
	ax.set_ylabel(ylabel, fontsize=19)
	ax.plot(xy_line, 'r--', label='y=x')
	ax.legend(loc="lower right")
	plt.savefig(name+'.png')
	plt.close()
