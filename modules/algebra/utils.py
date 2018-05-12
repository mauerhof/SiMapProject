import numpy as np


def solveTriangularMatrix(a, B):
	n = len(a)
	x = np.empty([n])
	for i in range(n-1, -1, -1):
		x[i] = a[i]/B[i,i]
		for j in range(n-1, i, -1):
			x[i] = x[i] - B[i,j]*x[j]/B[i,i]


	return x
	


# ~ B = np.empty([2,2])
# ~ B[0,0] = 1
# ~ B[0,1] = 2
# ~ B[1,0] = 0
# ~ B[1,1] = 3
# ~ print B

# ~ a = np.array([2,6])
# ~ print a

# ~ print solveTriangularMatrix(a,B)
