"""
This module provides methods to compute various controllability metrics of networks
and their control configurations.
"""

from numpy import zeros
from numpy.random import random

def control_idx(G,controls):
	controls_ = []
	for control in controls:
		controls_.append( tuple( [G.node_idx(c) for c in control]  ) )
	return controls_

def control_matrix_(n,m,controls_):
	B = zeros((n,m))
	for j,control_ in enumerate(controls_):
		for ci in control_:
			B[ci,j] = 1
	return B

def matrix_realization(M,a,b):
	n,m = M.shape
	R = M.copy()
	for i in range(n):
		for j in range(m):
			if R[i,j] != 0:
				R[i,j] = (b-a)*random() + a
	return R

def grammian(G,controls,T):
	return grammian_(G,control_idx(controls),T)

def grammian_(G,controls_,T):
	A = G.matrix().T
	n = len(A)
	B = B_(n,controls_)
	
	return
