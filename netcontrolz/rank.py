"""
This module provides implementations of rank-related functions.
"""
import zen
from zen import DiGraph, maximum_matching_
from zen.exceptions import type_check


def kalman_generic_rank(G,controls,repeats=100):
	"""
	Finds the reachability corresponding to a graph (adjacency matrix A) and its
	controls (a matrix B) by brute force computing the Kalman rank condition:
		rank [B AB A^2B A^3B ... A^(n-1)B].
	In order to compute the rank generically, we generate random entries for A and B,
	subject to their zero/non-zero sparsity patterns and compute the true rank. We
	repeat this "repeats" times (default is 100) and return the largest value.
	"""
	from numpy import array, zeros, nonzero, eye, dot
	from numpy.linalg import matrix_rank
	from numpy.random import rand as nprand
	
	rank = 0
	N = G.max_node_idx+1 # there could be some missing indexes in the graph (N > n)
	n = G.num_nodes
	A = G.matrix().T
	
	# build the B matrix
	m = len(controls)
	B = zeros((N,m))
	for i in range(m):
		for d_idx in controls[i]:
			B[d_idx,i] = 1
	
	nonzero_idxs = nonzero(A>0)
	num_nonzero_idxs = len(nonzero_idxs[0])
	for r in range(repeats):
		# create a randomized instance of A
		A1 = zeros((N,N))
		A1[nonzero_idxs] = nprand(num_nonzero_idxs)
		
		# compute the controllability matrix
		An = eye(N)
		C = zeros((N,m*N))
		C[:,0:m] = B
		for i in range(1,N):
			An = dot(An,A1)
			C[:,i*m:i*m+m] = dot(An,B)
			
		# generic rank is the max of all instance ranks
		new_rank = matrix_rank(C)
		if new_rank > rank:
			rank = new_rank
	
	return rank


def generic_rank(G):
	"""
	Returns the generic rank of the graph ``G``.
	"""
	type_check(G,DiGraph)
	
	matched_edges = maximum_matching_(G)
	
	return len(matched_edges)