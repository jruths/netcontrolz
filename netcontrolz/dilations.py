"""
This module provides tools to identify dilations in a directed network.
"""
from zen import DiGraph
from rank import generic_rank

__all__ = ['num_dilations','all_in_neighbors','all_in_neighbors_']

def num_dilations(G):
	"""
	Returns the number of dilations in the directed network ``G``. This is also 
	equivalent to the smallest number of controls that are required to structurally 
	control the graph ``G`` (except in the degenerate case when there are zero 
	dilations, but one control).
	"""
	type_check(G,DiGraph)
		
	return len(G) - generic_rank(G)

def all_in_neighbors(G,S):
	"""
	Returns the set of nodes (node objects) in the directed graph ``G`` that have 
	directed edges pointing to nodes in the set ``S``, a set of node objects.
	"""
    nbrs = set([])
    for nobj in S:
        nbrs.update(set(G.in_neighbors(nobj)))
    return nbrs
    
def all_in_neighbors_(G,S):
	"""
	Returns the set of nodes (node indices) in the directed graph ``G`` that have 
	directed edges pointing to nodes in the set ``S``, a set of node indices.
	"""
    nbrs = set([])
    for nidx in S:
        nbrs.update(set(G.in_neighbors_(nidx)))
    return nbrs
