"""
This module provides tools to identify dilations in a directed network.
"""
from zen import DiGraph
from rank import generic_rank

__all__ = ['num_dilations']

def num_dilations(G):
	"""
	Returns the number of dilations in the directed network ``G``. This is also 
	equivalent to the smallest number of controls that are required to structurally 
	control the graph ``G``.
	"""
	type_check(G,DiGraph)
		
	return max([len(G) - generic_rank(G),1])

