"""
This module provides implementations of rank-related functions.
"""
from zen import DiGraph, maximum_matching_
from zen.exceptions import type_check

from numpy import array, zeros, nonzero, eye, dot
from numpy.linalg import matrix_rank
from numpy.random import rand as nprand

__all__ = ['generic_rank']

def generic_rank(G):
    """
    Returns the generic rank of the graph ``G``. This is equivalent to the size of the maximum matching of the graph ``G``; also equivalent to the difference between the number of nodes in the network and the number of dilations, see method ``netcontrolz.num_dilations``.

    """
    type_check(G,DiGraph,'only directed graphs are supported')

    matched_edges = maximum_matching_(G)

    return len(matched_edges)
