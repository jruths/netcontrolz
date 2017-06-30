"""
This module provides convenience methods to work with the linear time-invariant (LTI)
dynamical system model for network control:
    x(t+1) = Ax(t) + Bu(t)
or the continuous time analog
    xdot(t) = Ax(t) + Bu(t)
"""
from numpy import zeros
from numpy.random import random

__all__ = ['controls_idx','input_matrix_','matrix_realization','fix_diagonals']

def controls_idx(G,controls):
    """
    Helper function to return the list of tuples of control indices of a network
    ``G`` corresponding to the list of tuples of node objects in ``controls``.
    """
    controls_ = []
    for control in controls:
        controls_.append( tuple( [G.node_idx(c) for c in control]  ) )
    return controls_

def input_matrix_(n,controls_,m=None):
    """
    Returns the ``n``x``m`` input matrix B corresponding to the list of tuples
    ``controls_`` of node indices that are driven by the controls. If ``m`` is
    not specified, then ``m`` is taken as the number of controls (i.e., the
    length of ``controls_``).
    """
    if m is None:
        m = len(controls_)
    B = zeros((n,m))
    for j,control_ in enumerate(controls_):
        for ci in control_:
            B[ci,j] = 1
    return B

def matrix_realization(M,a,b):
    """
    Returns a matrix with the same sparsity pattern as the matrix ``M``
    where nonzero elements of ``M`` are uniformly randomly selected from the
    interval [a,b].
    """
    n,m = M.shape
    R = M.copy()
    for i in range(n):
        for j in range(m):
            if R[i,j] != 0:
                R[i,j] = (b-a)*random() + a
    return R

def fix_diagonals(M,d):
    """
    Returns M with diagonal elements set to be:
        M_{ii} = -(d + \sum_{j=1}^n M_{ij})
    which is simply just the a perturbation d away from the row sum of M.
    Selecting d > 0 ensures the eigenvalues of a symmetric matrix are negative.
    """
    for i in range(len(M)):
        M[i,i] = -( d + M[i,:].sum() )
    return M
