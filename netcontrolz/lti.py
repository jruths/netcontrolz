"""
This module provides convenience methods to work with the linear time-invariant (LTI)
dynamical system model for network control:
    x(t+1) = Ax(t) + Bu(t)
or the continuous time analog
    xdot(t) = Ax(t) + Bu(t)
"""
from numpy import zeros
from numpy.random import random

__all__ = ['controls_idx','input_matrix','input_matrix_','matrix_realization','fix_diagonals','kalman_generic_rank','kalman_generic_rank_']

def controls_idx(G,controls):
    """
    Convenience function to return the list of tuples of control indices of a network
    ``G`` corresponding to the list of tuples of node objects in ``controls``.
    """
    return [tuple(G.node_idx(a) for a in c) for c in controls]

def input_matrix(G,controls,m=None):
    """
    Returns the ``n`` x ``m`` input matrix B corresponding to the list of tuples
    ``controls`` of node indices that are driven by the controls, where ``n`` is
    the number of nodes in the network ``G`` (if ``G`` is not compact, then the
    size ``n`` may be larger than the number of nodes). If ``m`` is
    not specified, then ``m`` is taken as the number of controls (i.e., the
    length of ``controls``).
    """
    return input_matrix_(G.max_node_idx+1,controls_idx(controls),m)

def input_matrix_(n,controls,m=None):
    """
    Returns the ``n`` x ``m`` input matrix B corresponding to the list of tuples
    ``controls`` of node indices that are driven by the controls. If ``m`` is
    not specified, then ``m`` is taken as the number of controls (i.e., the
    length of ``controls``).
    """
    if m is None:
        m = len(controls)
    B = zeros((n,m))
    for j,control in enumerate(controls):
        for ci in control:
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

    .. math::
        M_{ii} = -(d + \sum_{j=1}^n M_{ij})

    which is simply just a perturbation ``d`` away from the row sum of ``M``.
    Selecting ``d > 0`` ensures the eigenvalues of a symmetric matrix are negative.
    """
    for i in range(len(M)):
        M[i,i] = -( d + M[i,:].sum() )
    return M

def kalman_generic_rank(G,controls,repeats=None):
    """
    Finds the reachability corresponding to a graph ``G`` (adjacency matrix A) and its
    ``controls`` (a list of tuples of node objects) by brute force computing the Kalman rank condition:

    .. math::
        rank [B\ AB\ A^2B\ A^3B\ \cdots\ A^{(n-1)}B].

    In order to compute the rank generically, we generate random entries for A and B,
    subject to their zero/non-zero sparsity patterns and compute the true rank. We
    repeat this ``repeats`` times (default is 100) and return the largest value.
    """
    return kalman_generic_rank_(G,controls_idx(controls),repeats)

def kalman_generic_rank_(G,controls,repeats=None):
    """
    Finds the reachability corresponding to a graph ``G`` (adjacency matrix A) and its
    ``controls`` (a list of tuples of node indices) by brute force computing the Kalman rank condition:

    .. math::
        rank [B\ AB\ A^2B\ A^3B\ \cdots\ A^{(n-1)}B].

    In order to compute the rank generically, we generate random entries for A and B,
    subject to their zero/non-zero sparsity patterns and compute the true rank. We
    repeat this ``repeats`` times (default is 100) and return the largest value.
    """
    rank = 0
    N = G.max_node_idx+1 # there could be some missing indexes in the graph (N > n)
    n = G.num_nodes
    A = G.matrix().T

    if repeats == None:
        repeats = 100

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
