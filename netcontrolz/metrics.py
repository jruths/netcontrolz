"""
This module provides methods to compute various controllability metrics of networks
and their control configurations.
"""
from numpy import zeros, diag, dot, abs, concatenate, eye
from numpy.random import random
from numpy.linalg import cond
from scipy import allclose
from scipy.integrate import odeint
from scipy.linalg import expm, eig, inv, qr, eigh
from warnings import warn

__all__ = ['finite_horizon_gramian','infinite_horizon_gramian','finite_horizon_discrete_time_gramian']

# CURRENTLY ASSUMES SYMMETRIC MATRIX
def finite_horizon_gramian(A,B,T):
    """
    Returns the finite horizon (finite time) Gramian matrix when the system
    is driven to the origin:

    .. math::
        H = e^{-A*T} W e^{-A'*T},

    where :math:`'` stands for transpose

    .. math::
        W = \int_0^T e^{A*t} B B' e^{A'*t} dt

    Driving a system from :math:`x_0` to :math:`x_T` uses energy:

    .. math::
        \mathcal{E} = v' W^{-1} v,  \qquad v = x_T - e^{A*T}x_0

    If :math:`x_T = 0` (origin) then this simplifies to

    .. math::
        \mathcal{E} = x_0' H^{-1} x_0

    Bounds on the energy :math:`\mathcal{E}` can be established by the max and min eigenvalues
    of :math:`H` (:math:`\lambda_1 > \lambda_2 > \cdots > \lambda_n`):

    .. math::
        1/\lambda_1 \leq \mathcal{E} \leq 1/\lambda_n.

    Note that matrices ``A`` and ``B`` should be realizations, not sparsity patterns.
    """
    n = len(A)
    n2 = n**2
    d,V = eigh(A)
    if d.real.max() > 0 and d.real.min() > 0:
        print 'A is positive definite.'
    elif d.real.max() < 0 and d.real.min() < 0:
        print 'A is negative definite.'
    else:
        print 'A is neither positive or negative definite.'
    D = diag(d)

    if cond(V) > 10**12:
        warn("the condition number of the eigenvector matrix is large (>10^12).",RuntimeWarning)
    VBBV = dot(dot(V.T,dot(B,B.T)),V)

    def grammian_int(x,t):
        expDt = expm(D*t)
        Wdot = dot(dot(expDt,VBBV),expDt).reshape((n2,))
        return concatenate( (Wdot.real,Wdot.imag) )
    def grammian_int_imag(x,t):
        expDt = expm(D*t)
        return dot(dot(expDt,VBBV),expDt).reshape((n2,)).imag

    Wvec = odeint(grammian_int, grammian_int(0,0).reshape((2*n2,)), [0,T])[1,:]
    W = Wvec[0:n2].reshape((n,n)) + 1j * Wvec[n2:].reshape((n,n))

    expDT = expm(-D*T)
    H = dot( dot( dot(V,expDT), W), dot(expDT,V.T) )

    return H

# CURRENTLY ASSUMES SYMMETRIC MATRIX
def finite_horizon_discrete_time_gramian(A,B,T):
    """
    Returns the finite horizon (finite time, ``T``) Gramian matrix when the system
    is driven from the origin:

    .. math::
        W = \sum_{k=0}^{T-1} A^k B B' A^k

    Driving a system from the origin to :math:`x_T` uses energy:

    .. math::
        \mathcal{E} = x_T' W^{-1} x_T,

    Bounds on the energy :math:`E` can be established by the max and min eigenvalues
    of :math:`W` (:math:`\lambda_1 > \lambda_2 > \cdots > \lambda_n`):

    .. math::
        1/\lambda_1 \leq \mathcal{E} \leq 1/\lambda_n.

    Note that matrices ``A`` and ``B`` should be realizations, not sparsity patterns.
    """
    BB = dot(B,B.T)
    W = BB
    Ai = eye(len(A))
    for i in range(1,T):
        Ai = dot(Ai,A)
        W += dot(dot(Ai,BB),Ai)
    return W

# Attempt for directed networks
# def finite_horizon_gramian(A,B,T):
#     """
#     Returns the finite horizon (finite time) grammian matrix when the system
#     is driven to the origin:
#         H = e^{-A*T} W e^{-A'*T}, where ' stands for transpose
#         W = \int_0^T e^{A*t} B B' e^{A'*t} dt
#     Driving a system from x0 to xT uses energy:
#         E = v' W^{-1} v,  with v = xT - e^{A*T}x0
#     If xT = 0 (origin) then this simplifies to
#         E = x0' H^{-1} x0
#     Bounds on the energy E can be established by the max and min eigenvalues
#     of H:
#         1/eigmax(H) <= E <= 1/eigmin(H).
#     Note that matrices A and B should be realizations, not sparsity patterns.
#     """
#     n = len(A)
#     n2 = n**2
#     d,V = eig(A)
#     if d.real.max() > 0 and d.real.min() > 0:
#         print 'A is positive definite.'
#     elif d.real.max() < 0 and d.real.min() < 0:
#         print 'A is negative definite.'
#     else:
#         print 'A is neither positive or negative definite.'
#     D = diag(d)
#
#     if cond(V) > 10**12:
#         warn("the condition number of the eigenvector matrix is large (>10^12).",RuntimeWarning)
#     Vinv = inv(V)
#     diff = abs(dot(dot(V,D),Vinv) - A)
#     print diff.max().max()
#     print diag( dot(V,V.T) )
#     VBBV = dot(dot(Vinv,dot(B,B.T)),Vinv.T)
#
#     def grammian_int(x,t):
#         expDt = expm(D*t)
#         Wdot = dot(dot(expDt,VBBV),expDt).reshape((n2,))
#         return concatenate( (Wdot.real,Wdot.imag) )
#     def grammian_int_imag(x,t):
#         expDt = expm(D*t)
#         return dot(dot(expDt,VBBV),expDt).reshape((n2,)).imag
#
#     Wvec = odeint(grammian_int, grammian_int(0,0).reshape((2*n2,)), [0,T])[1,:]
#     W = Wvec[0:n2].reshape((n,n)) + 1j * Wvec[n2:].reshape((n,n))
#
#     expDT = expm(-D*T)
#     H = dot( dot( dot(V,expDT), W), dot(expDT,V.T) )
#
#     return H

# CURRENTLY ASSUMES SYMMETRIC MATRIX
def infinite_horizon_gramian(A,B):
    """
    Returns the infinite horizon (final time goes to infinity) Gramian matrix
    corresponiding to a symmetric stable matrix ``A``.
    """
    n = len(A)
    d, V = eigh(A)
    D = diag(d)

    VBBV = dot( V.T , dot( B , dot(B.T,V) ) )
    M = zeros((n,n))
    for i in range(n):
        for j in range(n):
            M[i,j] = VBBV[i,j] / (d[i]+d[j])
    G = dot( V , dot( M , V.T ) )
    return G
