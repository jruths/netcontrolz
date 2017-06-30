"""
This script reproduces the main results from the paper below.

The outstanding problem of controlling complex networks is relevant to many areas
of science and engineering, and has the potential to generate technological
breakthroughs as well. We address the physically important issue of the energy
required for achieving control by deriving and validating scaling laws for the
lower and upper energy bounds. These bounds represent a reasonable estimate of
the energy cost associated with control, and provide a step forward from the
current research on controllability toward ultimate control of complex networked
dynamical systems.

..seealso::

    Gang Yan, Jie Ren, Ying-Cheng Lai, Choy-Heng Lai, and Baowen Li. Controlling
    Complex Networks: How Much Energy Is Needed? Phys. Rev. Lett. 108, 218703
    (2012).

    url: http://link.aps.org/doi/10.1103/PhysRevLett.108.218703
==========================================
STILL BEING IMPLEMENTED
==========================================
"""
import zen
import netcontrolz
from numpy import abs, zeros, log10
from scipy.linalg import eigvals, inv, eigvalsh
import matplotlib.pyplot as plt
from numpy.random import choice

# 500 nodes, average degree of 20
n = 500
p = 20.0/n
G = zen.generating.erdos_renyi(n,p)

As = G.matrix().T
Bs = netcontrolz.input_matrix_(n,[(0,)])
# find realizations of each matrix
A = netcontrolz.matrix_realization(As,0,1)
A = netcontrolz.fix_diagonals(A,0.25)
B = Bs #netcontrolz.matrix_realization(Bs)

H = netcontrolz.finite_horizon_grammian(A,B,1)
d = eigvals(H)
eigmax = d.real.max()
eigmin = d.real.min()

print eigmin, eigmax

print 'Emin = %1.2f' % (1/eigmax)
