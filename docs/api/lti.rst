Linear Time-Invariant (LTI) Systems
===========
Much of network control assumes the state dynamics follow a linear time-invariant model.

.. math::
    x(t+1) = Ax(t) + Bu(t)

or the continuous time analog

.. math::
    \dot{x}(t) = Ax(t) + Bu(t).

These functions provide assistive functions to deal with both the transition from network to linear systems as well as generating specific realizations of and LTI system (with weights) from and structured network.

.. autofunction:: netcontrolz.kalman_generic_rank_(G,controls[,repeats=None])

.. autofunction:: netcontrolz.controls_idx(G,controls)

.. autofunction:: netcontrolz.input_matrix_(n,controls[,m=None])

.. autofunction:: netcontrolz.matrix_realization(M,a,b)

.. autofunction:: netcontrolz.fix_diagonals(M,d)
