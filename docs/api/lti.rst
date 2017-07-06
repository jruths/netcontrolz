Linear Time-Invariant (LTI) Systems
====================================

Much of network control assumes the state dynamics follow a linear time-invariant model.

.. math::
    x(t+1) = Ax(t) + Bu(t)

or the continuous time analog

.. math::
    \dot{x}(t) = Ax(t) + Bu(t).

These functions provide assistive functions to help with the transition from a network model to a linear systems model as well as facilitate generating specific realizations of a LTI system (with weights) from a structured network.

.. autofunction:: netcontrolz.kalman_generic_rank(G,controls[,repeats=None])

.. autofunction:: netcontrolz.kalman_generic_rank_(G,controls[,repeats=None])

.. autofunction:: netcontrolz.controls_idx(G,controls)

.. autofunction:: netcontrolz.input_matrix(G,controls[,m=None])

.. autofunction:: netcontrolz.input_matrix_(n,controls[,m=None])

.. autofunction:: netcontrolz.matrix_realization(M,a,b)

.. autofunction:: netcontrolz.fix_diagonals(M,d)
