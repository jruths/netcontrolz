Metrics for Quantifying Controllability
===========
Much of network control assumes the state dynamics follow a linear time-invariant model. These functions provide assistive functions to deal with both the transition from network to linear systems as well as generating specific realizations of and LTI system (with weights) from and structured network.

.. autofunction:: netcontrolz.finite_horizon_gramian(A,B,T)

.. autofunction:: netcontrolz.finite_horizon_discrete_time_gramian(A,B,T)

.. autofunction:: netcontrolz.infinite_horizon_gramian(A,B)
