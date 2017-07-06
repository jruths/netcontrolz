Metrics for Quantifying Controllability
==========================================

While controllability offers a binary quantification of control (are all state tranformations possible?), there are more fine-grained metrics that can be used. The amount of energy used to make a state transformation is contained in the Gramian matrix, and gives an assessment of how practical a control configuration might be, even if it is technically controllable.

.. autofunction:: netcontrolz.finite_horizon_gramian(A,B,T)

.. autofunction:: netcontrolz.finite_horizon_discrete_time_gramian(A,B,T)

.. autofunction:: netcontrolz.infinite_horizon_gramian(A,B)
