"""
The ``netcontrolz`` library provides functions and classes for analyzing the control structure of a network.  In particular, much of the
functionality currently implemented focuses on the perspective and machinery provided by structural controllability. ``netcontrolz`` is
built on top of the Zero-Effort Network (zen) Python network library.

Cacti
-------------
The ``py:class:Cacti`` class provides the heart of structural controllability analysis. A cactus is a collection of a single stem
and any number of buds, which are cycles that have a distinguished edge eminating from the stem or an existing bud. The ``cacti``
methods provide ways to compute the cacti that cover a directed graph with no controls provided, precise controls provided, or
only a requirement on the number of controls used. Once the cacti are formed, questions regarding the number of controls or
the number of nodes that are controllable can be answered. The ``py:class:Stem`` and ``py:class:Cycle`` classes provide handles to traverse the
cacti.

.. autofunction:: build_cacti(G)

.. autofunction:: build_cacti_free_controls(G, num_ctls)

.. autofunction:: build_cacti_fixed_controls(G, fixed_ctls[, with_cycles=False])

.. autofunction:: build_cacti_fixed_controls_(G, fixed_ctls[, with_cycles=False])

.. autofunction:: build_cacti_from_matching(G, fixed_ctls, matching[, roots=None])

.. autofunction:: build_cacti_from_matching_(G, fixed_ctls, matching[, roots=None])

Dilations
------------
Dilations play a central role in the structural control analysis of networks. At the most fundamental level, dilations
are the network structures the necessitate adding additional controls to render the network controllable. 

.. autofunction:: num_dilations(G)

.. autofunction:: all_in_neighbors(G,S)

.. autofunciton:: all_in_neighbors_(G,S)

Linear Time-Invariant (LTI) Systems
-----------------------------------
Much of network control assumes the state dynamics follow a linear time-invariant model. These functions provide
assistive functions to deal with both the transition from network to linear systems as well as generating specific
realizations of and LTI system (with weights) from and structured network.

.. autofunction:: controls_idx(G,controls)

.. autofunction:: input_matrix_(n,controls[,m=None])

.. autofunction:: matrix_realization(M,a,b)

.. autofunction:: fix_diagonals(M,d) 

Rank
------------

Under structural controllability, a set of directly controlled nodes can control a limited set of other nodes.  The
number of such nodes that can be controlled is called the *reachability* or *generic rank* of the network under those
controls.

.. autofunction:: kalman_generic_rank(G)
.. autofunction:: generic_rank(G)

Control Profiles
----------------

Control profiles were devised as a measure for quantifying the structures responsible for dictating how many
controls a network requires and where these controls much attach to the network.

.. seealso::

    J. Ruths and D. Ruths (2014). Control Profiles of Complex Networks. Science, 343(6177), 1373-1376.

.. autofunction:: profile(G,...)

Visualizing control profile plots can be particularly helpful when comparing the control profiles of different
networks.  Two functions are provided for this purpose.

.. autofunction:: profile_plot(G,...)

.. autofunction:: profile_heatmap(G,...)

"""

from cacti import *
from dilations import *
from lti import *
from metrics import *
from pplot import *
from profile import *
from rank import *
from robustness import *
import util
