"""
The ``netcontrolz`` library provides functions and classes for analyzing the control structure of a network.  In particular, much of the
functionality currently implemented focuses on the perspective and machinery provided by structural controllability. ``netcontrolz`` is
built on top of the Zero-Effort Network (zen) Python network library.


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
