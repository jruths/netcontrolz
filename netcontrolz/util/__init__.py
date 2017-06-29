"""
The util submodule provides subroutines used in the control package.

matching
------------

Under structural controllability, a set of directly controlled nodes can control a limited set of other nodes.  The
number of such nodes that can be controlled is called the *reachability* or *generic rank* of the network under those
controls.

.. autofunction:: maximum_weight_matching_(G,...)

"""

from matching import *
