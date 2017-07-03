(Linear) Structural Control
===========
Structural control is one of the major tools used to study control (and
specifically controllability) of directed networks.  Structural control is
a relaxation of conventional control methods in which the (edge) weights of the
interactions between nodes are omitted and properties are based only on the
structure (topology) of the interconnections. Most work with structural control
assumes a linear dynamics model of the underlying system (linear differential or
difference equation describing the state).

netcontrolz provides tools to build the :py:class:`Cacti`, which contains most
information about the structural control properties of a network; tools to compute
the control profile of a network, a statistic used to classify functional types
of networks; and tools to compute the robustness of the network control properties
of a graph under edge failure and attack.

To access the most basic structural control properties of a graph without
building a :py:class:`Cacti`, you may use the following complementary functions.

.. autofunction:: netcontrolz.num_dilations(G)

.. autofunction:: netcontrolz.generic_rank(G)


Constructing the Cacti of a Directed Graph
-------------------------
netcontrolz provides several handles to create the :py:class:`Cacti` for studying structural
control. A cactus is a collection of a single stem and any number of buds, which are cycles that have a distinguished edge emanating from the stem or an existing bud. The follow constructor functions provide ways to build the cacti that cover a directed graph with no controls provided, precise controls provided, or only a requirement on the number of controls used. Alternatively, you may provide a matching of your own choice and construct the :py:class:`Cacti` from it.

.. autofunction:: netcontrolz.build_cacti(G)

.. autofunction:: netcontrolz.build_cacti_free_controls(G, num_ctls)

.. autofunction:: netcontrolz.build_cacti_fixed_controls(G, fixed_ctls[, with_cycles=False])

.. autofunction:: netcontrolz.build_cacti_fixed_controls_(G, fixed_ctls[, with_cycles=False])

.. autofunction:: netcontrolz.build_cacti_from_matching(G, fixed_ctls, matching[, roots=None])

.. autofunction:: netcontrolz.build_cacti_from_matching_(G, fixed_ctls, matching[, roots=None])


Accessing the Cacti of a Directed Graph
------------------
Once the :py:class:`Cacti` are formed, questions regarding the number of controls or
the number of nodes that are controllable can be answered. The :py:class:`Stem` and :py:class:`Cycle` classes provide handles to traverse the cacti.

.. autoclass:: netcontrolz.cacti.Cacti()
    :members:

.. autoclass:: netcontrolz.cacti.Stem()
    :members:

.. autoclass:: netcontrolz.cacti.Cycle()
    :members:


Computing & Plotting the Control Profile
---------------------------------------
Control profiles were devised as a measure for quantifying the structures responsible for dictating how many
controls a network requires and where these controls much attach to the network.

.. seealso::

    J. Ruths and D. Ruths (2014). Control Profiles of Complex Networks. Science, 343(6177), 1373-1376.

.. autofunction:: netcontrolz.profile(G[,normalized=True])

Visualizing control profile plots can be particularly helpful when comparing the control profiles of different networks.  Two functions are provided for this purpose.

.. autofunction:: netcontrolz.profile_plot(G,...)

.. autofunction:: netcontrolz.profile_heatmap(G,...)

.. autofunction:: netcontrolz.profile_heatmap_weighted(G,...)

Robustness of Network Controllability
------------------------------------
Network controllability properties change as the graph is altered. An important
class of network alterations includes random failures and attacks.

.. autofunction:: netcontrolz.edge_percolation_(G,attack,...)
