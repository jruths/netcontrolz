"""
This module provides functionality for building, querying, and manipulating the cacti structure of a directed network.

For information cacti structures, see the following reference:
Commault, Dion, Van der Woude. Characterization of generic properties of linear structured systems for efficient computations. Kybernetika, 38(5):503-520, 2002.
"""
from sys import setrecursionlimit
from zen import DiGraph
from zen.algorithms.matching import maximum_matching_
from zen.exceptions import type_check
from netcontrolz.util import maximum_weight_matching_
from collections import OrderedDict
#from numpy import *
from exceptions import *
import itertools as it
setrecursionlimit(2**27)

__all__ = ['build_cacti','build_cacti_fixed_controls_','build_cacti_fixed_controls','build_cacti_free_controls','build_cacti_from_matching','build_cacti_from_matching_']

class Stem:
    """
    An instance of :py:class:`Stem` represents a single stem (and any associated buds) in a cacti structure.
    """

    def __init__(self,Gi,node_seq):
        """
            Creates a new :py:class:`Stem` from a sequence of nodes given by ``node_seq``.
        """
        self._node_dict = OrderedDict.fromkeys(node_seq)
        self._G = Gi

    def __contains__(self,x):
        """
        Return ``True`` if ``x`` is a node in the stem (containment doesn't extend to buds).
        """
        return x in self._node_dict

    def __len__(self):
        """
        Number of nodes in stem
        """
        return len(self._node_dict)

    def __iter__(self):
        """
        Iterate over the nodes in the stem itself.
        """
        return self._node_dict.iterkeys()

    def origin(self):
        """
        Return the origin or starting node (object) of the stem
        """
        return  self._G.node_object(self.origin_())

    def origin_(self):
        """
        Return the origin or starting node (index) of the stem
        """
        return next(self._node_dict.iterkeys(), None)

    def terminus(self):
        """
        Returns the terminus or end node (object) of the stem.
        """
        return self._G.node_object(self.terminus_())

    def terminus_(self):
        """
        Returns the terminus or end node (index) of the stem.
        """
        x = None
        for x in self._node_dict.iterkeys():
            pass
        return x

    def _extend(self, nodes):
        """
        Extend stem by adding nodes at the end. Needed when a bud is attached
        at terminus of the stem, in which case the bud is broken down and stem
        is extended.
        """
        for u in nodes:
            self._node_dict[u] = None

    def _set_bud(self, x, bud):
        """
        Private method to add bud to a node x in this stem.
        """
        self._node_dict[x] = bud

    def length_with_buds(self):
        """
         Returns length of stem including the buds attached to it.
        """
        l = len(self._node_dict)
        for bud in self._node_dict.itervalues():
            if bud is not None:
                l += len(bud)
        return l


    def buds(self):
        """
        Returns list of buds (Cycle objects) that are attached to this stem
        """
        return list(self.buds_iter())

    def buds_iter(self):
        """
        Returns iterator over buds (Cycle objects) that are attached to this stem
        """
        return (x for x in self._node_dict.itervalues() if x is not None)

class Cycle:
    """
    :py:class:Cycle can either be a bud that is connected to a stem by a
    distinguished node or an independent cycle (no associated stem).
    """
    def __init__(self, Gi, node_seq):
        """
            Creates a new :py:class:`Cycle` from a sequence of nodes given by
            ``node_seq``.
        """
        self._node_dict = OrderedDict.fromkeys(node_seq)
        self._G = Gi
        self._stem = None
        self._dnode = None

    def __contains__(self,x):
        """
        Return ``True`` if ``x`` is a node in the cycle
        """
        return x in self._node_dict

    def __len__(self):
        """
        Number of nodes in cycle
        """
        return len(self._node_dict)

    def __iter__(self):
        """
        Iterate over nodes in the cycle
        """
        return self._node_dict.iterkeys()

    def origin(self):
        """
        A starting node (object) of this cycle,  in case of bud this is the node where
        distinguished edge is pointing to. In case of non-bud cycle, returns None
        """
        res = self.origin_()
        return None if res == -1 else self._G.node_object(res)

    def origin_(self):
        """
        A starting node (index) of this cycle,  in case of bud this is the node where
        distinguished edge is pointing to. In case of non-bud cycle, returns -1
        """
        if not self.is_bud():
            return  -1

        return next(self._node_dict.iterkeys(), -1)

    def _reorder_origin(self, x):
        """
        Private method to make target node of the distinguished edge as the
        origin of this cycle when this cycle is a bud.
        """
        if x not in self._node_dict.iterkeys():
            return
        toremove = list(it.takewhile(lambda a: a!=x, self._node_dict.iterkeys()))
        for k in toremove:
            t = self._node_dict[k]
            del self._node_dict[k]
            self._node_dict[k] = t

    def _set_stem(self,stem,dnode):
        """
        Private method to set stem and distinguished node of this cycle when
        this cycle is a bud.
        """
        self._stem = stem
        self._dnode = dnode

    def stem(self):
        """
        Return stem to which this cycle is attached if it is a bud, otherwise
        return None
        """
        if not self.is_bud():
            return None

        return self._stem

    def dist_node(self):
        """
        Return distinguished node (object) of the stem to which this cycle is attached
        if it is a bud, otherwise return None
        """
        res = self.dist_node_()
        return None if res == -1 else self._G.node_object(res)

    def dist_node_(self):
        """
        Return distinguished node (index) of the stem to which this cycle is attached
        if it is a bud, otherwise return -1
        """
        if not self.is_bud():
            return -1

        return self._dnode

    def dist_edge(self):
        """
        Return distinguished edge (u,v) (node objects) where u is in the stem and v is in this cycle/bud of the stem.
        Returns None if this cycle is not a bud.
        """
        res = self.dist_edge_()
        return None if res == -1 else self._G.endpoints(res)

    def dist_edge_(self):
        """
        Return distinguished edge e (edge index) where src(e) is in the stem and tgt(e) is in this cycle/bud of the stem.
        Returns -1 if this cycle is not a bud.
        """
        if not self.is_bud():
            return -1

        return self._G.edge_idx_(self._dnode, self.origin_())

    def is_bud(self):
        """
        Returns ``True`` if this cycle is a bud
        """
        return self._stem != None

# helper function to construct stems and cycles from a given matching
def _stems_cycs_from_matching(Gi, matching, roots=None, origins=[]):
    outmap = {a:b for a,b in matching}
    vis = set()
    indegmap = {}
    for a,b in matching:
        indegmap.setdefault(a,0)
        indegmap[b] = indegmap.get(b, 0) + 1

    # if no roots are provided, pick all the unmatched nodes
    # isolated unmatched nodes would not be captured this way
    if roots is None:
        roots = [u for u in indegmap if indegmap[u] == 0]

    def recur(z,cur,stems,cycs):
    	# already visited this node
        if z in vis:
            cycs.append(Cycle(Gi, cur))
            return

        vis.add(z)
        cur.append(z)

        # reached the end of the stem
        if z not in outmap:
            stems.append(Stem(Gi, cur))
            return

        recur(outmap[z],cur,stems,cycs)

    stems, cycs = [],[]

    # starting at the roots, build stems & cycles
    for r in roots:
        if r not in vis:
            recur(r,[],stems,cycs)

    origins = set(origins)
    # find everything that wasn't reachable from the original roots
    roots = [x for x in outmap.iterkeys() if x not in vis]
    # give preference to any nodes who directly receive a control
    roots.sort(key=lambda x: 0 if x in origins else 1)
    # starting at the non-reachable nodes, sorted with controlled nodes first, build stems & cycles
    # these origins would include nodes controlled in cycles that are not buds
    for r in roots:
        if r not in vis:
           recur(r,[],stems,cycs)

    return stems, cycs

def build_cacti_from_matching(G, fixed_ctls, matching, roots=None):
    """
    This method constructs :py:class:`Cacti` for a given directed graph ``G``
    and a precomputed matching. This requires the matching, the roots of the matching
    and the subset of roots that are controlled (the fixed_ctrls).

    **Args**: see method build_cacti_from_matching_
              Only difference here is that fixed_ctls are node objs and not
              indices. The matching is passed as tuples of node pairs.

    Returns a py:class:`Cacti` object
    """
    type_check(G,DiGraph,'only directed graphs are supported')

    if fixed_ctls is None:
        raise ValueError, "fixed_ctls cannot be None."
    if matching is None:
        raise ValueError, "matching cannot be None."

    fixed_ctls_obj = [tuple(G.node_idx(a) for a in c) for c in fixed_ctls]
    matching_obj = [G.endpoints(eidx) for eidx in matching]
    roots_obj = None
    if roots:
        roots_obj = [G.node_idx(a) for a in roots]
    return build_cacti_from_matching_(G, fixed_ctls_obj, matching_obj, roots_obj)

def build_cacti_from_matching_(G, fixed_ctls, matching, roots=None):
    """
    This method constructs :py:class:`Cacti` for a given directed graph ``G``
    and a precomputed matching. This requires the matching, the roots of the matching
    and the subset of roots that are controlled (the fixed_ctrls).

    **Args**:

        *``fixed_ctls`` (``LIST_OF_TUPLES``)
            *``LIST_OF_TUPLES``: Represents control nodes that are
                attached to the nodes in G. e.g. [(1,),(3,)] represents two controls
                that are attached to node indices 1 and 3 in G. These indices should be
                a subset of the roots of the matching or part of a cycle.
        * ``matching'' (:py:class:`list')
            *:py:class:`list': The list of edge indices that indicate the edges belonging
                to a maximum-matching for the graph.
        * ``roots'' (:py:class:`list')
            *:py:class:`list': The roots of the matching (the indices of the nodes at the
                base of the stems of the matching); if none, the roots are identified automatically

    **Raises**:
        ``ValueError``: if fixed_ctls or matching is None

    Returns a py:class:`Cacti` object
    """
    type_check(G,DiGraph,'only directed graphs are supported')

    if fixed_ctls is None:
        raise ValueError, "fixed_ctls cannot be None."
    if matching is None:
        raise ValueError, "matching cannot be None."

    cact = Cacti(G)

    cact._from_matching_case(fixed_ctls, matching, roots)

    cact._build_cacti_from_stemscycles()
    return cact

def build_cacti_fixed_controls(G, fixed_ctls, **kwargs):
    """
    This method constructs :py:class:`Cacti` for a given directed graph ``G``
    and set of controls that are fixed to some nodes in G. Maximum perfect
    weighted matching is used to calculate reachable/controllable nodes given
    the fixed controls.

    **Args**: see method build_cacti_fixed_controls_
              Only difference here is that fixed_ctls are node objs and not
              indices.

    Returns a py:class:`Cacti` object
    """
    type_check(G,DiGraph,'only directed graphs are supported')

    if fixed_ctls is None:
        raise ValueError, "fixed_ctls cannot be None."

    ctls = [tuple(G.node_idx(a) for a in c) for c in fixed_ctls]
    return build_cacti_fixed_controls_(G, ctls, **kwargs)

def build_cacti_fixed_controls_(G, fixed_ctls, **kwargs):
    """
    This method constructs :py:class:`Cacti` for a given directed graph ``G``
    and set of controls that are fixed to some nodes in G. Maximum perfect
    weighted matching is used to calculate reachable/controllable nodes given
    the fixed controls.

    **Args**:

        *``fixed_ctls`` (``LIST_OF_TUPLES``)
            *``LIST_OF_TUPLES``: Represents control nodes that are
                attached to the nodes in G. e.g. [(1,),(3,)] represents two controls
                that are attached to node indices 1 and 3 in G.

    **KwArgs**:
        *``randomize[=False]`` (``Boolean``). Indicates whether the matching
            should be randomized
        *``with_cycles[=False]`` (``Boolean``). Indicates whether
            cycles not reachable from the ``fixed_ctls`` should be
            included in the matching/cacti

    **Raises**:
        ``ValueError``: if fixed_ctls is None

    Returns a py:class:`Cacti` object
    """
    type_check(G,DiGraph,'only directed graphs are supported')

    if fixed_ctls is None:
        raise ValueError, "fixed_ctls cannot be None."

    cact = Cacti(G)

    cact._fixed_control_case(fixed_ctls, **kwargs)

    cact._build_cacti_from_stemscycles()
    return cact


def build_cacti_free_controls(G, num_ctls, **kwargs):
    """
    This method constructs :py:class:`Cacti` for a given directed graph ``G``
    and a fixed number of controls that can be assigned to any nodes in G. Maximum perfect
    weighted matching is used to calculate reachable/controllable nodes given
    the number of controls.

    **Args**:

        * ``num_ctls`` (``int``): Number of control nodes that are
          to be attached to the nodes in G. If num_ctls is 0 a degenerate
          matching of only cycles will be found.

    **KwArgs**:
        * ``randomize[=False]`` (``Boolean``). Indicates whether the matching
          should be randomized

    Returns a py:class:`Cacti` object
    """
    type_check(G,DiGraph,'only directed graphs are supported')

    if num_ctls < 1:
        raise ValueError, "num_ctls cannot be less than one."

    cact = Cacti(G)

    cact._free_control_case(num_ctls, **kwargs)

    cact._build_cacti_from_stemscycles()
    return cact


def build_cacti(G):
    """
    This method constructs :py:class:`Cacti` for a given directed graph ``G``
    using maximum unweighted matching. The resulting matching forms a cacti
    which also gives the minimum number and locations of controls required for the full control.

    Returns a py:class:`Cacti` object
    """
    type_check(G,DiGraph,'only directed graphs are supported')

    cact = Cacti(G)

    cact._min_controls_case()

    cact._build_cacti_from_stemscycles()
    return cact

class Cacti:
    """
    The Cacti class represents the cacti  control structure of a given network. It can be used to find the location (and number) of inputs required to control the network as
    well as to access the underlying cacti structure.  It can also be used to find the extent to which a network can be controlled
    by a specified set of controls and the related constrained cacti.
    """

    def __init__(self,G):
        """
        Private constructor. Should not be used directly. build_cacti methods should be used, which return Cacti object.
        """
        self._G = G

        self._stems = []
        self._cycles = []
        self._matching = []
        self._controls = []
        self._controllable_nodes=set()

    def num_controls(self):
        """
        Return number of controls inputs nodes.
        """
        return len(self._controls)

    def stems(self):
        return list(self._stems)

    def cycles(self):
        return list(self._cycles)

    def non_bud_cycles(self):
        """
        Return the non-bud (independent cycles)
        """
        return filter(lambda x: not x.is_bud(), self._cycles)


    # helper function to build cacti from stems and cycles
    def _build_cacti_from_stemscycles(self):
        self._stems.sort(key=lambda x:len(x), reverse= True)
        self._cycles.sort(key=lambda x:len(x), reverse= True)

        # for every node in each stem, check if any of the nbrs are part of a cycle (that isn't already a bud)
        for stem in self._stems:
            for u in stem:
                for v in self._G.out_neighbors_(u):
                    for icyc,cyc in enumerate(self._cycles):
                        if cyc and not cyc.is_bud() and v in cyc:
                            cyc._reorder_origin(v)
                            if u == stem.terminus_():
                                stem._extend(cyc)
                                self._cycles[icyc] = None
                            else:
                                stem._set_bud(u, cyc)
                                cyc._set_stem(stem, u)

        # TODO: need to enable bud of bud...???
        # loop over all the nodes in the bud cycles to see if their out nbrs are part of a cycle (that isn't already a bud)
        # new_bud_cycles = [c for c in self._cycles if c.is_bud()]
        # while len(new_bud_cycles) > 0:
        # 	bud_cycles = new_bud_cycles
        # 	new_bud_cycles = []
        # 	for cycle in bud_cycles:
        # 		for u in cycle:
        # 			for v in self._G.out_neighbors_(u):
        # 				for icyc,cyc in enumerate(self._cycles):
        #                 if cyc and not cyc.is_bud() and v in cyc:
        #                     cyc._reorder_origin(v)
        #                     # ???NEED TO FIX STEM cycle._set_bud(u, cyc)
        #                     # ???NEED TO FIX CYCLE cyc._set_stem(cycle, u)
        #                     new_bud_cycles.append(cyc)

        # the remaining non-bud cycles
        self._cycles = [c for c in self._cycles if c is not None and len(c) > 0]

        # the location of the controls are the base of the stems
        controls = [ [stem.origin_()] for stem in self._stems ]

        # assign the remaining, non-bud cycles to the stems one at a time, looping through the stems as many times as needed
        # if there are no stems, then attach all cycles to a single control
        if controls:
            idx = 0
            for cyc in self._cycles:
                if not cyc.is_bud():
                    controls[idx].append(list(cyc)[0])
                    idx = (idx+1) % len(controls)
            controls = [tuple(a) for a in controls]
        else:
            controls = [ tuple(list(cyc)[0] for cyc in self._cycles) ]

        self._controls = controls

        ctlable_nodes = set()
        ctlable_nodes.update([a for s in self._stems for a in s])
        ctlable_nodes.update([a for c in self._cycles for a in c])
        self._controllable_nodes = ctlable_nodes


    # To find minimum controls required to fully control a network using
    # Unweighted maximum Matching
    def _min_controls_case(self):
        #TODO Randomize matching
        G = self._G
        matching = set(maximum_matching_(G))
        matching = [G.endpoints_(eidx) for eidx in matching]
        stems, cycles = _stems_cycs_from_matching(G, matching)
        # in the min controls case, any remaining unmatched nodes (isolated nodes), we need to add a stem and control
        for u in set(G.nodes_iter_()).difference(set(x for y in matching for x in y)) :
            stems.append(Stem(G, [u]))
            # ??? Do we need to add these nodes to the self._controllable_nodes set?
        self._matching, self._stems, self._cycles = matching, stems, cycles


    # To find cacti and hence controllable nodes given a fixed set of controls
    # (fixed_ctls) using maximum weighted matching
    def _fixed_control_case(self, fixed_ctls, **kwargs):
        kwargs['controls'] = fixed_ctls
        G = self._G
        origins = [a for b in fixed_ctls for a in b]
        matching,roots = maximum_weight_matching_(G, **kwargs)[1:3]
        matching = [G.endpoints_(eidx) for eidx in matching]
        self._stems, self._cycles = _stems_cycs_from_matching(G, matching, roots, origins)
        self._matching = matching


    # To find cacti and hence controllable nodes given a fixed number of
    # controls that can be assigned anywhere for maximum control possible
    def _free_control_case(self, num_ctls, **kwargs):
        G = self._G
        ctls = [ tuple(G.nodes_()) for x in xrange(num_ctls) ]
        kwargs['controls'] = ctls
        matching,roots = maximum_weight_matching_(G, **kwargs)[1:3]
        matching = [ G.endpoints_(eidx) for eidx in matching ]
        self._stems, self._cycles = _stems_cycs_from_matching(G, matching, roots=roots)
        self._matching = matching

    # Builds the cacti directly from the matching
    # fixed_ctls contains a subset of the roots, since some roots may not be attached to a control
    def _from_matching_case(self, fixed_ctls, matching, roots):
        origins = [a for b in fixed_ctls for a in b]
        matching = [self._G.endpoints_(eidx) for eidx in matching]
        self._stems, self._cycles = _stems_cycs_from_matching(self._G, matching, roots, origins)
        self._matching = matching

    def controls(self):
        """
        Returns list of nodes (objects) where controls should be attached.
        In case of unweighted matching, it is the minimum number
        of controls required for full control of the network. In case of
        weighted matching (when fixed set of controls are given), this method
        returns the controls that are sufficient for controlling the maximum
        possible nodes of the network.
        """
        ctls = self.controls_()
        return [tuple(self._G.node_object(x) for x in t) for t in ctls]

    def controls_(self):
        """
        Returns list of nodes (indices) where controls should be attached.
        In case of unweighted matching, it is the minimum number
        of controls required for full control of the network. In case of
        weighted matching (when fixed set of controls are given), this method
        returns the controls that are sufficient for controlling the maximum
        possible nodes of the network.
        """
        return list(self._controls)

    def controllable_nodes(self):
        """
        Returns set of nodes objects that are controllable
        """
        return set([self._G.node_object(a) for a in self._controllable_nodes])

    def controllable_nodes_(self):
        """
        Returns set of nodes indices that are controllable
        """
        return set(self._controllable_nodes)

    def num_controllable_nodes(self):
        """
        Return number of nodes that are controllable
        """
        return len(self._controllable_nodes)

    def matching(self):
        """
        Returns a  matching in form of a list of edges (u,v) where u and v are node objects, as calculated by the maximum matching algorithm (unweighted or weighted as the case maybe)
        """
        return [self._G.endpoints(e) for e in self.matching_()]

    def matching_(self):
        """
        Returns a  matching, a list of edge indices as calculated by the maximum matching algorithm (unweighted or weighted as the case maybe)
        """
        return [self._G.edge_idx_(u,v) for u,v in self._matching]
