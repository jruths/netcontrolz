"""
The ``zen.algorithms.matching`` module provides routines for computing `maximum-matchings <???>`_ on various types of graphs.

.. autofunction:: maximum_matching

.. autofunction:: maximum_matching_

.. autofunction:: hopcroft_karp_

Deven Parekh: added __max_weight_matching
    TODO: write max_weight_matching_ which takes input as  BipartiteGraph

"""

from zen.bipartite cimport BipartiteGraph
from zen.digraph cimport DiGraph
from zen.exceptions import *
import numpy as np, time, os
from numpy.random import randint
cimport numpy as np
from Queue import Queue
from libcpp.vector cimport vector
from libc.stdlib cimport malloc, free
from cython.operator cimport dereference as dref
from cython.operator cimport preincrement as incr


__all__ = ['maximum_weight_matching_']

# TODO(druths): Add support for matching undirected networks

cdef extern from "max_weight_matching.hpp":
    cdef struct edge_t:
        int src, tgt
        double w
        edge_t(int, int, double)
    cdef void mwb_matching(vector[vector[int]] &G, vector[edge_t] &E, \
        vector[int] &U, vector[int] &V, vector[int] &M)


#TODO write wrapper around mwb_matching for BipartiteGraph
#def max_weight_matching_(G):
#

def maximum_weight_matching_(Gi, **kwargs):
    """
        Find maximum weighted matching for given Directed Graph ``G``

        **KwArgs**:
            *``controls[=None]`` (LIST_OF_TUPLES)
                *``LIST_OF_TUPLES``: Representing control nodes that are
                attached to the nodes in G e.g. [(1,),(3,)] represents two controls
                that are attached to node indices 1 and 3.
                When controls is not given (or None), the result will consist
                of matching such that all cycles will be found (because cycles don't
                need any controls.)
            *``randomize[=False]`` (``Boolean``). Indicates whether the matching
                should be randomized
            *``with_cycles[=False]`` (``Boolean``). Indicates whether
                independent cycles not reachable from the ``controls`` should be
                included in the matching

        **Returns**:
            ``(n,M,R)``: where ``n`` (int) is number of matched nodes.
                        ``M`` is a list of edge indices representing the matching
                        ``R`` is a list of node indices representing origins of
                        the stems.
    """
    cdef:
        vector[int] *U
        vector[int] *V
        vector[edge_t] *E
        vector[vector[int]] *G
        vector[int].iterator it
        vector[edge_t].iterator ite
        unsigned int u,v,N,t,x,y,r,s,e
        double w,MAXW
        edge_t temp_edge
        vector[int] *M

    U = new vector[int]()
    V = new vector[int]()
    M = new vector[int]()
    E = new vector[edge_t]()

    controls = kwargs.pop('controls', None)
    doshuffle = kwargs.pop('randomize',False)
    indcyc = kwargs.pop('with_cycles',False)

    if controls is None:
        controls = []
    node2uv, uv2node = {}, {}

    if doshuffle:
        np.random.seed(int(time.time()+os.getpid()*np.random.random()*1000))

    #remove non reachable nodes
    if controls:
        vis = set()

        def dfs(r):
            vis.add(r)
            for v in Gi.out_neighbors_(r):
                if v not in vis:
                    dfs(v)

        for ctl in controls:
            for driver in ctl:
                if driver not in  vis:
                    dfs(driver)
    else:
        vis = set(Gi.nodes_())

    # first transform the graph itself
    N = 0
    #for t in vis:
    for t in Gi.nodes_iter_():
        if indcyc or t in vis:
            u, v = 2*N, 2*N+1
            node2uv[t] = N
            uv2node[N] = t
            N += 1
            U.push_back(u)
            V.push_back(v)
            if doshuffle:
                x, y = randint(U.size()), U.size()-1
                s = dref(U)[x]
                dref(U)[x] = dref(U)[y]
                dref(U)[y] = s
                x, y = randint(V.size()), V.size()-1
                s = dref(V)[x]
                dref(V)[x] = dref(V)[y]
                dref(V)[y] = s

    MAXW = sum(len(ctl) for ctl in controls) + Gi.size() + 100.0

    # add the edges
    for e in Gi.edges_iter_():
        x, y = Gi.endpoints_(e)
        #if x in vis and y in vis:
        if indcyc or (x in vis and y in vis):
            u,v = 2*node2uv[x], 2*node2uv[y]+1
            temp_edge.src = u
            temp_edge.tgt = v
            temp_edge.w = MAXW + 1.0
            E.push_back(temp_edge)

    # add control nodes and forward edges with weight 1
    START_CNODE = N
    for ctl in controls:
        cnode = N
        N += 1
        #if len(ctl) > 0:
            #d = ctl[0]
        for d in ctl:
            u, v = 2*cnode, 2*node2uv[d]+1
            temp_edge.src = u
            temp_edge.tgt = v
            temp_edge.w = MAXW + 1.0
            E.push_back(temp_edge)

        for x in vis:
            u, v = 2*node2uv[x], 2*cnode+1
            temp_edge.src = u
            temp_edge.tgt = v
            temp_edge.w = MAXW
            E.push_back(temp_edge)

        u, v = 2*cnode, 2*cnode+1
        U.push_back(u)
        V.push_back(v)
        if doshuffle:
            x, y = randint(U.size()), U.size()-1
            s = dref(U)[x]
            dref(U)[x] = dref(U)[y]
            dref(U)[y] = s
            x, y = randint(V.size()), V.size()-1
            s = dref(V)[x]
            dref(V)[x] = dref(V)[y]
            dref(V)[y] = s

    # add self loops with weight 0
    for s in xrange(N):
        if s >= START_CNODE or not Gi.has_edge_(uv2node[s],uv2node[s]):
            u, v = 2*s, 2*s+1
            temp_edge.src = u
            temp_edge.tgt = v
            temp_edge.w = MAXW
            E.push_back(temp_edge)

    G = new vector[vector[int]](2*N, vector[int]())

    for e in xrange(E.size()):
        u,v,w = dref(E)[e].src, dref(E)[e].tgt, dref(E)[e].w
        assert u%2==0 and v%2==1
        dref(G)[u].push_back(e)
        dref(G)[v].push_back(e)
        if doshuffle:
            x, y = randint(dref(G)[u].size()), dref(G)[u].size()-1
            s = dref(G)[u][x]
            dref(G)[u][x] = dref(G)[u][y]
            dref(G)[u][y] = s
            x, y = randint(dref(G)[v].size()), dref(G)[v].size()-1
            s = dref(G)[v][x]
            dref(G)[v][x] = dref(G)[v][y]
            dref(G)[v][y] = s

    #####
    # run the weighted bipartite matching

    mwb_matching(dref(G),dref(E),dref(U),dref(V),dref(M))

    result, roots = [], []
    num_matched = 0
    for x in xrange(M.size()):
        e = dref(M)[x]
        u,v,w = dref(E)[e].src, dref(E)[e].tgt, dref(E)[e].w
        assert u%2==0 and v%2==1
        if w > MAXW:
            num_matched += 1
            if (u/2 < START_CNODE):
                result.append(Gi.edge_idx_(uv2node[u/2],uv2node[v/2]))
            else:
                roots.append(uv2node[v/2])

    #  free memory
    del G, U, V, E, M

    return num_matched, result, roots

def directed_max_weight_matching(Gi, **kwargs):
    cdef:
        vector[int] *U
        vector[int] *V
        vector[edge_t] *E
        vector[vector[int]] *G
        vector[int].iterator it
        vector[edge_t].iterator ite
        unsigned int u,v,N,t,x,y,r,s,e
        double w
        edge_t temp_edge
        vector[int] *M

    U = new vector[int]()
    V = new vector[int]()
    M = new vector[int]()
    E = new vector[edge_t]()

    doshuffle = kwargs.pop('randomize',False)
    dologging = kwargs.pop('printlog',False)

    if dologging:
        print 'Logging is on.'

    node2uv, uv2node = {}, {}

    if doshuffle:
        np.random.seed(int(time.time()+os.getpid()*np.random.random()*1000))

    if dologging:
        print 'Transforming the graph'
    # first transform the graph itself
    N = 0
    for t in Gi.nodes_iter_():
        u, v = 2*N, 2*N+1
        node2uv[t] = N
        uv2node[N] = t
        N += 1
        U.push_back(u)
        V.push_back(v)
        if doshuffle:
            x, y = randint(U.size()), U.size()-1
            s = dref(U)[x]
            dref(U)[x] = dref(U)[y]
            dref(U)[y] = s
            x, y = randint(V.size()), V.size()-1
            s = dref(V)[x]
            dref(V)[x] = dref(V)[y]
            dref(V)[y] = s

    if dologging:
        print 'Adding original edges'
    # add the edges
    for e,w in Gi.edges_iter_(weight=True):
        x, y = Gi.endpoints_(e)
        u,v = 2*node2uv[x], 2*node2uv[y]+1
        temp_edge.src = u
        temp_edge.tgt = v
        temp_edge.w = w
        E.push_back(temp_edge)

    if dologging:
        print 'Building graph vector'
    # build the G vector: G[x] contains list of edge indices touching x
    G = new vector[vector[int]](2*N, vector[int]())
    for e in xrange(E.size()):
        u,v,w = dref(E)[e].src, dref(E)[e].tgt, dref(E)[e].w
        assert u%2==0 and v%2==1
        dref(G)[u].push_back(e)
        dref(G)[v].push_back(e)
        if doshuffle:
            x, y = randint(dref(G)[u].size()), dref(G)[u].size()-1
            s = dref(G)[u][x]
            dref(G)[u][x] = dref(G)[u][y]
            dref(G)[u][y] = s
            x, y = randint(dref(G)[v].size()), dref(G)[v].size()-1
            s = dref(G)[v][x]
            dref(G)[v][x] = dref(G)[v][y]
            dref(G)[v][y] = s

    if dologging:
        print "Starting call to external weighted maximum matching library"
    #####
    # run the weighted bipartite matching
    mwb_matching(dref(G),dref(E),dref(U),dref(V),dref(M))
    if dologging:
        print "Returned from external weighted maximum matching library"

    result = []
    num_matched = 0
    for x in xrange(M.size()):
        e = dref(M)[x]
        u,v,w = dref(E)[e].src, dref(E)[e].tgt, dref(E)[e].w
        assert u%2==0 and v%2==1
        num_matched += 1
        result.append(Gi.edge_idx_(uv2node[u/2],uv2node[v/2]))

    #  free memory
    del G, U, V, E, M

    return result

NO_COUNT_EDGE = '1234nocount'

def directed_fixed_roots_max_weight_matching(Gi, **kwargs):
    cdef:
        vector[int] *U
        vector[int] *V
        vector[edge_t] *E
        vector[vector[int]] *G
        vector[int].iterator it
        vector[edge_t].iterator ite
        unsigned int u,v,N,t,x,y,r,s,e
        double w,MAXW
        edge_t temp_edge
        vector[int] *M

    U = new vector[int]()
    V = new vector[int]()
    M = new vector[int]()
    E = new vector[edge_t]()

    roots = kwargs.pop('roots', None)
    doshuffle = kwargs.pop('randomize',False)
    indcyc = kwargs.pop('with_cycles',False)

    if roots is None:
        roots = []
    node2uv, uv2node = {}, {}

    if doshuffle:
        np.random.seed(int(time.time()+os.getpid()*np.random.random()*1000))

    #remove nodes not reachable from the roots
    if roots:
        vis = set()

        def dfs(r):
            vis.add(r)
            for v in Gi.out_neighbors_(r):
                if v not in vis:
                    dfs(v)

        for root in roots:
            for rt in root:
                if rt not in vis:
                    dfs(rt)
    else:  # if the roots are free, all nodes are potentially reachable
        vis = set(Gi.nodes_())

    # first transform the graph itself
    N = 0
    for t in Gi.nodes_iter_():
        if indcyc or t in vis:
            u, v = 2*N, 2*N+1
            node2uv[t] = N
            uv2node[N] = t
            N += 1
            U.push_back(u)
            V.push_back(v)
            if doshuffle:
                x, y = randint(U.size()), U.size()-1
                s = dref(U)[x]
                dref(U)[x] = dref(U)[y]
                dref(U)[y] = s
                x, y = randint(V.size()), V.size()-1
                s = dref(V)[x]
                dref(V)[x] = dref(V)[y]
                dref(V)[y] = s


    MAXW = ( sum(len(root) for root in roots) + sum(0 if d==NO_COUNT_EDGE else round(Gi.weight_(x)) for x,d in Gi.edges_iter_(data=True)) ) * 10.0

    # add the edges
    for e,w in Gi.edges_iter_(weight=True):
        x, y = Gi.endpoints_(e)
        if indcyc or (x in vis and y in vis):
            u,v = 2*node2uv[x], 2*node2uv[y]+1
            temp_edge.src = u
            temp_edge.tgt = v
            temp_edge.w = w + MAXW
            E.push_back(temp_edge)

    # add control nodes and forward edges with weight 1
    START_RNODE = N
    for root in roots:
        rnode = N
        N += 1
        #if len(ctl) > 0:
            #d = ctl[0]
        for d in root:
            u, v = 2*rnode, 2*node2uv[d]+1
            temp_edge.src = u
            temp_edge.tgt = v
            temp_edge.w = MAXW + 1.0
            E.push_back(temp_edge)

        for x in vis:
            u, v = 2*node2uv[x], 2*rnode+1
            temp_edge.src = u
            temp_edge.tgt = v
            temp_edge.w = MAXW
            E.push_back(temp_edge)

        u, v = 2*rnode, 2*rnode+1
        U.push_back(u)
        V.push_back(v)
        if doshuffle:
            x, y = randint(U.size()), U.size()-1
            s = dref(U)[x]
            dref(U)[x] = dref(U)[y]
            dref(U)[y] = s
            x, y = randint(V.size()), V.size()-1
            s = dref(V)[x]
            dref(V)[x] = dref(V)[y]
            dref(V)[y] = s

    # add self loops with weight 0
    for s in xrange(N):
        if s >= START_RNODE or not Gi.has_edge_(uv2node[s],uv2node[s]):
            u, v = 2*s, 2*s+1
            temp_edge.src = u
            temp_edge.tgt = v
            temp_edge.w = MAXW
            E.push_back(temp_edge)

    # build the G vector: G[x] contains list of edge indices touching x
    G = new vector[vector[int]](2*N, vector[int]())
    for e in xrange(E.size()):
        u,v,w = dref(E)[e].src, dref(E)[e].tgt, dref(E)[e].w
        assert u%2==0 and v%2==1
        dref(G)[u].push_back(e)
        dref(G)[v].push_back(e)
        if doshuffle:
            x, y = randint(dref(G)[u].size()), dref(G)[u].size()-1
            s = dref(G)[u][x]
            dref(G)[u][x] = dref(G)[u][y]
            dref(G)[u][y] = s
            x, y = randint(dref(G)[v].size()), dref(G)[v].size()-1
            s = dref(G)[v][x]
            dref(G)[v][x] = dref(G)[v][y]
            dref(G)[v][y] = s

    #####
    # run the weighted bipartite matching
    mwb_matching(dref(G),dref(E),dref(U),dref(V),dref(M))

    result, selected_roots = [], []
    num_matched = 0
    for x in xrange(M.size()):
        e = dref(M)[x]
        u,v,w = dref(E)[e].src, dref(E)[e].tgt, dref(E)[e].w
        assert u%2==0 and v%2==1
        if w > MAXW:
            num_matched += 1
            if (u/2 < START_RNODE):
                result.append(Gi.edge_idx_(uv2node[u/2],uv2node[v/2]))
            else:
                selected_roots.append(uv2node[v/2])

    #  free memory
    del G, U, V, E, M

    return num_matched, result, selected_roots

def directed_free_roots_max_weight_matching(Gi, **kwargs):

    num_roots = kwargs.pop('num_roots', 0)
    force_use_all = kwargs.pop('force_use_all', False)

    Gi2 = Gi.copy()
    extra_nodes = set(Gi2.add_nodes(num_roots))
    extra_edges = set([])
    kwargs['roots'] = [ (x,) for x in extra_nodes ]
    w = 1
    if force_use_all:
        w = sum(Gi.weight_(x) for x in Gi.edges_iter_()) * 10.0
    for u in extra_nodes:
        for v in Gi.nodes_iter_():
            eidx = Gi2.add_edge_(u,v,weight=w,data=NO_COUNT_EDGE)
            extra_edges.add(eidx)

    num_matched, result, selected_roots = directed_fixed_roots_max_weight_matching(Gi2, **kwargs)

    num_matched -= num_roots
    result = set(result)
    Gi_selected_roots = [ Gi2.endpoints_(eidx)[1] for eidx in extra_edges.intersection(result) ]
    result.difference_update(extra_edges)

    return num_matched, list(result), Gi_selected_roots
