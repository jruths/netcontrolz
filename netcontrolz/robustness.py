"""
This module provides a basic interface to compute control-based and
reachability-based robustness metrics for controllability of networks
under failure and attack.
"""
from cacti import build_cacti, build_cacti_fixed_controls_, build_cacti_free_controls
from lti import controls_idx
from zen import DiGraph
from zen.exceptions import type_check
from numpy import floor,zeros, linspace, argsort
from numpy.random import choice

__all__ = ['edge_percolation','edge_percolation_','ATTACK_RAND','EDGE_ATTACK_INOUT_DEG','EDGE_ATTACK_OUTIN_DEG','EDGE_ATTACK_ININ_DEG','EDGE_ATTACK_OUTOUT_DEG','EDGE_ATTACK_TOTAL_DEG','CONTROL_ROBUSTNESS','REACHABILITY_ROBUSTNESS_FIXED','REACHABILITY_ROBUSTNESS_FREE']

ATTACK_RAND = 'rand'
EDGE_ATTACK_INOUT_DEG = 'inout'
EDGE_ATTACK_OUTIN_DEG = 'outin'
EDGE_ATTACK_ININ_DEG = 'inin'
EDGE_ATTACK_OUTOUT_DEG = 'outout'
EDGE_ATTACK_TOTAL_DEG = 'total'

DEGREE_EDGE_ATTACKS = [EDGE_ATTACK_INOUT_DEG,
                        EDGE_ATTACK_OUTIN_DEG,
                        EDGE_ATTACK_ININ_DEG,
                        EDGE_ATTACK_OUTOUT_DEG,
                        EDGE_ATTACK_TOTAL_DEG]

CONTROL_ROBUSTNESS = 'control'
REACHABILITY_ROBUSTNESS_FIXED = 'reach_fixed'
REACHABILITY_ROBUSTNESS_FREE = 'reach_free'
METRICS_ALL = (CONTROL_ROBUSTNESS,REACHABILITY_ROBUSTNESS_FIXED,REACHABILITY_ROBUSTNESS_FREE)

def edge_percolation(G,attack,**kwargs):
    """
    Computes robustness metrics of controllability for the directed network ``G`` according to the edge
    selection method indicated by ``attack``.

    See method ``edge_percolation_``. The only difference is that the controls here are given as a list of tuples of node objects instead of node indices.
    """
    controls = kwargs.pop('controls',None)
    if controls:
        kwargs['controls'] = controls_idx(G,controls)

    return edge_percolation_(G,attack,kwargs)

def edge_percolation_(G,attack,**kwargs):
    """
    Computes robustness metrics of controllability for the directed network ``G`` according to the edge
    selection method indicated by ``attack``.

    **Args**
        * ``attack`` indicates the type of edge selection method to use. Supported values of ``attacks`` include:
            * ``netcontrolz.EDGE_ATTACK_ININ_DEG`` (degree-based attack)
            * ``netcontrolz.EDGE_ATTACK_INOUT_DEG`` (degree-based attack)
            * ``netcontrolz.EDGE_ATTACK_OUTIN_DEG`` (degree-based attack)
            * ``netcontrolz.EDGE_ATTACK_OUTOUT_DEG`` (degree-based attack)
            * ``netcontrolz.EDGE_ATTACK_TOTAL_DEG`` (degree-based attack)
            * ``netcontrolz.ATTACK_RAND`` (random failure)
            * callable function which takes in a :py:class:`zen.DiGraph` and returns the edge index that should be removed next

            For degree-based attacks, for example, INOUT ranks edges according to their source node's IN-degree and their target node's OUT-degree; OUTIN ranks edges according to their source node's OUT-degree and their target node's IN-degree; etc.

    **KwArgs**:

        * ``controls[=None]`` (``LIST_OF_TUPLES``)

            * ``LIST_OF_TUPLES``: Representing control nodes (indices) that are attached to the nodes in G e.g. [(1,),(3,)] represents two controls that are attached to node indices 1 and 3. When controls is not given (or None), a control set with minimal number of controls will be calculated and used.

        * ``frac_rm_edges [=0.5]`` (``float``). The fraction of edges to remove from the network. If ``num_steps`` does not divide this number of edges evenly, the actual fraction of edges removed may be slightly smaller than ``frac_rm_edges``.
        * ``num_steps [=10]`` (``int``). The number of steps of percolation.
        * ``metrics`` (``py:tuple``). The metrics to calculate and return. The options are:
            * ``netcontrolz.CONTROL_ROBUSTNESS`` which measures the increase in the minimum number of controls required to control the network.
            * ``netcontrolz.REACHABILITY_ROBUSTNESS_FIXED`` which measures the decrease in the number of nodes controllable by a fixed set of controls.
            * ``netcontrolz.REACHABILITY_ROBUSTNESS_FREE`` which measures the decrease in the number of nodes controllable by a fixed number of controls.


    **Returns**:
        * ``frac_edges_removed``. The fraction of edges remaining after each percolation steps, length ``num_steps+1``.
        * ``metrics_result``. A dictionary containing the chosen ``metrics`` at the percolation steps, each item in the dictionary has length ``num_steps+1``.

    """
    type_check(G,DiGraph,'only directed graphs are supported')
    if not G.is_compact():
        raise ValueError, "the original graph must be compact (no holes in the node/edge index arrays)."

    controls = kwargs.pop('controls',None)
    frac_rm_edges = kwargs.pop('frac_rm_edges',0.5)
    if frac_rm_edges > 1.0 or frac_rm_edges <= 0:
        raise ValueError, "the fraction of edges to remove must fall in the interval (0,1]."
    num_steps = kwargs.pop('num_steps',10)
    type_check(num_steps,int)
    metrics = kwargs.pop('metrics',METRICS_ALL)

    # make sure we don't destroy the original graph
    G = G.copy()

    # if no controls are provided, we generate a minimal set
    if controls is None:
        C = build_cacti(G)
        controls = C.controls_()

    # calculate the number of edges to remove at each step, round down if it doesn't divide evenly
    L = G.num_edges
    num_step_rm_edges = int(floor(frac_rm_edges*L/num_steps))
    frac_rm_edges = float(num_step_rm_edges*num_steps)/L
    frac_edges_removed = linspace(0,frac_rm_edges,num_steps+1)

    metrics_result = {}
    if CONTROL_ROBUSTNESS in metrics:
        metrics_result[CONTROL_ROBUSTNESS] = zeros(num_steps+1)
    if REACHABILITY_ROBUSTNESS_FIXED in metrics:
        metrics_result[REACHABILITY_ROBUSTNESS_FIXED] = zeros(num_steps+1)
    if REACHABILITY_ROBUSTNESS_FREE in metrics:
        metrics_result[REACHABILITY_ROBUSTNESS_FREE] = zeros(num_steps+1)
    if len(metrics_result) == 0:
        raise ValueError, "no recognized robustness metric provided."

    if attack == ATTACK_RAND:
        rm_edges = list(choice(L,L,replace=False))
    elif attack in DEGREE_EDGE_ATTACKS:
        # aggregate the in and our degrees of each edge
        s_indeg = zeros(L)
        s_outdeg = zeros(L)
        t_indeg = zeros(L)
        t_outdeg = zeros(L)
        for eidx in G.edges_iter_():
            u,v = G.endpoints_(eidx)
            s_indeg[eidx] = G.in_degree_(u)
            s_outdeg[eidx] = G.out_degree_(u)
            t_indeg[eidx] = G.in_degree_(v)
            t_outdeg[eidx] = G.out_degree_(v)
        # establish ranking according to the attack type
        order = []
        if attack == EDGE_ATTACK_ININ_DEG:
            order = s_indeg + t_indeg
        elif attack == EDGE_ATTACK_INOUT_DEG:
            order = s_indeg + t_outdeg
        elif attack == EDGE_ATTACK_OUTIN_DEG:
            order = s_outdeg + t_indeg
        elif attack == EDGE_ATTACK_OUTOUT_DEG:
            order = s_outdeg + t_outdeg
        elif attack == EDGE_ATTACK_TOTAL_DEG:
            order = s_indeg + s_outdeg + t_indeg + t_outdeg
        rm_edges = list(argsort(order))
    elif callable(attack):
        rm_edges = [attack(G)]
    else:
        raise ValueError, "unrecognized attack: %s." % str(attack)

    # do the percolation
    for step_i in range(num_steps+1):
        # compute robustness metrics
        if CONTROL_ROBUSTNESS in metrics:
            C = build_cacti(G)
            metrics_result[CONTROL_ROBUSTNESS][step_i] = C.num_controls()
        if REACHABILITY_ROBUSTNESS_FIXED in metrics:
            C = build_cacti_fixed_controls_(G,controls)
            metrics_result[REACHABILITY_ROBUSTNESS_FIXED][step_i] = C.num_controllable_nodes()
        if REACHABILITY_ROBUSTNESS_FREE in metrics:
            C = build_cacti_free_controls(G,len(controls))
            metrics_result[REACHABILITY_ROBUSTNESS_FREE][step_i] = C.num_controllable_nodes()

        if step_i == num_steps:
            break

        # remove edges according to attack
        for rm_i in range(num_step_rm_edges):
            # pop the first element
            rm_eidx = rm_edges.pop()
            try:
                u,v = G.endpoints_(rm_eidx)
                # remove the edge
                G.rm_edge_( rm_eidx )
            except:
                raise ValueError, "edge %i cannot be found. (step %i/%i; edge %i/%i)\n%s" % (rm_eidx,step_i,num_steps+1,rm_i,num_step_rm_edges)
            # update the ordering of the edges
            if attack in DEGREE_EDGE_ATTACKS:
                if attack == EDGE_ATTACK_ININ_DEG:
                    for eidx in list(G.in_edges_(v))+list(G.out_edges_(v)):
                        order[eidx] -= 1
                elif attack == EDGE_ATTACK_INOUT_DEG:
                    for eidx in list(G.in_edges_(u))+list(G.out_edges_(v)):
                        order[eidx] -= 1
                elif attack == EDGE_ATTACK_OUTIN_DEG:
                    for eidx in list(G.in_edges_(v))+list(G.out_edges_(u)):
                        order[eidx] -= 1
                elif attack == EDGE_ATTACK_OUTOUT_DEG:
                    for eidx in list(G.in_edges_(u))+list(G.out_edges_(u)):
                        order[eidx] -= 1
                elif attack == EDGE_ATTACK_TOTAL_DEG:
                    for eidx in list(G.in_edges_(u))+list(G.out_edges_(u))+list(G.in_edges_(v))+list(G.out_edges_(v)):
                        order[eidx] -= 1
                order[rm_eidx] = -1 # put the removed edge at the end of the list
                rm_edges = list(argsort(order))
            elif callable(attack):
                rm_edges = [attack(G)]

    return frac_edges_removed, metrics_result
