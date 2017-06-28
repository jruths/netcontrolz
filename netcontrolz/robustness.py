"""
This module provides a basic interface to compute control-based and
reachability-based robustness metrics for controllability of networks
under failure and attack.
"""
from cacti import build_cacti, build_cacti_fixed_controls_, build_cacti_free_controls
from zen import choose_edge_
from zen.exceptions import type_check
from numpy import floor,zeros, linspace, argsort

__all__ = ['edge_percolation_']

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

METRIC_CONTROL_ROBUSTNESS = 'control'
METRIC_REACHABILITY_ROBUSTNESS_FIXED = 'reach_fixed'
METRIC_REACHABILITY_ROBUSTNESS_FREE = 'reach_free'
METRICS_ALL = (METRIC_CONTROL_ROBUSTNESS,METRIC_REACHABILITY_ROBUSTNESS_FIXED,METRIC_REACHABILITY_ROBUSTNESS_FREE)

def edge_percolation_(G,attack,**kwargs):
    """
    Computes robustness metrics of controllability for the directed network ``G`` according to the edge
    selection method indicated by ``attack``.

    **Args**
        * ``attack`` (``string``). Indicates the type of edge selection method to use. Attack can also be a
            callable function which takes in a DiGraph and returns the edge index that should be removed next.

    **KwArgs**:
        *``controls[=None]`` (LIST_OF_TUPLES)
            *``LIST_OF_TUPLES``: Representing control nodes that are
            attached to the nodes in G e.g. [(1,),(3,)] represents two controls
            that are attached to node indices 1 and 3.
            When controls is not given (or None), a control set with minimal number of
            controls will be calculated and used.
        * ``frac_rm_edges [=0.5]`` (``float``). The fraction of edges to remove from the network. If
            ``num_steps`` does not divide this number of edges evenly, the actual fraction of edges removed
            may be slightly smaller than ``frac_rm_edges``.
        * ``num_steps [=10]`` (``int``). The number of steps of percolation.
        * ``metrics [=('control','reach_fixed','reach_free')]`` (``py:tuple``). The metrics to return:
            control-based robustness and/or reachability-based robustness (fixed and/or free).

    **Returns**:
        *``num_edges_removed``. The fraction of edges remaining after each percolation steps, length ``num_steps``+1.
        *``metrics_result``. A dictionary containing the chosen ``metrics`` at the percolation steps, each item
            in the dictionary has length ``num_steps``+1.
    """
    type_check(G,DiGraph,'only directed graphs are supported')

    controls = kwargs.pop('controls',None)
    frac_rm_edges = kwargs.pop('frac_rm_edges',0.5)
    type_check(frac_rm_edges,float)
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
    frac_rm_edges = float(num_step_rm_edges*10)/L
    num_edges_removed = linspace(0,frac_rm_edges,num_steps+1)

    metrics_result = {}
    if METRIC_CONTROL_ROBUSTNESS in metrics:
        metrics_result[METRIC_CONTROL_ROBUSTNESS] = zeros(num_steps+1)
    if METRIC_REACHABILITY_ROBUSTNESS_FIXED in metrics:
        metrics_result[METRIC_REACHABILITY_ROBUSTNESS_FIXED] = zeros(num_steps+1)
    if METRIC_REACHABILITY_ROBUSTNESS_FREE in metrics:
        metrics_result[METRIC_REACHABILITY_ROBUSTNESS_FREE] = zeros(num_steps+1)
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
        raise ValueError, "unrecognized attack."

    # do the percolation
    for step_i in range(num_steps+1):
        # compute robustness metrics
        if METRIC_CONTROL_ROBUSTNESS in metrics:
            C = build_cacti(G)
            metrics_result[METRIC_CONTROL_ROBUSTNESS][step_i] = C.num_controls()
        if METRIC_REACHABILITY_ROBUSTNESS_FIXED in metrics:
            C = build_cacti_fixed_controls_(G,controls)
            metrics_result[METRIC_REACHABILITY_ROBUSTNESS_FIXED][step_i] = C.num_controllable_nodes()
        if METRIC_REACHABILITY_ROBUSTNESS_FREE in metrics:
            C = build_cacti_free_controls_(G,len(controls))
            metrics_result[METRIC_REACHABILITY_ROBUSTNESS_FREE][step_i] = C.num_controllable_nodes()

        # remove edges according to attack
        for rm_i in range(num_step_rm_edges):
            # pop the first element
            rm_eidx = rm_edges.pop()
            u,v = G.endpoints_(rm_eidx)
            # update the ordering of the edges
            if attack in DEGREE_EDGE_ATTACKS:
                order = []
                if attack == EDGE_ATTACK_ININ_DEG:
                    for eidx in G.in_edges_(v)+G.out_edges_(v):
                        order[eidx] -= 1
                elif attack == EDGE_ATTACK_INOUT_DEG:
                    for eidx in G.in_edges_(u)+G.out_edges_(v):
                        order[eidx] -= 1
                elif attack == EDGE_ATTACK_OUTIN_DEG:
                    for eidx in G.in_edges_(v)+G.out_edges_(u):
                        order[eidx] -= 1
                elif attack == EDGE_ATTACK_OUTOUT_DEG:
                    for eidx in G.in_edges_(u)+G.out_edges_(u):
                        order[eidx] -= 1
                elif attack == EDGE_ATTACK_TOTAL_DEG:
                    for eidx in G.in_edges_(u)+G.out_edges_(u)+G.in_edges_(v)+G.out_edges_(v):
                        order[eidx] -= 1
                order[rm_eidx] = -1 # put the removed edge at the end of the list
                rm_edges = list(argsort(order))
            elif callable(attack):
                rm_edges = [attack(G)]

            # finally remove the edge
            G.rm_edge_( rm_eidx )

    return num_edges_removed, metrics_result

# from numpy import arange,trapz
# def mip_(G,controls=None,num_steps=10,num_repeats=10):
#     """
#     Returns R, R_std, percs reach reach_std
#     R: MIP robustness measure
#     R_std: MIT robustness measure standard deviation
#     percs: list of percolation steps, in absolute number of nodes
#     reach: list of the percent of the network reachable corresponding to percs
#     reach_std: list of the standard deviation in reachability (percent of network) corresponding to percs
#
#     User can specify number of controls (to restrict it less than m)
#     User can specify how many percolation steps to make and how many repitiions (repititions are used to calculate std)
#     """
#     type_check(G,DiGraph,'only directed graphs are supported')
#
#     if controls is None:
#         C = build_cacti(G)
#         controls = C.controls_()
#     else:
#         C = build_cacti_fixed_controls_(G,controls)
#
#     orig_rank = C.num_controllable_nodes()
#
#     # Planned percolations
#     num_steps += 1
#     percs = floor(linspace(0,G.size(),num_steps))
#
#     # Create matrix to hold reachability values
#     reach = zeros((num_repeats,num_steps))
#     for i in arange(num_repeats):
#         #print 'repeat %i/%i' % (i,num_repeats)
#         Gi = cp.copy(G)
#         reach[i,0] = orig_rank
#         for j in arange(1,num_steps):
#             #print 'step %i/%i' % (j,num_steps)
#             # Percolate mod_percs edges at a time
#             mod_percs = percs[j]-percs[j-1]
#             for k in arange(mod_percs):
#                 # Percolate a random edge
#                 eidx = choose_edge_(Gi)
#                 Gi.rm_edge_(eidx)
#
#             # Calculate the reachability
#             Ci = build_cacti_fixed_controls_(Gi,controls)
#             reach[i,j] = Ci.num_controllable_nodes()
#
#     if len(controls) >= len(G):
#         R = 1
#         R_std = 0
#     else:
#         R = trapz(reach.mean(axis=0),percs)/(G.size()*len(G))
#         R_std = trapz(reach.std(axis=0),percs)/(G.size()*len(G))
#
#     return R, R_std, percs, reach.mean(axis=0), reach.std(axis=0)
