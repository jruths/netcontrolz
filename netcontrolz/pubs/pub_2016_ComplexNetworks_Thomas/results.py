"""
Jijju Thomas, Supratim Ghosh, Deven Parek, Derek Ruths, Justin Ruths. *Robustness of Network Controllability to Degree-Based Edge Attacks*. Complex Networks & Their Applications V, 2016.

https://link.springer.com/chapter/10.1007/978-3-319-50901-3_42

This paper characterizes the robustness properties of Erdos-Renyi and Barabasi-Albert
random network models to various types of degree-based edge attacks as well as to
random edge failure.

**Abstract**: We analyze the tolerance of network controllability to degree-based edge
attacks as well as random edge failure. In particular, we leverage both
control-based and reachability-based robustness metrics to investigate the case
when a fixed number of controls are allowed to change locations following each attack.
This ability to change the locations of controls models the more realistic scenario
in which operators may have a fixed budget of resources but that these resources can
be redeployed in response to attacks on the system. We also identify that the most
potent targeted attack for network controllability selects edges (on average) based
on betweenness centrality.
"""

def highest_betweenness(G):
    """
    Returns the edge index that has the highest betweeness centrality.
    Currently we rely on the netowrkx to calculate this statistic.
    """
    Gnx = zen.nx.to_networkx(G)
    edge_btw = networkx.edge_betweenness_centrality(Gnx,normalized = False)
    high_btw_endpts = sorted(edge_btw.items(), key=lambda x: x[1],reverse = True)[0][0]
    eidx = G.edge_idx(high_btw_endpts[0],high_btw_endpts[1])
    return eidx

def calc_robustness_random_model(model,c,num_repeats):
    """
    Calculates and plots all robustness metrics for a specific ER/BA random ``model``
    with average degree ``c``. The plots average ``num_repeats`` number of percolation
    processes to attempt to get the expected behavior of the model.
    """
    p = float(c)/(N-1)
    q = c

    num_steps = 10
    # control robustness
    cR = {attack:zeros((num_repeats,num_steps+1)) for attack in attacks}
    # reachability robustness (fixed)
    rRfixed = {attack:zeros((num_repeats,num_steps+1)) for attack in attacks}
    # reachability robustness (free)
    rRfree = {attack:zeros((num_repeats,num_steps+1)) for attack in attacks}
    for r in range(num_repeats):
        if model == ER:
            G = zen.generating.erdos_renyi(N,p,directed=True)
        elif model == BA:
            G = zen.generating.barabasi_albert(N,q,directed=True)
        else:
            return None
        for attack in attacks:
            if attack == ATTACK_BETWEEN:
                a = highest_betweenness
            else:
                a = attack
            ell, M = netcontrolz.edge_percolation_(G,a,controls=None,frac_rm_edges=0.9,num_steps=10)
            cR[attack][r,:] = M[netcontrolz.CONTROL_ROBUSTNESS]/float(N)
            rRfixed[attack][r,:] = M[netcontrolz.REACHABILITY_ROBUSTNESS_FIXED]/float(N)
            rRfree[attack][r,:] = M[netcontrolz.REACHABILITY_ROBUSTNESS_FREE]/float(N)

    plt.figure(figsize=(15,4))
    for attack in attacks:
        plt.subplot(1,3,1)
        plt.errorbar(ell,cR[attack].mean(axis=0),yerr=cR[attack].std(axis=0), color=attack_colors[attack])
        plt.subplot(1,3,2)
        plt.errorbar(ell,rRfixed[attack].mean(axis=0),yerr=rRfixed[attack].std(axis=0), color=attack_colors[attack])
        plt.subplot(1,3,3)
        plt.errorbar(ell,rRfree[attack].mean(axis=0),yerr=rRfree[attack].std(axis=0), color=attack_colors[attack])

    plt.subplot(1,3,1)
    plt.xlabel('$\ell/L$')
    plt.ylabel('$n_c$')
    plt.subplot(1,3,2)
    plt.xlabel('$\ell/L$')
    plt.ylabel('$n_r$')
    plt.subplot(1,3,3)
    plt.xlabel('$\ell/L$')
    plt.ylabel('$n_f$')
    plt.show()

if __name__ == 'main':
    import zen
    import netcontrolz
    from numpy import zeros
    import matplotlib.pyplot as plt
    import networkx

    ER = 'er' # to indicate the Erdos Renyi model
    BA = 'ba' # to indicate the Barabasi-Albert model
    N = 100 # number of nodes

    # beyond the random and degree-based attacks, add a custom attack that
    # removes high-betweeness edges.
    ATTACK_BETWEEN = 'betweenness'
    attacks = [netcontrolz.ATTACK_RAND,
               netcontrolz.EDGE_ATTACK_INOUT_DEG,
               netcontrolz.EDGE_ATTACK_OUTIN_DEG,
               netcontrolz.EDGE_ATTACK_ININ_DEG,
               netcontrolz.EDGE_ATTACK_OUTOUT_DEG,
               netcontrolz.EDGE_ATTACK_TOTAL_DEG,
               ATTACK_BETWEEN]

    # colors to use in the plots
    attack_colors = {
        netcontrolz.ATTACK_RAND: 'm',
        netcontrolz.EDGE_ATTACK_INOUT_DEG: 'g',
        netcontrolz.EDGE_ATTACK_OUTIN_DEG: 'r',
        netcontrolz.EDGE_ATTACK_ININ_DEG: 'b',
        netcontrolz.EDGE_ATTACK_OUTOUT_DEG: 'c',
        netcontrolz.EDGE_ATTACK_TOTAL_DEG: 'y',
        ATTACK_BETWEEN: 'k'
    }

    # select a model, average degree, and number of graphs that should be averaged together
    calc_robustness_random_model(BA,2,10)
