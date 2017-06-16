import copy as cp

from numpy import floor,zeros,arange,trapz,linspace

from cacti import build_cacti, build_cacti_fixed_controls_
from zen import choose_edge_

def mip_(G,controls=None,num_steps=10,num_repeats=10):
	"""
	Returns R, R_std, percs reach reach_std
	R: MIP robustness measure
	R_std: MIT robustness measure standard deviation 
	percs: list of percolation steps, in absolute number of nodes
	reach: list of the percent of the network reachable corresponding to percs
	reach_std: list of the standard deviation in reachability (percent of network) corresponding to percs
	
	User can specify number of controls (to restrict it less than m) 
	User can specify how many percolation steps to make and how many repitiions (repititions are used to calculate std)
	"""
	if controls is None:
		C = build_cacti(G)
		controls = C.controls_()
	else:
		C = build_cacti_fixed_controls_(G,controls)
	
	orig_rank = C.num_controllable_nodes()

	# Planned percolations
	num_steps += 1
	percs = floor(linspace(0,G.size(),num_steps))

	# Create matrix to hold reachability values
	reach = zeros((num_repeats,num_steps))
	for i in arange(num_repeats):
		#print 'repeat %i/%i' % (i,num_repeats)
		Gi = cp.copy(G)
		reach[i,0] = orig_rank
		for j in arange(1,num_steps):
			#print 'step %i/%i' % (j,num_steps)
			# Percolate mod_percs edges at a time
			mod_percs = percs[j]-percs[j-1]
			for k in arange(mod_percs):
				# Percolate a random edge
				eidx = choose_edge_(Gi)
				Gi.rm_edge_(eidx)

			# Calculate the reachability
			Ci = build_cacti_fixed_controls_(Gi,controls)
			reach[i,j] = Ci.num_controllable_nodes()

	if len(controls) >= len(G):
		R = 1
		R_std = 0
	else:
		R = trapz(reach.mean(axis=0),percs)/(G.size()*len(G))
		R_std = trapz(reach.std(axis=0),percs)/(G.size()*len(G))

	return R, R_std, percs, reach.mean(axis=0), reach.std(axis=0)
