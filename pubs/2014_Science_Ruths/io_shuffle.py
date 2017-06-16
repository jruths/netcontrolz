from random import choice,randint
import zen

##
# Generates a shuffled version of the network in which the number of 
# source, sink, & isolated nodes are the same, but the rest is randomized

def count_nodes(G):
	src_nodes = set([])
	snk_nodes = set([])
	iso_nodes = set([])
	
	for n,data in G.nodes_iter(True):
		if G.in_degree(n) == 0:
			if G.out_degree(n) == 0:
				iso_nodes.add(n)
			else:
				src_nodes.add(n)
		elif G.out_degree(n) == 0:
			snk_nodes.add(n)
	
	#print '> src = %i, snk = %i, iso = %i' % (len(src_nodes),len(snk_nodes),len(iso_nodes))
	return src_nodes, snk_nodes, iso_nodes
	
def io_shuffle(G,self_loops=False):
	
	dG = zen.DiGraph()
	src_nodes, snk_nodes, iso_nodes = count_nodes(G)
	
	for n,data in G.nodes_iter(True):
		if not ((n in src_nodes) or (n in snk_nodes) or (n in iso_nodes)):
			dG.add_node(n,data)
	
	# shuffle edges
	all_nodes = dG.nodes()
	for e,data,weight in G.edges_iter_(-1,True,True):
		# pick two random nodes to connect
		x = choice(all_nodes)
		y = choice(all_nodes)
		while (not self_loops and x == y) or dG.has_edge(x,y):
			x = choice(all_nodes)
			y = choice(all_nodes)

		e = dG.add_edge(x,y,data)
		dG.set_weight_(e,weight)
	
	L = dG.size()
	# we could have created a few sinks, sources, or isolates along the way
	# move a few edges so that all nodes have in & out degree >= 1
	tmp_src_nodes, tmp_snk_nodes, tmp_iso_nodes = count_nodes(dG)
	
	# for each isolated node - pick two edges such that:
	# "from" node has >2 out degree and "to" node has >2 in degree
	for n in tmp_iso_nodes:
		e1 = randint(0,L-1)
		(x1,y1) = dG.endpoints(e1)
		while (dG.out_degree(x1) < 2) or (dG.in_degree(y1) < 2):
			e1 = randint(0,L-1)
			(x1,y1) = dG.endpoints(e1)
		
		e2 = randint(0,L-1)
		(x2,y2) = dG.endpoints(e2)
		while (dG.out_degree(x2) < 2) or (dG.in_degree(y2) < 2) or (e1 == e2):
			e2 = randint(0,L-1)
			(x2,y2) = dG.endpoints(e2)
		
		# Three mechanisms to bring in an isolated node (1,2 permutations are not important):
		#  1. n -> n  (remove just x1 -> y1)
		#  2. x1 -> n -> y1 (remove x1 -> y1 and x2 -> y2)
		#  3. x1 -> n -> y2 (remove x1 -> y1 and x2 -> y2)
		r = randint(1,3)
		
		data1 = dG.edge_data_(e1)
		weight1 = dG.weight_(e1)
		dG.rm_edge_(e1)
		
		if r == 1:
			dG.add_edge(n,n,data1,weight1)
		else:
			data2 = dG.edge_data_(e2)
			weight2 = dG.weight_(e2)
			dG.rm_edge_(e2)
			
			if r == 2:
				dG.add_edge(x1,n,data1,weight1)
				dG.add_edge(n,y1,data2,weight2)
			else:
				dG.add_edge(x1,n,data1,weight1)
				dG.add_edge(n,y2,data2,weight2)
	
	# for each source node - pick a node with 2 or more in degree
	for n in tmp_src_nodes:
		y = choice(all_nodes)
		while dG.in_degree(y) < 2:
			y = choice(all_nodes)
		# find an inbound neighbor of y
		x = choice(dG.in_neighbors(y))
		data = dG.edge_data(x,y)
		weight = dG.weight(x,y)
		dG.rm_edge(x,y)
		dG.add_edge(x,n,data,weight)
	
	# for each sink node - pick an node with 2 or more out degree
	for n in tmp_snk_nodes:
		x = choice(all_nodes)
		while dG.out_degree(x) < 2:
			x = choice(all_nodes)
		# find an outbound neighbor of x
		y = choice(dG.out_neighbors(x))
		data = dG.edge_data(x,y)
		weight = dG.weight(x,y)
		dG.rm_edge(x,y)
		dG.add_edge(n,y,data,weight)
	
	#count_nodes(dG)
	# now we have a network with no sources, sinks or isolated nodes
	
	# add isolated nodes
	for n in iso_nodes:
		dG.add_node(n,G.node_data(n))
	
	# add source nodes
	for n in src_nodes:
		dG.add_node(n,G.node_data(n))
		# pick a random edge
		e = randint(0,L-1)
		(x,y) = dG.endpoints(e)
		while (x in src_nodes) or (dG.out_degree(x) < 2):
			e = randint(0,L-1)
			(x,y) = dG.endpoints(e)
		# remove the edge and reattached with src node replacing the "from" node
		data = dG.edge_data_(e)
		weight = dG.weight_(e)
		dG.rm_edge_(e)
		dG.add_edge(n,y,data,weight)
		
	for n in snk_nodes:
		dG.add_node(n,G.node_data(n))
		# pick a random edge
		e = randint(0,L-1)
		(x,y) = dG.endpoints(e)
		while (y in snk_nodes) or (dG.in_degree(y) < 2):
			e = randint(0,L-1)
			(x,y) = dG.endpoints(e)
		# remove the edge and reattached with sink node replacing the "to" node
		data = dG.edge_data_(e)
		weight = dG.weight_(e)
		dG.rm_edge_(e)
		dG.add_edge(x,n,data,weight)
		
	return dG