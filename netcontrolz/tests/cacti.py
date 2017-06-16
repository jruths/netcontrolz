import unittest
from numpy.random import randint
from zen import *

class CactiUnweightedMatching(unittest.TestCase):

    def test_simple(self):

        G = DiGraph()
        G.add_edge('x','y')
        G.add_edge('x','z')
        G.add_edge('y','z')

        cacti = control.build_cacti(G)
        self.assertEqual(cacti.num_controls(),1)

        controls = cacti.controls_()
        self.assertEqual(controls[0][0],G.node_idx('x'))

        # change the graph
        G.rm_edge('y','z')
        cacti = control.build_cacti(G)
        self.assertEqual(cacti.num_controls(),2)
        controls = set(cacti.controls_())
        sol1 = set([(G.node_idx('x'),),(G.node_idx('z'),)])
        sol2 = set([(G.node_idx('x'),),(G.node_idx('y'),)])
        self.assertTrue(controls == sol1 or controls == sol2)

    def test_non_bud_cycles(self):
        # create a network consisting of three cycles
        G = DiGraph()
        G.add_edge(1,2)
        G.add_edge(2,1)
        c1 = set([G.node_idx(1),G.node_idx(2)])

        G.add_edge(3,4)
        G.add_edge(4,3)
        c2 = set([G.node_idx(3),G.node_idx(4)])

        G.add_edge(5,6)
        G.add_edge(6,7)
        G.add_edge(7,5)
        c3 = set([G.node_idx(5),G.node_idx(6),G.node_idx(7)])

        cacti = control.build_cacti(G)
        self.assertEqual(cacti.num_controls(),1)

        controls = set(cacti.controls_()[0])
        self.assertTrue(len(c1.intersection(controls)) == 1)
        self.assertTrue(len(c2.intersection(controls)) == 1)
        self.assertTrue(len(c3.intersection(controls)) == 1)

    def test_stem_cycle(self):
        G = DiGraph()
        G.add_edge(1,2)
        G.add_edge(2,3)
        G.add_edge(4,5)
        G.add_edge(5,4)

        cacti = control.build_cacti(G)
        self.assertEqual(cacti.num_controls(),1)

    def test_stem_bud(self):
        G = DiGraph()
        G.add_edge(1,2)
        G.add_edge(2,3)
        G.add_edge(2,4)
        G.add_edge(4,5)
        G.add_edge(5,4)

        cacti = control.build_cacti(G)
        self.assertEqual(cacti.num_controls(),1)


class CactiWeightedMatchingFixedControls(unittest.TestCase):
    '''
            graphs.append(G)
            G = barabasi_albert(N,k,directed=True,seed=int(time()+getpid()*1000))
            graphs.append(G)
            G = local_attachment(N,k,max(int(.25*k),1))
            graphs.append(G)
            G = local_attachment(N,k,max(int(.5*k),1))
            graphs.append(G)
            G = local_attachment(N,k,max(int(.75*k),1))
            graphs.append(G)
    '''
    def _test_with_kalman(self, G, ctls=None):
        if ctls is None:
            ctls = set( randint(len(G)) for _ in xrange(randint(1,len(G))) )
            ctls = [tuple([a]) for a in ctls]

        cact = control.build_cacti_fixed_controls_(G, ctls, randomize=True)

        NC = cact.num_controllable_nodes()
        K = control.reachability.kalman_generic_rank(G, ctls)

        self.assertEqual(NC, K)

        return cact


    def test_random_graphs_ER(self):
        N = 20
        for k in [2,4,6]:
            G = generating.erdos_renyi(N, float(k)/N, directed=True, seed=2**20)

            self._test_with_kalman(G)

    def test_random_graphs_BA(self):
        N = 20
        for k in [2,4,6]:
            G = generating.barabasi_albert(N,k,directed=True,seed=2**20)

            self._test_with_kalman(G)

    def test_stem_cycle1(self):
        G = DiGraph()
        G.add_nodes(5)
        G.add_edge_(0,1)
        G.add_edge_(1,2)
        G.add_edge_(3,4)
        G.add_edge_(4,3)

        ctls = [(0,)]
        cact = control.build_cacti_fixed_controls_(G, ctls, randomize=True, with_cycles=True)
        self.assertEqual(cact.num_controllable_nodes(),5)

    def test_stem_cycle2(self):
        G = DiGraph()
        G.add_nodes(5)
        G.add_edge_(0,1)
        G.add_edge_(1,2)
        G.add_edge_(3,4)
        G.add_edge_(4,3)

        ctls = [(0,)]
        cact = control.build_cacti_fixed_controls_(G, ctls, randomize=True)
        self.assertEqual(cact.num_controllable_nodes(),3)

    def test_stem_cycle3(self):
        G = DiGraph()
        G.add_nodes(5)
        G.add_edge_(0,1)
        G.add_edge_(1,2)
        G.add_edge_(3,4)
        G.add_edge_(4,3)

        ctls = [(3,)]
        cact = control.build_cacti_fixed_controls_(G, ctls, randomize=True)
        self.assertEqual(cact.num_controllable_nodes(),2)

    def test_stem_bud1(self):
        G = DiGraph()
        G.add_nodes(5)
        G.add_edge_(0,1)
        G.add_edge_(1,2)
        G.add_edge_(1,3)
        G.add_edge_(3,4)
        G.add_edge_(4,3)

        ctls = [(0,)]
        cact = control.build_cacti_fixed_controls_(G, ctls, randomize=True)
        self.assertEqual(cact.num_controllable_nodes(),5)

    def test_stem_bud2(self):
        G = DiGraph()
        G.add_nodes(5)
        G.add_edge_(0,1)
        G.add_edge_(1,2)
        G.add_edge_(1,3)
        G.add_edge_(3,4)
        G.add_edge_(4,3)

        ctls = [(1,)]
        cact = control.build_cacti_fixed_controls_(G, ctls, randomize=True)
        self.assertEqual(cact.num_controllable_nodes(),4)

    def test_cacti_form1(self):
        G = DiGraph()
        G.add_nodes(5)
        G.add_edge_(0,1)
        G.add_edge_(1,2)
        G.add_edge_(1,3)
        G.add_edge_(3,4)
        G.add_edge_(4,3)

        ctls = [(0,)]
        cact = control.build_cacti_fixed_controls_(G, ctls, randomize=True)
        stems, cycles = cact.stems(), cact.cycles()

        self.assertEqual(len(stems), 1)
        self.assertEqual(len(cycles), 1)
        self.assertEqual(stems[0].origin_(), 0)
        self.assertEqual(cycles[0].origin_(), 3)
        self.assertEqual(len(stems[0].buds()), 1)

        self.assertIsNotNone(cycles[0].stem())
        self.assertEqual(cycles[0].dist_node_(),1)

if __name__ == '__main__':
    unittest.main()
