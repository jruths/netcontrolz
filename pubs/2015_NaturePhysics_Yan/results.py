"""
Gang Yan, Georgios Tsekenis, Baruch Barzel, Jean-Jacques Slotine, Yang-Yu Liu, Albert-László Barabási. *Spectrum of controlling and observing complex networks*. **Nature Physics** 11, 779–786 (2015).

url: http://www.nature.com/nphys/journal/v11/n9/abs/nphys3422.html

This paper constructs the distribution of energies corresponding to the infinite
horizon controllability Gramian. It observes largely scale-free distributions
of energies when one or all nodes are driven by an external control(s). When
an intermediate fraction of nodes is controlled, multi-peak distributions are
observed and a gap between distinct bands of energies arises.

**Abstract**: Recent studies have made important advances in identifying sensor or
driver nodes, through which we can observe or control a complex system. But the
observational uncertainty induced by measurement noise and the energy required
for control continue to be significant challenges in practical applications.
Here we show that the variability of control energy and observational uncertainty
for different directions of the state space depend strongly on the number of
driver nodes. In particular, we find that if all nodes are directly driven,
control is energetically feasible, as the maximum energy increases sublinearly
with the system size. If, however, we aim to control a system through a single
node, control in some directions is energetically prohibitive, increasing
exponentially with the system size. For the cases in between, the maximum energy
decays exponentially when the number of driver nodes increases. We validate our
findings in several model and real networks, arriving at a series of fundamental
laws to describe the control energy that together deepen our understanding of
complex systems.
"""
import zen
import netcontrolz
from numpy import abs, zeros, log10
from scipy.linalg import eigvals, inv, eigvalsh
import matplotlib.pyplot as plt
from numpy.random import choice

def plot_pdf(data):
    """
    plots the averaged energy probability distribution given the spectra of
    several different realizations (weight choices) of the matrix A
    """
    pdfs, bins, patches = plt.hist(data, bins=25, normed=True, log=True)
    plt.cla()

    # calculate the center of the bins
    cbins = zeros(len(bins)-1)
    for i in range(len(cbins)):
        cbins[i] = 0.5*(bins[i+1] - bins[i]) + bins[i]

    # assemble all the pdfs into a 2D matrix
    pdf = zeros((len(pdfs),len(bins)-1))
    for i, pdfi in enumerate(pdfs):
        pdf[i,:] = pdfi

    plt.errorbar(cbins,pdf.mean(axis=0),yerr=pdf.std(axis=0),mec='b',mfc='w',marker='o',ms=4,lw=1)

# LOAD/BUILD A NETWORK
G = zen.generating.barabasi_albert(5000,5)

num_realizations = 10 # number of times to generate new weights
n = G.num_nodes
As = G.matrix().T # the sparsity pattern of the network

## ========================================================
# controlling every node
B = netcontrolz.input_matrix_(n,[(i,) for i in range(n)])
data = []
for r in range(num_realizations):
    A = netcontrolz.matrix_realization(As,0,1)
    A = netcontrolz.fix_diagonals(A,0.25)

    H = netcontrolz.infinite_horizon_gramian(A,B)
    d = abs(eigvalsh(H))
    data.append(1.0/d)

plt.figure()
plot_pdf(data)
ax = plt.gca()
ax.set_xscale("log")
ax.set_yscale("log")
plt.xlabel('$\mathcal{E}$')
plt.ylabel('$p(\mathcal{E})$')
plt.title('$N_D = N$')

## ========================================================
# controlling a single (random) node
data = []
for r in range(num_realizations):
    A = netcontrolz.matrix_realization(As,0,1)
    A = netcontrolz.fix_diagonals(A,0.25)
    # pick a random node to control
    B = netcontrolz.input_matrix_(n,[(zen.choose_node_(G),)])

    H = netcontrolz.infinite_horizon_gramian(A,B)
    d = abs(eigvalsh(H))
    data.append(1.0/d)

plt.figure()
plot_pdf(data)
ax = plt.gca()
ax.set_xscale("log")
ax.set_yscale("log")
plt.xlabel('$\mathcal{E}$')
plt.ylabel('$p(\mathcal{E})$')
plt.title('$N_D = 1$')

## ========================================================
# controlling a fraction of nodes
fraction = 0.2
data = []
for r in range(num_realizations):
    A = netcontrolz.matrix_realization(As,0,1)
    A = netcontrolz.fix_diagonals(A,0.25)
    B = netcontrolz.input_matrix_(n,[(i,) for i in choice(n,int(fraction*n),replace=False)])

    H = netcontrolz.infinite_horizon_gramian(A,B)
    d = abs(eigvalsh(H))

    data.append(log10(1.0/d))

plt.figure()
plot_pdf(data)
plt.xlabel('$\log(\mathcal{E})$')
plt.ylabel('$p[\log(\mathcal{E})]$')
plt.title('$N_D / N = %1.1f$' % fraction)

## ========================================================
# show figures
plt.show()
