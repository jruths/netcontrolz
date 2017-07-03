"""
Justin Ruths and Derek Ruths. *Control Profiles of Complex Networks*. **Science**, 343(6177), 1373-1376 (2014).

url: http://science.sciencemag.org/content/343/6177/1373

The paper shows that the correlation between the minimal number of controls
required to guarantee controllability and the degree distribution is mainly due
to the number of sources and sinks in a network. The paper also presents a new
statistic, the control profile, which categorizes networks into three groups,
which correspond to the three components of the control profile.

**Abstract**: Studying the control properties of complex networks provides insight
into how designers and engineers can influence these systems to achieve a desired
behavior. Topology of a network has been shown to strongly correlate with certain
control properties; here we uncover the fundamental structures that explain the
basis of this correlation. We develop the control profile, a statistic that
quantifies the different proportions of control-inducing structures present in a
network. We find that standard random network models do not reproduce the kinds
of control profiles that are observed in real-world networks. The profiles of
real networks form three well-defined clusters that provide insight into the
high-level organization and function of complex systems.
"""
import zen
import netcontrolz
import matplotlib.pyplot as plt
from numpy import array, arange
from numpy.random import random

ER,BA,LA,DD = 'Erdos-Renyi','Barabasi-Albert','Local Attachment','Duplication Divergence'
model_colors = {ER:'r',BA:'b',LA:'c',DD:'m'}
N = 100 # number of nodes
num_repeats = 100
profiles = {m:[] for m in [ER,BA,LA,DD]}
profiles_shuff = {m:[] for m in [ER,BA,LA,DD]}
for r in range(num_repeats):
    c = random()*19+1 # average degree between 2 and 20
    p = c/N
    for i,m in enumerate([ER,BA,LA,DD]):
        if m == ER:
            G = zen.generating.erdos_renyi(N,p, directed=True)
        elif m == BA:
            G = zen.generating.barabasi_albert(N,int(c), directed=True)
        elif m == LA:
            G = zen.generating.local_attachment(N,int(c), int(1+random()*(c-1)))
        elif m == DD:
            G = zen.generating.duplication_divergence_iky(N, 0.1*random()*0.8, directed=True)
        profiles[m].append(netcontrolz.profile(G, normalized=False))
        profiles_shuff[m].append(netcontrolz.profile(zen.randomize.shuffle(G,keep_degree=True), normalized=False))

# showing the correlation with degree distribution is mostly due to sources and sinks
def plot_num_dilations(profiles,sprofiles):
    plt.plot([0,1],[0,1],'k-')
    srcsnk_err = []
    shuff_err = []
    for i in range(len(profiles)):
        p = array(profiles[i])
        nc = p.sum()/float(N) # actual number of dilations
        nsrcsnk = p[:-1].sum()/float(N) # prediction based on number of source and sink
        nshuff = array(sprofiles[i]).sum()/float(N) # prediction based on the number of dilations in a shuffled network
        srcsnk_err.append( abs(nc-nsrcsnk) )
        shuff_err.append( abs(nc-nshuff) )
        plt.plot(nc,nsrcsnk,'bx',ms=4,mew=1)
        plt.plot(nc,nshuff,'ro',ms=4,mfc='none',mew=1)
    return array(srcsnk_err), array(shuff_err)

srcsnk_errors = {m:[] for m in [ER,BA,LA,DD]}
shuff_errors = {m:[] for m in [ER,BA,LA,DD]}
plt.figure(figsize=(10,6))
for i,m in enumerate([ER,BA,LA,DD]):
    plt.subplot(2,3,i+1)
    sse, sfe = plot_num_dilations(profiles[m],profiles_shuff[m])
    plt.text(0.0,0.9,m)
    if i==3:
        plt.xlabel('$n_c$')
    if i==0 or i==3:
        plt.ylabel('prediction')
    srcsnk_errors[m] = sse
    shuff_errors[m] = sfe

width = 0.3
sep = 0.1
plt.subplot(2,3,6)
for i,m in enumerate([ER,BA,LA,DD]):
    plt.bar([i+1], srcsnk_errors[m].mean(), width, yerr=srcsnk_errors[m].std(), color='b' )
    plt.bar([i+1+width+sep], shuff_errors[m].mean(), width, yerr=shuff_errors[m].std(), color='r' )
ax = plt.gca()
ax.set_xticks(arange(4) + 1 + (width+sep)/2)
ax.set_xticklabels(('ER','BA','LA','DD'))
plt.ylabel('|$n_c$-prediction|')

plt.show()


# Showing classification of networks by control profiles
plt.figure(figsize=(16,4))
for i,m in enumerate([ER,BA,LA,DD]):
    plt.subplot(1,4,i+1)
    netcontrolz.profile_heatmap(profiles[m],color=model_colors[m])
    plt.title(m)
plt.show()
