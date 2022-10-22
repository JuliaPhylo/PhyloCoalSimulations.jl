```@setup correlated
using PhyloNetworks, PhyloPlots, PhyloCoalSimulations
figpath = joinpath("..", "assets", "figures"); mkpath(figpath)
figname(x) = joinpath(figpath, x)
net = readTopology("((C:0.9,(B:0.2)#H1:0.7::0.6)I1:0.6,(#H1:0.6::0.4,A:1.0)I2:0.5)I3;");
```

# correlated inheritance

The typical assumption made by most network inference methods is that
of independent inheritance at a hybrid node. It means that if multiple
lineages (of a given locus) at present at a hybrid node, then each
one is inherited from each parent according to the γ inheritance probabilities
*independently of the other lineages*.

Another model might be that the lineages at a hybrid node are *all* inherited
from the same parent, still choosing a (common) parent according to the γ
inheritance probabilities. This model has full (and positive) correlation
between lineages (as in Meng & Kubatko [20xx]()).

The other extreme might be interesting for modelling allopolyploid events:
of 2 lineages of the same locus at a hybrid node, exactly 1 of them comes
from one parent, and the other comes from the other parent. This model
would have negative correlation between lineages.

By default, `simulatecoalescent` uses the traditional model, with independent
lineages.
It also has an option to simulate lineage inheritance with positive correlation.
