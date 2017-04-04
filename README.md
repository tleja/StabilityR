## StabilityR
A systems biology tool for network inference, partitioning and visualisation.

[//]: # (StabilityR is available at http://155.198.192.109:8080)

### Prerequisites

* LEMON graph library version 1.3.1: [http://lemon.cs.elte.hu/trac/lemon/wiki/Downloads](http://lemon.cs.elte.hu/trac/lemon/wiki/Downloads)
* Windows users, Rtools: [http://cran.r-project.org/bin/windows/Rtools](https://cran.r-project.org/bin/windows/Rtools)

### Installation 

From within the R console:

```coffee
require(devtools)
install_github('tleja/StabilityR')
```

### Example execution 

Load the library and an example graph

```coffee
library(StabilityR)
data(graph_rings)
mat = graph_rings
time = logspace(-3,1,61)
```

Execute the stability functions

```coffee
res1 = stability(mat=mat, time=time, full=TRUE, type='norm', M=100)
res2 = stability_evaluate(mat=mat, res=res1)
plotStability(res2)
plotStability_VI(res2)
plotStability_heatmap(res2)
```

Visualise the graph

```coffee
net = forceAtlas2(mat=mat, iterations=4000, plotstep=10, seed=4, ksmax=100)
plotNetwork(mat=mat, pos=net, group=res2$P[,55])
```
***
… more detailed documentation is coming soon.
