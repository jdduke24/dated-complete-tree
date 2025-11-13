import os
os.chdir('/Users/Jonathan/Library/CloudStorage/Dropbox/Imperial/Tree_of_Life/Open_Tree/python/dated-complete-tree/figures/figure1')

from pylab import *

import numpy as np

fin = open("bladj-v-me.csv", 'r')

spp = []
bladj = []
me = []
me_new = []

count = 0
for line in fin:
    if count == 0:
        count += 1
        continue

    parts = line.split(',')

    spp.append(int(parts[0]))
    bladj.append(float(parts[1]))
    me.append(float(parts[2]))
    me_new.append(float(parts[3].strip()))

    count += 1

fin.close()

fig = figure(figsize=(7,6))
ax1 = fig.add_subplot(111)

ax1.scatter(spp[:-3], bladj[:-3], label="Phylocom EQS-L implementation (BLADJ)")
ax1.scatter(spp[-3:], bladj[-3:], marker='o', facecolors="none", edgecolors="#1f77b4", label="Phylocom inferred run times")
ax1.scatter(spp, me, marker="^", label="Our EQS-L implementation")
ax1.scatter(spp, me_new, marker="+", label="Our EQS-L with parsimonious date cleaning", s=60)


ax1.plot([spp[0], spp[-1]], [bladj[0]*3, bladj[0]*spp[-1]**2/spp[0]**2], ':', label="Quadratic growth trend", zorder=-100)
ax1.plot([spp[0], spp[-1]], [me[0], me[0]*spp[-1]/spp[0]], '--', label="Linear growth trend", zorder=-100)

ax1.set_xscale("log")
ax1.set_yscale("log")

ax1.set_xlabel("Number of nodes (log scale)")
ax1.set_ylabel("Run time (minutes, log scale)")

ax1.legend(loc="upper left", frameon=False, fontsize="small")

ax1.text(spp[-1]*1.25,bladj[-1],"20 days", color="black", verticalalignment="center", size="small")
ax1.text(spp[-1]*1.25,me[-1],"22 seconds", color="black", verticalalignment="center", size="small")

x, y = ax1.get_xlim()
ax1.set_xlim(x, y*4)

fig.savefig("bladj_v_me_updated.pdf")
fig.savefig("bladj_v_me_updated.tif", dpi=300)
