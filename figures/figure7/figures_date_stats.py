import os
os.chdir('/Users/Jonathan/Library/CloudStorage/Dropbox/Imperial/Tree_of_Life/Open_Tree/python/dated-complete-tree')

import matplotlib.pyplot as plt
import numpy as np

indents = [
0,
0,
 1,
     2,
         4,
         4,
         4,
         4,
         4,
         4,
         4,
     2,
        3,
           4,
               5,
               5,
               5,
               5,
               5,
               # 5,
               # 5,
               5,
           4,
           4,
           4,
        3,
        3,
     2,
     2,
     2,
     2,
     2,
     2,
     2,
     2,
     2,
 1,
     2,
         4,
         4,
     2,
     2,
     2,
 1,
 1,
     2,
     2,
     2,
 1,
     2,
     2,
     2,
0,
0]

for i in range(1,len(indents)):
    indents[i] += 1

import csv

names = []
spp = []
dtnds = []
dtcov = []
dtsrc = []
phycov = []

with open("figures/figure7/date_stats.csv", newline='') as csvfile:
    rdr = csv.reader(csvfile)
    for idx, line in enumerate(rdr):
        if idx == 0:
            # first line has column headings
            continue
        names.append(line[0].strip())
        spp.append(line[1])
        dtnds.append(line[2])
        dtcov.append(line[3])
        dtsrc.append(line[4])
        phycov.append(line[5])

names.reverse()
indents.reverse()
spp.reverse()
dtnds.reverse()
dtcov.reverse()
dtsrc.reverse()
phycov.reverse()

################
# Plotting

fig = plt.figure(figsize=(5.2,7.5))

ax0 = fig.add_subplot(1,1,1)

green = "#DAF2D0"
blue = "#DAE9F8"
grey ="#A6A6A6"
lightgrey ="#E0E0E0"
moregrey = "#888888"

col_x = [0.01,
         0.39,
         0.54,
         0.68,
         0.83,
         0.99]

initial_indent = col_x[0]
indent_size = 0.017

for i in range(len(names)):
    ax0.text(col_x[0]+indents[i]*indent_size, i+0.95, names[i], horizontalalignment="left", verticalalignment="center", family=["Arial"], fontsize="x-small")
    ax0.text(col_x[1], i+0.95, spp[i], horizontalalignment="right", verticalalignment="center", family=["Arial"], fontsize="x-small")
    ax0.text(col_x[2], i+0.95, dtnds[i], horizontalalignment="right", verticalalignment="center", family=["Arial"], fontsize="x-small")
    ax0.text(col_x[3], i+0.95, dtcov[i], horizontalalignment="right", verticalalignment="center", family=["Arial"], fontsize="x-small")
    ax0.text(col_x[4], i+0.95, dtsrc[i], horizontalalignment="right", verticalalignment="center", family=["Arial"], fontsize="x-small")
    ax0.text(col_x[5], i+0.95, phycov[i], horizontalalignment="right", verticalalignment="center", family=["Arial"], fontsize="x-small")

    ax0.axhline(i+0.5, color=grey, lw=0.5, zorder=1000)

    if i >= 2 and i < len(names)-2:
        ax0.plot([initial_indent+(indents[i]-1)*indent_size+0.002, initial_indent+(indents[i]-1)*indent_size+0.009], [i+1, i+1], lw=0.4, color=moregrey)

        for j in range(1,indents[i]):
            if ((i > 40 and i < 48) or (i > 13 and i < 16)) and j == 3:
                continue

            if i < 5 and j == 1:
                continue

            if (indents[i-1] < indents[i] and j == indents[i]-1) or i == 5:
                ax0.plot([initial_indent+j*indent_size+0.002, initial_indent+j*indent_size+0.002], [i+1, i+1.5], lw=0.4, color=moregrey)
            else:
                ax0.plot([initial_indent+j*indent_size+0.002, initial_indent+j*indent_size+0.002], [i+0.5, i+1.5], lw=0.4, color=moregrey)


i = len(names)
ax0.text(col_x[0],i+0.65,"Clade", horizontalalignment="left", verticalalignment="bottom", family=["Arial"], weight="bold", fontsize="x-small")
ax0.text(col_x[1],i+0.65,"Species\nrichness", horizontalalignment="right", verticalalignment="bottom", family=["Arial"], weight="bold", fontsize="x-small")
ax0.text(col_x[2],i+0.65,"Dated\nnodes", horizontalalignment="right", verticalalignment="bottom", family=["Arial"], weight="bold", fontsize="x-small")
ax0.text(col_x[3],i+0.65,"Date\ncoverage", horizontalalignment="right", verticalalignment="bottom", family=["Arial"], weight="bold", fontsize="x-small")
ax0.text(col_x[4],i+0.655,"Date source\ntrees", horizontalalignment="right", verticalalignment="bottom", family=["Arial"], weight="bold", fontsize="x-small")
ax0.text(col_x[5],i+0.65,"Phylogenetic\ncoverage", horizontalalignment="right", verticalalignment="bottom", family=["Arial"], weight="bold", fontsize="x-small")
ax0.axhline(i+0.5, color=grey, lw=1)

ax0.get_yaxis().set_ticks([])
ax0.get_xaxis().set_ticks([])

ax0.axhspan(len(names)-0.5, len(names)+0.5, color=blue, alpha=1, zorder=-1000)

for i in range(len(names)-1):
    if i % 2 == 1:
        ax0.axhspan(i+1.5, i+0.5, color=lightgrey, alpha=0.35, zorder=-1000)

ax0.set_xlim(0,1)
ax0.set_ylim(0,len(names)+1.5)

ax0.set_frame_on(False)

fig.tight_layout()

import PIL
fig.savefig("figures/figure7/figure7_coverage_table.pdf")
fig.savefig("figures/figure7/figure7_coverage_table.tif", dpi=300)
