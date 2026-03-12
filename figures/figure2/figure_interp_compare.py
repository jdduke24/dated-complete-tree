import os
os.chdir('/Users/Jonathan/Library/CloudStorage/Dropbox/Imperial/Tree_of_Life/Open_Tree/python/dated-complete-tree/figures/figure2')

import pylab
import numpy as np

import pickle

fmc = 37
cmc = 0

with open('fullydated_avgs.pickle', 'rb') as f:
    fx, fy = pickle.load(f)

with open('classes_avgs.pickle', 'rb') as f:
    cx, cy = pickle.load(f)

algos = ["EQS-L", "EQS-S", "LnN", "BM", "EQS-LS"]

fig = pylab.figure(figsize=(7.5,7.5))

ax = fig.add_subplot(3,2,1)
for i, algo in enumerate(fy):
    ax.plot(fx[fmc:], fy[algo][fmc:], '-', label=algos[i])

ax.legend(loc="upper left", frameon=False, prop={"size":"small", "family":"Arial"})
ax.set_xlim((0,100))
ax.set_ylim((0,18))
ax.set_xticklabels(ax.get_xticklabels(), fontfamily="Arial")
ax.set_yticks([0, 2, 4, 6, 8, 10, 12, 14, 16, 18])
ax.set_yticklabels(["0%", "2%", "4%", "6%", "8%", "10%", "12%", "14%", "16%", "18%"], fontfamily="Arial")
ax.set_ylabel(r"$\mathbf{Date\ reproduction\ error}$" + "\n$abs$(orig. date - interp. date) / crown date", fontsize="x-small", fontfamily="Arial")

ymin, ymax = ax.get_ylim()
ax.text(95, ymax*0.95, "A1", fontsize="large", family=["Arial"], fontweight="bold", verticalalignment="top", horizontalalignment="right")

ax.text(50,19,'Fully-dated trees', family=["Arial"], horizontalalignment="center", fontsize="medium", fontweight="bold")

low = 65
high = 65.3

ax = fig.add_subplot(3,2,2)
for i, algo in enumerate(cy):
    ax.plot(cx[cmc:], cy[algo][cmc:], '-', label=algos[i])

ax.legend(loc="upper left", frameon=False, prop={"size":"small", "family":"Arial"})
ax.set_xlim((0,100))
ax.set_ylim((0,18))
ax.set_xticklabels(ax.get_xticklabels(), fontfamily="Arial")
ax.set_yticks([0, 2, 4, 6, 8, 10, 12, 14, 16, 18])
ax.set_yticklabels(["0%", "2%", "4%", "6%", "8%", "10%", "12%", "14%", "16%", "18%"], fontfamily="Arial")
ax.set_ylabel("  \n  ")

ymin, ymax = ax.get_ylim()
ax.text(95, ymax*0.95, "B1", fontsize="large", family=["Arial"], fontweight="bold", verticalalignment="top", horizontalalignment="right")

ax.text(50,19, 'Partially-dated clades', family=["Arial"], horizontalalignment="center", fontsize="medium", fontweight="bold")

########

# Pen means

with open('fullydated_pens.pickle', 'rb') as f:
    fx, fy = pickle.load(f)

with open('classes_pens.pickle', 'rb') as f:
    cx, cy = pickle.load(f)

ax = fig.add_subplot(3,2,3)
for i, algo in enumerate(fy):
    ax.plot(fx[fmc:], fy[algo][fmc:], '-', label=algos[i])

# ax.set_yscale("log", base=2)
# ax.set_title("Average across four 100%-dated clades", fontsize="medium")
ax.legend(loc="upper left", frameon=False, prop={"size":"small", "family":"Arial"})
ax.set_xlim((0,100))
ax.set_ylim((0,7))
ax.set_xticklabels(ax.get_xticklabels(), fontfamily="Arial")
ax.set_yticks([0, 1, 2, 3, 4, 5, 6, 7])
ax.set_xticklabels(ax.get_xticklabels(), fontfamily="Arial")
ax.set_yticklabels(["0", "1", "2", "3", "4", "5", "6", "7"], fontfamily="Arial")
ax.axhline(1, ls='--', lw=1, color="gray")
# ax.set_xlabel("Percentage of internal nodes interpolated")
ax.set_ylabel(r"$\mathbf{Mean\ pendant\ edge\ length}$" + "\n(relative to EQS-LS,\nMyr divided by initial value for EQS-LS)", fontsize="x-small", fontfamily="Arial")

ymin, ymax = ax.get_ylim()
ax.text(95, ymax*0.95, "A2", fontsize="large", family=["Arial"], fontweight="bold", verticalalignment="top", horizontalalignment="right")


ax = fig.add_subplot(3,2,4)
for i, algo in enumerate(cy):
    ax.plot(cx[cmc:], cy[algo][cmc:], '-', label=algos[i])

ax.legend(loc="upper left", frameon=False, prop={"size":"small", "family":"Arial"})
ax.set_xlim((0,100))
ax.set_ylim((0,7))
ax.set_xticklabels(ax.get_xticklabels(), fontfamily="Arial")
ax.set_yticks([0, 1, 2, 3, 4, 5, 6, 7])
ax.set_yticklabels(["0", "1", "2", "3", "4", "5", "6", "7"], fontfamily="Arial")
ax.axhline(1, ls='--', lw=1, color="gray")
ax.set_ylabel("  \n  ")

ymin, ymax = ax.get_ylim()
ax.text(95, ymax*0.95, "B2", fontsize="large", family=["Arial"], fontweight="bold", verticalalignment="top", horizontalalignment="right")


########

# Pen medians

with open('fullydated_penmeds.pickle', 'rb') as f:
    fx, fy = pickle.load(f)

with open('classes_penmeds.pickle', 'rb') as f:
    cx, cy = pickle.load(f)

ax = fig.add_subplot(3,2,5)
for i, algo in enumerate(fy):
    ax.plot(fx[fmc:], fy[algo][fmc:], '-', label=algos[i])

ax.legend(loc="upper left", frameon=False, prop={"size":"small", "family":"Arial"})
ax.set_xlim((0,100))
ax.set_ylim((0,11))
ax.set_xticklabels(ax.get_xticklabels(), fontfamily="Arial")
ax.set_yticks([0, 2, 4, 6, 8, 10])
ax.set_yticklabels(["0", "2", "4", "6", "8", "10"], fontfamily="Arial")
ax.axhline(1, ls='--', lw=1, color="gray")
ax.set_xlabel("Percentage of internal nodes interpolated", fontsize="medium", fontfamily="Arial")
ax.set_ylabel(r"$\mathbf{Median\ pendant\ edge\ length}$" + "\n(relative to EQS-LS,\nMyr divided by initial value for EQS-LS)", fontsize="x-small", fontfamily="Arial")

ymin, ymax = ax.get_ylim()
ax.text(95, ymax*0.95, "A3", fontsize="large", family=["Arial"], fontweight="bold", verticalalignment="top", horizontalalignment="right")


ax = fig.add_subplot(3,2,6)
for i, algo in enumerate(cy):
    ax.plot(cx[cmc:], cy[algo][cmc:], '-', label=algos[i])

ax.legend(loc="upper left", frameon=False, prop={"size":"small", "family":"Arial"})
ax.set_xlim((0,100))
ax.set_ylim((0,11))
ax.set_xticklabels(ax.get_xticklabels(), fontfamily="Arial")
ax.set_yticks([0, 2, 4, 6, 8, 10])
ax.set_yticklabels(["0", "2", "4", "6", "8", "10"], fontfamily="Arial")
ax.axhline(1, ls='--', lw=1, color="gray")
ax.set_xlabel("Percentage of internal nodes interpolated", fontsize="medium", fontfamily="Arial")
ax.set_ylabel("  \n  ")

ymin, ymax = ax.get_ylim()
ax.text(95, ymax*0.95, "B3", fontsize="large", family=["Arial"], fontweight="bold", verticalalignment="top", horizontalalignment="right")

fig.tight_layout()
fig.savefig("figure_interp_compare_narrow.svg")
fig.savefig("figure_interp_compare_narrow.tif", dpi=300)

pylab.close()
