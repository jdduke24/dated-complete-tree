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

fig = pylab.figure(figsize=(9,9))

ax = fig.add_subplot(3,2,1)
for i, algo in enumerate(fy):
    ax.plot(fx[fmc:], fy[algo][fmc:], '-', label=algos[i])

ax.legend(loc="upper left", frameon=False)
ax.set_xlim((0,100))
ax.set_ylim((0,18))
ax.set_yticks([0, 2, 4, 6, 8, 10, 12, 14, 16, 18])
ax.set_yticklabels(["0%", "2%", "4%", "6%", "8%", "10%", "12%", "14%", "16%", "18%"])
ax.set_ylabel("Date reproduction error:\n$abs$(orig. date - interp. date) / crown date", fontsize="small")

ax.text(-25, 17, "A", fontsize="x-large", family=["Arial"], fontweight="bold")

low = 65
high = 65.3

ax = fig.add_subplot(3,2,2)
for i, algo in enumerate(cy):
    ax.plot(cx[cmc:], cy[algo][cmc:], '-', label=algos[i])

ax.legend(loc="upper left", frameon=False)
ax.set_xlim((0,100))
ax.set_ylim((0,18))
ax.set_yticks([0, 2, 4, 6, 8, 10, 12, 14, 16, 18])
ax.set_yticklabels(["0%", "2%", "4%", "6%", "8%", "10%", "12%", "14%", "16%", "18%"])
ax.set_ylabel("  \n  ")

ax.text(-20, 17, "B", fontsize="x-large", family=["Arial"], fontweight="bold")


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
ax.legend(loc="upper left", frameon=False)
ax.set_xlim((0,100))
ax.set_ylim((0,7))
ax.set_yticks([0, 1, 2, 3, 4, 5, 6, 7])
ax.set_yticklabels(["0", "1", "2", "3", "4", "5", "6", "7"])
ax.axhline(1, ls='--', lw=1, color="gray")
# ax.set_xlabel("Percentage of internal nodes interpolated")
ax.set_ylabel("Mean pendant edge length relative to EQS-LS\n(Myr divided by initial value for EQS-LS)", fontsize="small")

ax.text(-25, 6.5, "C", fontsize="x-large", family=["Arial"], fontweight="bold")
# ax.axvspan(0, low, alpha=0.08, color="blue")
# ax.axvspan(high, 100, alpha=0.15, color="orange")


ax = fig.add_subplot(3,2,4)
for i, algo in enumerate(cy):
    ax.plot(cx[cmc:], cy[algo][cmc:], '-', label=algos[i])

ax.legend(loc="upper left", frameon=False)
ax.set_xlim((0,100))
ax.set_ylim((0,7))
ax.set_yticks([0, 1, 2, 3, 4, 5, 6, 7])
ax.set_yticklabels(["0", "1", "2", "3", "4", "5", "6", "7"])
ax.axhline(1, ls='--', lw=1, color="gray")
ax.set_ylabel("  \n  ")

ax.text(-16, 6.5, "D", fontsize="x-large", family=["Arial"], fontweight="bold")


########

# Pen medians

with open('fullydated_penmeds.pickle', 'rb') as f:
    fx, fy = pickle.load(f)

with open('classes_penmeds.pickle', 'rb') as f:
    cx, cy = pickle.load(f)

ax = fig.add_subplot(3,2,5)
for i, algo in enumerate(fy):
    ax.plot(fx[fmc:], fy[algo][fmc:], '-', label=algos[i])

ax.legend(loc="upper left", frameon=False)
ax.set_xlim((0,100))
ax.set_ylim((0,11))
ax.set_yticks([0, 2, 4, 6, 8, 10])
ax.set_yticklabels(["0", "2", "4", "6", "8", "10"])
ax.axhline(1, ls='--', lw=1, color="gray")
ax.set_xlabel("Percentage of internal nodes interpolated")
ax.set_ylabel("Median pendant edge length relative to EQS-LS\n(Myr divided by initial value for EQS-LS)", fontsize="small")

ax.text(-25, 10, "E", fontsize="x-large", family=["Arial"], fontweight="bold")


ax = fig.add_subplot(3,2,6)
for i, algo in enumerate(cy):
    ax.plot(cx[cmc:], cy[algo][cmc:], '-', label=algos[i])

ax.legend(loc="upper left", frameon=False)
ax.set_xlim((0,100))
ax.set_ylim((0,11))
ax.set_yticks([0, 2, 4, 6, 8, 10])
ax.set_yticklabels(["0", "2", "4", "6", "8", "10"])
ax.axhline(1, ls='--', lw=1, color="gray")
ax.set_xlabel("Percentage of internal nodes interpolated")
ax.set_ylabel("  \n  ")

ax.text(-16, 10, "F", fontsize="x-large", family=["Arial"], fontweight="bold")


fig.tight_layout()
fig.savefig("figure_interp_compare.svg")
fig.savefig("figure_interp_compare.tif", dpi=300)

pylab.close()
