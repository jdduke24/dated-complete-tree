# BSD 3-Clause License

# Copyright (c) 2025, Jonathan David Duke

# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:

# 1. Redistributions of source code must retain the above copyright notice, this
#    list of conditions and the following disclaimer.

# 2. Redistributions in binary form must reproduce the above copyright notice,
#    this list of conditions and the following disclaimer in the documentation
#    and/or other materials provided with the distribution.

# 3. Neither the name of the copyright holder nor the names of its
#    contributors may be used to endorse or promote products derived from
#    this software without specific prior written permission.

# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
# FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
# DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
# SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
# CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
# OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


import matplotlib.pyplot as plt


# read in comparison times (see benchmarking.py) from a csv and make the figure from them

fin = open("figures/figure6/bladj vs me.csv", 'r')

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

fig = plt.figure(figsize=(5.2,4.5))
ax1 = fig.add_subplot(111)

ax1.scatter(spp[:-3], bladj[:-3], label="Phylocom EQS-L implementation (BLADJ)")
ax1.scatter(spp[-3:], bladj[-3:], marker='o', facecolors="none", edgecolors="#1f77b4", label="Phylocom projected run times")
ax1.scatter(spp, me, marker="^", label="Our EQS-L implementation")
ax1.scatter(spp, me_new, marker="+", label="Our EQS-L with improved data cleaning", s=60)


ax1.plot([spp[0], spp[-1]], [bladj[0]*1.8, bladj[0]*1.8*spp[-1]**2/spp[0]**2], ':', label="Quadratic growth trend", zorder=-100)
ax1.plot([spp[0], spp[-1]], [me[0]*0.85, me[0]*0.85*spp[-1]/spp[0]], '--', label="Linear growth trend", zorder=-100, lw=1.2)
ax1.plot([spp[0], spp[-1]], [me_new[0], me_new[0]*spp[-1]/spp[0]], ls='dashdot', label="Linear growth trend", zorder=-100, lw=1)

ax1.set_xscale("log")
ax1.set_yscale("log")

ax1.set_xticklabels(labels=ax1.get_xticklabels(), fontfamily="Arial", fontsize="small")
ax1.set_yticklabels(labels=ax1.get_yticklabels(), fontfamily="Arial", fontsize="small")

ax1.set_xlabel("Number of nodes (log scale)", fontfamily="Arial", fontsize="small")
ax1.set_ylabel("Run time (minutes, log scale)", fontfamily="Arial", fontsize="small")

ax1.legend(loc="upper left", frameon=False, prop={"family":"Arial", "size":"small"})

ax1.text(spp[-1]*1.2,bladj[-1]*0.98,"14 days", color="black", verticalalignment="center", size="x-small", fontfamily="Arial")
ax1.text(spp[-1]*1.2,me[-1]*0.98,"26 seconds", color="black", verticalalignment="center", size="x-small", fontfamily="Arial")
ax1.text(spp[-1]*1.2,me_new[-1]*0.98,"51 seconds", color="black", verticalalignment="center", size="x-small", fontfamily="Arial")


xmin, xmax = ax1.get_xlim()
ax1.set_xlim(xmin, xmax*4)

ymin, ymax = ax1.get_ylim()

ax1.plot([spp[-1], spp[-1]], [ymin, ymin*10], lw=0.5, color="#888888", zorder=-200)
ax1.text(spp[-1]*1.06,2.3e-4,"Complete\ntree", color="#888888", verticalalignment="top", size="x-small", fontfamily="Arial")

ax1.plot([spp[-3], spp[-3]], [ymin, ymin*10], lw=0.5, color="#888888", zorder=-200)
ax1.text(spp[-3]*1.06,2.3e-4,"Flowering\nplants", color="#888888", verticalalignment="top", size="x-small", fontfamily="Arial")

ax1.plot([spp[-4], spp[-4]], [ymin, ymin*10], lw=0.5, color="#888888", zorder=-200)
ax1.text(spp[-4]*1.06,2.3e-4,"Gastropods", color="#888888", verticalalignment="top", size="x-small", fontfamily="Arial")

ax1.plot([spp[-7], spp[-7]], [ymin, ymin*10], lw=0.5, color="#888888", zorder=-200)
ax1.text(spp[-7]*1.06,2.3e-4,"Squamates", color="#888888", verticalalignment="top", size="x-small", fontfamily="Arial")

ax1.set_ylim(ymin, ymax)

fig.tight_layout()
fig.savefig("figures/figure6/bladj_v_me_updated.pdf")
fig.savefig("figures/figure6/Fig6.tif", dpi=300)

plt.close()
