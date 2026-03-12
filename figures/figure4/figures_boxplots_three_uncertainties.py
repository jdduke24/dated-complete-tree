import os
os.chdir('/Users/Jonathan/Library/CloudStorage/Dropbox/Imperial/Tree_of_Life/Open_Tree/python/dated-complete-tree')

import matplotlib.pyplot as plt
import numpy as np

clades = ['cellularorganisms_ott93302',
    'Eukaryota_ott304358',
 'Metazoa_ott691846',
     'Chordata_ott125642',
         'Actinopteri_ott285821',
         'Squamata_ott35888',
         'Aves_ott81461',
         'Amphibia_ott544595',
         'Mammalia_ott244265',
         'Chondrichthyes_ott278108',
         'Testudines_ott639666',
     'Arthropoda_ott632179',
        'Pancrustacea_ott985906',
           'Insecta_ott1062253',
               'Coleoptera_ott865243',
               'Lepidoptera_ott965954',
               'Diptera_ott661378',
               'Hymenoptera_ott753726',
               'Hemiptera_ott603650',
               # 'Orthoptera_ott1095594',
               # 'Trichoptera_ott457402',
               'Odonata_ott133665',
           'Malacostraca_ott212701',
           'Copepoda_ott461528',
           'Branchiopoda_ott632175',
        'Chelicerata_ott1041457',
        'Myriapoda_ott177526',
     'Mollusca_ott802117',
     'Platyhelminthes_ott555379',
     'Annelida_ott941620',
     'Nematoda_ott395057',
     'Cnidaria_ott641033',
     'Porifera_ott67819',
     'Echinodermata_ott451020',
     'Bryozoa_ott442934',
     # 'Rotifera_ott471706',
     'Tardigrada_ott111438',
 'Chloroplastida_ott361838',
     'Tracheophyta_ott10210',
         'Magnoliopsida_ott99252',
         'Polypodiopsida_ott166292',
     'Bryophyta_ott246594',
     'Marchantiophyta_ott56601',
     'Chlorophyta_ott979501',
 'Rhodophyta_ott878953',
 'Fungi_ott352914',
     'Ascomycota_ott439373',
     'Basidiomycota_ott634628',
     'Microsporidia_ott16113',
 'SAR_ott5246039',
     'Stramenopiles_ott266745',
     'Alveolata_ott266751',
     'Rhizaria_ott6929',
'Archaea_ott996421',
'Bacteria_ott844192']

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
     # 2,
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

############
# Data for topological variation

fin = open("figures/figure4/dated_tree_topo_pd_for_clades.txt", 'r')
pd_n = 500

phyla_pd_rel = []
phyla_pd = []
names = []
spp = []

skipped = []
for count, line in enumerate(fin):
    if count == 0:
        # ignore header row
        continue

    parts = line.split('\t')

    nm = parts[0]
    nm = nm[:parts[0].find("_")]
    names.append(nm)
    spp.append(int(parts[1]))

    phyla_pd_rel.append([])
    phyla_pd.append([])
    md = float(parts[7].strip())/1000
    start_idx = 12

    for i in range(start_idx,start_idx+pd_n):
        if len(parts) > i:
            phyla_pd_rel[-1].append((float(parts[i].strip())/1000)/md)
            phyla_pd[-1].append((float(parts[i].strip())/1000))


phyla_pd_rel_copy = phyla_pd_rel.copy()
phyla_pd_copy = phyla_pd.copy()
names_copy = names.copy()
spp_copy = spp.copy()

names_found = set()
for idx, clade in enumerate(clades):
    for i, name in enumerate(names_copy):
        if name in clade:
            phyla_pd_rel[idx] = phyla_pd_rel_copy[i]
            phyla_pd[idx] = int(np.median(phyla_pd_copy[i]))
            names[idx] = names_copy[i]
            spp[idx] = spp_copy[i]
            names_found.add(len(names)-i)

names_diff = len(names_copy) - len(names_found)
for i in range(names_diff):
    del names[-1]
    del phyla_pd_rel[-1]
    del phyla_pd[-1]

phyla_pd.reverse()
phyla_pd_rel.reverse()
names.reverse()
spp.reverse()
indents.reverse()

############
# Data for date variation

fin = open("figures/figure4/dated_tree_date_pd_for_clades.txt", 'r')
pd_n = 50

date_phyla_pd_rel = []
date_phyla_pd = []
dates_names = []

skipped = []
for count, line in enumerate(fin):
    if count == 0:
        # ignore header row
        continue

    parts = line.split('\t')

    nm = parts[0]
    nm = nm[:parts[0].find("_")]
    dates_names.append(nm)

    date_phyla_pd_rel.append([])
    date_phyla_pd.append([])
    md = float(parts[12].strip())/1000
    start_idx = 13

    for i in range(start_idx,start_idx+pd_n):
        if len(parts) > i:
            date_phyla_pd_rel[-1].append(round((float(parts[i].strip())/1000)/md,3))
            date_phyla_pd[-1].append((float(parts[i].strip())/1000))


date_phyla_pd_rel_copy = date_phyla_pd_rel.copy()
phyla_pd_copy = date_phyla_pd.copy()
dates_names_copy = dates_names.copy()

for idx, clade in enumerate(clades):
    for i, name in enumerate(dates_names):
        if name in clade:
            # print(name, np.percentile(phyla_pd_copy[i], 50), np.percentile(phyla_pd_copy[i], 2.5), np.percentile(phyla_pd_copy[i], 97.5))
            date_phyla_pd_rel[idx] = date_phyla_pd_rel_copy[i]
            date_phyla_pd[idx] = int(np.median(phyla_pd_copy[i]))
            dates_names[idx] = dates_names_copy[i]

for i in range(names_diff):
    del dates_names[-1]
    del date_phyla_pd_rel[-1]
    del date_phyla_pd[-1]

date_phyla_pd.reverse()
date_phyla_pd_rel.reverse()


############
# Data for median pds
fin = open("figures/figure4/dated_tree_both_pd_for_clades.txt", 'r')
pd_n = 500

both_phyla_pd_rel = []
both_phyla_pd = []
both_names = []
median_pd = []

skipped = []
for count, line in enumerate(fin):
    if count == 0:
        # ignore header row
        continue

    parts = line.split('\t')

    nm = parts[0]
    nm = nm[:parts[0].find("_")]
    both_names.append(nm)

    both_phyla_pd_rel.append([])
    both_phyla_pd.append([])
    md = float(parts[7].strip())/1000
    start_idx = 12

    for i in range(start_idx,start_idx+pd_n):
        if len(parts) > i:
            both_phyla_pd_rel[-1].append((float(parts[i].strip())/1000)/md)
            both_phyla_pd[-1].append((float(parts[i].strip())/1000))
            median_pd.append(int(round(float(parts[7])/1000, 0)))


both_phyla_pd_rel_copy = both_phyla_pd_rel.copy()
both_phyla_pd_copy = both_phyla_pd.copy()
both_names_copy = both_names.copy()

for idx, clade in enumerate(clades):
    for i, name in enumerate(both_names):
        if name in clade:
            print(name, np.percentile(both_phyla_pd_copy[i], 50), np.percentile(both_phyla_pd_copy[i], 0), np.percentile(both_phyla_pd_copy[i], 100))
            both_phyla_pd_rel[idx] = both_phyla_pd_rel_copy[i]
            both_phyla_pd[idx] = int(np.median(both_phyla_pd_copy[i]))
            both_names[idx] = both_names_copy[i]

for i in range(names_diff):
    del both_names[-1]
    del both_phyla_pd_rel[-1]
    del both_phyla_pd[-1]

both_phyla_pd.reverse()
both_phyla_pd_rel.reverse()
median_pd.reverse()
median_pd.append(0)


names[-1] = "All Life"

names.append("")
spp.append(0)
indents.append(0)


################
# Plotting

fig = plt.figure(figsize=(7.5,7.5))

spec = fig.add_gridspec(1, 11)
ax0 = fig.add_subplot(spec[0,:2])
ax1 = fig.add_subplot(spec[0,2:5])
ax2 = fig.add_subplot(spec[0,5:8])
ax3 = fig.add_subplot(spec[0,8:])

ax1.boxplot(phyla_pd_rel, vert=False, sym=".", boxprops={'lw':0.7}, whiskerprops={'lw':0.7}, flierprops={'markeredgewidth':0.4, 'markersize':4})

ax1.set_xscale("log")
ax1.set_xlabel("Distribution of PD relative to median", size="x-small", family=["Arial"])


ax2.boxplot(date_phyla_pd_rel, vert=False, sym=".", boxprops={'lw':0.7}, whiskerprops={'lw':0.7}, flierprops={'markeredgewidth':0.4, 'markersize':4})#, whis=(0,100))

ax2.set_xscale("log")
ax2.set_xlabel("Distribution of PD relative to baseline", size="x-small", family=["Arial"])

ax3.boxplot(both_phyla_pd_rel, vert=False, sym=".", boxprops={'lw':0.7}, whiskerprops={'lw':0.7}, flierprops={'markeredgewidth':0.4, 'markersize':4})

ax3.set_xscale("log")
ax3.set_xlabel("Distribution of PD relative to median", size="x-small", family=["Arial"])

green = "#DAF2D0"
blue = "#DAE9F8"
grey ="#A6A6A6"
lightgrey = "#E0E0E0"
moregrey = "#888888"

initial_indent = 0.02
indent_size = 0.07

for i in range(len(names)-1):
    ax0.text(initial_indent+indents[i]*indent_size, i+0.95, names[i], horizontalalignment="left", verticalalignment="center", family=["Arial"], fontsize="x-small")
    ax0.plot([0,1],[i+0.5,i+0.5], color=grey, lw=0.5, zorder=1000, clip_on=False)
    ax1.plot([0,2],[i+0.5,i+0.5], color=grey, lw=0.5, zorder=1000)
    ax2.plot([0,2],[i+0.5,i+0.5], color=grey, lw=0.5, zorder=1000)
    ax3.plot([0,2],[i+0.5,i+0.5], color=grey, lw=0.5, zorder=1000)

    if i >= 2 and i < len(names)-3:
        ax0.plot([initial_indent+(indents[i]-1)*indent_size+0.01, initial_indent+(indents[i]-1)*indent_size+0.04], [i+1, i+1], lw=0.4, color=moregrey)

        for j in range(1,indents[i]):
            if ((i > 40 and i < 48) or (i > 13 and i < 16)) and j == 3:
                continue

            if i < 5 and j == 1:
                continue

            if indents[i-1] < indents[i] and j == indents[i]-1 or i == 5:
                ax0.plot([initial_indent+j*indent_size+0.01, initial_indent+j*indent_size+0.01], [i+1, i+1.5], lw=0.4, color=moregrey)
            else:
                ax0.plot([initial_indent+j*indent_size+0.01, initial_indent+j*indent_size+0.01], [i+0.5, i+1.5], lw=0.4, color=moregrey)

# j = 1
# ax0.plot([initial_indent+j*indent_size, initial_indent+j*indent_size], [10, 50.5], lw = 0.5, color=grey, linestyle="dashed")
# ax0.plot([initial_indent+j*indent_size, initial_indent+j*indent_size+0.02], [10, 10], lw = 0.5, color=grey)

i = len(names)-1

ax0.axhline(i+0.5, color=grey, lw=1)
ax0.text(initial_indent,i+1,"Clade", horizontalalignment="left", verticalalignment="center", family=["Arial"], weight="bold", fontsize="x-small")

ax1.axhline(i+0.5, color=grey, lw=1)
ax1.text(1,i+1,"A: Variation in PD from topology", horizontalalignment="center", verticalalignment="center", family=["Arial"], weight="bold", fontsize="x-small")

ax2.axhline(i+0.5, color=grey, lw=1)
ax2.text(1,i+1,"B: Variation in PD from date sources", horizontalalignment="center", verticalalignment="center", family=["Arial"], weight="bold", fontsize="x-small")

ax3.axhline(i+0.5, color=grey, lw=1)
ax3.text(1,i+1,"C: Variation in PD from both", horizontalalignment="center", verticalalignment="center", family=["Arial"], weight="bold", fontsize="x-small")

ax0.set_xlim(0,1)

ax1.set_xticks([0.5,0.6,0.7,0.8,0.9, 1, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 2,3], labels=['0.5','0.6','0.7','0.8','0.9','1','', '1.2', '', '1.4', '', '1.6', '2','3'], family=["Arial"], size="x-small")
ax1.set_xlim(0.69,1.46)

ax2.set_xticks([0.5,0.6,0.7,0.8,0.9, 1, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 2,3], labels=['0.5','0.6','0.7','0.8','0.9','1','', '1.2', '', '1.4', '', '1.6', '2','3'], family=["Arial"], size="x-small")
ax2.set_xlim(0.69,1.46)

ax3.set_xticks([0.5,0.6,0.7,0.8,0.9, 1, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 2,3], labels=['0.5','0.6','0.7','0.8','0.9','1','', '1.2', '', '1.4', '', '1.6', '2','3'], family=["Arial"], size="x-small")
ax3.set_xlim(0.69,1.46)

ax0.get_yaxis().set_ticks([])
ax0.get_xaxis().set_ticks([])

ax1.get_yaxis().set_ticks([])
ax2.get_yaxis().set_ticks([])
ax3.get_yaxis().set_ticks([])

ax0.axhspan(i-0.5, i+0.5, color=blue, alpha=1, zorder=-1000)
ax1.axhspan(i-0.5, i+0.5, color=blue, alpha=1, zorder=-1000)
ax2.axhspan(i-0.5, i+0.5, color=blue, alpha=1, zorder=-1000)
ax3.axhspan(i-0.5, i+0.5, color=blue, alpha=1, zorder=-1000)

for i in range(len(names)-1):
    if i % 2 == 1:
        ax0.axhspan(i+1.5, i+0.5, color=lightgrey, alpha=0.35, zorder=-1000)
        ax1.axhspan(i+1.5, i+0.5, color=lightgrey, alpha=0.35, zorder=-1000)
        ax2.axhspan(i+1.5, i+0.5, color=lightgrey, alpha=0.35, zorder=-1000)
        ax3.axhspan(i+1.5, i+0.5, color=lightgrey, alpha=0.35, zorder=-1000)

# ax0.axhspan(len(names)-3.5, len(names)-11.5, color=blue, alpha=1)
# ax1.axhspan(len(names)-3.5, len(names)-11.5, color=blue, alpha=1)

# ax0.axhspan(len(names)-36.5, len(names)-40.5, color=green, alpha=1)
# ax1.axhspan(len(names)-36.5, len(names)-40.5, color=green, alpha=1)

ax0.set_frame_on(False)
ax1.set_frame_on(False)
ax2.set_frame_on(False)
ax3.set_frame_on(False)

ax0.set_ylim(ax1.get_ylim())
ax2.set_ylim(ax1.get_ylim())
ax3.set_ylim(ax1.get_ylim())

fig.tight_layout()

import PIL
fig.savefig("figures/figure4/boxplot_three_uncertainties_1500.pdf")
fig.savefig("figures/figure4/boxplot_three_uncertainties_1500.tif", dpi=300)
