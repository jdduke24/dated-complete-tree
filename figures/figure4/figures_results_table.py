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

phyla_pd_med = []
phyla_pd_low = []
phyla_pd_upp = []
names = []

skipped = []
for count, line in enumerate(fin):
    if count == 0:
        # ignore header row
        continue

    parts = line.split('\t')

    nm = parts[0]
    nm = nm[:parts[0].find("_")]
    names.append(nm)

    phyla_pd_med.append((float(parts[7].strip())/1000))
    phyla_pd_low.append((float(parts[4].strip())/1000))
    phyla_pd_upp.append((float(parts[10].strip())/1000))


phyla_pd_med_copy = phyla_pd_med.copy()
phyla_pd_low_copy = phyla_pd_low.copy()
phyla_pd_upp_copy = phyla_pd_upp.copy()
names_copy = names.copy()

names_found = set()
for idx, clade in enumerate(clades):
    for i, name in enumerate(names_copy):
        if name in clade:
            phyla_pd_med[idx] = phyla_pd_med_copy[i]
            phyla_pd_low[idx] = phyla_pd_low_copy[i]
            phyla_pd_upp[idx] = phyla_pd_upp_copy[i]
            names[idx] = names_copy[i]
            names_found.add(i)

names_diff = len(names_copy) - len(names_found)
for i in range(names_diff):
    del names[-1]
    del phyla_pd_med[-1]
    del phyla_pd_low[-1]
    del phyla_pd_upp[-1]

phyla_pd_med.reverse()
phyla_pd_low.reverse()
phyla_pd_upp.reverse()

############
# Data for date variation

fin = open("figures/figure4/dated_tree_date_pd_for_clades.txt", 'r')

dates_phyla_pd_med = []
dates_phyla_pd_low = []
dates_phyla_pd_upp = []
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

    dates_phyla_pd_med.append((float(parts[7].strip())/1000))
    dates_phyla_pd_low.append((float(parts[4].strip())/1000))
    dates_phyla_pd_upp.append((float(parts[10].strip())/1000))


dates_phyla_pd_med_copy = dates_phyla_pd_med.copy()
dates_phyla_pd_low_copy = dates_phyla_pd_low.copy()
dates_phyla_pd_upp_copy = dates_phyla_pd_upp.copy()
dates_names_copy = dates_names.copy()

n = set()
for idx, clade in enumerate(clades):
    for i, name in enumerate(dates_names):
        if name in clade:
            dates_phyla_pd_med[idx] = dates_phyla_pd_med_copy[i]
            dates_phyla_pd_low[idx] = dates_phyla_pd_low_copy[i]
            dates_phyla_pd_upp[idx] = dates_phyla_pd_upp_copy[i]
            dates_names[idx] = dates_names_copy[i]

for i in range(names_diff):
    del dates_names[-1]
    del dates_phyla_pd_med[-1]
    del dates_phyla_pd_low[-1]
    del dates_phyla_pd_upp[-1]

dates_phyla_pd_med.reverse()
dates_phyla_pd_low.reverse()
dates_phyla_pd_upp.reverse()


############
# Data for median pds
fin = open("figures/figure4/dated_tree_both_pd_for_clades.txt", 'r')

both_phyla_pd_med = []
both_phyla_pd_low = []
both_phyla_pd_upp = []
both_names = []

skipped = []
for count, line in enumerate(fin):
    if count == 0:
        # ignore header row
        continue

    parts = line.split('\t')

    nm = parts[0]
    nm = nm[:parts[0].find("_")]
    both_names.append(nm)

    both_phyla_pd_med.append((float(parts[7].strip())/1000))
    both_phyla_pd_low.append((float(parts[4].strip())/1000))
    both_phyla_pd_upp.append((float(parts[10].strip())/1000))


both_phyla_pd_med_copy = both_phyla_pd_med.copy()
both_phyla_pd_low_copy = both_phyla_pd_low.copy()
both_phyla_pd_upp_copy = both_phyla_pd_upp.copy()
both_names_copy = both_names.copy()

n = set()
for idx, clade in enumerate(clades):
    for i, name in enumerate(both_names):
        if name in clade:
            both_phyla_pd_med[idx] = both_phyla_pd_med_copy[i]
            both_phyla_pd_low[idx] = both_phyla_pd_low_copy[i]
            both_phyla_pd_upp[idx] = both_phyla_pd_upp_copy[i]
            both_names[idx] = both_names_copy[i]

for i in range(names_diff):
    del both_names[-1]
    del both_phyla_pd_med[-1]
    del both_phyla_pd_low[-1]
    del both_phyla_pd_upp[-1]

both_phyla_pd_med.reverse()
both_phyla_pd_low.reverse()
both_phyla_pd_upp.reverse()

################

names.reverse()
indents.reverse()

names[-1] = "All Life"

names.append("")
indents.append(0)


################
# Plotting

fig = plt.figure(figsize=(6.5,7.5))

ax0 = fig.add_subplot(1,1,1)

green = "#DAF2D0"
blue = "#DAE9F8"
grey ="#A6A6A6"
lightgrey ="#E0E0E0"
moregrey = "#888888"

col_x = [0.007,
         0.28,
         0.38,
         0.54,
         0.64,
         0.8,
         0.9]

initial_indent = col_x[0]
indent_size = 0.014

for i in range(len(names)-1):
    # ax0.text(0.02+indents[i]*0.02, i+1, names[i], horizontalalignment="left", verticalalignment="center", family=["Arial"], fontsize="x-small")
    ax0.text(col_x[0]+indents[i]*indent_size, i+0.95, names[i], horizontalalignment="left", verticalalignment="center", family=["Arial"], fontsize="x-small")
    ax0.text(col_x[1], i+0.95, "%d" % round(phyla_pd_med[i],0), horizontalalignment="right", verticalalignment="center", family=["Arial"], fontsize="x-small")
    ax0.text(col_x[2], i+0.95, "%d - %d" % (round(phyla_pd_low[i],0), round(phyla_pd_upp[i],0)), horizontalalignment="center", verticalalignment="center", family=["Arial"], fontsize="x-small")

    ax0.text(col_x[3], i+0.95, "%d" % round(dates_phyla_pd_med[i],0), horizontalalignment="right", verticalalignment="center", family=["Arial"], fontsize="x-small")
    ax0.text(col_x[4], i+0.95, "%d - %d" % (round(dates_phyla_pd_low[i],0), round(dates_phyla_pd_upp[i],0)), horizontalalignment="center", verticalalignment="center", family=["Arial"], fontsize="x-small")

    ax0.text(col_x[5], i+0.95, "%d" % round(both_phyla_pd_med[i],0), horizontalalignment="right", verticalalignment="center", family=["Arial"], fontsize="x-small")
    ax0.text(col_x[6], i+0.95, "%d - %d" % (round(both_phyla_pd_low[i],0), round(both_phyla_pd_upp[i],0)), horizontalalignment="center", verticalalignment="center", family=["Arial"], fontsize="x-small")

    ax0.axhline(i+0.5, color=grey, lw=0.5, zorder=1000)

    if i >= 2 and i < len(names)-3:
        ax0.plot([initial_indent+(indents[i]-1)*indent_size+0.0015, initial_indent+(indents[i]-1)*indent_size+0.007], [i+1, i+1], lw=0.4, color=moregrey)

        for j in range(1,indents[i]):
            if ((i > 40 and i < 48) or (i > 13 and i < 16)) and j == 3:
                continue

            if i < 5 and j == 1:
                continue

            if indents[i-1] < indents[i] and j == indents[i]-1 or i == 5:
                ax0.plot([initial_indent+j*indent_size+0.0015, initial_indent+j*indent_size+0.0015], [i+1, i+1.5], lw=0.4, color=moregrey)
            else:
                ax0.plot([initial_indent+j*indent_size+0.0015, initial_indent+j*indent_size+0.0015], [i+0.5, i+1.5], lw=0.4, color=moregrey)


i = len(names)-1
ax0.text(col_x[0],i+0.65,"Clade", horizontalalignment="left", verticalalignment="bottom", family=["Arial"], weight="bold", fontsize="x-small")
ax0.text(col_x[1]+0.03,i+0.65,"Median PD (Byr)", horizontalalignment="right", verticalalignment="bottom", family=["Arial"], weight="bold", fontsize="x-small")
ax0.text(col_x[2],i+0.65,"95% CI", horizontalalignment="center", verticalalignment="bottom", family=["Arial"], weight="bold", fontsize="x-small")

ax0.text(col_x[1]+(col_x[2]-col_x[1])/4,i+1.65,"A: Variation from topology", horizontalalignment="center", verticalalignment="bottom", family=["Arial"], weight="bold", fontsize="x-small")

ax0.text(col_x[3]+0.03,i+0.65,"Median PD (Byr)", horizontalalignment="right", verticalalignment="bottom", family=["Arial"], weight="bold", fontsize="x-small")
ax0.text(col_x[4],i+0.65,"95% CI", horizontalalignment="center", verticalalignment="bottom", family=["Arial"], weight="bold", fontsize="x-small")

ax0.text(col_x[3]+(col_x[4]-col_x[3])/4,i+1.65,"B: Variation from date sources", horizontalalignment="center", verticalalignment="bottom", family=["Arial"], weight="bold", fontsize="x-small")

ax0.text(col_x[5]+0.03,i+0.65,"Median PD (Byr)", horizontalalignment="right", verticalalignment="bottom", family=["Arial"], weight="bold", fontsize="x-small")
ax0.text(col_x[6],i+0.65,"95% CI", horizontalalignment="center", verticalalignment="bottom", family=["Arial"], weight="bold", fontsize="x-small")

ax0.text(col_x[5]+(col_x[6]-col_x[5])/4,i+1.65,"C: Variation from both", horizontalalignment="center", verticalalignment="bottom", family=["Arial"], weight="bold", fontsize="x-small")

ax0.axhline(i+0.5, color=grey, lw=1)

ax0.get_yaxis().set_ticks([])
ax0.get_xaxis().set_ticks([])

ax0.axhspan(len(names)-1.5, len(names)-0.5, color=blue, alpha=1, zorder=-1000)

for i in range(len(names)-1):
    if i % 2 == 1:
        ax0.axhspan(i+1.5, i+0.5, color=lightgrey, alpha=0.35, zorder=-1000)

ax0.set_xlim(0,1)
ax0.set_ylim(0,len(names)+1.5)

ax0.set_frame_on(False)

import PIL
fig.tight_layout()
fig.savefig("figures/figure4/pd_table_tight.pdf")
fig.savefig("figures/figure4/pd_table_tight.tif", dpi=300)
