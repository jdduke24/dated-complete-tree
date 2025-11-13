import os
os.chdir('/Users/Jonathan/Library/CloudStorage/Dropbox/Imperial/Tree_of_Life/Open_Tree/python/dated-complete-tree/figures/figure4')

import matplotlib.pyplot as plt
import numpy as np

clades = ['Eukaryota_ott304358',
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
               'Orthoptera_ott1095594',
               'Trichoptera_ott457402',
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
     'Rotifera_ott471706',
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

indents = [0,
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
               5,
               5,
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


fin = open("pd_by_clade_from_newicks.txt", 'r')

phyla_pd_rel = []
phyla_pd = []
names = []
spp = []

pd_n = 200

skipped = []
count = 0
for line in fin:
    if count == 0:
        count += 1
        continue
    parts = line.split('\t')

    nm = parts[0]
    nm = nm[:parts[0].find("_")]
    names.append(nm)
    spp.append(int(parts[1]))

    phyla_pd_rel.append([])
    phyla_pd.append([])
    md = float(parts[7].strip())/1000
    for i in range(12,12+pd_n):
        if len(parts) > i:
            phyla_pd_rel[-1].append((float(parts[i].strip())/1000)/md)
            phyla_pd[-1].append((float(parts[i].strip())/1000))

    count += 1

phyla_pd_rel_copy = phyla_pd_rel.copy()
phyla_pd_copy = phyla_pd.copy()
names_copy = names.copy()
spp_copy = spp.copy()

for idx, clade in enumerate(clades):
    for i, name in enumerate(names_copy):
        if 'Eukaryota' in name:
            euk_pd_rel = phyla_pd_rel_copy[i]
            euk_pd = phyla_pd_copy[i]
        elif name in clade:
            print(name, np.percentile(phyla_pd_copy[i], 50), np.percentile(phyla_pd_copy[i], 2.5), np.percentile(phyla_pd_copy[i], 97.5))
            phyla_pd_rel[idx] = phyla_pd_rel_copy[i]
            phyla_pd[idx] = int(np.median(phyla_pd_copy[i]))
            names[idx] = names_copy[i]
            spp[idx] = spp_copy[i]


del phyla_pd_rel[0]
del names[0]
del spp[0]
del phyla_pd[0]

phyla_pd.reverse()
phyla_pd_rel.reverse()
names.reverse()
spp.reverse()
indents.reverse()

phyla_pd.append(int(np.median(euk_pd)))
phyla_pd_rel.append(euk_pd_rel)
names.append("All Eukaryota")
spp.append(2222870)

phyla_pd_rel.append([])
names.append("")
spp.append(0)


fig = plt.figure(figsize=(7.3,12))

# ax = fig.add_subplot(111)

spec = fig.add_gridspec(1, 5)
ax0 = fig.add_subplot(spec[0,:3])

ax1 = fig.add_subplot(spec[0,3:])

ax1.boxplot(phyla_pd_rel, vert=False, sym=".")

ax1.set_xscale("log")
ax1.set_xlabel("Distribution of PD relative to median", size="small", family=["Arial"])

ax0.set_ylim(ax1.get_ylim())

green = "#DAF2D0"
blue = "#DAE9F8"
grey ="#A6A6A6"

for i in range(len(names)-1):
    ax0.text(0.02+indents[i]*0.02, i+1, names[i], horizontalalignment="left", verticalalignment="center", family=["Arial"], fontsize="small")
    ax0.text(0.62,i+1,format(spp[i], ",d"), horizontalalignment="right", verticalalignment="center", family=["Arial"], fontsize="small")
    ax0.text(0.98,i+1,format(phyla_pd[i], ",d"), horizontalalignment="right", verticalalignment="center", family=["Arial"], fontsize="small")
    ax0.axhline(i+0.5, color=grey, lw=1)
    ax1.axhline(i+0.5, color=grey, lw=1)

i = len(names)-1
ax0.text(0.02,i+1,"Phylum", horizontalalignment="left", verticalalignment="center", family=["Arial"], weight="bold", fontsize="small")
ax0.text(0.62,i+1,"Species richness", horizontalalignment="right", verticalalignment="center", family=["Arial"], weight="bold", fontsize="small")
ax0.text(0.98,i+1,"Median PD (Ga)", horizontalalignment="right", verticalalignment="center", family=["Arial"], weight="bold", fontsize="small")
ax0.axhline(i+0.5, color=grey, lw=1)

ax1.axhline(i+0.5, color=grey, lw=1)
ax1.text(1,i+1,"Topological variation", horizontalalignment="center", verticalalignment="center", family=["Arial"], weight="bold", fontsize="small")

ax1.set_xticks([0.5,0.6,0.7,0.8,0.9, 1, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 2,3], labels=['0.5','0.6','0.7','0.8','0.9','1','1.1', '1.2', '1.3', '1.4', '1.5', '1.6', '1.7', '2','3'], family=["Arial"], size="small")
ax1.set_xlim(0.7,1.43)


ax0.get_yaxis().set_ticks([])
ax0.get_xaxis().set_ticks([])
ax1.get_yaxis().set_ticks([])


ax0.axhspan(len(names)-2.5, len(names)-10.5, color=blue, alpha=1)
ax1.axhspan(len(names)-2.5, len(names)-10.5, color=blue, alpha=1)

ax0.axhspan(len(names)-36.5, len(names)-40.5, color=green, alpha=1)
ax1.axhspan(len(names)-36.5, len(names)-40.5, color=green, alpha=1)

ax0.set_frame_on(False)
ax1.set_frame_on(False)


import PIL
fig.savefig("boxplot_figure_3.pdf")
fig.savefig("boxplot_figure_3.tif", dpi=300, bbox_inches="tight", pad_inches=0)
