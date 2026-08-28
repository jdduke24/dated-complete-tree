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


import sys
import numpy as np

from dated_complete_tree import tree_loading
from dated_complete_tree import tree_labelling
from dated_complete_tree import tree_pruning
from dated_complete_tree import tree_topology
from dated_complete_tree import tree_dating

import logging
logger = logging.getLogger(__name__)
logging.basicConfig(filename="main.log", filemode="w", force=True, level=logging.ERROR)

sys.setrecursionlimit(10000)

rng = np.random.default_rng(seed=1)

#####################################################################################################################
# Load and prune tree

# Load metadata for tree: Open Tree taxonomy and phylogenetic annotations
phylogeny_nodes, taxa = tree_loading.load_metadata()

# Create ete4 tree structure for entire Open Tree of Life, with my annotations
whole_tre_unmodified = tree_loading.build_and_annotate_tree(phylogeny_nodes, taxa)

tree_pruning.remove_based_on_props(whole_tre_unmodified, ["extinct", "uncultured", "unidentified", "intergeneric", "hybrid"])

tree_pruning.strip_birds(whole_tre_unmodified)
tree_pruning.strip_turtles(whole_tre_unmodified)

tree_pruning.remove_subspecies(whole_tre_unmodified, rng)

tree_pruning.impute_species_into_empty_taxa(whole_tre_unmodified)

#####################################################################################################################
# Fix topology

tree_labelling.add_ancestral_ranks(whole_tre_unmodified)
tree_labelling.add_descendant_ranks(whole_tre_unmodified)

tree_labelling.fix_taxonomy_ordering(whole_tre_unmodified)

tree_topology.forced_taxa_moves(whole_tre_unmodified)

# Copy tree - we will change the copy, and keep the original unchanged so we can restore it next iteration without
# reloading everything
whole_tre = whole_tre_unmodified.copy()

# First, do labelling for steps 1-3:
#  - 1-2 are independent of each other; step 3 collects up nodes not labelled in 1-2.
#  - tree is only labelled at this stage; modifications are made in tree_fixing functions.
genus_dict = {}       # step 1, nodes below genus nodes
nmp_genus_dict = {}   # step 2, non-monophyletic genera
tree_labelling.populate_genus_dict(whole_tre, genus_dict, nmp_genus_dict, None)

tofix_dict = {}       # step 3, all other nodes from taxonomy (not phylogenies) to
                      # be moved to a suitable place in the tree, such that we
                      # generated a plausible hypothetical tree
tree_labelling.populate_tofix_dict(whole_tre, tofix_dict, nmp_genus_dict)

# Second, fix the topology based on the labels.
# Fix steps 1 and 2.
tree_topology.fix_polyphyly(genus_dict, rng)
tree_topology.fix_polyphyly(nmp_genus_dict, rng)

tree_topology.remove_nonspecies_leaves(whole_tre)

# Find and label backbone for step 3, after steps 1 an 2 already fixed.
tree_labelling.populate_tofix_bkb(whole_tre, tofix_dict, [])
fix_dict = tree_labelling.process_tofix_bkb(tofix_dict)

# Finally, fix step 3.
tree_topology.fix_polyphyly(fix_dict, rng, expand_parent_backbones=True)

tree_topology.remove_nonspecies_leaves(whole_tre)

# Last of all, polytomy resolution.
tree_topology.fix_remaining_polytomies(whole_tre, rng)

# Remove one-child nodes. Gives a fully bifurcating tree.
whole_tre = tree_topology.delete_one_child_nodes(whole_tre)

#####################################################################################################################
# Assign and interpolate dates

# Load dates from json
dates = tree_loading.load_dates()

# Assign dates
tree_dating.assign_dates(whole_tre, dates, num_sources=True)


clades = ['cellular_organisms_ott93302',
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



def label_pct_dates(parent):
    if parent.is_leaf:
        results = [0,0,1,0,set()]
    elif not parent.props["date"]:
        results = [1,0,0,0,set()]
    else:
        results = [1,1,0,0,parent.props["date_sourceids"]]

    if parent.props["ph_tx"] == "PH":
        results[3] += 1

    for child in parent.children:
        new_results = label_pct_dates(child)
        results[0] += new_results[0]
        results[1] += new_results[1]
        results[2] += new_results[2]
        results[3] += new_results[3]
        results[4] |= new_results[4]

    parent.add_prop("child_tree_size", results[0])
    parent.add_prop("num_dates", results[1])
    parent.add_prop("num_leaves", results[2])
    parent.add_prop("num_ph", results[3])
    parent.add_prop("num_date_sources", len(results[4]))

    return results

label_pct_dates(whole_tre)

def date_stats(parent, stats):
    if parent.name in clades:
        stats[parent.name] = (parent.props["tx_level"], parent.props["num_dates"], parent.props["num_dates"]/parent.props["child_tree_size"], parent.props["num_leaves"], parent.props["num_ph"]/(parent.props["child_tree_size"]+parent.props["num_leaves"]), parent.props["num_date_sources"])

    for child in parent.children:
        date_stats(child, stats)

stats = {}
date_stats(whole_tre, stats)

fout = open("figures/figure3/date_stats.csv", "wt")
fout.write("Clade,Species richness,Dated nodes,Date coverage,Date source trees,Phylogenetic coverage\n")
for clade in clades:
    parts = clade.split('_')
    if parts[0] == "cellular":
        name = "All Life"
    else:
        name = parts[0]

    # spaces = 36 - tx_levels[stats[clade][0]]
    # fout.write(' ' * spaces)

    fout.write(f'{name},"{stats[clade][3]:,}","{stats[clade][1]:,}",{100*stats[clade][2]:.1f}%,"{stats[clade][5]:,}",{100*stats[clade][4]:.1f}%\n')
fout.close()
