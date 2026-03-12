import os
os.chdir('/Users/Jonathan/Library/CloudStorage/Dropbox/Imperial/Tree_of_Life/Open_Tree/python/dated-complete-tree')

import sys
import numpy as np

import tree_loading
import tree_labelling
import tree_fixing
import tree_dating
import tree_metrics
import tree_checks

import logging
logger = logging.getLogger(__name__)
logging.basicConfig(filename="main.log", filemode="w", force=True, level=logging.WARNING)

sys.setrecursionlimit(10000)

#####################################################################################################################
# Load and prune tree

# Load metadata for tree from Open Tree and Chronosynth
dates, phylogeny_nodes, taxa = tree_loading.load_metadata()

# Create ETE3 tree structure for entire Open Tree of Life, with my annotations
whole_tre_unmodified = tree_loading.build_and_annotate_tree(phylogeny_nodes, taxa)

tree_fixing.strip_birds(whole_tre_unmodified)
tree_fixing.strip_turtles(whole_tre_unmodified)

rng = np.random.default_rng(seed=1)

tree_fixing.remove_subspecies(whole_tre_unmodified, rng)
tree_fixing.impute_species_into_empty_taxa(whole_tre_unmodified)

tree_fixing.fix_taxonomy_ordering(whole_tre_unmodified)

tree_labelling.add_anc_ranks(whole_tre_unmodified)
tree_labelling.add_desc_ranks(whole_tre_unmodified)

tree_fixing.forced_taxa_moves(whole_tre_unmodified)

# Copy tree - we will change the copy, and keep the original unchanged so we can restore it next iteration without
# reloading everything
whole_tre = whole_tre_unmodified.copy()

#####################################################################################################################
# Fix topology

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
tree_fixing.fix_polyphyly(genus_dict, rng)
tree_fixing.fix_polyphyly(nmp_genus_dict, rng)

tree_fixing.remove_nonspecies_leaves(whole_tre)

# Find and label backbone for step 3, after steps 1 an 2 already fixed.
tree_labelling.populate_tofix_bkb(whole_tre, tofix_dict, [])
fix_dict = tree_labelling.process_tofix_bkb(tofix_dict)

# Finally, fix step 3.
tree_fixing.fix_polyphyly(fix_dict, rng, expand_parent_backbones=True)

tree_fixing.remove_nonspecies_leaves(whole_tre)

# Last of all, polytomy resolution.
tree_fixing.fix_all_polytomies(whole_tre, rng)

# Remove one-child nodes. Gives a fully bifurcating tree.
whole_tre = tree_fixing.delete_one_child_nodes(whole_tre)

#####################################################################################################################
# Assign and interpolate dates

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


from taxonomy_utils import tx_levels

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

fout = open("date_stats_14.txt", "wt")
for clade in clades:
    spaces = 36 - tx_levels[stats[clade][0]]
    fout.write(' ' * spaces)
    fout.write("%s\t%s\t%d\t%f\t%d\t%d\t%f\n" % (clade, stats[clade][0], stats[clade][1], stats[clade][2], stats[clade][3], stats[clade][5], stats[clade][4]))
fout.close()
