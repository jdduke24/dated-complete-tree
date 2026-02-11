import os
os.chdir('/Users/Jonathan/Library/CloudStorage/Dropbox/Imperial/Tree_of_Life/Open_Tree/python/dated-complete-tree')

import sys
import gc

import tree_loading
import tree_labelling
import tree_fixing
import tree_dating
import tree_checks

from taxonomy_utils import tx_levels
import tree_plotting
import tree_metrics

import random
import numpy as np
from copy import deepcopy
from copy import copy


import logging
logger = logging.getLogger(__name__)
logging.basicConfig(filename="main.log", filemode="w", force=True, level=logging.ERROR)

sys.setrecursionlimit(10000)

rng = np.random.default_rng(seed=1)

# Load metadata for tree from Open Tree, Chronosynth and OneZoom
dates, phylogeny_nodes, taxa = tree_loading.load_metadata()

# Create ETE3 tree structure for entire Open Tree of Life, with my annotations
whole_tre_unmodified = tree_loading.build_and_annotate_tree(phylogeny_nodes, taxa)

tree_fixing.strip_birds(whole_tre_unmodified)
tree_fixing.strip_turtles(whole_tre_unmodified)

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

whole_tre = tree_fixing.delete_one_child_nodes(whole_tre)


# tree_dating.remove_inconsistent_dates(whole_tre, whole_tre.props["date"]+1)
tree_dating.assign_dates(whole_tre, dates)

count = 0
for node in whole_tre.traverse(strategy="preorder"):
    if node.props["date"] and node.props["date"] > 0:
        count += 1

# tree_dating.remove_inconsistent_dates(whole_tre)

# count = 0
# for node in whole_tre.traverse(strategy="preorder"):
#     if node.props["date"] and node.props["date"] > 0:
#         count += 1

# Date cleaning to ensure time consistency down the tree
tree_dating.label_older_descendants(whole_tre)
tree_dating.dq_date_removal(whole_tre)

count = 0
for node in whole_tre.traverse(strategy="preorder"):
    if node.props["date"] and node.props["date"] > 0:
        count += 1

def label_pct_dates(parent):
    if parent.is_leaf:
        results = [0,0,1]
    elif not parent.props["date"]:
        results = [1,0,0]
    else:
        results = [1,1,0]

    for child in parent.children:
        new_results = label_pct_dates(child)
        results[0] += new_results[0]
        results[1] += new_results[1]
        results[2] += new_results[2]

    parent.add_prop("child_tree_size", results[0])
    parent.add_prop("num_dates", results[1])
    parent.add_prop("num_leaves", results[2])

    return results

label_pct_dates(whole_tre)


clades = set(['Fagales_ott267709',
              'Mammalia_ott244265',
              'Aves_ott81461',
              'Squamata_ott35888',
              'Amphibia_ott544595',
              'Actinopteri_ott285821',
              'Polypodiales_ott93178'])

# date_pcts = [0.25, 0.5, 0.75] # pct to throw away
date_pcts = [0.2, 0.5, 0.75, 0.88, 0.95, 0.98, 1] # pct to throw away
# date_pcts = [0.25]

# params = [[1, 0, False, "Long", False],
#           [0, 0, False, "Short", False],
#           [1, 0, True, "LnN", False],
#           [1, 0, False, "Birth", True]]

params = [[1, 0, False, "Long", False],
          [0, 0, False, "Short", False],
          [1, 0, True, "LnN", False],
          [1, 0, False, "Birth", True],
          [0.25, 0, False, "LS", False]]

# params = [[1, 0, False, "Birth", True]]

repeats = 10

all_results = {}

rng_main = np.random.default_rng(seed=100)
dating_rng = np.random.default_rng(seed=100)

for node in whole_tre.traverse(strategy="preorder"):

    if node.name not in clades:
        continue

    clade = node.name
    clade_tre_orig = node.copy()

    root_date = clade_tre_orig.props["date"]

    # dists = []
    # for leaf in whole_tre.leaves():
    #     dists.append(leaf.dist)
    # mean_pen = np.mean(dists)
    # med_pen = np.median(dists)

    for date_pct in date_pcts:
        for rpt in range(repeats):

            diffs = [[] for p in params]
            lindiffs = [[] for p in params]
            pendants = [[] for p in params]
            gammas = []

            seed = rng_main.integers(1000)

            for i, param in enumerate(params):
                clade_tre = clade_tre_orig.copy()

                rng = np.random.default_rng(seed=seed)

                count = 0
                for node in clade_tre.traverse(strategy="preorder"):
                    node.add_prop("orig_date", None)
                    if node is clade_tre:
                        continue
                    if node.props["date"]:
                        # if node.props["date"] has a value and the value is > 0:
                        if rng.random() < date_pct:
                            node.props["orig_date"] = copy(node.props["date"])
                            node.props["date"] = None
                            count += 1

                print(clade, clade_tre.props["date"], param[3], count)

                tree_dating.date_labelling(clade_tre)

                tree_dating.impute_missing_dates(clade_tre, l=param[0], m=param[1], useLnN=param[2], useBirth=param[4], rng=dating_rng)
                tree_dating.compute_branch_lengths(clade_tre)

                for node in clade_tre.traverse(strategy="preorder"):
                    if node.props["orig_date"] is not None and node.props["imputed_date"] and node.props["imputation_type"] < 6:
                        diffs[i].append(np.abs(np.log(node.props["date"]/node.props["orig_date"])))
                        lindiffs[i].append(np.abs(node.props["date"]-node.props["orig_date"]))

                for leaf in clade_tre.leaves():
                    pendants[i].append(leaf.dist)

                gammas.append(tree_metrics.compute_gamma(clade_tre))

            for i in range(len(diffs)):
                print(len(diffs[i]))
                all_results[(clade, date_pct, params[i][3], rpt)] = (np.mean(diffs[i]), # 0
                                                                     clade_tre.props["num_leaves"], # 1
                                                                     clade_tre.props["child_tree_size"], # 2
                                                                     clade_tre_orig.props["num_dates"], # 3
                                                                     len(diffs[i]), # 4
                                                                     np.mean(pendants[i]), # 5
                                                                     gammas[i], # 6
                                                                     np.mean(lindiffs[i]), # 7
                                                                     np.median(pendants[i]), # 8
                                                                     root_date) # 9
                                                                     # mean_pen, # 10
                                                                     # med_pen) # 11


fout = open("figures/figure2/g_interp_class_medians.txt","w")

for date_pct in date_pcts:
    for clade in clades:
        fout.write("%s\t%f\t%d\t%d\t%d\t%f" % (clade,
                                           date_pct,
                                           all_results[(clade, date_pct, param[3], 0)][1],
                                           all_results[(clade, date_pct, param[3], 0)][2],
                                           all_results[(clade, date_pct, param[3], 0)][3],
                                           all_results[(clade, date_pct, param[3], 0)][9]))

        for param in params:
            means = []
            linmeans = []
            n = []
            pens = []
            penmeds = []
            gams = []
            for rpt in range(repeats):
                means.append(all_results[(clade, date_pct, param[3], rpt)][0])
                linmeans.append(all_results[(clade, date_pct, param[3], rpt)][7])
                n.append(all_results[(clade, date_pct, param[3], rpt)][4])
                pens.append(all_results[(clade, date_pct, param[3], rpt)][5])
                penmeds.append(all_results[(clade, date_pct, param[3], rpt)][8])
                gams.append(all_results[(clade, date_pct, param[3], rpt)][6])

            fout.write("\t%f\t%f\t%f\t%f\t%f" % (np.mean(means), np.mean(linmeans), np.mean(pens), np.mean(penmeds), np.mean(gams)))

        fout.write("\t%f" % (np.mean(n)))

        fout.write("\n")

fout.close()
