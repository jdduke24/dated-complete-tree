import os
os.chdir('/Users/Jonathan/Library/CloudStorage/Dropbox/Imperial/Tree_of_Life/Open_Tree/python/dated-complete-tree')

import sys
import numpy as np

import tree_loading
import tree_labelling
import tree_fixing
import tree_dating
import tree_metrics

import logging
logger = logging.getLogger(__name__)
logging.basicConfig(filename="main.log", filemode="w", force=True, level=logging.WARNING)

from copy import copy

#####################################################################################################################

sys.setrecursionlimit(10000)

# Load metadata for tree from Open Tree, Chronosynth and OneZoom
dates, phylogeny_nodes, taxa, descr_years = tree_loading.load_metadata()

# Create ETE3 tree structure for entire Open Tree of Life, with my annotations
whole_tre_unmodified = tree_loading.build_and_annotate_tree(dates, phylogeny_nodes, taxa, descr_years, tree_filename="frogs_dna_ott.tre", has_branch_lengths=True)

tree_dating.label_older_descendants(whole_tre_unmodified)
tree_dating.dq_date_removal(whole_tre_unmodified)

# Copy tree - we will change the copy, and keep the original unchanged so we can restore it next iteration without
# reloading everything
whole_tre = whole_tre_unmodified.copy()

# Remove anything below species level (this includes promoting some subspecies to species if they are the only
# example of their species)

rng = np.random.default_rng(seed=1)

tree_fixing.remove_subspecies(whole_tre, rng)
tree_fixing.impute_species_into_empty_taxa(whole_tre)

tree_labelling.add_anc_ranks(whole_tre)
tree_labelling.add_desc_ranks(whole_tre)

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

# Optional - mainly here because it's good for OneZoom: remove one-child nodes. Gives a fully bifurcating tree.
whole_tre = tree_fixing.delete_one_child_nodes(whole_tre, maintain_branch_lengths=True)

#####################################################################################################################

# Date imputation
tree_dating.remove_inconsistent_dates(whole_tre, whole_tre.date+1)
tree_dating.date_labelling(whole_tre)

counts = [0,0]
# tree_dating.impute_missing_dates(whole_tre, useLnN=True, counts=counts)

dating_rng = np.random.default_rng(seed=100)
tree_dating.impute_missing_dates(whole_tre, useLnN=True, rng=dating_rng, counts=counts)

tree_dating.compute_branch_lengths(whole_tre)

dists = []
for leaf in whole_tre.get_leaves():
    dists.append(leaf.dist)

print("Pendant mean:", np.mean(dists))
print("Pendant median:", np.median(dists))



# date_pcts = [0.25, 0.5] # pct to throw away
date_pcts = [0.1, 0.25, 0.33, 0.5, 0.6, 0.7, 0.8, 0.9, 0.93, 0.96, 0.98, 0.99, 0.996, 1] # pct to throw away

params = [[1, 0, False, "Long", False],
          [0, 0, False, "Short", False],
          [1, 0, True, "LnN", False],
          [1, 0, False, "Birth", True],
          [0.25, 0, False, "LS", False]]

# params = [[0.2, 0, False, "LS1", False],
#           [0.3, -1, False, "LS4", False],
#           [0.2, 0.25, False, "LS2", False],
#           [0.3, 0.25, False, "LS3", False],
#           [0.25, 0.3, False, "LS5", False]]

# params = [[1, 0, False, "Birth", True]]

repeats = 10

all_results = {}

rng_main = np.random.default_rng(seed=100)
dating_rng = np.random.default_rng(seed=100)


for date_pct in date_pcts:
    for rpt in range(repeats):

        diffs = [[] for p in params]
        lindiffs = [[] for p in params]
        pendants = [[] for p in params]
        gammas = []

        seed = rng_main.integers(1000)

        for i, param in enumerate(params):
            clade_tre = whole_tre.copy()

            rng = np.random.default_rng(seed=seed)

            count = 0
            for node in clade_tre.traverse():
                node.orig_date = None
                if node is clade_tre:
                    continue
                if node.date:
                    # if node.date has a value and the value is > 0:
                    if rng.random() < date_pct:
                        node.orig_date = copy(node.date)
                        node.date = None
                        count += 1

            print(clade_tre.date, param[3], count)

            tree_dating.date_labelling(clade_tre)

            tree_dating.impute_missing_dates(clade_tre, l=param[0], m=param[1], useLnN=param[2], useBirth=param[4], rng=dating_rng)
            tree_dating.compute_branch_lengths(clade_tre)

            for node in clade_tre.traverse():
                if node.orig_date is not None and node.imputed_date and node.imputation_type < 6:
                    diffs[i].append(np.abs(np.log(node.date/node.orig_date)))
                    lindiffs[i].append(np.abs(node.date-node.orig_date))

            for leaf in clade_tre.get_leaves():
                pendants[i].append(leaf.dist)

            gammas.append(tree_metrics.compute_gamma(clade_tre))

        for i in range(len(diffs)):
            print(len(diffs[i]))
            all_results[(date_pct, params[i][3], rpt)] = (np.mean(diffs[i]),
                                                            clade_tre.num_leaves,
                                                            clade_tre.child_tree_size,
                                                            whole_tre.num_dates,
                                                            len(diffs[i]),
                                                            np.mean(pendants[i]),
                                                            gammas[i],
                                                            np.mean(lindiffs[i]),
                                                            np.median(pendants[i]))


fout = open("frogs_interp_with_medians.txt","w")

clade = "Anura_ott991547"
for date_pct in date_pcts:
    fout.write("%s\t%f\t%d\t%d\t%d" % (clade,
                                       date_pct,
                                       all_results[(date_pct, param[3], 0)][1],
                                       all_results[(date_pct, param[3], 0)][2],
                                       all_results[(date_pct, param[3], 0)][3]))
    for param in params:
        means = []
        linmeans = []
        n = []
        pens = []
        penmeds = []
        gams = []
        for rpt in range(repeats):
            means.append(all_results[(date_pct, param[3], rpt)][0])
            linmeans.append(all_results[(date_pct, param[3], rpt)][7])
            n.append(all_results[(date_pct, param[3], rpt)][4])
            pens.append(all_results[(date_pct, param[3], rpt)][5])
            penmeds.append(all_results[(date_pct, param[3], rpt)][8])
            gams.append(all_results[(date_pct, param[3], rpt)][6])

        fout.write("\t%f\t%f\t%f\t%f\t%f" % (np.mean(means), np.mean(linmeans), np.mean(pens), np.mean(penmeds), np.mean(gams)))

    fout.write("\t%f" % (np.mean(n)))

    fout.write("\n")

fout.close()

# for date_pct in date_pcts:
#     fout.write("%f\t%f\t%f\t%f" % (date_pct,
#                                        all_results[(date_pct, param[3], 0)][1],
#                                        all_results[(date_pct, param[3], 0)][2],
#                                        all_results[(date_pct, param[3], 0)][3]))
#     for param in params:
#         means = []
#         n = []
#         for rpt in range(repeats):
#             means.append(all_results[(date_pct, param[3], rpt)][0])
#             n.append(all_results[(date_pct, param[3], rpt)][3])


#         fout.write("\t%f\t%f" % (np.mean(means), np.mean(n)))
#     fout.write('\n')

# fout.close()
