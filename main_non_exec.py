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
logging.basicConfig(filename="main.log", filemode="w", force=True, level=logging.ERROR)

#####################################################################################################################

sys.setrecursionlimit(10000)

# Load metadata for tree from Open Tree, Chronosynth and OneZoom
dates, phylogeny_nodes, taxa = tree_loading.load_metadata()

# Create ETE3 tree structure for entire Open Tree of Life, with my annotations
whole_tre_unmodified = tree_loading.build_and_annotate_tree(dates, phylogeny_nodes, taxa)

tree_fixing.strip_birds(whole_tre_unmodified)
tree_fixing.strip_turtles(whole_tre_unmodified)

tree_dating.label_older_descendants(whole_tre_unmodified)
tree_dating.dq_date_removal(whole_tre_unmodified)

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

# Date imputation
tree_dating.remove_inconsistent_dates(whole_tre, whole_tre.date+1)

tree_dating.date_labelling(whole_tre)
dating_rng = np.random.default_rng(seed=100)
tree_dating.impute_missing_dates(whole_tre, l=0.25, rng=dating_rng)

tree_dating.compute_branch_lengths(whole_tre)
tree_metrics.compute_pd(whole_tre)
