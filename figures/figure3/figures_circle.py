import os
os.chdir('/Users/Jonathan/Library/CloudStorage/Dropbox/Imperial/Tree_of_Life/Open_Tree/python/dated-complete-tree')

import sys
import gc

import tree_loading
import tree_labelling
import tree_fixing
import tree_dating

from taxonomy_utils import tx_levels
import tree_plotting
import tree_metrics

import random
import numpy as np
from copy import deepcopy

import logging
logger = logging.getLogger(__name__)
logging.basicConfig(filename="main.log", filemode="w", force=True, level=logging.ERROR)

sys.setrecursionlimit(10000)

rng = np.random.default_rng(seed=1)

# Load metadata for tree from Open Tree, Chronosynth and OneZoom
dates, phylogeny_nodes, taxa = tree_loading.load_metadata()

# Create ETE4 tree structure for entire Open Tree of Life, with my annotations
whole_tre_unmodified = tree_loading.build_and_annotate_tree(phylogeny_nodes,
                                                            taxa)

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

for node in whole_tre.search_nodes(name="Loricifera_ott199402"):
    lori = node
    break

tmp = lori.children[0].copy()

lori.add_child(tmp)

# Optional - mainly here because it's good for OneZoom: remove one-child nodes. Gives a fully bifurcating tree.
whole_tre = tree_fixing.delete_one_child_nodes(whole_tre)

# Date imputation
tree_dating.assign_dates(whole_tre, dates)

# Date cleaning to ensure time consistency down the tree
tree_dating.label_older_descendants(whole_tre)
tree_dating.dq_date_removal(whole_tre)

# Date imputation
tree_dating.date_labelling(whole_tre)
tree_dating.impute_missing_dates(whole_tre, l=0.25)

tree_dating.compute_branch_lengths(whole_tre)

lori.children[1].detach()

for node in whole_tre.traverse():
    if node is not whole_tre and node.dist < 0:
        print(node.name, node.dist)


# import tree_metrics
tree_metrics.compute_pd(whole_tre)


# ######### compute percentages of nodes with dates, cut to phylum or below, plot


def label_pct_dates(parent):
    if parent.is_leaf:
        results = [0,0,1]
    elif parent.props["imputed_date"]:
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


# whole_tre = whole_tre_orig.copy()

def get_nodes_for_trimming(parent, to_remove, keep_branches=True):
    if parent.props["num_leaves"] < 10:
        to_remove.append(parent)
        return
    elif parent.name == "Bacteria_ott844192" or parent.name == "Archaea_ott996421":
        to_remove.append(parent)
        return
    elif tx_levels[parent.props["tx_level"]] == tx_levels["phylum"] or tx_levels[parent.props["tx_level"]] == tx_levels["infrakingdom"]:
        to_remove.append(parent)
        return
    elif tx_levels[parent.props["desc_rank"]] < tx_levels["phylum"]:
        if keep_branches:
            if parent.props["tx_level"] != "mrca":
                to_remove.append(parent)
                return
        else:
            to_remove.append(parent)
            return

    for child in parent.children:
        get_nodes_for_trimming(child, to_remove, keep_branches)

to_remove = []
get_nodes_for_trimming(whole_tre, to_remove, keep_branches=False)

for node in to_remove:
    print(node.name, node.props["tx_level"], node.props["desc_rank"])

# whole_tre_orig = whole_tre.copy()

for node in to_remove:
    tree_fixing.remove_tree_below(node)

    if tx_levels[node.props["desc_rank"]] < tx_levels["phylum"] or tx_levels[node.props["desc_rank"]] > tx_levels["infrakingdom"] or node.props["num_leaves"] < 10:
        node.detach()
        del node

for node in whole_tre.traverse(strategy="preorder"):
    print(node.name, node.props["tx_level"], node.props["desc_rank"], node.props["ancestral_rank"])


for node in whole_tre.search_nodes(name="Eukaryota_ott304358"):
    euk_tre = node
    break

for node in euk_tre.traverse(strategy="preorder"):
    print(node.name, node.props["tx_level"], node.props["desc_rank"], node.props["ancestral_rank"])

for node in euk_tre.search_nodes(name="mrcaott4101ott21309"):
    node.name = "Rotifera_ott471706"
    node.props["tx_level"] = "phylum"
    break

# euk_tre = whole_tre.copy()

for node in euk_tre.traverse():
    if node.props["tx_level"] == "mrca" and len(node.children) == 0:
        node.detach()

tree_fixing.delete_one_child_nodes(euk_tre, maintain_branch_lengths=True)

euk_tre_orig = euk_tre.copy()

for leaf in euk_tre.leaves():
    print(leaf.name,
          leaf.props["tx_level"],
          leaf.props["ph_tx"],
          leaf.props["date"],
          leaf.props["imputed_date"],
          leaf.props["pd"],
          leaf.props["num_leaves"])

# tree_plotting.plot_dates_figure_outline(euk_tre_orig, "figures/figure3/euk_to_phylum_eqsls.svg", log_scale_dates=False, log_scale_branches=True, simple_label=True)
# tree_plotting.plot_dates_figure_outline(euk_tre_orig, "euk_to_phylum_eqsls.tif", log_scale_dates=False, log_scale_branches=True, simple_label=True)
