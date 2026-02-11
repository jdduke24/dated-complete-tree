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
logging.basicConfig(filename="main.log", filemode="w", force=True, level=logging.ERROR)

sys.setrecursionlimit(10000)

#####################################################################################################################
# Load and prune tree

# Load metadata for tree from Open Tree and Chronosynth
dates, phylogeny_nodes, taxa = tree_loading.load_metadata()

# Create ETE3 tree structure for entire Open Tree of Life, with my annotations
whole_tre_unmodified = tree_loading.build_and_annotate_tree(phylogeny_nodes, taxa)

tree_fixing.strip_birds(whole_tre_unmodified)
tree_fixing.strip_turtles(whole_tre_unmodified)

rng = np.random.default_rng(12345)

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
tree_dating.assign_dates(whole_tre, dates)

date_count = 0
for node in whole_tre.traverse():
    if node.props["date"]:
        date_count += 1

# Date cleaning to ensure time consistency down the tree
tree_dating.label_older_descendants(whole_tre)
tree_dating.dq_date_removal(whole_tre)

date_count = 0
for node in whole_tre.traverse():
    if node.props["date"]:
        date_count += 1

tree_dating.date_labelling(whole_tre)
tree_dating.impute_missing_dates(whole_tre, l=0.25)

tree_dating.compute_branch_lengths(whole_tre)
tree_metrics.compute_pd(whole_tre)

results = {}
base_pd = whole_tre.props["pd"]

date_count = 0
for node in whole_tre.traverse():
    if node.props["imputed_date"] == False and node.props["date"] > 0:
        date_count += 1

date_count = 0
for node in whole_tre.traverse():
    if node.props["imputed_date"] == True:
        date_count += 1


import datetime
itr_start = datetime.datetime.now()
num_itrs = 2000

for itr in range(num_itrs):
    for node in whole_tre.traverse():
        # reset all imputed dates
        if node.props["imputed_date"] == True:
            node.props["date"] = None

    for node in whole_tre.traverse():
        # change one date
        if "mrca" not in node.name and node.props["imputed_date"] == False and node.props["date"] > 0 and "sens_done" not in node.props:
            node.props["date"] *= 1.001
            node.add_prop("sens_done", True)
            break

    print(node.name, "iteration", itr+1, "/ projected end time:", itr_start + (num_itrs)*(datetime.datetime.now() - itr_start)/itr if itr > 0 else "")

    # Date imputation
    tree_dating.date_labelling(whole_tre)
    dating_rng = np.random.default_rng(seed=100)
    tree_dating.impute_missing_dates(whole_tre, l=0.25, rng=dating_rng)

    tree_dating.compute_branch_lengths(whole_tre)
    tree_metrics.compute_pd(whole_tre)

    results[node] = whole_tre.props["pd"]

    # change date back
    node.props["date"] /= 1.001

print(datetime.datetime.now())


def label_pct_dates(parent):
    if parent.is_leaf:
        results = [0,0,1,0]
    elif parent.props["imputed_date"]:
        results = [1,0,0,0]
    else:
        results = [1,1,0,0]

    if parent.props["ph_tx"] == "PH":
        results[3] += 1

    for child in parent.children:
        new_results = label_pct_dates(child)
        results[0] += new_results[0]
        results[1] += new_results[1]
        results[2] += new_results[2]
        results[3] += new_results[3]

    parent.add_prop("child_tree_size", results[0])
    parent.add_prop("num_dates", results[1])
    parent.add_prop("num_leaves", results[2])
    parent.add_prop("num_ph", results[3])

    return results

label_pct_dates(whole_tre)

fout = open("pd_date_sensitivity.txt", 'w')
for node in results:
    fout.write("%s\t%f\t%s\t%f\t%f\t%f\n" % (node.name,
                                         node.props["date"],
                                         node.props["tx_level"],
                                         results[node],
                                         node.props["child_tree_size"],
                                         node.props["num_dates"]))
fout.close()
