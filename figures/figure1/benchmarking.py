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
import gc

import tree_loading
import tree_labelling
import tree_fixing
import tree_dating

import random
import numpy as np

import logging
logger = logging.getLogger(__name__)
logging.basicConfig(filename="main.log", filemode="w", force=True, level=logging.ERROR)

sys.setrecursionlimit(10000)

#####################################################################################################################
# Load and prune tree

# Load metadata for tree from Open Tree and Chronosynth
dates, phylogeny_nodes, taxa = tree_loading.load_metadata()

# Create ete4 tree structure for entire Open Tree of Life, with my annotations
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

# trees:
# Columbiformes_ott363030: 685 nodes
# Muridae_ott816256: 2007
# mrcaott1822ott688506: 4783
# Jungermanniales_ott56621: 9069
# Passeriformes_ott1041547: 13179
# Aves_ott81461: 21843
# Teleostei_ott212201: 74731
# Gastropoda_ott409995: 140247
# Magnoliopsida_ott99252: 746891
# Metazoa_ott691846: 2950313
# cellular_organisms_ott93302: 4589559

os.chdir('/Users/Jonathan/Library/CloudStorage/Dropbox/Imperial/Tree_of_Life/Open_Tree/python/dated-complete-tree/figures/figure1')
import ete4

tree_dating.assign_dates(whole_tre, dates)

trees_to_write = [
    'Columbiformes_ott363030',
    'Muridae_ott816256',
    'mrcaott1822ott688506',
    'Jungermanniales_ott56621',
    'Passeriformes_ott1041547',
    'Aves_ott81461',
    'Squamata_ott35888',
    'Stramenopiles_ott266745',
    'Teleostei_ott212201',
    'Gastropoda_ott409995',
    'Magnoliopsida_ott99252',
    'Metazoa_ott691846',
    'cellular_organisms_ott93302']

for crown in trees_to_write:
    for node in whole_tre.search_nodes(name=crown):
        clade_tre = node
        break

    # fout = open("phylocom_tests/all_ages_and_sizes.txt",'wt')

    # for node in clade_tre.traverse():
    #     if node.props["date"]:
    #         fout.write("%s\t%f\t%d\n" % (node.name, node.props["date"], len(list(node.descendants()))+1))
    # fout.close()


    print(crown, len(list(clade_tre.descendants()))+1, clade_tre.props["date"])

    fout = open("phylocom_tests/%s.ages" % crown,'wt')
    for node in clade_tre.traverse():
        if node.props["date"]:
            fout.write("%s\t%f\n" % (node.name, node.props["date"]))
    fout.close()


    clade_tre.write("phylocom_tests/%s.phylo" % crown,
                    parser=1,
                    format_root_node=True)


##
# Read and write same tree files as Phylocom, for fair comparison

%%timeit -n1 -r2
clade_tre = ete4.Tree(open("phylocom_tests/phylo","r").read(), parser=1)
dates = {}
with open("phylocom_tests/ages", 'r') as fin:
    for line in fin:
        bits = line.split('\t')
        dates[bits[0]] = float(bits[1])

for node in clade_tre.traverse(strategy="preorder"):
    if node.name in dates:
        node.add_prop("date", dates[node.name])
    elif node.is_leaf:
        node.add_prop("date", 0.)
    else:
        node.add_prop("date", None)
    node.add_prop("imputed_date", False)
    node.add_prop("imputation_type", 0)

tree_dating.remove_inconsistent_dates(clade_tre)

# tree_dating.label_older_descendants(clade_tre)
# tree_dating.dq_date_removal(clade_tre)

tree_dating.date_labelling(clade_tre)
tree_dating.impute_missing_dates(clade_tre, l=1)

tree_dating.compute_branch_lengths(clade_tre)
clade_tre.write("phylocom_tests/my_dated_tree.phylo", parser=1)
