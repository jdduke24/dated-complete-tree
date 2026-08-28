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
from dated_complete_tree import tree_metrics

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
tree_dating.assign_dates(whole_tre, dates)


#####################################################################################################################

# Write out newick files for clades of the tree, each with an "ages" file listing the dates available in the phylesystem

# Crown nodes of the trees:

# Columbiformes_ott363030: 685 nodes
# Muridae_ott816256: 2007
# mrcaott1822ott688506: 4783
# Jungermanniales_ott56621: 9069
# Squamata_ott35888: 22441
# Stramenopiles_ott266745: 40175
# Teleostei_ott212201: 74731
# Gastropoda_ott409995: 140247
# Magnoliopsida_ott99252: 746891
# Metazoa_ott691846: 2950313
# cellular_organisms_ott93302: 4589559


trees_to_write = [
    'Columbiformes_ott363030',
    'Muridae_ott816256',
    'mrcaott1822ott688506',
    'Jungermanniales_ott56621',
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

    print(crown, len(list(clade_tre.descendants()))+1, clade_tre.props["date"])

    fout = open("phylocom_tests/%s.ages" % crown,'wt')
    for node in clade_tre.traverse():
        if node.props["date"]:
            fout.write("%s\t%f\n" % (node.name, node.props["date"]))
    fout.close()


    clade_tre.write("phylocom_tests/%s.phylo" % crown,
                    parser=1,
                    format_root_node=True)


#####################

# To be run separately in IPython: read in a Newick file called 'phylo' and a text file called 'ages' - as used
# by Phylocom - and write same output as Phylocom, for fair comparison

# use IPython timeit function
%%timeit -n1 -r2

# read 'phylo' Newick file
clade_tre = ete4.Tree(open("phylocom_tests/phylo","r").read(), parser=1)

# read 'ages' text file
dates = {}
with open("phylocom_tests/ages", 'r') as fin:
    for line in fin:
        bits = line.split('\t')
        dates[bits[0]] = float(bits[1])

# apply dates to nodes in tree
for node in clade_tre.traverse(strategy="preorder"):
    if node.name in dates:
        node.add_prop("date", dates[node.name])
    elif node.is_leaf:
        node.add_prop("date", 0.)
    else:
        node.add_prop("date", None)
    node.add_prop("imputed_date", False)
    node.add_prop("imputation_type", 0)

# either use Phylocom's date consistency fixing
tree_dating.remove_inconsistent_dates(clade_tre)

# or, use our new date consistency fixing
# tree_dating.dq_date_removal(clade_tre)

# interpolate with EQS-L
tree_dating.date_labelling(clade_tre)
tree_dating.impute_missing_dates(clade_tre, l=1)

tree_dating.compute_branch_lengths(clade_tre)

clade_tre.write("phylocom_tests/my_dated_tree.phylo", parser=1)
