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

# Date cleaning to ensure time consistency down the tree
tree_dating.dq_date_removal(whole_tre)

# Date imputation
dating_rng = np.random.default_rng(seed=100)
tree_dating.impute_missing_dates(whole_tre, l=0.25, rng=dating_rng)

tree_dating.compute_branch_lengths(whole_tre)
tree_metrics.compute_pd(whole_tre)
