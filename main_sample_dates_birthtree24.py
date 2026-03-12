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
dates, phylogeny_nodes, taxa = tree_loading.load_metadata(force_no_dates_refresh=True)

# Create ETE3 tree structure for entire Open Tree of Life, with my annotations
whole_tre = tree_loading.build_and_annotate_tree(phylogeny_nodes, taxa, tree_filename="trees_bm/birth_model_topo_sample_24.tre", has_branch_lengths=True)

tree_metrics.compute_pd(whole_tre)

base_pd = whole_tre.props["pd"]

# whole_tre.write("output/test_tre.tre", parser=1, format_root_node=True)



pd_clades = [cld.strip() for cld in list(open("pd_clades.txt"))]

pd_dict = {}
dates_dict = {}
spp_dict = {}
for clade in pd_clades:
    pd_dict[clade] = []
    dates_dict[clade] = []
    spp_dict[clade] = []

tree_metrics.save_pd_for_clades(whole_tre, pd_clades, pd_dict, dates_dict, spp_dict)


# Assign dates from varying sources
date_interpolation_rng = np.random.default_rng(seed=10)
date_source_rng = np.random.default_rng(seed=10)

results = []

import datetime
itr_start = datetime.datetime.now()
num_itrs = 101

for itr in range(num_itrs):
    print("iteration", itr+1, "/ projected end time:", itr_start + (num_itrs)*(datetime.datetime.now() - itr_start)/itr if itr > 0 else "")
    for node in whole_tre.traverse(strategy="preorder"):
        # reset all dates
        node.props["date"] = None
        node.props["imputed_date"] = False
        node.props["imputation_type"] = 0

    date_sources = tree_dating.assign_dates(whole_tre, dates, sample_dates=True, rng=date_source_rng)

    # Date cleaning to ensure time consistency down the tree
    tree_dating.label_older_descendants(whole_tre)
    tree_dating.dq_date_removal(whole_tre)

    # Date imputation
    tree_dating.date_labelling(whole_tre)
    tree_dating.impute_missing_dates(whole_tre, use_birth_model=True, rng=date_interpolation_rng)

    tree_dating.compute_branch_lengths(whole_tre)

    tree_dating.write_tree_with_branch_lengths(whole_tre, filename="output/birth_model_date_sample_%d.tre" % (itr+1))

    tree_metrics.compute_pd(whole_tre)

    results.append(whole_tre.props["pd"])
    tree_metrics.save_pd_for_clades(whole_tre, pd_clades, pd_dict, dates_dict, spp_dict)

tree_metrics.write_pd_dists("output/birth_model_date", pd_dict, dates_dict, spp_dict)
