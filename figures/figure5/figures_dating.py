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


from dated_complete_tree import tree_loading
from dated_complete_tree import tree_fixing_utils
from dated_complete_tree import tree_dating
from dated_complete_tree import tree_plotting

import random


# Load metadata for tree from Open Tree, Chronosynth and OneZoom
phylogeny_nodes, taxa = tree_loading.load_metadata()

test_tre = tree_loading.build_and_annotate_tree(phylogeny_nodes,
                                                taxa,
                                                tree_filename="figures/figure5/small_dating_tree.tre")


test_tre_unmod = test_tre.copy()

test_tre = test_tre_unmod.copy()

alphabet = ['A',
            'B',
            'C',
            'D',
            'E',
            'F',
            'G',
            'H',
            'I',
            'J',
            'K',
            'L',
            'M',
            'N',
            'O',
            'P',
            'Q',
            'R',
            'S',
            'T',
            'U']

test_tre = test_tre.children[1].detach()

test_tre.props["date"] = 5.19

test_tre.children[0].children[1].props["date"] = None
test_tre.children[1].children[1].props["date"] = None
test_tre.children[1].children[0].children[0].children[0].children[2].detach()

A = test_tre.children[0].children[0].children[0].children[0].children[0].children[0].children[0]
A.detach()
test_tre.children[0].children[0].children[0].children[0].children[0].add_child(A)
test_tre.children[0].children[0].children[0].children[0].children[0].children[0].detach()

C = test_tre.children[0].children[0].children[0].children[0].children[0].children[0].children[1].children[0]
C.detach()
C = test_tre.children[0].children[0].children[0].children[0].children[0].children[0].add_child(C)
test_tre.children[0].children[0].children[0].children[0].children[0].children[0].children[1].detach()

D = test_tre.children[0].children[0].children[0].children[0].children[1].children[0]
D.detach()
test_tre.children[0].children[0].children[0].children[0].add_child(D)
test_tre.children[0].children[0].children[0].children[0].children[1].detach()

C = test_tre.children[0].children[0].children[0].children[0].children[0].children[1]
C.detach()

new_node = tree_fixing_utils.create_node("mrca")

test_tre.children[0].children[0].children[0].children[0].children[0].add_child(new_node)
new_node.add_child(C)
new_node.add_child(C.copy())

test_tre.children[0].children[0].children[0].children[0].children[0].children[0].props["date"] = 0.52
test_tre.children[0].children[0].children[0].children[0].children[0].props["date"] = 1.14
test_tre.children[0].children[0].children[0].props["date"] = 2.44
test_tre.children[0].children[0].props["date"] = 3.39
test_tre.children[0].children[1].props["date"] = 3.04

print(test_tre.name, test_tre.props["date"])

i = 0
for node in test_tre.traverse(strategy='preorder'):
    if node.is_leaf:
        print(node.name, i)
        node.name = alphabet[i]
        i += 1
        node.add_prop("date", 0)

    if not "date" in node.props:
        node.add_prop("date", None)

    node.add_prop("imputed_date", False)
    node.add_prop("imputation_type", 0)

random.seed(1)

tree_dating.remove_inconsistent_dates(test_tre, test_tre.props["date"]+1)
tree_dating.date_labelling(test_tre)

tree_plotting.plot_dates_algo(test_tre, "figures/figure5/small_dating_tree_4c.svg", pinkblue=False, show_only_undated_paths=False, vs=0)
tree_plotting.plot_dates_algo(test_tre, "figures/figure5/small_dating_tree_4c.tif", pinkblue=False, show_only_undated_paths=False, vs=0)

tree_dating.impute_missing_dates(test_tre)

tree_dating.compute_branch_lengths(test_tre)

test_tre.dist = 0.7

scale = 100

for node in test_tre.traverse():
    if node.is_leaf:
        node.dist -= 10/scale
    elif node.props["imputed_date"]:
        node.dist -= 17/scale
    else:
        node.dist -= 18/scale

    if node.up and len(node.up.children) > 1:
        node.dist -= 1/scale
    node.props["adj"] = True

tree_plotting.plot_dates_algo(test_tre, "figures/figure5/small_dating_tree_5c.svg", pinkblue=False, show_paths=False, vs=7)
tree_plotting.plot_dates_algo(test_tre, "figures/figure5/small_dating_tree_5c.tif", pinkblue=False, show_paths=False, vs=7)

tree_dating.compute_branch_lengths(test_tre)

test_tre.dist = 0.7

for node in test_tre.traverse():
    if node.props["imputed_date"]:
        node.props["date"] = None

scale = 100

for node in test_tre.traverse():
    if node.is_leaf:
        node.dist -= 10/scale
    elif node.props["date"]:
        node.dist -= 18/scale
    else:
        node.dist -= 17/scale

    if node.up and len(node.up.children) > 1:
        node.dist -= 1/scale
    node.props["adj"] = True


tree_plotting.plot_dates_algo(test_tre, "figures/figure5/small_dating_tree_6c.svg", pinkblue=False, show_only_undated_paths=False, vs=0)
tree_plotting.plot_dates_algo(test_tre, "figures/figure5/small_dating_tree_6c.tif", pinkblue=False, show_only_undated_paths=False, vs=0)
