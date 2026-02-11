import os
os.chdir('/Users/Jonathan/Library/CloudStorage/Dropbox/Imperial/Tree_of_Life/Open_Tree/python/dated-complete-tree/figures/figure8')

import sys
import gc

import tree_loading
import tree_labelling
import tree_fixing
import tree_dating
import tree_plotting

import random
import numpy as np
from copy import deepcopy
import ete3

sys.setrecursionlimit(10000)


# Load metadata for tree from Open Tree, Chronosynth and OneZoom
dates, phylogeny_nodes, taxa = tree_loading.load_metadata()

test_tre = tree_loading.build_and_annotate_tree(dates,
                                                phylogeny_nodes,
                                                taxa,
                                                tree_filename="../../test_trees/small_dating_tree_2.tre")


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

test_tre = test_tre.children[1]

test_tre.date = 5.19

test_tre.children[0].children[1].date = None
test_tre.children[1].children[1].date = None
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

new_node = tree_fixing.create_node("mrca")

test_tre.children[0].children[0].children[0].children[0].children[0].add_child(new_node)
new_node.add_child(C)
new_node.add_child(C.copy())

test_tre.children[0].children[0].children[0].children[0].children[0].children[0].date = 0.52
test_tre.children[0].children[0].children[0].children[0].children[0].date = 1.14
test_tre.children[0].children[0].children[0].date = 2.44
test_tre.children[0].children[0].date = 3.39
test_tre.children[0].children[1].date = 3.04

print(test_tre.name, test_tre.date)

i = 0
for node in test_tre.traverse(strategy='preorder'):
    if node.is_leaf():
        print(node.name, i)
        node.name = alphabet[i]
        i += 1

random.seed(1)

tree_dating.remove_inconsistent_dates(test_tre, test_tre.date+1)
tree_dating.date_labelling(test_tre)

tree_plotting.plot_dates_algo(test_tre, "small_dating_tree_4c.svg", pinkblue=False, show_only_undated_paths=False, vs=0)
tree_plotting.plot_dates_algo(test_tre, "small_dating_tree_4c.tif", pinkblue=False, show_only_undated_paths=False, vs=0)

tree_dating.impute_missing_dates(test_tre)

tree_dating.compute_branch_lengths(test_tre)

test_tre.dist = 0.7

scale = 100

for node in test_tre.traverse():
    if node.is_leaf():
        node.dist -= 10/scale
    elif node.imputed_date:
        node.dist -= 17/scale
    else:
        node.dist -= 18/scale

    if len(node.up.children) > 1:
        node.dist -= 1/scale
    node.adj = True

tree_plotting.plot_dates_algo(test_tre, "small_dating_tree_5c.svg", pinkblue=False, show_paths=False, vs=7)
tree_plotting.plot_dates_algo(test_tre, "small_dating_tree_5c.tif", pinkblue=False, show_paths=False, vs=7)

tree_dating.compute_branch_lengths(test_tre)

test_tre.dist = 0.7

for node in test_tre.traverse():
    if node.imputed_date:
        node.date = None

scale = 100

for node in test_tre.traverse():
    if node.is_leaf():
        node.dist -= 10/scale
    elif node.date:
        node.dist -= 18/scale
    else:
        node.dist -= 17/scale

    if len(node.up.children) > 1:
        node.dist -= 1/scale
    node.adj = True


tree_plotting.plot_dates_algo(test_tre, "small_dating_tree_6c.svg", pinkblue=False, show_only_undated_paths=False, vs=0)
tree_plotting.plot_dates_algo(test_tre, "small_dating_tree_6c.tif", pinkblue=False, show_only_undated_paths=False, vs=0)
