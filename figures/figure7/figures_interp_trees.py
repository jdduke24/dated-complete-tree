import os
os.chdir('/Users/Jonathan/Library/CloudStorage/Dropbox/Imperial/Tree_of_Life/Open_Tree/python/dated-complete-tree/figures/figure7')

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
dates, phylogeny_nodes, taxa, descr_years = tree_loading.load_metadata()

test_tre = tree_loading.build_and_annotate_tree(dates,
                                                            phylogeny_nodes,
                                                            taxa,
                                                            descr_years,
                                                            tree_filename="test_trees/small_dating_tree_2.tre")


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



test_tre.date = 5.3827

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


test_tre.children[1].children[1].children[1].detach()
test_tre.children[1].children[1].children[0].detach()

test_tre.children[1].children[0].children[0].children[1].detach()
test_tre.children[1].children[0].children[0].children[0].detach()

test_tre.children[0].children[0].children[1].detach()
test_tre.children[0].children[0].children[0].detach()

i = 0
for node in test_tre.traverse(strategy='preorder'):
    if node.is_leaf():
        print(node.name, i)
        node.name = alphabet[i]
        i += 1


start_tre = test_tre.copy()

tree_plotting.plot_ultrametric(start_tre, "small_bm_test.svg")

for node in start_tre.traverse():
    if node.is_leaf():
        node.date = 0
    elif node is start_tre:
        node.date = 12.
    else:
        node.date = None

rng = np.random.default_rng()

tree_dating.date_labelling(start_tre)
tree_dating.impute_missing_dates(start_tre, useBirth=True, rng=rng)
tree_dating.compute_branch_lengths(start_tre)

tree_plotting.plot_ultrametric(start_tre, "small_bm_test_imputed.svg")


test_tre = start_tre.copy()

# long BLADJ
test_tre.children[0].dist = 4
test_tre.children[0].children[0].dist = 8

test_tre.children[0].children[1].dist = 4
test_tre.children[0].children[1].children[0].dist = 4
test_tre.children[0].children[1].children[1].dist = 4


test_tre.children[1].dist = 4
test_tre.children[1].children[1].dist = 8

test_tre.children[1].children[0].dist = 4
test_tre.children[1].children[0].children[0].dist = 4
test_tre.children[1].children[0].children[1].dist = 4

tree_plotting.plot_ultrametric(test_tre, "small_interp_tree_1.svg", color_long=True)


test_tre = start_tre.copy()

# short BLADJ
test_tre.children[0].dist = 6
test_tre.children[0].children[0].dist = 6

test_tre.children[0].children[1].dist = 3
test_tre.children[0].children[1].children[0].dist = 3
test_tre.children[0].children[1].children[1].dist = 3


test_tre.children[1].dist = 6
test_tre.children[1].children[1].dist = 6

test_tre.children[1].children[0].dist = 3
test_tre.children[1].children[0].children[0].dist = 3
test_tre.children[1].children[0].children[1].dist = 3

tree_plotting.plot_ultrametric(test_tre, "small_interp_tree_2.svg", color_short=True)


test_tre = start_tre.copy()

# log(N)
dt1 = 12 * np.log(3) / np.log(6)
dt2 = dt1 * np.log(2) / np.log(3)

print(dt1, dt2)

test_tre.children[0].dist = 12 - dt1
test_tre.children[0].children[0].dist = dt1

test_tre.children[0].children[1].dist = dt1 - dt2
test_tre.children[0].children[1].children[0].dist = dt2
test_tre.children[0].children[1].children[1].dist = dt2


test_tre.children[1].dist = 12 - dt1
test_tre.children[1].children[1].dist = dt1

test_tre.children[1].children[0].dist = dt1 - dt2
test_tre.children[1].children[0].children[0].dist = dt2
test_tre.children[1].children[0].children[1].dist = dt2

tree_plotting.plot_ultrametric(test_tre, "small_interp_tree_3.svg")


test_tre = start_tre.copy()

# BM
test_tre.children[0].dist = 12 - 7.36
test_tre.children[0].children[0].dist = 7.36

test_tre.children[0].children[1].dist = 7.36 - 1.22
test_tre.children[0].children[1].children[0].dist = 1.22
test_tre.children[0].children[1].children[1].dist = 1.22


test_tre.children[1].dist = 12 - 4.64
test_tre.children[1].children[1].dist = 4.64

test_tre.children[1].children[0].dist = 4.64 - 2.71
test_tre.children[1].children[0].children[0].dist = 2.71
test_tre.children[1].children[0].children[1].dist = 2.71

tree_plotting.plot_ultrametric(test_tre, "small_interp_tree_4.svg")


path_tre = ete3.TreeNode(name="root")

a_node = ete3.TreeNode(name="A")
path_tre.add_child(a_node)

int1 = ete3.TreeNode(name="int1")
int2 = ete3.TreeNode(name="int2")

sp1 = ete3.TreeNode(name="sp1")
sp2 = ete3.TreeNode(name="sp2")
sp3 = ete3.TreeNode(name="sp3")
sp4 = ete3.TreeNode(name="sp4")
sp5 = ete3.TreeNode(name="sp5")
sp5.dist = 4.2

path_tre.add_child(sp5)

a_node.add_child(int1)
a_node.add_child(int2)

int1.add_child(sp1)
int1.add_child(sp2)
int2.add_child(sp3)
int2.add_child(sp4)

tree_plotting.plot_ultrametric_interp(path_tre, "small_interp_tree_5.svg")


test_tre = start_tre.copy()

# EQS-LS
test_tre.children[0].dist = 5.5
test_tre.children[0].children[0].dist = 6.5

test_tre.children[0].children[1].dist = 3.25
test_tre.children[0].children[1].children[0].dist = 3.25
test_tre.children[0].children[1].children[1].dist = 3.25


test_tre.children[1].dist = 5.5
test_tre.children[1].children[1].dist = 6.5

test_tre.children[1].children[0].dist = 3.25
test_tre.children[1].children[0].children[0].dist = 3.25
test_tre.children[1].children[0].children[1].dist = 3.25

tree_plotting.plot_ultrametric(test_tre, "small_interp_tree_6.svg")
