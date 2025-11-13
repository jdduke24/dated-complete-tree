import os
os.chdir('/Users/Jonathan/Library/CloudStorage/Dropbox/Imperial/Tree_of_Life/Open_Tree/python/dated-complete-tree/figures/figure6')

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


root = tree_fixing.create_node("Passeriformes")
root.tx_level = "order"
root.ancestral_rank = "order"
root.desc_rank = "order"
root.ph_tx = "PH"
root.date = 100.

family_ph1 = tree_fixing.create_node("Alaudidae")
family_ph1.tx_level = "family"
family_ph1.ancestral_rank = "family"
family_ph1.desc_rank = "family"
family_ph1.ph_tx = "PH"

family_tx = tree_fixing.create_node("Colluricinclidae")
family_tx.tx_level = "family"
family_tx.ancestral_rank = "family"
family_tx.desc_rank = "family"
family_tx.ph_tx = "TX"

genus_tx9 = tree_fixing.create_node("Laphyctes")
genus_tx9.add_feature("genus_name", "Laphyctes")
genus_tx9.tx_level = "genus"
genus_tx9.ancestral_rank = "genus"
genus_tx9.desc_rank = "genus"
genus_tx9.ph_tx = "TX"

mrca1 = tree_fixing.create_node("mrca1")
mrca1.ancestral_rank = "order"
mrca1.desc_rank = "family"
mrca1.ph_tx = "PH"

root.add_child(family_ph1)
root.add_child(mrca1)
root.add_child(family_tx)
root.add_child(genus_tx9)



species_tx4 = tree_fixing.create_node("Laphyctes apolites")
species_tx4.add_feature("genus_name", "Laphyctes")
species_tx4.add_feature("species_name", "Laphyctes apolites")
species_tx4.tx_level = "species"
species_tx4.ancestral_rank = "species"
species_tx4.desc_rank = "species"
species_tx4.ph_tx = "TX"

genus_tx9.add_child(species_tx4)


species_tx5 = tree_fixing.create_node("Rhectes analogus")
species_tx5.add_feature("genus_name", "Rhectes")
species_tx5.add_feature("species_name", "Rhectes analogus")
species_tx5.tx_level = "species"
species_tx5.ancestral_rank = "species"
species_tx5.desc_rank = "species"
species_tx5.ph_tx = "TX"

species_tx6 = tree_fixing.create_node("Rhectes rubiensis")
species_tx6.add_feature("genus_name", "Rhectes")
species_tx6.add_feature("species_name", "Rhectes rubiensis")
species_tx6.tx_level = "species"
species_tx6.ancestral_rank = "species"
species_tx6.desc_rank = "species"
species_tx6.ph_tx = "TX"

species_tx7 = tree_fixing.create_node("Rhectes ferrugineus")
species_tx7.add_feature("genus_name", "Rhectes")
species_tx7.add_feature("species_name", "Rhectes ferrugineus")
species_tx7.tx_level = "species"
species_tx7.ancestral_rank = "species"
species_tx7.desc_rank = "species"
species_tx7.ph_tx = "TX"

family_tx.add_child(species_tx5)
family_tx.add_child(species_tx7)
family_tx.add_child(species_tx6)


mrca2 = tree_fixing.create_node("mrca2")
mrca2.ancestral_rank = "order"
mrca2.desc_rank = "family"
mrca2.ph_tx = "PH"

genus_ph1 = tree_fixing.create_node("Thamnomanes")
genus_ph1.add_feature("genus_name", "Thamnomanes")
genus_ph1.tx_level = "genus"
genus_ph1.ancestral_rank = "genus"
genus_ph1.desc_rank = "genus"
genus_ph1.ph_tx = "PH"

mrca1.add_child(mrca2)
mrca1.add_child(genus_ph1)


species_ph1 = tree_fixing.create_node("Thamnomanes ardesiacus")
species_ph1.add_feature("genus_name", "Thamnomanes")
species_ph1.add_feature("species_name", "Thamnomanes ardesiacus")
species_ph1.tx_level = "species"
species_ph1.ancestral_rank = "species"
species_ph1.desc_rank = "species"
species_ph1.ph_tx = "PH"

species_ph2 = tree_fixing.create_node("Thamnomanes saturninus")
species_ph2.add_feature("genus_name", "Thamnomanes")
species_ph2.add_feature("species_name", "Thamnomanes saturninus")
species_ph2.tx_level = "species"
species_ph2.ancestral_rank = "species"
species_ph2.desc_rank = "species"
species_ph2.ph_tx = "PH"

genus_ph1.add_child(species_ph1)
genus_ph1.add_child(species_ph2)


mrca3 = tree_fixing.create_node("mrca3")
mrca3.ancestral_rank = "order"
mrca3.desc_rank = "family"
mrca3.ph_tx = "PH"

species_ph3 = tree_fixing.create_node("Garrulax virgatus")
species_ph3.add_feature("genus_name", "Garrulax")
species_ph3.add_feature("species_name", "Garrulax virgatus")
species_ph3.tx_level = "species"
species_ph3.ancestral_rank = "species"
species_ph3.desc_rank = "species"
species_ph3.ph_tx = "PH"

mrca2.add_child(mrca3)
mrca2.add_child(species_ph3)


family_ph2 = tree_fixing.create_node("Irenidae")
family_ph2.tx_level = "family"
family_ph2.ancestral_rank = "family"
family_ph2.desc_rank = "family"
family_ph2.ph_tx = "PH"

genus_ph2 = tree_fixing.create_node("Cymbilaimus")
genus_ph1.add_feature("genus_name", "Cymbilaimus")
genus_ph2.tx_level = "genus"
genus_ph2.ancestral_rank = "genus"
genus_ph2.desc_rank = "genus"
genus_ph2.ph_tx = "PH"

mrca3.add_child(family_ph2)
mrca3.add_child(genus_ph2)


genus_ph3 = tree_fixing.create_node("Chloropsis")
genus_ph1.add_feature("genus_name", "Chloropsis")
genus_ph3.tx_level = "genus"
genus_ph3.ancestral_rank = "genus"
genus_ph3.desc_rank = "genus"
genus_ph3.ph_tx = "PH"

genus_ph4 = tree_fixing.create_node("Irena")
genus_ph1.add_feature("genus_name", "Irena")
genus_ph4.tx_level = "genus"
genus_ph4.ancestral_rank = "genus"
genus_ph4.desc_rank = "genus"
genus_ph4.ph_tx = "PH"

family_ph2.add_child(genus_ph3)
family_ph2.add_child(genus_ph4)


species_ph4 = tree_fixing.create_node("Chloropsis sonnerati")
species_ph4.add_feature("genus_name", "Chloropsis")
species_ph4.add_feature("species_name", "Chloropsis sonnerati")
species_ph4.tx_level = "species"
species_ph4.ancestral_rank = "species"
species_ph4.desc_rank = "species"
species_ph4.ph_tx = "PH"

species_ph5 = tree_fixing.create_node("Chloropsis aurifrons")
species_ph5.add_feature("genus_name", "Chloropsis")
species_ph5.add_feature("species_name", "Chloropsis aurifrons")
species_ph5.tx_level = "species"
species_ph5.ancestral_rank = "species"
species_ph5.desc_rank = "species"
species_ph5.ph_tx = "PH"

genus_ph3.add_child(species_ph4)
genus_ph3.add_child(species_ph5)


species_ph6 = tree_fixing.create_node("Irena cyanogastra")
species_ph6.add_feature("genus_name", "Irena")
species_ph6.add_feature("species_name", "Irena cyanogastra")
species_ph6.tx_level = "species"
species_ph6.ancestral_rank = "species"
species_ph6.desc_rank = "species"
species_ph6.ph_tx = "PH"

species_ph7 = tree_fixing.create_node("Irena puella")
species_ph7.add_feature("genus_name", "Irena")
species_ph7.add_feature("species_name", "Irena puella")
species_ph7.tx_level = "species"
species_ph7.ancestral_rank = "species"
species_ph7.desc_rank = "species"
species_ph7.ph_tx = "PH"

genus_ph4.add_child(species_ph6)
genus_ph4.add_child(species_ph7)


species_ph8 = tree_fixing.create_node("Cymbilaimus lineatus")
species_ph8.add_feature("genus_name", "Cymbilaimus")
species_ph8.add_feature("species_name", "Cymbilaimus lineatus")
species_ph8.tx_level = "species"
species_ph8.ancestral_rank = "species"
species_ph8.desc_rank = "species"
species_ph8.ph_tx = "PH"

species_ph9 = tree_fixing.create_node("Cymbilaimus sanctaemariae")
species_ph9.add_feature("genus_name", "Cymbilaimus")
species_ph9.add_feature("species_name", "Cymbilaimus sanctaemariae")
species_ph9.tx_level = "species"
species_ph9.ancestral_rank = "species"
species_ph9.desc_rank = "species"
species_ph9.ph_tx = "PH"

genus_ph2.add_child(species_ph8)
genus_ph2.add_child(species_ph9)


mrca4 = tree_fixing.create_node("mrca4")
mrca4.ancestral_rank = "family"
mrca4.desc_rank = "species"
mrca4.ph_tx = "PH"

genus_ph5 = tree_fixing.create_node("Calendulauda")
genus_ph5.add_feature("genus_name", "Calendulauda")
genus_ph5.tx_level = "genus"
genus_ph5.ancestral_rank = "genus"
genus_ph5.desc_rank = "genus"
genus_ph5.ph_tx = "PH"

genus_tx8 = tree_fixing.create_node("Otocoris")
genus_tx8.add_feature("genus_name", "Otocoris")
genus_tx8.tx_level = "genus"
genus_tx8.ancestral_rank = "genus"
genus_tx8.desc_rank = "genus"
genus_tx8.ph_tx = "TX"

family_ph1.add_child(genus_ph5)
family_ph1.add_child(mrca4)
family_ph1.add_child(genus_tx8)

species_tx3 = tree_fixing.create_node("Otocoris berlepschi")
species_tx3.add_feature("genus_name", "Otocoris")
species_tx3.add_feature("species_name", "Otocoris berlepschi")
species_tx3.tx_level = "species"
species_tx3.ancestral_rank = "species"
species_tx3.desc_rank = "species"
species_tx3.ph_tx = "TX"

genus_tx8.add_child(species_tx3)


species_ph10 = tree_fixing.create_node("Calendulauda burra")
species_ph10.add_feature("genus_name", "Calendulauda burra")
species_ph10.add_feature("species_name", "Calendulauda burra")
species_ph10.tx_level = "species"
species_ph10.ancestral_rank = "species"
species_ph10.desc_rank = "species"
species_ph10.ph_tx = "PH"

species_ph11 = tree_fixing.create_node("Calendulauda erythroclamys")
species_ph11.add_feature("genus_name", "Calendulauda")
species_ph11.add_feature("species_name", "Calendulauda erythroclamys")
species_ph11.tx_level = "species"
species_ph11.ancestral_rank = "species"
species_ph11.desc_rank = "species"
species_ph11.ph_tx = "PH"

species_tx1 = tree_fixing.create_node("Calendulauda albescens")
species_tx1.add_feature("genus_name", "Calendulauda")
species_tx1.add_feature("species_name", "Calendulauda albescens")
species_tx1.tx_level = "species"
species_tx1.ancestral_rank = "species"
species_tx1.desc_rank = "species"
species_tx1.ph_tx = "TX"

species_tx2 = tree_fixing.create_node("Calendulauda sabota")
species_tx2.add_feature("genus_name", "Calendulauda")
species_tx2.add_feature("species_name", "Calendulauda sabota")
species_tx2.tx_level = "species"
species_tx2.ancestral_rank = "species"
species_tx2.desc_rank = "species"
species_tx2.ph_tx = "TX"

genus_ph5.add_child(species_ph10)
genus_ph5.add_child(species_ph11)
genus_ph5.add_child(species_tx1)
genus_ph5.add_child(species_tx2)


mrca5 = tree_fixing.create_node("mrca5")
mrca5.ancestral_rank = "family"
mrca5.desc_rank = "species"
mrca5.ph_tx = "PH"

mrca6 = tree_fixing.create_node("mrca6")
mrca6.ancestral_rank = "family"
mrca6.desc_rank = "species"
mrca6.ph_tx = "PH"

species_tx3 = tree_fixing.create_node("Mirafra pulpa")
species_tx3.add_feature("genus_name", "Mirafra")
species_tx3.add_feature("species_name", "Mirafra pulpa")
species_tx3.tx_level = "species"
species_tx3.ancestral_rank = "species"
species_tx3.desc_rank = "species"
species_tx3.ph_tx = "TX"

mrca4.add_child(mrca5)
mrca4.add_child(mrca6)
mrca4.add_child(species_tx3)


species_ph12 = tree_fixing.create_node("Mirafra stresemanni")
species_ph12.add_feature("genus_name", "Mirafra")
species_ph12.add_feature("species_name", "Mirafra stresemanni")
species_ph12.tx_level = "species"
species_ph12.ancestral_rank = "species"
species_ph12.desc_rank = "species"
species_ph12.ph_tx = "PH"

species_ph13 = tree_fixing.create_node("Mirafra cantillans")
species_ph13.add_feature("genus_name", "Mirafra")
species_ph13.add_feature("species_name", "Mirafra cantillans")
species_ph13.tx_level = "species"
species_ph13.ancestral_rank = "species"
species_ph13.desc_rank = "species"
species_ph13.ph_tx = "PH"

mrca5.add_child(species_ph12)
mrca5.add_child(species_ph13)


species_ph14 = tree_fixing.create_node("Mirafra hova")
species_ph14.add_feature("genus_name", "Mirafra")
species_ph14.add_feature("species_name", "Mirafra hova")
species_ph14.tx_level = "species"
species_ph14.ancestral_rank = "species"
species_ph14.desc_rank = "species"
species_ph14.ph_tx = "PH"

species_ph15 = tree_fixing.create_node("Eremopterix verticalis")
species_ph15.add_feature("genus_name", "Eremopterix")
species_ph15.add_feature("species_name", "Eremopterix verticalis")
species_ph15.tx_level = "species"
species_ph15.ancestral_rank = "species"
species_ph15.desc_rank = "species"
species_ph15.ph_tx = "PH"

mrca6.add_child(species_ph14)
mrca6.add_child(species_ph15)


tree_plotting.plot_labels(root, "labels_test.svg")

for leaf in root.get_leaves():
    leaf.date = 0

tree_dating.date_labelling(root)
tree_dating.impute_missing_dates(root)
tree_dating.compute_branch_lengths(root)
root.dist = 20

for node in root.traverse():
    node.adj = False

scale = 5

for leaf in root.get_leaves():
    node = leaf
    path_length = 0
    while node.up:
        path_length += 1
        node = node.up
    print(leaf.name, path_length)

    node = leaf
    while node.up:
        if not node.adj:
            if node.ph_tx == "TX":
                node.dist -= 15/scale
            else:
                node.dist -= 10/scale
            if len(node.up.children) > 1:
                node.dist -= 1/scale
            node.adj = True

        node = node.up

rng = np.random.default_rng(seed=2)

before = root.copy()

genus_dict_before = {}
nmp_genus_dict_before = {}
tree_labelling.populate_genus_dict(before, genus_dict_before, nmp_genus_dict_before, None)

tofix_dict_before = {}
tree_labelling.populate_tofix_dict(before, tofix_dict_before, nmp_genus_dict_before)
# tree_plotting.plot_simple(root, "figure3_1.svg", name_mrcas=False)
# tree_plotting.plot_simple(root, "figure3_1a.svg", name_mrcas=False, info_colors=False)

tree_labelling.populate_tofix_bkb(before, tofix_dict_before, [])

tree_plotting.plot_figure_fixing_b(before, "figure3_12.svg", name_mrcas=False)
tree_plotting.plot_figure_fixing_b(before, "figure3_12.tif", name_mrcas=False)


genus_dict = {}
nmp_genus_dict = {}
tree_labelling.populate_genus_dict(root, genus_dict, nmp_genus_dict, None)

tofix_dict = {}
tree_labelling.populate_tofix_dict(root, tofix_dict, nmp_genus_dict)

tree_fixing.fix_polyphyly(genus_dict, rng)
tree_fixing.fix_polyphyly(nmp_genus_dict, rng)

tree_labelling.populate_tofix_bkb(root, tofix_dict, [])
fix_dict = tree_labelling.process_tofix_bkb(tofix_dict)
tree_fixing.fix_polyphyly(fix_dict, rng, expand_parent_backbones=True)

tree_fixing.fix_all_polytomies(root, rng)

for node in root.traverse():
    node.date = None
    node.imputed_date = False

root.date = 100

for leaf in root.get_leaves():
    leaf.date = 0

tree_dating.date_labelling(root)
tree_dating.impute_missing_dates(root)
tree_dating.compute_branch_lengths(root)
root.dist = 20

for node in root.traverse():
    node.adj = False

scale = 5

for leaf in root.get_leaves():
    node = leaf
    path_length = 0
    while node.up:
        path_length += 1
        node = node.up
    print(leaf.name, path_length)

    node = leaf
    while node.up:
        if not node.adj:
            if node.ph_tx == "TX":
                node.dist -= 15/scale
            else:
                node.dist -= 10/scale
            if len(node.up.children) > 1:
                node.dist -= 1/scale
            node.adj = True

        node = node.up

tree_plotting.plot_figure_fixing_b(root, "figure3_13.svg", name_mrcas=False, arrows=False)
tree_plotting.plot_figure_fixing_b(root, "figure3_13.tif", name_mrcas=False, arrows=False)
