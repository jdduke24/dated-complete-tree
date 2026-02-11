import os
os.chdir('/Users/Jonathan/Library/CloudStorage/Dropbox/Imperial/Tree_of_Life/Open_Tree/python/dated-complete-tree')

import sys
import gc

import tree_loading
import tree_labelling
import tree_fixing
import tree_dating
import tree_plotting

import numpy as np


root = tree_fixing.create_node("Passeriformes")
root.props["tx_level"] = "order"
root.props["ancestral_rank"] = "order"
root.props["desc_rank"] = "order"
root.props["ph_tx"] = "PH"
root.props["date"] = 100.

family_ph1 = tree_fixing.create_node("Alaudidae")
family_ph1.props["tx_level"] = "family"
family_ph1.props["ancestral_rank"] = "family"
family_ph1.props["desc_rank"] = "family"
family_ph1.props["ph_tx"] = "PH"

family_tx = tree_fixing.create_node("Colluricinclidae")
family_tx.props["tx_level"] = "family"
family_tx.props["ancestral_rank"] = "family"
family_tx.props["desc_rank"] = "family"
family_tx.props["ph_tx"] = "TX"

genus_tx9 = tree_fixing.create_node("Laphyctes")
genus_tx9.add_prop("genus_name", "Laphyctes")
genus_tx9.props["tx_level"] = "genus"
genus_tx9.props["ancestral_rank"] = "genus"
genus_tx9.props["desc_rank"] = "genus"
genus_tx9.props["ph_tx"] = "TX"

mrca1 = tree_fixing.create_node("mrca1")
mrca1.props["ancestral_rank"] = "order"
mrca1.props["desc_rank"] = "family"
mrca1.props["ph_tx"] = "PH"

root.add_child(family_ph1)
root.add_child(mrca1)
root.add_child(family_tx)
root.add_child(genus_tx9)



species_tx4 = tree_fixing.create_node("Laphyctes apolites")
species_tx4.add_prop("genus_name", "Laphyctes")
species_tx4.add_prop("species_name", "Laphyctes apolites")
species_tx4.props["tx_level"] = "species"
species_tx4.props["ancestral_rank"] = "species"
species_tx4.props["desc_rank"] = "species"
species_tx4.props["ph_tx"] = "TX"

genus_tx9.add_child(species_tx4)


species_tx5 = tree_fixing.create_node("Rhectes analogus")
species_tx5.add_prop("genus_name", "Rhectes")
species_tx5.add_prop("species_name", "Rhectes analogus")
species_tx5.props["tx_level"] = "species"
species_tx5.props["ancestral_rank"] = "species"
species_tx5.props["desc_rank"] = "species"
species_tx5.props["ph_tx"] = "TX"

species_tx6 = tree_fixing.create_node("Rhectes rubiensis")
species_tx6.add_prop("genus_name", "Rhectes")
species_tx6.add_prop("species_name", "Rhectes rubiensis")
species_tx6.props["tx_level"] = "species"
species_tx6.props["ancestral_rank"] = "species"
species_tx6.props["desc_rank"] = "species"
species_tx6.props["ph_tx"] = "TX"

species_tx7 = tree_fixing.create_node("Rhectes ferrugineus")
species_tx7.add_prop("genus_name", "Rhectes")
species_tx7.add_prop("species_name", "Rhectes ferrugineus")
species_tx7.props["tx_level"] = "species"
species_tx7.props["ancestral_rank"] = "species"
species_tx7.props["desc_rank"] = "species"
species_tx7.props["ph_tx"] = "TX"

family_tx.add_child(species_tx5)
family_tx.add_child(species_tx7)
family_tx.add_child(species_tx6)


mrca2 = tree_fixing.create_node("mrca2")
mrca2.props["ancestral_rank"] = "order"
mrca2.props["desc_rank"] = "family"
mrca2.props["ph_tx"] = "PH"

genus_ph1 = tree_fixing.create_node("Thamnomanes")
genus_ph1.add_prop("genus_name", "Thamnomanes")
genus_ph1.props["tx_level"] = "genus"
genus_ph1.props["ancestral_rank"] = "genus"
genus_ph1.props["desc_rank"] = "genus"
genus_ph1.props["ph_tx"] = "PH"

mrca1.add_child(mrca2)
mrca1.add_child(genus_ph1)


species_ph1 = tree_fixing.create_node("Thamnomanes ardesiacus")
species_ph1.add_prop("genus_name", "Thamnomanes")
species_ph1.add_prop("species_name", "Thamnomanes ardesiacus")
species_ph1.props["tx_level"] = "species"
species_ph1.props["ancestral_rank"] = "species"
species_ph1.props["desc_rank"] = "species"
species_ph1.props["ph_tx"] = "PH"

species_ph2 = tree_fixing.create_node("Thamnomanes saturninus")
species_ph2.add_prop("genus_name", "Thamnomanes")
species_ph2.add_prop("species_name", "Thamnomanes saturninus")
species_ph2.props["tx_level"] = "species"
species_ph2.props["ancestral_rank"] = "species"
species_ph2.props["desc_rank"] = "species"
species_ph2.props["ph_tx"] = "PH"

genus_ph1.add_child(species_ph1)
genus_ph1.add_child(species_ph2)


mrca3 = tree_fixing.create_node("mrca3")
mrca3.props["ancestral_rank"] = "order"
mrca3.props["desc_rank"] = "family"
mrca3.props["ph_tx"] = "PH"

species_ph3 = tree_fixing.create_node("Garrulax virgatus")
species_ph3.add_prop("genus_name", "Garrulax")
species_ph3.add_prop("species_name", "Garrulax virgatus")
species_ph3.props["tx_level"] = "species"
species_ph3.props["ancestral_rank"] = "species"
species_ph3.props["desc_rank"] = "species"
species_ph3.props["ph_tx"] = "PH"

mrca2.add_child(mrca3)
mrca2.add_child(species_ph3)


family_ph2 = tree_fixing.create_node("Irenidae")
family_ph2.props["tx_level"] = "family"
family_ph2.props["ancestral_rank"] = "family"
family_ph2.props["desc_rank"] = "family"
family_ph2.props["ph_tx"] = "PH"

genus_ph2 = tree_fixing.create_node("Cymbilaimus")
genus_ph1.add_prop("genus_name", "Cymbilaimus")
genus_ph2.props["tx_level"] = "genus"
genus_ph2.props["ancestral_rank"] = "genus"
genus_ph2.props["desc_rank"] = "genus"
genus_ph2.props["ph_tx"] = "PH"

mrca3.add_child(family_ph2)
mrca3.add_child(genus_ph2)


genus_ph3 = tree_fixing.create_node("Chloropsis")
genus_ph1.add_prop("genus_name", "Chloropsis")
genus_ph3.props["tx_level"] = "genus"
genus_ph3.props["ancestral_rank"] = "genus"
genus_ph3.props["desc_rank"] = "genus"
genus_ph3.props["ph_tx"] = "PH"

genus_ph4 = tree_fixing.create_node("Irena")
genus_ph1.add_prop("genus_name", "Irena")
genus_ph4.props["tx_level"] = "genus"
genus_ph4.props["ancestral_rank"] = "genus"
genus_ph4.props["desc_rank"] = "genus"
genus_ph4.props["ph_tx"] = "PH"

family_ph2.add_child(genus_ph3)
family_ph2.add_child(genus_ph4)


species_ph4 = tree_fixing.create_node("Chloropsis sonnerati")
species_ph4.add_prop("genus_name", "Chloropsis")
species_ph4.add_prop("species_name", "Chloropsis sonnerati")
species_ph4.props["tx_level"] = "species"
species_ph4.props["ancestral_rank"] = "species"
species_ph4.props["desc_rank"] = "species"
species_ph4.props["ph_tx"] = "PH"

species_ph5 = tree_fixing.create_node("Chloropsis aurifrons")
species_ph5.add_prop("genus_name", "Chloropsis")
species_ph5.add_prop("species_name", "Chloropsis aurifrons")
species_ph5.props["tx_level"] = "species"
species_ph5.props["ancestral_rank"] = "species"
species_ph5.props["desc_rank"] = "species"
species_ph5.props["ph_tx"] = "PH"

genus_ph3.add_child(species_ph4)
genus_ph3.add_child(species_ph5)


species_ph6 = tree_fixing.create_node("Irena cyanogastra")
species_ph6.add_prop("genus_name", "Irena")
species_ph6.add_prop("species_name", "Irena cyanogastra")
species_ph6.props["tx_level"] = "species"
species_ph6.props["ancestral_rank"] = "species"
species_ph6.props["desc_rank"] = "species"
species_ph6.props["ph_tx"] = "PH"

species_ph7 = tree_fixing.create_node("Irena puella")
species_ph7.add_prop("genus_name", "Irena")
species_ph7.add_prop("species_name", "Irena puella")
species_ph7.props["tx_level"] = "species"
species_ph7.props["ancestral_rank"] = "species"
species_ph7.props["desc_rank"] = "species"
species_ph7.props["ph_tx"] = "PH"

genus_ph4.add_child(species_ph6)
genus_ph4.add_child(species_ph7)


species_ph8 = tree_fixing.create_node("Cymbilaimus lineatus")
species_ph8.add_prop("genus_name", "Cymbilaimus")
species_ph8.add_prop("species_name", "Cymbilaimus lineatus")
species_ph8.props["tx_level"] = "species"
species_ph8.props["ancestral_rank"] = "species"
species_ph8.props["desc_rank"] = "species"
species_ph8.props["ph_tx"] = "PH"

species_ph9 = tree_fixing.create_node("Cymbilaimus sanctaemariae")
species_ph9.add_prop("genus_name", "Cymbilaimus")
species_ph9.add_prop("species_name", "Cymbilaimus sanctaemariae")
species_ph9.props["tx_level"] = "species"
species_ph9.props["ancestral_rank"] = "species"
species_ph9.props["desc_rank"] = "species"
species_ph9.props["ph_tx"] = "PH"

genus_ph2.add_child(species_ph8)
genus_ph2.add_child(species_ph9)


mrca4 = tree_fixing.create_node("mrca4")
mrca4.props["ancestral_rank"] = "family"
mrca4.props["desc_rank"] = "species"
mrca4.props["ph_tx"] = "PH"

genus_ph5 = tree_fixing.create_node("Calendulauda")
genus_ph5.add_prop("genus_name", "Calendulauda")
genus_ph5.props["tx_level"] = "genus"
genus_ph5.props["ancestral_rank"] = "genus"
genus_ph5.props["desc_rank"] = "genus"
genus_ph5.props["ph_tx"] = "PH"

genus_tx8 = tree_fixing.create_node("Otocoris")
genus_tx8.add_prop("genus_name", "Otocoris")
genus_tx8.props["tx_level"] = "genus"
genus_tx8.props["ancestral_rank"] = "genus"
genus_tx8.props["desc_rank"] = "genus"
genus_tx8.props["ph_tx"] = "TX"

family_ph1.add_child(genus_ph5)
family_ph1.add_child(mrca4)
family_ph1.add_child(genus_tx8)

species_tx3 = tree_fixing.create_node("Otocoris berlepschi")
species_tx3.add_prop("genus_name", "Otocoris")
species_tx3.add_prop("species_name", "Otocoris berlepschi")
species_tx3.props["tx_level"] = "species"
species_tx3.props["ancestral_rank"] = "species"
species_tx3.props["desc_rank"] = "species"
species_tx3.props["ph_tx"] = "TX"

genus_tx8.add_child(species_tx3)


species_ph10 = tree_fixing.create_node("Calendulauda burra")
species_ph10.add_prop("genus_name", "Calendulauda burra")
species_ph10.add_prop("species_name", "Calendulauda burra")
species_ph10.props["tx_level"] = "species"
species_ph10.props["ancestral_rank"] = "species"
species_ph10.props["desc_rank"] = "species"
species_ph10.props["ph_tx"] = "PH"

species_ph11 = tree_fixing.create_node("Calendulauda erythroclamys")
species_ph11.add_prop("genus_name", "Calendulauda")
species_ph11.add_prop("species_name", "Calendulauda erythroclamys")
species_ph11.props["tx_level"] = "species"
species_ph11.props["ancestral_rank"] = "species"
species_ph11.props["desc_rank"] = "species"
species_ph11.props["ph_tx"] = "PH"

species_tx1 = tree_fixing.create_node("Calendulauda albescens")
species_tx1.add_prop("genus_name", "Calendulauda")
species_tx1.add_prop("species_name", "Calendulauda albescens")
species_tx1.props["tx_level"] = "species"
species_tx1.props["ancestral_rank"] = "species"
species_tx1.props["desc_rank"] = "species"
species_tx1.props["ph_tx"] = "TX"

species_tx2 = tree_fixing.create_node("Calendulauda sabota")
species_tx2.add_prop("genus_name", "Calendulauda")
species_tx2.add_prop("species_name", "Calendulauda sabota")
species_tx2.props["tx_level"] = "species"
species_tx2.props["ancestral_rank"] = "species"
species_tx2.props["desc_rank"] = "species"
species_tx2.props["ph_tx"] = "TX"

genus_ph5.add_child(species_ph10)
genus_ph5.add_child(species_ph11)
genus_ph5.add_child(species_tx1)
genus_ph5.add_child(species_tx2)


mrca5 = tree_fixing.create_node("mrca5")
mrca5.props["ancestral_rank"] = "family"
mrca5.props["desc_rank"] = "species"
mrca5.props["ph_tx"] = "PH"

mrca6 = tree_fixing.create_node("mrca6")
mrca6.props["ancestral_rank"] = "family"
mrca6.props["desc_rank"] = "species"
mrca6.props["ph_tx"] = "PH"

species_tx3 = tree_fixing.create_node("Mirafra pulpa")
species_tx3.add_prop("genus_name", "Mirafra")
species_tx3.add_prop("species_name", "Mirafra pulpa")
species_tx3.props["tx_level"] = "species"
species_tx3.props["ancestral_rank"] = "species"
species_tx3.props["desc_rank"] = "species"
species_tx3.props["ph_tx"] = "TX"

mrca4.add_child(mrca5)
mrca4.add_child(mrca6)
mrca4.add_child(species_tx3)


species_ph12 = tree_fixing.create_node("Mirafra stresemanni")
species_ph12.add_prop("genus_name", "Mirafra")
species_ph12.add_prop("species_name", "Mirafra stresemanni")
species_ph12.props["tx_level"] = "species"
species_ph12.props["ancestral_rank"] = "species"
species_ph12.props["desc_rank"] = "species"
species_ph12.props["ph_tx"] = "PH"

species_ph13 = tree_fixing.create_node("Mirafra cantillans")
species_ph13.add_prop("genus_name", "Mirafra")
species_ph13.add_prop("species_name", "Mirafra cantillans")
species_ph13.props["tx_level"] = "species"
species_ph13.props["ancestral_rank"] = "species"
species_ph13.props["desc_rank"] = "species"
species_ph13.props["ph_tx"] = "PH"

mrca5.add_child(species_ph12)
mrca5.add_child(species_ph13)


species_ph14 = tree_fixing.create_node("Mirafra hova")
species_ph14.add_prop("genus_name", "Mirafra")
species_ph14.add_prop("species_name", "Mirafra hova")
species_ph14.props["tx_level"] = "species"
species_ph14.props["ancestral_rank"] = "species"
species_ph14.props["desc_rank"] = "species"
species_ph14.props["ph_tx"] = "PH"

species_ph15 = tree_fixing.create_node("Eremopterix verticalis")
species_ph15.add_prop("genus_name", "Eremopterix")
species_ph15.add_prop("species_name", "Eremopterix verticalis")
species_ph15.props["tx_level"] = "species"
species_ph15.props["ancestral_rank"] = "species"
species_ph15.props["desc_rank"] = "species"
species_ph15.props["ph_tx"] = "PH"

mrca6.add_child(species_ph14)
mrca6.add_child(species_ph15)


# tree_plotting.plot_labels(root, "labels_test.svg")

for leaf in root.leaves():
    leaf.props["date"] = 0

for node in root.traverse():
    if "date" not in node.props:
        node.add_prop("date", None)

tree_dating.date_labelling(root)
tree_dating.impute_missing_dates(root)
tree_dating.compute_branch_lengths(root)
root.dist = 20

for node in root.traverse():
    node.add_prop("adj", False)

scale = 6

before = root.copy()

for leaf in before.leaves():
    node = leaf
    path_length = 0
    while node.up:
        path_length += 1
        node = node.up
    print(leaf.name, path_length)

    node = leaf
    while node.up:
        if not node.props["adj"]:
            if node.props["ph_tx"] == "TX":
                node.dist -= 15/scale
            else:
                node.dist -= 10/scale
            if len(node.up.children) > 1:
                node.dist -= 1/scale
            node.props["adj"] = True

        node = node.up

rng = np.random.default_rng(seed=2)

genus_dict_before = {}
nmp_genus_dict_before = {}
tree_labelling.populate_genus_dict(before, genus_dict_before, nmp_genus_dict_before, None)

tofix_dict_before = {}
tree_labelling.populate_tofix_dict(before, tofix_dict_before, nmp_genus_dict_before)

tree_labelling.populate_tofix_bkb(before, tofix_dict_before, [])

tree_plotting.plot_figure_fixing_b(before, "figures/figure6/figure6_before.svg", name_mrcas=False, scale=scale)
tree_plotting.plot_figure_fixing_b(before, "figures/figure6/figure6_before.tif", name_mrcas=False, scale=scale)


##########################################################

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
    node.props["date"] = None
    node.props["imputed_date"] = False

root.props["date"] = 100

for leaf in root.leaves():
    leaf.props["date"] = 0

tree_dating.date_labelling(root)
tree_dating.impute_missing_dates(root)
tree_dating.compute_branch_lengths(root)
root.dist = 20

for node in root.traverse():
    node.props["adj"] = False

for leaf in root.leaves():
    node = leaf
    path_length = 0
    while node.up:
        path_length += 1
        node = node.up
    print(leaf.name, path_length)

    node = leaf
    while node.up:
        if not node.props["adj"]:
            if node.props["ph_tx"] == "TX":
                node.dist -= 15/scale
            else:
                node.dist -= 10/scale
            if len(node.up.children) > 1:
                node.dist -= 1/scale
            node.props["adj"] = True

        node = node.up

tree_plotting.plot_figure_fixing_b(root, "figures/figure6/figure6_after.svg", name_mrcas=False, arrows=False, scale=scale)
tree_plotting.plot_figure_fixing_b(root, "figures/figure6/figure6_after.tif", name_mrcas=False, arrows=False, scale=scale)
