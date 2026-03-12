#!/usr/bin/env python3

# import ete3

from PyQt6.QtGui import QFont
# Add the missing attribute that ETE is looking for
if not hasattr(QFont, "StyleItalic"):
    QFont.StyleItalic = QFont.Style.StyleItalic

import ete4
import numpy as np
from taxonomy_utils import tx_levels
import tree_fixing


# create string from a list of dates, rounding the dates to 1 d.p.
def date_tuple_to_str(date_tuple):
    rounded_dates = "%.1f" % date_tuple[2][0]
    for d in date_tuple[2][1:]:
        rounded_dates += ", %.1f" % d
    return "[%s]" % rounded_dates


def routes_to_str(routes):
    return_str = ""
    for age, dist in routes:
        if age == 0:
            return_str += "(0, %d) / " % (dist)
        else:
            return_str += "(%.1f, %d) / " % (age, dist)
    return return_str[:-3]


def oldest_path_to_str(oldest_path):
    age, dist = oldest_path
    if age == 0:
        return "(0, %d)" % (dist)
    else:
        return "(%.2f, %d)" % (age, dist)


def name_to_simple_name(name):
    parts = name.split("_")
    ret_str = parts[0]
    if len(parts) > 0:
        for p in parts[1:-1]:
            ret_str = ret_str + " " + p
    return ret_str.strip()


def plot_evoh_tree(input_tre, filename):
    tre = input_tre.copy()
    # plot tree
    for node in tre.traverse(strategy="preorder"):

        nstyle = ete4.treeview.NodeStyle()
        nstyle["vt_line_width"] = 1
        nstyle["hz_line_width"] = 1
        nstyle["vt_line_color"] = "gray"
        nstyle["hz_line_color"] = "gray"

        nstyle["size"] = 0
        node.set_style(nstyle)

        node.add_face(ete4.treeview.TextFace("p=%.3f" % (node.props["evoh_p"]), fsize=8), column=0, position="branch-top")
        node.add_face(ete4.treeview.TextFace("bl=%.3f" % (node.props["evoh_bl"]), fsize=8), column=0, position="branch-bottom")

        if node.is_leaf:
            node.add_face(ete4.treeview.TextFace("     %s, pext=%.2f" % (node.name, node.props["pext"]), fsize=8), column=0, position="branch-right")

    new_node = ete4.Tree()
    new_node.add_child(tre)
    new_node.dist = 3
    nstyle = ete4.treeview.NodeStyle()
    nstyle["hz_line_color"] = "white"
    nstyle["hz_line_width"] = 0
    new_node.set_style(nstyle)

    new_node.add_face(ete4.treeview.TextFace("rho=%.3f        " % (tre.props["rho"]), fsize=8), column=0, position="branch-top")
    new_node.add_face(ete4.treeview.TextFace("EvoH=%.3f        " % (tre.props["phi_rho"]), fsize=8), column=0, position="branch-bottom")

    ts = ete4.treeview.TreeStyle()
    ts.show_leaf_name = False
    ts.mode ="r"
    ts.branch_vertical_margin = 60
    ts.margin_right = 100
    ts.show_scale = False

    ts.scale = 50

    new_node.render(filename, tree_style=ts)


def plot_labels(input_tre, filename):
    tre = input_tre.copy()
    # plot tree
    for node in tre.traverse(strategy="preorder"):

        nstyle = ete4.treeview.NodeStyle()
        nstyle["vt_line_width"] = 1
        nstyle["hz_line_width"] = 1
        nstyle["vt_line_color"] = "gray"
        nstyle["hz_line_color"] = "gray"

        nstyle["size"] = 15
        node.set_style(nstyle)

        if node.props["info"]:
            if "GR FIX" in node.props["info"]:
                nstyle["fgcolor"] = "limegreen"
            elif "GR BKB" in node.props["info"]:
                nstyle["fgcolor"] = "greenyellow"
            elif "OTH FIX" in node.props["info"]:
                nstyle["fgcolor"] = "mediumvioletred"
                nstyle["shape"] = "sphere"
            elif "OTH BKB 1" in node.props["info"]:
                nstyle["fgcolor"] = "violet"
            elif "OTH BKB 2" in node.props["info"]:
                nstyle["fgcolor"] = "mediumorchid"
            elif "OTH BKB 3" in node.props["info"]:
                nstyle["fgcolor"] = "darkviolet"
            elif "OTH EX BKB" in node.props["info"]:
                nstyle["fgcolor"] = "darkmagenta"
            elif "NMP FIX" in node.props["info"]:
                nstyle["fgcolor"] = "deepskyblue"
            elif "NMP BKB" in node.props["info"]:
                nstyle["fgcolor"] = "turquoise"
            elif "NMP EX BKB" in node.props["info"]:
                nstyle["fgcolor"] = "royalblue"
        else:
            nstyle["fgcolor"] = "grey"

        if "date" in node.props:
            node.add_face(ete4.treeview.TextFace("Name: %s\nRank: %s, %s, %s" % (node.name, node.props["tx_level"], node.props["ph_tx"], str(node.props["date"]) if node.props["date"] and not node.props["imputed_date"] else "")), column=0, position="branch-top")
        else:
            node.add_face(ete4.treeview.TextFace("Name: %s\nRank: %s, %s" % (node.name, node.props["tx_level"], node.props["ph_tx"])), column=0, position="branch-top")
        node.add_face(ete4.treeview.TextFace("Anc: %s\nDesc: %s\n%s" % (node.props["ancestral_rank"], node.props["desc_rank"], node.props["info"])), column=0, position="branch-bottom")

    ts = ete4.treeview.TreeStyle()
    ts.show_leaf_name = False
    ts.mode ="r"
    ts.branch_vertical_margin = 8
    ts.show_scale = False

    ts.scale = 40

    tre.render(filename, tree_style=ts)


def plot_figure_fixing(input_tre, filename, name_mrcas=True, info_colors=True, arrows=True):
    tre = input_tre.copy()

    arrow0 = ete3.ImgFace("arrow0.png", width=10, height=15)
    arrow1 = ete3.ImgFace("arrow_new_1.png", width=10, height=15)
    arrow2 = ete3.ImgFace("arrow_new_2.png", width=10, height=15)
    arrow3 = ete3.ImgFace("arrow_new_3.png", width=10, height=15)
    arrow4 = ete3.ImgFace("arrow_new_4.png", width=10, height=15)
    arrow5 = ete3.ImgFace("arrow_new_5.png", width=10, height=15)

    # plot tree
    for node in tre.traverse(strategy="preorder"):

        nstyle = ete3.NodeStyle()
        nstyle["vt_line_width"] = 1
        nstyle["hz_line_width"] = 1
        nstyle["vt_line_color"] = "grey"
        nstyle["hz_line_color"] = "grey"

        if node.ph_tx == "TX":
            nstyle["shape"] = "square"
            nstyle["size"] = 15
        else:
            nstyle["shape"] = "circle"
            nstyle["size"] = 10

        if node.info and node.info != "OTH PARENT":
            if "GR FIX" in node.info and "Calendulauda" in node.name:
                nstyle["fgcolor"] = "#65B6F9"
            elif "OTH FIX" in node.info and "Otocoris" in node.name:
                nstyle["fgcolor"] = "#12083E"
            elif "OTH FIX" in node.info and "Colluricinclidae" in node.name:
                nstyle["fgcolor"] = "#FFD545"
            elif "OTH FIX" in node.info and "Laphyctes" in node.name:
                nstyle["fgcolor"] = "#EC803C"
            elif "NMP FIX" in node.info and "Rhectes" not in node.name:
                nstyle["fgcolor"] = "#DC2683"
            else:
                nstyle["fgcolor"] = "#E0E0E0"

            if arrows:
                if "GR BKB Calendulauda" in node.info:
                    node.add_face(arrow1, 0, position="branch-top")
                else:
                    node.add_face(arrow0, 0, position="branch-top")

                if ("OTH BKB" in node.info or "OTH EX BKB" in node.info) and "Alaudidae genus" in node.info:
                    node.add_face(arrow2, 1, position="branch-top")
                else:
                    node.add_face(arrow0, 1, position="branch-top")

                if ("NMP BKB" in node.info or "NMP EX BKB" in node.info) and "Mirafra" in node.info:
                    node.add_face(arrow3, 2, position="branch-top")
                else:
                    node.add_face(arrow0, 2, position="branch-top")

                if ("OTH BKB" in node.info or "OTH EX BKB" in node.info) and "Passeriformes family" in node.info:
                    node.add_face(arrow4, 3, position="branch-top")
                else:
                    node.add_face(arrow0, 3, position="branch-top")

                if ("OTH BKB" in node.info or "OTH EX BKB" in node.info) and "Passeriformes genus" in node.info:
                    node.add_face(arrow5, 4, position="branch-top")
                else:
                    node.add_face(arrow0, 4, position="branch-top")

                node.add_face(arrow0, 5, position="branch-top")
            else:
                for i in range(6):
                    node.add_face(arrow0, i, position="branch-top")

        else:
            nstyle["fgcolor"] = "#E0E0E0"


        node.set_style(nstyle)

        if not (tx_levels[node.tx_level] == tx_levels["mrca"] and not name_mrcas):
            if tx_levels[node.tx_level] == tx_levels["species"]:
                node.add_face(ete3.TextFace(name_to_simple_name(node.name), fstyle="italic", fsize=10), column=0, position="branch-right")
            else:
                # node.add_face(ete3.TextFace(name_to_simple_name(node.name), fsize=9), column=0, position="branch-top")
                node.add_face(ete3.TextFace("%s\n%s" % (name_to_simple_name(node.name), node.tx_level), fsize=9), column=0, position="branch-bottom")
                # node.add_face(ete3.TextFace("%s" % (node.tx_level), fsize=9), column=0, position="branch-bottom")


    ts = ete3.TreeStyle()
    ts.margin_right = 80
    ts.show_leaf_name = False
    ts.mode ="r"
    if arrows:
        ts.branch_vertical_margin = -5
    else:
        ts.branch_vertical_margin = -5
    ts.scale = 100
    ts.show_scale = False

    tre.render(filename, tree_style=ts)


def plot_figure_fixing_b(input_tre, filename, name_mrcas=True, info_colors=True, arrows=True, scale=5):
    tre = input_tre.copy()

    # plot tree
    for node in tre.traverse(strategy="preorder"):

        nstyle = ete4.treeview.NodeStyle()
        nstyle["vt_line_width"] = 1
        nstyle["hz_line_width"] = 1
        nstyle["vt_line_color"] = "grey"
        nstyle["hz_line_color"] = "grey"

        if arrows and node.props["ph_tx"] == "TX":
            nstyle["hz_line_type"] = 1
            nstyle["vt_line_type"] = 1
        else:
            nstyle["hz_line_type"] = 0
            nstyle["vt_line_type"] = 0

        if node.props["ph_tx"] == "TX":
            nstyle["shape"] = "square"
            nstyle["size"] = 15
        else:
            nstyle["shape"] = "circle"
            nstyle["size"] = 10

        colours = ["#65B6F9",
                   "#12083E",
                   "#FFA351",
                   "#DC2683",
                   "#FFD545"]

        spacing = 0
        if node.name == "Chloropsis":
            spacing = 2
        elif node.name == "Otocoris":
            spacing = 3
        elif node.name == "Irenidae":
            spacing = 4
        elif node.name == "Irena":
            spacing = 6


        if node is tre or (node.props["info"] and node.props["info"] != "OTH PARENT"):
            if "GR FIX" in node.props["info"] and "Calendulauda" in node.name:
                nstyle["fgcolor"] = colours[0]
            elif "OTH FIX" in node.props["info"] and "Otocoris" in node.name:
                nstyle["fgcolor"] = colours[1]
            elif "OTH FIX" in node.props["info"] and "Colluricinclidae" in node.name:
                nstyle["fgcolor"] = colours[4]
            elif "OTH FIX" in node.props["info"] and "Laphyctes" in node.name:
                nstyle["fgcolor"] = colours[2]
            elif "NMP FIX" in node.props["info"] and "Rhectes" not in node.name:
                nstyle["fgcolor"] = colours[3]
            else:
                nstyle["fgcolor"] = "#BBBBBB"

            if arrows:
                if "GR BKB Calendulauda" in node.props["info"] or node.name == "Calendulauda":
                    node.add_face(ete4.treeview.TextFace("B ", fgcolor=colours[0], fsize=10, bold=True), 1, position="branch-top")
                else:
                    node.add_face(ete4.treeview.TextFace("B ", fgcolor="white", fsize=10, bold=True), 1, position="branch-top")

                if (("OTH BKB" in node.props["info"] or "OTH EX BKB" in node.props["info"]) and "Alaudidae genus" in node.props["info"]) or node.name == "Alaudidae":
                    node.add_face(ete4.treeview.TextFace("C ", fgcolor=colours[1], fsize=10, bold=True), 2, position="branch-top")
                else:
                    node.add_face(ete4.treeview.TextFace("C ", fgcolor="white", fsize=10, bold=True), 2, position="branch-top")

                if ("NMP BKB" in node.props["info"] or "NMP EX BKB" in node.props["info"]) and "Mirafra" in node.props["info"]:
                    node.add_face(ete4.treeview.TextFace("A ", fgcolor=colours[3], fsize=10, bold=True), 0, position="branch-top")
                else:
                    node.add_face(ete4.treeview.TextFace("A ", fgcolor="white", fsize=10, bold=True), 0, position="branch-top")

                if (("OTH BKB" in node.props["info"] or "OTH EX BKB" in node.props["info"]) and "Passeriformes family" in node.props["info"]) or node.props["tx_level"] == "order":
                    node.add_face(ete4.treeview.TextFace("D ", fgcolor=colours[4], fsize=10, bold=True), 3, position="branch-top")
                else:
                    node.add_face(ete4.treeview.TextFace("D ", fgcolor="white", fsize=10, bold=True), 3, position="branch-top")

                if (("OTH BKB" in node.props["info"] or "OTH EX BKB" in node.props["info"]) and "Passeriformes genus" in node.props["info"]) or node.props["tx_level"] == "order":
                    node.add_face(ete4.treeview.TextFace("E ", fgcolor=colours[2], fsize=10, bold=True), 4, position="branch-top")
                else:
                    node.add_face(ete4.treeview.TextFace("E ", fgcolor="white", fsize=10, bold=True), 4, position="branch-top")

                if node.props["ph_tx"] == "TX":
                    if node.name[:12] == "Calendulauda":
                        node.add_face(ete4.treeview.TextFace("b ", fgcolor=colours[0], fsize=10, bold=True), 5, position="branch-top")
                    if node.name[:8] == "Otocoris":
                        node.add_face(ete4.treeview.TextFace("c ", fgcolor=colours[1], fsize=10, bold=True), 5, position="branch-top")
                    if node.name[:7] == "Mirafra":
                        node.add_face(ete4.treeview.TextFace("a ", fgcolor=colours[3], fsize=10, bold=True), 5, position="branch-top")
                    if node.name[:16] == "Colluricinclidae":
                        node.add_face(ete4.treeview.TextFace("d ", fgcolor=colours[4], fsize=10, bold=True), 5, position="branch-top")
                    if node.name[:9] == "Laphyctes":
                        node.add_face(ete4.treeview.TextFace("e ", fgcolor=colours[2], fsize=10, bold=True), 5, position="branch-top")

                else:
                    node.add_face(ete4.treeview.TextFace("A ", fgcolor="white", fsize=10, bold=True), 5, position="branch-top")

            else:
                for i in range(6):
                    node.add_face(ete4.treeview.TextFace("A ", fgcolor="white", fsize=10, bold=True), i, position="branch-top")

        else:
            nstyle["fgcolor"] = "#BBBBBB"


        node.set_style(nstyle)

        if not (tx_levels[node.props["tx_level"]] == tx_levels["mrca"] and not name_mrcas):
            if tx_levels[node.props["tx_level"]] == tx_levels["species"]:
                # node.add_face(ete4.treeview.TextFace(" " + name_to_simple_name(node.name), style={"fstyle": "italic", "fsize": 10}), column=0, position="branch-right")
                from ete4.treeview import TextFace
                # tf = TextFace(" " + name_to_simple_name(node.name), style={"fstyle": "italic", "fsize": 10})
                tf = TextFace(" " + name_to_simple_name(node.name), fstyle="italic")
                node.add_face(tf, column=0, position="branch-right")
            else:
                # node.add_face(ete4.treeview.TextFace(name_to_simple_name(node.name), fsize=9), column=0, position="branch-top")
                node.add_face(ete4.treeview.TextFace((' ' * spacing) +
                                            "%s\n" % name_to_simple_name(node.name) +
                                            (' ' * spacing) +
                                            "%s" % node.props["tx_level"], fsize=9), column=0, position="branch-bottom")
                # node.add_face(ete4.treeview.TextFace("%s" % (node.props["tx_level"]), fsize=9), column=0, position="branch-bottom")


    ts = ete4.treeview.TreeStyle()
    ts.margin_right = 80
    ts.show_leaf_name = False
    ts.mode ="r"
    if arrows:
        ts.branch_vertical_margin = -7
    else:
        ts.branch_vertical_margin = -7
    ts.scale = scale
    ts.show_scale = False

    tre.render(filename, tree_style=ts, units="in", w=6, dpi=300)


def plot_simple(input_tre, filename, name_mrcas=True, info_colors=True, color_GR=True, color_NMP=True, color_OTH=True, fix_node=None, bkb_label=None):
    tre = input_tre.copy()
    # plot tree
    for node in tre.traverse(strategy="preorder"):

        nstyle = ete3.NodeStyle()
        nstyle["vt_line_width"] = 1
        nstyle["hz_line_width"] = 1
        nstyle["vt_line_color"] = "grey"
        nstyle["hz_line_color"] = "grey"



        if info_colors:
            if node.info and node.info != "OTH PARENT":
                if "GR FIX" in node.info:
                    nstyle["fgcolor"] = "limegreen" if color_GR else "grey"
                elif "GR BKB" in node.info:
                    nstyle["fgcolor"] = "greenyellow" if color_GR else "grey"
                elif "OTH FIX" in node.info:
                    if fix_node:
                        if node.name == fix_node:
                            nstyle["fgcolor"] = "red" if color_OTH else "grey"
                            nstyle["shape"] = "sphere"
                        else:
                            nstyle["fgcolor"] = "grey"
                    else:
                        nstyle["fgcolor"] = "red" if color_OTH else "grey"
                        nstyle["shape"] = "sphere"
                elif "OTH BKB 1" in node.info:
                    if bkb_label:
                        if bkb_label in node.info:
                            nstyle["fgcolor"] = "violet" if color_OTH else "grey"
                        else:
                            nstyle["fgcolor"] = "grey"
                    else:
                         nstyle["fgcolor"] = "violet" if color_OTH else "grey"
                elif "OTH BKB 2" in node.info:
                    if bkb_label:
                        if bkb_label in node.info:
                            # nstyle["fgcolor"] = "mediumorchid" if color_OTH else "grey"
                            nstyle["fgcolor"] = "violet" if color_OTH else "grey"
                        else:
                            nstyle["fgcolor"] = "grey"
                    else:
                         nstyle["fgcolor"] = "violet" if color_OTH else "grey"
                elif "OTH BKB 3" in node.info:
                    if bkb_label:
                        if bkb_label in node.info:
                            # nstyle["fgcolor"] = "darkviolet" if color_OTH else "grey"
                            nstyle["fgcolor"] = "violet" if color_OTH else "grey"
                        else:
                            nstyle["fgcolor"] = "grey"
                    else:
                         nstyle["fgcolor"] = "violet" if color_OTH else "grey"
                elif "OTH EX BKB" in node.info:
                    if bkb_label:
                        if bkb_label in node.info:
                            # nstyle["fgcolor"] = "darkmagenta" if color_OTH else "grey"
                            nstyle["fgcolor"] = "violet" if color_OTH else "grey"
                        else:
                            nstyle["fgcolor"] = "grey"
                    else:
                         nstyle["fgcolor"] = "violet" if color_OTH else "grey"
                elif "NMP FIX" in node.info:
                    nstyle["fgcolor"] = "deepskyblue" if color_NMP else "grey"
                elif "NMP BKB" in node.info:
                    nstyle["fgcolor"] = "turquoise" if color_NMP else "grey"
                elif "NMP EX BKB" in node.info:
                    nstyle["fgcolor"] = "royalblue" if color_NMP else "grey"
            else:
                nstyle["fgcolor"] = "grey"
        else:
            if node.ph_tx == "TX":
                nstyle["fgcolor"] = "forestgreen"
            elif node.ph_tx == "PH":
                nstyle["fgcolor"] = "greenyellow"

        if nstyle["fgcolor"] == "grey":
            nstyle["size"] = 10
        else:
            nstyle["size"] = 20
        node.set_style(nstyle)

        if not (tx_levels[node.tx_level] == tx_levels["mrca"] and not name_mrcas):
            if tx_levels[node.tx_level] == tx_levels["species"]:
                node.add_face(ete3.TextFace(name_to_simple_name(node.name), fstyle="italic"), column=0, position="branch-right")
                node.add_face(ete3.TextFace(""), column=0, position="branch-bottom")
            else:
                node.add_face(ete3.TextFace(name_to_simple_name(node.name) + "      "), column=0, position="branch-top")
                node.add_face(ete3.TextFace("%s" % (node.tx_level)), column=0, position="branch-bottom")

    ts = ete3.TreeStyle()
    ts.margin_right=100
    ts.show_leaf_name = False
    ts.mode ="r"

    ts.scale = 20
    ts.show_scale = True

    tre.render(filename, tree_style=ts)


def plot_dates(input_tre, filename, show_ranks=False, show_dates=True, show_nodes=True, title=None):
    tre = input_tre.copy()

    # look for inconsistent dating, ie where a node has a later date than one of its ancestors
    # mrad = most recent ancestor date of the current node - if the current node is longer ago than this, that's bad
    def check_date_consistency(node, mrad):
        nstyle = ete3.NodeStyle()
        nstyle["vt_line_width"] = 1
        nstyle["hz_line_width"] = 1
        nstyle["vt_line_color"] = "gray"
        nstyle["hz_line_color"] = "gray"

        next_mrad = mrad
        # if we have a date, check it
        if node.date:
            if node.date > mrad:
                nstyle["fgcolor"] = "magenta"
                nstyle["size"] = 30
                node.set_style(nstyle)
            else:
                if node.imputed_date:
                    nstyle["fgcolor"] = "orange"
                    nstyle["size"] = 25
                    node.set_style(nstyle)
                else:
                    nstyle["fgcolor"] = "red"
                    nstyle["size"] = 20
                    node.set_style(nstyle)
                next_mrad = node.date
        else:
            if node.is_leaf():
                nstyle["fgcolor"] = "grey"
                nstyle["size"] = 12
                node.set_style(nstyle)
            else:
                nstyle["fgcolor"] = "blue"
                nstyle["size"] = 15
                node.set_style(nstyle)

        if not show_nodes:
            nstyle["size"] = 0
            # nstyle["vt_line_width"] = 2
            nstyle["hz_line_width"] = 2
            if node.imputed_date:
                if node.imputation_type == 5:
                    nstyle["hz_line_color"] = "green"
                    nstyle["vt_line_width"] = 1
                    nstyle["hz_line_width"] = 2
                elif node.imputation_type > 0:
                    nstyle["hz_line_color"] = "crimson"
                    nstyle["vt_line_width"] = 1
                    nstyle["hz_line_width"] = 2

                node.set_style(nstyle)

            else:
                if node.up and node.up.imputed_date:
                    if node.up.imputation_type == 5:
                        nstyle["hz_line_color"] = "green"
                        nstyle["vt_line_width"] = 1
                        nstyle["hz_line_width"] = 2
                    elif node.up.imputation_type > 0:
                        nstyle["hz_line_color"] = "crimson"
                        nstyle["vt_line_width"] = 1
                        nstyle["hz_line_width"] = 2

                    node.set_style(nstyle)

        for child in node.children:
            check_date_consistency(child, next_mrad)

    check_date_consistency(tre, 5000)

    # colors = ["#e6194b",
    #           "#f58231",
    #           "#3cb44b",
    #           "#4363d8",
    #           "#911eb4",
    #           "#f032e6",
    #           "#bfef45",
    #           "#42d4f4",
    #           "#800000",
    #           "#a9a9a9",
    #           "#808000",
    #           "#469990",
    #           "#000075",
    #           "#000000",
    #           "#9a6324",
    #           "#fabed4",
    #           "#aaffc3",
    #           "#dcbeff",
    #           "#ffd8b1",
    #           "#ffe119",
    #           "#ffd8b1",
    #           "#fffac8"]

    # sources_dict = {}
    # color_index = 0

    if show_dates:
        for node in tre.traverse(strategy='preorder'):
            if node.date is not None:
                # median / mean on top of branch
                #node.add_face(ete3.TextFace("Med %.1f / Mean %.1f" % (node.date[1], node.date[0])), column=0, position="branch-top")
                if node.imputed_date:
                    node.add_face(ete3.TextFace("Imp. %.2f" % (node.date)), column=0, position="branch-top")
                else:
                    if node.date == 0:
                        node.add_face(ete3.TextFace("tip"), column=0, position="branch-top")
                    else:
                        node.add_face(ete3.TextFace("%.2f " % (node.date)), column=0, position="branch-top")

                if show_ranks:
                    node.add_face(ete3.TextFace("%s, %s" % (node.name, node.tx_level), fgcolor="black"), column=0, position="branch-bottom")

                # list of dates on bottom of branch, coloured by the phylogeny they came from
                #for i in range(len(node.date[2])):
                #    date = node.date[2][i]
                    # source = node.date[3][i][1]
                    # if source not in sources_dict:
                    #     sources_dict[source] = colors[color_index]
                    #     color_index += 1
                    #     if color_index > 0 and color_index % len(colors) == 0:
                    #         color_index = 0
                    #node.add_face(ete3.TextFace("%.1f    " % date, fgcolor=sources_dict[source]), column=i, position="branch-bottom")
                #    node.add_face(ete3.TextFace("%.1f    " % date, fgcolor="black"), column=i, position="branch-bottom")
            else:
                node.add_face(ete3.TextFace(" ", fgcolor="black"), column=0, position="branch-top")
                if node.is_leaf():
                    node.add_face(ete3.TextFace("tip", fgcolor="black"), column=0, position="branch-bottom")
                else:
                    node.add_face(ete3.TextFace(" ", fgcolor="black"), column=0, position="branch-bottom")

    ts = ete3.TreeStyle()
    ts.show_leaf_name = False
    ts.mode ="r"
    ts.scale = 30

    if not show_dates:
        ts.min_leaf_separation = 10
    if not show_nodes:
        ts.complete_branch_lines_when_necessary = False

    if title:
        ts.title.add_face(ete3.faces.TextFace(title, fsize=16), 0)

    tre.render(filename, tree_style=ts)


def plot_dates_dq(input_tre, filename):
    tre = input_tre.copy()

    # look for inconsistent dating, ie where a node has a later date than one of its ancestors
    # mrad = most recent ancestor date of the current node - if the current node is longer ago than this, that's bad
    def check_date_consistency(node, mrad):
        nstyle = ete3.NodeStyle()
        nstyle["vt_line_width"] = 1
        nstyle["hz_line_width"] = 1
        nstyle["vt_line_color"] = "gray"
        nstyle["hz_line_color"] = "gray"

        next_mrad = mrad
        # if we have a date, check it
        if node.date:
            if np.median(node.date) > mrad:
                nstyle["fgcolor"] = "magenta"
                nstyle["size"] = 30
                node.set_style(nstyle)
            else:
                if node.imputed_date:
                    nstyle["fgcolor"] = "orange"
                    nstyle["size"] = 25
                    node.set_style(nstyle)
                else:
                    nstyle["fgcolor"] = "red"
                    nstyle["size"] = 20
                    node.set_style(nstyle)
                next_mrad = np.median(node.date)
        else:
            nstyle["fgcolor"] = "blue"
            nstyle["size"] = 15
            node.set_style(nstyle)

        for child in node.children:
            check_date_consistency(child, next_mrad)

    check_date_consistency(tre, 5000)

    colors = ["#e6194b",
              "#f58231",
              "#3cb44b",
              "#4363d8",
              "#911eb4",
              "#f032e6",
              "#bfef45",
              "#42d4f4",
              "#800000",
              "#a9a9a9",
              "#808000",
              "#469990",
              "#000075",
              "#000000",
              "#9a6324",
              "#fabed4",
              "#aaffc3",
              "#dcbeff",
              "#ffd8b1",
              "#ffe119",
              "#ffd8b1",
              "#fffac8"]

    sources_dict = {}
    color_index = 0
    for node in tre.traverse(strategy='preorder'):
        if node.date is not None:
            # median / mean on top of branch
            #node.add_face(ete3.TextFace("Med %.1f / Mean %.1f" % (node.date[1], node.date[0])), column=0, position="branch-top")

            if node.date == 0:
                node.date = [0]

            dt_str = "%.1f" % node.date[0]
            if node.date_sources:
                dt_str += " %s" % node.date_sources[0]
            for i in range(1,len(node.date)):
                dt_str += "\n%.1f" % node.date[i]
                if node.date_sources:
                    dt_str += " %s" % node.date_sources[i]

            node.add_face(ete3.TextFace("%s " % (dt_str)), column=0, position="branch-top")

            node.add_face(ete3.TextFace("%s\n%s" % (node.name, node.tx_level), fgcolor="black"), column=0, position="branch-bottom")

        else:
            #node.add_face(ete3.TextFace("", fgcolor="black"), column=0, position="branch-top")
            node.add_face(ete3.TextFace("%s\n%s" % (node.name, node.tx_level), fgcolor="black"), column=0, position="branch-bottom")

    ts = ete3.TreeStyle()
    ts.show_leaf_name = False
    ts.mode ="r"
    ts.scale = 100

    tre.render(filename, tree_style=ts)


def digit_to_hex(dig):
    """Convert an integer between 0 and 15 to a string representing a hex value between 0 and F."""
    if dig < 10:
        return str(dig)
    elif dig == 10:
        return "A"
    elif dig == 11:
        return "B"
    elif dig == 12:
        return "C"
    elif dig == 13:
        return "D"
    elif dig == 14:
        return "E"
    elif dig == 15:
        return "F"
    else:
        print(dig)
        raise ValueError


def pct_to_color_hex_str(pct):
    """Convert a number between 0 and 1 to a string representing a hex value between 00 and FF."""
    p = pct * 255
    first = int(p/16)
    second = int(p % 16)
    return digit_to_hex(first) + digit_to_hex(second)



def plot_dates_pct(input_tre, filename):
    tre = input_tre.copy()

    # look for inconsistent dating, ie where a node has a later date than one of its ancestors
    # mrad = most recent ancestor date of the current node - if the current node is longer ago than this, that's bad
    def check_date_consistency(node, mrad):
        nstyle = ete3.NodeStyle()
        nstyle["vt_line_width"] = 1
        nstyle["hz_line_width"] = 1
        nstyle["vt_line_color"] = "gray"
        nstyle["hz_line_color"] = "gray"

        next_mrad = mrad
        # # if we have a date, check it
        # if node.date:
        #     if np.median(node.date) > mrad:
        #         nstyle["fgcolor"] = "magenta"
        #         nstyle["size"] = 50
        #         node.set_style(nstyle)
        #     else:
        #         if node.imputed_date:
        #             nstyle["fgcolor"] = "orange"
        #             nstyle["size"] = 50
        #             node.set_style(nstyle)
        #         else:
        #             nstyle["fgcolor"] = "red"
        #             nstyle["size"] = 50
        #             node.set_style(nstyle)
        #         next_mrad = np.median(node.date)
        # else:
        #     nstyle["fgcolor"] = "blue"
        #     nstyle["size"] = 50
        #     node.set_style(nstyle)

        pct_child_dates = node.num_dates/node.child_tree_size
        pct_child_dates = min(0.8, 10*pct_child_dates)
        if pct_child_dates < 0.0000001:
            nstyle["fgcolor"] = "pink"
        else:
            nstyle["fgcolor"] = "#" + pct_to_color_hex_str(1-max(pct_child_dates,0.05)) + pct_to_color_hex_str(1-max(pct_child_dates,0.05)) + pct_to_color_hex_str(1-max(pct_child_dates,0.05))
        nstyle["size"] = 60 * (1-pct_child_dates)
        node.set_style(nstyle)

        for child in node.children:
            check_date_consistency(child, next_mrad)

    check_date_consistency(tre, 5000)

    for node in tre.traverse(strategy='preorder'):
        if node.date is not None and not node.imputed_date:
            # median / mean on top of branch
            #node.add_face(ete3.TextFace("Med %.1f / Mean %.1f" % (node.date[1], node.date[0])), column=0, position="branch-top")

            if node.date_dist == 0:
                node.date_dist = [0]

            dt_str = "%.1f" % node.date_dist[0]
            if node.date_sources:
                dt_str += " %s" % node.date_sources[0]
            for i in range(1,len(node.date_dist)):
                dt_str += "\n%.1f" % node.date_dist[i]
                if node.date_sources:
                    dt_str += " %s" % node.date_sources[i]

            node.add_face(ete3.TextFace("%s " % (dt_str)), column=0, position="branch-top")

            node.add_face(ete3.TextFace("%s\n%s\n%.2f%% / %d" % (node.name, node.tx_level, 100*node.num_dates/node.child_tree_size, node.child_tree_size), fgcolor="black"), column=0, position="branch-bottom")

        else:
            #node.add_face(ete3.TextFace("", fgcolor="black"), column=0, position="branch-top")
            node.add_face(ete3.TextFace("%s\n%s\n%.2f%% / %d" % (node.name, node.tx_level, 100*node.num_dates/node.child_tree_size, node.child_tree_size), fgcolor="black"), column=0, position="branch-bottom")

    ts = ete3.TreeStyle()
    ts.show_leaf_name = False
    ts.mode ="r"
    ts.scale = .1
    ts.branch_vertical_margin = 5

    tre.render(filename, tree_style=ts)


def plot_dates_compact(input_tre, filename):
    tre = input_tre.copy()

    # look for inconsistent dating, ie where a node has a later date than one of its ancestors
    # mrad = most recent ancestor date of the current node - if the current node is longer ago than this, that's bad
    def check_date_consistency(node, mrad):
        nstyle = ete3.NodeStyle()
        nstyle["vt_line_width"] = 1
        nstyle["hz_line_width"] = 1
        nstyle["vt_line_color"] = "gray"
        nstyle["hz_line_color"] = "gray"

        next_mrad = mrad
        # # if we have a date, check it
        # if node.date:
        #     if np.median(node.date) > mrad:
        #         nstyle["fgcolor"] = "magenta"
        #         nstyle["size"] = 50
        #         node.set_style(nstyle)
        #     else:
        #         if node.imputed_date:
        #             nstyle["fgcolor"] = "orange"
        #             nstyle["size"] = 50
        #             node.set_style(nstyle)
        #         else:
        #             nstyle["fgcolor"] = "red"
        #             nstyle["size"] = 50
        #             node.set_style(nstyle)
        #         next_mrad = np.median(node.date)
        # else:
        #     nstyle["fgcolor"] = "blue"
        #     nstyle["size"] = 50
        #     node.set_style(nstyle)

        pct_child_dates = node.num_dates/node.child_tree_size
        pct_child_dates = min(1, 10*pct_child_dates)
        if pct_child_dates < 0.0000001:
            nstyle["fgcolor"] = "pink"
        else:
            nstyle["fgcolor"] = "#" + pct_to_color_hex_str(1-max(pct_child_dates,0.05)) + pct_to_color_hex_str(1-max(pct_child_dates,0.05)) + pct_to_color_hex_str(1-max(pct_child_dates,0.05))
        #nstyle["size"] = 60 * (1-pct_child_dates)
        nstyle["size"] = 5
        node.set_style(nstyle)

        for child in node.children:
            check_date_consistency(child, next_mrad)

    check_date_consistency(tre, 5000)

    for node in tre.traverse(strategy='preorder'):
        if node.date is not None:
            # median / mean on top of branch
            #node.add_face(ete3.TextFace("Med %.1f / Mean %.1f" % (node.date[1], node.date[0])), column=0, position="branch-top")

            if node.date == 0:
                node.date = [0]

            dt_str = "%.1f" % node.date[0]
            if node.date_sources:
                dt_str += " %s" % node.date_sources[0]
            for i in range(1,len(node.date)):
                dt_str += "\n%.1f" % node.date[i]
                if node.date_sources:
                    dt_str += " %s" % node.date_sources[i]

            node.add_face(ete3.TextFace("%s " % (dt_str)), column=0, position="branch-top")

            node.add_face(ete3.TextFace("%s\n%s\n%.2f%% / %d" % (node.name, node.tx_level, 100*node.num_dates/node.child_tree_size, node.child_tree_size), fgcolor="black"), column=0, position="branch-bottom")

        else:
            #node.add_face(ete3.TextFace("", fgcolor="black"), column=0, position="branch-top")
            node.add_face(ete3.TextFace("%s\n%s\n%.2f%% / %d" % (node.name, node.tx_level, 100*node.num_dates/node.child_tree_size, node.child_tree_size), fgcolor="black"), column=0, position="branch-bottom")

    ts = ete3.TreeStyle()
    ts.show_leaf_name = False
    ts.mode ="r"
    ts.scale = 200
    ts.branch_vertical_margin = 0
    ts.show_scale = True

    tre.render(filename, tree_style=ts)


def plot_dates_circle(input_tre, filename, show_ranks=False, show_pct=False, not_circle=False):
    tre = input_tre.copy()

    # look for inconsistent dating, ie where a node has a later date than one of its ancestors
    # mrad = most recent ancestor date of the current node - if the current node is longer ago than this, that's bad
    def check_date_consistency(node, mrad):
        nstyle = ete3.NodeStyle()
        nstyle["vt_line_width"] = 1
        nstyle["hz_line_width"] = 1
        nstyle["vt_line_color"] = "gray"
        nstyle["hz_line_color"] = "gray"

        next_mrad = mrad
        # if we have a date, check it
        if node.date:
            if node.date > mrad:
                nstyle["fgcolor"] = "magenta"
                nstyle["size"] = 30
                node.set_style(nstyle)
            else:
                if node.imputed_date:
                    nstyle["fgcolor"] = "orange"
                    nstyle["size"] = 16
                    node.set_style(nstyle)
                else:
                    nstyle["fgcolor"] = "red"
                    nstyle["size"] = 16
                    node.set_style(nstyle)
                next_mrad = node.date
        else:
            nstyle["fgcolor"] = "blue"
            nstyle["size"] = 15
            node.set_style(nstyle)

        for child in node.children:
            check_date_consistency(child, next_mrad)

    check_date_consistency(tre, 5000)

    colors = ["#e6194b",
              "#f58231",
              "#3cb44b",
              "#4363d8",
              "#911eb4",
              "#f032e6",
              "#bfef45",
              "#42d4f4",
              "#800000",
              "#a9a9a9",
              "#808000",
              "#469990",
              "#000075",
              "#000000",
              "#9a6324",
              "#fabed4",
              "#aaffc3",
              "#dcbeff",
              "#ffd8b1",
              "#ffe119",
              "#ffd8b1",
              "#fffac8"]

    sources_dict = {}
    color_index = 0
    for node in tre.traverse(strategy='preorder'):
        if node.is_leaf():
            if node.tx_level != "mrca":
                node.add_face(ete3.TextFace("%s\n%s" % (name_to_simple_name(node.name), node.tx_level)), column=0, position="branch-right")
            if node.date is not None:
                node.add_face(ete3.TextFace("Date: %.1f mya\n%.2f%% of\n%d nodes dated\nPD: %.0f bn yrs" % (node.date, 100*node.num_dates/node.child_tree_size, node.child_tree_size, node.pd/1000), fgcolor="black"), column=0, position="branch-right")
            else:
                node.add_face(ete3.TextFace("% "), column=0, position="branch-right")
        else:
            if node.tx_level != "mrca":
                node.add_face(ete3.TextFace("%s\n%s" % (name_to_simple_name(node.name), node.tx_level)), column=0, position="branch-top")
            if node.date is not None:
                # node.add_face(ete3.TextFace("%.1f " % (node.date)), column=0, position="branch-bottom")
                node.add_face(ete3.TextFace("Date: %.1f mya\n%.2f%% of\n%d nodes dated\nPD: %.0f bn yrs" % (node.date, 100*node.num_dates/node.child_tree_size, node.child_tree_size, node.pd/1000), fgcolor="black"), column=0, position="branch-bottom")
            else:
                node.add_face(ete3.TextFace("% "), column=0, position="branch-bottom")


    ts = ete3.TreeStyle()
    ts.show_leaf_name = False
    if not_circle:
        ts.mode = "r"
    else:
        ts.mode = "c"
    ts.scale = .1
    ts.branch_vertical_margin = 8
    ts.show_scale = False
    ts.margin_right = 100

    tre.render(filename, tree_style=ts)


def plot_dates_figure(input_tre, filename, show_ranks=False, show_pct=False, log_scale=False):
    tre = input_tre.copy()

    tre.dist = 0.01

    if log_scale:
        for node in tre.traverse():
            node.dist = np.log(node.dist)

    # look for inconsistent dating, ie where a node has a later date than one of its ancestors
    # mrad = most recent ancestor date of the current node - if the current node is longer ago than this, that's bad
    def check_date_consistency(node, mrad):
        nstyle = ete3.NodeStyle()
        nstyle["vt_line_width"] = 3
        nstyle["hz_line_width"] = 3
        nstyle["vt_line_color"] = "black"
        nstyle["hz_line_color"] = "black"

        if "mrcaott2ott3973" in node.name:
            nstyle["bgcolor"] = "palegreen"
        elif "Metazoa" in node.name:
            nstyle["bgcolor"] = "lightblue"
        elif "Fungi" in node.name:
            nstyle["bgcolor"] = "#CFCFCF"

        next_mrad = mrad
        # if we have a date, check it
        if node.date:
            if node.date > mrad:
                nstyle["fgcolor"] = "magenta"
                nstyle["size"] = 30
                node.set_style(nstyle)
            else:
                if node.imputed_date:
                    nstyle["fgcolor"] = "orange"
                    nstyle["size"] = 16
                    node.set_style(nstyle)
                else:
                    nstyle["fgcolor"] = "red"
                    nstyle["size"] = 16
                    node.set_style(nstyle)
                next_mrad = node.date
        else:
            nstyle["fgcolor"] = "blue"
            nstyle["size"] = 15
            node.set_style(nstyle)

        for child in node.children:
            check_date_consistency(child, next_mrad)

    check_date_consistency(tre, 5000)

    colors = ["#e6194b",
              "#f58231",
              "#3cb44b",
              "#4363d8",
              "#911eb4",
              "#f032e6",
              "#bfef45",
              "#42d4f4",
              "#800000",
              "#a9a9a9",
              "#808000",
              "#469990",
              "#000075",
              "#000000",
              "#9a6324",
              "#fabed4",
              "#aaffc3",
              "#dcbeff",
              "#ffd8b1",
              "#ffe119",
              "#ffd8b1",
              "#fffac8"]

    sources_dict = {}
    color_index = 0
    for node in tre.traverse(strategy='preorder'):
        if node.is_leaf():
            if node.tx_level != "mrca":
                node.add_face(ete3.TextFace("%s" % (name_to_simple_name(node.name)),fsize=28), column=0, position="branch-right")
            if node.date is not None:
                node.add_face(ete3.TextFace("%.1f mya\n%.2f%% of\n%d nodes dated\nPD: %.0f bn yrs" % (
                    node.date, 100*node.num_dates/node.child_tree_size, node.child_tree_size, node.pd/1000), fsize=28),
                    column=0, position="branch-right")
            else:
                node.add_face(ete3.TextFace("% "), column=0, position="branch-right")
        else:
            if node.tx_level != "mrca":
                node.add_face(ete3.TextFace("%s" % (name_to_simple_name(node.name))), column=0, position="branch-top")
            if node.date is not None:
                # node.add_face(ete3.TextFace("%.1f " % (node.date)), column=0, position="branch-bottom")
                node.add_face(ete3.TextFace("%.1f mya" % (node.date), fgcolor="black"), column=0, position="branch-bottom")
            else:
                node.add_face(ete3.TextFace("% "), column=0, position="branch-bottom")


    ts = ete3.TreeStyle()
    ts.show_leaf_name = False
    ts.mode = "c"
    ts.scale = 30
    #ts.branch_vertical_margin = 8
    ts.show_scale = False
    ts.margin_right = 100
    ts.root_opening_factor = 0

    tre.render(filename, tree_style=ts)


def plot_big_tree(input_tre, filename, scale=0.5):
    import ete4
    tre = input_tre.copy()

    # orig_nodes = []
    # for node in tre.traverse(strategy="preorder"):
    #     orig_nodes.append(node)

    # for node in orig_nodes:
    #     # 1. Create dummy intermediate nodes
    #     # We give them a tiny distance to create the 'step' look
    #     if len(node.children) > 0:
    #         print(node.name, len(node.children))
    #         a = node.children[1].detach()
    #         b = node.children[0].detach()

    #         d1 = node.add_child(name="dummy_red", dist=0.1)
    #         d2 = node.add_child(name="dummy_blue", dist=0.1)

    #         d1.add_child(a)
    #         d2.add_child(b)
    #         d1.add_prop("pd", a.props["pd"])
    #         d2.add_prop("pd", b.props["pd"])

    import matplotlib
    cmap = matplotlib.colormaps['viridis']

    all_pds = []
    for node in tre.traverse(strategy="preorder"):
        all_pds.append(node.props["pd"])

    max_pd = np.log(np.max(all_pds))
    min_pd = np.log(np.min(all_pds))

    for node in tre.traverse(strategy="preorder"):
        nstyle = ete4.treeview.NodeStyle()

        pd_conversion = (np.log(node.props["pd"]) - min_pd) / (max_pd - min_pd)
        rgba = cmap(1-pd_conversion)
        color_str = "#"
        for i in range(3):
            color_str += pct_to_color_hex_str(rgba[i])

        print(node.name, node.props["pd"], pd_conversion, color_str)
        nstyle["vt_line_color"] = color_str
        nstyle["vt_line_type"] = 0
        nstyle["vt_line_width"] = 3
        nstyle["hz_line_color"] = color_str
        nstyle["hz_line_type"] = 0
        nstyle["hz_line_width"] = 3

        if node is tre:
            # colour the root node itself
            nstyle["size"] = 15
            nstyle["shape"] = "square"
            rgba = cmap(0)
            color_str = "#"
            for i in range(3):
                color_str += pct_to_color_hex_str(rgba[i])
            nstyle["fgcolor"] = color_str
        else:
            nstyle["size"] = 0

        if node.is_leaf:
            node.add_face(ete4.treeview.TextFace("%s" % (name_to_simple_name(node.name)), fsize=8), column=0, position="branch-right")

        node.set_style(nstyle)

    ts = ete4.treeview.TreeStyle()
    ts.mode = "r"
    ts.show_leaf_name = False

    ts.show_scale = False
    ts.margin_top = 0
    ts.margin_bottom = 0
    ts.margin_right = 30
    ts.margin_left = 0
    ts.root_opening_factor = 0
    ts.draw_guiding_lines = True
    ts.branch_vertical_margin = 5
    ts.draw_guiding_lines = True

    ts.scale = scale

    tre.render(filename, tree_style=ts)


def plot_dates_figure_outline(input_tre, filename, show_ranks=False, show_pct=False, log_scale_dates=False, log_scale_branches=False, simple_label=False, dpi=None):
    import ete4
    tre = input_tre.copy()

    import matplotlib

    cmap = matplotlib.cm.get_cmap('coolwarm')

    if log_scale_dates:
        for node in tre.traverse():
            node.props["date"] = np.log(node.props["date"])
        for node in tre.traverse():
            if node == tre:
                node.dist = 0.1
                continue
            node.dist = node.up.props["date"] - node.props["date"]
    elif log_scale_branches:
        for node in tre.traverse():
            # if node is not tre:
            node.dist = np.log(node.dist)

        edge_date = tre.props["date"] - tre.get_farthest_leaf()[1]
        farthest_dist = tre.get_farthest_leaf()[1]

        leaves = tre.leaves()

        for node in leaves:
            new_node = tree_fixing.create_node(node.name)
            node.add_child(new_node)
            new_node.props["date"] = node.props["date"]
            face_radius = 3 + node.props["num_leaves"]**(1/3)/4
            # face_radius = 0.7 * (5 + 18*np.sqrt(node.props["num_leaves"]) / 900) / 3.14
            # new_node.dist = farthest_dist - node.get_distance(node, tre) + farthest_dist/10 - face_radius
            new_node.dist = farthest_dist + farthest_dist/15 - node.get_distance(node, tre) - face_radius/4.5
            new_node.add_prop("num_leaves", node.props["num_leaves"])
            new_node.add_prop("child_tree_size", node.props["child_tree_size"])
            new_node.add_prop("num_dates", node.props["num_dates"])
            new_node.add_prop("pd", node.props["pd"])
            new_node.add_prop("imputed_date", node.props["imputed_date"])

    n_leaves = len(tre)

    # Map leaves to angles
    leaf_angles = {}
    for i, leaf in enumerate(tre.leaves()):
        angle = (360.0 * i) / n_leaves
        leaf_angles[leaf.name] = angle

    plants = False
    animals = False
    fungi = False

    cl = "black"
    for node in tre.traverse(strategy='preorder'):
        # if "Opisthokonta" in node.name:
        #     cl = "black"
        # if "mrcaott2ott3973" in node.name:
        #     cl = "#009E73"
        # elif "Metazoa" in node.name:
        #     cl = "#0072B2"
        # # elif "Fungi" in node.name:
        # #     cl = "red"
        # elif len(node.children) > 0 and "Microsporidia" in node.children[0].name:
        #     cl = "#E69F00"
        cl = "#333333"

        nstyle = ete4.treeview.NodeStyle()
        nstyle["vt_line_width"] = 3
        nstyle["hz_line_width"] = 3
        nstyle["vt_line_color"] = cl
        nstyle["hz_line_color"] = cl
        if log_scale_dates and "Bacillariophyta" in node.name:
            nstyle["hz_line_color"] = "grey"
            nstyle["hz_line_type"] = 1

        # if "Holozoa" in node.name:
        #     nstyle["bgcolor"] = "lightblue"
        # elif "mrcaott2ott3973" in node.name:
        #     nstyle["bgcolor"] = "lightgreen"
        # elif "Nucletmycea" in node.name:
        #     nstyle["bgcolor"] = "lightpink"

        if log_scale_dates and node.is_leaf:
            nstyle["vt_line_color"] = "white"
            nstyle["hz_line_color"] = "white"

        if log_scale_branches and (node.props["tx_level"] == "phylum" or node.props["tx_level"] == "infrakingdom"):
            # node_scaling = np.sqrt(node.props["num_leaves"]) / 1086
            # nstyle["size"] = 3 + 26*node_scaling
            node_size = node.props["num_leaves"]**(1/3) / 4
            nstyle["size"] = node_size

            if node.props["num_dates"] == 0:
                color_scale = 0
            else:
                # color_scale = 0.85 - 0.8*(np.log(1+100*node.num_dates/node.child_tree_size) / 3.22)
                color_scale = .15 + .85*(((100*node.props["num_dates"]/node.props["child_tree_size"])**0.25) / 2.2)

            rgba = cmap(color_scale)
            color_str = "#"
            for i in range(3):
                color_str += pct_to_color_hex_str(rgba[i])
            # nstyle["fgcolor"] = "#" + color_str + color_str + color_str
            # print(color_str)
            # nstyle["fgcolor"] = color_str
            nstyle["fgcolor"] = cl
            nstyle["shape"] = "sphere"

        else:
            nstyle["size"] = 0

        if log_scale_branches and node.is_leaf:
            nstyle["hz_line_color"] = "#CFCFCF"
            nstyle["hz_line_type"] = 1
            nstyle["hz_line_width"] = 1
            nstyle["shape"] = "square"

            nstyle["size"] = 19

            # pd_scale = ((node.pd/1000)**(1/6)-0.2) / 5.6
            # pd_scale = ((node.pd/1000)**(1/6)-0.745) / (5.77 - 0.745)
            # pd_scale = (np.log(node.pd/1000)+2.065) / (10.498+2.065)
            pd_scale = (np.log(node.props["pd"]/1000)-0.8788) / (9.44405685-0.8788)
            print('PD:', node.name, pd_scale, node.props["pd"], node.props["ph_tx"], node.props["imputed_date"])
            # print(pd_scale)
            rgba = cmap(pd_scale)
            color_str = "#"
            for i in range(3):
                color_str += pct_to_color_hex_str(rgba[i])
            nstyle["fgcolor"] = color_str

        node.set_style(nstyle)


        if node.is_leaf:
            if log_scale_dates:
                node.add_face(ete4.TextFace("%s" % (name_to_simple_name(node.name)), fsize=8),
                        column=0, position="branch-right")
            elif log_scale_branches:
                if node.props["num_leaves"] < 1000:
                    num_spp_str = "%.0f" % (node.props["num_leaves"])
                elif node.props["num_leaves"] < 10000:
                    num_spp_str = "%.2fk" % (node.props["num_leaves"]/1000)
                elif node.props["num_leaves"] < 100000:
                    num_spp_str = "%.1fk" % (node.props["num_leaves"]/1000)
                elif node.props["num_leaves"] < 1000000:
                    num_spp_str = "%.0fk" % (node.props["num_leaves"]/1000)
                else:
                    num_spp_str = "%.2fM" % (node.props["num_leaves"]/1000000)

                if node.props["num_dates"] < 0.000001:
                    pct_dated_str = "no"
                elif 100*node.props["num_dates"]/node.props["child_tree_size"] < 0.1:
                    pct_dated_str = "%.3f%%" % (100*node.props["num_dates"]/node.props["child_tree_size"])
                elif 100*node.props["num_dates"]/node.props["child_tree_size"] < 1:
                    pct_dated_str = "%.2f%%" % (100*node.props["num_dates"]/node.props["child_tree_size"])
                elif 100*node.props["num_dates"]/node.props["child_tree_size"] < 10:
                    pct_dated_str = "%.1f%%" % (100*node.props["num_dates"]/node.props["child_tree_size"])
                else:
                    pct_dated_str = "%.0f%%" % (100*node.props["num_dates"]/node.props["child_tree_size"])

                if node.props["pd"] < 1000:
                    pd_str = "%.0f Myr" % (node.props["pd"])
                elif node.props["pd"] < 10000:
                    pd_str = "%.2f Byr" % (node.props["pd"]/1000)
                elif node.props["pd"] < 100000:
                    pd_str = "%.1f Byr" % (node.props["pd"]/1000)
                elif node.props["pd"] < 1000000:
                    pd_str = "%.0f Byr" % (node.props["pd"]/1000)
                elif node.props["pd"] < 10000000:
                    pd_str = "%.2f Tyr" % (node.props["pd"]/1000000)
                elif node.props["pd"] < 100000000:
                    pd_str = "%.1f Tyr" % (node.props["pd"]/1000000)
                else:
                    pd_str = "%.0f Tyr" % (node.props["pd"]/1000000)

                # if simple_label:
                #     node.add_face(ete4.treeview.TextFace("%s" % (name_to_simple_name(node.name)), fsize=12),
                #         column=0, position="branch-right")
                #     node.add_face(ete4.treeview.TextFace("%.0f Mya, %s spp." % (node.props["date"], num_spp_str), fsize=10),
                #         column=0, position="branch-right")
                # else:
                #     node.add_face(ete4.treeview.TextFace("%s, %.0f Mya, %s\n%s spp., %s dated nodes" % (name_to_simple_name(node.name), node.props["date"], pd_str, num_spp_str, pct_dated_str), fsize=12),
                #         column=0, position="branch-right")


                if simple_label:
                    angle = leaf_angles.get(node.name, 0)
                    print(node.name, angle, node_size)
                    if 90 < angle < 271:
                        date_label_str = "%.0f Mya, %s spp. " % (node.props["date"], num_spp_str)
                        face1 = ete4.treeview.TextFace(date_label_str, fsize=11, ftype="Arial")

                        name_label = name_to_simple_name(node.name)
                        space = ""

                        for i in range(max(0,13-len(name_label))):
                            space += " "
                        for i in range(max(0,12-len(name_label))):
                            space += " "

                        if name_label == "Nematoda" or name_label == "Gnathostomulida" or name_label == "Rhombozoa" or name_label == "Platyhelminthes":
                            name_label += " "
                        if name_label == "Nematomorpha":
                            name_label += " "
                        if name_label == "Onychophora" or name_label == "Bryozoa":
                            name_label += " "
                        name_label += " "

                        face2 = ete4.treeview.TextFace("%s%s" % (space, name_label), fsize=14, ftype="Arial")
                        # face1.rotation = 181
                        # face2.rotation = 181

                        spacer = ete4.treeview.RectFace(width=32, height=1, fgcolor=None, bgcolor=None)
                        node.add_face(spacer, column=0, position="branch-right")

                    else:
                        face1 = ete4.treeview.TextFace(" %s" % (name_to_simple_name(node.name)), fsize=14, ftype="Arial")
                        face2 = ete4.treeview.TextFace(" %.0f Mya, %s spp." % (node.props["date"], num_spp_str), fsize=11, ftype="Arial")

                        spacer = ete4.treeview.RectFace(width=2, height=1, fgcolor=None, bgcolor=None)
                        node.add_face(spacer, column=0, position="branch-right")


                    node.add_face(face1,
                        column=1, position="branch-right")
                    node.add_face(face2,
                        column=1, position="branch-right")


                else:
                    node.add_face(ete4.treeview.TextFace("%s, %.0f Mya, %s\n%s spp., %s dated nodes" % (name_to_simple_name(node.name), node.props["date"], pd_str, num_spp_str, pct_dated_str), fsize=12),
                        column=0, position="branch-right")

                # print("%s\t%.0f\t%s\t%s\t%s\t%d" % (name_to_simple_name(node.name), node.props["date"], pd_str, num_spp_str, pct_dated_str, node.num_leaves))

        # else:
        #     if node.tx_level != "mrca":
        #         node.add_face(ete4.treeview.TextFace("%s" % (name_to_simple_name(node.name))), column=0, position="branch-top")
        #     if node.props["date"] is not None:
        #         # node.add_face(ete4.treeview.TextFace("%.1f " % (node.props["date"])), column=0, position="branch-bottom")
        #         node.add_face(ete4.treeview.TextFace("%.1f mya" % (node.props["date"]), fgcolor="black"), column=0, position="branch-bottom")
        #     else:
        #         node.add_face(ete4.treeview.TextFace("% "), column=0, position="branch-bottom")

    # def my_layout(node):
    #     if node.is_leaf:
    #         # Create a label face
    #         face = ete4.treeview.TextFace(node.name, fsize=10)

    #         # Calculate the angle of the node (used only in circular mode)
    #         # ETE will assign a node._angle attribute when using circular mode
    #         angle = getattr(face, "_angle", 0)
    #         print(node._angle)


    ts = ete4.treeview.TreeStyle()
    ts.mode = "c"
    # ts.layout_fn = my_layout
    ts.show_leaf_name = False

    #ts.branch_vertical_margin = 8
    ts.show_scale = False
    ts.margin_top = 0
    ts.margin_bottom = 0
    ts.margin_right = 0
    ts.margin_left = 0
    ts.root_opening_factor = 0
    # ts.complete_branch_lines_when_necessary = False
    ts.draw_guiding_lines = True
    ts.children_faces_on_top = False
    if log_scale_dates:
        ts.scale = 100
    elif log_scale_branches:
        ts.scale = 4.5
    ts.arc_start = 0 #235
    ts.arc_span = 359.99
    ts.min_leaf_separation = 0
    ts.allow_face_overlap = True
    ts.children_faces_on_top = False

    tre.render(filename, tree_style=ts)
    # import PIL
    # tre.render(filename + ".tif", tree_style=ts)


def plot_dates_algo(input_tre, filename, show_paths=True, show_only_undated_paths=True, mu=False, pinkblue=True, vs=3):
    tre = input_tre.copy()
    # tre.convert_to_ultrametric()

    # look for inconsistent dating, ie where a node has a later date than one of its ancestors
    # mrad = most recent ancestor date of the current node - if the current node is longer ago than this, that's bad
    def check_date_consistency(node, mrad):
        nstyle = ete3.NodeStyle()
        nstyle["vt_line_width"] = 1
        nstyle["hz_line_width"] = 1
        nstyle["vt_line_color"] = "gray"
        nstyle["hz_line_color"] = "gray"

        next_mrad = mrad
        # if we have a date, check it
        if node.is_leaf():
            nstyle["fgcolor"] = "#CFCFCF"
            nstyle["size"] = 10
            node.set_style(nstyle)
        else:
            if node.date and node.date > 0:
                if node.date > mrad:
                    nstyle["fgcolor"] = "magenta"
                    nstyle["size"] = 25
                    node.set_style(nstyle)
                else:
                    if node.imputed_date:
                        nstyle["size"] = 17
                        nstyle["shape"] = "square"
                    else:
                        nstyle["size"] = 18
                    if pinkblue:
                        nstyle["fgcolor"] = "#F1AEE8"
                    else:
                        nstyle["fgcolor"] = "#40B0A6"
                    node.set_style(nstyle)
                    next_mrad = node.date
            else:
                if pinkblue:
                    nstyle["fgcolor"] = "#768AE0"
                else:
                    nstyle["fgcolor"] = "#E1BE6A"
                nstyle["size"] = 17
                nstyle["shape"] = "square"
                node.set_style(nstyle)

        for child in node.children:
            check_date_consistency(child, next_mrad)

    check_date_consistency(tre, 5000)

    sources_dict = {}
    color_index = 0
    for node in tre.traverse(strategy='preorder'):
        if node.date is not None and node.date > 0:
            # median / mean on top of branch
            #node.add_face(ete3.TextFace("Med %.1f / Mean %.1f" % (node.date[1], node.date[0])), column=0, position="branch-top")
            if show_paths:
                node.add_face(ete3.TextFace("  %.2f " % (node.date), fgcolor="firebrick"), column=0, position="branch-top")
            else:
                node.add_face(ete3.TextFace("%.2f " % (node.date), fgcolor="firebrick"), column=0, position="branch-top")
            # if node.imputed_date and mu:
                # node.add_face(ete3.TextFace("%.2f  " % (node.date_long), fgcolor="blue"), column=0, position="branch-top")
                # node.add_face(ete3.TextFace("%.2f  " % (node.date_short), fgcolor="red"), column=0, position="branch-top")
                # node.add_face(ete3.TextFace("%.2f  " % (node.date_above_long), fgcolor="black"), column=0, position="branch-top")
                # node.add_face(ete3.TextFace("%.2f  " % (node.date_above_short), fgcolor="black"), column=0, position="branch-top")
                # node.add_face(ete3.TextFace(node.mu_spacing_long), column=0, position="branch-bottom")
                # node.add_face(ete3.TextFace(node.mu_spacing_short), column=0, position="branch-bottom")

        if node.is_leaf():
            node.add_face(ete3.TextFace(" "+node.name, fstyle="italic"), column=0, position="branch-right")

        if show_paths:
            if show_only_undated_paths:
                if node.date is None:
                    node.add_face(ete3.TextFace(oldest_path_to_str(node.oldest_path_long) + " "), column=0, position="branch-bottom")
                    # node.add_face(ete3.TextFace(oldest_path_to_str(node.oldest_path_short) + "  "), column=0, position="branch-bottom")
            else:
                node.add_face(ete3.TextFace(oldest_path_to_str(node.oldest_path_long) + " ", fsize=8), column=0, position="branch-bottom")
                # node.add_face(ete3.TextFace(oldest_path_to_str(node.oldest_path_short) + "  "), column=0, position="branch-bottom")

    ts = ete3.TreeStyle()
    ts.show_leaf_name = False
    ts.mode ="r"
    ts.scale = 100
    ts.branch_vertical_margin = vs
    ts.show_scale = False
    ts.margin_right = 10
    ts.complete_branch_lines_when_necessary = False

    tre.render(filename, tree_style=ts, units="in", w=6, dpi=300)


def plot_dates_algo_dq(input_tre, filename, show_paths=True, show_only_undated_paths=True):
    tre = input_tre.copy()

    # look for inconsistent dating, ie where a node has a later date than one of its ancestors
    # mrad = most recent ancestor date of the current node - if the current node is longer ago than this, that's bad
    def check_date_consistency(node, mrad):
        nstyle = ete3.NodeStyle()
        nstyle["vt_line_width"] = 1
        nstyle["hz_line_width"] = 1
        nstyle["vt_line_color"] = "gray"
        nstyle["hz_line_color"] = "gray"

        next_mrad = mrad
        # if we have a date, check it
        if node.is_leaf():
            nstyle["fgcolor"] = "#CFCFCF"
            nstyle["size"] = 10
            node.set_style(nstyle)
        else:
            if node.date and node.date > 0:
                if node.date > mrad:
                    nstyle["fgcolor"] = "magenta"
                    nstyle["size"] = 25
                    node.set_style(nstyle)
                else:
                    if node.imputed_date:
                        nstyle["size"] = 18
                        nstyle["shape"] = "square"
                    else:
                        nstyle["size"] = 18
                    nstyle["fgcolor"] = "lightpink"
                    node.set_style(nstyle)
                    next_mrad = node.date
            else:
                nstyle["fgcolor"] = "#80FFFF"
                nstyle["size"] = 18
                nstyle["shape"] = "square"
                node.set_style(nstyle)

        for child in node.children:
            check_date_consistency(child, next_mrad)

    check_date_consistency(tre, 5000)

    sources_dict = {}
    color_index = 0
    for node in tre.traverse(strategy='preorder'):
        if node.date is not None and node.date > 0 and not node.imputed_date:
            # median / mean on top of branch
            #node.add_face(ete3.TextFace("Med %.1f / Mean %.1f" % (node.date[1], node.date[0])), column=0, position="branch-top")
            node.add_face(ete3.TextFace("%.2f" % (node.date), fgcolor="red"), column=0, position="branch-top")
            if node.up:
                node.add_face(ete3.TextFace("%.2f" % (node.computed_date), fgcolor="blue"), column=0, position="branch-top")

        if node.is_leaf():
            node.add_face(ete3.TextFace(node.name, fstyle="italic"), column=0, position="branch-right")
        else:
            if show_paths:
                if show_only_undated_paths:
                    if node.date is None:
                        node.add_face(ete3.TextFace(oldest_path_to_str(node.oldest_path)), column=0, position="branch-bottom")
                else:
                    node.add_face(ete3.TextFace(oldest_path_to_str(node.oldest_path)), column=0, position="branch-bottom")
                    node.add_face(ete3.TextFace(oldest_path_to_str(node.oldest_path_dq)), column=0, position="branch-bottom")

    ts = ete3.TreeStyle()
    ts.show_leaf_name = False
    ts.mode ="r"
    ts.scale = 50
    ts.branch_vertical_margin = 6
    ts.show_scale = False
    ts.margin_right = 100

    tre.render(filename, tree_style=ts)


def plot_dates_algo_guo(input_tre, filename, show_paths=True, show_only_undated_paths=True, pinkblue=True):
    tre = input_tre.copy()

    # look for inconsistent dating, ie where a node has a later date than one of its ancestors
    # mrad = most recent ancestor date of the current node - if the current node is longer ago than this, that's bad
    def check_date_consistency(node, mrad):
        nstyle = ete3.NodeStyle()
        nstyle["vt_line_width"] = 1
        nstyle["hz_line_width"] = 1
        nstyle["vt_line_color"] = "gray"
        nstyle["hz_line_color"] = "gray"

        next_mrad = mrad
        # if we have a date, check it
        if node.is_leaf():
            nstyle["fgcolor"] = "#CFCFCF"
            nstyle["size"] = 10
            node.set_style(nstyle)
        else:
            if node.date and node.date > 0:
                if node.date > mrad:
                    nstyle["fgcolor"] = "magenta"
                    nstyle["size"] = 25
                    node.set_style(nstyle)
                else:
                    if node.imputed_date:
                        nstyle["size"] = 17
                        nstyle["shape"] = "square"
                    else:
                        nstyle["size"] = 18
                    if pinkblue:
                        nstyle["fgcolor"] = "#F1AEE8"
                    else:
                        nstyle["fgcolor"] = "#40B0A6"
                    node.set_style(nstyle)
                    next_mrad = node.date
            else:
                if pinkblue:
                    nstyle["fgcolor"] = "#768AE0"
                else:
                    nstyle["fgcolor"] = "#E1BE6A"
                nstyle["size"] = 17
                nstyle["shape"] = "square"
                node.set_style(nstyle)

        for child in node.children:
            check_date_consistency(child, next_mrad)

    check_date_consistency(tre, 5000)

    sources_dict = {}
    color_index = 0
    for node in tre.traverse(strategy='preorder'):
        if node.date is not None and node.date > 0:
            # median / mean on top of branch
            #node.add_face(ete3.TextFace("Med %.1f / Mean %.1f" % (node.date[1], node.date[0])), column=0, position="branch-top")
            node.add_face(ete3.TextFace("%.2f  " % (node.date), fgcolor="firebrick"), column=0, position="branch-top")

        if node.is_leaf():
            node.add_face(ete3.TextFace(node.name, fstyle="italic"), column=0, position="branch-right")
        else:
            if show_paths:
                if show_only_undated_paths:
                    if node.date is None:
                        node.add_face(ete3.TextFace(node.oldest_paths), column=0, position="branch-bottom")
                else:
                    if node.date is None:
                        node.add_face(ete3.TextFace(node.oldest_paths), column=0, position="branch-bottom")

    ts = ete3.TreeStyle()
    ts.show_leaf_name = False
    ts.mode ="r"
    ts.scale = 50
    ts.branch_vertical_margin = 6
    ts.show_scale = False
    ts.margin_right = 100

    tre.render(filename, tree_style=ts)


def plot_dates_algo_ed_guo(input_tre, filename, show_paths=True, show_only_undated_paths=True, pinkblue=True):
    tre = input_tre.copy()

    # look for inconsistent dating, ie where a node has a later date than one of its ancestors
    # mrad = most recent ancestor date of the current node - if the current node is longer ago than this, that's bad
    def check_date_consistency(node, mrad):
        nstyle = ete3.NodeStyle()
        nstyle["vt_line_width"] = 1
        nstyle["hz_line_width"] = 1
        nstyle["vt_line_color"] = "gray"
        nstyle["hz_line_color"] = "gray"

        next_mrad = mrad
        # if we have a date, check it
        if node.is_leaf():
            nstyle["fgcolor"] = "#CFCFCF"
            nstyle["size"] = 10
            node.set_style(nstyle)
        else:
            if node.date and node.date > 0:
                if node.date > mrad:
                    nstyle["fgcolor"] = "magenta"
                    nstyle["size"] = 25
                    node.set_style(nstyle)
                else:
                    if node.imputed_date:
                        nstyle["size"] = 17
                        nstyle["shape"] = "square"
                    else:
                        nstyle["size"] = 18
                    if pinkblue:
                        nstyle["fgcolor"] = "#F1AEE8"
                    else:
                        nstyle["fgcolor"] = "#40B0A6"
                    node.set_style(nstyle)
                    next_mrad = node.date
            else:
                if pinkblue:
                    nstyle["fgcolor"] = "#768AE0"
                else:
                    nstyle["fgcolor"] = "#E1BE6A"
                nstyle["size"] = 17
                nstyle["shape"] = "square"
                node.set_style(nstyle)

        for child in node.children:
            check_date_consistency(child, next_mrad)

    check_date_consistency(tre, 5000)

    sources_dict = {}
    color_index = 0
    for node in tre.traverse(strategy='preorder'):
        if node.date is not None and node.date > 0:
            # median / mean on top of branch
            #node.add_face(ete3.TextFace("Med %.1f / Mean %.1f" % (node.date[1], node.date[0])), column=0, position="branch-top")
            node.add_face(ete3.TextFace("%.2f  " % (node.date), fgcolor="firebrick"), column=0, position="branch-top")

        if node.is_leaf():
            node.add_face(ete3.TextFace(node.name, fstyle="italic"), column=0, position="branch-right")

        if show_paths:
            if show_only_undated_paths:
                if node.date is None:
                    node.add_face(ete3.TextFace(node.ed_scores_guo), column=0, position="branch-top")
                    node.add_face(ete3.TextFace(node.oldest_paths), column=0, position="branch-bottom")
            else:
                node.add_face(ete3.TextFace(node.ed_scores_guo), column=0, position="branch-top")
                node.add_face(ete3.TextFace(node.oldest_paths), column=0, position="branch-bottom")

    ts = ete3.TreeStyle()
    ts.show_leaf_name = False
    ts.mode ="r"
    ts.scale = 50
    ts.branch_vertical_margin = 6
    ts.show_scale = False
    ts.margin_right = 100

    tre.render(filename, tree_style=ts)



def plot_ed_scores(input_tre, filename):
    tre = input_tre.copy()

    for node in tre.traverse(strategy='preorder'):
        if node is not tre:
            node.add_face(ete3.TextFace("%.1f  " % node.dist, fsize=8, fgcolor="black"), column=0, position="branch-top")
            #if node.is_leaf():
            node.add_face(ete3.TextFace("%.2f  /  %d  " % (node.ed_score, node.desc_leaves), fsize=8, fgcolor="black"), column=0, position="branch-bottom")
            #else:
            #    node.add_face(ete3.TextFace("%d  " % node.desc_leaves, fsize=8, fgcolor="black"), column=0, position="branch-bottom")

    ts = ete3.TreeStyle()
    ts.show_leaf_name = True
    ts.branch_vertical_margin = 5
    ts.mode ="r"
    ts.scale = 80

    tre.render(filename, tree_style=ts)


def plot_bd(input_tre, filename, show_ranks=False, show_nodes=True):
    tre = input_tre.copy()

    # look for inconsistent dating, ie where a node has a later date than one of its ancestors
    # mrad = most recent ancestor date of the current node - if the current node is longer ago than this, that's bad
    def check_date_consistency(node, mrad):
        nstyle = ete3.NodeStyle()
        nstyle["vt_line_width"] = 1
        nstyle["hz_line_width"] = 1
        nstyle["vt_line_color"] = "gray"
        nstyle["hz_line_color"] = "gray"

        next_mrad = mrad

        # if we have a date, check it
        if node.date:
            if node.date > mrad:
                nstyle["fgcolor"] = "magenta"
                if show_nodes:
                    nstyle["size"] = 30
                else:
                    nstyle["size"] = 0
                node.set_style(nstyle)
            else:
                if node.extinct:
                    nstyle["fgcolor"] = "blue"
                    if show_nodes:
                        nstyle["size"] = 25
                    else:
                        nstyle["size"] = 0
                    node.set_style(nstyle)
                elif len(node.children) == 1:
                    nstyle["fgcolor"] = "cyan"
                    if show_nodes:
                        nstyle["size"] = 25
                    else:
                        nstyle["size"] = 0
                    node.set_style(nstyle)
                else:
                    nstyle["fgcolor"] = "red"
                    if show_nodes:
                        nstyle["size"] = 20
                    else:
                        nstyle["size"] = 0
                    node.set_style(nstyle)
                next_mrad = node.date
        else:
            if node.is_leaf():
                nstyle["fgcolor"] = "grey"
                nstyle["size"] = 12
                node.set_style(nstyle)
            else:
                nstyle["fgcolor"] = "blue"
                nstyle["size"] = 15
                node.set_style(nstyle)

        for child in node.children:
            check_date_consistency(child, next_mrad)

    check_date_consistency(tre, 5000)

    colors = ["#e6194b",
              "#f58231",
              "#3cb44b",
              "#4363d8",
              "#911eb4",
              "#f032e6",
              "#bfef45",
              "#42d4f4",
              "#800000",
              "#a9a9a9",
              "#808000",
              "#469990",
              "#000075",
              "#000000",
              "#9a6324",
              "#fabed4",
              "#aaffc3",
              "#dcbeff",
              "#ffd8b1",
              "#ffe119",
              "#ffd8b1",
              "#fffac8"]

    sources_dict = {}
    color_index = 0
    for node in tre.traverse(strategy='preorder'):
        if node.date is not None:
            # median / mean on top of branch
            #node.add_face(ete3.TextFace("Med %.1f / Mean %.1f" % (node.date[1], node.date[0])), column=0, position="branch-top")
            if node.imputed_date:
                node.add_face(ete3.TextFace("Imp. %.2f" % (node.date)), column=0, position="branch-top")
            else:
                node.add_face(ete3.TextFace("%.2f " % (node.date)), column=0, position="branch-top")
                node.add_face(ete3.TextFace("%.2f " % (node.dist)), column=0, position="branch-bottom")

            if show_ranks:
                node.add_face(ete3.TextFace("%s, %s" % (node.name, node.tx_level), fgcolor="black"), column=0, position="branch-bottom")

            # list of dates on bottom of branch, coloured by the phylogeny they came from
            #for i in range(len(node.date[2])):
            #    date = node.date[2][i]
                # source = node.date[3][i][1]
                # if source not in sources_dict:
                #     sources_dict[source] = colors[color_index]
                #     color_index += 1
                #     if color_index > 0 and color_index % len(colors) == 0:
                #         color_index = 0
                #node.add_face(ete3.TextFace("%.1f    " % date, fgcolor=sources_dict[source]), column=i, position="branch-bottom")
            #    node.add_face(ete3.TextFace("%.1f    " % date, fgcolor="black"), column=i, position="branch-bottom")
        else:
            node.add_face(ete3.TextFace(" ", fgcolor="black"), column=0, position="branch-top")
            if node.is_leaf():
                node.add_face(ete3.TextFace("tip", fgcolor="black"), column=0, position="branch-bottom")
            else:
                node.add_face(ete3.TextFace(" ", fgcolor="black"), column=0, position="branch-bottom")

    ts = ete3.TreeStyle()
    ts.show_leaf_name = False
    ts.mode ="r"
    ts.scale = 20
    ts.rotation = 90
    ts.complete_branch_lines_when_necessary = False

    tre.render(filename, tree_style=ts)


def plot_ultrametric(input_tre, filename, color_long=False, color_short=False):
    tre = input_tre.copy()

    # plot tree
    for node in tre.traverse(strategy="preorder"):

        nstyle = ete3.NodeStyle()
        nstyle["vt_line_width"] = 2
        nstyle["hz_line_width"] = 2
        nstyle["vt_line_color"] = "grey"
        nstyle["hz_line_color"] = "grey"

        nstyle["size"] = 0

        node.set_style(nstyle)

        # if node.is_leaf():
        #     node.add_face(ete3.TextFace(name_to_simple_name(node.name), fstyle="italic"), column=0, position="branch-right")
        #     node.add_face(ete3.TextFace(""), column=0, position="branch-bottom")

    if color_long:
        nstyle = ete3.NodeStyle()
        nstyle["vt_line_width"] = 2
        nstyle["vt_line_color"] = "grey"
        nstyle["hz_line_width"] = 2
        nstyle["hz_line_color"] = "crimson"
        nstyle["hz_line_type"] = 1
        nstyle["size"] = 0

        tre.children[0].set_style(nstyle)
        tre.children[0].children[1].set_style(nstyle)
        tre.children[0].children[1].children[1].set_style(nstyle)

        tre.children[1].set_style(nstyle)
        tre.children[1].children[0].set_style(nstyle)
        tre.children[1].children[0].children[0].set_style(nstyle)

    elif color_short:
        nstyle = ete3.NodeStyle()
        nstyle["vt_line_width"] = 2
        nstyle["vt_line_color"] = "grey"
        nstyle["hz_line_width"] = 2
        nstyle["hz_line_color"] = "crimson"
        nstyle["hz_line_type"] = 1
        nstyle["size"] = 0

        tre.children[0].set_style(nstyle)
        tre.children[0].children[0].set_style(nstyle)

        tre.children[1].set_style(nstyle)
        tre.children[1].children[1].set_style(nstyle)


    tre.dist = 0.01
    nstyle = ete3.NodeStyle()
    nstyle["vt_line_width"] = 2
    nstyle["hz_line_width"] = 2
    nstyle["vt_line_color"] = "grey"
    nstyle["hz_line_color"] = "grey"
    nstyle["size"] = 10
    tre.set_style(nstyle)

    ts = ete3.TreeStyle()
    ts.margin_right=100
    ts.show_leaf_name = False
    ts.mode ="r"
    ts.rotation = 90

    ts.scale = 20
    ts.show_scale = False
    # ts.show_branch_length = True
    ts.min_leaf_separation = 15

    tre.render(filename, tree_style=ts)


def plot_ultrametric_interp(input_tre, filename, color_long=False, color_short=False):
    tre = input_tre.copy()

    # plot tree
    for node in tre.traverse(strategy="preorder"):

        nstyle = ete3.NodeStyle()
        nstyle["vt_line_width"] = 2
        nstyle["hz_line_width"] = 2
        nstyle["vt_line_color"] = "grey"
        nstyle["hz_line_color"] = "grey"

        nstyle["size"] = 10

        node.set_style(nstyle)

        if node.name == "A":
            node.add_face(ete3.TextFace("A ", fgcolor="crimson"), column=0, position="branch-top")
            # node.add_face(ete3.TextFace("test"), column=0, position="branch-bottom")
        elif node.name == "root":
            node.add_face(ete3.TextFace("12 "), column=0, position="branch-top")
        elif node.name == "int1":
            node.add_face(ete3.TextFace("7"), column=0, position="branch-top")
        elif node.name == "int2":
            node.add_face(ete3.TextFace("11 "), column=0, position="branch-top")


    tre.dist = 0.01
    nstyle = ete3.NodeStyle()
    nstyle["vt_line_width"] = 2
    nstyle["hz_line_width"] = 2
    nstyle["vt_line_color"] = "grey"
    nstyle["hz_line_color"] = "grey"
    nstyle["size"] = 10
    tre.set_style(nstyle)

    ts = ete3.TreeStyle()
    ts.margin_right=100
    ts.show_leaf_name = False
    ts.mode ="r"
    # ts.rotation = 90

    ts.scale = 20
    ts.show_scale = False
    # ts.show_branch_length = True
    ts.min_leaf_separation = 20

    tre.render(filename, tree_style=ts)
