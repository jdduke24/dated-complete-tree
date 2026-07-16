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



import ete4
import numpy as np
from .taxonomy_utils import tx_levels
from . import tree_fixing


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


def plot_figure_polytomies(input_tre, filename, name_mrcas=True, info_colors=True, arrows=True, scale=5):
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
                from ete4.treeview import TextFace
                from PyQt6.QtGui import QFont
                # Add the missing attribute that ETE is looking for
                if not hasattr(QFont, "StyleItalic"):
                    QFont.StyleItalic = QFont.Style.StyleItalic
                tf = TextFace(" " + name_to_simple_name(node.name), fstyle="italic")
                node.add_face(tf, column=0, position="branch-right")
            else:
                node.add_face(ete4.treeview.TextFace((' ' * spacing) +
                                            "%s\n" % name_to_simple_name(node.name) +
                                            (' ' * spacing) +
                                            "%s" % node.props["tx_level"], fsize=9), column=0, position="branch-bottom")


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


def plot_dates(input_tre, filename, show_ranks=False, show_dates=True, show_nodes=True, title=None):
    tre = input_tre.copy()

    # look for inconsistent dating, ie where a node has a later date than one of its ancestors
    # mrad = most recent ancestor date of the current node - if the current node is longer ago than this, that's bad
    def check_date_consistency(node, mrad):
        nstyle = ete4.treeview.NodeStyle()
        nstyle["vt_line_width"] = 1
        nstyle["hz_line_width"] = 1
        nstyle["vt_line_color"] = "gray"
        nstyle["hz_line_color"] = "gray"

        next_mrad = mrad
        # if we have a date, check it
        if node.props["date"]:
            if node.props["date"] > mrad:
                nstyle["fgcolor"] = "magenta"
                nstyle["size"] = 30
                node.set_style(nstyle)
            else:
                if node.props["imputed_date"]:
                    nstyle["fgcolor"] = "orange"
                    nstyle["size"] = 25
                    node.set_style(nstyle)
                else:
                    nstyle["fgcolor"] = "red"
                    nstyle["size"] = 20
                    node.set_style(nstyle)
                next_mrad = node.props["date"]
        else:
            if node.is_leaf:
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
            if node.props["imputed_date"]:
                if node.props["imputation_type"] == 5:
                    nstyle["hz_line_color"] = "green"
                    nstyle["vt_line_width"] = 1
                    nstyle["hz_line_width"] = 2
                elif node.props["imputation_type"] > 0:
                    nstyle["hz_line_color"] = "crimson"
                    nstyle["vt_line_width"] = 1
                    nstyle["hz_line_width"] = 2

                node.set_style(nstyle)

            else:
                if node.up and node.up.props["imputed_date"]:
                    if node.up.props["imputation_type"] == 5:
                        nstyle["hz_line_color"] = "green"
                        nstyle["vt_line_width"] = 1
                        nstyle["hz_line_width"] = 2
                    elif node.up.props["imputation_type"] > 0:
                        nstyle["hz_line_color"] = "crimson"
                        nstyle["vt_line_width"] = 1
                        nstyle["hz_line_width"] = 2

                    node.set_style(nstyle)

        for child in node.children:
            check_date_consistency(child, next_mrad)

    check_date_consistency(tre, 5000)

    if show_dates:
        for node in tre.traverse(strategy='preorder'):
            if node.props["date"] is not None:
                if node.props["imputed_date"]:
                    node.add_face(ete4.treeview.TextFace("Imp. %.2f" % (node.props["date"])), column=0, position="branch-top")
                else:
                    if node.props["date"] == 0:
                        # node.add_face(ete4.treeview.TextFace("tip"), column=0, position="branch-top")
                        node.add_face(ete4.treeview.TextFace(node.name), column=0, position="branch-top")
                    else:
                        node.add_face(ete4.treeview.TextFace("%.2f " % (node.props["date"])), column=0, position="branch-top")

                if show_ranks:
                    node.add_face(ete4.treeview.TextFace("%s, %s" % (node.name, node.props["tx_level"]), fgcolor="black"), column=0, position="branch-bottom")
            else:
                node.add_face(ete4.treeview.TextFace(" ", fgcolor="black"), column=0, position="branch-top")
                if node.is_leaf:
                    # node.add_face(ete4.treeview.TextFace("tip", fgcolor="black"), column=0, position="branch-bottom")
                    node.add_face(ete4.treeview.TextFace(node.name, fgcolor="black"), column=0, position="branch-bottom")
                else:
                    node.add_face(ete4.treeview.TextFace(" ", fgcolor="black"), column=0, position="branch-bottom")

    ts = ete4.treeview.TreeStyle()
    ts.show_leaf_name = False
    ts.mode ="r"
    ts.scale = 30

    if not show_dates:
        ts.min_leaf_separation = 10
    if not show_nodes:
        ts.complete_branch_lines_when_necessary = False

    if title:
        ts.title.add_face(ete4.treeview.faces.TextFace(title, fsize=16), 0)

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


def plot_dates_algo(input_tre, filename, show_paths=True, show_only_undated_paths=True, mu=False, pinkblue=True, vs=3):
    """For Figure 5 of the paper."""

    tre = input_tre.copy()

    # look for inconsistent dating, ie where a node has a later date than one of its ancestors
    # mrad = most recent ancestor date of the current node - if the current node is longer ago than this, that's bad
    def check_date_consistency(node, mrad):
        nstyle = ete4.treeview.NodeStyle()
        nstyle["vt_line_width"] = 1
        nstyle["hz_line_width"] = 1
        nstyle["vt_line_color"] = "gray"
        nstyle["hz_line_color"] = "gray"

        next_mrad = mrad
        # if we have a date, check it
        if node.is_leaf:
            nstyle["fgcolor"] = "#CFCFCF"
            nstyle["size"] = 10
            node.set_style(nstyle)
        else:
            if node.props["date"] and node.props["date"] > 0:
                if node.props["date"] > mrad:
                    nstyle["fgcolor"] = "magenta"
                    nstyle["size"] = 25
                    node.set_style(nstyle)
                else:
                    if node.props["imputed_date"]:
                        nstyle["size"] = 17
                        nstyle["shape"] = "square"
                    else:
                        nstyle["size"] = 18
                    if pinkblue:
                        nstyle["fgcolor"] = "#F1AEE8"
                    else:
                        nstyle["fgcolor"] = "#40B0A6"
                    node.set_style(nstyle)
                    next_mrad = node.props["date"]
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

    for node in tre.traverse(strategy='preorder'):
        if node.props["date"] is not None and node.props["date"] > 0:
            if show_paths:
                node.add_face(ete4.treeview.TextFace("  %.2f " % (node.props["date"]), fgcolor="firebrick"), column=0, position="branch-top")
            else:
                node.add_face(ete4.treeview.TextFace("%.2f " % (node.props["date"]), fgcolor="firebrick"), column=0, position="branch-top")

        if node.is_leaf:
            from PyQt6.QtGui import QFont
            # Add the missing attribute that ETE is looking for
            if not hasattr(QFont, "StyleItalic"):
                QFont.StyleItalic = QFont.Style.StyleItalic

            node.add_face(ete4.treeview.TextFace(" "+node.name, fstyle="italic"), column=0, position="branch-right")

        if show_paths:
            if show_only_undated_paths:
                if node.props["date"] is None:
                    node.add_face(ete4.treeview.TextFace(oldest_path_to_str(node.props["oldest_path_long"]) + " "), column=0, position="branch-bottom")
            else:
                node.add_face(ete4.treeview.TextFace(oldest_path_to_str(node.props["oldest_path_long"]) + " ", fsize=8), column=0, position="branch-bottom")

    ts = ete4.treeview.TreeStyle()
    ts.show_leaf_name = False
    ts.mode ="r"
    ts.scale = 100
    ts.branch_vertical_margin = vs
    ts.show_scale = False
    ts.margin_right = 10
    ts.complete_branch_lines_when_necessary = False

    tre.render(filename, tree_style=ts, units="in", w=6, dpi=300)


def plot_undated_tree(input_tre, filename=None, vs=3, w=800, to_highlight=[], hl_polytomies=False, inline=False):
    """Draw a cladogram. Nodes named in to_highlight will be labelled in bold. Nodes named in hl_polytomies will be
    highlighted red. Nodes with props["ph_tx"] == 'IN' will be highlighted red. Use inline=True for display in IPython."""
    tre = input_tre.copy()

    # strip OTT stuff from names
    for node in tre.traverse():
        name_parts = node.name.split('_')
        if len(name_parts) == 2:
            node.name = name_parts[0]
        elif len(name_parts) == 3:
            node.name = name_parts[0] + " " + name_parts[1]

    # strip out OTT stuff from names of nodes to highlight
    hl = []
    for name in to_highlight:
        name_parts = name.split('_')
        if len(name_parts) == 2:
            hl.append(name_parts[0])
        elif len(name_parts) == 3:
            hl.append(name_parts[0] + " " + name_parts[1])

    # to draw the tree nicely we want to space out the branch lengths - we can use the date interpolation for this
    for node in tre.traverse():
        if node.is_leaf:
            node.add_prop("date", 0)
        else:
            node.add_prop("date", None)
        node.add_prop("imputed_date", False)
        node.add_prop("imputation_type", 0)

    tre.props["date"] = 4

    from dated_complete_tree import tree_dating
    tree_dating.date_labelling(tre)
    tree_dating.impute_missing_dates(tre)

    tree_dating.compute_branch_lengths(tre)

    tre.dist = 1

    scale = 100

    # do some rescaling of branch lengths because we want to draw nodes, and nodes are
    # added on top of the scaled branch lengths by ete. This means the leaf nodes don't
    # line up neatly on the right - we fix that here.
    for node in tre.traverse():
        if node.is_leaf:
            node.dist -= 10/scale
        else:
            node.dist -= 12/scale

        if node.up and len(node.up.children) > 1:
            node.dist -= 1/scale
        node.props["adj"] = True

    # look for inconsistent dating, ie where a node has a later date than one of its ancestors
    # mrad = most recent ancestor date of the current node - if the current node is longer ago than this, that's bad
    def add_node_styles(node):
        nstyle = ete4.treeview.NodeStyle()
        nstyle["vt_line_width"] = 1
        nstyle["hz_line_width"] = 1
        nstyle["vt_line_color"] = "gray"
        nstyle["hz_line_color"] = "gray"

        if node.is_leaf:
            nstyle["fgcolor"] = "#AAAAAA"
            nstyle["size"] = 10
            if node.props["ph_tx"] == "TX":
                nstyle["shape"] = "square"
            node.set_style(nstyle)
        elif hl_polytomies and node.name == "mrcapoly":
            nstyle["size"] = 12
            nstyle["shape"] = "circle"
            nstyle["fgcolor"] = "#D81B60"
            node.set_style(nstyle)
        else:
            nstyle["size"] = 12
            if node.props["ph_tx"] == "TX":
                nstyle["shape"] = "square"
                nstyle["fgcolor"] = "#40B0A6"
                node.set_style(nstyle)
            elif node.props["ph_tx"] == "PH":
                nstyle["shape"] = "circle"
                nstyle["fgcolor"] = "#40B0A6"
                node.set_style(nstyle)
            elif node.props["ph_tx"] == "IN_old":
                nstyle["shape"] = "circle"
                nstyle["fgcolor"] = "#E1BE6A"
                node.set_style(nstyle)
            else:
                nstyle["shape"] = "circle"
                nstyle["fgcolor"] = "#D81B60"
                node.set_style(nstyle)

        for child in node.children:
            add_node_styles(child)

    add_node_styles(tre)

    for node in tre.traverse(strategy='preorder'):
        if not node.is_leaf and node.name[:4] != "mrca":
            node.add_face(ete4.treeview.TextFace(node.name, bold=(node.name in hl), fsize=9), column=0, position="float")
            node.add_face(ete4.treeview.TextFace(node.props["tx_level"], bold=(node.name in hl), fsize=9), column=0, position="float")

        if node.is_leaf:
            node.add_face(ete4.treeview.TextFace(" "+node.name, bold=(node.name in hl), fsize=9), column=0, position="branch-right")

    ts = ete4.treeview.TreeStyle()
    ts.show_leaf_name = False
    ts.mode ="r"
    ts.scale = scale
    ts.branch_vertical_margin = vs
    ts.show_scale = False
    ts.margin_right = 10
    ts.complete_branch_lines_when_necessary = False
    ts.allow_face_overlap = True

    if inline:
        from IPython.display import display
        display(tre.render(file_name="%%inline", tree_style=ts, w=w))
    else:
        tre.render(filename, tree_style=ts, units="in", w=8, dpi=300)

    # mark the newly inserted nodes as having been drawn once already
    for node in input_tre.traverse():
        if node.props["ph_tx"] == "IN":
            node.props["ph_tx"] = "IN_old"


def plot_dated_tree(input_tre, filename=None, vs=2, leaves_to_highlight=[], true_branch_lengths=False, inline=False, w=800):
    """Plot a partially- or fully-dated tree. If true_branch_lengths==True then branch lengths will be scaled to time. Otherwise
    a cladogram is drawn with any dates present labelled. Nodes named in to_highlight will be labelled in bold.
    Nodes named in hl_polytomies will be highlighted red. Nodes with props["ph_tx"] == 'IN' will be highlighted red.
    Use inline=True for display in IPython."""

    tre = input_tre.copy()

    # strip OTT stuff from names
    for node in tre.traverse():
        name_parts = node.name.split('_')
        if len(name_parts) == 2:
            node.name = name_parts[0]
        elif len(name_parts) == 3:
            node.name = name_parts[0] + " " + name_parts[1]

    # strip out OTT stuff from names of nodes to highlight
    hl = []
    for name in leaves_to_highlight:
        name_parts = name.split('_')
        if len(name_parts) == 2:
            hl.append(name_parts[0])
        elif len(name_parts) == 3:
            hl.append(name_parts[0] + " " + name_parts[1])


    # look for inconsistent dating, ie where a node has a later date than one of its ancestors
    # mrad = most recent ancestor date of the current node - if the current node is longer ago than this, that's bad
    def add_node_styles(node):
        nstyle = ete4.treeview.NodeStyle()
        nstyle["vt_line_width"] = 1
        nstyle["hz_line_width"] = 1
        nstyle["vt_line_color"] = "gray"
        nstyle["hz_line_color"] = "gray"
        nstyle["size"] = 1
        node.set_style(nstyle)

        for child in node.children:
            add_node_styles(child)

    add_node_styles(tre)

    for node in tre.traverse(strategy='preorder'):
        if node.props["date"] is not None and node.props["date"] > 0:
            node.add_face(ete4.treeview.TextFace("%.2f" % (node.props["date"]), fgcolor="firebrick", fsize=8), column=0, position="float")
            node.add_face(ete4.treeview.TextFace("" % (node.props["date"]), fgcolor="firebrick"), column=0, position="float")

        if node.is_leaf:
            node.add_face(ete4.treeview.TextFace(" "+node.name, fsize=9), column=0, position="branch-right")

    tre.dist = 1

    from dated_complete_tree import tree_dating

    if not true_branch_lengths:
        scale = 100
        for node in tre.traverse():
            if node.is_leaf:
                node.add_prop("date", 0)
            else:
                node.add_prop("date", None)
            node.add_prop("imputed_date", False)
            node.add_prop("imputation_type", 0)

        tre.props["date"] = 4

        tree_dating.date_labelling(tre)
        tree_dating.impute_missing_dates(tre)

        tree_dating.compute_branch_lengths(tre)
    else:
        scale = 100
        scaling_factor = 4./tre.props["date"]
        for node in tre.traverse():
            node.props["date"] *= scaling_factor

    tree_dating.compute_branch_lengths(tre)

    most_nodes = 0
    for leaf in tre.leaves():
        path = 0
        nd = leaf
        while nd.up != tre:
            path += 1
            nd = nd.up
        if path > most_nodes:
            most_nodes = path

    for leaf in tre.leaves():
        path = 0
        nd = leaf
        while nd.up != tre:
            path += 1
            nd = nd.up

        leaf.dist += (most_nodes-path) * 1/scale
        leaf.dist += (most_nodes-path) * 1/scale

    ts = ete4.treeview.TreeStyle()
    ts.show_leaf_name = False
    ts.mode ="r"
    ts.scale = scale
    ts.branch_vertical_margin = vs
    ts.show_scale = False
    ts.margin_right = 10
    ts.complete_branch_lines_when_necessary = False
    ts.allow_face_overlap = True

    if inline:
        from IPython.display import display
        display(tre.render(file_name="%%inline", tree_style=ts, w=w))
    else:
        tre.render(filename, tree_style=ts, units="in", w=8, dpi=300)
